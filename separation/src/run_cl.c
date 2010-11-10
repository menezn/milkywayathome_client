/* Copyright 2010 Matthew Arsenault, Travis Desell, Boleslaw
Szymanski, Heidi Newberg, Carlos Varela, Malik Magdon-Ismail and
Rensselaer Polytechnic Institute.

This file is part of Milkway@Home.

Milkyway@Home is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Milkyway@Home is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "separation_types.h"
#include "milkyway_util.h"
#include "show_cl_types.h"
#include "setup_cl.h"
#include "milkyway_cl.h"
#include "separation_cl_buffers.h"
#include "separation_cl_defs.h"
#include "calculated_constants.h"
#include "run_cl.h"
#include "r_points.h"
#include "integrals_common.h"

typedef struct
{
    cl_event outMap;  /* reading main output buffer */
    cl_event probMap; /* reading stream probs buffer */

    cl_event outUnmap;
    cl_event probUnmap;

    cl_event endTmp;   /* end of the NDRange writing to the temporary output buffers */
} SeparationCLEvents;

static void printSeparationEventTimes(SeparationCLEvents* evs)
{
    printf("NDrange:    %.15g s\n"
           "Read out:   %.15g s\n"
           "Read probs: %.15g s\n",
           mwEventTime(evs->endTmp),
           mwEventTime(evs->outMap),
           mwEventTime(evs->probMap));
}

static real* mapIntegralResults(CLInfo* ci,
                                SeparationCLMem* cm,
                                SeparationCLEvents* evs,
                                size_t resultsSize)
{
    cl_int err;
    real* mapOutMu;

    mapOutMu = (real*) clEnqueueMapBuffer(ci->queue,
                                          cm->outMu,
                                          CL_TRUE, CL_MAP_READ,
                                          0, resultsSize,
                                          0, NULL,
                                          NULL,
                                          &err);
    if (err != CL_SUCCESS)
        warn("Error mapping integral result buffer: %s\n", showCLInt(err));

    return mapOutMu;
}

static real* mapProbsResults(CLInfo* ci,
                             SeparationCLMem* cm,
                             SeparationCLEvents* evs,
                             size_t probsResultsSize)
{
    cl_int err;
    real* mapOutProbs;

    mapOutProbs = (real*) clEnqueueMapBuffer(ci->queue,
                                             cm->outProbs,
                                             CL_TRUE, CL_MAP_READ,
                                             0, probsResultsSize,
                                             0, NULL,
                                             NULL,
                                             &err);
    if (err != CL_SUCCESS)
        warn("Error mapping probs result buffer: %s\n", showCLInt(err));

    return mapOutProbs;
}

static void sumStreamResults(KAHAN* probs_results,
                             const real* probs_V_reff_xr_rp3,
                             const cl_uint number_streams)
{
    cl_uint i;

    for (i = 0; i < number_streams; ++i)
        KAHAN_ADD(probs_results[i], probs_V_reff_xr_rp3[i]);
}

static void sumProbsResults(real* probs_results,
                            const real* st_probs_V_reff_xr_rp3_mu_r,
                            const cl_uint mu_steps,
                            const cl_uint r_steps,
                            const cl_uint number_streams)
{
    cl_uint i, j, idx;
    KAHAN* probs_sum;

    probs_sum = (KAHAN*) callocSafe(number_streams, sizeof(KAHAN));

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            idx = (i * r_steps * number_streams) + (j * number_streams);
            sumStreamResults(probs_sum, &st_probs_V_reff_xr_rp3_mu_r[idx], number_streams);
        }
    }

    for (i = 0; i < number_streams; ++i)
        probs_results[i] = probs_sum[i].sum + probs_sum[i].correction;

    free(probs_sum);
}

/* Find the smallest integer i >= x that is divisible by n */
static size_t nextMultiple(const size_t x, const size_t n)
{
    return (x % n == 0) ? x : x + (n - x % n);
}

static cl_int cpuGroupSizes(const INTEGRAL_AREA* ia, size_t global[], size_t local[])
{
    global[0] = 1;
    global[1] = ia->mu_steps;
    global[2] = ia->r_steps;
    local[0] = 1;
    local[1] = 1;
    local[2] = 1;

    return CL_SUCCESS;
}

static cl_int gpuGroupSizes(const INTEGRAL_AREA* ia, size_t numChunks, size_t global[], size_t local[])
{
    size_t groupSize = 64;
    size_t chunkSize = ia->r_steps / numChunks;

    if (ia->r_steps % numChunks != 0)
    {
        warn("Error: r_steps (%zu) not divisible by number of chunks (%zu)\n", ia->r_steps, numChunks);
        return -1;
    }

    /* Ideally these are already nicely divisible by 32/64/128, otherwise
     * round up a bit. */
    global[0] = 1;
    global[1] = ia->mu_steps;
    global[2] = chunkSize;

    /* Bias towards mu steps seems to be better */
    local[0] = 1;
    local[1] = groupSize;
    local[2] = 1;

    return CL_SUCCESS;
}

/* TODO: Do we need hadware specific workgroup sizes? */
static cl_bool findWorkGroupSizes(CLInfo* ci,
                                  const INTEGRAL_AREA* ia,
                                  size_t numChunks,
                                  size_t global[],
                                  size_t local[])
{
    cl_int err;
    WGInfo wgi;

    err = getWorkGroupInfo(ci, &wgi);
    if (err != CL_SUCCESS)
        warn("Failed to get work group info: %s\n", showCLInt(err));
    else
        printWorkGroupInfo(&wgi);

    if (ci->devType == CL_DEVICE_TYPE_CPU)
        err = cpuGroupSizes(ia, global, local);
    else
        err = gpuGroupSizes(ia, numChunks, global, local);

    if (err != CL_SUCCESS)
    {
        warn("Failed to choose group size\n");
        return CL_TRUE;
    }

    size_t localSize = local[0] * local[1] * local[2];

    warn("Range is { nu_steps = %zu, mu_steps = %zu, r_steps = %zu }\n"
         "Chunked range is  { nu_steps = %zu, mu_steps = %zu, r_steps = %zu }\n"
         "Attempting to use a workgroup size of { %zu, %zu, %zu } = %zu \n",
         ia->nu_steps, ia->mu_steps, ia->r_steps,
         global[0], global[1], global[2],
         local[0], local[1], local[2],
         localSize);

    /* Sanity checks */
    /* TODO: also check against CL_DEVICE_MAX_WORK_ITEM_SIZES */
    return (global[0] % local[0] || global[1] % local[1] || global[2] % local[2] || localSize > wgi.wgs);
}

static cl_int enqueueIntegralKernel(CLInfo* ci,
                                    SeparationCLEvents* evs,
                                    const size_t offset[],
                                    const size_t global[],
                                    const size_t local[])
{
    cl_int err;

    err = clEnqueueNDRangeKernel(ci->queue,
                                 ci->kern,
                                 3,
                                 offset, global, local,
                                 0, NULL, &evs->endTmp);
    if (err != CL_SUCCESS)
    {
        warn("Error enqueueing integral kernel execution: %s\n", showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

/* The mu and r steps for a nu step were done in parallel. After that
 * we need to add the result to the running total + correction from
 * all the nu steps */
static real sumMuResults(const real* mu_results,
                         const cl_uint mu_steps,
                         const cl_uint r_steps)
{
    cl_uint i, j;
    KAHAN bg_prob = ZERO_KAHAN;

    for (i = 0; i < mu_steps; ++i)
    {
        for (j = 0; j < r_steps; ++j)
        {
            KAHAN_ADD(bg_prob, mu_results[r_steps * i + j]);
        }
    }

    return bg_prob.sum + bg_prob.correction;
}

static cl_int setNuKernelArgs(CLInfo* ci, const INTEGRAL_AREA* ia, const cl_uint nu_step)
{
    cl_int err;
    NU_ID nuid;

    /* Avoid doing any trig in the broken ATI math. Also trig seems to
     * be more expensive there. Not doing the coordinate conversion
     * there also halves number of required registers, which prevents
     * enough threads to hide the horrible latency of the other
     * required reads. */
    nuid = calc_nu_step(ia, nu_step);
    err = clSetKernelArg(ci->kern, 9, sizeof(real), &nuid.id);
    if (err != CL_SUCCESS)
    {
        warn("Error setting nu_id argument for step %u: %s\n", nu_step, showCLInt(err));
        return err;
    }

    return CL_SUCCESS;
}

static inline real readKernelResults(CLInfo* ci,
                                     SeparationCLMem* cm,
                                     SeparationCLEvents* evs,
                                     real* probs_results,
                                     const cl_uint mu_steps,
                                     const cl_uint r_steps,
                                     const cl_uint number_streams)
{
    cl_int err;
    real* mu_results;
    real* probs_tmp;
    real result;

    size_t resultSize = sizeof(real) * mu_steps * r_steps;
    mu_results = mapIntegralResults(ci, cm, evs, resultSize);
    if (!mu_results)
    {
        warn("Failed to map integral results\n");
        return NAN;
    }

    result = sumMuResults(mu_results, mu_steps, r_steps);

    err = clEnqueueUnmapMemObject(ci->queue, cm->outMu, mu_results, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        warn("Failed to unmap results buffer: %s\n", showCLInt(err));
        return NAN;
    }

    size_t probsSize = sizeof(real) * mu_steps * r_steps * number_streams;
    probs_tmp = mapProbsResults(ci, cm, evs, probsSize);
    if (!probs_tmp)
    {
        warn("Failed to map probs results\n");
        return NAN;
    }

    sumProbsResults(probs_results, probs_tmp, mu_steps, r_steps, number_streams);

    err = clEnqueueUnmapMemObject(ci->queue, cm->outProbs, probs_tmp, 0, NULL, &evs->probUnmap);
    if (err != CL_SUCCESS)
    {
        warn("Failed to unmap probs buffer: %s\n", showCLInt(err));
        return NAN;
    }

    return result;
}

static cl_int runNuStep(CLInfo* ci,
                        SeparationCLMem* cm,
                        SeparationCLEvents* evs,

                        const INTEGRAL_AREA* ia,
                        const size_t numChunks,
                        const size_t global[],
                        const size_t local[],
                        const cl_uint number_streams,
                        const cl_uint nu_step)
{
    cl_int err;
    size_t offset[2];
    unsigned int i;

    printf("Nu step %u:\n", nu_step);
    //printSeparationEventTimes(evs);

    err = setNuKernelArgs(ci, ia, nu_step);
    if (err != CL_SUCCESS)
    {
        warn("Failed to set nu kernel argument\n");
        return err;
    }

    size_t chunkSize = ia->r_steps / numChunks;

    offset[0] = nu_step;
    offset[1] = 0;
    for (i = 0; i < numChunks; ++i)
    {
        offset[2] = i * chunkSize;

        err = enqueueIntegralKernel(ci, evs, offset, global, local);
        if (err != CL_SUCCESS)
        {
            warn("Failed to enqueue integral kernel: %s\n", showCLInt(err));
            return err;
        }

        #if 0
        clWaitForEvents(1, &evs->endTmp);
        printf("Step took %f\n", mwEventTime(evs->endTmp) * 1000.0);

        err = mwWaitReleaseEvent(&evs->endTmp);
        if (err != CL_SUCCESS)
        {
            warn("Failed to wait/release NDRange event: %s\n", showCLInt(err));
            return err;
        }
        #endif
    }

    return CL_SUCCESS;
}

static real runIntegral(CLInfo* ci,
                        SeparationCLMem* cm,
                        real* probs_results,
                        const ASTRONOMY_PARAMETERS* ap,
                        const INTEGRAL_AREA* ia)
{
    cl_uint i;
    cl_int err;
    real result;
    SeparationCLEvents evs;

    err = mwEnableProfiling(ci);
    if (err != CL_SUCCESS)
        warn("Failed to enable profiling: %s\n", showCLInt(err));

    err = separationSetOutputBuffers(ci, cm);
    if (err != CL_SUCCESS)
    {
        warn("Failed to set output buffer arguments: %s\n", showCLInt(err));
        return NAN;
    }

    size_t global[3];
    size_t local[3];

    /* TODO: Figure out a number based on GPU speed to keep execution
     * times under 100ms to prevent unusable system */
    size_t numChunks = 1;

    if (findWorkGroupSizes(ci, ia, numChunks, global, local))
    {
        warn("Failed to calculate acceptable work group sizes\n");
        return NAN;
    }

    printf("Sizes: ap = %zu ia = %zu, sc = %u * %zu, rcs = %u * %zu, sg_dx = %u * %zu total = %zu\n",
           sizeof(ASTRONOMY_PARAMETERS), sizeof(INTEGRAL_AREA),
           ap->number_streams, sizeof(STREAM_CONSTANTS), ia->r_steps, sizeof(R_CONSTS),
           ia->r_steps, sizeof(real),
           sizeof(ASTRONOMY_PARAMETERS) + sizeof(INTEGRAL_AREA) + ap->number_streams * sizeof(STREAM_CONSTANTS) + ia->r_steps * sizeof(R_CONSTS) + ia->r_steps * sizeof(real)

        );

    for (i = 0; i < ia->nu_steps; ++i)
    {
        double t1 = mwGetTimeMilli();
        err = runNuStep(ci, cm, &evs,
                        ia, numChunks, global, local,
                        ap->number_streams, i);
        if (err != CL_SUCCESS)
        {
            warn("Failed to run nu step: %s\n", showCLInt(err));
            return NAN;
        }

        err = clWaitForEvents(1, &evs.endTmp);
        if (err != CL_SUCCESS)
        {
            warn("Failed to wait for event: %s\n", showCLInt(err));
            return NAN;
        }

        double t2 = mwGetTimeMilli();
        printf("Loop time: %f ms\n", t2 - t1);
    }

    clFinish(ci->queue);

    /* Read results from final step */
    result = readKernelResults(ci, cm, &evs, probs_results, ia->mu_steps, ia->r_steps, ap->number_streams);
    if (isnan(result))
        warn("Failed to read final kernel results: %s\n", showCLInt(err));

    return result;
}

/* FIXME: This can only work right now for 1 integral */
real integrateCL(const ASTRONOMY_PARAMETERS* ap,
                 const INTEGRAL_AREA* ia,
                 const STREAM_CONSTANTS* sc,
                 const STREAM_GAUSS sg,
                 real* probs_results,
                 const CLRequest* clr)
{
    real result = NAN;
    CLInfo ci = EMPTY_CL_INFO;
    SeparationCLMem cm = EMPTY_SEPARATION_CL_MEM;

    if (setupSeparationCL(&ci, &cm, ap, ia, sc, sg, clr) != CL_SUCCESS)
        warn("Failed to setup up CL\n");
    else
        result = runIntegral(&ci, &cm, probs_results, ap, ia);

    destroyCLInfo(&ci);
    releaseSeparationBuffers(&cm);

    return result;
}

