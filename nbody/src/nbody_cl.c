/*
 * Copyright (c) 2011 Matthew Arsenault
 * Copyright (c) 2011 Rensselaer Polytechnic Institute
 *
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "milkyway_cl.h"
#include "milkyway_util.h"
#include "nbody_cl.h"
#include "nbody_show.h"
#include "nbody_util.h"
#include "nbody_curses.h"

/* CHECKME: Padding between these fields might be a good idea */
typedef struct NBODY_ALIGN
{
    real radius;
    int bottom;
    int maxDepth;
    int errorCode;
    unsigned int blkCnt;
} TreeStatus;

typedef struct
{
    real f[32];
    int i[64];
} Debug;


static void printNBodyWorkSizes(const NBodyWorkSizes* ws)
{
    mw_printf("\n"
              "Kernel launch sizes:\n"
              "  Bounding box kernel:  "ZU", "ZU"\n"
              "  Tree build kernel:    "ZU", "ZU"\n"
              "  Summarization kernel: "ZU", "ZU"\n"
              "  Sort kernel:          "ZU", "ZU"\n"
              "  Force kernel:         "ZU", "ZU"\n"
              "  Integration kernel:   "ZU", "ZU"\n"
              "\n",
              ws->global[0], ws->local[0],
              ws->global[1], ws->local[1],
              ws->global[2], ws->local[2],
              ws->global[3], ws->local[3],
              ws->global[4], ws->local[4],
              ws->global[5], ws->local[5]
        );
}

cl_bool setWorkSizes(NBodyWorkSizes* ws, const DevInfo* di)
{
    cl_uint blocks = di->maxCompUnits;

    ws->global[0] = ws->threads[0] * ws->factors[0] * blocks;
    ws->local[0] = ws->threads[0];

    ws->global[1] = ws->threads[1] * ws->factors[1] * blocks;
    ws->local[1] = ws->threads[1];

    ws->global[2] = ws->threads[2] * ws->factors[2] * blocks;
    ws->local[2] = ws->threads[2];

    ws->global[3] = ws->threads[3] * ws->factors[3] * blocks;
    ws->local[3] = ws->threads[3];

    ws->global[4] = ws->threads[4] * ws->factors[4] * blocks;
    ws->local[4] = ws->threads[4];

    ws->global[5] = ws->threads[5] * ws->factors[5] * blocks;
    ws->local[5] = ws->threads[5];

    return CL_FALSE;
}

/* Return CL_TRUE if some error */
cl_bool setThreadCounts(NBodyWorkSizes* ws, const DevInfo* di)
{
    ws->factors[0] = 1;
    ws->factors[1] = 2;
    ws->factors[2] = 1;
    ws->factors[3] = 1;
    ws->factors[4] = 1;
    ws->factors[5] = 1;

    if (di->devType == CL_DEVICE_TYPE_CPU)
    {
        ws->threads[0] = 1;
        ws->threads[1] = 1;
        ws->threads[2] = 1;
        ws->threads[3] = 1;
        ws->threads[4] = 1;
        ws->threads[5] = 1;
    }
    else if (computeCapabilityIs(di, 1, 3))
    {
        ws->threads[0] = 256;
        ws->threads[1] = 288;
        ws->threads[2] = 256;
        ws->threads[3] = 512;
        ws->threads[4] = 384;
        ws->threads[5] = 512;
    }
    else if (minComputeCapabilityCheck(di, 2, 0))
    {
        ws->threads[0] = 512;
        ws->threads[1] = 1024;
        ws->threads[2] = 1024;
        ws->threads[3] = 256;
        ws->threads[4] = 256;
        ws->threads[5] = 512;
    }
    else
    {
        ws->threads[0] = 256;
        ws->threads[1] = 256;
        ws->threads[2] = 256;
        ws->threads[3] = 256;
        ws->threads[4] = 256;
        ws->threads[5] = 256;
    }

    return CL_FALSE;
}


static void* mapBuffer(CLInfo* ci, cl_mem mem, cl_map_flags flags, size_t size)
{
    return clEnqueueMapBuffer(ci->queue, mem, CL_TRUE, flags, 0,
                              size,
                              0, NULL, NULL, NULL);

}

static void printDebug(const Debug* d)
{
    int i;
    for (i = 0; i < 32; ++i)
    {
        mw_printf("Debug.int[%d] = %d\n", i, d->i[i]);
    }

    for (i = 0; i < 32; ++i)
    {
        mw_printf("Debug.float[%d] = %.15f\n", i, d->f[i]);
    }
}

static void printTreeStatus(const TreeStatus* ts)
{
    mw_printf("TreeStatus = {\n"
              "  radius    = %.15f\n"
              "  bottom    = %d\n"
              "  maxDepth  = %d\n"
              "  errorCode = %d\n"
              "  blckCnt   = %u\n"
              "}\n",
              ts->radius,
              ts->bottom,
              ts->maxDepth,
              ts->errorCode,
              ts->blkCnt);
}

/* Set arguments (idx, idx + 1, idx + 2) to the buffers in mem[3] */
static cl_int setMemArrayArgs(cl_kernel kern, cl_mem mem[3], cl_uint idx)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;

    for (i = 0; i < 3; ++i)
    {
        err |= clSetKernelArg(kern, idx + i, sizeof(cl_mem), &mem[i]);
    }

    return err;
}

static cl_int nbodySetKernelArguments(cl_kernel kern, NBodyBuffers* nbb, cl_bool exact)
{
    cl_int err = CL_SUCCESS;
    cl_int step = 0;

    err |= setMemArrayArgs(kern, nbb->pos, 0);
    err |= setMemArrayArgs(kern, nbb->vel, 3);
    err |= setMemArrayArgs(kern, nbb->acc, 6);

    if (!exact)
    {
        err |= setMemArrayArgs(kern, nbb->max, 9);
        err |= setMemArrayArgs(kern, nbb->min, 12);

        err |= clSetKernelArg(kern, 15, sizeof(cl_mem), &nbb->masses);

        err |= clSetKernelArg(kern, 16, sizeof(cl_mem), &nbb->start);
        err |= clSetKernelArg(kern, 17, sizeof(cl_mem), &nbb->count);
        err |= clSetKernelArg(kern, 18, sizeof(cl_mem), &nbb->child);
        err |= clSetKernelArg(kern, 19, sizeof(cl_mem), &nbb->sort);
        err |= clSetKernelArg(kern, 20, sizeof(cl_mem), &nbb->treeStatus);


        if (nbb->critRadii)
        {
            err |= clSetKernelArg(kern, 23, sizeof(cl_mem), &nbb->critRadii);
        }
        else
        {
            /* Just be anything. Won't be used, it just needs to be not null to not error */
            err |= clSetKernelArg(kern, 23, sizeof(cl_mem), &nbb->treeStatus);
        }
    }
    else
    {
        /* For exact kernel just set whatever. These won't actually be used */
        err |= setMemArrayArgs(kern, nbb->pos, 9);
        err |= setMemArrayArgs(kern, nbb->pos, 12);
        err |= clSetKernelArg(kern, 15, sizeof(cl_mem), &nbb->pos[0]);
        err |= clSetKernelArg(kern, 16, sizeof(cl_mem), &nbb->pos[0]);
        err |= clSetKernelArg(kern, 17, sizeof(cl_mem), &nbb->pos[0]);
        err |= clSetKernelArg(kern, 18, sizeof(cl_mem), &nbb->pos[0]);
        err |= clSetKernelArg(kern, 19, sizeof(cl_mem), &nbb->pos[0]);
        err |= clSetKernelArg(kern, 20, sizeof(cl_mem), &nbb->pos[0]);

        err |= clSetKernelArg(kern, 23, sizeof(cl_mem), &nbb->pos[0]);
    }


    /* Set the step to be 0 here. Only the force calculation kernel
     * actually cares. We'll set it before then. */
    err |= clSetKernelArg(kern, 21, sizeof(cl_int), &step);

    /* Value also doesn't matter */
    err |= clSetKernelArg(kern, 22, sizeof(cl_int), &step);


    err |= clSetKernelArg(kern, 24, sizeof(cl_mem), &nbb->debug);

    return err;
}

cl_int nbodySetAllKernelArguments(NBodyState* st)
{
    cl_int err = CL_SUCCESS;
    NBodyKernels* k = st->kernels;
    cl_bool exact = st->usesExact;

    if (!exact)
    {
        err |= nbodySetKernelArguments(k->boundingBox, st->nbb, exact);
        err |= nbodySetKernelArguments(k->buildTree, st->nbb, exact);
        err |= nbodySetKernelArguments(k->summarization, st->nbb, exact);
        err |= nbodySetKernelArguments(k->sort, st->nbb, exact);
        err |= nbodySetKernelArguments(k->forceCalculation, st->nbb, exact);
    }
    else
    {
        err |= nbodySetKernelArguments(k->forceCalculation_Exact, st->nbb, exact);
    }


    err |= nbodySetKernelArguments(k->integration, st->nbb, exact);

    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error setting kernel arguments", err);
    }

    return err;
}

cl_int nbodyCreateKernels(NBodyState* st)
{
    cl_int err = CL_SUCCESS;
    NBodyKernels* kernels = st->kernels;
    CLInfo* ci = st->ci;

    err |= mwCreateKernel(&kernels->boundingBox, ci, "boundingBox");
    err |= mwCreateKernel(&kernels->buildTree, ci, "buildTree");
    err |= mwCreateKernel(&kernels->summarization, ci, "summarization");
    err |= mwCreateKernel(&kernels->sort, ci, "sort");
    err |= mwCreateKernel(&kernels->forceCalculation, ci, "forceCalculation");
    err |= mwCreateKernel(&kernels->integration, ci, "integration");

    err |= mwCreateKernel(&kernels->forceCalculation_Exact, ci, "forceCalculation_Exact");

    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error creating kernels", err);
    }

    return err;
}

static cl_int clReleaseKernel_quiet(cl_kernel kern)
{
    if (!kern)
        return CL_SUCCESS;

    return clReleaseKernel(kern);
}

cl_int nbodyReleaseKernels(NBodyState* st)
{
    cl_int err = CL_SUCCESS;
    NBodyKernels* kernels = st->kernels;

    err |= clReleaseKernel_quiet(kernels->boundingBox);
    err |= clReleaseKernel_quiet(kernels->buildTree);
    err |= clReleaseKernel_quiet(kernels->summarization);
    err |= clReleaseKernel_quiet(kernels->sort);
    err |= clReleaseKernel_quiet(kernels->forceCalculation);
    err |= clReleaseKernel_quiet(kernels->integration);

    if (err != CL_SUCCESS)
        mwCLWarn("Error releasing kernels", err);

    return err;
}

static cl_uint findNNode(const DevInfo* di, cl_int nbody)
{
    cl_uint nNode = 2 * nbody;

    if (nNode < 1024 * di->maxCompUnits)
        nNode = 1024 * di->maxCompUnits;
    while ((nNode & (di->warpSize - 1)) != 0)
        ++nNode;

    return nNode - 1;
}

static char* getCompileFlags(const NBodyCtx* ctx, const NBodyState* st, const DevInfo* di)
{
    char* buf;
    const NBodyWorkSizes* ws = st->workSizes;
    const Potential* p = &ctx->pot;

    /* Put a space between the -D. if the thing begins with D, there's
     * an Apple OpenCL compiler bug where the D will be
     * stripped. -DDOUBLEPREC=1 will actually define OUBLEPREC */
    if (asprintf(&buf,
                 "-D DOUBLEPREC=%d "
               #if !DOUBLEPREC
                 "-cl-single-precision-constant "
               #endif

                 "-D NBODY=%d "
                 "-D NNODE=%u "
                 "-D WARPSIZE=%u "

                 "-D NOSORT=%d "

                 "-D THREADS1="ZU" "
                 "-D THREADS2="ZU" "
                 "-D THREADS3="ZU" "
                 "-D THREADS4="ZU" "
                 "-D THREADS5="ZU" "
                 "-D THREADS6="ZU" "
                 "-D MAXDEPTH=%u "

                 "-D TIMESTEP=%a "
                 "-D EPS2=%a "
                 "-D THETA=%a "

                 "-D NEWCRITERION=%d "
                 "-D SW93=%d "
                 "-D BH86=%d "
                 "-D EXACT=%d "

                 /* Potential */
                 "-D USE_EXTERNAL_POTENTIAL=%d "
                 "-D MIYAMOTO_NAGAI_DISK=%d "
                 "-D EXPONENTIAL_DISK=%d "
                 "-D LOG_HALO=%d "
                 "-D NFW_HALO=%d "
                 "-D TRIAXIAL_HALO=%d "

                 /* Spherical constants */
                 "-DSPHERICAL_MASS=%a "
                 "-DSPHERICAL_SCALE=%a "

                 /* Disk constants */
                 "-DDISK_MASS=%a "
                 "-DDISK_SCALE_LENGTH=%a "
                 "-DDISK_SCALE_HEIGHT=%a "

                 /* Halo constants */
                 "-DHALO_VHALO=%a "
                 "-DHALO_SCALE_LENGTH=%a "
                 "-DHALO_FLATTEN_Z=%a "
                 "-DHALO_FLATTEN_Y=%a "
                 "-DHALO_FLATTEN_X=%a "
                 "-DHALO_TRIAX_ANGLE=%a "
                 "-DHALO_C1=%a "
                 "-DHALO_C2=%a "
                 "-DHALO_C3=%a "

                 "%s ",
                 DOUBLEPREC,

                 st->nbody,
                 findNNode(di, st->nbody),
                 di->warpSize,

                 (di->devType == CL_DEVICE_TYPE_CPU),

                 ws->threads[0],
                 ws->threads[1],
                 ws->threads[2],
                 ws->threads[3],
                 ws->threads[4],
                 ws->threads[5],
                 NBODY_MAXDEPTH,

                 ctx->timestep,
                 ctx->eps2,
                 ctx->theta,

                 /* SEt criterion */
                 ctx->criterion == NewCriterion,
                 ctx->criterion == SW93,
                 ctx->criterion == BH86,
                 ctx->criterion == Exact,


                 /* Set potential */
                 ctx->potentialType == EXTERNAL_POTENTIAL_DEFAULT,

                 p->disk.type == MiyamotoNagaiDisk,
                 p->disk.type == ExponentialDisk,
                 p->halo.type == LogarithmicHalo,
                 p->halo.type == NFWHalo,
                 p->halo.type == TriaxialHalo,

                 /* Set potential constants */
                 /* Spherical constants */
                 p->sphere[0].mass,
                 p->sphere[0].scale,

                 /* Disk constants */
                 p->disk.mass,
                 p->disk.scaleLength,
                 p->disk.scaleHeight,

                 /* Halo constants */
                 p->halo.vhalo,
                 p->halo.scaleLength,
                 p->halo.flattenZ,
                 p->halo.flattenY,
                 p->halo.flattenX,
                 p->halo.triaxAngle,
                 p->halo.c1,
                 p->halo.c2,
                 p->halo.c3,

                 /* Misc. other stuff */
                 hasNvidiaCompilerFlags(di) ? "-cl-nv-verbose" : ""
            ) < 1)
    {
        mw_printf("Error getting compile flags\n");
        return NULL;
    }

    return buf;
}

cl_int nbodyLoadKernels(const NBodyCtx* ctx, NBodyState* st)
{
    CLInfo* ci = st->ci;
    cl_int err = CL_SUCCESS;
    char* src = NULL;
    char* compileFlags = NULL;

    src = mwReadFile("kernels/nbody_kernels.cl");
    if (!src)
    {
        mw_printf("Failed to read N-body kernel source\n");
        return MW_CL_ERROR;
    }

    compileFlags = getCompileFlags(ctx, st, &ci->di);
    assert(compileFlags);
    err = mwSetProgramFromSrc(ci, (const char**) &src, 1, compileFlags);
    if (err != CL_SUCCESS)
    {
        mw_printf("Failed build flags: %s\n", compileFlags);
    }

    free(src);
    free(compileFlags);

    return err;
}

/* Return CL_FALSE if device isn't capable of running this */
cl_bool nbodyCheckDevCapabilities(const DevInfo* di, const NBodyCtx* ctx, NBodyState* st)
{
    (void) di;
    (void) ctx;
    (void) st;
    return CL_TRUE;
}

static cl_int debug(CLInfo* ci, NBodyBuffers* nbb)
{
    cl_int err;
    Debug d;

    memset(&d, 0, sizeof(d));

    err = clEnqueueReadBuffer(ci->queue,
                              nbb->debug,
                              CL_TRUE,
                              0, sizeof(d), &d,
                              0, NULL, NULL);

    mw_printf("Debug:\n");
    printDebug(&d);
    return err;
}

static cl_int readTreeStatus(TreeStatus* tc, CLInfo* ci, NBodyBuffers* nbb)
{
    assert(tc);

    return clEnqueueReadBuffer(ci->queue,
                               nbb->treeStatus,
                               CL_TRUE,
                               0, sizeof(*tc), tc,
                               0, NULL, NULL);
}

static cl_int printBuffer(CLInfo* ci, cl_mem mem, size_t n, const char* name, int type)
{
    size_t i;
    void* p;

    p = mapBuffer(ci, mem, CL_MAP_READ, n * (type == 0 ? sizeof(real) : sizeof(int)));
    if (!p)
        return MW_CL_ERROR;

    if (type == 0)
    {
        const real* pr = (const real*) p;
        for (i = 0; i < n; ++i)
        {
            mw_printf("%s["ZU"] = %.15f\n", name, i, pr[i]);
        }
    }
    else
    {
        const int* pi = (const int*) p;
        for (i = 0; i < n; ++i)
        {
            mw_printf("%s["ZU"] = %d\n", name, i, pi[i]);
        }
    }

    return clEnqueueUnmapMemObject(ci->queue, mem, p, 0, NULL, NULL);
}

static void stdDebugPrint(CLInfo* ci, NBodyState* st, NBodyBuffers* nbb)
{
    mw_printf("--------------------------------------------------------------------------------\n");
    cl_int err;
    cl_uint nNode = findNNode(&ci->di, st->nbody);
    mw_printf("BEGIN CHILD\n");
    printBuffer(ci, nbb->child, NSUB * (nNode + 1), "child", 1);
    mw_printf("BEGIN START\n");
    printBuffer(ci, nbb->start, nNode, "start", 1);

    mw_printf("BEGIN MASS\n");
    printBuffer(ci, nbb->masses, nNode + 1, "mass", 0);

    {
        TreeStatus tc;
        memset(&tc, 0, sizeof(tc));
        err = readTreeStatus(&tc, ci, nbb);
        if (err != CL_SUCCESS)
            mwCLWarn("Reading tree status failed\n", err);
        else
            printTreeStatus(&tc);
    }

    debug(ci, nbb);
    mw_printf("--------------------------------------------------------------------------------\n");
}

/* Check the error code */
static cl_bool checkKernelErrorCode(CLInfo* ci, NBodyBuffers* nbb)
{
    cl_int err;
    TreeStatus ts;

    err = readTreeStatus(&ts, ci, nbb);
    if (err != CL_SUCCESS)
        return CL_TRUE;

    if (ts.errorCode != 0)
    {
        mw_printf("Kernel reported error: %d (%s)\n", ts.errorCode, showNBodyKernelError(ts.errorCode));
        return CL_TRUE;
    }

    return CL_FALSE;
}

static cl_double waitReleaseEventWithTime(cl_event ev)
{
    cl_double t;
    cl_int err;

    err = clWaitForEvents(1, &ev);
    if (err != CL_SUCCESS)
        return 0.0;

    t = mwEventTimeMS(ev);

    err = clReleaseEvent(ev);
    if (err != CL_SUCCESS)
        return 0.0;

    return t;
}

/* If kernels are taking too long, try running with fewer threads if
 * possible. Return CL_TRUE if anything changed
 */
static cl_bool reevaluateWorkSizes(NBodyWorkSizes* ws)
{
    cl_double buildTreeThreshold = 33.0;
    cl_double forceThreshold = 33.0;

    if (ws->chunkTimings[1] > buildTreeThreshold)
    {
        /* Try fewer threads if possible */
    }

    if (ws->chunkTimings[4] > forceThreshold)
    {
        /* Try fewer threads if possible */
    }

    return CL_FALSE;
}

static void reportProgressWithTimings(const NBodyCtx* ctx, const NBodyState* st)
{
    NBodyWorkSizes* ws = st->workSizes;

    mw_mvprintw(0, 0,
                "Step %d (%f%%):\n"
                "  boundingBox:      %15f ms\n"
                "  buildTree:        %15f ms%15f ms\n"
                "  summarization:    %15f ms\n"
                "  sort:             %15f ms\n"
                "  forceCalculation: %15f ms%15f ms\n"
                "  integration:      %15f ms\n"
                "\n",
                st->step,
                100.0 * st->tnow / ctx->timeEvolve,
                ws->timings[0],
                ws->timings[1], ws->chunkTimings[1],
                ws->timings[2],
                ws->timings[3],
                ws->timings[4], ws->chunkTimings[4],
                ws->timings[5]);
    mw_refresh();
}

/* Run kernels used only by tree versions. */
static cl_int executeTreeConstruction(NBodyState* st)
{
    cl_int err;
    size_t chunk;
    size_t nChunk;
    int upperBound;
    size_t offset[1];
    CLInfo* ci = st->ci;
    NBodyWorkSizes* ws = st->workSizes;
    NBodyKernels* kernels = st->kernels;
    cl_event boxEv, sumEv, sortEv;

    err = clEnqueueNDRangeKernel(ci->queue, kernels->boundingBox, 1,
                                 NULL, &ws->global[0], &ws->local[0],
                                 0, NULL, &boxEv);
    if (err != CL_SUCCESS)
        return err;

    ws->timings[0] += waitReleaseEventWithTime(boxEv);

    nChunk     = st->ignoreResponsive ?         1 : mwDivRoundup((size_t) st->nbody, ws->global[1]);
    upperBound = st->ignoreResponsive ? st->nbody : (int) ws->global[1];
    offset[0] = 0;
    for (chunk = 0; chunk < nChunk; ++chunk)
    {
        cl_event ev;

        if (upperBound > st->nbody)
            upperBound = st->nbody;

        err = clSetKernelArg(kernels->buildTree, 22, sizeof(int), &upperBound);
        if (err != CL_SUCCESS)
            return err;

        err = clEnqueueNDRangeKernel(ci->queue, kernels->buildTree, 1,
                                     offset, &ws->global[1], &ws->local[1],
                                     0, NULL, &ev);
        if (err != CL_SUCCESS)
            return err;

        ws->timings[1] += waitReleaseEventWithTime(ev);

        upperBound += (int) ws->global[1];
        offset[0] += ws->global[1];
    }
    ws->chunkTimings[1] = ws->timings[1] / (double) nChunk;

    err = clEnqueueNDRangeKernel(ci->queue, kernels->summarization, 1,
                                 NULL, &ws->global[2], &ws->local[2],
                                 0, NULL, &sumEv);
    if (err != CL_SUCCESS)
        return err;

    ws->timings[2] += waitReleaseEventWithTime(sumEv);

    /* FIXME: This does not work unless ALL of the threads are
     * launched at once. This may be bad when we need
     * responsiveness. This also means it will always hang with
     * CPUs */
    err = clEnqueueNDRangeKernel(ci->queue, kernels->sort, 1,
                                 NULL, &ws->global[3], &ws->local[3],
                                 0, NULL, &sortEv);
    if (err != CL_SUCCESS)
        return err;

    ws->timings[3] += waitReleaseEventWithTime(sortEv);

    return CL_SUCCESS;
}

/* Run force calculation and integration kernels */
static cl_int executeForceKernels(NBodyState* st)
{
    cl_int err;
    size_t chunk = 0;
    CLInfo* ci = st->ci;
    NBodyWorkSizes* ws = st->workSizes;
    NBodyKernels* kernels = st->kernels;
    size_t nChunk = st->ignoreResponsive ?          1 : mwDivRoundup((size_t) st->nbody, ws->global[4]);
    cl_int upperBound = st->ignoreResponsive ? st->nbody : (cl_int) ws->global[4];
    size_t offset[1] = { 0 };
    cl_event integrateEv;
    cl_kernel forceKern = st->usesExact ? kernels->forceCalculation_Exact : kernels->forceCalculation;


    for (chunk = 0; chunk < nChunk; ++chunk)
    {
        cl_event ev;
        double t;

        if (upperBound > st->nbody)
            upperBound = st->nbody;

        err = clSetKernelArg(forceKern, 22, sizeof(cl_int), &upperBound);
        if (err != CL_SUCCESS)
            return err;

        err = clEnqueueNDRangeKernel(ci->queue, forceKern, 1,
                                     offset, &ws->global[4], &ws->local[4],
                                     0, NULL, &ev);

        if (err != CL_SUCCESS)
            return err;

        t = waitReleaseEventWithTime(ev);

        ws->timings[4] += t;
        upperBound += (cl_int) ws->global[4];
        offset[0] += ws->global[4];
    }
    ws->chunkTimings[4] = ws->timings[4] / (double) nChunk;


    err = clEnqueueNDRangeKernel(ci->queue, kernels->integration, 1,
                                 NULL, &ws->global[5], &ws->local[5],
                                 0, NULL, &integrateEv);
    if (err != CL_SUCCESS)
        return err;

    ws->timings[5] += waitReleaseEventWithTime(integrateEv);


    return CL_SUCCESS;
}

static cl_int stepSystemCL(const NBodyCtx* ctx, NBodyState* st)
{
    cl_int err;
    cl_uint i;
    CLInfo* ci = st->ci;
    NBodyWorkSizes* ws = st->workSizes;
    NBodyKernels* kernels = st->kernels;
    cl_int step = st->step;

    memset(ws->timings, 0, sizeof(ws->timings));

    if (!st->usesExact)
    {
        err = clSetKernelArg(kernels->forceCalculation, 21, sizeof(cl_int), &step);
        if (err != CL_SUCCESS)
            return err;

        err = executeTreeConstruction(st);
        if (err != CL_SUCCESS)
            return err;
    }
    else
    {
        err = clSetKernelArg(kernels->forceCalculation_Exact, 21, sizeof(cl_int), &step);
        if (err != CL_SUCCESS)
            return err;
    }

    err = executeForceKernels(st);
    if (err != CL_SUCCESS)
        return err;

    err = clFinish(ci->queue);
    if (err != CL_SUCCESS)
        return err;

    if (st->reportProgress)
    {
        reportProgressWithTimings(ctx, st);
    }

    for (i = 0; i < 6; ++i) /* Add timings to running totals */
    {
        ws->kernelTimings[i] += ws->timings[i];
    }

    return err;
}

static cl_int nbodyMainLoop(const NBodyCtx* ctx, NBodyState* st)
{
    cl_int err = CL_SUCCESS;
    const real tstop = ctx->timeEvolve - ctx->timestep / 1024.0;

    st->tnow = 0;
    st->step = 0;
    while (err == CL_SUCCESS && st->tnow < tstop)
    {
        st->dirty = TRUE;
        if (!st->usesExact && checkKernelErrorCode(st->ci, st->nbb))
        {
            err = MW_CL_ERROR;
            break;
        }

        err = stepSystemCL(ctx, st);

        st->tnow += ctx->timestep;
        st->step++;
    }

    return err;
}

/* This is dumb and errors if mem isn't set */
static cl_int clReleaseMemObject_quiet(cl_mem mem)
{
    if (!mem)
        return CL_SUCCESS;

    return clReleaseMemObject(mem);
}

cl_int nbodyReleaseBuffers(NBodyState* st)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;
    NBodyBuffers* nbb = st->nbb;

    for (i = 0; i < 3; ++i)
    {
        err |= clReleaseMemObject_quiet(nbb->pos[i]);
        err |= clReleaseMemObject_quiet(nbb->vel[i]);
        err |= clReleaseMemObject_quiet(nbb->acc[i]);
        err |= clReleaseMemObject_quiet(nbb->max[i]);
        err |= clReleaseMemObject_quiet(nbb->min[i]);
    }

    err |= clReleaseMemObject_quiet(nbb->masses);
    err |= clReleaseMemObject_quiet(nbb->treeStatus);
    err |= clReleaseMemObject_quiet(nbb->start);
    err |= clReleaseMemObject_quiet(nbb->count);
    err |= clReleaseMemObject_quiet(nbb->child);
    err |= clReleaseMemObject_quiet(nbb->sort);

    if (err != CL_SUCCESS)
        mwCLWarn("Error releasing buffers", err);

    return err;
}


cl_int nbodySetInitialTreeStatus(NBodyState* st)
{
    cl_int err;
    TreeStatus iniTreeStatus;
    CLInfo* ci = st->ci;
    NBodyBuffers* nbb = st->nbb;

    memset(&iniTreeStatus, 0, sizeof(iniTreeStatus));

    iniTreeStatus.radius = 0;
    iniTreeStatus.bottom = 0;
    iniTreeStatus.maxDepth = 1;
    iniTreeStatus.errorCode = 0;
    iniTreeStatus.blkCnt = 0;

    err = clEnqueueWriteBuffer(ci->queue,
                               nbb->treeStatus,
                               CL_TRUE,
                               0, sizeof(TreeStatus), &iniTreeStatus,
                               0, NULL, NULL);

    if (err != CL_SUCCESS)
    {
        mwCLWarn("Error writing initial tree status", err);
    }

    return err;
}

static cl_uint findInc(cl_uint warpSize, cl_uint nbody)
{
    return (nbody + warpSize - 1) & (-warpSize);
}

cl_int nbodyCreateBuffers(const NBodyCtx* ctx, NBodyState* st)
{
    cl_uint i;
    CLInfo* ci = st->ci;
    NBodyBuffers* nbb = st->nbb;
    size_t massSize;
    cl_uint nNode = findNNode(&ci->di, st->nbody);

    for (i = 0; i < 3; ++i)
    {
        nbb->pos[i] = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
        nbb->vel[i] = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(real));
        nbb->acc[i] = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(real));

        if (!nbb->pos[i] || !nbb->vel[i] || !nbb->acc[i])
        {
            return MW_CL_ERROR;
        }

        if (ctx->criterion != Exact)
        {
            nbb->min[i] = mwCreateZeroReadWriteBuffer(ci, ci->di.maxCompUnits * sizeof(real));
            nbb->max[i] = mwCreateZeroReadWriteBuffer(ci, ci->di.maxCompUnits * sizeof(real));
            if (!nbb->min[i] || !nbb->max[i])
            {
                return MW_CL_ERROR;
            }
        }
    }

    nbb->debug = mwCreateZeroReadWriteBuffer(ci, sizeof(Debug));
    if (!nbb->debug)
    {
        return MW_CL_ERROR;
    }

    massSize = st->usesExact ? st->nbody * sizeof(real) : (nNode + 1) * sizeof(real);
    nbb->masses = mwCreateZeroReadWriteBuffer(ci, massSize);
    if (!nbb->masses)
    {
        return MW_CL_ERROR;
    }

    /* If we are doing an exact Nbody, we don't need the rest */
    if (ctx->criterion != Exact)
    {
        nbb->treeStatus = mwCreateZeroReadWriteBuffer(ci, sizeof(TreeStatus));
        nbb->start = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(cl_int));
        nbb->count = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(cl_int));
        nbb->sort = mwCreateZeroReadWriteBuffer(ci, st->nbody * sizeof(cl_int));
        nbb->child = mwCreateZeroReadWriteBuffer(ci, NSUB * (nNode + 1) * sizeof(cl_int));

        if (!nbb->treeStatus || !nbb->start || !nbb->count || !nbb->sort || !nbb->child)
            return MW_CL_ERROR;

        if (ctx->criterion == SW93 || ctx->criterion == NewCriterion)
        {
            /* This only is for cells, so we could subtract nbody if we wanted */
            nbb->critRadii = mwCreateZeroReadWriteBuffer(ci, (nNode + 1) * sizeof(real));
            if (!nbb->critRadii)
                return MW_CL_ERROR;
        }
    }

    return CL_SUCCESS;
}

static cl_int mapBodies(real* pos[3], real* vel[3], real** mass, NBodyBuffers* nbb, CLInfo* ci, cl_map_flags flags, NBodyState* st)
{
    cl_uint i;

    for (i = 0; i < 3; ++i)
    {
        pos[i] = (real*) mapBuffer(ci, nbb->pos[i], flags, st->nbody * sizeof(real));
        vel[i] = (real*) mapBuffer(ci, nbb->vel[i], flags, st->nbody * sizeof(real));
        if (!pos[i] || !vel[i])
        {
            return MW_CL_ERROR;
        }
    }

    *mass = (real*) mapBuffer(ci, nbb->masses, flags, st->nbody * sizeof(real));

    return CL_SUCCESS;
}

static cl_int unmapBodies(real* pos[3], real* vel[3], real* mass, NBodyBuffers* nbb, CLInfo* ci)
{
    cl_uint i;
    cl_int err = CL_SUCCESS;

    for (i = 0; i < 3; ++i)
    {
        if (pos[i])
            err |= clEnqueueUnmapMemObject(ci->queue, nbb->pos[i], pos[i], 0, NULL, NULL);

        if (vel[i])
            err |= clEnqueueUnmapMemObject(ci->queue, nbb->vel[i], vel[i], 0, NULL, NULL);
    }

    if (mass)
        err |= clEnqueueUnmapMemObject(ci->queue, nbb->masses, mass, 0, NULL, NULL);

    return err;
}

/* If last parameter is true, copy to the buffers. if false, copy from the buffers */
cl_int nbodyMarshalBodies(NBodyState* st, cl_bool marshalIn)
{
    cl_int i;
    cl_int err = CL_SUCCESS;
    const Body* b;
    real* pos[3] = { NULL, NULL, NULL };
    real* vel[3] = { NULL, NULL, NULL };
    real* mass = NULL;
    CLInfo* ci = st->ci;
    NBodyBuffers* nbb = st->nbb;
    cl_map_flags flags = marshalIn ? CL_MAP_WRITE : CL_MAP_READ;

    err = mapBodies(pos, vel, &mass, nbb, ci, flags, st);
    if (err != CL_SUCCESS)
    {
        unmapBodies(pos, vel, mass, nbb, ci);
        return err;
    }

    if (marshalIn)
    {
        for (i = 0, b = st->bodytab; b < st->bodytab + st->nbody; ++i, ++b)
        {
            pos[0][i] = X(Pos(b));
            pos[1][i] = Y(Pos(b));
            pos[2][i] = Z(Pos(b));

            vel[0][i] = X(Vel(b));
            vel[1][i] = Y(Vel(b));
            vel[2][i] = Z(Vel(b));

            mass[i] = Mass(b);
        }

        st->dirty = FALSE;
    }
    else
    {
        for (i = 0, b = st->bodytab; b < st->bodytab + st->nbody; ++i, ++b)
        {
            X(Pos(b)) = pos[0][i];
            Y(Pos(b)) = pos[1][i];
            Z(Pos(b)) = pos[2][i];

            X(Vel(b)) = vel[0][i];
            Y(Vel(b)) = vel[1][i];
            Z(Vel(b)) = vel[2][i];

            Mass(b) = mass[i];
        }
    }

    return unmapBodies(pos, vel, mass, nbb, ci);
}

void nbodyPrintKernelTimings(const NBodyState* st)
{
    double totalTime = 0.0;
    cl_uint i;
    double nStep = (double) st->step;
    const double* kernelTimings = st->workSizes->kernelTimings;

    for (i = 0; i < 6; ++i)
    {
        totalTime += kernelTimings[i];
    }

    mw_printf("\n--------------------------------------------------------------------------------\n"
              "Total timing over %d steps:\n"
              "                         Average             Total            Fraction\n"
              "                    ----------------   ----------------   ----------------\n"
              "  boundingBox:      %16f   %16f   %15.4f%%\n"
              "  buildTree:        %16f   %16f   %15.4f%%\n"
              "  summarization:    %16f   %16f   %15.4f%%\n"
              "  sort:             %16f   %16f   %15.4f%%\n"
              "  forceCalculation: %16f   %16f   %15.4f%%\n"
              "  integration:      %16f   %16f   %15.4f%%\n"
              "\n--------------------------------------------------------------------------------\n"
              "\n",
              st->step,
              kernelTimings[0] / nStep, kernelTimings[0], 100.0 * kernelTimings[0] / totalTime,
              kernelTimings[1] / nStep, kernelTimings[1], 100.0 * kernelTimings[1] / totalTime,
              kernelTimings[2] / nStep, kernelTimings[2], 100.0 * kernelTimings[2] / totalTime,
              kernelTimings[3] / nStep, kernelTimings[3], 100.0 * kernelTimings[3] / totalTime,
              kernelTimings[4] / nStep, kernelTimings[4], 100.0 * kernelTimings[4] / totalTime,
              kernelTimings[5] / nStep, kernelTimings[5], 100.0 * kernelTimings[5] / totalTime
        );
}


NBodyStatus runSystemCL(const NBodyCtx* ctx, NBodyState* st, const NBodyFlags* nbf)
{
    cl_int err = CL_SUCCESS;

    (void) nbf;

    err = nbodyMainLoop(ctx, st);
    if (err != CL_SUCCESS)
        return NBODY_CL_ERROR;

    err = nbodyMarshalBodies(st, CL_FALSE);
    if (err != CL_SUCCESS)
        return NBODY_CL_ERROR;

    return NBODY_SUCCESS;
}
