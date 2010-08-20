/*
Copyright 2008-2010 Travis Desell, Dave Przybylo, Nathan Cole, Matthew
Arsenault, Boleslaw Szymanski, Heidi Newberg, Carlos Varela, Malik
Magdon-Ismail and Rensselaer Polytechnic Institute.

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

#include <stdlib.h>
#include <stdio.h>

#include "separation_types.h"
#include "evaluation_state.h"
#include "milkyway_util.h"
#include "calculated_constants.h"
#include "integrals.h"
#include "likelihood.h"
#include "star_points.h"
#include "run_cl.h"

static void final_stream_integrals(FINAL_STREAM_INTEGRALS* fsi,
                                   const EVALUATION_STATE* es,
                                   const unsigned int number_streams,
                                   const unsigned int number_integrals)
{
    unsigned int i, j;

    fsi->stream_integrals = callocSafe(number_streams, sizeof(double));

    fsi->background_integral = es->integrals[0].background_integral;
    for (i = 0; i < number_streams; ++i)
        fsi->stream_integrals[i] = es->integrals[0].stream_integrals[i];

    for (i = 1; i < number_integrals; ++i)
    {
        fsi->background_integral -= es->integrals[i].background_integral;
        for (j = 0; j < number_streams; j++)
            fsi->stream_integrals[j] -= es->integrals[i].stream_integrals[j];
    }

}

static void free_final_stream_integrals(FINAL_STREAM_INTEGRALS* fsi)
{
    free(fsi->stream_integrals);
}

static void print_stream_integrals(const FINAL_STREAM_INTEGRALS* fsi, const unsigned int number_streams)
{
    unsigned int i;
    fprintf(stderr, "<background_integral> %.20lf </background_integral>\n", fsi->background_integral);
    fprintf(stderr, "<stream_integrals>");
    for (i = 0; i < number_streams; i++)
        fprintf(stderr, " %.20lf", fsi->stream_integrals[i]);
    fprintf(stderr, " </stream_integrals>\n");
}

inline static void calculate_stream_integrals(const ST_PROBS* probs,
                                              double* stream_integrals,
                                              const unsigned int number_streams)
{
    unsigned int i;

    for (i = 0; i < number_streams; ++i)
        stream_integrals[i] = probs[i].st_prob_int + probs[i].st_prob_int_c;
}

/* Add up completed integrals for progress reporting */
inline static double completed_integral_progress(const ASTRONOMY_PARAMETERS* ap,
                                                 const EVALUATION_STATE* es)
{
    INTEGRAL_AREA* ia;
    unsigned int i, current_calc_probs = 0;

    for (i = 0; i < es->current_integral; ++i)
    {
        ia = &ap->integral[i];
        current_calc_probs += ia->r_steps * ia->mu_steps * ia->nu_steps;
    }

    return current_calc_probs;
}

/* returns background integral */
static double integrate(const ASTRONOMY_PARAMETERS* ap,
                        const INTEGRAL_AREA* ia,
                        const STREAM_CONSTANTS* sc,
                        const STREAM_GAUSS* sg,
                        ST_PROBS* probs,
                        EVALUATION_STATE* es)
{
    double result;

    NU_CONSTANTS* nu_consts = prepare_nu_constants(ia->nu_steps, ia->nu_step_size, ia->nu_min);
    R_POINTS* r_pts = mallocSafe(sizeof(R_POINTS) * ap->convolve);
    vector* xyz = mallocSafe(sizeof(vector) * ap->convolve);

    result = r_sum(ap, ia, sc, sg, nu_consts, r_pts, probs, xyz, es);
    es->r_step = 0;

    free(nu_consts);
    free(r_pts);
    free(xyz);

    return result;
}

static void calculate_integrals(const ASTRONOMY_PARAMETERS* ap,
                                const STREAM_CONSTANTS* sc,
                                const STREAM_GAUSS* sg,
                                EVALUATION_STATE* es)
{
    INTEGRAL* integral;
    INTEGRAL_AREA* ia;

    double t1, t2;

    for (; es->current_integral < ap->number_integrals; es->current_integral++)
    {
        integral = &es->integrals[es->current_integral];
        ia = &ap->integral[es->current_integral];
        es->current_calc_probs = completed_integral_progress(ap, es);

        t1 = get_time();
        //integral->background_integral = integrate(ap, ia, sc, sg, integral->probs, es);
        integral->background_integral = integrateCL(ap, ia, sc, sg, integral->probs);
        t2 = get_time();

        printf("Time = %.20g\n", t2 - t1);

        calculate_stream_integrals(integral->probs, integral->stream_integrals, ap->number_streams);

        CLEAR_BG_PROB(es->r_acc);
    }

}

double evaluate(const ASTRONOMY_PARAMETERS* ap,
                const STREAMS* streams,
                const STREAM_CONSTANTS* sc,
                const char* star_points_file)
{
    double likelihood_val;
    EVALUATION_STATE es = EMPTY_EVALUATION_STATE;
    STREAM_GAUSS* sg;
    FINAL_STREAM_INTEGRALS fsi;
    STAR_POINTS sp = EMPTY_STAR_POINTS;

    initialize_state(ap, &es);
    sg = get_stream_gauss(ap->convolve);

  #if BOINC_APPLICATION && !SEPARATION_OPENCL
    if (boinc_file_exists(CHECKPOINT_FILE))
    {
        warn("Checkpoint exists. Attempting to resume from it\n");

        if (read_checkpoint(&es))
        {
            warn("Reading checkpoint failed\n");
            boinc_delete_file(CHECKPOINT_FILE);
            mw_finish(EXIT_FAILURE);
        }
    }
  #endif

    calculate_integrals(ap, sc, sg, &es);

  #if BOINC_APPLICATION && !SEPARATION_OPENCL
    /* Final checkpoint. */
    write_checkpoint(&es);
  #endif

    final_stream_integrals(&fsi, &es, ap->number_streams, ap->number_integrals);
    print_stream_integrals(&fsi, ap->number_streams);
    free_evaluation_state(&es);

    read_star_points(&sp, star_points_file);

    likelihood_val = likelihood(ap, &sp, sc, streams, &fsi, sg);

    free_star_points(&sp);
    free_final_stream_integrals(&fsi);
    free(sg);

  #if BOINC_APPLICATION && !SEPARATION_OPENCL
    boinc_delete_file(CHECKPOINT_FILE);
  #endif

    return likelihood_val;
}


