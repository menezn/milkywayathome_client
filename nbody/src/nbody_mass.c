/*
 * Copyright (c) 2012 Rensselaer Polytechnic Institute
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

#include "nbody_mass.h"

#include "milkyway_math.h"
#include "nbody_types.h"


/*In order to decrease the size of the numbers
 * computed all these functions are
 * calculated in log space*/

static real factorial(int n){
     int counter;
     real result = 0;

     for (counter = n; counter >= 1; counter--)
       {
          result += mw_log(counter);
       }

     return result;
  }


static real choose(int n, int c)
{
    unsigned int i;
    real result = 0;

    /* This for loop calulates log(n!/(n-c)!) */
    for (i = n - c + 1; i <= (unsigned int) n; ++i)
    {
        result += mw_log(i);
    }

    result -= factorial(c);
    return result;
}

real probability_match(int n, int k, real pobs)
{
    real result = 0;

    //result +=  (real) mw_log(choose(n, k));                                                                                                                                                        
    //result += mw_log(pobs) * (real) k;                                                                                                                                                             
    //result += mw_log(1.0 - pobs) * (real)(n - k);                                                                                                                                                  

    //The previous calculation does not return the right values.  Furthermore, we need a zeroed metric.                                                                                              
    result = (real)  choose(n,k) + k * mw_log(pobs) + (n-k) * mw_log(1 - pobs);
    //mw_printf("Coeff = %10.10f\n",choose(n,k));
    //mw_printf("Term 1 = %10.10f\n",k * mw_log(pobs));
    //mw_printf("Term 2 = %10.10f\n",(n-k) * mw_log(1 - pobs));

    return(mw_pow(10,result));
}


