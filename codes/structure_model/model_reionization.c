/*
 * model_reionization.c
 *
 *  Created on: Apr 10, 2013
 *      Author: luyu
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_sf_erf.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"

double do_reionization(double mvir, double Zcurr)
{
  double alpha, a, a0, ar, f_of_a, a_on_a0, a_on_ar, mass, Mfiltering, Mjeans, Mchar, mass_to_use, modifier;
  double Tvir, Vchar, deltacritZ, HubbleZ;

  /*  we use the fitting formulae given by Kravtsov et al. (2004) Appendix B,
   * the model is described in Gnedin (2000) */

  mass = mvir*xhubble*1e-10; // the rest of the code is written in units of 1e10Msun/h

  /*  here are two parameters that Kravtsov et al keep fixed. */
  /*  alpha gives the best fit to the Gnedin data */
  alpha = 6.0;
  Tvir = 1e4;

  /*  calculate the filtering mass */
  a = 1.0 / (1.0 + Zcurr);
  a0 = 1./(1.+Par.Reionization_z0);
  ar = 1./(1.+Par.Reionization_zr);
  a_on_a0 = a / a0;
  a_on_ar = a / ar;

  if(a <= a0)
    f_of_a = 3.0 * a / ((2.0 * alpha) * (5.0 + 2.0 * alpha)) * pow(a_on_a0, alpha);
  else if((a > a0) && (a < ar))
    f_of_a =
      (3.0 / a) * a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
      a * a / 10.0 - (a0 * a0 / 10.0) * (5.0 - 4.0 * pow(a_on_a0, -0.5));
  else
    f_of_a =
      (3.0 / a) * (a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
                   (ar * ar / 10.0) * (5.0 - 4.0 * pow(a_on_ar, -0.5)) - (a0 * a0 / 10.0) * (5.0 -
                                                                                             4.0 *
                                                                                             pow(a_on_a0,
                                                                                                 -0.5)) +
                   a * ar / 3.0 - (ar * ar / 3.0) * (3.0 - 2.0 * pow(a_on_ar, -0.5)));

  /*  this is in units of 10^10Msun/h, note mu=0.59 and mu^-1.5 = 2.21 */
  Mjeans = 25.0 * pow(Omega_M0, -0.5) * 2.21;
  Mfiltering = Mjeans * pow(f_of_a, 1.5);


  /*  calculate the characteristic mass coresponding to a halo temperature of 10^4K */
  Vchar = sqrt(Tvir / 36.0);
  //double omegaZ = Omega * (pow(1.0 + Zcurr, 3.0) / (Omega * pow(1.0 + Zcurr, 3.0) + Omega_Lambda0));
  //double xZ = omegaZ - 1.0;
  //deltacritZ = 18.0 * M_PI * M_PI + 82.0 * xZ - 39.0 * xZ * xZ;
  deltacritZ = Delta_vir(Zcurr); // use written function to replace above 3 lines
  HubbleZ = H0 * sqrt(Omega_M0 * pow(1.0 + Zcurr, 3.0) + Omega_Lambda0);

  Mchar = 1e-10 * Vchar * Vchar * Vchar / (G * HubbleZ * sqrt(0.5 * deltacritZ));


  /*  we use the maximum of Mfiltering and Mchar */
  mass_to_use = dmax(Mfiltering, Mchar);
  modifier = 1.0 / pow(1.0 + 0.26 * (mass_to_use / mass), 3.0);

  //printf("%g %g %g %g\n", mass, mass_to_use, Mfiltering, Mchar);
  return modifier;

}
