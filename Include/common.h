#ifndef __common__
#define __common__

/*****************************************************************************************************************************************************************************************
******************************************************************************************************************************************************************************************
*
*
* This module contains all the constant used through the computation. All are in natural units : hbar = c = k = 1.
* Our energy unit is the MeV.
*
*****************************************************************************************************************************************************************************************
*****************************************************************************************************************************************************************************************/

#define  ALPHA    (1./137.)                          /* The fine structure constant */
#define  pi       (3.14159)
#define  m_e      (0.511)                            /* The mass of an electron */
#define  r_e      (ALPHA)/(m_e)                      /* The classical electron radius */
#define  a_0      r_e/((ALPHA)*(ALPHA))              /* The bohr radius */
#define  T_0      (2.7255*0.862*pow(10.,-10))        /* Temperature of the CMB today */
#define  sigma_T  8*(pi)*(r_e)*(r_e)/3.              /* Thompson cross-section */
#define  Y        0.25                               /* The primordial helium abundance */
#define  n_y_0    (3.154*pow(10.,-30))               /* The number density of CMB photons today */
#define  eta      (6.05*pow(10.,-10))                /* The baryon-to-photon ratio */
#define  n_e      (eta*n_y_0*(1-(Y)/(2)))                /* The number density of electrons today */
#define  E_x_0    (0.261121/(80*T_0))                /* Threshold above which photon-photon scattering starts to dominate today */
#define  E_c_0    (0.261121/(22*T_0))                /* Pair-production threshold today */
#define  H_0      (2.187*pow(10,-18))                /* Hubble rate today */
#define  H_r      (H_0*pow((1+7./8*pow(4./11,4./3)*3.046)*5.46*pow(10,-5),0.5))  /* Effective Hubble parameter today */


// #define Gamma_Table_Size    200                  /* Number of points in the gamma spectrum table */
// #define Electron_Table_Size 200                  /* Number of points in the electron spectrum table */
// #define E_min               1                    /* Minimal energy of the table, the maximum will be set by the mass of the decaying particle */

/*****************************************************************************************************************************************************************************************/
/************************************************************************Observationnal values of the nuclei abundances ******************************************************************/
#define omega_b  0.02225
#define tau_n  880.3
#define tau_n_0  880.3
#define Y_4He_0  (0.248*pow(omega_b/0.02273,0.39)*pow(tau_n/tau_n_0,0.72))
#define Y_4He_Min  0.2368
#define Y_4He_Max  1.
#define Y_3He_0  (1.02*pow(10,-5)*pow(omega_b/0.02273,-0.59)*pow(tau_n/tau_n_0,0.15))
#define Y_3He_Min  0.
#define Y_3He_Max  (1.5*pow(10,-5))
#define Y_2H_0  (2.53*pow(10,-5)*pow(omega_b/0.02273,-1.62)*pow(tau_n/tau_n_0,0.41))
#define Y_2H_Min  (2.56*pow(10,-5))
#define Y_2H_Max  (3.48*pow(10,-5))
#define Y_7Li_0  (1e-10)
#define Y_7Li_Min  0.
#define Y_7Li_Max  1.
#define Y_7Be_0  (1e-10)
#define Y_7Be_Min  0.
#define Y_7Be_Max  1.
/*****************************************************************************************************************************************************************************************/
/*****************************************************************************************************************************************************************************************/

// #define verbose 2                               /* Define the level of text printing during computation. 0 = no text, 1 = only main steps, 2 = maximal text.*/
#endif
