irfuproj_JUICE_noise
====================

Code descriptions:
====================
eo.noise.m
--------------------
Compute and plot the output voltage spectral densities of all the
noise components in the small signal integrated noise model, together with
the impedance of the QTN antenna model (Z_A), the IU-curve and the
associated resistances (R_s, R_ph).
   noise(n_e, T_e, V_bias) computes and plots as described above with the
   electron density n_e [cm^-3] and temperature T_e [eV], plus a bias
   potential of V_bias [V].
   
Calls on eo.qtnmod.m to evaluate the QTN voltage spectral density and antenna impedance
for the given plasma parameters n_e and T_e.
Calls on eo.oml.m and eo.grard.m to evaluate the current-voltage characteristics (IU-curve) of 
the probes for the given plasma parameters n_e and T_e.
Differentiates the probe currents at the operating point set by V_bias to obtain the
effective resistances between the probes (R_s and R_ph).
Computes and plots the contributions of the various noise sources to the output noise
voltage spectral density at the antenna terminal. This is done for each noise source by
filtering it through the potential divider made up by all the impedances and
resistances of the small signal integrated noise model, with all the other sources
turned off (superposition principle).

eo.qtn_plot_n_e.m
--------------------
This script plots the quasi-thermal noise (QTN) and associated antenna
impedance and capacitance for a number of electron densities at
constant temperature.

Calls on eo.qtnmod.m to evaluate the QTN voltage spectral densities and antenna impedances,
which are then plotted in a systematic way for the various electron densities.

eo.qtn_plot_T_e.m
--------------------
This script plots the quasi-thermal noise (QTN) and associated antenna
impedance and capacitance for a number of electron temperatures at
constant density.

Calls on eo.qtnmod.m to evaluate the QTN voltage spectral densities and antenna impedances,
which are then plotted in a systematic way for the various electron temperatures.

eo.qtnmod.m
--------------------
Compute quasi-thermal noise and associated antenna impedance for kappa-distributed electrons.
   qtnmod(n_e, T_e) computes the quasi-thermal noise in a double-sphere 
   dipole antenna in a plasma where the electron velocities obey a "Kappa-
   distribution". n_e is the electron density in cm^-3 and T_e is the
   electron temperature in eV.
   
Matlab's built-in numerical quadrature function quadgk.m is used to compute the integrals in the QTN and 
antenna impedance formulas for kappa-distributed plasmas.

eo.oml.m
--------------------
Computes current from spherical Langmuir probe due to single particle species under OML approximation.
   OML(n_j, T_j, v_Dj, {q_j, m_j, V_B, r_p}) computes and plots the particle
   current from given plasma parameters (from equations (7.16) - (7.21) in
   H?ymork et al) particle density n_j, temperature T_j, drift velocity
   v_Dj, charge q_j, mass m_j and probe radius r_p at bias voltage points
   in vector V_B.

eo.grard.m
--------------------
Computes photoelectron current vs voltage using expression in Grard
(1973).
   grard(V, V_p, J_ao, R, A) computes the photoelectron current as
   predicted by Grard (1973) for a spherical probe at the voltage points 
   contained in vector V, for photoelectron e-folding energy V_p, saturated
   current density J_ao (A/m^2 at 1 AU, material property), distance from the sun R (AU) 
   and sunward projected probe area A (m^2).

eo.debye.m
--------------------
Compute the Debye length (in meters) from n_e and T_e.

eo.plasmafreq.m
--------------------
Compute the plasma frequency (in Hertz) from n_e.

   
   
