# verif_coupling
Verification code for the PoP paper "Simulating a pulsed power-driven plasma with ideal MHD"

They realize 1D problem with pulsed power generator driving a thin shell of material into internal magnetic field in RZ (_r) as well as cartesian (_z) geometry.
This alternatively known as flux compression (internal flux is compressed by external forces).

The \*.cpp files are "problem" files for Athena++ MHD code
The \*.txt files are "input" files for Athena++
The \*.py files are ODE solvers that realize the same problem
