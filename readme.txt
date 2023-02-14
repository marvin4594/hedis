Hedis solves the equations of fluiddynamics specialized to
axisymmetric, thin and isothermal systems. Additionally the
equation for disk evolution proposed by Pringle (1981) is
solved. The time integration is performed using an explicit
scheme 4th order or an implicit scheme 2nd order.
For a detailed description of the physics, equations and
the numerical schemes see the publication from M. Blank 2010
"Das Wachstum Schwarzer Löcher in aktiven galaktischen Kernen"

Hedis needs the software libraries BLAS and LAPACK. Type
"make PARAMS=noh-problem.f90" to compile the code with the
parameterfile for the Noh Problem, or choose another ex-
ample from the folder "examples". Then type "./hedis" to
run the code. The results will appear in the "results"
folder.
