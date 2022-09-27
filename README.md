# ConstraintSolverACMCSingleKerr
Constraintsolver for single Kerr-Black-Hole on ACMC-Slices.

This code solves the constraint equations for a single perturbed Kerr-Black-Hole as described in https://arxiv.org/abs/1301.6984. 
Afterwards the apparent horizon is searched. At the horizon the Bondi-Mass and angular momentum, as well as some inequalities are calculated.

# Compilation

Run _make_ in _src/_ to execute the Makefile.

# Configuration

- **Newton_itmin** the minimum number of steps in the Newton-Raphson-Method
- **Newton_itmax** the maximum number of steps in the Newton-Raphson-Method
- **Newton_tol** the satisfing residuum of the Newton-Raphson-Method
- **Newton_verb** sets the verbosity of the output, 1 for verbose, 0 for silent

- **which_data(Directory)** directory for the precalculated coefficients of the constraint equation
- **ns_domain_0** and **ns_domain_0** represents the number of grid-points in radial direction. Domain 0 is located at the horizon, domain 1 at Scri+.
- **nt** is the number of radial gridpoints in both domains

- **Goal_value_sig0** is the final value of the radial coordinate sigma for domain 0.
- **Goal_value_sig1** is the final value of the radial coordinate sigma for domain 1.

- **Number_of_Solutions_in_the_Sequence** sets the number of intermediate steps

- **Create_ID** if 1 creates new initial data for the Newton-Raphson-Method and exits. The new initial data is stored in _run/Created_InitialData_

# Basic Usage

For j=0.5 a set of Coefficient functions is provided within the repository.

Edit the _run/Config_ file, to your needs. For pertubations set a value of **Goal_value_sig1** unequal 1. Start the code with by executing _startACMC.sh_ in _run/_.


# Advanced Usage

To calculate pertubations for different starting j another set of precalculated coefficients has to be generated. 
The numerical precision of the data can be modified by editing the used datatype in _src/ftype.h.
The default ist _long double_ by setting _double_ by setting mpfr/boost datatype the precision can be drastically improved.

