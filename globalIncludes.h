# pragma once

//# define PARTIO
# include <string>
# include <map>
# include <iostream>
# include <vector>
# ifdef PARTIO
# include "Partio.h"
# endif
# include "Eigen/Dense"
# include "Eigen/Core"
# include <cmath>
# include "Eigen/IterativeLinearSolvers"
//# include <unsupported/Eigen/IterativeSolvers>

# define M_PI           3.14159265358979323846  /* pi */
# define M_2PI          6.28318530717958647692  /* 2pi */
# define M_hPI          1.57079632679489661923  /* pi / 2 */

// to switch between double and float
// ConjugateGradient method requires double
typedef double fReal;

enum gridType { SMOKE, SOLID };