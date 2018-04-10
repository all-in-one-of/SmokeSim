# pragma once

# include <string>
# include <map>
# include <iostream>
# include <vector>
# include "Partio.h"
# include "Eigen/Dense"
# include "Eigen/Core"
# include <cmath>
# include <Eigen/IterativeLinearSolvers>

# define M_PI           3.14159265358979323846  /* pi */
# define M_2PI          6.28318530717958647692  /* 2pi */
# define M_hPI          1.57079632679489661923  /* pi / 2 */

// to switch between double and float
typedef float fReal;

enum gridType { SMOKE, SOLID };

// Linear interpolation
// TODO cubic interpolation
fReal Lerp(const fReal &fromEndPoint, const fReal &toEndPoint, fReal ratio)
{
    return (1.0 - ratio) * fromEndPoint + ratio * toEndPoint;
}