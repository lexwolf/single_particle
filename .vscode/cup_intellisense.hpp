#pragma once

// Ensure CUP headers are parsed in the same backend mode used by this project.
#if !defined(CUP_BACKEND_QUASI_STATIC) && !defined(CUP_BACKEND_MIE)
#define CUP_BACKEND_QUASI_STATIC
#endif

#include <cup/cup.hpp>
