// debugLevels.h --- Define different levels of debugging. Determines how assert works.
// Copyright (C) 2006-2008 Wouter Boomsma
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

#ifndef DEBUGLEVELS_H
#define DEBUGLEVELS_H

// DEFINE DEBUG LEVELS:
// 0: assert->No checks, strong_assert->No checks
// 1: assert->No checks, strong_assert->Warnings only
// 2: assert->Warnings only, strong_assert->Warnings only
// 3: assert->core dump, strong_assert->core dump
#ifndef DEBUGLEVEL
#define DEBUGLEVEL 3
#endif

#if DEBUGLEVEL == 0
#define NDEBUG
#include <cassert>
#undef assert
#define assert(x) ((void)0)
#define strong_assert assert
#elif DEBUGLEVEL == 1
#define NDEBUG
#include <cassert>
#undef assert
#define assert(x) ((void)0)
#define strong_assert(x)				\
if (! (x)) \
{ \
     std::cerr << "Warning: Assertion " << #x << " failed at " << __FILE__ << ":" << __LINE__ << "\n"; \
}
#elif DEBUGLEVEL == 2
#undef assert
#define assert(x)				\
if (! (x)) \
{ \
     std::cerr << "Warning: Assertion " << #x << " failed at " << __FILE__ << ":" << __LINE__ << "\n"; \
}
#define strong_assert assert
#else
// use normal behaviour: assertion -> core dump
#include <cassert>
#define strong_assert assert
#endif


#endif
