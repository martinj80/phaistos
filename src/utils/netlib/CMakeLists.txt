# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(integrate
  slatec/dqags.f 
  slatec/dqagse.f 
  slatec/d1mach.f
  slatec/dqelg.f
  slatec/dqk21.f
  slatec/dqpsrt.f
  slatec/fdump.f
  slatec/i1mach.f
  slatec/j4save.f
  slatec/xercnt.f
  slatec/xerhlt.f
  slatec/xermsg.f
  slatec/xerprn.f
  slatec/xersve.f
  slatec/xgetua.f
  integrate.cpp
  integrate.h
)