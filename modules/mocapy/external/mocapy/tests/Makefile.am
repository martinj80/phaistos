AM_CXXFLAGS = -ggdb -Wall -I$(top_srcdir)/src $(BOOST_CPPFLAGS)
AM_LDFLAGS =  $(BOOST_LDFLAGS) $(LAPACK_LIBS) $(BLAS_LIBS) $(LIBS) 
LDADD = $(top_srcdir)/src/libMocapy.la $(BOOST_SERIALIZATION_LIB)

noinst_PROGRAMS = infenginehmm_test infenginemm_test test_mdarray_permute test_mdarray

infenginehmm_test_SOURCES = infenginehmm_test.cpp

infenginemm_test_SOURCES = infenginemm_test.cpp 

test_mdarray_permute_SOURCES = test_mdarray_permute.cpp 

test_mdarray_SOURCES = test_mdarray.cpp 
