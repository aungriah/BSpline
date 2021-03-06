cmake_minimum_required(VERSION 3.0.0)

project(BSplineCompiler)

# FITPACK FORTRAN code
enable_language(Fortran)

if(${CMAKE_Fortran_COMPILER} MATCHES "ifort.*")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -fltconsistency")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -fltconsistency")
endif()

set(fitpack_SRC
    fitpack/bispev.f
    fitpack/curfit.f
    fitpack/fpback.f
    fitpack/fpbisp.f
    fitpack/fpchec.f
    fitpack/fpcurf.f
    fitpack/fpbspl.f
    fitpack/fpdisc.f
    fitpack/fpgivs.f
    fitpack/fpknot.f
    fitpack/fporde.f
    fitpack/fprank.f
    fitpack/fprati.f
    fitpack/fprota.f
    fitpack/fpsurf.f
    fitpack/parder.f
    fitpack/splder.f
    fitpack/splev.f
    fitpack/surfit.f
    fitpack/parcur.f
    fitpack/fppara.f
    fitpack/curev.f 
    fitpack/cualde.f
    fitpack/fpader.f
)

add_library(fitpack STATIC 
            ${fitpack_SRC}
)

enable_language(CXX)

# C++ interface
include(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE ""
                         SYMBOLS curfit splev splder surfit bispev parder parcur curev cualde)
FortranCInterface_VERIFY(CXX)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(fitpackpp_SRC
	fitpackpp/BSplineCurve.cpp
    fitpackpp/BSplineSurface.cpp
    fitpackpp/BSpline.cpp)

add_library(fitpackpp STATIC
	        ${fitpackpp_SRC})

target_link_libraries(fitpackpp fitpack)