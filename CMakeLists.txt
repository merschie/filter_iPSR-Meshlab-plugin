# Copyright 2019-2020, Collabora, Ltd.
# SPDX-License-Identifier: BSL-1.0


set(SOURCES filter_iPSR.cpp)

set(HEADERS
    filter_iPSR.h
    src/Allocator.h
    src/Array.h
    src/BinaryNode.h
    src/BSplineData.h
    src/CmdLineParser.h
    src/Factor.h
    src/FEMTree.h
    src/FunctionData.h
    src/Geometry.h
    src/kdtree.h
    src/LinearSolvers.h
    src/MarchingCubes.h
    src/MAT.h
    src/MyMiscellany.h
    src/Ply.h
    src/PointStream.h
    src/PointStreamData.h
    src/PoissonRecon.h
    src/Polynomial.h
    src/PPolynomial.h
    src/RegularTree.h
    src/SparseMatrix.h
    src/SparseMatrixInterface.h
    src/utility.h
    src/Window.h
    )

add_meshlab_plugin(filter_iPSR ${SOURCES} ${HEADERS})
