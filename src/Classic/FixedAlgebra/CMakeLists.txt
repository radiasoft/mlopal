set (_SRCS
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources (${_SRCS})

set (HDRS
    FArray1D.h
    FArray2D.h
    FDoubleEigen.h
    FLieGenerator.h
    FLieGenerator.hpp
    FLUMatrix.h
    FMatrix.h
    FMonomial.h
    FSlice.h
    FTpsData.h
    FTps.h
    FTps.hpp
    FTpsMath.h
    FVector.h
    FVps.h
    FVps.hpp
    LinearFun.h
    LinearFun.hpp
    LinearMap.h
    LinearMap.hpp
    Taylor.h
    Taylor.hpp
    TransportFun.h
    TransportFun.hpp
    TransportMap.h
    TransportMap.hpp
    )

install (FILES ${HDRS}  DESTINATION "${CMAKE_INSTALL_PREFIX}/include/FixedAlgebra")
# Some source files included in header files
install (FILES ${_SRCS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/FixedAlgebra")
