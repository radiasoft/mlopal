set (_SRCS
    OpalSample.cpp
    RNGStream.cpp
    SampleCmd.cpp
    Sampler.cpp
    )

include_directories (
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_opal_sources(${_SRCS})

set (HDRS
    FromFile.h
    LatinHyperCube.h
    Normal.h
    OpalSample.h
    RNGStream.h
    SampleCmd.h
    SampleGaussianSequence.h
    SampleIndividual.h
    SamplePilot.h
    Sampler.h
    SampleRandomizedSequence.h
    SampleSequence.h
    SampleWorker.h
    SamplingMethod.h
    Uniform.h
    )

install (FILES ${HDRS} DESTINATION "${CMAKE_INSTALL_PREFIX}/include/Sample")
