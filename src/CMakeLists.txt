SUBDIRS(libseq)

INCLUDE_DIRECTORIES(libseq)

SET(LIB_SOURCES
  Cigar.cpp Genbank.cpp Genome.cpp
  NTSequence6AA.cpp SimpleScorer.cpp SubstitutionMatrix.cpp
  SearchRange.cpp
)  

ADD_LIBRARY(agalib ${LIB_SOURCES})

ADD_EXECUTABLE(aga Aga.cpp)
TARGET_LINK_LIBRARIES(aga agalib seq)

INSTALL_TARGETS(/lib agalib)
INSTALL_TARGETS(/bin aga)
INSTALL_FILES(/include .*\\.h)

IF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/geval)
  SUBDIRS(geval)
ENDIF (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/geval)
