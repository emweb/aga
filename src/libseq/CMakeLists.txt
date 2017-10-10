SET(SOURCES
  Nucleotide.cpp AminoAcid.cpp NTSequence.cpp
  AASequence.cpp Codon.cpp CodingSequence.cpp
)

ADD_LIBRARY(seq ${SOURCES})

INSTALL_TARGETS(/lib seq)
INSTALL_FILES(/include .*\\.h)
