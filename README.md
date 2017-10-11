# Advanced Genomic Aligner (AGA)

AGA computes the optimal pairwise alignment of a nucleic acid sequence
against a reference genome, taking into account all CDS annotations of
the reference to simultaneously score nucleic acid and amino acid
similarity, minimize frameshifts, and avoid codon misaligned gaps.

It outputs the optimal nucleic acid sequence alignment and all CDS and
protein amino acid sequence alignments.

## Web version

A web version is available at
[http://www.genomedetective.com/app/aga](http://www.genomedetective.com/app/aga)

## Build

To compile AGA from source, you need a standard-compliant C++11
compiler, and CMake.

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

## Usage

As input, AGA requires a reference genome and a query sequence.

The reference genome is provided as a GENBANK flatfile record. For
example, go to [https://www.ncbi.nlm.nih.gov/nuccore/NC_001802.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001802.1), choose Send to:, File, Format Genbank.

The query sequence is provided in a FASTA file.

```
  aga {OPTIONS} [REFERENCE.GB] [QUERY.FASTA] [ALIGNMENT.FASTA]

    This is AGA, a Genomic Aligner, (c) Emweb bvba
    See http://github.com/emweb/aga/LICENSE.txt for terms of use.

  OPTIONS:

      --help                            Display this help menu
      --version                         Display the version
      Alignment mode, specify one of:
        --global                          Global alignment
        --local                           Local alignment
      Nucleic Acid Score options
        --nt-weight=[WEIGHT]              Weight for NT score fraction
                                          (default=1)
        --nt-gap-open=[COST]              Nucleotide Gap Open penalty
                                          (default=-10)
        --nt-gap-extend=[COST]            Nucleotide Gap Extension penalty
                                          (default=-1)
        --nt-match=[SCORE]                Score for a nucleotide match
                                          (default=2)
        --nt-mismatch=[COST]              Penalty for a nucleotide mismatch
                                          (default=-2)
      Amino Acid Score options
        --aa-weight=[WEIGHT]              Total weight for AA score fraction
                                          (default=1)
        --aa-gap-open=[COST]              Amino Acid Gap Open penalty
                                          (default=-6)
        --aa-gap-extend=[COST]            Amino Acid Gap Extension penalty
                                          (default=-2)
        --aa-matrix=[MATRIX]              Substitution matrix for amino acid
                                          matches: BLOSUM62 or BLOSUM30
                                          (default=BLOSUM30)
        --aa-frameshift=[COST]            Frameshift penalty (default=-100)
        --aa-misalign=[COST]              Codon misalignment penalty
                                          (default=-20)
      Amino acid alignments output
        --cds-aa-alignments=[ALIGNMENT.FASTA]
                                          Amino acid alignments output file of
                                          CDS (FASTA)
        --cds-nt-alignments=[ALIGNMENT.FASTA]
                                          Nucleic acid CDS alignments output
                                          file of CDS (FASTA)
        --protein-aa-alignments=[ALIGNMENT.FASTA]
                                          Amino acid alignments output file of
                                          Protein Products (FASTA)
        --protein-nt-alignments=[ALIGNMENT.FASTA]
                                          Nucleic acid CDS alignments output
                                          file of Protein Products (FASTA)
      REFERENCE.GB                      Annotated reference (Genbank Record)
      QUERY.FASTA                       FASTA file with nucleic acid query
                                        sequence
      ALIGNMENT.FASTA                   Nucleic acid alignment output file
                                        (FASTA)
      "--" can be used to terminate flag options and force all following
      arguments to be treated as positional options

    AGA will compute the optimal pairwise alignment of a nucleic acid query
    sequence (QUERY.FASTA) against a reference genome (REFERENCE.GB), taking
    into account CDS annotations in the genbank record to include in the
    alignment score all amino acid alignments and minimizing frameshifts within
    these open reading frames. It writes the resulting alignment to
    ALIGNMENT.FASTA
```

## Example

```
aga --global NC_001802.gb query.fasta alignment.fasta
```

## License

This project is licensed under the Emweb Non-Commercial Public
License. See the [LICENSE.txt](LICENSE.txt) file for details.
