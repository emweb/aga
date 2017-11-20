// This may look like C code, but it's really -*- C++ -*-
/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */
#ifndef GENBANK_H_
#define GENBANK_H_

#include <map>
#include <string>
#include <vector>
#include "NTSequence.h"

class Genome;

struct GenbankRecord
{
  typedef std::map<std::string, std::string> string_map;
  struct Feature {
    std::string type;
    std::string location;
    string_map qualifiers;
  };

  std::string locus, definition, accession, version, dbLink, keywords, source, organism, comment;
  seq::NTSequence sequence;
  std::vector<string_map> references;
  std::vector<Feature> features;
};

class CdsFeature;

extern std::istream& operator>>(std::istream& s, GenbankRecord& record);

extern Genome getGenome(const GenbankRecord& record);

extern std::string removeNewLines(const std::string& input);
extern std::string makeValidId(const std::string& input);

extern std::vector<CdsFeature> getProteins(const Genome& genome, const GenbankRecord& record);

#endif // GENBANK_H_
