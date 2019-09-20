/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include "Genbank.h"
#include "Genome.h"

void expect(bool value, const std::string& msg)
{
  if (!value)
    throw std::runtime_error("Parse error, expecting " + msg);
}

std::string nextLine(std::istream& i)
{
  std::string line;
  do {
    std::getline(i, line);
    if (!i)
      return std::string();
  } while (line.empty());

  return line;
}

struct ParsedLine {
  std::string key;
  std::string value;
  int valueCol;  
};

ParsedLine parseKeyValue(const std::string& line)
{
  ParsedLine result;

  std::size_t pos = 0;

  while (pos < line.length() && std::isspace(line[pos]))
    ++pos;

  if (pos == line.length())
    throw std::runtime_error("Error parsing line: " + line + ", pos " + std::to_string(pos));
 
  while (pos < line.length() && !std::isspace(line[pos])) {
    result.key += line[pos];
    ++pos;
  }

  if (pos == line.length()) {
    result.valueCol = -1;
    return result;
  }

  while (pos < line.length() && std::isspace(line[pos]))
    ++pos;
  
  if (pos == line.length()) {
    result.valueCol = -1;
    return result;
  }

  result.valueCol = pos;
  result.value += line.substr(pos);

  return result;
}

bool isValueContinuation(const std::string& line, int valuePos)
{
  if (line.length() < valuePos)
    return false;

  for (int i = 0; i < valuePos; ++i) {
    if (!std::isspace(line[i]))
      return false;
  }

  return true;
}

std::istream& operator>>(std::istream& i, GenbankRecord& record)
{
  std::string line;

  for (;;) {
    if (line.empty())
      line = nextLine(i);

    if (!i)
      return i;

    expect(!std::isspace(line[0]), line + " does not start with a space");

    if (line == "//")
      return i;
    
    GenbankRecord::string_map reference;
    std::string *v = 0;
    
    ParsedLine entry = parseKeyValue(line);
    if (entry.key == "LOCUS")
      v = &record.locus;
    else if (entry.key == "DEFINITION")
      v = &record.definition;
    else if (entry.key == "ACCESSION")
      v = &record.accession;
    else if (entry.key == "VERSION")
      v = &record.version;
    else if (entry.key == "DBLINK")
      v = &record.dbLink;
    else if (entry.key == "KEYWORDS")
      v = &record.keywords;
    else if (entry.key == "SOURCE")
      v = &record.source;
    else if (entry.key == "COMMENT")
      v = &record.comment;
    else if (entry.key == "REFERENCE")
      v = &reference["KEY"];
    else if (entry.key == "FEATURES") {
      // see below
    } else if (entry.key == "ORIGIN") {
      // see below
    } else {
      std::cerr << "Warning: ignoring entry key: " << entry.key << std::endl;
    }

    line = std::string();

    if (v)
      *v = entry.value;

    GenbankRecord::Feature *feature = 0;
    
    for (;;) {
      line = nextLine(i);
      if (!i)
	return i;

      if (isValueContinuation(line, entry.valueCol)) {
	std::string value = line.substr(entry.valueCol);
	if (feature) {
	  if (value[0] == '/') {
	    std::size_t e = value.find('=');
	    std::string key = value.substr(1, e - 1);
	    v = &feature->qualifiers[key];
	    value = value.substr(e + 1);
	  }
	}
	if (v) {
	  if (!v->empty())
	    *v += "\n";
	  *v += value;
	}
      } else if (std::isspace(line[0])) {
	ParsedLine subEntry = parseKeyValue(line);
	if (entry.key == "SOURCE" && subEntry.key == "ORGANISM") {
	  v = &record.organism;
	  *v = subEntry.value;
	} else if (!reference.empty()) {
	  v = &reference[subEntry.key];
	  *v = subEntry.value;
	} else if (entry.key == "FEATURES") {
	  if (feature) {
	    for (auto& q : feature->qualifiers) {
	      if (q.second.length() >= 2) {
		q.second.erase(0, 1);
		q.second.erase(q.second.length() - 1, 1);
	      }
	    }
	  }
	  record.features.push_back(GenbankRecord::Feature());
	  feature = &record.features.back();
	  feature->type = subEntry.key;
	  v = &feature->location;
	  *v = subEntry.value;
	} else if (entry.key == "ORIGIN") {
	  std::string& seq = subEntry.value;
	  for (int i = 0; i < seq.length(); ++i)
	    if (!std::isspace(seq[i]))
	      record.sequence.push_back(seq::Nucleotide(seq[i]));
	} else {
	  std::cerr << "Warning: ignoring sub entry key: " << subEntry.key << " for " << entry.key << std::endl;
	  v = 0;
	}
      } else
	break;
    }

    if (!reference.empty())
      record.references.push_back(reference);
    else if (feature) {
      for (auto& q : feature->qualifiers) {
	if (q.second.length() >= 2) {
	  q.second.erase(0, 1);
	  q.second.erase(q.second.length() - 1, 1);
	}
      }
    }
  }

  return i;
}

std::string removeNewLines(const std::string& input)
{
  std::string result;

  for (std::size_t i = 0; i < input.size(); ++i) {
    if (input[i] == '\n' || input[i] == '\r')
      result += ' ';
    else
      result += input[i];
  }

  return result;
}

std::string makeValidId(const std::string& input)
{
  std::string result;

  for (std::size_t i = 0; i < input.size(); ++i) {
    if (input[i] == ' ')
      result += '_';
    else
      result += input[i];
  }

  return result;
}

Genome getGenome(const GenbankRecord& record)
{
  Genome::Geometry geometry = Genome::Geometry::Linear;
  if (record.locus.find("circular") != std::string::npos)
    geometry = Genome::Geometry::Circular;
  
  Genome result(record.sequence, geometry);
  result.setName(removeNewLines(record.version));
  result.setDescription(removeNewLines(record.definition));
  result.sampleAmbiguities();

  for (const auto& f : record.features) {
    if (f.type == "CDS") {
      std::string name, description;
      if (f.qualifiers.count("gene") > 0)
	name = f.qualifiers.at("gene");
      if (name.empty())
	if (f.qualifiers.count("locus_tag") > 0)	
	  name = f.qualifiers.at("locus_tag");
      if (f.qualifiers.count("product") > 0)
	if (name.empty())
	  name = f.qualifiers.at("product");
	else
	  description = f.qualifiers.at("product");
      if (name.empty())
      	if (f.qualifiers.count("note") > 0)
	  name = f.qualifiers.at("note");

      int i = result.cdsFeatures().size() + 1;

      name = makeValidId(removeNewLines(name));

      result.addCdsFeature(CdsFeature(std::to_string(i) + "_" + name, f.location, description));
    }
  }

  /*
  for (const auto& f : record.features) {
    if (f.type == "gene") {
      std::string name;
      if (f.qualifiers.count("gene") > 0)
	name = f.qualifiers.at("gene");
      if (name.empty())
	if (f.qualifiers.count("locus_tag") > 0)	
	  name = f.qualifiers.at("locus_tag");

      CdsFeature feat(name, f.location);

      bool contained = false;
      for (const auto& g : result.cdsFeatures()) {
	if (g.contains(feat)) {
	  contained = true;
	  break;
	}
      }

      if (!contained)
	result.addCdsFeature(feat);
    }
  }
  */

  return result;
}

std::vector<CdsFeature> getProteins(const Genome& genome, const GenbankRecord& record)
{
  std::vector<CdsFeature> result;

  for (const auto& f : record.features) {
    if (f.qualifiers.count("protein_id") > 0) {
      std::string name = f.qualifiers.at("protein_id");
      if (f.qualifiers.count("product") > 0)
	name += " " + f.qualifiers.at("product");
      else if (f.qualifiers.count("gene") > 0)	
	name += " " + f.qualifiers.at("gene");
      CdsFeature cdsFeature(name, f.location);
      if (genome.processCdsFeature(cdsFeature))
	result.push_back(cdsFeature);
    }
  }

  return result;
}
