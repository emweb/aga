/*
 * Copyright Emweb BVBA, 3020 Herent, Belgium
 *
 * See LICENSE.txt for terms of use.
 */

#include <limits>

#include "SimpleScorer.h"

AlignmentScoreVector::AlignmentScoreVector()
  : begin(std::numeric_limits<int>::max()),
    end(std::numeric_limits<int>::max())
{ }

std::ostream& operator<< (std::ostream& o, const AlignmentStats& stats)
{
  double alignLength = stats.matchCount + stats.insertCount + stats.deleteCount;

  if (alignLength == 0)
    return o << " N/A";
  
  o << "begin: " << (stats.begin + 1)
    << " end: " << stats.end
    << " coverage: " << 100.0 * stats.coverage / stats.refLength << "%"
    << " score: " << stats.score
    << " quality: " << (double)stats.score / stats.coverage
    << " matches: " << stats.matchCount
    << " (" << (stats.matchCount / alignLength * 100) << "%)"
    << " identities: " << stats.identityCount
    << " (" << (stats.identityCount / alignLength * 100) << "%)"
    << " inserts: " << stats.insertCount << " deletes: " << stats.deleteCount
    << " misaligned: " << stats.misaligned << " frameshifts: " << stats.frameShifts;

  return o;
}

void asJson(std::ostream& o, const std::string& id, const AlignmentStats& stats,
	    const std::string& mutationStr, const std::string& cds, int cdsBegin, int cdsEnd)
{
  double alignLength = stats.matchCount + stats.insertCount + stats.deleteCount;

  o << "{ \"id\" : \"" << id << "\", \"alignLength\" : " << alignLength << ", "
    << "\"cds\" : \"" << cds << "\", "
    << "\"cdsBegin\" : " << cdsBegin << ", "
    << "\"cdsEnd\" : " << cdsEnd;

  if (alignLength != 0) {
    o << ", "
      << "\"begin\" : " << (stats.begin + 1) << ", "
      << "\"end\" : " << stats.end << ", "
      << "\"coverage\" : " <<  100.0 * stats.coverage / stats.refLength << ", "
      << "\"score\" : " << stats.score << ", "
      << "\"quality\" : " << (double)stats.score / stats.coverage << ", "
      << "\"matches\" : " << stats.matchCount << ", "
      << "\"identities\" : " << stats.identityCount << ", "
      << "\"inserts\" : " << stats.insertCount << ", "
      << "\"deletes\" : " << stats.deleteCount << ", "
      << "\"misaligned\" : " << stats.misaligned << ", "
      << "\"frameshifts\" : " << stats.frameShifts << ", "
      << "\"ambiguities\" : " << stats.ambiguities << ", "
      << "\"stopCodons\" : " << stats.stopCodons << ", "
      << "\"mutations\" : \"" << mutationStr << "\"";
  }

  o << " }";
}

