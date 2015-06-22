//
//  EpigeneticScoring.h
//  Filehandling
//
//  Created by Manuel Tuschen on 10.11.14.
//  Copyright (c) 2014 Manuel Tuschen. All rights reserved.
//

#ifndef Filehandling_EpigeneticScoring_h
#define Filehandling_EpigeneticScoring_h

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm> 

#include <fs_formats.h>

/**
 *  This class is specifically designed to score specific genetic regions
 *  according to one epigenetic feature. To do so first a list of genetic
 *  locations in with gff objectd is required. It must at least contain the
 *  chromosome in the seqname field and the start and end position. To score
 *  each region a list with scores in bedgraph/wig format is also required.
 *  The output is a gff file with the new score entries. Be aware that this 
 *  class works on the original data given!
 */
class ScoreGenomicROI {
protected:
  
  /**
   *  A boolean indicating if negarive values should be ignored during scoring.
   */
  bool _ignoreNegative;
  
  /**
   *  A boolean indicating if zero values should be ignored during the scoring.
   */
  bool _ignoreZero;
  
  /**
   *  A boolean indicating if the number of CpGs should be counted.
   */
  bool _countBins;
  
  /**
   *  A list pointer to gff objects which contain the regions to score.
   */
  std::list<fs::GffFormat>* _positions;
  
  /**
   *  A list pointer to bedgraph objects which contain the scores.
   */
  std::list<fs::BedgraphFormat>* _values;
  
  /**
   *  A hash map used to store the starting positions of each chromosomein the 
   *  bedgrah file. This makes it much faster as it it not necessay to go to 
   *  through the whole file.
   */
  std::map<std::string, std::list<fs::BedgraphFormat>::iterator>* _startChromosome;
  
  
  
protected:
  
  /**
   *  Create a hash list with iterators to specific chromosomes in the value list for faster processing.
   */
  void setHashList();
  
  
  
public:
  
  /**
   *  The constructor. It makes sure, that all values necessay are set.
   */
  ScoreGenomicROI(std::list<fs::GffFormat>& positions, std::list<fs::BedgraphFormat>& values, bool ignoreNegative = false, bool ignoreZero = false, bool countBins = false);
  
  /**
   *  The destructor.
   */
  virtual ~ScoreGenomicROI();
  
  
  /**
   *  Set the boolean if negative values should be ignored.
   *
   *  @param ignoreNegative True or False.
   */
  void setIgnoreNegative(bool ignoreNegative);
  
  /**
   *  Set the boolean if zeros should be ignored.
   *
   *  @param ignoreZero True or False.
   */
  void setIgnoreZero(bool ignoreZero);
  
  
  /**
   *  The main function to score each region.
   */
  void scoreRegions();
  
  /**
   *  Transform continous scores into discrete ones based on a vector of upper
   *  bounds.
   *
   *  @param bounds A vector containing the upper bound for one discrete bin.
   */
  void boundDiscretization(std::vector<double> const & bounds);
  
  /**
   *  Transform continous scores into discrete ones based on numb quantiles.
   *
   *  @param numb The number of quantiles.
   */
  void quantileDiscretization(unsigned int numb);
  
  /**
   *  Transform continous scores into discrete ones based on numb intervalls.
   *
   *  @param numb The number of intervalls.
   */
  void intervallDiscretization(unsigned int numb);
  
  
};

#endif


