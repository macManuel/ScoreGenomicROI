//
//  EpigeneticScoring.cpp
//  Filehandling
//
//  Created by Manuel Tuschen on 10.11.14.
//  Copyright (c) 2014 Manuel Tuschen. All rights reserved.
//

#include <fs_statistics.h>

#include "ScoreGenomicROI.h"


ScoreGenomicROI::ScoreGenomicROI(std::list<fs::GffFormat>& positions, std::list<fs::BedgraphFormat>& values, bool ignoreNegative , bool ignoreZero, bool countBins) {
  
  this->_positions = &positions;
  this->_values = &values;
  this->_ignoreNegative = ignoreNegative;
  this->_ignoreZero = ignoreZero;
  this->_countBins = countBins;
  
  this->_startChromosome = new std::map<std::string, std::list<fs::BedgraphFormat>::iterator>;
}

ScoreGenomicROI::~ScoreGenomicROI() {
  
  delete this->_startChromosome;
  this->_startChromosome = nullptr;
}


void ScoreGenomicROI::setIgnoreNegative(bool ignoreNegative) {
    this->_ignoreNegative = ignoreNegative;
}


void ScoreGenomicROI::setIgnoreZero(bool ignoreZero) {
    this->_ignoreZero = ignoreZero;
}



void ScoreGenomicROI::setHashList() {
  
  // sort lists first to be sure of the correct ordering
  this->_values->sort();
  this->_positions->sort();
  
  // hash the starting of a chromosome in the score file
  std::list<fs::BedgraphFormat>::iterator iteratorAtPosition = this->_values->begin();
  std::string oldChromosome = this->_values->begin()->chrom();
  
  this->_startChromosome->emplace(oldChromosome, iteratorAtPosition);
  
  for (auto it = this->_values->begin(); it != this->_values->end(); ++it) {
    std::string chomosome = it->chrom();
    
    if (oldChromosome != chomosome) {
      oldChromosome = chomosome;
      iteratorAtPosition = it;
      
      this->_startChromosome->emplace(oldChromosome, iteratorAtPosition);
    }
  }
}



void ScoreGenomicROI::scoreRegions() {
  
  // make sure the hash list is set
  this->setHashList();
  
  // first go through the genomicRegions to considere
  for (auto it = this->_positions->begin(); it != this->_positions->cend(); ++it) {
    
    // first get a frame to mark our region of interest
    unsigned int frameStart = it->start();
    unsigned int frameEnd = it->end();
    // Gff is a 1 based and closed intervall !
    unsigned int frameLength = frameEnd - frameStart + 1;
    unsigned int binCounter = 0;

    
    // It is not necessary to considere the strand as only the region matters
        
    // Now look in the value list for the first position matching the frame region
    // Choose correct iterator from hash list first
    double score = 0;
    unsigned int length = 0;

    std::list<fs::BedgraphFormat>::iterator jt = this->_startChromosome->at(it->seqname());
    
    // now iterate over all scored genomicRegions in the bedgraph file als long as those regions start not behind the frame
    while (jt->chromStart() < frameEnd) {
      
      // Be aware that the bedgraph format is zero-based hand the end value is not included in the intervall
      unsigned int binStart = jt->chromStart() + 1;
      unsigned int binEnd = jt->chromEnd();
      
      // if the frame is totally embraced by one bin of the bedgraph
      if (binStart <= frameStart &&  binEnd >= frameEnd) {
        
        // the score is the score * the frame size
        score = jt->score() * frameLength;
        jt++;
        continue;
      } else {
        
        // if there is a first partial overlap of the frame and the bin
        if (binStart < frameStart && binEnd >= frameStart) {
          if (this->_ignoreNegative == true && jt->score() > 0) {
            score += (binEnd - frameStart + 1) * jt->score();
          }
          if (this->_ignoreNegative == false) {
            score += (binEnd - frameStart + 1) * jt->score();
          }
          length += (binEnd - frameStart + 1);
          ++binCounter;
          jt++;
          continue;
        }
        
        // if a bin is completely embraced by the frame
        if (binStart >= frameStart && binEnd <= frameEnd) {
          if (this->_ignoreNegative == true && jt->score() > 0) {
            score += (binEnd - binStart + 1) * jt->score();
          }
          if (this->_ignoreNegative == false) {
            score += (binEnd - binStart + 1) * jt->score();
          }
          length += (binEnd - binStart + 1);
          ++binCounter;
          jt++;
          continue;
        }
        
        // if there is a partial overlapp at the end of the frame with a bin
        if (binStart <= frameEnd && binEnd > frameEnd) {
          if (this->_ignoreNegative == true && jt->score() > 0) {
            score += (frameEnd - binStart +1) * jt->score();
          }
          if (this->_ignoreNegative == false) {
            score += (frameEnd - binStart +1) * jt->score();
          }
          length += (frameEnd - binStart +1);
          ++binCounter;
          jt++;
          continue;
        }
        std::cout << "Error in scoring! Some unhandeld case occured. Please report a bug!" << std::endl;
        jt++;
      }
    }
    
    // perform a small test if the code works
    if (length > (frameLength)) {
      std::cout << "Error in scoring! Unvalid length calculated. Please report a bug!" << std::endl;
    }
    
    // save score
    it->setScore(score/ (double) frameLength);
    
    // save number of CpGs
    if (this->_countBins == true && it->hasAttribute_string() == true) {
      std::string tmp_str = it->attribute_string();
      it->setAttribute_string(tmp_str + "; CpGs "+ std::to_string(binCounter));
    } else {
      it->addAttribute("CpGs", std::to_string(binCounter));
    }
  }
}


void ScoreGenomicROI::boundDiscretization(std::vector<double> const & bounds) {
    // calculate number of classes to sort in
    unsigned long bins = bounds.size();

    // calculate the discrete scores
    for (auto it = this->_positions->begin(); it != this->_positions->end() ; ++it) {
        bool scored = false;
        
        // first check if the score is already zero and eventually need not to be taken into account.
        if (this->_ignoreZero == true && it->score() == 0) {
            scored = true;
        }
        // first check if the score is smaller zero and eventually need not to be taken into account.
        if (_ignoreNegative == true && it->score() < 0) {
            scored = true;
            it->setScore(0);
        }
        // set score according to border
        for (unsigned long i = 1; i <= bins; i++) {
            if (scored == false && it->score() <= bounds.at(i-1)) {
                it->setScore(i);
                scored = true;
            }
        }
        
        if (scored == false) {
            it->setScore(bins+1);
        }
    }
}


void ScoreGenomicROI::quantileDiscretization(unsigned int numb) {
  
  // set up a distribution from which the quantile borders can be determined.
  fs::Distribution dist;
  
  // get all scores and detemine the quantiles depending on the ignore conditions
  for (auto it = this->_positions->begin(); it != this->_positions->end(); ++it) {
    if (it->score() > 0) {
      dist.addElement(it->score());
    }
    if (this->_ignoreZero == false && it->score() == 0) {
      dist.addElement(it->score());
    }
    if (this->_ignoreNegative == false && it->score() < 0) {
      dist.addElement(it->score());
    }
  }
  

  
  std::vector<double> quantiles;
  // determine the quantiles
  for (unsigned int i = 1; i <= numb; ++i) {
      quantiles.push_back(dist.quantile((double) i / (double) numb));
  }
    
    this->boundDiscretization(quantiles);
}


void ScoreGenomicROI::intervallDiscretization(unsigned int numb) {
    std::vector<double> all_scores;
    
    // get all scores and detemine the intervall depending on the ignore conditions
    for (auto it = this->_positions->begin(); it != this->_positions->end(); ++it) {
        if (it->score() > 0) {
            all_scores.push_back(it->score());
        }
        if (this->_ignoreZero == false && it->score() == 0) {
            all_scores.push_back(it->score());
        }
        if (this->_ignoreNegative == false && it->score() < 0) {
            all_scores.push_back(it->score());
        }
    }
    
    std::sort(all_scores.begin(), all_scores.end());
    
    // determine the intervalls
    double min = all_scores.at(0);
    double max = all_scores.at(all_scores.size()-1);
    double intervall_size = (max - min) / numb;
    
    std::vector<double> intervalls;
    
    for (unsigned int i = 1; i <= numb; ++i) {
        intervalls.push_back(min + intervall_size * i);
    }
    this->boundDiscretization(intervalls);
}


