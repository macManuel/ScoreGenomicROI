//
//  main.cpp
//  EpigeneticScoring
//
//  Created by Manuel Tuschen on 12.11.14.
//  Copyright (c) 2014 Manuel Tuschen. All rights reserved.
//

#include <iostream>

#include <fs_clap.h>
#include <fs_formats.h>
#include <fs_io.h>
#include <fs_exception.h>

#include "ScoreGenomicROI.h"

// This programm is designed to score a genomic region according to scores given in a .wig, .bedgraph, .bed file. Regions must be provided in a .gff, .gtf or .gff3 file.
//
//
// Parameters:
// # --regions <filename>:  A .gff file (version 1-3) including the path with those genomic regrions to extract from.
// # --scores <filename>: .bedgraph or .wig file containing a score for the whole genome.
// # --promoter <-bp +bp>: If this flag is set the start and stop values are set for the promoter size. -bp and +bp reflect the start and end position of the promoter relative to the TSS

int main(int argc, const char * argv[]) {    
  
  // First set up all the command line arguments to be parsed
  fs::CommandlineParser myParser("1.0", "Score a given genomic region.");
  
  fs::ClapSingleArgument<std::string> regionFile(myParser, "regions", "<filename>", "A .gff file (version 1-3) including the path with those genomic regrions to score.");
  fs::ClapSingleArgument<std::string> scoreFile(myParser, "scores", "<filename>", ".bedgraph or .wig file containing a score for the whole genome.");
  
  fs::ClapSingleArgument<bool> ignoreByDiscretization(myParser, 'd', "", "A flag indicating if the -n and -n flag should also be taken into account during discretization when given. This means, that when -d is set, only values ≥ or > are considered. All other scores were set to 0 and are therfore smaller than the first disrete class 1.", false);
  fs::ClapSingleArgument<bool> ignoreNegative(myParser, 'n', "", "A flag to indicate if negative values should be ignored during the scoring.", false);
  fs::ClapSingleArgument<bool> ignoreZero(myParser, 'z', "", "A flag to indicate if zeros should be ignored during the scoring.", false);
  fs::ClapSingleArgument<bool> countCpG(myParser, 'c', "", "A flag to indicate if the number of CpGs should be counted. Setting this flag makes only sense if one score per CpG is set in the score file!", false);
  
  fs::ClapSingleArgument<double> quantileDiscetize(myParser, "quantile", "<int> ", "One integer to determine in how many cases the scores should be classified in based on the quantile they are located instarting with 1.", false);
  fs::ClapSingleArgument<double> intervalDiscetize(myParser, "interval", "<int>", "One integer to determine in how many cases the scores should be classified in based on the intervall they are located in starting with 1.", false);
  fs::ClapMultiArgument<double> boundsDiscetize(myParser, "bounds", "<double> <double> ...", "A list of upper bounds for discetization. According to these bounds, the scores are classified into cases starting with 1.", false);

  // Now perform the parsing and extract the values
  myParser.parse(argc, argv);
  

  
  // As we now have everything, the program may start its computations
  
  // We need a list to store our gff objects in:
  std::list<fs::GffFormat> gff_list;
  
  // Read in the gff file:
  fs::GffReader gff_input(regionFile.value());
  gff_input.readFromFile(gff_list);

  // Next we need the scores
  std::list<fs::BedgraphFormat> bedgraph_list;

  try {
    fs::FileChecker checkSuffix;
    checkSuffix.check(scoreFile.value());
    
    if (checkSuffix.isBedgraph() == false && checkSuffix.isWig() == false ) {
      throw fs::Exception("No valid scoring file was given!");
    }
  } catch (fs::Exception &e) {
    std::cout << std::endl;
    std::cout << e.message() << std::endl;
    exit(1);
  }
  
  // Now score the genomic regions
  ScoreGenomicROI scoring(gff_list, bedgraph_list, ignoreNegative.value(), ignoreZero.value(), countCpG.value());
  scoring.scoreRegions();
  
  // write the output file
  fs::GffWriter gff_output(gff_input.path() + gff_input.file() + "_scored" + gff_input.suffix());
  gff_output.writeIntoFile(gff_list);
  
  // Next performe the disretizations
  // Therefore we need a copyof the scored regions
  std::list<fs::GffFormat> gff_scored = gff_list;

  // should the negative or zero values ignored by discretization
  bool ignore_zero = false;
  bool ignore_negative = false;
  
  if (ignoreByDiscretization .value()== true && ignoreZero.value() == true) {
    ignore_zero = true;
  }
  if (ignoreByDiscretization .value()== true && ignoreNegative.value() == true) {
    ignore_negative = true;
  }
  
  ScoreGenomicROI discretizing(gff_scored, bedgraph_list, ignore_negative, ignore_zero, countCpG.value());
  
  // last perform all disretizations necessary.
  if (boundsDiscetize.isSet() == true) {
    gff_scored = gff_list;
    discretizing.boundDiscretization(boundsDiscetize.value());
    
    // write the output file
    fs::GffWriter gff_output(gff_input.path() + gff_input.file() + "_scored_discrete_by_bounds" + gff_input.suffix());
    gff_output.writeIntoFile(gff_scored);
  }
  
  if (intervalDiscetize.isSet() == true) {
    gff_scored = gff_list;
    discretizing.intervallDiscretization(intervalDiscetize.value());
    
    // write the output file
    fs::GffWriter gff_output(gff_input.path() + gff_input.file() + "_scored_discrete_by_intervall" + gff_input.suffix());
    gff_output.writeIntoFile(gff_scored);
  }
  
  if (quantileDiscetize.isSet() == true) {
    gff_scored = gff_list;
    discretizing.intervallDiscretization(quantileDiscetize.value());
    
    // write the output file
    fs::GffWriter gff_output(gff_input.path() + gff_input.file() + "_scored_discrete_by_quantile" + gff_input.suffix());
    gff_output.writeIntoFile(gff_scored);
  }


  std::cout << "Finished!" << std::endl;
  return 0;
}
