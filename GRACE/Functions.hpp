//
//  Functions.hpp
//  GRACE
//
//  Created by HA BUI from Hartgerink Lab on 10/18/24.
//


#ifndef Functions_hpp
#define Functions_hpp

#include <stdio.h>

#include "SCEPTTr.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <string>
#include <cstdlib>
#include <vector>

#pragma once

using namespace std;

bool Gly_repetition(std::vector<char> sequence);
bool findGlyAtEveryThird(int seqLength, std::vector<char>& userSequence, int& firstGly);

// // // // // // // // // // // /// // // // // // // // IMPORTANT PARAMETERS // // // // // // // // // // // /// // // // // // // //


struct GA_Parameters {
    
    int maxHelices; // the array size of SCEPTTr library
    bool AAB; // if true, generate AAB ht instead of ABC ht
    
    // cannonical amino acids  abbreviations
    vector<char> alphabet;
  
    int alphabetSize;
    
    int populationSize;
    int numAA;
    int numPep;
    int GlyPos;
    vector<vector<vector<char> > > Helices;
    
    double crossoverRate;
    double mutationRate;
    vector<double> MutationRateXaa;
    vector<double> MutationRateYaa;
    
    
    bool excludeXaa;
    bool excludeYaa;
    int numExcludedXaa;
    int numExcludedYaa;
    vector<char> excludedXaaList;
    vector<char> excludedYaaList;

    
    bool haveMotif;
    int motifLength;
    int randomSeqLength;
    vector<vector<char> > MotifSequences;

    

    // Constructor to initialize the mutation rate vectors
    GA_Parameters();
    
    // Set mutation rate of excluded amino acids to 0 (no mutation)
    void setXaaMutationRateToZero(char aminoAcid);
     
    void setYaaMutationRateToZero(char aminoAcid);
    


    
    
};


//
//struct GA_Parameters;

GA_Parameters initialPopulationGenerator(GA_Parameters GAparameters);

GA_Parameters CrossOver_withMotif(GA_Parameters parents);
GA_Parameters CrossOver (GA_Parameters Parents);

GA_Parameters Mutation_withMotif(GA_Parameters parents);
GA_Parameters Mutation (GA_Parameters parents);

vector<double> fitnessScore (GA_Parameters GAparameters,
                             TripleHelix * Library,
                             parameterType helixParameters);

GA_Parameters Selection (GA_Parameters Parents, const vector<double>& fitnessScores);

#endif /* Functions_hpp */
