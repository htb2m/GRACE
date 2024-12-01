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

bool _Gly_repetition(std::vector<char> sequence);
bool findGlyAtEveryThird(int seqLength, std::vector<char>& userSequence, int& firstGly);

// // // // // // // // // // // /// // // // // // // // IMPORTANT PARAMETERS // // // // // // // // // // // /// // // // // // // //


struct GA_Parameters {
    
    int maxHelices = 1250; // the array size of SCEPTTr library
    bool AAB = false; // if true, generate AAB ht instead of ABC ht
    
    // cannonical amino acids and hydroxyproline abbreviations
    vector<char> alphabet = {'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'};
  
    int alphabetSize = alphabet.size();
    
    int populationSize = 100;
    int numAA = 30;
    int numPep = 3;
    int GlyPos = 0;
    vector<vector<vector<char>>> Helices;
    
    double crossoverRate = 0;
    double mutationRate = 0.3;
    vector<double> MutationRateXaa;
    vector<double> MutationRateYaa;
    
    
    bool excludeXaa = false;
    bool excludeYaa = false;
    int numExcludedXaa = 0;
    int numExcludedYaa = 0;
    vector<char> excludedXaaList;
    vector<char> excludedYaaList;

    
    bool haveMotif = false;
    int motifLength = 0;
    int randomSeqLength = 15;
    vector<vector<char>> MotifSequences;

    

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
