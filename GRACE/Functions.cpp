//
//  Functions.cpp
//  GRACE
//
//  Created by Hartgerink Lab on 11/30/24.
//

#include "Functions.hpp"
#include "SCEPTTr.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <unordered_set>
#include <algorithm> // For std::shuffle
#include <random>



using namespace std;


    
// GA_Parameters constructor to initialize parameters
GA_Parameters::GA_Parameters() {
    // Initialize alphabet with canonical amino acids and hydroxyproline abbreviations
    alphabet = {'a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y'};
    alphabetSize = alphabet.size();
    
    maxHelices = 1250; // the array size of SCEPTTr library
   
    AAB = false; // if true, generate AAB ht instead of ABC ht
    
    populationSize = 500;
    numAA = 30;
    numPep = 3;
    GlyPos = 0;

    crossoverRate = 0.6;
    mutationRate = 0.2;
    
    excludeXaa = false;
    excludeYaa = false;
    numExcludedXaa = 0;
    numExcludedYaa = 0;
    
    haveMotif = false;
    motifLength = 0;
    randomSeqLength = 15;
    
    // Initialize alphabetSize based on the alphabet vector size
    alphabetSize = alphabet.size();

    MutationRateXaa.resize(alphabetSize, mutationRate);
    MutationRateYaa.resize(alphabetSize, mutationRate);
}

// GA_Parameters function to turn Mutation rate of exluded amino acids to 0
void GA_Parameters::setXaaMutationRateToZero(char aminoAcid) {
    auto it = std::find(alphabet.begin(), alphabet.end(), aminoAcid);
    if (it != alphabet.end()) {
        int index = std::distance(alphabet.begin(), it);
        MutationRateXaa[index] = 0;
    }
}

void GA_Parameters::setYaaMutationRateToZero(char aminoAcid) {
    auto it = std::find(alphabet.begin(), alphabet.end(), aminoAcid);
    if (it != alphabet.end()) {
        int index = std::distance(alphabet.begin(), it);
        MutationRateYaa[index] = 0;
    }
}

bool Gly_repetition(vector<char> sequence) {
    short x, xCount, yCount, zCount;
    bool goodPeptide;
    
    xCount = 0;
    yCount = 0;
    zCount = 0;
    goodPeptide = false;
    
    int numAA = sequence.size();
    for (x=0; x < numAA; x++) {
        if ( (x%3 == 0) && ( (sequence[x] == 'G') or (sequence[x] == 'g') ) ) xCount++;
        if ( (x%3 == 1) && ( (sequence[x] == 'G') or (sequence[x] == 'g') ) ) yCount++;
        if ( (x%3 == 2) && ( (sequence[x] == 'G') or (sequence[x] == 'g') ) ) zCount++;
    }
    
    if (xCount >= (numAA/3 - 2)) {goodPeptide = true;}
    if (yCount >= (numAA/3 - 2)) {goodPeptide = true;}
    if (zCount >= (numAA/3 - 2)) {goodPeptide = true;}
    
    return goodPeptide;
}

// A global functions to ensure Gly at every third positions and return the position of first Gly
bool findGlyAtEveryThird(int seqLength, vector<char>& userSequence, int& firstGly) {
   
    
    bool GlyAtEveryThird = false;
    firstGly = -1;
    
   
    for (int i = 0; i < 3; ++i) {

        if (tolower(userSequence[i]) == 'g') {
            firstGly = i; // find first Gly
        } else {
            continue;
        }
        
        if (seqLength == 3 and firstGly != -1) {
            break;
        } else if (seqLength > 3) {
            bool correctGlyPos = true;
            for (int j = firstGly; j < seqLength; j+= 3) {
                if (tolower(userSequence[j]) != 'g') {
                    correctGlyPos = false;
                    break;
                }
            }
            if (!correctGlyPos) {
                continue;
            } else {
                GlyAtEveryThird = true;
                break;
            }
        }
    }
    
    
    return GlyAtEveryThird;
    
}


// // // // // // // // // // // /// // // // // // // // INITIAL POPULATION // // // // // // // // // // // /// // // // // // // //


GA_Parameters initialPopulationGenerator(GA_Parameters GAparameters) {
    
    vector<char> alphabet = GAparameters.alphabet; // vector to store the amino acid abbreviations

    short randomIndex = -1; // Initialize randomIndex of alphabet
    vector<vector<char>> helix; // Create a vector of vectors to store the sequences
    helix.reserve(GAparameters.numPep);
    vector<char> sequence; // Create a vector to store amino acids
    sequence.reserve(GAparameters.numAA);

    GAparameters.Helices.clear();
    

    for (int helixCount = 0; helixCount < GAparameters.populationSize; helixCount++) {

        helix.clear(); // Reinitialize the sequence vector
        helix.resize(GAparameters.numPep);
        sequence.clear();
        sequence.resize(GAparameters.numAA);

    
        if (GAparameters.haveMotif) {
            // -----random-seq-1-----user-seq--------random-seq-2-----
            // -----15-residues-----user-length------15-residues------
            
          
            
            GAparameters.numAA = GAparameters.randomSeqLength*2 + GAparameters.motifLength;
            // Generate 2 or 3 random sequences
            helix.resize(GAparameters.numPep);
            
            for (int i = 0; i < GAparameters.numPep; i++) {
                sequence.clear();
                sequence.resize(GAparameters.numAA);
                
                // generate the first 12 random residues
                for (int j = 0; j < GAparameters.randomSeqLength; j++) {
                    // Generate a random index within the alphabet vector
                    randomIndex = rand() % alphabet.size();
                    // Assign the randomly selected letter to the sequence
                    sequence[j]= alphabet[randomIndex];
                }
                
                
                // insert user sequences
                for (int j = 0; j < GAparameters.motifLength; j++) {
                    sequence[j+GAparameters.randomSeqLength] = GAparameters.MotifSequences[i][j];
                }
                
                // generate random sequence 2
                for (int j = (GAparameters.randomSeqLength + GAparameters.motifLength); j < GAparameters.numAA; j++) {
                    // Generate a random index within the alphabet vector
                    randomIndex = rand() % alphabet.size();
                    // Assign the randomly selected letter to the sequence
                    sequence[j]= alphabet[randomIndex];
                   
                }
            
                
                helix[i] = sequence;
            }
            
            
            // Insert Gly to every third position
                for (int i = 0; i < GAparameters.numPep; i++) {
                    for (int j = GAparameters.GlyPos; j < GAparameters.numAA; j+=3) {
                    //cout << " j " << j << " ";
                    helix[i][j] = 'G';
                }
              
            }
            
            
            // Delete excluded amino acids
            int Xaa = -1 , Yaa = -1;
                int GlyPos = GAparameters.GlyPos;
            if (GlyPos == 0) {Xaa = 1; Yaa = 2;}
            if (GlyPos == 1) {Xaa = 2; Yaa = 0;}
            if (GlyPos == 2) {Xaa = 0; Yaa = 1;}
            
            
                if (GAparameters.excludeXaa) {
                for (int i = 0; i < GAparameters.numPep; i++) {
                    for (int j = Xaa; j < GAparameters.numAA; j+=3) {
                        if ((j < GAparameters.randomSeqLength) or
                            (j > (GAparameters.randomSeqLength + GAparameters.motifLength))) {
                            while (find(GAparameters.excludedXaaList.begin(), GAparameters.excludedXaaList.end(), helix[i][j]) != GAparameters.excludedXaaList.end()) {
                                helix[i][j] = alphabet[rand() % alphabet.size()];
                            }
                        }
                    }
                }
            }
            
                if (GAparameters.excludeYaa) {
                    for (int i = 0; i < GAparameters.numPep; i++) {
                        for (int j = Yaa; j < GAparameters.numAA; j+=3) {
                            if ((j < GAparameters.randomSeqLength) or
                                (j > (GAparameters.randomSeqLength + GAparameters.motifLength))) {
                                while (find(GAparameters.excludedYaaList.begin(), GAparameters.excludedYaaList.end(), helix[i][j]) != GAparameters.excludedYaaList.end()) {
                                helix[i][j] = alphabet[rand() % alphabet.size()];
                            }
                        }
                    }
                }
            }
            
            
            
            // Check if Gly at everythird position
            bool GlyCheck = true;
                for (int i=0; i< GAparameters.numPep; i++) {
                GlyCheck = Gly_repetition(helix[i]);
                if (!GlyCheck) break;
            }
            if (GlyCheck) {
                GAparameters.Helices.push_back(helix);
            }
            
        }
        

    else { // IF NO USER MOTIF

            // GENERATE 3 (if ABC) or 2 (if AAB) RANDOM SEQUENCES
            for (int i = 0; i < GAparameters.numPep; i++ ) {
                sequence.clear();
                sequence.resize(GAparameters.numAA);
               
                for (int j = 0; j < GAparameters.numAA; j++) {
                   // Generate a random index within the alphabet vector
                    randomIndex = rand() % alphabet.size();
                    // Assign the randomly selected letter to the sequence
                    sequence[j] = alphabet[randomIndex];
                    
                }
                helix[i] = sequence;
            }
            
            
            // Randomly generate frameshift, GlyPos can be 0 (Gxy), 1 (yGx), or 2 (xyG)
            //GAparameters.GlyPos = rand() % 3;
            GAparameters.GlyPos = 0;

            for (int i = 0; i < GAparameters.numPep; i++) { // INSERT Gly to EVERY THIRD POS
                for (int j = GAparameters.GlyPos; j < GAparameters.numAA; j+=3) {
                    helix[i][j] = 'G';
                }
            }

            // DELETE EXCLUDED AMINO AICDS IF ANY
            int Xaa = -1 , Yaa = -1;

            if (GAparameters.GlyPos == 0) {Xaa = 1; Yaa = 2;}
            if (GAparameters.GlyPos == 1) {Xaa = 2; Yaa = 0;}
            if (GAparameters.GlyPos == 2) {Xaa = 0; Yaa = 1;}

            if (GAparameters.excludeXaa) {
                //cout << "XAA " << endl;
                for (int j = 0; j < helix.size(); j++) {
                    for (int i = Xaa; i < GAparameters.numAA; i += 3) {
                        while (find(GAparameters.excludedXaaList.begin(), GAparameters.excludedXaaList.end(), helix[j][i]) != GAparameters.excludedXaaList.end()) {
                            helix[j][i] = alphabet[rand() % alphabet.size()];
                        }
                    }
                }
            }


            if (GAparameters.excludeYaa) {
                //cout << "YAA" << endl;
                for (int j = 0; j < helix.size(); j++) {
                    for (int i = Yaa; i < GAparameters.numAA; i += 3) {
                        while (find(GAparameters.excludedYaaList.begin(), GAparameters.excludedYaaList.end(), helix[j][i]) != GAparameters.excludedYaaList.end()) {
                            helix[j][i] = alphabet[rand() % alphabet.size()];
                        }
                    }
                }
            }



            // CHECK IF Gly AT EVERY THIRD POS
            bool GlyCheck = true;
            for (int i=0; i<GAparameters.numPep; i++) {
                GlyCheck = Gly_repetition(helix[i]);
                if (!GlyCheck) break;
            }
            if (GlyCheck) {
                GAparameters.Helices.push_back(helix);
                //GAparameters.GlyPos.push_back(GlyPos);
            }


        }
    }

    return GAparameters;
}

// // // // // // // // // // // /// // // // // // // // CROSSOVER // // // // // // // // // // // /// // // // // // // //

// Function to find indices of 'G' at every third position starting from the first 'G'
vector<int> findGIndices(const vector<char>& randomSeq) {
    vector<int> indices;

    // Find the first occurrence of 'G' within the first 3 residues
    size_t start = randomSeq.size();
    for (size_t i = 0; i < 3; ++i) {
        if (randomSeq[i] == 'G') {
            start = i;
            break;
        }
    }

    // If no 'G' is found, return an empty vector and issue a warning
    if (start == randomSeq.size()) {
        cerr << "No 'G' found in the first 3 residues." << endl;
        return indices;
    }

    // Check every third position starting from the first 'G'
    for (size_t i = start; i < randomSeq.size(); i += 3) {
        if (randomSeq[i] == 'G') {
            indices.push_back(i);
        }
    }

    return indices;
}


// Function to process vectors and return indices of randomly selected G in left and right
pair<int, int> chooseCrossingPoints (const vector<char>& left, const vector<char>& right) {
    // Find G indices in both vectors
    vector<int> leftGIndices = findGIndices(left);
    vector<int> rightGIndices = findGIndices(right);

    // Randomly select a G from left
    if (!leftGIndices.empty()) {
        int randomIndex = rand() % leftGIndices.size(); // Randomly pick an index
        int selectedGIndexLeft = leftGIndices[randomIndex];

        // Find corresponding index in right for the same G count
        int selectedGIndexRight = -1;
        if (randomIndex < rightGIndices.size()) {
            selectedGIndexRight = rightGIndices[randomIndex];
        }
        // Return the indices
        return {selectedGIndexLeft, selectedGIndexRight};
    }

    // If no G found in left, return -1 for both
    return {-1, -1};
}

// Function to concatenate to vector <char> into one
vector<char> concatenate(const vector<char>& vec1, const vector<char>& vec2) {
    vector<char> result(vec1);
    result.insert(result.end(), vec2.begin(), vec2.end());
    return result;
}


// Converts a helix (vector of vector of chars) into a hashable string
string hashHelix(const std::vector<std::vector<char>>& helix) {
    ostringstream oss;
    for (const auto& sequence : helix) {
        for (char c : sequence) {
            oss << c; // Add each character
        }
        oss << '|'; // Use a delimiter to separate sequences
    }
    return oss.str(); // Return as a single string
}


GA_Parameters CrossOver_withMotif(GA_Parameters parents) {
    GA_Parameters Offsprings = parents;
    short numPep = parents.numPep;

    // Generate a random number from 0 to 1
    double randomNumber = static_cast<double>(rand()) / RAND_MAX;

    if (randomNumber < parents.crossoverRate) {
        int popSize = parents.populationSize;
        int randomSeqLen = parents.randomSeqLength;
        int motifLen = parents.motifLength;

        // Pre-reserve memory
        vector<vector<char>> randomSeqLeft, randomSeqRight;
        randomSeqLeft.reserve(popSize * numPep);
        randomSeqRight.reserve(popSize * numPep);

        // Collect random sequences on the left and right of the motif
        for (int i = 0; i < popSize; i++) {
            for (int j = 0; j < numPep; j++) {
                const std::vector<char>& seq = parents.Helices[i][j];
                randomSeqLeft.emplace_back(seq.begin(), seq.begin() + randomSeqLen);
                randomSeqRight.emplace_back(seq.begin() + randomSeqLen + motifLen, seq.end());
            }
        }

        // Generate crossover points
        auto crossingPoints = chooseCrossingPoints(randomSeqLeft[0], randomSeqRight[0]);

        // Perform crossover for left and right sequences
        auto performCrossover = [&](const std::vector<std::vector<char>>& seqs, int crossingPoint) {
            vector<vector<char>> crossedSeqs;
            crossedSeqs.reserve(seqs.size() * seqs.size());
            for (const auto& seq1 : seqs) {
                for (const auto& seq2 : seqs) {
                    if (seq1 == seq2) continue; // Skip identical pairs
                    crossedSeqs.emplace_back(concatenate(vector<char>(seq1.begin(), seq1.begin() + crossingPoint + 1),
                                                          vector<char>(seq2.begin() + crossingPoint + 1, seq2.end())));
                }
            }
            return crossedSeqs;
        };

        auto newRandomSeqLeft = performCrossover(randomSeqLeft, crossingPoints.first);
        auto newRandomSeqRight = performCrossover(randomSeqRight, crossingPoints.second);

        // Assemble new helices directly while avoiding full combination generation
        vector<vector<vector<char>>> newHelices;
        unordered_set<string> helixSet; // For uniqueness check

        for (size_t i = 0; i < newRandomSeqLeft.size() && newHelices.size() < parents.maxHelices; ++i) {
            for (size_t j = 0; j < newRandomSeqRight.size(); ++j) {
                vector<vector<char>> newHelix;
                bool validHelix = true;
                for (int k = 0; k < numPep; ++k) {
                    auto leftSeq = newRandomSeqLeft[i];
                    auto rightSeq = newRandomSeqRight[j];
                    auto motif = parents.MotifSequences[k];

                    auto helix = concatenate(concatenate(leftSeq, motif), rightSeq);

                    if (!Gly_repetition(helix)) {
                        validHelix = false;
                        break;
                    }
                    newHelix.push_back(helix);
                }

                if (validHelix) {
                    string helixHash = hashHelix(newHelix); // Hash for uniqueness
                    if (helixSet.find(helixHash) == helixSet.end()) {
                        helixSet.insert(helixHash);
                        newHelices.push_back(newHelix);
                    }
                }
            }
        }

        // Add new helices to offspring struct
        int totalExistingHelices = parents.Helices.size();
        int newHelicesToAdd = std::min(static_cast<int>(newHelices.size()), (parents.maxHelices - totalExistingHelices));
        Offsprings.Helices.insert(Offsprings.Helices.end(), newHelices.begin(), newHelices.begin() + newHelicesToAdd);
        Offsprings.populationSize = Offsprings.Helices.size();
    }
    
    return Offsprings;
}


GA_Parameters CrossOver (GA_Parameters Parents) {
    short i, j, k;
    GA_Parameters Offsprings = Parents;
    short numPep = Parents.numPep;
    

    vector<vector<vector<char> > > totalOffsprings;
    vector <int> OffspringsGlyPos;
    
    // generate a random number from 0 to 1
    double randomNumber = static_cast<double>(rand()) / RAND_MAX;
    
    if (randomNumber < Parents.crossoverRate) {
        
        // Select a random crossoverPoint
        int crossoverPoint = rand() % Parents.numAA;  // pick a random crossover point
        

        vector<vector<vector<char>>> totalFragment1;
        vector<vector<vector<char>>> totalFragment2;
        
        vector<vector<char>> P1_Fragment1;
        vector<vector<char>> P1_Fragment2;
        
        vector<vector<char>> P2_Fragment1;
        vector<vector<char>> P2_Fragment2;
        
        for (i = 0; i < Parents.populationSize; i++) {
            vector<vector<char>> Fragment1;
            vector<vector<char>> Fragment2;
            for (j=0; j<numPep; j++) {
                vector <char> frag1;
                
                for (k=0; k <= crossoverPoint; k++) {
                    frag1.push_back(Parents.Helices[i][j][k]);
                }
             
                Fragment1.push_back(frag1);
                frag1.clear();
                vector <char> frag2;
                
                for (k = (crossoverPoint + 1); k < Parents.numAA; k++) {
                    frag2.push_back(Parents.Helices[i][j][k]);
                }
            
                Fragment2.push_back(frag2);
                frag2.clear();
            }
            
            totalFragment1.push_back(Fragment1);
            totalFragment2.push_back(Fragment2);
            Fragment1.clear();
            Fragment2.clear();
        }
        
        
       
        
//        Generate all posible offspring combination
//        possibles offspring = {F1,F2} {F2,F1} {F1,F4} {F4,F1} {F3,F4} {F4,F3} {F3,F2} {F2,F3}
        vector<vector<char>> offspring1;
        vector<vector<char>> offspring2;
        for (i = 0; i < totalFragment1.size(); i++) { // fragment length should be 2
            for (j = 0; j < totalFragment2.size(); j++) {
                vector <char> strand1;
                vector <char> strand2;
                
                for (k = 0; k < Parents.numPep; k ++) {
                    strand1 = totalFragment1[i][k];
                    strand1.insert(strand1.end(), totalFragment2[j][k].begin(), totalFragment2[j][k].end());
                    offspring1.push_back(strand1);
                    strand1.clear();
                    
                    
                    strand2 = totalFragment2[i][k];
                    strand2.insert(strand2.end(), totalFragment1[j][k].begin(), totalFragment1[j][k].end());
                    offspring2.push_back(strand2);
                    strand2.clear();
                }
                totalOffsprings.push_back(offspring1);
                offspring1.clear();
                totalOffsprings.push_back(offspring2);
                offspring2.clear();
            }
        }
        
        
        // Check if there is Gly at every third position and correct if there is not
        for (i = 0; i < totalOffsprings.size(); i++) {
            int GlyPos = -1;
            for (j = 0; j < totalOffsprings[i].size(); j++) {
               // get Gly position
                for (k = 0; k < 2; k++) { // looping through first 3 position 0,1,2
                    if (totalOffsprings[i][j][k] == 'G' or totalOffsprings[i][j][k] == 'g') {
                        GlyPos = k;
                        break; // Found a Gly, stop here
                    }
                }
            }
            
            if (GlyPos < 0) GlyPos = 0; // if cannot find G in the first 3 residues, assign the GPO frame
            
            for (j = 0; j< totalOffsprings[i].size(); j++) {
                for (k = GlyPos; k < totalOffsprings[i][j].size(); k += 3) {
                    totalOffsprings[i][j][k] = 'G';
                }

            }
                
            
        }
  
        
  
        
            Offsprings.populationSize = totalOffsprings.size();
            Offsprings.Helices = totalOffsprings;
    }
        
    
    return Offsprings;
}


// // // // // // // // // // // /// // // // // // // // MUTATION  // // // // // // // // // // // /// // // // // // // //

vector<char> mutateSequence(vector<char> sequence, GA_Parameters parents) {
    vector<char> mutatedSequence;
    double randomNumber;
    int randomIndex;
    int numAA = sequence.size();
    int k;
    vector<char> alphabet = parents.alphabet;
    int Xaa = -1 , Yaa = -1;
    
    int GlyPos = -1;
    
    findGlyAtEveryThird(numAA, sequence, GlyPos);
    
    if (GlyPos == 0) {Xaa = 1; Yaa = 2;}
    if (GlyPos == 1) {Xaa = 2; Yaa = 0;}
    if (GlyPos == 2) {Xaa = 0; Yaa = 1;}
    
    
    // Mutate XaaPos
    for (k=Xaa; k < numAA; k+=3) {
        // Generate a random number between 0 and 1 to compare with the mutation rate
        randomNumber = static_cast<double>(rand()) / RAND_MAX;
        // Generate random index of the alphabet to mutate
        randomIndex = rand() % alphabet.size();
        if (randomNumber < parents.MutationRateXaa[randomIndex]) {
            //parents.Helices[helicesIndex][sequenceIndex][k] = alphabet[randomIndex];}
            sequence[k] = alphabet[randomIndex];}
    }
    // Mutate YaaPos
    for (k=Yaa; k < numAA; k+=3) {
        // Generate a random number between 0 and 1 to compare with the mutation rate
        randomNumber = static_cast<double>(rand()) / RAND_MAX;
        // Generate random index of the alphabet to mutate
        randomIndex = rand() % alphabet.size();
        if (randomNumber < parents.MutationRateYaa[randomIndex]) {
            sequence[k] = alphabet[randomIndex];}
    }
        
    mutatedSequence = sequence;
    
    return mutatedSequence;
}


// Function to check whether the chosen helix has already presented in the population before adding it to the population.
bool isUnique(const vector<vector<char>>& offspring, const vector<vector<vector<char>>>& existingHelices, int& indexOfOffspring) {
    
    for (int helix = 0; helix < existingHelices.size(); helix++) {
        if (helix == indexOfOffspring) continue; // Skip the helix's index
        
        bool isSame = true;
        for (int strand = 0; strand < existingHelices[helix].size(); strand++) {
            if (offspring[strand] != existingHelices[helix][strand]) {
                isSame = false;
                break;
            }
        }
        
        if (isSame) return false;
    }
    return true;
}

// Generate new combination of 3 sequences
void generateCombinations(vector<vector<char>>& mutatedAndOriginal, int start, vector<vector<char>>& current, vector<vector<vector<char>>>& result, short numPep) {
    // Add the current combination to the result
    if (current.size() == numPep) result.push_back(current);

    // Iterate through the remaining elements to generate combinations
    for (int i = start; i < mutatedAndOriginal.size(); ++i) {
        // Choose the current element
        current.push_back(mutatedAndOriginal[i]);

        // Recursively generate combinations with the remaining elements
        generateCombinations(mutatedAndOriginal, i + 1, current, result, numPep);

        // Backtrack - remove the current element to try other combinations
        current.pop_back();
    }
}

// Generate new combination of 3 sequences
vector<vector<vector<char>>> findCombinations(std::vector<vector<char>>& mutatedAndOriginal, const short& numPep) {
    vector<vector<vector<char>>> result;
    vector<vector<char>> current;
    
    // Start the process from the beginning of the array
    generateCombinations(mutatedAndOriginal, 0, current, result, numPep);

    return result;
}



GA_Parameters Mutation_withMotif(GA_Parameters parents) {
    int i;
    short numPep = parents.numPep;


    GA_Parameters mutatedOffsprings = parents;

    vector <double> MutationRate = parents.MutationRateXaa;
    vector<char> alphabet = parents.alphabet;


    int userSeqStart = parents.randomSeqLength;
    int userSeqEnd = parents.randomSeqLength + parents.motifLength;

    vector<vector<vector<char>>> newPopulation = parents.Helices; // vector to store the final results
    vector<int> newGlyPos;

    for (i=0; i< parents.populationSize; i++) {

        // Store mutated sequences in a vector
        vector<vector<char>> combined = parents.Helices[i];
        vector<vector<char>> mutatedSequences;
        mutatedSequences.clear();
        int GlyPos = parents.GlyPos;

        for (int j = 0; j < numPep; ++j) {
            vector<char> mutated;
            mutated = mutateSequence(parents.Helices[i][j], parents);
            if (parents.Helices[i][j] != mutated) mutatedSequences.push_back(mutated);
        }

        // push the mutated sequences to the combined vector that already includes the original sequences
        for (vector<char> seq: mutatedSequences) {
            combined.push_back(seq);

        }

        // Generate all possibles combination of mutated and original sequences
        vector<vector<vector<char>>> combinations = findCombinations(combined, numPep);
        
      
        for (int helixIndex = 0; helixIndex < combinations.size(); helixIndex++) {
            // Insert user motif to mutated sequences
            for (int pep = 0; pep < combinations[helixIndex].size(); pep++) {
                for (int aa = 0; aa < parents.motifLength; aa++) {
                    combinations[helixIndex][pep][aa+ parents.randomSeqLength] = parents.MotifSequences[pep][aa];
                }
            }
            if (isUnique(combinations[helixIndex], combinations, helixIndex)){
                newPopulation.push_back(combinations[helixIndex]);
                newGlyPos.push_back(GlyPos);
            }
        }


    }

    // update the mutatedOffspring struct
   
        mutatedOffsprings.populationSize = newPopulation.size();
        mutatedOffsprings.Helices = newPopulation;
        //mutatedOffsprings.GlyPosition = newGlyPos;



    return mutatedOffsprings;
}



GA_Parameters Mutation (GA_Parameters parents) {
    int i;
    short numPep = parents.numPep;
    
    GA_Parameters mutatedOffsprings = parents;

    vector<vector<vector<char>>> newPopulation; // vector to store the final results
    vector<int> newGlyPos;
    
    for (i=0; i< parents.populationSize; i++) {
       
        // Store mutated sequences in a vector
        vector<vector<char>> combined = parents.Helices[i];
        vector<vector<char>> mutatedSequences;
        mutatedSequences.clear();
        int GlyPos = parents.GlyPos;
        
        for (int j = 0; j < numPep; ++j) {
            vector<char> mutated;
            mutated = mutateSequence(parents.Helices[i][j], parents);
            if (parents.Helices[i][j] != mutated and Gly_repetition(mutated)) mutatedSequences.push_back(mutated);
        }
        
        // push the mutated sequences to the combined vector that already includes the original sequences
        for (vector<char> seq: mutatedSequences) {
            combined.push_back(seq);
            
        }
        
        // Generate all possibles combination of mutated and original sequences
    
        vector<vector<vector<char>>> combinations = findCombinations(combined, numPep);
        
        
        for (int helixIndex = 0; helixIndex < combinations.size(); helixIndex++) {
            if (isUnique(combinations[helixIndex], combinations, helixIndex)) {
                newPopulation.push_back(combinations[helixIndex]);
                newGlyPos.push_back(GlyPos);
            }
            
        }
        
    }
    
    // update the mutatedOffspring struct
        mutatedOffsprings.populationSize = newPopulation.size();
        mutatedOffsprings.Helices = newPopulation;

    
    
    return mutatedOffsprings;
}


// // // // // // // // // // // /// // // // // // // // FITNESS FUNCTION // // // // // // // // // // // /// // // // // // // //

void SCEPTTr12(TripleHelix * Library,
               parameterType helixParameters,
               GA_Parameters GAparameters,
               vector<double>& outTm,
               vector<double>& outSpecificity,
               vector<vector<int>>& outBestReg){
    short numPep = GAparameters.numPep;
   
    parameterType parameters;
    parameters = ReadParameters();
    
    // Read Parent Helices data into the TripleHelix struct
    short i,j;
    short x,y;
    
    // Initialize the Library
    for (int i=0; i<GAparameters.maxHelices; i++) {
        Library[i].initializeAll();
    }
    
    
    // Input values to the Library
    for (int i=0; i<GAparameters.populationSize; i++) {
        Library[i].initializeAll();
        Library[i].Nterm = "ac";
        Library[i].Cterm = "am";
        Library[i].numPep = numPep;
        Library[i].numAA = GAparameters.numAA;
        for (x=0; x< numPep; x++){
            for (y=0; y<GAparameters.numAA; y++){
                Library[i].sequences[x][y] = toupper(GAparameters.Helices[i][x][y]);
            }
        }
        Library[i].determine_reptition();
    }
    
    
    for (i=0;i<GAparameters.populationSize;i++)
        {
            Library[i] = ScoreHelix(parameters, Library[i]);
        }
    
    
    
    outTm.reserve(GAparameters.populationSize);
    outSpecificity.reserve(GAparameters.populationSize);
    vector<int> bestReg;
    bestReg.reserve(3);
    
    for (i=0; i< GAparameters.populationSize; i++) {
        
        bestReg.clear();
        for (j=0; j<3; j++) {
            bestReg.push_back(Library[i].bestRegister[j]);
        }
        outBestReg.push_back(bestReg);
        
        double CCTm = Library[i].CCTm;
        double spec = Library[i].specificity;
        
        if ((Library[i].bestRegister[0] == Library[i].bestRegister[1]) and (Library[i].bestRegister[0] == Library[i].bestRegister[2])) {
            spec = -10; // homotrimer
            
        }
        
        outSpecificity.push_back(spec);
        outTm.push_back(CCTm);
    }
    
}

vector<double> fitnessScore (GA_Parameters GAparameters,
                             TripleHelix * Library,
                            parameterType helixParameters,
                             vector<double>& MeltingTemp,
                             vector<double>& Spec) {
    vector<double> FitnessScore;
    FitnessScore.reserve(GAparameters.populationSize);
    FitnessScore.clear();
    
    
    parameterType parameters;
    parameters = ReadParameters();
    
    vector<double> Tm;
    Tm.reserve(GAparameters.populationSize);
    vector<double> specificity;
    specificity.reserve(GAparameters.populationSize);
    vector<vector<int>> bestRegister;
   
    
    // Score the library
    SCEPTTr12(Library, helixParameters, GAparameters, Tm, specificity, bestRegister);
    
 
   

    vector<int> CCRegisterScore;
   
    for (int i=0; i< bestRegister.size(); i++) {
        if ((bestRegister[i][0] == 0) and (bestRegister[i][1] == 1) and (bestRegister[i][2] == 2)) {
            CCRegisterScore.push_back(50);
            
        } else {
            CCRegisterScore.push_back(-50);
            
        }
        
    }
    
    
    FitnessScore.assign(GAparameters.populationSize, 0.0); // Fill with zeros
    
    double a = 0.5; // spec
    double b = 0.5; // melting
    double c = 0; // CCReg
    if (GAparameters.haveMotif) {
        a = 0.4;
        b = 0.4;
        c = 0.2;
    }
        

    for (size_t i = 0; i < GAparameters.populationSize; i++) {
        double fs = a * specificity[i] + b * Tm[i] + c * CCRegisterScore[i];
        FitnessScore[i] = fs;
       
    }
    
    MeltingTemp = Tm;
    Spec = specificity;
    
    
    return FitnessScore;
}

// // // // // // // // // // // /// // // // // // // // SELECTION // // // // // // // // // // // /// // // // // // // //

vector<int> findIndicesOfTwoHighest(GA_Parameters parameters, const vector<double>& fitnessScores) {
    vector <int> Indices(2);
    if (parameters.populationSize < 2) {
        // Handle the case where the vector has less than 2 elements
        // Return an invalid pair of indices
        Indices[0] = -1;
        Indices[1] = -1;
        return Indices;
    }

    int highestIndex = 0;
    int secondHighestIndex = 1;

    if (fitnessScores[1] > fitnessScores[0]) {
        highestIndex = 1;
        secondHighestIndex = 0;
    }

    for (int i = 2; i < parameters.populationSize; ++i) {
        if (fitnessScores[i] > fitnessScores[highestIndex]) {
            secondHighestIndex = highestIndex;
            highestIndex = i;
        } else if (fitnessScores[i] > fitnessScores[secondHighestIndex]) {
            secondHighestIndex = i;
        }
    }
    
    Indices[0] = highestIndex;
    Indices[1] = secondHighestIndex;

    return Indices;
}

GA_Parameters Selection (GA_Parameters Parents, const vector<double>& fitnessScores, vector<int>& bestHelixID) {
    
    vector<int> bestHelicesIndex = findIndicesOfTwoHighest(Parents,  fitnessScores);
    bestHelixID = bestHelicesIndex;
  
    GA_Parameters bestParents = Parents;
    // Redefine parameters of bestParents
    bestParents.populationSize = 2;

    int j,k;
    bestParents.Helices.clear();
    
    vector<vector<char>> BestHelix;
    vector<vector<char>> secBestHelix;
    vector<char> seq1(Parents.numAA);
    vector<char> seq2(Parents.numAA);
    for (j = 0; j < Parents.numPep; j++) {
        seq1.clear();
        seq2.clear();
        for (k = 0; k < Parents.numAA; k++) {
            seq1.push_back(Parents.Helices[bestHelicesIndex[0]][j][k]);
            seq2.push_back(Parents.Helices[bestHelicesIndex[1]][j][k]);
        }
        BestHelix.push_back(seq1);
        secBestHelix.push_back(seq2);
    }
    
    bestParents.Helices.push_back(BestHelix);
    bestParents.Helices.push_back(secBestHelix);
    
    return bestParents;
}





