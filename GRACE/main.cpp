//
//  main.cpp
//  GRACE
//
//  Created by HA BUI from Hartgerink Lab on 10/18/23.
//


#include "Functions.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <string>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <unordered_set>
#include <set>
#include <chrono>
#include <algorithm>

using namespace std;

int main(int argc, const char * argv[]) {
   
 
    short i;
    
    srand(static_cast<unsigned>(time(0))); // Initialize random seed

    
    
// // // // // // // // // // // // // // // // USER INTERFACE //  // // // // // // // // // // // // // //
    cout << "---------------------------------------------------------" << endl;
    cout << "Genetically Refining Algorithm for Collagen Engineering" << endl;
    cout << "                      GRACE v1.0                       " << endl;
    cout << "    by T.HA BUI from Hartgerink Lab @ Rice University  " << endl << endl;
    cout << "---------------------------------------------------------" << endl;
    
    cout << endl;
    
    // Read the AminoAcid.csv file
    ifstream file("AminoAcids.csv");
    vector<int> Xaa;
    vector<int> Yaa;
    vector<char> validAA;

    string line;
    getline(file, line); //skip the table header
    
    // Check if file is opened successfully
    if (!file.is_open()) {
        cerr << "Unable to open file" << endl;
        return 1; // Return with error code
    }
    


    while (getline(file, line)) {
        stringstream ss(line);
        string item;
        int value;
        char AA;

        // Read the amino acid
        getline(ss, item, ',');
        if (stringstream(item) >> AA) {
            validAA.push_back(tolower(AA)); // Store amino acid as lowercase
        }

        // Get Xaa value
        getline(ss, item, ',');
        if (stringstream(item) >> value) {
            if (value != 0 && value != 1) {
                cerr << "Error: Invalid Xaa value '" << value << "' found in AminoAcids.csv file AminoAcids.csv file. Only '0' and '1' are allowed." << endl;
                return 1; // Exit with error code
            }
            Xaa.push_back(value);
        } else {
            cerr << "Error: Invalid Xaa value found. Must be a number (0 or 1)." << endl;
            return 1;
        }

        // Get Yaa value
        getline(ss, item, ',');
        if (stringstream(item) >> value) {
            if (value != 0 && value != 1) {
                cerr << "Error: Invalid Yaa value '" << value << "' found in AminoAcids.csv file. Only '0' and '1 are allowed." << endl;
                return 1; // Exit with error code
            }
            Yaa.push_back(value);
        }
        else {
            cerr << "Error: Invalid Yaa value found AminoAcids.csv file. Must be a number (0 or 1)." << endl;
            return 1;
        }
    }

    file.close();
    
    bool excludeXaa = false;
    int numXaaExcluded = 0;
    vector<char> excludedXaaList;
    excludedXaaList.clear();

    bool excludeYaa = false;
    int numYaaExcluded = 0;
    vector<char> excludedYaaList;
    excludedYaaList.clear();
  
    for (i = 0; i < validAA.size(); i++) {
        if (Xaa[i] == 0) {
            numXaaExcluded += 1;
            excludedXaaList.push_back(validAA[i]);
        }
        if (Yaa[i] == 0) {
            numYaaExcluded += 1;
            excludedYaaList.push_back(validAA[i]);
        }
    }
    
    if (numXaaExcluded != 0) excludeXaa = true;
    if (numYaaExcluded != 0) excludeYaa = true;
  
    
    // ABC or AAB heterotrimer
    
    int tempInput = -1;
    bool userAAB = false;
    int userLength = -1;
    int userTm = -1;
    int userSpec = -1;
    
   
    bool haveMotif = false;

    cout << "Do you want to generate a heterotrimer or a heterotrimer including recognition epitope?" << endl;
    cout << "Please enter '0' for novel heterotrimer or '1' heterotrimer including recognition epitope. " << endl;
        while (!(cin >> tempInput) || (tempInput != 1 and tempInput != 0)) {
           cout << "Please enter either '0' or '1': " << endl;
           cin.clear();  // Clear error flag
           cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
       }

    if (tempInput == 1) haveMotif = true;

 
    vector<char> Seq1;
    vector<char> Seq2;
    vector<char> Seq3;
    
    int firstGPosition = -1;
    int motifLength = -1;
    

    
    if (haveMotif) {
        cout << "Please specify the number of amino acids present in your sequences" << endl;
        cout << " (more than 3 and less than 15) " << endl;
        motifLength = -1;
        while ( !(cin >> motifLength) || (motifLength < 3) || motifLength > 15) {
            cout << "Please enter a valid number" << endl;
            cin.clear();  // Clear error flag
            cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
        }
        
        
        Seq1.clear();
        int firstG1 = -1;
        while (true) {
            Seq1.clear();
            string input1;
            cout << "Please enter your sequence 1: " << endl;
            cin >> input1;

            // Check if the input matches the pre-determined length
            if (input1.length() != motifLength) {
                cout << "Invalid input length. Please try again." << endl;
                continue; // Skip the remaining code and start the next loop iteration
            }
            // Check if there are only cannonical amino acids + Hyp
           
            bool allCharsValid = true;
            for (char c: input1) {
                if (std::find(validAA.begin(), validAA.end(), tolower(c)) == validAA.end()) {
                   allCharsValid = false;
                   break; // No need to check the rest
                }
            }
            
            if (allCharsValid) {
                for (char c: input1) {
                    Seq1.push_back(toupper(c));
                }
                //break; // Break out of the loop if all characters are valid
            } else {
                cout << "Invalid amino acid found. Please make sure the sequence only include canonical amino acids." << std::endl;
                continue;
            }

            // Check Gly at every third position
            bool GlyAtEveryThird = false;
            firstG1 = -1;

            GlyAtEveryThird = findGlyAtEveryThird(motifLength, Seq1, firstG1);
            
            string checkGly = " ";
            if (GlyAtEveryThird) {
                checkGly = "Yes";
            } else {checkGly = "No";}
            

            cout << "User input is a collagen-like peptide ? " << checkGly << endl;
            if (GlyAtEveryThird) {
                break;
            } else {
                cout << "Please make sure the input sequence is a collagen-like sequence. " << endl;
            }
        }
        

        while (true) {
            Seq2.clear();

            string input2;
            cout << "Please enter your sequence 2: " << endl;
            cin >> input2;

            // Check if the input length
            if (input2.length() != motifLength) {
                cout << "Invalid input length. Please try again." << endl;
                continue; // Skip the remaining code and start the next loop iteration
            }
            bool allCharsValid = true;
            for (char c: input2) {
                if (find(validAA.begin(), validAA.end(), tolower(c)) == validAA.end()) {
                   allCharsValid = false;
                   break; // No need to check the rest
                }
            }
            if (allCharsValid) {
                for (char c: input2) {
                    Seq2.push_back(toupper(c));
                }
                //continue; // Break out of the loop if all characters are valid
            } else {
                cout << "Invalid amino acid found. Please make sure the sequence only include canonical amino acids." << std::endl;
                continue; // Skip the remaining code and start the next loop iteration
            }

            // Check Gly at every third position
            bool GlyAtEveryThird = false;
            int firstG2 = -1;

            GlyAtEveryThird = findGlyAtEveryThird(motifLength, Seq2, firstG2);

            //cout << "check Gly " << GlyAtEveryThird << " First G2 " << firstG2 << endl;

            if (GlyAtEveryThird == false) {
                cout << "Please make sure the input sequence is a collagen-like sequence. " << endl;
                continue;
            }


            if (firstG2 == firstG1) {
                break;
            } else {
                cout << "Glycine positions do not matched. Please try again." << endl;
                continue;
            }

        }
        while (true) {
            Seq3.clear();

            string input3;
            cout << "Please enter your sequence 3: " << endl;
            cin >> input3;

            // Check if the input length
            if (input3.length() != motifLength) {
                cout << "Invalid input length. Please try again." << endl;
                continue; // Skip the remaining code and start the next loop iteration
            }
            bool allCharsValid = true;
            for (char c: input3) {
                if (find(validAA.begin(), validAA.end(), tolower(c)) == validAA.end()) {
                   allCharsValid = false;
                   break; // No need to check the rest
                }
            }
            if (allCharsValid) {
                for (char c: input3) {
                    Seq3.push_back(toupper(c));
                }
                //continue; // Break out of the loop if all characters are valid
            } else {
                cout << "Invalid amino acid found. Please make sure the sequence only include canonical amino acids." << std::endl;
                continue; // Skip the remaining code and start the next loop iteration
            }

            // Check Gly at every third position
            bool GlyAtEveryThird = false;
            int firstG3 = -1;

            GlyAtEveryThird = findGlyAtEveryThird(motifLength, Seq3, firstG3);

            //cout << "check Gly " << GlyAtEveryThird << " First G3 " << firstG3 << endl;

            if (GlyAtEveryThird == false) {
                cout << "Please make sure the input sequence is a collagen-like sequence. " << endl;
                continue;
            }


            if (firstG3 == firstG1) {
                break;
            } else {
                cout << "Glycine positions do not matched. Please try again." << endl;
                continue;
            }

        }
            
        firstGPosition = firstG1;
            
    
    } else {
        
        cout << "Do you want to generate AAB or ABC heterotrimer?" << endl;
        cout << "Please enter '0' for AAB or '1' for ABC." << endl;
        while (!(cin >> tempInput) || (tempInput != 0 and tempInput != 1)) {
            cout << "Please enter either '0' for AAB or '1' for ABC: " << endl;
            cin.clear();  // Clear error flag
            cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
        }
        
        if (tempInput == 0) userAAB = true;
      
        
        // Input sequence length
        cout << "Please enter a number between 24 and 40 as peptide Length: " << endl;
        while (!(cin >> userLength) || userLength < 21 || userLength > 40) {
            cout << "Please enter a number between 21 and 40: " << endl;
            cin.clear();  // Clear error flag
            cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
        }
    }
    
    
    cout << "Please enter a number between 40 to 55 as target Tm: " << endl;
    while (!(cin >> userTm) || userTm < 30 || userTm > 70) {
        cout << "Please enter a number between 40 to 55: " << endl;
        cin.clear();  // Clear error flag
        cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
    }


    cout << "Please enter a number between 10 to 25 as target Specificity: " << endl;
    while (!(cin >> userSpec) || userSpec < 10 || userSpec > 35) {
        cout << "Please enter a number between 10 to 25: " << endl;
        cin.clear();  // Clear error flag
        cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
    }
    
   
   
    
   
// // // // // // // // // // // // // // ///   END OF UI  // // // // // // // // // // // // // // // //
    

    cout << endl << endl;
    //cout << "Sequence Length: " << userLength << endl;
    cout << "Target Tm = " << userTm << "; Target Specificity = " << userSpec << endl;
    
    
    GA_Parameters Parents;
    
    
    TripleHelix*Library = new TripleHelix[Parents.maxHelices];
    parameterType parameters;
    parameters = ReadParameters();
   


   
    // // // // // // // // // // //  INITIALIZE Parents STRUCT WITH UI // // // // // // // // // // //
    
    Parents.AAB = userAAB; // ABC or AAB; default is ABC
    if (userAAB) Parents.numPep = 2;
   
 
    // Excluded amino acid
    Parents.excludeXaa = excludeXaa;
    Parents.excludedXaaList = excludedXaaList;
    Parents.numExcludedXaa = numXaaExcluded;
    Parents.excludedXaaList.clear();
    for (char aminoAcid : excludedXaaList) {
        Parents.setXaaMutationRateToZero(aminoAcid);
        Parents.excludedXaaList.push_back(aminoAcid);
    }
    
    Parents.excludeYaa = excludeYaa;
    Parents.excludedYaaList = excludedYaaList;
    Parents.numExcludedYaa = numYaaExcluded;
    Parents.excludedYaaList.clear();
    for (char aminoAcid : excludedYaaList) {
        Parents.setYaaMutationRateToZero(aminoAcid);
        Parents.excludedYaaList.push_back(aminoAcid);
    }
    
    Parents.numAA = userLength;  // Initialize sequence length (default = 30)
    
    Parents.haveMotif = haveMotif;
    
    if (Parents.haveMotif) {
        Parents.haveMotif = haveMotif;
        Parents.motifLength = motifLength;
        Parents.MotifSequences.resize(Parents.numPep);
        Parents.MotifSequences[0] = Seq1;
        Parents.MotifSequences[1] = Seq2;
        if (not Parents.AAB) {
            Parents.MotifSequences[2] = Seq3;
        }
        Parents.GlyPos = firstGPosition;
        Parents.numAA = motifLength + Parents.randomSeqLength*2;
    }
    

        
        
    cout << endl;
    cout << "Generating random population..." << endl << endl;
    
    
    
 

    // Holder for fitness lanscape tracking
    vector <double> FSHolder;
    vector <int> GenerationHolder;
    vector <double> TmHolder;
    vector<double> SpecHolder;
    vector<double> TimeElapsed;
    vector<double> MeltingTemp;
    vector<double> Spec;
    auto start_time = std::chrono::high_resolution_clock::now();

    // // // // // // // // // // //   GENERATE POPULATION  // // // // // // // // // // // // //
    

    // Generate ParentHelices
    Parents = initialPopulationGenerator(Parents); // random sequence
    

   
    // // // // // // // // // // //   SCORE POPULATION // // // // // // // // // // // // //
    
    MeltingTemp.clear();
    Spec.clear();
    vector<double> FNScore = fitnessScore(Parents, Library, parameters, MeltingTemp, Spec);

    
    // // // // // // // // // // // // //  SELECTION   // // // // // // // // // // // // //
   
    vector<int> bestHelixID;
    //GA_Parameters bestParents = Selection(Parents, FNScore);
    GA_Parameters bestParents = Selection(Parents, FNScore, bestHelixID);
    
    int bestID = bestHelixID[0];
    FSHolder.push_back(FNScore[bestID]);
    TmHolder.push_back(MeltingTemp[bestID]);
    SpecHolder.push_back(Spec[bestID]);
    GenerationHolder.push_back(0);
    bestHelixID.clear();
    // set the clock to track run time
    //auto start_time = std::chrono::high_resolution_clock::now();
    
    auto now = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = now - start_time;
    TimeElapsed.push_back(elapsed.count());
   
    Parents = bestParents; // 2 helices with best fitness score are now Parents for the next generation
   


    
    // // // // // // // // // // // //  CROSSOVER  // // // // // // // // // // // // //
    
    GA_Parameters CrossoverOS;
    if (Parents.haveMotif) {
        CrossoverOS = CrossOver_withMotif(Parents);
        Parents = CrossoverOS;
        MeltingTemp.clear();
        Spec.clear();
        FNScore = fitnessScore(Parents, Library, parameters, MeltingTemp, Spec);
          // // // // // // // // // // // // //  SELECTION   // // // // // // // // // // // // //
          bestParents = Selection(Parents, FNScore, bestHelixID);

          Parents = bestParents; // 2 helices with best fitness score are now Parents for the next generation
    } else {
        CrossoverOS = CrossOver(Parents);
        Parents = CrossoverOS;
    }


     // // // // // // // // // // //  MUTATION  // // // // // // // // // // // // //
   
   
   
      GA_Parameters MutatedOS;
      if (!haveMotif) {
          MutatedOS = Mutation(Parents);
      } else {
          MutatedOS = Mutation_withMotif(Parents);
      }
   
       Parents = MutatedOS;
   
     
    // // // // // // // // // // //   SCORE POPULATION // // // // // // // // // // // // //
    
    MeltingTemp.clear();
    Spec.clear();
    FNScore = fitnessScore(Parents, Library, parameters, MeltingTemp, Spec);
    
   

   
    // // // // // // // // // // // // //  SELECTION   // // // // // // // // // // // // //
   
    
    bestParents = Selection(Parents, FNScore, bestHelixID);
    
    Parents = bestParents; // 2 helices with best fitness score are now Parents for the next generation
   
    bestID = bestHelixID[0];
    FSHolder.push_back(FNScore[bestID]);
    GenerationHolder.push_back(1);
    TmHolder.push_back(MeltingTemp[bestID]);
    SpecHolder.push_back(Spec[bestID]);
    bestHelixID.clear();
    
    now = chrono::high_resolution_clock::now();
    elapsed = now - start_time;
    TimeElapsed.push_back(elapsed.count());
   
    // // // // // // // // // // // // //  EVOLUTION LOOP  // // // // // // // // // // // // //
    

    int rounds = 0;
    bool done = false;
    
    while (not done) {
        
        
        // // // // // // // // // // // DISPLAY SEARCHING PROGRESS  // // // // // // // // // // //
        if ((rounds % 100 == 0) and ((rounds +1) > 1)) {
            cout << endl ;
            cout << "-----------------------------------------------------------------------------------------------------" << endl;
            cout << "GENERATION NUMBER " << rounds << endl;
            cout << "High Tm = Propensity + PairWise = " << Library[0].HighTm << " = " << Library[0].bestPropensity << " + " << Library[0].bestPairwise << endl;
            cout << "BestRegister " << Library[0].bestRegister[0] << Library[0].bestRegister[1] << Library[0].bestRegister[2] << ". Specificity = " << Library[0].specificity << endl << endl;
            
            for (int a = 0; a < Parents.numPep; a++) {
                for (int b = 0; b < Parents.numAA; b++) {
                    cout << Parents.Helices[0][a][b];
                }
                cout << endl;
            }
            
           
            cout << "-----------------------------------------------------------------------------------------------------" << endl;
            cout << endl << endl;
            
        }
       
        cout << "Generation number " << rounds + 1 << " - ";
        
        // // // // // // // // // // // // //  CROSSOVER  // // // // // // // // // // // // //
        //CrossoverOS;
        if (Parents.haveMotif) {
           
            CrossoverOS = CrossOver_withMotif(Parents);
            Parents = CrossoverOS;
            MeltingTemp.clear();
            Spec.clear();
            FNScore = fitnessScore(Parents, Library, parameters, MeltingTemp, Spec);
            // // // // // // // // // // // // //  SELECTION   // // // // // // // // // // // // //
            bestParents = Selection(Parents, FNScore, bestHelixID);
            Parents = bestParents;
        } else {
            CrossoverOS = CrossOver(Parents);
            Parents = CrossoverOS;
        }
        
       
        cout << "Crossover..";
        
           
        // // // // // // // // // // // // //   MUTATION  // // // // // // // // // // // // //
        
        if (Parents.haveMotif) {
            MutatedOS = Mutation_withMotif(Parents);
        } else {
            MutatedOS = Mutation(Parents);
        }
        Parents = MutatedOS;
        cout << "Mutation...";
        
    
        
        // // // // // // // // // // // // FITNESS COMPUTING  // // // // // // // // // // // // //
        MeltingTemp.clear();
        Spec.clear();
        FNScore = fitnessScore(Parents, Library, parameters, MeltingTemp, Spec);
        
        cout << "Fitness computing..." ;
       
        // // // // // // // // // // // BEST PARENTS SELECTION  // // // // // // // // // // // //
        bestParents = Selection(Parents, FNScore, bestHelixID);
        
        Parents = bestParents; // bestParents are now Parent for the next generation
        
        cout << "Selecting parent helices...";
        
        cout << endl;
       
        int bestID = bestHelixID[0];
        FSHolder.push_back(FNScore[bestID]);
        TmHolder.push_back(MeltingTemp[bestID]);
        SpecHolder.push_back(Spec[bestID]);
        GenerationHolder.push_back(rounds + 2);
        bestHelixID.clear();
        
        now = chrono::high_resolution_clock::now();
        elapsed = now - start_time;
        TimeElapsed.push_back(elapsed.count());
        
        // // // // // // // // // // // LOOP-BREAKING CONDITIONS // // // // // // // // // // // //
        
        if (((Library[0].HighTm >= userTm ) and (Library[0].specificity >= userSpec)  and  ((Library[0].bestRegister[0] != Library[0].bestRegister[1]) or (Library[0].bestRegister[0] != Library[0].bestRegister[2]) or (Library[0].bestRegister[1] != Library[0].bestRegister[2])))
            or
            (rounds == 500000)) // return the closest result if no satisfied sequences were found after 500000 generations
        {
            done = true;
        }
       
        rounds += 1;
       
    }
    
    // // // // // // // // // // // DISPLAY FINAL RESULTS // // // // // // // // // // // //
    cout << endl << endl << endl;
    cout << "-------------------------------------------------------------" << endl;
    if (rounds == 500001) {
        cout << "Cannot come up with sequences that satisfy all the given conditions. Displaying the best result found: " << endl << endl;
    }
    
    cout << "STOP AT GENERATION NUMBER " << rounds << endl << endl;
    
    cout << "FITNESS SCORE: ";
    
    cout << FNScore[0] << endl;
    
    if (Parents.haveMotif) {
        cout << "User inputs: " << endl;
        for (int i = 0; i < Parents.MotifSequences.size(); i++) {
            cout << "Sequence " << i << ": ";
            for (int j = 0; j < Parents.MotifSequences[i].size(); j++) {
                cout << Parents.MotifSequences[i][j];
            }
            cout << endl;
        }
    }

    
    cout << endl;
    cout << "-----------------------------------" << endl << endl;
    cout << "Sequences were generated by pararameter sets curated on " << parameters.Date << endl;
    cout << endl;
    
    for (i=0; i < 1; i++) {
        cout << "High Tm = Propensity + PairWise = " << Library[i].HighTm << " = " << Library[i].bestPropensity << " + " << Library[i].bestPairwise << endl;
        Library[i].dissect();
        Library[i].userOutput();
    
    
    }
    

    
 
    // Check size consistency
    if (FSHolder.size() != GenerationHolder.size()) {
        cerr << "Error: Mismatched vector sizes.\n";
        return 1;
    }

    if (FSHolder.size() != TmHolder.size()) {
        cerr << "Error: Mismatched vector sizes Tm.\n";
        return 1;
    }
    
    if (FSHolder.size() != SpecHolder.size()) {
        cerr << "Error: Mismatched vector sizes Spec.\n";
        return 1;
    }
    

    
    // Open CSV file
    std::ofstream outFile("FitnessLandscape.csv");
    if (!outFile.is_open()) {
        cerr << "Failed to open file for writing.\n";
        return 1;
    }

    // Write header
    outFile << "Generation,TimeElapsed,FitnessScore,Tm,Spec\n";
   
    // Write data rows
    for (size_t i = 0; i < FSHolder.size(); ++i) {
        outFile << GenerationHolder[i] << "," << TimeElapsed[i] << "," << FSHolder[i] << "," << TmHolder[i] << "," << SpecHolder[i] << "\n";
    }

    outFile.close();
    cout << "Fitness score to FitnessLandscape.csv successfully.\n";
    
    
    
    delete [] Library;
    return 0;
}
