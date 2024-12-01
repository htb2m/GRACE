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

using namespace std;


int main(int argc, const char * argv[]) {
   
 
    short i;
    
    srand(static_cast<unsigned>(time(0))); // Initialize random seed

    
    
// // // // // // // // // // // // // // // // USER INPUT //  // // // // // // // // // // // // // //
    
    // Read the AminoAcid.csv file
    ifstream file("AminoAcids.csv");
    vector<int> Xaa;
    vector<int> Yaa;
    vector<char> validAA;

    string line;
    getline(file, line); //skip the table header
    
    // Check if file is opened successfully
    if (!file.is_open()) {
        cerr << "Unable to open file" << std::endl;
        return 1; // Return with error code
    }

    while (getline(file, line)) {
        stringstream ss(line);
        string item;
        int value;
        char AA;

        getline(ss, item, ',');
        if (stringstream(item) >> AA) validAA.push_back(tolower(AA));
        
        // Get Xaa value
        getline(ss, item, ',');
        if (stringstream(item) >> value) Xaa.push_back(value);
        
        // Get Yaa value
        getline(ss, item, ',');
        if (stringstream(item) >> value) Yaa.push_back(value);
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

    cout << "Do you want to generate a novel heterotrimer or a heterotrimer including recognition epitope? Please enter '0' for novel heterotrimer or '1' heterotrimer including recognition epitope. " << endl; // Sounds bad, will change this prompt later
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
        cout << "Please specify the number of amino acids present in your sequences (more than 3 and less than 15) " << endl;
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
            cout << "Please enter your flank sequence 1: " << endl;
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
            cout << "Please enter your flank sequence 2: " << endl;
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
            cout << "Please enter your flank sequence 3: " << endl;
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
        
        cout << "Do you want to generate AAB or ABC heterotrimer? Please enter '0' for AAB or '1' for ABC." << endl;
        while (!(cin >> tempInput) || (tempInput != 0 and tempInput != 1)) {
            cout << "Please enter either '0' for AAB or '1' for ABC: " << endl;
            cin.clear();  // Clear error flag
            cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
        }
        
        if (tempInput == 0) userAAB = true;
      
        
        // Input sequence length
        cout << "Please enter a number between 21 and 40 as peptide Length: " << endl;
        while (!(cin >> userLength) || userLength < 21 || userLength > 40) {
            cout << "Please enter a number between 21 and 40: " << endl;
            cin.clear();  // Clear error flag
            cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
        }
    }
    
    
    cout << "Please enter a number between 30 to 50 as target Tm: " << endl;
    while (!(cin >> userTm) || userTm < 30 || userTm > 70) {
        cout << "Please enter a number between 30 to 50: " << endl;
        cin.clear();  // Clear error flag
        cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
    }


    cout << "Please enter a number between 10 to 20 as target Specificity: " << endl;
    while (!(cin >> userSpec) || userSpec < 10 || userSpec > 35) {
        cout << "Please enter a number between 10 to 20: " << endl;
        cin.clear();  // Clear error flag
        cin.ignore(numeric_limits<streamsize>::max(), '\n');  // Ignore invalid input
    }
    
   
   
    
    
// // // // // // // // // // // // // // ///   END OF UI  // // // // // // // // // // // // // // // //
    

        cout << endl << endl;
        //cout << "Sequence Length: " << userLength << endl;
        cout << "Target Tm = " << userTm << "; Target Specificity = " << userSpec << endl;
        cout << endl;
        cout << "GENERATING A HETEROTRIMER" << endl << endl;
        
        GA_Parameters Parents;
        
        
        TripleHelix Library[Parents.maxHelices];
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
        }
        
        
        Parents.numAA = userLength;  // Initialize sequence length (default = 30)
        if (Parents.haveMotif) Parents.numAA = motifLength + Parents.randomSeqLength*2;


        // // // // // // // // // // //   GENERATE POPULATION  // // // // // // // // // // // // //
        

        // Generate ParentHelices
        Parents = initialPopulationGenerator(Parents); // random sequence
     
        
        // // // // // // // // // // //   SCORE POPULATION // // // // // // // // // // // // //
        
        
        vector<double> FNScore = fitnessScore(Parents, Library, parameters);

        
        // // // // // // // // // // // // //  SELECTION   // // // // // // // // // // // // //
       
        
        //GA_Parameters bestParents = Selection(Parents, FNScore);
        GA_Parameters bestParents = Selection(Parents, FNScore);
        
       
        Parents = bestParents; // 2 helices with best fitness score are now Parents for the next generation
       

        // // // // // // // // // // // //  MUTATION  // // // // // // // // // // // // //
     
       
     
       GA_Parameters MutatedOS;
       if (!haveMotif) {
           MutatedOS = Mutation(Parents);
       } else {
           MutatedOS = Mutation_withMotif(Parents);
       }
    
        Parents = MutatedOS;
       
        
        if (Parents.haveMotif) {
                // // // // // // // // // // //   SCORE POPULATION // // // // // // // // // // // // //

                FNScore = fitnessScore(Parents, Library, parameters);


                // // // // // // // // // // // // //  SELECTION   // // // // // // // // // // // // //

                bestParents = Selection(Parents, FNScore);

                Parents = bestParents; // 2 helices with best fitness score are now Parents for the next generation

        }


        
        // // // // // // // // // // // //  CROSSOVER  // // // // // // // // // // // // //
        
        GA_Parameters CrossoverOS;
        if (Parents.haveMotif) {
            CrossoverOS = CrossOver_withMotif(Parents);
        } else {
            CrossoverOS = CrossOver(Parents);
        }
        
        
         
        // // // // // // // // // // //   SCORE POPULATION // // // // // // // // // // // // //
        
        
        FNScore = fitnessScore(Parents, Library, parameters);
       
        
        // // // // // // // // // // // // //  SELECTION   // // // // // // // // // // // // //
       
        
        bestParents = Selection(Parents, FNScore);
        
        Parents = bestParents; // 2 helices with best fitness score are now Parents for the next generation
       
       
       
        // // // // // // // // // // // // //  EVOLUTION LOOP  // // // // // // // // // // // // //
        

        int rounds = 0;
        bool done = false;
        
        while (not done) {
            
            
            // // // // // // // // // // // DISPLAY SEARCHING PROGRESS  // // // // // // // // // // //
            if ((rounds % 50 == 0) and ((rounds +1) > 1)) {
                cout << endl << endl;
                cout << "--------------------------------------------------------" << endl;
                cout << "GENERATION NUMBER " << rounds + 1 << endl;
                cout << "High Tm = Propensity + PairWise = " << Library[0].HighTm << " = " << Library[0].bestPropensity << " + " << Library[0].bestPairwise << endl;
                cout << "BestRegister " << Library[0].bestRegister[0] << Library[0].bestRegister[1] << Library[0].bestRegister[2] << ". Specificity = " << Library[0].specificity << endl << endl;
                
                
                // Display the ParentHelices
                cout << "parents" << endl;
                for (int x = 0; x < Parents.populationSize; x++) {
                    cout << "Helix : " << x +1 << endl;
                    for (int y = 0; y < Parents.numPep; y++) {
                        for (int z = 0; z < Parents.numAA; z++) {
                            cout << Parents.Helices[x][y][z];
                        }
                        cout << endl;
                    }
                }
               
                cout << endl;
                cout <<" Fitness score" << endl;
                for (int x = 0; x < Parents.populationSize; x++) {
                    cout << FNScore[x] << " ";
                }
                cout << endl;
                
            }
           
            cout << "Generation number " << rounds + 1 << " - ";
            
            // // // // // // // // // // // // //   MUTATION  // // // // // // // // // // // // //
            
            if (Parents.haveMotif) {
                MutatedOS = Mutation_withMotif(Parents);
            } else {
                MutatedOS = Mutation(Parents);
            }
            Parents = MutatedOS;
            cout << "Mutation..";
            
            
            if (Parents.haveMotif) {
                // // // // // // // // // // // // FITNESS COMPUTING  // // // // // // // // // // // // //
                
                FNScore = fitnessScore(Parents, Library, parameters);
                
                cout << "Fitness computing..." ;
               
                // // // // // // // // // // // BEST PARENTS SELECTION  // // // // // // // // // // // //
                bestParents = Selection(Parents, FNScore);
                
                Parents = bestParents; // bestParents are now Parent for the next generation
                
                cout << "Selecting parent helices...";
                
            }

            
            // // // // // // // // // // // // //  CROSSOVER  // // // // // // // // // // // // //
            
            GA_Parameters CrossoverOS;
            if (Parents.haveMotif) {
                CrossoverOS = CrossOver_withMotif(Parents);
            } else {
                CrossoverOS = CrossOver(Parents);
            }
            
            Parents = CrossoverOS;
            cout << "Crossover..";
            

            
           
            
            
            
            
            
            // // // // // // // // // // // // FITNESS COMPUTING  // // // // // // // // // // // // //
            
            FNScore = fitnessScore(Parents, Library, parameters);
            
            cout << "Fitness computing..." ;
           
            // // // // // // // // // // // BEST PARENTS SELECTION  // // // // // // // // // // // //
            bestParents = Selection(Parents, FNScore);
            
            Parents = bestParents; // bestParents are now Parent for the next generation
            
            cout << "Selecting parent helices...";
            
            cout << endl;
           
            
            // // // // // // // // // // // LOOP-BREAKING CONDITIONS // // // // // // // // // // // //
            
            if (((Library[0].HighTm >= userTm ) and (Library[0].specificity >= userSpec)  and  ((Library[0].bestRegister[0] != Library[0].bestRegister[1]) or (Library[0].bestRegister[0] != Library[0].bestRegister[2]) or (Library[0].bestRegister[1] != Library[0].bestRegister[2])))
                or
                (rounds == 500000)) // return the closest result
            {
                done = true;
            }
           
            rounds += 1;
           
        }
        
        // // // // // // // // // // // DISPLAY FINAL RESULTS // // // // // // // // // // // //
        cout << endl << endl << endl;
        cout << "-------------------------------------------------------------" << endl;
        if (rounds == 500000) {
            cout << "Cannot come up with sequences that satisfy all the given conditions. Displaying the best result found: " << endl << endl;
        }
        
        cout << "STOP AT GENERATION NUMBER " << rounds << endl << endl;
        
        cout << "FITNESS SCORE: " << endl;
        for (i=0; i<Parents.populationSize; i++){
                cout << FNScore[i] << " ";
        }
        
        cout << endl;
        cout << "------------------------" << endl << endl;
        
        for (i=0; i < 1; i++) {
            cout << "High Tm = Propensity + PairWise = " << Library[i].HighTm << " = " << Library[i].bestPropensity << " + " << Library[i].bestPairwise << endl;
            Library[i].dissect();
            Library[i].userOutput();
        
        
        }
        
        cout << endl;
        
     
        
    
    
    
   
    return 0;
}
