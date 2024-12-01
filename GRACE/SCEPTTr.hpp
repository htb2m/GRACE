//
//  SCEPTTr.hpp
//  GRACE
//
//  Created by Hartgerink Lab on 11/30/24.
//

#ifndef SCEPTTr_hpp
#define SCEPTTr_hpp

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <ctime>

using namespace std;

struct parameterType
{
    // These are all the parameters that can be read in from file or optimized
    double axial[27][27];
    double lateral[27][27];
    double propensityX[27];
    double propensityY[27];
    // Added JuL 01 // //
    double FrameShift[6][6];
    // // // // // // //
    double  A, B, C; // A + Bx + Cx^2 for the length parameter. Starting values of -70, 6.55, -0.065. We are not currently optimizing these.
    double charge;
    double Nterm;
    double Cterm;
    
    // these are all booleans that set flags to determin if that parameter will be optimized or not
    bool    optLength;
    bool    optPropX[27];
    bool    optPropY[27];
    bool    optAxial[27][27];
    bool    optLat[27][27];
    bool    optFrameShift[6][6];
    
    // These are the "initial" values set from experimentation (usually). We use these values during optimization to make sure we do not deviate from them by "too much". Too much is usually set to "by more than 2".
    double exAxial[27][27];
    double exLateral[27][27];
    double exPropensityX[27];
    double exPropensityY[27];
    // Added Jul 1 // //
    double exFrameShift[6][6];
    // // // // // // //
    double exA, exB, exC; // A + Bx + Cx^2 for the length parameter. Starting values of -70, 6.55, -0.065.
    double exCharge;
    double exNterm;
    double exCterm;
    
    string Nterm_type[6];
    string Cterm_type[6];
};

struct TripleHelix
{
    short   numPep;
    short   numAA;
    char    sequences[3][100];
    string  Nterm, Cterm;
    //short   n_term, c_term;
    
    double  expTm;
    double  CCTm; // Tm of the "correct" composition. For A2B systems this requires at least one of each peptide (no homotrimers). For ABC systems this requires one of each peptide (no homotrimers and no A2B systems).
    
    double  HighTm;
    double  deviation;
    double  secTm; // melting temperature of the second best register/composition
    double  specificity; // difference in melting temperature between the best and second best composition / register
    
    short   bestRegister[4];
    short   secRegister[4];
    short   CCRegister[4];
    
    
    string year;
    string ref;
    double SCEPTTr_10;
    double SCEPTTr_11;
    int index;
    
    
    // These arrays contain the following:
    // [0] which peptide is in the leading position (peptide numbers are 0-2)
    // [1] which peptide is in the middle position
    // [2] which peptide is in the trailing position
    // [3] This is the value of the Offset with 9 possible values as follows:
    // Values range from 0-8. Nine Offsets or Staggers can be considered.
    // 0) {012} This is the canonical offset ** mostly we consider only this **
    // --------------
    // 1) {015} -1 triplet.     trailing strand offset by an additional 3 amino acids
    // 2) {042} -1 triplet.     middle strand offset by an additional 3 amino acids
    // 3) {045} -1 triplet.     middle and trailing strand offset by an additional 3 amino acids
    // --------------
    // 4) {018} -2 triplets.    trailing strand offset by an additional 6 amino acids
    // 5) {048} -2 triplets.    middle stand offset by 3, trailing by 6
    // 6) {072} -2 triplets.    middle strand offset by 6
    // 7) {075} -2 triplets.    middle strand offset by 6, trailing by 3
    // 8) {078} -2 triplets.    middle and trailing strand offset by 6
    // Scoring will be handled by trimming the helix down to a canonical triple helix and then scoring based on this new smaller helix.
        
    double  bestPropensity;
    double  bestPairwise;
    double  Propensity[3][3][3][9];
    double  PairWise[3][3][3][9];
    double  Tm[3][3][3][9];

    
    short   XaaPos; // The position of the first Xaa amino acid.
    
    
    
    void initializeAll(void)
    {
        short a,b,c,d;
        
        numPep = 0;
        numAA = 0;
        for (a=0;a<3;a++)for(b=0;b<100;b++) sequences[a][b] = '.';
        Nterm = "initial";
        Cterm = "initial";
        //c_term = -1;
        //n_term = -1;
        expTm = 0;
        CCTm = 0;
        HighTm = 0;
        deviation = 0;
        secTm = 0;
        specificity = 0;
        for (a=0;a<4;a++)
        {
            bestRegister[a] = 0;
            secRegister[a] = 0;
            CCRegister[a] = 0;
        }
        
        bestPropensity = 0;
        bestPairwise = 0;
        
        for (a=0;a<3;a++) for(b=0;b<3;b++) for(c=0;c<3;c++) for(d=0;d<9;d++)
        {
            Propensity[a][b][c][d] = 0;
            PairWise[a][b][c][d] = 0;
            Tm[a][b][c][d] = 0;
        }
        
        XaaPos = -1;    //Should be 0, 1 or 2 depending on (POG), (GPO), (OGP) respectively
                        // Is eventually set by "determine_reptition"
        
    };
    
    bool isXaa(short position)
    {
        if (abs(position + (3 - this->XaaPos)) % 3 == 0) return true;
        return false;
    };
    
    bool isYaa(short position)
    {
        if (abs(position + (3 - this->XaaPos)) % 3 == 1) return true;
        return false;
    };
    
    bool isGly(short position)
    {
        if (abs(position + (3 - this->XaaPos)) % 3 == 2) return true;
        return false;
    };
    
    void dissect(void)
    {
        short x, y;
        cout << "numPep = " << numPep << endl;
        cout << "numAA =  " << numAA << endl;
        for (x=0;x<numPep;x++)
        {
            for (y=0;y<numAA;y++)
            {
                cout << sequences[x][y];
            }
            cout << endl;
        }
        cout << "termination: " << Nterm << " " << Cterm << endl;
        cout << "XaaPos = " << XaaPos << endl;
        cout << "expTm = " << expTm << ". CCTm = " << CCTm << ". Deviation = " << deviation << endl;
        cout << "CCregister = " << CCRegister[0] << "," << CCRegister[1] << "," << CCRegister[2] << "." << CCRegister[3] << endl;
        cout << "High Tm = " << HighTm << " = " << bestPropensity << " + " << bestPairwise << endl;
        cout << "Best register = " << bestRegister[0] << "," << bestRegister[1] << "," << bestRegister[2] << "." << bestRegister[3] <<endl;
        cout << "Second highest Tm = " << secTm << endl;
        cout << "Second Best register = " << secRegister[0] << "," << secRegister[1] << "," << secRegister[2] << "." << secRegister[3] << endl;
        cout << "Specificity = " << specificity << "." << endl;
        cout << endl;
    };
    
    void userOutput(void)
    {
        short a,b,c,x,y;
        cout << endl;
        cout << "----------------------------------------------------------" << endl;
        cout << "The most stable register/composition is {" << bestRegister[0] <<  bestRegister[1]  << bestRegister[2] << "}. Tm = " << HighTm << "." << endl;
        if (numPep == 2)
        {
            if ((bestRegister[0] != bestRegister[1]) || (bestRegister[0] != bestRegister[2]) || (bestRegister[1] != bestRegister[2]))
            {
                //cout << "This matches the input diversity of peptides (2)." << endl;
            }
            else
            {
                cout << "WARNING: The most stable register/composition does not include all the peptides you input." << endl;
            }
        }
        if (numPep == 3)
        {
            if ((bestRegister[0] != bestRegister[1]) && (bestRegister[0] != bestRegister[2]) && (bestRegister[1] != bestRegister[2]))
            {
                //cout << "This matches the input diversity of peptides (3)." << endl;
            }
            else
            {
                cout << "WARNING: The most stable register/composition does not include all the peptides you input." << endl;
            }
        }
        x = bestRegister[0];
        cout << bestRegister[0] << ": ";
        for (y=0;y<numAA;y++)
        {
            if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
            if (sequences[x][y] == 'R') cout << "\x1b[34m";
            if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
            if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
            if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
            cout << sequences[x][y];
            cout << "\x1b[0m";
        }
        cout << endl;
        x = bestRegister[1];
        cout << bestRegister[1] << ":  ";
        for (y=0;y<numAA;y++)
        {
            if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
            if (sequences[x][y] == 'R') cout << "\x1b[34m";
            if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
            if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
            if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
            cout << sequences[x][y];
            cout << "\x1b[0m";
        }
        cout << endl;
        x = bestRegister[2];
        cout << bestRegister[2] << ":   ";
        for (y=0;y<numAA;y++)
        {
            if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
            if (sequences[x][y] == 'R') cout << "\x1b[34m";
            if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
            if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
            if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
            cout << sequences[x][y];
            cout << "\x1b[0m";
        }
        cout << endl;
        
        if (numPep != 1)
        {
            cout << endl;
            cout << "The second most stable register/composition is {" << secRegister[0] << secRegister[1] << secRegister[2] << "}. Tm = " << secTm << "." << endl;
            x = secRegister[0];
            cout << secRegister[0] << ": ";
            for (y=0;y<numAA;y++)
            {
                if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
                if (sequences[x][y] == 'R') cout << "\x1b[34m";
                if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
                if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
                if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
                cout << sequences[x][y];
                cout << "\x1b[0m";
            }
            cout << endl;
            x = secRegister[1];
            cout << secRegister[1] << ":  ";
            for (y=0;y<numAA;y++)
            {
                if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
                if (sequences[x][y] == 'R') cout << "\x1b[34m";
                if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
                if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
                if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
                cout << sequences[x][y];
                cout << "\x1b[0m";
            }
            cout << endl;
            x = secRegister[2];
            cout << secRegister[2] << ":   ";
            for (y=0;y<numAA;y++)
            {
                if (sequences[x][y] == 'K') cout << "\x1b[1m\x1b[34m";
                if (sequences[x][y] == 'R') cout << "\x1b[34m";
                if ((sequences[x][y] == 'E') || (sequences[x][y] == 'D')) cout << "\x1b[1m\x1b[31m";
                if ((sequences[x][y] == 'F') || (sequences[x][y] == 'Y') || (sequences[x][y] == 'W')) cout << "\x1b[1m";
                if (sequences[x][y] == 'Q') cout << "\x1b[1m\x1b[32m";
                cout << sequences[x][y];
                cout << "\x1b[0m";
            }
            cout << endl;
            cout << endl;
            cout << "The specificity is = " << specificity << "." << endl;
            cout << endl;
        }
        
        cout << "Melting temperatures of all canonical registers." << endl;
        cout << "Best in blue, second best in red." << endl;
        cout << "Tm < 10C faded to indicate experimentally unreliable (frequently will not fold)." << endl;
        for (a=0;a<numPep;a++)
        {
            for (b=0;b<numPep;b++)
            {
                for (c=0;c<numPep;c++)
                {
                    cout << "\x1b[0m";
                    if ((a == bestRegister[0]) && (b == bestRegister[1]) && (c == bestRegister[2])) cout << "\x1b[1m\x1b[34m";
                    if ((a == secRegister[0]) && (b == secRegister[1]) && (c == secRegister[2])) cout << "\x1b[1m\x1b[31m";
                    if (Tm[a][b][c][0] < 10) cout << "\x1b[2m";
                    cout << "{" << a << b << c << "} = " << Tm[a][b][c][0] << endl;
                    cout << "\x1b[0m";
                }
            }
            cout << endl;
        }
    };
    
    void determine_reptition(void)
    {
        short x, xCount, yCount, zCount;
        bool goodPeptide;
        
        xCount = 0;
        yCount = 0;
        zCount = 0;
        goodPeptide = false;
        
        for (x=0;x<numAA;x++)
        {
            if ((x%3 == 0) && (sequences[0][x] == 'G')) xCount++;
            if ((x%3 == 1) && (sequences[0][x] == 'G')) yCount++;
            if ((x%3 == 2) && (sequences[0][x] == 'G')) zCount++;
        }
        if (xCount >= (numAA)/3)
        {
            // we are ok
            // cout << n << " Gly was found in the first position and in every subsequent i+3 position. xCount = " << xCount << endl;
            goodPeptide = true;
            XaaPos = 1;
        }
        if (yCount >= (numAA)/3)
        {
            // we are ok
            // cout << n << " Gly was found in the second position and in every subsequent i+3 position. yCount = " << yCount << endl;
            goodPeptide = true;
            XaaPos = 2;
        }
        if (zCount >= (numAA)/3)
        {
            // we are ok
            // cout << n << " Gly was found in the third position and in every subsequent i+3 position. zCount = " << zCount << endl;
            goodPeptide = true;
            XaaPos = 0;
        }
        if (not goodPeptide)
        {
            cout << "This peptide does not appear to have a Gly every third residue!" << endl;
            dissect();
            // error handling would be better if we had determine_reptition return an invalid result and kill app
        }
    };
    
};


// This will read in values from parameters.txt, parameters_exp.txt, and opt_list.txt.
parameterType ReadParameters(void);

void DisplayParameters(parameterType parameters);
void WriteParameters(parameterType parameters);

double PairWiseCalc (double XPW[], double LPW[], short currentPair, short lastPair, short previousPWType, double currentPWSum, double BestPWTotal);

TripleHelix ScoreHelix (parameterType parameters, TripleHelix theHelix);

short readLibrary (TripleHelix * Lib, string Lib_Name);
void writeLibrary(TripleHelix*Lib, int TotalOutput, string filename);


#endif /* SCEPTTr_hpp */
