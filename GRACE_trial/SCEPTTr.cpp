//
//  SCEPTTr.cpp
//  GRACE
//
//  Created by Hartgerink Lab on 11/30/24.
//

#include "SCEPTTr.hpp"

#include <vector>
#include <cmath>

//
//  SCEPTTrC 1.20
//
//  Created by Jeffrey Hartgerink, February 2022
//
//  Reads parameters from parameters.txt.



#include <iostream>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <ctime>
using namespace std;


// This will read in values from parameters.txt, parameters_exp.txt, and opt_list.txt.
parameterType ReadParameters(void)
{
    // Read parameters from file.
    parameterType parameters;
    short x, y;
    char SingleAminoAcid, parameterChar;
    short optValue;
    string StringLine;

    // Zero initial values.
    for (x=0; x<6; x++) {
        for (y=0; y<6; y++) {
            parameters.FrameShift[x][y] = 0;
        }
    }
    for (x=0;x<27;x++)
    {
        parameters.propensityX[x] = 0;
        parameters.propensityY[x] = 0;
        
        parameters.exPropensityX[x] = 0;
        parameters.exPropensityY[x] = 0;
        
        parameters.optPropX[x] = false;
        parameters.optPropY[x] = false;
    
        for (y=0;y<27;y++)
        {
            parameters.axial[x][y] = 0;
            parameters.lateral[x][y] = 0;
            
            parameters.exAxial[x][y] = 0;
            parameters.exLateral[x][y] = 0;
            
            parameters.optAxial[x][y] = false;
            parameters.optLat[x][y] = false;
        }
    }
    
    ifstream parameterFile("parameters.txt");
    if (!parameterFile.is_open()) {
        cout << "We couldn't open the parameters.txt file." << endl;
        exit(1);
    }
    
    // Reading parameter file is entirely dependant on the structure of the parameter text file being exactly accurate.
    while (!parameterFile.eof())
    {
        getline(parameterFile, StringLine);
        if (StringLine == "Length")
        {
            //cout << "we found Length" << endl;
            parameterFile >> parameters.A;
            parameterFile >> parameters.B;
            parameterFile >> parameters.C;
            
        }
        
        // // // // // //
        // added 07/01 //
        // // // // /// /
        if (StringLine == "FrameShift") {
            //cout << "we found FrameShift" << endl;
            string temp = "";
            for (x=0; x<6; x++) {
                parameterFile >> parameters.Cterm_type[x];
            }
            for (x=0; x<6; x++) {
                parameterFile >> parameters.Nterm_type[x];
                for (y=0; y<6; y++){
                    parameterFile >> parameters.FrameShift[x][y];
                }
            }
        }
        
        // /// // // // // // // //
        
        if (StringLine == "XaaPropensity")
        {
            //cout << "we found XaaPropensity" << endl;
            for (x=1;x<27;x++)
            {
                parameterFile >> SingleAminoAcid;
                parameterFile >> parameters.propensityX[(short)SingleAminoAcid-64];
            }
        }
        if (StringLine == "YaaPropensity")
        {
            //cout << "we found YaaPropensity" << endl;
            for (x=1;x<27;x++)
            {
                parameterFile >> SingleAminoAcid;
                parameterFile >> parameters.propensityY[(short)SingleAminoAcid-64];
            }
        }
        if (StringLine == "PairwiseLateral")
        {
            //cout << "we found PairwiseLateral" << endl;
            getline(parameterFile, StringLine);
            for (y=1; y<27; y++)
            {
                parameterFile >> parameterChar;
                for (x=1; x<27; x++)
                {
                    parameterFile >> parameters.lateral[y][x];
                }
            }
        }
        
        if (StringLine == "PairwiseAxial")
        {
            //cout << "we found PairwiseAxial" << endl;
            getline(parameterFile, StringLine);
            for (y=1; y<27; y++)
            {
                parameterFile >> parameterChar;
                for (x=1; x<27; x++)
                {
                    parameterFile >> parameters.axial[y][x];
                }
            }
        }
        if (StringLine == "EOF")
        {
            //cout << "we found EOF" << endl;
        }
    }
    parameterFile.close();
    

    
    
    return parameters;
}

void DisplayParameters(parameterType parameters)
{
    short x, y;
  
    cout << endl;
    cout << "OptXaaPropensity" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.optPropX[x] << endl;
    cout << "OptYaaPropensity" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.optPropY[x] << endl;
    cout << "OptPairwiseLateral" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.optLat[x][y] << "\t";
        }
        cout << endl;
    }
    
    cout << "OptPairwiseAxial" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.optAxial[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "EOF" << endl;
    cout << endl;
    cout << "----------------------------------------" << endl;
    cout << endl;
    cout << "LengthEx" << endl;
    cout << parameters.exA << endl;
    cout << parameters.exB << endl;
    cout << parameters.exC << endl;
    cout << "XaaPropensityEx" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.exPropensityX[x] << endl;
    cout << "YaaPropensityEx" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.exPropensityY[x] << endl;
    cout << "PairwiseLateralEX" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.exLateral[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "PairwiseAxial" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.exAxial[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "EOF" << endl;
    cout << endl;
    cout << "----------------------------------------" << endl;
    cout << endl;
    cout << "Length" << endl;
    cout << parameters.A << endl;
    cout << parameters.B << endl;
    cout << parameters.C << endl;
    cout << "XaaPropensity" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.propensityX[x] << endl;
    cout << "YaaPropensity" << endl;
    for (x=1;x<27;x++) cout << char(64+x) << "\t" << parameters.propensityY[x] << endl;
    cout << "PairwiseLateral" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.lateral[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "PairwiseAxial" << endl;
    for (x=1;x<27;x++) cout << "\t" << char(64+x);
    cout << endl;
    for (x=1;x<27;x++)
    {
        cout << char(64+x) << "\t";
        for (y=1;y<27;y++)
        {
            cout << parameters.axial[x][y] << "\t";
        }
        cout << endl;
    }
    cout << "EOF" << endl;
    cout << endl;
}



// This is a recursive function!
// It returns a double containing the maximum possible change in Tm from pairwise
// interactions from that start point of the call through to the end of the peptide strand.
// The final return should be the total maximum pairwise interactions from that strand.
double PairWiseCalc (double XPW[], double LPW[], short currentPair, short lastPair, short previousPWType, double currentPWSum, double BestPWTotal)
{
    bool NoPairT = false;
    bool AxPairT = false;
    bool LatPairT = false;
        
    if (currentPair > lastPair)
    {
        return BestPWTotal;
    }
    
    if ((not NoPairT) && (previousPWType != 0))
    {
        currentPWSum += 0;
        if (currentPair == lastPair)
        {
            if (currentPWSum > BestPWTotal)
            {
                BestPWTotal = currentPWSum;
            }
            return BestPWTotal;
        }
        NoPairT = true;
                
        // "We need to go deeper!"
        BestPWTotal = PairWiseCalc(XPW, LPW, (currentPair +1), lastPair, 0 /* the previous is no pair or zero */, currentPWSum, BestPWTotal);
    }
        
    if ((not LatPairT) && (previousPWType != 2))
    {
        // only sum stabilizing interactions. All possible destabilizing intereactions will be accounted for elsewhere.
        if (LPW[currentPair] > 0) currentPWSum += LPW[currentPair];
        if (currentPair == lastPair)
         {
             if (currentPWSum > BestPWTotal)
             {
                 BestPWTotal = currentPWSum;
             }
             return BestPWTotal;
         }
         LatPairT = true;
                  
     // "We need to go deeper!"
     BestPWTotal = PairWiseCalc(XPW, LPW, (currentPair +1), lastPair, 1 /* the previous is lateral or one */, currentPWSum, BestPWTotal);
    }
            
    if (not AxPairT)
    {
        // only sum stabilizing interactions. All possible destabilizing intereactions will be accounted for elsewhere.
        if (XPW[currentPair] > 0) currentPWSum += XPW[currentPair];
        if (currentPair == lastPair)
        {
            if (currentPWSum > BestPWTotal)
            {
                BestPWTotal = currentPWSum;
            }
            return BestPWTotal;
        }
        AxPairT = true;
        
        // "We need to go deeper!"
        BestPWTotal = PairWiseCalc(XPW, LPW, (currentPair +1), lastPair, 2 /* the previous is axial or two */, currentPWSum, BestPWTotal);
    }
    //cout << "terminal return" << endl;
    return BestPWTotal;
}

// Would be better to pass a pointer to parameters rather than the entire struct.
// This change should be incorporated into future version, but this works for now.
TripleHelix ScoreHelix (parameterType parameters, TripleHelix theHelix)
{
    short a, b, c, d, x, y, i;
    
    for (a=0;a<3;a++)for(b=0;b<3;b++)for(c=0;c<3;c++)for(d=0;d<9;d++)
    {
        theHelix.Propensity[a][b][c][d] = 0;
        theHelix.PairWise[a][b][c][d] = 0;
        theHelix.Tm[a][b][c][d] = 0;
        
    }
    
    double XinteractionThread[20];
    double LinteractionThread[20];
    for (x=0; x<20; x++)
    {
        XinteractionThread[x] = 0;
        LinteractionThread[x] = 0;
    }
    
    // this triple helix will hold the temporary values for all trimmed non-canonical offsets.
    TripleHelix trimmedHelix;
    
    short numYaa = 0;
    numYaa = theHelix.numAA / 3;
    
    double maxTm, secondBestTm;
    double bestReg[4], secondBestReg[4];
    double lengthBasis = 0;
    secondBestTm = -2000;
    secondBestReg[0] = 5;
    secondBestReg[1] = 5;
    secondBestReg[2] = 5;
    secondBestReg[3] = 10;
    
    maxTm = -1000;
    bestReg[0] = 6;
    bestReg[1] = 6;
    bestReg[2] = 6;
    bestReg[3] = 11;
    
    theHelix.CCTm = -1500;
    theHelix.CCRegister[0] = 7;
    theHelix.CCRegister[1] = 7;
    theHelix.CCRegister[2] = 7;
    theHelix.CCRegister[3] = 12;
    
    // NOTE: d loop is being short circuited to only look at canonical registers here! //
    for (a=0; a<theHelix.numPep; a++) for (b=0; b<theHelix.numPep; b++) for (c=0; c<theHelix.numPep; c++) for (d=0;d<1;d++)
    {
        // // // // // // // //
        // Offset Switch       //
        // // // // // // // //
        trimmedHelix.numPep = theHelix.numPep;
        trimmedHelix.Nterm = theHelix.Nterm;
        trimmedHelix.Cterm = theHelix.Cterm;
        trimmedHelix.XaaPos = theHelix.XaaPos;
        //trimmedHelix.n_term = theHelix.n_term;
        //trimmedHelix.c_term = theHelix.c_term;
        // Trim triple helix for offset scoring.
        switch(d)
        {
            case 0:
                // 0 is {012}, canonical and therefore needs no trimming.
                // Still need to move thisHelix into trimmedHelix for scoring below
                for(y=0;y<theHelix.numAA;y++) trimmedHelix.sequences[a][y] = theHelix.sequences[a][y];
                for(y=0;y<theHelix.numAA;y++) trimmedHelix.sequences[b][y] = theHelix.sequences[b][y];
                for(y=0;y<theHelix.numAA;y++) trimmedHelix.sequences[c][y] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA;
                break;
                
            case 1:
                // 1 is {015}
                // trailing strand is offset 3 and therefore needs to be cut at C term while leading and middle strands need to lose first triplet
                for(y=3;y<theHelix.numAA;y++) trimmedHelix.sequences[a][y-3] = theHelix.sequences[a][y];
                for(y=3;y<theHelix.numAA;y++) trimmedHelix.sequences[b][y-3] = theHelix.sequences[b][y];
                for(y=0;y<theHelix.numAA-3;y++) trimmedHelix.sequences[c][y] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA - 3;
                break;

            case 2:
                // 3 is {042}
                // middle strand is offset 3, others are not shifted. Cut lead and trail at N-term, mid at C term
                for(y=3;y<theHelix.numAA-0;y++) trimmedHelix.sequences[a][y-3] = theHelix.sequences[a][y];
                for(y=0;y<theHelix.numAA-3;y++) trimmedHelix.sequences[b][y-0] = theHelix.sequences[b][y];
                for(y=3;y<theHelix.numAA-0;y++) trimmedHelix.sequences[c][y-3] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA - 3;
                break;
                
            case 3:
                // 4 is {045}
                // mid and trail are offset 3. Cut lead at N-term, mid and trail at C term
                for(y=3;y<theHelix.numAA-0;y++) trimmedHelix.sequences[a][y-3] = theHelix.sequences[a][y];
                for(y=0;y<theHelix.numAA-3;y++) trimmedHelix.sequences[b][y-0] = theHelix.sequences[b][y];
                for(y=0;y<theHelix.numAA-3;y++) trimmedHelix.sequences[c][y-0] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA - 3;
                break;
                
            case 4:
                // 2 is {018}
                for(y=6;y<theHelix.numAA;y++) trimmedHelix.sequences[a][y-6] = theHelix.sequences[a][y];
                for(y=6;y<theHelix.numAA;y++) trimmedHelix.sequences[b][y-6] = theHelix.sequences[b][y];
                for(y=0;y<theHelix.numAA-6;y++) trimmedHelix.sequences[c][y] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA - 6;
                break;
                
            case 5:
                // 5 is {048}
                // mid by 3, trail by 6.
                for(y=6;y<theHelix.numAA-0;y++) trimmedHelix.sequences[a][y-6] = theHelix.sequences[a][y];
                for(y=3;y<theHelix.numAA-3;y++) trimmedHelix.sequences[b][y-3] = theHelix.sequences[b][y];
                for(y=0;y<theHelix.numAA-6;y++) trimmedHelix.sequences[c][y-0] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA - 6;
                break;
            case 6:
                // 6 is {072}
                // mid by 6, trail by 3
                for(y=6;y<theHelix.numAA-0;y++) trimmedHelix.sequences[a][y-6] = theHelix.sequences[a][y];
                for(y=0;y<theHelix.numAA-6;y++) trimmedHelix.sequences[b][y-0] = theHelix.sequences[b][y];
                for(y=3;y<theHelix.numAA-6;y++) trimmedHelix.sequences[c][y-3] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA - 6;
                break;
            case 7:
                // 7 is {075}
                // mid by 6, trail by 3
                for(y=6;y<theHelix.numAA-0;y++) trimmedHelix.sequences[a][y-6] = theHelix.sequences[a][y];
                for(y=0;y<theHelix.numAA-6;y++) trimmedHelix.sequences[b][y-0] = theHelix.sequences[b][y];
                for(y=3;y<theHelix.numAA-3;y++) trimmedHelix.sequences[c][y-3] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA - 6;
                break;
            case 8:
                // 8 is {078}
                // mid by 6, trail by 6
                for(y=6;y<theHelix.numAA-0;y++) trimmedHelix.sequences[a][y-6] = theHelix.sequences[a][y];
                for(y=0;y<theHelix.numAA-6;y++) trimmedHelix.sequences[b][y-0] = theHelix.sequences[b][y];
                for(y=0;y<theHelix.numAA-6;y++) trimmedHelix.sequences[c][y-0] = theHelix.sequences[c][y];
                trimmedHelix.numAA = theHelix.numAA - 6;
                break;
            default:
                cout << "Value of d out of bounds. d = " << d << endl;
                break;
                
        }
        
        // // // // //
        // Length   //
        // // // // //
        if (trimmedHelix.numAA >50) lengthBasis = parameters.A + (parameters.B*50) + (parameters.C*50*50);
        else lengthBasis = parameters.A + (parameters.B*trimmedHelix.numAA) + (parameters.C*theHelix.numAA*trimmedHelix.numAA);
        
        //cout << "lengthBasis = " << lengthBasis << endl;
        theHelix.Propensity[a][b][c][d] = lengthBasis;
    
    
        // // // // // //
        // TERMINATION //
        // // // // // //
        
        // // // Added 07/01 // // //
        // // // // // // // // // //
        
        // 0: capped XYG
        // 1: capped YGX
        // 2: capped GXY
        // 3: uncapped XYG
        // 4: uncapped YGX
        // 5: uncapped GXY
        
        short nterm = -1;
        short cterm = -1;
       
        if (trimmedHelix.Nterm == "ac") {
            if (trimmedHelix.isXaa(0)) {nterm = 0;}
            if (trimmedHelix.isYaa(0)) {nterm = 1;}
            if (trimmedHelix.isGly(0)) {nterm = 2;}
        }
        else {
            if (trimmedHelix.isXaa(0)) {nterm = 3;}
            if (trimmedHelix.isYaa(0)) {nterm = 4;}
            if (trimmedHelix.isGly(0)) {nterm = 5;}
        }
        
        if (trimmedHelix.Cterm == "am") {
            if (trimmedHelix.isGly(trimmedHelix.numAA-1)) {cterm = 0;}
            if (trimmedHelix.isXaa(trimmedHelix.numAA-1)) {cterm = 1;}
            if (trimmedHelix.isYaa(trimmedHelix.numAA-1)) {cterm = 2;}
        }
        else {
            if (trimmedHelix.isGly(trimmedHelix.numAA-1)) {cterm = 3;}
            if (trimmedHelix.isXaa(trimmedHelix.numAA-1)) {cterm = 4;}
            if (trimmedHelix.isYaa(trimmedHelix.numAA-1)) {cterm = 5;}
        }
        
        theHelix.Propensity[a][b][c][d] += parameters.FrameShift[nterm][cterm];
       
        //cout << "Terminal " << theHelix.Propensity[a][b][c][d] << endl;
        
        /*
        if (trimmedHelix.Nterm == "n") {
            theHelix.Propensity[a][b][c][d] -= 1.8;
        }
        if (trimmedHelix.Cterm == "c") {
            theHelix.Propensity[a][b][c][d] -= 1.8;
        }
         */
        
        // // // // // // // // //
        // Tyrosine and Tryptophan Termination //
        // // // // // // // // //
        if ((trimmedHelix.sequences[a][0] == 'Y') && (trimmedHelix.sequences[b][0] == 'Y') && (trimmedHelix.sequences[c][0] == 'Y'))
        {
            theHelix.Propensity[a][b][c][d] += 3;
        }
        if ((trimmedHelix.sequences[a][trimmedHelix.numAA-1] == 'Y') && (trimmedHelix.sequences[b][trimmedHelix.numAA-1] == 'Y') && (trimmedHelix.sequences[c][trimmedHelix.numAA-1] == 'Y'))
        {
            theHelix.Propensity[a][b][c][d] += 3;
        }
        
        if ((trimmedHelix.sequences[a][0] == 'W') && (trimmedHelix.sequences[b][0] == 'W') && (trimmedHelix.sequences[c][0] == 'W'))
        {
            theHelix.Propensity[a][b][c][d] += 2.8;
        }
        if ((trimmedHelix.sequences[a][trimmedHelix.numAA-1] == 'W') && (trimmedHelix.sequences[b][trimmedHelix.numAA-1] == 'W') && (trimmedHelix.sequences[c][trimmedHelix.numAA-1] == 'W'))
        {
            theHelix.Propensity[a][b][c][d] += 2.8;
        }
        
        /*
        //cout << "capping mod = " << Propensity[a][b][c] << endl;
        
        // // // // // // // // // // //
        // Terminal Hydrogen Bonding  //
        // // // // // // // // // // //
        // if (not theHelix.isXYG) Propensity[a][b][c][d] -= 3.6;
        if (not trimmedHelix.isXaa(0)) {
            theHelix.Propensity[a][b][c][d] -= 1.8;
        }
        if (not trimmedHelix.isGly(theHelix.numAA-1)) {
            theHelix.Propensity[a][b][c][d] -= 1.8;
        }
        
        //cout << "terminal H-bond mod = " << Propensity[a][b][c] << endl;
        */
        
        // // // // // // // //
        // Single AA Score   //
        // // // // // // // //
            
        for (x=0;x<theHelix.numAA;x++)
        {
            
            if ((x>2) && (x<(theHelix.numAA-2))) // not the tips
            {
            // cout << trimmedHelix.sequences[a][x] << " " << (short)trimmedHelix.sequences[a][x]-64 << endl;
            if (trimmedHelix.isXaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[a][x]-64];
            if (trimmedHelix.isYaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[a][x]-64];
            
            if (trimmedHelix.isXaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[b][x]-64];
            if (trimmedHelix.isYaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[b][x]-64];
            
            if (trimmedHelix.isXaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[c][x]-64];
            if (trimmedHelix.isYaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[c][x]-64];
            }
            else // the tips
            {
                if (trimmedHelix.isXaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[a][x]-64]/3;
                if (trimmedHelix.isYaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[a][x]-64]/3;
                
                if (trimmedHelix.isXaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[b][x]-64]/3;
                if (trimmedHelix.isYaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[b][x]-64]/3;
                
                if (trimmedHelix.isXaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityX[(short)trimmedHelix.sequences[c][x]-64]/3;
                if (trimmedHelix.isYaa(x)) theHelix.Propensity[a][b][c][d] += parameters.propensityY[(short)trimmedHelix.sequences[c][x]-64]/3;
            }
        }
        

        // // // // // // // //
        // set pairwise Tm   //
        // // // // // // // //
        
        // FIRST THREAD
        i = 0;
        for (x=0;x<theHelix.numAA;x++)
        {
            if (theHelix.isYaa(x))
            {
                // Populate all values for the First Interaction Thread (both axial and lateral options)
                // a, b & c are the peptide number of this particular composition / registration
                // i is tracking the number of Yaa's
                // x is tracking the amino acid position of the peptide
                if ((x+2) < trimmedHelix.numAA) XinteractionThread[i] = parameters.axial[(short)trimmedHelix.sequences[a][x]-64][(short)trimmedHelix.sequences[b][x+2]-64];
                    else XinteractionThread[i] = 0;
                if ((x-1) >= 0) LinteractionThread[i] = parameters.lateral[(short)trimmedHelix.sequences[a][x]-64][(short)trimmedHelix.sequences[b][x-1]-64];
                    else LinteractionThread[i] = 0;
                i++;
            }
        }
                    
        // find best combination of stabilizing interactions
        theHelix.PairWise[a][b][c][d] = PairWiseCalc (XinteractionThread, LinteractionThread, 0, numYaa, 9, 0, 0);
        
        // force *ALL* possible destabilizing interactions
        for (x=0;x<numYaa;x++)
        {
            //cout << "Axial Interaction Thread = " << XinteractionThread[i] << " ";
            if (XinteractionThread[x] < 0) theHelix.PairWise[a][b][c][d] += XinteractionThread[x];
            if (LinteractionThread[x] < 0) theHelix.PairWise[a][b][c][d] += LinteractionThread[x];
        }
        
        // SECOND THREAD
        i = 0;
        for (x=0;x<trimmedHelix.numAA;x++)
        {
            if (trimmedHelix.isYaa(x))
            {
                // Populate all values for the Second Interaction Thread (both axial and lateral options)
                if ((x+2) < theHelix.numAA)
                    XinteractionThread[i] = parameters.axial[(short)trimmedHelix.sequences[b][x]-64][(short)trimmedHelix.sequences[c][x+2]-64];
                else
                    XinteractionThread[i] = 0;
                if ((x-1) >= 0)
                    LinteractionThread[i] = parameters.lateral[(short)trimmedHelix.sequences[b][x]-64][(short)trimmedHelix.sequences[c][x-1]-64];
                else
                    LinteractionThread[i] = 0;
                i++;
            }
        }
        
        // find best combination of stabilizing interactions
        theHelix.PairWise[a][b][c][d] += PairWiseCalc (XinteractionThread, LinteractionThread, 0, numYaa, 9, 0, 0);
        
        // force *ALL* possible destabilizing interactions
        for (x=0;x<numYaa;x++)
        {
            if (XinteractionThread[x] < 0) theHelix.PairWise[a][b][c][d] += XinteractionThread[x];
            if (LinteractionThread[x] < 0) theHelix.PairWise[a][b][c][d] += LinteractionThread[x];
        }
        
        // THIRD THREAD
        i = 0;
        for (x=0;x<trimmedHelix.numAA;x++)
        {
            if (trimmedHelix.isYaa(x))
            {
                // Populate all values for the Third Interaction Thread (both axial and lateral options)
                if ((x+5) < theHelix.numAA)
                    XinteractionThread[i] = parameters.axial[(short)trimmedHelix.sequences[c][x]-64][(short)trimmedHelix.sequences[a][x+5]-64];
                else
                    XinteractionThread[i] = 0;
                if ((x+2) < theHelix.numAA)
                    LinteractionThread[i] = parameters.lateral[(short)trimmedHelix.sequences[c][x]-64][(short)trimmedHelix.sequences[a][x+2]-64];
                else
                    LinteractionThread[i] = 0;
                i++;
            }
        }
        
        // find best combination of stabilizing interactions
        theHelix.PairWise[a][b][c][d] += PairWiseCalc (XinteractionThread, LinteractionThread, 0, numYaa, 9, 0, 0);
        // force *ALL* possible destabilizing interactions
        for (x=0;x<numYaa;x++)
        {
            if (XinteractionThread[x] < 0) theHelix.PairWise[a][b][c][d] += XinteractionThread[x];
            if (LinteractionThread[x] < 0) theHelix.PairWise[a][b][c][d] += LinteractionThread[x];
        }
        
        //cout << "Pairwise Mod = " << PairWise[a][b][c] << endl;
        
        // Calculate the Tm for this composition / register determined.
        theHelix.Tm[a][b][c][d] = theHelix.Propensity[a][b][c][d] + theHelix.PairWise[a][b][c][d];
        
        //cout << "Tm = " << Tm[a][b][c] << " = " << Propensity[a][b][c] << " + " << PairWise[a][b][c] << endl;
        
        // Determine if this is the best and remember the best & second best registers.
        if (theHelix.Tm[a][b][c][d] >= maxTm)
        {
            secondBestTm = maxTm;
            secondBestReg[0] = bestReg[0];
            secondBestReg[1] = bestReg[1];
            secondBestReg[2] = bestReg[2];
            secondBestReg[3] = bestReg[3];
            
            maxTm = theHelix.Tm[a][b][c][d];
            bestReg[0] = a;
            bestReg[1] = b;
            bestReg[2] = c;
            bestReg[3] = d;
        }
        else
        {
            if (theHelix.Tm[a][b][c][d] >= secondBestTm)
            {
                secondBestTm = theHelix.Tm[a][b][c][d];
                secondBestReg[0] = a;
                secondBestReg[1] = b;
                secondBestReg[2] = c;
                secondBestReg[3] = d;
            }
        }
        
        if (theHelix.numPep == 1)
        {
            if (theHelix.Tm[a][b][c][d] >= theHelix.CCTm)
            {
                //cout << "Setting CCTm for homotrimer" << endl;
                theHelix.CCTm = theHelix.Tm[a][b][c][d];
                theHelix.CCRegister[0] = a;
                theHelix.CCRegister[1] = b;
                theHelix.CCRegister[2] = c;
            }
        }
        
        if ((theHelix.numPep == 2) && ((a != b) || (a != c) || (b != c)))
        {
            
            if (theHelix.Tm[a][b][c][d] >= theHelix.CCTm)
            {
                //cout << "This is an A2B register: " << n << " " <<bestReg[0]<<bestReg[1]<<bestReg[2] << endl;
                theHelix.CCTm = theHelix.Tm[a][b][c][d];
                theHelix.CCRegister[0] = a;
                theHelix.CCRegister[1] = b;
                theHelix.CCRegister[2] = c;
                theHelix.CCRegister[3] = d;
            }
        }
            
        if ((theHelix.numPep == 3) && ((a != b) && (a != c) && (b != c)))
        {
            
            if (theHelix.Tm[a][b][c][d] >= theHelix.CCTm)
            {
                //cout << "This is an ABC register: " << n << " " <<bestReg[0]<<bestReg[1]<<bestReg[2] << endl;
                theHelix.CCTm = theHelix.Tm[a][b][c][d];
                theHelix.CCRegister[0] = a;
                theHelix.CCRegister[1] = b;
                theHelix.CCRegister[2] = c;
                theHelix.CCRegister[3] = d;
            }
        }
        
        
    } // end for abcd
    
    
    // move all values to the helix parameters...
    
    
    theHelix.HighTm = maxTm;
    theHelix.bestRegister[0] = bestReg[0];
    theHelix.bestRegister[1] = bestReg[1];
    theHelix.bestRegister[2] = bestReg[2];
    theHelix.bestRegister[3] = bestReg[3];
    theHelix.bestPropensity = theHelix.Propensity[theHelix.bestRegister[0]][theHelix.bestRegister[1]][theHelix.bestRegister[2]][theHelix.bestRegister[3]];
    theHelix.bestPairwise = theHelix.PairWise[theHelix.bestRegister[0]][theHelix.bestRegister[1]][theHelix.bestRegister[2]][theHelix.bestRegister[3]];
    theHelix.secTm = secondBestTm;
    theHelix.secRegister[0] = secondBestReg[0];
    theHelix.secRegister[1] = secondBestReg[1];
    theHelix.secRegister[2] = secondBestReg[2];
    theHelix.secRegister[3] = secondBestReg[3];
    theHelix.specificity = maxTm - secondBestTm;
    theHelix.deviation = theHelix.CCTm - theHelix.expTm;
    
    
    // // /// /// /// //

    // an experimental Tm of -10 indicates that no transition was observed experimentally. Only penalize Tm's higher than 10C.
    if (theHelix.expTm == -10)
    {
        if (theHelix.CCTm <= 10) theHelix.deviation = 0; // if it is below 10 it is effectively unknown and good enough.
        else theHelix.deviation = theHelix.CCTm - 10; // if it is above 10, penalize it for each point above 10.
    }
    else
        theHelix.deviation = theHelix.CCTm - theHelix.expTm;
        
    return theHelix;
}


