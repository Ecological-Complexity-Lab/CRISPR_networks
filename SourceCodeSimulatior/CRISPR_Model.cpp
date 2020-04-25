
//  Created by Tong Wang on 5/24/18.
//  Modified by Sergio A. Alcala-Corona during 2018-2019
//  Creative Commons License: CC BY-NC-SA 4.0 (2019)

#include <iostream>
#include <math.h>
#include <stdio.h>
// #include <cstdlib>
#include <cmath>
#include <fstream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <random>
#include <iomanip>      // std::setprecision

#include <bits/stdc++.h>

using namespace std;

// ##########################################################################################################################
// ########################################## PDI Functions  ################################################################
// ##########################################################################################################################


// This function returns the union of two sets (of ProtoSpacers and Spacers) to compute PDI.
vector<int> Union(int arr1[], int arr2[], int n1, int n2)
{
    set<int> hs;
 
    for (int i = 0; i < n1; i++){
        if (arr1[i] != -1)
            hs.insert(arr1[i]);
    }
 
    for (int i = 0; i < n2; i++){
        if (arr1[i] != -1)
            hs.insert(arr2[i]);
    }
 
    vector<int> v(hs.size());
    copy(hs.begin(), hs.end(), v.begin());
    
    return v;    
}
 

// This function returns the intersection of two sets (of ProtoSpacers and Spacers) to compute PDI.
vector<int> Intersection(int arr1[], int arr2[],int n1, int n2) {
    
    set<int> hs;
    vector<int> v;
    
    for (int i = 0; i < n1; i++){
        if (arr1[i] != -1)
            hs.insert(arr1[i]);
    }
        

    for (int i = 0; i < n2; i++)
        if (hs.find(arr2[i]) != hs.end()){
            if (arr2[i] != -1){
                v.push_back(arr2[i]);
            }
        }
    
    return v;   
}

// This function computes the sigma term (based on the modified version) for the  PDI.
double Sigma(int SpacerSet1[], int SpacerSet2[], int ProtospacerSet[], int NumOfSpacers, int NumOfProtospacers){
    
    double S=0;
    double J=0;
    int Cn=1;
    int Cu=1;
    vector<int> R1;
    vector<int> R2;
    vector<int> UN;
    vector<int> IN;

   
    R1 = Intersection(SpacerSet1, ProtospacerSet, NumOfSpacers, NumOfProtospacers);
    R2 = Intersection(SpacerSet2, ProtospacerSet, NumOfSpacers, NumOfProtospacers);
    
    
    int R1a[R1.size()];
    int R2a[R2.size()];    
    
    copy(R1.begin(), R1.end(), R1a);
    copy(R2.begin(), R2.end(), R2a);
    
    
    int r1 = R1.size();
    int r2 = R2.size();
    

    if (R1.size() == 0 || R2.size()== 0){
        Cn=1;
        Cu=1;
    } else {        
        IN=Intersection(R1a,R2a,r1,r2);
        UN=Union(R1a,R2a,r1,r2);
        Cn=IN.size();
        Cu=UN.size();
    }
        
   
        
    J = (float)Cn/(float)Cu;
    S = (1 - J);
    
    return S;
    
}

// ##########################################################################################################################
// ##########################################################################################################################



// Set up an uniform distribution between 0 and 1 and fix the seed for the random generator:
// std::mt19937 generator(1729);
// std::uniform_real_distribution<double> distribution(0.0,1.0);

// CRISPRImmuneFunction is the function checking whether oneSpacerSet[] and oneProtospacerSet[] shares elements:
int CRISPRImmuneFunction(int oneSpacerSet[],int oneProtospacerSet[], int NumOfSpacers, int NumOfProtospacers); 

// SpacerProtospacerOverlap is the function checking how many elements are shared by oneSpacerSet[] and oneProtospacerSet[]:
int SpacerProtospacerOverlap(int oneSpacerSet[],int oneProtospacerSet[], int NumOfSpacers, int NumOfProtospacers);


// int AllPAM_mutated(int oneProtospacerSet[], int NumOfProtospacers);  This not matters

int main(int argc, const char * argv[]) {
    
    int Seed = stoi(argv[6]);    // Here we collect the seed in order to run different stochastic simulations.
 
    
    std::mt19937 generator(Seed);
    std::uniform_real_distribution<double> distribution(0.0,1.0);    
    
   
    string stringSeed = to_string(Seed);
    
    // Parameters (can be compared with Childs et.al. model. The names of those parameters should be self-explanatory):
    
    double CRISPRfailure = 1.0e-5;           // CRISPR failure rate even if the bacteria have the resistance.  //p in Childs et.al. model
    double SpacerAcquisition = 1.0e-5;       // Spacer acquisition rate per bacterium phage encounter.          //q in Childs et.al. model
    double carryingCapacity = 1e5*sqrt(10);   // Carrying capacity for the bacterial population.                // K in Childs et.al. model
    int burstSize = 50;                     // burst size for phages.                                           // beta in Childs et.al. model
    double adsorptionRate = 1.0e-7;          // adsorption rate.                                                // phi in Childs et.al. model
    double decayRate = 0.1;                 // decay rate for phages.                                           \\ m in Childs et.al. model
    double MutationRate = 5e-7;             // mutation rate for phages per phage generation per protospacer.   \\ mu in thChilds et.al. model
    string stringMutationRate = "5e-7";     // A string stringMutationRate only used for generated file name.
    stringMutationRate = argv[2];
    MutationRate = stod(stringMutationRate);
    double densityCutoff = 0.1;            // densityCutoff is the lower bound below which a single bacterial strain and phage strain is considered to be extinct.  \\ \pho_c
    
    int NumOfProtospacers = 10;           // Number of protospacer.                                               
    int NumOfSpacers = 8;                 // Number of spacers.
    
    string stringNumOfProtospacers = argv[4];
    string stringNumOfSpacers = argv[5];
    NumOfProtospacers = stoi(stringNumOfProtospacers);          // Number of protospacer.                                               
    NumOfSpacers = stoi(stringNumOfSpacers);                 // Number of spacers.        
    
    
    // Because the protospacers and spacers are written as numbers 0,1,2... So the newest mutation of a protospacer will simply generate a new number. LargestIntegerUsed is the largest number used in the protospacer sets.
    int LargestIntegerUsed = NumOfProtospacers-1;
    double timeInterval = 1.0/100;          // time steps for calculating the ODEs in units of "hour".   // DT  for Euler
    double time = 0.0;                     // time is used to save the time in the coevolutionary simulation.  // Reord initial time
    int timesOfRecord = 1;                // Not the information of every single time point is written to the file. Only information every 1hour is written to the file in the simulation.
    double nextTimeToWriteFile=0;         // nextTimeToWriteFile defines the next time to write information to the file. nextTimeToWriteFile is added a fixed time separation like 1hour to make sure the information is written to the file every 1 hour.
    
    
    int runningTime=10000;                  // controlling the running time as hours for the simulation.      
    string stringrunningTime = argv[3];
    runningTime = stoi(stringrunningTime);
    double evolulationRate = 0.0;
    double BacteriaEvolutionRate = 0.0;    // total bacterial evolutionary rate.
    double PhageEvolutionRate = 0.0;       // total phage evolutionary rate.
    double sumOfBacteria = 0;             // summation of all bacterial population.
    double sumOfPhage = 0;                // summation of all phage strains.
    double effDb;                         // effective diversity for bacteria denoted by 1/(\sum(p_i^2)) where p_i is the proportion of strains.
    double effDp;                         // effective diversity for phage denoted by 1/(\sum(p_i^2)) where p_i is the proportion of strains.
    
    
    int Dp = 1;                           // Db is the number of phage strains existing. For a cold start, Dp = 1.
    int Db = 1;                           // Db is the number of bacterial strains existing.
    
#define count 50000
    double* Bold = new double[count];         // Bold saves the bacterial abundance.   Number of Bacteria per millilitre
    double* Pold = new double[count];         // Pold saves the phage abundance.        Number of Phages per millilitre
    double* Bnew = new double[count];         // Bnew temporarily saves the bacterial abundance calculated by ODEs.
    double* Pnew = new double[count];         // Pnew temporarily saves the phage abundance calculated by ODEs.
    int* Blabel = new int[count];            // Blablel saves the labelling for bacterial strains.
    int* Plabel = new int[count];            // Plablel saves the labelling for phage strains.
    int BlabelLargest=0;                   // BlabelLargest saves the largest used labelling for bacterial strains.
    int PlabelLargest=0;                   // PlabelLargest saves the largest used labelling for phage strains.
   
    int** Spacers = new int*[count];         // Spacers save the spacer set for each corresponded bacterial strain.
    for(int i=0;i<count;i++){
        Spacers[i] = new int[NumOfSpacers];
    }
    
    
    int** Protospacers = new int*[count];    // Protospacers save the spacer set for each corresponded phage strain.
    for(int i=0;i<count;i++){
        Protospacers[i] = new int[NumOfProtospacers];
    }
    
    double* SpacerAcquisitionRate = new double[count]; //SpacerAcquisitionRate is the calculated spacer acquisition rate for each bacterial strain.    
    double* ProtospacerMutationRate = new double[count]; //ProtospacerMutationRate is the calculated mutation rate for each phage strain.
    
    
    int** infectionMatrix = new int*[count];
    for(int i=0;i<count;i++){
        infectionMatrix[i] = new int[count];      // infectionMatrix is the infection matrix/network between bacterial strains and phages. The entry of the matrix equals 1/0 if the phage strain can/cannot infect the bacterial strain.
    }
    
   
    
    // Initial conditions:
    time = 0.0;
    timesOfRecord = 1;
    nextTimeToWriteFile=0;

    string stringDp = argv[1];          
    Dp = stod(stringDp);
    
    Db = 1;                             
    for(int j=0;j<Db;j++){
        Bold[j] = 1e5;                // Assign the initial bacterial population.
        Blabel[j] = j+1;             // Assign a unique positive number to each strain.
        BlabelLargest = j+1;         // BlabelLargest is the largest label number used by the bacterial population.
        for(int i=0;i<NumOfSpacers;i++){
            Spacers[j][i] = -1;      // spacer set initially is set as [-1,-1,...,-1]/
        }
    }
    for(int j=0;j<Dp;j++){
        Pold[j] = 1e5;               // Assign the initial phage population.
        Plabel[j] = j+1;            // Assign a unique positive number to each strain.
        PlabelLargest = j+1;        // PlabelLargest is the largest label number used by the phage population.
        for(int i=0;i<NumOfProtospacers;i++){
            Protospacers[j][i] = i + j*NumOfProtospacers;      // Protospacer set initially is set as [0,1,..,9], [10,11,...,19]... to make sure there is no matching between bacterial spacer set and phage protospacer sets.
            LargestIntegerUsed = Protospacers[j][i];           // LargestIntegerUsed saves the largest integer used for protospace sets.
        }
    }
    
    
    
    
// #############################################################################################################################################################   
    
    
// Write the file name: (You need to change the path according to your system.)

    
    ofstream fileForBacteria("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_data-bact"+".txt");     // data1 records data of population from 1000h to 2500h per 1h
    ofstream fileForPhage("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_data-phage"+".txt");
    ofstream fileForTimeSeriesData("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_time-series-data"+".txt");
    ofstream fileForBacteriaAbundance("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_Bacteria-abundance"+".txt");
    ofstream fileForBacteriaAnimation("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_Bacteria-animation"+".txt");
    ofstream fileForPhageAbundance("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_Phage-abundance"+".txt");
    ofstream fileForPhageAnimation("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_Phage-animation"+".txt");
    
    ofstream fileForPhageTree("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_Phage-TREE"+".txt"); //   intento de Arbol!!!!!
    ofstream fileForBacteriaTree("mu"+stringMutationRate+"_initialDiffDp"+stringDp+"_S"+to_string(NumOfSpacers)+"P"+to_string(NumOfProtospacers)+"_R-"+stringSeed+"_Bacteria-TREE"+".txt"); //   intento de Arbol!!!!!
        

    
    fileForBacteria << "timesOfRecord time label density";
    for(int i=0;i<NumOfSpacers;i++){
        fileForBacteria << " Spacer" << to_string(i+1);
    }

    fileForBacteria << endl;
    fileForPhage << "timesOfRecord time label density";
    for(int i=0;i<NumOfProtospacers;i++){
        fileForPhage << " Protospacer" << to_string(i+1);
    }
    fileForPhage << endl;
//     fileForTimeSeriesData << "timesOfRecord time Db Dp Bdensity Pdensity effDb effDp" << endl;
    fileForTimeSeriesData << "timesOfRecord time Db Dp Bdensity Pdensity effDb effDp PDI" << endl;

    fileForBacteriaAbundance << "timesOfRecord time label Bdensity ";// << endl;
    fileForBacteriaAbundance << "SpacersOccupied"  << endl;  // Aqui Agregue esta linea para cotar los spacers ocupados!!!!!!!
    
    fileForBacteriaAnimation << "timesOfRecord time label Bdensity Bfit" << endl;
    fileForPhageAbundance << "timesOfRecord time label Pdensity" << endl;
    fileForPhageAnimation << "timesOfRecord time label Pdensity Pfit" << endl;
    
// Code for generate phylogeny tree files

    for(int j=0;j<Dp;j++){
        fileForPhageTree << timesOfRecord << "\t" << Plabel[j] << "\t" << "NA" << "\t" << time << endl;     
    }  
    
    
    
// ################################################################################################################################################################################ 
    
    
    // Start the simulation:
    while(time<runningTime){     // Continue the simulation when the time is less than the ruuningTime.
        // Eliminate the strains below the densityCutoff (i.e. keep the strains above the densityCutoff):
        int tempDp = 0;
        int tempDb = 0;
        effDb = 0;
        effDp = 0;
        
     
        
        // ########### caulculation of Effective Diversity   ###############
        
        for(int i=0;i<Db;i++){            
            sumOfBacteria = sumOfBacteria + Bold[i];
            effDb = effDb + pow(Bold[i],2);
        }
        
        
        effDb = pow(sumOfBacteria,2) / effDb;
                
        for(int i=0;i<Dp;i++){
            sumOfPhage = sumOfPhage + Pold[i];
            effDp = effDp + pow(Pold[i],2);
        }
        effDp = pow(sumOfPhage,2) / effDp;
        
        
        
        // ############ extinction part #######################
        
        for(int i=0;i<Db;i++){
            if(Bold[i]>densityCutoff){
                Bold[tempDb] = Bold[i];
                Blabel[tempDb] = Blabel[i];
                for(int j=0;j<NumOfSpacers;j++){
                    Spacers[tempDb][j] = Spacers[i][j];
                }
                tempDb++;
            }
        }
        for(int i=0;i<Dp;i++){
            if(Pold[i]>densityCutoff){
                Pold[tempDp] = Pold[i];
                Plabel[tempDp] = Plabel[i];
                for(int j=0;j<NumOfProtospacers;j++){
                    Protospacers[tempDp][j] = Protospacers[i][j];
                }
                tempDp++;
            }
            
        }
        Db = tempDb;
        Dp = tempDp;
        
        // ##################################################################
        
        
        // Write the information out when time is larger than nextTimeToWriteFile. nextTimeToWriteFile is added by 1 hour each time. So the information is written out every 1 hour.
        if(time > nextTimeToWriteFile){
            // Calculate the infection_matrix:
            for(int i=0;i<Db;i++){
                for(int j=0;j<Dp;j++){
                    infectionMatrix[i][j] = 1.0 - CRISPRImmuneFunction(Spacers[i], Protospacers[j], NumOfSpacers, NumOfProtospacers);
                }
            }
//              ##################################################################!!!!!!!!!!
            for(int i=0;i<Db;i++){
                fileForBacteria << timesOfRecord << " " << time << " " << Blabel[i] << " " << Bold[i];
                for(int j=0;j<NumOfSpacers;j++){
                    fileForBacteria << " " << Spacers[i][j];
/                }
                fileForBacteria << endl;
            }
            
            
            for(int i=0;i<Dp;i++){
                fileForPhage << timesOfRecord << " " << time << " " << Plabel[i] << " " << Pold[i];
                                for(int j=0;j<NumOfProtospacers;j++){
                    fileForPhage << " " << Protospacers[i][j];
                }
                fileForPhage << endl;
            }
            



            double Ni = 0;       
            double Nj = 0;       
            double Vk = 0; 
            double PDI = 0;     
            double Nmax=0;
            double sumOfBacteria2 = 0;
            double sumOfPhage2 = 0;

            
            for(int i=0;i<Db;i++){
                sumOfBacteria2 = sumOfBacteria2 + Bold[i];
                if(Bold[i]>Nmax){
                    Nmax=Bold[i];
                }
            }
            
            Nmax=Nmax/sumOfBacteria2;
            
            for(int i=0;i<Dp;i++){
                sumOfPhage2 = sumOfPhage2 + Pold[i];
            }            

            

            for(int i=0;i<Db;i++){
                for(int j=0;j<Db;j++){
                    for(int k=0;k<Dp;k++){
                        Ni = Bold[i]/sumOfBacteria2;
                        Nj = Bold[j]/sumOfBacteria2;
                        Vk = Pold[k]/sumOfPhage2;
                        

//   ################### PDI computation  #########################                       
                        
                        if (Blabel[i] != Blabel[j]) {

                            double S = Sigma(Spacers[i],Spacers[j], Protospacers[k], NumOfSpacers, NumOfProtospacers);
                            PDI = PDI + (1-((Ni-Nj)/Nmax))*S*(Ni*Nj*Vk);
                            
                        }                        
                    }                    
                }                
            }
            
            fileForTimeSeriesData << timesOfRecord << " " << time << " " << Db << " " << Dp << " " << sumOfBacteria2 << " " << sumOfPhage2 << " " << effDb << " " << effDp << " " << PDI <<  endl;
            

            for(int i=0;i<Db;i++){  // Abundant bacteria population
                int tmpSPcounter = 0; // number of used spacers
                for(int j=0;j<NumOfSpacers;j++){
                    if (Spacers[i][j] != -1) {tmpSPcounter++;} // Aqui Agregue esta linea para cotar los spacers ocupados!!!!!!!
                }                    
                    
                    fileForBacteriaAbundance << timesOfRecord << " " << time << " " << Blabel[i] << " " << Bold[i] << " " << tmpSPcounter << endl;
            }
            
            
            double Bfit;
            for(int i=0;i<Db;i++){  // bacteria animation
                Bfit = 0;
                for(int j=0;j<Dp;j++){
                    Bfit = Bfit + (1.0 - infectionMatrix[i][j]) * Pold[j]/sumOfPhage;
                }
                fileForBacteriaAnimation << timesOfRecord << " " << time << " " << Blabel[i] << " " << Bold[i] << " " << Bfit << endl;
            }
            
       
            for(int i=0;i<Dp;i++){  // Abundant phage population
                    fileForPhageAbundance << timesOfRecord << " " << time << " " << Plabel[i] << " " << Pold[i] << endl;
              
            }
            double Pfit;
            for(int i=0;i<Dp;i++){  // phage animation
                Pfit = 0;
                for(int j=0;j<Db;j++){
                    Pfit = Pfit + infectionMatrix[j][i] * Bold[j]/sumOfBacteria;
                }
                fileForPhageAnimation << timesOfRecord << " " << time << " " << Plabel[i] << " " << Pold[i] << " " << Pfit << endl;
            }
            
            timesOfRecord = timesOfRecord+1;
            nextTimeToWriteFile = nextTimeToWriteFile + 1.0;
        }
        
        
        // Check when the next mutation event happens:
        evolulationRate = 0.0;
        BacteriaEvolutionRate = 0.0;
        PhageEvolutionRate = 0.0;
        sumOfPhage = 0;
        sumOfBacteria = 0;
        for(int i=0;i<Dp;i++){
            sumOfPhage = sumOfPhage + Pold[i];
        }
        for(int i=0;i<Db;i++){
            sumOfBacteria = sumOfBacteria + Bold[i];
        }
        // Calculate the infection_matrix:
        for(int i=0;i<Db;i++){
            for(int j=0;j<Dp;j++){
                infectionMatrix[i][j] = 1.0 - CRISPRImmuneFunction(Spacers[i], Protospacers[j], NumOfSpacers, NumOfProtospacers);
            }
        }
        
        // Bacteria evolution rate:
        for(int i=0;i<Db;i++){
            SpacerAcquisitionRate[i] = SpacerAcquisition * adsorptionRate * Bold[i] * sumOfPhage;
            evolulationRate = evolulationRate + SpacerAcquisitionRate[i];
            BacteriaEvolutionRate = BacteriaEvolutionRate + SpacerAcquisitionRate[i];
        }
        // Phage evolution rate:
        for(int i=0;i<Dp;i++){
            ProtospacerMutationRate[i] = 0.0;
            for(int j=0;j<Db;j++){
                ProtospacerMutationRate[i] = ProtospacerMutationRate[i] + NumOfProtospacers * (1.0 - SpacerAcquisition) * MutationRate * burstSize * adsorptionRate * Pold[i] * Bold[j] * infectionMatrix[j][i] + NumOfProtospacers * CRISPRfailure * MutationRate * burstSize * adsorptionRate * Pold[i] * Bold[j] * (1.0 - infectionMatrix[j][i]);
            }
            evolulationRate = evolulationRate + ProtospacerMutationRate[i];
            PhageEvolutionRate = PhageEvolutionRate + ProtospacerMutationRate[i];
        }
        // Use Gillespie algorithm to calculate the time to the next evolutionary event:
        double timeToEvolution = -(1.0/evolulationRate) * log(distribution(generator));

// #########################################################################################################################################     
        
        // Calculate the ODEs before the evolution:
        sumOfBacteria = 0.0;
        for(int i=0;i<Db;i++){
            sumOfBacteria = sumOfBacteria + Bold[i];
        }
        if(timeToEvolution<=0.1){     // If the timeToEvolution(time to the next evolutionary event) is less than 0.1, timeToEvolution is used as the time interval to calculate the ODEs.
            timeInterval = timeToEvolution;
            // Calculate the change of bacteria population:
            for(int i=0;i<Db;i++){
                Bnew[i] = Bold[i] + timeInterval * (1-sumOfBacteria/carryingCapacity) * Bold[i];
                for(int j=0;j<Dp;j++){
                    Bnew[i] = Bnew[i] - timeInterval * ((1.0 - SpacerAcquisition) * adsorptionRate * Bold[i] * Pold[j] * infectionMatrix[i][j] - CRISPRfailure * adsorptionRate * Bold[i] * Pold[j] * (1.0 - infectionMatrix[i][j]));
                }
            }
            // Calculate the change of phage population:
            for(int i=0;i<Dp;i++){
                Pnew[i] = Pold[i] - timeInterval * decayRate * Pold[i];
                for(int j=0;j<Db;j++){
                    Pnew[i] = Pnew[i] + timeInterval * ((1.0 - SpacerAcquisition) * burstSize * adsorptionRate * Pold[i] * Bold[j] * infectionMatrix[j][i] + CRISPRfailure * burstSize * adsorptionRate * Pold[i] * Bold[j] * (1.0 - infectionMatrix[j][i]) - adsorptionRate * Pold[i] * Bold[j]);
                }
            }
            time = time + timeInterval;
            for(int i=0;i<Db;i++){
                Bold[i] = Bnew[i];            // Transfer the value from Bnew back to Bold. Bnew is just used for temporary saving.
            }
            for(int i=0;i<Dp;i++){
                Pold[i] = Pnew[i];            // Transfer the value from Pnew back to Pold. Pnew is just used for temporary saving.
            }
        }
        else if(timeToEvolution<=1){     // If the timeToEvolution(time to the next evolutionary event) is less than 0.1, timeInterval=0.1 is used as the time interval to calculate the ODEs.
            int CalcTimes = floor(timeToEvolution/0.1) + 1;
            timeInterval = timeToEvolution/CalcTimes;
            for(int k=0;k<CalcTimes;k++){
                // Calculate the change of bacteria population:
                for(int i=0;i<Db;i++){
                    Bnew[i] = Bold[i] + timeInterval * (1-sumOfBacteria/carryingCapacity) * Bold[i];
                    for(int j=0;j<Dp;j++){
                        Bnew[i] = Bnew[i] - timeInterval * ((1.0 - SpacerAcquisition) * adsorptionRate * Bold[i] * Pold[j] * infectionMatrix[i][j] - CRISPRfailure * adsorptionRate * Bold[i] * Pold[j] * (1.0 - infectionMatrix[i][j]));
                    }
                }
                // Calculate the change of phage population:
                for(int i=0;i<Dp;i++){
                    Pnew[i] = Pold[i] - timeInterval * decayRate * Pold[i];
                    for(int j=0;j<Db;j++){
                        Pnew[i] = Pnew[i] + timeInterval * ((1.0 - SpacerAcquisition) * burstSize * adsorptionRate * Pold[i] * Bold[j] * infectionMatrix[j][i] + CRISPRfailure * burstSize * adsorptionRate * Pold[i] * Bold[j] * (1.0 - infectionMatrix[j][i]) - adsorptionRate * Pold[i] * Bold[j]);
                    }
                }
                time = time + timeInterval;
                for(int i=0;i<Db;i++){
                    Bold[i] = Bnew[i];
                }
                for(int i=0;i<Dp;i++){
                    Pold[i] = Pnew[i];
                }
            }
        }
        else{          // If the timeToEvolution(time to the next evolutionary event) is less than 0.1, timeInterval=0.1 is used as the time interval to calculate the ODEs.
            int CalcTimes = 10;
            timeInterval = 0.1;
            for(int k=0;k<CalcTimes;k++){
                // Calculate the change of bacteria population:
                for(int i=0;i<Db;i++){
                    Bnew[i] = Bold[i] + timeInterval * (1-sumOfBacteria/carryingCapacity) * Bold[i];
                    for(int j=0;j<Dp;j++){
                        Bnew[i] = Bnew[i] - timeInterval * ((1.0 - SpacerAcquisition) * adsorptionRate * Bold[i] * Pold[j] * infectionMatrix[i][j] - CRISPRfailure * adsorptionRate * Bold[i] * Pold[j] * (1.0 - infectionMatrix[i][j]));
                    }
                }
                // Calculate the change of phage population:
                for(int i=0;i<Dp;i++){
                    Pnew[i] = Pold[i] - timeInterval * decayRate * Pold[i];
                    for(int j=0;j<Db;j++){
                        Pnew[i] = Pnew[i] + timeInterval * ((1.0 - SpacerAcquisition) * burstSize * adsorptionRate * Pold[i] * Bold[j] * infectionMatrix[j][i] + CRISPRfailure * burstSize * adsorptionRate * Pold[i] * Bold[j] * (1.0 - infectionMatrix[j][i]) - adsorptionRate * Pold[i] * Bold[j]);
                    }
                }
                time = time + timeInterval;
                for(int i=0;i<Db;i++){
                    Bold[i] = Bnew[i];
                }
                for(int i=0;i<Dp;i++){
                    Pold[i] = Pnew[i];
                }
            }
        }
        
        
        // Perform the evolution event only when the timeToEvolution is less than 1 hour:
        if(timeToEvolution<=1){
            double randomNumber = distribution(generator);
            // Bacteria evolve (i.e. acquire a spacer):
            if(BacteriaEvolutionRate/evolulationRate>=randomNumber){
                double randomNumber1 = distribution(generator);   // randomNumber1 is used to determine which bacterial strain will acquire a spacer.
                double cumulative = 0.0;
                int needToBreak = 0;
                for(int i=0;i<Db;i++){
                    cumulative = cumulative + SpacerAcquisitionRate[i]/BacteriaEvolutionRate;
                    if(cumulative>randomNumber1){   // It is when the bacteria evolution happens.
                        double totalOfPhage = 0.0;
                        for(int j=0;j<Dp;j++){
                            totalOfPhage = totalOfPhage + Pold[j];
                        }
                        double randomNumber2 = distribution(generator);
                        double cumulative2 = 0.0;
                        for(int j=0;j<Dp;j++){
                            cumulative2 = cumulative2 + Pold[j]/totalOfPhage;  // randomNumber2 is used to determine from which phage strain the bacterial strain is going to take the protospace.
                            if(cumulative2>randomNumber2){
                                int whichProtospacerToChoose = floor(distribution(generator)*NumOfProtospacers); // whichProtospacerToChoose determines which protospacer from the phage strain to choose as the newly acquired spacer.
                                for(int k=1;k<NumOfSpacers;k++){
                                    Spacers[Db][k] = Spacers[i][k-1];  // <<<---------
                                }
                                Spacers[Db][0] = Protospacers[j][whichProtospacerToChoose];
                                // Check if this evolved bacterium strain already existed:
                                int ifMutatedStrainOfBacteriaAlreadyExisted=0;
                                for(int k=0;k<Db;k++){
                                    for(int l=0;l<NumOfSpacers;l++){
                                        if(Spacers[k][l]==Spacers[Db][l]){
                                            ifMutatedStrainOfBacteriaAlreadyExisted = 1;
                                        }
                                        else{
                                            ifMutatedStrainOfBacteriaAlreadyExisted = 0;
                                            break;
                                        }
                                    }
                                    if(ifMutatedStrainOfBacteriaAlreadyExisted==1){
                                        break;
                                    }
                                }
                                if(ifMutatedStrainOfBacteriaAlreadyExisted==1){
                                    needToBreak=1;
                                    break;
                                }
                                else{

                                    Bold[Db] = 1.1 * densityCutoff;
                                    Blabel[Db] = BlabelLargest+1;
                                    Db = Db + 1;
                                    BlabelLargest = BlabelLargest+1;
                                    needToBreak=1;

                                    break;
                                }
                            }
                        }
                    }
                    if(needToBreak==1){
                        break;
                    }
                }
            }
            // Phages evolve (i.e. mutates):
            else{
                double randomNumber1 = distribution(generator);  // randomNumber1 is used to determine which phage strain will mutate its protospacer.
                double cumulative = 0.0;
                for(int i=0;i<Dp;i++){
                    cumulative = cumulative + ProtospacerMutationRate[i]/PhageEvolutionRate;
                    if(cumulative>randomNumber1){
                        for(int j=0;j<NumOfProtospacers;j++){
                            Protospacers[Dp][j] = Protospacers[i][j];
                        }
//                         cout << timesOfRecord << " " << time << " " << Plabel[i] << " " << endl;  // Intento de Arbol!!!!!
                        int whichProtospacerToMutate = floor(distribution(generator)*NumOfProtospacers);
                        Protospacers[Dp][whichProtospacerToMutate] = LargestIntegerUsed + 1;
                        Pold[Dp] = 1.1 * densityCutoff;
                        Plabel[Dp] = PlabelLargest+1;
                        LargestIntegerUsed = LargestIntegerUsed + 1;
                        Dp = Dp + 1;
                        PlabelLargest = PlabelLargest+1;

                        fileForPhageTree << timesOfRecord << "\t" << PlabelLargest << "\t" << Plabel[i] << "\t" << time << endl;  // IPhage phylogeny tree
                        break;
                    }
                }
            }
        }

        
    }
    
    
    fileForBacteria.close();
    fileForPhage.close();
    fileForTimeSeriesData.close();
    fileForBacteriaAbundance.close();
    fileForBacteriaAnimation.close();
    fileForPhageAbundance.close();
    fileForPhageAnimation.close();
    
    fileForPhageTree.close();
    
    return 0;
    
    
}


// ################################################################################################################################################################################

// CRISPRImmuneFunction is the function checking whether oneSpacerSet[] and oneProtospacerSet[] shares elements:
// CRISPRImmuneFunction compares protospacers with spacers to check if they are matched or not and return 1 for TRUE and 0 as FALSE:
int CRISPRImmuneFunction(int oneSpacerSet[],int oneProtospacerSet[], int NumOfSpacers, int NumOfProtospacers){
    for(int i=0;i<NumOfSpacers;i++){
        for(int j=0;j<NumOfProtospacers;j++){
            if(oneSpacerSet[i]==oneProtospacerSet[j]){
                return 1;
            }
        }
    }
    return 0;
}


// SpacerProtospacerOverlap is the function checking how many elements are shared by oneSpacerSet[] and oneProtospacerSet[]:
int SpacerProtospacerOverlap(int oneSpacerSet[],int oneProtospacerSet[], int NumOfSpacers, int NumOfProtospacers){
    int overlap = 0;
    for(int i=0;i<NumOfSpacers;i++){
        for(int j=0;j<NumOfProtospacers;j++){
            if(oneSpacerSet[i]==oneProtospacerSet[j]){
                overlap = overlap+1;
            }
        }
    }
    return overlap;
}
