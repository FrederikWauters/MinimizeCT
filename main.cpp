/* c++ version of my ct minimizer. 

Combines the minimization, construction of pdf surface and confidence intervals/contours)

dof = 2 -> scan over ct and cs, and minimize Vud and Ca 

dof = 4 -> scan over ct+ctp and ct-ctp 

*/

#include <cstdio>
#include <unistd.h>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "DataManager.h"
#include "TData.h"
#include "Analyzer.h"
#include "TDatime.h"
#include "TTimeStamp.h"

//global variables
char *fout_name;
char *fin_name;
double stepsize, spread_x,spread_y;
int dof;
const int nArgc =13;

using namespace std;

int analyze_command_line (int argc, char **argv);

int main(int argc, char* argv[])
{
  std::cout << "Running tensor minimizer ... " << std::endl;
    
  //Input arguments
  int ret = analyze_command_line(argc, argv);
  if(ret) { return 1; }

  //Set up i/o
  DataManager* io = new DataManager(fin_name,fout_name,dof);
  io->PrintData();
  io->InitHistos(-spread_x,-spread_y,spread_x,spread_y,stepsize);

  //The actual analysis
  Analyzer* ana = new Analyzer(dof,io,"Minuit","Combined");
  ana->TestRun(0.00,0.00,dof);
  ana->Run(-spread_x,-spread_y,spread_y,spread_y,stepsize,dof);

  //Wrap things up
  //io->Plot();
  io->MakeCLContours(1000,10);
  io->Make1DContours(0.9,1000);
  io->GetCorrelation();
  io->WriteOutput();


  return 0;
}

int analyze_command_line(int argc, char **argv)
{
if(argc==nArgc)
  {
    for(int i=1; i<argc; /*increment in loop*/)
    {
      if(argv[i][0] != '-'){ printf("ERROR: Wrong argument %s\n",argv[i]); }
      else
      {
        if(strlen(&argv[i][1]) == 1)
        {
	  switch(argv[i][1])
          {
	    case 'o':
	      if(i+1 < argc)
              {
	        fout_name = argv[i+1];
	        i+=2;
	      }
              else
              {
	        printf("ERROR: No argument for output file specified\n");
	      }
              break;

            case 'i':
	      if(i+1 < argc)
              {
	        fin_name = argv[i+1];
	        i+=2;
	      }
              else
              {
	        printf("ERROR: No argument for input file specified\n");
	      }
              break;

            case 's':
	      if(i+1 < argc)
              {
	        stepsize = atof(argv[i+1]);
	        i+=2;
	      }
              else
              {
	        printf("ERROR: No argument for step size specified\n");
	      }
              break;  

            case 'x':
	      if(i+1 < argc)
              {
	        spread_x = atof(argv[i+1]);
	        i+=2;
	      }
              else
              {
	        printf("ERROR: No argument for plot range specified\n");
	      }
              break;
            
            case 'y':
	      if(i+1 < argc)
              {
	        spread_y = atof(argv[i+1]);
	        i+=2;
	      }
              else
              {
	        printf("ERROR: No argument for plot range specified\n");
	      }
              break;   

            case 'd':
	      if(i+1 < argc)
              {
	        dof = atoi(argv[i+1]);
	        i+=2;
	      }
              else
              {
	        printf("ERROR: No argument for plot range specified\n");
	      }
              break;       

             default:
	     printf("Argument %s not recognized\n",argv[i]);
             return 1;
          }
        }
      }
    }
  }
  else
  {
    cout << " use as ./Minimize -i inputfile -o outputfile -s step -x Xrange -y Yrange -d dof" << endl;
    return 1;
  }
  return 0;
}

