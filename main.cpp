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
int stepsize;
const int nArgc = 7;

using namespace std;

int analyze_command_line (int argc, char **argv);

int main(int argc, char* argv[])
{
  std::cout << "Running Delayed Event Tree analysis ... " << std::endl;
    
  //Input arguments
  int ret = analyze_command_line(argc, argv);
  if(ret) { return 1; }

  //Set up i/o
  DataManager* io = new DataManager(fin_name,fout_name);
  io->PrintData();

  //The actual analysis
  //TTimeStamp *time = new TTimeStamp;
  int dof = 2;
  Analyzer* ana = new Analyzer(2,io,"Minuit","Combined");
  
  



  //Wrap things up

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
	        printf("ERROR: No argument for input file specified\n");
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
	        stepsize = atoi(argv[i+1]);
	        i+=2;
	      }
              else
              {
	        printf("ERROR: No argument for input file specified\n");
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
    cout << " use as ./Minimize -i inputfile -o outputfile -s step" << endl;
    return 1;
  }
  return 0;
}

