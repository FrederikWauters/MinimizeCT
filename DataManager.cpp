/**

Process input

**/

//Std includes
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cmath>


#include "TData.h"
#include "DataManager.h"

using std::cout;
using std::endl;


//constructor
DataManager::DataManager(char* fin, char* fout)
{
  cout << "Construct data manager ... " << endl;

  data.clear();
  
  std::ifstream *inputfile = new std::ifstream(fin);
  if(inputfile->is_open()) cout << fin << " open ... " << endl;
    else cout << "failed to open inputfile ... " << endl; 
  if(!ProcessInputFile(inputfile)) cout << "input not properly processed" << endl; 
  inputfile->close();
 
  if(SetOutputFile(fout)==1) cout << fout << " open ... " << endl;
}

//Destructor
DataManager::~DataManager()
{
  cout << "Delete data manager ... " << endl;
}

int DataManager::ProcessInputFile(std::ifstream* file)
{
  if(file==NULL) return 0;
  std::string line;
  int nInput = 0;
  
  while(std::getline(*file,line))
  {
    if(line[0]=='#' || line.size()<2) continue;
    
    else  
    {
      std::stringstream ss;
      ss << line;
      TData* newData = new TData();
      newData->SetData(ss);
      if(newData->Use()) data.push_back(newData);
  
      nInput++;
    }
  }

  cout<< "size " << data.size() << endl;
  return 1;
}

int DataManager::SetOutputFile(char* outfile)
{
  fout = new TFile ( outfile, "ReCreate");

  if ( fout->IsZombie() ) {
    std::cerr << "**ERROR! Cannot open file [" << outfile << "]" << endl;
    return 0;
  }
  fout->cd();
  return 1;
}

void DataManager::PrintData()
{
  cout << endl;
  cout << "***************************************" << endl;
  cout << "************* Input Data **************" << endl;
  cout << "***************************************" << endl << endl;
  if(data.size()<1) return;
cout << "isotope" << "	" << "par" << "	" << "jInit"  << "	"<< "jFinal" << "	" << "zDaught" << "	" << "sign" << "	"<< "mF" << "	" << "mGT" << "	"<< "mE" << "	"<< "expVal" << "	" << "error" << "	" << "use" << endl;
  for(int i = 0; i<data.size(); i++) data.at(i)->Print();

  cout << endl;
   
}
