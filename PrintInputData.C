//Open output file and execute

//gSystem->Load("/home/frederik/Analysis/CTMin/New/MinimizeCT/libDataClasses.so");

{  
   
   //tree_file0
   int nentries = InputData->GetEntries();
   TBranch* br = InputData->GetBranch("datapoint");
   cout << nentries << " datapoints " << endl;
   TData *data = new TData(); 
   br->SetAddress(&data);
   
   cout << "isotope" << "	" << "par" << "	" << "jInit"  << "	"<< "jFinal" << "	" << "zDaught" << "	" << "b +/-" << "	"<< "mF" << "	" << "mGT" << "	"<< "mE" << "	"<< "expVal" << "	" << "error" << "	" << "use" << endl<< endl;
   
  //for(int i = 0; i<data.size(); i++) data.at(i)->Print
   for(int i = 1; i<=nentries; i++)
   {
     InputData->GetEntry(i);
     data->Print();
   }

}
