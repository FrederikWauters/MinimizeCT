//Open output file and execute

{  
   gSystem->Load("libDataClasses.so");
   tree_file0
   int nentries = InputData->GetEntries();
   TBranch* br = InputData->GetBranch("datapoint");
   
   TData *data = new TData(); 
   /*br->SetAddress(&data);
   
   for(int i = 1; i<=entries; i++)
   {
     InputData->GetEntry(i);
   }*/

}
