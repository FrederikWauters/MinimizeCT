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
   for(Int_t i = 0; i<nentries; i++)
   {
     InputData->GetEntry(i);
     data->Print();
   }
   
   cout <<" Start loop in CT C'T space for 6He" << endl;
   
   InputData->GetEntry(0);
   
   double cs = 0;
   double cpt = 0;
   double ct = 0.;
   double cps = 0.;
   double ca = -1.2723;
   double vud = 0.97425;
   double cv = 1;
   
   double value;
   
   value = data->a(cv, ca, cs, cps, ct, cpt);
   
   cout << " value = " << value << endl;
   
   int N = 1200;
   double min = -0.6;
   double max = 0.6;
   double step = (max-min)/(1.*N); 
   
   TH2F* h = new TH2F("h","h",N,min,max,N,min,max);
   
//   ca = ca*0.9;
//   value = data->a(cv, ca, cs, cps, ct, cpt);
//   cout << " value = " << value << endl;
//   

//   ca = ca*0.9;
//   value = data->a(cv, ca, cs, cps, ct, cpt);
//   cout << " value = " << value << endl;
//   ca = ca*0.9;
//   value = data->a(cv, ca, cs, cps, ct, cpt);
//   cout << " value = " << value << endl;
//   ca = ca*0.9;
//   value = data->a(cv, ca, cs, cps, ct, cpt);
//   cout << " value = " << value << endl;
//   ca = ca*0.9;
//   value = data->a(cv, ca, cs, cps, ct, cpt);
//   cout << " value = " << value << endl;
//   ca = ca*0.9;
//   value = data->a(cv, ca, cs, cps, ct, cpt);
//   cout << " value = " << value << endl;
//   
   
   
   for(int i = 0; i <= N; i++)
   {
     for(int j = 0; j <= N;j++)
     {
       double x = min+i*step;
       double y = min+j*step;
       ct = ca*(x + y)/2.;
       cpt = ca*(x - y)/2.;
       value = data->a(cv, ca, cs, cps, ct, cpt);
       h->SetBinContent(i,j,value);       
     }
   }

}
