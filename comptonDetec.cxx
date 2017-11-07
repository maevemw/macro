int comptonDetec(Int_t e = 100, Int_t events = 1000000){

  if(e%20 != 0){
    return 2;
  }

  TString comp_name; comp_name.Form("data/G3He/%d_%d.root",e,events);
  TString total_name; total_name.Form("dataMC/G3He/hth%d_%d.root",e,events);


  //Gathering the root files
  TFile* compton = new TFile(comp_name);
  TFile* total = new TFile(total_name);


  //Getting the missing mass histograms
  TH2F* compton_2h = (TH2F*) compton->GetObjectChecked("Compton_CDegCompP_v_EmCompTotP;1","TH2F");
  TH1D* total_1h = (TH1D*) total->GetObjectChecked("hTheta;1","TH1D");

  compton_2h->RebinY();

  //Finding number of events in peak range
  Double_t xhigh = compton_2h->GetXaxis()->FindFixBin(1);
  Double_t xlow = compton_2h->GetXaxis()->FindFixBin(-3);
  Double_t thetatotal = 0.0;
  Double_t yhigh = compton_2h->GetYaxis()->FindFixBin(169);
  Double_t ylow = compton_2h->GetYaxis()->FindFixBin(10);
  Double_t yhighish = total_1h->FindFixBin(169);
  Double_t ylowish = total_1h->FindFixBin(10);
  Double_t totaltrue = total_1h->Integral(ylowish, yhighish);
  Double_t comptotaltrue = compton_2h->Integral(xlow, xhigh, ylow, yhigh);
  Double_t comptotal = 0.0;
    
  for( Int_t x=20; x < 180; x = x+20){

    Double_t yhighc = compton_2h->GetYaxis()->FindFixBin(x+9);
    Double_t ylowc = compton_2h->GetYaxis()->FindFixBin(x-10);
    Double_t yhight = total_1h->FindFixBin(x+9);
    Double_t ylowt = total_1h->FindFixBin(x-10);

    Double_t total_peak = total_1h->Integral(ylowt,yhight);
    Double_t comp_peak = compton_2h->Integral(xlow,xhigh,ylowc,yhighc);
    comptotal = comptotal + comp_peak;
    thetatotal = thetatotal + total_peak;

    cout<<"At "<<x<<" degrees:"<<endl;
    cout<<"Number of events in total: "<<total_peak<<endl;
    cout<<"Number of events in Compton peak: "<<comp_peak<<endl;
    cout<<"Detection Efficiency : "<<comp_peak/total_peak<<endl;
  }

  cout<<"To check for double counting this should be "<<totaltrue <<" : "<< thetatotal << endl;
  cout<<"this should be " << comptotaltrue<< ": "<<comptotal<<endl;
  return 0;
}
