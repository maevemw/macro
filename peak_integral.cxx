int peak_integral(Int_t e = 100, Int_t events = 100000){

  if(e%20 != 0){
    return 2;
  }

  //Setting path names
  TString comp_name; comp_name.Form("data/G3He/%d_%d.root",e,events);
  TString GP2H_name; GP2H_name.Form("data/GP2H_QF/%d_%d.root",e,events);
  TString PPN_name; PPN_name.Form("data/PPN/%d_%d.root",e,events);
  if(e >=140){
    TString HePi0_name; HePi0_name.Form("data/3HePi0/%d_%d.root",e,events);
    TString PPi02H_name; PPi02H_name.Form("data/PPi02H_QF/%d_%d.root",e,events);
    TString P2HPi0_name; P2HPi0_name.Form("data/P2HPi0_QF/%d_%d.root",e,events);
  }
    /*
  TString comp_name; comp_name.Form("LowLim/G3He/%d_%d.root",e,events);
  TString GP2H_name; GP2H_name.Form("LowLim/GP2H_QF/%d_%d.root",e,events);
  TString PPN_name; PPN_name.Form("LowLim/PPN/%d_%d.root",e,events);
  if(e >=140){
    TString HePi0_name; HePi0_name.Form("LowLim/3HePi0/%d_%d.root",e,events);
    TString PPi02H_name; PPi02H_name.Form("LowLim/PPi02H_QF/%d_%d.root",e,events);
    TString P2HPi0_name; P2HPi0_name.Form("LowLim/P2HPi0_QF/%d_%d.root",e,events);
    }*/

  //Gathering the root files
  TFile* compton = new TFile(comp_name);
  TFile* GP2H = new TFile(GP2H_name);
  TFile* PPN = new TFile(PPN_name);
  if(e >= 140){
    TFile* HePi0 = new TFile(HePi0_name);
    TFile* PPi02H = new TFile(PPi02H_name);
    TFile* P2HPi0 = new TFile(P2HPi0_name);
  }

  //Getting the missing mass histograms
  TH2F* compton_2h = (TH2F*) compton->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");
  TH2F* GP2H_2h = (TH2F*) GP2H->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");
  TH2F* PPN_2h = (TH2F*) PPN->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");
  if(e >= 140){
    TH2F* HePi0_2h = (TH2F*) HePi0->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");
    TH2F* PPi02H_2h = (TH2F*) PPi02H->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");
    TH2F* P2HPi0_2h = (TH2F*) P2HPi0->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");
  }

  //Add Histograms
  TH2F* total_2h=new TH2F(*compton_2h);
  total_2h->Add(PPN_2h);
  total_2h->Add(GP2H_2h);
  if(e >= 140){
    total_2h->Add(HePi0_2h);
    total_2h->Add(PPi02H_2h);
    total_2h->Add(P2HPi0_2h);
  }

  //Finding number of events in peak range
  Double_t xhigh = compton_2h->GetXaxis()->FindFixBin(1);
  Double_t xlow = compton_2h->GetXaxis()->FindFixBin(-3);
  
  for( Int_t x=0; x < 180; x = x+20){

    Double_t theta_min = ((Double_t)x * TMath::Pi())/180;
    Double_t theta_max = ((Double_t)(x+20) * TMath::Pi())/180;
    Double_t yhigh = compton_2h->GetYaxis()->FindFixBin(TMath::Cos(theta_min));
    Double_t ylow = compton_2h->GetYaxis()->FindFixBin(TMath::Cos(theta_max));

    Double_t total_peak = total_2h->Integral(xlow,xhigh,ylow,yhigh);
    Double_t comp_peak = compton_2h->Integral(xlow,xhigh,ylow,yhigh);
    Double_t GP2H_peak = GP2H_2h->Integral(xlow,xhigh,ylow,yhigh);
    Double_t PPN_peak = PPN_2h->Integral(xlow,xhigh,ylow,yhigh);
    if(e >= 140){
      Double_t HePi0_peak = HePi0_2h->Integral(xlow,xhigh,ylow,yhigh);
      Double_t PPi02H_peak = PPi02H_2h->Integral(xlow,xhigh,ylow,yhigh);
      Double_t P2HPi0_peak = P2HPi0_2h->Integral(xlow,xhigh,ylow,yhigh);
    }


    cout<<"At "<<x<<" degrees:"<<endl;
    cout<<"Number of events in total: "<<total_peak<<endl;
    cout<<"Number of events in Compton peak: "<<comp_peak<<endl;
    cout<<"Number of events in GP2H peak: "<<GP2H_peak<<endl;
    cout<<"Number of events in PPN peak: "<<PPN_peak<<endl;
    if(e >= 140){
      cout<<"Number of events in HePi0 peak: "<<HePi0_peak<<endl;
      cout<<"Number of events in PPi02H peak: "<<PPi02H_peak<<endl;
      cout<<"Number of events in P2HPi0 peak: "<<P2HPi0_peak<<endl;
    }
  }
  return 0;
}
