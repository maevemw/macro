int mmass_theta_1D(Int_t E = 100, Int_t theta = 80, Int_t events = 100000){
  
  if(E%20 != 0){
    return 2;
  }
  if(theta%20 != 0){
    return 3;
  }

//Setting up paths to files
  TString comp_name; comp_name.Form("data/G3He/%d_%d.root",E,events);
  TString GP2H_name; GP2H_name.Form("data/GP2H_QF/%d_%d.root",E,events);
  TString PPN_name; PPN_name.Form("data/PPN/%d_%d.root",E,events);
  if(E >= 140){
    TString HePi0_name; HePi0_name.Form("data/G3HePi0/%d_%d.root",E,events);
    TString PPi02H_name; PPi02H_name.Form("data/PPi02H_QF/%d_%d.root",E,events);
    TString P2HPi0_name; P2HPi0_name.Form("data/P2HPi0_QF/%d_%d.root",E,events);
  }
  
//Import data from files
  TFile* compton = new TFile(comp_name);
  TFile* GP2H = new TFile(GP2H_name);
  TFile* PPN = new TFile(PPN_name);
  if(E >= 140){
    TFile* HePi0 = new TFile(HePi0_name);
    TFile* PPi02H = new TFile(PPi02H_name);
    TFile* P2HPi0 = new TFile(P2HPi0_name);
  }

//Convert data to a histogram
  TH2F *compton_2h = (TH2F*) compton->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");
  
  TH2F *GP2H_2h = (TH2F*) GP2H->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");

  TH2F *PPN_2h = (TH2F*) PPN->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");

  if(E >=140){
    TH2F* HePi0_2h = (TH2F*) HePi0->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");

    TH2F* PPi02H_2h = (TH2F*) PPi02H->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");

    TH2F* P2HPi0_2h = (TH2F*) P2HPi0->GetObjectChecked("Compton_CThCompP_v_EmCompTotP;1","TH2F");
  }

//Scale histograms
  if(E == 100){
    GP2H_2h->Scale(2.13);
    PPN_2h->Scale(200);
  }
  else if(E == 120){
    GP2H_2h->Scale(1.60);
    PPN_2h->Scale(100);
  }
  else if(E == 200){
    HePi0_2h->Scale(58.18);
    GP2H_2h->Scale(0.79);
    PPi02H_2h->Scale(40.38);
    P2HPi0_2h->Scale(192.3);
    PPN_2h->Scale(9.62);
  }

//Project the theta range onto the X axis
    Double_t theta_min = ((Double_t)theta * TMath::Pi())/180;
    Double_t theta_max = ((Double_t)(theta+20) * TMath::Pi())/180;
    Double_t yhigh = compton_2h->GetYaxis()->FindFixBin(TMath::Cos(theta_min));
    Double_t ylow = compton_2h->GetYaxis()->FindFixBin(TMath::Cos(theta_max));

    TH1D* compton_h = compton_2h->ProjectionX("compton_1h", ylow, yhigh);
    TH1D* GP2H_h = GP2H_2h->ProjectionX("GP2H_h", ylow, yhigh);
    TH1D* PPN_h = PPN_2h->ProjectionX("PPN_h", ylow, yhigh);
    if(E >= 140){
      TH1D* HePi0_h = HePi0_2h->ProjectionX("HePi0_h", ylow, yhigh);
      TH1D* PPi02H_h = PPi02H_2h->ProjectionX("PPi02H_h", ylow, yhigh);
      TH1D* P2HPi0_h = P2HPi0_2h->ProjectionX("P2HPi0_h", ylow, yhigh);
    }
      
//Format and Draw Histograms
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas();

  TString title; title.Form("E=%dMeV: Missing Energy",E);
  compton_h->SetTitle(title);
  compton_h->SetXTitle("Missing Energy");
  compton_h->SetYTitle("Number of Events");
  compton_h->SetAxisRange(-10,10,"X");
  compton_h->SetLineColor(1);
  compton_h->DrawClone("SameHist");

  GP2H_h->SetLineColor(3);
  GP2H_h->DrawClone("SameHist");
  
  PPN_h->SetLineColor(7);
  PPN_h->DrawClone("Samehist");

  if(E >=140){
    HePi0_h->SetLineColor(2);
    HePi0_h->DrawClone("SameHist");

    PPi02H_h->SetLineColor(4);
    PPi02H_h->DrawClone("SameHist");

    P2HPi0_h->SetLineColor(6);
    P2HPi0_h->DrawClone("SameHist");
  }

  //Add Histograms
  TH1D *total_h=new TH1D(*compton_h);
  total_h->Add(PPN_h);
  total_h->Add(GP2H_h);
  if(E >= 140){
    total_h->Add(HePi0_h);
    total_h->Add(PPi02H_h);
    total_h->Add(P2HPi0_h);
  }
  
//Calculate integral
  Double_t lowbin, highbin;
  lowbin = total_h->FindFixBin(-3);
  highbin = total_h->FindFixBin(1);
  Double_t total_peak = total_h->Integral(lowbin,highbin);
  Double_t comp_peak = compton_h->Integral(lowbin,highbin);
  Double_t PPN_peak = PPN_h->Integral(lowbin,highbin);
  Double_t GP2H_peak = GP2H_h->Integral(lowbin,highbin);
  if(E >=140){
    Double_t HePi0_peak = HePi0_h->Integral(lowbin,highbin);
    Double_t PPi02H_peak = PPi02H_h->Integral(lowbin,highbin);
    Double_t P2HPi0_peak = P2HPi0_h->Integral(lowbin,highbin);
  }
  
//Determine contamination
  Double_t comp_cont = comp_peak/total_peak * 100;
  Double_t GP2H_cont = GP2H_peak/total_peak * 100;
  Double_t PPN_cont = PPN_peak/total_peak * 100;
  if(E >=140){
    Double_t HePi0_cont = HePi0_peak/total_peak * 100;
    Double_t PPi02H_cont = PPi02H_peak/total_peak * 100;
    Double_t P2HPi0_cont = P2HPi0_peak/total_peak * 100;
  }

//Print out contaminations
  cout<<"Compton: "<<comp_cont<<"%"<<endl;
  cout<<"GP2H_QF: "<<GP2H_cont<<"%"<<endl;
  cout<<"PPN: "<<PPN_cont<<"%"<<endl;
  if(E >=140){
    cout<<"3HePi0: "<<HePi0_cont<<"%"<<endl;
    cout<<"PPi02H_QF: "<<PPi02H_cont<<"%"<<endl;
    cout<<"P2HPi0_QF: "<<P2HPi0_cont<<"%"<<endl;
  }

  return 0;
}
