int mmass_theta(Int_t E = 100, Int_t events = 100000){
  
  if(E%20 != 0){
    return 2;
  }
  
//Setting up the Palette
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

//Setting up paths to files
  
  TString comp_name; comp_name.Form("data/G3He/%d_%d.root",E,events);
  TString GP2H_name; GP2H_name.Form("data/GP2H_QF/%d_%d.root",E,events);
  TString PPN_name; PPN_name.Form("data/PPN/%d_%d.root",E,events);
  if(E >= 140){
    TString HePi0_name; HePi0_name.Form("data/3HePi0/%d_%d.root",E,events);
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

  gStyle->SetOptStat(0);
  
//Format and Draw Histograms
  TString title; title.Form("E=%dMeV: Missing Energy",E);
  TCanvas* c = new TCanvas(title,title);
  if(E < 140){
    c->Divide(3);
  }
  else{
    c->Divide(3,2);
  }

  c->cd(1);
  compton_2h->SetTitle("Compton");
  compton_2h->SetXTitle("Missing Energy (MeV)");
  compton_2h->SetYTitle("Cos(#theta)");
  compton_2h->SetTitleOffset(0.8, "XY");
  compton_2h->SetTitleSize(0.05, "XY");
  compton_2h->SetAxisRange(-10,10,"X");
  compton_2h->DrawClone("col");

  c->cd(2);
  GP2H_2h->SetTitle("QF from Proton");
  GP2H_2h->SetXTitle("Missing Energy (MeV)");
  GP2H_2h->SetYTitle("Cos(#theta)");
  GP2H_2h->SetTitleOffset(0.8, "XY");
  GP2H_2h->SetTitleSize(0.05, "XY");
  GP2H_2h->SetAxisRange(-10,10,"X");
  GP2H_2h->DrawClone("col");
  
  c->cd(3);
  PPN_2h->SetTitle("Pure Breakup");
  PPN_2h->SetXTitle("Missing Energy (MeV)");
  PPN_2h->SetYTitle("Cos(#theta)");
  PPN_2h->SetTitleOffset(0.8, "XY");
  PPN_2h->SetTitleSize(0.05, "XY");
  PPN_2h->SetAxisRange(-10,10,"X");
  PPN_2h->DrawClone("col");

  if(E >=140){
    c->cd(4);
    HePi0_2h->SetTitle("Pi0 production");
    HePi0_2h->SetXTitle("Missing Energy (MeV)");
    HePi0_2h->SetYTitle("Cos(#theta)");
    HePi0_2h->SetTitleOffset(0.8, "XY");
    HePi0_2h->SetTitleSize(0.05, "XY");
    HePi0_2h->SetAxisRange(-10,10,"X");
    HePi0_2h->DrawClone("col");

    c->cd(5);
    PPi02H_2h->SetTitle("QF Pi0 from proton");
    PPi02H_2h->SetXTitle("Missing Energy (MeV)");
    PPi02H_2h->SetYTitle("Cos(#theta)");
    PPi02H_2h->SetTitleOffset(0.8, "XY");
    PPi02H_2h->SetTitleSize(0.05, "XY");
    PPi02H_2h->SetAxisRange(-10,10,"X");
    PPi02H_2h->DrawClone("col");

    c->cd(6);
    P2HPi0_2h->SetTitle("QF Pi0 from deuteron");
    P2HPi0_2h->SetXTitle("Missing Energy (MeV)");
    P2HPi0_2h->SetYTitle("Cos(#theta)");
    P2HPi0_2h->SetTitleOffset(0.8, "XY");
    P2HPi0_2h->SetTitleSize(0.05, "XY");
    P2HPi0_2h->SetAxisRange(-10,10,"X");
    P2HPi0_2h->DrawClone("col");
  }
  	
//Add Histograms
  TH2F* total_2h=new TH2F(*compton_2h);
  total_2h->Add(PPN_2h);
  total_2h->Add(GP2H_2h);
  if(E >= 140){
    total_2h->Add(HePi0_2h);
    total_2h->Add(PPi02H_2h);
    total_2h->Add(P2HPi0_2h);
  }
  
//Integrate over 20 degree intervals
  Double_t xlow, xhigh;
  xlow = compton_2h->GetXaxis()->FindFixBin(-3);
  xhigh = total_2h->GetXaxis()->FindFixBin(1);
  for(Int_t x=0; x < 180; x = x+20){

    Double_t theta_min = ((Double_t)x * TMath::Pi())/180;
    Double_t theta_max = ((Double_t)(x+20) * TMath::Pi())/180;
    Double_t yhigh = compton_2h->GetYaxis()->FindFixBin(TMath::Cos(theta_min));
    Double_t ylow = compton_2h->GetYaxis()->FindFixBin(TMath::Cos(theta_max)+0.02);
  
    Double_t total_peak = total_2h->Integral(xlow,xhigh,ylow,yhigh);
    Double_t comp_peak = compton_2h->Integral(xlow,xhigh,ylow,yhigh);
    Double_t PPN_peak = PPN_2h->Integral(xlow,xhigh,ylow,yhigh);
    Double_t GP2H_peak = GP2H_2h->Integral(xlow,xhigh,ylow,yhigh);
    if(E >=140){
      Double_t HePi0_peak = HePi0_2h->Integral(xlow,xhigh,ylow,yhigh);
      Double_t PPi02H_peak = PPi02H_2h->Integral(xlow,xhigh,ylow,yhigh);
      Double_t P2HPi0_peak = P2HPi0_2h->Integral(xlow,xhigh,ylow,yhigh);
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
    cout<<"Contaminaiton at "<<x+20<<" degrees scattering angle"<<endl;
    cout<<"Compton: "<<comp_cont<<"%"<<endl;
    cout<<"GP2H_QF: "<<GP2H_cont<<"%"<<endl;
    cout<<"PPN: "<<PPN_cont<<"%"<<endl;
    if(E >=140){
      cout<<"3HePi0: "<<HePi0_cont<<"%"<<endl;
      cout<<"PPi02H_QF: "<<PPi02H_cont<<"%"<<endl;
      cout<<"P2HPi0_QF: "<<P2HPi0_cont<<"%"<<endl;
    }
  }
}

