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
  TH1F* compton_mm = (TH1F*) compton->GetObjectChecked("Compton_EmCompTotP;1","TH1F");
  TH1F* GP2H_mm = (TH1F*) GP2H->GetObjectChecked("Compton_EmCompTotP;1","TH1F");
  TH1F* PPN_mm = (TH1F*) PPN->GetObjectChecked("Compton_EmCompTotP;1","TH1F");
  if(e >= 140){
    TH1F* HePi0_mm = (TH1F*) HePi0->GetObjectChecked("Compton_EmCompTotP;1","TH1F");
    TH1F* PPi02H_mm = (TH1F*) PPi02H->GetObjectChecked("Compton_EmCompTotP;1","TH1F");
    TH1F* P2HPi0_mm = (TH1F*) P2HPi0->GetObjectChecked("Compton_EmCompTotP;1","TH1F");
  }

  //Getting the theta histograms
  TH1F* compton_th = (TH1F*) compton->GetObjectChecked("Compton_CThCompP;1","TH1F");
  TH1F* GP2H_th = (TH1F*) GP2H->GetObjectChecked("Compton_CThCompP;1","TH1F");
  TH1F* PPN_th = (TH1F*) PPN->GetObjectChecked("Compton_CThCompP;1","TH1F");
  if(e >= 140){
    TH1F* HePi0_th = (TH1F*) HePi0->GetObjectChecked("Compton_CThCompP;1","TH1F");
    TH1F* PPi02H_th = (TH1F*) PPi02H->GetObjectChecked("Compton_CThCompP;1","TH1F");
    TH1F* P2HPi0_th = (TH1F*) P2HPi0->GetObjectChecked("Compton_CThCompP;1","TH1F");
  }

  //Creating 2D Histograms
  Int_t xbin = compton_mm->FindFixBin(25);
  Int_t ybin = compton_th->FindFixBin(1);
  TH2F *compton_2h = new TH2F("compton_2h","compton",xbin,-25,25,ybin,-1,1);
  TH2F *GP2H_2h = new TH2F("GP2H_2h","GP2H",xbin,-25,25,ybin,-1,1);
  TH2F *PPN_2h = new TH2F("PPN_2h","PPN",xbin,-25,25,ybin,-1,1);
  if( e >= 140){
    TH2F *HePi0_2h = new TH2F("HePi0_2h","G3HePi0",xbin,-25,25,ybin,-1,1);
    TH2F *PPi02H_2h = new TH2F("PPi02H_2h","PPi02H",xbin,-25,25,ybin,-1,1);
    TH2F *P2HPi0_2h = new TH2F("P2HPi0_2h","P2HPi0",xbin,-25,25,ybin,-1,1);
  }

  //Filling the histogram
  for(int i=0; i<xbin; i++){
    for(int j=0; j<ybin; j++){
       Double_t compt_xcont = compton_mm->GetBinContent(i);
       Double_t compt_ycont = compton_th->GetBinContent(j); 
       Double_t GP2H_xcont = GP2H_mm->GetBinContent(i);
       Double_t GP2H_ycont = GP2H_th->GetBinContent(j);
       Double_t PPN_xcont = PPN_mm->GetBinContent(i);
       Double_t PPN_ycont = PPN_th->GetBinContent(j);

       if(e >=140){
	 Double_t HePi0_xcont = HePi0_mm->GetBinContent(i);
	 Double_t HePi0_ycont = HePi0_th->GetBinContent(j);
	 Double_t PPi02H_xcont = PPi02H_mm->GetBinContent(i);
	 Double_t PPi02H_ycont = PPi02H_th->GetBinContent(j);
	 Double_t P2HPi0_xcont = P2HPi0_mm->GetBinContent(i);
	 Double_t P2HPi0_ycont = P2HPi0_th->GetBinContent(j);
       }

       Double_t compt_cont = TMath::Sqrt((compt_xcont+compt_ycont));
       Double_t GP2H_cont = TMath::Sqrt((GP2H_xcont+GP2H_ycont));
       Double_t PPN_cont = TMath::Sqrt((PPN_xcont+PPN_ycont));
       if(e >= 140){
	 Double_t HePi0_cont = (HePi0_xcont+HePi0_ycont);
	 Double_t PPi02H_cont = (PPi02H_xcont+PPi02H_ycont);
	 Double_t P2HPi0_cont = (P2HPi0_xcont+P2HPi0_ycont);
       }

       compton_2h->SetBinContent(i,j,compt_cont);
       GP2H_2h->SetBinContent(i,j,GP2H_cont);
       PPN_2h->SetBinContent(i,j,PPN_cont);
       if(e >= 140){
	 HePi0_2h->SetBinContent(i,j,HePi0_cont);
	 PPi02H_2h->SetBinContent(i,j,PPi02H_cont);
	 P2HPi0_2h->SetBinContent(i,j,P2HPi0_cont);
       }
    }
  }

  //Scale histograms
  if(e == 100){
    GP2H_2h->Scale(2.13);
    PPN_2h->Scale(200);
  }
  else if(e == 120){
    GP2H_2h->Scale(1.60);
    PPN_2h->Scale(100);
  }
  else if(e == 200){
    HePi0_2h->Scale(58.18);
    GP2H_2h->Scale(0.79);
    PPi02H_2h->Scale(40.38);
    P2HPi0_2h->Scale(192.3);
    PPN_2h->Scale(9.62);
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
