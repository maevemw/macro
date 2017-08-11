gROOT->Reset();

#include "/home/maeve/Code/lib/physics.h"
#include "/home/maeve/Code/lib/functions.h"

void DiffXSecError(Double_t e = 120)
{
        // Tagging Efficiency
        Double_t tageff = 0.5;

	// # tagger channels binned together
	Double_t ebin = 10;

	// Flux per channel (Hz/MeV)
	Double_t tagflux = 1000000;

	// Photon Flux (Hz/MeV)
	Double_t photonflux = tageff*tagflux*ebin;

	// Data Acqusition Live Time
	Double_t livet = 0.7;

	// Target Thickness (nuclei/nb)
	Double_t targt = 1.128*pow(10, -11);

	// Beam Time (seconds)
	Double_t beamt = 200*3600;

//Alpha values//////////////////////////////////////////////////////////////////
	
	// Solid angle, Yield , Yield error
	Double_t adomega, ayield, aerror;

	// Photon energy, Recoil angle, Differential cross section -2,
	// Differential cross section +2, Varied Standard Unit *2 ,Chi squared
	Double_t aegamma, atheta, adiffxssm, adiffxsbg, adiffxs, adxs, achi_sq, achi_sq_sm, achi_sq_bg;

	// Line values
	// Theta, Theta error,
	// Differential cross section (-2),
	// Differential cross section (-2) error,
	// Differential cross section (+2),
	// Differential cross section (+2) error,
	// Differential cross section,
	// Differential cross section error
	Double_t athet1[20], adthet1[20], adixssm1[20], addixssm1[20], adixsbg1[20], addixsbg1[20], adixs1[20], addixs1[20];

	// Point values
	// Theta, Theta error,
	// Differential cross section (-2),
	// Differential cross section (-2) error,
	// Differential cross section (+2),
	// Differential cross section (+2) error,
	// Differential cross section,
	// Differential cross section error	
	Double_t athet2[20], adthet2[20], adixssm2[20], addixssm2[20], adixsbg2[20], addixsbg2[20], adixs2[20], addixs2[20];

	// Smeared DXS values alpha
	Double_t adixsran[20];

	Int_t i, j;

	// Data file to be read in
	TString nalpha = "3HeXsec_neutron_alphavaried.dat";

	// Sets Detection efficiency based on Energy
	Double_t deteff = 1;
	Int_t lowDXS = 1, highDXS = 100;
	Int_t lowX = 5, highX = 100;
	if (e == 100) 
	{
	       deteff = 0.326; //My value
	       lowDXS = 15;
	       highDXS = 50;
	}	
	else if (e == 120) 
	{
	        deteff = 0.818; //Meg's Value
		lowDXS = 10;
		highDXS = 25;
	}

	// Ensures file nalpha exists and can be opened
	ifstream alph(nalpha);
	if ( !alph.is_open()) 
	{
		cout << "Error opening file ";
		cout << nalpha;
		cout << endl;
		break;
	}

	TRandom3 r(0);

	// Set all DXS values to be within range of DXS high, low
	TF1 *f1 = new TF1("f1", "gaus", lowDXS, highDXS);

	// Random amplitude
	f1->SetParameter(0,10);

	// Zero all variables to prevent errors
	i = 0;
	achi_sq = 0;
	achi_sq_sm = 0;
	achi_sq_bg = 0;
	j = 0;

	// Read through data until the end of file
	while(!alph.eof())
	{
	        // Read each line into appropriate variables
		alph >> aegamma >> atheta >> adiffxssm >> adiffxs >> adiffxsbg >> adxs;

		// Examine lines pertaining to the set energy
		if(aegamma == e)
		{
		        // Error for all the following values are 0
		        // Set x data = theta		        
			athet1[i] = atheta;    
			adthet1[i] = 0;

			// Set y data line 3 = DXS -2
			adixssm1[i] = adiffxssm;	
			addixssm1[i] = 0;

			// Set y data line 2 = DXS +2
			adixsbg1[i] = adiffxsbg;	
			addixsbg1[i] = 0;	    	

			// Set y data line 4 = DXS
			adixs1[i] = adiffxs;   
			addixs1[i] = 0;	  

			
			j = i/2;

			// For angles 20-160 incremented by 20
			if ( i%2 == 0 && j!=0 && j!=9)
			{
			        // Set x data = theta
			        // Theta error = 0
				athet2[j] = atheta;	
				adthet2[j] = 0;	       
				
				// Set y data line 3 = DXS -2
				adixssm2[j] = adiffxssm;	
				
				// Set y data line 2 = DXS +2
				adixsbg2[j] = adiffxsbg;	
				
				// Set y data line 4 = DXS
				adixs2[j] = adiffxs;       
				
				// Solid angle
				adomega = 2*kPI*(cos((athet2[j]-10)*kPI/180)-cos((athet2[j]+10)*kPI/180));
				
				// Yield and Error
				ayield = adixs2[j]*photonflux*deteff*livet*targt*adomega*beamt;
				aerror = sqrt(ayield);

				// DXS uncertainty
				addixs2[j] = adixs2[j]/aerror;
			      	addixssm2[j] = adixssm2[j]/aerror;
				addixsbg2[j] = adixsbg2[j]/aerror;
				
				// Smear DXS values
				// Mean = central DXS
				f1->SetParameter(1,adixs2[j]);
				
				// St. Dev. = DXS error
				f1->SetParameter(2,addixs2[j]); 

				// Randomize point location
				adixsran[j] = f1->GetRandom();

				// Print theta, DXS, Randomize point location
				cout << athet2[j] << " " << adixsran[j] << " " << addixs2[j]<< endl;

				// Calculates the chi-sq for this point
				achi_sq += pow((adixsran[j] - adixs2[j])/addixs2[j], 2);
			       	achi_sq_sm += pow((adixsran[j] - adixssm2[j])/addixssm2[j], 2);
				achi_sq_bg += pow((adixsran[j] - adixsbg2[j])/addixsbg2[j], 2);

			}
			i++;	
		}
	}
	alph.close();

	// Calculate chi squared
	if (j > 1){
	        achi_sq /= j-2;
		achi_sq_sm /= j-2;
		achi_sq_bg /= j-2;
	}

	// For unacceptable values set Chi-squared to zero to prevent error
	else if ( (j == 1) || ( j == 0)){
	        achi_sq = 0;
		achi_sq_sm = 0;
		achi_sq_bg = 0;
	}

	// Print Chi-squared for alpha
	cout << "Alpha Reduced Chi Squared = " << achi_sq << endl;
	//cout << "Alpha Reduced Chi Squared -2 = " << achi_sq_sm << endl;
	//cout << "Alpha Reduced Chi Squared +2 = " << achi_sq_bg << endl;

	// Plot trend line for alpha unvaried
	TGraphErrors *gr1 = new TGraphErrors( i, athet1, adixs1, adthet1, addixs1);
	TString titlea = Form("Alpha Varied Differential Cross Sections for E_{#gamma}=%3.0f MeV", e);
	gr1->SetTitle(titlea);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gr1->SetLineWidth(3);
	gr1->SetLineStyle(1);
	gr1->GetXaxis()->SetTitleOffset(1.1);
	gr1->GetYaxis()->SetTitleOffset(1.0);
	gr1->GetXaxis()->SetTitle("#theta (deg)");
	gr1->GetYaxis()->SetTitle("d#sigma/d#Omega (nb/sr)");
	gr1->GetXaxis()->SetLabelSize(0.03);
	gr1->GetYaxis()->SetLabelSize(0.03);
	gr1->GetXaxis()->CenterTitle();
	gr1->GetYaxis()->CenterTitle();
	gr1->GetXaxis()->SetRangeUser(0,180);
	gr1->GetYaxis()->SetRangeUser( lowDXS, highDXS);

	// Plot points with error bars for alpha unvaried, every 20 deg ** smeared **
	TGraphErrors *gr2 = new TGraphErrors( i , athet2, adixsran, adthet2, addixs2);
	gr2->SetTitle();
	gr2->SetMarkerColor(4);
	gr2->SetMarkerStyle(20);
	gr2->SetMarkerSize(1.1);
	gr2->SetLineWidth(2);

	// Plot trend line for alpha varied by +2
	TGraphErrors *gr3 = new TGraphErrors( i, athet1, adixsbg1, adthet1, addixsbg1);
	gr3->SetTitle();
	gr3->SetLineWidth(3);
	gr3->SetLineStyle(9);

	// Plot trend line for alpha varied by -2
	TGraphErrors *gr4 = new TGraphErrors( i, athet1, adixssm1, adthet1, addixssm1);
	gr4->SetTitle();
	gr4->SetLineWidth(3);
	gr4->SetLineStyle(2);

	TCanvas *c1 = new TCanvas("c1", "Alpha and Beta Varied Differential Cross Sections", 0, 0, 900, 900);	
	c1->Divide(1,2);
	c1->cd(1);

	gr1->Draw("AL");
	gr2->Draw("same" "P");
	gr3->Draw("same");	
	gr4->Draw("same");

	pt = new TLegend(0.6, 0.7, 0.85, 0.85);
	pt->SetFillColor(0);
	pt->SetBorderSize(0);
	pt->SetTextSize(0.04);

	pt->AddEntry(gr2, "Alpha", "lp");
	pt->AddEntry(gr3, "Alpha+2", "l");
	pt->AddEntry(gr4, "Alpha-2", "l");

	pt->Draw();

//Beta values//////////////////////////////////////////////////////////////////

	// Solid angle, Yield, Yield error
	Double_t bdomega, byield, berror;

	// Photon energy, Recoil angle, Differential cross section -2,
	// Differential cross section +2, Varied Standard Unit *2 ,Chi squared
	Double_t begamma, btheta, bdiffxssm, bdiffxsbg, bdiffxs, bdxs, bchi_sq, bchi_sq_sm, bchi_sq_bg;

	// Theta for line, Theta error for line, Theta for points,
	// Theta error for points, Differential cross section (-2),
	// Differential cross section (-2) error,
	// Differential cross section (+2),
	// Differential cross section (+2) error,
	// Differential cross section for points,
	// Differential cross section error for points,
	// Differential cross section for line,
	// Differential cross section error for line	
	Double_t bthet1[20], bdthet1[20], bdixssm1[20], bddixssm1[20], bdixsbg1[20], bddixsbg1[20], bdixs1[20], bddixs1[20];

	// Point values
	// Theta, Theta error,
	// Differential cross section (-2),
	// Differential cross section (-2) error,
	// Differential cross section (+2),
	// Differential cross section (+2) error,
	// Differential cross section,
	// Differential cross section error	
	Double_t bthet2[20], bdthet2[20], bdixssm2[20], bddixssm2[20], bdixsbg2[20], bddixsbg2[20], bdixs2[20], bddixs2[20];

	// Smeared DXS values
	Double_t bdixsran[20];

	Int_t m, n;

	// Data file to be read in
	TString nbeta = "3HeXsec_neutron_betavaried.dat";	
	
 	Double_t deteff = 1;
	if (e == 100) 
	{
	        deteff = 0.326; //My values
	}	
	else if (e == 120) 
	{
	        deteff = 0.818; //Meg's values
	}

	// Ensure data file is open and can be read in 
	ifstream bet(nbeta);
	if ( !bet.is_open()) 
	{
		cout << "Error opening file ";
		cout << nbeta;
		cout << endl;
		break;
	}

	// Set all DXS values to be within range of DXS low and high
	TF1 *f2 = new TF1("f2", "gaus", lowDXS, highDXS);

	// Random amplitude
	f2->SetParameter(0,10);

	// Zero all variables to prevent errors
	bchi_sq = 0;
	m = 0;
	n = 0;

	// Read through data until the end of file
	while(!bet.eof())
	{
	        // Read each line into appropriate variables
		bet >> begamma >> btheta >> bdiffxssm >> bdiffxs >> bdiffxsbg >> bdxs;

		// Examine lines pertaining to the set energy
		if(begamma == e)
		{

		        // Error for all the following values are 0
		        // Set x data = theta
			bthet1[m] = btheta;
			bdthet1[m] = 0;

			// Set y data line 3 = DXS -2
			bdixssm1[m] = bdiffxssm;
			bddixssm1[m] = 0;

			// Set y data line 2 = DXS +2
			bdixsbg1[m] = bdiffxsbg;
			bddixsbg1[m] = 0;

			// Set y data line 4 = DXS
			bdixs1[m] = bdiffxs;
			bddixs1[m] = 0;

			n = m/2;

			// For angles 20 -160 incremented by 20
			if ( m%2 == 0 && n!=0 && n!=9)
			{

			        // Set x data = theta
			        // Theta error = 0
				bthet2[n] = btheta;
				bdthet2[n] = 0;

				// Set y data line 3 = DXS -2
				bdixssm2[n] = bdiffxssm;
				
				// Set y data line 2 = DXS +2
				bdixsbg2[n] = bdiffxsbg;

				// Set y data line 1 = DXS
				bdixs2[n] = bdiffxs;
				
				// Solid angle
				bdomega = 2*kPI*(cos((bthet2[n]-10)*kPI/180)-cos((bthet2[n]+10)*kPI/180));
				
				// Yield and Error
				byield = bdixs2[n]*photonflux*deteff*livet*targt*bdomega*beamt;
				berror = sqrt(byield);

				//DXS uncertainty
				bddixs2[n] = bdixs2[n]/berror;
				bddixssm2[n] = bdixssm2[n]/berror;
				bddixsbg2[n] = bdixsbg2[n]/berror;

				// Smear DXS values
				// Main = central DXS
				f2->SetParameter(1,bdixs2[n]);

				// St. Dev. = DXS error
				f2->SetParameter(2,bddixs2[n]);

				// Randomize point location
				bdixsran[n] = f2->GetRandom();

				cout << bthet2[n] << " " << bdixs2[n] << " " << bdixsran[n] << endl;	       

				//Calculates the chi-sq for this point
				bchi_sq += pow((bdixsran[n] - bdixs2[n])/bddixs2[n], 2);
				bchi_sq_sm += pow((bdixsran[n] - bdixssm2[n])/bddixssm2[n], 2);
				bchi_sq_bg += pow((bdixsran[n] - bdixsbg2[n])/bddixsbg2[n], 2);
			     

			}	
			m++;	
		}
	}
	bet.close();

	// Calculate chi squared
	if (n > 1){
	        bchi_sq /= n-2;
		bchi_sq_sm /= n-2;
		bchi_sq_bg /= n-2;
	}
	  
	// For unacceptable values set Chi-squared to zero to prevent error
	else if ( (n == 1) || ( n == 0)){
	        bchi_sq = 0;
		bchi_sq_sm = 0;
		bchi_sq_bg = 0;
	}

	// Print Chi-squared for beta
	cout << "Beta Reduced Chi Squared = " << bchi_sq << endl;
	cout << "Beta Reduced Chi Squared -2 = " << bchi_sq_sm << endl;
	cout << "Beta Reduced Chi Squared +2 = " << bchi_sq_bg << endl;

	// Plot trend line for beta unvaried
	TGraphErrors *gr5 = new TGraphErrors( m, bthet1, bdixs1, bdthet1, bddixs1);
	TString titleb = Form("Beta Varied Differential Cross Sections for E_{#gamma}=%3.0f MeV", e);
	gr5->SetTitle(titleb);
	gStyle->SetTitleX(0.5);
	gStyle->SetTitleAlign(23);
	gr5->SetLineWidth(3);
	gr5->SetLineStyle(1);
	gr5->GetXaxis()->SetTitleOffset(1.1);
	gr5->GetYaxis()->SetTitleOffset(1.0);
	gr5->GetXaxis()->SetTitle("#theta (deg)");
	gr5->GetYaxis()->SetTitle("d#sigma/d#Omega (nb/sr)");
	gr5->GetXaxis()->SetLabelSize(0.03);
	gr5->GetYaxis()->SetLabelSize(0.03);
	gr5->GetXaxis()->CenterTitle();
	gr5->GetYaxis()->CenterTitle();
	gr5->GetXaxis()->SetRangeUser(0,180);
	gr5->GetYaxis()->SetRangeUser( lowDXS, highDXS);

	// Plot points with error bars for beta unvaried, every 20 deg ** smeared ** 
	TGraphErrors *gr6 = new TGraphErrors( m, bthet2, bdixsran, bdthet2, bddixs2);
	gr6->SetTitle();
	gr6->SetMarkerColor(2);
	gr6->SetMarkerStyle(20);
	gr6->SetMarkerSize(1.1);
	gr6->SetLineWidth(2);

	// Plot trend line for beta varied by +2
	TGraphErrors *gr7 = new TGraphErrors( m, bthet1, bdixsbg1, bdthet1, bddixsbg1);
	gr7->SetTitle();
	gr7->SetLineWidth(3);
	gr7->SetLineStyle(9);

	// Plot trend line for beta varied by -2
	TGraphErrors *gr8 = new TGraphErrors( m, bthet1, bdixssm1, bdthet1, bddixssm1);
	gr8->SetTitle();
	gr8->SetLineWidth(3);
	gr8->SetLineStyle(2);
	c1->cd(2);

	gr5->Draw("AL");
	gr6->Draw("same" "P");
	gr7->Draw("same");	
	gr8->Draw("same");

	pt2 = new TLegend(0.6, 0.7, 0.85, 0.85);
	pt2->SetFillColor(0);
	pt2->SetBorderSize(0);
	pt2->SetTextSize(0.04);

	pt2->AddEntry(gr6, "Beta", "lp");
	pt2->AddEntry(gr7, "Beta+2", "l");
	pt2->AddEntry(gr8, "Beta-2", "l");

	pt2->Draw();

}


