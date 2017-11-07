gROOT->Reset();

#include "/home/maeve/Code/lib/physics.h"
#include "/home/maeve/Code/lib/functions.h"

void AlphaBeta(Double_t e = 120)
{

        Int_t i, j, k, l;
	
	// Theta, differential cross section data,
	// differential cross section error
	Double_t astheta, as, aserror;
	
	// Energy, differential cross section -2,
	// Differential cross section +2, Differential cross section,
	// Amount varied to get +/-2
	Double_t aegamma, adiffxssm, adiffxsbg, adiffxs, adxs;
	
	//Angle
	Int_t atheta;
	
	// Accepted alpha value, varied alpha value,
	// chi-squared, error on alpha
	Double_t  alpha0, alpha, achi_sq, alpherr;
	
	// Differential cross section smeared value,
	// Differential cross section error
	Double_t asmear[20], aserr[20];

	// Change in Cross section, New alpha value,
	// Central Differential Cross section 
	Double_t adeltaxs[20], axsnew[20], adixs[20];

	// Chi-Squared array, Alpha array
	Double_t achi_sqa[100], alphaa[100];

	// Simulated cross section data
	TString nalpha = "3HeXsec_neutron_alphavaried.dat";

	// Smeared data
	TString alphfile = Form("alpha_smear_%3.0f.dat", e);
	
	// Check if file (nalpha) is open and readable
	ifstream alph(nalpha);
	if ( !alph.is_open()) 
	{
		cout << "Error opening file ";
		cout << nalpha;
		cout << endl;
		break;
	}

	// Check if file (alphfile) is open and readable
	ifstream asmearf(alphfile);
	if ( !asmearf.is_open())
	{
		cout << "Error opening file ";
		cout << alphfile;
		cout << endl;
		break;
	}

	// Read through smeared data file, placing variables into arrays
	k = 0;
	while(!asmearf.eof())
	{
		asmearf >> astheta >> as >> aserror;
		asmear[k] = as;
		aserr[k] = aserror;
		k++;
	}
	asmearf.close();

	// Read through theory DXS data file
	i = 0;
	while(!alph.eof())
	{
		alph >> aegamma >> atheta >> adiffxssm >> adiffxs >> adiffxsbg >> adxs;

		// Place variables from lines into arrays
		if(aegamma == e && atheta%20 == 0 && atheta != 0)
		{

			// y data = DXS
			adixs[i] = adiffxs;

			// amount DXS is varied +/-2
			adeltaxs[i] = adxs;
			i++;
		}			
	}
	alph.close();

	achi_sq = 0;

	// Choose central value of accepted alpha amount
	alpha0 = 12.5;
	l = 0;
	
	// Determine Chi-Squared for alpha values +/- 1.5 from alpha0
	for(alpha = 9.0; alpha <= 15.0; alpha = alpha + 0.1)
	{

	        alphaa[l] = alpha;
		
	        // Calculate Chi-squared for each angle
		for (j = 0; j < 8; j++)
		{

		        axsnew[j] = adixs[j] + ((alpha-alpha0)/2)*adeltaxs[j];
		
			// Calculates the chi-sq for this point
			achi_sq += pow((axsnew[j] - asmear[j])/aserr[j], 2);
		}
		
		// Calculate chi squared for all points
		if (j > 1)
		{
		        achi_sqa[l] = achi_sq/(j-2);
		}	
		else if ( (j == 1) || ( j == 0))
		{
			achi_sqa[l] = 0;
		}

		// Need to check this value
	       	alpherr = alpha0 - alpha;

		cout << "Alpha = " << alpha << endl;
	       	cout << "Alpha Error = " << alpherr << endl;
		cout << "Alpha Reduced Chi Squared  = " << achi_sqa[l] << endl;
		
		achi_sq = 0;
		l++;
		
	}
	
	TString titlea = Form("Alpha Varied Chi-Squared Values for 100 MeV", e);
	
	TGraph *gr1 = new TGraph(l, alphaa, achi_sqa);
	gr1->SetTitle(titlea);
	gr1->SetLineWidth(3);
	gr1->SetLineColor(kRed);
	gr1->GetXaxis()->SetTitle("Alpha (fm^{3})");
	gr1->GetYaxis()->SetTitle("#chi^{2}_{red}");
	gr1->GetXaxis()->SetTitleSize(0.05);
	gr1->GetYaxis()->SetTitleSize(0.05);
	gr1->GetXaxis()->SetLabelSize(0.03);
	gr1->GetYaxis()->SetLabelSize(0.03);
	gr1->GetXaxis()->SetTitleOffset(0.8);
	gr1->GetYaxis()->SetTitleOffset(1.0);
	gr1->GetXaxis()->CenterTitle();
	gr1->GetYaxis()->CenterTitle();
	gr1->GetXaxis()->SetRangeUser( 9, 15);
	gr1->GetYaxis()->SetRangeUser( 1, 2.5);

	TLine *line1 = new TLine(9, 2, 15, 2);
	line1->SetLineWidth(2);
	
	TCanvas *c1 = new TCanvas("c1", "Errors", 0, 0, 800,400);
       	c1->Divide(2);
	c1->cd(1);
	
	gr1->Draw("AL");
	line1->Draw("clone");
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Int_t m, n, o, p;

	// Theta, differential cross section data,
	//differential cross section error
	Double_t bstheta, bs, bserror;

	// Energy, Differential cross section -2,
	// Differential cross section +2, Differential cross section,
	// Amount varied to get +/- 2
	Double_t begamma, bdiffxssm, bdiffxsbg, bdiffxs, bdxs;

	// Angle
	Int_t btheta;

	// Accepted alpha value, varied alpha value
	// Chi-squared, error on alpha
	Double_t beta0, beta, bchi_sq, beterr;

	// Differential cross section ameared value,
	// Differential cross section error
	Double_t bsmear[20], bserr[20];

	// Change is cross section, new alpha value,
	// Central Differential Cross section
	Double_t bdeltaxs[20], bxsnew[20], bdixs[20];

	// Chi-Squared array, beta array
	Double_t bchi_sqa[100], betaa[100];

	// Simulated cross section data
	TString nbeta = "3HeXsec_neutron_betavaried.dat";

	// Smeared data
	TString betafile = Form("beta_smear_%3.0f.dat", e);  

	// Check if file (nbeta) is open and readable
	ifstream bet(nbeta);
	if ( !bet.is_open()) 
	{
		cout << "Error opening file ";
		cout << nbeta;
		cout << endl;
		break;
	}

	// Check if file (betafile) is open and readable
	ifstream bsmearf(betafile);
	if ( !bsmearf.is_open())
	{
		cout << "Error opening file ";
		cout << betafile;
		cout << endl;
		break;
	}

	// Read through smeared data file, placing variables into arrays
	m = 0;
	while(!bsmearf.eof())
	{
		bsmearf >> bstheta >> bs >> bserror;
		bsmear[m] = bs;
		bserr[m] = bserror;
		m++;
	}
	bsmearf.close();

	// Read through theory DXS data file
	n = 0;
	while(!bet.eof())
	{
		bet >> begamma >> btheta >> bdiffxssm >> bdiffxs >> bdiffxsbg >> bdxs;

		// Place variables from lines into arrays
		if(begamma == e && btheta%20 == 0 && btheta != 0)
		{

		        // ydata = DXS
		        bdixs[n] = bdiffxs;

			// amount DXS is varied +/- 2
			bdeltaxs[n] = bdxs;   
			n++;
		}			
	}
	bet.close();

	bchi_sq = 0;

	// Choose central value of accepted beta amount
	beta0 = 1.25;
	p=0;

	// Determine Chi-squared for alpha values +/- 1.5 from beta0
	for (beta = -1.5; beta <= 4.4; beta = beta + 0.15)
	{

	        betaa[p] = beta;
	  
	        // Calculate Chi-squared for each angle
		for (o = 0; o < 8; o++)
			{
				bxsnew[o] = bdixs[o] + ((beta-beta0)*bdeltaxs[o])/2;
				
				//Calculates the chi-sq for this point
				bchi_sq += pow((bxsnew[o] - bsmear[o])/bserr[o], 2);

			}

		// Calculate chi squared
		if (o > 1)
		{
		        bchi_sqa[p] = bchi_sq/(o-2);
		}	
		else if ( (o == 1) || ( o == 0))
		{
			bchi_sqa[p] = 0;
		}

		beterr = beta0 - beta;

		cout << "Beta = " << beta << endl;
		cout << "Beta Error = " << beterr << endl;
		cout << "Beta Reduced Chi Squared  = " << bchi_sqa[p] << endl;
		
		bchi_sq = 0;
		p++;
	}

	TString titlea = Form("Beta Varied Chi-Squared Values for 100 MeV", e);
	
	TGraph *gr2 = new TGraph(p, betaa, bchi_sqa);
	gr2->SetTitle(titlea);
	gr2->SetLineWidth(3);
	gr2->SetLineColor(2);
	gr2->GetXaxis()->SetTitle("Beta (fm^{3})");
	gr2->GetYaxis()->SetTitle("#chi^{2}_{red}");
	gr2->GetXaxis()->SetLabelSize(0.03);
	gr2->GetYaxis()->SetLabelSize(0.03);
	gr2->GetXaxis()->SetTitleOffset(0.8);
	gr2->GetYaxis()->SetTitleOffset(1.0);
	gr2->GetXaxis()->SetTitleSize(0.05);
	gr2->GetYaxis()->SetTitleSize(0.05);
	gr2->GetXaxis()->CenterTitle();
	gr2->GetYaxis()->CenterTitle();
	gr2->GetXaxis()->SetRangeUser( -1.5, 4);
	gr2->GetYaxis()->SetRangeUser( 1, 2.5);

	TLine *line2 = new TLine(-1.5, 2, 4, 2);
	line2->SetLineWidth(2);
     
	c1->cd(2);
	gr2->Draw("AL");
	line2->Draw("clone");
	
}
