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
	
	// Differential cross section, Differential cross section error, Theta
	Double_t asmear[20], aserr[20], asthet[20];

	// Change in Cross section, New alpha value
	Double_t adeltaxs[20], axsnew[20];

	// Angle, Chi-Squared
	Int_t athet[20];

	Double_t achi_sqa[40], alphaa[40];

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
		asthet[k] = astheta;
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
		if(aegamma == e && atheta%20 == 0)
		{

		        // x data = theta
			athet[i] = atheta;

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
	for (alpha = 11.0; alpha <= 14.0; alpha = alpha + 0.1)
	{

	        alphaa[l] = alpha;
		
	        // Calculate Chi-squared for each angle
		for (j = 0; j < 8; j++)
		{

		        //I don't think these are right
		        axsnew[j] = asmear[j] + ((alpha-alpha0)/2)*adeltaxs[j];
		
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
		l++;
	}

	TGraphErrors *gr1 = new TGraphErrors(l, alphaa, achi_sqa);
	gr1->GetXaxis()->SetRangeUser( 11, 14);
	gr1->GetYaxis()->SetRangeUser( 0, 6);

	TCanvas *c1 = new TCanvas("c1", "Errors", 0, 0, 500,500);

	gr1->Draw("AL");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Int_t l, m, n;

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
	// Chi-squared, error on alphs
	Double_t beta0, beta, bchi_sq, beterr;

	// Differential cross section, Differential cross section error, Theta
	Double_t bsmear[20], bsthet[20], bserr[20];

	// Change is cross section, new alpha value
	Double_t bdeltaxs[20], bxsnew[20];

	// Angle
	Int_t bthet[20];

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
	l = 0;
	while(!bsmearf.eof())
	{
		bsmearf >> bstheta >> bs >> bserror;
		bsthet[l] = bstheta;
		bsmear[l] = bs;
		bserr[l] = bserror;
		l++;
	}
	bsmearf.close();

	// Read through theory DXS data file
	m = 0;
	while(!bet.eof())
	{
		bet >> begamma >> btheta >> bdiffxssm >> bdiffxs >> bdiffxsbg >> bdxs;

		// Place variables from lines into arrays
		if(begamma == e && btheta%20 == 0)
		{

		        // x data = theta
			bthet[m] = btheta;

			// amount DXS is varied +/- 2
			bdeltaxs[m] = bdxs;   
			m++;
		}			
	}
	bet.close();

	bchi_sq = 0;

	// Choose central value of accepted beta amount
	beta0 = 1.25;

	// Determine Chi-squared for alpha values +/- 1.5 from beta0
	for (beta = 0.00; beta <= 2.50; beta = beta + 0.15)
	{

	        // Calculate Chi-squared for each angle
		for (n = 0; n < 8; n++)
			{
				bxsnew[n] = bsmear[n] + ((beta-beta0)*bdeltaxs[n])/2;
				
				//Calculates the chi-sq for this point
				bchi_sq += pow((bxsnew[n] - bsmear[n])/bserr[n], 2);

			}

		// Calculate chi squared
		if (n > 1)
		{
			bchi_sq /= n-2;
		}	
		else if ( (n == 1) || ( n == 0))
		{
			bchi_sq = 0;
		}

		beterr = beta0 - beta;

		cout << "Beta = " << beta << endl;
		cout << "Beta Error = " << beterr << endl;
		cout << "Beta Reduced Chi Squared  = " << bchi_sq << endl;

	}

*/
}
