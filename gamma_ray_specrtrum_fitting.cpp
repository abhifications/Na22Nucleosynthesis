#include <iostream>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TMath.h"
#include "Math/SpecFuncMathMore.h"

using namespace std;

void plot_GePD_new() 
{
    // Prompt for filename.
    cout << "Enter name of file: ";
    string filename;
    cin >> filename;

    // Create canvas and set style options.
    TCanvas* c1 = new TCanvas("c1", "plots", 1600, 1000);
    gStyle->SetOptStat(1);

    // Open the ROOT file and retrieve the original histogram.
    TFile* file = new TFile(filename.c_str());
    TH1D* hOriginal = (TH1D*)file->Get("GePD");

    // Create a fresh histogram for the raw spectrum (16384 bins).
    TH1D* GePD = new TH1D("GePD", "GePD", 16384, 0, 16384);
    for(int bin = 1; bin <= 16384; bin++) {
        double binVal = hOriginal->GetBinContent(bin);
        GePD->SetBinContent(bin, binVal);
    }

    // Apply a linear calibration with slope & intercept, then produce a calibrated histogram.
    double intercept = -49.5375;
    double slope     =  0.70384;
    int nbins = GePD->GetXaxis()->GetNbins();
    double xMin = GePD->GetXaxis()->GetBinCenter(1);
    double xMax = GePD->GetXaxis()->GetBinCenter(nbins);
    double calMin = intercept + slope * xMin;
    double calMax = intercept + slope * xMax;

    TH1D* hCalibrated = new TH1D("hnew", "calib", nbins, calMin, calMax);
    for(int i = 1; i <= nbins; i++) {
        double content = GePD->GetBinContent(i);
        double ch      = GePD->GetXaxis()->GetBinCenter(i);
        double energy  = intercept + slope * ch;
        // Fill the calibrated histogram at the corresponding bin:
        // (Note: for strictly correct binning you would map 'energy' to bin index, 
        // but here the code copies the content as if same bin index.)
        hCalibrated->SetBinContent(i, content);
    }

    // Log-scale for Y-axis.
    c1->SetLogy();

    // Label axes and set titles for both raw & calibrated.
    hCalibrated->GetXaxis()->SetTitle("Energy_{#gamma} [keV]");
    hCalibrated->GetYaxis()->SetTitle("Counts");
    hCalibrated->GetXaxis()->CenterTitle();
    hCalibrated->GetYaxis()->CenterTitle();
    hCalibrated->SetTitle("GePD (Calibrated)");

    GePD->GetXaxis()->SetTitle("Channel number");
    GePD->GetYaxis()->SetTitle("Counts");
    GePD->GetXaxis()->CenterTitle();
    GePD->GetYaxis()->CenterTitle();
    GePD->SetTitle("GePD (Raw)");

    // Draw the calibrated histogram by default (comment/uncomment as desired).
    hCalibrated->Draw();

    // -----------------------------
    // Example: Distorted Gaussian fit around 2287 keV region.
    // -----------------------------
    {
        // Define a 6-parameter Distorted Gaussian + linear background.
        Double_t par[6];
        TF1* fitDistGauss = new TF1("g1",
            "([0]/(2*[3]))*exp((x-[1])/[3] + [2]^2/(2*[3]^2)) * erfc( ( (x-[1])/[2] + [2]/[3] )/sqrt(2) ) + [4] + [5]*x",
            3308, 3329);
        fitDistGauss->SetParameter(0, 1400.0);   // net area
        fitDistGauss->SetParLimits(0, 1, 1e7);
        fitDistGauss->SetParameter(1, 3319.0);   // mean
        fitDistGauss->SetParLimits(1, 3319-5, 3319+5);
        fitDistGauss->SetParameter(2, 4.0);      // sigma
        fitDistGauss->SetParLimits(2, 2, 8);
        fitDistGauss->SetParameter(3, 2.0);      // tail parameter
        fitDistGauss->SetParLimits(3, 1, 8);
        fitDistGauss->SetParameter(4, 40.0);     // constant background
        fitDistGauss->SetParameter(5, 0.0);      // linear slope of background

        fitDistGauss->SetLineColor(2);
        GePD->Fit(fitDistGauss, "R"); // Fit the raw spectrum in that region (channel domain).

        fitDistGauss->GetParameters(&par[0]);
        double chisq = fitDistGauss->GetChisquare();
        double ndf   = fitDistGauss->GetNDF();
        cout << "[Distorted Gaussian Fit] Chi2/ndf = " << chisq/ndf << endl;
    }

    // -----------------------------
    // Background-subtraction method (Gilmore) around peak near 2287 keV
    // -----------------------------
    {
        // Window definitions in channel space:
        int E1L = 3292; // left bkg region start
        int E2L = 3308; // left bkg region end
        int EL  = 3312; // peak region start
        int ER  = 3327; // peak region end
        int E1R = 3330; // right bkg region start
        int E2R = 3345; // right bkg region end

        // Convert to bin indices:
        TAxis* xax   = GePD->GetXaxis();
        int binL     = xax->FindBin(EL);
        int binR     = xax->FindBin(ER);
        int bin1L    = xax->FindBin(E1L);
        int bin2L    = xax->FindBin(E2L);
        int bin1R    = xax->FindBin(E1R);
        int bin2R    = xax->FindBin(E2R);

        // Summations.
        double rawPeakSum = GePD->Integral(binL, binR);
        double errRawPeak = sqrt(rawPeakSum);

        double bckLeftSum  = GePD->Integral(bin1L, bin2L);
        double bckRightSum = GePD->Integral(bin1R, bin2R);

        int N_L = (bin2L - bin1L + 1);
        int N_R = (bin2R - bin1R + 1);
        int N_P = (binR - binL + 1);

        // Centroids to determine weighting factors.
        double peakCentroid = 0;
        double totalPeak    = 0;
        for(int iBin = binL; iBin <= binR; iBin++){
            double content = GePD->GetBinContent(iBin);
            peakCentroid  += content * iBin;
            totalPeak     += content;
        }
        if(totalPeak > 0) {
            peakCentroid /= totalPeak;
        }

        double leftBckCentroid = 0;
        double totalLeft       = 0;
        for(int iBin = bin1L; iBin <= bin2L; iBin++){
            double content = GePD->GetBinContent(iBin);
            leftBckCentroid += content * iBin;
            totalLeft       += content;
        }
        if(totalLeft > 0) {
            leftBckCentroid /= totalLeft;
        }

        double rightBckCentroid = 0;
        double totalRight       = 0;
        for(int iBin = bin1R; iBin <= bin2R; iBin++){
            double content = GePD->GetBinContent(iBin);
            rightBckCentroid += content * iBin;
            totalRight       += content;
        }
        if(totalRight > 0) {
            rightBckCentroid /= totalRight;
        }

        // Weighted average background (Gilmore).
        double distLeft  = fabs(peakCentroid - leftBckCentroid);
        double distRight = fabs(peakCentroid - rightBckCentroid);

        double wL = distRight / (distLeft + distRight);
        double wR = distLeft  / (distLeft + distRight);

        double bckPerBin = (bckLeftSum / N_L) * wL + (bckRightSum / N_R) * wR;
        double totalBck  = bckPerBin * N_P;
        double errBck    = sqrt(totalBck);

        // Net area & error in that region.
        double netPeak    = rawPeakSum - totalBck;
        double errNetPeak = sqrt(pow(errRawPeak, 2) + pow(errBck, 2));

        // Print results:
        cout << "Peak Raw Counts: " << rawPeakSum 
             << " +/- " << errRawPeak << endl;
        cout << "Net Peak Counts: " << netPeak 
             << " +/- " << errNetPeak << endl;
    }

    // Optionally draw the raw spectrum on top if needed:
    // GePD->Draw("same");
}
