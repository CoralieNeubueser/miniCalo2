#include "TSystem.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TEventList.h"
#include "TStyle.h"
#include "THStack.h"
#include "TPad.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TLatex.h"
#include "TObject.h"

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>


using namespace std;

double showDecimals(const double& x, const int&numDecimals){
  int y=x;
  double z=x-y;
  double m=pow(10,numDecimals);
  double q=z*m;
  double r=round(q);

  return static_cast<double>(y)+(1.0/m)*r;
}

void fitGauss(TH1 *hErec, double &fitResponse, double &fitResponseErr, double &fitReso, double &fitResoErr, bool addFitFuncToHisto) {
  string fitOption = addFitFuncToHisto?"Q":"QN";
 
  TF1 *fGauss = new TF1("fGauss","gaus",0,100);
  double constantGuess = hErec->GetBinContent(hErec->GetMaximumBin()); //guess for gauss constant is content of maximum bin
  double meanGuess = hErec->GetMean(1); //guess for gauss mean is mean of histo
  double sigmaGuess = hErec->GetRMS(1); //guess for gauss sigma is RMS of histo
  fGauss->SetParameters(constantGuess, meanGuess, sigmaGuess);

  hErec->Fit(fGauss,fitOption.c_str(),"",meanGuess-2*sigmaGuess, meanGuess + 2*sigmaGuess);

  meanGuess = fGauss->GetParameter(1);
  sigmaGuess= fGauss->GetParameter(2);
  hErec->Fit(fGauss,fitOption.c_str(),"",meanGuess-2*sigmaGuess, meanGuess + 2*sigmaGuess);

  meanGuess = fGauss->GetParameter(1);
  sigmaGuess= fGauss->GetParameter(2);
  hErec->Fit(fGauss,fitOption.c_str(),"",meanGuess-2*sigmaGuess, meanGuess + 2*sigmaGuess);

  double mean = fGauss->GetParameter(1);
  double meanErr = fGauss->GetParError(1);
  double sigma = fGauss->GetParameter(2);
  double sigmaErr = fGauss->GetParError(2);

  if (addFitFuncToHisto){
    TF1* fGaussFull=new TF1(*fGauss);
    fGaussFull->SetName("fitNovoFull");
    fGaussFull->SetLineStyle(2);
    fGaussFull->SetRange(mean-4*sigma, mean+4*sigma);
    hErec->GetListOfFunctions()->Add(fGaussFull);
  }
  //hErec->GetListOfFunctions()->Add(fGauss);
  fitReso = sigma/mean * 100.;
  fitResoErr = 1/mean * TMath::Sqrt(TMath::Power(sigmaErr,2)+TMath::Power((fitReso/100.*meanErr),2)) * 100.;
  fitResponse = mean;
  fitResponseErr = meanErr;
}

void resolution(string name) {

  double eBeam[10]={10,20,30,40,50,60,70,80,90,100};
  double eBeamErr[10]={0};

  double Mean[10] = {0};
  double MeanErr[10] = {0};
  double Sigma[10] = {0};
  double SigmaErr[10] = {0};

  THStack *allHists = new THStack();

  std::stringstream path;
  path << "/eos/user/c/cneubuse/miniCalo2/results/"<< name; 
  
  std::stringstream inputFileName;
  inputFileName << path.str() << "/response.root";
  
  TFile* inputFile = TFile::Open(inputFileName.str().c_str(), "update");
  
  for (unsigned int i=0; i<10; i++){
    std::cout << "Energy: " << eBeam[i] << std::endl;
    
    std::stringstream energy;
    energy << eBeam[i] << " GeV";
    
    TH1D* E = (TH1D*)inputFile->Get(energy.str().c_str());
    E->Sumw2();
    E->Scale(1./E->Integral());
    
    fitGauss(E, Mean[i], MeanErr[i], Sigma[i], SigmaErr[i], true);
    
    E->GetFunction("fitNovoFull")->SetLineColor(E->GetLineColor());
    E->GetFunction("fGauss")->SetLineColor(E->GetLineColor());
    
    E->Write(E->GetName(),TObject::kOverwrite);
    
    std::stringstream energyTitle;
    energyTitle << eBeam[i] << "GeV #pi^{-} @ #eta=0.36";
    
    std::stringstream printName; 
    printName << "~/root_plots/FccHcal/Clustering/simProd_v03_topoClusters_calibratedPions_eta036_" << eBeam[i] << "GeV.pdf";     

    TLatex *Tright = new TLatex(0.2, 0.95, energyTitle.str().c_str()); 
    Tright->SetNDC(kTRUE) ;
    Tright->SetTextSize(0.044);
    Tright->SetTextFont(132);

    TLegend *legende = new TLegend(0.22,0.75,0.55,0.9);
    legende->SetTextFont(132) ;
    legende->AddEntry(E, energyTitle.str().c_str(),"l");
    
    allHists->Add(E, "hist l");
    
  }

  TGraphErrors* lin = new TGraphErrors(10, eBeam, Mean, eBeamErr, MeanErr);
  lin->SetTitle("linearity");
  lin->SetMarkerStyle(22);  
  lin->SetMarkerSize(1.5);
  lin->SetMarkerColor(1);

  lin->Write(lin->GetTitle(),TObject::kOverwrite);

  TGraphErrors* reso = new TGraphErrors(10, eBeam, Sigma, eBeamErr, SigmaErr);
  reso->SetTitle("resolution");
  reso->SetMarkerStyle(22);
  reso->SetMarkerSize(1.5);
  reso->SetMarkerColor(1);
  reso->SetLineColor(1);

  reso->Write(reso->GetTitle(),TObject::kOverwrite);

  TF1 *fitReso = new TF1("fitReso", "sqrt(pow([0]/sqrt(x),2)+pow([1],2))", 20,90);
  fitReso->SetLineColor(1);
  reso->Fit("fitReso","R");

  TMultiGraph *multi = new TMultiGraph();
  multi->Add(reso);
  TMultiGraph *multiLin = new TMultiGraph();
  multiLin->Add(lin);

  std::cout << "Resolution for pions (PbSpacer) : a = " << fitReso->GetParameter(0) << " , c = " << fitReso->GetParameter(1) << std::endl;

  std::stringstream fitResults0;
  fitResults0 << "#frac{" << showDecimals(fitReso->GetParameter(0),0) << "%#sqrt{GeV}}{#sqrt{E}} #oplus "<< showDecimals(fitReso->GetParameter(1),1) << "%";

  TLegend *legendeNew = new TLegend(.4,.6,.95,.95);
  legendeNew->SetFillColor(0);
  legendeNew->SetTextFont(132);
  legendeNew->SetTextSize(0.05);
  legendeNew->AddEntry(reso,     fitResults0.str().c_str(),"pl");
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TCanvas *c2 = new TCanvas("c2","resolutions+linearity",600,700);
  TPad *pad3 = new TPad("pad3","pad3",0,0,1,0.66);
  TPad *pad4 = new TPad("pad4","pad4",0,0.66,1,1);
  pad4->SetBottomMargin(0.01);
  pad3->SetBorderMode(0);
  pad3->SetTopMargin(0.01);
  pad3->SetBottomMargin(0.15);
  pad4->SetBorderMode(0);
  pad3->SetTickx(1);
  pad4->SetTickx(1);
  pad3->SetTicky(1);
  pad4->SetTicky(1);
  pad3->Draw();
  pad4->Draw();

  pad3->cd();
  multi->Draw("AP");
  multi->GetXaxis()->SetTitle("E [GeV]");
  multi->GetYaxis()->SetTitle("#sigma_{E_{rec}}/#LT E_{rec}#GT [%]");
  multi->GetXaxis()->SetLimits(0,100);
  multi->SetMinimum(0);
  multi->SetMaximum(20);
  multi->Draw("AP");
  legendeNew->Draw();

  pad4->cd();
  pad4->SetGridy();
  multiLin -> Draw("AP");
  multiLin -> GetXaxis() -> SetTitle ("E [GeV]");
  multiLin -> GetYaxis() -> SetTitle ("#LTE_{rec}/E_{true}#GT");
  multiLin -> GetYaxis() -> SetLabelSize(0.09);
  multiLin -> GetYaxis() -> SetTitleOffset(0.9);
  multiLin -> GetYaxis() -> SetTickLength(0.04);
  multiLin -> GetYaxis() -> SetTitleSize(0.11);
  multiLin -> GetXaxis() -> SetTickLength(0.04);
  multiLin -> GetXaxis() -> SetTitleOffset(1.);
  multiLin -> GetXaxis() -> SetLabelSize(0.1);
  multiLin -> GetXaxis() -> SetTitleSize(0.09);
  multiLin -> GetXaxis() -> SetLimits(0,100);
  multiLin -> GetYaxis() -> SetRangeUser(0.95,1.05);
  multiLin -> Draw("AP");

  c2->Write(c2->GetName(),TObject::kOverwrite);
  std::stringstream outfile;
  outfile << path.str() << "/resolution.pdf";
  c2->Print(outfile.str().c_str());
}
