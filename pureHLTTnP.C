//root -b -q 'draw_effS.C()'

#include <TROOT.h>                        // access to gROOT, entry point to ROOT system
#include <TSystem.h>                      // interface to OS
#include <TFile.h>                        // file handle class
#include <TTree.h>                        // class to access ntuples
#include <TBenchmark.h>                   // class to track macro running statistics
#include <TH1D.h>                         // histogram class
#include <vector>                         // STL vector class
#include <iomanip>                        // functions to format standard I/O
#include <string>                         // C++ string class
#include <sstream>                        // class for parsing strings
#include <TRandom3.h>
#include <TGaxis.h>
#include "TLorentzVector.h"           // 4-vector class
#include "TGraphAsymmErrors.h"

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF1.h>

#include "RooFit.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"

#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooCMSShape.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooPlot.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "TText.h"
#include "RooFitResult.h"
#include "RooRealBinding.h"
#include "RooBrentRootFinder.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TPad.h"

#include <iostream>
#include <fstream>

using namespace RooFit ;

int layer(int info,int hitNr)
{
  int layersInfo = (info >> 4)&0xFF;
  int hitCounter=0;
  for(int bitNr=0;bitNr<8;bitNr++){
    int bit = 0x1 << bitNr;
    if((layersInfo&bit)!=0){
      if(hitCounter==hitNr) return bitNr+1;
      else hitCounter++;
    }
  }
  return -1;

}

int nrLayers(int info)
{
  return (info >> 12)&0xF;
}

int nrHits(int info)
{
  int hitCounter=0;
  int layersInfo = (info >> 4)&0xFF;
  for(int bitNr=0;bitNr<8;bitNr++){
    int bit = 0x1 << bitNr;
    if((layersInfo&bit)!=0) hitCounter++;
  }
  return hitCounter;
}

// Cut and cout
TGraph* draw_eff(TTree* tree, TString title, TString xTitle, int nrBins,float xMin,float xMax,const TString& var,const TString& sampleCuts, TString& cuts){

 TH1* pass = new TH1D("passTemp","pass",nrBins,xMin,xMax);
 TH1* all = new TH1D("allTemp","all",nrBins,xMin,xMax);
 pass->Sumw2();
 all->Sumw2();

 float nrPass = tree->Draw((var+">>passTemp"),sampleCuts+" && "+cuts,"goff");
 float nrAll = tree->Draw(var+">>allTemp",sampleCuts,"goff");

 all->SetDirectory(0);
 pass->SetDirectory(0);

 TGraphAsymmErrors* graph = new TGraphAsymmErrors(pass,all,"B");
 //if(!plotXErrOnAsymEff_){ 
 //  for(int pointNr=0;pointNr<graph->GetN();pointNr++){
 //    graph->SetPointEXhigh(pointNr,0);
 //    graph->SetPointEXlow(pointNr,0);
 //  }
 //}
 
 delete all;
 delete pass;
 
 std::cout <<"nrPass "<<nrPass<< " / "<<nrAll<<std::endl;
 
 graph->SetTitle(title+";"+xTitle+";Efficiency");
 return graph;
}

TH1* MakePassHist(TTree* tree, TString title, TString xTitle, int nrBins,float xMin,float xMax,const TString& var,const TString& sampleCuts, TString& cuts){

 TH1* pass = new TH1D("passTemp","pass",nrBins,xMin,xMax);
 pass->Sumw2();

 float nrPass = tree->Draw((var+">>passTemp"),sampleCuts+" && "+cuts,"goff");

 pass->SetDirectory(0);

 return pass;
}

TH1* MakeFailHist(TTree* tree, TString title, TString xTitle, int nrBins,float xMin,float xMax,const TString& var,const TString& sampleCuts, TString& cuts){
 
 TH1* fail = new TH1D("failTemp","fail",nrBins,xMin,xMax);
 fail->Sumw2();
 
 float nrPass = tree->Draw((var+">>failTemp"),sampleCuts+" && "+cuts,"goff");
 
 fail->SetDirectory(0);
 
 return fail;
}

void pureHLTTnP()
{

 TFile* infile1 = NULL;
 TString filename1 = "/Volumes/Samsung_T3/2018_HLT_Performances/PixelMatchStudy/CMSSW_10_1_5/ntuple/ntuple_CMSSW_10_1_5_BadPixModuleFix.root";
 infile1 = new TFile(filename1);

 TTree *tree=0;

 // for pixel matching filter
 TString var = "mass";
 TString sampleCuts = " tagHLTRegion==0 && mass>60 && mass<120 && probeHLT.et > 30 && probeHLT.et < 60 && probeHLTRegion==0 && (evtTrigs[0]&0x4000000)!=0 && (tagTrigs[2]&0x400000)!=0 &&  (probeTrigs[2]&0x8000)!=0";
 TString passingProbe = "(probeTrigs[2]&0x10000)!=0"; // pixel matching filter bit
 TString failingProbe = "(probeTrigs[2]&0x10000)==0"; // pixel matching filter bit

 TH1* passHist;
 TH1* failHist;
   
 tree = (TTree*)infile1->Get("hltTPTree");
 passHist = MakePassHist(tree, "Run 315322, Eff. w.r.t. HCalIso filter", "Mass", 60, 60., 120., var, sampleCuts, passingProbe); 
 failHist = MakePassHist(tree, "Run 315322, Eff. w.r.t. HCalIso filter", "Mass", 60, 60., 120., var, sampleCuts, failingProbe); 

 RooRealVar massP("massP","m_{ee} [GeV]", 60., 120.);
 RooRealVar massF("massF","m_{ee} [GeV]", 60., 120.);
 RooDataHist dsDataP("dsDataP","dsDataP",RooArgSet(massP),Import(*passHist));
 RooDataHist dsDataF("dsDataF","dsDataF",RooArgSet(massF),Import(*failHist));

 RooPlot* frameP = massP.frame(Name("passing"), Title("passing")) ;
 RooPlot* frameF = massF.frame(Name("failing"), Title("failing")) ;
 frameP->SetMaximum(1.2*passHist->GetMaximum());
 frameF->SetMaximum(1.2*failHist->GetMaximum());
 dsDataP.plotOn(frameP, MarkerColor(kBlack), MarkerStyle(21), XErrorSize(0), MarkerSize(1.2)) ;
 dsDataF.plotOn(frameF, MarkerColor(kBlack), MarkerStyle(21), XErrorSize(0), MarkerSize(1.2)) ;

 // Crystal ball function for data
 RooRealVar cbMeanDataP("cbMeanDataP","cbMean", -1., -4., 4.);
 RooRealVar cbSigmaDataP("cbSigmaDataP","cbSigma", 2., 0.5, 4.);
 RooRealVar cbAlphaDataP("cbAlphaDataP","cbAlpha", 1., 0.1, 50);
 RooRealVar cbNDataP("cbNDataP","cbN", 2., 0.2, 50);

 RooRealVar cbMeanDataF("cbMeanDataF","cbMean", -1., -4., 4.);
 RooRealVar cbSigmaDataF("cbSigmaDataF","cbSigma", 2., 0.5, 4.);
 RooRealVar cbAlphaDataF("cbAlphaDataF","cbAlpha", 1., 0.1, 50);
 RooRealVar cbNDataF("cbNDataF","cbN", 2., 0.2, 50);

 // BreitWigner 
 RooRealVar meanP("meanP","meanP",91.2);
 RooRealVar sigmaP("sigmaP","sigmaP",2.5);

 RooRealVar meanF("meanF","meanF",91.2);
 RooRealVar sigmaF("sigmaF","sigmaF",2.5);

 // RooCMSShape
 RooRealVar alphacmsP("alphacmsP","alphacmsP", 60., 50., 80.);
 RooRealVar betacmsP("betacmsP","betacmsP", 0.05,0.01,0.2);
 RooRealVar gammacmsP("gammacmsP","gammacmsP", 0.1, 0, 1);
 RooRealVar peakcmsP("peakcmsP","peakcmsP", 90.0);

 RooRealVar alphacmsF("alphacmsF","alphacmsF", 60., 50., 80.);
 RooRealVar betacmsF("betacmsF","betacmsF", 0.05,0.01,0.2);
 RooRealVar gammacmsF("gammacmsF","gammacmsF", 0.1, 0, 1);
 RooRealVar peakcmsF("peakcmsF","peakcmsF", 90.0);

 RooRealVar nsigP("nsigP","signal events1",10,0.,100000.);
 RooRealVar nbkgP("nbkgP","signal background events1",10,0.,100000.);

 RooRealVar nsigF("nsigF","signal events1",10,0.,100000.);
 RooRealVar nbkgF("nbkgF","signal background events1",10,0.,100000.);

 RooBreitWigner bwP("bwP","bwP", massP, meanP, sigmaP);
 RooCBShape cbDataP("cbDataP","cb", massP, cbMeanDataP, cbSigmaDataP, cbAlphaDataP, cbNDataP);
 RooFFTConvPdf BWxCBDataP("BWxCBDataP", "BW (x) CB", massP, bwP, cbDataP);
 RooCMSShape bgP("bgP", "bgP", massP, alphacmsP, betacmsP, gammacmsP, peakcmsP);

 RooBreitWigner bwF("bwF","bwF", massF, meanF, sigmaF);
 RooCBShape cbDataF("cbDataF","cb", massF, cbMeanDataF, cbSigmaDataF, cbAlphaDataF, cbNDataF);
 RooFFTConvPdf BWxCBDataF("BWxCBDataF", "BW (x) CB", massF, bwF, cbDataF);
 RooCMSShape bgF("bgF", "bgF", massF, alphacmsF, betacmsF, gammacmsF, peakcmsF);

 RooAddPdf dataModelP("dataModelP","dataModelP",RooArgList(BWxCBDataP,bgP), RooArgList(nsigP,nbkgP)) ;
 RooAddPdf dataModelF("dataModelF","dataModelF",RooArgList(BWxCBDataF,bgF), RooArgList(nsigF,nbkgF)) ;

 RooFitResult* dataFitP = dataModelP.fitTo(dsDataP,Save(),  SumW2Error(kTRUE)) ;
 RooFitResult* dataFitF = dataModelF.fitTo(dsDataF,Save(),  SumW2Error(kTRUE)) ;

 dataModelP.plotOn(frameP,LineColor(kRed)) ;
 dataModelP.plotOn(frameP,Components(BWxCBDataP) ,LineColor(kBlue)) ;
 dataModelP.plotOn(frameP,Components(bgP) ,LineColor(kCyan)) ;

 dataModelF.plotOn(frameF,LineColor(kRed)) ;
 dataModelF.plotOn(frameF,Components(BWxCBDataF) ,LineColor(kBlue)) ;
 dataModelF.plotOn(frameF,Components(bgF) ,LineColor(kCyan)) ;

 Double_t chi2P = frameP->chiSquare("dataModelP", "dsDataP", 9);
 Double_t chi2F = frameP->chiSquare("dataModelF", "dsDataF", 9);
 std::cout<<"Chi Square=:"<<chi2P<<std::endl;
 std::cout<<"Chi Square=:"<<chi2F<<std::endl;

 massP.setRange("signal",81,101) ;
 RooAbsReal* fracSigP = BWxCBDataP.createIntegral(massP, massP, "signal") ;
 Double_t nsig_fracP  = nsigP.getVal() * fracSigP->getVal();
 RooAbsReal* fracBkgP = bgP.createIntegral(massP, massP, "background") ;
 Double_t nbkg_fracP  = nbkgP.getVal() * fracBkgP->getVal();

 cout<< "Pass Signal = " << nsig_fracP << endl;
 cout<< "Pass BKGD = " << nbkg_fracP << endl;

 massF.setRange("signal",81,101) ;
 RooAbsReal* fracSigF = BWxCBDataF.createIntegral(massF, massF, "signal") ;
 Double_t nsig_fracF  = nsigF.getVal() * fracSigF->getVal();
 RooAbsReal* fracBkgF = bgF.createIntegral(massF, massF, "background") ;
 Double_t nbkg_fracF  = nbkgF.getVal() * fracBkgF->getVal();

 cout<< "Fail Signal = " << nsig_fracF << endl;
 cout<< "Fail BKGD = " << nbkg_fracF << endl;

 TCanvas *c1 = new TCanvas("c1","c1",2000,900);
 c1->Divide(2,1,0.001,0.001);
 gStyle->SetOptStat(0);

 float yTitleOffset = 1.0;
 float yTitleSize = 0.05;
 float yLabelSize = 0.05;
 float xTitleOffset = 1.;
 float xTitleSize = 0.05;
 float xLabelSize = 0.05;

 c1->cd(1);
 gPad->SetTicky(1);
 gPad->SetTickx(1);
 gPad->SetBottomMargin(0.12);
 gPad->SetRightMargin(0.05);
 frameP->GetYaxis()->SetTitleOffset(yTitleOffset);
 frameP->GetYaxis()->SetTitleSize(yTitleSize);
 frameP->GetYaxis()->SetLabelSize(yLabelSize);
 frameP->GetYaxis()->SetDecimals(2);
 frameP->GetXaxis()->SetTitleOffset(xTitleOffset);
 frameP->GetXaxis()->SetLabelSize(xLabelSize);
 frameP->GetXaxis()->SetTitleSize(xTitleSize);
 frameP->SetMarkerStyle(20);
 //passHist->SetMinimum(0.5);
 //passHist->SetMaximum(1.01);
 frameP->SetMarkerColor(kRed);
 frameP->SetLineColor(kRed);
 frameP->Draw("pe");

 TLegend* leg = new TLegend(0.25, 0.20, 0.5, 0.4,"","brNDC");
 leg->SetTextSize(0.05);
 leg->AddEntry(passHist, "PM filter, +EE", "pl");
 leg->SetBorderSize(0);
 //leg->Draw();

 c1->cd(2);
 gPad->SetTicky(1);
 gPad->SetTickx(1);
 gPad->SetBottomMargin(0.12);
 gPad->SetRightMargin(0.05);
 frameF->GetYaxis()->SetTitleOffset(yTitleOffset);
 frameF->GetYaxis()->SetTitleSize(yTitleSize);
 frameF->GetYaxis()->SetLabelSize(yLabelSize);
 frameF->GetYaxis()->SetDecimals(2);
 frameF->GetXaxis()->SetTitleOffset(xTitleOffset);
 frameF->GetXaxis()->SetLabelSize(xLabelSize);
 frameF->GetXaxis()->SetTitleSize(xTitleSize);
 frameF->SetMarkerStyle(20);
 //passHist->SetMinimum(0.5);
 //passHist->SetMaximum(1.01);
 frameF->SetMarkerColor(kRed);
 frameF->SetLineColor(kRed);
 frameF->Draw("pe");

 c1->SaveAs("pass_fail.png");
}
