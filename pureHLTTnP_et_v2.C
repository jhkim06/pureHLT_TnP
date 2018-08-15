#include <TROOT.h>                    // access to gROOT, entry point to ROOT system
#include <TSystem.h>                  // interface to OS
#include <TFile.h>                    // file handle class
#include <TTree.h>                    // class to access ntuples
#include <TBenchmark.h>               // class to track macro running statistics
#include <TH1D.h>                     // histogram class
#include <vector>                     // STL vector class
#include <iomanip>                    // functions to format standard I/O
#include <string>                     // C++ string class
#include <sstream>                    // class for parsing strings
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
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "TText.h"
#include "RooFitResult.h"
#include "RooRealBinding.h"
#include "RooBrentRootFinder.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TPad.h"

/// include pdfs
#include "RooCBExGaussShape.h"
#include "RooCMSShape.h"

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
 //TString filename1 = "/Volumes/Samsung_T3/2018_HLT_Performances/PixelMatchStudy/CMSSW_10_1_5/ntuple/ntuple_CMSSW_10_1_5_BadPixModuleFix.root";
 //TString filename1 = "/Volumes/Samsung_T3/2018_HLT_Performances/PixelMatchStudy/CMSSW_10_1_2_patch2/ntuple/ntuple_CMSSW_10_1_2_patch.root";
 TString filename1 = "/Volumes/Samsung_T3/2018_HLT_Performances/PixelMatchStudy/TnPsamples/2018Bv1/2018Bv1.root";
 infile1 = new TFile(filename1);

 TTree *tree=0;
 tree = (TTree*)infile1->Get("hltTPTree");

 //double etBins[] = {25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0, 80.0, 100.0, 150.0, 200.0, 250.0};  
 double etBins[] = {35.0, 40.0};  
 const int NetBins = sizeof(etBins)/sizeof(double); 

 vector<TString> etbinNames;
 for( int i = 0; i < NetBins-1; i++){
      TString lowerBound, upperBound;
      lowerBound.Form("%f", etBins[i]);
      upperBound.Form("%f", etBins[i+1]);
      etbinNames.push_back("et"+lowerBound+"et"+upperBound);
      cout << "check bin names: " << "et"+lowerBound+"et"+upperBound <<  endl;
 }

 TH1D* htotalSig = new TH1D("htotalSig","htotalSig", NetBins-1, etBins);
 TH1D* hPtotalSig = new TH1D("hPtotalSig","hPtotalSig", NetBins-1, etBins);

 for( int ibin = 0; ibin < NetBins-1; ibin++){

    //if(ibin != 8) continue;
    TString lowerBound, upperBound;
    lowerBound.Form("%f", etBins[ibin]);
    upperBound.Form("%f", etBins[ibin+1]);

    // for pixel matching filter
    TString var = "mass";
    //TString sampleCuts = " tagHLTRegion==0 && mass>60 && mass<120 && probeHLT.et > " + lowerBound + "&& probeHLT.et < " + upperBound + " && probeHLT.eta > 0 && probeHLTRegion==1 && probeHLT.et > 30 && (evtTrigs[0]&0x4000000)!=0 && (tagTrigs[2]&0x400000)!=0 &&  (probeTrigs[2]&0x8000)!=0";
    TString sampleCuts = " tagHLTRegion==0 && mass>50 && mass<130 && probeHLT.et > " + lowerBound + "&& probeHLT.et < " + upperBound + " && (tagTrigs[2]&0x800)!=0 &&  (probeTrigs[2]&0x10)!=0";

    cout << sampleCuts << endl;

    //TString passingProbe = "(probeTrigs[2]&0x10000)!=0"; // pixel matching filter bit
    //TString failingProbe = "(probeTrigs[2]&0x10000)==0"; // pixel matching filter bit

    TString passingProbe = "(probeTrigs[2]&0x20)!=0"; // pixel matching filter bit
    TString failingProbe = "(probeTrigs[2]&0x20)==0"; // pixel matching filter bit

    TH1* passHist;
    TH1* failHist;

    bool _useMinos = false;
    double _xFitMin = 60.;
    double _xFitMax = 120.;
      
    //passHist = MakePassHist(tree, "Run 315322, Eff. w.r.t. HCalIso filter", "Mass", 60, 60., 120., var, sampleCuts, passingProbe); 
    //failHist = MakePassHist(tree, "Run 315322, Eff. w.r.t. HCalIso filter", "Mass", 60, 60., 120., var, sampleCuts, failingProbe); 

    passHist = MakePassHist(tree, "2018 Run Bv1, Eff. w.r.t. HCalIso filter", "Mass", 80, 50., 130., var, sampleCuts, passingProbe); 
    failHist = MakePassHist(tree, "2018 Run Bv1, Eff. w.r.t. HCalIso filter", "Mass", 80, 50., 130., var, sampleCuts, failingProbe); 

    // Set up Extended binned Maximum likelihood fit

    // total number of probes in the histogram
    double nTotP = passHist->Integral();
    double nTotF = failHist->Integral();

    RooWorkspace *_work;
    _work = new RooWorkspace("w") ;

    _work->factory("massP[50,130]");
    _work->factory("massF[50,130]");

    RooDataHist dsDataP("dsDataP","dsDataP",*_work->var("massP"),passHist);
     _work->import(dsDataP);
    RooDataHist dsDataF("dsDataF","dsDataF",*_work->var("massF"),failHist);
     _work->import(dsDataF);

    RooPlot* frameP = _work->var("massP")->frame(60,120);
    frameP->SetTitle("passing probe");
    RooPlot* frameF = _work->var("massF")->frame(60,120);
    frameF->SetTitle("failing probe");

     // set initial fitting variables for passing pdf
    _work->factory("meanP[0.2,-5.0,5.0]");
    _work->factory("sigmaP[1,0.2,6.0]");
    _work->factory("alphaP[0.3,0.0,3.5]");
    _work->factory("nP[3,-5,10]");
    _work->factory("sigmaP_2[1.5,0.,6.0]");
    _work->factory("sosP[1,0.,5.0]");
    _work->factory("tailLeft[-1]");

    if( etBins[ibin+1] <= 35 ) _work->factory("tailLeft[1]");

    _work->factory("acmsP[60.,50.,150.]");
    _work->factory("betaP[0.04,0.01,0.06]");
    _work->factory("gammaP[0.1, 0.005, 1]");
    _work->factory("peakP[90.0]");

    _work->factory("mean1P[85., 80., 100.]");
    _work->factory("sigma1P[2.0., 1.0, 3.]");

    _work->factory("mean2P[91.2,85.,95.]");
    _work->factory("sigma2P[0.1,0.,0.5]");

    // set initial fitting variables for failing pdf
    _work->factory("meanF[-0.0,-5.0,5.0]");
    _work->factory("sigmaF[2,0.2,15.0]");
    _work->factory("alphaF[2.0,1.2,3.5]");
    _work->factory("nF[3,-5,5]");
    _work->factory("sigmaF_2[2.0,0.,6.0]");
    _work->factory("sosF[1,0.,5.0]");

    _work->factory("acmsF[60.,50.,75.]");
    _work->factory("betaF[0.04,0.01,0.06]");
    _work->factory("gammaF[0.1, 0.005, 1]");
    _work->factory("peakF[90.0]");

    _work->factory("mean1F[91.2]");
    _work->factory("sigma1F[2.5]");

    _work->factory("mean2F[91.2,85.,95.]");
    _work->factory("sigma2F[0.1,0.,0.5]");

    // pdf for passing probes
    _work->factory("RooCBExGaussShape::sigResPass(massP,meanP,expr('sqrt(sigmaP*sigmaP+sosP*sosP)',{sigmaP,sosP}),alphaP,nP, expr('sqrt(sigmaP_2*sigmaP_2+sosP*sosP)',{sigmaP_2,sosP}),tailLeft)");
    _work->factory("BreitWigner::bwP(massP,mean1P,sigma1P)");
    _work->factory("Gaussian::g1P(massP,mean1P,sigma1P)");
    _work->factory("Gaussian::g2P(massP,mean2P,sigma2P)");
    _work->factory(TString::Format("ng1P[%f,0.5,%f]",nTotP*0.5,nTotP));
    _work->factory(TString::Format("ng2P[%f,0.5,%f]",nTotP*0.1,nTotP));
    _work->factory("SUM::gaussP(ng1P*g1P,ng2P*g2P)");
    _work->factory("FCONV::sigPass(massP, bwP, sigResPass)");
    _work->factory("RooCMSShape::bkgPass(massP, acmsP, betaP, gammaP, peakP)");
    _work->factory(TString::Format("nSigP[%f,0.5,%f]",nTotP*0.9,nTotP*1.5));
    _work->factory(TString::Format("nBkgP[%f,0.5,%f]",nTotP*0.1,nTotP*1.5));
    //_work->factory("RooExtendPdf::esigPass(sigPass, nSigP)");
    //_work->factory("RooExtendPdf::ebkgPass(bkgPass, nBkgP)");
    _work->var("massP")->setRange("window",81.,101.);
    RooExtendPdf esigPass("esigPass", "esigPass", *_work->pdf("sigPass"), *_work->var("nSigP"), "window");
    RooExtendPdf ebkgPass("ebkgPass", "ebkgPass", *_work->pdf("bkgPass"), *_work->var("nBkgP"), "window");
     //_work->import(esigPass);
     //_work->import(ebkgPass);
    //_work->factory("SUM::pdfPass(nSigP*esigPass,nBkgP*ebkgPass)");
     RooAddPdf _pdfPass("pdfPass","pdfPass",RooArgList(esigPass,ebkgPass)) ;
     //_work->import(_pdfPass);
    //_work->factory("SUM::pdfPass(esigPass,ebkgPass)");

    // pdf for failing probes
    _work->factory("RooCBExGaussShape::sigResFail(massF,meanF,expr('sqrt(sigmaF*sigmaF+sosF*sosF)',{sigmaF,sosF}),alphaF,nF, expr('sqrt(sigmaF_2*sigmaF_2+sosF*sosF)',{sigmaF_2,sosF}),tailLeft)");
    _work->factory("BreitWigner::bwF(massF,mean1F,sigma1F)");
    _work->factory("Gaussian::g1F(massF,mean1F,sigma1F)");
    _work->factory("Gaussian::g2F(massF,mean2F,sigma2F)");
    _work->factory(TString::Format("ng1F[%f,0.5,%f]",nTotF*0.5,nTotF));
    _work->factory(TString::Format("ng2F[%f,0.5,%f]",nTotF*0.1,nTotF));
    _work->factory("SUM::gaussF(ng1F*g1F,ng2F*g2F)");
    _work->factory("FCONV::sigFail(massF, bwF, sigResFail)");
    _work->factory("RooCMSShape::bkgFail(massF, acmsF, betaF, gammaF, peakF)");
    _work->factory(TString::Format("nSigF[%f,0.5,%f]",nTotF*0.9,nTotF*1.5));
    _work->factory(TString::Format("nBkgF[%f,0.5,%f]",nTotF*0.1,nTotF*1.5));
    _work->var("massF")->setRange("window",81.,101.);
    RooExtendPdf esigFail("esigFail", "esigFail", *_work->pdf("sigFail"), *_work->var("nSigF"), "window");
    RooExtendPdf ebkgFail("ebkgFail", "ebkgFail", *_work->pdf("bkgFail"), *_work->var("nBkgF"), "window");
    //_work->factory("SUM::pdfFail(nSigF*sigFail,nBkgF*bkgFail)");
    RooAddPdf _pdfFail("pdfFail","pdfFail",RooArgList(esigFail,ebkgFail)) ;

    _work->var("massP")->setRange("fitMassRange", 50., 130.);
    _work->var("massF")->setRange("fitMassRange", 50., 130.);

    //RooFitResult* resPass = pdfPass->fitTo(*_work->data("dsDataP"), Minos(_useMinos), Range("fitMassRange"));
    RooFitResult* resPass = _pdfPass.fitTo(*_work->data("dsDataP"), Save(),Extended(), SumCoefRange("fitMassRange"));
    //RooFitResult* resFail = pdfFail->fitTo(*_work->data("dsDataF"), Minos(_useMinos), Range("fitMassRange"));
    RooFitResult* resFail = _pdfFail.fitTo(*_work->data("dsDataF"),Save(),Extended(), SumCoefRange("fitMassRange"));

    RooPlot *pPass = _work->var("massP")->frame(60.,120.);
    pPass->SetTitle("passing probe");

    dsDataP.plotOn(pPass, Name("dsDataP"), MarkerColor(kBlack), MarkerStyle(21), MarkerSize(1.2)) ;
    _pdfPass.plotOn(pPass,Name("pdfPass"),LineColor(kRed)) ;
    _pdfPass.plotOn(pPass,Components(esigPass) ,LineColor(kBlack), LineStyle(kDashed)) ;
    _pdfPass.plotOn(pPass,Components(ebkgPass) ,LineColor(kBlue), LineStyle(kDashed)) ;
    
    //_work->data("dsDataP")->plotOn( pPass );
    //_work->pdf("pdfPass")->plotOn( pPass, LineColor(kRed) );
    //_work->pdf("pdfPass")->plotOn( pPass, Components("esigPass"),LineColor(kBlack),LineStyle(kDashed));
    //_work->pdf("pdfPass")->plotOn( pPass, Components("ebkgPass"),LineColor(kBlue),LineStyle(kDashed));

    RooPlot *pFail = _work->var("massF")->frame(60.,120.);
    pFail->SetTitle("Failing probe");

    dsDataF.plotOn(pFail, Name("dsDataF"), MarkerColor(kBlack), MarkerStyle(21), MarkerSize(1.2)) ;
    _pdfFail.plotOn(pFail,Name("pdfFail"),LineColor(kRed)) ;
    _pdfFail.plotOn(pFail,Components(esigFail) ,LineColor(kBlack), LineStyle(kDashed)) ;
    _pdfFail.plotOn(pFail,Components(ebkgFail) ,LineColor(kBlue), LineStyle(kDashed)) ;

    //_work->data("dsDataF")->plotOn( pFail );
    //_work->pdf("pdfFail")->plotOn( pFail, LineColor(kRed) );
    //_work->pdf("pdfFail")->plotOn( pFail, Components("sigFail"),LineColor(kBlack),LineStyle(kDashed));
    //_work->pdf("pdfFail")->plotOn( pFail, Components("bkgFail"),LineColor(kBlue),LineStyle(kDashed));
    
    Double_t chi2P = pPass->chiSquare("pdfPass", "dsDataP", 13);
    Double_t chi2F = pFail->chiSquare("pdfFail", "dsDataF", 13);
    std::cout<<"Chi Square=:"<<chi2P<<std::endl;
    std::cout<<"Chi Square=:"<<chi2F<<std::endl;

    cout << "Passing nsig: " << _work->var("nSigP")->getVal() << " error: " << _work->var("nSigP")->getError()<< endl;
    cout << "Failing nsig: " << _work->var("nSigF")->getVal() << " error: " << _work->var("nSigF")->getError()<< endl;

    //dataFitP->Print();

    // fill total signal and passing signal histograms 
    htotalSig->SetBinContent(ibin+1, _work->var("nSigP")->getVal() + _work->var("nSigF")->getVal());
    htotalSig->SetBinError(ibin+1, _work->var("nSigP")->getError() + _work->var("nSigF")->getError());
    hPtotalSig->SetBinContent(ibin+1, _work->var("nSigP")->getVal());
    hPtotalSig->SetBinError(ibin+1, _work->var("nSigP")->getError());

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
    pPass->SetMaximum(1.2*passHist->GetMaximum());
    pPass->GetYaxis()->SetTitleOffset(yTitleOffset);
    pPass->GetYaxis()->SetTitleSize(yTitleSize);
    pPass->GetYaxis()->SetLabelSize(yLabelSize);
    pPass->GetYaxis()->SetDecimals(2);
    pPass->GetXaxis()->SetTitleOffset(xTitleOffset);
    pPass->GetXaxis()->SetLabelSize(xLabelSize);
    pPass->GetXaxis()->SetTitleSize(xTitleSize);
    pPass->SetMarkerStyle(20);
    //passHist->SetMinimum(0.5);
    //passHist->SetMaximum(1.01);
    pPass->SetMarkerColor(kRed);
    pPass->SetLineColor(kRed);
    pPass->Draw("pe");

    TLatex passFitInfo;
    passFitInfo.SetNDC();                                                                                  
    passFitInfo.SetTextFont(42);                                                                           
    passFitInfo.SetTextSize(0.03);  
    
    TString chiPass;                                                                                     
    chiPass.Form("chi^{2}/NDF = %.3f", chi2P);  
    passFitInfo.DrawLatex(0.65,0.85, chiPass);

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
    pFail->SetMaximum(1.2*failHist->GetMaximum());
    pFail->GetYaxis()->SetTitleOffset(yTitleOffset);
    pFail->GetYaxis()->SetTitleSize(yTitleSize);
    pFail->GetYaxis()->SetLabelSize(yLabelSize);
    pFail->GetYaxis()->SetDecimals(2);
    pFail->GetXaxis()->SetTitleOffset(xTitleOffset);
    pFail->GetXaxis()->SetLabelSize(xLabelSize);
    pFail->GetXaxis()->SetTitleSize(xTitleSize);
    pFail->SetMarkerStyle(20);
    //passHist->SetMinimum(0.5);
    //passHist->SetMaximum(1.01);
    pFail->SetMarkerColor(kRed);
    pFail->SetLineColor(kRed);
    pFail->Draw("pe");

    TLatex failFitInfo;
    failFitInfo.SetNDC();
    failFitInfo.SetTextFont(42);
    failFitInfo.SetTextSize(0.03);

    TString chiFail;
    chiFail.Form("chi^{2}/NDF = %.3f", chi2F);
    failFitInfo.DrawLatex(0.65,0.85, chiFail);

    c1->SaveAs("pass_fail_" + etbinNames.at(ibin) + ".png");
    delete passHist;
    delete failHist;
    delete c1;
 }

 //TGraphAsymmErrors* eff = new TGraphAsymmErrors(hPtotalSig,htotalSig,"B");
 TGraphAsymmErrors* eff = new TGraphAsymmErrors(hPtotalSig,htotalSig);

 TCanvas *c1 = new TCanvas("c1","c1",1000,800);
 gStyle->SetOptStat(0);

 float yTitleOffset = 1.0;
 float yTitleSize = 0.05;
 float yLabelSize = 0.05;
 float xTitleOffset = 1.;
 float xTitleSize = 0.05;
 float xLabelSize = 0.05;

 c1->cd();
 gPad->SetGridy();
 gPad->SetTicky(1);
 gPad->SetTickx(1);
 gPad->SetBottomMargin(0.12);
 gPad->SetRightMargin(0.05);
 //eff->SetTitle("Run 315322, Eff. w.r.t. HCalIso filter");
 eff->SetTitle("Pixel match filter Eff. w.r.t. HCal Iso filter");
 eff->GetYaxis()->SetTitleOffset(yTitleOffset);
 eff->GetYaxis()->SetTitleSize(yTitleSize);
 eff->GetYaxis()->SetLabelSize(yLabelSize);
 eff->GetYaxis()->SetDecimals(2);
 eff->GetXaxis()->SetTitleOffset(xTitleOffset);
 eff->GetXaxis()->SetLabelSize(xLabelSize);
 eff->GetXaxis()->SetTitleSize(xTitleSize);
 eff->GetXaxis()->SetRangeUser(20.,250.);
 eff->GetXaxis()->SetTitle("et_{HLT}");
 eff->GetYaxis()->SetTitle("Efficiency");
 eff->SetMarkerStyle(20);
 eff->SetMinimum(0.5);
 eff->SetMaximum(1.01);
 eff->SetMarkerColor(kRed);
 eff->SetLineColor(kRed);
 eff->Draw("ape");

 TLatex lumi(210., 1.02, "2.7 fb^{-1} (13 TeV)");
 lumi.SetTextSize(0.03);
 lumi.Draw();

 c1->SaveAs("eff_et.png");

}
