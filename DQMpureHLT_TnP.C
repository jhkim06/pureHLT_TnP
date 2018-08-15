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

/// include pdfs
#include "RooCBExGaussShape.h"
#include "RooCMSShape.h"

#include <iostream>
#include <fstream>

using namespace RooFit ;

TH1* MakeFailHist(TH1* passed, TH1* total){
  TH1F *h_temp = (TH1F*)total->Clone("total");
  TH1F *h_passed = (TH1F*)passed->Clone("passed");

  h_temp->Add(h_passed, -1);

  return h_temp; 
}

void DQMpureHLT_TnP()
{

 TFile* infile1 = NULL;
 infile1 = new TFile("./DQM_V0001_R000320804__HLT__TrigObjTnpSource__All.root");

 // FIX ME: update later to allow navigation to find target 2D histograms
 TH2F* hist_total_2d = dynamic_cast<TH2F*> (( infile1->Get("DQMData/Run 320804/HLT/Run summary/EGM/TrigObjTnP/stdTag_hltEG32L1SingleEGOrEtFilter_phi"))->Clone());
 TH2F* hist_passed_2d = dynamic_cast<TH2F*> (( infile1->Get("DQMData/Run 320804/HLT/Run summary/EGM/TrigObjTnP/stdTag_hltEle32WPTightHcalIsoFilter_phi"))->Clone());

 //TH2F* hist_total_2d = dynamic_cast<TH2F*> (( infile1->Get("DQMData/Run 320804/HLT/Run summary/EGM/TrigObjTnP/stdTag_hltEG32L1SingleEGOrEtFilter_phi"))->Clone());
 //TH2F* hist_passed_2d = dynamic_cast<TH2F*> (( infile1->Get("DQMData/Run 320804/HLT/Run summary/EGM/TrigObjTnP/stdTag_hltEle32WPTightGsfTrackIsoFilter_phi"))->Clone());

 //TH2F* hist_total_2d = dynamic_cast<TH2F*> (( infile1->Get("DQMData/Run 320804/HLT/Run summary/EGM/TrigObjTnP/stdTag_hltEle32WPTightHcalIsoFilter_phi"))->Clone());
 //TH2F* hist_passed_2d = dynamic_cast<TH2F*> (( infile1->Get("DQMData/Run 320804/HLT/Run summary/EGM/TrigObjTnP/stdTag_hltEle32WPTightPixelMatchFilter_phi"))->Clone());

 //cout << "x min: " << hist_total_2d->GetXaxis()->GetXmin() << " max: " << hist_total_2d->GetXaxis()->GetXmax() << endl;
 const double *xaxis = (hist_total_2d->GetXaxis()->GetXbins())->GetArray();

 int nbinx = hist_total_2d->GetNbinsX();

 vector<TString> binNames;
 for( int i = 0; i < nbinx; i++){
      TString lowerBound, upperBound;
      lowerBound.Form("%f", xaxis[i]);
      upperBound.Form("%f", xaxis[i+1]);
      binNames.push_back("phi"+lowerBound+"phi"+upperBound);
      cout << "check bin names: " << "phi"+lowerBound+"phi"+upperBound <<  endl;
 }

 TH1D* htotalSig = new TH1D("htotalSig","htotalSig", nbinx, xaxis);
 TH1D* hPtotalSig = new TH1D("hPtotalSig","hPtotalSig", nbinx, xaxis);

 for(int ibin = 1; ibin < nbinx + 1; ibin++){
 
  TH1* passHist;
  TH1* failHist;
 
  passHist = hist_passed_2d->ProjectionY("passed plot",ibin, ibin,"");
  failHist = MakeFailHist(hist_passed_2d->ProjectionY("passed plot",ibin, ibin,""), hist_total_2d->ProjectionY("total plot",ibin, ibin,""));;
 
  // total number of probes in the histogram
  double nTotP = passHist->Integral();
  double nTotF = failHist->Integral();
 
  RooWorkspace *_work;
  _work = new RooWorkspace("w") ;
 
  _work->factory("massP[60,120]");
  _work->factory("massF[60,120]");
 
  RooDataHist dsDataP("dsDataP","dsDataP",*_work->var("massP"),passHist);
   _work->import(dsDataP);
  RooDataHist dsDataF("dsDataF","dsDataF",*_work->var("massF"),failHist);
   _work->import(dsDataF);
 
  RooPlot* frameP = _work->var("massP")->frame(60,120);
  frameP->SetTitle("passing probe");
  RooPlot* frameF = _work->var("massF")->frame(60,120);
  frameF->SetTitle("failing probe");
 
   // set initial fitting variables for passing pdf
  _work->factory("meanP[-0.0,-5.0,5.0]");
  _work->factory("sigmaP[1,0.2,6.0]");
  _work->factory("alphaP[2.0,1.2,3.5]");
  _work->factory("nP[3,-5,5]");
  _work->factory("sigmaP_2[1.5,0.,6.0]");
  _work->factory("sosP[1,0.,5.0]");
  _work->factory("tailLeft[-1]");
 
  _work->factory("acmsP[60.,50.,150.]");
  _work->factory("betaP[0.04,0.01,0.06]");
  _work->factory("gammaP[0.1, 0.005, 1]");
  _work->factory("peakP[90.0]");
 
  _work->factory("mean1P[91.2,85.,95.]");
  _work->factory("sigma1P[2.5,0.,10.0]");
 
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
 
  _work->factory("mean1F[91.2,85.,95.]");
  _work->factory("sigma1F[2.5,0.,10.0]");
 
  _work->factory("mean2F[91.2,85.,95.]");
  _work->factory("sigma2F[0.1,0.,0.5]");
 
  // pdf for passing probes
  _work->factory("RooCBExGaussShape::sigResPass(massP,meanP,expr('sqrt(sigmaP*sigmaP+sosP*sosP)',{sigmaP,sosP}),alphaP,nP, expr('sqrt(sigmaP_2*sigmaP_2+sosP*sosP)',{sigmaP_2,sosP}),tailLeft)");
  _work->factory("BreitWigner::bwP(massP,mean1P,sigma1P)");
  _work->factory(TString::Format("ng1P[%f,0.5,%f]",nTotP*0.5,nTotP));
  _work->factory(TString::Format("ng2P[%f,0.5,%f]",nTotP*0.1,nTotP));
  _work->factory("FCONV::sigPass(massP, bwP, sigResPass)"); // try later template fit: https://github.com/michelif/egm_tnp_analysis/blob/egm_tnp_Moriond18_v3.0/libPython/fitUtils.py#L117
  _work->factory("RooCMSShape::bkgPass(massP, acmsP, betaP, gammaP, peakP)");
  _work->factory(TString::Format("nSigP[%f,0.5,%f]",nTotP*0.9,nTotP*1.5));
  _work->factory(TString::Format("nBkgP[%f,0.5,%f]",nTotP*0.1,nTotP*1.5));
  _work->var("massP")->setRange("window",81.,101.);
  RooExtendPdf esigPass("esigPass", "esigPass", *_work->pdf("sigPass"), *_work->var("nSigP"), "window");
  RooExtendPdf ebkgPass("ebkgPass", "ebkgPass", *_work->pdf("bkgPass"), *_work->var("nBkgP"), "window");
  RooAddPdf _pdfPass("pdfPass","pdfPass",RooArgList(esigPass,ebkgPass)) ;
 
  // pdf for failing probes
  _work->factory("RooCBExGaussShape::sigResFail(massF,meanF,expr('sqrt(sigmaF*sigmaF+sosF*sosF)',{sigmaF,sosF}),alphaF,nF, expr('sqrt(sigmaF_2*sigmaF_2+sosF*sosF)',{sigmaF_2,sosF}),tailLeft)");
  _work->factory("BreitWigner::bwF(massF,mean1F,sigma1F)");
  _work->factory(TString::Format("ng1F[%f,0.5,%f]",nTotF*0.5,nTotF));
  _work->factory(TString::Format("ng2F[%f,0.5,%f]",nTotF*0.1,nTotF));
  _work->factory("FCONV::sigFail(massF, bwF, sigResFail)");
  _work->factory("RooCMSShape::bkgFail(massF, acmsF, betaF, gammaF, peakF)");
  _work->factory(TString::Format("nSigF[%f,0.5,%f]",nTotF*0.9,nTotF*1.5));
  _work->factory(TString::Format("nBkgF[%f,0.5,%f]",nTotF*0.1,nTotF*1.5));
  _work->var("massF")->setRange("window",81.,101.);
  RooExtendPdf esigFail("esigFail", "esigFail", *_work->pdf("sigFail"), *_work->var("nSigF"), "window");
  RooExtendPdf ebkgFail("ebkgFail", "ebkgFail", *_work->pdf("bkgFail"), *_work->var("nBkgF"), "window");
  RooAddPdf _pdfFail("pdfFail","pdfFail",RooArgList(esigFail,ebkgFail)) ;
 
  _work->var("massP")->setRange("fitMassRange", 60., 120.);
  _work->var("massF")->setRange("fitMassRange", 60., 120.);
 
 
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
 
  RooPlot *pFail = _work->var("massF")->frame(60.,120.);
  pFail->SetTitle("Failing probe");
 
  dsDataF.plotOn(pFail, Name("dsDataF"), MarkerColor(kBlack), MarkerStyle(21), MarkerSize(1.2)) ;
  _pdfFail.plotOn(pFail,Name("pdfFail"),LineColor(kRed)) ;
  _pdfFail.plotOn(pFail,Components(esigFail) ,LineColor(kBlack), LineStyle(kDashed)) ;
  _pdfFail.plotOn(pFail,Components(ebkgFail) ,LineColor(kBlue), LineStyle(kDashed)) ;
 
  Double_t chi2P = pPass->chiSquare("pdfPass", "dsDataP", 13);
  Double_t chi2F = pFail->chiSquare("pdfFail", "dsDataF", 13);
  std::cout<<"Chi Square=:"<<chi2P<<std::endl;
  std::cout<<"Chi Square=:"<<chi2F<<std::endl;
 
  cout << "Passing nsig: " << _work->var("nSigP")->getVal() << " error: " << _work->var("nSigP")->getError()<< endl;
  cout << "Failing nsig: " << _work->var("nSigF")->getVal() << " error: " << _work->var("nSigF")->getError()<< endl;
 
  // fill total signal and passing signal histograms 
  htotalSig->SetBinContent(ibin, _work->var("nSigP")->getVal() + _work->var("nSigF")->getVal());
  htotalSig->SetBinError(ibin, _work->var("nSigP")->getError() + _work->var("nSigF")->getError());
  hPtotalSig->SetBinContent(ibin, _work->var("nSigP")->getVal());
  hPtotalSig->SetBinError(ibin, _work->var("nSigP")->getError());
 
 
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
  gPad->SetLogy();
  gPad->SetTicky(1);
  gPad->SetTickx(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.05);
  pPass->SetMaximum(1.2*pPass->GetMaximum());
  pPass->SetMinimum(1e-2);
  pPass->GetYaxis()->SetTitleOffset(yTitleOffset);
  pPass->GetYaxis()->SetTitleSize(yTitleSize);
  pPass->GetYaxis()->SetLabelSize(yLabelSize);
  pPass->GetYaxis()->SetDecimals(2);
  pPass->GetXaxis()->SetTitleOffset(xTitleOffset);
  pPass->GetXaxis()->SetLabelSize(xLabelSize);
  pPass->GetXaxis()->SetTitleSize(xTitleSize);
  pPass->SetMarkerStyle(20);
  //pPass->SetMinimum(0.5);
  //pPass->SetMaximum(1.01);
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
  leg->AddEntry(pPass, "PM filter, +EE", "pl");
  leg->SetBorderSize(0);
  //leg->Draw();
  c1->cd(2);
  gPad->SetLogy();
  gPad->SetTicky(1);
  gPad->SetTickx(1);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.05);
  pFail->SetMaximum(1.2*pFail->GetMaximum());
  pFail->SetMinimum(1e-2);
  pFail->GetYaxis()->SetTitleOffset(yTitleOffset);
  pFail->GetYaxis()->SetTitleSize(yTitleSize);
  pFail->GetYaxis()->SetLabelSize(yLabelSize);
  pFail->GetYaxis()->SetDecimals(2);
  pFail->GetXaxis()->SetTitleOffset(xTitleOffset);
  pFail->GetXaxis()->SetLabelSize(xLabelSize);
  pFail->GetXaxis()->SetTitleSize(xTitleSize);
  pFail->SetMarkerStyle(20);
  //pPass->SetMinimum(0.5);
  //pPass->SetMaximum(1.01);
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
 
  c1->SaveAs("pass_fail_" + binNames.at(ibin-1) + "_barrel.png");
  //c1->SaveAs("pass_fail.png");
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
 eff->SetTitle("HCal Iso filter Eff. w.r.t. Et filter");
 //eff->SetTitle("Pixel match filter Eff. w.r.t. HCal Iso filter");
 //eff->SetTitle("Track Iso filter Eff. w.r.t. Et Filter");
 eff->GetYaxis()->SetTitleOffset(yTitleOffset);
 eff->GetYaxis()->SetTitleSize(yTitleSize);
 eff->GetYaxis()->SetLabelSize(yLabelSize);
 eff->GetYaxis()->SetDecimals(2);
 eff->GetXaxis()->SetTitleOffset(xTitleOffset);
 eff->GetXaxis()->SetLabelSize(xLabelSize);
 eff->GetXaxis()->SetTitleSize(xTitleSize);
 eff->GetXaxis()->SetRangeUser(-3.15,3.15);
 eff->GetXaxis()->SetTitle("#phi_{HLT}");
 eff->GetYaxis()->SetTitle("Efficiency");
 eff->SetMarkerStyle(20);
 eff->SetMinimum(0.5);
 eff->SetMaximum(1.01);
 eff->SetMarkerColor(kRed);
 eff->SetLineColor(kRed);
 eff->Draw("ape");

 TLatex et_cut(1.7,0.8,"SC E_{T} > 50 GeV");
 et_cut.SetTextSize(0.035);
 //et_cut.Draw();

 //TLatex lumi(1.8, 1.02, "2.7 fb^{-1} (13 TeV)");
 TLatex lumi(1.8, 1.02, "Run 320804");
 lumi.SetTextSize(0.03);
 lumi.Draw();

 c1->SaveAs("eff_phi_barrel.png");
}
