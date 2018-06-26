//root -b -q 'draw_effS.C()'

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
 //TString filename1 = "/Volumes/Samsung_T3/2018_HLT_Performances/PixelMatchStudy/CMSSW_10_1_5/ntuple/ntuple_CMSSW_10_1_5_BadPixModuleFix.root";
 //TString filename1 = "/Volumes/Samsung_T3/2018_HLT_Performances/PixelMatchStudy/CMSSW_10_1_2_patch2/ntuple/ntuple_CMSSW_10_1_2_patch.root";
 TString filename1 = "/Volumes/Samsung_T3/2018_HLT_Performances/PixelMatchStudy/TnPsamples/2018Av3/2018Av3.root";
 infile1 = new TFile(filename1);

 TTree *tree=0;
 tree = (TTree*)infile1->Get("hltTPTree");

 double phiBins[] = {-3.15, -2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4, 3.15};  
 const int NphiBins = sizeof(phiBins)/sizeof(double); 

 vector<TString> etbinNames;
 for( int i = 0; i < NphiBins-1; i++){
      TString lowerBound, upperBound;
      lowerBound.Form("%f", phiBins[i]);
      upperBound.Form("%f", phiBins[i+1]);
      etbinNames.push_back("phi"+lowerBound+"phi"+upperBound);
      cout << "check bin names: " << "phi"+lowerBound+"phi"+upperBound <<  endl;
 }

 TH1D* htotalSig = new TH1D("htotalSig","htotalSig", NphiBins-1, phiBins);
 TH1D* hPtotalSig = new TH1D("hPtotalSig","hPtotalSig", NphiBins-1, phiBins);

 // variables in the input tnp tree
 struct HLTEgammaStruct {
     float et,nrgy,eta,phi,hadem,sigmaIEtaIEta,dEtaIn,dPhiIn,nrMissHits,nrClus,seedClusEFrac,pmDPhi1,pmDPhi2,pmDPhi3,pmDPhi4,pmDPhi1Info,pmDPhi2Info,pmDPhi3Info,pmDPhi4Info,pmDRZ1,pmDRZ2,pmDRZ3,pmDRZ4,pmDRZ1Info,pmDRZ2Info,pmDRZ3Info,pmDRZ4Info,phiWidth,etaWidth,s2,dPhi1BestS2,dPhi2BestS2,dPhi3BestS2,dRZ2BestS2,dRZ3BestS2,ecalIso,hcalIso,trkIso,trkChi2,invEOInvP,trkIso2016;
 };
 HLTEgammaStruct probeHLT;
 HLTEgammaStruct tagHLT;
 int tagTrigs[4];
 int probeTrigs[4];
 int tagHLTRegion;
 int probeHLTRegion;
 float mass;
 tree->SetBranchAddress("tagHLTRegion",&tagHLTRegion);
 tree->SetBranchAddress("probeHLTRegion",&probeHLTRegion);
 tree->SetBranchAddress("probeHLT",&probeHLT);
 tree->SetBranchAddress("tagHLT",&tagHLT);
 tree->SetBranchAddress("tagTrigs",&tagTrigs);
 tree->SetBranchAddress("probeTrigs",&probeTrigs);
 tree->SetBranchAddress("mass",&mass);

 //Long64_t nentries = tree->GetEntries();
 //  for (Long64_t i=0;i<100;i++) {
 //    tree->GetEntry(i);
 //    cout << "tagHLTRegion: " << tagHLTRegion << " eta: " << tagHLT.eta << " phi: " << tagHLT.phi << endl;   
 //}

 for( int ibin = 0; ibin < NphiBins-1; ibin++){

    if(ibin != 1) continue;
    TString lowerBound, upperBound;
    lowerBound.Form("%f", phiBins[ibin]);
    upperBound.Form("%f", phiBins[ibin+1]);

    // for pixel matching filter
    TString var = "mass";
    //TString sampleCuts = " tagHLTRegion==0 && mass>60 && mass<120 && probeHLT.phi > " + lowerBound + "&& probeHLT.phi < " + upperBound + " && probeHLT.eta > 0 && probeHLTRegion==1 && probeHLT.et > 30 && (evtTrigs[0]&0x4000000)!=0 && (tagTrigs[2]&0x400000)!=0 &&  (probeTrigs[2]&0x8000)!=0";
    TString sampleCuts = " tagHLTRegion==0 && mass>60 && mass<120 && probeHLT.phi > " + lowerBound + "&& probeHLT.phi < " + upperBound + " && probeHLT.et > 50  && (tagTrigs[2]&0x800)!=0 &&  (probeTrigs[2]&0x10)!=0";
    cout << sampleCuts << endl;
    //TString passingProbe = "(probeTrigs[2]&0x10000)!=0"; // pixel matching filter bit
    //TString failingProbe = "(probeTrigs[2]&0x10000)==0"; // pixel matching filter bit

    TString passingProbe = "(probeTrigs[2]&0x20)!=0"; // pixel matching filter bit
    TString failingProbe = "(probeTrigs[2]&0x20)==0"; // pixel matching filter bit

    TH1* passHist;
    TH1* failHist;
      
    //passHist = MakePassHist(tree, "Run 315322, Eff. w.r.t. HCalIso filter", "Mass", 60, 60., 120., var, sampleCuts, passingProbe); 
    //failHist = MakePassHist(tree, "Run 315322, Eff. w.r.t. HCalIso filter", "Mass", 60, 60., 120., var, sampleCuts, failingProbe); 

    passHist = MakePassHist(tree, "2018 RunBv1, Eff. w.r.t. HCalIso filter", "Mass", 60, 60., 120., var, sampleCuts, passingProbe); 
    failHist = MakePassHist(tree, "2018 RunBv1, Eff. w.r.t. HCalIso filter", "Mass", 60, 60., 120., var, sampleCuts, failingProbe); 


    TFile *temp_file = new TFile("temp.root","RECREATE");
    TTree* tnp_tree = new TTree("tnpTree","tnpTree");

    temp_file->cd();
    float passBranch; tnp_tree->Branch("massP", &passBranch);
    float failBranch; tnp_tree->Branch("massF", &failBranch);

    // fill tree
    Long64_t nentries = tree->GetEntries();
      for (Long64_t i=0;i<nentries;i++) {
        tree->GetEntry(i);

        passBranch = -999.;
        failBranch = -999.;

        // probes
        if( mass > 60 && mass < 120 && tagHLTRegion==0 && probeHLT.phi > phiBins[ibin] && probeHLT.phi < phiBins[ibin+1] && probeHLT.et > 50 && (tagTrigs[2]&0x800)!=0 &&  (probeTrigs[2]&0x10)!=0){
           // passing probe
           if((probeTrigs[2]&0x20)!=0){ passBranch = mass; tnp_tree->Fill();} 
           // failing probes
           if((probeTrigs[2]&0x20)==0){ failBranch = mass; tnp_tree->Fill();}

        }
        //tnp_tree->Fill(); 
    }

    infile1->cd();
    // Set up Extended binned Maximum likelihood fit

    // total number of probes in the histogram
    double nToP = passHist->Integral();
    double nToF = failHist->Integral();

    RooRealVar massP("massP","m_{ee} [GeV]", 60., 120.);
    RooRealVar massF("massF","m_{ee} [GeV]", 60., 120.);

    //RooDataHist dsDataP("dsDataP","dsDataP",massP,passHist);
    //RooDataHist dsDataF("dsDataF","dsDataF",massF,failHist);

    RooDataSet dsDataP("dsDataP","dsDataP",tnp_tree,massP);
    RooDataSet dsDataF("dsDataF","dsDataF",tnp_tree,massF);

    RooPlot* frameP = massP.frame(Name("passing"), Title("passing")) ;
    RooPlot* frameF = massF.frame(Name("failing"), Title("failing")) ;
    //frameP->SetMaximum(1.2*passHist->GetMaximum());
    //frameF->SetMaximum(1.2*failHist->GetMaximum());

    // plot histograms on frames
    dsDataP.plotOn(frameP, Name("dsDataP"), MarkerColor(kBlack), MarkerStyle(21), MarkerSize(1.2)) ;
    dsDataF.plotOn(frameF, Name("dsDataF"), MarkerColor(kBlack), MarkerStyle(21), MarkerSize(1.2)) ;

    // Crystal ball function
    RooRealVar cbMeanDataP("cbMeanDataP","cbMean", -2., -4., 4.);
    RooRealVar cbSigmaDataP("cbSigmaDataP","cbSigma", 2., 0.5, 4.);
    RooRealVar cbAlphaDataP("cbAlphaDataP","cbAlpha", 20., 0.1, 50);
    RooRealVar cbNDataP("cbNDataP","cbN", 30., 0.2, 50);

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
    RooRealVar alphacmsP("alphacmsP","alphacmsP", 60., 50., 100.);
    RooRealVar betacmsP("betacmsP","betacmsP", 0.05,0.01,0.4);
    RooRealVar gammacmsP("gammacmsP","gammacmsP", 0.1, 0, 1);
    RooRealVar peakcmsP("peakcmsP","peakcmsP", 90.0);

    // values from egm_tnp_analysis
    //RooRealVar alphacmsP("alphacmsP","alphacmsP", 90., 50., 100.);
    //RooRealVar betacmsP("betacmsP","betacmsP", 0.05,0.01,0.5);
    //RooRealVar gammacmsP("gammacmsP","gammacmsP", 0.1, 0., 1.);
    //RooRealVar peakcmsP("peakcmsP","peakcmsP", 90.0, 85., 95.);

    RooRealVar alphacmsF("alphacmsF","alphacmsF", 60., 50., 80.);
    RooRealVar betacmsF("betacmsF","betacmsF", 0.05,0.01,0.2);
    RooRealVar gammacmsF("gammacmsF","gammacmsF", 0.1, 0, 1);
    RooRealVar peakcmsF("peakcmsF","peakcmsF", 90.0);

    //RooRealVar alphacmsF("alphacmsF","alphacmsF", 80., 50., 120.);
    //RooRealVar betacmsF("betacmsF","betacmsF", 0.05,0.01,0.2);
    //RooRealVar gammacmsF("gammacmsF","gammacmsF", 0.1, -2., 2.);
    //RooRealVar peakcmsF("peakcmsF","peakcmsF", 90.0);

    RooRealVar nsigP("nsigP","signal events1", nToP*0.99,0.,nToP); // https://github.com/michelif/egm_tnp_analysis/blob/egm_tnp_Moriond18_v3.0/libCpp/histFitter.C#L121
    RooRealVar nbkgP("nbkgP","signal background events1",nToP*0.01,0.,nToP);

    RooRealVar nsigF("nsigF","signal events1",nToF*0.9,0.,nToF);
    RooRealVar nbkgF("nbkgF","signal background events1",nToF*0.1,0.,nToF);

    // create passing pdfs
    RooBreitWigner bwP("bwP","bwP", massP, meanP, sigmaP);
    RooCBShape cbDataP("cbDataP","cb", massP, cbMeanDataP, cbSigmaDataP, cbAlphaDataP, cbNDataP);
    RooFFTConvPdf BWxCBDataP("BWxCBDataP", "BW (x) CB", massP, bwP, cbDataP);
    RooCMSShape bgP("bgP", "bgP", massP, alphacmsP, betacmsP, gammacmsP, peakcmsP);

    // use RooExtendPdf
    massP.setRange("signal",81,101) ;
    RooExtendPdf eBWxCBDataP("esigP", "esigP", BWxCBDataP, nsigP, "signal");
    RooExtendPdf ebgP("ebkgP", "ebkgP", bgP, nbkgP,"signal");

    // create failing pdfs
    RooBreitWigner bwF("bwF","bwF", massF, meanF, sigmaF);
    RooCBShape cbDataF("cbDataF","cb", massF, cbMeanDataF, cbSigmaDataF, cbAlphaDataF, cbNDataF);
    RooFFTConvPdf BWxCBDataF("BWxCBDataF", "BW (x) CB", massF, bwF, cbDataF);
    RooCMSShape bgF("bgF", "bgF", massF, alphacmsF, betacmsF, gammacmsF, peakcmsF);

    // use RooExtendPdf
    massF.setRange("signal",81,101) ;
    RooExtendPdf eBWxCBDataF("esigF", "esigF", BWxCBDataF, nsigF, "signal");
    RooExtendPdf ebgF("ebkgF", "ebkgF", bgF, nbkgF,"signal");

    RooAddPdf dataModelP("dataModelP","dataModelP",RooArgList(eBWxCBDataP,ebgP)) ;
    RooAddPdf dataModelF("dataModelF","dataModelF",RooArgList(eBWxCBDataF,ebgF)) ;

    // fit passing probes
    //RooFitResult* dataFitP = dataModelP.fitTo(dsDataP,Save(),Extended(), SumW2Error(kTRUE)) ;
    RooFitResult* dataFitP = dataModelP.fitTo(dsDataP,Save(),Extended()) ;
    dataModelP.plotOn(frameP,Name("dataModelP"), LineColor(kBlack)) ;
    dataModelP.plotOn(frameP,Components(BWxCBDataP) ,LineColor(kRed), LineStyle(kDashed)) ;
    dataModelP.plotOn(frameP,Components(bgP) ,LineColor(kBlue), LineStyle(kDashed)) ;

    // fit failing probes
    //RooFitResult* dataFitF = dataModelF.fitTo(dsDataF,Save(),Extended(), SumW2Error(kTRUE)) ;
    RooFitResult* dataFitF = dataModelF.fitTo(dsDataF,Save(),Extended()) ;
    dataModelF.plotOn(frameF,Name("dataModelF"),LineColor(kBlack)) ;
    dataModelF.plotOn(frameF,Components(BWxCBDataF) ,LineColor(kRed), LineStyle(kDashed)) ;
    dataModelF.plotOn(frameF,Components(bgF) ,LineColor(kBlue), LineStyle(kDashed)) ;

    Double_t chi2P = frameP->chiSquare("dataModelP", "dsDataP", 10);
    Double_t chi2F = frameF->chiSquare("dataModelF", "dsDataF", 10);
    std::cout<<"Chi Square=:"<<chi2P<<std::endl;
    std::cout<<"Chi Square=:"<<chi2F<<std::endl;

    cout << "Passing nsig: " << nsigP.getVal() << " error: " << nsigP.getError()<< endl;
    cout << "Failing nsig: " << nsigF.getVal() << " error: " << nsigF.getError()<< endl;

    dataFitP->Print();

    // fill total signal and passing signal histograms 
    htotalSig->SetBinContent(ibin+1, nsigP.getVal()+nsigF.getVal());
    htotalSig->SetBinError(ibin+1, nsigP.getError()+nsigF.getError());
    hPtotalSig->SetBinContent(ibin+1, nsigP.getVal());
    hPtotalSig->SetBinError(ibin+1, nsigP.getError());

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
    delete tnp_tree;
    delete temp_file;
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
 eff->SetTitle("2018 RunBv1, Eff. w.r.t. HCalIso filter");
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

 c1->SaveAs("eff_phi.png");

}
