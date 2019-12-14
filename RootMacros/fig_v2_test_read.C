//
//macro to calculate v2(pt) for charged particles and compared them with STAR data: 
//J.~Adams {\it et al.} [STAR Collaboration], Phys. Rev. {\bf C 72} (2005) 014904.  

#include "stdlib.h"
#include "stdio.h"

void fig_v2_test_read()
{
  Double_t v2pi[100], v1pi[100], v2K[100], v2p[100], ev2pi[100], ev2K[100], ev2p[100],pt, ptt[100], ept[100];
  Int_t n2pi[100], n2K[100], n2p[100];

  std::ifstream out1("v2_STAR/v2star0-5"); 
  Double_t ptstar[100],ptlow[100],pthight[100],v2star[100],ev2star[100],epts[100],x1,x2;
  Int_t i=0;

  TH1D *hv2res = new TH1D("hv2res", "hv2res", 100, 0.0, 10.);
   
   while(out1)
   {
   out1>>ptstar[i]>>ptlow[i]>>pthight[i]>>v2star[i]>>ev2star[i]>>x1>>x2;
   epts[i]= 0.01;//(-ptlow[i]+pthight[i])/ptstar[i];           
  // printf("%f %f %3f\n",ptstar[i], v2star[i],ev2star[i]);
   i++;       
   }
   out1.close();

   TFile *f = new TFile("../RunOutputHisto.root");

   gStyle->SetOptStat(10000000);
   gStyle->SetStatBorderSize(0);

   TCanvas *c1 = new TCanvas("c1", "c1",364,44,699,499);
   gStyle->SetOptStat(0);
   c1->Range(-24.9362,1.71228,25.0213,4.77534);
   c1->SetFillColor(10);
   c1->SetFillStyle(4000);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameFillStyle(0);
    
   TGraph *gr1= new TGraphErrors(20,ptt,v2pi,ept,ev2pi);
   gr1->SetMarkerColor(1);
   gr1->SetLineColor(2);
   gr1->SetMarkerStyle(20);
	  
   TGraph *gr5= new TGraphErrors(20,ptstar,v2star,epts,ev2star);
   gr5->SetMarkerColor(4);
   gr5->SetLineColor(4);
   gr5->SetMarkerStyle(20);
   
   TMultiGraph *mg=new TMultiGraph();
   mg->SetMaximum(0.3);
   mg->SetMinimum(0.);
   mg->Add(gr5);          
   mg->Add(gr1);
   mg->Draw("AP");

   hv2res->SetLineWidth(2.);
   
   hv2res->Divide(hv2,hv0,1.,1.);
   hv2res->SetLineWidth(2);
   hv2res->Draw("histo:same");


   TLegend *legend=new TLegend(0.6, 0.6, 0.9, 0.9);   
   legend->AddEntry(gr5, " STAR: all charged ", "p"); 
   legend->AddEntry(hv2res, " HYDJET++ ", "l");
   legend->Draw();

   mg->GetXaxis()->SetTitleSize(0.06);
   mg->GetYaxis()->SetTitleSize(0.06); 
   mg->GetXaxis()->SetTitleOffset(0.7);
   mg->GetYaxis()->SetTitleOffset(0.7);     
   mg->GetXaxis()->SetTitle("p_{t} (GeV/c)");
   mg->GetYaxis()->SetTitle(" v_{2} ");

                 
 }
