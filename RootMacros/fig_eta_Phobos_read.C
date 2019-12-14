#include "stdlib.h"
#include "stdio.h"

//
//macro to calculate
//pseudorapidity distributions for charged particles and compared them with PHOBOS data: dN/deta 
//B. B. Back et al. [PHOBOS Collaboration], Nucl.Phys. A757 (2005) 28



void fig_eta_Phobos_read()
{

  Int_t nevents = 10;
 
//star
   std::ifstream out1("eta_PHOBOS/eta_0-6"); 
   Double_t etaPHOBOS[100],eetaPHOBOS[100],dNdetaPHOBOS[100], edNPHOBOSP[100],edNPHOBOSM[100], epts[100] ;
   Int_t i=0;

   while(out1)
   {
   out1>>etaPHOBOS[i]>>eetaPHOBOS[i]>>dNdetaPHOBOS[i]>>edNPHOBOSP[i]>>edNPHOBOSM[i];
   epts[i]= 0.01;//(-ptlow[i]+pthight[i])/ptstar[i];           
   printf("%f %f %f %f %4f\n",etaPHOBOS[i],eetaPHOBOS[i],dNdetaPHOBOS[i],edNPHOBOSM[i],edNPHOBOSP[i]);
   i++;       
   }
   out1.close();
 
   TFile *f = new TFile("../RunOutputHisto.root");

   Int_t nbx=hy->GetNbinsX();
   Double_t binmin=hy->GetXaxis()->GetXmin();
   Double_t binmax=hy->GetXaxis()->GetXmax();
   Double_t delta=binmax-binmin;

   hy->Scale(nbx/(nevents*delta));
   hyjets->Scale(nbx/ (nevents*delta));
   hyhydro->Scale(nbx/ (nevents*delta));
    
   TCanvas *c2 = new TCanvas("c2", "c2",364,44,699,499);
   c2->Range(-24.9362,1.71228,25.0213,4.77534);
   c2->SetFillColor(10);
   c2->SetFillStyle(4000);
   c2->SetBorderSize(2);
   c2->SetFrameFillColor(0);
   c2->SetFrameFillStyle(0);

   gStyle->SetOptStat(0);
   gStyle->SetOptStat(10000000);
   gStyle->SetStatBorderSize(0);
 
   TGraph *gr1= new TGraphErrors(53,etaPHOBOS,dNdetaPHOBOS,epts,0);
   gr1->SetMarkerColor(2);
   gr1->SetLineColor(2);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(0.4);
    
   TMultiGraph *mg=new TMultiGraph();
   mg->Add(gr1);     
   mg->Draw("AP");
   mg->GetXaxis()->SetTitle("#eta");
   mg->GetYaxis()->SetTitle("dN/d#eta");
    
   hy->SetLineWidth(2.);
   hy->SetLineStyle(1);
   hyjets->SetLineWidth(2.);
   hyjets->SetLineStyle(2);
   hyhydro->SetLineWidth(2.);
   hyhydro->SetLineStyle(3);
    
   hy->Draw("same:histo");
   hyjets->Draw("SAME:histo");
   hyhydro->Draw("SAME:histo");

   TLegend *legend=new TLegend(0.6, 0.6, 0.9, 0.9);       
   legend->AddEntry(hy, " HYDJET++: all charged  ", "l");
   legend->AddEntry(hyjets, " HYDJET++: jet part ", "l");
   legend->AddEntry(hyhydro, " HYDJET++: hydro part ", "l");
   legend->AddEntry(gr1, " PHOBOS: all charged ", "p");
   legend->Draw();

 }
