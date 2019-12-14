//
//macro to calculate pt-spectra for pi+ and compared them with STAR data: 
//B.I.~Abelev {\it et al.} [STAR Collaboration],
//Phys. Rev. Lett. {\bf 97} (2006) 152301.


 #include "stdlib.h"
#include "stdio.h"

void fig_PTH_STAR()
{
 
//star
   std::ifstream out1("HPT_STAR/pipPHTSTAR_0-12"); 
   Double_t mtstarpi[100],dNdmtpi[100],edNdmtpi[100], epts[100];
   Int_t i=0;
   
   while(out1)
   {
   out1>>mtstarpi[i]>>dNdmtpi[i]>>edNdmtpi[i];
   epts[i]= 0.01;//(-ptlow[i]+pthight[i])/ptstar[i];   
        
  // printf("%f %f %3f\n",ptstar[i], v2star[i],ev2star[i]);
   i++;       
   }
   out1.close();

 
   TFile *f = new TFile("../RunOutput.root");
   const Int_t kMax = 20000; 

   TTree *td = (TTree*)f->Get("td");
   Int_t nevents = td->GetEntries();  
   Info("mult.C", "Nevents %d ", nevents);  
  
 
   TH1D *hpt1 = new TH1D("hpt1", "hpt1", 100, 0., 20.);
   TH1D *hpt1j = new TH1D("hpt1j", "hpt1j", 100, 0., 20.);
   TH1D *hpt1h = new TH1D("hpt1h", "hpt1h", 100, 0., 20.);

   Int_t nbx=hpt1->GetNbinsX();
   Double_t binmin=hpt1->GetXaxis()->GetXmin();
   Double_t binmax=hpt1->GetXaxis()->GetXmax();
   Double_t delta=binmax-binmin;
   Double_t delta_y=2; //[-1;1]
   

   td->Draw("sqrt(Px*Px+Py*Py)>>hpt1","(1.0/(sqrt(Px*Px+Py*Py)))*(final==1 && pdg==211 && abs(0.5*log((E+Pz)/(E-Pz)))<1.)", "goff");
   hpt1->Scale(nbx/ (2.0*TMath::Pi()*nevents*delta_y*delta));
   td->Draw("sqrt(Px*Px+Py*Py)>>hpt1j","(1.0/(sqrt(Px*Px+Py*Py)))*(final==1 && type==1 && pdg==211 && abs(0.5*log((E+Pz)/(E-Pz)))<1.)", "goff");
   hpt1j->Scale(nbx/ (2.0*TMath::Pi()*nevents*delta_y*delta)); 
   td->Draw("sqrt(Px*Px+Py*Py)>>hpt1h","(1.0/(sqrt(Px*Px+Py*Py)))*(final==1 && type==0 && pdg==211 && abs(0.5*log((E+Pz)/(E-Pz)))<1.)", "goff");
   hpt1h->Scale(nbx/ (2.0*TMath::Pi()*nevents*delta_y*delta)); 
    
   gStyle->SetOptStat(10000000);
   gStyle->SetStatBorderSize(0);

   TCanvas *c2 = new TCanvas("c2", "c2",364,44,699,499);
   gStyle->SetOptStat(0);
   c2->Range(-24.9362,1.71228,25.0213,4.77534);
   c2->SetFillColor(10);
   c2->SetFillStyle(4000);
   c2->SetBorderSize(2);
   c2->SetLogy();
   c2->SetFrameFillColor(0);
   c2->SetFrameFillStyle(0);
                                                                                              
   TGraph *gr1= new TGraphErrors(31,mtstarpi,dNdmtpi,epts,edNdmtpi);
   gr1->SetMarkerColor(1);
   gr1->SetLineColor(2);
   gr1->SetMarkerStyle(20);
 
   TMultiGraph *mg=new TMultiGraph();
   mg->Add(gr1);     
   mg->Draw("AP");
   mg->GetXaxis()->SetTitleSize(0.06);
   mg->GetYaxis()->SetTitleSize(0.06); 
   mg->GetXaxis()->SetTitleOffset(0.7); 
   mg->GetYaxis()->SetTitleOffset(0.7); 
   mg->GetXaxis()->SetTitle("p_{t} (GeV/c)");
   mg->GetYaxis()->SetTitle("1/(2 #pi) d^{2}N/ N p_{t} dp_{t} dY, c^{2}/GeV^{2}");
      
   hpt1->SetLineWidth(2.);
   hpt1j->SetLineWidth(2.);
   hpt1h->SetLineWidth(2.);
 
   hpt1->Draw("same::hist");
   hpt1j->SetLineStyle(2);
   hpt1h->SetLineStyle(3);
   hpt1h->SetLineColor(2);       
   hpt1j->Draw("same::hist");
   hpt1h->Draw("same::hist");

   TLegend *legend=new TLegend(0.6, 0.6, 0.9, 0.9);   
   legend->AddEntry(gr1, " STAR: #pi^{+} ", "p"); 
   legend->AddEntry(hpt1, " HYDJET++: all #pi^{+}  ", "l");
   legend->AddEntry(hpt1j, " HYDJET++: jet part ", "l");
   legend->AddEntry(hpt1h, " HYDJET++: hydro part ", "l");
 
   legend->Draw();


 }
