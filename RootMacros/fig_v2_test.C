//
//macro to calculate v2(pt) for charged particles and compared them with STAR data: 
//J.~Adams {\it et al.} [STAR Collaboration], Phys. Rev. {\bf C 72} (2005) 014904.  


 #include "stdlib.h"
#include "stdio.h"

void fig_v2_test()
{
  Double_t v2pi[100], v1pi[100], v2K[100], v2p[100], ev2pi[100], ev2K[100], ev2p[100],pt, ptt[100], ept[100];
  Int_t n2pi[100], n2K[100], n2p[100];

   std::ifstream out1("v2_STAR/v2star0-5"); 
   Double_t ptstar[100],ptlow[100],pthight[100],v2star[100],ev2star[100],epts[100],x1,x2;
   Int_t i=0;
   
   while(out1)
   {
   out1>>ptstar[i]>>ptlow[i]>>pthight[i]>>v2star[i]>>ev2star[i]>>x1>>x2;
   epts[i]= 0.01;//(-ptlow[i]+pthight[i])/ptstar[i];   
        
   printf("%f %f %3f\n",ptstar[i], v2star[i],ev2star[i]);
   i++;       
   }
   out1.close();
 
   TH1D *hv2 = new TH1D("hv2", "hv2", 100, 0.0, 10.);
   TH1D *hv0 = new TH1D("hv0", "hv0", 100, 0.0, 10.);
   TH1D *hv2res1 = new TH1D("hv2res1", "hv2res1", 100, 0.0, 10.);
   
   Int_t nev = 0;
   Float_t v2, phi;
 
   TFile *f = new TFile("../RunOutput.root");
   const Int_t kMax = 200000; 
 
   Int_t Ntot;
   Int_t   pdg[kMax];
   Int_t   Mpdg[kMax];
   Float_t Px[kMax];
   Float_t Py[kMax];
   Float_t Pz[kMax];
   Float_t E[kMax];   
   Float_t X[kMax];
   Float_t Y[kMax];
   Float_t Z[kMax];
   Float_t T[kMax]; 
   Int_t type[kMax]; 
   Int_t NDaughters[kMax]; 
   Int_t final[kMax];

  TTree *td = (TTree*)f->Get("td");
  Int_t nevents = td->GetEntries();  
  Info("mult.C", "Nevents %d ", nevents);  
 
  td->SetBranchAddress("Ntot",&Ntot);
  td->SetBranchAddress("Px",Px);
  td->SetBranchAddress("Py",Py);
  td->SetBranchAddress("Pz",Pz);
  td->SetBranchAddress("E",E);
  td->SetBranchAddress("pdg",pdg);
  td->SetBranchAddress("Mpdg",Mpdg);
  td->SetBranchAddress("NDaughters",NDaughters);
  td->SetBranchAddress("final",final);
  
  for (Int_t k=0;k<nevents;k++) 
  {
            
      td->GetEntry(k);  
              
     for (Int_t i=0;i<Ntot;i++) 
	   {

//all charged
      if(((pdg[i]==211)||(abs(pdg[i])==321)||(abs(pdg[i])==2212)) 
      && (abs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i]+0.0001)))<1.0) && final[i]==1){

      pt = TMath::Sqrt(Px[i]*Px[i]+Py[i]*Py[i]);      
      phi = TMath::ATan2(Py[i],Px[i]);
  
      v2 = TMath::Cos(2*phi); 
      
      hv2->Fill(pt,v2);
      hv0->Fill(pt,1.);
      
    }     
  } 
} 


    TCanvas *c1= new TCanvas("c1","c1",200,10,800,600);
                                                                                              
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
    mg->GetXaxis()->SetTitle("pt, GeV/c");
    mg->GetYaxis()->SetTitle(" v2 ");
    hv2res1->Divide(hv2,hv0,1.,1.);
    hv2res1->SetLineWidth(2);
    hv2res1->Draw("histo:same");


   TLegend *legend=new TLegend(0.6, 0.6, 0.9, 0.9);   
   legend->AddEntry(gr5, " STAR: all charged ", "p"); 
   legend->AddEntry(hv2res1, " HYDJET++ ", "l");
   legend->Draw();


                 
 }
