/*******************************************************************************
 *                                                                             *
 *    HYDJET++ , event generator under the ROOT FRAMEWORK for simulation of    *
 *    relativistic heavy ion AA collisions as a superposition of the soft,     *
 *    hydro-type state and the hard, multi-parton state.                       *
 *                                                                             *
 *     The main routine is written in the object-oriented C++ language         *        
 *     under the ROOT environment. The hard, multi-partonic part of            *  
 *     HYDJET++ event is identical to the hard part of Fortran-written         *
 *     HYDJET (PYTHIA6.4xx + PYQUEN1.5) and is included in the generator       *
 *     structure as the separate directory. The soft part of HYDJET++          * 
 *     event represents the thermal hadronic state obtained with the           *
 *     Bjorken-like parameterization of freeze-out hypersurface and            *
 *     includes longitudinal, radial and elliptic flow effects and             *
 *     decays of hadronic resonances. The corresponding fast                   * 
 *     Monte-Carlo simulation procedure (C++ code) FAST MC is adapted.         *
 *     --------------------------------------------------------------          *
 *     Web-page:                                                               *
 *     http://cern.ch/lokhtin/hydjet++                                         *   
 *     --------------------------------------------------------------          *  
 *                                                                             *
 *                                                                             *
 * This program is a free software; one can use and redistribute it freely.    *  
 * Any publication of results obtained using this code must reference          * 
 *                                                                             *
 *                                                                             * 
 *                                                                             *
 *     Main reference for HYDJET++:                                            *
 *     I.P. Lokhtin, L.V. Malinina, S.V. Petrushanko, A.M. Snigirev,           *
 *     I. Arsene, K. Tywoniuk, Comp. Phys. Comm. 180 (2009) 779;               *
 *     http://cern.ch/lokhtin/hydjet++.                                        *
 *                                                                             * 
 *     Reference for HYDJET and PYQUEN:                                        *
 *     I.P. Lokhtin, A.M. Snigirev, Eur. Phys. J. C 46 (2006) 211;             *
 *     http://cern.ch/lokhtin/hydro/hydjet.html                                * 
 *     http://cern.ch/lokhtin/pyquen.                                          *  
 *                                                                             *    
 *     Reference for PYTHIA6.4:                                                *
 *     T.Sjostrand, S. Mrenna and P. Skands, JHEP05 (2006) 026;                *
 *     http://home.thep.lu.se/~torbjorn/Pythia.html.                           * 
 *                                                                             * 
 *     References for FAST MC:                                                 *  
 *     N.S. Amelin, R. Lednicky, T.A. Pocheptsov, I.P. Lokhtin,                * 
 *     L.V. Malinina, A.M. Snigirev, Iu.A. Karpenko and Yu.M. Sinyukov,        * 
 *     Phys. Rev. C 74 (2006) 064901;                                          *
 *     N.S. Amelin, I. Arsene, L. Bravina, Iu.A. Karpenko, R. Lednicky,        *  
 *     I.P. Lokhtin, L.V. Malinina, A.M. Snigirev and Yu.M. Sinyukov,          *  
 *     Phys. Rev. C 77 (2008) 014903;                                          *
 *     http://uhkm.jinr.ru.                                                    *   
 *                                                                             *
 *     Reference for nuclear shadowing model:                                  *
 *     K. Tywoniuk, I.C. Arsene, L. Bravina, A. Kaidalov and                   *
 *     E. Zabrodin, Phys. Lett. B 657 (2007) 170.                              *
 *                                                                             * 
 *     Version 2.4:                                                            *
 *                                                                             *
 *     Igor Lokhtin, SINP MSU, Moscow, RU                                      *
 *     e-mail: Igor.Lokhtin@cern.ch                                            *
 *                                                                             *
 *     Ludmila Malinina, SINP MSU, Moscow, RU                                  *   
 *     e-mail: malinina@lav01.sinp.msu.ru                                      * 
 *                                                                             *                                                              *
 *******************************************************************************/ 
#include <iostream> 
#include <fstream>
#include <vector>
#include <time.h>

#include <TNtuple.h>
#include <TError.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TFile.h>

#include "InitialState.h"
#include "InitialStateHydjet.h"


#include <TRandom.h>

#include "Particle.h"
#include "HYJET_COMMONS.h"
//extern SERVICECommon SERVICE;

extern SERVICEEVCommon SERVICEEV;

//Main program:
//reads input parameters from file "RunInputBjorken" or "RunInputHubble";
//calculates particle densities and average initial multiplicities and writes them
//in output file "multiplicities.txt";
//creates trees (tree with direct hadrons and hadrons after resonance decays)
//with space-time and momentum-energy information of produced hadrons;
//writes trees in file "RunOutput.root".

Int_t main(Int_t argc, Char_t** argv) {

  clock_t start;
  start = clock();

 
//new
  time_t  now;
  struct tm  *ts;
  char       buf[80];
         
 // Get the current time
   time(&now);
              
 // Format and print the time, "ddd yyyy-mm-dd hh:mm:ss zzz"
    ts = localtime(&now);
    strftime(buf, sizeof(buf), "%a %Y-%m-%d %H:%M:%S %Z", ts);
    printf("%s\n", buf);
 

std::cout<<" *******************************************************************************"<<std::endl;
std::cout<<" *                         HYDJET++ version 2.4:                               *"<<std::endl;
std::cout<<" *                      Last revision: 03-DEC-2018                             *"<<std::endl;
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" *    event generator under the ROOT FRAMEWORK for simulation of               *"<<std::endl;
std::cout<<" *    relativistic heavy ion AA collisions as a superposition of the           *"<<std::endl;
std::cout<<" *    soft, hydro-type state and hard, multi-parton state.                     *"<<std::endl;
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" *     The main routine is written in the object-oriented C++ language         *"<<std::endl;        
std::cout<<" *     under the ROOT environment. The hard, multi-partonic part of            *"<<std::endl;  
std::cout<<" *     HYDJET++ event is identical to the hard part of Fortran-written         *"<<std::endl;
std::cout<<" *     HYDJET (PYTHIA6.4xx + PYQUEN1.5) and is included in the generator       *"<<std::endl;
std::cout<<" *     structure as the separate directory. The soft part of HYDJET++          *"<<std::endl; 
std::cout<<" *     event represents the thermal hadronic state obtained with the           *"<<std::endl;
std::cout<<" *     Bjorken-like parameterization of freeze-out hypersurface and            *"<<std::endl;
std::cout<<" *     includes longitudinal, radial and elliptic flow effects and             *"<<std::endl;
std::cout<<" *     decays of hadronic resonances. The corresponding fast                   *"<<std::endl; 
std::cout<<" *     Monte-Carlo simulation procedure (C++ code) FAST MC is adapted.         *"<<std::endl;
std::cout<<" *     --------------------------------------------------------------          *"<<std::endl;
std::cout<<" *     Web-page:                                                               *"<<std::endl;
std::cout<<" *     http://cern.ch/lokhtin/hydjet++                                         *"<<std::endl;   
std::cout<<" *     --------------------------------------------------------------          *"<<std::endl;
std::cout<<" *                                                                             *"<<std::endl;                                                                            
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" * This program is a free software, one can use and redistribute it freely.     *"<<std::endl;  
std::cout<<" * Any publication of results obtained using this code must reference          *"<<std::endl; 
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" *                                                                             *"<<std::endl; 
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" *     Main reference for HYDJET++:                                            *"<<std::endl;
std::cout<<" *     I.P. Lokhtin, L.V. Malinina, S.V. Petrushanko, A.M. Snigirev,           *"<<std::endl;
std::cout<<" *     I. Arsene, K. Tywoniuk, Comp. Phys. Comm. 180 (2009) 779                *"<<std::endl;
std::cout<<" *     http://cern.ch/lokhtin/hydjet++.                                        *"<<std::endl;
std::cout<<" *                                                                             *"<<std::endl; 
std::cout<<" *     Reference for HYDJET and PYQUEN:                                        *"<<std::endl;
std::cout<<" *     I.P. Lokhtin, A.M. Snigirev, Eur. Phys. J. C 46 (2006) 211              *"<<std::endl;
std::cout<<" *     http://cern.ch/lokhtin/hydro/hydjet.html                                *"<<std::endl; 
std::cout<<" *     http://cern.ch/lokhtin/pyquen.                                          *"<<std::endl;  
std::cout<<" *                                                                             *"<<std::endl;    
std::cout<<" *     Reference for PYTHIA6.4:                                                *"<<std::endl;
std::cout<<" *     T.Sjostrand, S. Mrenna and P. Skands, JHEP05 (2006) 026;                *"<<std::endl;
std::cout<<" *     http://home.thep.lu.se/~torbjorn/Pythia.html.                           *"<<std::endl; 
std::cout<<" *                                                                             *"<<std::endl; 
std::cout<<" *     References for FAST MC:                                                 *"<<std::endl;  
std::cout<<" *     N.S. Amelin, R. Lednicky, T.A. Pocheptsov, I.P. Lokhtin,                *"<<std::endl; 
std::cout<<" *     L.V. Malinina, A.M. Snigirev, Iu.A. Karpenko and Yu.M. Sinyukov,        *"<<std::endl; 
std::cout<<" *     Phys. Rev. C 74 (2006) 064901;                                          *"<<std::endl;
std::cout<<" *     N.S. Amelin, I. Arsene, L. Bravina, Iu.A. Karpenko, R. Lednicky,        *"<<std::endl;  
std::cout<<" *     I.P. Lokhtin, L.V. Malinina, A.M. Snigirev and Yu.M. Sinyukov,          *"<<std::endl;  
std::cout<<" *     Phys. Rev. C 77 (2008) 014903                                           *"<<std::endl;
std::cout<<" *     http://uhkm.jinr.ru.                                                    *"<<std::endl;   
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" *     Reference for nuclear shadowing model:                                  *"<<std::endl;
std::cout<<" *     K. Tywoniuk, I.C. Arsene, L. Bravina, A. Kaidalov and                   *"<<std::endl;
std::cout<<" *     E. Zabrodin, Phys. Lett. B 657 (2007) 170.                              *"<<std::endl;
std::cout<<" *                                                                             *"<<std::endl; 
std::cout<<" *     Version 2.4:                                                          *"<<std::endl;
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" *     Igor Lokhtin, SINP MSU, Moscow, RU                                      *"<<std::endl;
std::cout<<" *     e-mail: Igor.Lokhtin@cern.ch                                            *"<<std::endl;
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" *     Ludmila Malinina, SINP MSU, Moscow, RU                                  *"<<std::endl;   
std::cout<<" *     e-mail: malinina@lav01.sinp.msu.ru                                      *"<<std::endl; 
std::cout<<" *                                                                             *"<<std::endl;
std::cout<<" *******************************************************************************"<<std::endl; 



 
 
  TString outputFilename;
  if(argc>1) outputFilename = argv[1];
  else outputFilename = "RunOutputHisto.root";
  TFile *outputFile=new TFile(outputFilename.Data(), "RECREATE"); 
  if(outputFile->IsZombie()) {
    std::cout << "Error: The output file " << outputFilename.Data() << " could not be opened !!!" << std::endl;
    return 1;
  }

  //SET MAXIMAl VALUE OF PARTICLE MULTIPLICITY!!!
  const Int_t kMax = 500000; 
  //define event number
  Int_t nev;
  //define hadron characteristic vectors
  std::vector<Int_t> pdg(kMax); //pdg encodings
  std::vector<Int_t> Mpdg(kMax);//pdg encodings for mother hadrons
  std::vector<Int_t> type(kMax);//type: 0-from hydro or decay, 1 from jets
  std::vector<Int_t> pythiaStatus(kMax); // pythia status code
  std::vector<Float_t> Px(kMax);//x-hadron momentum component,[GeV/c]
  std::vector<Float_t> Py(kMax);//y-hadron momentum component,[GeV/c]
  std::vector<Float_t> Pz(kMax);//z-hadron momentum component,[GeV/c]
  std::vector<Float_t> E(kMax); //hadron total energy,[GeV]  
  std::vector<Float_t> X(kMax);//x-hadron coordinate component,[fm]
  std::vector<Float_t> Y(kMax);//y-hadron coordinate component,[fm]
  std::vector<Float_t> Z(kMax);//z-hadron coordinate component,[fm]
  std::vector<Float_t> T(kMax);//hadron time,[fm/c] 

  TH1D *hpt1 = new TH1D("hpt1", "hpt1", 100, 0., 20.);
  TH1D *hpt1j = new TH1D("hpt1j", "hpt1j", 100, 0., 20.);
  TH1D *hpt1h = new TH1D("hpt1h", "hpt1h", 100, 0., 20.);

  
  TH1D *hv2 = new TH1D("hv2", "hv2", 100, 0.0, 10.);
  TH1D *hv3 = new TH1D("hv3", "hv3", 100, 0.0, 10.);
  TH1D *hv4 = new TH1D("hv4", "hv4", 100, 0.0, 10.);
  TH1D *hv5 = new TH1D("hv5", "hv5", 100, 0.0, 10.);
  TH1D *hv6 = new TH1D("hv6", "hv6", 100, 0.0, 10.);
  TH1D *hv0 = new TH1D("hv0", "hv0", 100, 0.0, 10.);

  TH1D *hv3ps = new TH1D("hv3ps", "hv3ps", 100, 0.0, 10.);
  TH1D *hv5ps = new TH1D("hv5ps", "hv5ps", 100, 0.0, 10.);

  TH1D *hv2e = new TH1D("hv2e", "hv2e", 50, -5, 5.);
  TH1D *hv0e = new TH1D("hv0e", "hv0e", 50, -5, 5.);


  TH1D *hv2h = new TH1D("hv2h", "hv2h", 100, 0.0, 10.);
  TH1D *hv0h = new TH1D("hv0h", "hv0h", 100, 0.0, 10.);

  TH1D *hv2eh = new TH1D("hv2eh", "hv2eh", 50, -5, 5.);
  TH1D *hv0eh = new TH1D("hv0eh", "hv0eh", 50, -5, 5.);


  TH1D *hv2j = new TH1D("hv2j", "hv2j", 100, 0.0, 10.);
  TH1D *hv0j = new TH1D("hv0j", "hv0j", 100, 0.0, 10.);

  TH1D *hv2ej = new TH1D("hv2ej", "hv2ej", 50, -5, 5.);
  TH1D *hv0ej = new TH1D("hv0ej", "hv0ej", 50, -5, 5.);



  
  TH1D *hy = new TH1D("hy", "hy", 51, -5.1, 5.1);
  TH1D *hyjets = new TH1D("hyjets", "hyjets", 51, -5.1, 5.1);
  TH1D *hyhydro = new TH1D("hyhydro", "hyhydro", 51, -5.1, 5.1);

  double pdg1, Mpdg1, Px1, Py1, E1, Pz1, pt, phi, v2, v3, v4, v5, v6, eta, v3ps, v3psi, v5ps;
  int type1;

  InitialState *FASTMC;

  FASTMC = new InitialStateHydjet();
   
  if(!FASTMC->ReadParams()) {
    Error("RunHadronSource::main", "No initial model parameters found!!\n");
    return 0;
  }

  if(!FASTMC->MultIni()) {
    Error("RunHadronSource::main", "Initial multiplicities are zero!!\n");
    return 0;
  }

  ParticleAllocator allocator;
  List_t particleList;
  std::cout << "Generating " << FASTMC->GetNev() << " events" << std::endl;
  std::cout << "Starting the event loop" << std::endl;
  
  
  double nnnn = 0.;
  double vvvv = 0.;
    
  // Loop over events  
  for(Int_t ev = 0; ev < FASTMC->GetNev(); ++ev) {
    nev = ev;
    // Initialize the source
    FASTMC->Initialize(particleList, allocator);

    if(particleList.empty()){
      if(/*TODO: Not allow empty events prarameter here*/false){
        --ev;
        continue;
      }
      std::cout<<" Warning! Empty event! "<<std::endl;  
    }else if(FASTMC->RunDecays()) FASTMC->Evolve(particleList, allocator, FASTMC->GetWeakDecayLimit());
   
    std::cout << "event #" << ev << "\r" << std::flush;
    LPIT_t it;
    LPIT_t e;
    if(!particleList.empty())      
     for(it = particleList.begin(), e = particleList.end(); it != e; ++it) {
      TVector3 pos(it->Pos().Vect());
      TVector3 mom(it->Mom().Vect());
      Float_t m1 = it->TableMass();
      pdg1 = it->Encoding();
      Mpdg1 = it->GetLastMotherPdg();
      Px1 = mom[0];
      Py1 = mom[1];
      Pz1 = mom[2];
      E1 =  TMath::Sqrt(mom.Mag2() + m1*m1);
      type1 = it->GetType();  // 0-hydro; 1-jet
      Int_t final = 0;
      if(type1==0 && it->GetNDaughters()==0) // hydro
	final=1;
      if(type1==1 && it->GetNDaughters()==0 && it->GetPythiaStatusCode()==1) // jet
	final=1;
      if(final==1 && pdg1==211 && TMath::Abs(0.5*log((E1+Pz1)/(E1-Pz1)))<1.) {
	hpt1->Fill(sqrt(Px1*Px1+Py1*Py1),1./sqrt(Px1*Px1+Py1*Py1));
      }
      
     double etaeta=0.5*TMath::Log((sqrt(Px1*Px1+Py1*Py1+Pz1*Pz1)+Pz1)/(sqrt(Px1*Px1+Py1*Py1+Pz1*Pz1)-Pz1));

      if(final==1 && pdg1==211 && TMath::Abs(0.5*log((E1+Pz1)/(E1-Pz1)))<0.8 && type1==0) 
        hpt1h->Fill(sqrt(Px1*Px1+Py1*Py1),1./sqrt(Px1*Px1+Py1*Py1));
      if(final==1 && pdg1==211 && TMath::Abs(0.5*log((E1+Pz1)/(E1-Pz1)))<0.8 && type1==1) 
        hpt1j->Fill(sqrt(Px1*Px1+Py1*Py1),1./sqrt(Px1*Px1+Py1*Py1));

      if(final==1 && ((TMath::Abs(pdg1)==211)||(TMath::Abs(pdg1)==321)||(TMath::Abs(pdg1)==2212)) ) {

	pt = TMath::Sqrt(Px1*Px1+Py1*Py1);      
	phi = TMath::ATan2(Py1,Px1);
	v2 = TMath::Cos(2*phi);       
	v3 = TMath::Cos(3*phi);       
	v4 = TMath::Cos(4*phi);       
	v5 = TMath::Cos(5*phi);       
	v6 = TMath::Cos(6*phi);       
       
        v3psi = SERVICEEV.psiv3;
        v3ps = TMath::Cos(3*(phi-v3psi));
        v5ps = TMath::Cos(5*(phi-v3psi));
        
//  std::cout << " v3psi " << v3psi <<  " " << phi <<std::endl;

// v2 all

	if (TMath::Abs(etaeta)<0.8) {
	  hv2->Fill(pt,v2);
	  hv3->Fill(pt,v3);
	  hv4->Fill(pt,v4);
	  hv5->Fill(pt,v5);
	  hv6->Fill(pt,v6);
          hv0->Fill(pt,1.);
          hv3ps->Fill(pt,v3ps);
          hv5ps->Fill(pt,v5ps);
	 }

        if (pt>0.3 && pt<3) hv2e->Fill(etaeta,v2);
        if (pt>0.3 && pt<3) hv0e->Fill(etaeta,1.);


// hydro

        if (TMath::Abs(etaeta)<0.8 && type1==0) hv2h->Fill(pt,v2);
        if (TMath::Abs(etaeta)<0.8 && type1==0) hv0h->Fill(pt,1.);

        if (pt>0.3 && pt<3 && type1==0) hv2eh->Fill(etaeta,v2);
        if (pt>0.3 && pt<3 && type1==0) hv0eh->Fill(etaeta,1.);

// jets

        if (TMath::Abs(etaeta)<0.8 && type1==1) hv2j->Fill(pt,v2);
        if (TMath::Abs(etaeta)<0.8 && type1==1) hv0j->Fill(pt,1.);

        if (pt>0.3 && pt<3 && type1==1) hv2ej->Fill(etaeta,v2);
        if (pt>0.3 && pt<3 && type1==1) hv0ej->Fill(etaeta,1.);


       
        if (pt>0.3 && pt<3 && TMath::Abs(etaeta)<0.8) 
        {nnnn = nnnn + 1;
         vvvv = vvvv + v2;
        }

      }
       
      if(final==1 && ((TMath::Abs(pdg1)==211)||(TMath::Abs(pdg1)==321)||(TMath::Abs(pdg1)==2212))){    
	eta=0.5*TMath::Log((sqrt(Px1*Px1+Py1*Py1+Pz1*Pz1)+Pz1)/(sqrt(Px1*Px1+Py1*Py1+Pz1*Pz1)-Pz1));
	if(type1==1)hyjets->Fill(eta);
	if(type1==0)hyhydro->Fill(eta);
	hy->Fill(eta);
      }
    }
 
    allocator.FreeList(particleList); 
  }
  
  hpt1->Write();
  hpt1h->Write();
  hpt1j->Write();

  hv2->Write();
  hv3->Write();
  hv4->Write();
  hv5->Write();
  hv6->Write();
  
  hv3ps->Write();
  hv5ps->Write();

  hv0->Write();
  hv2e->Write();
  hv0e->Write();

  hv2h->Write();
  hv0h->Write();
  hv2eh->Write();
  hv0eh->Write();

  hv2j->Write();
  hv0j->Write();
  hv2ej->Write();
  hv0ej->Write();

  hyhydro->Write();
  hyjets->Write();
  hy->Write();
  
  outputFile->Close();
  
  clock_t stop;
  stop = clock();
  
  std::cout << " v2 integral " << vvvv/nnnn << "  "<< nnnn << std::endl;
  
  std::cout << "*********************************************" << std::endl;
  std::cout << "Execution time: " << (stop - start)/CLOCKS_PER_SEC << " seconds" << std::endl;
  std::cout << "*********************************************" << std::endl;

  time_t  now1;
  struct tm  *ts1;
  char       buf1[80];
         
  // Get the current time
  time(&now1);
              
  // Format and print the time, "ddd yyyy-mm-dd hh:mm:ss zzz"
  ts1 = localtime(&now1);
  strftime(buf1, sizeof(buf1), "%a %Y-%m-%d %H:%M:%S %Z", ts1);
  printf("%s\n", buf1);
  
  return 0;
}
