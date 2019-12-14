/*                                                                           
         HYDJET++ 
         version 2.4:  
         InitialStateHydjet is the modified InitialStateBjorken 
         The high-pt part related with PYTHIA-PYQUEN is included       
         InitialStateBjorken (FASTMC) was used.


         
         InitialStateBjorken           
         version 2.0: 
         Ludmila Malinina  malinina@lav01.sinp.msu.ru,   SINP MSU/Moscow and JINR/Dubna
         Ionut Arsene  i.c.arsene@fys.uio.no,            Oslo University                                                
                     June 2007
        
         version 1.0:                                                               
         Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
         amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005 

                     
*/


//expanding localy equilibated fireball with volume hadron radiation
//thermal part: Blast wave model, Bjorken-like parametrization
//high-pt: PYTHIA + jet quenching model PYQUEN


#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TMath.h>

#include "InitialStateHydjet.h"
#include "RandArrayFunction.h"
#include "HadronDecayer.h"
#include "GrandCanonical.h"
#include "StrangePotential.h"
#include "EquationSolver.h"
#include "Particle.h"
#include "ParticlePDG.h"
#include "UKUtility.h"
#include <iostream> 
#include <fstream>
#include "HYJET_COMMONS.h"
extern "C" void  hyevnt_();
extern "C" void  myini_();
//extern "C" void  hyinit_();

extern HYIPARCommon HYIPAR;
extern HYFPARCommon HYFPAR;
extern HYJPARCommon HYJPAR;
extern HYPARTCommon HYPART;
extern SERVICECommon SERVICE;
extern SERVICEEVCommon SERVICEEV;

using std::cout;
using std::endl;

// definition of the static member fLastIndex
Int_t Particle::fLastIndex;

void InitialStateHydjet::Initialize(List_t &source, ParticleAllocator & allocator) {
  // Initialize the static "last index variable"
  Particle::InitIndexing();
 
  //----- high-pt part------------------------------
  TLorentzVector partJMom, partJPos, zeroVec;
  
  //generate non-equilibrated part event
  hyevnt_();

  if(fParams.fNhsel != 0){   
    //get number of particles in jets
    Int_t numbJetPart = HYPART.njp;
    
    for(Int_t i = 0; i <numbJetPart; ++i) {
      Int_t pythiaStatus    = Int_t(HYPART.ppart[i][0]);   // PYTHIA status code
      Int_t pdg             = Int_t(HYPART.ppart[i][1]);   // PYTHIA species code
      Double_t px           = HYPART.ppart[i][2];          // px
      Double_t py           = HYPART.ppart[i][3];          // py
      Double_t pz           = HYPART.ppart[i][4];          // pz
      Double_t e            = HYPART.ppart[i][5];          // E
      Double_t vx           = HYPART.ppart[i][6];          // x
      Double_t vy           = HYPART.ppart[i][7];          // y
      Double_t vz           = HYPART.ppart[i][8];          // z
      Double_t vt           = HYPART.ppart[i][9];          // t
      // particle line number in pythia are 1 based while we use a 0 based numbering
      Int_t mother_index    = Int_t(HYPART.ppart[i][10])-1;  //line number of parent particle
      Int_t daughter_index1 = Int_t(HYPART.ppart[i][11])-1;  //line number of first daughter
      Int_t daughter_index2 = Int_t(HYPART.ppart[i][12])-1;  //line number of last daughter

      // For status codes 3, 13, 14 the first and last daughter indexes have a different meaning
      // used for color flow in PYTHIA. So these indexes will be reset to zero.
      if(TMath::Abs(daughter_index1)>numbJetPart || TMath::Abs(daughter_index2)>numbJetPart ||
	 TMath::Abs(daughter_index1)>TMath::Abs(daughter_index2)) {
	daughter_index1 = -1;
	daughter_index2 = -1;
      }
      
      ParticlePDG *partDef = fDatabase->GetPDGParticle(pdg);
            
      Int_t type=1; //from jet
      if(partDef) {
	Int_t motherPdg = Int_t(HYPART.ppart[mother_index][1]);
	if(motherPdg==0) motherPdg = -1;
	partJMom.SetXYZT(px, py, pz, e); 
	partJPos.SetXYZT(vx, vy, vz, vt);
	Particle particle(partDef, partJPos, partJMom, 0, 0, type, motherPdg, zeroVec, zeroVec);
	Int_t index = particle.SetIndex();
	if(index!=i) {
	  cout << "InitialStateHydjet::Initialize(): Allocated HYDJET++ index is not synchronized with the PYTHIA index!" << endl
	       << "                              Collision history information is destroyed! It happens when a PYTHIA code is not" << endl
	       << "                              implemented in HYDJET++ particle list particles.data! Check it out!" << endl;
	}
	particle.SetPythiaStatusCode(pythiaStatus);
	particle.SetMother(mother_index);
	particle.SetFirstDaughterIndex(daughter_index1);
	particle.SetLastDaughterIndex(daughter_index2);
	if(pythiaStatus!=1) particle.SetDecayed();
	allocator.AddParticle(particle, source);
      }
      else {
	cout << "InitialStateHydjet::Initialize(): PYTHIA particle of specie " << pdg << " is not in HYDJET++ particle list" << endl
	     <<                                  "Please define it in particles.data, otherwise the history information will be de-synchronized and lost!" << endl;
      }
    }
  } //nhsel !=0 not only hydro!        

  //----------HYDRO part------------------------------------------------

  // if fluctuations are included fv2!=0 and fv3!=0 

  Double_t kv2 = fParams.fv2; //0.34 or 0
  Double_t kv3 = fParams.fv3; //0.52

    //-------------------------------------
    // get impact parameter    
   Double_t impactParameter = HYFPAR.bgen;
   Double_t e3_old,  Eps3Fluc, Epsilon_old, Eps2Fluc;

   Double_t psiforv3 = 0.;

   Double_t e3 = (0.2/5.5)*TMath::Power(impactParameter,1./3.);
 
    if(kv3!=0){ 
    e3_old = e3;
    Eps3Fluc = e3_old*kv3;
    e3 =gRandom->Gaus(e3_old,Eps3Fluc);
    }

   psiforv3 = TMath::TwoPi() *  (-0.5 + gRandom->Rndm())  / 3.;
   SERVICEEV.psiv3 = -psiforv3;

  if(fParams.fNhsel < 3){    
    const Double_t  weightMax = 2*TMath::CosH(fParams.fUmax);
    const Int_t nBins = 100;
    Double_t probList[nBins];
    RandArrayFunction arrayFunctDistE(nBins);
    RandArrayFunction arrayFunctDistR(nBins);
    
    TLorentzVector partPos, partMom, n1, p0;
    TVector3 vec3;
    const TLorentzVector zeroVec;
    //set maximal hadron energy
    const Double_t eMax = 5.;  
    
    Double_t Delta = fParams.fDelta;
    Double_t Epsilon = fParams.fEpsilon; //our standart 

    // if fluctuatins are requested
    if(kv2!=0){
    Epsilon_old = Epsilon;
    Eps2Fluc = Epsilon_old*kv2; 
    Epsilon =gRandom->Gaus(Epsilon_old, Eps2Fluc);
    }

    
    if(fParams.fIfDeltaEpsilon>0){
      Double_t Epsilon0 = 0.5*impactParameter; //e0=b/2Ra


     Double_t coeff;
    if(kv2==0 && kv3==0){
       coeff = (HYIPAR.RA/fParams.fR)/12.;//phenomenological coefficient
     }
     else
       {
       coeff = (HYIPAR.RA/fParams.fR)/9.9;//phenomenological coefficient
       }


      Epsilon = Epsilon0 * coeff; 

       if(kv2!=0){
          Epsilon =  Epsilon0 * coeff;   
          Epsilon_old = Epsilon;
          Eps2Fluc = Epsilon_old*kv2; 
          Epsilon =gRandom->Gaus(Epsilon_old, Eps2Fluc);
           }

      Double_t C=5.6; 
      Double_t A = C*Epsilon*(1-Epsilon*Epsilon);

      if(TMath::Abs(Epsilon)<0.0001 || TMath::Abs(A)<0.0001 )Delta=0.0;
      if(TMath::Abs(Epsilon)>0.0001 && TMath::Abs(A)>0.0001)Delta = 0.5*(TMath::Sqrt(1+4*A*(Epsilon+A))-1)/A; 
     
    }
    
    //effective volume for central     
    double dYl= 2 * fParams.fYlmax; //uniform distr. [-Ylmax; Ylmax]  
    if (fParams.fEtaType >0) dYl = TMath::Sqrt(2 * TMath::Pi()) * fParams.fYlmax ;  //Gaussian distr.                                                                            
    Double_t VolEffcent = 2 * TMath::Pi() * fParams.fTau * dYl * (fParams.fR * fParams.fR)/TMath::Power((fParams.fUmax),2) * 
      ((fParams.fUmax)*TMath::SinH((fParams.fUmax))-TMath::CosH((fParams.fUmax))+ 1);
    
    //effective volume for non-central Simpson2 
    Double_t VolEffnoncent = fParams.fTau * dYl * SimpsonIntegrator2(0., 2.*TMath::Pi(), Epsilon, Delta);
    
    fVolEff = VolEffcent * HYFPAR.npart/HYFPAR.npart0;
    
    Double_t coeff_RB = TMath::Sqrt(VolEffcent * HYFPAR.npart/HYFPAR.npart0/VolEffnoncent);
    Double_t coeff_R1 = HYFPAR.npart/HYFPAR.npart0;
    coeff_R1 = TMath::Power(coeff_R1, 0.333333);
    
    double Veff=fVolEff;
    
    
    //------------------------------------
    //cycle on particles types

    Double_t Nbcol = HYFPAR.nbcol;
    Double_t NccNN = SERVICE.charm;
    Double_t Ncc = Nbcol * NccNN/dYl;
    Double_t Nocth = fParams.fNocth;
    Double_t NJPsith = fParams.fNccth;
   
    Double_t gammaC=1.0;
    if(fParams.fCorrC<=0){
     gammaC=CharmEnhancementFactor(Ncc, Nocth, NJPsith, 0.001);
     }
     else
       {
        gammaC=fParams.fCorrC; 
       }


    for(Int_t i = 0; i < fParams.fNPartTypes; ++i) {
      double Mparam = fParams.fPartMult[2 * i] * Veff;
      const Int_t encoding = fParams.fPartEnc[i];

	ParticlePDG *partDef0 = fDatabase->GetPDGParticle(encoding);
       
	if(partDef0->GetCharmQNumber()!=0 || partDef0->GetCharmAQNumber()!=0)Mparam = Mparam * gammaC;
        if(abs(encoding)==443)Mparam = Mparam * gammaC;  
    
   
        Int_t multiplicity = gRandom->Poisson(Mparam);

      //      cout << "specie: " << encoding << "; average mult: = " << Mparam << "; multiplicity = " << multiplicity << endl;
      
      if (multiplicity > 0) {
	ParticlePDG *partDef = fDatabase->GetPDGParticle(encoding);
	if(!partDef) {
	  Error("InitialStateHydjet::Initialize", "No particle with encoding %d", encoding);
	  continue;
	}
	
	 
	if(fParams.fCharmProd<=0 && (partDef->GetCharmQNumber()!=0 || partDef->GetCharmAQNumber()!=0)){
	  //cout<<"statistical charmed particle not allowed ! "<<encoding<<endl;
	  continue;
	}
	
	//compute chemical potential for single f.o. mu==mu_ch
	//compute chemical potential for thermal f.o.                
	Double_t mu = fParams.fPartMu[2 * i];
               
	//choose Bose-Einstein or Fermi-Dirac statistics
	const Double_t d    = !(Int_t(2*partDef->GetSpin()) & 1) ? -1 : 1;
	const Double_t mass = partDef->GetMass();                	 
	
	//prepare histogram to sample hadron energy: 
	Double_t h = (eMax - mass) / nBins;
	Double_t x = mass + 0.5 * h;
	Int_t i;        
	for(i = 0; i < nBins; ++i) {
	  if(x>=mu && fParams.fThFO>0)probList[i] = x * TMath::Sqrt(x * x - mass * mass) / 
					(TMath::Exp((x - mu) / (fParams.fThFO)) + d);
	  if(x>=mu && fParams.fThFO<=0)probList[i] = x * TMath::Sqrt(x * x - mass * mass) / 
					 (TMath::Exp((x - mu) / (fParams.fT)) + d);								 
	  if(x<mu)probList[i] = 0.; 
	  x += h;
	}
	arrayFunctDistE.PrepareTable(probList);
	
	//prepare histogram to sample hadron transverse radius: 
	h = (fParams.fR) / nBins;
	x =  0.5 * h;
	Double_t param = (fParams.fUmax) / (fParams.fR);
	for (i = 0; i < nBins; ++i) {
	  probList[i] = x * TMath::CosH(param*x);
	  x += h;
	}
	arrayFunctDistR.PrepareTable(probList);
	
	//loop over hadrons, assign hadron coordinates and momenta
	Double_t weight = 0., yy = 0., px0 = 0., py0 = 0., pz0 = 0.;
	Double_t e = 0., x0 = 0., y0 = 0., z0 = 0., t0 = 0., etaF = 0.; 
	Double_t r, RB, phiF;

        RB = fParams.fR * coeff_RB * coeff_R1 * TMath::Sqrt((1+TMath::Abs(e3))/(1-TMath::Abs(e3)));
        


	for(Int_t j = 0; j < multiplicity; ++j) {               
	  do {

	M1: fParams.fEtaType <=0 ? etaF = fParams.fYlmax * (2. * gRandom->Rndm() - 1.) 
	      : etaF = (fParams.fYlmax) * (gRandom->Gaus());                                                                              
	    n1.SetXYZT(0.,0.,TMath::SinH(etaF),TMath::CosH(etaF));  
	    
	    /////////ML!!!	    
	    //  if(TMath::Abs(etaF)>5.)continue;

	    if(TMath::Abs(etaF)>5)goto M1;
	       
	    
     M2:    Double_t rho = TMath::Sqrt(gRandom->Rndm());
	    Double_t phi = TMath::TwoPi() * gRandom->Rndm();
	    Double_t Rx =  TMath::Sqrt(1-Epsilon)*RB; 
	    Double_t Ry =  TMath::Sqrt(1+Epsilon)*RB;
	    
	    x0 = Rx * rho * TMath::Cos(phi);
	    y0 = Ry * rho * TMath::Sin(phi);
	    r = TMath::Sqrt(x0*x0+y0*y0);
	    phiF = TMath::Abs(TMath::ATan(y0/x0));
	    
	    if(x0<0&&y0>0)phiF = TMath::Pi()-phiF;
	    if(x0<0&&y0<0)phiF = TMath::Pi()+phiF;
	    if(x0>0&&y0<0)phiF = 2.*TMath::Pi()-phiF;

	    //new Nov2012 AS-ML	    
	    // if(r>RB*(1+e3*TMath::Cos(3*(phiF+psiforv3)))/(1+e3))continue;

	    ///////////////  
            if(r>=RB*(1+e3*TMath::Cos(3*phiF+3*psiforv3))/(1+TMath::Abs(e3)))goto M2;
	    
	    //proper time with emission duration                                                               
	    Double_t tau = coeff_R1 * fParams.fTau +  sqrt(2.) * fParams.fSigmaTau * coeff_R1 * (gRandom->Gaus());                                            	  
	    z0 = tau  * TMath::SinH(etaF);                                                                             
	    t0 = tau  * TMath::CosH(etaF);
	    
	    Double_t rhou = fParams.fUmax * r / RB;          
	    //Double_t rhou0 = rhou;  

 	    
	    //Double_t rhou3 = 0.063*TMath::Sqrt((0.5*impactParameter)/0.67);
	    //Double_t rhou4 = 0.023*((0.5*impactParameter)/0.67);

	   Double_t rhou2 = fParams.frhou2;
           Double_t rhou3 = fParams.frhou3*TMath::Sqrt((0.5*impactParameter)/0.67);
           //Double_t rhou4 = fParams.frhou4*((0.5*impactParameter)/0.67);

//           Double_t rhou4 = 0.023*TMath::Sqrt((0.5*impactParameter)/0.67);

//           Double_t rrcoeff = TMath::Sqrt(1+Epsilon*TMath::Cos(2*phiF))/sqrt(1-Epsilon*Epsilon);
//           Double_t rrcoeff = 1.;

           Double_t rrcoeff = 1./TMath::Sqrt(1. + Delta*TMath::Cos(2*phiF));



	   //last s rho:
	   //	   rhou = rhou * (1 + rrcoeff*rhou3*TMath::Cos(3*(phiF+psiforv3)) + rrcoeff*rhou4*TMath::Cos(4*phiF) );

	   //last 22.06.16  
           rhou = rhou * (1 + rhou2*TMath::Cos(2*phiF) + rrcoeff*rhou3*TMath::Cos(3*(phiF+psiforv3)));


          Double_t delta1 = 0.;
          Delta = Delta * (1.0 + delta1 * TMath::Cos(phiF) - delta1 * TMath::Cos(3*phiF));

	    Double_t uxf = TMath::SinH(rhou)*TMath::Sqrt(1+Delta)*TMath::Cos(phiF); 
	    Double_t uyf = TMath::SinH(rhou)*TMath::Sqrt(1-Delta)*TMath::Sin(phiF);
	    Double_t utf = TMath::CosH(etaF) * TMath::CosH(rhou) * 
	      TMath::Sqrt(1+Delta*TMath::Cos(2*phiF)*TMath::TanH(rhou)*TMath::TanH(rhou));
	    Double_t uzf = TMath::SinH(etaF) * TMath::CosH(rhou) * 
	      TMath::Sqrt(1+Delta*TMath::Cos(2*phiF)*TMath::TanH(rhou)*TMath::TanH(rhou));
	    
	    vec3.SetXYZ(uxf / utf, uyf / utf, uzf / utf);
	    n1.Boost(-vec3); 
	    
	    yy = weightMax * gRandom->Rndm();        
               
	    Double_t php0 = TMath::TwoPi() * gRandom->Rndm();
	    Double_t ctp0 = 2. * gRandom->Rndm() - 1.;
	    Double_t stp0 = TMath::Sqrt(1. - ctp0 * ctp0); 
	    e = mass + (eMax - mass) * arrayFunctDistE(); 
	    Double_t pp0 = TMath::Sqrt(e * e - mass * mass);
	    px0 = pp0 * stp0 * TMath::Sin(php0); 
	    py0 = pp0 * stp0 * TMath::Cos(php0);
	    pz0 = pp0 * ctp0;
	    p0.SetXYZT(px0, py0, pz0, e);
	    
	    //weight for rdr          
	    weight = (n1 * p0) /e;  // weight for rdr gammar: weight = (n1 * p0) / n1[3] / e; 
	                                                           
	  } while(yy >= weight); 
          
	  //	  if(abs(z0)>1000 || abs(x0)>1000)
	  //std::cout<<"====== etaF==== "<<etaF<<std::endl;
          
	  partMom.SetXYZT(px0, py0, pz0, e);
	  partPos.SetXYZT(x0, y0, z0, t0);
	  partMom.Boost(vec3);
	  
	  Int_t type =0; //hydro
	  Particle particle(partDef, partPos, partMom, 0., 0, type, -1, zeroVec, zeroVec);
	  particle.SetIndex();
	  allocator.AddParticle(particle, source);
	} //nhsel==4 , no hydro part
      }
    }
  }
}

Bool_t InitialStateHydjet::ReadParams() {     
  Float_t par[200] = {0.};
  Int_t i = 0; 
  std::string s(40,' '); 
  std::ifstream input("RunInputHydjet");
  if (!input) {
    Error("Ukm::ReadParams", "Cannot open RunInputHydjet");
    return kFALSE;
  }
      
  while (std::getline(input, s)) {
    input>>par[i];
    if (i < 140) 
      std::cout<<s<<"     =  "<<par[i]<<std::endl;
    ++i;
    std::getline(input,s);
  }

  std::cout<<"\nFor output use the files RunOutput.root  \n\n"<< std::endl; 
   
  fParams.fNevnt  = Int_t(par[0]); //number of events
  fParams.fSqrtS  = par[1];        //cms energy per nucleon
  fParams.fAw     = par[2];        // atomic number of colliding nuclei
  fParams.fIfb    = Int_t(par[3]);      // flag of type of centrality generation (=0 is fixed by fBfix, not 0 
                                   //impact parameter is generated in each event between fBfmin 
                                   //and fBmax according with Glauber model (f-la 30)
   fParams.fBmin = par[4];         //minimum impact parameter in units of nuclear radius RA 
   fParams.fBmax = par[5];         //maximum impact parameter in units of nuclear radius RA
   fParams.fBfix = par[6];         //fix impact parameter in units of nuclear radius RA

   fParams.fSeed = Int_t(par[7]);         //parameter to set the random nuber seed (=0 the current time is used
                                   //to set the random generator seed, !=0 the value fSeed is 
                                   //used to set the random generator seed and then the state of random
                                   //number generator in PYTHIA MRPY(1)=fSeed
       
   fParams.fT         = par[8];     //chemical freeze-out temperature in GeV    
   fParams.fMuB       = par[9];     //baryon potential 
   fParams.fMuS       = par[10];    //strangeness potential 
   fParams.fMuC       = par[11];    //charm potential 
   fParams.fMuI3      = par[12];    //isospin potential   
   fParams.fThFO      = par[13];    //thermal freeze-out temperature T^th in GeV
   fParams.fMu_th_pip = par[14];    // effective chemical potential of positivly charged pions at thermal in GeV 

       
   fParams.fTau       = par[15];     //proper time value
   fParams.fSigmaTau  = par[16];     //its standart deviation (emission duration)
   fParams.fR         = par[17];     //maximal transverse radius 
   fParams.fYlmax     = par[18];     //maximal longitudinal rapidity 
   fParams.fUmax      = par[19];     //maximal transverse velocity multiplaed on \gamma_r 

   fParams.frhou2= par[20]; //parameter to swich ON/OFF(0) rhou2 
   fParams.frhou3= par[21];  //parameter to swich ON/OFF(0) rhou3  
   fParams.frhou4= par[22];  //parameter to swich ON/OFF(0) rhou4  


   fParams.fDelta     = par[23];     //momentum asymmetry parameter
   fParams.fEpsilon   = par[24];     //coordinate asymmetry parameter

   fParams.fv2= par[25]; //parameter to swich ON/OFF(0) epsilon2 fluctuations 
   fParams.fv3= par[26];  //parameter to swich ON/OFF(0) epsilon3 fluctuations 

   
   fParams.fIfDeltaEpsilon= Int_t(par[27]); // Flag to specify fDelta and fEpsilon values 
  
   fParams.fDecay      = Int_t(par[28]);    // flag to switch on/off hadron decays =0: decays off, >=1: decays on, (default: 1)
   fParams.fWeakDecay  = par[29];         //flag to set a decay width threshold; particles with lower decay widths will not be decayed; 
                                          //can be used e.g. to switch on/off weak hadron decays
  
   fParams.fEtaType   = Int_t(par[30]);     // flag to choose rapidity distribution, if fEtaType<=0, 
                                     //then uniform rapidity distribution in [-fYlmax,fYlmax] if fEtaType>0,
                                     //then Gaussian with dispertion = fYlmax 
  
   fParams.fTMuType   = Int_t(par[31]);     // flag to use calculated chemical freeze-out temperature,
                                     //baryon potential and strangeness potential as a function of fSqrtS 

   fParams.fCorrS     = par[32];            // flag and value to include strangeness supression factor    
   fParams.fCharmProd = Int_t(par[33]);     // flag to include charm production
   fParams.fCorrC     = par[34];            // charm supression factor

   fParams.fNhsel = Int_t(par[35]);         //flag to switch on/off jet and hydro-state production (0: jet
                                     // production off and hydro on, 1: jet production on and jet quenching
                                     // off and hydro on, 2: jet production on and jet quenching on and
                                     // hydro on, 3: jet production on and jet quenching off and hydro
                                     // off, 4: jet production on and jet quenching on and hydro off

   fParams.fPyhist = Int_t(par[36]);  // suppress PYTHIA partonic history (=1 only final state particle; =0 include particle history)
   fParams.fIshad= Int_t(par[37]);         //flag to switch on/off impact parameter dependent nuclear
                                    // shadowing for gluons and light sea quarks (u,d,s) (0: shadowing off,
                                    // 1: shadowing on for fAw=207, 197, 110, 40, default: 1
  
   fParams.fPtmin = par[38];       //minimal transverse momentum transfer p_T of hard
                                   // parton-parton scatterings in GeV (the PYTHIA parameter ckin(3)=fPtmin)
   
//  PYQUEN energy loss model parameters:
 
   fParams.fT0 = par[39];          // initial temperature (in GeV) of QGP for
                                   //central Pb+Pb collisions at mid-rapidity (initial temperature for other
                                  //centralities and atomic numbers will be calculated automatically) (allowed range is 0.2<fT0<2) 
  
   fParams.fTau0= par[40];        //proper QGP formation time in fm/c (0.01<fTau0<10)
   fParams.fNf= Int_t(par[41]);          //number of active quark flavours N_f in QGP fNf=0, 1,2 or 3 
   fParams.fIenglu= Int_t(par[42]);      // flag to fix type of in-medium partonic energy loss 
                                  //(0: radiative and collisional loss, 1: radiative loss only, 2:
                                  //collisional loss only) (default: 0);
   fParams.fIanglu= Int_t(par[43]);      //flag to fix type of angular distribution of in-medium emitted
                                  // gluons (0: small-angular, 1: wide-angular, 2:collinear) (default: 0).

  


   //PYTHIA parameters:
   Int_t jj;
   for (Int_t j = 0; j <25; ++j) {
     jj= 45+j;
     SERVICE.parPYTH[j]=par[jj];
   } 

   // Set Random Number seed 
        
   gRandom->SetSeed(fParams.fSeed); //Set 0 to use the current time
//to send seed in PYTHIA
   SERVICE.iseed_fromC=gRandom->GetSeed(); 
   std::cout<<"Seed for random number generation= "<<gRandom->GetSeed()<<std::endl;  

   fParams.fNPartTypes = 0;         //counter of hadron species

  
   return kTRUE; 
}



Bool_t InitialStateHydjet::MultIni() {
  //check and redefine input parameters
  if(fParams.fTMuType>0 &&  fParams.fSqrtS > 2.24) {
    if(fParams.fSqrtS < 2.24){
      Error("InitialStateHydjet::MultIni", "SqrtS<2.24 not allowed with fParams.fTMuType>0");
      return 0;
    }
    
    //sqrt(s) = 2.24 ==> T_kin = 0.8 GeV
    //see J. Cleymans, H. Oeschler, K. Redlich,S. Wheaton, Phys Rev. C73 034905 (2006)
    fParams.fMuB = 1.308/(1. + fParams.fSqrtS*0.273);
    fParams.fT = 0.166 - 0.139*fParams.fMuB*fParams.fMuB - 0.053*fParams.fMuB*fParams.fMuB*
      fParams.fMuB*fParams.fMuB;
    fParams.fMuI3 = 0.;
    fParams.fMuS = 0.;
    //create strange potential object and set strangeness density 0
    NAStrangePotential* psp = new NAStrangePotential(0., fDatabase);
    psp->SetBaryonPotential(fParams.fMuB);
    psp->SetTemperature(fParams.fT);
    //compute strangeness potential
    if(fParams.fMuB > 0.01)
      fParams.fMuS = psp->CalculateStrangePotential();
    cout << "fMuS = " << fParams.fMuS << endl;  

    //if user choose fYlmax larger then allowed by kinematics at the specified beam energy sqrt(s)     
    if(fParams.fYlmax > TMath::Log(fParams.fSqrtS/0.94)){
      Error("InitialStateHydjet::MultIni", "fParams.fYlmax > TMath::Log(fParams.fSqrtS/0.94)!!! ");
      return 0;
    }
      
    if(fParams.fCorrS <= 0.) {
      //see F. Becattini, J. Mannien, M. Gazdzicki, Phys Rev. C73 044905 (2006)
      fParams.fCorrS = 1. - 0.386* TMath::Exp(-1.23*fParams.fT/fParams.fMuB);
      std::
      cout<<"The phenomenological f-la F. Becattini et al. PRC73 044905 (2006) for CorrS was used." << std::endl;
      std::cout<<"Strangeness suppression parameter = "<<fParams.fCorrS << std::endl;
    }
    std::cout<<"The phenomenological f-la J. Cleymans et al. PRC73 034905 (2006) for Tch mu_B was used." << std::endl;
    std::cout<<"The simulation will be done with the calculated parameters:" << std::endl;
    std::cout<<"Baryon chemical potential = "<<fParams.fMuB<< " [GeV]" << std::endl;
    std::cout<<"Strangeness chemical potential = "<<fParams.fMuS<< " [GeV]" << std::endl;
    std::cout<<"Isospin chemical potential = "<<fParams.fMuI3<< " [GeV]" << std::endl;
    std::cout<<"Strangeness suppression parameter = "<<fParams.fCorrS << std::endl;
    std::cout<<"Eta_max = "<<fParams.fYlmax<<  std::endl;
    std::cout << std::endl;
  }
  
  std::cout<<"Used eta_max = "<<fParams.fYlmax<<  std::endl;
  std::cout<<"maximal allowed eta_max TMath::Log(fParams.fSqrtS/0.94)=  "<<TMath::Log(fParams.fSqrtS/0.94)<<std::endl;

  //initialisation of high-pt part 
                                                                                                     
  HYJPAR.nhsel = fParams.fNhsel;
  HYJPAR.ptmin = fParams.fPtmin;
  HYJPAR.ishad = fParams.fIshad;
  HYJPAR.iPyhist = fParams.fPyhist;
  HYIPAR.bminh = fParams.fBmin;
  HYIPAR.bmaxh = fParams.fBmax;
  HYIPAR.AW = fParams.fAw;
  
  HYPYIN.ifb = fParams.fIfb;
  HYPYIN.bfix = fParams.fBfix;
  HYPYIN.ene = fParams.fSqrtS;

  PYQPAR.T0 = fParams.fT0;
  PYQPAR.tau0 = fParams.fTau0;
  PYQPAR.nf = fParams.fNf;
  PYQPAR.ienglu = fParams.fIenglu;
  PYQPAR.ianglu = fParams.fIanglu;
  
  myini_();
 

  // calculation of  multiplicities of different particle species
  // according to the grand canonical approach
  GrandCanonical gc(15, fParams.fT, fParams.fMuB, fParams.fMuS, fParams.fMuI3, fParams.fMuC);
  GrandCanonical gc_ch(15, fParams.fT, fParams.fMuB, fParams.fMuS, fParams.fMuI3, fParams.fMuC);
  GrandCanonical gc_pi_th(15, fParams.fThFO, 0., 0., fParams.fMu_th_pip, fParams.fMuC);
  GrandCanonical gc_th_0(15, fParams.fThFO, 0., 0., 0., 0.);
  
   // std::ofstream outMult("densities.txt");
//    outMult<<"encoding    particle density      chemical potential "<<std::endl;


   Double_t Nocth=0; //open charm
   Double_t NJPsith=0; //JPsi
 
    
  //effective volume for central     
  double dYl= 2 * fParams.fYlmax; //uniform distr. [-Ylmax; Ylmax]  
  if (fParams.fEtaType >0) dYl = TMath::Sqrt(2 * TMath::Pi()) * fParams.fYlmax ;  //Gaussian distr.                                                                            
  fVolEff = 2 * TMath::Pi() * fParams.fTau * dYl * (fParams.fR * fParams.fR)/TMath::Power((fParams.fUmax),2) * 
  ((fParams.fUmax)*TMath::SinH((fParams.fUmax))-TMath::CosH((fParams.fUmax))+ 1);
  std::cout <<"central Effective volume = " << fVolEff << " [fm^3]" << std::endl;
  
  Double_t particleDensity_pi_ch=0;
  Double_t particleDensity_pi_th=0;
  //  Double_t particleDensity_th_0=0;

  if(fParams.fThFO != fParams.fT && fParams.fThFO > 0){
    GrandCanonical gc_ch(15, fParams.fT, fParams.fMuB, fParams.fMuS, fParams.fMuI3, fParams.fMuC);
    GrandCanonical gc_pi_th(15, fParams.fThFO, 0., 0., fParams.fMu_th_pip, fParams.fMuC);
    GrandCanonical gc_th_0(15, fParams.fThFO, 0., 0., 0., 0.);
    particleDensity_pi_ch = gc_ch.ParticleNumberDensity(fDatabase->GetPDGParticle(211));
    particleDensity_pi_th = gc_pi_th.ParticleNumberDensity(fDatabase->GetPDGParticle(211));
  }

  for(Int_t particleIndex = 0; particleIndex < fDatabase->GetNParticles(); particleIndex++) {
    ParticlePDG *currParticle = fDatabase->GetPDGParticleByIndex(particleIndex);
    Int_t encoding = currParticle->GetPDG();

    //strangeness supression
    Double_t gammaS = 1;
    Int_t S = Int_t(currParticle->GetStrangeness());
    if(encoding == 333)S = 2;
    if(fParams.fCorrS < 1. && S != 0)gammaS = TMath::Power(fParams.fCorrS,-TMath::Abs(S));


    //average densities      
    Double_t particleDensity = gc.ParticleNumberDensity(currParticle)/gammaS;
    
    //compute chemical potential for single f.o. mu==mu_ch
    Double_t mu = fParams.fMuB  * Int_t(currParticle->GetBaryonNumber()) + 
      fParams.fMuS  * Int_t(currParticle->GetStrangeness()) +
      fParams.fMuI3 * Int_t(currParticle->GetElectricCharge()) +
      fParams.fMuC * Int_t(currParticle->GetCharmness());

    //thermal f.o.
    if(fParams.fThFO != fParams.fT && fParams.fThFO > 0){
      Double_t particleDensity_ch = gc_ch.ParticleNumberDensity(currParticle);
      Double_t particleDensity_th_0 = gc_th_0.ParticleNumberDensity(currParticle);
      Double_t numb_dens_bolt = particleDensity_pi_th*particleDensity_ch/particleDensity_pi_ch;               
      mu = fParams.fThFO*TMath::Log(numb_dens_bolt/particleDensity_th_0);
      if(abs(encoding)==211 || encoding==111)mu= fParams.fMu_th_pip; 
      particleDensity = numb_dens_bolt;         
    }
    
    // set particle densities to zero for some particle codes
    // pythia quark codes
    if(abs(encoding)<=9) {
      particleDensity=0;        
    }
    // leptons
    if(abs(encoding)>10 && abs(encoding)<19) {
      particleDensity=0;   
    }
    // exchange bosons
    if(abs(encoding)>20 && abs(encoding)<30) {
      particleDensity=0;
    }
    // pythia special codes (e.g. strings, clusters ...)
    if(abs(encoding)>80 && abs(encoding)<100) {
      particleDensity=0;
    }
    // pythia di-quark codes
    // Note: in PYTHIA all diquark codes have the tens digits equal to zero
    if(abs(encoding)>1000 && abs(encoding)<6000) {
      Int_t tens = ((abs(encoding)-(abs(encoding)%10))/10)%10;
      if(tens==0) {             // its a diquark;
	particleDensity=0;
      }
    }
    // K0S and K0L
    if(abs(encoding)==130 || abs(encoding)==310) {
      particleDensity=0;
    }
    // charmed particles
     

    if(encoding==443)NJPsith=particleDensity*fVolEff/dYl;
 	
    
    // We generate thermo-statistically only J/psi(443), D_+(411), D_-(-411), D_0(421), 
    //Dbar_0(-421), D1_+(413), D1_-(-413), D1_0(423), D1bar_0(-423)
    //Dcs(431) Lambdac(4122)
    if(currParticle->GetCharmQNumber()!=0 || currParticle->GetCharmAQNumber()!=0) {
//ml if(abs(encoding)!=443 &&  
//ml	 abs(encoding)!=411 && abs(encoding)!=421 && 
//ml	 abs(encoding)!=413 && abs(encoding)!=423 && abs(encoding)!=4122 && abs(encoding)!=431) {
//ml	particleDensity=0; } 

     if(abs(encoding)==441 ||  
	 abs(encoding)==10441 || abs(encoding)==10443 || 
	 abs(encoding)==20443 || abs(encoding)==445 || abs(encoding)==4232 || abs(encoding)==4322 ||
	 abs(encoding)==4132 || abs(encoding)==4312 || abs(encoding)==4324 || abs(encoding)==4314 ||
	 abs(encoding)==4332 || abs(encoding)==4334 
	 ) {
	particleDensity=0; } 
	else
      {
      if(abs(encoding)!=443){ //only open charm
      Nocth=Nocth+particleDensity*fVolEff/dYl;       
//      cout<<encoding<<" Nochth "<<Nocth<<endl;
//      particleDensity=particleDensity*fParams.fCorrC;
//      if(abs(encoding)==443)particleDensity=particleDensity*fParams.fCorrC;
       }
      }
      
    }
    

    // bottom mesons
    if((abs(encoding)>500 && abs(encoding)<600) ||
       (abs(encoding)>10500 && abs(encoding)<10600) ||
       (abs(encoding)>20500 && abs(encoding)<20600) ||
       (abs(encoding)>100500 && abs(encoding)<100600)) {
      particleDensity=0;
    }
    // bottom baryons
    if(abs(encoding)>5000 && abs(encoding)<6000) {
      particleDensity=0;
    }
    ////////////////////////////////////////////////////////////////////////////////////////


    if(particleDensity > 0.) {
      fParams.fPartEnc[fParams.fNPartTypes] = encoding;
      fParams.fPartMult[2 * fParams.fNPartTypes] = particleDensity;
      fParams.fPartMu[2 * fParams.fNPartTypes] = mu;
      ++fParams.fNPartTypes;
      if(fParams.fNPartTypes > 1000)
	Error("in Bool_t MultIni:", "fNPartTypes is too large %d", fParams.fNPartTypes);

    //outMult<<encoding<<" "<<particleDensity*fVolEff/dYl <<" "<<mu<<std::endl;

    }
  }

//put open charm number and cc number in Params
     fParams.fNocth = Nocth;
     fParams.fNccth = NJPsith;


  return kTRUE;
}


Double_t InitialStateHydjet::SimpsonIntegrator2(Double_t a, Double_t b, Double_t Epsilon, Double_t Delta) {

 // std::cout<<"in SimpsonIntegrator2: epsilon"<<Epsilon<<"delta"<<Delta<<std::endl; 
  Int_t nsubIntervals=10000;
  Double_t h = (b - a)/nsubIntervals; //-1-pi, phi
  Double_t s=0;
//  Double_t h2 = (fParams.fR)/nsubIntervals; //0-R maximal RB ?

  Double_t x = 0; //phi
  for(Int_t j = 1; j < nsubIntervals; j++) {
    x += h; // phi
    Double_t e = Epsilon;
    Double_t RsB = fParams.fR; //test: podstavit' *coefff_RB
    Double_t RB = RsB *(TMath::Sqrt(1-e*e)/TMath::Sqrt(1+e*TMath::Cos(2*x))); //f-la7 RB    
    Double_t sr = SimpsonIntegrator(0,RB,x,Delta);
    s += sr;
  }
  return s*h;

}

Double_t InitialStateHydjet::SimpsonIntegrator(Double_t a, Double_t b, Double_t phi, Double_t Delta) {
//  std::cout<<"in SimpsonIntegrator"<<"delta"<<Delta<<std::endl; 
  Int_t nsubIntervals=100;
  Double_t h = (b - a)/nsubIntervals;
  Double_t s = f2(phi,a + 0.5*h,Delta);
  Double_t t = 0.5*(f2(phi,a,Delta) + f2(phi,b,Delta));
  Double_t x = a;
  Double_t y = a + 0.5*h;
  for(Int_t i = 1; i < nsubIntervals; i++) {
    x += h;
    y += h;
    s += f2(phi,y,Delta);
    t += f2(phi,x,Delta);
  }	
  t += 2.0*s;
  return t*h/3.0;
}


//f2=f(phi,r)
Double_t InitialStateHydjet::f2(Double_t x, Double_t y, Double_t Delta) {
  //std::cout<<"in f2: "<<"delta"<<Delta<<std::endl; 
  Double_t RsB = fParams.fR; //test: podstavit' *coefff_RB
  Double_t rhou =  fParams.fUmax * y / RsB;
  Double_t ff = y*TMath::CosH(rhou)*
    TMath::Sqrt(1+Delta*TMath::Cos(2*x)*TMath::TanH(rhou)*TMath::TanH(rhou));
//n_mu u^mu f-la 20
  return ff;
}


Double_t InitialStateHydjet::MidpointIntegrator2(Double_t a, Double_t b, Double_t Delta, Double_t Epsilon) {

  Int_t nsubIntervals=2000; 
  Int_t nsubIntervals2=1; 
  Double_t h = (b - a)/nsubIntervals; //0-pi , phi
  Double_t h2 = (fParams.fR)/nsubIntervals; //0-R maximal RB ?

  Double_t x = a + 0.5*h;
  Double_t y = 0;
      
  Double_t t = f2(x,y,Delta);                    
 
  Double_t e = Epsilon;

  for(Int_t j = 1; j < nsubIntervals; j++) {
    x += h; // integr  phi

    Double_t RsB = fParams.fR; //test: podstavit' *coefff_RB
    Double_t  RB = RsB *(TMath::Sqrt(1-e*e)/TMath::Sqrt(1+e*TMath::Cos(2*x))); //f-la7 RB

    nsubIntervals2 = Int_t(RB / h2)+1;
    // integr R 
    y=0;
    for(Int_t i = 1; i < nsubIntervals2; i++) 
      t += f2(x,(y += h2),Delta);
  }
  return t*h*h2;
}

Double_t InitialStateHydjet::CharmEnhancementFactor(Double_t Ncc, Double_t Ndth, Double_t NJPsith, Double_t Epsilon) {
  
  Double_t gammaC=100.;
  Double_t x1 = gammaC*Ndth; 
  Double_t var1 = Ncc-0.5*gammaC*Ndth*TMath::BesselI1(x1)/TMath::BesselI0(x1)-gammaC*gammaC*NJPsith;
 // cout<<"gammaC 20"<<" var "<<var1<<endl;
  gammaC=1.;
  Double_t x0 = gammaC*Ndth; 
  Double_t var0 = Ncc-0.5*gammaC*Ndth*TMath::BesselI1(x0)/TMath::BesselI0(x0)-gammaC*gammaC*NJPsith;
 // cout<<"gammaC 1"<<" var "<<var0<<endl;
  
  for(Int_t i=1; i<1000; i++){ 
    if(var1 * var0<0){
     gammaC=gammaC+0.01*i;
     Double_t x = gammaC*Ndth;  
     var0 = Ncc-0.5*gammaC*Ndth*TMath::BesselI1(x)/TMath::BesselI0(x)-gammaC*gammaC*NJPsith;
      }
      else
      {
//      cout<<"gammaC "<<gammaC<<" var0 "<<var0<<endl;
      return gammaC;
      } 

    }
//      cout<<"gammaC not found ? "<<gammaC<<" var0 "<<var0<<endl;
  return -100;
  }

  
  
