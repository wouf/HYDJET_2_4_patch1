#ifndef INITIALSTATEHYDJET_INCLUDED
#define INITIALSTATEHYDJET_INCLUDED

/*                                                                            
                                                                            
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/

#include "DatabasePDG.h"
#include "Particle.h"
#include "InitialState.h"

struct InitialParamsHydjet_t {

  Int_t fNevnt;                      // number of events
  Double_t fSqrtS;                   // cms energy per nucleon
  Double_t fAw;                      // atomic number of colliding nuclei
  Int_t fIfb;                        // flag of type of centrality generation (=0 is fixed by fBfix, not 0 
                                     // impact parameter is generated in each event between fBfmin 
                                     // and fBmax according with Glauber model (f-la 30)
  Double_t fBmin;                    // minimum impact parameter in units of nuclear radius RA 
  Double_t fBmax;                    // maximum impact parameter in units of nuclear radius RA
  Double_t fBfix;                    // fix impact parameter in units of nuclear radius RA
  
  Int_t fSeed;                       // parameter to set the random nuber seed (=0 the current time is used
                                     // to set the random generator seed, !=0 the value fSeed is 
                                     // used to set the random generator seed and then the state of random
                                     // number generator in PYTHIA MRPY(1)=fSeed
       
  Double_t fT;                       // chemical freeze-out temperature in GeV    
  Double_t fMuB;                     // baryon potential 
  Double_t fMuS;                     // strangeness potential 
  Double_t fMuC;                     // charm potential 
  Double_t fMuI3;                    // isospin potential   
  Double_t fThFO;                    // thermal freeze-out temperature T^th in GeV
  Double_t fMu_th_pip;               // effective chemical potential of positivly charged pions at thermal in GeV 

       
  Double_t fTau;                     // proper time value
  Double_t fSigmaTau;                // its standart deviation (emission duration)
  Double_t fR;                       // maximal transverse radius 
  Double_t fYlmax;                   // maximal longitudinal rapidity 
  Double_t fUmax;                    // maximal transverse velocity multiplaed on \gamma_r 

  Double_t frhou2;                  //parameter to swich ON/OFF(0) rhou2 
  Double_t frhou3;                  //parameter to swich ON/OFF(0) rhou3  
  Double_t frhou4;                  //parameter to swich ON/OFF(0) rhou4  

  Double_t fDelta;                   // momentum asymmetry parameter
  Double_t fEpsilon;                 // coordinate asymmetry parameter

  Double_t  fv2;             //parameter to swich ON/OFF(0) epsilon2 fluctuations 
  Double_t  fv3;            //parameter to swich ON/OFF(0) epsilon3 fluctuations 

  Int_t    fIfDeltaEpsilon;          // flag to specify fDelta and fEpsilon values(=0 user's ones, >=1 parametrized)
  
  Int_t fDecay;                      // flag to switch on/off hadron decays (=0 decays off, >=1 decays on), (default: 1)
  Double_t fWeakDecay;               // flag to switch on/off weak hadron decays <0: decays off, >0: decays on, (default: 0)
  Int_t fPythDecay;                  // Flag to choose how to decay resonances in high-pt part, fPythDecay: 0 by PYTHIA decayer, 
                                     // 1 by FASTMC decayer(mstj(21)=0)  
  
  Int_t fEtaType;                    // flag to choose rapidity distribution, if fEtaType<=0, 
                                     // then uniform rapidity distribution in [-fYlmax,fYlmax] if fEtaType>0,
                                     // then Gaussian with dispertion = fYlmax 
  
  Int_t fTMuType;                    // flag to use calculated chemical freeze-out temperature,
                                     // baryon potential and strangeness potential as a function of fSqrtS 

  Double_t fCorrS;                   // flag and value to include strangeness supression factor    
  Int_t fCharmProd;                  // flag to include statistical charm production    
  Double_t fCorrC;                   // flag and value to include charmness supression factor

  Int_t fNhsel;                      // flag to switch on/off jet and hydro-state production (0: jet
                                     // production off and hydro on, 1: jet production on and jet quenching
                                     // off and hydro on, 2: jet production on and jet quenching on and
                                     // hydro on, 3: jet production on and jet quenching off and hydro
                                     // off, 4: jet production on and jet quenching on and hydro off
  Int_t fPyhist;                     // Suppress PYTHIA particle history (=1 only final state particles from hard part; =0 include full particle history) 
  Int_t fIshad;                      // flag to switch on/off impact parameter dependent nuclear
                                     // shadowing for gluons and light sea quarks (u,d,s) (0: shadowing off,
                                     // 1: shadowing on for fAw=207, 197, 110, 40, default: 1
  
  Double_t fPtmin;                   // minimal transverse momentum transfer p_T of hard
                                     // parton-parton scatterings in GeV (the PYTHIA parameter ckin(3)=fPtmin)
   
  //  PYQUEN energy loss model parameters:
 
  Double_t fT0;                      // initial temperature (in GeV) of QGP for
                                     // central Pb+Pb collisions at mid-rapidity (initial temperature for other
                                     // centralities and atomic numbers will be calculated automatically) (allowed range is 0.2<fT0<2) 
  
  Double_t fTau0;                    // proper QGP formation time in fm/c (0.01<fTau0<10)
  Int_t fNf;                         // number of active quark flavours N_f in QGP fNf=0, 1,2 or 3 
  Int_t fIenglu;                     // flag to fix type of in-medium partonic energy loss 
                                     // (0: radiative and collisional loss, 1: radiative loss only, 2:
                                     // collisional loss only) (default: 0);
  Int_t fIanglu;                     // flag to fix type of angular distribution of in-medium emitted
                                     // gluons (0: small-angular, 1: wide-angular, 2:collinear) (default: 0).
  
  
  Int_t    fNPartTypes;              //counter of hadron species  
  Int_t    fPartEnc[1000];           //Hadron encodings. Maximal number of hadron species is 1000!!!
  Double_t fPartMult[2000];          //Multiplicities of hadron species
  Double_t fPartMu[2000];            //Chemical potentials of hadron species

  Double_t fMuTh[1000];              //Chemical potentials at thermal freezeout of hadron species

//open charm
  
  Double_t fNocth;
  Double_t fNccth;
  

};

class InitialStateHydjet : public InitialState {
 private:
  Double_t fVolEff;                           // the effective volume
  InitialParamsHydjet_t fParams;             // the list of initial state parameters
  
 public:
  InitialStateHydjet() {};
  ~InitialStateHydjet() {};
  
  void SetVolEff(Double_t value) {fVolEff = value;}
  Double_t GetVolEff() {return fVolEff;}
  virtual Bool_t RunDecays() {return (fParams.fDecay>0 ? kTRUE : kFALSE);}
  virtual Int_t GetNev() {return fParams.fNevnt;}
  virtual Double_t GetWeakDecayLimit() {return fParams.fWeakDecay;}  
  
  virtual void Initialize(List_t &source, ParticleAllocator &allocator);
  virtual Bool_t ReadParams();
  virtual Bool_t MultIni();
  Bool_t IniOfThFreezeoutParameters();
 
  Double_t f(Double_t);
  Double_t f2(Double_t, Double_t, Double_t);
   
  Double_t SimpsonIntegrator(Double_t, Double_t, Double_t, Double_t);
  Double_t SimpsonIntegrator2(Double_t, Double_t, Double_t, Double_t);
  Double_t MidpointIntegrator2(Double_t, Double_t, Double_t, Double_t);
  
  Double_t CharmEnhancementFactor(Double_t, Double_t, Double_t, Double_t);

  
};

#endif
