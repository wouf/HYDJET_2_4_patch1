/*
  Ludmila Malinina  malinina@lav01.sinp.msu.ru,   SINP MSU/Moscow and JINR/Dubna
  Ionut Arsene  i.c.arsene@fys.uio.no,            Oslo University and ISS-Bucharest
  Date        : 2007/05/30
*/

#ifndef INITIALSTATE
#define INITIALSTATE
// Virtual class for the initial state classes
// Include here common methods, but always declare them as virtual
#include "Particle.h"
#include "DatabasePDG.h"

class InitialState {
 protected:
   DatabasePDG *fDatabase;
 public:
  InitialState() {
    fDatabase = new DatabasePDG();
    fDatabase->LoadData();
    fDatabase->SetMassRange(0.0, 200.);
    fDatabase->SetWidthRange(0., 10.);
  };
  virtual ~InitialState() {
    delete fDatabase;
  };
  
  virtual void Initialize(List_t &source, ParticleAllocator &allocator) = 0;
  virtual Bool_t ReadParams() = 0;
  virtual Bool_t MultIni() = 0;
  virtual Bool_t RunDecays() = 0;
  virtual Int_t GetNev() = 0;
  virtual Double_t GetWeakDecayLimit() = 0;
  
  virtual void Evolve(List_t &secondaries, ParticleAllocator &allocator, Double_t weakDecayLimit);
};

#endif
