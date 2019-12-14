#ifndef HADRONDECAYER_INCLUDED
#define HADRONDECAYER_INCLUDED

/*                                                                            
                                                                           
        Nikolai Amelin, Ludmila Malinina, Timur Pocheptsov (C) JINR/Dubna
      amelin@sunhe.jinr.ru, malinina@sunhe.jinr.ru, pocheptsov@sunhe.jinr.ru 
                           November. 2, 2005                                

*/
//#include "Rtypes.h"
#include "DatabasePDG.h"
#include "Particle.h"

Double_t GetDecayTime(const Particle &p, Double_t weakDecayLimit);
void Decay(List_t &output, Particle &p, ParticleAllocator &allocator, DatabasePDG* database);

#endif
