#ifndef MOLECULE_H
#define MOLECULE_H

#include <set>
#include <bond.hpp>

class Bond;

class Molecule
{
   private:
      std::set<Bond> bondList;
   public:
      void add (Bond newBond) {bondList.insert(newBond);}
};

#endif
