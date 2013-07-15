#ifndef INTERACTION_H
#define INTERACTION_H

#include <set>
#include <bond.hpp>

class Bond;

class Interaction
{
   private:
      std::set<Bond> bondList;
      std::set<Bond>::const_iterator it;
   public:
      void add (Bond newBond) {bondList.insert(newBond); if(bondList.size()==1)it=bondList.begin();}
      const Bond get () {if(it==bondList.end())it=bondList.begin(); return *it++;}
      int size() const {return bondList.size();}
      void restartIter() {it = bondList.begin();}
};

#endif
