#ifndef BOND_H
#define BOND_H

#include <atom.hpp>

class Atom;

class Bond
{
   private:
      Atom * a;
      Atom * b;
   public:
      Bond (Atom*, Atom*);
      bool operator== (const Bond&) const;
      bool operator< (const Bond&) const;
      int index(int) const;
};

#endif
