#ifndef BOND_H
#define BOND_H

#include <atom.hpp>

class Atom;
class Angle;

class Bond
{
   private:
      Atom * a;
      Atom * b;
   protected:
      friend class Angle;
      Atom * otherAtom (Atom * c) const {if (c == a) return b; else return a;}
   public:
      Bond (Atom*, Atom*);
      bool operator== (const Bond&) const;
      bool operator< (const Bond&) const;
      int operator[] (unsigned int) const;
};

#endif
