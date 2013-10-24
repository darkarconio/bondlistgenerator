#ifndef BOND_H
#define BOND_H

#include <atom.hpp>

class Atom;
class Angle;

class Bond
{
   private:
      const Atom * a;
      const Atom * b;
      static const int MIN_BOND = 2;
      void copy (const Bond&);
   protected:
      friend class Angle;
      const Atom * otherAtom (const Atom * c) const {if (c == a) return b; else return a;}
   public:
      Bond (const Atom*, const Atom*);
      Bond (const Bond& other) {copy(other);}
      bool operator== (const Bond&) const;
      bool operator< (const Bond&) const;
      int operator[] (unsigned int) const;
      Bond& operator= (const Bond& other) { copy(other); return *this; }
      bool offCandidate() const; 
};

#endif
