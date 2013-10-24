#include <atom.hpp>
#include <bond.hpp>

Bond::Bond (const Atom* c, const Atom* d)
{
   if (c > d)
   {
      a = d;
      b = c;
   }
   else
   {
      a = c;
      b = d;
   }
}

void Bond::copy (const Bond& target)
{
   a = target.a;
   b = target.b;
}

bool Bond::operator== (const Bond & other) const 
{
   return (a == other.a && b == other.b) || (a == other.b && b == other.a);
}

bool Bond::operator< (const Bond & other) const
{
   return a < other.a || ( b < other.b && a == other.a );
}

int Bond::operator[] (unsigned int i) const 
{
   if (i==0) return a->getIndex();
   return b->getIndex(); 
}

bool Bond::offCandidate() const
{
   return (a->getNumBonds() > MIN_BOND && b->getNumBonds() > MIN_BOND);
}
