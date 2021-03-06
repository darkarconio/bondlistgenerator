#include <atom.hpp>
#include <bond.hpp>
#include <fxn.hpp>
#include <iostream>

#define debug(x) std::cout << __LINE__ << ' ' << x << std::endl; std::cout.flush();

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

bool Bond::offCandidate(unsigned int guideDel) const
{
   return guideCriteria(guideDel);
}

bool Bond::offCandidate (unsigned int guideDel, double distPercent) const
{
   Point loc = location();
   Point mid = a->cellInfo.midpoint();
   double dist = distPercent/200*a->cellInfo.len().average();

   if (a->cellInfo.getRealDiff(loc, mid).distance() > dist && distPercent != 0)
      return (false);
   
   return guideCriteria(guideDel);
}


Point Bond::location() const
{
   Point locA = a->getPos();
   Point locB = b->getPos();
   return a->cellInfo.getRealAvg(locA, locB);
}

bool Bond::guideCriteria(unsigned int guideDel) const
{
   if (guideDel == 0)
   {
      return (a->getNumBonds() > MIN_BOND && b->getNumBonds() > MIN_BOND);
   }
   else if (guideDel == 1)
   {
      return (a->getNumBonds() > MIN_NEIGH_BOND && b->getNumBonds() > MIN_BOND) ||
             (b->getNumBonds() > MIN_NEIGH_BOND && a->getNumBonds() > MIN_BOND);
   }
   else if (guideDel == 2)
   {
      return (a->getNumBonds() > MIN_NEIGH_BOND && b->getNumBonds() > MIN_NEIGH_BOND);
   }
   else if (guideDel == 3)
   {
      if (Atom::getNumCoordX(FULL_BOND) > Atom::getNumCoordX(MIN_NEIGH_BOND))
         return (a->getNumBonds() > MIN_NEIGH_BOND && b->getNumBonds() > MIN_NEIGH_BOND);
      else
         return (a->getNumBonds() > MIN_BOND && b->getNumBonds() > MIN_BOND);
   }
   else if (guideDel == 4)
   {
      if (Atom::getNumCoordX(FULL_BOND) > Atom::getNumCoordX(MIN_NEIGH_BOND)*5)
         return (a->getNumBonds() > MIN_NEIGH_BOND && b->getNumBonds() > MIN_NEIGH_BOND);
      else
         return (a->getNumBonds() > MIN_BOND && b->getNumBonds() > MIN_BOND);
   }
   else if (guideDel == 5)
   {
      if (Atom::getNumCoordX(FULL_BOND)*5 > Atom::getNumCoordX(MIN_NEIGH_BOND))
         return (a->getNumBonds() > MIN_NEIGH_BOND && b->getNumBonds() > MIN_NEIGH_BOND);
      else
         return (a->getNumBonds() > MIN_BOND && b->getNumBonds() > MIN_BOND);
   }
   else //if (guideDel == 6)
   {
      if (Atom::getNumCoordX(FULL_BOND) > Atom::getNumCoordX(MIN_NEIGH_BOND)*5)
         return (a->getNumBonds() > MIN_NEIGH_BOND && b->getNumBonds() > MIN_NEIGH_BOND);
      else
         return (a->getNumBonds() > MIN_NEIGH_BOND && b->getNumBonds() > MIN_BOND) ||
                (b->getNumBonds() > MIN_NEIGH_BOND && a->getNumBonds() > MIN_BOND);
   }
}
