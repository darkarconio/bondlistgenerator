#include <atom.hpp>
#include <point.hpp>
#include <cstdlib>
#include <ctime>
#include <matrix3.hpp>
#include <minexcept.hpp>
#include <iostream>
#include <vector>

#define debug(x) std::cout << __LINE__ << ' ' << x << std::endl; std::cout.flush();

Atom::Atom (Parameters *info)
{
   pos.setAll(0);
   cellInfo = info;
}

Atom::Atom(const Point& target, int index, Parameters* info)
{
   cellInfo = info;
   setPos(target);
   atomIndex = index;
}
      
void Atom::copy (const Atom& target)
{  
   cellInfo = target.cellInfo;
   neighbours = target.neighbours;
   pos = target.pos;
   atomIndex = target.atomIndex;
}

Atom::Atom (const Atom & target)
{
   copy (target);
}

void Atom::clearNeighbours()
{
   for (int i=neighbours.size()-1; i>=0; i--) delNeighbour(i);
}

void Atom::setPos(const Point & coord)
{
   pos = coord.changeBasis(cellInfo->dim());
}

void Atom::setNeighbour(Atom * pTarget)
{
   bool alreadyNeighbour = false;
   for(unsigned int i=0;i<neighbours.size();i++)
   {
      if (neighbours[i] == pTarget) alreadyNeighbour = true;
   }
   if (!alreadyNeighbour)
   {
      neighbours.push_back(pTarget);
      pTarget->neighbours.push_back(this);
   }
}

void Atom::setParam (Parameters *info)
{
   cellInfo = info;
}

Atom * Atom::getNeighbour(int num) const
{
   return neighbours[num];
}

void Atom::delNeighbour(int num)
{
   Atom * that = neighbours[num];
   for(unsigned int i = 0; i < that->neighbours.size(); i++)
   {
      if (that->neighbours[i] == this) that->neighbours.erase(that->neighbours.begin()+i);
   }
   
   neighbours.erase(neighbours.begin()+num);
}

void Atom::delRandNeighbour()
{
   srand ( time(NULL) );
   delNeighbour(rand() % neighbours.size());
}

double Atom::getPos(char dimen) const
{
   Point posXYZ = pos.xyzBasis(cellInfo->dim());
   switch (dimen)
   {
      case 'x':
         return posXYZ.x();
      case 'y':
         return posXYZ.y();
      case 'z':
         return posXYZ.z();
      default:
         throw BadIndex();
   }
}

double Atom::getPos(int i) const
{
   Point posXYZ = pos.xyzBasis(cellInfo->dim());
   return posXYZ.coord(i);
}

Point Atom::getPos() const
{
   return pos.xyzBasis(cellInfo->dim());
}

void Atom::setIndex(int num)
{
   atomIndex = num;
}

int Atom::getIndex() const
{
   return atomIndex;
}

int Atom::getNeighbourIndex(int num) const
{
   return neighbours[num]->atomIndex;
}
