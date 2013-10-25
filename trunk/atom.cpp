#include <atom.hpp>
#include <point.hpp>
#include <matrix3.hpp>
#include <minexcept.hpp>
#include <iostream>
#include <fstream>

#include <cstdlib>
#include <ctime>

#include <vector>
#include <string>

#define debug(x) std::cout << __LINE__ << ' ' << x << std::endl; std::cout.flush();

using namespace std;

Atom::Atom ()
{
   pos.setAll(0);
}

Atom::Atom(const Point& target, int index)
{
   setPos(target);
   atomIndex = index;
}
      
void Atom::copy (const Atom& target)
{  
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
   pos = coord.changeBasis(cellInfo.dim());
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

void Atom::setParam (Parameters info)
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

bool Atom::delRandBond()
{
   Bond candidate(*(cellInfo.moffCandidates.begin()));
   int i = 1;
   bool success = false;
   for (set<Bond>::iterator it=cellInfo.moffCandidates.begin(); it!=cellInfo.moffCandidates.end();)
   {
      if (!it->offCandidate())
      {
         cellInfo.moffCandidates.erase(*it++);
         continue;
      }
      if (rand() % i == 0)
      {
         candidate = *it;
         success = true;
      }
      it++;
      i++;
   }
   if (success)
   {
      cellInfo.moffCandidates.erase(candidate);
      bondOff(candidate);
   }
   return success;
}

int Atom::getNumBonds() const
{
   int numBond = 0;
   for (int i=0; i<adjMatrix.columns(); i++)
   {
      if (adjMatrix.get(atomIndex,i) == 1)
         numBond += 1;
   }
   return numBond;
}

double Atom::getPos(char dimen) const
{
   Point posXYZ = pos.xyzBasis(cellInfo.dim());
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
   Point posXYZ = pos.xyzBasis(cellInfo.dim());
   return posXYZ.coord(i);
}

Point Atom::getPos() const
{
   return pos.xyzBasis(cellInfo.dim());
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

void Atom::readAtoms(string fileName)
{
   string buffer;
   fstream file (fileName.c_str(), fstream::in | fstream::out);
   double x,y,z;
   int i,n;
   Point tmp;

   file >> buffer >> buffer;
   for (i=0;i<3;i++)
   {
      file >> x >> y >> z;
      tmp.x(x);
      tmp.y(y);
      tmp.z(z);
      cellInfo.setDim(tmp, i);
   }
   
   file >> buffer >> n >> buffer;
   cellInfo.pnt(n);
   for (i=0;i<cellInfo.pnt();i++)
   {
      file >> n >> x >> y >> z;
      tmp.x(x);
      tmp.y(y);
      tmp.z(z);
      Atom newAtom(tmp, i);
      atomList.push_back(newAtom);
   }
   file.close();
}

//Assumes atom connections will be generated afterwards
void Atom::multiplyCell(Point n)
{
   Parameters newInfo = cellInfo;
   Parameters oldInfo = cellInfo;
   std::vector<Atom> oldList = atomList;
   
   int i,j,k,m;
   unsigned int index;
   Matrix3 factor;
   factor.setDiag(n);

   newInfo.pnt( oldInfo.pnt() * (n.x() * n.y() * n.z()) );
   newInfo.setDim( oldInfo.dim()*factor );

   for (i=0;i<n.x();i++)
   {
      for (j=0;j<n.y();j++)
      {
         for (k=0;k<n.z();k++)
         {
            for (m=0;m<oldInfo.pnt();m++)
            {
               cellInfo = oldInfo;
               Point pos(oldList[m].getPos());
               cellInfo = newInfo;
               index = m + k*oldInfo.pnt() + j*n.z()*oldInfo.pnt() + i*n.y()*n.z()*oldInfo.pnt();
               pos += (oldInfo.dim(0)*i);
               pos += (oldInfo.dim(1)*j);
               pos += (oldInfo.dim(2)*k);
               if (index < atomList.size())
               {
                  atomList[index].setPos(pos);
                  atomList[index].setIndex(index);
                  atomList[index].clearNeighbours();
               }
               else
               {
                  Atom newAtom(pos, index);
                  atomList.push_back(newAtom);
               }
            }
         }
      }
   }
   cellInfo = newInfo;
}

void Atom::connectAtoms(int target)
{
   int size = atomList.size();
   adjMatrix = Matrix(size);
   Point a,b;
   int nBond = 0;
   cellInfo.cxn(0);
   double cellDist;
   double epsilon = 1e-8;
   
   cellInfo.dist( fabs( (atomList[0].getPos() - atomList[target].getPos()).distance() ) );
   
   for (int i=0; i<size; i++)
   {
      for (int j=0; j<size; j++)
      {
         if (i == j)
            continue;
         a = atomList[i].getPos();
         b = atomList[j].getPos();
         cellDist = fabs( (cellInfo.getRealDiff(a,b)).distance());

         if ( fabs(cellDist - cellInfo.dist()) < epsilon)
         {
            atomList[i].setNeighbour(&atomList[j]);
            adjMatrix.setSym(i,j,1);
            nBond++;
         }
      }
      if (nBond > cellInfo.cxn())
         cellInfo.cxn(nBond);
      nBond = 0;
   }
}

void Atom::outputAtoms(string fileName)
{
   Parameters p = cellInfo;
   int i;
   set<Bond>::const_iterator it;
   set<Angle>::const_iterator it2;
   fstream file ( (fileName.append(".out")).c_str(), fstream::in | fstream::out | fstream::trunc);
   
   file << fileName << endl << endl;

   file << cellInfo.pnt() << " atoms" << endl;
   file << cellInfo.nBonds() << " bonds" << endl;
   file << cellInfo.nAngles() << " angles" << endl;
   file << "0 dihedrals" << endl;
   file << "0 impropers" << endl << endl;

   file << "1 atom types" << endl;
   file << "1 bond types" << endl;
   if (cellInfo.angles().size() != 0) file << "1 angle types" << endl;
   else file << "0 angle types" << endl;
   file << "0 dihedral types" << endl;
   file << "0 improper types" << endl << endl;
   
   file << "0 " << p.dim(0).x() << " xlo xhi" << endl;
   file << "0 " << p.dim(1).y() << " ylo yhi" << endl;
   file << "0 " << p.dim(2).z() << " zlo zhi" << endl;
   file << cellInfo.dim(1).x() << ' ' << cellInfo.dim(2).x() << ' ' << cellInfo.dim(2).y() 
        << " xy xz yz" << endl << endl;

   file << "Atoms" << endl << endl;

   for (i=0;i<cellInfo.pnt();i++)
   {
      file << i+1 << " 1 1 " << atomList[i].getPos('x') << " " << atomList[i].getPos('y') << " " 
           << atomList[i].getPos('z') << endl;
   }
   
   file << endl << "Bonds" << endl << endl;
   
   i = 0;
   for (it=cellInfo.bonds().begin();it!=cellInfo.bonds().end();++it)
   {
      file << ++i << " 1 " << (*it)[0]+1 << ' ' << (*it)[1]+1 << endl;
   }

   if (cellInfo.angles().size() != 0)
   {
      file << endl << "Angles" << endl << endl;
      i = 0;
      for (it2=cellInfo.angles().begin();it2!=cellInfo.angles().end();++it2)
      {
         file << ++i << " 1 " << (*it2)[0]+1 << ' ' << (*it2)[1]+1 << ' ' << (*it2)[2]+1 << endl;
      }
   }

   file.close();
}

void Atom::genBondList()
{
   cellInfo.mbonds.clear();
   unsigned int i,j;

   for (i=0;i<atomList.size();i++)
   {
      for (j=i+1;j<atomList.size();j++)
      {
         if (adjMatrix.get(i,j) == 1) cellInfo.mbonds.insert(Bond(&atomList[i], &atomList[j]) );
      }
   }
}

void Atom::genAngleList()
{
   cellInfo.mangles.clear();
   unsigned int i,j,k;

   for (i=0;i<atomList.size();i++)
   {
      for (j=0;j<atomList.size();j++)
      {
         for (k=j+1;k<atomList.size();k++)
         {
            if (adjMatrix.get(i,j) == 1 && adjMatrix.get(i,k) == 1)
            {
               try { cellInfo.mangles.insert( Angle( &atomList[i], &atomList[j], &atomList[k]) ); }
               catch (BadStructureException e) { continue; }
            }
         }
      }
   }
}

void Atom::genDelList()
{
   cellInfo.moffCandidates = cellInfo.mbonds;
}

void Atom::bondOff(Bond target)
{
   int i = target[0];
   int j = target[1];
   if (adjMatrix.get(i,j) != 0)
   {
      adjMatrix.setSym(i,j,-1);
   }

}

void Atom::bondOn(Bond target)
{
   int i = target[0];
   int j = target[1];
   if (adjMatrix.get(i,j) != 0)
   {
      adjMatrix.setSym(i,j,1);
   }
}
