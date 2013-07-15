#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include <cmath>
#include <cstdlib>

#include <atom.hpp>
#include <point.hpp>
#include <fxn.hpp>
#include <bond.hpp>
#include <interaction.hpp>

#define debug(x) cout << __LINE__ << ' ' << x << endl; cout.flush();

using namespace std;

void readAtoms(vector<Atom>&, string, Parameters&);
void multiplyCell(vector<Atom>&, Point&, Parameters&, Parameters&);
void connectAtoms(vector<Atom>&, int, Parameters&);
void outputAtoms(vector<Atom>&, string, Parameters&);

int mod(int, int); 
double degRad(double);
double pi();

void changeAtomParam(vector<Atom>&, Parameters&);
//double energy (vector<Atom>&, Parameters&, double, double);

int main (int argc, char* argv[])
{
   string inputFile = argv[1];
   vector<Atom> atomList;
   Parameters param, unitCellParam, strainParam;
   
   int exBond = atoi(argv[2]);
   int cellFactor = atoi(argv[3]);
   Point cellMultiply(cellFactor);
  
   readAtoms(atomList, inputFile, unitCellParam);
   param = unitCellParam;
   multiplyCell(atomList, cellMultiply, unitCellParam, param);
   connectAtoms(atomList, exBond, param);
//   param.delRandBond(atomList);
   param.genBondList(atomList);

   outputAtoms(atomList, inputFile, param);

   return 0;
}

/*double energy (vector<Atom>& atomList, Parameters& param, double stepSize, double tolerance)
{
   double potential;
   int i;

   changeAtomParam(atomList, param);
   vector<Point> positions(param.pnt());
   for (i=0; i<param.pnt(); i++)
      positions[i] = atomList[i].getPos();

   param.genBondList(atomList);

   potential = optimizer(positions, param, stepSize, tolerance);
   for (i=0; i<param.pnt(); i++)
   {
      atomList[i].setPos(positions[i]);
   }

   return potential;
}*/

int mod(int x, int m) 
{
   return (x%m + m)%m;
}

double degRad(double deg)
{
   return (deg * pi() / 180);
}

double pi()
{
   return 3.14159265359;
}

void readAtoms(vector<Atom> & atoms, string fileName, Parameters & p)
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
      p.setDim(tmp, i);
   }
   Atom newAtom(&p);
   
   file >> buffer >> n >> buffer;
   p.pnt(n);
   for (i=0;i<p.pnt();i++)
   {
      file >> n >> x >> y >> z;
      tmp.x(x);
      tmp.y(y);
      tmp.z(z);
      newAtom.setPos(tmp);
      newAtom.setIndex(i);
      atoms.push_back(newAtom);
   }
   file.close();
}

//Assumes atom connections will be generated afterwards
void multiplyCell(vector<Atom> & atoms, Point & n, Parameters & op, Parameters & np)
{
   int i,j,k,m;
   unsigned int index;
   Atom newAtom(&np);
   Matrix3 factor;
   factor.setDiag(n);

   np = op;
   np.pnt( op.pnt() * (n.x() * n.y() * n.z()) );
   np.setDim( op.dim()*factor );

   for (i=0;i<n.x();i++)
   {
      for (j=0;j<n.y();j++)
      {
         for (k=0;k<n.z();k++)
         {
            for (m=0;m<op.pnt();m++)
            {
               Point pos;
               pos = atoms[m].getPos();
               index = m + k*op.pnt() + j*n.z()*op.pnt() + i*n.y()*n.z()*op.pnt();
               pos += (op.dim(0)*i);
               pos += (op.dim(1)*j);
               pos += (op.dim(2)*k);
               if (index < atoms.size())
               {
                  atoms[index].setParam(&np);
                  atoms[index].setPos(pos);
                  atoms[index].setIndex(index);
                  atoms[index].clearNeighbours();
               }
               else
               {
                  newAtom.setPos(pos);
                  newAtom.setIndex(index);
                  atoms.push_back(newAtom);
               }
            }
         }
      }
   }
}

void connectAtoms(vector<Atom> & atoms, int target, Parameters & p)
{
   int size = atoms.size();
   Point a,b;
   int nBond = 0;
   p.cxn(0);
   double cellDist;
   double epsilon = 1e-8;

   p.dist( fabs( (atoms[0].getPos() - atoms[target].getPos()).distance() ) );
   
   for (int i=0; i<size; i++)
   {
      for (int j=0; j<size; j++)
      {
         if (i == j)
            continue;
         a = atoms[i].getPos();
         b = atoms[j].getPos();
         cellDist = fabs( (p.getRealDiff(a,b)).distance());

         if ( fabs(cellDist - p.dist()) < epsilon)
         {
            atoms[i].setNeighbour(&atoms[j]);
            nBond++;
         }
      }
      if (nBond > p.cxn())
         p.cxn(nBond);
      nBond = 0;
   }
}

void outputAtoms(vector<Atom> & atoms, string fileName, Parameters & p)
{
   int i;
   fstream file ( (fileName.append(".out")).c_str(), fstream::in | fstream::out | fstream::trunc);
   
   file << fileName << endl << endl;

   file << p.pnt() << " atoms" << endl;
   file << p.pnt()*p.cxn()/2 << " bonds" << endl; //does not support deleting bonds
   file << "0 angles" << endl;
   file << "0 dihedrals" << endl;
   file << "0 impropers" << endl << endl;

   file << "1 atom types" << endl;
   file << "1 bond types" << endl;
   file << "0 angle types" << endl;
   file << "0 dihedral types" << endl;
   file << "0 improper types" << endl << endl;
   
   file << "0 " << p.dim(0).x() << " xlo xhi" << endl;
   file << "0 " << p.dim(1).y() << " ylo yhi" << endl;
   file << "0 " << p.dim(2).z() << " zlo zhi" << endl;
   file << p.dim(1).x() << ' ' << p.dim(2).x() << ' ' << p.dim(2).y() << " xy xz yz" << endl << endl;

   file << "Atoms" << endl << endl;

   for (i=0;i<p.pnt();i++)
   {
      file << i+1 << " 1 1 " << atoms[i].getPos('x') << " " << atoms[i].getPos('y') << " " 
           << atoms[i].getPos('z') << endl;
   }
   
   file << endl << "Bonds" << endl << endl;
   
   for (i=0;i<p.bonds().size();i++)
   {
      Bond outBond = p.bonds().get();
      file << i+1 << " 1 " << outBond.index(0) << ' ' << outBond.index(1) << endl;
   }

   file.close();
}

void changeAtomParam(vector<Atom> & atoms, Parameters & p)
{
   int size = atoms.size();
   for (int i=0; i<size; i++)
   {
      atoms[i].setParam(&p);
   }
}
