#include <iostream>
#include <fstream>

#include <vector>
#include <string>

#include <cmath>
#include <cstdlib>
#include <ctime>

#include <atom.hpp>
#include <point.hpp>
#include <fxn.hpp>
#include <bond.hpp>
//#include <interaction.hpp>

#define debug(x) cout << __LINE__ << ' ' << x << endl; cout.flush();

using namespace std;

int mod(int, int); 
double degRad(double);

const long double PI = 3.141592653589793238;

vector<Atom> Atom::atomList;
Matrix Atom::adjMatrix = Matrix();
Parameters Atom::cellInfo = Parameters();

int main (int argc, char* argv[])
{
   string inputFile = argv[1];
   int exBond = atoi(argv[2]);
   int cellFactor = atoi(argv[3]);
   int numDelBonds = atoi(argv[4]);
   srand ( time(NULL) );

   Atom::readAtoms(inputFile);
   Atom::multiplyCell(Point(cellFactor));
   Atom::connectAtoms(exBond);
  
   Atom::genBondList();
   Atom::genDelList();
   
   for (int i=0; i<numDelBonds;)
   {
      if (Atom::delRandBond())
         i++;
      if (Atom::cellInfo.nCandidates() == 0)
         break;
   }
   
   Atom::genBondList();
   Atom::genAngleList();
   
   Atom::outputAtoms(inputFile);
   
   cout << (double)Atom::cellInfo.nBonds()/(double)Atom::atomList.size() << endl;//Bond density
   cout << Atom::cellInfo.nCandidates() << endl;

   return 0;
}

int mod(int x, int m) 
{
   return (x%m + m)%m;
}

double degRad(double deg)
{
   return (deg * PI / 180);
}
