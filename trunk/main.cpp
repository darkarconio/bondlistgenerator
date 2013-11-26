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

#define debug(x) cout << __LINE__ << ' ' << x << endl; cout.flush();

using namespace std;

vector<Atom> Atom::atomList;
Matrix Atom::adjMatrix = Matrix();
Parameters Atom::cellInfo = Parameters();

int main (int argc, char* argv[])
{
   string inputFile = argv[1];
   int exBond = atoi(argv[2]);
   int cellFactor = atoi(argv[3]);
   double percentDelBonds = atof(argv[4]);
   bool atomNotBond = !(bool)(strcmp(argv[5], "atom"));
   unsigned int modDelete = abs(atoi(argv[6]));
   srand ( time(NULL) );
   time_t start = time(NULL);

   Atom::readAtoms(inputFile);
   Atom::multiplyCell(Point(cellFactor));
   
   cout << "Connecting atoms..." << endl;
   Atom::connectAtoms(exBond);
   time_t connected = time(NULL);
   cout << "Connecting time: " << difftime(time(NULL), start) << endl;
  
   cout << "Deleting Atoms..." << endl;
   Atom::genBondList();
   Atom::genAtomDelList();
   Atom::genBondDelList();
   if (atomNotBond)
   {
      Atom::delPercentAtom(percentDelBonds);
   }
   else
   {
      Atom::delPercentBond(percentDelBonds, modDelete);
   }
   Atom::genBondList();
   Atom::genBondDelList();
   time_t deletions = time(NULL);
   cout << "Deletion time: " << difftime(time(NULL), connected) << endl;
   
   cout << "Generating angle list..." << endl;
   Atom::genAngleList();
   cout << "Angle list time: " << difftime(time(NULL), deletions) << endl;
   
   Atom::outputAtoms(inputFile);
   Atom::cellInfo.writeBondLoc();
   
   cout << "Average Coordination: " << 2*(double)Atom::cellInfo.nBonds()/(double)Atom::cellInfo.nAtomCandidates() << endl;
   cout << "Density: " << Atom::cellInfo.nAtomCandidates()/Atom::cellInfo.volume() << endl;
   cout << "Num Coord 2: " << Atom::getNumCoordX(2) << endl;
   cout << "Num Coord 3: " << Atom::getNumCoordX(3) << endl;
   cout << "Num Coord 4: " << Atom::getNumCoordX(4) << endl;
   cout << "Deletable Bonds: " << Atom::cellInfo.nBondCandidates() << endl;
   cout << "Deletable Atoms: " << Atom::cellInfo.nAtomCandidates() << endl;
   cout << "Total time: " << difftime(time(NULL), start) << endl;

   return 0;
}
