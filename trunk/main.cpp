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
   double percentDelBonds = atof(argv[4]);
   bool modDelete = !(bool)(strcmp(argv[5], "guide"));
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
   Atom::genDelList();
   Atom::delPercentBond(percentDelBonds, modDelete);
   Atom::genBondList();
   time_t deletions = time(NULL);
   cout << "Deletion time: " << difftime(time(NULL), connected) << endl;
   
   cout << "Generating angle list..." << endl;
   Atom::genAngleList();
   cout << "Angle list time: " << difftime(time(NULL), deletions) << endl;
   
   Atom::outputAtoms(inputFile);
   
   cout << "Average Coordination: " << 2*(double)Atom::cellInfo.nBonds()/(double)Atom::atomList.size() << endl;
   cout << "Num Coord 2: " << Atom::getNumCoordX(2) << endl;
   cout << "Num Coord 3: " << Atom::getNumCoordX(3) << endl;
   cout << "Num Coord 4: " << Atom::getNumCoordX(4) << endl;
   cout << "Deletable Atoms: " << Atom::cellInfo.nCandidates() << endl;
   cout << "Total time: " << difftime(time(NULL), start) << endl;

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
