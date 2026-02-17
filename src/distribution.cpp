#include <iostream>
#include <fstream>
#include <math.h>
#include "macros.hpp"
#include "base.hpp"
#include <omp.h>


using namespace std;

string dots = "\033[1;34m:: \033[0m";
string reset = "\033[0m";
string bold  = "\033[1m";

ofstream fout;

/*--------------------------------------------------------------------*/

bool askToProceed() {
  string input;
  cout << dots << bold << "Proceed with scan? [Y/n] " << reset;
  getline(cin, input);
  if (input.empty()||(input=="Y")||(input == "y")) return true;
  if ((input == "N")||(input == "n")) return false;
  return false;
}

int main() {
  using namespace itp;

  init(10,10.,3.);

  if (!askToProceed()) return 0;

  return 0;
}
