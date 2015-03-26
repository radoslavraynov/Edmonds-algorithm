// Radoslav Zlatev, rzlatev@uni-bonn.de
// Programming Assignment, 26.11.2008

#include <vector>
#include <iostream>
#include <fstream>

#include "classes2.h"

using namespace std;

int main (int argc, char** argv) {
	if (argc == 1) {
		cout << "WARNING: NO SPECIFICATION FILE" << endl;
		return 0;
	}
		
	Graph GG;
		GG.step1(argv[1]);
		GG.step2();

	return 0;
}
