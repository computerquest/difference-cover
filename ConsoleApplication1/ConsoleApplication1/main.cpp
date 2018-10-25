// Ankur Gupta
// Mathematics Honors Thesis
// Difference Covers
// 8 May 2000

#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <string>
#include <algorithm>
#include <mpi.h>

using namespace std;

void cover(const int p, ofstream &out);
void print(const int p, string * differenceCover, ofstream& out);
inline int size(const int * cover, const int p);
int isCover(string cover, const int p);

string comb(const int P, const int K)
{
	std::string bitmask(K, 1); // K leading 1's
	bitmask.resize(P, 0); // N-K trailing 0's
						  // print integers and permute bitmask
	do {
		if (!bitmask[1]) {
			return "";
		}

		if (isCover(bitmask, P)) {
			return bitmask;
		}
	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));

	return "";
}

int id; //the id of this process
int ierr;
int p; //number of connected nodes

/*COMMANDS TO RUN/COMPILE
 * COMPILE: mpicxx \-o main main.cpp
 * RUN: mpiexec \-n [NUMBER OF INSTANCES] main [PARAMETERS]
 * EXAMPLE: mpiexec \-n 4 main stuff.txt 100
 */
int main(int argc, char *argv[]) {
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &p);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (argc == 1) {
        cout << "This program computes difference covers of a range specified by" << endl;
        cout << "the user. This program then saves the covers to a file that you" << endl;
        cout << "specify. The following format is used to execute this command." << endl;
    }
    if (argc < 3) {
        cout << "Usage: mhonors <filename> <startValue> [numberToCompute]" << endl << endl;
        cout << "<filename> - enter a filename." << endl << endl;
        cout << "<startValue> - enter an integer number. This is the first cover" << endl;
        cout << "               the program will attempt to find." << endl << endl;
        cout << "[numberToCompute] - enter an integer number. This field is optional, and" << endl;
        cout << "                    will tell the program to compute this many covers," << endl;
        cout << "                    starting at <startValue>. If you choose not to provide" << endl;
        cout << "                    this value, the default value is 10." << endl;
        return 0;
    }

    int startValue = atoi(argv[2]);
    int numberToCompute = 10;

    if (argc == 4)
        numberToCompute = atoi(argv[3]);
    ofstream out;
    out.open(argv[1], ios::out);
    if (!out) {
        cerr << "File could not be opened." << endl;
        exit(1);
    } // end if
    cout << setw(4) << "p" << setw(9) << "f(p)" << setw(37) << "difference cover" << endl;
    out << setw(4) << "p" << setw(9) << "f(p)" << setw(37) << "difference cover" << endl;

    int instanceStart = numberToCompute/p;
    if(id < numberToCompute%p) {
        instanceStart += 1;

        if(id != 0) {
            startValue += id*instanceStart;
        }
    } else {
        startValue += numberToCompute%p+id*instanceStart;
    }

    cout << "Thread: " << id << " starting on: " << startValue << " ending on: " << startValue+instanceStart-1 << endl;
    for (int x = startValue; x < startValue + instanceStart; x++) {
        cover(x, out);
    } // end for
    out.close();
    return 0;
} // end main

void cover(const int p, ofstream &out) {
	// The min is the smallest possible difference cover of a set of size p.
	int min = 1;
	for (int x = 2; x < p; x++) {
		if (x * (x - 1) + 1 >= p) {
			min = x;
			break;
		} // end if
	}
	// The max is the largest possible difference cover of a set of size p.
	int max = p;
	for (int x = 2; x < p / 2; x++) {
		if (x * x > p) {
			max = 3 * (p + x - 1) / (2 * (x - 1)) + 4;
			break;
		} // end if
	}

	string ans = "";
	for (int x = min; x <= max; x++) {
		ans = comb(p, x);
		if (ans != "") {
			print(p, ans, out);
			return;
		} // end if
	}
} // end cove

int isCover(string cover, const int p) {
	int * testCover = new int[p];

	//this does all the mod operations to find which numbers are generated
	testCover[0] = 1;
	for (int x = 0; x < p - 1; x++)
		if (cover[x])
			for (int y = x + 1; y < p; y++)
				if (cover[y]) {
					// testCover[(x - y + p) % p] = 1;
					// testCover[(y - x + p) % p] = 1;
					// Use below instead of the above, since y > x, (y-x) < p, and (y-x)%p + (x-y)%p = p x=0 y=2 n=100 2+98 so works
					int q = y - x;
					testCover[q] = testCover[p - q] = 1;
					// testCover[p - q] = 1;
				} // end if

	if (size(testCover, p) == p) { //could optimize so it can exit early if there is a zero
		delete[] testCover;
		return 1;
	}

	delete[] testCover;
	return 0;
}

void print(const int p, string differenceCover, ofstream &out) {
    cout << setw(4) << p;
    out << setw(4) << p;
    int f = differenceCover.length();
    cout << setw(9) << f;
    out << setw(9) << f;
    cout << setw(8);
    out << setw(8);
    for (int x = 0; x < p; x++)
        if (differenceCover[x]) {
            cout << x << setw(3);
            out << x << setw(3);
        } // end if
    cout << endl;
    out << endl;
} // end print

int size(const int *cover, const int p) {
    int f = 0;
    for (int x = 0; x < p; x++)
        f += cover[x];
    return f;
} // end size
