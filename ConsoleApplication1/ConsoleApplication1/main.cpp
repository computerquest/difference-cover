// Ankur Gupta
// Mathematics Honors Thesis
// Difference Covers
// 8 May 2000

#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <string>
#include <sys/stat.h>
#include <sstream>

using namespace std;

void cover(const int p, ofstream &out);

int coverOfSize(const int p, int& dSize, int *differenceCover, int *testCover);

int choose(const int p, int dSize, int *pattern, int &beginning, int init);

int isCover(const int p, const int *differenceCover, int *testCover);

void print(const int p, const int *differenceCover, ofstream &out);

inline int size(const int *cover, const int p);


int id; //the id of this process
int ierr;
int p; //number of connected nodes
string pFile = "";

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
    if (!out) {
        cerr << "File could not be opened." << endl;
        exit(1);
    } // end if
    cout << setw(4) << "p" << setw(9) << "f(p)" << setw(37) << "difference cover" << endl;
    //out << setw(4) << "p" << setw(9) << "f(p)" << setw(37) << "difference cover" << endl;

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

		std::stringstream ss;
		ss << x;
		pFile = ss.str()+".txt";

		out.open(argv[1], ios::out);
        cover(x, out);
		out.close();
    } // end for

	MPI_Finalize();

    return 0;
} // end main

void cover(const int p, ofstream &out) {
    // The min is the smallest possible difference cover of a set of size p.
    int min = 1;
    for (int x = 2; x < p; x++)
        if (x * (x - 1) + 1 >= p) {
            min = x;
            break;
        } // end if

    // The max is the largest possible difference cover of a set of size p.
    int max = p;
    for (int x = 2; x < p / 2; x++)
        if (x * x > p) {
            max = 3 * (p + x - 1) / (2 * (x - 1)) + 4;
            break;
        } // end if

    int *differenceCover = new int[p];
    int *testCover = new int[p];
    for (int x = min; x <= max; x++)
        if (coverOfSize(p, x, differenceCover, testCover)) {
            print(p, differenceCover, out);
            break;
        } // end if
    delete[] differenceCover;
    delete[] testCover;
} // end cover

int coverOfSize(const int p, int& dSize, int *differenceCover, int *testCover) {
    // differenceCover is an array that stores the current pattern generated

	struct stat buffer;
	if (stat(pFile.c_str(), &buffer) == 0) {
		ifstream infile(pFile.c_str());
		string line;
		getline(infile, line);
		infile.close();

		for (int i = 0; i < line.length(); i++) {
			differenceCover[i] = line[i]-48; //the character value of the 0 or 1 -48 gives the actual value
		}

		dSize = size(differenceCover, p);

		int beginning = 0;
		while (choose(p, dSize, differenceCover, beginning, p-1))
			if (isCover(p, differenceCover, testCover))
				return 1;
	}
	else {
		int beginning = 1;
		while (choose(p, dSize, differenceCover, beginning, 0))
			if (isCover(p, differenceCover, testCover))
				return 1;
	}

    return 0;
} // end coverOfSize

int choose(const int p, int dSize, int *pattern, int &beginning, int init) {
	static int index = init;

	int z = 0;
    if (beginning) {
        for (int x = 0; x < p; x++)
            pattern[x] = 0;
        for (int x = 0; x < dSize - 1; x++)
            pattern[x] = 1;
        pattern[p - 1] = 1;
        index = p - 1;
        beginning = 0;
        return 1;
    }

    if (pattern[index - 1] == 0) {
        pattern[index] = 0;
        pattern[--index] = 1;
        return 1;
    } // end if
    else if (index != p - 1) { // && pattern[index - 1] == 1)
        pattern[--index] = 0;
        pattern[p - 1] = 1;
        index = p - 1;

		if (z % 10000 == 0) {
			//THIS IS WHERE TO WRITE TO THE FILE
			ofstream myfile;
			myfile.open(pFile.c_str(), ios::trunc);
			for (int i = 0; i < p; i++) {
				myfile << pattern[i];
			}
			myfile.close();
		} 
		z++;

        return 1;
    } // end else if
    else { // pattern[index - 1] == 1 && index == p - 1        
		int flipCount = 0;
        while (pattern[--index]) {
            pattern[index] = 0;
            flipCount++;
        } // end while
        while (index > 0 && !pattern[--index]);
		if (index == 1) {
			return 0;
		}
        else {
            pattern[index] = 0;
            flipCount++;
            for (int y = 1; y <= flipCount; y++)
                pattern[index + y] = 1;
            index = p - 1;
            return 1;
        } // end else
    } // end else
} // end choose

int isCover(const int p, const int *differenceCover, int *testCover) {
    // Takes half the time, unless the double pointer casting takes too long. (Check this.)
    double *temp = (double *) testCover;
    for (int x = 0; x < p / 2; x++)
        temp[x] = 0;
    // Faster to not check if it is even or not.
    testCover[p - 1] = 0;

    testCover[0] = 1;
    for (int x = 0; x < p - 1; x++)
        if (differenceCover[x])
            for (int y = x + 1; y < p; y++)
                if (differenceCover[y]) {
                    // testCover[(x - y + p) % p] = 1;
                    // testCover[(y - x + p) % p] = 1;
                    // Use below instead of the above, since y > x, (y-x) < p, and (y-x)%p + (x-y)%p = p
                    int q = y - x;
                    testCover[q] = testCover[p - q] = 1;
                    // testCover[p - q] = 1;
                } // end if

    if (size(testCover, p) == p)
        return 1;
    return 0;
} // end isCover

void print(const int p, const int *differenceCover, ofstream &out) {
    cout << setw(4) << p;
    out << setw(4) << p;
    int f = size(differenceCover, p);
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
