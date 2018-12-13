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
#include <vector>
#include <limits> 
/*
import operator as op
def nCr(n, r):
	r = min(r, n-r)
	numer = reduce(op.mul, xrange(n, n-r, -1), 1)
	denom = reduce(op.mul, xrange(1, r+1), 1)
	return numer//denom

k = combination number, l = numbers to choose from, r = slots available
def kthCombination(k, l, r):

	if r == 0:
		return []
	elif len(l) == r:
		return l
	else:
		i=nCr(len(l)-1, r-1)
		if k < i:
			return l[0:1] + kthCombination(k, l[1:], r-1)
		else:
			return kthCombination(k-i, l[1:], r)

#kthCombination(3000, [0,1,2,3, 4,5,6,7,8,9,10,11,12,13,14], 5)
for k in range(1,3000):
global counter
counter = 0
print kthCombination(k, [0,1,2,3, 4,5,6,7,8,9,10,11,12,13,14], 5), counter

*/
namespace patch
{
	template < typename T > std::string stringMaker(const T& n)
	{
		std::ostringstream stm;
		stm << n;
		return stm.str();
	}
}

using namespace std;

unsigned long long
gcd(unsigned long long x, unsigned long long y)
{
	while (y != 0)
	{
		unsigned long long t = x % y;
		x = y;
		y = t;
	}
	return x;
}

unsigned long long
nChoosek(unsigned long long n, unsigned long long k)
{
	if (k > n)
		throw std::invalid_argument("invalid argument in choose");
	unsigned long long r = 1;
	for (unsigned long long d = 1; d <= k; ++d, --n)
	{
		unsigned long long g = gcd(r, d);
		r /= g;
		unsigned long long t = n / (d / g);
		if (r > std::numeric_limits<unsigned long long>::max() / t)
			throw std::overflow_error("overflow in choose");
		r *= t;
	}
	return r;
}

vector<int> kthCombination(unsigned long long k, vector<int> l, int r) {
	if (r == 0) {
		vector<int> ans;
		return ans;
	}
	else if (l.size() == r) {
		return l;
	}
	else {
		unsigned long long i = nChoosek(l.size() - 1, r - 1);//calculare number of combinations
		if (k < i) {
			vector<int> ans;
			ans.push_back(l[0]);

			vector<int> tertiary(l.begin() + 1, l.end());

			vector<int> secondary = kthCombination(k, tertiary, r - 1);

			ans.insert(ans.end(), secondary.begin(), secondary.end());

			return ans;
		}
		else {
			vector<int> ans;

			vector<int> tertiary(l.begin() + 1, l.end());

			vector<int> secondary = kthCombination(k - i, tertiary, r);

			return secondary;
		}
	}
}

void cover(string out);

int coverOfSize(int *differenceCover, int *testCover);

int choose(int *pattern);

int isCover(const int *differenceCover, int *testCover);

bool check();

void print(const int *differenceCover, string file);

void quicksave(unsigned long long pos, const int *differenceCover);

inline int size(const int *cover, const int p);

int id; //the id of this process
int ierr;
int nn; //number of connected nodes
string pFile = "";
unsigned long long batchSize = 0;
unsigned long long lastComb = 0;
int p = -1;
int dSize = -1;

/*COMMANDS TO RUN/COMPILE
 * COMPILE: mpicxx \-o main main.cpp
 * RUN: mpiexec \-n [NUMBER OF INSTANCES] main [PARAMETERS]
 * EXAMPLE: mpiexec \-n 4 main stuff.txt 100
 */
int main(int argc, char *argv[]) {
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nn);
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
		cout << "[batchSize] - enter an integer number. This field sets the batch size. If not set the program will run continiously." << endl;
		return 0;
	}

	int startValue = atoi(argv[2]);
	p = startValue;
	int numberToCompute = 10;

	if (argc >= 4)
		numberToCompute = atoi(argv[3]);
	if (argc >= 5)
		batchSize = atoi(argv[4]);
	/*ofstream out;
	if (!out) {
		cerr << "File could not be opened." << endl;
		exit(1);
	} // end if*/

	string outFile = argv[1];

	//cout << setw(4) << "p" << setw(9) << "f(p)" << setw(37) << "difference cover" << endl;
	//out << setw(4) << "p" << setw(9) << "f(p)" << setw(37) << "difference cover" << endl;

	cout << "mpi has started with " << nn << " and a batch size of " << batchSize << endl;

	while (p < startValue + numberToCompute) {
		pFile = patch::stringMaker(p);

		if (check()) {
			p++;
			continue;
		}
		cout << "Thread: " << id << " on: " << p << endl;

		cover(outFile);

		if (id == 0 && batchSize != 0) {
			unsigned long long numCombo = nChoosek(p - 2, dSize - 2);

			cout << "Thread: " << id << " batch completed." << endl;
			ofstream myfilea;
			myfilea.open((pFile + "_batch.txt").c_str(), ios::trunc);

			if (nn * batchSize + lastComb <= numCombo) {
				myfilea << nn * batchSize + lastComb << endl;
				myfilea << dSize << endl;
			}
			else {
				myfilea << 0 << endl;
				myfilea << dSize + 1 << endl;
			}
			myfilea.close();

			break;
		}

		p++;
	} // end for 

	cout << "Thread: " << id << " is finished" << endl;

	MPI_Finalize();

	cout << "done" << endl;
	return 0;
} // end main

void cover(string out) {
	// The min is the smallest possible difference cover of a set of size p.
	int min = 1;
	for (int x = 2; x < p; x++)
		if (x * (x - 1) + 1 >= p) {
			min = x;
			break;
		} // end if
	dSize = min;

	// The max is the largest possible difference cover of a set of size p.
	int max = p;
	for (int x = 2; x < p / 2; x++)
		if (x * x > p) {
			max = 3 * (p + x - 1) / (2 * (x - 1)) + 4;
			break;
		} // end if

	//these can and should be minimized in the future
	int *differenceCover = new int[max];
	int *testCover = new int[p];

	if (batchSize != 0) {
		struct stat buffer;
		if (stat((pFile + "_batch.txt").c_str(), &buffer) == 0) {
			ifstream infile((pFile + "_batch.txt").c_str());
			string line;
			getline(infile, line);

			string linea;
			getline(infile, linea);

			infile.close();

			lastComb = strtoull(line.c_str(), 0, 16);
			dSize = atoi(linea.c_str());

			cout << "reading yeilded: " << dSize << " " << lastComb << endl;
		}
	}

	while (dSize <= max) {
		if (coverOfSize(differenceCover, testCover)) {
			if (isCover(differenceCover, testCover)) {
				ofstream myfilea;
				myfilea.open((pFile + ".txt").c_str(), ios::trunc);
				for (int a = 0; a < dSize - 1; a++) {
					myfilea << differenceCover[a] << " ";
				}
				myfilea << differenceCover[dSize - 1] << endl;
				myfilea.close();

				print(differenceCover, out);
			}
		}

		cout << "Thread " << id << " is waiting for the rest" << endl;

		MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave

		if (check() || batchSize != 0) {
			return;
		}

		dSize++;
	}

	delete[] testCover;
	delete[] differenceCover;
} // end cover

int coverOfSize(int *differenceCover, int *testCover) {
	unsigned long long startValue = 0;
	unsigned long long instanceStart = 0;

	unsigned long long numCombo = nChoosek(p - 2, dSize - 2);

	if (batchSize == 0 || nn * batchSize + lastComb > numCombo) {
		cout << "irregular init " << nn * batchSize+lastComb << " against " << numCombo << endl;
		instanceStart = numCombo / nn;
		if (id < numCombo%nn) {
			instanceStart += 1;

			if (id != 0) {
				startValue += id * instanceStart;
			}
		}
		else {
			startValue += numCombo % nn + id * instanceStart;
		}
	}
	else if (batchSize != 0) {
		startValue = lastComb + batchSize * id;
		instanceStart = batchSize;
	}

	unsigned long long startingIndex = startValue + 1;

	struct stat buffer;
	bool reset = true; //this is to make sure when the computer goes to the next dSize it still works
	if (batchSize == 0 && stat((pFile + "_" + patch::stringMaker(id) + ".txt").c_str(), &buffer) == 0) {
		ifstream infile((pFile + "_" + patch::stringMaker(id) + ".txt").c_str());
		string line;
		getline(infile, line);

		string linea;
		getline(infile, linea);

		infile.close();

		int value = 0;
		std::stringstream  lineStream(line);
		int count = 0;
		for (int a = 0; a < dSize && lineStream >> value; a++)
		{
			// Add the integers from a line to a 1D array (vector)
			differenceCover[a] = value;
			count++;
		}

		if (count >= dSize) {
			dSize = count;

			std::stringstream  lineaStream(linea);
			lineaStream >> startingIndex;

			reset = false;
		}
	}


	if (reset) {
		vector<int> sc;
		for (int i = 0; i < p; i++) {
			sc.push_back(i);
		}
		vector<int> starter = kthCombination(startValue, sc, dSize);

		for (int i = 0; i < starter.size(); i++) {
			differenceCover[i] = starter[i];
		}
	}

	for (int i = 0; i < p; i++) {
		testCover[i] = 0;
	}

	cout << "Thread " << id << " starting cover for: " << p << " is: [" << startingIndex << "] ";
	for (int x = 0; x < dSize; x++) {
		cout << differenceCover[x] << " ";
	}

	{
		cout << " ending cover is: [" << startValue + instanceStart - 1<< "] ";
		{
			vector<int> sc;
			for (int i = 0; i < p; i++) {
				sc.push_back(i);
			}
			vector<int> starter = kthCombination(startValue + instanceStart - 1, sc, dSize);

			for (int i = 0; i < starter.size(); i++) {
				cout << starter[i] << " ";
			}
		}
	}

	cout << " out of " << numCombo << endl;

	if (isCover(differenceCover, testCover)) {
		//cout << "we found " << p << "//////////////////////////////////////////////" << endl;

		quicksave(0, differenceCover);

		return 1;
	}

	unsigned long long writeTime = 47000000;

	cout << "Thread " << id << " starting Index: " << startingIndex << " ending condition " << startValue + instanceStart << endl;
	for (unsigned long long z = startingIndex; z < startValue + instanceStart && choose(differenceCover); z++) {
		if (isCover(differenceCover, testCover)) {
			//cout << p << " is done /////////////////////////////////////////////////////////////" << endl;

			quicksave(z, differenceCover);

			return 1;
		}
		else if (z % writeTime == 0) {
			if (check()) {
				return 1;
			}

			cout << "Thread " << id << " writing for " << p << " complete " << (double)(z) / (startValue+instanceStart) << endl;


			quicksave(z, differenceCover);
		}

		/*cout << "Thread " << id << " z "<< z << " out of " << (startValue+instanceStart) << " " << z % (int)(.01*(startValue + instanceStart)) << " ";
		for (int i = 0; i < dSize; i++) {
			cout << differenceCover[i] << " ";
		}
		cout << endl;*/
	}

	return 0;
} // end coverOfSize

int choose(int *pattern) {
	if (pattern[dSize - 1] < p - 1) {
		pattern[dSize - 1]++;
		return 1;
	}
	else if (pattern[dSize - 1] == p - 1) {
		for (int a = 1; dSize - a >= 0; a++) {
			pattern[dSize - a]++;
			if (pattern[dSize - a] <= p - a) {
				for (int z = 1; z + dSize - a < dSize; z++) {
					pattern[z + dSize - a] = pattern[dSize - a] + z;
				}

				break;
			}
		}

		return 1;
	} // end else if
//	cout << "something is super wrong ////////////////////////////////////////////////////////// " << pattern[dSize-1]  << " " << p << endl;
} // end choose

int isCover( const int *differenceCover, int *testCover) {
	double *temp = (double *)testCover;
	for (int x = 0; x < p / 2; x++)
		temp[x] = 0;
	testCover[p - 1] = 0;
	testCover[0] = 1;

	for (int x = 0; x < dSize - 1; x++) {
		for (int y = x + 1; y < dSize; y++) {
			int xa = differenceCover[x];
			int ya = differenceCover[y];

			int q = ya - xa;
			testCover[q] = testCover[p - q] = 1;
		}
	}

	if (size(testCover, p) == p)
		return 1;
	return 0;
} // end isCover

void quicksave( unsigned long long pos, const int *differenceCover) {
	if (batchSize != 0) {
		return;
	}

	ofstream myfile;
	myfile.open((pFile + "_" + patch::stringMaker(id) + ".txt").c_str(), ios::trunc);
	for (int a = 0; a < dSize - 1; a++) {
		myfile << differenceCover[a] << " ";
	}
	myfile << differenceCover[dSize - 1] << endl;
	myfile << pos << endl;
	myfile.close();
}

void print(const int *differenceCover, string file) {
	ofstream out;
	out.open(file.c_str(), ios::app);

	cout << setw(4) << p;
	out << setw(4) << p;
	int f = dSize;
	cout << setw(9) << f;
	out << setw(9) << f;
	cout << setw(8);
	out << setw(8);
	for (int x = 0; x < dSize; x++) {
		cout << differenceCover[x] << setw(3);
		out << differenceCover[x] << setw(3);
	} // end if
	out << endl;
	cout << endl;
	out.close();
} // end print

bool check() {
	struct stat b;
	if (stat((pFile + ".txt").c_str(), &b) == 0) {
		//cout << "someone else found it //////////////////////////////////////" << endl;
		return 1;
	}

	return 0;
}

int size(const int *cover, const int p) {
	int f = 0;
	for (int x = 0; x < p; x++)
		f += cover[x];
	return f;
} // end size
