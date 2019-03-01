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
namespace patch {
	template<typename T>
	std::string stringMaker(const T &n) {
		std::ostringstream stm;
		stm << n;
		return stm.str();
	}
}

using namespace std;

unsigned long long
gcd(unsigned long long x, unsigned long long y) {
	while (y != 0) {
		unsigned long long t = x % y;
		x = y;
		y = t;
	}
	return x;
}

unsigned long long
nChoosek(unsigned long long n, unsigned long long k) {
	if (k > n)
		throw std::invalid_argument("invalid argument in choose");
	unsigned long long r = 1;
	for (unsigned long long d = 1; d <= k; ++d, --n) {
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
		cout << "returning " << ans.size() << endl;
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

			cout << "returning " << ans.size() << endl;

			return ans;
		}
		else {
			vector<int> ans;

			vector<int> tertiary(l.begin() + 1, l.end());

			vector<int> secondary = kthCombination(k - i, tertiary, r);

			cout << "returning " << ans.size() << endl;
			return secondary;
		}
	}
}

void cover(string out);

int choose(int *pattern);

int choose(int *pattern, int localp, int localdSize, int safe);

int isCover(const int *differenceCover, int *testCover);

bool check();

void print(const int *differenceCover, string file);

void quicksave(unsigned long long pos, int startingThird, int *differenceCover);

inline int size(const int *cover, const int p);

int recursiveLock(int *differenceCover, int *testCover, int localp, int localdSize, int localThird, vector<int> starting);

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
		cout
			<< "[batchSize] - enter an integer number. This field sets the batch size. If not set the program will run continiously."
			<< endl;
		return 0;
	}

	int startValue = atoi(argv[2]);
	p = startValue;
	int numberToCompute = 10;

	if (argc >= 4)
		numberToCompute = atoi(argv[3]);
	if (argc >= 5)
		batchSize = atoi(argv[4]);

	string outFile = argv[1];

	if (id == 0) {
		cout << "mpi has started with " << nn << " and a batch size of " << batchSize << endl;
	}

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
			myfilea.close();;
		}

		if (batchSize != 0) {
			break;
		}

		p++;
	}

	cout << "Thread: " << id << " is finished" << endl;

	MPI_Finalize();

	return 0;
} // end main

void cover(string out) {
	// The min is the smallest possible difference cover of a set of size p.
	int min = 1;
	for (int x = 2; x < p; x++)
		if (x * (x - 1) + 1 >= p) {
			min = x;
			break;
		}
	dSize = min;

	// The max is the largest possible difference cover of a set of size p.
	int max = p;
	for (int x = 2; x < p / 2; x++)
		if (x * x > p) {
			max = 3 * (p + x - 1) / (2 * (x - 1)) + 4;
			break;
		}

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
		}
	}

	vector<int> starting;
	starting.push_back(0);
	starting.push_back(1);
	while (dSize <= max) {
		for (int i = int((p + 1) / 2) + 1; i > 1; i--) {
			cout << "the new starting third is " << i << endl;
			if (recursiveLock(differenceCover, testCover, p, dSize-2, i, starting)) {
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

			MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave

			if (check() || batchSize != 0) {
				//return;
			}
		}
		if (!check()) {
			dSize++;
		}
		else {
			return;
		}
	}

	delete[] testCover;
	delete[] differenceCover;
} // end cover

int recursiveLock(int *differenceCover, int *testCover, int localp, int localdSize, int localThird, vector<int> starting) {
	cout << "new step down " << localp << " dsize: " << localdSize << " localThird: " << localThird << " dc: ";

	for (int i = 0; i < dSize; i++) {
		cout << " " << differenceCover[i];
	}
	
	cout << " | starting ";
	for (int i = 0; i < starting.size(); i++) {
		cout << " " << starting[i];
	}
	cout << endl;

	for (int i = 0; i < localp; i++) {
		testCover[i] = 0;
	}

	int ans = 0;

	if (localThird < int((p + 1) / 2)) { //this is for below the half
		cout << "this is the below the half stuff" << endl;

		for (unsigned long long lock = p + 1 - localThird; lock < localp; lock++) {
            cout << "this is the first for lock: " << lock << " for " << localThird << " for localp " << localp << endl;
            for (int i = 0; i < localdSize+starting.size(); i++) { //the plus is to add the starting numbers
				if (i < starting.size()) {
					differenceCover[i] = starting[i];
				}
				else if (i == starting.size()) {
					differenceCover[i] = localThird;
				}
				else {
					differenceCover[i] = differenceCover[i - 1] + 1;
				}
			}

			differenceCover[starting.size()+localdSize-1] = lock; //the higher ups do the same thing so no need to worry about what comes after

			if (differenceCover[localdSize - 2] >= differenceCover[localdSize - 1]) {
				continue;
			}

			if (localdSize-2 > 1) { //-1 for the current third -1 for the lock				
				cout << "we are sending this down a level " << localdSize - 2 << endl;
				starting.push_back(localThird);

				for(int i = lock-(localdSize-2); i > localThird; i--) {
				    cout << "going through next recursion " << i << endl;
                    recursiveLock(differenceCover, testCover, lock, localdSize - 2, i, starting);
                }
				cout << "out of here bois" << localp << " " << localdSize << " " << localThird << endl;
				starting.erase(starting.begin() + starting.size() - 1);
				continue;
			}

			cout << endl;
			for (int i = 0; i < dSize; i++) {
				cout << differenceCover[i] << " ";
			}
			cout << endl;

			if (isCover(differenceCover, testCover)) {
				cout << "cover found" << endl;
				quicksave(0, localThird, differenceCover);

				ofstream myfilea;
				myfilea.open((pFile + ".txt").c_str(), ios::trunc);
				for (int a = 0; a < localdSize - 1; a++) {
					myfilea << differenceCover[a] << " ";
				}
				myfilea << differenceCover[localdSize - 1] << endl;
				myfilea.close();

				print(differenceCover, "testing.txt");

				//return 1;
			}

			cout << "going to do the general checking now; starting: " << starting.size()-1 << " localdsize " << localdSize << endl;
			if(localdSize-2 == 0) {
			    continue;
			}

			for (unsigned long long count = 0; choose(differenceCover, lock, localdSize-2, starting.size()); count++) { //this compensates for adding the num we check and subtracting 2 dsize
				cout << endl;
				cout << "lock  at " << lock << " for " << localThird << " for localp " << localp << " with localdsize " << localdSize << " with starting " << starting.size()-1 << endl;
				for (int i = 0; i < dSize; i++) {
					cout << differenceCover[i] << " ";
				}
				cout << endl;

				if (isCover(differenceCover, testCover)) {
					cout << "cover found" << endl;
					quicksave(count, localThird, differenceCover);

					ofstream myfilea;
					myfilea.open((pFile + ".txt").c_str(), ios::trunc);
					for (int a = 0; a < localdSize - 1; a++) {
						myfilea << differenceCover[a] << " ";
					}
					myfilea << differenceCover[localdSize - 1] << endl;
					myfilea.close();

					print(differenceCover, "testing.txt");

					//return 1;
				}
			}
		}
	}
	else { //this is for half or above
		unsigned long long startingIndex = 0;
		unsigned long long startValue = 0;
		
		cout << "regular above half " << starting.size() << " localdsize " << localdSize << " dsize is " << dSize << " the third is " << localThird << endl;
		for (int i = 0; i < localdSize+starting.size(); i++) {
			if (i < starting.size()) {
				differenceCover[i] = starting[i];
			}
			else if (i == starting.size()) {
				differenceCover[i] = localThird;
			}
			else {
				differenceCover[i] = differenceCover[i - 1] + 1;
			}
		}
		if (isCover(differenceCover, testCover)) {
			ofstream myfilea;
			myfilea.open((pFile + ".txt").c_str(), ios::trunc);
			for (int a = 0; a < localdSize - 1; a++) {
				myfilea << differenceCover[a] << " ";
			}
			myfilea << differenceCover[localdSize - 1] << endl;
			myfilea.close();

			print(differenceCover, "testing.txt");
			//return 1;
		}

		quicksave(startValue, localThird, differenceCover);
		unsigned long long writeTime = 47000000;

		cout << "this is the starting cover" << endl;
		for (int i = 0; i < dSize; i++) {
			cout << differenceCover[i] << " ";
		}
		cout << endl;

		cout << "this is regular " << localp << " " << localdSize-1 << " " << starting.size() << endl;
		for (unsigned long long z = startingIndex; choose(differenceCover, localp, localdSize-1, starting.size()); z++) {
			cout << endl;
			for (int i = 0; i < dSize; i++) {
				cout << differenceCover[i] << " ";
			}
			cout << endl;

			if (isCover(differenceCover, testCover)) {
				cout << "cover found" << endl;
				quicksave(z, localThird, differenceCover);

				ofstream myfilea;
				myfilea.open((pFile + ".txt").c_str(), ios::trunc);
				for (int a = 0; a < localdSize - 1; a++) {
					myfilea << differenceCover[a] << " ";
				}
				myfilea << differenceCover[localdSize - 1] << endl;
				myfilea.close();

				print(differenceCover, "testing.txt");

				//return 1;
			}
			else if (z % writeTime == 0) {
				if (check()) {
					cout << "someone else found the cover" << endl;

					ofstream myfilea;
					myfilea.open((pFile + ".txt").c_str(), ios::trunc);
					for (int a = 0; a < localdSize - 1; a++) {
						myfilea << differenceCover[a] << " ";
					}
					myfilea << differenceCover[localdSize - 1] << endl;
					myfilea.close();

					print(differenceCover, "testing.txt");

					//return 1;
				}

				quicksave(z, localThird, differenceCover);
			}
		}
	}
	return ans;
}


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
	}
} // end choose
int choose(int *pattern, int localp, int localdSize, int safe) {
	localdSize += safe+1;

	if (pattern[localdSize - 1] < localp - 1) {
		pattern[localdSize - 1]++;
		return 1;
	}
	else if (pattern[localdSize - 1] == localp - 1) {
		bool ret = 0;

		for (int a = 1; localdSize - a >= 0; a++) {
			if (localdSize - a <= safe) {
				cout << "returning 0 on " << localdSize - a << " value: " << pattern[localdSize-a] << endl;
				ret = 1;
			}
			pattern[localdSize - a]++;


			if (pattern[localdSize - a] <= localp - a) {
				for (int z = 1; z + localdSize - a < localdSize; z++) {
					pattern[z + localdSize - a] = pattern[localdSize - a] + z;
				}

				break;
			}
		}

		if (ret) {
			return 0;
		}
		return 1;
	}
} // end choose

int isCover(const int *differenceCover, int *testCover) {
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
} // end isCover //set the p per method and add a cut off

void quicksave(unsigned long long pos, int startingThird, int *differenceCover) {
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
	myfile << startingThird << endl;
	myfile.close();
}

void print(const int *differenceCover, string file) {
	ofstream out;
	out.open(file.c_str(), ios::app);

	out << setw(4) << p;
	int f = dSize;
	out << setw(9) << f;
	out << setw(8);
	for (int x = 0; x < dSize; x++) {
		out << differenceCover[x] << setw(3);
	}
	out << endl;
	out.close();
} // end print

bool check() {
	struct stat b;
	if (stat((pFile + ".txt").c_str(), &b) == 0) {
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
