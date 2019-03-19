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
        throw std::invalid_argument("invalid argument in nChoosek");
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d, --n) {
        unsigned long long g = gcd(r, d);
        r /= g;
        unsigned long long t = n / (d / g);
        if (r > std::numeric_limits<unsigned long long>::max() / t)
            throw std::overflow_error("overflow in nChoosek");
        r *= t;
    }
    return r;
}

vector<int> kthCombination(unsigned long long k, vector<int> l, int r) {
    if (r == 0) {
        vector<int> ans;
        return ans;
    } else if (l.size() == r) {
        return l;
    } else {
        unsigned long long i = nChoosek(l.size() - 1, r - 1);//calculare number of combinations
        if (k < i) {
            vector<int> ans;
            ans.push_back(l[0]);

            vector<int> tertiary(l.begin() + 1, l.end());

            vector<int> secondary = kthCombination(k, tertiary, r - 1);

            ans.insert(ans.end(), secondary.begin(), secondary.end());

            return ans;
        } else {
            vector<int> ans;

            vector<int> tertiary(l.begin() + 1, l.end());

            vector<int> secondary = kthCombination(k - i, tertiary, r);

            return secondary;
        }
    }
}

void cover(string out);

int choose(int *pattern, int localp, int localdSize, int safe);

int isCover(const int *differenceCover, int *testCover);

bool check();

void print(const int *differenceCover, string file);

inline int size(const int *cover, const int p);

int recursiveLock(int *differenceCover, int *testCover, int localp, int localdSize, int localThird, vector<int> starting);

void calcBounds(unsigned long long numNum, unsigned long long &iters, unsigned long long &starting);

int id; //the id of this process
int ierr;
int nn; //number of connected nodes
string pFile = "";
int p = -1;
int dSize = -1;
int groupNodes;
int groupid;
string combinedOut = "testing.txt";

void calcBounds(unsigned long long numNum, unsigned long long &iters, unsigned long long &starting) {
    int preGroup = groupNodes;
    int preId = groupid;

    groupNodes = int(groupNodes / numNum);

    if (groupNodes == 0) {
        groupNodes = 1;
    }

    if (groupid >= numNum) {
        groupid = groupid % numNum;
    }

    iters = numNum / preGroup;

    if (groupid < preGroup % numNum && preGroup > numNum) {
        groupNodes++;//+= nn / numNum;
    }

    if (groupid < numNum % preGroup) {
        iters += 1;

        if (groupid != 0) {
            starting += groupid * iters;
        }
    } else {
        starting += numNum % preGroup + groupid * iters;
    }

    if (groupNodes == 1) {
        groupid = 0;
    } else {
        groupid = preId / numNum;
    }
}

bool checkWrite(int *differenceCover, int *testCover) {
    if (isCover(differenceCover, testCover)) {
        ofstream myfilea;
        myfilea.open((pFile + ".txt").c_str(), ios::trunc);
        for (int a = 0; a < dSize - 1; a++) {
            myfilea << differenceCover[a] << " ";
        }
        myfilea << differenceCover[dSize - 1] << endl;
        myfilea.close();

        print(differenceCover, combinedOut);

        return true;
    }

    return false;
}

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

        return 0;
    }

    int startValue = atoi(argv[2]);
    p = startValue;
    int numberToCompute = 10;

    if (argc >= 4)
        numberToCompute = atoi(argv[3]);

    string outFile = argv[1];
    combinedOut = outFile;

    groupNodes = nn;
    groupid = id;

    while (p < startValue + numberToCompute) {
        pFile = patch::stringMaker(p);

        if (check()) {
            cout << "Thread: " << id << " found " << p << " before starting" << endl;
            p++;
            continue;
        }

        MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave

        cout << "Thread: " << id << " on: " << p << endl;

        cover(outFile);

        p++;
    }

    cout << "Thread: " << id << " is finished. Checked: " << totalCheck << endl;

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

    int *differenceCover = new int[max];
    int *testCover = new int[p];

    vector<int> starting;
    starting.push_back(0);
    starting.push_back(1);


    unsigned long long totalCombo = nChoosek(p - 2, dSize - 2) / 2;
    int startingThird = int((p + 1) / 2) + 1;
    struct stat buffer;
    if (stat((pFile + "_0.txt").c_str(), &buffer) == 0) {
        ifstream infile((pFile + "_0.txt").c_str());
        string line;
        getline(infile, line);

        string linea;
        getline(infile, linea);

        infile.close();

        {
            std::stringstream lineaStream(line);
            lineaStream >> dSize;
        }

        {
            std::stringstream lineaStream(linea);
            lineaStream >> startingThird;
        }

        cout << "Thread: " << id << " read " << dSize << " " << startingThird << endl;
    }

    while (dSize <= max) {
        cout << "dSize is now " << dSize << endl;
        for (int i = startingThird; i > 1; i--) {
            if (id == 0) {
                cout << "Thread: " << id << " is writing " << dSize << " " << i << endl;
                ofstream myfile;
                myfile.open((pFile + "_0.txt").c_str(), ios::trunc);
                myfile << dSize << endl;
                myfile << i << endl;
                myfile.close();
            }

            //this probably isn't neededd but is a safety for low cost
            groupid = id;
            groupNodes = nn;

            cout << "Thread: " << id << " the new starting third is " << i << endl;

            recursiveLock(differenceCover, testCover, p, dSize - 2, i, starting);

            cout << "Thread: " << id << " is done checking startingThird: " << i << endl;


            cout << "Thread: " << id << " is waiting for the rest" << endl;
            MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave
            if (check()) {
                groupid = id;
                groupNodes = nn;

                delete[] testCover;
                delete[] differenceCover;

                return;
            }
        }

        startingThird = int((p + 1) / 2) + 1;
        dSize++;
    }

    delete[] testCover;
    delete[] differenceCover;
} // end cover

int recursiveLock(int *differenceCover, int *testCover, int localp, int localdSize, int localThird, vector<int> starting) { //should starting be a reference and have a reserved size???
    //is this needed????
    for (int i = 0; i < localp; i++) {
        testCover[i] = 0;
    }

    if (localThird < int((p + 1) / 2)) { //this is for below the half
        int preGlobalGroup = groupNodes;
        int preGlobalId = groupid;
        unsigned long long startingGlobalLock = p + 1 - localThird;
        unsigned long long globalIters = 0;
        int numGlobalNum = localp - startingGlobalLock;

        calcBounds(numGlobalNum, globalIters, startingGlobalLock);
        unsigned long long globalUpperBound = startingGlobalLock + globalIters;

        for (unsigned long long lock = startingGlobalLock; lock < globalUpperBound; lock++) {
            for (int i = 0; i < localdSize + starting.size(); i++) { //the plus is to add the starting numbers
                if (i < starting.size()) {
                    differenceCover[i] = starting[i];
                } else if (i == starting.size()) {
                    differenceCover[i] = localThird;
                } else {
                    differenceCover[i] = differenceCover[i - 1] + 1;
                }
            }

            differenceCover[starting.size() + localdSize - 1] = lock; //the higher ups do the same thing so no need to worry about what comes after

            bool perfectRef = true;

            if (2 > starting.size()) {
                perfectRef = false;
            } else {
                for (int i = 2; i <= starting.size(); i++) {
                    if (p + 1 - differenceCover[i] != differenceCover[dSize + 1 - i]) {
                        perfectRef = false;
                        break;
                    }
                }
            }

            if (localdSize - 2 == 0) {
                if (checkWrite(differenceCover, testCover)) {
                    groupid = preGlobalId;
                    groupNodes = preGlobalGroup;

                    return 1;
                }

                continue;
            }

            if (localdSize - 2 > 1 && perfectRef) {
                starting.push_back(localThird);

                int preGroup = groupNodes;
                int preId = groupid;
                int numNum = lock - (localdSize - 2) - localThird - 1;
                unsigned long long startingLock = localThird + 1;
                unsigned long long iters = 0;

                if (numNum == 0) {
                    starting.erase(starting.begin() + starting.size() - 1);
                    continue;
                }

                calcBounds(numNum, iters, startingLock);

                for (int i = startingLock + iters; i >= startingLock; i--) {
                    if (recursiveLock(differenceCover, testCover, lock, localdSize - 2, i, starting) || check()) { //added check here in case none of the recursives have the time to check
                        groupid = preGlobalId;
                        groupNodes = preGlobalGroup;

                        return 1;
                    }
                }

                groupid = preId;
                groupNodes = preGroup;

                starting.erase(starting.begin() + starting.size() - 1);
                continue;
            }

            int localStartLock = localThird + 1;

            if (perfectRef && localThird < int(p / 2) + 1 && localdSize - 2 == 1) {
                localStartLock = int(p / 2) + 1;
            }

            vector<int> sc;
            sc.reserve(lock - localStartLock);

            for (int i = localStartLock; i < lock; i++) {
                sc.push_back(i);
            }

            if (sc.size() < localdSize - 2) {
                continue; //this needs to be continue because lock icreases each loop so this condition might not be true next time around
            }

            int preGroup = groupNodes;
            int preId = groupid;
            unsigned long long numNum = nChoosek(sc.size(), localdSize - 2);
            unsigned long long startValue = 0;
            unsigned long long instanceStart = 0;

            calcBounds(numNum, instanceStart, startValue);
            unsigned long long upperBound = instanceStart + startValue;

            vector<int> localStart = kthCombination(startValue, sc, localdSize - 2);
            for (int i = 0; i < localdSize + starting.size(); i++) {
                if (i < starting.size()) {
                    differenceCover[i] = starting[i];
                } else if (i == starting.size()) {
                    differenceCover[i] = localThird;
                } else {
                    differenceCover[i] = localStart[i - (1 + starting.size())];
                }
            }

            differenceCover[localdSize + starting.size() - 1] = lock;

            unsigned long long startingIndex = startValue + 1;

            if (checkWrite(differenceCover, testCover)) {
                groupid = preGlobalId;
                groupNodes = preGlobalGroup;

                return 1;
            }

            for (unsigned long long count = startingIndex; count < upperBound && choose(differenceCover, lock, localdSize - 2, starting.size()); count++) { //this compensates for adding the num we check and subtracting 2 dsize
                if (checkWrite(differenceCover, testCover)) {
                    groupid = preGlobalId;
                    groupNodes = preGlobalGroup;

                    return 1;
                } else if ((count-startingIndex) % 47000000 == 0 && check()) { //count-startingIndex makes it so that count is treated like 0
                    groupid = preId;
                    groupNodes = preGroup;

                    return 1;
                }
            }

            groupid = preId; //this can be the pre because it will loop through and we want it to be normal so it can be looped and assigned as before
            groupNodes = preGroup; //this is the equivalent of setting it to the resulting assigned global (normal again)
        }

        groupid = preGlobalId;
        groupNodes = preGlobalGroup;
    } else { //this is for half or above
        int preGroup = groupNodes;
        int preId = groupid;

        if (localp - localThird - 1 < localdSize - 1) {
            groupid = preGroup;
            groupNodes = preId;

            return 0; //we return here because there is no hope for change
        }

        unsigned long long numNum = nChoosek(localp - localThird - 1, localdSize - 1);
        unsigned long long startValue = 0;
        unsigned long long instanceStart = 0;

        calcBounds(numNum, instanceStart, startValue);
        unsigned long long upperBound = startValue + instanceStart;

        vector<int> localStart;
        {
            vector<int> sc;
            sc.reserve(localp - (localThird + 1));

            for (int i = localThird + 1; i < localp; i++) {
                sc.push_back(i);
            }

            localStart = kthCombination(startValue, sc, localdSize - 1);
        }

        for (int i = 0; i < localdSize + starting.size(); i++) {
            if (i < starting.size()) {
                differenceCover[i] = starting[i];
            } else if (i == starting.size()) {
                differenceCover[i] = localThird;
            } else {
                differenceCover[i] = localStart[i - (1 + starting.size())];
            }
        }

        unsigned long long startingIndex = startValue + 1;

        if (checkWrite(differenceCover, testCover)) {
            groupNodes = preGroup;
            groupid = preId;

            return 1;
        }

        for (unsigned long long z = startingIndex; z < upperBound && choose(differenceCover, localp, localdSize - 1, starting.size()); z++) {
            if (checkWrite(differenceCover, testCover)) {
                groupid = preId;
                groupNodes = preGroup;

                return 1;
            } else if ((z-startingIndex) % 47000000 == 0 && check()) {
                groupid = preId;
                groupNodes = preGroup;

                return 1;
            }
        }

        groupid = preId;
        groupNodes = preGroup;
    }

    //the ids should already be reset and we just need to say we didn't find anything
    return 0;
}


int choose(int *pattern, int localp, int localdSize, int safe) {
    localdSize += safe + 1;

    if (pattern[localdSize - 1] < localp - 1) {
        pattern[localdSize - 1]++;
        return 1;
    } else if (pattern[localdSize - 1] == localp - 1) {
        bool ret = 0;

        for (int a = 1; localdSize - a >= 0; a++) {
            if (localdSize - a <= safe) {
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
            return 0; //TODO: RETURN SOONER (CHECK)
        }
        return 1;
    }
} // end choose

int isCover(const int *differenceCover, int *testCover) {
    checkCount++;
    double *temp = (double *) testCover;
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
