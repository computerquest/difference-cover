//
// Created by jstigter on 6/10/19.
//
#include <string>
#include "main.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <bits/stat.h>

using namespace std;

///////////////////////////////////////////////////GLOBALS///////////////////////////
int p;
int dSize;
vector<int> differenceCover; //this is a vector but it acts as a stack
int *testCover;
string pFile;

///////////////////mpi
vector<int> groupNodes; //this is a stack
vector<int> groupid; //this is a stack
int id; //the id of this process
int ierr;
int nn; //number of connected nodes

///////////////////////////////////////////////////RESOURCE FUNCTIONS///////////////////////////
void calcBounds(unsigned long long numNum, unsigned long long &iters, unsigned long long &starting) {
    int groupNodes = groupNodes.back();
    int groupid = groupid.back();

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

    //TODO make sure this isn't buggy af
    this->groupNodes.push_back(groupNodes);
    this->groupid.push_back(groupid);
}

unsigned long long gcd(unsigned long long x, unsigned long long y) {
    while (y != 0) {
        unsigned long long t = x % y;
        x = y;
        y = t;
    }
    return x;
}

unsigned long long nChoosek(unsigned long long n, unsigned long long k) {
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

/////////////////////////////
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

    groupNodes.push_back(nn);
    groupid.push_back(id);

    while (p < startValue + numberToCompute) {
        pFile = patch::stringMaker(p);

        if (check()) {
            cout << "Thread: " << id << " found " << p << " before starting" << endl;
            p++;
            continue;
        }

        MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave

        cout << "Thread: " << id << " on: " << p << endl;

        startSearch();

        p++;
    }

    cout << "Thread: " << id << " is finished." << endl;

    MPI_Finalize();

    return 0;
} // end main

void startSearch() {
    //calcs min
    int min = 1;
    for (int x = 2; x < p; x++) {
        if (x * (x - 1) + 1 >= p) {
            min = x;
            break;
        }
    }

    //calcs max
    int max = p;
    for (int x = 2; x < p / 2; x++) {
        if (x * x > p) {
            max = 3 * (p + x - 1) / (2 * (x - 1)) + 4;
            break;
        }
    }

    differenceCover.reserve(max);
    testCover = new int[p];
    dSize = min;


    unsigned long long totalCombo = nChoosek(p - 2, dSize - 2) / 2;
    int startingThird = int((p + 1) / 2) + 1;

    //TODO update and clean up this copy paste
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

    differenceCover.push_back(0);
    differenceCover.push_back(1);

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

            cout << "Thread: " << id << " the new starting third is " << i << endl;

            //TODO this section will be subject to change because it makes the start uneven
            searchCovers(startingThird, dSize - 2, true);

            cout << "Thread: " << id << " is done checking startingThird: " << i << endl;


            cout << "Thread: " << id << " is waiting for the rest" << endl;
            MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave
            if (check()) {
                delete[] testCover;

                return;
            }
        }

        startingThird = int((p + 1) / 2) + 1;
        dSize++;
    }

    delete[] testCover;
}

int searchCovers(int localThird, int localdSize, bool perfect) {
    //TODO smooth this out and get rid of it (this is a good temp fix)
    int localp = differenceCover.back();

    if (localp == 1) {
        localp = p;
    }

    differenceCover.push_back(localThird);

    //TODO any returns need to have difference cover and the mpi id stuff pop before hand
    if (localThird < int((p + 1) / 2)) { //this is for below the half //TODO have ending function?
        //this is all init setup for splitting
        unsigned long long startingGlobalLock = p + 1 - localThird;
        unsigned long long globalIters = 0;
        int numGlobalNum = localp - startingGlobalLock;

        calcBounds(numGlobalNum, globalIters, startingGlobalLock);
        unsigned long long globalUpperBound = startingGlobalLock + globalIters;


        for (unsigned long long lock = startingGlobalLock; lock < globalUpperBound; lock++) {
            differenceCover.push_back(lock);

            //TODO update this (fully)
            //this is for if there is nothing to change after assigning the lock
            if (differenceCover.size() == dSize) {
                if (checkWrite(differenceCover, testCover)) {
                    groupid.pop_back();
                    groupNodes.pop_back();

                    return 1;
                }

                continue;
            }

            if (localdSize - 2 > 1 && (perfect = perfect && lock == startingGlobalLock)) { //TODO make sure that this is correct. I want it to be that the lock is the perfect reflection of the starting third
                int numNum = lock - (localdSize - 2) - localThird - 1;
                unsigned long long startingLock = localThird + 1;
                unsigned long long iters = 0;

                if (numNum == 0) {
                    continue;
                }

                calcBounds(numNum, iters, startingLock);

                for (int i = startingLock + iters; i >= startingLock; i--) {
                    if (recursiveLock(differenceCover, testCover, lock, localdSize - 2, i, starting) || check()) { //added check here in case none of the recursives have the time to check
                        groupid.pop_back();
                        groupNodes.pop_back();

                        return 1;
                    }
                }

                //shouldn't need any id reset because the preceeding layer handles

                starting.erase(starting.begin() + starting.size() - 1);
                continue;
            }

            int localStartLock = localThird + 1;

            if (perfect && localThird < int(p / 2) + 1 && localdSize - 2 == 1) {
                localStartLock = int(p / 2) + 1;
            }

            //calcs all this based on the segment that has yet to be determined for just a plugin (not bad)
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

            if (numNum == 0) {
                continue;
            }

            calcBounds(numNum, instanceStart, startValue);
            unsigned long long upperBound = instanceStart + startValue;

            //this needs to populate differenceCover the rest of the way with values
            vector<int> localStart = kthCombination(startValue, sc, localdSize - 2);
            for(int i = differenceCover.size(), count = 0; i < dSize; i++, count++) {
                differenceCover.push_back(localStart[count]);
            }

            //TODO this might go because of the way choose is changed to
            unsigned long long startingIndex = startValue + 1;

            if (checkWrite(differenceCover, testCover)) {
                groupid.pop_back();
                groupNodes.pop_back();

                return 1;
            }

            //TODO update this to reflect
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

            groupid.pop_back();
            groupNodes.pop_back();
            //TODO this will need to pop all the way back to origin (double check)
            differenceCover.pop_back();
        }
    } else { //this is for half or above
    }
}