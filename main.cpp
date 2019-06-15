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
#include <limits>
#include <mpi.h>
#include <iomanip>

using namespace std;

///////////////////////////////////////////////////GLOBALS///////////////////////////
int p;
int dSize;
vector<int> differenceCover; //this is a vector but it acts as a stack
int *testCover;
string pFile;
string globalFile;

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

void popLayer() {
    //pop();
    groupNodes.pop_back();
    groupid.pop_back();
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

bool check() {
    struct stat b;
    if (stat((pFile + ".txt").c_str(), &b) == 0) {
        return 1;
    }

    return 0;
}

int push(int num) {
    updateTest(num);

    differenceCover.push_back(num);

    //this is to make sure that the change we are making can still yield a cover
    //TODO make sure that this doesn't fuck with adding 0 and 1
    if (isCover()) {
        return 1;
    } else {
        pop(num); //popping here takes a lot of control away
        return 0;
    }
}

int pop() {
    int num = differenceCover.back();

    differenceCover.pop_back();

    undoTest(num);

    return num;
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
    globalFile = outFile;

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


    //TODO might need to do these manually (hopefully shouldn't)
    push(0);
    push(1);

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
            if(searchCovers(startingThird, dSize - 2, true) && testSize() == p) {
                //this is the individual write
                ofstream myfilea;
                myfilea.open((pFile + ".txt").c_str(), ios::trunc);
                for (int a = 0; a < dSize - 1; a++) {
                    myfilea << differenceCover[a] << " ";
                }
                myfilea << differenceCover[dSize - 1] << endl;
                myfilea.close();


                //this is the all together write
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
            }

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

    if(!push(localThird)) {
        return 0;
    }

    //TODO any returns need to have difference cover and the mpi id stuff pop before hand
    if (localThird < int((p + 1) / 2)) { //this is for below the half //TODO have ending function?
        //this is all init setup for splitting
        unsigned long long startingGlobalLock = p + 1 - localThird;
        unsigned long long globalIters = 0;
        int numGlobalNum = localp - startingGlobalLock;

        calcBounds(numGlobalNum, globalIters, startingGlobalLock);
        unsigned long long globalUpperBound = startingGlobalLock + globalIters;

        localdSize -= 2;

        for (unsigned long long lock = startingGlobalLock; lock < globalUpperBound; lock++) {
            if(!push(lock)) {
                continue;
            } else if(differenceCover.size() == dSize) {
                popLayer();

                return 1;
            }

            if (localdSize > 1 && (perfect = perfect && lock == startingGlobalLock)) { //TODO make sure that this is correct. I want it to be that the lock is the perfect reflection of the starting third
                int numNum = lock - (localdSize) - localThird - 1;
                unsigned long long startingLock = localThird + 1;
                unsigned long long iters = 0;

                if (numNum == 0) {
                    pop();
                    continue;
                }

                calcBounds(numNum, iters, startingLock);

                for (int i = startingLock + iters; i >= startingLock; i--) {
                    if (searchCovers(localdSize, i, perfect) || check()) { //added check here in case none of the recursives have the time to check
                        popLayer();

                        return 1;
                    }
                }

                //shouldn't need any id reset because the preceeding layer handles

                pop();
                continue;
            }

            if (check()) { //added check here in case none of the recursives have the time to check
                popLayer();

                return 1;
            }

            int localStartLock = localThird + 1;

            if (perfect && localThird < int(p / 2) + 1 && localdSize == 1) {
                localStartLock = int(p / 2) + 1;
            }

            //calcs all this based on the segment that has yet to be determined for just a plugin (not bad)
            vector<int> sc; //this is all the possible numbers between the startingThird and the lock
            sc.reserve(lock - localStartLock);

            for (int i = localStartLock; i < lock; i++) {
                sc.push_back(i);
            }

            //this condition asks if there are more spots than numbers to choose from
            if (sc.size() < localdSize) {
                continue; //this needs to be continue because lock icreases each loop so this condition might not be true next time around
            }

            unsigned long long numNum = nChoosek(sc.size(), localdSize);
            unsigned long long startValue = 0;
            unsigned long long instanceStart = 0;

            if (numNum == 0) {
                continue;
            }

            calcBounds(numNum, instanceStart, startValue);
            unsigned long long upperBound = instanceStart + startValue;

            //TODO add the reassignments for the check (might need to get rid of the check and do manually for these)
            /*
             * thinking that you could use generate cover after manually adding one less than the first value and updating test because it will immidiately try out you number by adding one
             * then it will procees to go through the rest of the checks
             * this wouldn't allow the accuracy of the kth combo but it would make stuff a lot easier at the cost of a couple hundred or thousand of covers
             * it would also probably be faster
             */
            //this needs to populate differenceCover the rest of the way with values
            vector<int> localStart = kthCombination(startValue, sc, localdSize);

            differenceCover.push_back(localStart.front());
            updateTest(differenceCover.back());

            vector<int> endCover = kthCombination(upperBound, sc, localdSize);

            if(generateCover(differenceCover[differenceCover.size() - 1 - localdSize], dSize-localStart.size(), endCover.front()+1)) {
                popLayer();

                return 1;
            }

            popLayer();

            //TODO this will need to pop all the way back to origin (double check) (might need to pop more off), but does it really matter? This branch would have been explored already
            //this is the pop back for any other number that needed to be filled
            while(differenceCover.size() > dSize-localStart.size()) {
                pop();
            }

            pop(); //this is the popback for the lock
        }

    } else { //this is for half or above
        localdSize--;

        if(!push(localThird)) {
            return 0;
        }

        if (localp - localThird - 1 < localdSize) {
            return 0; //we return here because there is no hope for change
        }

        unsigned long long numNum = nChoosek(localp - localThird - 1, localdSize);
        unsigned long long startValue = 0;
        unsigned long long instanceStart = 0;

        calcBounds(numNum, instanceStart, startValue);
        unsigned long long upperBound = startValue + instanceStart;

        vector<int> sc;
        sc.reserve(localp - (localThird + 1));

        for (int i = localThird + 1; i < localp; i++) {
            sc.push_back(i);
        }

        vector<int> localStart = kthCombination(startValue, sc, localdSize);


        differenceCover.push_back(localStart.front());
        updateTest(differenceCover.back());

        vector<int> endCover = kthCombination(upperBound, sc, localdSize);

        if(generateCover(differenceCover[differenceCover.size() - 1 - localdSize], dSize-localStart.size(), endCover.front()+1)) {
            popLayer();

            return 1;
        }

        popLayer();
    }

    //this is the popback for the starting third
    pop();
}

/*(
 * the expectation is that you run this
 * it will increment the las number and see if it works
 * when it doesn't work it starts popping off
 * eventually this is done leaving only what is in bounds
 * TODO test
 */
int generateCover(int localp, int minSize, int stop) {
    int pval = pop();
    for(int i = 1; differenceCover.back() < localp - (dSize - differenceCover.size())-1; i++) { //the -1 should make it so that when another one is added it is still under the pval the other piece is to adjust localp for position (startingThird max value of 3 overall
        if(push(pval+i)) {
            return 1;
        }
    }

    int lastPop = 0;
    //we can't pop off the last number because it has to increment to build back up again
    while (differenceCover.back() >= localp - (dSize - differenceCover.size())) { //the one before the edit space
        if(differenceCover.size() > minSize) {
            lastPop = pop();
        } else {
            return 0;
        }
    }

    push(lastPop); //there shouldn't be an issue with pushing this again because it was certified once and that is how it got onto the stack

    while(differenceCover.size() < dSize) {
        if(!generateCover(localp, minSize, stop)) {
            return 0;
        }
    }

    return 1;
}

void updateTest(int num) {
    for (int i = 0; i < differenceCover.size(); i++) {
        int q = 0;
        if (num < differenceCover[i]) {
            q = differenceCover[i] - num;
        } else {
            q = num - differenceCover[i];
        }

        testCover[q]++;
        testCover[p - q]++;
    }

    //TODO might not need because the number wouldn't have been added yet
    testCover[0] = 1;
}

//exact same as isCover just the reverse
void undoTest(int num) {
    for (int i = 0; i < differenceCover.size(); i++) {
        int q = 0;
        if (num < differenceCover[i]) {
            q = differenceCover[i] - num;
        } else {
            q = num - differenceCover[i];
        }

        testCover[q]--;
        testCover[p - q]--;
    }

    //TODO might not need because the number would not have been added yet
    testCover[0] = 1;
}

int isCover() {
    //TODO make sure this actually says if you got a whole cover
    if ((p - 1) - differenceCover.size() * 2 <= testSize()) {
        return 1;
    }

    return 0;
}

int testSize() {
    int size = 0;
    for (int i = 0; i < p; i++) {
        if (testCover[i] != 0) {
            size++;
        }
    }

    return size;
}