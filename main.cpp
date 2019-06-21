//
// Created by jstigter on 6/10/19.
//
//TODO add all the mpi stuff back
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
#include <math.h>

using namespace std;

///////////////////////////////////////////////////GLOBALS///////////////////////////
int p;
int dSize;
vector<int> differenceCover; //this is a vector but it acts as a stack
vector<int> testCover;
string pFile;
string globalFile;

///////////////////mpi
vector<int> groupNodes; //this is a stack
vector<int> groupid;    //this is a stack
int id;                 //the id of this process
int ierr;
int nn; //number of connected nodes

///////////////////////////////////////////////////RESOURCE FUNCTIONS///////////////////////////
namespace patch
{
template <typename T>
std::string stringMaker(const T &n)
{
    std::ostringstream stm;
    stm << n;
    return stm.str();
}
} // namespace patch

void calcBounds(unsigned long long numNum, unsigned long long &iters, unsigned long long &starting)
{
    int localNodes = groupNodes.back();
    int localId = groupid.back();

    int preGroup = localNodes;
    int preId = localId;

    localNodes = int(localNodes / numNum);

    if (localNodes == 0)
    {
        localNodes = 1;
    }

    if (localId >= numNum)
    {
        localId = localId % numNum;
    }

    iters = numNum / preGroup;

    if (localId < preGroup % numNum && preGroup > numNum)
    {
        localNodes++; //+= nn / numNum;
    }

    if (localId < numNum % preGroup)
    {
        iters += 1;

        if (localId != 0)
        {
            starting += localId * iters;
        }
    }
    else
    {
        starting += numNum % preGroup + localId * iters;
    }

    if (localNodes == 1)
    {
        localId = 0;
    }
    else
    {
        localId = preId / numNum;
    }

    groupNodes.push_back(localNodes);
    groupid.push_back(localId);
}

void popLayer()
{
    groupNodes.pop_back();
    groupid.pop_back();
}

unsigned long long gcd(unsigned long long x, unsigned long long y)
{
    while (y != 0)
    {
        unsigned long long t = x % y;
        x = y;
        y = t;
    }
    return x;
}

unsigned long long nChoosek(unsigned long long n, unsigned long long k)
{
    if (k > n)
        throw std::invalid_argument("invalid argument in nChoosek");
    unsigned long long r = 1;
    for (unsigned long long d = 1; d <= k; ++d, --n)
    {
        unsigned long long g = gcd(r, d);
        r /= g;
        unsigned long long t = n / (d / g);
        if (r > std::numeric_limits<unsigned long long>::max() / t)
            throw std::overflow_error("overflow in nChoosek");
        r *= t;
    }
    return r;
}

vector<int> kthCombination(unsigned long long k, vector<int> l, int r)
{
    if (r == 0)
    {
        vector<int> ans;
        return ans;
    }
    else if (l.size() == r)
    {
        return l;
    }
    else
    {
        unsigned long long i = nChoosek(l.size() - 1, r - 1); //calculare number of combinations
        if (k < i)
        {
            vector<int> ans;
            ans.push_back(l[0]);

            vector<int> tertiary(l.begin() + 1, l.end());

            vector<int> secondary = kthCombination(k, tertiary, r - 1);

            ans.insert(ans.end(), secondary.begin(), secondary.end());

            return ans;
        }
        else
        {
            vector<int> ans;

            vector<int> tertiary(l.begin() + 1, l.end());

            vector<int> secondary = kthCombination(k - i, tertiary, r);

            return secondary;
        }
    }
}

bool check()
{
    struct stat b;
    if (stat((pFile + ".txt").c_str(), &b) == 0)
    {
        return 1;
    }

    return 0;
}

int push(int num)
{
    updateTest(num);
    differenceCover.push_back(num);

    if (isCover())
    {
        return 1;
    }
    else
    {
        pop(); //popping here takes a lot of control away
        return 0;
    }
}

int pop()
{
    int num = differenceCover.back();

    differenceCover.pop_back();

    undoTest(num);

    return num;
}

/////////////////////////////
int main(int argc, char *argv[])
{
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &nn);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &id);

    if (argc == 1)
    {
        cout << "This program computes difference covers of a range specified by" << endl;
        cout << "the user. This program then saves the covers to a file that you" << endl;
        cout << "specify. The following format is used to execute this command." << endl;
    }
    if (argc < 3)
    {
        cout << "Usage: mhonors <filename> <startValue> [numberToCompute]" << endl
             << endl;
        cout << "<filename> - enter a filename." << endl
             << endl;
        cout << "<startValue> - enter an integer number. This is the first cover" << endl;
        cout << "               the program will attempt to find." << endl
             << endl;
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

    testCover = vector<int>(startValue + numberToCompute-1, 0);

    while (p < startValue + numberToCompute)
    {
        pFile = patch::stringMaker(p);

        if (check())
        {
            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " found " << p << " before starting" << endl;
            p++;
            continue;
        }

        MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave

        cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " on: " << p << endl;

        for (int i = 0; i < p; i++)
        {
            testCover[i] = 0;
        }

        startSearch();

        p++;
    }

    cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " is finished." << endl;

    MPI_Finalize();

    cout << "mpi was finalized" << endl;

    return 0;
} // end main

void startSearch()
{
    //calcs min
    int min = 1;
    for (int x = 2; x < p; x++)
    {
        if (x * (x - 1) + 1 >= p)
        {
            min = x;
            break;
        }
    }

    //calcs max
    int max = p;
    for (int x = 2; x < p / 2; x++)
    {
        if (x * x > p)
        {
            max = 3 * (p + x - 1) / (2 * (x - 1)) + 4;
            break;
        }
    }

    differenceCover.reserve(max);

    dSize = min;

    unsigned long long totalCombo = nChoosek(p - 2, dSize - 2) / 2;
    int startingThird = int((p + 1) / 2) + 1;

    struct stat buffer;
    if (stat((pFile + "_0.txt").c_str(), &buffer) == 0)
    {
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

        cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " read " << dSize << " " << startingThird << endl;
    }

    push(0);
    push(1);

    while (dSize <= max)
    {
        cout << "dSize is now " << dSize << endl;
        for (int i = startingThird; i > 1; i--)
        {
            if (id == 0)
            {
                cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " is writing " << dSize << " " << i << endl;
                ofstream myfile;
                myfile.open((pFile + "_0.txt").c_str(), ios::trunc);
                myfile << dSize << endl;
                myfile << i << endl;
                myfile.close();
            }

            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " the new starting third is " << i << endl;

            if (searchCovers(i, dSize - 2, true) && testSize() == p)
            {
                cout << "winning cover is: ";
                {
                    cout << "differenceCover: ";
                    for (int i = 0; i < differenceCover.size(); i++)
                    {
                        cout << differenceCover[i] << " ";
                    }

                    cout << endl;
                }

                {
                    cout << "testCover: ";
                    for (int i = 0; i < p; i++)
                    {
                        if (testCover[i] != 0)
                        {
                            cout << "(" << testCover[i] << "," << i << ")";
                        }
                    }
                    cout << endl;
                }

                //this is the individual write
                ofstream indivOut;
                indivOut.open((pFile + ".txt").c_str(), ios::trunc);
                for (int i = 0; i < differenceCover.size() - 1; i++)
                {
                    indivOut << differenceCover[i] << "   ";
                }
                indivOut << differenceCover.back();
                indivOut << endl;
                indivOut.flush();
                indivOut.close();

                //this is the all together write
                ofstream out;
                out.open(globalFile.c_str(), ios::app);
                out << p;
                out << "    " << differenceCover.size();
                out << "            ";
                for (int i = 0; i < differenceCover.size() - 1; i++)
                {
                    out << differenceCover[i] << "   ";
                }
                out << differenceCover.back();
                out << endl;
                out.flush();
                out.close();
            }

            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " is done checking startingThird: " << i << endl;

            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " is waiting for the rest" << endl;
            MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave
            if (check())
            {
                cout << "we are returning: " << endl;

                for (int x = 0; x < differenceCover.size(); x++)
                {
                    cout << differenceCover[x] << " ";
                }
                cout << endl;

                return;
            }
        }

        startingThird = int((p + 1) / 2) + 1;
        dSize++;
    }
}

int searchCovers(int localThird, int localdSize, bool perfect)
{
    //TODO smooth this out and get rid of it (this is a good temp fix)
    int localp = differenceCover.back();

    if (localp == 1)
    {
        localp = p;
    }

    if (!push(localThird))
    {
        return 0;
    }

    //TODO any returns need to have difference cover and the mpi id stuff pop before hand
    if (localThird < int((p + 1) / 2))
    { //this is for below the half //TODO have ending function?
        cout << "we are below the half" << endl;
        //this is all init setup for splitting
        unsigned long long startingGlobalLock = p + 1 - localThird;
        unsigned long long globalIters = 0;
        int numGlobalNum = localp - startingGlobalLock;

        calcBounds(numGlobalNum, globalIters, startingGlobalLock);

        unsigned long long globalUpperBound = startingGlobalLock + globalIters;

        localdSize -= 2;

        cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " starting lock is: " << startingGlobalLock << " globalUpperBound is: " << globalUpperBound << endl;
        for (unsigned long long lock = startingGlobalLock; lock < globalUpperBound; lock++)
        {
            if (!push(lock))
            {
                continue;
            }
            else if (differenceCover.size() == dSize)
            {
                popLayer();

                return 1;
            }

            if (localdSize > 1 && (perfect = (perfect && lock == p + 1 - localThird)))
            { //TODO make sure that this is correct. I want it to be that the lock is the perfect reflection of the starting third
                int numNum = lock - (localdSize)-localThird - 1;
                unsigned long long startingLock = localThird + 1;
                unsigned long long iters = 0;

                if (numNum == 0)
                {
                    pop();
                    continue;
                }

                calcBounds(numNum, iters, startingLock);

                for (int i = startingLock + iters; i >= startingLock; i--)
                {
                    cout << "going another layer deep " << i << " " << localdSize << " perfect: " << perfect << " actual " << lock << " v " << p + 1 - localThird << endl;
                    if (searchCovers(i, localdSize, perfect) || check())
                    {               //added check here in case none of the recursives have the time to check //TODO add back || check()
                        popLayer(); //this is to get rid of the recursive
                        popLayer(); //this is to get rid of this functions

                        return 1;
                    }
                }

                popLayer();
                pop(); //gets rid of the lock
                continue;
            }

            if (check())
            { //added check here in case none of the recursives have the time to check
                cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " someone else found something" << endl;
                popLayer();

                return 1;
            }

            int localStartLock = localThird + 1;

            if (perfect && localThird < int(p / 2) + 1 && localdSize == 1)
            {
                localStartLock = int(p / 2) + 1;
            }

            //calcs all this based on the segment that has yet to be determined for just a plugin (not bad)
            vector<int> sc; //this is all the possible numbers between the startingThird and the lock
            sc.reserve(lock - localStartLock);

            for (int i = localStartLock; i < lock; i++)
            {
                sc.push_back(i);
            }

            //this condition asks if there are more spots than numbers to choose from
            if (sc.size() < localdSize)
            {
                pop();
                continue; //this needs to be continue because lock icreases each loop so this condition might not be true next time around
            }

            /*unsigned long long numNum = nChoosek(sc.size(), localdSize);
            unsigned long long startValue = 0;
            unsigned long long instanceStart = 0;

            if (numNum == 0)
            {
                pop();
                continue;
            }

            calcBounds(numNum, instanceStart, startValue);*/

            unsigned long long numNum = differenceCover.back()-localdSize-differenceCover[differenceCover.size()-2]-1;
            unsigned long long startValue = differenceCover[differenceCover.size()-2]+1;
            unsigned long long instanceStart = 0;

            if (numNum == 0)
            {
                pop();
                continue;
            }

            calcBounds(numNum, instanceStart, startValue);
            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " " << numNum << " " << instanceStart << " " << startValue << endl;
            unsigned long long upperBound = instanceStart + startValue;

            //TODO add the reassignments for the check (might need to get rid of the check and do manually for these)
            /*
             * thinking that you could use generate cover after manually adding one less than the first value and updating test because it will immidiately try out you number by adding one
             * then it will procees to go through the rest of the checks
             * this wouldn't allow the accuracy of the kth combo but it would make stuff a lot easier at the cost of a couple hundred or thousand of covers
             * it would also probably be faster
             */
            //this needs to populate differenceCover the rest of the way with values
            //vector<int> localStart = kthCombination(startValue, sc, localdSize);

            differenceCover.push_back(startValue-1);
            updateTest(differenceCover.back());

            //TODO update this to pass the entire vector because it is more percise and wont have a large impact
            //vector<int> endCover = kthCombination(upperBound, sc, localdSize);

            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " starting is: " << differenceCover.back()+1 << " ending is: " <<  upperBound << " | ";

            for (int i = 0; i < differenceCover.size(); i++)
            {
                cout << " " << differenceCover[i];
            }
            cout << endl;

            if (generateCover(differenceCover[differenceCover.size() - 2], dSize - localdSize, upperBound))
            { //the lock becomes the localp
                cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " we found something" << endl;

                popLayer();

                return 1;
            }

            popLayer(); 

            //TODO this will need to pop all the way back to origin (double check) (might need to pop more off), but does it really matter? This branch would have been explored already
            //this is the pop back for any other number that needed to be filled
            while (differenceCover.size() > dSize - localdSize - 1)
            { //teh minux 1 is to get rid of the lock
                pop();
            }

            //pop(); //this is the popback for the lock
        }
    }
    else
    { //this is for half or above
        cout << "we are above the half" << endl;
        localdSize--;

        if (localp - localThird - 1 < localdSize)
        {
            return 0; //we return here because there is no hope for change
        }

        /*unsigned long long numNum = nChoosek(localp - localThird - 1, localdSize);
        unsigned long long startValue = 0;
        unsigned long long instanceStart = 0;

        calcBounds(numNum, instanceStart, startValue);*/
        unsigned long long numNum = localp - localdSize-differenceCover.back()-1;
        unsigned long long startValue = differenceCover.back()+1; //needs to be the same so it can increment to +1
        unsigned long long instanceStart = 0;

        if(numNum == 0) {
            pop();
            return 0;
        }

        calcBounds(numNum, instanceStart, startValue);
            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " " << numNum << " " << instanceStart << " " << startValue << endl;

        unsigned long long upperBound = startValue + instanceStart;

        /*vector<int> sc;
        sc.reserve(localp - (localThird + 1));

        for (int i = localThird + 1; i < localp; i++)
        {
            sc.push_back(i);
        }

        vector<int> localStart = kthCombination(startValue, sc, localdSize);*/

        differenceCover.push_back(startValue-1); // the minus one is so that it is immidiately incremented to be normal
        updateTest(differenceCover.back());

        //vector<int> endCover = kthCombination(upperBound, sc, localdSize);

        cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " starting is (+): " << differenceCover.back()+1 << " ending is: " << upperBound << " | ";

        for (int i = 0; i < differenceCover.size(); i++)
        {
            cout << " " << differenceCover[i];
        }
        cout << endl;

        //you need to adjust for having numbers at the end which affects the localp value max for each position
        if (generateCover(localp, dSize - localdSize, upperBound))
        { //TODO make sure that setting the localp to p doesn't screw the recursion
            popLayer();
            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " we found something (above half)" << endl;

            return 1;
        }

        while (differenceCover.size() > dSize - localdSize)
        {
            pop();
        }
    }

    //this is the popback for the starting third
    pop();
    popLayer();

    return 0;
}

/*(
 * the expectation is that you run this
 * it will increment the las number and see if it works
 * when it doesn't work it starts popping off
 * eventually this is done leaving only what is in bounds
 * TODO test
 */
int generateCover(int localp, int minSize, int stop)
{
    while (differenceCover.size() < dSize)
    {
        /*cout << "Generating Cover: ";
        {
            cout << "differenceCover: ";
            for (int i = 0; i < differenceCover.size(); i++)
            {
                cout << differenceCover[i] << " ";
            }

            cout << endl;
        }*/

        if (differenceCover.size() <= minSize)
        {
            return 0;
        }

        int pval = pop();

        //this is the part that fills the cover because it doesn't break after pushing it keeps boing getting higher and higher
        for (int i = 1; pval + i <= localp - (dSize - (differenceCover.size() + 1)) - 1; i++)
        { //the -1 should make it so that when another one is added it is still under the pval the other piece is to adjust localp for position (startingThird max value of 3 overall
            if (push(pval + i))
            { //this is to make sure that this function does not exit before filling everything
                if (differenceCover.size() == dSize)
                {
                    return 1;
                }
                else if (differenceCover.size() == minSize + 1 && differenceCover.back() == stop)
                {
                    return 0;
                }
            }
        }

        int lastPop = 0;
        //we can't pop off the last number because it has to increment to build back up again
        while (differenceCover.back() > localp - (dSize - differenceCover.size()) - 1)
        { //the one before the edit space
            if (differenceCover.size() - 1 >= minSize)
            { //need sto be less than because you add one back
                lastPop = pop();
            }
            else
            {
                return 0;
            }
        }

        if (differenceCover.size() == minSize)
        {
            return 0;
        }

        if (lastPop != 0)
        {
            push(lastPop); //there shouldn't be an issue with pushing this again because it was certified once and that is how it got onto the stack
        }
    }
    return 1;
}

void updateTest(int num)
{
    for (int i = 0; i < differenceCover.size(); i++)
    {
        int q = 0;
        if (num < differenceCover[i])
        {
            q = differenceCover[i] - num;
        }
        else
        {
            q = num - differenceCover[i];
        }

        testCover[q]++;
        testCover[p - q]++;
    }

    //TODO might not need because the number wouldn't have been added yet
    testCover[0] = 1;
}

//exact same as isCover just the reverse
void undoTest(int num)
{
    for (int i = 0; i < differenceCover.size(); i++)
    {
        int q = 0;
        if (num < differenceCover[i])
        {
            q = differenceCover[i] - num;
        }
        else
        {
            q = num - differenceCover[i];
        }

        testCover[q]--;
        testCover[p - q]--;
    }

    //TODO might not need because the number would not have been added yet
    testCover[0] = 1;
}

int isCover()
{
    int dif = dSize - differenceCover.size();
    //TODO make sure this actually says if you got a whole cover
    /*{
        cout << "differenceCover: ";
        for (int i = 0; i < differenceCover.size(); i++) {
            cout << differenceCover[i] << " ";
        }

        cout << endl;
    }

    {
        cout << "testCover: ";
        for (int i = 0; i < p; i++) {
            if (testCover[i] != 0) {
                cout << "(" << testCover[i] << "," << i << ")";
            }
        }
        cout << endl;
    }*/

    //TODO this might still bewrong
    //got rid of the -1 for p because you get zero for free in the count so you would get the discount twice with the -1
    //cout << "dif: " << dif << " need: " << (p - 1) - (pow(dif, 2) + dif + dif*2*(differenceCover.size()-1)) << " testSize: " << testSize() << " " << differenceCover.size() << endl;
    if ((p) - (pow(dif, 2) + dif + dif * 2 * (differenceCover.size() - 1)) <= testSize())
    {
        return 1;
    }

    /*cout << "this cover isn't possible "  << " needed: " << (p - 1) - (pow(dif, 2) + dif + dif*2*(differenceCover.size()-1)) << " had: " << testSize() << " other: " << dif << endl;
    {
        cout << "differenceCover: ";
        int lastMin = -1;
        for(int x = 0; x < differenceCover.size(); x++) {
            int min = 100000;
            for (int i = 0; i < differenceCover.size(); i++) {
                if(min > differenceCover[i] && lastMin < differenceCover[i]) {
                    min = differenceCover[i];
                }
                //cout << differenceCover[i] << " ";
            }

            lastMin = min;
            cout << min << " ";
        }

        cout << "    |     ";
        for(int i = 0; i < differenceCover.size(); i++) {
            cout << differenceCover[i] << " ";
        }

        cout << endl;
    }*/

    /*{
        cout << "testCover: ";
        for (int i = 0; i < p; i++) {
            if (testCover[i] != 0) {
                cout << "(" << testCover[i] << "," << i << ")";
            }
        }
        cout << endl;
    }*/
    return 0;
}

int testSize()
{
    int size = 0;
    for (int i = 0; i < p; i++)
    {
        if (testCover[i] != 0)
        {
            size++;
        }
    }

    return size;
}