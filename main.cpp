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
#include <stdlib.h>

using namespace std;

//todo delete
unsigned long long totalCover = 0;
unsigned long long totalCheck = 0;
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

    testCover = vector<int>(startValue + numberToCompute, 0);

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

        MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave

        unsigned long long *allChecked = NULL;
        unsigned long long *allCover = NULL;

        if (id == 0)
        {
            allChecked = new unsigned long long[nn];
            allCover = new unsigned long long[nn];
        }

        MPI_Gather(&totalCheck, 1, MPI::UNSIGNED_LONG_LONG, allChecked, 1, MPI::UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
        MPI_Gather(&totalCover, 1, MPI::UNSIGNED_LONG_LONG, allCover, 1, MPI::UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

        if (id == 0)
        {
            cout << "TOTALS: " << endl;
            unsigned long long tCheck = 0;
            cout << "checked:";
            for (int i = 0; i < nn; i++)
            {
                cout << " " << allChecked[i];
                tCheck += allChecked[i];
            }
            cout << " total: " << tCheck << endl;

            unsigned long long tCover = 0;
            cout << "cover:";
            for (int i = 0; i < nn; i++)
            {
                cout << " " << allCover[i];
                tCover += allCover[i];
            }
            cout << " total: " << tCover << endl;

            delete[] allChecked;
            delete[] allCover;
        }

        totalCheck = 0;
        totalCover = 0;

        p++;
    }

    cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " is finished." << endl;

    MPI_Finalize();

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

    differenceCover = vector<int>();
    differenceCover.reserve(max);

    dSize = min;

    //todo do i need this line?????
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

    differenceCover.push_back(0);
    updateTest(0);
    push(1);

    while (dSize <= max)
    {
        if (id == 0)
        {
            cout << "dSize is now " << dSize << endl;
        }
        for (int i = startingThird; i > 1; i--)
        {
            if (id == 0)
            {
                cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " the new starting third is " << i << endl;
                ofstream myfile;
                myfile.open((pFile + "_0.txt").c_str(), ios::trunc);
                myfile << dSize << endl;
                myfile << i << endl;
                myfile.close();
            }

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

                //todo add back
                //this is the individual write
                /*ofstream indivOut;
                indivOut.open((pFile + ".txt").c_str(), ios::trunc);
                for (int i = 0; i < differenceCover.size() - 1; i++)
                {
                    indivOut << differenceCover[i] << "   ";
                }
                indivOut << differenceCover.back();
                indivOut << endl;
                indivOut.flush();
                indivOut.close();

                {
                    //this is the all together write
                    ofstream out;
                    out.open(globalFile.c_str(), ios::app);
                    out << p;
                    out << "    " << differenceCover.size();
                    out << "            ";

                    {
                        int lastMin = -1;
                        for(int x = 0; x < differenceCover.size(); x++) {
                            int min = 100000;
                            for (int i = 0; i < differenceCover.size(); i++) {
                                if(min > differenceCover[i] && lastMin < differenceCover[i]) {
                                    min = differenceCover[i];
                                }
                            }
                            lastMin = min;
                            out << min << "   ";
                        }
                    }
                    out << endl;
                    out.flush();
                    out.close();
                }*/
            }

            cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " is done checking startingThird: " << i << endl;

            MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave
            //todo add back
            /*if (check())
            {
                cout << "we are returning: " << endl;

                for (int x = 0; x < differenceCover.size(); x++)
                {
                    cout << differenceCover[x] << " ";
                }
                cout << endl;

                return;
            }*/
        }

        MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave

        if (check())
        {
            return;
        }

        MPI_Barrier(MPI_COMM_WORLD); //this is to sync all the processes for the next wave

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
    else if (dSize == differenceCover.size())
    {
        {
            //this is the all together write
            ofstream out;
            out.open(globalFile.c_str(), ios::app);
            out << p;
            out << "    " << differenceCover.size();
            out << "            ";

            {
                int lastMin = -1;
                for (int x = 0; x < differenceCover.size(); x++)
                {
                    int min = 100000;
                    for (int i = 0; i < differenceCover.size(); i++)
                    {
                        if (min > differenceCover[i] && lastMin < differenceCover[i])
                        {
                            min = differenceCover[i];
                        }
                    }
                    lastMin = min;
                    out << min << "   ";
                }
            }
            out << endl;
            out.flush();
            out.close();
        }

        //todo remove
        pop();
        return 0; //todo change this to return 1;
    }

    localdSize--;

    //TODO any returns need to have difference cover and the mpi id stuff pop before hand
    if (localThird < int((p + 1) / 2))
    { //this is for below the half //TODO have ending function?
        //this is all init setup for splitting
        unsigned long long startingGlobalLock = p + 1 - localThird;
        unsigned long long globalIters = 0;
        int numGlobalNum = localp - startingGlobalLock;

        if (numGlobalNum <= 0)
        {
            return 0;
        }

        calcBounds(numGlobalNum, globalIters, startingGlobalLock);

        unsigned long long globalUpperBound = startingGlobalLock + globalIters;

        localdSize--;

        //cout << "Thread: " << id << " gid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " the lock is start: " << startingGlobalLock << " end: " << globalUpperBound << " num: " << numGlobalNum << " localp: " << localp << " third: " << localThird << endl;
        for (unsigned long long lock = startingGlobalLock; lock < globalUpperBound; lock++)
        {
            if (!push(lock))
            {
                continue;
            }
            else if (differenceCover.size() == dSize)
            {

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

                {
                    //this is the all together write
                    ofstream out;
                    out.open(globalFile.c_str(), ios::app);
                    out << p;
                    out << "    " << differenceCover.size();
                    out << "            ";

                    {
                        int lastMin = -1;
                        for (int x = 0; x < differenceCover.size(); x++)
                        {
                            int min = 100000;
                            for (int i = 0; i < differenceCover.size(); i++)
                            {
                                if (min > differenceCover[i] && lastMin < differenceCover[i])
                                {
                                    min = differenceCover[i];
                                }
                            }
                            lastMin = min;
                            out << min << "   ";
                        }
                    }
                    out << endl;
                    out.flush();
                    out.close();
                }

                pop();
                continue;
                //todo replace
                //                pop();
                //                popLayer();

                //return 0;
            }
            //localdSize > 1 &&
            if (localdSize > 1 & (perfect = (perfect && lock == p + 1 - localThird)))
            { //TODO make sure that this is correct. I want it to be that the lock is the perfect reflection of the starting third
                int numNum = lock - (localdSize)-localThird;
                unsigned long long startingLock = localThird + 1;
                unsigned long long iters = 0;

                if (numNum <= 0)
                {
                    pop();
                    continue;
                }

                calcBounds(numNum, iters, startingLock);

                //cout << "Thread: " << id << " gid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " the recursive bounds are start: " << startingLock << " end: " << startingLock + iters << " num: " << numNum << endl;
                for (unsigned long long i = startingLock; i < startingLock + iters; i++)
                {
                    if (searchCovers(i, localdSize, perfect)) //TODO add back || check())
                    {                                         //added check here in case none of the recursives have the time to check //TODO add back || check()
                        popLayer();                           //this is to get rid of the recursive
                        popLayer();                           //this is to get rid of this functions

                        return 1;
                    }
                }

                popLayer();
                pop(); //gets rid of the lock
                continue;
            }

            //TODO add back
            /*if (check())
            { //added check here in case none of the recursives have the time to check
                cout << "Thread: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " someone else found something" << endl;
                popLayer();

                return 1;
            }*/

            int localStartLock = differenceCover[differenceCover.size() - 2];

            if (perfect && localThird < int(p / 2) + 1 && localdSize == 1)
            {
                localStartLock = int(p / 2) + 1;
            }

            exhaustiveSearch(localStartLock - 1, differenceCover.back(), localdSize);

            pop(); //this is the popback for the lock
        }

        popLayer();
    }
    else
    { //this is for half or above
        exhaustiveSearch(localThird, localp, localdSize);
    }

    //this is the popback for the starting third
    pop();

    return 0;
}

//TODO fix this it is causing a lot of covers to go missing (100% sure this is cause no matter what)
int exhaustiveSearch(int floor, int localp, int localdSize)
{
    if (localp - floor - 1 < localdSize)
    {
        return 0; //we return here because there is no hope for change
    }

    unsigned long long numNum = localp - localdSize - floor; //-1;
    unsigned long long startValue = floor + 1;               //needs to be the same so it can increment to +1
    unsigned long long instanceStart = 0;

    if (numNum <= 0)
    {
        pop();
        return 0;
    }

    calcBounds(numNum, instanceStart, startValue);

    unsigned long long upperBound = startValue + instanceStart;

    //TODO switch order
    updateTest(startValue - 1);
    differenceCover.push_back(startValue - 1); // the minus one is so that it is immidiately incremented to be normal

    //TODO add back
    //you need to adjust for having numbers at the end which affects the localp value max for each position

    /*cout << "Thead: " << id << " groupid: " << groupid.back() << " groupNodes: " << groupNodes.back() << " start: " << differenceCover.back() + 1 << " end: " << upperBound << " |";
    for (int i = 0; i < differenceCover.size(); i++)
    {
        cout << " " << differenceCover[i];
    }
    cout << endl;*/

    if (generateCover(localp, dSize - localdSize, upperBound)) //todo is the minsize differenceCover.size()-1?
    {
        //popLayer();

        //return 1;
    }

    while (differenceCover.size() > dSize - localdSize)
    {
        pop();
    }

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
    //todo remove redundant staements
    do
    {
        //todo needed ?
        if (differenceCover.size() <= minSize)
        {
            return 0;
        }

        int pval = pop();

        //this is the part that fills the cover because it doesn't break after pushing it keeps boing getting higher and higher
        for (int i = 1; pval + i <= localp - (dSize - (differenceCover.size() + 1)) - 1; i++)
        { //the -1 should make it so that when another one is added it is still under the pval the other piece is to adjust localp for position (startingThird max value of 3 overall
            if (differenceCover.size() == minSize && pval + i == stop)
            {
                return 0;
            }

            if (push(pval + i) && differenceCover.size() == dSize)
            { //this is to make sure that this function does not exit before filling everything
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

                {
                    //this is the all together write
                    ofstream out;
                    out.open(globalFile.c_str(), ios::app);
                    out << p;
                    out << "    " << differenceCover.size();
                    out << "            ";

                    {
                        int lastMin = -1;
                        for (int x = 0; x < differenceCover.size(); x++)
                        {
                            int min = 100000;
                            for (int i = 0; i < differenceCover.size(); i++)
                            {
                                if (min > differenceCover[i] && lastMin < differenceCover[i])
                                {
                                    min = differenceCover[i];
                                }
                            }
                            lastMin = min;
                            out << min << "   ";
                        }
                    }
                    out << endl;
                    out.flush();
                    out.close();
                }

                pop();

                //return 1;
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
    } while (differenceCover.size() < dSize);

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

        testCover.at(q)++;
        testCover.at(p - q)++;
    }
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

        testCover.at(q)--;
        testCover.at(p - q)--;
    }
}

int isCover()
{
    if (differenceCover.size() == dSize)
    {
        totalCover++;
    }

    totalCheck++;

    int dif = dSize - differenceCover.size();

    //got rid of the -1 for p because you get zero for free in the count so you would get the discount twice with the -1
    if ((p) - (pow(dif, 2) + dif + dif * 2 * (differenceCover.size() - 1)) <= testSize())
    {
        return 1;
    }

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
