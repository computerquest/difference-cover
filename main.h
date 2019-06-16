//
// Created by jstigter on 6/10/19.
//
#include <vector>
using namespace std;
#ifndef DIFFERENCE_COVER_MAIN_H
#define DIFFERENCE_COVER_MAIN_H


void updateTest(int num);

int isCover();

void undoTest(int num);

int push(int num);

int pop();

int testSize();

int searchCovers(int localThird, int localdSize, bool perfect);

void startSearch();

void popLayer();

void calcBounds(unsigned long long numNum, unsigned long long &iters, unsigned long long &starting);

unsigned long long gcd(unsigned long long x, unsigned long long y);

unsigned long long nChoosek(unsigned long long n, unsigned long long k);

int generateCover(int localp, int minSize, int stop);

vector<int> kthCombination(unsigned long long k, vector<int> l, int r);

bool check();

#endif //DIFFERENCE_COVER_MAIN_H
