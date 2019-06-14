//
// Created by jstigter on 6/10/19.
//

#ifndef DIFFERENCE_COVER_MAIN_H
#define DIFFERENCE_COVER_MAIN_H


void updateTest();

void isCover();

void undoTest();
void quickCheck();

void nextCover();

int searchCovers(int localThird, int localdSize, bool perfect);
void startSearch();

void calcBounds(unsigned long long numNum, unsigned long long &iters, unsigned long long &starting);

unsigned long long gcd(unsigned long long x, unsigned long long y);

unsigned long long nChoosek(unsigned long long n, unsigned long long k);

#endif //DIFFERENCE_COVER_MAIN_H
