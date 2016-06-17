//
//
//                      Microtubule dynamics simulation
//
//
// The goal was to simulate the growth of polymers called microtubules within
// living plant cells, and investigate if there are conditions under which a set
// of randomly growing microtubules might form a stable ordered state.
//
// This code simulates the microtubules as a number of rigid rods. These rods
// grow from the front end at a rate 'vplus' (um/min) and shrink from the rear
// end at a rate 'vmin'. New rods are generated on the screen at a rate of 'rInj'
// in units (um^-2 * sec^-1).
//
// Rod dynamics (my collission / pass-through algorithm):
//   As rods grow, they can interact by either (i) colliding with each other and
//   halting their growth as long as their front end is obstructed by another
//   rod, or (ii) passing through another rod. The probabilities of (i) and (ii)
//   are controlled through variables 'lBundlingProb' and 'rBundlingProb' which
//   set the probabilities of rods passing through other rods if they encounter
//   it from left or right, respectively.
//
// Rod dynamics (bundling):
//   If 'bundling' is true, whenever a rod hits another rod, it will not stop
//   growing. Instead it will continue growing 'alongside' the encountered rod,
//   forming an obtuse angle with its original direction, e.g. \__ .
//
// Initial conditions:
//   Rods are generated in random directions if startOriented == false. If this
//   is true then the rods are generated using Normal distribution with the mean
//
// Boundary conditions:
//   This simulation uses exclusively square periodic boundary bonditions.
//
// Quantification of the order of the entire system:
//   The level of order in this system is quantified using a coarse-grained
//   ordered parameters (see e.g. the statistical physics textbook from Jim
//   Sethna 'Statistical physics, order paramters and complexity'). The entire
//   system is first coarse-grained into a grid of 2^fftPowerOfTwo squares, the
//   order parameter is defined as the average of r*exp(i*theta) in each square.
//   We plot the first and second moment of this order parameter.


static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

#include <cstdio>
#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include "Ran.h"
#include "deviates.h"
#include "erf.h"
using namespace std;

// constants and input arguments
#define PI 3.141592653589793
#define TWOPI 6.283185307179586
#define PIHALF 1.570796326794897
#define DEBUG false // setting this to TRUE will display rod growth numbers in stdout
#define startOriented true
#define bundling true

// global variables
map <unsigned int, vector<vector<double> > > rods;
map <vector<int>, vector<double> > grid;
map <unsigned int, double> phiMoments; // phi[m] = \sum_{k=1...N} L_k exp(i m theta_k)

double screenUnit = 0.033333333333333; // 1 screen unit in um
double boxSize = 10;	  // um; was 10
double winSize = boxSize/screenUnit; // 300 screen unit
double precision = 1e-5;

// Declare functions used in main()
void 	snapshotPs(map <unsigned int, vector<vector<double> > > rods, const char * fileName, unsigned int k);
void	outputParams(bool toFile);
void 	rodSegmentCreate(unsigned int rIdx, double x0, double x1, double y0, double y1, double phi, double status, double clock);
double	rodDistance(vector<double> rodsJ, vector<double> rodsL);
double	rodHitAngle(vector<double> rodsJ, vector<double> rodsL);
void 	rodGrow(unsigned int rIdx, unsigned int m, double len);
void 	rodShrink(unsigned int rIdx, double len);
double	rodLength(unsigned int rIdx);
double  	rodSegLength(unsigned int rIdx, unsigned int sIdx);
double  	rodAngle(unsigned int rIdx, unsigned int sIdx);
double  	rodAngle2(vector<double> rod);
void 	gridInit(unsigned int fftPowerOfTwo, double len);
double	complexAbs(vector<double>);
map<unsigned int, double> orderMoments(unsigned int m); // calculate order phiMoment_m


int main() {
	//
	//  Main function, change simulation parameters here for now, later implement it to load from a file
	//

	// simulation parameters
	double dt = 0.1;			// Time step in min; default 0.1
	unsigned int steps = 10000;	// Number of steps to run the simulation
	unsigned int dsstep = 100;	// Status display step
	unsigned int psstep = 50;	// print .ps file every psstep-th step
	unsigned int pMstep = 50;	// calculate phiMoment every pMstep-th step
	double startTheta = 1.5708;	// starting orientation (if startOriented == true) (1.5708 = PI/2)
	double startThetaSpread = 0.05;	// Gaussian spread around the startTheta.
	unsigned int kSwitch = 3000;	// step at which we switch from ordered to random generation
									//  in case we start oriented=true
	unsigned int phi_m = 2;
	char phi_mFname[23] = "phi_moments.txt";
	double pfix = 1;		// distance between collided rods (to easier recalculate distance
									//  for rechecking when other rod is shrinking)
	unsigned int fftPowerOfTwo = 8;		// 2^fftPowerOfTwo points for FFT

	// open files for writing statistical moments of the order parameter
	ofstream phiMomentsFile;
	phiMomentsFile.open(phi_mFname);

	// model parameters
	double rInj	= 0.0005;			// rate of injenction per sec per um^2 (def. 1)
	double rc = rInj*dt*pow(screenUnit*winSize,2); // Rate of creation of new rods/min/entire area
	double vmin = 0.3;			// Rate of depolymerization of - end (in um/min)
	double vplus = 1.0;			// Rate of polymerization of + end (in um/min)
	double bundlingAngle = PI/6;			// 30 degrees; from [1], p. 110
	double bundlingDistance	= 1;				// arbitrary, has to be small
	double lBundlingProb = 1;				// lBundling
	double rBundlingProb = 0;

	// derived model parameters
	double lContinueProb = 1-lBundlingProb; // probability of "going through" from left
	double rContinueProb = 1-rBundlingProb;

	// [1] D. W. Ehrhardt, \Straighten up and
	// y right-microtubule dynamics and organization of non-
	// centrosomal arrays in higher plants", Curr. Opinion Cell Biol. 20 (Feb 2008), 107-116
	// 300 screen units = 10 um
	// rate of injection = 100 /s/um^2

	//    um/min [Walker 1998]
	// storage arrays
	// using vector of vectors (resizable array)
	// storage of rods' coords, direction and status (-1 del, 0 not growing, 1 growing)
	map <unsigned int, vector<unsigned int> > inters;
	map <unsigned int, unsigned int> rodsB;	// for each rod (key), write which one blocked it (value)

	// Initialize random number generator
	int randSeed = 5345;
	Ran random(randSeed);
	Normaldev ranGauss(startTheta, startThetaSpread, randSeed);
	//srand(randSeed); // debug random seed

	double dr   = vplus*dt;

	// make first rod
	double rnd = random.doub();
	double phi, x1, y1, x2, y2, rndLR, bundlingNewAngle;
	if (startOriented) {
		phi = ranGauss.dev()*TWOPI;              // starting growth direction
		if (random.doub()<0.5) {
			if (phi < PI) phi += PI;
			else phi -= PI;
		}
	}
	else {
		phi = random.doub()*TWOPI;              // starting growth direction
	}
	x1 = random.doub()*winSize;
	y1 = random.doub()*winSize;
	x2 = x1 + dr*cos(phi);
	y2 = y1 + dr*sin(phi);
	unsigned int rodsI = 0;
	rodSegmentCreate(rodsI, x1, x1, y1, y1, phi, 1, 0);

	double Dx, Dy, d;
	double shLen;
	double newrodsRand;
	double bNx1, bNy1;

	unsigned int i,j,k,l,vv;
	map<unsigned int, vector<vector<double> > >::iterator jit, lit;
	map<unsigned int, unsigned int>::iterator nit;
	vector<unsigned int>::iterator iit, mit;
	vector<vector<double> >::iterator vvit;
	vector<vector<double> > vvec;

	char fname[23];
	int tmp = 0;
	vector<double> tmp2, vec;
	// Initialize coarse-grained grid for calculating order parameter
	map<vector<int>, vector<double> > gridTh;
	map<vector<int>, vector<double> >::iterator git;

	cout << "Microtubule simulation: Starting V6.2..." << endl;

	// main loop
	// each k is a different time point
	for (k=0; k<steps; k++) {

		// Display progress bar
		if (k % dsstep == 0) cout << ((double)k/(double)steps)*100 << "% completed." << endl;

		// Output an eps graphic file every 'psstep'
		if (k % psstep == 0) {
			tmp = k / psstep;
			sprintf(fname, "gfx/snapshot_%0.5d.eps", tmp); // k == step
			snapshotPs(rods,fname,k);
		}

		newrodsRand = random.doub();
		if (newrodsRand < rc*dt) {
			if (k<=kSwitch && startOriented) {
				phi = ranGauss.dev()*TWOPI;
				if (random.doub()<0.5) {
					if (phi < PI) phi += PI;
					else phi -= PI;
				}
			}
			else {
				phi = random.doub()*TWOPI;
			}
			x1  = random.doub()*winSize;
			y1  = random.doub()*winSize;
			rodsI++;
			rodSegmentCreate(rodsI, x1, x1, y1, y1, phi, 1, 0);
		}

		// for each rod, check whether it will hit another rod as it grows
		for (jit=rods.begin(); jit!=rods.end(); ++jit) {

			j = jit->first;
			// here check only last member of the rods[j] vector, because this will be the final
			// segment of the rod that (potentially) crossed the bounding box.
			unsigned int m = rods[j].size()-1; // replace with rods[j].back; final segment

			if (rods[j][m][5] == -1) continue; // rod was deleted
			Dx   = dr*cos(rodAngle(j,m));
			Dy   = dr*sin(rodAngle(j,m));
			// periodic boundary conditions are automatically accounted for in rodGrow

			if (rods[j][m][5] == 1) { // only check rods that are growing

				// data for the rod that is the closest one.
				int lMin = -1;
				int mMin = -1;
				double dMin = -1;

				if (DEBUG) cout << "step: " << k << ", rod j: " << j << " (" << rods[j][m][0] << ","
					<< rods[j][m][2] << "),(" << rods[j][m][1] << "," << rods[j][m][3] << ")" << endl;

				// check the distance from rod j to other rods (l)
				for (lit=rods.begin(); lit!=rods.end(); ++lit) {

					l = lit->first;
					if (l==j) continue; // don't check vs itself
					vv=0;

					// check distance to each rod segment
					for (vvit=rods[l].begin(); vvit!=rods[l].end(); ++vvit) {	// go through rods vector, with each element specifying rod segment
						vec = *vvit;
						if (DEBUG) cout << " .. checking " << j << " vs " << l << "(" << vv << "): [" << vec[0] << ", " << vec[2] << "]-[" <<
							vec[1] << ", " << vec[3] << "]@ " << vec[4] << endl;
						d = rodDistance(rods[j][m], vec);

						if (d == -1.0) {
							vv++;
							continue; // they won't intersect
						}
						else { // if it will intersect in some time
							// d != -1 here
							if (lMin == -1) { // ... and its first one
								lMin = l;
								mMin = vv;
								dMin = d;
							}
							else { // if there are other possibly intersecting
								if (dMin>d) {
									lMin = l;
									mMin = vv;
									dMin = d;
								}
							}
						}
						vv++;
					}
				} // end for 'l' loop

				// after previous loop i should have closest rod data in *Mins
				// if there is some rods in front that can block growing rod.
				// otherwise i'll just have -1 everywhere.

				if (dMin < dr && dMin != -1) { // === this is where rod is blocked ===

					// first check from which side is the rod (j?) hitting
					int leftright = 0;
					// to do: correct L-R check; use leftright +-1
					rndLR = random.doub();
					cout << "random number is: " << rndLR << endl;
					if (rods[j][m][4]<=PIHALF && rods[j][m][4]>3*PI/4) { // hitting from left
						if (rndLR <= lBundlingProb) { // make it bundle
							if (bundling) { //rodBundle(j, m, dMin, lMin, mMin, pfix);
								// rod growth function when bundling is enabled
								// double bundlingNewAngle, bNx1, bNy1;
								rodGrow(j, m, dMin);
								rods[j][m][1] -= pfix*cos(rods[j][m][4]); // x_final ***
								rods[j][m][3] -= pfix*sin(rods[j][m][4]); // y_final ***
								if (rodHitAngle(rods[j][m], rods[lMin][mMin]) <= pfix) { // used to be bundlingDistance instead of pfix
									// start growing from rods[j][m][1] (x2) rods[j][m][3] (y2) at direction rods[lMin][mMin][4]
									//
									//if ( fabs(rods[lMin][mMin][4]-rods[j][m][4]) > PIHALF ) {
									double phiDif = max(rods[lMin][mMin][4],rods[j][m][4]) - min(rods[lMin][mMin][4],rods[j][m][4]);
									if (phiDif > PIHALF) {
										// change angle, rod is hitting at a narrow angle but towards the tail of the 'l' rod
										bundlingNewAngle = rods[lMin][mMin][4]>PI ? rods[lMin][mMin][4]-PI : PI+rods[lMin][mMin][4];
									}
									else {
										bundlingNewAngle = rods[lMin][mMin][4];
									}
									bNx1 = rods[j][m][1];
									bNy1 = rods[j][m][3];
									rodSegmentCreate(j, bNx1, bNx1, bNy1, bNy1, bundlingNewAngle, 1, 0);
								}
								else {
									for (unsigned int kk=0; kk<rods[j].size(); ++kk) {
										rods[j][kk][5] = 0;
									}
									inters[lMin].push_back(j);
									rodsB[j] = lMin; // j was blocked by l
								}
							}
							else { // if bundling is disabled, just stop the rod
								// set status of all rod segments to 0
								for (unsigned int kk=0; kk<rods[j].size(); ++kk) {
									rods[j][kk][5] = 0;
								}

								inters[lMin].push_back(j);
								rodsB[j] = lMin; // j was blocked by l
								rodGrow(j, m, dMin); // grow only up to rod that will block it.
								// set the coordinates manually in case it just grew 'too close'
								// make it grow up to 0.01 of the
								rods[j][m][1] -= pfix*cos(rods[j][m][4]); // x_final ***
								rods[j][m][3] -= pfix*sin(rods[j][m][4]); // y_final ***
							}
						}
						else { // make it go through
							// probably I don't need anything here, just leave rods[j][m][5] bit to zero

						}
					} // end if hitting from the left
					else { // hitting from the right
						if (rndLR <= rBundlingProb) {
							if (bundling) { //rodBundle(j, m, dMin, lMin, mMin, pfix);
								// rod growth function when bundling is enabled
								//double bundlingNewAngle, bNx1, bNy1, inters;
								rodGrow(j, m, dMin);
								rods[j][m][1] -= pfix*cos(rods[j][m][4]); // x_final ***
								rods[j][m][3] -= pfix*sin(rods[j][m][4]); // y_final ***
								if (rodHitAngle(rods[j][m], rods[lMin][mMin]) <= pfix) { // used to be bundlingDistance instead of pfix
									// start growing from rods[j][m][1] (x2) rods[j][m][3] (y2) at direction rods[lMin][mMin][4]
									//
									//if ( fabs(rods[lMin][mMin][4]-rods[j][m][4]) > PIHALF ) {
									double phiDif = max(rods[lMin][mMin][4],rods[j][m][4]) - min(rods[lMin][mMin][4],rods[j][m][4]);
									if (phiDif > PIHALF) {
										// change angle, rod is hitting at a narrow angle but towards the tail of the 'l' rod
										bundlingNewAngle = rods[lMin][mMin][4]>PI ? rods[lMin][mMin][4]-PI : PI+rods[lMin][mMin][4];
									}
									else {
										bundlingNewAngle = rods[lMin][mMin][4];
									}
									bNx1 = rods[j][m][1];
									bNy1 = rods[j][m][3];
									rodSegmentCreate(j, bNx1, bNx1, bNy1, bNy1, bundlingNewAngle, 1, 0);
								}
								else {
									for (unsigned int kk=0; kk<rods[j].size(); ++kk) {
										rods[j][kk][5] = 0;
									}
									inters[lMin].push_back(j);
									rodsB[j] = lMin; // j was blocked by l
								}
							}
							else {
								// set status of all rod segments to 0
								for (unsigned int kk=0; kk<rods[j].size(); ++kk) {
									rods[j][kk][5] = 0;
								}

								inters[lMin].push_back(j);
								rodsB[j] = lMin; // j was blocked by l
								rodGrow(j, m, dMin); // grow only up to rod that will block it.
								// set the coordinates manually in case it just grew 'too close'
								// make it grow up to 0.01 of the
								rods[j][m][1] -= pfix*cos(rods[j][m][4]); // x_final ***
								rods[j][m][3] -= pfix*sin(rods[j][m][4]); // y_final ***
							}
						}

					} // end if hitting from the right
					// continue;
				} // if d<dr

			} // end if growth == 1

			if (rods[j][m][5] == 1) { // if growth is still 1, then grow; if the rod bumped into another, its 0 now.
				// growth with length d
				rodGrow(j, m, dr);
	            rods[j][m][6]++; // update rod's internal "clock"
			}

			// ----------------------------------------------------------------------
			//                            Shrink rods here
			// ----------------------------------------------------------------------
			//  * every time we shrink a rod, check whether it was blocking any other rod
			//  * if it was, check if it is still blocking all of them
			//  *  if any of them is unblocked change its state [5] = 1;
			shLen = vmin*dt;
			rodShrink(j, shLen);

			if (inters[j].size() != 0) { 	// if this rod blocks/was blocking any other rods
                                          	// they can now resume growth.
				for (iit=inters[j].begin(); iit!=inters[j].end(); ++iit) { // go through all others that its blocking
					// check if it is still blocking each one
					i = *iit;
					// for each other rod that was being blocked, calculate the distance from l tip to j rod
					// check all segments and pick the minimum distance.
					// if this min is less than one growth step (make it 'vplus*dt') its still being blocked
					// so don't do anything with it.

					// since only [0] segment of the j rod gets shrunk, check only the distance between growing segments [last]
					// of all the rods it was blocking and the [0] segment of j
					double ddMin = -1;
					vv = 0;
					for (vvit=rods[j].begin(); vvit!=rods[j].end(); ++vvit) {
						vec = *vvit;
						d = rodDistance(rods[i][rods[i].size()-1], vec);
						if (d!=-1) {
							if (ddMin!=-1) {
								if (d<ddMin) {
									ddMin = d;
								}
							}
							else { // ddMin=-1; d!=-1
								ddMin = d;
							}
						}
						vv++;
					}
					// not good enough to check only [0]:
					// i might have hit j like this
					//      /        /
					//     / i[last]/
					//  j0/    ----/
					//   /        /
					//  /        /j1

					if (ddMin == -1) {
						// rod rods[i][rods[i].size()-1] just became free;
						// change this and make sure all segments[5] are updated.
						rods[i][rods[i].size()-1][5] = 1;
						inters[j].erase(iit); // this invalidates for inters iterator, so we have to break
						break;
					}

				}

			} // end if this rod blocks other rods

		} // end for j loop

		// -----------------------------------------------------------
		//            Write order parameter moments to file
		//------------------------------------------------------------
		if (k % pMstep == 0) {
			tmp = k/pMstep;
			map<unsigned int, double> kMoments = orderMoments(phi_m);
			phiMomentsFile << k;
			for (unsigned int ii=0; ii<=phi_m; ++ii) {
				phiMomentsFile << "\t" << kMoments[ii];
			}
			phiMomentsFile << endl;

		}
	} // end main k loop

	return 0;
}


double rodLength(unsigned int rIdx) {
	// Calculate the rod length of rod with id 'rIdx' (sum over all segments)
	double totLen = 0.0;
	vector<vector<double> >::iterator vvit;
	vector<double> vec;
	for (vvit=rods[rIdx].begin(); vvit!=rods[rIdx].end(); ++vvit) {
		vec = *vvit;
		totLen += sqrt( pow(vec[1]-vec[0],2) + pow(vec[3]-vec[2],2) );
	}
	return totLen;
}

double rodSegLength(unsigned int rIdx, unsigned int sIdx) {
	// Calculate the length of the segment 'sIdx' for the rod 'rIdx'
	return sqrt( pow(rods[rIdx][sIdx][1]-rods[rIdx][sIdx][0],2) + pow(rods[rIdx][sIdx][3]-rods[rIdx][sIdx][2],2) );
}

void rodSegmentCreate(unsigned int rIdx, double x0, double x1, double y0, double y1, double phi, double status, double clock) {
	// Create a new rod segment
	vector<double> rodSegment;
	rodSegment.push_back(x0);
	rodSegment.push_back(x1);
	rodSegment.push_back(y0);
	rodSegment.push_back(y1);
	rodSegment.push_back(phi);
	rodSegment.push_back(status);
	rodSegment.push_back(clock);
	rods[rIdx].push_back(rodSegment);
	rodSegment.clear();

}

void rodShrink(unsigned int rIdx, double len) {
  //
	// Shrink the rod with the id 'rIdx', by length 'len'
  //
	// sIdx=0, because 0th rod is always being shortened (this is the "tail")
	double Dx = 0;
	double Dy = 0;
	double angle = 0;
	double sLen = rodSegLength(rIdx,0);
	double dLen = 0; // total length deleted (need to deduce it when i'm shrinking final rod segment)

	// try to shrink from (x0,y0) at the beginning of the first rod segment
	// check whether the first segment is longer then 'len'
	// if it is, then remove it. then check again new 1st segment, etc.
	// this will be while() loop
	// be sure to check if there are any rod segments left after the loop is over.
	while (sLen < len) {
		if (rods[rIdx].size() == 1) {
			rods[rIdx][0][5] = -1.0;
			break;
		}
		rods[rIdx].erase(rods[rIdx].begin());
		dLen += sLen;
		sLen  = rodSegLength(rIdx,0);
	}

	// update angle
	angle = rodAngle(rIdx,0);

	if (rods[rIdx][0][5] != -1) {
		Dx	= (len-dLen)*cos(angle);
		Dy	= (len-dLen)*sin(angle);
		rods[rIdx][0][0] += Dx;
		rods[rIdx][0][2] += Dy;
	}
}

void rodGrow(unsigned int rIdx, unsigned int m, double len) {
	//
	//  Rod grow function. Try to grow a rod 'rIdx' by 'len' if there's nothing
	//  obstructing it's path. If there's another rod in the path, either collide,
	//  pass-through, or bundle depending on simulation parameters.
	//
	double phi  = 0;

	if (rods[rIdx][m][0] == rods[rIdx][m][1]) phi = rods[rIdx][m][4];
	else phi  = rodAngle(rIdx,m);

	double Dx	= len*cos(phi);
	double Dy	= len*sin(phi);
	double newX = rods[rIdx][m][1]+Dx;
	double newY = rods[rIdx][m][3]+Dy;

	if (newX > winSize) {
		// x is going out of the bounding box at the right side
		// - make new rod segment and close the current one with the
		//   boundaries of the bbox
		rods[rIdx][m][1] = winSize;
		rodSegmentCreate(rIdx, newX - winSize, newX - winSize, newY, newY, rods[rIdx][m][4],
			rods[rIdx][m][5], rods[rIdx][m][6]);
		//rods[rIdx][1] = newX - winSize;
	}
	else if (newX < 0) {
		// x is going out of the bounding box at the left side
		rods[rIdx][m][1] = 0;
		rodSegmentCreate(rIdx, newX + winSize, newX + winSize, newY, newY, rods[rIdx][m][4],
			rods[rIdx][m][5], rods[rIdx][m][6]);
	}
	else {
		// x stays inside
		rods[rIdx][m][1] = newX;
	}

	if (newY > winSize) {
		rods[rIdx][m][3] = winSize;
		rodSegmentCreate(rIdx, newX, newX, newY - winSize, newY - winSize, rods[rIdx][m][4],
			rods[rIdx][m][5], rods[rIdx][m][6]);
	}
	else if (newY < 0) {
		rods[rIdx][m][3] = 0;
		rodSegmentCreate(rIdx, newX, newX, newY + winSize, newY + winSize, rods[rIdx][m][4],
			rods[rIdx][m][5], rods[rIdx][m][6]);
	}
	else {
		rods[rIdx][m][3] = newY;
	}
		// 0 ... x1,  1 ... x2,  2 ... y1,  3 ... y2,  4 ... phi
}

double rodAngle(unsigned int rIdx, unsigned int sIdx) {
	// Extract the smaller angle between two rods.
	double phi;
	if (rods[rIdx][sIdx][1]!=rods[rIdx][sIdx][0]) phi = atan2(rods[rIdx][sIdx][3]-rods[rIdx][sIdx][2], rods[rIdx][sIdx][1]-rods[rIdx][sIdx][0]);
	else phi = rods[rIdx][sIdx][4];
	if (phi < 0) phi = 2*PI + phi;
	return phi;
}

double rodAngle2(vector<double> rod) {
	// Calculate the actual rod angle using atan2
	double phi;
	if (rod[0]==rod[1]) phi=rod[4];
	else phi=atan2(rod[3]-rod[2],rod[1]-rod[0]);
	return phi;
}


void snapshotPs(map <unsigned int, vector<vector<double> > > rods, const char * fileName, unsigned int k) {
	//
	//  Write the instantaneous state of the rods (microtubules) to the EPS file.
	//
	ofstream outFile;
	outFile.open(fileName);
	// output .ps header
	outFile << "%!PS-Adobe-2.0 EPSF-2.0" << endl;
	outFile << "%%BoundingBox: 0 0 " << winSize << " " << winSize << endl;
	outFile << "%%HiResBoundingBox: 0 0 " << winSize << " " << winSize << endl;
	outFile << "%%CropBox: 0 0 " << winSize << " " << winSize << endl;
	outFile << "%%EndComments" << endl;
	outFile << "%%BeginProlog" << endl;
	outFile << "newpath 0 0 moveto 0 " << winSize << " lineto " << winSize << " " << winSize << " lineto " << winSize << " 0 lineto closepath clip" << endl;
	outFile << "save" << endl;
	outFile << "countdictstack" << endl;
	outFile << "mark" << endl;
	outFile << "/showpage {} def" << endl;
	outFile << "/setpagedevice {pop} def" << endl;
	outFile << "/winSize " << winSize << " def" << endl;
	outFile << "%%EndProlog" << endl;
	outFile << "%%Page: (1) 1" << endl;
	outFile << "% draw a bounding box of winSize" << endl;
	// white background
	outFile << "/Times-Roman findfont" << endl;
	outFile << "14 scalefont" << endl;
	outFile << "setfont" << endl;
	outFile << "newpath" << endl;
	outFile << "0 0 moveto" << endl;
	outFile << "winSize 0 lineto" << endl;
	outFile << "winSize winSize lineto" << endl;
	outFile << "0 winSize lineto" << endl;
	outFile << "closepath" << endl;
	outFile << "gsave" << endl;
	outFile << "1.0 setgray" << endl;
	outFile << "fill" << endl;
	outFile << "grestore" << endl;
	outFile << "0.5 setlinewidth" << endl;
	outFile << "newpath" << endl;
	outFile << "5 5 moveto" << endl;
	outFile << "(" << k << ") show" << endl;
	outFile << "stroke" << endl;
	// draw all the lines
  for (map<unsigned int, vector<vector<double> > >::iterator jit=rods.begin(); jit!=rods.end(); ++jit) {
		int i = jit->first;
		for (vector<vector<double > >::iterator vvit=rods[i].begin(); vvit != rods[i].end(); ++vvit) {
			vector<double> rod = *vvit;
			if (rod[5] == -1) continue;
			outFile << rod[0] << " " << rod[2] << " moveto " << rod[1]
					<< " " << rod[3] << " lineto" << " 0.1 setlinewidth stroke" << endl;
		}
  }
	outFile << "showpage" << endl;
	outFile << "%%Trailer" << endl;
	outFile << "cleartomark" << endl;
	outFile << "countdictstack" << endl;
	outFile << "exch sub { end } repeat" << endl;
	outFile << "restore" << endl;
	outFile << "%%EOF" << endl;

}


double rodDistance(vector<double> rodsJ, vector<double> rodsL) {
	//
	// Calcualte distance between two rods if they will intersect as rodJ is growing from its end.
	// If they won't intersect, this function returns -1.0
	//
	// rodsL is a vector containing possible multiple elements (rod segments)
	// so we need to check whether it can hit either of those.
	double h1x,h1y,h2x,h2y,phi1,phi2,phii,k1,l1,k2,l2,Dk,intX,intY,d;
	// coords of the 1st edge
	h1x  = rodsL[0] - rodsJ[1];
	h1y  = rodsL[2] - rodsJ[3];
	// coords of the 2nd edge
	h2x  = rodsL[1] - rodsJ[1];
	h2y  = rodsL[3] - rodsJ[3];

	// angles of the two helper vectors h1, h2
	// in [0, 2PI].
	phi1 = atan2(h1y,h1x);
	phi2 = atan2(h2y,h2x);

	if (rodsJ[0] == rodsJ[1]) phii = rodsJ[4]; // if the length is 0 (x1==x2) just take initial angle
	else phii = atan2(rodsJ[3]-rodsJ[2], rodsJ[1]-rodsJ[0]); // atan( (y2-y1)/(x2-x1) )

	// make both phi12 defined [0,2PI]
	if (phii<0) phii = 2*PI + phii;
	if (phi1<0.0) phi1 = 2*PI + phi1;
	if (phi2<0.0) phi2 = 2*PI + phi2;

	// check whether rod 'l' is on a way of growing rod j
		// This is done by checking whether the angle of growth
		// direction is between two angles spanning the rod j
		// from the point of the growth of rod j (O).
		// -> is the direction of growth of j (its not necessarily
		// perpendicular to rod l).
		//            /|
		//         h1/ | r
		//  rod j   /  | o
		// --------O-> | d
		//          \  |
		//         h2\ | l
		//            \|
		// If the rod l isn't on collision with j, we skip this
		// loop iteration (continue).

	// check for collision (don't caluclate distance if rods aren't colliding)
	if (fabs(fabs(phi2-phi1)-PI) < precision) {
		// this means phi1 and phi2 span angle of PI up to high accuracy (precision)
		// and this means that the incoming rod L hit J
		// just return zero distance
		return 0.0;
	}
	//if (DEBUG) cout << " .. .. .. phi1: " << phi1 << ", phi2: " << phi2 << ", phii: "
	//	<< phii << ":: dif: " << phi2-phi1 << endl;

	if (phi2 > phi1) {
		//cout << phi2-phi1 << " vs " << PI << endl;
		if ((phi2-phi1)<=PI) {
			//if (DEBUG) cout << " .. ... .. case 1" << endl;
			if (phii<phi1 || phii>phi2) return -1.0;
		}
		else {
			//if (DEBUG) cout << " .. ... .. case 2" << endl;
			if (phii<phi2 && phii>phi1) return -1.0;
		}
	}
	else { // phi1>phi2
		if ((phi1-phi2)<=PI) {
			//if (DEBUG) cout << " .. ... .. case 3" << endl;
			if (phii<=phi2 || phii>=phi1) return -1.0;
		}
		else {
			//if (!(phi1<=phii || phii<=phi2)) return -1.0;
			//if (DEBUG) cout << " .. ... .. case 4" << endl;
			if (phii>phi2 && phii<=phi1) return -1.0;
		}
	}
	//}

	//else {
		// h1 and h2 are at ~ 180 degrees, which
		// means the rodJ that is incoming almost hit the rodL
		// its unlikely in that case that the rodL is "behind"
	//}

		// We can only get at this point, if it's shown that two rods
		// will eventually collide.

	// Find the intersection point, then calc distance between
	// it and growing point of j.

		// We will find intersection between (inf) lines of growing
		// and the one rod j sits on.
		// eqs of lines; f_{12}(x) = k_{12}x + l_{12}

	// make function that returns true if two rods j and l will intersect
	k1   = tan(phii);
	l1   = rodsJ[3] - k1*rodsJ[1];
	k2   = (rodsL[3] - rodsL[2]) / (rodsL[1] - rodsL[0]);
	l2   = rodsL[3] - k2*rodsL[1];
	// intersection
	Dk   = k1-k2;
	intX = (l2-l1)/Dk;
	intY = (k1*l2 - l1*k2)/Dk;
	// distance
	d    = sqrt( pow(intX-rodsJ[1],2) + pow(intY-rodsJ[3],2) );
	return d;
}



double rodHitAngle(vector<double> rodsJ, vector<double> rodsL) {
	// calculates hit angle; returns angle between 0 and PI/2
	double th1 = fabs(rodAngle2(rodsJ)-rodAngle2(rodsL));
	double th2 = fabs(PI-th1);
	return min(th1, th2);
}


void gridInit(unsigned int fftPowerOfTwo, double len) {
	//
	// Initialize the coarse-grained grid for calculating the order parameter.
	//
	map <vector<int>, vector<double> >::iterator git; // grid iterator
	unsigned int n = (unsigned int)pow(2.0, (double)fftPowerOfTwo);
	double		dn = n/len;
	int			in = 0;
	vector<int> coords;
	vector<double> thetas;
	thetas.clear();

	for (unsigned int i=0; i<n; ++i) {
		coords.clear();
		in = i*n;
		coords.push_back(in);
		coords.push_back(in+n);
		grid[coords] = thetas;
	}

}

void gridSnapshot() {
	//
	//  Calculate the current state of the coarse-grained grid.
	//
	map <vector<int>, vector<double> >::iterator git;
	map<unsigned int, vector<vector<double> > >::iterator jit;
	vector<vector<double> >::iterator vvit;
	vector<double> vec;
	unsigned int l;
	//double x1,y1,x2,y2;
	int nx1,ny1,nx2,ny2;
	// go through all of the rods
	for (jit=rods.begin(); jit!=rods.end(); ++jit) {
		// go through each rod segment and treat it like a separate rod
		l = jit->first;
		for (vvit=rods[l].begin(); vvit!=rods[l].end(); ++vvit) {
			vec = *vvit;
			nx1 = (int)floor(vec[0]);
			nx2 = (int)floor(vec[1]);
			ny1 = (int)floor(vec[2]);
			ny2 = (int)floor(vec[3]);
		}
	}

}

double complexAbs(map<unsigned int, double> z) {
	//
	// Return the absolute value of a complex number z
	//
	double a = z[0];
	double b = z[1];
	double aa = fabs(a);
	double ab = fabs(b);
	if (aa < ab) return ab*sqrt(1+pow(a/b,2));
	else return aa*sqrt(1+pow(b/a,2));
}

map<unsigned int, double> orderMoments(unsigned int m) {
	//
	// Calculate the order parameter moments up to m-th.
	//
	map<unsigned int, vector<vector<double> > >::iterator jit;
	vector<vector<double> >::iterator vvit;
	map<unsigned int, double> ls; // total length of rod k
	map<unsigned int, double> thetas; // rod k orientation
	//map<unsigned int, vector<double> > phiM;
	map<unsigned int, map<unsigned int, double> > sumPhi;
	map<unsigned int, double> res;
	vector<double> vec;
	double lsk, thk; // rod segment length

	for (unsigned int i=0; i<=m; ++i) {
		sumPhi[i][0] = 0.0;
		sumPhi[i][1] = 0.0;
	}
	unsigned int j,k;
	for (jit=rods.begin(); jit!=rods.end(); ++jit) { // go through all rods
		k = jit->first;
		if (rods[k][0][5] == -1) continue; // deleted rod; skip it
		for (vvit=rods[k].begin(); vvit!=rods[k].end(); ++vvit) { // go through all rod segments
			vec = *vvit; // each vec is one rod segment
			lsk = sqrt(pow(vec[1]-vec[0],2)+pow(vec[3]-vec[2],2));
			thk = rodAngle2(vec);
			for (unsigned int i=0; i<=m; ++i) {
				if (i!=0) {
					sumPhi[i][0]  += lsk*cos(i*thk);
					sumPhi[i][1]  += lsk*sin(i*thk);
				}
				else {
					sumPhi[i][0]  += lsk;
				}
			}
		}
	}
	for (j=0; j<=m; ++j) {
		res[j] = complexAbs(sumPhi[j]);
	}
	return res;

}
