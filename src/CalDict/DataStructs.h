/*
	Data structures used by the AnasenCal program.
	Many of these are given to ROOT for dictionary generation, therefore
	it may be necessary to port some or all of this file to projects which use
	data made by AnasenCal.

	Written by Gordon McCann Nov 2021
*/
#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <vector>
#include "Math/Point3D.h"


struct DetectorHit
{
	int globalChannel;
	double energy;
	double timestamp;
};

struct FQQQDetector
{
	std::vector<DetectorHit> rings;
	std::vector<DetectorHit> wedges;
};

struct BarrelDetector
{
	std::vector<DetectorHit> frontsUp;
	std::vector<DetectorHit> frontsDown;
	std::vector<DetectorHit> backs;
};

struct BarcDetector
{
	std::vector<DetectorHit> fronts;
	std::vector<DetectorHit> backs;
};

struct CoincEvent
{
	BarrelDetector barrel[12];
	FQQQDetector fqqq[4];
	BarcDetector barcUp[6];
	BarcDetector barcDown[6];
};

struct GraphData
{
	std::vector<double> xvals;
	std::vector<double> yvals;
};

struct CalParams
{
	double slope = 0;
	double intercept = 0;
};

struct CalibratedSX3Hit
{
	double frontUpEnergyAdc = -1;
	double frontDownEnergyAdc = -1;
	double backEnergy = -1;
	int frontUpGlobalChannel = -1;
	int frontUpLocalChannel = -1;
	int frontDownGlobalChannel = -1;
	int frontDownLocalChannel = -1;
	int backGlobalChannel = -1;
	int backLocalChannel = -1;
	int detectorIndex = -1;

	ROOT::Math::XYZPoint coordinates;
};

struct CalibratedQQQHit
{
	double ringEnergy = -1;
	double wedgeEnergy = -1;
	int ringGlobalChannel = -1;
	int ringLocalChannel = -1;
	int wedgeGlobalChannel = -1;
	int wedgeLocalChannel = -1;
	int detectorIndex = -1;
	ROOT::Math::XYZPoint coordinates;
};

struct CalibratedBarcHit
{
	double frontEnergy = -1.0;
	int frontGlobalChannel = -1;
	int frontLocalChannel = -1;
	int detectorIndex = -1;
	ROOT::Math::XYZPoint coordinates;
};

struct CalibratedEvent
{
	std::vector<CalibratedSX3Hit> barrel;
	std::vector<CalibratedQQQHit> fqqq;
	std::vector<CalibratedBarcHit> barcUp;
	std::vector<CalibratedBarcHit> barcDown;
};


bool EnforceDictionaryLinked();
#endif
