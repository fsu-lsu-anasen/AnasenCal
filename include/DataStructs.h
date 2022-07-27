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
	double frontup_energy_adc = -1;
	double frontdown_energy_adc = -1;
	double back_energy = -1;
	int frontup_gchan = -1;
	int frontdown_gchan = -1;
	int back_gchan = -1;
	int detector_index = -1;
};

struct CalibratedQQQHit
{
	double ring_energy = -1;
	double wedge_energy = -1;
	int ring_gchan = -1;
	int wedge_gchan = -1;
	int detector_index = -1;
};

struct CalibratedBarcHit
{
	double front_energy = -1.0;
	int front_gchan = -1;
	int detector_index = -1;
};

struct CalibratedEvent
{
	std::vector<CalibratedSX3Hit> barrel;
	std::vector<CalibratedQQQHit> fqqq;
	std::vector<CalibratedBarcHit> barcUp;
	std::vector<CalibratedBarcHit> barcDown;
};

#endif
