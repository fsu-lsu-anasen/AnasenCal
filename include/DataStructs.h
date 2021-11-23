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

struct SiliconHit
{
	int global_chan = -1;
	int local_chan = -1;
	double energy = -1;
};

struct SX3Data
{
	std::vector<SiliconHit> fronts_up;
	std::vector<SiliconHit> fronts_down;
	std::vector<SiliconHit> backs;
};

struct QQQData
{
	std::vector<SiliconHit> rings;
	std::vector<SiliconHit> wedges;
};

struct AnasenEvent
{
	SX3Data barrel1[12];
	SX3Data barrel2[12];
	QQQData fqqq[4];
	QQQData bqqq[4];
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

struct CalibratedEvent
{
	std::vector<CalibratedSX3Hit> barrel1;
	std::vector<CalibratedSX3Hit> barrel2;
	std::vector<CalibratedQQQHit> fqqq;
	std::vector<CalibratedQQQHit> bqqq;
};

#endif