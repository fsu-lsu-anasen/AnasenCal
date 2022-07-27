/*
	EnergyCalibrator
	Class designed to perform energy calibrations. Designed for use with gain-matched results and source data.
	Main method is Run.

	Written by Gordon McCann Nov 2021
*/

#ifndef ENERGYCALIBRATOR_H
#define ENERGYCALIBRATOR_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <THashTable.h>
#include <TSpectrum.h>
#include "ChannelMap.h"
#include "ParameterMap.h"
#include "DataStructs.h"


class EnergyCalibrator {
  
public:
	EnergyCalibrator(const std::string& channelfile, const std::string& backmatch, const std::string& updownmatch, const std::string& frontbackmatch);
	~EnergyCalibrator();
	void Run(const std::string& inputname, const std::string& plotname, const std::string& outputname);

private:
	CalParams CalibrateEnergy(THashTable* table, const std::string& name, const GraphData& data);
	void FillHistogram(THashTable* table, const std::string& name, const std::string& title, int bins, double minx, double maxx, double value);
	GraphData GetPoints(THashTable* table, int gchan, const std::string& name);

	TSpectrum spec;

	ChannelMap cmap;
	ParameterMap bmap, udmap, fbmap;

	double sigma, threshold;
	const int nchannels = 640;

	std::vector<double> energyValues = {5.155, 5.486, 5.805}; //Will need modified for each experiment.

};

#endif
