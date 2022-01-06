/*
	ZeroCalibrator
	Class for generating zero-offset calibrations for ANASEN silicon detectors. ASICS do not fix the zero point on the ADC scale,
	therefore it is necessary to determine the zero point after the experiment. Pulser data is fed in, and TSpectrum is used to
	identify the pulser peaks. The offset is then determined from a linear fit.

	Written by Gordon McCann Nov 2021
*/
#ifndef ZEROCALIBRATOR_H
#define ZEROCALIBRATOR_H

#include <string>
#include <vector>
#include <THashTable.h>
#include <TSpectrum.h>
#include "ChannelMap.h"
#include "DataStructs.h"

class ZeroCalibrator
{

public:
	ZeroCalibrator(const std::string& channelfile);
	~ZeroCalibrator();
	void Run(const std::string& inputname, const std::string& plotname, const std::string& outputname);
	void RecoverOffsets(const std::string& inputname, const std::string& plotname, const std::string& outputname);

private:
	void FillHistogram(THashTable* table, const std::string& name, const std::string& title, int bins, double minx, double maxx, double value);
	void FillHistogram(THashTable* table, const std::string& name, const std::string& title, int binsx, double minx, double maxx, double valuex,
																							int binsy, double miny, double maxy, double valuey);
	GraphData GetPoints(THashTable* table, int gchan, const std::string& histoname);
	GraphData GetPointsAlphas(THashTable* table, int gchan, const std::string& histoname);
	double MakeGraph(THashTable* table, const std::string& name, const GraphData& data);

	TSpectrum spec;

	double sigma, threshold;

	ChannelMap cmap;

	/****Experiment parameters****/
	std::vector<double> frontPulseValues = {1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 10.0}; //Values from experiment, should be adjusted each data set
	std::vector<double> backPulseValues = {0.25, 0.5, 0.75, 1.0, 1.25};
	std::vector<double> energyValues = {5.155, 5.486, 5.805}; //Will need modified for each experiment.

};

#endif