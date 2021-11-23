/*
	DataCalibrator
	Class which applies all of the calibration results to a data set, performs front-back hit assignment (kind of)
	and saves to a condensed data format (CalibratedEvent) for further analysis. This essentially a single function
	with some class-global parameter maps.

	Written by Gordon McCann Nov 2021
*/
#ifndef DATACALIBRATOR_H
#define DATACALIBRATOR_H

#include <string>
#include "ChannelMap.h"
#include "ZeroCalMap.h"
#include "ParameterMap.h"

class DataCalibrator
{
public:
	DataCalibrator(const std::string& channelfile, const std::string& zerofile, const std::string& backmatch, const std::string& updownmatch, 
					const std::string& frontbackmatch, const std::string& energyfile);
	~DataCalibrator();
	void Run(const std::string& inputname, const std::string& outputname);

private:
	ChannelMap channel_map;
	ZeroCalMap zero_map;
	ParameterMap back_map, updown_map, frontback_map, energy_map;
	const int updown_list[8] = {1, 0, 3, 2, 5, 4, 7, 6}; //matches index -> index of up/down pair for SX3 fronts
};

#endif