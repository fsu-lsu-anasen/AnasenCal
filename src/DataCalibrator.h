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
#include "ParameterMap.h"
#include "Geometry/AnasenArray.h"

class DataCalibrator
{
public:
	DataCalibrator(const std::string& channelfile, const std::string& backmatch, const std::string& updownmatch, 
					const std::string& frontbackmatch, const std::string& energyfile);
	~DataCalibrator();
	void Run(const std::string& inputname, const std::string& outputname);

private:
	ChannelMap m_channelMap;
	ParameterMap m_backGainMap, m_sx3UpDownGainMap, m_frontBackGainMap, m_energyCalMap;
	static constexpr int s_sx3UpDownMatch[8] = {1, 0, 3, 2, 5, 4, 7, 6}; //matches index -> index of up/down pair for SX3 fronts
};

#endif
