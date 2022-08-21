/*
	MapChecker
	Simple class which reads in various calibration files and checks the statistics of the different
	calibration stages. Mostly a debug concept.

	Written by Gordon McCann Nov 2021
*/
#ifndef MAPCHECKER_H
#define MAPCHECKER_H

#include <string>
#include "ChannelMap.h"

class MapChecker
{
public:
	MapChecker(const std::string& channelfile);
	~MapChecker();
	void CheckBackGainMatch(const std::string& filename);
	void CheckUpDownGainMatch(const std::string& filename);
	void CheckFrontBackGainMatch(const std::string& filename);

private:
	ChannelMap m_channelMap;
	static constexpr int s_nchannels = 640;
};

#endif
