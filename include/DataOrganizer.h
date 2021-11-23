/*
	DataOrganizer
	Class to convert data from the raw motherboard-channel-based structure of evt2root to
	a more practical detector-based structure, AnasenEvent. This requires that a channel map
	exists, otherwise the conversion is not possible.

	Written by Gordon McCann Nov. 2021
*/
#ifndef DATAORGANIZER_H
#define DATAORGANIZER_H

#include <string>
#include "DataStructs.h"
#include "ChannelMap.h"
#include <TRandom3.h>

class DataOrganizer
{
public:
	DataOrganizer(const std::string& channelfile);
	~DataOrganizer();

	void Run(const std::string& inputname, const std::string& outputname);
private:
	void FillEvent(AnasenEvent& event, int gchan, int energy);
	//When switching from integers to floating point, need to smear within the bin.
	inline double ConvertInt2Double(int value) { return value + generator->Uniform(0.0, 1.0); }


	ChannelMap cmap;
	TRandom3* generator;
};


#endif