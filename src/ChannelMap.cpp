/*
	ChannelMap
	Wrapper around std::unordered_map which gives acces to channel information using a global channel number key.
	Global channel number is given by chipboard_num*32 + channel (with offsets to account for changes in MB).
	This particular version is for the solid-target ANASEN detector array,

	Key method is FindChannel. This returns an iterator, which can be checked against the end of the map for validity.
	If the iterator is equal to the end, the key does not link to a member of the map.
	
	Written by Gordon McCann Nov 2021
*/

#include "ChannelMap.h"
#include <fstream>
#include <iostream>

ChannelMap::ChannelMap(const std::string& filename) :
	valid_flag(false), name(filename)
{
	FillMap(filename);
}

ChannelMap::~ChannelMap() {}

void ChannelMap::FillMap(const std::string& filename)
{
	std::ifstream input(filename);
	if(!input.is_open())
	{
		valid_flag=false;
		std::cerr<<"Bad file at ChannelMap::FillMap!"<<std::endl;
		return;
	}

	ChannelData data;
	int gchan;
	while(input>>gchan)
	{
		input>>data.detectorType>>data.detectorID>>data.detectorComponent>>data.detectorDirection>>data.channel;
		cmap[gchan] = data;
	}

	valid_flag = true;
}

int ChannelMap::ConvertSX3Name2Index(const std::string& detectorType, const std::string& detectorID)
{
	if(detectorType == "BARREL1A" || detectorType == "BARREL2A")
	{
		if(detectorID == "A")
		{
			return 0;
		}
		else if(detectorID == "B")
		{
			return 1;
		}
		else if(detectorID == "C")
		{
			return 2;
		}
		else if(detectorID == "D")
		{
			return 3;
		}
		else if(detectorID == "E")
		{
			return 4;
		}
		else if(detectorID == "F")
		{
			return 5;
		}
	} 
	else if(detectorType == "BARREL1B" || detectorType == "BARREL2B")
	{
		if(detectorID == "A")
		{
			return 6;
		}
		else if(detectorID == "B")
		{
			return 7;
		}
		else if(detectorID == "C")
		{
			return 8;
		}
		else if(detectorID == "D")
		{
			return 9;
		}
		else if(detectorID == "E")
		{
			return 10;
		}
		else if(detectorID == "F")
		{
			return 11;
		}
	}
	
	std::cerr<<"Invalid detectorType "<<detectorType<<" or detectorID "<<detectorID<<" at ChannelMap::ConvertSX3Name2Index! Returning -1."<<std::endl;
	return -1;
}

int ChannelMap::ConvertQQQName2Index(const std::string& detectorID)
{
	if(detectorID == "0")
		return 0;
	else if (detectorID == "1")
		return 1;
	else if (detectorID == "2")
		return 2;
	else if (detectorID == "3")
		return 3;

	std::cerr<<"Invalid detectorID "<<detectorID<<" at ChannelMap::CovnertQQQName2Index! Returning -1."<<std::endl;
	return -1;
}

int ChannelMap::InverseFindChannel(const ChannelData& data)
{
	for(auto& channel : cmap)
	{
		ChannelData& this_data = channel.second;
		if(this_data == data)
			return channel.first;
	}
	return -1;
}
