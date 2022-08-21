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
	m_isValid(false), m_name(filename)
{
	FillMap(filename);
}

ChannelMap::~ChannelMap() {}

void ChannelMap::FillMap(const std::string& filename)
{
	std::ifstream input(filename);
	if(!input.is_open())
	{
		m_isValid=false;
		std::cerr<<"Bad file at ChannelMap::FillMap!"<<std::endl;
		return;
	}

	ChannelData data;
	int gchan;
	while(input>>gchan)
	{
		input>>data.detectorType>>data.detectorID>>data.detectorComponent>>data.localChannel;
		m_cmap[gchan] = data;
	}

	m_isValid = true;
}

int ChannelMap::InverseFindChannel(const ChannelData& data)
{
	for(auto& channel : m_cmap)
	{
		ChannelData& this_data = channel.second;
		if(this_data == data)
			return channel.first;
	}
	return -1;
}
