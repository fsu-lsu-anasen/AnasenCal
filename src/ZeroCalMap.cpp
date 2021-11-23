/*
	ZeroCalMap
	Wrapper around std::unordered_map which gives acces to an offset parameter using a global channel key.
	See ChannelMap documentation for more information on global channels.

	Key method is FindOffset. This returns an iterator, which can be checked against the end of the map for validity.
	If the iterator is equal to the end, the key does not link to a member of the map.

	Written by Gordon McCann Nov 2021
*/
#include "ZeroCalMap.h"
#include <fstream>

ZeroCalMap::ZeroCalMap(const std::string& filename) :
	valid_flag(false)
{
	FillMap(filename);
}

ZeroCalMap::~ZeroCalMap() {}

void ZeroCalMap::FillMap(const std::string& filename)
{
	std::ifstream input(filename);
	if(!input.is_open())
	{
		valid_flag = false;
		return;
	}

	int gchan;
	double offset;

	while(input>>gchan)
	{
		input>>offset;
		zmap[gchan] = offset;
	}

	input.close();
	valid_flag = true;
}