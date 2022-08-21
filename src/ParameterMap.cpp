/*
	ParameterMap
	Wrapper around std::unordered_map which gives acces to fitting parameters based on a global channel key.
	See ChannelMap documentation for more information on global channels.

	Key method is FindParameters. This returns an iterator, which can be checked against the end of the map for validity.
	If the iterator is equal to the end, the key does not link to a member of the map.

	Written by Gordon McCann Nov 2021
*/
#include "ParameterMap.h"
#include <fstream>

ParameterMap::ParameterMap(const std::string& filename) :
	m_isValid(false)
{
	FillMap(filename);
}

ParameterMap::~ParameterMap() {}

void ParameterMap::FillMap(const std::string& filename)
{
	std::ifstream input(filename);
	if(!input.is_open())
	{
		m_isValid = false;
		return;
	}

	int gchan;
	CalParams params;
	while(input>>gchan)
	{
		input>>params.intercept>>params.slope;
		m_map[gchan] = params;
	}
	input.close();

	m_isValid = true;
}