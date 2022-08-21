/*
	ParameterMap
	Wrapper around std::unordered_map which gives acces to fitting parameters based on a global channel key.
	See ChannelMap documentation for more information on global channels.

	Key method is FindParameters. This returns an iterator, which can be checked against the end of the map for validity.
	If the iterator is equal to the end, the key does not link to a member of the map.

	Written by Gordon McCann Nov 2021
*/
#ifndef PARAMETERMAP_H
#define PARAMETERMAP_H

#include <string>
#include <unordered_map>
#include "DataStructs.h"

class ParameterMap
{
public:
	typedef std::unordered_map<int, CalParams>::iterator Iter; 
	ParameterMap(const std::string& filename);
	~ParameterMap();
	inline Iter FindParameters(int gchan) { return map.find(gchan); }
	inline Iter End() { return map.end(); }
	inline const bool IsValid() const { return valid_flag; }

private:
	void FillMap(const std::string& filename);
	
	bool valid_flag;
	std::unordered_map<int, CalParams> map;
};

#endif