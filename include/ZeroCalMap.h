/*
	ZeroCalMap
	Wrapper around std::unordered_map which gives acces to an offset parameter using a global channel key.
	See ChannelMap documentation for more information on global channels.

	Key method is FindOffset. This returns an iterator, which can be checked against the end of the map for validity.
	If the iterator is equal to the end, the key does not link to a member of the map.

	Written by Gordon McCann Nov 2021
*/
#ifndef ZEROCALMAP_H
#define ZEROCALMAP_H

#include <unordered_map>
#include <string>

class ZeroCalMap
{

public:
	typedef std::unordered_map<int, double>::iterator Iter;
	
	ZeroCalMap(const std::string& filename);
	~ZeroCalMap();
	inline Iter FindOffset(int gchan) { return zmap.find(gchan); }
	inline Iter End() { return zmap.end(); }
	inline const bool IsValid() const { return valid_flag; }

private:
	void FillMap(const std::string& filename);

	std::unordered_map<int, double> zmap;
	bool valid_flag;
};

#endif