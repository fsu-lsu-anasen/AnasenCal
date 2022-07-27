/*
	ChannelMap
	Wrapper around std::unordered_map which gives acces to channel information using a global channel number key.
	Global channel number is given by chipboard_num*32 + channel (with offsets to account for changes in MB).
	This particular version is for the solid-target ANASEN detector array,

	Key method is FindChannel. This returns an iterator, which can be checked against the end of the map for validity.
	If the iterator is equal to the end, the key does not link to a member of the map.

	Written by Gordon McCann Nov 2021
*/
#ifndef CHANNELMAP_H
#define CHANNELMAP_H

#include <string>
#include <unordered_map>

struct ChannelData
{
	std::string detectorType="None"; //Identifier such as FQQQ, BARREL1A, etc.
	int detectorID = -1;
	std::string detectorComponent="None"; //RING/WEDGE or FRONT/FRONTUP/FRONTDOWN/BACK
	int localChannel=0; //Local channel in detector

	bool operator==(const ChannelData& rhs)
	{
		if(this->detectorType == rhs.detectorType && this->detectorID == rhs.detectorID && this->detectorComponent == rhs.detectorComponent
			&& this->localChannel == rhs.localChannel)
		{
			return true;
		}
		return false;
	}
};

class ChannelMap
{
public:
	typedef std::unordered_map<int, ChannelData>::iterator Iterator;

	ChannelMap(const std::string& filename);
	~ChannelMap();
	inline const bool IsValid() const { return valid_flag; }
	inline Iterator FindChannel(int gchan) { return cmap.find(gchan); }
	inline Iterator End() { return cmap.end(); }
	int InverseFindChannel(const ChannelData& data);

private:
	void FillMap(const std::string& filename);
	std::unordered_map<int, ChannelData> cmap;
	bool valid_flag;
	std::string name;
};

#endif
