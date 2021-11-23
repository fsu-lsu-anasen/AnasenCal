/*
	MapChecker
	Simple class which reads in various calibration files and checks the statistics of the different
	calibration stages. Mostly a debug concept.

	Written by Gordon McCann Nov 2021
*/
#include "MapChecker.h"
#include <iostream>
#include "ZeroCalMap.h"
#include "ParameterMap.h"

MapChecker::MapChecker(const std::string& channelfile) :
	cmap(channelfile)
{
}

MapChecker::~MapChecker() {}


void MapChecker::CheckZOffset(const std::string& filename)
{
	std::cout<<"Loading map "<<filename<<"..."<<std::endl;
	ZeroCalMap zmap(filename);

	//Counters
	int frontups=0, frontdowns=0, backs=0, rings=0, wedges=0;
	int max_fups=24*4, max_fdowns=24*4, max_backs=24*4, max_rings=8*16, max_wedges=8*16;

	std::vector<int> missing_channels;
	std::cout<<"Checking channels..."<<std::endl;
	for(int i=0; i<nchannels; i++)
	{
		auto channel = cmap.FindChannel(i);
		auto offset = zmap.FindOffset(i);
		if(offset == zmap.End())
		{
			missing_channels.push_back(i);
			if(channel->second.detectorComponent == "FRONT")
			{
				if(channel->second.detectorDirection == "DOWN")
					frontdowns++;
				else if(channel->second.detectorDirection == "UP")
					frontups++;
			}
			else if(channel->second.detectorComponent == "BACK")
				backs++;
			else if(channel->second.detectorComponent == "RING")
				rings++;
			else if(channel->second.detectorComponent == "WEDGE")
				wedges++;
		}
	}
	std::cout<<"Results..."<<std::endl;
	std::cout<<"Total number of missing channels: "<<missing_channels.size()<<" (percentage: "<<((double)missing_channels.size())/nchannels*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing SX3 Upstream Fronts: "<<frontups<<" (percentage: "<<((double)frontups)/max_fups*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing SX3 Downstream Fronts: "<<frontdowns<<" (percentage: "<<(double)frontdowns/max_fdowns*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing SX3 Backs: "<<backs<<" (percentage: "<<((double)backs)/max_backs*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing QQQ Rings: "<<rings<<" (percentage: "<<((double)rings)/max_rings*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing QQQ Wedges: "<<wedges<<" (percentage: "<<((double)wedges)/max_wedges*100.0<<"%)"<<std::endl;
	std::cout<<"List of missing channels..."<<std::endl;
	for(auto& i : missing_channels)
	{
		std::cout<<i<<std::endl;
	}

}

void MapChecker::CheckBackGainMatch(const std::string& filename)
{
	std::cout<<"Loading map "<<filename<<"..."<<std::endl;
	ParameterMap pmap(filename);

	//Counters
	int backs=0, wedges=0;
	int  max_backs=24*4, max_wedges=8*16;

	std::vector<int> missing_channels;
	std::cout<<"Checking channels..."<<std::endl;
	for(int i=0; i<nchannels; i++)
	{
		auto channel = cmap.FindChannel(i);
		auto offset = pmap.FindParameters(i);
		if(offset == pmap.End())
		{
			if(channel->second.detectorComponent == "BACK")
			{
				backs++;
				missing_channels.push_back(i);
			}
			else if(channel->second.detectorComponent == "WEDGE")
			{
				wedges++;
				missing_channels.push_back(i);
			}
		}
	}
	std::cout<<"Results..."<<std::endl;
	std::cout<<"Total number of missing back channels: "<<missing_channels.size()<<" (percentage: "<<((double)missing_channels.size())/(max_backs+max_wedges)*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing SX3 Backs: "<<backs<<" (percentage: "<<((double)backs)/max_backs*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing QQQ Wedges: "<<wedges<<" (percentage: "<<((double)wedges)/max_wedges*100.0<<"%)"<<std::endl;
	std::cout<<"List of missing channels..."<<std::endl;
	for(auto& i : missing_channels)
	{
		std::cout<<i<<std::endl;
	}

}

void MapChecker::CheckUpDownGainMatch(const std::string& filename)
{
	std::cout<<"Loading map "<<filename<<"..."<<std::endl;
	ParameterMap pmap(filename);

	//Counters
	int frontups=0;
	int max_fups=24*4;

	std::vector<int> missing_channels;
	std::cout<<"Checking channels..."<<std::endl;
	for(int i=0; i<nchannels; i++)
	{
		auto channel = cmap.FindChannel(i);
		auto offset = pmap.FindParameters(i);
		if(offset == pmap.End())
		{
			if(channel->second.detectorComponent == "FRONT" && channel->second.detectorDirection == "UP")
			{
				frontups++;
				missing_channels.push_back(i);
			}
		}
	}
	std::cout<<"Results..."<<std::endl;
	std::cout<<"Number of missing SX3 Upstream Fronts: "<<frontups<<" (percentage: "<<((double)frontups)/max_fups*100.0<<"%)"<<std::endl;
	std::cout<<"List of missing channels..."<<std::endl;
	for(auto& i : missing_channels)
	{
		std::cout<<i<<std::endl;
	}

}

void MapChecker::CheckFrontBackGainMatch(const std::string& filename)
{
	std::cout<<"Loading map "<<filename<<"..."<<std::endl;
	ParameterMap pmap(filename);

	//Counters
	int frontups=0, rings=0;
	int max_fups=24*4, max_rings=8*16;

	std::vector<int> missing_channels;
	std::cout<<"Checking channels..."<<std::endl;
	for(int i=0; i<nchannels; i++)
	{
		auto channel = cmap.FindChannel(i);
		auto offset = pmap.FindParameters(i);
		if(offset == pmap.End())
		{
			if(channel->second.detectorComponent == "FRONT")
			{
				if(channel->second.detectorDirection == "UP")
				{
					frontups++;
					missing_channels.push_back(i);
				}
			}
			else if(channel->second.detectorComponent == "RING")
			{
				rings++;
				missing_channels.push_back(i);
			}
		}
	}
	std::cout<<"Results..."<<std::endl;
	std::cout<<"Total number of missing channels: "<<missing_channels.size()<<" (percentage: "<<((double)missing_channels.size())/(max_fups+max_rings)*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing SX3 Upstream Fronts: "<<frontups<<" (percentage: "<<((double)frontups)/max_fups*100.0<<"%)"<<std::endl;
	std::cout<<"Number of missing QQQ Rings: "<<rings<<" (percentage: "<<((double)rings)/max_rings*100.0<<"%)"<<std::endl;
	std::cout<<"List of missing channels..."<<std::endl;
	for(auto& i : missing_channels)
	{
		std::cout<<i<<std::endl;
	}

}