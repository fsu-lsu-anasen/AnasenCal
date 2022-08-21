/*
	MapChecker
	Simple class which reads in various calibration files and checks the statistics of the different
	calibration stages. Mostly a debug concept.

	Written by Gordon McCann Nov 2021
*/
#include "MapChecker.h"
#include <iostream>
#include "ParameterMap.h"

MapChecker::MapChecker(const std::string& channelfile) :
	m_channelMap(channelfile)
{
}

MapChecker::~MapChecker() {}


void MapChecker::CheckBackGainMatch(const std::string& filename)
{
	std::cout<<"Loading map "<<filename<<"..."<<std::endl;
	ParameterMap paramMap(filename);

	//Counters
	int backs=0, wedges=0;
	int  max_backs=12*4, max_wedges=4*16;

	std::vector<int> missing_channels;
	std::cout<<"Checking channels..."<<std::endl;
	for(int i=0; i<s_nchannels; i++)
	{
		auto channel = m_channelMap.FindChannel(i);
		auto offset = paramMap.FindParameters(i);
		if(offset == paramMap.End())
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
	ParameterMap paramMap(filename);

	//Counters
	int frontups=0;
	int max_fups=12*4;

	std::vector<int> missing_channels;
	std::cout<<"Checking channels..."<<std::endl;
	for(int i=0; i<s_nchannels; i++)
	{
		auto channel = m_channelMap.FindChannel(i);
		auto offset = paramMap.FindParameters(i);
		if(offset == paramMap.End())
		{
			if(channel->second.detectorComponent == "FRONTUP")
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
	ParameterMap paramMap(filename);

	//Counters
	int frontups=0, rings=0;
	int max_fups=12*4, max_rings=4*16;

	std::vector<int> missing_channels;
	std::cout<<"Checking channels..."<<std::endl;
	for(int i=0; i<s_nchannels; i++)
	{
		auto channel = m_channelMap.FindChannel(i);
		auto offset = paramMap.FindParameters(i);
		if(offset == paramMap.End())
		{
			if(channel->second.detectorComponent == "FRONTUP")
			{
				frontups++;
				missing_channels.push_back(i);
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
