#include "DeadChannelMap.h"
#include "ChannelMap.h"
#include "ParameterMap.h"
#include <iostream>
#include <fstream>
#include <vector>

void GenerateDeadChannelMap(const std::string& backmatch, const std::string& updownmatch, const std::string& frontbackmatch, const std::string& energycal,
		            const std::string& channelfile, const std::string& deadname)
{

	ChannelMap cmap(channelfile);
	ParameterMap backmap(backmatch), updownmap(updownmatch), frontbackmap(frontbackmatch), energymap(energycal);

	if(!cmap.IsValid() || !backmap.IsValid() || !updownmap.IsValid() || !frontbackmap.IsValid() || !energymap.IsValid())
	{
		std::cerr<<"Bad maps at GenerateDeadChannelMap(). Exiting."<<std::endl;
		return;
	}

	const int nchannels = 640;
	const int updown_list[8] = {1, 0, 3, 2, 5, 4, 7, 6}; //matches index -> index of up/down pair for SX3 fronts

	std::ofstream output(deadname);

	std::vector<bool> dead_flags;
	dead_flags.resize(nchannels);
	for(int i=0; i<nchannels; i++)
		dead_flags[i] = true;

	for(int i=0; i<nchannels; i++)
	{
		auto channel = cmap.FindChannel(i);
		auto backs = backmap.FindParameters(i);
		auto updowns = updownmap.FindParameters(i);
		auto frontbacks = frontbackmap.FindParameters(i);
		auto energy = energymap.FindParameters(i);

		if(channel == cmap.End())
		{
			std::cerr<<"Bad channel "<<i<<" found at GenerateDeadChannelMap!"<<std::endl;
			continue;
		}

		if(channel->second.detectorComponent == "FRONTUP" && updowns != updownmap.End() && frontbacks != frontbackmap.End())
		{
			dead_flags[i] = false;
			ChannelData downdata;
			downdata = channel->second;
			downdata.localChannel = updown_list[channel->second.localChannel];
			downdata.detectorComponent = "FRONTDOWN";
			int down_gchan = cmap.InverseFindChannel(downdata);
			if(down_gchan == -1)
			{
				std::cerr<<"Blergh"<<std::endl;
			}
			dead_flags[down_gchan] = false;
		}
		else if((channel->second.detectorComponent == "BACK" || channel->second.detectorComponent == "WEDGE")
				&& backs != backmap.End() && energy != energymap.End())
		{
			dead_flags[i] = false;
		}
		else if(channel->second.detectorComponent == "RING" && frontbacks != frontbackmap.End() && energy != energymap.End())
		{
			dead_flags[i] = false;
		}
		else if(channel->second.detectorComponent == "FRONT" && energy != energymap.End()) //Barcelonas
		{
			dead_flags[i] = false;
		}
	}

	for(int i=0; i<nchannels; i++)
	{
		if(dead_flags[i])
		{
			auto channel = cmap.FindChannel(i);
			output<<i<<" "<<channel->second.detectorType<<" "<<channel->second.detectorID<<" "<<channel->second.detectorComponent<<" ";
			output<<" "<<channel->second.localChannel<<std::endl;
		}
	}
	output.close();
}
