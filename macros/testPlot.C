#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include "../include/DataStructs.h"

R__LOAD_LIBRARY(../objs/libAnasenEvent_dict.so)

void testPlot()
{
	TFile* file = TFile::Open("/data1/gwm17/7BeNov2021/merged/alphaData_updown.root", "READ");
	TTree* tree = (TTree*) file->Get("EventTree");

	AnasenEvent* event = new AnasenEvent();

	tree->SetBranchAddress("event", &event);

	TH2F* histo = new TH2F("channel_correlations","Channel Correlations;global channel;global channel",545,-1,544,545,-1,544);

	std::vector<int> fired_channels;
	std::vector<std::vector<int>> channel_matrix;
	channel_matrix.resize(544);
	for(auto& list : channel_matrix)
		list.resize(545);

	for(auto& list : channel_matrix)
		for(auto& entry : list)
			entry=0;

	int small, large;

	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);

		fired_channels.clear();

		for(int j=0; j<12; j++)
		{
			for(auto& hit : event->barrel1[j].fronts_up)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
			for(auto& hit : event->barrel1[j].fronts_down)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
			for(auto& hit : event->barrel1[j].backs)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
			for(auto& hit : event->barrel2[j].fronts_up)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
			for(auto& hit : event->barrel2[j].fronts_down)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
			for(auto& hit : event->barrel2[j].backs)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
		}

		for(int j=0; j<4; j++)
		{
			for(auto& hit : event->fqqq[j].rings)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
			for(auto& hit : event->fqqq[j].wedges)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
			for(auto& hit : event->bqqq[j].rings)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
			for(auto& hit : event->bqqq[j].wedges)
				if(hit.energy > 700)
					fired_channels.push_back(hit.global_chan);
		}

		for(unsigned int j=0; j<fired_channels.size(); j++)
		{
			for(unsigned int k=(j+1); k<fired_channels.size(); k++)
			{
				small = std::min(fired_channels[j], fired_channels[k]);
				large = std::max(fired_channels[j], fired_channels[k]);
				histo->Fill(small, large);
				channel_matrix[small][large]++;
			}
		}

		if(fired_channels.size() == 1)
		{
			histo->Fill(fired_channels[0], -1);
			channel_matrix[fired_channels[0]][544]++;
		}

	}

	std::cout<<"Channel Matrix"<<std::endl;
	for(unsigned int i=0; i<545; i++)
	{
		std::cout<<"| "<<i<<" ";
	}
	std::cout<<std::endl;
	for(unsigned int i=0; i<channel_matrix.size(); i++)
	{
		std::cout<<std::setw(3)<<i<<" | ";
		for(unsigned int j=0; j<channel_matrix[i].size(); j++)
		{
			std::cout<<std::setw(5)<<channel_matrix[i][j]<<" | ";
		}
		std::cout<<std::endl;
	}

	histo->Draw("colz");
}
