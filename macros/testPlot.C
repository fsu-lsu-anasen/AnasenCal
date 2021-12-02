#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <THashTable.h>
#include <TH1.h>
#include <TH2.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "../include/DataStructs.h"

#include "../include/ChannelMap.h"

R__LOAD_LIBRARY(../objs/libAnasenEvent_dict.so)

ChannelMap::ChannelMap(const std::string& filename) :
	valid_flag(false), name(filename)
{
	FillMap(filename);
}

ChannelMap::~ChannelMap() {}

void ChannelMap::FillMap(const std::string& filename)
{
	std::ifstream input(filename);
	if(!input.is_open())
	{
		valid_flag=false;
		std::cerr<<"Bad file at ChannelMap::FillMap!"<<std::endl;
		return;
	}

	ChannelData data;
	int gchan;
	while(input>>gchan)
	{
		input>>data.detectorType>>data.detectorID>>data.detectorComponent>>data.detectorDirection>>data.channel;
		cmap[gchan] = data;
	}

	valid_flag = true;
}

int ChannelMap::ConvertSX3Name2Index(const std::string& detectorType, const std::string& detectorID)
{
	if(detectorType == "BARREL1A" || detectorType == "BARREL2A")
	{
		if(detectorID == "A")
		{
			return 0;
		}
		else if(detectorID == "B")
		{
			return 1;
		}
		else if(detectorID == "C")
		{
			return 2;
		}
		else if(detectorID == "D")
		{
			return 3;
		}
		else if(detectorID == "E")
		{
			return 4;
		}
		else if(detectorID == "F")
		{
			return 5;
		}
	} 
	else if(detectorType == "BARREL1B" || detectorType == "BARREL2B")
	{
		if(detectorID == "A")
		{
			return 6;
		}
		else if(detectorID == "B")
		{
			return 7;
		}
		else if(detectorID == "C")
		{
			return 8;
		}
		else if(detectorID == "D")
		{
			return 9;
		}
		else if(detectorID == "E")
		{
			return 10;
		}
		else if(detectorID == "F")
		{
			return 11;
		}
	}
	
	std::cerr<<"Invalid detectorType "<<detectorType<<" or detectorID "<<detectorID<<" at ChannelMap::ConvertSX3Name2Index! Returning -1."<<std::endl;
	return -1;
}

int ChannelMap::ConvertQQQName2Index(const std::string& detectorID)
{
	if(detectorID == "0")
		return 0;
	else if (detectorID == "1")
		return 1;
	else if (detectorID == "2")
		return 2;
	else if (detectorID == "3")
		return 3;

	std::cerr<<"Invalid detectorID "<<detectorID<<" at ChannelMap::CovnertQQQName2Index! Returning -1."<<std::endl;
	return -1;
}

int ChannelMap::InverseFindChannel(const ChannelData& data)
{
	for(auto& channel : cmap)
	{
		ChannelData& this_data = channel.second;
		if(this_data == data)
			return channel.first;
	}
	return -1;
}

void MyFill(THashTable* table, const std::string& name, const std::string& title, int bins, double minx, double maxx, double value)
{
	TH1* h = (TH1*) table->FindObject(name.c_str());
	if(h == nullptr)
	{	
		TH1F* histo = new TH1F(name.c_str(), title.c_str(), bins, minx, maxx);
		histo->Fill(value);
		table->Add(histo);
	}
	else
	{
		h->Fill(value);
	}
}

void MyFill(THashTable* table, const std::string& name, const std::string& title, int binsx, double minx, double maxx, double valuex,
			int binsy, double miny, double maxy, double valuey)
{
	TH2* h = (TH2*) table->FindObject(name.c_str());
	if(h == nullptr)
	{	
		TH2F* histo = new TH2F(name.c_str(), title.c_str(), binsx, minx, maxx, binsy, miny, maxy);
		histo->Fill(valuex, valuey);
		table->Add(histo);
	}
	else
	{
		h->Fill(valuex, valuey);
	}
}

void testPlot()
{
	TFile* file = TFile::Open("/data1/gwm17/7BeNov2021/merged/alphaData_updown.root", "READ");
	TTree* tree = (TTree*) file->Get("EventTree");

	AnasenEvent* event = new AnasenEvent();

	tree->SetBranchAddress("event", &event);

	TFile* outfile = TFile::Open("/data1/gwm17/7BeNov2021/histograms/ChannelChecking.root", "RECREATE");

	THashTable* table = new THashTable();

	ChannelMap cmap("../etc/AnasenChannelMap_fixedOrientation.txt");

	std::vector<int> fired_channels;
	std::vector<int> total_counts;
	std::vector<std::vector<int>> channel_matrix;
	std::vector<std::vector<int>> improper_channel_matrix;
	channel_matrix.resize(544);
	total_counts.resize(544);
	for(auto& list : channel_matrix)
		list.resize(545);

	for(auto& list : channel_matrix)
		for(auto& entry : list)
			entry=0;

	improper_channel_matrix = channel_matrix;

	for(auto& chan : total_counts)
		chan = 0;

	int small, large;

	std::string name;
	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);

		fired_channels.clear();

		for(int j=0; j<12; j++)
		{
			for(auto& hit : event->barrel1[j].fronts_up)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table,name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
			for(auto& hit : event->barrel1[j].fronts_down)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table,name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
			for(auto& hit : event->barrel1[j].backs)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table,name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
			for(auto& hit : event->barrel2[j].fronts_up)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table,name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
			for(auto& hit : event->barrel2[j].fronts_down)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table,name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
			for(auto& hit : event->barrel2[j].backs)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table, name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
		}

		for(int j=0; j<4; j++)
		{
			for(auto& hit : event->fqqq[j].rings)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table, name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
			for(auto& hit : event->fqqq[j].wedges)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table, name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
			for(auto& hit : event->bqqq[j].rings)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table, name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
			for(auto& hit : event->bqqq[j].wedges)
			{
				if(hit.energy > 700 && hit.energy < 8000)
				{
					name = "channel_"+std::to_string(hit.global_chan);
					fired_channels.push_back(hit.global_chan);
					total_counts[hit.global_chan]++;
					MyFill(table, name.c_str(), name.c_str(), 16384,0,16384,hit.energy);
				}
			}
		}

		for(unsigned int j=0; j<fired_channels.size(); j++)
		{
			for(unsigned int k=(j+1); k<fired_channels.size(); k++)
			{
				small = std::min(fired_channels[j], fired_channels[k]);
				large = std::max(fired_channels[j], fired_channels[k]);
				MyFill(table, "channel_correlations","Channel Correlations;global channel;global channel",545,-1,544,small,545,-1,544,large);
				channel_matrix[small][large]++;

				auto channel1 = cmap.FindChannel(small);
				auto channel2 = cmap.FindChannel(large);

				if(channel1 == cmap.End() || channel2 == cmap.End())
				{
					std::cerr<<"Bad channel at fixed chan analysis...."<<std::endl;
					continue;
				}

				if(channel1->second.detectorType == channel2->second.detectorType && channel1->second.detectorID == channel2->second.detectorID)
					continue;
				improper_channel_matrix[small][large]++;

			}
		}

		if(fired_channels.size() == 1)
		{
			MyFill(table, "channel_correlations","Channel Correlations;global channel;global channel",545,-1,544,fired_channels[0],545,-1,544,-1);
			channel_matrix[fired_channels[0]][544]++;
			improper_channel_matrix[fired_channels[0]][544]++;
		}

	}

	std::ofstream output("ChannelMatrix.txt");
	for(unsigned int i=0; i<channel_matrix.size(); i++)
	{
		output<<"----------------------------------------------------------"<<std::endl;
		output<<std::setw(3)<<i<<" | ";
		for(unsigned int j=0; j<channel_matrix[i].size(); j++)
		{
			if(channel_matrix[i][j] < 50)
				continue;
			output<<j<<":"<<channel_matrix[i][j]<<", ";
		}
		output<<"total counts:"<<total_counts[i]<<std::endl;
		output<<"----------------------------------------------------------"<<std::endl;
	}
	output.close();

	std::ofstream imoutput("ImproperChannelMatrix.txt");
	for(unsigned int i=0; i<improper_channel_matrix.size(); i++)
	{
		imoutput<<"----------------------------------------------------------"<<std::endl;
		imoutput<<std::setw(3)<<i<<" | ";
		for(unsigned int j=0; j<improper_channel_matrix[i].size(); j++)
		{
			if(improper_channel_matrix[i][j] < 1000)
				continue;
			imoutput<<j<<":"<<improper_channel_matrix[i][j]<<", ";
		}
		imoutput<<"total counts:"<<total_counts[i]<<std::endl;
		imoutput<<"----------------------------------------------------------"<<std::endl;
	}
	imoutput.close();
	file->Close();

	outfile->cd();
	table->Write();
	outfile->Close();

}
