/*
	DataOrganizer
	Class to convert data from the raw motherboard-channel-based structure of evt2root to
	a more practical detector-based structure, AnasenEvent. This requires that a channel map
	exists, otherwise the conversion is not possible.

	Written by Gordon McCann Nov. 2021
*/
#include "DataOrganizer.h"

#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include "ChannelMap.h"
#include "DataStructs.h"

DataOrganizer::DataOrganizer(const std::string& channelfile) :
	cmap(channelfile), generator(new TRandom3())
{
	generator->SetSeed(0);
}

DataOrganizer::~DataOrganizer() 
{
	delete generator;
}

/*
	Method which actually fills the AnasenEvent. The event is passed by reference.
*/
void DataOrganizer::FillEvent(AnasenEvent& event, int gchan, int energy, int time)
{
	if(energy == -1.0)
		return;
	SiliconHit hit;
	auto channel = cmap.FindChannel(gchan);
	if(channel == cmap.End())
	{
		std::cerr<<"Bad global channel "<<gchan<<" at DataOrganizer::FillEvent(). Skipping hit."<<std::endl;
		return;
	}
	hit.global_chan = gchan;
	hit.local_chan = channel->second.channel;
	hit.energy = ConvertInt2Double(energy);
	hit.time = ConvertInt2Double(time);

	if(channel->second.detectorType == "BARREL1A" || channel->second.detectorType == "BARREL1B")
	{
		int detIndex = cmap.ConvertSX3Name2Index(channel->second.detectorType, channel->second.detectorID);
		if(channel->second.detectorComponent == "FRONT" && channel->second.detectorDirection == "UP")
			event.barrel1[detIndex].fronts_up.push_back(hit);
		else if(channel->second.detectorComponent == "FRONT" && channel->second.detectorDirection == "DOWN")
			event.barrel1[detIndex].fronts_down.push_back(hit);
		else if(channel->second.detectorComponent == "BACK")
			event.barrel1[detIndex].backs.push_back(hit);
	}
	else if(channel->second.detectorType == "BARREL2A" || channel->second.detectorType == "BARREL2B")
	{
		int detIndex = cmap.ConvertSX3Name2Index(channel->second.detectorType, channel->second.detectorID);
		if(channel->second.detectorComponent == "FRONT" && channel->second.detectorDirection == "UP")
			event.barrel2[detIndex].fronts_up.push_back(hit);
		else if(channel->second.detectorComponent == "FRONT" && channel->second.detectorDirection == "DOWN")
			event.barrel2[detIndex].fronts_down.push_back(hit);
		else if(channel->second.detectorComponent == "BACK")
			event.barrel2[detIndex].backs.push_back(hit);
	}
	else if(channel->second.detectorType == "FQQQ")
	{
		int detIndex = cmap.ConvertQQQName2Index(channel->second.detectorID);
		if(channel->second.detectorComponent == "RING")
			event.fqqq[detIndex].rings.push_back(hit);
		else if(channel->second.detectorComponent == "WEDGE")
			event.fqqq[detIndex].wedges.push_back(hit);
	}
	else if(channel->second.detectorType == "BQQQ")
	{
		int detIndex = cmap.ConvertQQQName2Index(channel->second.detectorID);
		if(channel->second.detectorComponent == "RING")
			event.bqqq[detIndex].rings.push_back(hit);
		else if(channel->second.detectorComponent == "WEDGE")
			event.bqqq[detIndex].wedges.push_back(hit);
	}
}

/*
	Main loop function. Takes in an input file name and an outputfile name
*/
void DataOrganizer::Run(const std::string& inputname, const std::string& outputname)
{
	if(!cmap.IsValid())
	{
		std::cerr<<"Bad channel map at DataOrganizer::Run()! Exiting."<<std::endl;
		return;
	}


	TFile* input = TFile::Open(inputname.c_str(), "READ");
	TTree* intree = (TTree*) input->Get("DataTree");

	int mb1_energy[9][32];
	int mb2_energy[9][32];
	int mb1_time[9][32];
	int mb2_time[9][32];

	intree->SetBranchAddress("mb1_energy", &mb1_energy);
	intree->SetBranchAddress("mb1_time", &mb1_time);
	intree->SetBranchAddress("mb2_energy", &mb2_energy);
	intree->SetBranchAddress("mb2_time", &mb2_time);

	TFile* output = TFile::Open(outputname.c_str(), "RECREATE");
	TTree* outtree = new TTree("EventTree", "EventTree");

	AnasenEvent event, blank;
	int gchan, mb2_gchan_offset = 9*32;

	outtree->Branch("event", &event);

	int nentries = intree->GetEntries();
	int count=0, flush_count=0, flush_val=0.01*nentries;

	std::cout<<"Orgainizing data into detector structures... Total number of entries: "<<nentries<<std::endl;

	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			count=0;
			flush_count++;
			std::cout<<"\rPercent of data formated: "<<0.01*flush_count*100.0<<"%"<<std::flush;
		}

		event = blank;

		for(int j=0; j<9; j++)
		{
			for(int k=0; k<32; k++)
			{
				gchan = j*32 + k;
				FillEvent(event, gchan, mb1_energy[j][k], mb1_time[j][k]);
			}
		}

		for(int j=0; j<8; j++)
		{
			for(int k=0; k<32; k++)
			{
				gchan = mb2_gchan_offset + j*32 + k;
				FillEvent(event, gchan, mb2_energy[j][k], mb2_time[j][k]);
			}
		}

		outtree->Fill();
	}
	std::cout<<std::endl;

	input->Close();
	output->cd();
	outtree->Write(outtree->GetName(), TObject::kOverwrite);
	output->Close();
}