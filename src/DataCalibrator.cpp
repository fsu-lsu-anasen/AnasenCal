/*
	DataCalibrator
	Class which applies all of the calibration results to a data set, performs front-back hit assignment (kind of)
	and saves to a condensed data format (CalibratedEvent) for further analysis. This essentially a single function
	with some class-global parameter maps.

	Written by Gordon McCann Nov 2021
*/
#include "DataCalibrator.h"
#include "DataStructs.h"
#include <iostream>

#include <TFile.h>
#include <TTree.h>

//Requires a file from each calibration stage
DataCalibrator::DataCalibrator(const std::string& channelfile, const std::string& backmatch, const std::string& updownmatch,
				const std::string& frontbackmatch, const std::string& energyfile) :
	channel_map(channelfile), back_map(backmatch), updown_map(updownmatch), frontback_map(frontbackmatch), energy_map(energyfile)
{
}

DataCalibrator::~DataCalibrator() {}

/*
	Main loop. Takes in an input data fle, and an output data file. These should both be ROOT formated, where input data should be of AnasenEvent
	type, and the output will be saved as CalibratedEvent data.
*/
void DataCalibrator::Run(const std::string& inputname, const std::string& outputname)
{
	if(!channel_map.IsValid() || !back_map.IsValid() || !updown_map.IsValid() || !frontback_map.IsValid() || !energy_map.IsValid())
	{
		std::cerr<<"Bad maps at DataCalibrator::Run()! Exiting."<<std::endl;
		return;
	}	

	TFile* input = TFile::Open(inputname.c_str(), "READ");
	if(!input->IsOpen())
	{
		std::cerr<<"Unable to open input file "<<inputname<<" at DataCalibrator::Run()! Exiting."<<std::endl;
		return;
	}
	TTree* intree = (TTree*) input->Get("EventTree");

	CoincEvent* event = new CoincEvent();
	intree->SetBranchAddress("event", &event);

	TFile* output = TFile::Open(outputname.c_str(), "RECREATE");
	if(!output->IsOpen())
	{
		std::cerr<<"Unable to open output file "<<outputname<<" at DataCalibrator::Run()! Exiting."<<std::endl;
		return;
	}
	TTree* outtree = new TTree("CalTree", "CalTree");

	CalibratedEvent calevent, blank_event;
	outtree->Branch("event", &calevent);

	long nentries = intree->GetEntries();

	long count=0, flush_count=0, flush_val = 0.01*nentries;

	CalibratedSX3Hit sx3hit, blank_sx3;
	CalibratedQQQHit qqqhit, blank_qqq;
	CalibratedBarcHit barchit, blank_barc;
	double cal_back, cal_up_energy, cal_down_energy, cal_sum;
	for(long i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			count=0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count<<"%"<<std::flush;
		}

		calevent = blank_event;

		/*
			Where the fun happens. We want to convert to a condensed format consisting of assiciated information
			(front & back data) that makes up a particle hit. To this end: first look for a good back hit, then search
			for a front (either QQQ ring or SX3 up-down pair as appropriate) and write to a CalibratedEvent. Looks scarier
			than it is

			NOTE: As currently implemented a front is NOT required to make a good hit. This is due primarily to the 
			poor SX3 front efficiency
		*/

		for(int j=0; j<12; j++)
		{
			for(auto& backhit : event->barrel[j].backs)
			{
				sx3hit = blank_sx3;
				auto backgains = back_map.FindParameters(backhit.globalChannel);
				auto backcal = energy_map.FindParameters(backhit.globalChannel);
				if(backgains == back_map.End() || backcal == energy_map.End())
					continue;
				sx3hit.back_energy = backcal->second.slope*(backgains->second.slope*(backhit.energy) + backgains->second.intercept) + backcal->second.intercept;
				sx3hit.back_gchan = backhit.globalChannel;
				sx3hit.detector_index = j;				 
				for(auto& fuphit : event->barrel[j].frontsUp)
				{
					auto fup_channel = channel_map.FindChannel(fuphit.globalChannel);
					for(auto& fdownhit : event->barrel[j].frontsDown)
					{
						auto fdown_channel = channel_map.FindChannel(fdownhit.globalChannel);
						if(fup_channel->second.localChannel != updown_list[fdown_channel->second.localChannel])
							continue;
						auto upgains = updown_map.FindParameters(fuphit.globalChannel);
						if(upgains == updown_map.End())
							continue;
						auto frontbackgains = frontback_map.FindParameters(fuphit.globalChannel);
						if(frontbackgains == frontback_map.End())
							continue;
		
						cal_back = backgains->second.slope*(backhit.energy) + backgains->second.intercept;
						cal_up_energy = cal_back - upgains->second.slope*(fuphit.energy) - upgains->second.intercept*cal_back;
						cal_down_energy = fdownhit.energy;
						cal_sum = frontbackgains->second.slope*(cal_down_energy+cal_up_energy) + frontbackgains->second.intercept;
						if(cal_sum/cal_back > 1.2 || cal_sum/cal_back < 0.8)
						{
							continue;
						}
						else
						{
							sx3hit.frontup_energy_adc = cal_up_energy;
							sx3hit.frontdown_energy_adc = cal_down_energy;
							sx3hit.frontup_gchan = fuphit.globalChannel;
							sx3hit.frontdown_gchan = fdownhit.globalChannel;
							break;
						}
					}
					calevent.barrel.push_back(sx3hit);
				}
			}
		}

		for(int j=0; j<4; j++)
		{
			for(auto& wedgehit : event->fqqq[j].wedges)
			{
				qqqhit = blank_qqq;
				auto backgains = back_map.FindParameters(wedgehit.globalChannel);
				auto backcal = energy_map.FindParameters(wedgehit.globalChannel);
				if(backgains == back_map.End() || backcal == energy_map.End())
					continue;
				
				qqqhit.wedge_energy = backcal->second.slope*(backgains->second.slope*(wedgehit.energy) + backgains->second.intercept) + backcal->second.intercept;
				qqqhit.wedge_gchan = wedgehit.globalChannel;
				qqqhit.detector_index = j;
				for(auto& ringhit : event->fqqq[j].rings)
				{
					auto frontgains = frontback_map.FindParameters(ringhit.globalChannel);
					auto frontcal = energy_map.FindParameters(ringhit.globalChannel);
					if(frontgains == frontback_map.End() || frontcal == energy_map.End())
						continue;
					cal_up_energy = frontcal->second.slope*(frontgains->second.slope*(ringhit.energy) + frontgains->second.intercept) + frontcal->second.intercept;
					if(cal_up_energy/qqqhit.wedge_energy > 1.2 || cal_up_energy/qqqhit.wedge_energy < 0.8)
						continue;
					else
					{
						qqqhit.ring_energy = cal_up_energy;
						qqqhit.ring_gchan = ringhit.globalChannel;
						break;
					} 
				}
				calevent.fqqq.push_back(qqqhit);
			}
		}

		for(int j=0; j<6; j++)
		{
			for(auto& fronthit : event->barcUp[j].fronts)
			{
				barchit = blank_barc;
				auto frontcal = energy_map.FindParameters(fronthit.globalChannel);
				if(frontcal == energy_map.End())
					continue;

				barchit.front_energy = frontcal->second.slope*fronthit.energy + frontcal->second.intercept;
				barchit.front_gchan = fronthit.globalChannel;
				barchit.detector_index = j;
				calevent.barcUp.push_back(barchit);
			}
			for(auto& fronthit : event->barcDown[j].fronts)
			{
				barchit = blank_barc;
				auto frontcal = energy_map.FindParameters(fronthit.globalChannel);
				if(frontcal == energy_map.End())
					continue;

				barchit.front_energy = frontcal->second.slope*fronthit.energy + frontcal->second.intercept;
				barchit.front_gchan = fronthit.globalChannel;
				barchit.detector_index = j;
				calevent.barcDown.push_back(barchit);
			}
		}

		if(calevent.fqqq.size() + calevent.barrel.size() + calevent.barcUp.size() + calevent.barcDown.size() > 0)
			outtree->Fill();

	}
	std::cout<<std::endl;

	input->Close();
	output->cd();
	outtree->Write(outtree->GetName(), TObject::kOverwrite);
	output->Close();
	delete event;
}
