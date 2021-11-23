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
DataCalibrator::DataCalibrator(const std::string& channelfile, const std::string& zerofile, const std::string& backmatch, const std::string& updownmatch,
								const std::string& frontbackmatch, const std::string& energyfile) :
	channel_map(channelfile), zero_map(zerofile), back_map(backmatch), updown_map(updownmatch), frontback_map(frontbackmatch), energy_map(energyfile)
{
}

DataCalibrator::~DataCalibrator() {}

/*
	Main loop. Takes in an input data fle, and an output data file. These should both be ROOT formated, where input data should be of AnasenEvent
	type, and the output will be saved as CalibratedEvent data.
*/
void DataCalibrator::Run(const std::string& inputname, const std::string& outputname)
{
	if(!channel_map.IsValid() || !zero_map.IsValid() || !back_map.IsValid() || !updown_map.IsValid() || !frontback_map.IsValid() || !energy_map.IsValid())
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

	AnasenEvent* event = new AnasenEvent();
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
			for(auto& backhit : event->barrel1[j].backs)
			{
				sx3hit = blank_sx3;
				auto backoffset = zero_map.FindOffset(backhit.global_chan);
				auto backgains = back_map.FindParameters(backhit.global_chan);
				auto backcal = energy_map.FindParameters(backhit.global_chan);
				if(backgains == back_map.End() || backcal == energy_map.End() || backoffset == zero_map.End())
					continue;
				sx3hit.back_energy = backcal->second.slope*(backgains->second.slope*(backhit.energy - backoffset->second) + backgains->second.intercept) + backcal->second.intercept;
				sx3hit.back_gchan = backhit.global_chan;
				sx3hit.detector_index = j;				 
				for(auto& fuphit : event->barrel1[j].fronts_up)
				{
					for(auto& fdownhit : event->barrel1[j].fronts_down)
					{
						if(fuphit.local_chan != updown_list[fdownhit.local_chan])
							continue;
						auto upgains = updown_map.FindParameters(fuphit.global_chan);
						if(upgains == updown_map.End())
							continue;
						auto fupzero_offset = zero_map.FindOffset(fuphit.global_chan);
						auto fdownzero_offset = zero_map.FindOffset(fdownhit.global_chan);
						if(fupzero_offset == zero_map.End() || fdownzero_offset == zero_map.End())
							continue;
						auto frontbackgains = frontback_map.FindParameters(fuphit.global_chan);
						if(frontbackgains == frontback_map.End())
							continue;
		
						cal_back = backgains->second.slope*(backhit.energy - backoffset->second) + backgains->second.intercept;
						cal_up_energy = cal_back - upgains->second.slope*(fuphit.energy - fupzero_offset->second) - upgains->second.intercept*cal_back;
						cal_down_energy = fdownhit.energy - fdownzero_offset->second;
						cal_sum = frontbackgains->second.slope*(cal_down_energy+cal_up_energy) + frontbackgains->second.intercept;
						if(cal_sum/cal_back > 1.2 || cal_sum/cal_back < 0.8)
						{
							continue;
						}
						else
						{
							sx3hit.frontup_energy_adc = cal_up_energy;
							sx3hit.frontdown_energy_adc = cal_down_energy;
							sx3hit.frontup_gchan = fuphit.global_chan;
							sx3hit.frontdown_gchan = fdownhit.global_chan;
							break;
						}
					}
					if(sx3hit.frontdown_gchan != -1 && sx3hit.frontup_gchan !=  -1)
					{
						break;
					}
				}
				calevent.barrel1.push_back(sx3hit);
			}

			for(auto& backhit : event->barrel2[j].backs)
			{
				sx3hit = blank_sx3;
				auto backoffset = zero_map.FindOffset(backhit.global_chan);
				auto backgains = back_map.FindParameters(backhit.global_chan);
				auto backcal = energy_map.FindParameters(backhit.global_chan);
				if(backgains == back_map.End() || backcal == energy_map.End() || backoffset == zero_map.End())
					continue;
				sx3hit.back_energy = backcal->second.slope*(backgains->second.slope*(backhit.energy - backoffset->second) + backgains->second.intercept) + backcal->second.intercept;
				sx3hit.back_gchan = backhit.global_chan;
				sx3hit.detector_index = j;				 
				for(auto& fuphit : event->barrel2[j].fronts_up)
				{
					for(auto& fdownhit : event->barrel2[j].fronts_down)
					{
						if(fuphit.local_chan != updown_list[fdownhit.local_chan])
							continue;
						auto upgains = updown_map.FindParameters(fuphit.global_chan);
						if(upgains == updown_map.End())
							continue;
						auto fupzero_offset = zero_map.FindOffset(fuphit.global_chan);
						auto fdownzero_offset = zero_map.FindOffset(fdownhit.global_chan);
						if(fupzero_offset == zero_map.End() || fdownzero_offset == zero_map.End())
							continue;
						auto frontbackgains = frontback_map.FindParameters(fuphit.global_chan);
						if(frontbackgains == frontback_map.End())
							continue;
		
						cal_back = backgains->second.slope*(backhit.energy - backoffset->second) + backgains->second.intercept;
						cal_up_energy = cal_back - upgains->second.slope*(fuphit.energy - fupzero_offset->second) - upgains->second.intercept*cal_back;
						cal_down_energy = fdownhit.energy - fdownzero_offset->second;
						cal_sum = frontbackgains->second.slope*(cal_down_energy+cal_up_energy) + frontbackgains->second.intercept;
						if(cal_sum/cal_back > 1.2 || cal_sum/cal_back < 0.8)
						{
							continue;
						}
						else
						{
							sx3hit.frontup_energy_adc = cal_up_energy;
							sx3hit.frontdown_energy_adc = cal_down_energy;
							sx3hit.frontup_gchan = fuphit.global_chan;
							sx3hit.frontdown_gchan = fdownhit.global_chan;
							break;
						}
					}
					if(sx3hit.frontdown_gchan != -1 && sx3hit.frontup_gchan !=  -1)
					{
						break;
					}
				}
				calevent.barrel2.push_back(sx3hit);
			}
		}

		for(int j=0; j<4; j++)
		{
			for(auto& wedgehit : event->fqqq[j].wedges)
			{
				qqqhit = blank_qqq;
				auto backoffset = zero_map.FindOffset(wedgehit.global_chan);
				auto backgains = back_map.FindParameters(wedgehit.global_chan);
				auto backcal = energy_map.FindParameters(wedgehit.global_chan);
				if(backgains == back_map.End() || backcal == energy_map.End() || backoffset == zero_map.End())
					continue;
				qqqhit.wedge_energy = backcal->second.slope*(backgains->second.slope*(wedgehit.energy - backoffset->second) + backgains->second.intercept) + backcal->second.intercept;
				qqqhit.wedge_gchan = wedgehit.global_chan;
				qqqhit.detector_index = j;
				for(auto& ringhit : event->fqqq[j].rings)
				{
					auto frontoffset = zero_map.FindOffset(ringhit.global_chan);
					auto frontgains = frontback_map.FindParameters(ringhit.global_chan);
					if(frontoffset == zero_map.End() || frontgains == frontback_map.End())
						continue;
					cal_up_energy = frontgains->second.slope*(ringhit.energy - frontoffset->second) + frontgains->second.intercept;
					if(cal_up_energy/qqqhit.wedge_energy > 1.2 || cal_up_energy/qqqhit.wedge_energy < 0.8)
						continue;
					else
					{
						qqqhit.ring_energy = cal_up_energy;
						qqqhit.ring_gchan = ringhit.global_chan;
						break;
					} 
				}
				calevent.fqqq.push_back(qqqhit);
			}

			for(auto& wedgehit : event->bqqq[j].wedges)
			{
				qqqhit = blank_qqq;
				auto backoffset = zero_map.FindOffset(wedgehit.global_chan);
				auto backgains = back_map.FindParameters(wedgehit.global_chan);
				auto backcal = energy_map.FindParameters(wedgehit.global_chan);
				if(backgains == back_map.End() || backcal == energy_map.End() || backoffset == zero_map.End())
					continue;
				qqqhit.wedge_energy = backcal->second.slope*(backgains->second.slope*(wedgehit.energy - backoffset->second) + backgains->second.intercept) + backcal->second.intercept;
				qqqhit.wedge_gchan = wedgehit.global_chan;
				qqqhit.detector_index = j;
				for(auto& ringhit : event->bqqq[j].rings)
				{
					auto frontoffset = zero_map.FindOffset(ringhit.global_chan);
					auto frontgains = frontback_map.FindParameters(ringhit.global_chan);
					if(frontoffset == zero_map.End() || frontgains == frontback_map.End())
						continue;
					cal_up_energy = frontgains->second.slope*(ringhit.energy - frontoffset->second) + frontgains->second.intercept;
					if(cal_up_energy/qqqhit.wedge_energy > 1.2 || cal_up_energy/qqqhit.wedge_energy < 0.8)
						continue;
					else
					{
						qqqhit.ring_energy = cal_up_energy;
						qqqhit.ring_gchan = ringhit.global_chan;
						break;
					} 
				}
				calevent.bqqq.push_back(qqqhit);
			}
		}

		if(calevent.bqqq.size() + calevent.fqqq.size() + calevent.barrel1.size() + calevent.barrel2.size() > 0)
			outtree->Fill();

	}
	std::cout<<std::endl;

	input->Close();
	output->cd();
	outtree->Write(outtree->GetName(), TObject::kOverwrite);
	output->Close();
	delete event;
}
