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
	int nbqqq_ws=0, nbqqq_ws_matched=0;
	int nbqqq0_ws=0, nbqqq0_ws_matched=0;
	int nbqqq1_ws=0, nbqqq1_ws_matched=0;
	int nbqqq2_ws=0, nbqqq2_ws_matched=0;
	int nbqqq3_ws=0, nbqqq3_ws_matched=0;
	int nfqqq_ws=0, nfqqq_ws_matched=0;
	int nfqqq0_ws=0, nfqqq0_ws_noRings=0, nfqqq0_ws_manyRings=0, nfqqq0_ws_oneRing=0, nfqqq0_ws_matched=0;
	int nfqqq1_ws=0, nfqqq1_ws_noRings=0, nfqqq1_ws_manyRings=0, nfqqq1_ws_oneRing=0, nfqqq1_ws_matched=0;
	int nfqqq2_ws=0, nfqqq2_ws_noRings=0, nfqqq2_ws_manyRings=0, nfqqq2_ws_oneRing=0, nfqqq2_ws_matched=0;
	int nfqqq3_ws=0, nfqqq3_ws_noRings=0, nfqqq3_ws_manyRings=0, nfqqq3_ws_oneRing=0, nfqqq3_ws_matched=0;
	int nfqqq0_ringsOnly=0;
	int nfqqq1_ringsOnly=0;
	int nfqqq2_ringsOnly=0;
	int nfqqq3_ringsOnly=0;
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
					/* This should not have been here
					if(sx3hit.frontdown_gchan != -1 && sx3hit.frontup_gchan !=  -1)
					{
						break;
					}
					*/
					calevent.barrel1.push_back(sx3hit);
				}
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

					/* This should not have been there
					if(sx3hit.frontdown_gchan != -1 && sx3hit.frontup_gchan !=  -1)
					{
						break;
					}
					*/
					calevent.barrel2.push_back(sx3hit); //This should be here....
				}
			}
		}

		for(int j=0; j<4; j++)
		{
			/*
			if(event->fqqq[j].wedges.size() != 0)
				std::cout<<"looking over fqqq"<<j<<std::endl;
			*/
			if(event->fqqq[j].wedges.size() == 0 && event->fqqq[j].rings.size() != 0)
			{
				switch(j)
				{
					case 0: nfqqq0_ringsOnly++; break;
					case 1: nfqqq1_ringsOnly++; break;
					case 2: nfqqq2_ringsOnly++; break;
					case 3: nfqqq3_ringsOnly++; break;
				}
			} 
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
				if(qqqhit.wedge_energy < 2.8)
					continue;
				nfqqq_ws++;
				switch(j)
				{
					case 0: nfqqq0_ws++; break;
					case 1: nfqqq1_ws++; break;
					case 2: nfqqq2_ws++; break;
					case 3: nfqqq3_ws++; break;
				}
				/*
				if(j==1)
				{
					std::cout<<"attempting to match fqqq"<<j<<" wedge with energy "<<qqqhit.wedge_energy<<" time "<<wedgehit.time<<" and gchan "<<qqqhit.wedge_gchan<<std::endl;
					std::cout<<"number of candidates "<<event->fqqq[j].rings.size()<<std::endl;
				}
				*/
				for(auto& ringhit : event->fqqq[j].rings)
				{
					auto frontoffset = zero_map.FindOffset(ringhit.global_chan);
					auto frontgains = frontback_map.FindParameters(ringhit.global_chan);
					auto frontcal = energy_map.FindParameters(ringhit.global_chan);
					if(frontoffset == zero_map.End() || frontgains == frontback_map.End() || frontcal == energy_map.End())
						continue;
					cal_up_energy = frontcal->second.slope*(frontgains->second.slope*(ringhit.energy - frontoffset->second) + frontgains->second.intercept) + frontcal->second.intercept;
					//if(j == 1)std::cout<<"ring candidate energy "<<cal_up_energy<<" time "<<ringhit.time<<" gchan "<<ringhit.global_chan<<std::endl;
					if(cal_up_energy/qqqhit.wedge_energy > 1.2 || cal_up_energy/qqqhit.wedge_energy < 0.8)
						continue;
					else
					{
						//if(j==1)std::cout<<"matched!"<<std::endl;
						nfqqq_ws_matched++;
						switch(j)
						{
							case 0: nfqqq0_ws_matched++; break;
							case 1: nfqqq1_ws_matched++; break;
							case 2: nfqqq2_ws_matched++; break;
							case 3: nfqqq3_ws_matched++; break;
						}
						qqqhit.ring_energy = cal_up_energy;
						qqqhit.ring_gchan = ringhit.global_chan;
						break;
					} 
				}
				calevent.fqqq.push_back(qqqhit);
				if(qqqhit.ring_gchan == -1)
				{
					switch(j)
					{
						case 0:
						{
							if(event->fqqq[j].rings.size() == 0)
								nfqqq0_ws_noRings++;
							else if(event->fqqq[j].rings.size() == 1)
								nfqqq0_ws_oneRing++;
							else if(event->fqqq[j].rings.size() > 1)
								nfqqq0_ws_manyRings++;
							break;
						}
						case 1:
						{
							if(event->fqqq[j].rings.size() == 0)
								nfqqq1_ws_noRings++;
							else if(event->fqqq[j].rings.size() == 1)
								nfqqq1_ws_oneRing++;
							else if(event->fqqq[j].rings.size() > 1)
								nfqqq1_ws_manyRings++;
							break;
						}
						case 2:
						{
							if(event->fqqq[j].rings.size() == 0)
								nfqqq2_ws_noRings++;
							else if(event->fqqq[j].rings.size() == 1)
								nfqqq2_ws_oneRing++;
							else if(event->fqqq[j].rings.size() > 1)
								nfqqq2_ws_manyRings++;
							break;
						}
						case 3:
						{
							if(event->fqqq[j].rings.size() == 0)
								nfqqq3_ws_noRings++;
							else if(event->fqqq[j].rings.size() == 1)
								nfqqq3_ws_oneRing++;
							else if(event->fqqq[j].rings.size() > 1)
								nfqqq3_ws_manyRings++;
							break;
						}
					}
				}
			}

			for(auto& wedgehit : event->bqqq[j].wedges)
			{
				qqqhit = blank_qqq;
				auto backoffset = zero_map.FindOffset(wedgehit.global_chan);
				auto backgains = back_map.FindParameters(wedgehit.global_chan);
				auto backcal = energy_map.FindParameters(wedgehit.global_chan);
				if(backgains == back_map.End() || backcal == energy_map.End() || backoffset == zero_map.End())
					continue;
				nbqqq_ws++;
				switch(j)
				{
					case 0: nbqqq0_ws++; break;
					case 1: nbqqq1_ws++; break;
					case 2: nbqqq2_ws++; break;
					case 3: nbqqq3_ws++; break;
				}
				qqqhit.wedge_energy = backcal->second.slope*(backgains->second.slope*(wedgehit.energy - backoffset->second) + backgains->second.intercept) + backcal->second.intercept;
				qqqhit.wedge_gchan = wedgehit.global_chan;
				qqqhit.detector_index = j;
				for(auto& ringhit : event->bqqq[j].rings)
				{
					auto frontoffset = zero_map.FindOffset(ringhit.global_chan);
					auto frontgains = frontback_map.FindParameters(ringhit.global_chan);
					auto frontcal = energy_map.FindParameters(ringhit.global_chan);
					if(frontoffset == zero_map.End() || frontgains == frontback_map.End() || frontcal == energy_map.End())
						continue;
					cal_up_energy = frontcal->second.slope*(frontgains->second.slope*(ringhit.energy - frontoffset->second) + frontgains->second.intercept) + frontcal->second.intercept;
					if(cal_up_energy/qqqhit.wedge_energy > 1.2 || cal_up_energy/qqqhit.wedge_energy < 0.8)
						continue;
					else
					{
						nbqqq_ws_matched++;
						switch(j)
						{
							case 0: nbqqq0_ws_matched++; break;
							case 1: nbqqq1_ws_matched++; break;
							case 2: nbqqq2_ws_matched++; break;
							case 3: nbqqq3_ws_matched++; break;
						}
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

	std::cout<<"nbqqq_ws: "<<nbqqq_ws<<" matched: "<<nbqqq_ws_matched<<std::endl;
	std::cout<<"nbqqq0_ws: "<<nbqqq0_ws<<" matched: "<<nbqqq0_ws_matched<<std::endl;
	std::cout<<"nbqqq1_ws: "<<nbqqq1_ws<<" matched: "<<nbqqq1_ws_matched<<std::endl;
	std::cout<<"nbqqq2_ws: "<<nbqqq2_ws<<" matched: "<<nbqqq2_ws_matched<<std::endl;
	std::cout<<"nbqqq3_ws: "<<nbqqq3_ws<<" matched: "<<nbqqq3_ws_matched<<std::endl;
	std::cout<<"nfqqq_ws: "<<nfqqq_ws<<" matched: "<<nfqqq_ws_matched<<std::endl;
	std::cout<<"nfqqq0_ws: "<<nfqqq0_ws<<" matched: "<<nfqqq0_ws_matched<<" no rings present: "<<nfqqq0_ws_noRings;
	std::cout<<" many rings present: "<<nfqqq0_ws_manyRings<<" one ring present: "<<nfqqq0_ws_oneRing<<std::endl;
	std::cout<<"nfqqq1_ws: "<<nfqqq1_ws<<" matched: "<<nfqqq1_ws_matched<<" no rings present: "<<nfqqq1_ws_noRings;
	std::cout<<" many rings present: "<<nfqqq1_ws_manyRings<<" one ring present: "<<nfqqq1_ws_oneRing<<std::endl;
	std::cout<<"nfqqq2_ws: "<<nfqqq2_ws<<" matched: "<<nfqqq2_ws_matched<<" no rings present: "<<nfqqq2_ws_noRings;
	std::cout<<" many rings present: "<<nfqqq2_ws_manyRings<<" one ring present: "<<nfqqq2_ws_oneRing<<std::endl;
	std::cout<<"nfqqq3_ws: "<<nfqqq3_ws<<" matched: "<<nfqqq3_ws_matched<<" no rings present: "<<nfqqq3_ws_noRings;
	std::cout<<" many rings present: "<<nfqqq3_ws_manyRings<<" one ring present: "<<nfqqq3_ws_oneRing<<std::endl;
	std::cout<<"nfqqq0_ringsOnly: "<<nfqqq0_ringsOnly<<std::endl;
	std::cout<<"nfqqq1_ringsOnly: "<<nfqqq1_ringsOnly<<std::endl;
	std::cout<<"nfqqq2_ringsOnly: "<<nfqqq2_ringsOnly<<std::endl;
	std::cout<<"nfqqq3_ringsOnly: "<<nfqqq3_ringsOnly<<std::endl;

	input->Close();
	output->cd();
	outtree->Write(outtree->GetName(), TObject::kOverwrite);
	output->Close();
	delete event;
}
