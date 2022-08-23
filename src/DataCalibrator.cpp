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
	m_channelMap(channelfile), m_backGainMap(backmatch), m_sx3UpDownGainMap(updownmatch), m_frontBackGainMap(frontbackmatch), m_energyCalMap(energyfile)
{
}

DataCalibrator::~DataCalibrator() {}

/*
	Main loop. Takes in an input data fle, and an output data file. These should both be ROOT formated, where input data should be of AnasenEvent
	type, and the output will be saved as CalibratedEvent data.
*/
void DataCalibrator::Run(const std::string& inputname, const std::string& outputname)
{
	if(!m_channelMap.IsValid() || !m_backGainMap.IsValid() || !m_sx3UpDownGainMap.IsValid() || !m_frontBackGainMap.IsValid() || !m_energyCalMap.IsValid())
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
				auto backchannel = m_channelMap.FindChannel(backhit.globalChannel);
				auto backgains = m_backGainMap.FindParameters(backhit.globalChannel);
				auto backcal = m_energyCalMap.FindParameters(backhit.globalChannel);
				if(backgains == m_backGainMap.End() || backcal == m_energyCalMap.End())
					continue;
				sx3hit.backEnergy = backcal->second.slope*(backgains->second.slope*(backhit.energy) + backgains->second.intercept) +
									 backcal->second.intercept;
				sx3hit.backGlobalChannel = backhit.globalChannel;
				sx3hit.backLocalChannel = backchannel->second.localChannel;
				sx3hit.detectorIndex = j;				 
				for(auto& fuphit : event->barrel[j].frontsUp)
				{
					auto fup_channel = m_channelMap.FindChannel(fuphit.globalChannel);
					for(auto& fdownhit : event->barrel[j].frontsDown)
					{
						auto fdown_channel = m_channelMap.FindChannel(fdownhit.globalChannel);
						if(fup_channel->second.localChannel != s_sx3UpDownMatch[fdown_channel->second.localChannel])
							continue;
						auto upgains = m_sx3UpDownGainMap.FindParameters(fuphit.globalChannel);
						if(upgains == m_sx3UpDownGainMap.End())
							continue;
						auto frontbackgains = m_frontBackGainMap.FindParameters(fuphit.globalChannel);
						if(frontbackgains == m_frontBackGainMap.End())
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
							sx3hit.frontUpEnergyAdc = cal_up_energy;
							sx3hit.frontDownEnergyAdc = cal_down_energy;
							sx3hit.frontUpGlobalChannel = fuphit.globalChannel;
							sx3hit.frontDownGlobalChannel = fdownhit.globalChannel;
							sx3hit.frontUpLocalChannel = fup_channel->second.localChannel;
							sx3hit.frontDownLocalChannel = fdown_channel->second.localChannel;
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
				auto backchannel = m_channelMap.FindChannel(wedgehit.globalChannel);
				auto backgains = m_backGainMap.FindParameters(wedgehit.globalChannel);
				auto backcal = m_energyCalMap.FindParameters(wedgehit.globalChannel);
				if(backgains == m_backGainMap.End() || backcal == m_energyCalMap.End())
					continue;
				
				qqqhit.wedgeEnergy = backcal->second.slope*(backgains->second.slope*(wedgehit.energy) + backgains->second.intercept) +
									  backcal->second.intercept;
				qqqhit.wedgeGlobalChannel = wedgehit.globalChannel;
				qqqhit.wedgeLocalChannel = backchannel->second.localChannel;
				qqqhit.detectorIndex = j;
				for(auto& ringhit : event->fqqq[j].rings)
				{
					auto frontchannel = m_channelMap.FindChannel(ringhit.globalChannel);
					auto frontgains = m_frontBackGainMap.FindParameters(ringhit.globalChannel);
					auto frontcal = m_energyCalMap.FindParameters(ringhit.globalChannel);
					if(frontgains == m_frontBackGainMap.End() || frontcal == m_energyCalMap.End())
						continue;
					cal_up_energy = frontcal->second.slope*(frontgains->second.slope*(ringhit.energy) + frontgains->second.intercept) +
									frontcal->second.intercept;
					if(cal_up_energy/qqqhit.wedgeEnergy > 1.2 || cal_up_energy/qqqhit.wedgeEnergy < 0.8)
						continue;
					else
					{
						qqqhit.ringEnergy = cal_up_energy;
						qqqhit.ringGlobalChannel = ringhit.globalChannel;
						qqqhit.ringLocalChannel = frontchannel->second.localChannel;
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
				auto frontchannel = m_channelMap.FindChannel(fronthit.globalChannel);
				auto frontcal = m_energyCalMap.FindParameters(fronthit.globalChannel);
				if(frontcal == m_energyCalMap.End())
					continue;

				barchit.frontEnergy = frontcal->second.slope*fronthit.energy + frontcal->second.intercept;
				barchit.frontGlobalChannel = fronthit.globalChannel;
				barchit.frontLocalChannel = frontchannel->second.localChannel;
				barchit.detectorIndex = j;
				calevent.barcUp.push_back(barchit);
			}
			for(auto& fronthit : event->barcDown[j].fronts)
			{
				barchit = blank_barc;
				auto frontchannel = m_channelMap.FindChannel(fronthit.globalChannel);
				auto frontcal = m_energyCalMap.FindParameters(fronthit.globalChannel);
				if(frontcal == m_energyCalMap.End())
					continue;

				barchit.frontEnergy = frontcal->second.slope*fronthit.energy + frontcal->second.intercept;
				barchit.frontGlobalChannel = fronthit.globalChannel;
				barchit.frontLocalChannel = frontchannel->second.localChannel;
				barchit.detectorIndex = j;
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
