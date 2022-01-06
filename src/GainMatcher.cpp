/*
	GainMatcher
	Class oriented around performing silicon detector gain matching in ANASEN. Two types of detectors: SX3's and QQQ's.
	There are three gain matching methods: MatchBacks, MatchSX3UpDown, and MatchFrontBack. They should be performed in that order,
	as each step relies upon the previous results.

	Written by Gordon McCann Nov 2021
*/
#include "GainMatcher.h"
#include "ParameterMap.h"
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TF1.h>
#include <TTree.h>
#include <fstream>
#include <iostream>

//For use with std::sort
bool SortGainData(const double i, const double j)
{
	return i<j;
}

GainMatcher::GainMatcher(const std::string& channelfile, const std::string& zerofile) :
	cmap(channelfile), zmap(zerofile)
{
}

GainMatcher::~GainMatcher() {}

//Wrappers around histogram creation, storage, and filling
void GainMatcher::MyFill(THashTable* table, const std::string& name, const std::string& title, int bins, double minx, double maxx, double value)
{
	TH1* histo = (TH1*) table->FindObject(name.c_str());
	if(histo == nullptr)
	{
		TH1F* h = new TH1F(name.c_str(), title.c_str(), bins, minx, maxx);
		h->Fill(value);
		table->Add(h);
	}
	else
	{
		histo->Fill(value);
	}
}

void GainMatcher::MyFill(THashTable* table, const std::string& name, const std::string& title, int binsx, double minx, double maxx, double valuex,
																								int binsy, double miny, double maxy, double valuey)
{
	TH2* histo = (TH2*) table->FindObject(name.c_str());
	if(histo == nullptr)
	{
		TH2F* h = new TH2F(name.c_str(), title.c_str(), binsx, minx, maxx, binsy, miny, maxy);
		h->Fill(valuex, valuey);
		table->Add(h);
	}
	else
	{
		histo->Fill(valuex, valuey);
	}
}

/*
	Method which calls TSpectrum to obtain the number of peaks in a histogram. The peak locations are returned as the
	x-coordinate values of GraphData.
*/
GraphData GainMatcher::GetPoints(THashTable* table, const std::string& name)
{
	GraphData data;
	TH1* histo = (TH1*) table->FindObject(name.c_str());
	if(histo == nullptr)
	{
		//std::cerr<<"Histogram named "<<name<<" not found at ZeroCalibrator::GetPoints! Returning empty data."<<std::endl;
		return data;
	}

	if(histo->Integral() < 1000.0)
	{
		//std::cerr<<"Spectrum "<<name<<" has less than a thousand total counts. Rejected."<<std::endl;
		return data;
	}

	int npeaks = spec.Search(histo, sigma, "nobackground", threshold);
	if(npeaks == 0)
	{
		std::cerr<<"No peaks found in spectrum "<<name<<" returning empty data."<<std::endl;
		return data;
	}
	else if(npeaks != 3)
	{
		std::cerr<<"Less than 3 peaks found in spectrum "<<name<<"! returning empty."<<std::endl;
		return data;
	}

	for(int i=0; i<npeaks; i++)
	{
		double x = spec.GetPositionX()[i];
		data.xvals.push_back(x);
	}

	std::sort(data.xvals.begin(), data.xvals.end(), SortGainData);

	return data;
}

//Wrapper around graph creation from std::vector data and fitting. Returns fit parameters.
CalParams GainMatcher::MakeGraph(THashTable* table, int gchan, const GraphData& data)
{
	std::string name = "channel_"+std::to_string(gchan)+"_graph";
	TGraph* graph = new TGraph(data.xvals.size(), &(data.xvals[0]), &(data.yvals[0]));
	TF1* func = new TF1("linear","pol1",0,16384);
	graph->SetName(name.c_str());
	graph->SetTitle(name.c_str());
	graph->Fit(func, "ROB|Q+");
	table->Add(graph);
	CalParams params;
	params.slope = func->GetParameter(1);
	params.intercept = func->GetParameter(0);

	return params;
}

/*
	Main loop for gain-matching all of the backs within each detector. Takes in an input data file, which should contain source calibration
	data, and two output files: a ROOT file which will contain all of the graphs and histograms, and a text file which will contain all of the
	calibration parameters. Additionally, takes in a detector channel number for both SX3s and QQQs; this channel number indicates which back channel
	will be the "fixed" channel to which all other backs are matched. Trial and error is best for chosing this.
*/
void GainMatcher::MatchBacks(const std::string& inputname, const std::string& graphname, const std::string& outputname, int sx3match, int qqqmatch)
{
	if(!cmap.IsValid() || !zmap.IsValid())
	{
		std::cerr<<"Bad map files at GainMatcher::Run! Exiting."<<std::endl;
		return;
	}

	TFile* input = TFile::Open(inputname.c_str(), "READ");
	TTree* intree = (TTree*) input->Get("EventTree");

	AnasenEvent* event = new AnasenEvent();

	intree->SetBranchAddress("event", &event);

	TFile* graphoutput = TFile::Open(graphname.c_str(), "RECREATE");

	THashTable* graph_table = new THashTable();
	THashTable* histo_table = new THashTable();

	std::ofstream output(outputname);

	std::vector<GraphData> gain_data;
	std::vector<int> match_channel; //List of channels which are to be matched against
	gain_data.resize(max_chan);
	match_channel.resize(max_chan);

	ChannelData match;
	//Create a quick lookup of global channel to match against
	for(int i=0; i<max_chan; i++)
	{
		auto channel = cmap.FindChannel(i);
		if(channel->second.detectorComponent == "FRONT" || channel->second.detectorComponent == "RING")
			continue;

		match = channel->second;
		if(match.detectorComponent == "BACK")
			match.channel = sx3match;
		else if (match.detectorComponent == "WEDGE")
			match.channel = qqqmatch;
		else
			std::cerr<<"weird channel at GainMatcher::MatchBacks()"<<std::endl;

		match_channel[i] = cmap.InverseFindChannel(match);
		if(match_channel[i] == -1)
			std::cerr<<"Found a zero to match against for GainMatcher::MatchBacks() gchan: "<<i<<" trying to match to detector channel: "<<match.channel<<std::endl;
	}

	int nentries = intree->GetEntries();
	int count=0, flush_count=0, flush_val=nentries*0.01;

	std::string name;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			count=0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*0.01*100.0<<"%"<<std::flush;
		}

		/*
			For each back (SX3 back, QQQ wedge), generate the energy spectrum from which
			peaks will be extracted
		*/
		for(int j=0; j<12; j++)
		{
			for(auto& hit : event->barrel1[j].backs)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				if(zero_offset == zmap.End())
					continue;
				MyFill(histo_table, name.c_str(), name.c_str(), 875, 1000.0, 8000.0, hit.energy - zero_offset->second);
			}
			for(auto& hit : event->barrel2[j].backs)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				if(zero_offset == zmap.End())
					continue;
				MyFill(histo_table, name.c_str(), name.c_str(), 875, 1000.0, 8000.0, hit.energy - zero_offset->second);
			}
		}

		for(int j=0; j<4; j++)
		{
			for(auto& hit : event->fqqq[j].wedges)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				if(zero_offset == zmap.End())
					continue;
				MyFill(histo_table, name.c_str(), name.c_str(), 875, 1000.0, 8000.0, hit.energy - zero_offset->second);
			}
			for(auto& hit : event->bqqq[j].wedges)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				if(zero_offset == zmap.End())
					continue;
				MyFill(histo_table, name.c_str(), name.c_str(), 875, 1000.0, 8000.0, hit.energy - zero_offset->second);
			}
		}
	}

	//Find the peaks from the energy spectra and store in an array.
	for(int i=0; i<max_chan; i++)
	{
		if(cmap.FindChannel(i)->second.detectorComponent == "FRONT" || cmap.FindChannel(i)->second.detectorComponent == "RING")
			continue;
		name = "channel_"+std::to_string(i);
		gain_data[i] = GetPoints(histo_table, name);
	}

	//Assign the data to match against, generate the graph, and obtain the fit parameters
	CalParams params;
	for(int i=0; i<max_chan; i++)
	{
		if(gain_data[i].xvals.size() == 0)
			continue;

		gain_data[i].yvals = gain_data[match_channel[i]].xvals;
		if(gain_data[i].yvals.size() == 0)
		{
			std::cout<<"Bad match condition for channel "<<i<<" matching to channel "<<match_channel[i]<<std::endl;
			continue;
		}
		params = MakeGraph(graph_table, i, gain_data[i]);
		if (i == match_channel[i])
		{
			output<<i<<"\t"<<0<<"\t"<<1<<std::endl;
			continue;
		}
		output<<i<<"\t"<<params.intercept<<"\t"<<params.slope<<std::endl;
	}

	/*
		Test the results
	*/
	ParameterMap backmap(outputname);
	if(!backmap.IsValid())
	{
		std::cerr<<"Unable to load back-gain-matching data in GainMatcher::MatchBacks()!"<<std::endl;
	}

	count=0;
	flush_count=0;
	std::cout<<"Generating test plots for back channel gain-matching..."<<std::endl;

	std::string before_name, after_name;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			count=0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*0.01*100.0<<"%"<<std::flush;
		}

		for(int j=0; j<12; j++)
		{
			before_name = "detector_barrel1_"+std::to_string(j)+"_before";
			after_name = "detector_barrel1_"+std::to_string(j)+"_after";
			for(auto& hit : event->barrel1[j].backs)
			{
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto backgains = backmap.FindParameters(hit.global_chan);
				if(zero_offset == zmap.End())
					continue;
				MyFill(histo_table, before_name.c_str(), before_name.c_str(), 875, 1000.0, 8000.0, hit.energy - zero_offset->second);
				if(backgains == backmap.End())
					continue;
				MyFill(histo_table, after_name.c_str(), after_name.c_str(), 875, 1000.0, 8000.0, backgains->second.slope*(hit.energy - zero_offset->second)+backgains->second.intercept);
			}
			before_name = "detector_barrel2_"+std::to_string(j)+"_before";
			after_name = "detector_barrel2_"+std::to_string(j)+"_after";
			for(auto& hit : event->barrel2[j].backs)
			{
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto backgains = backmap.FindParameters(hit.global_chan);
				if(zero_offset == zmap.End())
					continue;
				MyFill(histo_table, before_name.c_str(), before_name.c_str(), 875, 1000.0, 8000.0, hit.energy - zero_offset->second);
				if(backgains == backmap.End())
					continue;
				MyFill(histo_table, after_name.c_str(), after_name.c_str(), 875, 1000.0, 8000.0, backgains->second.slope*(hit.energy - zero_offset->second)+backgains->second.intercept);
			}
		}

		for(int j=0; j<4; j++)
		{
			before_name = "detector_fqqq_"+std::to_string(j)+"_before";
			after_name = "detector_fqqq_"+std::to_string(j)+"_after";
			for(auto& hit : event->fqqq[j].wedges)
			{
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto backgains = backmap.FindParameters(hit.global_chan);
				if(zero_offset == zmap.End())
					continue;
				MyFill(histo_table, before_name.c_str(), before_name.c_str(), 875, 1000.0, 8000.0, hit.energy - zero_offset->second);
				if(backgains == backmap.End())
					continue;
				MyFill(histo_table, after_name.c_str(), after_name.c_str(), 875, 1000.0, 8000.0, backgains->second.slope*(hit.energy - zero_offset->second)+backgains->second.intercept);
			}
			before_name = "detector_bqqq_"+std::to_string(j)+"_before";
			after_name = "detector_bqqq_"+std::to_string(j)+"_after";
			for(auto& hit : event->bqqq[j].wedges)
			{
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto backgains = backmap.FindParameters(hit.global_chan);
				if(zero_offset == zmap.End())
					continue;
				MyFill(histo_table, before_name.c_str(), before_name.c_str(), 875, 1000.0, 8000.0, hit.energy - zero_offset->second);
				if(backgains == backmap.End())
					continue;
				MyFill(histo_table, after_name.c_str(), after_name.c_str(), 875, 1000.0, 8000.0, backgains->second.slope*(hit.energy - zero_offset->second)+backgains->second.intercept);
			}
		}
	}
	std::cout<<std::endl;


	input->Close();
	graphoutput->cd();
	histo_table->Write();
	graph_table->Write();
	graphoutput->Close();
	output.close();
	delete event;
}

/*
	Method which gain-matches the upstream and downstream SX3 front channels. SX3 detector front strips are resistive strips -- charge is distriubted
	based on the location of the hit, which allows for positional information to be extracted by comparing the energy signals from each end of the
	strip. Gain matching is done by taking upstream and downstream data, normalizing each to the back energy, and plotting the correlation. The
	data is then fit and corrected such that it follows
		down/back = -1.0*up/back + 1.0
	The function takes in input data, which can be either source or total run data, and two output files: one is a ROOT file for storing the graphs
	and the other is a text file for storing calibration parameters. It also takes in the name of a file containg results from MatchBack, as all
	back channels need to be gain-matched prior to this analysis.
*/
void GainMatcher::MatchSX3UpDown(const std::string& inputname, const std::string& graphname, const std::string& outputname, const std::string& backmatchname)
{
	if(!cmap.IsValid() || !zmap.IsValid())
	{
		std::cerr<<"Bad map files at GainMatcher::Run! Exiting."<<std::endl;
		return;
	}

	ParameterMap backmap(backmatchname);
	if(!backmap.IsValid())
	{
		std::cerr<<"Back back-only gain-matching map at GainMatcher::MatchSX3UpDown(). Exiting."<<std::endl;
		return;
	}

	TFile* input = TFile::Open(inputname.c_str(), "READ");
	TTree* intree = (TTree*) input->Get("EventTree");

	AnasenEvent* event = new AnasenEvent();

	intree->SetBranchAddress("event", &event);

	TFile* graphoutput = TFile::Open(graphname.c_str(), "RECREATE");

	THashTable* graph_table = new THashTable();
	THashTable* histo_table = new THashTable();

	std::ofstream output(outputname);

	std::vector<GraphData> gain_data;
	gain_data.resize(max_chan);

	int nentries = intree->GetEntries();
	int count=0, flush_count=0, flush_val=nentries*0.01;


	std::string name;
	double cal_back, up_rel_energy, down_rel_energy;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			count=0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*0.01*100.0<<"%"<<std::flush;
		}

		/*
			Need to associate several chunks of data. A back and an upstream and downstream hit must all be gathered
			together. Iterating over the arrays, toss away combinations that do not meet anti-noise conditions.
			Note that there is no front-back hit assignment; rely on robust fitting to eliminate choices of back-fronts
			that don't actually come from the same hit.
		*/
		for(int j=0; j<12; j++)
		{
			if(event->barrel1[j].fronts_up.size() > 0 && event->barrel1[j].fronts_down.size() > 0 && event->barrel1[j].backs.size() > 0)
			{
				for(auto& backhit : event->barrel1[j].backs)
				{
					if(backhit.energy < 1000.0)
						continue;
					for(auto& fuphit : event->barrel1[j].fronts_up)
					{
						for(auto& fdownhit : event->barrel1[j].fronts_down)
						{
							if(fuphit.local_chan != updown_list[fdownhit.local_chan])
							{
								continue;
							}

							auto backgains = backmap.FindParameters(backhit.global_chan);
							if(backgains == backmap.End())
								continue;
							auto backzero_offset = zmap.FindOffset(backhit.global_chan);
							auto fupzero_offset = zmap.FindOffset(fuphit.global_chan);
							auto fdownzero_offset = zmap.FindOffset(fdownhit.global_chan);
							if(backzero_offset == zmap.End() || fupzero_offset == zmap.End() || fdownzero_offset == zmap.End())
								continue;
		
							cal_back = backgains->second.slope*(backhit.energy - backzero_offset->second) + backgains->second.intercept;
							up_rel_energy = (fuphit.energy - fupzero_offset->second)/(cal_back);
							down_rel_energy = (fdownhit.energy - fdownzero_offset->second)/cal_back;
							if(up_rel_energy > 1.3 || down_rel_energy > 1.3 || cal_back < 0 || up_rel_energy < 0 || down_rel_energy < 0
								|| (up_rel_energy+down_rel_energy) < 0.5 || (up_rel_energy + down_rel_energy)>1.5)
								continue;
							gain_data[fuphit.global_chan].xvals.push_back(up_rel_energy);
							gain_data[fuphit.global_chan].yvals.push_back(down_rel_energy);
						}
					}
				}
			}

			if(event->barrel2[j].fronts_up.size() > 0 && event->barrel2[j].fronts_down.size() > 0 && event->barrel2[j].backs.size() > 0)
			{
				for(auto& backhit : event->barrel2[j].backs)
				{
					if(backhit.energy < 1000.0)
						continue;
					for(auto& fuphit : event->barrel2[j].fronts_up)
					{
						for(auto& fdownhit : event->barrel2[j].fronts_down)
						{
							if(fuphit.local_chan != updown_list[fdownhit.local_chan])
							{
								continue;
							}

							auto backgains = backmap.FindParameters(backhit.global_chan);
							if(backgains == backmap.End())
								continue;
							auto backzero_offset = zmap.FindOffset(backhit.global_chan);
							auto fupzero_offset = zmap.FindOffset(fuphit.global_chan);
							auto fdownzero_offset = zmap.FindOffset(fdownhit.global_chan);
							if(backzero_offset == zmap.End() || fupzero_offset == zmap.End() || fdownzero_offset == zmap.End())
								continue;
		
							cal_back = backgains->second.slope*(backhit.energy - backzero_offset->second) + backgains->second.intercept;
							up_rel_energy = (fuphit.energy - fupzero_offset->second)/(cal_back);
							down_rel_energy = (fdownhit.energy - fdownzero_offset->second)/cal_back;
							if(up_rel_energy > 1.5 || down_rel_energy > 1.5 || cal_back < 0 || up_rel_energy < 0 || down_rel_energy < 0
								|| (up_rel_energy+down_rel_energy) < 0.5 || (up_rel_energy + down_rel_energy)>1.5)
								continue;
							gain_data[fuphit.global_chan].xvals.push_back(up_rel_energy);
							gain_data[fuphit.global_chan].yvals.push_back(down_rel_energy);
						}
					}
				}
			}
		}
	}
	std::cout<<std::endl;

	CalParams params;
	//Fit the data and write the parameters
	for(int i=0; i<max_chan; i++)
	{
		if(gain_data[i].xvals.size() < 50 || gain_data[i].yvals.size() < 50)
			continue;
		params = MakeGraph(graph_table, i, gain_data[i]);
		output<<i<<"\t"<<params.intercept<<"\t"<<params.slope<<std::endl;
	}
	output.close();

	/*
		Testing
	*/

	ParameterMap updownmap(outputname);
	if(!updownmap.IsValid())
	{
		std::cerr<<"Unable to open up-down gain-matching map at GainMatcher::MatchSX3UpDown()!"<<std::endl;
	}

	std::cout<<"Generating test plots for up-down gain-matching..."<<std::endl;
	count=0;
	flush_count=0;
	std::string before_name, after_name;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			count=0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*0.01*100.0<<"%"<<std::flush;
		}

		for(int j=0; j<12; j++)
		{
			if(event->barrel1[j].fronts_up.size() > 0 && event->barrel1[j].fronts_down.size() > 0 && event->barrel1[j].backs.size() > 0)
			{
				for(auto& backhit : event->barrel1[j].backs)
				{
					if(backhit.energy < 1000.0)
						continue;
					for(auto& fuphit : event->barrel1[j].fronts_up)
					{
						for(auto& fdownhit : event->barrel1[j].fronts_down)
						{
							if(fuphit.local_chan != updown_list[fdownhit.local_chan])
							{
								continue;
							}

							auto backgains = backmap.FindParameters(backhit.global_chan);
							if(backgains == backmap.End())
								continue;
							auto backzero_offset = zmap.FindOffset(backhit.global_chan);
							auto fupzero_offset = zmap.FindOffset(fuphit.global_chan);
							auto fdownzero_offset = zmap.FindOffset(fdownhit.global_chan);
							if(backzero_offset == zmap.End() || fupzero_offset == zmap.End() || fdownzero_offset == zmap.End())
								continue;
		
							cal_back = backgains->second.slope*(backhit.energy - backzero_offset->second) + backgains->second.intercept;
							up_rel_energy = (fuphit.energy - fupzero_offset->second)/(cal_back);
							down_rel_energy = (fdownhit.energy - fdownzero_offset->second)/cal_back;
							if(up_rel_energy > 1.5 || down_rel_energy > 1.5 || cal_back < 0 || up_rel_energy < 0 || down_rel_energy < 0)
								continue;
							before_name = "detector_barrel1_"+std::to_string(j)+"_before_channels_"+std::to_string(fuphit.local_chan)+"_"+std::to_string(fdownhit.local_chan);
							MyFill(histo_table, before_name, ";Up;Down", 1000.0, 0.0, 1.0, up_rel_energy, 1000.0, 0.0, 1.0, down_rel_energy);
							auto upgains = updownmap.FindParameters(fuphit.global_chan);
							if(upgains == updownmap.End())
								continue;
							after_name = "detector_barrel1_"+std::to_string(j)+"_after_channels_"+std::to_string(fuphit.local_chan)+"_"+std::to_string(fdownhit.local_chan);
							MyFill(histo_table, after_name, ";Up;Down", 1000.0, 0.0, 1.0, 1.0 - upgains->second.slope*up_rel_energy-upgains->second.intercept, 1000.0, 0.0, 1.0, down_rel_energy);
						}
					}
				}
			}

			if(event->barrel2[j].fronts_up.size() > 0 && event->barrel2[j].fronts_down.size() > 0 && event->barrel2[j].backs.size() > 0)
			{
				for(auto& backhit : event->barrel2[j].backs)
				{
					if(backhit.energy < 1000.0)
						continue;
					for(auto& fuphit : event->barrel2[j].fronts_up)
					{
						for(auto& fdownhit : event->barrel2[j].fronts_down)
						{
							if(fuphit.local_chan != updown_list[fdownhit.local_chan])
							{
								continue;
							}

							auto backgains = backmap.FindParameters(backhit.global_chan);
							if(backgains == backmap.End())
								continue;
							auto backzero_offset = zmap.FindOffset(backhit.global_chan);
							auto fupzero_offset = zmap.FindOffset(fuphit.global_chan);
							auto fdownzero_offset = zmap.FindOffset(fdownhit.global_chan);
							if(backzero_offset == zmap.End() || fupzero_offset == zmap.End() || fdownzero_offset == zmap.End())
								continue;
		
							cal_back = backgains->second.slope*(backhit.energy - backzero_offset->second) + backgains->second.intercept;
							up_rel_energy = (fuphit.energy - fupzero_offset->second)/(cal_back);
							down_rel_energy = (fdownhit.energy - fdownzero_offset->second)/cal_back;
							if(up_rel_energy > 1.5 || down_rel_energy > 1.5 || cal_back < 0 || up_rel_energy < 0 || down_rel_energy < 0)
								continue;
							before_name = "detector_barrel2_"+std::to_string(j)+"_before_channels_"+std::to_string(fuphit.local_chan)+"_"+std::to_string(fdownhit.local_chan);
							MyFill(histo_table, before_name, ";Up;Down", 1000.0, 0.0, 1.0, up_rel_energy, 1000.0, 0.0, 1.0, down_rel_energy);
							auto upgains = updownmap.FindParameters(fuphit.global_chan);
							if(upgains == updownmap.End())
								continue;
							after_name = "detector_barrel2_"+std::to_string(j)+"_after_channels_"+std::to_string(fuphit.local_chan)+"_"+std::to_string(fdownhit.local_chan);
							MyFill(histo_table, after_name, ";Up;Down", 1000.0, 0.0, 1.0, upgains->second.slope*up_rel_energy+upgains->second.intercept, 1000.0, 0.0, 1.0, down_rel_energy);
						}
					}
				}
			}
		}
	}
	std::cout<<std::endl;


	input->Close();
	graphoutput->cd();
	graph_table->Write();
	histo_table->Write();
	graphoutput->Close();
	delete event;
}

/*
	Method which gain-matches front channels to back channels. After matching backs to each other, and then correcting for SX3 up-down effects, can
	now gain-match all front signals to all back signals (front: SX3 sum of up-down, QQQ ring). Takes an input data file, which should be a decently
	large run to cover as much of the dynamic range as possible and two outputs: a ROOT file for graph storage and a text file for calibration
	results. Also requires a file contaning the results of MatchBack and MatchSX3UpDown as they are necessary to perform this step.
*/
void GainMatcher::MatchFrontBack(const std::string& inputname, const std::string& graphname, const std::string& outputname, const std::string& backmatchname, const std::string& updownmatchname)
{
	if(!cmap.IsValid() || !zmap.IsValid())
	{
		std::cerr<<"Bad map files at GainMatcher::Run! Exiting."<<std::endl;
		return;
	}

	ParameterMap backmap(backmatchname);
	ParameterMap updownmap(updownmatchname);
	if(!backmap.IsValid() || !updownmap.IsValid())
	{
		std::cerr<<"Bad back and up-down gain-matching files at GainMatcher::MatchFrontBack(). Exiting."<<std::endl;
		return;
	}

	TFile* input = TFile::Open(inputname.c_str(), "READ");
	TTree* intree = (TTree*) input->Get("EventTree");

	AnasenEvent* event = new AnasenEvent();

	intree->SetBranchAddress("event", &event);

	TFile* graphoutput = TFile::Open(graphname.c_str(), "RECREATE");

	THashTable* graph_table = new THashTable();
	THashTable* histo_table = new THashTable();

	std::ofstream output(outputname);

	std::vector<GraphData> gain_data;
	gain_data.resize(max_chan);

	int nentries = intree->GetEntries();
	int count=0, flush_count=0, flush_val=nentries*0.01;

	double cal_back, cal_up_energy, cal_down_energy;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			count=0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*0.01*100.0<<"%"<<std::flush;
		}

		/*
			Loop over all front-back combinations, using anti-noise conditions to reject. Note that again we do not
			make a front-back hit assignment. All valid (non-noise) combinations are made and robust fitting is used
			to reject any mismatched data.
		*/
		for(int j=0; j<12; j++)
		{
			if(event->barrel1[j].fronts_up.size() > 0 && event->barrel1[j].fronts_down.size() > 0 && event->barrel1[j].backs.size() > 0)
			{
				for(auto& backhit : event->barrel1[j].backs)
				{
					for(auto& fuphit : event->barrel1[j].fronts_up)
					{
						for(auto& fdownhit : event->barrel1[j].fronts_down)
						{
							if(fuphit.local_chan != updown_list[fdownhit.local_chan])
								continue;
							auto backgains = backmap.FindParameters(backhit.global_chan);
							if(backgains == backmap.End())
								continue;
							auto upgains = updownmap.FindParameters(fuphit.global_chan);
							if(upgains == updownmap.End())
								continue;
							auto backzero_offset = zmap.FindOffset(backhit.global_chan);
							auto fupzero_offset = zmap.FindOffset(fuphit.global_chan);
							auto fdownzero_offset = zmap.FindOffset(fdownhit.global_chan);
							if(backzero_offset == zmap.End() || fupzero_offset == zmap.End() || fdownzero_offset == zmap.End())
								continue;
		
							cal_back = backgains->second.slope*(backhit.energy - backzero_offset->second) + backgains->second.intercept;
							cal_up_energy = cal_back - upgains->second.slope*(fuphit.energy - fupzero_offset->second) - upgains->second.intercept*cal_back;
							cal_down_energy = fdownhit.energy - fdownzero_offset->second;
							if(cal_back < 100 || cal_up_energy < 100 || cal_down_energy < 100 || (cal_up_energy+cal_down_energy)/cal_back > 1.2 || (cal_up_energy+cal_down_energy)/cal_back < 0.8)
								continue;
							gain_data[fuphit.global_chan].xvals.push_back(cal_up_energy+cal_down_energy);
							gain_data[fuphit.global_chan].yvals.push_back(cal_back);
						}
					}
				}
			}

			if(event->barrel2[j].fronts_up.size() > 0 && event->barrel2[j].fronts_down.size() > 0 && event->barrel2[j].backs.size() > 0)
			{
				for(auto& backhit : event->barrel2[j].backs)
				{
					for(auto& fuphit : event->barrel2[j].fronts_up)
					{
						for(auto& fdownhit : event->barrel2[j].fronts_down)
						{
							if(fuphit.local_chan != updown_list[fdownhit.local_chan])
								continue;
							auto backgains = backmap.FindParameters(backhit.global_chan);
							if(backgains == backmap.End())
								continue;
							auto upgains = updownmap.FindParameters(fuphit.global_chan);
							if(upgains == updownmap.End())
								continue;
							auto backzero_offset = zmap.FindOffset(backhit.global_chan);
							auto fupzero_offset = zmap.FindOffset(fuphit.global_chan);
							auto fdownzero_offset = zmap.FindOffset(fdownhit.global_chan);
							if(backzero_offset == zmap.End() || fupzero_offset == zmap.End() || fdownzero_offset == zmap.End())
								continue;
		
							cal_back = backgains->second.slope*(backhit.energy - backzero_offset->second) + backgains->second.intercept;
							cal_up_energy = cal_back - upgains->second.slope*(fuphit.energy - fupzero_offset->second) - upgains->second.intercept*cal_back;
							cal_down_energy = fdownhit.energy - fdownzero_offset->second;
							if(cal_back < 100 || cal_up_energy < 100 || cal_down_energy < 100 || (cal_up_energy+cal_down_energy)/cal_back > 1.2 || (cal_up_energy+cal_down_energy)/cal_back < 0.8)
								continue;
							gain_data[fuphit.global_chan].xvals.push_back(cal_up_energy+cal_down_energy);
							gain_data[fuphit.global_chan].yvals.push_back(cal_back);
						}
					}
				}
			}
		}

		for(int j=0; j<4; j++)
		{
			if(event->fqqq[j].rings.size() > 0 && event->fqqq[j].wedges.size() > 0)
			{
				for(auto& wedgehit : event->fqqq[j].wedges)
				{
					for(auto& ringhit : event->fqqq[j].rings)
					{
						auto wedgegains = backmap.FindParameters(wedgehit.global_chan);
						if(wedgegains == backmap.End())
							continue;
						auto wedgezero_offset = zmap.FindOffset(wedgehit.global_chan);
						auto ringzero_offset = zmap.FindOffset(ringhit.global_chan);
						if(wedgezero_offset == zmap.End() || ringzero_offset == zmap.End())
							continue;
	
						cal_back = wedgegains->second.slope*(wedgehit.energy - wedgezero_offset->second) + wedgegains->second.intercept;
						cal_up_energy = ringhit.energy - ringzero_offset->second;
						if(cal_back < 0 || cal_up_energy < 0 || cal_up_energy/cal_back > 1.2 || cal_up_energy/cal_back < 0.8)
							continue;
						gain_data[ringhit.global_chan].xvals.push_back(cal_up_energy);
						gain_data[ringhit.global_chan].yvals.push_back(cal_back);
					}
				}
			}
			if(event->bqqq[j].rings.size() > 0 && event->bqqq[j].wedges.size() > 0)
			{
				for(auto& wedgehit : event->bqqq[j].wedges)
				{
					for(auto& ringhit : event->bqqq[j].rings)
					{
						auto wedgegains = backmap.FindParameters(wedgehit.global_chan);
						if(wedgegains == backmap.End())
							continue;
						auto wedgezero_offset = zmap.FindOffset(wedgehit.global_chan);
						auto ringzero_offset = zmap.FindOffset(ringhit.global_chan);
						if(wedgezero_offset == zmap.End() || ringzero_offset == zmap.End())
							continue;
	
						cal_back = wedgegains->second.slope*(wedgehit.energy - wedgezero_offset->second) + wedgegains->second.intercept;
						cal_up_energy = ringhit.energy - ringzero_offset->second;
						if(cal_back < 0 || cal_up_energy < 0 || cal_up_energy/cal_back > 1.2 || cal_up_energy/cal_back < 0.8)
							continue;
						gain_data[ringhit.global_chan].xvals.push_back(cal_up_energy);
						gain_data[ringhit.global_chan].yvals.push_back(cal_back);
					}
				}
			}
		}
	}
	std::cout<<std::endl;

	//Generate graphs, obtain and write fit data
	CalParams params;
	for(int i=0; i<max_chan; i++)
	{
		if(gain_data[i].xvals.size() < 10|| gain_data[i].yvals.size() < 10)
			continue;
		params = MakeGraph(graph_table, i, gain_data[i]);
		output<<i<<"\t"<<params.intercept<<"\t"<<params.slope<<std::endl;
	}
	output.close();

	/*
		Testing
	*/
	ParameterMap frontbackmap(outputname);
	if(!frontbackmap.IsValid())
	{
		std::cerr<<"Unable to open front-back gain-matching map at GainMatcher::MatchFrontBack()!"<<std::endl;
	}

	count=0;
	flush_count=0;
	std::string before_name, after_name;
	std::cout<<"Generating front-back gain-matching test plots..."<<std::endl;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			count=0;
			flush_count++;
			std::cout<<"\rPercent of data processed: "<<flush_count*0.01*100.0<<"%"<<std::flush;
		}

		for(int j=0; j<12; j++)
		{
			if(event->barrel1[j].fronts_up.size() > 0 && event->barrel1[j].fronts_down.size() > 0 && event->barrel1[j].backs.size() > 0)
			{
				for(auto& backhit : event->barrel1[j].backs)
				{
					for(auto& fuphit : event->barrel1[j].fronts_up)
					{
						for(auto& fdownhit : event->barrel1[j].fronts_down)
						{
							if(fuphit.local_chan != updown_list[fdownhit.local_chan])
								continue;
							auto backgains = backmap.FindParameters(backhit.global_chan);
							if(backgains == backmap.End())
								continue;
							auto upgains = updownmap.FindParameters(fuphit.global_chan);
							if(upgains == updownmap.End())
								continue;
							auto backzero_offset = zmap.FindOffset(backhit.global_chan);
							auto fupzero_offset = zmap.FindOffset(fuphit.global_chan);
							auto fdownzero_offset = zmap.FindOffset(fdownhit.global_chan);
							if(backzero_offset == zmap.End() || fupzero_offset == zmap.End() || fdownzero_offset == zmap.End())
								continue;
		
							cal_back = backgains->second.slope*(backhit.energy - backzero_offset->second) + backgains->second.intercept;
							cal_up_energy = cal_back - upgains->second.slope*(fuphit.energy - fupzero_offset->second) - upgains->second.intercept*cal_back;
							cal_down_energy = fdownhit.energy - fdownzero_offset->second;
							if(cal_back < 100 || cal_up_energy < 100 || cal_down_energy < 100 || (cal_up_energy+cal_down_energy)/cal_back > 1.2 || (cal_up_energy+cal_down_energy)/cal_back < 0.8)
								continue;
							before_name = "channel_"+std::to_string(fuphit.global_chan)+"_before";
							MyFill(histo_table, before_name,";Front;Back",1024,0.0,16384,cal_up_energy+cal_down_energy,1024,0,16384,cal_back);
							auto frontbackgains = frontbackmap.FindParameters(fuphit.global_chan);
							if(frontbackgains == frontbackmap.End())
								continue;
							after_name = "channel_"+std::to_string(fuphit.global_chan)+"_after";
							MyFill(histo_table, after_name, ";Front;Back",1024,0,16384,frontbackgains->second.slope*(cal_up_energy+cal_down_energy)+frontbackgains->second.intercept,1024,0,16384,cal_back);
						}
					}
				}
			}

			if(event->barrel2[j].fronts_up.size() > 0 && event->barrel2[j].fronts_down.size() > 0 && event->barrel2[j].backs.size() > 0)
			{
				for(auto& backhit : event->barrel2[j].backs)
				{
					for(auto& fuphit : event->barrel2[j].fronts_up)
					{
						for(auto& fdownhit : event->barrel2[j].fronts_down)
						{
							if(fuphit.local_chan != updown_list[fdownhit.local_chan])
								continue;
							auto backgains = backmap.FindParameters(backhit.global_chan);
							if(backgains == backmap.End())
								continue;
							auto upgains = updownmap.FindParameters(fuphit.global_chan);
							if(upgains == updownmap.End())
								continue;
							auto backzero_offset = zmap.FindOffset(backhit.global_chan);
							auto fupzero_offset = zmap.FindOffset(fuphit.global_chan);
							auto fdownzero_offset = zmap.FindOffset(fdownhit.global_chan);
							if(backzero_offset == zmap.End() || fupzero_offset == zmap.End() || fdownzero_offset == zmap.End())
								continue;
		
							cal_back = backgains->second.slope*(backhit.energy - backzero_offset->second) + backgains->second.intercept;
							cal_up_energy = cal_back - upgains->second.slope*(fuphit.energy - fupzero_offset->second) - upgains->second.intercept*cal_back;
							cal_down_energy = fdownhit.energy - fdownzero_offset->second;
							if(cal_back < 100 || cal_up_energy < 100 || cal_down_energy < 100 || (cal_up_energy+cal_down_energy)/cal_back > 1.2 || (cal_up_energy+cal_down_energy)/cal_back < 0.8)
								continue;
							before_name = "channel_"+std::to_string(fuphit.global_chan)+"_before";
							MyFill(histo_table, before_name,";Front;Back",1024,0.0,16384,cal_up_energy+cal_down_energy,1024,0,16384,cal_back);
							auto frontbackgains = frontbackmap.FindParameters(fuphit.global_chan);
							if(frontbackgains == frontbackmap.End())
								continue;
							after_name = "channel_"+std::to_string(fuphit.global_chan)+"_after";
							MyFill(histo_table, after_name, ";Front;Back",1024,0,16384,frontbackgains->second.slope*(cal_up_energy+cal_down_energy)+frontbackgains->second.intercept,1024,0,16384,cal_back);
						}
					}
				}
			}
		}

		for(int j=0; j<4; j++)
		{
			if(event->fqqq[j].rings.size() > 0 && event->fqqq[j].wedges.size() > 0)
			{
				for(auto& wedgehit : event->fqqq[j].wedges)
				{
					for(auto& ringhit : event->fqqq[j].rings)
					{
						auto wedgegains = backmap.FindParameters(wedgehit.global_chan);
						if(wedgegains == backmap.End())
							continue;
						auto wedgezero_offset = zmap.FindOffset(wedgehit.global_chan);
						auto ringzero_offset = zmap.FindOffset(ringhit.global_chan);
						if(wedgezero_offset == zmap.End() || ringzero_offset == zmap.End())
							continue;
	
						cal_back = wedgegains->second.slope*(wedgehit.energy - wedgezero_offset->second) + wedgegains->second.intercept;
						cal_up_energy = ringhit.energy - ringzero_offset->second;
						if(cal_back < 0 || cal_up_energy < 0 || cal_up_energy/cal_back > 1.2 || cal_up_energy/cal_back < 0.8)
							continue;
						before_name = "channel_"+std::to_string(ringhit.global_chan)+"_before";
						MyFill(histo_table, before_name,";Front;Back",1024,0.0,16384,cal_up_energy,1024,0,16384,cal_back);
						auto frontbackgains = frontbackmap.FindParameters(ringhit.global_chan);
						if(frontbackgains == frontbackmap.End())
							continue;
						after_name = "channel_"+std::to_string(ringhit.global_chan)+"_after";
						MyFill(histo_table, after_name, ";Front;Back",1024,0,16384,frontbackgains->second.slope*cal_up_energy+frontbackgains->second.intercept,1024,0,16384,cal_back);
					}
				}
			}
			if(event->bqqq[j].rings.size() > 0 && event->bqqq[j].wedges.size() > 0)
			{
				for(auto& wedgehit : event->bqqq[j].wedges)
				{
					for(auto& ringhit : event->bqqq[j].rings)
					{
						auto wedgegains = backmap.FindParameters(wedgehit.global_chan);
						if(wedgegains == backmap.End())
							continue;
						auto wedgezero_offset = zmap.FindOffset(wedgehit.global_chan);
						auto ringzero_offset = zmap.FindOffset(ringhit.global_chan);
						if(wedgezero_offset == zmap.End() || ringzero_offset == zmap.End())
							continue;
	
						cal_back = wedgegains->second.slope*(wedgehit.energy - wedgezero_offset->second) + wedgegains->second.intercept;
						cal_up_energy = ringhit.energy - ringzero_offset->second;
						if(cal_back < 0 || cal_up_energy < 0 || cal_up_energy/cal_back > 1.2 || cal_up_energy/cal_back < 0.8)
							continue;
						before_name = "channel_"+std::to_string(ringhit.global_chan)+"_before";
						MyFill(histo_table, before_name,";Front;Back",1024,0.0,16384,cal_up_energy,1024,0,16384,cal_back);
						auto frontbackgains = frontbackmap.FindParameters(ringhit.global_chan);
						if(frontbackgains == frontbackmap.End())
							continue;
						after_name = "channel_"+std::to_string(ringhit.global_chan)+"_after";
						MyFill(histo_table, after_name, ";Front;Back",1024,0,16384,frontbackgains->second.slope*cal_up_energy+frontbackgains->second.intercept,1024,0,16384,cal_back);
					}
				}
			}
		}
	}
	std::cout<<std::endl;

	input->Close();
	graphoutput->cd();
	graph_table->Write();
	histo_table->Write();
	graphoutput->Close();
	delete event;
}