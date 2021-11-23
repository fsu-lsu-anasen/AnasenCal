/*
	ZeroCalibrator
	Class for generating zero-offset calibrations for ANASEN silicon detectors. ASICS do not fix the zero point on the ADC scale,
	therefore it is necessary to determine the zero point after the experiment. Pulser data is fed in, and TSpectrum is used to
	identify the pulser peaks. The offset is then determined from a linear fit.

	Written by Gordon McCann Nov 2021
*/
#include "ZeroCalibrator.h"
#include "ZeroCalMap.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>

//for use with std::sort
bool SortZeroData(const double i, const double j) 
{
	return i<j;
}


/*
	Sigma is the width TSpetrum uses and threshold is the percentage less than the max peak height used as a cutoff by TSpectrum.
	These may need to be adjusted on an experiment by experiment basis.
*/
ZeroCalibrator::ZeroCalibrator(const std::string& channelfile) :
	sigma(25.0), threshold(0.15), cmap(channelfile)
{
}

ZeroCalibrator::~ZeroCalibrator() {}

//Wrappers which make, store, and fill histograms
void ZeroCalibrator::FillHistogram(THashTable* table, const std::string& name, const std::string& title, int bins, double minx, double maxx, double value)
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

void ZeroCalibrator::FillHistogram(THashTable* table, const std::string& name, const std::string& title, int binsx, double minx, double maxx, double valuex,
																										int binsy, double miny, double maxy, double valuey)
{
	TH1* histo = (TH1*) table->FindObject(name.c_str());
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
	Method which calls TSpectrum to find peaks; retrieves the histogram from the table, checks channel info to deterime
	how many peaks are expected 
*/
GraphData ZeroCalibrator::GetPoints(THashTable* table, int gchan, const std::string& name)
{
	GraphData data;

	TH1* histo = (TH1*) table->FindObject(name.c_str());
	if(histo == nullptr)
	{
		//std::cerr<<"Histogram named "<<name<<" not found at ZeroCalibrator::GetPoints! Returning empty data."<<std::endl;
		return data;
	} 
	else if(histo->Integral() < 1000.0)
	{
		//std::cerr<<"Histogram "<<name<<" with less than one thousand total counts found. Rejected."<<std::endl;
		return data;
	}

	auto channel = cmap.FindChannel(gchan);
	if(channel->second.detectorComponent == "RING" || channel->second.detectorComponent ==  "FRONT")
		threshold = 0.025;
	else if(channel->second.detectorComponent == "BACK" || channel->second.detectorComponent == "WEDGE")
		threshold = 0.15;

	int npeaks = spec.Search(histo, sigma, "nobackground", threshold);
	int nfrontpeaks = frontPulseValues.size();
	int nbackpeaks = backPulseValues.size();
	if(npeaks == 0)
	{
		std::cerr<<"No peaks found in spectrum "<<name<<" returning empty data."<<std::endl;
		return data;
	}
	else if ((channel->second.detectorComponent == "FRONT" || channel->second.detectorComponent == "RING") && npeaks != nfrontpeaks && (gchan < 224 || gchan > 287))
	{
		std::cerr<<"At ZeroCalibrator::GetPoints front/ring channel (gchan="<<gchan<<") does not have "<<nfrontpeaks<<" peaks (found "<<npeaks<<"). Returning empty data."<<std::endl;
		return data;
	}
	else if(gchan>223 && gchan <287 && (npeaks != nfrontpeaks-1 && npeaks != nfrontpeaks))
	{
		std::cerr<<"At ZeroCalibrator::GetPoints front/ring channel (gchan="<<gchan<<") does not have "<<nfrontpeaks-1<<" peaks (found "<<npeaks<<"). Returning empty data."<<std::endl;
		return data;
	}
	else if ((channel->second.detectorComponent == "BACK" || channel->second.detectorComponent == "WEDGE") && npeaks != nbackpeaks)
	{
		std::cerr<<"At ZeroCalibrator::GetPoints back/wedge channel (gchan="<<gchan<<") does not have "<<nbackpeaks<<" peaks (found "<<npeaks<<"). Returning empty data."<<std::endl;
		return data;
	}

	for(int i=0; i<npeaks; i++)
	{
		double x = spec.GetPositionX()[i];
		data.xvals.push_back(x);
	}

	

	std::sort(data.xvals.begin(), data.xvals.end(), SortZeroData);

	if((channel->second.detectorComponent == "FRONT" || channel->second.detectorComponent == "RING"))
	{
		data.yvals = frontPulseValues;
	}
	else if ((channel->second.detectorComponent == "BACK" || channel->second.detectorComponent == "WEDGE"))
	{
		data.yvals = backPulseValues;
	}

	return data;
}

//Wrapper which generates graphs from vectors, fits, and reports fit parameters
double ZeroCalibrator::MakeGraph(THashTable* table, const std::string& name, const GraphData& data)
{
	TGraph* graph = new TGraph(data.xvals.size(), &(data.xvals[0]), &(data.yvals[0]));
	graph->SetName(name.c_str());

	TF1* func = new TF1("linear", "pol1", 0, 16000);
	graph->Fit(func, "Q|ROB+");

	double offset = func->GetParameter(0);
	double slope = func->GetParameter(1);

	table->Add(graph);

	return -offset/slope; //convert from y-intercept to x-intercept
}

/*
	Main loop, takes in an input file and two output files: one output is a ROOT file for plots, the other a textfile
	for the calibration results
*/
void ZeroCalibrator::Run(const std::string& inputname, const std::string& plotname, const std::string& outputname)
{

	TFile* input = TFile::Open(inputname.c_str(), "READ");
	if(!input->IsOpen())
	{
		std::cerr<<"Unable to open input datafile "<<inputname<<"! Quitting."<<std::endl;
		return;
	}
	TTree* intree = (TTree*) input->Get("EventTree");

	AnasenEvent* event = new AnasenEvent();

	intree->SetBranchAddress("event", &event);

	TFile* graphoutput = TFile::Open(plotname.c_str(), "RECREATE");
	if(!graphoutput->IsOpen())
	{
		input->Close();
		std::cerr<<"Unable to create output graph file "<<plotname<<"! Quitting."<<std::endl;
		return;
	}

	THashTable* histo_table = new THashTable();
	THashTable* graph_table = new THashTable();

	std::ofstream output(outputname);
	if(!output.is_open())
	{
		input->Close();
		graphoutput->Close();
		std::cerr<<"Unable to open output file "<<outputname<<". Quitting."<<std::endl;
		return;
	}

	int nentries = intree->GetEntries();
	int count=0, flush_count=0, flush_val = 0.05*nentries;

	std::string name;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			flush_count++;
			count=0;
			std::cout<<"\rPercent of data histogrammed: "<<flush_count*5<<"%"<<std::flush;
		}

		/*
			In the case of zero-offset calibrations, each channel should be calibrated
			independently.
		*/
		for(int j=0; j<12; j++)
		{
			for(auto& hit : event->barrel1[j].fronts_up)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}
			for(auto& hit : event->barrel1[j].fronts_down)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}
			for(auto& hit : event->barrel1[j].backs)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}

			for(auto& hit : event->barrel2[j].fronts_up)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}
			for(auto& hit : event->barrel2[j].fronts_down)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}
			for(auto& hit : event->barrel2[j].backs)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}
		}

		for(int j=0; j<4; j++)
		{
			for(auto& hit : event->fqqq[j].rings)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}
			for(auto& hit : event->fqqq[j].wedges)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}

			for(auto& hit : event->bqqq[j].rings)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}
			for(auto& hit : event->bqqq[j].wedges)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				FillHistogram(histo_table, name, name, 16384,1400,16384, hit.energy);
			}
		}
	}

	

	int nchannels = 544;
	GraphData data;
	double offset;
	for(int i=0; i<nchannels; i++)
	{
		name = "channel_"+std::to_string(i);
		data = GetPoints(histo_table, i, name);

		if(data.xvals.size() ==  0)
			continue;

		name += "_graph";
		offset = MakeGraph(graph_table, name, data);

		output<<i<<"\t"<<offset<<std::endl;
	}
	output.close();

	ZeroCalMap zmap(outputname);
	if(!zmap.IsValid())
	{
		std::cerr<<"Unable to open a map after creating calibrations in ZeroCalibrator::Run()."<<std::endl;
	}

	count=0;
	flush_count=0;
	std::string before_name="before_offset_cal";
	std::string before_title="before_offset_cal;Channel;Energy(arb)";
	std::string after_name="after_offset_cal";
	std::string after_title="after_offset_cal;Channel;Energy(arb)";

	std::cout<<"Generating zero-offset calibration test plots..."<<std::endl;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			flush_count++;
			count=0;
			std::cout<<"\rPercent of data histogrammed: "<<flush_count*5<<"%"<<std::flush;
		}

		for(int j=0; j<12; j++)
		{
			for(auto& hit : event->barrel1[j].fronts_up)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}
			for(auto& hit : event->barrel1[j].fronts_down)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}
			for(auto& hit : event->barrel1[j].backs)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}

			for(auto& hit : event->barrel2[j].fronts_up)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}
			for(auto& hit : event->barrel2[j].fronts_down)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}
			for(auto& hit : event->barrel2[j].backs)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}
		}

		for(int j=0; j<4; j++)
		{
			for(auto& hit : event->fqqq[j].rings)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}
			for(auto& hit : event->fqqq[j].wedges)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}

			for(auto& hit : event->bqqq[j].rings)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}
			for(auto& hit : event->bqqq[j].wedges)
			{
				FillHistogram(histo_table, before_name, before_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy);
				auto offset_data = zmap.FindOffset(hit.global_chan);
				if(offset_data == zmap.End())
					continue;
				FillHistogram(histo_table, after_name, after_title, nchannels, 0, nchannels, hit.global_chan, 16384,1400,16384, hit.energy-offset_data->second);
			}
		}
	}
	std::cout<<std::endl;

	input->Close();

	graphoutput->cd();
	histo_table->Write();
	graph_table->Write();
	graphoutput->Close();
	delete event;
}