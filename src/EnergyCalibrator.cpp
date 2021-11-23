/*
	EnergyCalibrator
	Class designed to perform energy calibrations. Designed for use with gain-matched results and source data.
	Main method is Run.

	Written by Gordon McCann Nov 2021
*/

#include "EnergyCalibrator.h"
#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH1.h>
#include <TF1.h>
#include <algorithm>

//For use with std::sort
bool SortEnergyData(const double i, const double j) 
{
	return i<j;
}

/*
	Currently requires entire suite of gain-matching results. In principle this could be relaxed... but
	probably shouldn't be.
	Note sigma, threshold: parameters passed to TSpectrum for peak searching, sigma is width and threshold is fraction
	of maximum peak height. These may need adjusted for each experiment.
*/
EnergyCalibrator::EnergyCalibrator(const std::string& channelfile, const std::string& zerofile, const std::string& backmatch, const std::string& updownmatch, const std::string& frontbackmatch) :
	cmap(channelfile), zmap(zerofile), bmap(backmatch), udmap(updownmatch), fbmap(frontbackmatch), sigma(1.0), threshold(0.4)
{
}

EnergyCalibrator::~EnergyCalibrator() {}

//Wrapper on histogram creation, storage, and filling.
void EnergyCalibrator::FillHistogram(THashTable* table, const std::string& name, const std::string& title, int bins, double minx, double maxx, double value)
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

/*
	Function which calls TSpectrum to obtain peak locations from histograms.
	Peak locations are then returned as x-coordinates of GraphData, with associated
	peak energy values as y-coordinates.
*/
GraphData EnergyCalibrator::GetPoints(THashTable* table, int gchan, const std::string& name)
{
	GraphData data;

	TH1* histo = (TH1*) table->FindObject(name.c_str());
	if(histo == nullptr)
	{
		//std::cerr<<"Histogram named "<<name<<" not found at EnergyCalibrator::GetPoints! Returning empty data."<<std::endl;
		return data;
	}

	int npeaks = spec.Search(histo, sigma, "nobackground", threshold);
	int nepeaks = energyValues.size();
	if(npeaks != nepeaks)
	{
		std::cerr<<"Npeaks not equal to "<<nepeaks<<" in spectrum "<<name<<" returning empty data."<<std::endl;
		return data;
	}

	for(int i=0; i<npeaks; i++)
	{
		double x = spec.GetPositionX()[i];
		data.xvals.push_back(x);
	}

	std::sort(data.xvals.begin(), data.xvals.end(), SortEnergyData);

	data.yvals = energyValues;

	return data;
}

/*
	Method which generates a graph, fits it, and returns the fit parameters.
*/
CalParams EnergyCalibrator::CalibrateEnergy(THashTable* table, const std::string& name, const GraphData& data) 
{
	CalParams parameters;
	std::string linename = name+"_fit";
	TGraph* graph = new TGraph(data.xvals.size(), &(data.xvals[0]), &(data.yvals[0]));
	graph->SetTitle(name.c_str());
	graph->SetName(name.c_str());
	TF1* func = new TF1(linename.c_str(),"pol1",0.0,16384.0);
	graph->Fit(func, "Q|ROB+");

	parameters.intercept = func->GetParameter(0);
	parameters.slope = func->GetParameter(1);

	table->Add(graph);
	return parameters;
}

/*
	Main loop. Takes in input data file, which should contain source calibration data, and two output files: one which is 
	a ROOT file for storing graphs, and a text file for storing calibraton parameters.
*/
void EnergyCalibrator::Run(const std::string& inputname, const std::string& plotname, const std::string& outputname) 
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
	int count=0, flush_count=0, flush_val = 0.01*nentries;

	std::string name;
	double cal_energy;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			flush_count++;
			count=0;
			std::cout<<"\rPercent of data histogrammed: "<<flush_count*0.01*100.0<<"%"<<std::flush;
		}

		/*
			Generate a calibration spectrum for each channel, excluding SX3 fronts (only used for positional data)
		*/
		for(int j=0; j<12; j++)
		{
			for(auto& hit : event->barrel1[j].backs)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = bmap.FindParameters(hit.global_chan);
				if(gains == bmap.End() || zero_offset == zmap.End())
					continue;
				cal_energy = gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept;
				FillHistogram(histo_table, name, name, 925,600.0,8000.0, cal_energy);
			}

			for(auto& hit : event->barrel2[j].backs)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = bmap.FindParameters(hit.global_chan);
				if(gains == bmap.End() || zero_offset == zmap.End())
					continue;
				cal_energy = gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept;
				FillHistogram(histo_table, name, name, 925,600.0,8000.0, cal_energy);
			}
		}

		for(int j=0; j<4; j++)
		{
			for(auto& hit : event->fqqq[j].rings)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = fbmap.FindParameters(hit.global_chan);
				if(gains == fbmap.End() || zero_offset == zmap.End())
					continue;
				cal_energy = gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept;
				FillHistogram(histo_table, name, name, 925,600.0,8000.0, cal_energy);
			}
			for(auto& hit : event->fqqq[j].wedges)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = bmap.FindParameters(hit.global_chan);
				if(gains == bmap.End() || zero_offset == zmap.End())
					continue;
				cal_energy = gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept;
				FillHistogram(histo_table, name, name, 925,600.0,8000.0, cal_energy);
			}
			for(auto& hit : event->bqqq[j].rings)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = fbmap.FindParameters(hit.global_chan);
				if(gains == fbmap.End() || zero_offset == zmap.End())
					continue;
				cal_energy = gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept;
				FillHistogram(histo_table, name, name, 925,600.0,8000.0, cal_energy);
			}
			for(auto& hit : event->bqqq[j].wedges)
			{
				name = "channel_"+std::to_string(hit.global_chan);
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = bmap.FindParameters(hit.global_chan);
				if(gains == bmap.End() || zero_offset == zmap.End())
					continue;
				cal_energy = gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept;
				FillHistogram(histo_table, name, name, 925,600.0,8000.0, cal_energy);
			}
		}
	}
	std::cout<<std::endl;

	//Generate graphs, obtain fit parameters
	GraphData data;
	CalParams parameters;
	for(int i=0; i<nchannels; i++)
	{
		name = "channel_"+std::to_string(i);
		data = GetPoints(histo_table, i, name);

		if(data.xvals.size() ==  0)
			continue;

		name += "_graph";
		parameters = CalibrateEnergy(graph_table, name, data);

		output<<i<<"\t"<<parameters.intercept<<"\t"<<parameters.slope<<std::endl;
	}
	output.close();

	/*
		Testing
	*/
	ParameterMap energymap(outputname);
	if(!energymap.IsValid())
	{
		std::cerr<<"Energy calibration map is not valid at EnergyCalibrator::Run()!"<<std::endl;
	}

	std::cout<<"Generating energy calibration test plots..."<<std::endl;
	count=0;
	flush_count=0;
	for(int i=0; i<nentries; i++)
	{
		intree->GetEntry(i);
		count++;
		if(count == flush_val)
		{
			flush_count++;
			count=0;
			std::cout<<"\rPercent of data histogrammed: "<<flush_count*0.01*100.0<<"%"<<std::flush;
		}

		for(int j=0; j<12; j++)
		{
			for(auto& hit : event->barrel1[j].backs)
			{
				name = "channel_"+std::to_string(hit.global_chan)+"_calibrated";
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = bmap.FindParameters(hit.global_chan);
				auto ecal = energymap.FindParameters(hit.global_chan);
				if(gains == bmap.End() || zero_offset == zmap.End() || ecal == energymap.End())
					continue;
				cal_energy = ecal->second.slope*(gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept) + ecal->second.intercept;
				FillHistogram(histo_table, name, name, 1000.0,0.0,10.0, cal_energy);
			}

			for(auto& hit : event->barrel2[j].backs)
			{
				name = "channel_"+std::to_string(hit.global_chan)+"_calibrated";
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = bmap.FindParameters(hit.global_chan);
				auto ecal = energymap.FindParameters(hit.global_chan);
				if(gains == bmap.End() || zero_offset == zmap.End() || ecal == energymap.End())
					continue;
				cal_energy = ecal->second.slope*(gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept) + ecal->second.intercept;
				FillHistogram(histo_table, name, name, 1000.0,0.0,10.0, cal_energy);
			}
		}

		for(int j=0; j<4; j++)
		{
			for(auto& hit : event->fqqq[j].rings)
			{
				name = "channel_"+std::to_string(hit.global_chan)+"_calibrated";
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = fbmap.FindParameters(hit.global_chan);
				auto ecal = energymap.FindParameters(hit.global_chan);
				if(gains == fbmap.End() || zero_offset == zmap.End() || ecal == energymap.End())
					continue;
				cal_energy = ecal->second.slope*(gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept) + ecal->second.intercept;
				FillHistogram(histo_table, name, name, 1000.0,0.0,10.0, cal_energy);
			}
			for(auto& hit : event->fqqq[j].wedges)
			{
				name = "channel_"+std::to_string(hit.global_chan)+"_calibrated";
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = bmap.FindParameters(hit.global_chan);
				auto ecal = energymap.FindParameters(hit.global_chan);
				if(gains == bmap.End() || zero_offset == zmap.End() || ecal == energymap.End())
					continue;
				cal_energy = ecal->second.slope*(gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept) + ecal->second.intercept;
				FillHistogram(histo_table, name, name, 1000.0,0.0,10.0, cal_energy);
			}
			for(auto& hit : event->bqqq[j].rings)
			{
				name = "channel_"+std::to_string(hit.global_chan)+"_calibrated";
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = fbmap.FindParameters(hit.global_chan);
				auto ecal = energymap.FindParameters(hit.global_chan);
				if(gains == fbmap.End() || zero_offset == zmap.End() || ecal == energymap.End())
					continue;
				cal_energy = ecal->second.slope*(gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept) + ecal->second.intercept;
				FillHistogram(histo_table, name, name, 1000.0,0.0,10.0, cal_energy);
			}
			for(auto& hit : event->bqqq[j].wedges)
			{
				name = "channel_"+std::to_string(hit.global_chan)+"_calibrated";
				auto zero_offset = zmap.FindOffset(hit.global_chan);
				auto gains = bmap.FindParameters(hit.global_chan);
				auto ecal = energymap.FindParameters(hit.global_chan);
				if(gains == bmap.End() || zero_offset == zmap.End() || ecal == energymap.End())
					continue;
				cal_energy = ecal->second.slope*(gains->second.slope*(hit.energy - zero_offset->second) + gains->second.intercept) + ecal->second.intercept;
				FillHistogram(histo_table, name, name, 1000.0,0.0,10.0, cal_energy);
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

