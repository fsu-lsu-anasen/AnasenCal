#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <vector>
#include <iostream>
#include "../include/DataStructs.h"

R__LOAD_LIBRARY(../objs/libAnasenEvent_dict.so)

void chainFiles()
{
	int start = 2448;
	int stop = 2477;
	std::string mod = "_130keV_timed";
	std::string outname = "/data1/gwm17/7BeNov2021/merged/runs-"+std::to_string(start)+"-"+std::to_string(stop)+mod+".root";
	TFile* output = TFile::Open(outname.c_str(), "RECREATE");

	std::string searchdir = "/data1/gwm17/7BeNov2021/runs/run-";

	std::vector<int> runs;
	for(int i=start; i<(stop+1); i++)
		runs.push_back(i); 

	TChain* chainer = new TChain("EventTree", "EventTree");

	std::string filename;
	std::cout<<"Merging ROOT files..."<<std::endl;
	for(auto run : runs)
	{
		filename = searchdir + std::to_string(run) + ".root";
		chainer->Add(filename.c_str());
		std::cout<<"Added file "<<filename<<std::endl;
	}

	std::cout<<"Performing merge to "<<outname<<"..."<<std::endl;
	chainer->Merge(output, 0, "fast");
	std::cout<<"Finished."<<std::endl;
	
}
