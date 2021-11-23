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
	std::string outname = "/media/gordon/GordonData/gwm17/7BeNov2021/merged/runs-2448-2477_130keV.root";
	TFile* output = TFile::Open(outname.c_str(), "RECREATE");

	std::string searchdir = "/media/gordon/GordonData/gwm17/7BeNov2021/runs/run-";

	std::vector<int> runs;
	for(int i=2448; i<2478; i++)
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
