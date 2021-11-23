#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <iostream>
#include "../include/DataStructs.h"

R__LOAD_LIBRARY(../objs/libAnasenEvent_dict.so)

void testPlot()
{
	TFile* file = TFile::Open("/media/gordon/GordonData/gwm17/7BeNov2021/merged/pulserData_bothSides.root");
	TTree* tree = (TTree*) file->Get("EventTree");

	AnasenEvent* event = new AnasenEvent();

	tree->SetBranchAddress("event", &event);

	TH1F* histo = new TH1F("fqqq0_0_energy", "fqqq0_0_energy",16384,0, 16384);

	for(int i=0; i<tree->GetEntries(); i++)
	{
		tree->GetEntry(i);

		if(event->fqqq[0].rings.size() != 0)
		{
				histo->Fill(event->fqqq[0].rings[0].energy);
		}
	}

	histo->Draw();
}
