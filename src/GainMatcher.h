/*
	GainMatcher
	Class oriented around performing silicon detector gain matching in ANASEN. Two types of detectors: SX3's and QQQ's.
	There are three gain matching methods: MatchBacks, MatchSX3UpDown, and MatchFrontBack. They should be performed in that order,
	as each step relies upon the previous results.

	Written by Gordon McCann Nov 2021

	Modified for use with the Barcelona-ANASEN scheme
*/
#ifndef GAINMATCHER_H
#define GAINMATCHER_H

#include <string>
#include <vector>
#include <THashTable.h>
#include "ChannelMap.h"
#include "DataStructs.h"
#include "TSpectrum.h"

class GainMatcher
{
public:
	GainMatcher(const std::string& channelfile);
	~GainMatcher();
	void MatchBacks(const std::string& inputname, const std::string& plotname, const std::string& outputname, int sx3match, int qqqmatch);
	void MatchSX3UpDown(const std::string& inputname, const std::string& plotname, const std::string& outputname, const std::string& backmatchname);
	void MatchFrontBack(const std::string& inputname, const std::string& graphname, const std::string& outputname, const std::string& backmatchname, const std::string& updownmatchname);

private:
	void MyFill(THashTable* table, const std::string& name, const std::string& title, int bins, double minx, double maxx, double value);
	void MyFill(THashTable* table, const std::string& name, const std::string& title, int binsx, double minx, double maxx, double valuex,
																						int binsy, double miny, double maxy, double valuey);
	GraphData GetPoints(THashTable* table, const std::string& name);
	CalParams MakeGraph(THashTable* table, int gchan, const GraphData& data);
	
	ChannelMap cmap;
	TSpectrum spec;
	const int max_chan=640; //May need modified if ANASEN is modified
	const double sigma = 1.0, threshold=0.4; //May need modified for each experiment

	const int updown_list[8] = {1, 0, 3, 2, 5, 4, 7, 6}; //matches index -> index of up/down pair for SX3 fronts
};

#endif
