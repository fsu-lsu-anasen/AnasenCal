#ifndef DEADCHANNELMAP_H
#define DEADCHANNELMAP_H

#include <string>

void GenerateDeadChannelMap(const std::string& zoffset, const std::string& backmatch, const std::string& updownmatch, 
									const std::string& frontbackmatch, const std::string& energycal, const std::string& channelfile, 
									const std::string& deadname);

#endif