#include <string>
#include <fstream>
#include <vector>
#include <iostream>

void makeMap() {
	std::string outputname = "./etc/AnasenChannelMap_fixedOrientation.txt";
	std::ofstream output(outputname.c_str());

	int gchan;
	std::string detType, detId, channelType;
	int channel;

	int mb2_gchan_offset = 9*32;
	int local_chan;
	for(int i=0; i<9; i++)
	{
		std::cout<<"Running MB1 chipboard "<<i<<std::endl;
		for(int j=0; j<32; j++)
		{
			gchan = i*32+j;
			if(gchan>=0 && gchan<16)
			{
				output<<gchan<<" "<<"FQQQ"<<" "<<"1"<<" "<<"RING NONE "<<j%16<<std::endl;
			}
			else if(gchan>=16 && gchan<32)
			{
				output<<gchan<<" "<<"FQQQ"<<" "<<"0"<<" "<<"RING NONE "<<j%16<<std::endl;
			}
			else if(gchan>=32 && gchan<48)
			{
				output<<gchan<<" "<<"FQQQ"<<" "<<"1"<<" "<<"WEDGE NONE "<<15-j%16<<std::endl;
			}
			else if(gchan>=48 && gchan<64)
			{
				output<<gchan<<" "<<"FQQQ"<<" "<<"0"<<" "<<"WEDGE NONE "<<15-j%16<<std::endl;
			}
			else if(gchan>=64 && gchan<80)
			{
				output<<gchan<<" "<<"FQQQ"<<" "<<"2"<<" "<<"RING NONE "<<j%16<<std::endl;
			}
			else if(gchan>=80 && gchan<96)
			{
				output<<gchan<<" "<<"FQQQ"<<" "<<"3"<<" "<<"RING NONE "<<j%16<<std::endl;
			}
			else if(gchan>=96 && gchan<112)
			{
				output<<gchan<<" "<<"FQQQ"<<" "<<"2"<<" "<<"WEDGE NONE "<<15-j%16<<std::endl;
			}
			else if(gchan>=112 && gchan<128)
			{
				output<<gchan<<" "<<"FQQQ"<<" "<<3<<" "<<"WEDGE NONE "<<15-j%16<<std::endl;
			}
			else if(gchan>=128 && gchan<136)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"B"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=136 && gchan<144)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"A"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=144 && gchan<152)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"C"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=152 && gchan<160)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"D"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=160 && gchan<168)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"F"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=168 && gchan<176)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"F"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=176 && gchan<184)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"A"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=184 && gchan<192)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"B"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=192 && gchan<208)
			{
				output<<gchan<<" "<<"BQQQ"<<" "<<"3"<<" "<<"RING NONE "<<j%16<<std::endl;
			}
			else if(gchan>=208 && gchan<224)
			{
				output<<gchan<<" "<<"BQQQ"<<" "<<"2"<<" "<<"RING NONE "<<j%16<<std::endl;
			}
			else if(gchan>=224 && gchan<232)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"C"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=232 && gchan<240)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"D"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=240 && gchan<248)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"E"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=248 && gchan<256)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"F"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=256 && gchan<272)
			{
				output<<gchan<<" "<<"BQQQ"<<" "<<"1"<<" "<<"RING NONE "<<j%16<<std::endl;
			}
			else if(gchan>=272 && gchan<288)
			{
				output<<gchan<<" "<<"BQQQ"<<" "<<"0"<<" "<<"RING NONE "<<j%16<<std::endl;
			}
		}
	}

	for(int i=0; i<8; i++)
	{
		std::cout<<"Running MB2 chipboard "<<i<<std::endl;
		for(int j=0; j<32; j++)
		{
			gchan = mb2_gchan_offset+i*32+j;
			if(gchan>=288 && gchan<296)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"C"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=296 && gchan<304)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"D"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=304 && gchan<312)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"B"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=312 && gchan<320)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"A"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=320 && gchan<324)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"C"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=324 && gchan<328)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"D"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=328 && gchan<332)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"E"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=332 && gchan<336)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"F"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=336 && gchan<340)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"C"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=340 && gchan<344)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"D"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=344 && gchan<348)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"E"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=348 && gchan<352)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"F"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=352 && gchan<360)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"F"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=360 && gchan<368)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"E"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=368 && gchan<376)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"B"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=376 && gchan<384)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"A"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=384 && gchan<388)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"C"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=388 && gchan<392)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"D"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=392 && gchan<396)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"E"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=396 && gchan<400)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"F"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=400 && gchan<404)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"C"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=404 && gchan<408)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"D"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=408 && gchan<412)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"E"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=412 && gchan<416)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"F"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=416 && gchan<424)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"F"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=424 && gchan<432)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"E"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=432 && gchan<440)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"C"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=440 && gchan<448)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"D"<<" "<<"FRONT";
				local_chan = (7-j%8);
				if((7-j%8) < 4)
				{
					output<<(local_chan%2 == 1 ? " UP " : " DOWN ");
				}
				else
				{
					output<<(local_chan%2 == 0 ? " UP " : " DOWN ");
				}
				output<<7-j%8<<std::endl;
			}
			else if(gchan>=448 && gchan<452)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"A"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=452 && gchan<456)
			{
				output<<gchan<<" "<<"BARREL2A"<<" "<<"B"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=456 && gchan<460)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"A"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=460 && gchan<464)
			{
				output<<gchan<<" "<<"BARREL1A"<<" "<<"B"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=464 && gchan<468)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"A"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=468 && gchan<472)
			{
				output<<gchan<<" "<<"BARREL2B"<<" "<<"B"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=472 && gchan<476)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"A"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=476 && gchan<480)
			{
				output<<gchan<<" "<<"BARREL1B"<<" "<<"B"<<" "<<"BACK NONE "<<3-j%4<<std::endl;
			}
			else if(gchan>=480 && gchan<496)
			{
				output<<gchan<<" "<<"BQQQ"<<" "<<"1"<<" "<<"WEDGE NONE "<<15-j%16<<std::endl;
			}
			else if(gchan>=496 && gchan<512)
			{
				output<<gchan<<" "<<"BQQQ"<<" "<<"0"<<" "<<"WEDGE NONE "<<15-j%16<<std::endl;
			}
			else if(gchan>=512 && gchan<528)
			{
				output<<gchan<<" "<<"BQQQ"<<" "<<"2"<<" "<<"WEDGE NONE "<<15-j%16<<std::endl;
			}
			else if(gchan>=528 && gchan<544)
			{
				output<<gchan<<" "<<"BQQQ"<<" "<<"3"<<" "<<"WEDGE NONE "<<15-j%16<<std::endl;
			}
		}
	}
	output.close();
}
