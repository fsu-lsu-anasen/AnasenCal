#include "AnasenArray.h"
#include <fstream>

AnasenArray::AnasenArray()
{
    //Generate our geometry
    for(int i=0; i<s_nSX3; i++)
        m_sx3Barrel.emplace_back(s_sx3BarrelPhiList[i], s_sx3BarrelZ, s_sx3BarrelRhoList[i]);

    for(int i=0; i<s_nBarcDownstream; i++)
        m_barcDownstreamBarrel.emplace_back(s_barcDownstreamBarrelPhiList[i], s_barcDownstreamBarrelZ, s_barcDownstreamBarrelRhoList[i], true);

    for(int i=0; i<s_nBarcUpstream; i++)
        m_barcUpstreamBarrel.emplace_back(s_barcUpstreamBarrelPhiList[i], s_barcUpstreamBarrelZ, s_barcUpstreamBarrelRhoList[i]);
    
    for(int i=0; i<s_nQQQDownstream; i++)
        m_qqqDownstreamCap.emplace_back(s_qqqDownstreamCapPhiList[i], s_qqqDownstreamCapZ);
}

AnasenArray::~AnasenArray() {}

//Will need to change experiment to experiment
void AnasenArray::FillCoordinates(CalibratedEvent& event)
{
    double sx3FrontRatio;
    for(CalibratedSX3Hit& hit : event.barrel)
    {
        //Avoid using both front energies, as in some extreme cases can result in unreliable analysis
        if(hit.frontDownEnergyAdc > hit.frontUpEnergyAdc)
            sx3FrontRatio = 2.0*hit.frontDownEnergyAdc/hit.backEnergy - 1.0;
        else
            sx3FrontRatio = 1.0 - 2.0*hit.frontUpEnergyAdc/hit.backEnergy;
        hit.coordinates = m_sx3Barrel[hit.detectorIndex].GetHitCoordinates(s_sx3FrontStripConversion[hit.frontUpLocalChannel], sx3FrontRatio);
    }

    for(CalibratedQQQHit& hit : event.fqqq)
        hit.coordinates = m_qqqDownstreamCap[hit.detectorIndex].GetHitCoordinates(hit.ringLocalChannel, hit.wedgeLocalChannel);

    for(CalibratedBarcHit& hit : event.barcDown)
        hit.coordinates = m_barcDownstreamBarrel[hit.detectorIndex].GetHitCoordinates(hit.frontLocalChannel);

    for(CalibratedBarcHit& hit : event.barcUp)
        hit.coordinates = m_barcDownstreamBarrel[hit.detectorIndex].GetHitCoordinates(hit.frontLocalChannel);
}

void AnasenArray::WriteDetectorArray(const std::string& filename)
{
    std::ofstream output(filename);
    if(!output.is_open())
    {
        std::cerr<<"Unable to open output file "<<filename<<" at AnasenArray::WriteDetectorArray!"<<std::endl;
        return;
    }

    //Format is X | Y | Z
    ROOT::Math::XYZPoint point;
    for(auto& sx3 : m_sx3Barrel)
    {
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                point = sx3.GetRotatedFrontStripCoordinates(i, j);
                output << point.X() << " " << point.Y() << " " << point.Z() << std::endl;
                point = sx3.GetRotatedBackStripCoordinates(i, j);
                output << point.X() << " " << point.Y() << " " << point.Z() << std::endl;
            }
        }
    }

    for(auto& qqq : m_qqqDownstreamCap)
    {
        for(int i=0; i<16; i++)
        {
            for(int j=0; j<4; j++)
            {
                point = qqq.GetRingCoordinates(i, j);
                output << point.X() << " " << point.Y() << " " << point.Z() << std::endl;
                point = qqq.GetWedgeCoordinates(i, j);
                output << point.X() << " " << point.Y() << " " << point.Z() << std::endl;
            }
        }
    }

    for(auto& barc : m_barcDownstreamBarrel)
    {
        for(int i=0; i<32; i++)
        {
            for(int j=0; j<4; j++)
            {
                point = barc.GetRotatedStripCoordinates(i, j);
                output << point.X() << " " << point.Y() << " " << point.Z() << std::endl;
            }
        }
    }

    for(auto& barc : m_barcUpstreamBarrel)
    {
        for(int i=0; i<32; i++)
        {
            for(int j=0; j<4; j++)
            {
                point = barc.GetRotatedStripCoordinates(i, j);
                output << point.X() << " " << point.Y() << " " << point.Z() << std::endl;
            }
        }
    }

    output.close();
}