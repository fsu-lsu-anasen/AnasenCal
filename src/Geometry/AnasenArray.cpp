#include "AnasenArray.h"

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

//Will need to change experiment to experiment
void AnasenArray::FillCoordinates(CalibratedEvent& event)
{
    double sx3FrontRatio;
    for(CalibratedSX3Hit& hit : event.barrel)
    {
        if(hit.frontDownEnergyAdc > hit.frontUpEnergyAdc)
        {
            //Idk some updown shit, depends on where origin is
        }
        else
        {
            //Idk some updown shit, depends on where origin is
        }
        hit.coordinates = m_sx3Barrel[hit.detectorIndex].GetHitCoordinates(s_sx3FrontStripConversion[hit.frontUpLocalChannel], sx3FrontRatio);
    }

    for(CalibratedQQQHit& hit : event.fqqq)
        hit.coordinates = m_qqqDownstreamCap[hit.detectorIndex].GetHitCoordinates(hit.ringLocalChannel, hit.wedgeLocalChannel);

    for(CalibratedBarcHit& hit : event.barcDown)
        hit.coordinates = m_barcDownstreamBarrel[hit.detectorIndex].GetHitCoordinates(hit.frontLocalChannel);

    for(CalibratedBarcHit& hit : event.barcUp)
        hit.coordinates = m_barcDownstreamBarrel[hit.detectorIndex].GetHitCoordinates(hit.frontLocalChannel);
}