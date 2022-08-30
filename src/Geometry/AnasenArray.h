#ifndef ANASEN_ARRAY_H
#define ANASEN_ARRAY_H

#include "SX3Detector.h"
#include "QQQDetector.h"
#include "BarcelonaDetector.h"
#include "CalDict/DataStructs.h"

/*
    For each experiment, please state here where you are defining the origin in the ANASEN array:
    (For example: solid target makes most sense with origin at target, active target with origin at either end of the gas volume, etc.)

    Aug 2022 18Ne(a,p) origin:

*/

class AnasenArray
{
public:
    AnasenArray();
    ~AnasenArray();

    //Takes in CalibratedEvent ref, and fills all detector positions in the event
    void FillCoordinates(CalibratedEvent& event);

    void WriteDetectorArray(const std::string& filename); //So that the array can be drawn

private:
    std::vector<SX3Detector> m_sx3Barrel;
    std::vector<BarcelonaDetector> m_barcDownstreamBarrel;
    std::vector<BarcelonaDetector> m_barcUpstreamBarrel;
    std::vector<QQQDetector> m_qqqDownstreamCap;

    /*** Geometry definitions ***/
    static constexpr int s_nSX3 = 12;
    static constexpr int s_nBarcDownstream = 5;
    static constexpr int s_nBarcUpstream = 6;
    static constexpr int s_nQQQDownstream = 4;

    //All distances in meters, all angles in degrees (will be converted to radians in detector classes)
    static constexpr double s_sx3BarrelZ = 0.1;
    static constexpr double s_barcDownstreamBarrelZ = 0.1;
    static constexpr double s_barcUpstreamBarrelZ = -0.1;
    static constexpr double s_qqqDownstreamCapZ = 0.2;

    static constexpr double s_sx3BarrelPhiList[s_nSX3] = {285.0, 315.0, 345.0, 15.0, 45.0, 75.0,
                                                          105.0, 135.0, 165.0, 195.0, 225.0, 255.0};
    static constexpr double s_sx3BarrelRhoList[s_nSX3] = {0.0890601, 0.0889871, 0.0890354, 0.0890247, 0.0890354, 0.0890354,
                                                          0.0890247, 0.0890354, 0.0890354, 0.0890247, 0.0890354, 0.0890354};
    //Need to convert SX3 front up/down channel to actual strip number on the detector.
    static constexpr int s_sx3FrontStripConversion[8] = {0, 0, 1, 1, 2, 2, 3, 3}; //idk if this is right order (could be flipped)

    //Barc geometry probably not right (especially assignment of det ID to geometry)
    //Rho just guessed as by eye matching hole of QQQs
    static constexpr double s_barcUpstreamBarrelPhiList[s_nBarcUpstream] = {0., 60., 120., 180., 240., 300.};
    static constexpr double s_barcUpstreamBarrelRhoList[s_nBarcUpstream] = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03};

    static constexpr double s_barcDownstreamBarrelPhiList[s_nBarcDownstream] = {60., 120., 180., 240., 300.};
    static constexpr double s_barcDownstreamBarrelRhoList[s_nBarcDownstream] = {0.03, 0.03, 0.03, 0.03, 0.03};


    static constexpr double s_qqqDownstreamCapPhiList[s_nQQQDownstream] = {315.0, 45.0, 135.0, 225.0};
};

#endif