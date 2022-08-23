#ifndef BARCELONA_DETECTOR_H
#define BARCELONA_DETECTOR_H

#include "RandomGenerator.h"

#include <cmath>
#include <vector>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/RotationZ.h"

class BarcelonaDetector
{
public:
    BarcelonaDetector(double length, double width, double centerPhi, double centerZ, double centerRho);
    ~BarcelonaDetector();

    const ROOT::Math::XYZPoint& GetStripCoordinates(int stripch, int corner) { return m_stripCoords[stripch][corner]; }
    const ROOT::Math::XYZPoint& GetRotatedStripCoordinates(int stripch, int corner) { return m_rotStripCoords[stripch][corner]; }
    ROOT::Math::XYZVector GetRotatedNorm() { return m_zRotation * m_norm; }

    void SetPixelSmearing(bool isSmearing) { m_isSmearing = isSmearing; }

    ROOT::Math::XYZPoint GetHitCoordinates(int stripch);
    int GetHitChannel(double theta, double phi);

private:
    bool IsChannelValid(int channel) { return (channel >= 0 && channel < s_nStrips); }
    void CalculateCorners();

    double m_centerRho;
    double m_centerPhi;
    double m_centerZ;

    double m_stripLength;
    double m_stripWidth;
    double m_totalLength;

    bool m_isSmearing;
    std::vector<std::vector<ROOT::Math::XYZPoint>> m_stripCoords;
    std::vector<std::vector<ROOT::Math::XYZPoint>> m_rotStripCoords;

    ROOT::Math::RotationZ m_zRotation;
    ROOT::Math::XYZVector m_norm;

	std::uniform_real_distribution<double> m_uniformFraction;

    static constexpr int s_nStrips = 32;
    static constexpr int s_nCorners = 4;
};

#endif