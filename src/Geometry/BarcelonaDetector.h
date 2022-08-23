#ifndef BARCELONA_DETECTOR_H
#define BARCELONA_DETECTOR_H

#include "RandomGenerator.h"

#include <cmath>
#include <vector>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/RotationZ.h"
#include "Math/RotationY.h"
#include "Math/Translation3D.h"

class BarcelonaDetector
{
public:
    BarcelonaDetector(double centerPhi, double centerZ, double centerRho, bool isFlippedX = false);
    ~BarcelonaDetector();

    const ROOT::Math::XYZPoint& GetStripCoordinates(int stripch, int corner) { return m_stripCoords[stripch][corner]; }
    const ROOT::Math::XYZPoint& GetRotatedStripCoordinates(int stripch, int corner) { return m_rotStripCoords[stripch][corner]; }
    ROOT::Math::XYZVector GetRotatedNorm() { return m_zRotation * m_norm; }

    void SetPixelSmearing(bool isSmearing) { m_isSmearing = isSmearing; }

    ROOT::Math::XYZPoint GetHitCoordinates(int stripch);
    int GetHitChannel(double theta, double phi);

private:
    bool IsChannelValid(int channel) { return (channel >= 0 && channel < s_nStrips); }
    ROOT::Math::XYZPoint Transform(const ROOT::Math::XYZPoint& point)
    { 
        return m_isFlippedY ? m_zRotation * ( m_translation * (m_yRotation * point)) : m_zRotation * (m_translation * point);
    }
    void CalculateCorners();

    double m_centerRho;
    double m_centerPhi;
    double m_centerZ;

    bool m_isFlippedY;
    bool m_isSmearing;
    std::vector<std::vector<ROOT::Math::XYZPoint>> m_stripCoords;
    std::vector<std::vector<ROOT::Math::XYZPoint>> m_rotStripCoords;

    ROOT::Math::RotationZ m_zRotation;
    ROOT::Math::RotationY m_yRotation;
    ROOT::Math::Translation3D m_translation;
    ROOT::Math::XYZVector m_norm;

	std::uniform_real_distribution<double> m_uniformFraction;

    static constexpr int s_nStrips = 32;
    static constexpr int s_nCorners = 4;

    //All distances defined in meters
    static constexpr double s_totalLength = 0.075;
    static constexpr double s_stripLength = s_totalLength/s_nStrips;
    static constexpr double s_stripWidth = 0.02;
    static constexpr double s_deg2rad = M_PI/180.0;
};

#endif