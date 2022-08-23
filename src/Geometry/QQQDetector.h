/*
	QQQDetector.h
	Class implementing geometry for QQQDetector where the detector is perpendicular to the beam axis.
	Detector is first generated centered on the x-axis (phi=0)
	Coordinate convention : +z is downstream, -z is upstream. +y is vertically down in the lab.
*/
#ifndef QQQDETECTOR_H
#define QQQDETECTOR_H

#include "RandomGenerator.h"

#include <cmath>
#include <vector>

#include "Math/Point3D.h"
#include "Math/Vector3D.h"
#include "Math/RotationZ.h"
#include "Math/Translation3D.h"

class QQQDetector
{
public:
	QQQDetector(double R_in, double R_out, double deltaPhi, double phiCentral, double z, double x=0, double y=0);
	~QQQDetector();
	const ROOT::Math::XYZPoint& GetRingCoordinates(int ringch, int corner) { return m_ringCoords[ringch][corner]; }
	const ROOT::Math::XYZPoint& GetWedgeCoordinates(int wedgech, int corner) { return m_wedgeCoords[wedgech][corner]; }
	const ROOT::Math::XYZVector& GetNorm() { return m_norm; }
	ROOT::Math::XYZPoint GetTrajectoryCoordinates(double theta, double phi);
	std::pair<int, int> GetTrajectoryRingWedge(double theta, double phi);
	ROOT::Math::XYZPoint GetHitCoordinates(int ringch, int wedgech);
	
	void SetSmearing(bool isSmearing) { m_isSmearing = isSmearing; }

	int GetNumberOfRings() { return s_nRings; }
	int GetNumberOfWedges() { return s_nWedges; }

private:

	bool CheckChannel(int ch) { return (ch >=0 && ch < s_nRings); }
	bool CheckCorner(int corner) { return (corner >=0 && corner < 4); }

	void CalculateCorners();
	ROOT::Math::XYZPoint TransformCoordinates(ROOT::Math::XYZPoint& vector) { return m_translation * (m_zRotation * vector) ; }

	double m_innerR;
	double m_outerR;
	double m_deltaR;
	double m_deltaPhi;
	double m_deltaPhiWedge;
	double m_centralPhi;

	std::vector<std::vector<ROOT::Math::XYZPoint>> m_ringCoords, m_wedgeCoords;
	ROOT::Math::Translation3D m_translation;
	ROOT::Math::XYZVector m_norm;
	ROOT::Math::RotationZ m_zRotation;

	std::uniform_real_distribution<double> m_uniformFraction;
	bool m_isSmearing;

	static constexpr int s_nRings = 16;
	static constexpr int s_nWedges = 16;
	static constexpr double s_deg2rad  = M_PI/180.0;
};

#endif