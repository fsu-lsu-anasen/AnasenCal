
#include "SX3Detector.h"
#include <iostream>

/*
  Corner layout for each strip in the un-rotated frame
  0--------------------------1
  |                          |              
  |                          |              
  |                          |              
  |                          |              
  |                          |              x
  2--------------------------3    z<--------X
																						|
																						|
																						|
																						y
*/

//TODO: Implement rng
SX3Detector::SX3Detector(double centerPhi, double centerZ, double centerRho) :
	m_centerPhi(centerPhi), m_centerZ(centerZ), m_centerRho(std::abs(centerRho)), m_norm(1.0,0.0,0.0), m_uniformFraction(0.0, 1.0), m_isSmearing(false)
{
	m_frontStripWidth = s_totalWidth/s_nStrips;
	m_backStripLength = s_stripLength/s_nStrips;

	m_zRotation.SetAngle(m_centerPhi);

	m_frontStripCoords.resize(s_nStrips);
	m_backStripCoords.resize(s_nStrips);
	m_rotFrontStripCoords.resize(s_nStrips);
	m_rotBackStripCoords.resize(s_nStrips);
	for(int i=0; i<s_nStrips; i++)
	{
		m_frontStripCoords[i].resize(s_nCorners);
		m_backStripCoords[i].resize(s_nCorners);
		m_rotFrontStripCoords[i].resize(s_nCorners);
		m_rotBackStripCoords[i].resize(s_nCorners);
	}
	CalculateCorners();
	
}

SX3Detector::~SX3Detector() {}

void SX3Detector::CalculateCorners()
{
	double y_min, y_max, z_min, z_max; 
	for (int s=0; s<s_nStrips; s++)
	{
		y_max = s_totalWidth/2.0 - m_frontStripWidth*s;
		y_min = s_totalWidth/2.0 - m_frontStripWidth*(s+1);
		z_max = m_centerZ + s_stripLength/2.0;
		z_min = m_centerZ - s_stripLength/2.0;
		m_frontStripCoords[s][2].SetXYZ(m_centerRho, y_max, z_max);
		m_frontStripCoords[s][3].SetXYZ(m_centerRho, y_max, z_min);
		m_frontStripCoords[s][0].SetXYZ(m_centerRho, y_min, z_max);
		m_frontStripCoords[s][1].SetXYZ(m_centerRho, y_min, z_min);

		z_max = (m_centerZ - s_stripLength/2.0) + (s+1)*m_backStripLength;
		z_min = (m_centerZ - s_stripLength/2.0) + (s)*m_backStripLength;
		y_max = s_totalWidth/2.0;
		y_min = -s_totalWidth/2.0;
		m_backStripCoords[s][2].SetXYZ(m_centerRho, y_max, z_max);
		m_backStripCoords[s][3].SetXYZ(m_centerRho, y_max, z_min);
		m_backStripCoords[s][0].SetXYZ(m_centerRho, y_min, z_max);
		m_backStripCoords[s][1].SetXYZ(m_centerRho, y_min, z_min);
	}

	for(int s=0; s<s_nStrips; s++)
	{
		m_rotFrontStripCoords[s][0] = m_zRotation*m_frontStripCoords[s][0];
		m_rotFrontStripCoords[s][1] = m_zRotation*m_frontStripCoords[s][1];
		m_rotFrontStripCoords[s][2] = m_zRotation*m_frontStripCoords[s][2];
		m_rotFrontStripCoords[s][3] = m_zRotation*m_frontStripCoords[s][3];

		m_rotBackStripCoords[s][0] = m_zRotation*m_backStripCoords[s][0];
		m_rotBackStripCoords[s][1] = m_zRotation*m_backStripCoords[s][1];
		m_rotBackStripCoords[s][2] = m_zRotation*m_backStripCoords[s][2];
		m_rotBackStripCoords[s][3] = m_zRotation*m_backStripCoords[s][3];
	}
}

ROOT::Math::XYZPoint SX3Detector::GetHitCoordinates(int front_stripch, double front_strip_ratio)
{

	if (!ValidChannel(front_stripch) || !ValidRatio(front_strip_ratio))
		return ROOT::Math::XYZPoint(0,0,0);

	double y;
	if(m_isSmearing)
		y = -s_totalWidth*0.5 + (front_stripch + m_uniformFraction(RandomGenerator::GetInstance().GetGenerator()))*m_frontStripWidth;
	else
		y = -s_totalWidth*0.5 + (front_stripch+0.5)*m_frontStripWidth;

	//recall we're still assuming phi=0 det:
	ROOT::Math::XYZPoint coords(m_centerRho, y, front_strip_ratio*(s_stripLength*0.5) + m_centerZ);
  
	//NOW rotate by appropriate phi
	return m_zRotation*coords;

}

SX3Hit SX3Detector::GetChannelRatio(double theta, double phi)
{

	SX3Hit hit;
	while (phi < 0)
		phi += 2*M_PI;

	//to make the math easier (and the same for each det), rotate the input phi
	//BACKWARD by the phi of the det, s.t. we are looking along x-axis
	phi -= m_centerPhi;

	if (phi > M_PI)
		phi -= 2*M_PI;

	//then we can check easily whether it even hit the detector in phi
	double det_max_phi = std::atan2(s_totalWidth/2,m_centerRho);
	double det_min_phi = -det_max_phi;
  
	if (phi < det_min_phi || phi > det_max_phi)
		return hit;

	//for theta it's not so simple, so we have to go through the typical plane-intersect method
	//first thing's first: we have a fixed x for the entire detector plane:
	double xhit = m_centerRho;
	//thus we find the corresponding y and z for that fixed x, given the input theta and phi:
	double yhit = xhit*tan(phi);
	double zhit = sqrt(xhit*xhit+yhit*yhit)/tan(theta);

	for (int s=0; s<s_nStrips; s++) {
		if (xhit >=m_frontStripCoords[s][0].X() && xhit <=m_frontStripCoords[s][0].X() && //Check min and max x (constant in flat)
			yhit >=m_frontStripCoords[s][1].Y() && yhit <=m_frontStripCoords[s][2].Y() && //Check min and max y
			zhit >=m_frontStripCoords[s][1].Z() && zhit <=m_frontStripCoords[s][0].Z()) //Check min and max z
		{
			hit.frontStripIndex = s;
			hit.frontRatio = (zhit-m_centerZ)/(s_stripLength/2);
			break;
		}
	}

	for (int s=0; s<s_nStrips; s++) {
		if (xhit >= m_backStripCoords[s][0].X() && xhit <= m_backStripCoords[s][0].X() && //Check min and max x (constant in flat)
			yhit >= m_backStripCoords[s][1].Y() && yhit <= m_backStripCoords[s][2].Y() && //Check min and max y
			zhit >= m_backStripCoords[s][1].Z() && zhit <= m_backStripCoords[s][0].Z()) //Check min and max z
		{
			hit.backStripIndex = s;
			break;
		}
	}

	return hit;
}