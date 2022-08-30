#include "BarcelonaDetector.h"
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
BarcelonaDetector::BarcelonaDetector(double centerPhi, double centerZ, double centerRho, bool isFlippedY) :
    m_centerRho(std::abs(centerRho)), m_centerPhi(centerPhi*s_deg2rad), m_centerZ(centerZ), m_isFlippedY(isFlippedY), 
    m_norm(1.0, 0.0, 0.0), m_uniformFraction(0.0, 1.0)
{
    m_zRotation.SetAngle(m_centerPhi);
    m_yRotation.SetAngle(M_PI);
    m_translation.SetXYZ(m_centerRho, 0., m_centerZ);

    m_stripCoords.resize(s_nStrips);
    m_rotStripCoords.resize(s_nStrips);
    for(int i=0; i<s_nStrips; i++)
    {
        m_stripCoords[i].resize(s_nCorners);
        m_rotStripCoords[i].resize(s_nCorners);
    }

    CalculateCorners();
}

BarcelonaDetector::~BarcelonaDetector() {}

void BarcelonaDetector::CalculateCorners()
{
    double z_min, z_max, y_min, y_max;
    for(int s=0; s<s_nStrips; s++)
    {
        z_max = (-s_stripLength*0.5) + (s+1)*s_stripLength;
		z_min = (-s_stripLength*0.5) + (s)*s_stripLength;
		y_max = s_stripWidth/0.5;
		y_min = -s_stripWidth/0.5;
		m_stripCoords[s][2].SetXYZ(0.0, y_max, z_max);
		m_stripCoords[s][3].SetXYZ(0.0, y_max, z_min);
		m_stripCoords[s][0].SetXYZ(0.0, y_min, z_max);
		m_stripCoords[s][1].SetXYZ(0.0, y_min, z_min);
    }

    for(int i=0; i<s_nStrips; i++)
    {
        for(int j=0; j<s_nCorners; j++)
            m_rotStripCoords[i][j] = Transform(m_stripCoords[i][j]);
    }
}

ROOT::Math::XYZPoint BarcelonaDetector::GetHitCoordinates(int stripch)
{
    //Calculate first in un-rotated det
    double y, z;
    if(m_isSmearing)
    {
        RandomGenerator& rng = RandomGenerator::GetInstance();
        y = -s_stripWidth * 0.5 + s_stripWidth * m_uniformFraction(rng.GetGenerator());
        z = -s_totalLength * 0.5 + s_stripLength * (stripch + m_uniformFraction(rng.GetGenerator()));
    }
    else
    {
        y = 0.0;
        z = -s_totalLength*0.5 + s_stripLength * (stripch + 0.5);
    }
    
    ROOT::Math::XYZPoint coords(m_centerRho, y, z);
    return Transform(coords);
}

int BarcelonaDetector::GetHitChannel(double theta, double phi)
{
    int hitChannel = -1;
	while (phi < 0)
		phi += 2*M_PI;

	//to make the math easier (and the same for each det), rotate the input phi
	//BACKWARD by the phi of the det, s.t. we are looking along x-axis
	phi -= m_centerPhi;

	if (phi > M_PI)
		phi -= 2*M_PI;

	//then we can check easily whether it even hit the detector in phi
	double det_max_phi = std::atan2(s_stripWidth*0.5, m_centerRho);
	double det_min_phi = -det_max_phi;
  
	if (phi < det_min_phi || phi > det_max_phi)
		return hitChannel;

	//for theta it's not so simple, so we have to go through the typical plane-intersect method
	//first thing's first: we have a fixed x for the entire detector plane:
	double xhit = m_centerRho;
	//thus we find the corresponding y and z for that fixed x, given the input theta and phi:
	double yhit = xhit*std::tan(phi);
	double zhit = std::sqrt(xhit*xhit+yhit*yhit)/std::tan(theta);

	for (int s=0; s<s_nStrips; s++)
    {
		if (xhit >=m_stripCoords[s][0].X() && xhit <=m_stripCoords[s][0].X() && //Check min and max x (constant in flat)
			yhit >=m_stripCoords[s][1].Y() && yhit <=m_stripCoords[s][2].Y() && //Check min and max y
			zhit >=m_stripCoords[s][1].Z() && zhit <=m_stripCoords[s][0].Z()) //Check min and max z
		{
			hitChannel = s;
			break;
		}
	}

	return hitChannel;
}