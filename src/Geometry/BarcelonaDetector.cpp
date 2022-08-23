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
BarcelonaDetector::BarcelonaDetector(double length, double width, double centerPhi, double centerZ, double centerRho) :
    m_centerRho(centerRho), m_centerPhi(centerPhi), m_centerZ(centerZ), m_stripWidth(width), m_totalLength(length),
    m_norm(1.0, 0.0, 0.0), m_uniformFraction(0.0, 1.0)
{
    m_stripLength = m_totalLength/s_nStrips;

    if(m_centerPhi < 0)
        m_centerPhi += 2.0*M_PI;
    if(m_centerRho < 0)
        m_centerRho *= -1.0;
    
    m_zRotation.SetAngle(m_centerPhi);

    m_stripCoords.resize(s_nStrips);
    m_rotStripCoords.resize(s_nCorners);
    for(int i=0; i<s_nStrips; i++)
    {
        m_stripCoords.resize(s_nCorners);
        m_rotStripCoords.resize(s_nCorners);
    }

    CalculateCorners();
}

void BarcelonaDetector::CalculateCorners()
{
    double z_min, z_max, y_min, y_max;
    for(int s=0; s<s_nStrips; s++)
    {
        z_max = (m_centerZ - m_stripLength/2.0) + (s+1)*m_stripLength;
		z_min = (m_centerZ - m_stripLength/2.0) + (s)*m_stripLength;
		y_max = m_stripWidth/2.0;
		y_min = -m_stripWidth/2.0;
		m_stripCoords[s][2].SetXYZ(m_centerRho, y_max, z_max);
		m_stripCoords[s][3].SetXYZ(m_centerRho, y_max, z_min);
		m_stripCoords[s][0].SetXYZ(m_centerRho, y_min, z_max);
		m_stripCoords[s][1].SetXYZ(m_centerRho, y_min, z_min);
    }

    for(int i=0; i<s_nStrips; i++)
    {
        for(int j=0; j<s_nCorners; j++)
            m_rotStripCoords[i][j] = m_zRotation * m_stripCoords[i][j];
    }
}

ROOT::Math::XYZPoint BarcelonaDetector::GetHitCoordinates(int stripch)
{
    //Calculate first in un-rotated det
    double y, z;
    if(m_isSmearing)
    {
        RandomGenerator& rng = RandomGenerator::GetInstance();
        y = -m_stripWidth * 0.5 + m_stripWidth * m_uniformFraction(rng.GetGenerator());
        z = -m_totalLength * 0.5 + m_stripLength * (stripch + m_uniformFraction(rng.GetGenerator()));
    }
    else
    {
        y = 0.0;
        z = -m_totalLength*0.5 + m_stripLength * (stripch + 0.5);
    }
    
    ROOT::Math::XYZPoint coords(m_centerRho, y, z);
    return m_zRotation * coords;
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
	double det_max_phi = std::atan2(m_stripWidth*0.5, m_centerRho);
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