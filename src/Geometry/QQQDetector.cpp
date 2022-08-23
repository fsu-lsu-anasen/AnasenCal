
#include "QQQDetector.h"

QQQDetector::QQQDetector(double R_in, double R_out, double deltaPhi, double phiCentral, double z, double x, double y) :
	m_innerR(R_in), m_outerR(R_out), m_deltaPhi(deltaPhi), m_centralPhi(phiCentral), m_translation(x,y,z), m_norm(0.0,0.0,1.0),
	m_uniformFraction(0.0, 1.0), m_isSmearing(false)
{
	m_deltaR = (m_outerR - m_innerR)/s_nRings;
	m_deltaPhiWedge = m_deltaPhi/s_nWedges;
	m_zRotation.SetAngle(m_centralPhi);
	m_ringCoords.resize(s_nRings);
	m_wedgeCoords.resize(s_nWedges);
	for(auto& ring : m_ringCoords)
		ring.resize(4);
	for(auto& wedge : m_wedgeCoords)
		wedge.resize(4);

	CalculateCorners();
}

QQQDetector::~QQQDetector() {}

void QQQDetector::CalculateCorners()
{
	double x0, x1, x2, x3;
	double y0, y1, y2, y3;
	double z0, z1, z2, z3;

	//Generate flat ring corner coordinates
	for(int i=0; i<s_nRings; i++)
	{
		x0 = (m_innerR + m_deltaR*(i+1))*std::cos(-m_deltaPhi/2.0);
		y0 = (m_innerR + m_deltaR*(i+1))*std::sin(-m_deltaPhi/2.0);
		z0 = 0.0;
		m_ringCoords[i][0].SetXYZ(x0, y0, z0);

		x1 = (m_innerR + m_deltaR*(i))*std::cos(-m_deltaPhi/2.0);
		y1 = (m_innerR + m_deltaR*(i))*std::sin(-m_deltaPhi/2.0);
		z1 = 0.0;
		m_ringCoords[i][1].SetXYZ(x1, y1, z1);

		x2 = (m_innerR + m_deltaR*(i))*std::cos(m_deltaPhi/2.0);
		y2 = (m_innerR + m_deltaR*(i))*std::sin(m_deltaPhi/2.0);
		z2 = 0.0;
		m_ringCoords[i][2].SetXYZ(x2, y2, z2);

		x3 = (m_innerR + m_deltaR*(i+1))*std::cos(m_deltaPhi/2.0);
		y3 = (m_innerR + m_deltaR*(i+1))*std::sin(m_deltaPhi/2.0);
		z3 = 0.0;
		m_ringCoords[i][3].SetXYZ(x3, y3, z3);
	}

	//Generate flat wedge corner coordinates
	for(int i=0; i<s_nWedges; i++)
	{
		x0 = m_outerR * std::cos(-m_deltaPhi/2.0 + i*m_deltaPhiWedge);
		y0 = m_outerR * std::sin(-m_deltaPhi/2.0 + i*m_deltaPhiWedge);
		z0 = 0.0;
		m_wedgeCoords[i][0].SetXYZ(x0, y0, z0);

		x1 = m_innerR * std::cos(-m_deltaPhi/2.0 + i*m_deltaPhiWedge);
		y1 = m_innerR * std::sin(-m_deltaPhi/2.0 + i*m_deltaPhiWedge);
		z1 = 0.0;
		m_wedgeCoords[i][1].SetXYZ(x1, y1, z1);

		x2 = m_innerR * std::cos(-m_deltaPhi/2.0 + (i+1)*m_deltaPhiWedge);
		y2 = m_innerR * std::sin(-m_deltaPhi/2.0 + (i+1)*m_deltaPhiWedge);
		z2 = 0.0;
		m_wedgeCoords[i][2].SetXYZ(x2, y2, z2);

		x3 = m_outerR * std::cos(-m_deltaPhi/2.0 + (i+1)*m_deltaPhiWedge);
		y3 = m_outerR * std::sin(-m_deltaPhi/2.0 + (i+1)*m_deltaPhiWedge);
		z3 = 0.0;
		m_wedgeCoords[i][3].SetXYZ(x3, y3, z3);
	}

	for(int i=0; i<s_nRings; i++)
	{
		for(int j=0; j<4; j++)
			m_ringCoords[i][j] = TransformCoordinates(m_ringCoords[i][j]);
	}

	for(int i=0; i<s_nWedges; i++)
	{
		for(int j=0; j<4; j++)
			m_wedgeCoords[i][j] = TransformCoordinates(m_wedgeCoords[i][j]);
	}

}

ROOT::Math::XYZPoint QQQDetector::GetTrajectoryCoordinates(double theta, double phi)
{
	double z_to_detector = m_translation.Vect().Z();
	double rho_traj = z_to_detector*std::tan(theta);
	double r_traj = std::sqrt(rho_traj*rho_traj + z_to_detector*z_to_detector);
	double min_rho, max_rho, min_phi, max_phi;

	ROOT::Math::XYZPoint result;

	for(auto& ring : m_ringCoords)
	{
		min_rho = ring[1].Rho();
		max_rho = ring[0].Rho();
		if(rho_traj < max_rho && rho_traj > min_rho)
		{
			for(auto& wedge : m_wedgeCoords)
			{
				min_phi = wedge[0].Phi();
				max_phi = wedge[3].Phi();
				if(phi < min_phi && phi < max_phi)
				{
					result.SetXYZ(std::sin(theta)*std::cos(phi)*r_traj, 
								  std::sin(theta)*std::sin(phi)*r_traj, 
								  std::cos(theta)*r_traj);
					break;
				}
			}
		}
	}
	
	return result;
}

std::pair<int,int> QQQDetector::GetTrajectoryRingWedge(double theta, double phi)
{
	double z_to_detector = m_translation.Vect().Z();
	double rho_traj = z_to_detector*std::tan(theta);
	double min_rho, max_rho, min_phi, max_phi;


	for(int r=0; r<s_nRings; r++)
	{
		auto& ring = m_ringCoords[r];
		min_rho = ring[1].Rho();
		max_rho = ring[0].Rho();
		if(rho_traj < max_rho && rho_traj > min_rho)
		{
			for(int w=0; w<s_nWedges; w++)
			{
				auto& wedge = m_wedgeCoords[w];
				min_phi = wedge[0].Phi();
				max_phi = wedge[3].Phi();
				if(phi > min_phi && phi < max_phi)
					return std::make_pair(r, w);
			}
		}
	}
	
	return std::make_pair(-1, -1);
}

ROOT::Math::XYZPoint QQQDetector::GetHitCoordinates(int ringch, int wedgech)
{
	if(!CheckChannel(ringch) || !CheckChannel(wedgech))
		return ROOT::Math::XYZPoint();

	double r_center, phi_center;

	if(m_isSmearing)
	{
		RandomGenerator& rng  = RandomGenerator::GetInstance();
		r_center  = m_innerR + (m_uniformFraction(rng.GetGenerator())+ringch)*m_deltaR;
		phi_center = -m_deltaPhi/2.0 + (m_uniformFraction(rng.GetGenerator())+wedgech)*m_deltaPhiWedge;
	}
	else
	{
		r_center  = m_innerR + (0.5+ringch)*m_deltaR;
		phi_center = -m_deltaPhi/2.0 + (0.5+wedgech)*m_deltaPhiWedge;
	}
	double x = r_center*std::cos(phi_center);
	double y = r_center*std::sin(phi_center);
	double z = 0;

	ROOT::Math::XYZPoint hit(x, y, z);

	return TransformCoordinates(hit);
}