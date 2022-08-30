#include <TROOT.h>
#include <TGraph2D.h>
#include <TCanvas.h>

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

void drawGeometry()
{
    std::ifstream input("../etc/AnasenGeo.txt");

    std::vector<double> xvec, yvec, zvec;

    if(!input.is_open())
    {
        std::cout << "Cannot open geometry file, please create it." << std::endl;
        return;
    }

    double x, y , z;

    while(input >> x)
    {
        input >> y >> z;
        xvec.push_back(x);
        yvec.push_back(y);
        zvec.push_back(z);
    }

    input.close();

    TGraph2D* graph = new TGraph2D(xvec.size(), xvec.data(), yvec.data(), zvec.data());
    graph->SetTitle("Anasen Geometry");
    graph->SetMarkerStyle(2);
    graph->GetXaxis()->SetTitle("X");
    graph->GetYaxis()->SetTitle("Y");
    graph->GetZaxis()->SetTitle("Z");

    TCanvas* c1 = new TCanvas();

    c1->cd();
    graph->Draw("P");
}