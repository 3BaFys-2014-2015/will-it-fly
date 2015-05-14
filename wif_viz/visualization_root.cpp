#include "visualization_root.hpp"
#include <vector>
#include "TROOT.h"
#include "TArrow.h"
#include "THistPainter.h"
#include "TCanvas.h"
#include "TH2.h"
#include <iostream>
#include "TStyle.h"
#include "TF2.h"
#include "TAttLine.h"
#include "TLegend.h"
namespace wif_viz
{


visualization_root_c::visualization_root_c(std::shared_ptr<flow_c> flow, const vector_2d_c & min_range, const vector_2d_c & max_range) :
	visualization_c(flow, min_range, max_range)
{
	Flow = flow;
}


visualization_root_c::~visualization_root_c()
{
	//
}

void visualization_root_c::draw(const std::string & filename)
{

	/*double contours[contour_locations.size()];

	for(unsigned int i = 0; i < contour_locations.size(); i++)
	{
		contours[i] = contour_locations[i];
	}*/

	FillBins();
	fillbinStagnatie();

	if(velocity_bins.x != 0 && velocity_bins.y != 0)
	{
		std::string filename1 = filename + "velocity.pdf";
		const char * savename = filename1.c_str();



		TCanvas * c = new TCanvas("c", "c", 1000, 1000);
		gStyle->SetOptStat(0);
		//velocity->SetContour(number_of_contours, contour_locations);

		//velocity->Draw();
		stag->Draw("CONT1");
		addArrows();
		velocity->GetXaxis()->SetTitle("x");
		velocity->GetYaxis()->SetTitle("y");
		velocity->GetXaxis()->CenterTitle();
		velocity->GetYaxis()->CenterTitle();
		velocity->SetTitle("velocity");
		c->SaveAs(savename);
		c->Destructor();
	}

	if(psi_bins.x != 0 && psi_bins.y != 0)
	{
		std::string filename2 = filename + "psi.pdf";
		const char * savename = filename2.c_str();
		TCanvas * c = new TCanvas("c", "c", 1000, 1000);
		gStyle->SetOptStat(0);
		gStyle->SetPalette(1);

		/*if(contour_locations.size() >= 1)
		{
			psi->SetContour(contour_locations.size(), contours);
		}*/

		psi->Draw("CONT1Z");
		psi->GetXaxis()->SetTitle("x");
		psi->GetYaxis()->SetTitle("y");
		psi->GetXaxis()->CenterTitle();
		psi->GetYaxis()->CenterTitle();
		psi->SetTitle("psi");
		c->SaveAs(savename);
		c->Destructor();
	}

	if(phi_bins.x != 0 && phi_bins.y != 0)
	{
		std::string filename3 = filename + "phi.pdf";
		const char * savename = filename3.c_str();
		TCanvas * c = new TCanvas("c", "c", 1000, 1000);
		gStyle->SetOptStat(0);
		gStyle->SetPalette(1);
		/*if(contour_locations.size()!=0)
		{
			phi->SetContour(contour_locations.size(), contours);
		}*/
		phi->Draw("CONT1Z");
		phi->GetXaxis()->SetTitle("x");
		phi->GetYaxis()->SetTitle("y");
		phi->GetXaxis()->CenterTitle();
		phi->GetYaxis()->CenterTitle();
		phi->SetTitle("contour plot phi");
		c->SaveAs(savename);
		c->Destructor();
	}

}

void visualization_root_c::FillBins()
{
	if(velocity_bins.x != 0 && velocity_bins.y != 0)
	{
		velocity = new TH2F("harr", "test", velocity_bins.x, min_range.x, max_range.x, velocity_bins.y, min_range.y, max_range.y);
	}

	if(psi_bins.x != 0 && psi_bins.y != 0)
	{
		psi = new TH2F("harr1", "test", psi_bins.x, min_range.x, max_range.x, psi_bins.y, min_range.y, max_range.y);
	}

	if(phi_bins.x != 0 && phi_bins.y != 0)
	{
		phi = new TH2F("harr2", "test", phi_bins.x, min_range.x, max_range.x, phi_bins.y, min_range.y, max_range.y);
	}

	double clip_max = 10;
	double clip_min = -10;
	std::vector<wif_core::vector_2d_c> flows(3);
	flows[0] = velocity_bins;
	flows[1] = psi_bins;
	flows[2] = phi_bins;
	wif_core::vector_2d_c pos;
	wif_core::vector_2d_c vel;

	for(int k = 0; k <= 2; k++)
	{
		double stepsx = (max_range.x - min_range.x) / flows[k].x;
		double stepsy = (max_range.y - min_range.y) / flows[k].y;
		double bincontent;

		for(int i = 1; i <= flows[k].x; i++)
		{
			//std::string harrname= "harr"+(i+4);
			//std::cout<<harrname<<std::endl;
			for(int j = 1; j <= flows[k].y; j++)
			{
				double evaluateX = min_range.x + i * stepsx - stepsx / 2;
				double evaluateY = min_range.y + j * stepsy - stepsy / 2;
				pos.x = evaluateX;
				pos.y = evaluateY;

				if(k == 0 && velocity_bins.x != 0 && velocity_bins.y != 0)
				{
					vel = Flow->get_velocity(pos);

					bincontent = sqrt(pow(vel.x, 2) + pow(vel.y, 2));

					if(bincontent > clip_max)
					{
						bincontent = clip_max;
					}
					else if(bincontent < clip_min)
					{
						bincontent = clip_min;
					}

					velocity->SetBinContent(i, j, bincontent);
				}
				else if(k == 1 && psi_bins.x != 0 && psi_bins.y != 0)
				{
					bincontent = Flow->get_psi(pos);

					if(bincontent > clip_max)
					{
						bincontent = clip_max;
					}
					else if(bincontent < clip_min)
					{
						bincontent = clip_min;
					}

					psi->SetBinContent(i, j, bincontent);
				}
				else if(k == 2 && psi_bins.x != 0 && psi_bins.y != 0)
				{
					bincontent = Flow->get_phi(pos);

					if(bincontent > clip_max)
					{
						bincontent = clip_max;
					}
					else if(bincontent < clip_min)
					{
						bincontent = clip_min;
					}

					phi->SetBinContent(i, j, bincontent);

				}

			}
		}



	}

}

void visualization_root_c::addArrows()
{
	int number_of_arrowbins = 31;
	double stepsx = (max_range.x - min_range.x) / number_of_arrowbins;
	double stepsy = (max_range.y - min_range.y) / number_of_arrowbins;
	wif_core::vector_2d_c pos;
	wif_core::vector_2d_c vel;


	TArrow ARR(0, 0, 0, 0, 0, ">");

	for(int i = 1; i < number_of_arrowbins + 1; i++)
	{
		for(int j = 1; j < number_of_arrowbins + 1; j++)
		{
			double evaluateX = min_range.x + i * stepsx - stepsx / 2;
			double evaluateY = min_range.y + j * stepsy - stepsy / 2;
			pos.x = evaluateX;
			pos.y = evaluateY;

			vel = Flow->get_velocity(pos);
			double Vx = evaluateX + 0.2 * (vel.x);
			double Vy = evaluateY + 0.2 * (vel.y);

			if(Vx > max_range.x || Vy > max_range.y || Vx < min_range.x || Vy < min_range.y || sqrt(pow(vel.x, 2) + pow(vel.y, 2)) <= 0.05)
			{}
			else
			{
				ARR.DrawArrow(evaluateX, evaluateY, Vx, Vy, 0.01, ">");
			}
		}
	}
}

void visualization_root_c::fillbinStagnatie()
{
	double stepsx = (max_range.x - min_range.x) / 51;
	double stepsy = (max_range.y - min_range.y) / 51;
	wif_core::vector_2d_c pos;
	wif_core::vector_2d_c vel;

	if(psi_bins.x != 0 && psi_bins.y != 0)
	{
		stag = new TH2F("harr4", "test", 51, min_range.x, max_range.x, 51, min_range.y, max_range.y);

		for(int i = 1; i < 51 + 1; i++)
		{
			for(int j = 1; j < 51 + 1; j++)
			{
				double evaluateX = min_range.x + i * stepsx - stepsx / 2;
				double evaluateY = min_range.y + j * stepsy - stepsy / 2;
				pos.x = evaluateX;
				pos.y = evaluateY;

				vel = Flow->get_velocity(pos);
				double velstag = sqrt(pow(vel.x, 2) + pow(vel.y, 2));


				if(velstag <= 0.05)
				{
					stag->SetBinContent(i, j, 10);
				}
			}
		}
	}

	Double_t contours[1];
	contours[0] = 10;
	stag->SetContour(1, contours);

}
void visualization_root_c::plotVectors(std::vector<std::vector<double>> yVector,std::vector<double> xVector,std::vector<std::string> legend)
{
	
	unsigned int aantalTH2F=yVector.size();
	TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
	leg->SetHeader("The Legend Title");
	TCanvas * c = new TCanvas("c", "c", 1000, 1000);
	for(unsigned int i=0;i<aantalTH2F;i++)
	{
		std::string number = std::to_string(i+4);
		std::string harrname= "harr"+number;
		const char * harrname2 = harrname.c_str();
		double xBegin=xVector[0];
		double xEnd=xVector[xVector.size()];
		double yBegin=yVector[i][0];
		double yEnd=yVector[i][yVector[i].size()];
		TH2F* harr=new TH2F(harrname2,"test",xVector.size(), xBegin, xEnd, yVector[i].size(), yBegin, yEnd);
		
		
		for(unsigned int j=0;j<yVector[i].size();j++)
		{
		 	for(unsigned int k=0;k<xVector.size();k++)
			{
				harr->SetBinContent(xVector[k],yVector[i][j]);
			}
		}
		
		harr->Draw("scat=1 same");
		std::string legende=legend[i];
		const char * legendname = legende.c_str();
		leg->AddEntry("harr",legendname,"f");
	}
	
	leg->Draw();
	c->SaveAs("test.pdf");
	c->Destructor();
	
}
}// namespace wif_viz
