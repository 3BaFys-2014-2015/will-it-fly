#ifndef __VISUALISATION_ROOT_HPP_INCLUDED__
#define __VISUALISATION_ROOT_HPP_INCLUDED__

#include "visualization.hpp"
#include <TROOT.h>
#include <THistPainter.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TStyle.h>
#include <TMarker.h>
#include <TArrow.h>
#include <TGraph.h>
namespace wif_viz
{


class visualization_root_c : public visualization_c
{
public:

	
	visualization_root_c(std::shared_ptr<flow_c> flow, const vector_2d_c & min_range, const vector_2d_c & max_range);

	virtual ~visualization_root_c();

	virtual void draw(const std::string & filename = "");

	void FillBins();

	void addArrows();

	void fillbinStagnatie();
	
	void plotVectors(std::vector<std::vector<double>>,std::vector<double>,std::vector<std::string>);

protected:
	TH2F * velocity, *psi, *phi;
	TH2F * stag;
	TH2F * harr;
	std::shared_ptr<flow_c> Flow;
	//
};


} // namespace wif_viz

#endif // __VISUALISATION_ROOT_HPP_INCLUDED__

