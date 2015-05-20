#ifndef __VISUALISATION_VTK_HPP_INCLUDED__
#define __VISUALISATION_VTK_HPP_INCLUDED__

#include "visualization.hpp"

//#include <vtkSmartPointer>

#include <iostream>

#include "vtkCubeAxesActor.h"

#include "vtkVersion.h"
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkGlyph3D.h>
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkMath.h"
#include "vtkSphereSource.h"
#include "vtkStructuredGrid.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkArrowSource.h"

#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkInteractorStyleImage.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include "vtkPlaneSource.h"

namespace wif_viz
{


class visualization_vtk_c : public visualization_c
{
public:
	visualization_vtk_c(std::shared_ptr<flow_c> flow, const vector_2d_c & min_range, const vector_2d_c & max_range);

	virtual ~visualization_vtk_c();

	virtual void draw(const std::string & name = "");

private:
	vtkSmartPointer<vtkPoints> construct_points(const vector_2d_c & binning) const;

	vtkSmartPointer<vtkDoubleArray> construct_field(const vector_2d_c & binning, bool scalar) const;

	vtkSmartPointer<vtkStructuredGrid> combine_grid(const vector_2d_c & binning, vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> field) const;

	vtkSmartPointer<vtkStructuredGrid> construct_psi_grid() const;

	vtkSmartPointer<vtkStructuredGrid> construct_phi_grid() const;

	vtkSmartPointer<vtkStructuredGrid> construct_velocity_grid() const;

	vtkSmartPointer<vtkActor> separating_streamlines(vtkSmartPointer<vtkPlaneSource> plane) const;

	vtkSmartPointer<vtkActor> geef_actor_lijnen(std::vector<wif_core::line_2d_c>);

	vtkSmartPointer<vtkActor> geef_actor_punten(std::vector<wif_core::vector_2d_c>);

	void get_data_range(double & min, double & max, std::vector<double> & contours, vtkSmartPointer<vtkStructuredGrid> grid) const;
};


} // namespace wif_viz

#endif // __VISUALISATION_VTK_HPP_INCLUDED__

