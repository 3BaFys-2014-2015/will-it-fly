#include "visualization_vtk.hpp"
#include "visualization.hpp"

#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkAppendPolyData.h>
#include "vtkStructuredGridWriter.h"
#include <vtkClipPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkContourFilter.h>

#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkScalarsToColors.h>
#include <vtkLookupTable.h>
#include "vtkColorTransferFunction.h"
#include <vtkPolyDataMapper.h>
#include "vtkDataSetMapper.h"
#include <vtkProperty.h>
#include <vtkActor.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPlaneSource.h>
#include <vtkDoubleArray.h>

#include <vector>
#include <algorithm>
#include <iostream>

namespace wif_viz
{


visualization_vtk_c::visualization_vtk_c(std::shared_ptr<flow_c> flow, const vector_2d_c & min_range, const vector_2d_c & max_range) :
	visualization_c(flow, min_range, max_range)/*, actors()*/
{
	//this-> actors = vtkSmartPointer<vtkActor>::New();
}


visualization_vtk_c::~visualization_vtk_c()
{
	//
}

uint32_t round_abs(double x)
{
	return std::round(std::abs(x));
}

bool one_is_zero(const vector_2d_c & vec)
{
	return (round_abs(vec.x) * round_abs(vec.y)) != 0;
}

void write_to_file(vtkSmartPointer<vtkStructuredGrid> grid, const std::string & name)
{
	vtkSmartPointer<vtkStructuredGridWriter> writer = vtkStructuredGridWriter::New();
	writer->SetFileName(name.c_str());
	writer->SetInput(grid);
	writer->Write();
}

float red(uint32_t col)
{
	return ((double)((col & 0xff0000) >> 16)) / 255.0;
}

float green(uint32_t col)
{
	return ((double)((col & 0xff00) >> 8)) / 255.0;
}

float blue(uint32_t col)
{
	return ((double)(col & 0xff)) / 255.0;
}

void visualization_vtk_c::draw(const std::string & filename)
{
	//
	if((filename.size() > 1) && (filename.find('.') == std::string::npos))
	{
		/// Enkel uitvoeren als de filename geen extensie heeft.

		if(one_is_zero(psi_bins))
		{
			write_to_file(construct_psi_grid(), filename + "-psi.vtk");
		}

		if(one_is_zero(phi_bins))
		{
			write_to_file(construct_phi_grid(), filename + "-phi.vtk");
		}

		if(one_is_zero(velocity_bins))
		{
			write_to_file(construct_velocity_grid(), filename + "-velocity.vtk");
		}
	}
	else
	{
		vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();

		vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
		renWin->AddRenderer(renderer);
		renWin->SetSize(1000 , 1000);

		vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();

		vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

		iren->SetInteractorStyle(style);
		iren->SetRenderWindow(renWin);

		renderer->SetBackground(0, 0, 0);     // stelt de kleur van de achtergrond in
		renderer->ResetCamera();

		renWin->AddRenderer(renderer);
		iren->SetRenderWindow(renWin);

		if(one_is_zero(psi_bins))
		{
			vtkSmartPointer<vtkStructuredGrid> psi_grid = construct_psi_grid();

			vtkSmartPointer<vtkDataSetMapper> mapperveld = vtkDataSetMapper::New();
			mapperveld->SetInput(psi_grid);

			//vtkSmartPointer<vtkColorTransferFunction> lut = vtkColorTransferFunction::New();

			/*
			double values[2];

			psi_grid->GetPointData()->GetScalars()->GetRange(values);

			std::cout << values[0] << std::endl;

			lut->AddRGBPoint(values[0], 0, 0, 1);
			lut->AddRGBPoint(values[1], 1, 0, 0);
			*/

			double min_v = this->clip_min;
			double max_v = this->clip_max;

			if(clip_min == clip_max)
			{
				double values[2];

				psi_grid->GetPointData()->GetScalars()->GetRange(values);

				min_v = values[0];
				max_v = values[1];
			}
			else
			{
				min_v = this->clip_min;
				max_v = this->clip_max;
			}

			double delta = max_v - min_v;

			double r1 = red(0x00134e5e);
			double g1 = green(0x00134e5e);
			double b1 = blue(0x00134e5e);


			double r2 = red(0x0071b280);
			double g2 = green(0x0071b280);
			double b2 = blue(0x0071b280);

			//

			vtkSmartPointer<vtkColorTransferFunction> lut = vtkColorTransferFunction::New();


			lut->SetNanColor(0, 0, 1);

			lut->AddRGBPoint(min_v, r1, g1, b1);

			lut->AddRGBPoint(std::nextafter(min_v, -std::numeric_limits<double>::infinity()), 1, 0, 0);

			lut->AddRGBPoint(max_v, r2, g2, b2);

			lut->AddRGBPoint(std::nextafter(max_v, std::numeric_limits<double>::infinity()), 0, 0, 1);

			//

			mapperveld->SetLookupTable(lut);

			vtkSmartPointer<vtkActor> actorveld = vtkActor::New();
			actorveld->SetMapper(mapperveld);

			renderer->AddActor(actorveld);
		}

		if(one_is_zero(phi_bins))
		{
			vtkSmartPointer<vtkStructuredGrid> phi_grid = construct_phi_grid();

			vtkSmartPointer<vtkDataSetMapper> mapperveld = vtkDataSetMapper::New();
			mapperveld->SetInput(phi_grid);

			double min_v = this->clip_min;
			double max_v = this->clip_max;

			if(clip_min == clip_max)
			{
				double values[2];

				phi_grid->GetPointData()->GetScalars()->GetRange(values);

				min_v = values[0];
				max_v = values[1];
			}
			else
			{
				min_v = this->clip_min;
				max_v = this->clip_max;
			}

			double delta = max_v - min_v;

			double r1 = red(0x00134e5e);
			double g1 = green(0x00134e5e);
			double b1 = blue(0x00134e5e);


			double r2 = red(0x0071b280);
			double g2 = green(0x0071b280);
			double b2 = blue(0x0071b280);

			//

			vtkSmartPointer<vtkColorTransferFunction> lut = vtkColorTransferFunction::New();


			lut->SetNanColor(0, 0, 1);

			lut->AddRGBPoint(min_v, r1, g1, b1);

			lut->AddRGBPoint(std::nextafter(min_v, -std::numeric_limits<double>::infinity()), 1, 0, 0);

			lut->AddRGBPoint(max_v, r2, g2, b2);

			lut->AddRGBPoint(std::nextafter(max_v, std::numeric_limits<double>::infinity()), 0, 0, 1);


			mapperveld->SetLookupTable(lut);

			vtkSmartPointer<vtkActor> actorveld = vtkActor::New();
			actorveld->SetMapper(mapperveld);

			renderer->AddActor(actorveld);
		}

		if(one_is_zero(velocity_bins))
		{
			vtkSmartPointer<vtkStructuredGrid> velocity_grid = construct_velocity_grid();

			vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();

			vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
			glyph3D->SetSourceConnection(arrowSource->GetOutputPort());
			//glyph3D->SetVectorModeToUseVector();
#if VTK_MAJOR_VERSION <= 5
			glyph3D->SetInput(velocity_grid);
#else
			glyph3D->SetInputData(input);
#endif

			glyph3D->SetColorMode(2);
			//glyph3D->SetScaleModeToScaleByVector();
			//glyph3D->OrientOff();
			glyph3D->SetScaleFactor(.05);
			glyph3D->Update();



			vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			mapper->SetInputConnection(glyph3D->GetOutputPort());

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(mapper);
			actor->GetProperty()->SetRepresentationToSurface();

			renderer->AddActor(actor);
		}

		iren->Initialize();
		renWin->Render();
		iren->Start();

		iren->TerminateApp();
	}

	return;

	std::cout << "hier" << std::endl;

	std::vector<double> conts;

	vtkSmartPointer<vtkPlaneSource> phi_plane = construct_phi_plane();
	vtkSmartPointer<vtkPlaneSource> psi_plane = construct_psi_plane();

	double phiRange[2], psiRange[2];
	phi_plane->GetOutput()->GetPointData()->GetScalars()->GetRange(phiRange);
	psi_plane->GetOutput()->GetPointData()->GetScalars()->GetRange(psiRange);

	std::vector<double> contvec_phi, contvec_psi;

	double delta_phi = std::abs((phiRange[1] - phiRange[0]) / (20));
	double delta_psi = std::abs((psiRange[1] - psiRange[0]) / (20));

	for(int i = 0; i < 20; ++i)
	{
		contvec_phi.push_back(phiRange[0] + delta_phi * i);
		std::cout << contvec_phi[i] << std::endl;
		contvec_psi.push_back(psiRange[0] + delta_psi * i);
		std::cout << contvec_psi[i] << std::endl;
	}

	std::cout << "hier" << std::endl;
	contour_plot(phi_plane, contvec_phi);
	std::cout << "hier" << std::endl;
	contour_plot(psi_plane, contvec_psi);
	std::cout << "hier" << std::endl;

	return;

}

vtkSmartPointer<vtkPoints> visualization_vtk_c::construct_points(const vector_2d_c & binning) const
{
	uint32_t bin_x = round_abs(binning.x);
	uint32_t bin_y = round_abs(binning.y);

	uint32_t s = (bin_x + 1) * (bin_y + 1);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->Allocate(s);

	const vector_2d_c delta = max_range - min_range;
	const vector_2d_c delta_it(delta.x / bin_x, delta.y / bin_y);

	for(uint32_t i = 0; i <= bin_x; i++)
	{
		for(uint32_t j = 0; j <= bin_y; j++)
		{
			vector_2d_c pos(i * delta_it.x, j * delta_it.y);

			pos += min_range;

			points->InsertNextPoint(pos.x, pos.y, 0.0);
		}
	}

	return points;
}

vtkSmartPointer<vtkDoubleArray> visualization_vtk_c::construct_field(const vector_2d_c & binning, bool scalar) const
{
	//uint32_t bin_x = round_abs(binning.x);
	//uint32_t bin_y = round_abs(binning.y);

	vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();
	//vectors->Allocate((bin_x + 1) * (bin_y + 1));

	if(scalar)
	{
		vectors->SetNumberOfComponents(1);
		vectors->SetName("ScalarField");
	}
	else
	{
		vectors->SetNumberOfComponents(3);
		vectors->SetName("VelocityField");
	}

	return vectors;
}

vtkSmartPointer<vtkStructuredGrid> visualization_vtk_c::combine_grid(const vector_2d_c & binning, vtkSmartPointer<vtkPoints> points, vtkSmartPointer<vtkDoubleArray> field) const
{
	uint32_t bin_x = round_abs(binning.x);
	uint32_t bin_y = round_abs(binning.y);

	vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();

	grid->SetDimensions(bin_x + 1, bin_y + 1, 1);
	grid->SetPoints(points);

	if(field->GetNumberOfComponents() > 1)
	{
		grid->GetPointData()->SetVectors(field);
	}
	else
	{
		grid->GetPointData()->SetScalars(field);
	}

	return grid;
}

vtkSmartPointer<vtkStructuredGrid> visualization_vtk_c::construct_psi_grid() const
{
	vtkSmartPointer<vtkPoints> points = construct_points(psi_bins);
	vtkSmartPointer<vtkDoubleArray> field = construct_field(psi_bins, true);

	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);

		double t = clip_value(flow->get_psi(pos));

		//field->InsertNextTuple(&t);
		field->InsertNextValue(t);
	}

	return combine_grid(psi_bins, points, field);
}

vtkSmartPointer<vtkStructuredGrid> visualization_vtk_c::construct_phi_grid() const
{
	vtkSmartPointer<vtkPoints> points = construct_points(phi_bins);
	vtkSmartPointer<vtkDoubleArray> field = construct_field(phi_bins, true);

	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);

		double t = clip_value(flow->get_phi(pos));

		//field->InsertNextTuple1(t);
		field->InsertNextValue(t);
	}

	return combine_grid(phi_bins, points, field);
}

vtkSmartPointer<vtkStructuredGrid> visualization_vtk_c::construct_velocity_grid() const
{
	vtkSmartPointer<vtkPoints> points = construct_points(velocity_bins);
	vtkSmartPointer<vtkDoubleArray> field = construct_field(velocity_bins, false);

	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);
		const vector_2d_c velo = flow->get_velocity(pos);
		double v[3] = {velo.x, velo.y, 0.0};

		field->InsertNextTuple(v);
	}

	return combine_grid(velocity_bins, points, field);
}

vtkSmartPointer<vtkPlaneSource> visualization_vtk_c::construct_phi_plane() const
{
	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();

	plane->SetXResolution(round_abs(phi_bins.x));
	plane->SetYResolution(round_abs(phi_bins.y));
	plane->SetOrigin(min_range.x, min_range.y, 0);
	plane->SetPoint1(max_range.x, min_range.y, 0);
	plane->SetPoint2(min_range.x, max_range.y, 0);
	plane->Update();
	vtkSmartPointer<vtkDoubleArray> field = construct_field(phi_bins, true);
	vtkSmartPointer<vtkPoints> points = plane->GetOutput()->GetPoints();

	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);

		double t = clip_value(flow->get_phi(pos));

		//field->InsertNextTuple1(t);
		field->InsertNextValue(t);
	}

	plane->GetOutput()->GetPointData()->SetScalars(field);
	return plane;
}

vtkSmartPointer<vtkPlaneSource> visualization_vtk_c::construct_psi_plane() const
{
	vtkSmartPointer<vtkPlaneSource> plane = vtkSmartPointer<vtkPlaneSource>::New();

	plane->SetXResolution(round_abs(psi_bins.x));
	plane->SetYResolution(round_abs(psi_bins.y));
	plane->SetOrigin(min_range.x, min_range.y, 0);
	plane->SetPoint1(max_range.x, min_range.y, 0);
	plane->SetPoint2(min_range.x, max_range.y, 0);
	plane->Update();
	vtkSmartPointer<vtkDoubleArray> field = construct_field(psi_bins, true);
	vtkSmartPointer<vtkPoints> points = plane->GetOutput()->GetPoints();

	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);

		double t = clip_value(flow->get_psi(pos));

		if(t > vtkMax)
		{
			t = vtkMax;
		}

		if(t < -vtkMax)
		{
			t = -vtkMax;
		}

		//field->InsertNextTuple1(t);
		field->InsertNextValue(t);
	}

	plane->GetOutput()->GetPointData()->SetScalars(field);
	return plane;
}


void visualization_vtk_c::contour_plot(vtkSmartPointer<vtkPlaneSource> plane, std::vector<double> contlvls) const //int ncont
{
	std::sort(contlvls.begin(), contlvls.end());

	//int Nvec = contlvls.size();

	/*while (*(contlvls.end()) >= vtkMax && contlvls.size() > 2)
	{
		contlvls.pop_back();
		std::cout << "removed" << std::endl;
	}*/

	double scalarRange[2];
	double planeRange[2];
	//plane->GetOutput()->GetPointData()->GetScalars()->GetRange(scalarRange);
	plane->GetOutput()->GetPointData()->GetScalars()->GetRange(planeRange);
	scalarRange[0] = *(contlvls.begin());
	scalarRange[1] = *(contlvls.end());

	vtkSmartPointer<vtkAppendPolyData> appendFilledContours = vtkSmartPointer<vtkAppendPolyData>::New();

	//int numberOfContours = ncont;
	int numberOfContours = contlvls.size();
	contlvls.push_back(planeRange[1] + 0.5);

	//double delta =(scalarRange[1] - scalarRange[0]) /static_cast<double> (numberOfContours - 1);

	// Keep the clippers alive
	std::vector<vtkSmartPointer<vtkClipPolyData> > clippersLo;
	std::vector<vtkSmartPointer<vtkClipPolyData> > clippersHi;

	for(int i = 0; i < numberOfContours; i++)
	{
		//double valueLo = scalarRange[0] + static_cast<double> (i) * delta;
		double valueLo = contlvls[i];
		//double valueHi = scalarRange[0] + static_cast<double> (i + 1) * delta;
		double valueHi = contlvls[i + 1];
		clippersLo.push_back(vtkSmartPointer<vtkClipPolyData>::New());
		clippersLo[i]->SetValue(valueLo);

		if(i == 0)
		{
			clippersLo[i]->SetInputConnection(plane->GetOutputPort());
		}
		else
		{
			clippersLo[i]->SetInputConnection(clippersHi[i - 1]->GetOutputPort(1));
		}

		clippersLo[i]->InsideOutOff();
		clippersLo[i]->Update();

		clippersHi.push_back(vtkSmartPointer<vtkClipPolyData>::New());
		clippersHi[i]->SetValue(valueHi);
		clippersHi[i]->SetInputConnection(clippersLo[i]->GetOutputPort());
		clippersHi[i]->GenerateClippedOutputOn();
		clippersHi[i]->InsideOutOn();
		clippersHi[i]->Update();

		if(clippersHi[i]->GetOutput()->GetNumberOfCells() == 0)
		{
			continue;
		}

		vtkSmartPointer<vtkFloatArray> cd = vtkSmartPointer<vtkFloatArray>::New();
		cd->SetNumberOfComponents(1);
		cd->SetNumberOfTuples(clippersHi[i]->GetOutput()->GetNumberOfCells());
		cd->FillComponent(0, valueLo);

		clippersHi[i]->GetOutput()->GetCellData()->SetScalars(cd);
		appendFilledContours->AddInputConnection(clippersHi[i]->GetOutputPort());
	}

	vtkSmartPointer<vtkCleanPolyData> filledContours = vtkSmartPointer<vtkCleanPolyData>::New();
	filledContours->SetInputConnection(appendFilledContours->GetOutputPort());

	vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
	lut->SetNumberOfTableValues(numberOfContours + 1);
	lut->Build();
	vtkSmartPointer<vtkPolyDataMapper> contourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	contourMapper->SetInputConnection(filledContours->GetOutputPort());
	contourMapper->SetScalarRange(scalarRange[0], scalarRange[1]);
	contourMapper->SetScalarModeToUseCellData();
	contourMapper->SetLookupTable(lut);

	vtkSmartPointer<vtkActor> contourActor =
	    vtkSmartPointer<vtkActor>::New();
	contourActor->SetMapper(contourMapper);
	contourActor->GetProperty()->SetInterpolationToFlat();

	vtkSmartPointer<vtkContourFilter> contours =
	    vtkSmartPointer<vtkContourFilter>::New();
	contours->SetInputConnection(filledContours->GetOutputPort());
	contours->GenerateValues(numberOfContours, scalarRange[0], scalarRange[1]);

	vtkSmartPointer<vtkPolyDataMapper> contourLineMapperer =
	    vtkSmartPointer<vtkPolyDataMapper>::New();
	contourLineMapperer->SetInputConnection(contours->GetOutputPort());
	contourLineMapperer->SetScalarRange(scalarRange[0], scalarRange[1]);
	contourLineMapperer->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> contourLineActor =
	    vtkSmartPointer<vtkActor>::New();
	contourLineActor->SetMapper(contourLineMapperer);
	contourLineActor->GetProperty()->SetLineWidth(2);

	// The usual renderer, render window and interactor
	vtkSmartPointer<vtkRenderer> ren1 =
	    vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin =
	    vtkSmartPointer<vtkRenderWindow>::New();
	vtkSmartPointer<vtkRenderWindowInteractor>
	iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

	ren1->SetBackground(.1, .2, .3);
	renWin->AddRenderer(ren1);
	iren->SetRenderWindow(renWin);

	// Add the actors
	ren1->AddActor(contourActor);
	ren1->AddActor(contourLineActor);

	// Begin interaction
	renWin->Render();
	iren->Start();
}


} // namespace wif_viz
