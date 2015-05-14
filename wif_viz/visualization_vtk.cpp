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
#include <vtkScalarBarActor.h>

#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPlaneSource.h>
#include <vtkDoubleArray.h>

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
#include <vtkStreamLine.h>


#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

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

void visualization_vtk_c::draw_ivo(const std::string & filename)
{
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

		if(this->airfoil)
		{
			renderer->AddActor(geef_actor_lijnen(this->airfoil->get_lines()));
		}

		iren->Initialize();
		renWin->Render();
		iren->Start();

		iren->TerminateApp();
	}
}

void visualization_vtk_c::draw(const std::string & filename)
{
	//


	//return;

	std::cout << "hier" << std::endl;

	std::vector<double> psi_pot_vec;

	vtkSmartPointer<vtkPlaneSource> psi_plane = construct_psi_plane();

	double psiRange[2];
	psi_plane->GetOutput()->GetPointData()->GetScalars()->GetRange(psiRange);

	double delta_psi = std::abs((psiRange[1] - psiRange[0]) / (20));

	//std::cout << phiRange[0] << ", " << phiRange[1] << std::endl;
	for(int i = 0; i < 20; ++i)
	{
		psi_pot_vec.push_back(psiRange[0] + (delta_psi * i)) ;
		//std::cout << contvec_psi[i] << std::endl;
	}

	//contour_plot(phi_plane, 20);//contvec_phi);
	contour_plot(psi_plane, psi_pot_vec);//contvec_psi);



	streamlines_plot(construct_velocity_grid(), 20);

	arrow_plot();




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

		if(velo.get_length_sq() < stagnation_tolerance)
		{
			this->stagnation_point.push_back(vector_2d_c(pos));
		}

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
	vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkPoints> points = plane->GetOutput()->GetPoints();

	double icout = points->GetNumberOfPoints() * 0.01;


	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);

		double t = clip_value(flow->get_phi(pos));

		if(i > icout)
		{
			std::cout << i << ": " << pos.x << ", " << pos.y << "// ->" << t <<  std::endl;
			icout = icout + (points->GetNumberOfPoints() / 100);
		}

		/*if (i > icout)
		{
			std::cout << i <<": " << pos.x << ", " << pos.y << "// ->" << t <<  std::endl;
			icout = icout + (points->GetNumberOfPoints()/100);
		}*/
		/*if(t > vtkMax)
		{
			t = vtkMax;
		}

		if(t < -vtkMax)
		{
			t = -vtkMax;
		}*/

		//std::cout << t << std::endl;
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
	vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkPoints> points = plane->GetOutput()->GetPoints();

	double icout = points->GetNumberOfPoints() * 0.01;

	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);


		double t = clip_value(flow->get_psi(pos));
		//std::cout <<  i <<": " << pos.x << ", " << pos.y << "// ->" << t << std::endl;


		if(i > icout)
		{
			std::cout << i << ": " << pos.x << ", " << pos.y << "// ->" << t <<  std::endl;
			icout = icout + (points->GetNumberOfPoints() * 0.01);
		}

		/*if(t > vtkMax)
		{
			t = vtkMax;
		}

		if(t < -vtkMax)
		{
			t = -vtkMax;
		}*/

		//std::cout << t << std::endl;
		//field->InsertNextTuple1(t);
		field->InsertNextValue(t);
	}

	plane->GetOutput()->GetPointData()->SetScalars(field);
	return plane;
}


void visualization_vtk_c::contour_plot(vtkSmartPointer<vtkPlaneSource> plane, std::vector<double> contlvls) const //int ncont
{
	std::sort(contlvls.begin(), contlvls.end());

	double scalarRange[2];
	double planeRange[2];
	//plane->GetOutput()->GetPointData()->GetScalars()->GetRange(scalarRange);
	plane->GetOutput()->GetPointData()->GetScalars()->GetRange(planeRange);
	scalarRange[0] = contlvls[contlvls.size() - 1];
	scalarRange[1] = contlvls[0];

	int Nvec = contlvls.size();

	while(contlvls[Nvec - 1] >= planeRange[1] && Nvec > 2)
	{
		contlvls.pop_back();
		Nvec = contlvls.size();
		std::cout << "removed" << std::endl;
	}



	vtkSmartPointer<vtkAppendPolyData> appendFilledContours = vtkSmartPointer<vtkAppendPolyData>::New();

	//int numberOfContours = ncont;
	int numberOfContours = contlvls.size();
	contlvls.push_back(planeRange[1] + 0.5);

	//double delta = (scalarRange[1] - scalarRange[0]) / static_cast<double>(numberOfContours - 1);

	// Keep the clippers alive
	std::vector<vtkSmartPointer<vtkClipPolyData> > clippersLo;
	std::vector<vtkSmartPointer<vtkClipPolyData> > clippersHi;

	for(int i = 0; i < numberOfContours; i++)
	{
		//double valueLo = scalarRange[0] + static_cast<double>(i) * delta;
		double valueLo = contlvls[i];
		//double valueHi = scalarRange[0] + static_cast<double>(i + 1) * delta;
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
	contourMapper->SetScalarRange(planeRange[0], planeRange[1]);
	contourMapper->SetScalarModeToUseCellData();
	contourMapper->SetLookupTable(lut);

	//schaal bar
	vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar->SetLookupTable(lut);
	//scalarBar->GetLabelTextProperty()->SetColor(0,0,0); //werkt blijkbaar niet

	vtkSmartPointer<vtkActor> contourActor = vtkSmartPointer<vtkActor>::New();
	contourActor->SetMapper(contourMapper);
	contourActor->GetProperty()->SetInterpolationToFlat();

	vtkSmartPointer<vtkContourFilter> contours = vtkSmartPointer<vtkContourFilter>::New();
	contours->SetInputConnection(filledContours->GetOutputPort());
	contours->GenerateValues(numberOfContours, scalarRange[0], scalarRange[1]);

	vtkSmartPointer<vtkPolyDataMapper> contourLineMapperer = vtkSmartPointer<vtkPolyDataMapper>::New();
	contourLineMapperer->SetInputConnection(contours->GetOutputPort());
	contourLineMapperer->SetScalarRange(scalarRange[0], scalarRange[1]);
	contourLineMapperer->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> contourLineActor = vtkSmartPointer<vtkActor>::New();
	contourLineActor->SetMapper(contourLineMapperer);
	contourLineActor->GetProperty()->SetLineWidth(2);

	// The usual renderer, render window and interactor
	vtkSmartPointer<vtkRenderer> ren1 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	vtkSmartPointer<vtkRenderWindowInteractor>
	iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

	ren1->SetBackground(1, 1, 1);
	renWin->AddRenderer(ren1);
	iren->SetRenderWindow(renWin);

	// assen
	vtkSmartPointer<vtkCubeAxesActor> assen = axis(plane, ren1);
	assen->GetProperty()->SetColor(0, 0, 0);
	assen->SetXTitle("x");
	assen->SetYTitle("y");

	// Add the actors
	ren1->AddActor2D(scalarBar);
	ren1->AddActor(contourActor);
	ren1->AddActor(contourLineActor);
	ren1->AddActor(assen);
	ren1->ResetCamera();

	// Begin interaction
	renWin->Render();
	iren->Start();
}


vtkSmartPointer<vtkActor> visualization_vtk_c::geef_actor_lijnen(std::vector<wif_core::line_2d_c> mylines)
{
	unsigned char black[3] = {0, 0, 0};
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

	int N = mylines.size() ;                                                 // aantal lijnen dat getekend moet worden.

	for(int j = 0; j < N; j = j + 1)
	{
		pts->InsertNextPoint(mylines[j].begin.x, mylines[j].begin.y, 0);    // beginpunt lijn
		pts->InsertNextPoint(mylines[j].end.x, mylines[j].end.y, 0);			// eindpunt lijn     opslagen in vtk object
		colors->InsertNextTupleValue(black);
	}



	std::vector<vtkSmartPointer<vtkLine> > line;

// maak cell array om de lijnen in te steken
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();

	for(int j = 0; j < N; j = j + 1)
	{
		line.push_back(vtkSmartPointer<vtkLine>::New());

		line[j]->GetPointIds()->SetId(0, 2 * j); //het tweede getal is de index van O1 in de vtkPoints
		line[j]->GetPointIds()->SetId(1, (2 * j) + 1); //het tweede getal is de index van P1 in de vtkPoints

		lines->InsertNextCell(line[j]);
	}


// maak dataset
	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();

// voeg punten aan dataset toe
	linesPolyData->SetPoints(pts);

// voeg lijnen aan dataset toe
	linesPolyData->SetLines(lines);

// kleur de lijnen door elk component aan pollydata te koppelen aan colors
	linesPolyData->GetCellData()->SetScalars(colors);



	vtkSmartPointer<vtkPolyDataMapper> mapperlijn = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapperlijn->SetInput(linesPolyData);
#else
	mapperlijn->SetInputData(linesPolyData);
#endif

	vtkSmartPointer<vtkActor> actorlijn = vtkSmartPointer<vtkActor>::New();
	actorlijn->SetMapper(mapperlijn);
	actorlijn->GetProperty()->SetLineWidth(5);    // de dikte van de lijn aanpassen
	return actorlijn ;


}

vtkSmartPointer<vtkActor> visualization_vtk_c::geef_actor_punten(std::vector<wif_core::vector_2d_c> mypoints)
{


	vtkSmartPointer<vtkPoints> points =  vtkSmartPointer<vtkPoints>::New();

	unsigned char red[3] = {255, 0, 0};
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");


	int n = mypoints.size() ;

	for(int j = 0; j < n; j = j + 1)
	{
		points->InsertNextPoint(mypoints[j].x, mypoints[j].y, 0);    // beginpunt lijn
		colors->InsertNextTupleValue(red);
	}

	vtkSmartPointer<vtkPolyData> pointsPolydata = vtkSmartPointer<vtkPolyData>::New();
	pointsPolydata->SetPoints(points);


	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter = vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
	vertexFilter->SetInputConnection(pointsPolydata->GetProducerPort());
#else
	vertexFilter->SetInputData(pointsPolydata);
#endif
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->ShallowCopy(vertexFilter->GetOutput());
	polydata->GetPointData()->SetScalars(colors);


// VISUALISATIE
	vtkSmartPointer<vtkPolyDataMapper> mapperpunt = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapperpunt->SetInputConnection(polydata->GetProducerPort());
#else
	mapperpunt->SetInputData(polydata);
#endif

	vtkSmartPointer<vtkActor> actorpunt = vtkSmartPointer<vtkActor>::New();
	actorpunt->SetMapper(mapperpunt);
	actorpunt->GetProperty()->SetPointSize(5);

	return actorpunt;
}
void visualization_vtk_c::arrow_plot() const
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

	vtkSmartPointer<vtkRenderer> renderer =
	    vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(1, 1, 1);

	vtkSmartPointer<vtkPolyDataMapper> mapper =
	    vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(glyph3D->GetOutputPort());

	vtkSmartPointer<vtkActor> actor =
	    vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetRepresentationToSurface();

	vtkSmartPointer<vtkActor> superquadricActor =
	    vtkSmartPointer<vtkActor>::New();
	superquadricActor->SetMapper(mapper);
	vtkSmartPointer<vtkCubeAxesActor> cubeAxesActor = axis(glyph3D, renderer);

	cubeAxesActor->GetProperty()->SetColor(0, 0, 0);
	cubeAxesActor->SetXTitle("x");
	cubeAxesActor->SetYTitle("y");

	/*renderer->AddActor(actor);
	renderer->AddActor(cubeAxesActor );
	renderer->ResetCamera();*/

	renderer->AddActor(cubeAxesActor);
	renderer->AddActor(superquadricActor);
	renderer->ResetCamera();

	//renderer->GetActiveCamera()->Azimuth(0);
	//renderer->GetActiveCamera()->Elevation(0);

	vtkSmartPointer<vtkRenderWindow> renderWindow =
	    vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkInteractorStyleImage> imageStyle =
	    vtkSmartPointer<vtkInteractorStyleImage>::New();
	renderWindow->GetInteractor()->SetInteractorStyle(imageStyle);

	renderWindow->Render();
	//

	renderWindowInteractor->Start();
}

void visualization_vtk_c::streamlines_plot(vtkSmartPointer<vtkStructuredGrid> sgrid, uint32_t number_of_streamlines) const
{
	// Source of the streamlines

	vtkSmartPointer<vtkLineSource> seeds = vtkSmartPointer<vtkLineSource>::New();
	seeds->SetResolution(number_of_streamlines);
	seeds->SetPoint1(min_range.x, max_range.y, 0);
	seeds->SetPoint2(min_range.x, min_range.y, 0);

	// Setting Streamline properties
	vtkSmartPointer<vtkStreamLine> streamLine = vtkSmartPointer<vtkStreamLine>::New();

	streamLine->SetInput(sgrid);
	streamLine->SetSource(seeds->GetOutput());

	// Integration properties

	streamLine->SetMaximumPropagationTime(200);
	streamLine->SetIntegrationStepLength(.2);
	streamLine->SetStepLength(.001);
	streamLine->SetNumberOfThreads(1);
	streamLine->SetIntegrationDirectionToForward();
	//streamLine->VorticityOn();

	// Creating a Mapper and Actor

	vtkSmartPointer<vtkPolyDataMapper> streamLineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	streamLineMapper->SetInputConnection(streamLine->GetOutputPort());

	vtkSmartPointer<vtkActor> streamLineActor = vtkSmartPointer<vtkActor>::New();
	streamLineActor->SetMapper(streamLineMapper);
	streamLineActor->GetProperty()->SetColor(0, 0.2, 0.5);
	streamLineActor->GetProperty()->SetLineWidth(3);
	streamLineActor->VisibilityOn();

	//return streamLineActor;

	vtkSmartPointer<vtkRenderer> renderer24 = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer24);

	vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	interactor->SetRenderWindow(renderWindow);

	vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
	interactor->SetInteractorStyle(style);



	vtkSmartPointer<vtkCubeAxesActor> cubeAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
	cubeAxesActor->SetBounds(min_range.x, max_range.x, min_range.y, max_range.y, 0, 0);
	cubeAxesActor->SetCamera(renderer24->GetActiveCamera());

	cubeAxesActor->DrawXGridlinesOff();
	cubeAxesActor->DrawYGridlinesOff();
	cubeAxesActor->DrawZGridlinesOff();
	cubeAxesActor->ZAxisLabelVisibilityOff();
#if VTK_MAJOR_VERSION > 5
	cubeAxesActor->SetGridLineLocation(VTK_GRID_LINES_FURTHEST);
#endif

	cubeAxesActor->XAxisMinorTickVisibilityOff();
	cubeAxesActor->YAxisMinorTickVisibilityOff();
	cubeAxesActor->ZAxisMinorTickVisibilityOff();

	cubeAxesActor->GetProperty()->SetColor(0, 0, 0);
	cubeAxesActor->SetXTitle("x");
	cubeAxesActor->SetYTitle("y");
	cubeAxesActor->SetCamera(renderer24->GetActiveCamera());

	renderer24->AddActor(streamLineActor);
	renderer24->AddActor(cubeAxesActor);
	renderer24->ResetCamera();

	// Add the actors to the renderer, set the background and size
	renderer24->SetBackground(1, 1, 1);
	renderWindow->SetSize(300, 300);
	interactor->Initialize();
	std::cout << "test3" << std::endl;
	renderWindow->Render();

	interactor->Start();
}

vtkSmartPointer<vtkCubeAxesActor> visualization_vtk_c::axis(vtkSmartPointer<vtkPlaneSource> object, vtkSmartPointer<vtkRenderer> renderer) const
{
	vtkSmartPointer<vtkCubeAxesActor> cubeAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
	cubeAxesActor->SetBounds(object->GetOutput()->GetBounds());
	cubeAxesActor->SetCamera(renderer->GetActiveCamera());

	cubeAxesActor->DrawXGridlinesOff();
	cubeAxesActor->DrawYGridlinesOff();
	cubeAxesActor->DrawZGridlinesOff();
	cubeAxesActor->ZAxisLabelVisibilityOff();
#if VTK_MAJOR_VERSION > 5
	cubeAxesActor->SetGridLineLocation(VTK_GRID_LINES_FURTHEST);
#endif

	cubeAxesActor->XAxisMinorTickVisibilityOff();
	cubeAxesActor->YAxisMinorTickVisibilityOff();
	cubeAxesActor->ZAxisMinorTickVisibilityOff();
	return cubeAxesActor;
}

vtkSmartPointer<vtkCubeAxesActor> visualization_vtk_c::axis(vtkSmartPointer<vtkGlyph3D> object, vtkSmartPointer<vtkRenderer> renderer) const
{
	vtkSmartPointer<vtkCubeAxesActor> cubeAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
	cubeAxesActor->SetBounds(object->GetOutput()->GetBounds());
	cubeAxesActor->SetCamera(renderer->GetActiveCamera());

	cubeAxesActor->DrawXGridlinesOff();
	cubeAxesActor->DrawYGridlinesOff();
	cubeAxesActor->DrawZGridlinesOff();
	cubeAxesActor->ZAxisLabelVisibilityOff();
#if VTK_MAJOR_VERSION > 5
	cubeAxesActor->SetGridLineLocation(VTK_GRID_LINES_FURTHEST);
#endif

	cubeAxesActor->XAxisMinorTickVisibilityOff();
	cubeAxesActor->YAxisMinorTickVisibilityOff();
	cubeAxesActor->ZAxisMinorTickVisibilityOff();

	cubeAxesActor->GetProperty()->SetColor(0, 0, 0);
	cubeAxesActor->SetXTitle("x");
	cubeAxesActor->SetYTitle("y");
	cubeAxesActor->SetCamera(renderer->GetActiveCamera());

	return cubeAxesActor;
}

} // namespace wif_viz
