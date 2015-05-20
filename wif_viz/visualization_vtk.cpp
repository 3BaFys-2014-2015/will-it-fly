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
#include <vtkGeometryFilter.h>
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
	visualization_c(flow, min_range, max_range)
{
	//
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

vtkSmartPointer<vtkCleanPolyData> clip_into_pieces(vtkSmartPointer<vtkStructuredGrid> grid, const std::vector<double> & levels)
{
	// This is my last resort

	vtkSmartPointer<vtkGeometryFilter> filter = vtkSmartPointer<vtkGeometryFilter>::New();

	filter->SetInput(grid);
	filter->Update();

	vtkSmartPointer<vtkPolyData> data = filter->GetOutput();

	//

	std::vector<double> lvls = levels;

	std::sort(lvls.begin(), lvls.end());

	vtkSmartPointer<vtkAppendPolyData> append_poly_data = vtkSmartPointer<vtkAppendPolyData>::New();

	std::vector<vtkSmartPointer<vtkClipPolyData>> clippers_low;
	std::vector<vtkSmartPointer<vtkClipPolyData>> clippers_high;

	for(int i = 0; i < (levels.size() - 1); i++)
	{
		const double low = lvls[i];
		const double high = lvls[i + 1];
		clippers_low.push_back(vtkSmartPointer<vtkClipPolyData>::New());
		clippers_low[i]->SetValue(low);

		if(i == 0)
		{
			clippers_low[i]->SetInput(data);
		}
		else
		{
			clippers_low[i]->SetInput(clippers_high[i - 1]->GetOutput(1));
		}

		clippers_low[i]->InsideOutOff();
		clippers_low[i]->Update();

		clippers_high.push_back(vtkSmartPointer<vtkClipPolyData>::New());
		clippers_high[i]->SetValue(high);
		clippers_high[i]->SetInput(clippers_low[i]->GetOutput());
		clippers_high[i]->GenerateClippedOutputOn();
		clippers_high[i]->InsideOutOn();
		clippers_high[i]->Update();

		if(clippers_high[i]->GetOutput()->GetNumberOfCells() == 0)
		{
			continue;
		}

		vtkSmartPointer<vtkDoubleArray> midvalue = vtkSmartPointer<vtkDoubleArray>::New();
		midvalue->SetNumberOfComponents(1);
		midvalue->SetNumberOfTuples(clippers_high[i]->GetOutput()->GetNumberOfCells());
		midvalue->FillComponent(0, (low + high) * 0.5);

		clippers_high[i]->GetOutput()->GetCellData()->SetScalars(midvalue);
		append_poly_data->AddInput(clippers_high[i]->GetOutput());
	}

	append_poly_data->Update();

	vtkSmartPointer<vtkCleanPolyData> pieces = vtkSmartPointer<vtkCleanPolyData>::New();
	pieces->SetInput(append_poly_data->GetOutput());
	pieces->Update();

	return pieces;
}
/*
void create_contour_actors(vtkSmartPointer<vtkStructuredGrid> grid, const std::vector<double> & levels, vtkSmartPointer<vtkActor> contour_levels, vtkSmartPointer<vtkActor> contour_lines)
{
	vtkSmartPointer<vtkCleanPolyData> clips = clip_into_pieces(grid, levels);

	if(contour_levels)
	{
		vtkSmartPointer<vtkLookupTable> lut2 = vtkSmartPointer<vtkLookupTable>::New();
		lut2->SetNumberOfTableValues(levels.size() + 1);
		lut2->Build();

		vtkSmartPointer<vtkPolyDataMapper> contourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		contourMapper->SetInput(clips->GetOutput());
		contourMapper->SetScalarRange(levels.front(), levels.back());
		contourMapper->SetScalarModeToUseCellData();
		contourMapper->SetLookupTable(lut2);
		contourMapper->Update();

		//schaal bar
		vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
		scalarBar->SetLookupTable(lut2);

		contour_levels->SetMapper(contourMapper);
		contour_levels->GetProperty()->SetInterpolationToFlat();
	}


	if(contour_lines)
	{
		vtkSmartPointer<vtkContourFilter> contours3 = vtkSmartPointer<vtkContourFilter>::New();
		contours3->SetInput(clips->GetOutput());
		contours3->GenerateValues(levels.size(), levels.front(), levels.back());

		vtkSmartPointer<vtkPolyDataMapper> contourLineMapperer = vtkSmartPointer<vtkPolyDataMapper>::New();
		contourLineMapperer->SetInput(contours3->GetOutput());
		contourLineMapperer->SetScalarRange(levels.front(), levels.back());
		contourLineMapperer->ScalarVisibilityOff();

		contour_lines->SetMapper(contourLineMapperer);
		contour_lines->GetProperty()->SetLineWidth(2);
	}
}

*/

void data_range_two(double & min_v, double & max_v, vtkSmartPointer<vtkStructuredGrid> grid)
{
	double values[2];

	grid->GetPointData()->GetScalars()->GetRange(values);

	min_v = values[0];
	max_v = values[1];
}

void fix_contour_range(double min_v, double max_v, std::vector<double> & contours)
{
	if(contours.size() == 0)
	{
		contours.resize(10, 0);
	}

	double min_c = *std::min_element(contours.begin(), contours.end());
	double max_c = *std::max_element(contours.begin(), contours.end());

	if((min_c < min_v) || (max_c > max_v) || (min_c == max_c))
	{
		std::cout << "Rescaling contours" << std::endl;

		double delta = (max_v - min_v) / ((double)contours.size() - 1);

		for(uint32_t i = 0; i < contours.size(); i++)
		{
			contours[i] = min_v + i * delta;
		}
	}

	std::sort(contours.begin(), contours.end());
}

void visualize_scalar_field(const scalar_field_c & field, vtkSmartPointer<vtkStructuredGrid> grid, vtkSmartPointer<vtkRenderer> renderer)
{
	if((field.style & ESDS_GRADIENT) && (field.style & ESDS_DISCRETE))
	{
		std::cout << "Can't draw both a gradient scalar field, and a discrete one." << std::endl;
		return;
	}

	double min_v = 0.0;
	double max_v = 0.0;

	data_range_two(min_v, max_v, grid);

	if(field.style & ESDS_GRADIENT)
	{
		vtkSmartPointer<vtkDataSetMapper> data_set_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
		data_set_mapper->SetInput(grid);
		data_set_mapper->SetScalarRange(min_v, max_v);
		data_set_mapper->InterpolateScalarsBeforeMappingOn();

        vtkSmartPointer<vtkScalarBarActor> scalar_bar = vtkSmartPointer<vtkScalarBarActor>::New();
        scalar_bar->SetLookupTable(data_set_mapper->GetLookupTable());
        renderer->AddActor(scalar_bar);

		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(data_set_mapper);
		actor->GetProperty()->SetInterpolationToFlat();

		renderer->AddActor(actor);
	}

	if((field.style & ESDS_DISCRETE) || (field.style & ESDS_CONTOURS))
	{
		std::vector<double> corrected_contours = field.contour_locations;

		fix_contour_range(min_v, max_v, corrected_contours);

		vtkSmartPointer<vtkCleanPolyData> filled_contours = clip_into_pieces(grid, corrected_contours);

		if(field.style & ESDS_DISCRETE)
		{
			vtkSmartPointer<vtkLookupTable> lut2 = vtkSmartPointer<vtkLookupTable>::New();
			lut2->SetNumberOfTableValues(corrected_contours.size() + 1);
			lut2->Build();

			vtkSmartPointer<vtkPolyDataMapper> contourMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			contourMapper->SetInput(filled_contours->GetOutput());
			contourMapper->SetScalarRange(min_v, max_v);
			contourMapper->SetScalarModeToUseCellData();
			contourMapper->SetLookupTable(lut2);
			contourMapper->Update();

			vtkSmartPointer<vtkScalarBarActor> scalar_bar = vtkSmartPointer<vtkScalarBarActor>::New();
			scalar_bar->SetLookupTable(lut2);
			renderer->AddActor(scalar_bar);

			vtkSmartPointer<vtkActor> contourActor = vtkSmartPointer<vtkActor>::New();
			contourActor->SetMapper(contourMapper);
			contourActor->GetProperty()->SetInterpolationToFlat();
			renderer->AddActor(contourActor);
		}

		if(field.style & ESDS_CONTOURS)
		{
			vtkSmartPointer<vtkContourFilter> contour_filter = vtkSmartPointer<vtkContourFilter>::New();
			contour_filter->SetInput(filled_contours->GetOutput());
			contour_filter->GenerateValues(corrected_contours.size(), min_v, max_v);

			vtkSmartPointer<vtkPolyDataMapper> poly_data_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			poly_data_mapper->SetInput(contour_filter->GetOutput());
			poly_data_mapper->SetScalarRange(min_v, max_v);
			poly_data_mapper->ScalarVisibilityOff();

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(poly_data_mapper);
			actor->GetProperty()->SetLineWidth(5);

			renderer->AddActor(actor);
		}
	}
}

void visualize_vector_field(const vector_field_c & field, vtkSmartPointer<vtkStructuredGrid> grid, vtkSmartPointer<vtkRenderer> renderer)
{
	if(field.style & EVDS_ARROWS)
	{
		vtkSmartPointer<vtkArrowSource> arrow_source = vtkSmartPointer<vtkArrowSource>::New();

		vtkSmartPointer<vtkGlyph3D> glyph_3d = vtkSmartPointer<vtkGlyph3D>::New();
		glyph_3d->SetSourceConnection(arrow_source->GetOutputPort());
#if VTK_MAJOR_VERSION <= 5
		glyph_3d->SetInput(grid);
#else
		glyph_3d->SetInputData(grid);
#endif
		glyph_3d->SetColorMode(2);
		glyph_3d->SetScaleFactor(field.arrow_scale);
		glyph_3d->Update();

		vtkSmartPointer<vtkPolyDataMapper> poly_data_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		poly_data_mapper->SetInputConnection(glyph_3d->GetOutputPort());

		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		actor->SetMapper(poly_data_mapper);
		actor->GetProperty()->SetRepresentationToSurface();

		renderer->AddActor(actor);
	}

	if(field.style & EVDS_STREAMLINES)
	{
		if(field.streamline_resolution == 0)
		{
			std::cout << "Streamline resolution is zero, nothing will be drawn." << std::endl;
		}
		else
		{
			vtkSmartPointer<vtkLineSource> seeds = vtkSmartPointer<vtkLineSource>::New();
			seeds->SetResolution(field.streamline_resolution);

			if(field.streamline_seeds.begin == field.streamline_seeds.end)
			{
				seeds->SetPoint1(field.streamline_seeds.begin.x, field.streamline_seeds.begin.y, 0);
				seeds->SetPoint2(field.streamline_seeds.end.x, field.streamline_seeds.end.y, 0);
			}
			else
			{
				seeds->SetPoint1(field.streamline_seeds.begin.x, field.streamline_seeds.begin.y, 0);
				seeds->SetPoint2(field.streamline_seeds.end.x, field.streamline_seeds.end.y, 0);
			}

			vtkSmartPointer<vtkStreamLine> stream_line = vtkSmartPointer<vtkStreamLine>::New();
			stream_line->SetInput(grid);
			stream_line->SetSource(seeds->GetOutput());
			stream_line->SetMaximumPropagationTime(20);
			stream_line->SetIntegrationStepLength(.2);
			stream_line->SetStepLength(.001);
			stream_line->SetNumberOfThreads(1);
			stream_line->SetIntegrationDirectionToForward();

			vtkSmartPointer<vtkPolyDataMapper> poly_data_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			poly_data_mapper->SetInputConnection(stream_line->GetOutputPort());

			vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
			actor->SetMapper(poly_data_mapper);
			actor->GetProperty()->SetColor(0, 0.2, 0.5);
			actor->GetProperty()->SetLineWidth(3);
			actor->VisibilityOn();

			renderer->AddActor(actor);
		}
	}
}

void visualization_vtk_c::draw(const std::string & name)
{
	if(this->global_settings.output_to_file)
	{
		if(one_is_zero(this->psi_field.bins))
		{
			write_to_file(construct_psi_grid(), name + "-psi.vtk");
		}

		if(one_is_zero(this->phi_field.bins))
		{
			write_to_file(construct_phi_grid(), name + "-phi.vtk");
		}

		if(one_is_zero(this->v_field.bins))
		{
			write_to_file(construct_velocity_grid(), name + "-velocity.vtk");
		}
	}
	else
	{
		vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
		vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();

		renWin->AddRenderer(renderer);
		renWin->SetSize(1000, 1000);
		renWin->SetWindowName(name.c_str());

		renderer->SetBackground(1, 1, 1);
		renderer->ResetCamera();

		if(this->global_settings.draw_scale)
		{
			vtkSmartPointer<vtkCubeAxesActor> axes = vtkSmartPointer<vtkCubeAxesActor>::New();
			axes->SetBounds(min_range.x, max_range.x, min_range.y, max_range.y, 0, 0);

			axes->DrawXGridlinesOff();
			axes->DrawYGridlinesOff();
			axes->DrawZGridlinesOff();
			axes->ZAxisLabelVisibilityOff();
#if VTK_MAJOR_VERSION > 5
			axes->SetGridLineLocation(VTK_GRID_LINES_FURTHEST);
#endif
			axes->XAxisMinorTickVisibilityOff();
			axes->YAxisMinorTickVisibilityOff();
			axes->ZAxisMinorTickVisibilityOff();

			axes->GetProperty()->SetColor(0, 0, 0);

			/// VTK wilt dit niet verschuiven, zodat de cijfers op de assen niet overlappen.
			axes->SetXTitle("");
			axes->SetYTitle("");
			axes->SetCamera(renderer->GetActiveCamera());

			renderer->AddActor(axes);
		}

		if(one_is_zero(this->psi_field.bins))
		{
			visualize_scalar_field(this->psi_field, construct_psi_grid(), renderer);
		}

		if(one_is_zero(this->phi_field.bins))
		{
			visualize_scalar_field(this->phi_field, construct_phi_grid(), renderer);
		}

		if(one_is_zero(this->v_field.bins))
		{
			visualize_vector_field(this->v_field, construct_velocity_grid(), renderer);
		}

		if(this->airfoil)
		{
			renderer->AddActor(geef_actor_lijnen(this->airfoil->get_lines()));
			renderer->AddActor(geef_actor_punten(this->airfoil->get_points()));
		}

		renderer->ResetCamera();

		vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
		vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();

		iren->SetInteractorStyle(style);
		iren->SetRenderWindow(renWin);

		iren->Initialize();
		renWin->Render();
		iren->Start();

		iren->TerminateApp();
	}
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
	vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();

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
	vtkSmartPointer<vtkPoints> points = construct_points(this->psi_field.bins);
	vtkSmartPointer<vtkDoubleArray> field = construct_field(this->psi_field.bins, true);

	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);

		field->InsertNextValue(flow->get_psi(pos));
	}

	return combine_grid(this->psi_field.bins, points, field);
}

vtkSmartPointer<vtkStructuredGrid> visualization_vtk_c::construct_phi_grid() const
{
	vtkSmartPointer<vtkPoints> points = construct_points(this->phi_field.bins);
	vtkSmartPointer<vtkDoubleArray> field = construct_field(this->phi_field.bins, true);

	for(int i = 0; i < points->GetNumberOfPoints(); i++)
	{
		double x[3];

		points->GetPoint(i, x);

		const vector_2d_c pos(x[0], x[1]);

		field->InsertNextValue(flow->get_phi(pos))
		;
	}

	return combine_grid(this->phi_field.bins, points, field);
}

vtkSmartPointer<vtkStructuredGrid> visualization_vtk_c::construct_velocity_grid() const
{
	vtkSmartPointer<vtkPoints> points = construct_points(this->v_field.bins);
	vtkSmartPointer<vtkDoubleArray> field = construct_field(this->v_field.bins, false);

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

	return combine_grid(this->v_field.bins, points, field);
}

vtkSmartPointer<vtkActor> visualization_vtk_c::separating_streamlines(vtkSmartPointer<vtkPlaneSource> plane) const
{
	vtkSmartPointer<vtkContourFilter> contours = vtkSmartPointer<vtkContourFilter>::New();
	contours->SetInputConnection(plane->GetOutputPort());
	vector<double> values;

	for(int i = 0; i < stagnation_point.size(); ++i)
	{
		values[i] = flow->get_phi(stagnation_point[i]);
	}

	std::sort(values.begin(), values.end());

	for(int i = 0; i < stagnation_point.size(); ++i)
	{
		contours->SetValue(i, values[i]);
	}

	vtkSmartPointer<vtkPolyDataMapper> contourLineMapperer = vtkSmartPointer<vtkPolyDataMapper>::New();
	contourLineMapperer->SetInputConnection(contours->GetOutputPort());
	contourLineMapperer->SetScalarRange(values[0], values[values.size() - 1]);
	contourLineMapperer->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> contourLineActor = vtkSmartPointer<vtkActor>::New();
	contourLineActor->SetMapper(contourLineMapperer);
	contourLineActor->GetProperty()->SetLineWidth(2);

	return contourLineActor;

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

} // namespace wif_viz
