#include "vtkDoubleArray.h"
#include "vtkImageData.h"
#include "vtkPoints.h"
#include "vtkStructuredGrid.h"
#include "vtkStructuredGridWriter.h"
#include "vtkImageWriter.h"
#include "vtkPointData.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkColorTransferFunction.h"
#include "vtkWriter.h"
#include "vtkCamera.h"
#include "vtkArrowSource.h"

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

#include <vtkActor2D.h>
#include <vtkIntArray.h>
#include <vtkLabelPlacementMapper.h>
#include <vtkPointSetToLabelHierarchy.h>
#include <vtkPointSource.h>
#include <vtkStringArray.h>







#include<iostream>

vtkStructuredGrid * createGrid()
{

	int nx = 11;
	int ny = 11;
	int nz = 1;

	vtkPoints * points = vtkPoints::New();

	vtkDoubleArray * data1 = vtkDoubleArray::New();
	data1->SetNumberOfComponents(1);
	data1->SetName("myscalars");

	for(int iz = 0; iz < nz; iz++)
	{
		for(int iy = 0; iy < ny; iy++)
		{
			for(int ix = 0; ix < nx; ix++)
			{
				double x = ix - 5;
				double y = iy - 5;
				double f = x * x + y * y;
				//double vx = 2 * x ;
				//double vy = 2 * y ;

				points->InsertNextPoint(x, y, iz);
				data1->InsertNextValue(f);

			}
		}
	}

	vtkStructuredGrid * grid = vtkStructuredGrid::New();
	grid->SetDimensions(nx, ny, nz);
	grid->SetPoints(points);
	grid->GetPointData()->SetScalars(data1);

	points->Delete();
	data1->Delete();

	return grid;
}


int main()
{

// Maakt de punten
	double O1[3] = {0.0, 0.0, 0.0};   //oorsprong 1
	double Q1[3] = {4.0, -2.0, 0.0};	//punt 1
	double O2[3] = { -1.0, 3.0, 0.0};	//oorsprong 2
	double Q2[3] = { -3.0, 0.0, 0.0};	//punt 2

// Maak een vtk object om de punten in op te slagen
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	pts->InsertNextPoint(O1);
	pts->InsertNextPoint(Q1);
	pts->InsertNextPoint(O2);
	pts->InsertNextPoint(Q2);

// kleuren aanmaken
	//unsigned char red[3] = {255, 0, 0};
	//unsigned char green[3] = {0, 255, 0};
	unsigned char zwart[3] = {0, 0, 0};

// Setup the colors array
	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

// Add the colors we created to the colors array
	colors->InsertNextTupleValue(zwart);
	colors->InsertNextTupleValue(zwart);

// maak eerste lijn (tussen O1 en P1)
	vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
	line0->GetPointIds()->SetId(0, 0); //het tweede getal is de index van O1 in de vtkPoints
	line0->GetPointIds()->SetId(1, 1); //het tweede getal is de index van P1 in de vtkPoints

// Maak tweede lijn
	vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
	line1->GetPointIds()->SetId(0, 2); //het tweede getal is de index van O2 in de vtkPoints
	line1->GetPointIds()->SetId(1, 3); //het tweede getal is de index van P2 in de vtkPoints

// Create a cell array to store the lines in and add the lines to it
	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(line0);
	lines->InsertNextCell(line1);

// Create a polydata to store everything in
	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();

// Add the points to the dataset
	linesPolyData->SetPoints(pts);

// Add the lines to the dataset
	linesPolyData->SetLines(lines);

// Color the lines - associate the first component (red) of the
// colors array with the first component of the cell array (line 0)
// and the second component (green) of the colors array with the
// second component of the cell array (line 1)
	linesPolyData->GetCellData()->SetScalars(colors);




// PUNT
// definieer de plaats van het punt
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	const float P[3] = {1.0, 2.0, 0.0};

// Create the topology of the point (a vertex)
	vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType pid[1];
	pid[0] = points->InsertNextPoint(P);
	vertices->InsertNextCell(1, pid);

// Create a polydata object
	vtkSmartPointer<vtkPolyData> point = vtkSmartPointer<vtkPolyData>::New();

// Set the points and vertices we created as the geometry and topology of the polydata
	point->SetPoints(points);
	point->SetVerts(vertices);










// VISUALISEREN VELD

	vtkStructuredGrid * sgrid = createGrid();

// opslagen in vtk bestand
	vtkStructuredGridWriter * writer = vtkStructuredGridWriter::New();
	writer->SetFileName("output.vtk");
	writer->SetInput(sgrid);
	writer->Write();

// min en max kleuren instellen
	vtkColorTransferFunction * lut = vtkColorTransferFunction::New();
	double minvalue = 0;
	double maxvalue = 40;
	lut->AddRGBPoint(minvalue, 0, 0, 1);
	lut->AddRGBPoint(maxvalue, 1, 0, 0);



	vtkDataSetMapper * mapperveld = vtkDataSetMapper::New();
	mapperveld->SetInput(sgrid);
	mapperveld->SetLookupTable(lut);


	vtkActor * actorveld = vtkActor::New();
	actorveld->SetMapper(mapperveld);

	vtkRenderer * renderer = vtkRenderer::New();


	vtkRenderWindow * renWin = vtkRenderWindow::New();
	renWin->AddRenderer(renderer);
	renWin->SetSize(1000 , 1000);

	vtkInteractorStyleTrackballCamera * style = vtkInteractorStyleTrackballCamera::New();

	vtkRenderWindowInteractor * iren = vtkRenderWindowInteractor::New();

	iren->SetInteractorStyle(style);
	iren->SetRenderWindow(renWin);

	renderer->AddActor(actorveld);
	renderer->SetBackground(0  , 0 , 0);     // stelt de kleur van de achtergrond in
	renderer->ResetCamera();



//  VISUALISEREN LIJN
	vtkSmartPointer<vtkPolyDataMapper> mapperlijn = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapperlijn->SetInput(linesPolyData);
#else
	mapperlijn->SetInputData(linesPolyData);
#endif

	vtkSmartPointer<vtkActor> actorlijn = vtkSmartPointer<vtkActor>::New();
	actorlijn->SetMapper(mapperlijn);
	actorlijn->GetProperty()->SetLineWidth(5);    // de dikte van de lijn aanpassen


	renWin->AddRenderer(renderer);
	iren->SetRenderWindow(renWin);
	renderer->AddActor(actorlijn);


	// VISUALISATIE PUNT
	vtkSmartPointer<vtkPolyDataMapper> mapperpunt = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapperpunt->SetInput(point);
#else
	mapperpunt->SetInputData(point);
#endif

	vtkSmartPointer<vtkActor> actorpunt = vtkSmartPointer<vtkActor>::New();
	actorpunt->SetMapper(mapperpunt);
	actorpunt->GetProperty()->SetPointSize(10);


	renderer->AddActor(actorpunt);











	renWin->Render();
	iren->Initialize();
	iren->Start();

	iren->Delete();
	style->Delete();
	renWin->Delete();
	renderer->Delete();
	actorveld->Delete();
	actorlijn->Delete();
	mapperveld->Delete();
	mapperlijn->Delete();
	lut->Delete();
	writer->Delete();
	sgrid->Delete();


	return 0;
}
