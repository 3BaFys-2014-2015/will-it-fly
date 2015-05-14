#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkRenderer.h>
#include <vtkColorTransferFunction.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkXMLPolyDataReader.h>
#include <vtkContourFilter.h>
#include <vtkStripper.h>
#include <vtkPolyDataMapper.h>
#include <vtkLabeledDataMapper.h>
#include <vtkPlaneSource.h>
#include <vtkLookupTable.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkActor2D.h>

#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>

#include <vtkMath.h>

#include <vector>
#include <algorithm>

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

void print_lut(vtkSmartPointer<vtkColorTransferFunction> lut, double min, double max, uint32_t n)
{
	double delta = (max - min) / n;

	for(uint32_t i = 0; i <= n; i++)
	{
		const double v = min + delta * i;

		std::cout << v << " " << lut->GetColor(v)[0] << std::endl;
	}
}

int main()
{
	double clip_max = 1.0;
	double clip_min = -1.0;

	vtkSmartPointer<vtkColorTransferFunction> lut = vtkColorTransferFunction::New();

	//double delta = clip_max - clip_min;

	lut->AddRGBPoint(std::nextafter(clip_min, -std::numeric_limits<double>::infinity()), 0, 0, 1);
	lut->SetNanColor(0, 0, 1);

	lut->AddRGBPoint(clip_min, 0.4, 0, 0);
	lut->AddRGBPoint(clip_max, 0.5, 0, 0);

	lut->AddRGBPoint(std::nextafter(clip_max, std::numeric_limits<double>::infinity()), 0, 0, 1);

	print_lut(lut, -1.1, -1.0, 10);
	return 0;
}
