#include <stdio.h>
#include <cmath>
#include <ctime>
#include <gsl/gsl_integration.h>

using namespace std;

double f(double x, void * params)
{
	double alpha = *(double *) params;
	double f = log(alpha * x) / sqrt(x);
	return f;
}

struct my_ff_params
{
	double a;
	double b;
	double c;
};

double ff(double x, void * p)
{
	struct my_ff_params * params = (struct my_ff_params *)p;
	double alpha = (params->a);
	double beta = (params->b);
	double gamma = (params->c);
	double ff = alpha * pow(x, 2) + beta * x + gamma;
	return ff;
}

struct my_velocity_params
{
	double d;
	double e;
	double f;
};

double velocityx(double XI, void * p)
{
	struct my_velocity_params * params = (struct my_velocity_params *)p;
	double sigma = (params->d);
	double x = (params->e);
	double y = (params->f);
	double velocityx = sigma / 2 / 3.14 * x / (pow(x, 2) + pow(y - XI, 2));
	return velocityx;
}

double velocityy(double XI, void * p)
{
	struct my_velocity_params * params = (struct my_velocity_params *)p;
	double sigma = (params->d);
	double x = (params->e);
	double y = (params->f);
	double velocityy = sigma / 2 / 3.14 * (y - XI) / (pow(x, 2) + pow(y - XI, 2));
	return velocityy;
}

int main()
{

	clock_t start;
	double duration;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);  //Maakt plaats vrij voor intervallen

	double result, error;
	double alpha = 1.0;
	double int_bound_min = 0;
	double int_bound_max = 1;

	gsl_function F;
	F.function = &f;
	F.params = &alpha;

	start = clock();

	for(int i = 0; i < 1000; i++)
	{
		gsl_integration_qags(&F, int_bound_min, int_bound_max, 0, 1e-7, 1000, w, &result, &error);
	}

	duration = (clock() - start) / (double) CLOCKS_PER_SEC;

	printf("result          = % .18f\n", result);
	printf("estimated error = % .18f\n", error);
	printf("intervals =  %zd\n", w->size);
	printf("duration = %.18f\n", duration);

	gsl_integration_workspace_free(w);

	gsl_integration_workspace * ww = gsl_integration_workspace_alloc(1000);

	struct my_ff_params alpha2 = {1., 2., 3.};
	gsl_function FF;
	FF.function = &ff;
	FF.params = &alpha2;

	start = clock();

	for(int i = 0; i < 1000; i++)
	{
		gsl_integration_qags(&FF, int_bound_min, int_bound_max, 0, 1e-7, 1000, w, &result, &error);
	}

	duration = (clock() - start) / (double) CLOCKS_PER_SEC;

	printf("result          = % .18f\n", result);
	printf("estimated error = % .18f\n", error);
	printf("intervals =  %zd\n", ww->size);
	printf("duration = %.18f\n", duration);

	gsl_integration_workspace_free(ww);

	gsl_integration_workspace * w3 = gsl_integration_workspace_alloc(1000);

	struct my_velocity_params alpha3 = {1., 1., 1.};
	gsl_function VELOCITYX;
	VELOCITYX.function = &velocityx;
	VELOCITYX.params = &alpha3;

	start = clock();

	for(int i = 0; i < 1000; i++)
	{
		gsl_integration_qags(&VELOCITYX, int_bound_min, int_bound_max, 0, 1e-7, 1000, w3, &result, &error);
	}

	duration = (clock() - start) / (double) CLOCKS_PER_SEC;

	double sigma = (alpha3.d);
	double x = (alpha3.e);
	double y = (alpha3.f);
	double Act_X = sigma * x / 2 / 3.14 * log(((pow(x, 2) + (y - int_bound_min)) / ((pow(x, 2) + (y - int_bound_max)))));

	printf("result velocity x = % .18f\n", result);
	printf("estimated error = % .18f\n", error);
	printf("actual error =%.18f\n", result - Act_X);
	printf("intervals =  %zd\n", w3->size);
	printf("duration = %.18f\n", duration);


	gsl_function VELOCITYY;
	VELOCITYY.function = &velocityy;
	VELOCITYY.params = &alpha3;

	start = clock();

	for(int i = 0; i < 1000; i++)
	{
		gsl_integration_qags(&VELOCITYY, int_bound_min, int_bound_max, 0, 1e-7, 1000, w3, &result, &error);
	}

	duration = (clock() - start) / (double) CLOCKS_PER_SEC;

	double Act_Y = sigma / 3.14 * log((pow(x, 2) + pow(x - int_bound_min, 2)) / ((pow(x, 2) + pow(y - int_bound_max, 2))));

	printf("result velocity y = % .18f\n", result);
	printf("estimated error = % .18f\n", error);
	printf("actual error =%.18f\n", result - Act_Y);
	printf("intervals =  %zd\n", w3->size);
	printf("duration = %.18f\n", duration);

	gsl_integration_workspace_free(w3);


	return 0;
}
