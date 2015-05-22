#include "airfoil_c.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <assert.h>

#include <iterator>
#include <sstream>



namespace wif_core
{


airfoil_c::airfoil_c() :
	name("")
{
	//
}

bool is_empty_line(const std::string & line)
{
	return line.find_first_not_of(" \t\n\r") == std::string::npos;
}

bool contains_non_digits(const std::string & line)
{
	return line.find_first_not_of("0123456789.- \t\n\r") != std::string::npos;
}

bool string_to_double_doubles(const std::string & line, double & a, double & b)
{
	std::istringstream is(line);

	if(is >> a)
	{
		if(is >> b)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}
}

airfoil_c::airfoil_c(const std::string & filename) :
	name("")
{
	std::ifstream file(filename);

	if(!file.is_open())
	{
		std::cout << "Warning: file " << filename << "Couldn't be loaded, using this airfoil might give segmentation faults." << std::endl;
		return;
	}

	std::string line;
	bool ok = false;

	/// Skip whitespace before the header.
	while(!file.eof())
	{
		std::getline(file, line);

		if(!is_empty_line(line))
		{
			ok = true;
			break;
		}
	}

	if(!ok)
	{
		std::cout << "ERROR READING AIRFOIL: FILE CLOSED TOO SOON." << std::endl;
		return;
	}

	ok = false;
	name.append(line);

	while(!file.eof())
	{
		std::getline(file, line);

		// Ignore empty lines.
		if(is_empty_line(line))
		{
			std::cout << "SKIPPING WHITESPACE IN HEADER: THIS IS WEIRD, BUT OK." << std::endl;
			continue;
		}

		// Add all text to the name.
		if(contains_non_digits(line))
		{
			name.append(line);
		}
		else
		{
			ok = true;
			break;
		}
	}

	if(!ok)
	{
		std::cout << "ERROR READING AIRFOIL: FILE CLOSED TOO SOON." << std::endl;
		return;
	}

	ok = false;
	double x = 0.0;
	double y = 0.0;

	if(!string_to_double_doubles(line, x, y))
	{
		std::cout << "ERROR READING AIRFOIL: COULDN'T EXTRACT LINE." << std::endl;
		return;
	}

	const double lednicer_cutoff = 1.2;

	if((x < lednicer_cutoff) && (y < lednicer_cutoff))
	{
		/// Selig

		this->points.emplace_back(x, y);

		std::cout << x << y << std::endl;

		while(!file.eof())
		{
			std::cout << line << std::endl;
			std::getline(file, line);

			if(is_empty_line(line))
			{
				std::cout << "SKIPPING WHITESPACE IN SELIGER FORMAT: THIS IS WEIRD, BUT OK. IT MIGHT BE A TRAILING NEWLINE." << std::endl;
				std::cout << line << std::endl;
				continue;
			}
			else if(contains_non_digits(line))
			{
				std::cout << "ERROR READING AIRFOIL: TEXT IN NUMBER ONLY AREA." << std::endl;
				return;
			}
			else if(string_to_double_doubles(line, x, y))
			{
				this->points.emplace_back(x, y);
			}
			else
			{
				std::cout << "ERROR READING AIRFOIL: COULDN'T EXTRACT DOUBLES FROM LINE." << std::endl;
				return;
			}
		}
	}
	else
	{
		/// LEDNICER

		bool topside = true;
		bool started_reading = false;

		while(!file.eof())
		{
			std::getline(file, line);

			if(is_empty_line(line))
			{
				if(topside)
				{
					if(started_reading)
					{
						if(this->points.size() == 0)
						{
							std::cout << "ERROR: WEIRD AIRFOIL, NO POINTS ON TOP. THIS SHOULD NEVER HAPPEN." << std::endl;
							return;
						}

						topside = false;

						std::reverse(this->points.begin(), this->points.end());
						continue;
					}
					else
					{
						continue;
					}
				}
				else
				{
					std::cout << "WARNING: MORE THAN ONE EMPTY LINE IN LEDNICER FORMAT IN DATA PART. THIS MIGHT BE A TRAILING NEWSPACE AND NOT AN ISSUE." << std::endl;
				}
			}
			else if(contains_non_digits(line))
			{
				std::cout << "ERROR READING AIRFOIL: TEXT IN NUMBER ONLY AREA." << std::endl;
				return;
			}
			else if(string_to_double_doubles(line, x, y))
			{
				if(!started_reading)
				{
					started_reading = true;
				}

				if(topside)
				{
					this->points.emplace_back(x, y);
				}
				else
				{
					const vector_2d_c p(x, y);

					if((p - this->points.back()).get_length_sq() > (0.0001))
					{
						this->points.push_back(p);
					}
				}
			}
			else
			{
				std::cout << "ERROR READING AIRFOIL: COULDN'T EXTRACT DOUBLES FROM LINE." << std::endl;
				return;
			}
		}
	}

	//
	/* DO NOT USE
		std::ifstream detect(filename);

		if(!detect.is_open())
		{
			return; //just give up if file does not open
		}

		std::string line1;
		std::string data_pit;
		std::getline(detect, line1);
		double testval;
		detect >> testval;
		//std::cout << "testvalue :" << testval << std::endl;
		detect.close();
		std::ifstream data(filename);

		if(!data.good())
		{
			std::cout << "Could not open file." << filename << std::endl;
			return; //just give up if file does not open
		}


		if(testval > 1)
		{

			std::getline(data, this->name);
			std::getline(data, data_pit);
			double x;
			double y;
			data >> x >> y;


			while(data.good())
			{
				double x;
				double y;
				data >> x >> y;
				this->points.emplace_back(x, y);
				//std::cout << x << "\t" << y << std::endl;
			}

			unsigned int half_size = this->points.size() / 2;
			std::reverse(this->points.begin(), this->points.end() - half_size - 1);
		}
		else
		{

			//selig format
			std::getline(data, this->name);

			while(data.good())
			{
				double x;
				double y;
				data >> x >> y;
				this->points.emplace_back(x, y);
				//std::cout << x << "\t" << y << std::endl;
			}
		}

		this->points.pop_back(); //reads last line double. Not any more
	*/
}

airfoil_c::airfoil_c(const std::vector<vector_2d_c> & points, const std::string & name) :
	name(name),
	points(points)
{
	//
}


airfoil_c::airfoil_c(const vector_2d_c & midpoint, double radius, unsigned int corners) :
	name("circle"),
	points(corners + 1, vector_2d_c(0.0, 0.0))
{
	for(unsigned int i = 0; i <= corners; i++)
	{
		points[i] = (vector_2d_radian(radius, (2.0 * M_PI * i) / corners) + midpoint);
	}
}

bool airfoil_c::check_lengths() const
{
	if(!is_valid())
	{
		return false;
	}

	std::vector<line_2d_c> lines = get_lines();

	for(const line_2d_c & line : lines)
	{
		if(line.get_length() == 0.0)
		{
			return false;
		}
	}

	return true;
}

std::vector<line_2d_c> airfoil_c::get_lines() const
{
	if(!is_valid())
	{
		std::cout << "WARNING: TRYING TO GET LINES FROM AN INVALID AIRFOIL" << std::endl;
	}

	std::vector<line_2d_c> ret;
	ret.resize(this->points.size() - 1, line_2d_c(0.0, 0.0, 0.0, 0.0));

	for(unsigned int index = 0; index < this->points.size() - 1; index++)
	{

		line_2d_c r(this->points[index], this->points[index + 1]);
		ret[index] = r;

	}

	return ret;
}

std::vector<line_2d_c> airfoil_c::get_upper_panels() const
{
	std::vector<line_2d_c> lines = get_lines();
	std::vector<line_2d_c> upper;

	for(const line_2d_c & line : lines)
	{
		if(line.get_difference().x < 0.0)
		{
			upper.push_back(line);
		}
		else
		{
			break;
		}
	}

	std::sort(upper.begin(), upper.end(), [](const line_2d_c & lhs, const line_2d_c & rhs)
	{
		return lhs.get_center_point().x < rhs.get_center_point().x;
	});

	return upper;
}

std::vector<double> airfoil_c::get_upper_x() const
{
	std::vector<line_2d_c> upper = get_upper_panels();
	std::vector<double> x;

	for(const line_2d_c & line : upper)
	{
		x.push_back(line.get_center_point().x);
	}

	return x;
}

std::vector<line_2d_c> airfoil_c::get_lower_panels() const
{
	std::vector<line_2d_c> lines = get_lines();
	std::vector<line_2d_c> lower;

	std::reverse(lines.begin(), lines.end());

	for(const line_2d_c & line : lines)
	{
		if(line.get_difference().x > 0.0)
		{
			lower.push_back(line);
		}
		else
		{
			break;
		}
	}

	std::sort(lower.begin(), lower.end(), [](const line_2d_c & lhs, const line_2d_c & rhs)
	{
		return lhs.get_center_point().x < rhs.get_center_point().x;
	});

	return lower;
}

std::vector<double> airfoil_c::get_lower_x() const
{
	std::vector<line_2d_c> lower = get_lower_panels();
	std::vector<double> x;

	for(const line_2d_c & line : lower)
	{
		x.push_back(line.get_center_point().x);
	}

	return x;
}

std::vector<double> airfoil_c::select_upper_data(const std::vector<double> & input) const
{
	std::vector<line_2d_c> upper_panels = get_upper_panels();

	std::vector<double> output(input.begin(), input.begin() + upper_panels.size());

	std::reverse(output.begin(), output.end());

	return output;
}

std::vector<double> airfoil_c::select_lower_data(const std::vector<double> & input) const
{
	std::vector<line_2d_c> lower_panels = get_lower_panels();
	std::vector<double> reverse_input = input;

	std::reverse(reverse_input.begin(), reverse_input.end());

	std::vector<double> output(reverse_input.begin(), reverse_input.begin() + lower_panels.size());

	std::reverse(output.begin(), output.end());

	return output;
}

const std::vector<vector_2d_c> & airfoil_c::get_points() const
{
	if(!is_valid())
	{
		std::cout << "WARNING: TRYING TO GET POINTS FROM AN INVALID AIRFOIL" << std::endl;
	}

	return points;
}


std::vector<line_2d_c> airfoil_c::get_lines_reversed() const
{
	std::vector<line_2d_c> ret = this->get_lines();
	std::reverse(ret.begin(), ret.end());
	return ret;
}


vector_2d_c airfoil_c::get_intersection_first(const line_2d_c & line) const
{
	for(const line_2d_c & l : this->get_lines())
	{
		vector_2d_c intersect(0, 0);
		E_INTERSECTION intersect_type = line.get_intersection(l, intersect);

		if(intersect_type == EI_SEGMENT)
		{
			return(intersect);
		}
	}

	return vector_2d_c(INFINITY, INFINITY);
}


vector_2d_c airfoil_c::get_intersection_last(const line_2d_c & line) const
{
	for(const line_2d_c & l : this->get_lines_reversed())
	{
		vector_2d_c intersect(0, 0);
		E_INTERSECTION intersect_type = line.get_intersection(l, intersect);

		if(intersect_type == EI_SEGMENT)
		{
			return(intersect);
		}
	}

	return vector_2d_c(INFINITY, INFINITY);
}


airfoil_c airfoil_c::get_circle_projection(uint32_t n, const vector_2d_c & projection_center, double radius, double epsilon) const
{
	if(!is_valid())
	{
		std::cout << "WARNING: TRYING TO GET PROJECT ON AN INVALID AIRFOIL" << std::endl;
		return airfoil_c();
	}

	std::vector<vector_2d_c> newpoints;
	std::stringstream newname;

	airfoil_c temp = closed_merge(0.0001);

	for(unsigned int i = 0; i < n; i++)
	{
		vector_2d_c circle_point = vector_2d_radian(radius, (M_PI * i) / n) + projection_center;
		vector_2d_c top_point = vector_2d_c(circle_point.x, 1);
		vector_2d_c inverse_point = vector_2d_c(circle_point.x, -1);
		line_2d_c vert_line(top_point, inverse_point);
		vector_2d_c intersect = temp.get_intersection_first(vert_line);

		if(intersect.x != INFINITY)
		{
			newpoints.push_back(intersect);
		}
	}

	for(unsigned int i = n; i < 2 * n; i++)
	{
		vector_2d_c circle_point = vector_2d_radian(radius, (M_PI * i) / n) + projection_center;
		vector_2d_c top_point = vector_2d_c(circle_point.x, 1);
		vector_2d_c inverse_point = vector_2d_c(circle_point.x, -1);
		line_2d_c vert_line(top_point, inverse_point);
		vector_2d_c intersect = temp.get_intersection_last(vert_line);

		if(intersect.x != INFINITY)
		{
			newpoints.push_back(intersect);
		}
	}

	newpoints.push_back(newpoints.front());
	newname << this->name << " circle projected with " << n << " subdivisions centered on " << projection_center;


	return airfoil_c(newpoints, newname.str());
}

airfoil_c airfoil_c::get_circle_projection(uint32_t n, double epsilon) const
{
	if(!is_valid())
	{
		std::cout << "WARNING: TRYING TO GET PROJECT ON AN INVALID AIRFOIL" << std::endl;

		return airfoil_c();
	}

	double x_max = (*std::max_element(this->points.begin(), this->points.end(), [](const vector_2d_c & lhs, const vector_2d_c & rhs)
	{
		return lhs.x < rhs.x;
	})).x;
	double x_min = (*std::min_element(this->points.begin(), this->points.end(), [](const vector_2d_c & lhs, const vector_2d_c & rhs)
	{
		return lhs.x < rhs.x;
	})).x;

	return get_circle_projection(n, {0.5 * (x_max - x_min), 0.0}, 0.5 * (x_max - x_min), epsilon);
}


bool airfoil_c::is_closed(double epsilon) const
{
	double lsq = (this->points.front() - this->points.back()).get_length_sq();

	return !is_valid() or (lsq <= (epsilon * epsilon));
}


airfoil_c airfoil_c::closed_merge(double epsilon) const
{
	if(this->is_closed(epsilon))
	{
		return *this;
	}

	vector_2d_c endpoint = (points.front() + points.back()) * 0.5;
	std::vector<vector_2d_c> newpoints = this->points;
	newpoints.front() = endpoint;
	newpoints.back() = endpoint;

	std::stringstream newname;
	newname << this->name << " closed";
	return airfoil_c(newpoints, newname.str());
}

airfoil_c airfoil_c::closed_intersect(double epsilon) const
{
	if(this->is_closed(epsilon))
	{
		return *this;
	}

	vector_2d_c endpoint;

	this->get_lines().front().get_intersection(this->get_lines().back(), endpoint, 0);

	std::vector<vector_2d_c> newpoints;
	newpoints.push_back(endpoint);
	newpoints.insert(newpoints.end(), this->points.begin(), this->points.end());
	newpoints.push_back(endpoint);

	std::stringstream newname;
	newname << this->name << " closed";
	return airfoil_c(newpoints, newname.str());
}


bool airfoil_c::is_valid() const
{
	return this->points.size() > 1;
}


std::string airfoil_c::get_name() const
{
	return this->name;
}


void airfoil_c::set_name(const std::string & new_name)
{
	this->name = new_name;
}


std::ostream & operator << (std::ostream & output, const airfoil_c & airfoil)
{
	output << airfoil.name << std::endl;

	for(const vector_2d_c & v : airfoil.points)
	{
		output << v.x << "\t" << v.y << std::endl;
	}

	return output;
}


} //namespace wif_core
