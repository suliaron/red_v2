#include <algorithm>
#include <cstring>
#include <iostream>

#include "parameter.h"
#include "redutil2.h"

#include "constants.h"

using namespace std;
using namespace redutil2;

parameter::parameter(string& dir, string& filename, bool verbose) :
	filename(filename),
	verbose(verbose)
{
	create_default();

	string path;
	if (0 < file::get_directory(filename).length())
	{
		path = filename;
	}
	else
	{
		path = file::combine_path(dir, filename);
	}
	file::load_ascii_file(path, data);
	parse();
}

parameter::~parameter() 
{
}

void parameter::create_default()
{
	adaptive           = false;
	error_check_for_tp = false;
	output_data_rep    = DATA_REPRESENTATION_ASCII;
	int_type           = INTEGRATOR_RUNGEKUTTA7;
	tolerance          = 1.0e-10;
    eta                = 0.02;
	simulation_length  = 0.0;
	output_interval    = 0.0;
	info_dt            = 5.0;      // [sec]
	dump_dt            = 3600.0;   // [sec]

	memset(threshold, 0, THRESHOLD_N * sizeof(var_t));
}

void parameter::parse()
{
	// instantiate Tokenizer classes
	Tokenizer data_tokenizer;
	Tokenizer line_tokenizer;
	string line;

	data_tokenizer.set(data, "\r\n");

	while ((line = data_tokenizer.next()) != "")
	{
		line_tokenizer.set(line, "=");
		string token;
		int tokenCounter = 1;

		string key; 
		string value;
		while ((token = line_tokenizer.next()) != "" && tokenCounter <= 2)
		{
			if (1 == tokenCounter)
			{
				key = token;
			}
			if (2 == tokenCounter)
			{
				value = token;
			}
			tokenCounter++;
		}
		if (2 < tokenCounter)
		{
			set_param(key, value);
		}
		else
		{
			throw string("Invalid key/value pair: " + line + ".");
		}
	}
}

void parameter::set_param(string& key, string& value)
{
	static char n_call = 0;

	n_call++;
	key = tools::trim(key);
	value = tools::trim(value);
	transform(key.begin(), key.end(), key.begin(), ::tolower);

	if (     key == "name")
	{
		simulation_name = value;
    } 
    else if (key == "description")
	{
		simulation_desc = value;
    }
    else if (key == "output data representation" || key == "odr")
	{
		simulation_desc = value;
		transform(value.begin(), value.end(), value.begin(), ::tolower);
		if (     value == "t" || value == "text")
		{
			output_data_rep = DATA_REPRESENTATION_ASCII;
		}
		else if (value == "b" || value == "binary")
		{
			output_data_rep = DATA_REPRESENTATION_BINARY;
		}
		else
		{
			throw string("Invalid data type: " + value + ".");
		}
    }
    else if (key == "integrator")
	{
		transform(value.begin(), value.end(), value.begin(), ::tolower);
		if (     value == "e"   || value == "euler")
		{
			int_type = INTEGRATOR_EULER;
		}
		else if (value == "rk2" || value == "rungekutta2")
		{
			int_type = INTEGRATOR_RUNGEKUTTA2;
		}
		else if (value == "rk4" || value == "rungekutta4")
		{
			int_type = INTEGRATOR_RUNGEKUTTA4;
		}
		else if (value == "rk5" || value == "rungekutta5")
		{
			int_type = INTEGRATOR_RUNGEKUTTA5;
		}
		else if (value == "rk7" || value == "rungekutta7")
		{
			int_type = INTEGRATOR_RUNGEKUTTA7;
		}			
		else if (value == "h4" || value == "hermite4")
		{
			int_type = INTEGRATOR_HERMITE4;
		}			
		else if (value == "h4b" || value == "hermite4b")
		{
			int_type = INTEGRATOR_HERMITE4B;
		}			
		else
		{
			throw string("Invalid integrator type: " + value + ".");
		}
	}
    else if (key == "tolerance")
	{
		if (!tools::is_number(value))
		{
			throw string("Invalid number at: " + key + ".");
		}
		adaptive = true;
		tolerance = atof(value.c_str());
	}
    else if (key == "eta")
	{
		if (!tools::is_number(value))
		{
			throw string("Invalid number at: " + key + ".");
		}
		adaptive = true;
		eta = atof(value.c_str());
	}
	else if (key == "error check for tp")
	{
		if      (value == "true")
		{
			error_check_for_tp = true;
		}
		else if (value == "false")
		{
			error_check_for_tp = false;
		}
		else
		{
			throw string("Invalid value at: " + key + ".");
		}
	}
    else if (key == "length")
	{
		if (!tools::is_number(value))
		{
			throw string("Invalid number at: " + key + ".");
		}
		simulation_length = atof(value.c_str()) * constants::YearToDay;
	}
    else if (key == "output interval")
	{
		if (!tools::is_number(value))
		{
			throw string("Invalid number at: " + key + ".");
		}
		output_interval = atof(value.c_str()) * constants::YearToDay;
	}
	else if (key == "info dt")
	{
		if (!tools::is_number(value)) 
		{
			throw string("Invalid number at: " + key + ".");
		}
		info_dt = atof(value.c_str());
	}
	else if (key == "dump dt")
	{
		if (!tools::is_number(value)) 
		{
			throw string("Invalid number at: " + key + ".");
		}
		dump_dt = atof(value.c_str());
	}
    else if (key == "ejection")
    {
        if (!tools::is_number(value))
        {
            throw string("Invalid number at: " + key);
        }
        threshold[THRESHOLD_EJECTION_DISTANCE] = atof(value.c_str());
    }
    else if (key == "hit centrum")
    {
        if (!tools::is_number(value))
        {
            throw string("Invalid number at: " + key);
        }
        threshold[THRESHOLD_HIT_CENTRUM_DISTANCE] = atof(value.c_str());
    }
    else if (key == "radii enhance factor" || key == "collision factor")
    {
        if (!tools::is_number(value))
        {
            throw string("Invalid number at: " + key);
        }
        threshold[THRESHOLD_RADII_ENHANCE_FACTOR] = atof(value.c_str());
    }
    else
	{
		throw string("Invalid parameter: " + key + ".");
	}

	if (verbose)
	{
		if (n_call == 1)
		{
			cout << "The following parameters are setted (from file :'" << filename << "'): " << endl;
		}
		cout << "\t'" << key << "' was assigned to '" << value << "'" << endl;
	}
}

ostream& operator<<(ostream& stream, const parameter* p)
{
	const char* integrator_name[] = 
		{
			"INTEGRATOR_EULER"
			"INTEGRATOR_RUNGEKUTTA2",
			"INTEGRATOR_RUNGEKUTTA4",
			"INTEGRATOR_RUNGEKUTTA5",
			"INTEGRATOR_RUNGEKUTTA7"
		};

	stream << "simulation name           : " << p->simulation_name << endl;
	stream << "simulation description    : " << p->simulation_desc << endl;
	stream << "output data representation: " << (DATA_REPRESENTATION_ASCII == p->output_data_rep ? "text" : "binary") << endl;
	stream << "simulation integrator     : " << integrator_name[p->int_type] << endl;
	stream << "simulation tolerance      : " << p->tolerance << endl;
	stream << "simulation adaptive       : " << (p->adaptive ? "true" : "false") << endl;
	stream << "simulation length         : " << p->simulation_length << endl;
	stream << "simulation output interval: " << p->output_interval << endl;

	return stream;
}
