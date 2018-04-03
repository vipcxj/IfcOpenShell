/********************************************************************************
 *                                                                              *
 * This file is part of IfcOpenShell.                                           *
 *                                                                              *
 * IfcOpenShell is free software: you can redistribute it and/or modify         *
 * it under the terms of the Lesser GNU General Public License as published by  *
 * the Free Software Foundation, either version 3.0 of the License, or          *
 * (at your option) any later version.                                          *
 *                                                                              *
 * IfcOpenShell is distributed in the hope that it will be useful,              *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of               *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                 *
 * Lesser GNU General Public License for more details.                          *
 *                                                                              *
 * You should have received a copy of the Lesser GNU General Public License     *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                              *
 ********************************************************************************/

/********************************************************************************
 *                                                                              *
 * This started as a brief example of how IfcOpenShell can be interfaced from   * 
 * within a C++ context, it has since then evolved into a fullfledged command   *
 * line application that is able to convert geometry in an IFC files into       *
 * several tesselated and topological output formats.                           *
 *                                                                              *
 ********************************************************************************/

// #include "../ifcconvert/ColladaSerializer.h"
// #include "../ifcconvert/IgesSerializer.h"
// #include "../ifcconvert/StepSerializer.h"
// #include "../ifcconvert/WavefrontObjSerializer.h"
// #include "../ifcconvert/XmlSerializer.h"
// #include "../ifcconvert/SvgSerializer.h"

// #include "../ifcgeom/IfcGeomIterator.h"

// #include <IGESControl_Controller.hxx>
#include <Standard_Version.hxx>
#include <map>
#include <set>
#include <vector>
#include <limits>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include <boost/program_options.hpp>
#include "../ifcparse/IfcFile.h"
#include "../ifcgeom/IfcGeom.h"
#include "../ifcgeom/IfcGeomElement.h"
#include "../ifcgeom/IfcGeomMaterial.h"
#include "../ifcgeom/IfcGeomIteratorSettings.h"
#include "../ifcgeom/IfcRepresentationShapeItem.h"
#include "../ifcgeom/IfcGeomFilter.h"

#include <fstream>
#include <sstream>
#include <set>
#include <time.h>

namespace po = boost::program_options;

void print_version()
{
    std::cout << "IfcOpenShell "
              << IfcSchema::Identifier
              << " IfcConvert "
              << IFCOPENSHELL_VERSION
              << " (OCC "
              << OCC_VERSION_STRING_EXT
              << ")\n";
}

void print_usage(bool t_suggest_help = true)
{
  std::cout << "Usage: IfcInfo [options] <input.ifc>\n"
            << "\n"
            << "Displays the geometry in an IFC file.\n";
  if (t_suggest_help) {
    std::cout << "\nRun 'IfcInfo --help' for more information.";
  }
  std::cout << std::endl;
}

/// @todo Add help for single option
void print_options(const po::options_description& options)
{
    std::cout << "\n" << options;
    std::cout << std::endl;
}

bool file_exists(const std::string& filename)
{
    /// @todo Windows Unicode support
    std::ifstream file(filename.c_str());
    return file.good();
}

// bool init_input_file(const std::string& filename, IfcParse::IfcFile& ifc_file, bool no_progress, bool mmap);

int main(int ac, char** av)
{
  
  print_version();
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "produce help message")
      ("y", "yes on all")
      ("input-file", po::value<std::string>(), "input file");

  po::positional_options_description positional_options;
  positional_options.add("input-file", -1);

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).
            options(desc).positional(positional_options).run(), vm);
  
  po::notify(vm);

  if(vm.count("help"))
  {
    std::cout << "help: \n"
              << desc << std::endl;
    return EXIT_SUCCESS;
  }else if (!vm.count("input-file"))
  {
    std::cerr << "[Error] Input file not specified" << std::endl;
    std::cout << desc << std::endl;
    return EXIT_FAILURE;
  }

  const std::string input_filename = vm["input-file"].as<std::string>();
  if (!file_exists(input_filename)) {
    std::cerr << "[Error] Input file '" << input_filename << "' does not exist" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "input file: " << input_filename << std::endl;
  IfcParse::IfcFile ifc_file;
  
  
  if (!ifc_file.Init(input_filename)) {
    Logger::Error("Unable to parse input file '" + input_filename + "'");
    return EXIT_FAILURE;
  }
  //if (no_progress) { Logger::SetOutput(&std::cout, &log_stream); }


  
  return 1;
}
