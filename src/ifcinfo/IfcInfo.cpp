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

// #include <IGESControl_Controller.hxx>
#include <Standard_Version.hxx>
#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "../ifcgeom/IfcGeom.h"
#include "../ifcgeom/IfcGeomElement.h"
#include "../ifcgeom/IfcGeomFilter.h"
#include "../ifcgeom/IfcGeomIterator.h"
#include "../ifcgeom/IfcGeomIteratorSettings.h"
#include "../ifcgeom/IfcGeomMaterial.h"
#include "../ifcgeom/IfcRepresentationShapeItem.h"
#include "../ifcparse/IfcFile.h"

#include <fstream>
#include <set>
#include <sstream>
#include <time.h>

#if USE_VLD
#include <vld.h>
#endif

using namespace IfcSchema;
int main(int argc, char** argv) 
{

  if (argc != 2)
  {
    std::cout << "usage: IfcParseExamples <filename.ifc>" << std::endl;
    return 1;
  }

  // Redirect the output (both progress and log) to stdout
  Logger::SetOutput(&std::cout, &std::cout);

  // Parse the IFC file provided in argv[1]
  IfcParse::IfcFile file;
  if (!file.Init(argv[1]))
  {
    std::cout << "Unable to parse .ifc file" << std::endl;
    return 1;
  }

  // Lets get a list of IfcBuildingElements, this is the parent
  // type of things like walls, windows and doors.
  // entitiesByType is a templated function and returns a
  // templated class that behaves like a std::vector.
  // Note that the return types are all typedef'ed as members of
  // the generated classes, ::list for the templated vector class,
  // ::ptr for a shared pointer and ::it for an iterator.
  // We will simply iterate over the vector and print a string
  // representation of the entity to stdout.
  //
  // Secondly, lets find out which of them are IfcWindows.
  // In order to access the additional properties that windows
  // have on top af the properties of building elements,
  // we need to cast them to IfcWindows. Since these properties
  // are optional we need to make sure the properties are
  // defined for the window in question before accessing them.

  // IfcBuildingElement::list::ptr elements = file.entitiesByType<IfcBuildingElement>();

  // std::cout << "Found " << elements->size() << " elements in " << argv[1] << ":" << std::endl;

  // for (IfcBuildingElement::list::it it = elements->begin(); it != elements->end(); ++it)
  // {
  //   const IfcBuildingElement *element = *it;
  //   std::cout << element->entity->toString() << std::endl;
  //   if (element->is(IfcWindow::Class()))
  //   {
  //     const IfcWindow *window = (IfcWindow *)element;
  //     if (window->hasOverallWidth() && window->hasOverallHeight())
  //     {
  //       const double area = window->OverallWidth() * window->OverallHeight();
  //       std::cout << "The area of this window is " << area << std::endl;
  //     }
  //   }
  // }

  IfcProduct::list::ptr products = file.entitiesByType<IfcProduct>();

  std::cout << "Found " << products->size() << " products in " << argv[1] << ":" << std::endl;

  for (IfcProduct::list::it it = products->begin(); it != products->end(); ++it)
  {
    const IfcProduct *product = *it;
    std::cout << product->entity->toString() << std::endl;
    // if (element->is(IfcWindow::Class()))
    // {
    //   const IfcWindow *window = (IfcWindow *)element;
    //   if (window->hasOverallWidth() && window->hasOverallHeight())
    //   {
    //     const double area = window->OverallWidth() * window->OverallHeight();
    //     std::cout << "The area of this window is " << area << std::endl;
    //   }
    // }
  }


  return 1;
}
