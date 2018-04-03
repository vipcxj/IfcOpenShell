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
 * Geometrical data in an IFC file consists of shapes (IfcShapeRepresentation)  *
 * and instances (SUBTYPE OF IfcBuildingElement e.g. IfcWindow).                *
 *                                                                              *
 * IfcGeom::Representation::Triangulation is a class that represents a          *
 * triangulated IfcShapeRepresentation.                                         *
 *   Triangulation.verts is a 1 dimensional vector of float defining the        *
 *      cartesian coordinates of the vertices of the triangulated shape in the  *
 *      format of [x1,y1,z1,..,xn,yn,zn]                                        *
 *   Triangulation.faces is a 1 dimensional vector of int containing the        *
 *     indices of the triangles referencing positions in Triangulation.verts    *
 *   Triangulation.edges is a 1 dimensional vector of int in {0,1} that dictates*
 *	   the visibility of the edges that span the faces in Triangulation.faces   *
 *                                                                              *
 * IfcGeom::Element represents the actual IfcBuildingElements.                  *
 *   IfcGeomObject.name is the GUID of the element                              *
 *   IfcGeomObject.type is the datatype of the element e.g. IfcWindow           *
 *   IfcGeomObject.mesh is a pointer to an IfcMesh                              *
 *   IfcGeomObject.transformation.matrix is a 4x3 matrix that defines the       *
 *     orientation and translation of the mesh in relation to the world origin  *
 *                                                                              *
 * IfcGeom::Iterator::initialize()                                              *
 *   finds the most suitable representation contexts. Returns true iff          *
 *   at least a single representation will process successfully                 *
 *                                                                              *
 * IfcGeom::Iterator::get()                                                     *
 *   returns a pointer to the current IfcGeom::Element                          *
 *                                                                              * 
 * IfcGeom::Iterator::next()                                                    *
 *   returns true iff a following entity is available for a successive call to  *
 *   IfcGeom::Iterator::get()                                                   *
 *                                                                              *
 * IfcGeom::Iterator::progress()                                                *
 *   returns an int in [0..100] that indicates the overall progress             *
 *                                                                              *
 ********************************************************************************/

#ifndef IFCGEOMITERATORMAP_H
#define IFCGEOMITERATORMAP_H

#include <map>
#include <set>
#include <vector>
#include <limits>
#include <algorithm>

#include <boost/algorithm/string.hpp>

#include <gp_Mat.hxx>
#include <gp_Mat2d.hxx>
#include <gp_GTrsf.hxx>
#include <gp_GTrsf2d.hxx>
#include <gp_Trsf.hxx>
#include <gp_Trsf2d.hxx>

#include "../ifcparse/IfcFile.h"

#include "../ifcgeom/IfcGeom.h"
#include "../ifcgeom/IfcGeomElement.h"
#include "../ifcgeom/IfcGeomMaterial.h"
#include "../ifcgeom/IfcGeomIteratorSettings.h"
#include "../ifcgeom/IfcRepresentationShapeItem.h"
#include "../ifcgeom/IfcGeomFilter.h"

// The infamous min & max Win32 #defines can leak here from OCE depending on the build configuration
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

// style freedom
// http://www.ivanism.com/Articles/CodingStandards.html
// https://lefticus.gitbooks.io/cpp-best-practices/content/03-Style.html
// prefix m_ for private data (member/private data and functions)
// prefix t_ for function parameters
// forget starting with _, possibilitiy for reserved compiler and stl implementations

namespace IfcGeom {
	
template <typename P>
class IteratorMap {
 private:
  Kernel m_kernel;
  IteratorSettings m_settings;
  IfcParse::IfcFile* m_ifc_file;
  IfcSchema::IfcRepresentation::list::ptr m_representations;
  IfcSchema::IfcRepresentation::list::it m_representation_iterator;
  TriangulationElement<P>* m_current_triangulation;
  BRepElement<P>* m_current_shape_model;
  SerializedElement<P>* m_current_serialization;
  IfcSchema::IfcProduct::list::ptr m_ifcproducts;
  IfcSchema::IfcProduct::list::it m_ifcproduct_iterator;
  int m_done;
  int m_total;
  std::string m_unit_name;
  double m_unit_magnitude;
  gp_XYZ m_bounds_min;
  gp_XYZ m_bounds_max;
  std::vector<filter_t> m_filters;

  struct filter_match
  {
    filter_match (IfcSchema::IfcProduct *prod) : product(prod) {}
    bool operator() (const filter_t& filter) const { return filter(product);  }
    IfcSchema::IfcProduct* product;
  };

  void initUnits() {
    IfcSchema::IfcProject::list::ptr projects = ifc_file->entitiesByType<IfcSchema::IfcProject>();
    if (projects->size() == 1) {
      IfcSchema::IfcProject* project = *projects->begin();
      std::pair<std::string, double> length_unit = kernel.initializeUnits(project->UnitsInContext());
      unit_name = length_unit.first;
      unit_magnitude = length_unit.second;
    }
  }
  void m_initialize() {
    current_triangulation = 0;
    current_shape_model = 0;
    current_serialization = 0;

    unit_name = "METER";
    unit_magnitude = 1.f;

    kernel.setValue(IfcGeom::Kernel::GV_MAX_FACES_TO_SEW, settings.get(IteratorSettings::SEW_SHELLS) ? 1000 : -1);
    kernel.setValue(IfcGeom::Kernel::GV_DIMENSIONALITY, (settings.get(IteratorSettings::INCLUDE_CURVES)
                                                         ? (settings.get(IteratorSettings::EXCLUDE_SOLIDS_AND_SURFACES) ? -1. : 0.) : +1.));
    if (settings.get(IteratorSettings::BUILDING_LOCAL_PLACEMENT)) {
      if (settings.get(IteratorSettings::SITE_LOCAL_PLACEMENT)) {
        Logger::Message(Logger::LOG_WARNING, "building-local-placement takes precedence over site-local-placement");
      }
      kernel.set_conversion_placement_rel_to(IfcSchema::Type::IfcBuilding);
    } else if (settings.get(IteratorSettings::SITE_LOCAL_PLACEMENT)) {
      kernel.set_conversion_placement_rel_to(IfcSchema::Type::IfcSite);
    }
  }

  bool owns_ifc_file;

 public:
  IteratorMap(
      const IteratorSettings& settings,
      IfcParse::IfcFile* file,
      std::vector<IfcGeom::filter_t>& filters
              ) : settings(settings), ifc_file(file), owns_ifc_file(false), filters_(filters)
  {
    m_initialize();
  }
  
};
}

#endif

