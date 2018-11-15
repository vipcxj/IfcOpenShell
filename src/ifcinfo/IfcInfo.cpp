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
// #include "../ifcconvert/SvgSerializer.h"
// #include "../ifcconvert/XmlSerializer.h"
#include "../ifcconvert/GeometrySerializer.h"

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

bool reuse_ok_(SerializerSettings settings, const IfcSchema::IfcProduct::list::ptr &products)
{
  IfcGeom::Kernel kernel;

  // With world coords enabled, object transformations are directly applied to
  // the BRep. There is no way to re-use the geometry for multiple products.
  if (settings.get(IfcGeom::IteratorSettings::USE_WORLD_COORDS))
  {
    return false;
  }

  std::set<const IfcSchema::IfcMaterial *> associated_single_materials;

  for (IfcSchema::IfcProduct::list::it it = products->begin(); it != products->end(); ++it)
  {
    IfcSchema::IfcProduct *product = *it;
    if (!settings.get(IfcGeom::IteratorSettings::DISABLE_OPENING_SUBTRACTIONS) &&
        kernel.find_openings(product)->size())
    {
      return false;
    }
    if (settings.get(IfcGeom::IteratorSettings::APPLY_LAYERSETS))
    {
      IfcSchema::IfcRelAssociates::list::ptr associations = product->HasAssociations();
      for (IfcSchema::IfcRelAssociates::list::it jt = associations->begin();
           jt != associations->end(); ++jt)
      {
        IfcSchema::IfcRelAssociatesMaterial *assoc =
            (*jt)->as<IfcSchema::IfcRelAssociatesMaterial>();
        if (assoc)
        {
          if (assoc->RelatingMaterial()->is(IfcSchema::Type::IfcMaterialLayerSetUsage))
          {
            // TODO: Check whether single layer?
            return false;
          }
        }
      }
    }
    // Note that this can be a nullptr (!), but the fact that set size should be one still holds
    associated_single_materials.insert(kernel.get_single_material_association(product));
    if (associated_single_materials.size() > 1)
      return false;
  }
  return associated_single_materials.size() == 1;
}

int main(int argc, char **argv)
{
  std::chrono::time_point<std::chrono::system_clock> start, end;
  if (argc != 2)
  {
    std::cout << "usage: IfcParseExamples <filename.ifc>" << std::endl;
    return 1;
  }
  // Redirect the output (both progress and log) to stdout
  Logger::SetOutput(&std::cout, &std::cout);
  // Parse the IFC file provided in argv[1]
  Logger::Status("Ifc_File.Init");
  start = std::chrono::system_clock::now();
  IfcParse::IfcFile ifc_file;
  if (!ifc_file.Init(argv[1]))
  {
    std::cout << "Unable to parse .ifc file" << std::endl;
    return 1;
  }
  end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  Logger::Status("Ifc_File.Init duration: " + std::to_string(elapsed_seconds.count()));

  SerializerSettings settings;
  /// @todo Make APPLY_DEFAULT_MATERIALS configurable? Quickly tested setting this to false and
  /// using obj exporter caused the program to crash and burn.
  settings.set(IfcGeom::IteratorSettings::APPLY_DEFAULT_MATERIALS, true);
  settings.set(IfcGeom::IteratorSettings::USE_WORLD_COORDS, true);
  settings.set(IfcGeom::IteratorSettings::WELD_VERTICES, true);
  settings.set(IfcGeom::IteratorSettings::SEW_SHELLS, true);
  settings.set(IfcGeom::IteratorSettings::CONVERT_BACK_UNITS, true);
#if OCC_VERSION_HEX < 0x60900
  settings.set(IfcGeom::IteratorSettings::FASTER_BOOLEANS, true);
#endif
  settings.set(IfcGeom::IteratorSettings::DISABLE_OPENING_SUBTRACTIONS, true);
  settings.set(IfcGeom::IteratorSettings::INCLUDE_CURVES, false);
  settings.set(IfcGeom::IteratorSettings::EXCLUDE_SOLIDS_AND_SURFACES, true);
  settings.set(IfcGeom::IteratorSettings::APPLY_LAYERSETS, true);
  settings.set(IfcGeom::IteratorSettings::NO_NORMALS, false);
  settings.set(IfcGeom::IteratorSettings::GENERATE_UVS, true);
  settings.set(IfcGeom::IteratorSettings::SEARCH_FLOOR, true);
  settings.set(IfcGeom::IteratorSettings::SITE_LOCAL_PLACEMENT, true);
  settings.set(IfcGeom::IteratorSettings::BUILDING_LOCAL_PLACEMENT, true);
  settings.set(SerializerSettings::USE_ELEMENT_NAMES, true);
  settings.set(SerializerSettings::USE_ELEMENT_GUIDS, true);
  settings.set(SerializerSettings::USE_MATERIAL_NAMES, true);
  settings.set(SerializerSettings::USE_ELEMENT_TYPES, true);
  settings.set(SerializerSettings::USE_ELEMENT_HIERARCHY, true);
  settings.set_deflection_tolerance(0.001);
  settings.precision = 7;

  std::set<std::string> allowed_context_types;
  allowed_context_types.insert("model");
  allowed_context_types.insert("plan");
  allowed_context_types.insert("notdefined");
  std::set<std::string> context_types;
  if (!settings.get(IfcGeom::IteratorSettings::EXCLUDE_SOLIDS_AND_SURFACES))
  {
    context_types.insert("model");
    context_types.insert("design");
    context_types.insert("model view");
    context_types.insert("detail view");
  }
  if (settings.get(IfcGeom::IteratorSettings::INCLUDE_CURVES))
  {
    context_types.insert("plan");
  }

  double lowest_precision_encountered = std::numeric_limits<double>::infinity();
  bool any_precision_encountered = false;

  IfcSchema::IfcRepresentation::list::ptr representations =
      IfcSchema::IfcRepresentation::list::ptr(new IfcSchema::IfcRepresentation::list);
  IfcSchema::IfcRepresentation::list::it representation_iterator;
  IfcSchema::IfcRepresentation::list::ptr ok_mapped_representations;
  IfcSchema::IfcGeometricRepresentationContext::list::ptr contexts =
      ifc_file.entitiesByType<IfcSchema::IfcGeometricRepresentationContext>();
  IfcSchema::IfcGeometricRepresentationContext::list::ptr filtered_contexts(
      new IfcSchema::IfcGeometricRepresentationContext::list);
  IfcSchema::IfcGeometricRepresentationContext::list::it it;
  IfcSchema::IfcGeometricRepresentationSubContext::list::it jt;

  ////////////////////////////////////////////////////////////
  ///////////////////////// filter contexts for representations
  ////////////////////////////////////////////////////////////
  Logger::Status("go through contexts");
  start = std::chrono::system_clock::now();
  for (it = contexts->begin(); it != contexts->end(); ++it)
  {
    IfcSchema::IfcGeometricRepresentationContext *context = *it;
    if (context->is(IfcSchema::Type::IfcGeometricRepresentationSubContext))
    {
      // Continue, as the list of subcontexts will be considered
      // by the parent's context inverse attributes.
      continue;
    }
    try
    {
      if (context->hasContextType())
      {
        std::string context_type = context->ContextType();
        boost::to_lower(context_type);
        if (allowed_context_types.find(context_type) == allowed_context_types.end())
        {
          Logger::Message(Logger::LOG_ERROR,
                          std::string("ContextType '") + context->ContextType() + "' not allowed:",
                          context->entity);
        }
        if (context_types.find(context_type) != context_types.end())
        {
          filtered_contexts->push(context);
        }
      }
    }
    catch (const std::exception &e)
    {
      Logger::Error(e);
    }
  }

  // In case no contexts are identified based on their ContextType, all contexts are
  // considered. Note that sub contexts are excluded as they are considered later on.
  if (filtered_contexts->size() == 0)
  {
    Logger::Status("no filtered contexts, all contexts are considered..");
    for (it = contexts->begin(); it != contexts->end(); ++it)
    {
      IfcSchema::IfcGeometricRepresentationContext *context = *it;
      if (!context->is(IfcSchema::Type::IfcGeometricRepresentationSubContext))
      {
        filtered_contexts->push(context);
      }
    }
  }
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  Logger::Status("context duration: " + std::to_string(elapsed_seconds.count()));

  Logger::Status("Iterate over filtered_contexts ... ");
  start = std::chrono::system_clock::now();
  for (it = filtered_contexts->begin(); it != filtered_contexts->end(); ++it)
  {
    IfcSchema::IfcGeometricRepresentationContext *context = *it;
    representations->push(context->RepresentationsInContext());
    try
    {
      if (context->hasPrecision() && context->Precision() < lowest_precision_encountered)
      {
        lowest_precision_encountered = context->Precision();
        any_precision_encountered = true;
      }
    }
    catch (const std::exception &e)
    {
      Logger::Error(e);
    }
    IfcSchema::IfcGeometricRepresentationSubContext::list::ptr sub_contexts =
        context->HasSubContexts();
    for (jt = sub_contexts->begin(); jt != sub_contexts->end(); ++jt)
    {
      representations->push((*jt)->RepresentationsInContext());
    }
    // There is no need for full recursion as the following is governed by the schema:
    // WR31: The parent context shall not be another geometric representation sub context.
  }
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  Logger::Status("filtered_context duration: " + std::to_string(elapsed_seconds.count()));

  ////////////////////////////////////////////////////////////
  /////// representations filled, now set units and precision
  ////////////////////////////////////////////////////////////
  IfcGeom::Kernel kernel;
  // initUnits()
  IfcSchema::IfcProject::list::ptr projects = ifc_file.entitiesByType<IfcSchema::IfcProject>();
  std::string unit_name;
  double unit_magnitude;
  unit_name = "METER";
  unit_magnitude = 1.f;
  if (projects->size() == 1)
  {
    IfcSchema::IfcProject *project = *projects->begin();
    std::pair<std::string, double> length_unit = kernel.initializeUnits(project->UnitsInContext());
    unit_name = length_unit.first;
    unit_magnitude = length_unit.second;
  }
  else
  {
    Logger::Error("A single IfcProject is expected (encountered " +
                  boost::lexical_cast<std::string>(projects->size()) +
                  "); unable to read unit information.");
  }
  if (any_precision_encountered)
  {
    // Some arbitrary factor that has proven to work better for the models in the set of test
    // files.
    lowest_precision_encountered *= 10.;
    lowest_precision_encountered *= unit_magnitude;
    if (lowest_precision_encountered < 1.e-7)
    {
      Logger::Message(Logger::LOG_WARNING, "Precision lower than 0.0000001 meter not enforced");
      kernel.setValue(IfcGeom::Kernel::GV_PRECISION, 1.e-7);
    }
    else
    {
      kernel.setValue(IfcGeom::Kernel::GV_PRECISION, lowest_precision_encountered);
    }
  }
  else
  {
    kernel.setValue(IfcGeom::Kernel::GV_PRECISION, 1.e-5);
  }

  if (representations->size() == 0)
  {
    Logger::Message(Logger::LOG_ERROR, "No geometries found");
    return 0;
  }

  ////////////////////////////////////////////////////////////
  /////// find products associated with representations (... ?)
  ////////////////////////////////////////////////////////////
  IfcSchema::IfcProduct::list::ptr ifcproducts;
  IfcSchema::IfcProduct::list::it ifcproduct_iterator;
  std::vector<IfcGeom::filter_t> filters_;

  IfcGeom::layer_filter layer_filter;
  IfcGeom::entity_filter entity_filter;
  IfcGeom::string_arg_filter guid_filter(IfcSchema::Type::IfcRoot, 0);
  IfcGeom::string_arg_filter name_filter(IfcSchema::Type::IfcRoot, 2);
  IfcGeom::string_arg_filter desc_filter(IfcSchema::Type::IfcRoot, 3);
  IfcGeom::string_arg_filter tag_filter(IfcSchema::Type::IfcProxy, 8, IfcSchema::Type::IfcElement,
                                        7);
  filters_.emplace_back(boost::ref(layer_filter));
  filters_.emplace_back(boost::ref(entity_filter));
  filters_.emplace_back(boost::ref(guid_filter));
  filters_.emplace_back(boost::ref(name_filter));
  filters_.emplace_back(boost::ref(desc_filter));
  filters_.emplace_back(boost::ref(tag_filter));
  bool geometry_reuse_ok_for_current_representation_;

  // functor
  struct filter_match
  {
    filter_match(IfcSchema::IfcProduct *prod) : product(prod) {}
    bool operator()(const IfcGeom::filter_t &filter) const { return filter(product); }
    IfcSchema::IfcProduct *product;
  };
  
struct IfcproductRepresentation
{
  int index;
  IfcSchema::IfcRepresentation *representation;
  IfcSchema::IfcProduct *product;
};

std::vector<IfcproductRepresentation> IfcproductRepresentations;

  start = std::chrono::system_clock::now();
  int index_count=0;
  for (representation_iterator = representations->begin();
       representation_iterator != representations->end(); representation_iterator++)
  {
    IfcSchema::IfcRepresentation *representation;
    representation = *representation_iterator;
    ifcproducts.reset();
    ifcproducts = IfcSchema::IfcProduct::list::ptr(new IfcSchema::IfcProduct::list);
    IfcSchema::IfcProduct::list::ptr unfiltered_products =
        kernel.products_represented_by(representation);
    geometry_reuse_ok_for_current_representation_ = reuse_ok_(settings, unfiltered_products);
    IfcSchema::IfcRepresentationMap::list::ptr maps = representation->RepresentationMap();
    if (!geometry_reuse_ok_for_current_representation_ && maps->size() == 1)
    {
      // unfiltered_products contains products represented by this representation by means of
      // mapped items. For example because of openings applied to products, reuse might not be
      // acceptable and then the products will be processed by means of their immediate
      // representation and not the mapped representation.

      // IfcRepresentationMaps are also used for IfcTypeProducts, so an additional check is
      // performed whether the map is indeed used by IfcMappedItems.
      IfcSchema::IfcRepresentationMap *map = *maps->begin();
      if (map->MapUsage()->size() > 0)
      {
        // _nextShape();
        // continue;
        // NOTE(sander): is this equivalent to _nextShape() ?
        continue;
      }
    }

    bool representation_processed_as_mapped_item = false;
    IfcSchema::IfcRepresentation *representation_mapped_to =
        kernel.representation_mapped_to(representation);
    if (representation_mapped_to)
    {
      // Check if this representation has (or will be) processed as part its mapped
      // representation
      representation_processed_as_mapped_item =
          ok_mapped_representations->contains(representation_mapped_to) ||
          reuse_ok_(settings, kernel.products_represented_by(representation_mapped_to));
    }
    if (representation_processed_as_mapped_item)
    {
      ok_mapped_representations->push(representation_mapped_to);
      // _nextShape();
      // continue;
      
      continue;
    }

    // Filter the products based on the set of entities and/or names being included or excluded
    // for processing.
    for (IfcSchema::IfcProduct::list::it jt = unfiltered_products->begin();
         jt != unfiltered_products->end(); ++jt)
    {
      IfcSchema::IfcProduct *prod = *jt;
      if (boost::all(filters_, filter_match(prod)))
      {
        ifcproducts->push(prod);
      }
    } // end for unfiltered_products

    for (ifcproduct_iterator = ifcproducts->begin(); ifcproduct_iterator != ifcproducts->end();
         ifcproduct_iterator++)
    {
      IfcproductRepresentation ir;
      ir.index=index_count;
      ir.product=*ifcproduct_iterator;
      ir.representation=representation;
      IfcproductRepresentations.push_back(ir);
      index_count++;
    } // end for ifcproducts

  } // for representation in representations
  end = std::chrono::system_clock::now();
  elapsed_seconds = end - start;
  Logger::Status("iterated over representations: " + std::to_string(elapsed_seconds.count()));

  return 1;
}
