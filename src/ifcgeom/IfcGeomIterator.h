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

#ifndef IFCGEOMITERATOR_H
#define IFCGEOMITERATOR_H

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

namespace IfcGeom {
	
template <typename P>
class Iterator {

 public:
  Iterator(const IteratorSettings& settings, IfcParse::IfcFile* file, std::vector<IfcGeom::filter_t>& filters)
      : settings_(settings)
      , ifc_file_(file)
      , owns_ifc_file_(false)
      , filters_(filters)
  {
    _initialize();
  }
  Iterator(const IteratorSettings& settings, IfcParse::IfcFile* file)
      : settings_(settings)
      , ifc_file_(file)
      , owns_ifc_file_(false)
  {
    _initialize();
  }
  Iterator(const IteratorSettings& settings, const std::string& filename)
      : settings_(settings)
      , ifc_file_(new IfcParse::IfcFile)
      , owns_ifc_file_(true)
  {
    ifc_file_->Init(filename);
    _initialize();
  }
  Iterator(const IteratorSettings& settings, void* data, int length)
      : settings_(settings)
      , ifc_file_(new IfcParse::IfcFile)
      , owns_ifc_file_(true)
  {
    ifc_file_->Init(data, length);
    _initialize();
  }
  Iterator(const IteratorSettings& settings, std::istream& filestream, int length)
      : settings_(settings)
      , ifc_file_(new IfcParse::IfcFile)
      , owns_ifc_file_(true)
  {
    ifc_file_->Init(filestream, length);
    _initialize();
  }
  ~Iterator() {
    if (owns_ifc_file_) {
      delete ifc_file_;
    }
    _free_shapes();
  }

  bool initialize() {
    try {
      _initUnits();
    } catch (const std::exception& e) {
      Logger::Error(e);
    }

    std::set<std::string> allowed_context_types;
    allowed_context_types.insert("model");
    allowed_context_types.insert("plan");
    allowed_context_types.insert("notdefined");
    std::set<std::string> context_types;

    if (!settings_.get(IteratorSettings::EXCLUDE_SOLIDS_AND_SURFACES)) {
      // Really this should only be 'Model', as per 
      // the standard 'Design' is deprecated. So,
      // just for backwards compatibility:
      context_types.insert("model");
      context_types.insert("design");
      // Some earlier (?) versions DDS-CAD output their own ContextTypes
      context_types.insert("model view");
      context_types.insert("detail view");
    }
    if (settings_.get(IteratorSettings::INCLUDE_CURVES)) {
      context_types.insert("plan");
    }			

    double lowest_precision_encountered = std::numeric_limits<double>::infinity();
    bool any_precision_encountered = false;

    representations_ = IfcSchema::IfcRepresentation::list::ptr(new IfcSchema::IfcRepresentation::list);

    IfcSchema::IfcGeometricRepresentationContext::list::it it;
    IfcSchema::IfcGeometricRepresentationSubContext::list::it jt;
    IfcSchema::IfcGeometricRepresentationContext::list::ptr contexts = 
        ifc_file_->entitiesByType<IfcSchema::IfcGeometricRepresentationContext>();

    IfcSchema::IfcGeometricRepresentationContext::list::ptr filtered_contexts (new IfcSchema::IfcGeometricRepresentationContext::list);

    for (it = contexts->begin(); it != contexts->end(); ++it) {
      IfcSchema::IfcGeometricRepresentationContext* context = *it;
      if (context->is(IfcSchema::Type::IfcGeometricRepresentationSubContext)) {
        // Continue, as the list of subcontexts will be considered
        // by the parent's context inverse attributes.
        continue;
      }
      try {
        if (context->hasContextType()) {
          std::string context_type = context->ContextType();
          boost::to_lower(context_type);

          if (allowed_context_types.find(context_type) == allowed_context_types.end()) {
            Logger::Message(Logger::LOG_ERROR, std::string("ContextType '") + context->ContextType() + "' not allowed:", context->entity);
          }
          if (context_types.find(context_type) != context_types.end()) {
            filtered_contexts->push(context);
          }
        }
      } catch (const std::exception& e) {
        Logger::Error(e);
      }
    }

    // In case no contexts are identified based on their ContextType, all contexts are
    // considered. Note that sub contexts are excluded as they are considered later on.
    if (filtered_contexts->size() == 0) {
      for (it = contexts->begin(); it != contexts->end(); ++it) {
        IfcSchema::IfcGeometricRepresentationContext* context = *it;
        if (!context->is(IfcSchema::Type::IfcGeometricRepresentationSubContext)) {
          filtered_contexts->push(context);
        }
      }
    }

    for (it = filtered_contexts->begin(); it != filtered_contexts->end(); ++it) {
      IfcSchema::IfcGeometricRepresentationContext* context = *it;
      representations_->push(context->RepresentationsInContext());
      try {
        if (context->hasPrecision() && context->Precision() < lowest_precision_encountered) {
          lowest_precision_encountered = context->Precision();
          any_precision_encountered = true;
        }
      } catch (const std::exception& e) {
        Logger::Error(e);
      }
      IfcSchema::IfcGeometricRepresentationSubContext::list::ptr sub_contexts = context->HasSubContexts();
      for (jt = sub_contexts->begin(); jt != sub_contexts->end(); ++jt) {
        representations_->push((*jt)->RepresentationsInContext());
      }
      // There is no need for full recursion as the following is governed by the schema:
      // WR31: The parent context shall not be another geometric representation sub context. 
    }
    
    if (any_precision_encountered) {
      // Some arbitrary factor that has proven to work better for the models in the set of test files.
      lowest_precision_encountered *= 10.;
      lowest_precision_encountered *= unit_magnitude_;
      if (lowest_precision_encountered < 1.e-7) {
        Logger::Message(Logger::LOG_WARNING, "Precision lower than 0.0000001 meter not enforced");
        kernel_.setValue(IfcGeom::Kernel::GV_PRECISION, 1.e-7);
      } else {
        kernel_.setValue(IfcGeom::Kernel::GV_PRECISION, lowest_precision_encountered);
      }
    } else {
      kernel_.setValue(IfcGeom::Kernel::GV_PRECISION, 1.e-5);
    }
    
    if (representations_->size() == 0) {
      Logger::Message(Logger::LOG_ERROR, "No geometries found");
      return false;
    }
    
    representation_iterator_ = representations_->begin();
    ifcproducts_.reset();
    if (!create()) {
      return false;
    }
    
    done_ = 0;
    total_ = representations_->size();
    for (int i = 1; i < 4; ++i) {
      bounds_min_.SetCoord(i, std::numeric_limits<double>::infinity());
      bounds_max_.SetCoord(i, -std::numeric_limits<double>::infinity());
    }
    
    IfcSchema::IfcProduct::list::ptr products = ifc_file_->entitiesByType<IfcSchema::IfcProduct>();
    for (IfcSchema::IfcProduct::list::it iter = products->begin(); iter != products->end(); ++iter) {
      IfcSchema::IfcProduct* product = *iter;
      if (product->hasObjectPlacement()) {
        // Use a fresh trsf every time in order to prevent the result to be concatenated
        gp_Trsf trsf; 
        bool success = false;
        try {
          success = kernel_.convert(product->ObjectPlacement(), trsf);
        } catch (const std::exception& e) {
          Logger::Error(e);
        } catch (...) {
          Logger::Error("Failed to construct placement");
        }
        if (!success) {
          continue;
        }
        const gp_XYZ& pos = trsf.TranslationPart();
        bounds_min_.SetX(std::min(bounds_min_.X(), pos.X()));
        bounds_min_.SetY(std::min(bounds_min_.Y(), pos.Y()));
        bounds_min_.SetZ(std::min(bounds_min_.Z(), pos.Z()));
        bounds_max_.SetX(std::max(bounds_max_.X(), pos.X()));
        bounds_max_.SetY(std::max(bounds_max_.Y(), pos.Y()));
        bounds_max_.SetZ(std::max(bounds_max_.Z(), pos.Z()));
      }
    }
    return true;
  } // end Iterator::initialize

  int progress() const { return 100 * done_ / total_; }

  const std::string& getUnitName() const { return unit_name_; }

  /// @note Double always as per IFC specification.
  double getUnitMagnitude() const { return unit_magnitude_; }
	
  std::string getLog() const { return Logger::GetLog(); }

  IfcParse::IfcFile* getFile() const { return ifc_file_; }

  const std::vector<IfcGeom::filter_t>& filters() const { return filters_; }
  std::vector<IfcGeom::filter_t>& filters() { return filters_; }

  const gp_XYZ& bounds_min() const { return bounds_min_; }

  const gp_XYZ& bounds_max() const { return bounds_max_; }

  /// Returns what would be the product for the next shape representation
  /// @todo Double-check and test the impl.
  //IfcSchema::IfcProduct* peek_next() const
  //{
  //    if (ifcproducts && ifcproduct_iterator_ + 1 != ifcproducts_->end()){
  //        return *(ifcproduct_iterator_ + 1);
  //    } else {
  //        return 0;
  //    }
  //}

  /// @todo Would this be as simple as the following code?
  //void skip_next() { if (ifcproducts) { ++ifcproduct_iterator_; } }

  /// Moves to the next shape representation, create its geometry, and returns the associated product.
  /// Use get() to retrieve the created geometry.
  IfcSchema::IfcProduct* next() {
    // Increment the iterator over the list of products using the current
    // shape representation
    if (ifcproducts_) {
      ++ifcproduct_iterator_;
    }

    return create();
  }

  /// Gets the representation of the current geometrical entity.
  Element<P>* get()
  {
    // TODO: Test settings and throw
    Element<P>* ret = 0;
    if (current_triangulation_) { ret = current_triangulation_; }
    else if (current_serialization_) { ret = current_serialization_; }
    else if (current_shape_model_) { ret = current_shape_model_; }
    // If we want to organize the element considering their hierarchy
    if (settings_.get(IteratorSettings::SEARCH_FLOOR))
    {
      // We are going to build a vector with the element parents.
      // First, create the parent vector
      std::vector<const IfcGeom::Element<P>*> parents;
      // if the element has a parent
      if (ret->parent_id() != -1)
      {
        const IfcGeom::Element<P>* parent_object = NULL;
        bool hasParent = true;
        // get the parent 
        try {
          parent_object = getObject(ret->parent_id());
        } catch (const std::exception& e) {
          Logger::Error(e);
          hasParent = false;
        }
        // Add the previously found parent to the vector
        if (hasParent) parents.insert(parents.begin(), parent_object);
        // We need to find all the parents
        while (parent_object != NULL && hasParent && parent_object->parent_id() != -1)
        {
          // Find the next parent
          try {
            parent_object = getObject(parent_object->parent_id());
          } catch (const std::exception& e) {
            Logger::Error(e);
            hasParent = false;
          }
          // Add the previously found parent to the vector
          if (hasParent) parents.insert(parents.begin(), parent_object);
          hasParent = hasParent && parent_object->parent_id() != -1;
        }
        // when done_ push the parent list in the Element object
        ret->SetParents(parents);
      }
    }
    return ret;
  }

  /// Gets the native (Open Cascade) representation of the current geometrical entity.
  BRepElement<P>* get_native()
  {
    // TODO: Test settings and throw
    return current_shape_model_;
  }

  const Element<P>* getObject(int id) {
    gp_Trsf trsf;
    int parent_id = -1;
    std::string instance_type, product_name, product_guid;
    IfcSchema::IfcProduct* ifc_product = 0;

    try {
      IfcUtil::IfcBaseClass* ifc_entity = ifc_file_->entityById(id);
      instance_type = IfcSchema::Type::ToString(ifc_entity->type());

      if (ifc_entity->is(IfcSchema::Type::IfcRoot)) {
        IfcSchema::IfcRoot* ifc_root = ifc_entity->as<IfcSchema::IfcRoot>();
        product_guid = ifc_root->GlobalId();
        product_name = ifc_root->hasName() ? ifc_root->Name() : "";
      }

      if (ifc_entity->is(IfcSchema::Type::IfcProduct)) {
        ifc_product = ifc_entity->as<IfcSchema::IfcProduct>();
        parent_id = -1;
        try {
          IfcSchema::IfcObjectDefinition* parent_object = kernel_.get_decomposing_entity(ifc_product);
          if (parent_object) {
            parent_id = parent_object->entity->id();
          }
        } catch (const std::exception& e) {
          Logger::Error(e);
        } catch (...) {
          Logger::Error("Failed to find decomposing entity");
        }

        try {
          kernel_.convert(ifc_product->ObjectPlacement(), trsf);
        } catch (const std::exception& e) {
          Logger::Error(e);
        } catch (...) {
          Logger::Error("Failed to construct placement");
        }
      }
    } catch (const std::exception& e) {
      Logger::Error(e);
    } catch (const Standard_Failure& e) {
      if (e.GetMessageString() && strlen(e.GetMessageString())) {
        Logger::Error(e.GetMessageString());
      } else {
        Logger::Error("Unknown error returning product");
      }
    } catch (...) {
      Logger::Error("Unknown error returning product");
    }

    ElementSettings element_settings(settings_, unit_magnitude_, instance_type);

    Element<P>* ifc_object = new Element<P>(element_settings, id, parent_id, product_name, instance_type, product_guid, "", trsf, ifc_product);
    return ifc_object;
  }

  IfcSchema::IfcProduct* create() {
    IfcGeom::BRepElement<P>* next_shape_model = 0;
    IfcGeom::SerializedElement<P>* next_serialization = 0;
    IfcGeom::TriangulationElement<P>* next_triangulation = 0;

    try {
      next_shape_model = _create_shape_model_for_next_entity();
    } catch (const std::exception& e) {
      Logger::Error(e);
    } catch (const Standard_Failure& e) {
      if (e.GetMessageString() && strlen(e.GetMessageString())) {
        Logger::Error(e.GetMessageString());
      } else {
        Logger::Error("Unknown error creating geometry");
      }
    } catch (...) {
      Logger::Error("Unknown error creating geometry");
    }

    if (next_shape_model) {
      if (settings_.get(IteratorSettings::USE_BREP_DATA)) {
        try {
          next_serialization = new SerializedElement<P>(*next_shape_model);
        } catch (...) {
          Logger::Message(Logger::LOG_ERROR, "Getting a serialized element from model failed.");
        }
      } else if (!settings_.get(IteratorSettings::DISABLE_TRIANGULATION)) {
        try {
          if (ifcproduct_iterator_ == ifcproducts_->begin() || !geometry_reuse_ok_for_current_representation_) {
            next_triangulation = new TriangulationElement<P>(*next_shape_model);
          } else {
            next_triangulation = new TriangulationElement<P>(*next_shape_model, current_triangulation_->geometry_pointer());
          }
        } catch (...) {
          Logger::Message(Logger::LOG_ERROR, "Getting a triangulation element from model failed.");
        }
      }
    }

    _free_shapes();

    current_shape_model_ = next_shape_model;
    current_serialization_ = next_serialization;
    current_triangulation_ = next_triangulation;

    return next_shape_model ? next_shape_model->product() : 0;
  }

 private:
  IteratorSettings settings_;
  Iterator(const Iterator&); // N/I
  Iterator& operator=(const Iterator&); // N/I
  Kernel kernel_;
  IfcParse::IfcFile* ifc_file_;
  // A container and iterator for IfcRepresentations
  IfcSchema::IfcRepresentation::list::ptr representations_;
  IfcSchema::IfcRepresentation::list::it representation_iterator_;
  // The object is fetched beforehand to be sure that get() returns a valid element
  TriangulationElement<P>* current_triangulation_;
  BRepElement<P>* current_shape_model_;
  SerializedElement<P>* current_serialization_;
  // A container and iterator for IfcBuildingElements for the current IfcRepresentation referenced by *representation_iterator_
  IfcSchema::IfcProduct::list::ptr ifcproducts_;
  IfcSchema::IfcProduct::list::it ifcproduct_iterator_;

  int done_;
  int total_;
  bool owns_ifc_file_;
  std::string unit_name_;
  double unit_magnitude_;
  gp_XYZ bounds_min_;
  gp_XYZ bounds_max_;
  std::vector<filter_t> filters_;
  
  bool geometry_reuse_ok_for_current_representation_;
  
  struct filter_match_
  {
    filter_match_(IfcSchema::IfcProduct *prod) : product(prod) {}
    bool operator()(const filter_t& filter) const { return filter(product);  }

    IfcSchema::IfcProduct* product;
  };

  void _initUnits() {
    IfcSchema::IfcProject::list::ptr projects = ifc_file_->entitiesByType<IfcSchema::IfcProject>();
    if (projects->size() == 1) {
      IfcSchema::IfcProject* project = *projects->begin();
      std::pair<std::string, double> length_unit = kernel_.initializeUnits(project->UnitsInContext());
      unit_name_ = length_unit.first;
      unit_magnitude_ = length_unit.second;
    }
  }

  // Move to the next IfcRepresentation
  void _nextShape() {
    // In order to conserve memory and reduce cache insertion times, the cache is
    // cleared after an arbitary number of processed representations_. This has been
    // benchmarked extensively: https://github.com/IfcOpenShell/IfcOpenShell/pull/47
    static const int clear_interval = 64;
    if (done_ % clear_interval == clear_interval - 1) {
      kernel_.purge_cache();
    }
    ifcproducts_.reset();
    ++ representation_iterator_;
    ++ done_;
  }

  bool _reuse_ok(const IfcSchema::IfcProduct::list::ptr& products) {
    // With world coords enabled, object transformations are directly applied to
    // the BRep. There is no way to re-use the geometry for multiple products.
    if (settings_.get(IteratorSettings::USE_WORLD_COORDS)) {
      return false;
    }

    std::set<const IfcSchema::IfcMaterial*> associated_single_materials;

    for (IfcSchema::IfcProduct::list::it it = products->begin(); it != products->end(); ++it) {
      IfcSchema::IfcProduct* product = *it;

      if (!settings_.get(IteratorSettings::DISABLE_OPENING_SUBTRACTIONS) && kernel_.find_openings(product)->size()) {
        return false;
      }

      if (settings_.get(IteratorSettings::APPLY_LAYERSETS)) {
        IfcSchema::IfcRelAssociates::list::ptr associations = product->HasAssociations();
        for (IfcSchema::IfcRelAssociates::list::it jt = associations->begin(); jt != associations->end(); ++jt) {
          IfcSchema::IfcRelAssociatesMaterial* assoc = (*jt)->as<IfcSchema::IfcRelAssociatesMaterial>();
          if (assoc) {
            if (assoc->RelatingMaterial()->is(IfcSchema::Type::IfcMaterialLayerSetUsage)) {
              // TODO: Check whether single layer? 
              return false;
            }
          }
        }
      }

      // Note that this can be a nullptr (!), but the fact that set size should be one still holds
      associated_single_materials.insert(kernel_.get_single_material_association(product));
    }

    return associated_single_materials.size() == 1;
  }

  BRepElement<P>* _create_shape_model_for_next_entity() {
    for (;;) {
      IfcSchema::IfcRepresentation* representation;

      // Have we reached the end of our list of representations_?
      if ( representation_iterator_ == representations_->end() ) {
        representations_.reset();
        return 0;
      }
      representation = *representation_iterator_;

      // Has the list of IfcProducts for this representation been initialized?
      if (!ifcproducts_) {
        ifcproducts_ = IfcSchema::IfcProduct::list::ptr(new IfcSchema::IfcProduct::list);
        IfcSchema::IfcProduct::list::ptr unfiltered_products = kernel_.products_represented_by(representation);

        geometry_reuse_ok_for_current_representation_ = _reuse_ok(unfiltered_products);

        if (!geometry_reuse_ok_for_current_representation_ && representation->RepresentationMap()->size() == 1) {
          // unfiltered_products contains products represented by this representation by means of mapped items.
          // For example because of openings applied to products, reuse might not be acceptable and then the
          // products will be processed by means of their immediate representation and not the mapped representation.
          _nextShape();
          continue;
        }

        bool representation_processed_as_mapped_item = false;

        IfcSchema::IfcRepresentation* representation_mapped_to = kernel_.representation_mapped_to(representation);
        if (representation_mapped_to) {
          // Check if this represenation has (or will be) processed as part its mapped representation
          representation_processed_as_mapped_item = _reuse_ok(kernel_.products_represented_by(representation_mapped_to));
        }

        if (representation_processed_as_mapped_item) {
          _nextShape();
          continue;
        }

        // Filter the products based on the set of entities and/or names being included or excluded for processing.
        for (IfcSchema::IfcProduct::list::it jt = unfiltered_products->begin(); jt != unfiltered_products->end(); ++jt) {
          IfcSchema::IfcProduct* prod = *jt;
          if (boost::all(filters_, filter_match_(prod))) {
            ifcproducts_->push(prod);
          }
        }

        ifcproduct_iterator_ = ifcproducts_->begin();
      }

      // Have we reached the end of our list of IfcProducts?
      if ( ifcproduct_iterator_ == ifcproducts_->end() ) {
        _nextShape();
        continue;
      }

      IfcSchema::IfcProduct* product = *ifcproduct_iterator_;
      Logger::SetProduct(product);

      BRepElement<P>* element;
      if (ifcproduct_iterator_ == ifcproducts_->begin() || !geometry_reuse_ok_for_current_representation_) {
        element = kernel_.create_brep_for_representation_and_product<P>(settings_, representation, product);
      } else {
        element = kernel_.create_brep_for_processed_representation(settings_, representation, product, current_shape_model_);
      }

      Logger::SetProduct(boost::none);

      if (!element) {
        _nextShape();
        continue;
      }

      return element;
    }
  }

  void _free_shapes() {
    // Free all possible representations_ of the current geometrical entity
    delete current_triangulation_;
    current_triangulation_ = 0;
    delete current_serialization_;
    current_serialization_ = 0;
    delete current_shape_model_;
    current_shape_model_ = 0;
  }

  void _initialize() {
    current_triangulation_ = 0;
    current_shape_model_ = 0;
    current_serialization_ = 0;
    unit_name_ = "METER";
    unit_magnitude_ = 1.f;
    kernel_.setValue(IfcGeom::Kernel::GV_MAX_FACES_TO_SEW,
                    settings_.get(IteratorSettings::SEW_SHELLS) ? 1000 : -1);
    kernel_.setValue(IfcGeom::Kernel::GV_DIMENSIONALITY,
                    (settings_.get(IteratorSettings::INCLUDE_CURVES)
                                                         ? (settings_.get(IteratorSettings::EXCLUDE_SOLIDS_AND_SURFACES) ? -1. : 0.) : +1.));
    if (settings_.get(IteratorSettings::SITE_LOCAL_PLACEMENT)) {
      kernel_.set_conversion_placement_rel_to(IfcSchema::Type::IfcSite);
    }
  } //end Iterator::_initialize
  
}; //end class Iterator
} // end namespace IfcGeom

#endif
