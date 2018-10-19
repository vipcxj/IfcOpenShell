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

#include "IfcGeom.h"

#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>
#include <gp_Dir.hxx>
#include <gp_GTrsf.hxx>
#include <gp_GTrsf2d.hxx>
#include <gp_Mat.hxx>
#include <gp_Mat2d.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Trsf2d.hxx>

#include "../ifcparse/IfcBaseClass.h"
#include "../ifcparse/IfcParse.h"


// #include "IfcGeomShapeType.h" //moved definition to ifcgeom.h, it defines ST_SHAPELIST

using namespace IfcSchema;
using namespace IfcUtil;

bool IfcGeom::Kernel::convert_shapes(const IfcBaseClass *l, IfcRepresentationShapeItems &r)
{
  if (shape_type(l) != ST_SHAPELIST)
  {
    TopoDS_Shape shp;
    if (convert_shape(l, shp))
    {
      r.push_back(IfcGeom::IfcRepresentationShapeItem(
          shp, get_style(l->as<IfcSchema::IfcRepresentationItem>())));
      return true;
    }
    return false;
  }

#define SHAPES(T)                                                                                  \
  if (l->is(T::Class()))                                                                           \
  {                                                                                                \
    try                                                                                            \
    {                                                                                              \
      return convert((T *)l, r);                                                                   \
    }                                                                                              \
    catch (const std::exception &e)                                                                \
    {                                                                                              \
      Logger::Message(Logger::LOG_ERROR,                                                           \
                      std::string(e.what()) + "\nFailed to convert:", l->entity);                  \
    }                                                                                              \
    catch (const Standard_Failure &f)                                                              \
    {                                                                                              \
      if (f.GetMessageString())                                                                    \
        Logger::Message(                                                                           \
            Logger::LOG_ERROR,                                                                     \
            std::string("Error in: ") + f.GetMessageString() + "\nFailed to convert:", l->entity); \
      else                                                                                         \
        Logger::Message(Logger::LOG_ERROR, "Failed to convert:", l->entity);                       \
    }                                                                                              \
    return false;                                                                                  \
  }

  SHAPES(IfcShellBasedSurfaceModel);
  SHAPES(IfcFaceBasedSurfaceModel);
  SHAPES(IfcRepresentation);
  SHAPES(IfcMappedItem);
  // IfcFacetedBrep included
  // IfcAdvancedBrep included
  // IfcFacetedBrepWithVoids included
  // IfcAdvancedBrepWithVoids included
  SHAPES(IfcManifoldSolidBrep);
  SHAPES(IfcGeometricSet);
  
  Logger::Message(Logger::LOG_ERROR, "No operation defined for:", l->entity);
  return false;
}

IfcGeom::ShapeType IfcGeom::Kernel::shape_type(const IfcBaseClass *l)
{
  /////////////////////////////////////////////////////////
  // continue:

  // #define SHAPES(T)                                                                                  \
//   if (l->is(T::Class()))                                                                           \
//     return ST_SHAPELIST;

  std::vector<Type::Enum> shapes_type_list = {
      Type::IfcShellBasedSurfaceModel, Type::IfcFaceBasedSurfaceModel,
      Type::IfcRepresentation,         Type::IfcMappedItem,
      Type::IfcManifoldSolidBrep,      Type::IfcGeometricSet};

  for (Type::Enum T : shapes_type_list)
  {
    if (l->is(T))
      return ST_SHAPELIST;
  }

  // #define SHAPE(T)                                                                                   \
//   if (l->is(T::Class()))                                                                           \
//     return ST_SHAPE;

  std::vector<Type::Enum> shape_list = {Type::IfcPlane,
                                        Type::IfcExtrudedAreaSolid,
                                        Type::IfcRevolvedAreaSolid,
                                        Type::IfcConnectedFaceSet,
                                        Type::IfcBooleanResult,
                                        Type::IfcPolygonalBoundedHalfSpace,
                                        Type::IfcHalfSpaceSolid,
                                        Type::IfcSurfaceOfLinearExtrusion,
                                        Type::IfcSurfaceOfRevolution,
                                        Type::IfcBlock,
                                        Type::IfcRectangularPyramid,
                                        Type::IfcRightCircularCylinder,
                                        Type::IfcRightCircularCone,
                                        Type::IfcSphere,
                                        Type::IfcCsgSolid,
                                        Type::IfcCurveBoundedPlane,
                                        Type::IfcRectangularTrimmedSurface,
                                        Type::IfcSurfaceCurveSweptAreaSolid,
                                        Type::IfcSweptDiskSolid};

  std::vector<Type::Enum> shape_list_ifc4 = {
      Type::IfcCylindricalSurface, Type::IfcAdvancedBrep, Type::IfcBSplineSurfaceWithKnots,
      Type::IfcTriangulatedFaceSet, Type::IfcExtrudedAreaSolidTapered};

#ifdef USE_IFC4
  shape_list.insert(shape_list.end(), shape_list_ifc4.begin(), shape_list_ifc4.end());
#endif

  for (Type::Enum T : shape_list)
  {
    if (l->is(T))
      return ST_SHAPE;
  }

  // #define FACE(T)                                                                                    \
//   if (l->is(T::Class()))                                                                           \
//     return ST_FACE;

  std::vector<Type::Enum> face_list = {Type::IfcArbitraryProfileDefWithVoids,
                                       Type::IfcArbitraryClosedProfileDef,
                                       Type::IfcRoundedRectangleProfileDef,
                                       Type::IfcRectangleHollowProfileDef,
                                       Type::IfcRectangleProfileDef,
                                       Type::IfcTrapeziumProfileDef,
                                       Type::IfcCShapeProfileDef,
                                       Type::IfcIShapeProfileDef,
                                       Type::IfcLShapeProfileDef,
                                       Type::IfcTShapeProfileDef,
                                       Type::IfcUShapeProfileDef,
                                       Type::IfcZShapeProfileDef,
                                       Type::IfcCircleHollowProfileDef,
                                       Type::IfcCircleProfileDef,
                                       Type::IfcEllipseProfileDef,
                                       Type::IfcCenterLineProfileDef,
                                       Type::IfcCompositeProfileDef,
                                       Type::IfcDerivedProfileDef,
                                       Type::IfcFace};

  for (Type::Enum T : face_list)
  {
    if (l->is(T))
      return ST_FACE;
  }

  // #define WIRE(T)                                                                                    \
//   if (l->is(T::Class()))                                                                           \
//     return ST_WIRE;
  std::vector<Type::Enum> wire_list = {Type::IfcEdgeCurve,    Type::IfcSubedge,
                                       Type::IfcOrientedEdge, Type::IfcEdge,
                                       Type::IfcEdgeLoop,     Type::IfcPolyline,
                                       Type::IfcPolyLoop,     Type::IfcCompositeCurve,
                                       Type::IfcTrimmedCurve, Type::IfcArbitraryOpenProfileDef};
#ifdef USE_IFC4
  wire_list.emplace_back(Type::IfcIndexedPolyCurve);
#endif
  for (Type::Enum T : wire_list)
  {
    if (l->is(T))
      return ST_WIRE;
  }

  // #define CURVE(T)                                                                                   \
//   if (l->is(T::Class()))                                                                           \
//     return ST_CURVE;

  std::vector<Type::Enum> curve_list{Type::IfcCircle, Type::IfcEllipse, Type::IfcLine};
#ifdef USE_IFC4
  curve_list.emplace_back(Type::IfcBSplineCurveWithKnots);
#endif
  for (Type::Enum T : curve_list)
  {
    if (l->is(T))
      return ST_CURVE;
  }

  return ST_OTHER;
}

bool IfcGeom::Kernel::convert_shape(const IfcBaseClass *l, TopoDS_Shape &r)
{
  const unsigned int id = l->entity->id();
  bool success = false;
  bool processed = false;
  bool ignored = false;

  // #ifndef NO_CACHE                                                                                \
    //   std::map<int, TopoDS_Shape>::const_iterator it = cache.Shape.find(id);                    \
    //   if (it != cache.Shape.end())                                                              \
    //   {                                                                                         \
    //     r = it->second;                                                                         \
    //     return true;                                                                            \
    //   }                                                                                         \
    // #endif

  const bool include_curves = getValue(GV_DIMENSIONALITY) != +1;
  const bool include_solids_and_surfaces = getValue(GV_DIMENSIONALITY) != -1;

  IfcGeom::ShapeType st = shape_type(l);
  ignored = (!include_solids_and_surfaces && (st == ST_SHAPE || st == ST_FACE)) ||
            (!include_curves && (st == ST_WIRE || st == ST_CURVE));
  if (st == ST_SHAPELIST)
  {
    processed = true;
    IfcRepresentationShapeItems items;
    success = convert_shapes(l, items) && flatten_shape_list(items, r, false);
  }
  else if (st == ST_SHAPE && include_solids_and_surfaces)
  {
    //#include "IfcRegisterConvertShape.h" // Dude, no.
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
#define SHAPE(T)                                                                                   \
  if (!processed && l->is(T::Class()))                                                             \
  {                                                                                                \
    processed = true;                                                                              \
    try                                                                                            \
    {                                                                                              \
      if (convert((T *)l, r))                                                                      \
      {                                                                                            \
        success = true;                                                                            \
      }                                                                                            \
    }                                                                                              \
    catch (const std::exception &e)                                                                \
    {                                                                                              \
      Logger::Message(Logger::LOG_ERROR,                                                           \
                      std::string(e.what()) + "\nFailed to convert:", l->entity);                  \
      return false;                                                                                \
    }                                                                                              \
    catch (const Standard_Failure &f)                                                              \
    {                                                                                              \
      if (f.GetMessageString() && strlen(f.GetMessageString()))                                    \
        Logger::Message(                                                                           \
            Logger::LOG_ERROR,                                                                     \
            std::string("Error in: ") + f.GetMessageString() + "\nFailed to convert:", l->entity); \
      else                                                                                         \
        Logger::Message(Logger::LOG_ERROR, "Failed to convert:", l->entity);                       \
      return false;                                                                                \
    }                                                                                              \
    if (!success)                                                                                  \
    {                                                                                              \
      Logger::Message(Logger::LOG_ERROR, "Failed to convert:", l->entity);                         \
      return false;                                                                                \
    }                                                                                              \
  }
#ifdef USE_IFC4
    SHAPE(IfcCylindricalSurface);
    SHAPE(IfcAdvancedBrep);
    // FIXME: Surfaces should have a shape type of their own
    SHAPE(IfcBSplineSurfaceWithKnots);
    SHAPE(IfcTriangulatedFaceSet);
    SHAPE(IfcExtrudedAreaSolidTapered);
#endif
    SHAPE(IfcPlane);
    SHAPE(IfcExtrudedAreaSolid);
    SHAPE(IfcRevolvedAreaSolid);
    SHAPE(IfcConnectedFaceSet);
    SHAPE(IfcBooleanResult);
    SHAPE(IfcPolygonalBoundedHalfSpace);
    SHAPE(IfcHalfSpaceSolid);
    // FIXME: Surfaces should have a shape type of their own
    SHAPE(IfcSurfaceOfLinearExtrusion);
    SHAPE(IfcSurfaceOfRevolution);
    SHAPE(IfcBlock);
    SHAPE(IfcRectangularPyramid);
    SHAPE(IfcRightCircularCylinder);
    SHAPE(IfcRightCircularCone);
    SHAPE(IfcSphere);
    SHAPE(IfcCsgSolid);
    SHAPE(IfcCurveBoundedPlane);
    SHAPE(IfcRectangularTrimmedSurface);
    SHAPE(IfcSurfaceCurveSweptAreaSolid);
    SHAPE(IfcSweptDiskSolid);
    /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////
  }
  else if (st == ST_FACE && include_solids_and_surfaces)
  {
    processed = true;
    success = convert_face(l, r);
  }
  else if (st == ST_WIRE && include_curves)
  {
    processed = true;
    TopoDS_Wire w;
    success = convert_wire(l, w);
    if (success)
    {
      r = w;
    }
  }
  else if (st == ST_CURVE && include_curves)
  {
    processed = true;
    Handle(Geom_Curve) crv;
    TopoDS_Wire w;
    success = convert_curve(l, crv) && convert_curve_to_wire(crv, w);
    if (success)
    {
      r = w;
    }
  }

  if (processed && success)
  {
    const double precision = getValue(GV_PRECISION);
    apply_tolerance(r, precision);
    // #ifndef NO_CACHE
    //     cache.Shape[id] = r;
    // #endif
  }
  else if (!ignored)
  {
    const char *const msg = processed ? "Failed to convert:" : "No operation defined for:";
    Logger::Message(Logger::LOG_ERROR, msg, l->entity);
  }
  return success;
}

bool IfcGeom::Kernel::convert_wire(const IfcBaseClass *l, TopoDS_Wire &r)
{
  //#include "IfcRegisterConvertWire.h"
#define WIRE(T)                                                                                    \
  if (l->is(T::Class()))                                                                           \
    return convert((T *)l, r);

  WIRE(IfcEdgeCurve);
  WIRE(IfcSubedge);
  WIRE(IfcOrientedEdge);
  WIRE(IfcEdge);
  WIRE(IfcEdgeLoop);
  WIRE(IfcPolyline);
  WIRE(IfcPolyLoop);
  WIRE(IfcCompositeCurve);
  WIRE(IfcTrimmedCurve);
  WIRE(IfcArbitraryOpenProfileDef);
#ifdef USE_IFC4
  WIRE(IfcIndexedPolyCurve)
#endif

  Handle(Geom_Curve) curve;
  if (IfcGeom::Kernel::convert_curve(l, curve))
  {
    return IfcGeom::Kernel::convert_curve_to_wire(curve, r);
  }
  Logger::Message(Logger::LOG_ERROR, "No operation defined for:", l->entity);
  return false;
}

bool IfcGeom::Kernel::convert_face(const IfcBaseClass *l, TopoDS_Shape &r)
{

#define FACE(T)                                                                                    \
  if (l->is(T::Class()))                                                                           \
    return convert((T *)l, r);

  FACE(IfcArbitraryProfileDefWithVoids);
  FACE(IfcArbitraryClosedProfileDef);
  FACE(IfcRoundedRectangleProfileDef);
  FACE(IfcRectangleHollowProfileDef);
  FACE(IfcRectangleProfileDef);
  FACE(IfcTrapeziumProfileDef)
  FACE(IfcCShapeProfileDef);
  // IfcAsymmetricIShapeProfileDef included
  FACE(IfcIShapeProfileDef);
  FACE(IfcLShapeProfileDef);
  FACE(IfcTShapeProfileDef);
  FACE(IfcUShapeProfileDef);
  FACE(IfcZShapeProfileDef);
  FACE(IfcCircleHollowProfileDef);
  FACE(IfcCircleProfileDef);
  FACE(IfcEllipseProfileDef);
  FACE(IfcCenterLineProfileDef);
  FACE(IfcCompositeProfileDef);
  FACE(IfcDerivedProfileDef);
  // IfcFaceSurface included
  // IfcAdvancedFace included in case of IFC4
  FACE(IfcFace);

  Logger::Message(Logger::LOG_ERROR, "No operation defined for:", l->entity);
  return false;
}

bool IfcGeom::Kernel::convert_curve(const IfcBaseClass *l, Handle(Geom_Curve) & r)
{
  //#include "IfcRegisterConvertCurve.h"
  /*
#include "IfcRegisterUndef.h"
#define CURVE(T) \
        if ( l->is(T::Class()) ) return convert((T*)l,r);
#include "IfcRegisterDef.h"

#include "IfcRegister.h"
*/
#define CURVE(T)                                                                                   \
  if (l->is(T::Class()))                                                                           \
    return convert((T *)l, r);

  CURVE(IfcCircle);
  CURVE(IfcEllipse);
  CURVE(IfcLine);
#ifdef USE_IFC4
  // IfcRationalBSplineCurveWithKnots included
  CURVE(IfcBSplineCurveWithKnots);
#endif

  Logger::Message(Logger::LOG_ERROR, "No operation defined for:", l->entity);
  return false;
}
