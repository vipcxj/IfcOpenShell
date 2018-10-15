#ifndef IFCGEOMCURVES_H

#define IFCGEOMCURVES_H
namespace IfcGeom
{
class IFC_GEOM_API Kernel
{
public:
bool convert(const IfcSchema::IfcCircle *l, Handle(Geom_Curve) & curve);
bool convert(const IfcSchema::IfcEllipse *l, Handle(Geom_Curve) & curve);
bool convert(const IfcSchema::IfcLine *l, Handle(Geom_Curve) & curve);

#ifdef USE_IFC4
  bool convert(const IfcSchema::IfcBSplineCurveWithKnots *l, Handle(Geom_Curve) & curve);
#endif

}; // end class
}

#endif
