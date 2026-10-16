//! Compile only the exact selected immutable claims. Zero strain is a law-
//! selection sample, not a claim that the evolving specimen is unstrained.
use super::{StrainCurve, TabulatedMembrane};
use crate::sheet::{SheetError, MEMBRANE_DIMS, AREAL_DENSITY, THICKNESS, property};
use crate::state_point::ResolvedMaterialStatePoint;
use fs_matdb::{MaterialCard, InterpolationPolicy, PropertyValue};
use fs_qty::{Dims, QuantitySpec};

/// Constitutive lookup coordinate, dimensionless Green strain.
pub const STRAIN_AXIS: &str = "sheet_green_strain";
/// Work-conjugate stress tables; keys identify direction AND stress measure.
pub const CURVE_PROPERTIES: [&str;3] = [
    "sheet_second_piola_warp", "sheet_second_piola_weft", "sheet_second_piola_shear",
];
/// Constant cross derivative dS11/dE22 = dS22/dE11 [N/m].
pub const COUPLING_PROPERTY: &str = "sheet_green_strain_coupling";

impl TabulatedMembrane {
    /// Compile the selected original curves without refitting. Only explicit
    /// linear-inside interpolation is admitted, never tabulated-only samples.
    pub fn from_material_source(material:&ResolvedMaterialStatePoint, card:&MaterialCard)
        -> Result<Self,SheetError> {
        if material.card_identity()!=card.content_hash()
            || material.query_point().iter().find(|(n,_)|n==STRAIN_AXIS).map(|(_,v)|*v)!=Some(0.0)
            || material.axis_quantities().get(STRAIN_AXIS)
                .is_some_and(|q|*q!=QuantitySpec::dimensional(Dims::NONE)) {
            return Err(SheetError::Invalid("nonlinear sheet needs its original card and dimensionless zero-strain selection sample"));
        }
        for key in [AREAL_DENSITY, THICKNESS, COUPLING_PROPERTY] {
            let p=material.property(key).ok_or(SheetError::Invalid(key))?;
            let c=card.claims().claim(p.answer().receipt.selected).ok_or(SheetError::Invalid("missing selected sheet claim"))?;
            if c.validity.bounds().contains_key(STRAIN_AXIS) {
                return Err(SheetError::Invalid("reference mass, thickness and constant coupling cannot depend on the evolving strain"));
            }
        }
        let coupling=property(material,COUPLING_PROPERTY,MEMBRANE_DIMS)?;
        let compile=|key:&'static str| -> Result<StrainCurve,SheetError> {
            let p=material.property(key).ok_or(SheetError::Invalid(key))?;
            if p.requirement().quantity()!=QuantitySpec::dimensional(MEMBRANE_DIMS) {
                return Err(SheetError::Invalid("second-Piola sheet resultant must have N/m units"));
            }
            let claim=card.claims().claim(p.answer().receipt.selected)
                .ok_or(SheetError::Invalid("missing selected sheet curve"))?;
            let (abscissa,abscissa_dims,knots,dims)=match &claim.value {
                PropertyValue::Curve {abscissa,abscissa_dims,knots,dims}=>(abscissa,abscissa_dims,knots,dims),
                _=>return Err(SheetError::Invalid("nonlinear sheet requires actual stress/strain curves")),
            };
            if abscissa!=STRAIN_AXIS || *abscissa_dims!=Dims::NONE || *dims!=MEMBRANE_DIMS
                || claim.interpolation!=InterpolationPolicy::LinearInside {
                return Err(SheetError::Invalid("nonlinear sheet requires linear-inside Green-strain curves, not sampled or temperature curves"));
            }
            let mut curve=StrainCurve::new(knots)?;
            if let Some(&(lo,hi))=claim.validity.bounds().get(STRAIN_AXIS) {
                curve.support=[curve.support[0].max(lo),curve.support[1].min(hi)];
                if !lo.is_finite() || !hi.is_finite() || curve.support[0]>=0.0 || curve.support[1]<=0.0 {
                    return Err(SheetError::Invalid("source support must contain a neighborhood of the stress-free origin"));
                }
            }
            Ok(curve)
        };
        Self::new([compile(CURVE_PROPERTIES[0])?,compile(CURVE_PROPERTIES[1])?,compile(CURVE_PROPERTIES[2])?],coupling)
    }
}
