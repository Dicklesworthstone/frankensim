//! Original phase-dependent transport on one fixed-geometry lumped body.
use crate::ConductionError;
use crate::radiation::{MaterialGraySurfaceLaw, MaterialGraySurfacePath};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::PropertyPathSlice;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::thermal::conduction::phase::PreparedPhaseConductivity;

/// Borrowed original conductivity and gray surface law on the SAME caloric
/// object. Preparation admits source support once; accepted cooling continues
/// to use the native radiation law at the actual endpoint enthalpy and phase.
/// No material data, thermal storage, geometry evolution or solver is added.
#[derive(Debug)]
pub struct LumpedPhaseThermalTransport<'a> {
    conductivity: &'a PreparedPhaseConductivity<'a>,
    radiation: &'a MaterialGraySurfaceLaw<'a>,
    radiation_path: MaterialGraySurfacePath,
    minimum_conductivity_w_per_m_k: f64,
    identity: ContentHash,
}

impl PartialEq for LumpedPhaseThermalTransport<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.identity == other.identity
    }
}
impl Eq for LumpedPhaseThermalTransport<'_> {}

impl<'a> LumpedPhaseThermalTransport<'a> {
    pub fn new(
        conductivity: &'a PreparedPhaseConductivity<'a>,
        radiation: &'a MaterialGraySurfaceLaw<'a>,
        mut checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<Self, ConductionError> {
        checkpoint()?;
        let invalid = |what| super::lumped_error("phase thermal transport", what);
        let curve = conductivity.caloric_curve();
        let knots = curve.knots();
        if !core::ptr::eq(curve, radiation.caloric_curve())
            || conductivity.enthalpy_bounds_j_kg()
                != [
                    knots[0].specific_enthalpy_j_kg,
                    knots[knots.len() - 1].specific_enthalpy_j_kg,
                ]
        {
            return Err(invalid(
                "phase conductivity and radiation must share the whole original caloric object"
                    .to_owned(),
            ));
        }
        let material = radiation.material_point();
        let context = &conductivity.source_paths()[0].start;
        if conductivity.temperature_axis() != material.temperature_axis
            || conductivity.liquid_fraction_axis() != material.liquid_fraction_axis
            || context.axes().keys().ne(material.context.axes().keys())
            || context.axis_quantities() != material.context.axis_quantities()
            || context.axes().iter().any(|(axis, value)| {
                axis != material.temperature_axis
                    && Some(axis.as_str()) != material.liquid_fraction_axis
                    && material.context.axes().get(axis) != Some(value)
            })
        {
            return Err(invalid(
                "phase conductivity and emissivity require the same typed source context"
                    .to_owned(),
            ));
        }
        let minimum_conductivity_w_per_m_k = conductivity
            .knots()
            .iter()
            .map(|point| point[1])
            .fold(f64::INFINITY, f64::min);
        if !minimum_conductivity_w_per_m_k.is_finite() || minimum_conductivity_w_per_m_k <= 0.0 {
            return Err(invalid(
                "an exposed lumped body needs positive conductivity throughout its chart"
                    .to_owned(),
            ));
        }
        let radiation_path = radiation.prepare_caloric_path(&mut checkpoint)?;
        let mut identity =
            DomainHasher::new("org.frankensim.fs-conduction.lumped-phase-thermal-transport.v1");
        identity.update(curve.identity().as_bytes());
        identity.update(conductivity.identity().as_bytes());
        identity.update(radiation_path.identity().as_bytes());
        checkpoint()?;
        Ok(Self {
            conductivity,
            radiation,
            radiation_path,
            minimum_conductivity_w_per_m_k,
            identity: identity.finalize(),
        })
    }

    pub fn caloric_curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.conductivity.caloric_curve()
    }
    pub const fn conductivity(&self) -> &'a PreparedPhaseConductivity<'a> {
        self.conductivity
    }
    pub const fn radiation(&self) -> &'a MaterialGraySurfaceLaw<'a> {
        self.radiation
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn minimum_conductivity_w_per_m_k(&self) -> f64 {
        self.minimum_conductivity_w_per_m_k
    }
    pub const fn maximum_emissivity(&self) -> f64 {
        self.radiation_path.maximum_emissivity()
    }
    pub fn emissivity_source_paths(&self) -> &[PropertyPathSlice] {
        self.radiation_path.source_paths()
    }

    /// Bound for inward radiation against specific enthalpy [W/(J/kg)]. The
    /// native thermal owner divides by mass and combines it with internal power.
    pub fn maximum_positive_radiative_feedback_w_per_j_kg(
        &self,
        area_m2: f64,
        reservoir_k: f64,
        checkpoint: impl FnMut() -> Result<(), ConductionError>,
    ) -> Result<f64, ConductionError> {
        self.radiation_path
            .maximum_positive_radiative_feedback_w_per_j_kg(area_m2, reservoir_k, checkpoint)
    }
}
