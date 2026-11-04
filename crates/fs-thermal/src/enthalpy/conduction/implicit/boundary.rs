//! Complete conductive surface fluxes on the original implicit heat graph.

use super::{ImplicitConductionChain, ImplicitConductionError};
use crate::enthalpy::{CellId, DepositionError, DepositionSimulation, ThermalLaw};
use std::collections::BTreeMap;

/// Original centre-to-face geometry and the external film/contact resistance.
/// Area is the actual reservoir coverage or exposed ambient area. Resistance
/// is area-specific; a material owner converts it to total K/W exactly once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitBoundaryGeometry {
    pub area_m2: f64,
    pub half_distance_m: f64,
    pub surface_resistance_m2_k_w: f64,
}

/// Passive outward heat flux and its complete original-enthalpy derivative.
/// The native owner applies `Q_out = G * (T_cell - T_bath)` and retains its
/// original heat and entropy accounts. The derivative differentiates Q_out,
/// not G, with the bath temperature and geometry held fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitBoundaryState {
    pub conductance_w_k: f64,
    pub power_enthalpy_derivative_w_kg_j: f64,
}

impl ImplicitBoundaryState {
    pub(super) fn admit(self) -> Result<(), DepositionError> {
        if !self.conductance_w_k.is_finite()
            || self.conductance_w_k < 0.0
            || !self.power_enthalpy_derivative_w_kg_j.is_finite()
        {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

/// Immutable complete material-to-bath flux on the original caloric object.
///
/// A source owner may integrate k(T) from the cell trace to its physical face
/// and close the film/contact jump to the bath. Only material temperatures
/// need material-source support; the bath is not another material sample.
/// Each trial must enforce its actual source domain and return a consistent
/// secant and analytic tangent without changing physical history. This law
/// adds no face storage, extra boundary heat or separate time integrator.
pub trait ImplicitBoundaryLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw;
    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError>;
}

impl<'a> ImplicitConductionChain<'a> {
    /// Replace selected reservoir and ambient bulk-series approximations with
    /// complete original-law fluxes. All unselected channels keep their exact
    /// previous arithmetic. Each selected face must already have a geometric
    /// conduction binding and use the timeline's original caloric object.
    ///
    /// One cell may have both channels because their original covered/exposed
    /// areas are distinct. Duplicate ownership within a channel, repeated
    /// attachment and absent faces refuse. The initial source query is checked
    /// before returning; no thermal state or authored boundary is changed.
    pub fn with_boundary_laws(
        mut self,
        thermal: &DepositionSimulation<'_>,
        reservoirs: &[(CellId, &'a dyn ImplicitBoundaryLaw)],
        ambient: &[(CellId, &'a dyn ImplicitBoundaryLaw)],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !self.reservoir_laws.is_empty()
            || !self.ambient_laws.is_empty()
            || (reservoirs.is_empty() && ambient.is_empty())
            || reservoirs.len() > self.conduction.reservoirs.len()
            || ambient.len() > self.conduction.ambient_faces.len()
        {
            return Err(ImplicitConductionError::Invalid(
                "one bounded nonempty implicit boundary attachment is required",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let mut reservoir_laws = BTreeMap::new();
        let mut ambient_laws = BTreeMap::new();
        for (bindings, selected, reservoir) in [
            (&mut reservoir_laws, reservoirs, true),
            (&mut ambient_laws, ambient, false),
        ] {
            for &(id, law) in selected {
                checkpoint()?;
                let present = if reservoir {
                    self.conduction.reservoirs.contains_key(&id)
                } else {
                    self.conduction.ambient_faces.contains_key(&id)
                };
                if !present || bindings.insert(id, law).is_some() {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit boundary flux needs one owner on an original conducting face",
                    ));
                }
                let cell = &thermal.cells[&id];
                if !thermal.uses_material_law(cell.material_id, law.caloric_law()) {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit boundary flux requires the original caloric object",
                    ));
                }
            }
        }
        self.reservoir_laws = reservoir_laws;
        self.ambient_laws = ambient_laws;
        self.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }
}
