//! Complete conductive surface fluxes on the original implicit heat graph.

use super::{ImplicitConductionChain, ImplicitConductionError};
use crate::enthalpy::backward_euler::PrescribedBathHeat;
use crate::enthalpy::{
    CellId, DepositionError, DepositionSimulation, ThermalLaw, finite, positive,
};
use std::collections::BTreeMap;

/// Continuous prescribed temperatures on the original boundary channels.
#[path = "boundary/history.rs"]
mod history;
pub use history::{BathTemperatureKnot, BathTemperatureSchedule};
pub(super) use history::{BoundBathHistory, temperature_at};

/// Gray radiation and convection through the original finite half-cell.
#[path = "boundary/radiation.rs"]
mod radiation;
pub use radiation::{GrayAmbientBoundary, GrayAmbientControls, GrayAmbientResponse};

/// Separate convective and radiative baths on one shared exterior surface.
#[path = "boundary/two_bath.rs"]
mod two_bath;
pub use two_bath::{GrayTwoBathBoundary, GrayTwoBathResponse};

/// Original enthalpy-dependent gray emissivity and convective surface response.
#[path = "boundary/surface.rs"]
mod surface;
pub(super) use surface::BoundGraySurface;
pub use surface::{ImplicitGraySurfaceLaw, ImplicitGraySurfaceState};

/// Original energy and actual exterior-temperature gray source response.
#[path = "boundary/exterior.rs"]
mod exterior;
pub use exterior::{
    ImplicitGrayExteriorFluxState, ImplicitGrayExteriorLaw, ImplicitGrayExteriorState,
};

/// Original centre-to-face geometry and the external film/contact resistance.
/// Area is the actual reservoir coverage or exposed ambient area. Resistance
/// is area-specific; a material owner converts it to total K/W exactly once.
/// Explicitly bound current mechanical measures replace reference area and
/// distance, held fixed during a solve; they never replace reference inventory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitBoundaryGeometry {
    pub area_m2: f64,
    pub half_distance_m: f64,
    pub surface_resistance_m2_k_w: f64,
}

/// Passive outward heat flux and its complete original-enthalpy derivative.
/// The native owner applies `Q_out = G * (T_cell - T_bath)` and retains its
/// original heat and entropy accounts. The derivative differentiates Q_out,
/// not G, with the bath temperature held fixed. The ordinary query also holds
/// geometry fixed; the geometry-aware query includes its half-distance chain rule.
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

    /// Opt in to current mechanical area/distance in the supplied geometry.
    /// Both flux and tangent must use these actual measures, subject to the
    /// law's geometric domain. The default protects reference-pinned adapters.
    fn supports_current_geometry(&self) -> bool {
        false
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError>;

    /// Original flux and COMPLETE dQ_out/dh for a density-driven half-distance.
    /// The supplied derivative is (dd/dh)/d [kg/J]; geometry contains the
    /// current positive distance. Bath temperature, area and coating/film
    /// resistance are held fixed. Resolve the same original face balance once,
    /// including its geometry chain rule without adding heat or pressure work.
    /// None leaves moving half-distance geometry explicitly unsupported.
    fn state_at_specific_enthalpy_with_geometry_tangent(
        &self,
        _specific_h_j_kg: f64,
        _bath_temperature_k: f64,
        _geometry: ImplicitBoundaryGeometry,
        _half_distance_relative_enthalpy_derivative_kg_j: f64,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitBoundaryState>, DepositionError> {
        Ok(None)
    }
}

/// Up to two signed bath exchanges through one original physical boundary face.
///
/// Powers are positive outward. Their sum is the heat removed from the cell;
/// their separate temperatures determine the reservoir entropy exchange. They
/// can oppose each other and remain nonzero when the cell receives no net heat.
/// The derivative differentiates their SUM with respect to the original h,
/// with both prescribed temperatures fixed. The ordinary query holds geometry
/// fixed; its geometry-aware counterpart includes the shared half-distance
/// response. This is not two independently conducted faces or an extra heat source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitBoundaryExchangeState {
    pub baths: [PrescribedBathHeat; 2],
    pub power_enthalpy_derivative_w_kg_j: f64,
}

impl ImplicitBoundaryExchangeState {
    pub(super) fn validate(self) -> Result<Self, DepositionError> {
        if !self.power_enthalpy_derivative_w_kg_j.is_finite()
            || self
                .baths
                .iter()
                .any(|bath| !bath.outward_power_w.is_finite() || !positive(bath.temperature_k))
        {
            return Err(DepositionError::MaterialDomain);
        }
        self.outward_power_w()?;
        Ok(self)
    }

    pub(super) fn outward_power_w(&self) -> Result<f64, DepositionError> {
        finite(self.baths[0].outward_power_w + self.baths[1].outward_power_w)
    }
}

/// Complete exchange with two prescribed baths on the original caloric object.
///
/// The selected channel supplies its actual stage bath temperature. An
/// implementation must explicitly bind the second bath and return both actual
/// temperatures. Each query is immutable and keeps one shared half-cell and
/// surface balance. The native heat and entropy owners consume the same two
/// powers; no effective single reservoir or separate source is inferred.
pub trait ImplicitBoundaryExchangeLaw {
    fn caloric_law(&self) -> &dyn ThermalLaw;

    /// Opt in to current mechanical area/distance for both bath powers and
    /// their total tangent. Source geometry admission remains authoritative;
    /// the default prevents a fixed-metric law from using stale dimensions.
    fn supports_current_geometry(&self) -> bool {
        false
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryExchangeState, DepositionError>;

    /// The same exchange at the physical endpoint or SDIRK stage time.
    /// Prescribed baths owned by this law may vary continuously with time;
    /// the returned powers and actual temperatures share that same instant.
    /// Autonomous laws retain their existing query and arithmetic.
    fn state_at_specific_enthalpy_and_time(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryExchangeState, DepositionError> {
        checkpoint()?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(DepositionError::InvalidInput(
                "invalid boundary exchange time",
            ));
        }
        self.state_at_specific_enthalpy(specific_h_j_kg, bath_temperature_k, geometry, checkpoint)
    }

    /// The same two actual bath powers and their TOTAL d(Q_0+Q_1)/dh when
    /// the shared half-distance evolves. The supplied derivative is
    /// (dd/dh)/d [kg/J], and geometry contains its current value. Both bath
    /// temperatures, area and coating resistance remain fixed; their separate
    /// entropy accounts and original shared surface solve are retained.
    /// None leaves this density/geometry coupling explicitly unsupported.
    fn state_at_specific_enthalpy_with_geometry_tangent(
        &self,
        _specific_h_j_kg: f64,
        _bath_temperature_k: f64,
        _geometry: ImplicitBoundaryGeometry,
        _half_distance_relative_enthalpy_derivative_kg_j: f64,
        _checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitBoundaryExchangeState>, DepositionError> {
        Ok(None)
    }

    /// Physical-time counterpart retaining the complete moving-distance tangent.
    /// Time and both prescribed bath temperatures remain fixed in dQ/dh.
    fn state_at_specific_enthalpy_with_geometry_tangent_and_time(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        half_distance_relative_enthalpy_derivative_kg_j: f64,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitBoundaryExchangeState>, DepositionError> {
        checkpoint()?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(DepositionError::InvalidInput(
                "invalid boundary exchange time",
            ));
        }
        self.state_at_specific_enthalpy_with_geometry_tangent(
            specific_h_j_kg,
            bath_temperature_k,
            geometry,
            half_distance_relative_enthalpy_derivative_kg_j,
            checkpoint,
        )
    }

    /// Next finite, strictly future slope change of a continuous owned bath.
    /// Fixed and adaptive native intervals stop here. Discontinuous forcing
    /// requires an explicit event law; a stage query does not infer jump sides.
    fn next_transition_time_s(&self, _after_s: f64) -> Option<f64> {
        None
    }
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
                if self.conduction.current_geometry.is_some() && !law.supports_current_geometry() {
                    return Err(ImplicitConductionError::Invalid(
                        "current mechanical boundary measures require an explicitly compatible flux law",
                    ));
                }
                let present = if reservoir {
                    self.conduction.reservoirs.contains_key(&id)
                } else {
                    self.conduction.ambient_faces.contains_key(&id)
                };
                let exchanged = if reservoir {
                    self.reservoir_exchanges.contains_key(&id)
                } else {
                    self.ambient_exchanges.contains_key(&id)
                };
                if !present || exchanged || bindings.insert(id, law).is_some() {
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

    /// Replace selected original channels with one composed two-bath exchange.
    ///
    /// The original covered/exposed area and half-cell geometry remain owned
    /// by the conduction binding. A selected channel cannot also have a legacy
    /// boundary-law owner. Separate channels on the same cell remain distinct;
    /// each must use the timeline's original caloric object. Admission checks
    /// the initial complete endpoint before returning any prepared solver.
    pub fn with_boundary_exchanges(
        mut self,
        thermal: &DepositionSimulation<'_>,
        reservoirs: &[(CellId, &'a dyn ImplicitBoundaryExchangeLaw)],
        ambient: &[(CellId, &'a dyn ImplicitBoundaryExchangeLaw)],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, ImplicitConductionError> {
        checkpoint()?;
        if !self.reservoir_exchanges.is_empty()
            || !self.ambient_exchanges.is_empty()
            || (reservoirs.is_empty() && ambient.is_empty())
            || reservoirs.len() > self.conduction.reservoirs.len()
            || ambient.len() > self.conduction.ambient_faces.len()
        {
            return Err(ImplicitConductionError::Invalid(
                "one bounded nonempty implicit boundary-exchange attachment is required",
            ));
        }
        self.admit(thermal, &mut checkpoint)?;
        let mut reservoir_exchanges = BTreeMap::new();
        let mut ambient_exchanges = BTreeMap::new();
        for (bindings, selected, reservoir) in [
            (&mut reservoir_exchanges, reservoirs, true),
            (&mut ambient_exchanges, ambient, false),
        ] {
            for &(id, law) in selected {
                checkpoint()?;
                if self.conduction.current_geometry.is_some() && !law.supports_current_geometry() {
                    return Err(ImplicitConductionError::Invalid(
                        "current mechanical boundary measures require an explicitly compatible exchange law",
                    ));
                }
                let (present, occupied) = if reservoir {
                    (
                        self.conduction.reservoirs.contains_key(&id),
                        self.reservoir_laws.contains_key(&id),
                    )
                } else {
                    (
                        self.conduction.ambient_faces.contains_key(&id),
                        self.ambient_laws.contains_key(&id),
                    )
                };
                if !present || occupied || bindings.insert(id, law).is_some() {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit boundary exchange needs sole ownership of an original conducting face",
                    ));
                }
                let cell = &thermal.cells[&id];
                if !thermal.uses_material_law(cell.material_id, law.caloric_law()) {
                    return Err(ImplicitConductionError::Invalid(
                        "implicit boundary exchange requires the original caloric object",
                    ));
                }
            }
        }
        self.reservoir_exchanges = reservoir_exchanges;
        self.ambient_exchanges = ambient_exchanges;
        self.endpoint(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(self)
    }
}
