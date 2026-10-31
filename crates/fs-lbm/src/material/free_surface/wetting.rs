//! Source-pinned hysteresis drives the existing planar capillary boundary.
//!
//! This explicit geometric closure estimates contact-line speed from successive
//! accepted wall-adjacent VOF fields. It does not identify wall-fluid velocity
//! with a no-slip contact-line velocity. The estimate is first order in space
//! at the wall and backward in time; its angle drives the NEXT native step.
//! In the deadband, the measured angle is projected onto the sourced interval.
//! That is a bounded hysteresis closure, not a proof of exact pinning, calibrated
//! dynamic wetting, or closed surface/mechanical/thermal energy conservation.

use std::collections::BTreeMap;
use std::sync::Arc;

use super::*;
use crate::d3q19::freesurface3::wetting::WallFill3;
use crate::freesurface::{
    ADVANCING_CONTACT_ANGLE_PROPERTY, ContactLineRegime2, DYNAMIC_WETTING_SPEED_AXIS,
    DYNAMIC_WETTING_SPEED_DIMS, RECEDING_CONTACT_ANGLE_PROPERTY,
};
use fs_matdb::{InterfaceSystemCard, MatDbError};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_qty::{Dims, QuantitySpec};

mod cox_voinov;
pub use cox_voinov::{CoxVoinovControls, CoxVoinovLaw, CoxVoinovSample};

/// Local liquid properties supplied by the existing enthalpy source owner.
#[derive(Clone, Copy)]
pub(super) struct WettingThermalPoint {
    pub temperature_k: f64,
    pub dynamic_viscosity_pa_s: f64,
}

#[derive(Clone)]
struct HydrodynamicWetting {
    law: CoxVoinovLaw,
    reference_viscosity_pa_s: f64,
    surface_tension_n_m: f64,
    outer_gas: Option<Arc<ResolvedThermochemicalGasState>>,
}

/// Explicit model/numerical inputs; the deadband has no invented source receipt.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceWettingControls {
    /// Symmetric physical speed deadband [m/s], finite and nonnegative.
    pub pinning_speed_m_s: f64,
    /// Maximum active wall-adjacent cells, including gas history.
    pub maximum_wall_cells: usize,
}

/// One accepted contact estimate and the two ORIGINAL source endpoint claims.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceWettingSample {
    /// Positive is advancing liquid along its outward in-wall normal [m/s].
    pub signed_speed_m_s: f64,
    /// Direction/deadband selection; `Pinned` does not certify zero line motion.
    pub regime: ContactLineRegime2,
    /// Geometric angle through liquid from the accepted first-cell fill [rad].
    pub geometric_angle_radians: f64,
    /// Angle used by the next native curvature reconstruction [rad].
    pub applied_angle_radians: f64,
    /// Optional explicit hydrodynamic prediction, separate from static source claims.
    pub hydrodynamic: Option<CoxVoinovSample>,
    /// Explicit outer-phase transport source, including its own temperature,
    /// pressure and composition. Shared by this wall field, never substituted
    /// for the liquid/interface temperature. None for a declared fixed ratio.
    pub outer_gas: Option<Arc<ResolvedThermochemicalGasState>>,
    /// Both pinned endpoint claims, ordered surfaces, source point and receipts.
    pub material: ResolvedInterfaceStatePoint,
}

/// A source, geometry, numerical or cancellation refusal; no state is published.
#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceWettingError {
    Material(MaterialStatePointError),
    Query(MatDbError),
    Flow(MaterialLbmError),
    Invalid(&'static str),
}

impl core::fmt::Display for FreeSurfaceWettingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => write!(f, "free-surface wetting material: {e}"),
            Self::Query(e) => write!(f, "free-surface wetting query: {e}"),
            Self::Flow(e) => write!(f, "free-surface wetting: {e}"),
            Self::Invalid(message) => write!(f, "free-surface wetting: {message}"),
        }
    }
}
impl std::error::Error for FreeSurfaceWettingError {}
impl From<MaterialStatePointError> for FreeSurfaceWettingError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<MatDbError> for FreeSurfaceWettingError {
    fn from(e: MatDbError) -> Self {
        Self::Query(e)
    }
}
impl From<crate::d3q19::freesurface3::FreeSurfaceError3> for FreeSurfaceWettingError {
    fn from(e: crate::d3q19::freesurface3::FreeSurfaceError3) -> Self {
        Self::Flow(e.into())
    }
}

type Error = FreeSurfaceWettingError;

/// One frozen isothermal liquid and its evolving planar contact-line history.
///
/// Bulk density, viscosity and sigma remain the original liquid's coefficients.
/// Surface A is the assigned stationary substrate, surface B must be this exact
/// liquid phase/process/revision. The original full liquid query is retained;
/// only the explicit physical contact-line speed coordinate evolves. Named gas
/// and surface/history context belongs to the interface card, not a guessed EOS.
pub struct MaterialFreeSurfaceWetting3<'s, 'm> {
    flow: MaterialFreeSurface3<'s>,
    binding: WettingBinding<'m>,
}

// Frozen and thermal flow share one source/geometry/history owner. Thermal
// callers supply temperatures; the geometric speed and angle arithmetic stays
// identical to the original isothermal closure.
#[derive(Clone)]
pub(super) struct WettingBinding<'m> {
    card: &'m InterfaceSystemCard,
    source: ResolvedInterfaceStatePoint,
    wall: PlanarWall3,
    controls: FreeSurfaceWettingControls,
    geometry: BTreeMap<[usize; 3], WallFill3>,
    samples: BTreeMap<[usize; 3], FreeSurfaceWettingSample>,
    hydrodynamic: Option<HydrodynamicWetting>,
}

impl<'s> MaterialFreeSurface3<'s> {
    /// Bind a rate-dependent advancing/receding interface to one planar wall.
    /// Initial history is stationary: both source claims must support speed zero.
    /// The caller assigns substrate and named environment to the scene. Existing
    /// static wetting and the legacy Wetting ghost cannot be combined with this
    /// law. Other domain walls and internal solid tiles remain neutral.
    pub fn with_dynamic_planar_wetting<'m>(
        mut self,
        wall: PlanarWall3,
        card: &'m InterfaceSystemCard,
        selection: MaterialPropertySelection,
        controls: FreeSurfaceWettingControls,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceWetting3<'s, 'm>, Error> {
        let binding = WettingBinding::try_new(
            &mut self, wall, card, selection, controls, gate,
        )?;
        binding.install(&mut self.lattice, gate)?;
        checkpoint(gate)?;
        Ok(MaterialFreeSurfaceWetting3 { flow: self, binding })
    }
}

impl<'m> WettingBinding<'m> {
    pub(super) fn try_new(
        flow: &mut MaterialFreeSurface3<'_>,
        wall: PlanarWall3,
        card: &'m InterfaceSystemCard,
        selection: MaterialPropertySelection,
        controls: FreeSurfaceWettingControls,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        checkpoint(gate)?;
        if !controls.pinning_speed_m_s.is_finite()
            || controls.pinning_speed_m_s < 0.0
            || controls.maximum_wall_cells == 0
            || controls.maximum_wall_cells > 262_144
        {
            return Err(Error::Invalid(
                "invalid speed deadband or bounded wall-cell budget",
            ));
        }
        if flow.scaling.planar_wetting.is_some() || !flow.lattice.admits_dynamic_planar_wetting() {
            return Err(Error::Invalid(
                "dynamic planar wetting requires otherwise neutral walls",
            ));
        }
        let liquid = flow.scaling.material.material();
        if &card.surface_b().material != liquid.material()
            || liquid
                .query_point()
                .iter()
                .any(|(axis, _)| axis == DYNAMIC_WETTING_SPEED_AXIS)
        {
            return Err(Error::Invalid(
                "wetting requires this liquid identity and a separate evolving speed axis",
            ));
        }
        let point = liquid.query().with_quantity(
            DYNAMIC_WETTING_SPEED_AXIS,
            QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS),
            0.0,
        )?;
        let requirements = [
            RECEDING_CONTACT_ANGLE_PROPERTY,
            ADVANCING_CONTACT_ANGLE_PROPERTY,
        ]
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                Dims::NONE,
                ScalarAdmissibility::OpenInterval {
                    lower: 0.0,
                    upper: core::f64::consts::PI,
                },
            )
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
        let source = resolve_interface_state_point(card, &point, &requirements, selection)?;
        require_owned_axes(&source)?;
        bracket(&source)?;
        let geometry = flow
            .lattice
            .planar_wall_fill(wall, controls.maximum_wall_cells, gate)?;
        let mut candidate = Self {
            card,
            source,
            wall,
            controls,
            geometry,
            samples: BTreeMap::new(),
            hydrodynamic: None,
        };
        let samples = candidate.resolve_samples(
            &candidate.geometry, None, &flow.scaling, None, gate,
        )?;
        checkpoint(gate)?;
        candidate.samples = samples;
        Ok(candidate)
    }
}

impl MaterialFreeSurfaceWetting3<'_, '_> {
    /// Add the full-angle Cox law to STATIC advancing/receding endpoints.
    /// Any finite outer/liquid viscosity ratio is an explicit fixed model input;
    /// this boundary does not evolve a second-fluid momentum field.
    /// Source curves with a contact-speed axis are refused to avoid counting
    /// measured rate dependence twice. Source receipts remain at speed zero;
    /// actual speed, viscosity and model-predicted angle are separate samples.
    /// The apparent angle feeds the existing native capillary reconstruction.
    pub fn with_cox_voinov(mut self, law: CoxVoinovLaw, gate: &CancelGate) -> Result<Self, Error> {
        let binding = self
            .binding
            .with_cox_voinov(law, &self.flow.scaling, gate)?;
        binding.install(&mut self.flow.lattice, gate)?;
        checkpoint(gate)?;
        self.binding = binding;
        Ok(self)
    }

    /// Explicit scales and limits of the optional hydrodynamic closure.
    #[must_use]
    pub fn cox_voinov_law(&self) -> Option<&CoxVoinovLaw> {
        self.binding.cox_voinov_law()
    }

    /// Read-only actual flow, native populations, mass ledger and physical units.
    #[must_use]
    pub const fn flow(&self) -> &MaterialFreeSurface3<'_> {
        &self.flow
    }

    /// Original two-claim binding at stationary initialization.
    #[must_use]
    pub const fn source(&self) -> &ResolvedInterfaceStatePoint {
        &self.binding.source
    }

    /// Accepted endpoint queries and next-step angles by wall-cell coordinate.
    #[must_use]
    pub const fn samples(&self) -> &BTreeMap<[usize; 3], FreeSurfaceWettingSample> {
        &self.binding.samples
    }

    /// Advance native momentum/mass and atomically refresh the contact history.
    /// Endpoint source-domain, degeneracy, budget and cancellation refusals keep
    /// populations, fill, carry, conversions, angle field and receipts unchanged.
    pub fn step_with_gate(&mut self, gate: &CancelGate) -> Result<(), Error> {
        checkpoint(gate)?;
        let mut lattice = self.flow.lattice.clone();
        lattice.step_with_gate(gate)?;
        let binding = self.binding.at_endpoint(&mut lattice, &self.flow.scaling, None, gate)?;
        checkpoint(gate)?;
        self.flow.lattice = lattice;
        self.binding = binding;
        Ok(())
    }
}

impl WettingBinding<'_> {
    pub(super) fn cox_voinov_law(&self) -> Option<&CoxVoinovLaw> {
        self.hydrodynamic.as_ref().map(|binding| &binding.law)
    }

    pub(super) fn with_cox_voinov(
        mut self,
        law: CoxVoinovLaw,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        checkpoint(gate)?;
        if self.hydrodynamic.is_some() {
            return Err(Error::Invalid("contact line already owns a Cox model"));
        }
        if scaling.material.rheology_axes().is_some() {
            return Err(Error::Invalid(
                "Cox wetting requires Newtonian liquid viscosity",
            ));
        }
        for property in self.source.properties() {
            if property
                .source_axes()
                .iter()
                .any(|(axis, _)| axis == DYNAMIC_WETTING_SPEED_AXIS)
            {
                return Err(Error::Invalid(
                    "Cox wetting cannot add rate dependence to a measured dynamic-angle curve",
                ));
            }
            law.admit_angle(property.value_si())?;
        }
        let surface_tension_n_m = scaling
            .material
            .surface_tension_n_m()
            .ok_or(Error::Invalid(
                "Cox wetting requires sourced positive surface tension",
            ))?;
        self.hydrodynamic = Some(HydrodynamicWetting {
            law,
            reference_viscosity_pa_s: scaling.material.dynamic_viscosity_pa_s(),
            surface_tension_n_m,
            outer_gas: None,
        });
        let mut samples = BTreeMap::new();
        for (&position, current) in &self.geometry {
            checkpoint(gate)?;
            let speed = self.samples[&position].signed_speed_m_s;
            samples.insert(position, self.sample(current, speed, None, gate)?);
        }
        self.samples = samples;
        Ok(self)
    }

    pub(super) fn with_cox_voinov_gas(
        self,
        controls: CoxVoinovControls,
        gas: &ResolvedThermochemicalGasState,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let ratio = gas.state().dynamic_viscosity / scaling.material.dynamic_viscosity_pa_s();
        if !ratio.is_finite() || ratio <= 0.0 {
            return Err(Error::Invalid("Cox outer gas requires positive representable source viscosities"));
        }
        let law = CoxVoinovLaw::try_new_with_viscosity_ratio(controls, ratio, gate)?;
        let mut candidate = self.with_cox_voinov(law, scaling, gate)?;
        let source = Arc::new(gas.clone());
        candidate.hydrodynamic.as_mut().expect("admitted Cox law").outer_gas = Some(source.clone());
        for sample in candidate.samples.values_mut() {
            checkpoint(gate)?;
            sample.outer_gas = Some(source.clone());
        }
        Ok(candidate)
    }

    pub(super) fn outer_gas(&self) -> Option<&ResolvedThermochemicalGasState> {
        self.hydrodynamic.as_ref().and_then(|binding| binding.outer_gas.as_deref())
    }

    /// The finite-atmosphere timeline admits original transport sources before
    /// passing this prepared snapshot. Existing fixed-ratio bindings stay fixed.
    pub(super) fn at_prepared_outer_gas(&self, gas: &ResolvedThermochemicalGasState) -> Self {
        let mut candidate = self.clone();
        if let Some(binding) = &mut candidate.hydrodynamic
            && binding.outer_gas.is_some()
            && binding.outer_gas.as_deref() != Some(gas)
        {
            binding.outer_gas = Some(Arc::new(gas.clone()));
        }
        candidate
    }

    pub(super) fn source(&self) -> &ResolvedInterfaceStatePoint {
        &self.source
    }

    pub(super) fn samples(&self) -> &BTreeMap<[usize; 3], FreeSurfaceWettingSample> {
        &self.samples
    }

    pub(super) fn install(
        &self,
        lattice: &mut FreeSurface3,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        lattice.set_planar_contact_field(
            self.wall,
            self.samples.iter().map(|(&p, s)| (p, s.applied_angle_radians)).collect(),
            gate,
        )?;
        Ok(())
    }

    pub(super) fn at_endpoint(
        &self,
        lattice: &mut FreeSurface3,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        temperatures: Option<(&str, &BTreeMap<[usize; 3], WettingThermalPoint>)>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let geometry = lattice.planar_wall_fill(self.wall, self.controls.maximum_wall_cells, gate)?;
        let samples = self.resolve_samples(
            &geometry, Some(&self.geometry), scaling, temperatures, gate,
        )?;
        let candidate = Self { geometry, samples, ..self.clone() };
        candidate.install(lattice, gate)?;
        Ok(candidate)
    }

    // Temperature changes before motion hold the last accepted geometric
    // speed, rather than resetting the hysteresis branch to stationary.
    pub(super) fn at_temperatures(
        &self,
        axis: &str,
        temperatures: &BTreeMap<[usize; 3], WettingThermalPoint>,
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        let mut samples = BTreeMap::new();
        for (&position, current) in &self.geometry {
            checkpoint(gate)?;
            let speed = self.samples[&position].signed_speed_m_s;
            let temperature = temperatures.get(&position).map(|&t| (axis, t));
            samples.insert(position, self.sample(current, speed, temperature, gate)?);
        }
        Ok(Self { samples, ..self.clone() })
    }

    pub(super) fn admit_temperature(
        &self,
        position: [usize; 3],
        axis: &str,
        temperature_k: f64,
        dynamic_viscosity_pa_s: f64,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        if let Some(current) = self.geometry.get(&position) {
            self.sample(
                current,
                self.samples[&position].signed_speed_m_s,
                Some((axis, WettingThermalPoint { temperature_k, dynamic_viscosity_pa_s })),
                gate,
            )?;
        }
        Ok(())
    }

    fn resolve_samples(
        &self,
        geometry: &BTreeMap<[usize; 3], WallFill3>,
        previous: Option<&BTreeMap<[usize; 3], WallFill3>>,
        scaling: &MaterialFreeSurfaceScaling<'_>,
        temperatures: Option<(&str, &BTreeMap<[usize; 3], WettingThermalPoint>)>,
        gate: &CancelGate,
    ) -> Result<BTreeMap<[usize; 3], FreeSurfaceWettingSample>, Error> {
        let mut result = BTreeMap::new();
        for (&position, current) in geometry {
            checkpoint(gate)?;
            let old_fill = previous.map_or(current.fill, |map| {
                map.get(&position).map_or(0.0, |g| g.fill)
            });
            let speed = contact_speed(
                current.fill - old_fill,
                current.tangent_gradient,
                scaling.cell_size_m(),
                scaling.timestep_s(),
            )?;
            let temperature = temperatures
                .and_then(|(axis, field)| field.get(&position).map(|&t| (axis, t)));
            result.insert(position, self.sample(current, speed, temperature, gate)?);
        }
        Ok(result)
    }

    fn sample(
        &self,
        current: &WallFill3,
        speed: f64,
        temperature: Option<(&str, WettingThermalPoint)>,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceWettingSample, Error> {
        // A hydrodynamic prediction uses zero-speed endpoint measurements.
        // Its rate dependence must never be laundered into a material receipt.
        let query_speed = if self.hydrodynamic.is_some() {
            0.0
        } else {
            speed.abs()
        };
        let mut point = self
            .source
            .query()
            .with(DYNAMIC_WETTING_SPEED_AXIS, query_speed)?;
        if let Some((axis, value)) = temperature {
            point = point.with(axis, value.temperature_k)?;
        }
        let material = self.source.at_query_point(self.card, &point)?;
        let (receding, advancing) = bracket(&material)?;
        let (regime, mut applied_angle_radians) = select_angle(
            speed, self.controls.pinning_speed_m_s, current.angle_radians, receding, advancing,
        );
        let hydrodynamic = if let Some(binding) = &self.hydrodynamic {
            binding.law.admit_angle(receding)?;
            binding.law.admit_angle(advancing)?;
            let viscosity = temperature.map_or(binding.reference_viscosity_pa_s, |(_, state)| {
                state.dynamic_viscosity_pa_s
            });
            if !viscosity.is_finite() || viscosity <= 0.0 {
                return Err(Error::Invalid(
                    "Cox wetting needs finite positive source viscosity",
                ));
            }
            let capillary = (viscosity * speed) / binding.surface_tension_n_m;
            let ratio = binding.outer_gas.as_ref().map_or(
                binding.law.outer_to_liquid_viscosity_ratio(),
                |gas| gas.state().dynamic_viscosity / viscosity,
            );
            if binding.outer_gas.is_some() && ratio <= 0.0 {
                return Err(Error::Invalid("Cox source viscosity ratio is not representable"));
            }
            binding.law.admit_capillary_number_at_ratio(capillary, ratio)?;
            if regime != ContactLineRegime2::Pinned {
                applied_angle_radians = binding
                    .law
                    .apparent_angle_at_ratio(applied_angle_radians, capillary, ratio, gate)?;
            }
            Some(CoxVoinovSample {
                signed_capillary_number: capillary,
                dynamic_viscosity_pa_s: viscosity,
                surface_tension_n_m: binding.surface_tension_n_m,
                outer_to_liquid_viscosity_ratio: ratio,
            })
        } else {
            None
        };
        Ok(FreeSurfaceWettingSample {
            signed_speed_m_s: speed,
            regime,
            geometric_angle_radians: current.angle_radians,
            applied_angle_radians,
            hydrodynamic,
            outer_gas: self.hydrodynamic.as_ref().and_then(|binding| binding.outer_gas.clone()),
            material,
        })
    }
}

fn checkpoint(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(crate::d3q19::freesurface3::FreeSurfaceError3::Cancelled.into())
    } else {
        Ok(())
    }
}

// At a fixed point of a translating fill profile, d(phi)/dt = V |grad_t(phi)|
// for outward liquid normal -grad_t(phi)/|grad_t(phi)|. Gradient here is per
// LATTICE CELL, so dx/dt converts its displacement to physical metres/second.
fn contact_speed(
    delta_fill: f64,
    gradient_per_cell: f64,
    dx_m: f64,
    dt_s: f64,
) -> Result<f64, Error> {
    if gradient_per_cell == 0.0 {
        return if delta_fill == 0.0 {
            Ok(0.0)
        } else {
            Err(Error::Invalid(
                "moving contact geometry has no resolved tangential gradient",
            ))
        };
    }
    let speed = (delta_fill / gradient_per_cell) * (dx_m / dt_s);
    if !speed.is_finite() || delta_fill != 0.0 && speed == 0.0 {
        return Err(Error::Invalid("contact-line speed is unrepresentable"));
    }
    Ok(speed)
}

fn select_angle(
    speed: f64,
    deadband: f64,
    geometric: f64,
    receding: f64,
    advancing: f64,
) -> (ContactLineRegime2, f64) {
    if speed > deadband {
        (ContactLineRegime2::Advancing, advancing)
    } else if speed < -deadband {
        (ContactLineRegime2::Receding, receding)
    } else {
        (
            ContactLineRegime2::Pinned,
            geometric.clamp(receding, advancing),
        )
    }
}

fn bracket(material: &ResolvedInterfaceStatePoint) -> Result<(f64, f64), Error> {
    let angle = |name| {
        material
            .property(name)
            .map(|p| p.value_si())
            .ok_or(Error::Invalid(
                "both original contact-angle endpoint claims are required",
            ))
    };
    let receding = angle(RECEDING_CONTACT_ANGLE_PROPERTY)?;
    let advancing = angle(ADVANCING_CONTACT_ANGLE_PROPERTY)?;
    if receding > advancing {
        return Err(Error::Invalid("receding angle exceeds advancing angle"));
    }
    Ok((receding, advancing))
}

// The contact history owns this speed only. Temperature and fixed environment
// remain the frozen liquid context; frequency, strain, stress, or another rate
// cannot silently become constants merely because an initial sample exists.
fn require_owned_axes(source: &ResolvedInterfaceStatePoint) -> Result<(), Error> {
    for property in source.properties() {
        if !property.has_only_owned_and_fixed_external_conditions(
            DYNAMIC_WETTING_SPEED_AXIS,
            QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS),
        ) {
            return Err(Error::Invalid(
                "wetting source needs its owned speed and canonical fixed external conditions",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g1_translating_fill_has_physical_speed_and_refinement_invariance() {
        // Independent physical plane phi(x,t)=1/2 - a*(x - V*t), a [1/m].
        // Its lattice derivative is a*dx and its fixed-cell change is a*V*dt.
        for speed in [-0.021, 0.0, 0.037] {
            for dx in [0.001, 0.0005] {
                for dt in [0.0001, 0.000025] {
                    let slope_per_m = 31.0;
                    let actual =
                        contact_speed(slope_per_m * speed * dt, slope_per_m * dx, dx, dt).unwrap();
                    assert!((actual - speed).abs() < 1e-16);
                }
            }
        }
        assert_eq!(contact_speed(0.0, 0.0, 0.001, 0.0001).unwrap(), 0.0);
        assert!(contact_speed(0.001, 0.0, 0.001, 0.0001).is_err());
        assert!(contact_speed(1.0, f64::MIN_POSITIVE, 1.0, 1e-4).is_err());
    }

    #[test]
    fn g1_direction_and_closed_deadband_use_the_sourced_hysteresis_bracket() {
        let receding = 0.9;
        let advancing = 2.1;
        assert_eq!(
            select_angle(0.02, 0.01, 1.5, receding, advancing),
            (ContactLineRegime2::Advancing, advancing)
        );
        assert_eq!(
            select_angle(-0.02, 0.01, 1.5, receding, advancing),
            (ContactLineRegime2::Receding, receding)
        );
        for speed in [-0.01, 0.0, 0.01] {
            assert_eq!(
                select_angle(speed, 0.01, 1.5, receding, advancing),
                (ContactLineRegime2::Pinned, 1.5)
            );
            assert_eq!(
                select_angle(speed, 0.01, 0.3, receding, advancing),
                (ContactLineRegime2::Pinned, receding)
            );
            assert_eq!(
                select_angle(speed, 0.01, 2.9, receding, advancing),
                (ContactLineRegime2::Pinned, advancing)
            );
        }
    }
}
