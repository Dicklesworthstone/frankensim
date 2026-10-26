//! Moving-liquid extensive enthalpy on the original VOF mass routes.
//!
//! Each step applies declared heat to the accepted liquid, refreshes the SAME
//! source claims at its caloric temperature, then runs native collision with
//! that cell's viscosity. Backward-Euler donor transport carries `mass*h` over
//! the native face and conversion/carry routes. The existing sparse Krylov
//! owner solves that transport; no velocity-derived second mass flux exists.
//!
//! This is the single-liquid, constant-density, stationary-wall rung. Capillary
//! energy, phase mobility, viscous heat, pressure work and conduction are not
//! supplied. In particular, h is enthalpy, not internal or total energy. The
//! incompressible lattice pressure is not a thermodynamic EOS. A negative native
//! accepted interface/carry mass refuses rather than becoming thermal storage.
//! A source-domain, numerical-budget or cancellation refusal publishes nothing.

use std::collections::{BTreeMap, BTreeSet};

use fs_exec::CancelGate;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid::surface::LiquidCapillarity;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError, SolidLiquidPhase};
use fs_material::state_point::MaterialStatePointError;
use fs_qty::Temperature;

use super::{MaterialFreeSurface3, MaterialLbmError};
use crate::d3q19::freesurface3::transport::MassLocation3;
use crate::d3q19::freesurface3::{FreeSurface3, FreeSurfaceError3};
use crate::material::plan_momentum_scaling;

mod transport;

/// Explicit bounds for the sparse companion transport.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceEnthalpyControls {
    /// Maximum wet/storage/route nodes, including the conversion carry.
    pub max_cells: usize,
    /// Maximum nonzero native pair and conversion transfers per step.
    pub max_mass_routes: usize,
    /// Shared sparse Krylov iteration bound.
    pub max_iterations: usize,
    /// True scaled transport residual and extensive balance allowance.
    pub relative_tolerance: f64,
}

impl Default for FreeSurfaceEnthalpyControls {
    fn default() -> Self {
        Self {
            max_cells: 16_384,
            max_mass_routes: 400_000,
            max_iterations: 2_000,
            relative_tolerance: 1e-11,
        }
    }
}

/// Heat supplied during one step to an already occupied cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceCellHeat {
    /// Global Cartesian cell coordinates at the beginning of the step.
    pub coordinates: [usize; 3],
    /// Signed external heat [J]. This is not an inferred mechanical loss.
    pub energy_j: f64,
}

/// Actual conservation checks from one accepted step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceEnthalpyReport {
    /// Original VOF mass before and after the step [kg], including carry.
    pub mass_before_kg: f64,
    /// Original VOF mass after the step [kg], including carry.
    pub mass_after_kg: f64,
    /// Extensive liquid enthalpy before the step [J], including carry.
    pub enthalpy_before_j: f64,
    /// Extensive liquid enthalpy after the step [J], including carry.
    pub enthalpy_after_j: f64,
    /// Exactly the sum of the caller's admitted cell heat [J].
    pub external_heat_j: f64,
    /// `H_after - H_before - external_heat` [J]. No balancing heat is added.
    pub balance_defect_j: f64,
    /// True residual returned by the existing sparse solver.
    pub relative_transport_residual: f64,
    /// Actual sparse transport iterations.
    pub transport_iterations: usize,
    /// Count of native nonzero mass routes, including conversion redistribution.
    pub mass_routes: usize,
}

/// Typed source, native-flow, numerical or physical-domain refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceEnthalpyError {
    /// The original material/free-surface scaling or flow owner refused.
    Flow(MaterialLbmError),
    /// Original source query refusal, including pinned-claim support.
    Material(MaterialStatePointError),
    /// Original caloric chart refusal; no temperature extrapolation occurs.
    Phase(PhaseStateError),
    /// The existing sparse transport solver refused.
    Solve(crate::conjugate::ChtError),
    /// Missing physical closure, invalid input or exhausted explicit budget.
    Invalid(&'static str),
}

impl core::fmt::Display for FreeSurfaceEnthalpyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Flow(e) => write!(f, "free-surface enthalpy: {e}"),
            Self::Material(e) => write!(f, "free-surface enthalpy source: {e}"),
            Self::Phase(e) => write!(f, "free-surface enthalpy chart: {e}"),
            Self::Solve(e) => write!(f, "free-surface enthalpy solve: {e}"),
            Self::Invalid(e) => write!(f, "free-surface enthalpy: {e}"),
        }
    }
}
impl std::error::Error for FreeSurfaceEnthalpyError {}
impl From<MaterialLbmError> for FreeSurfaceEnthalpyError {
    fn from(e: MaterialLbmError) -> Self {
        Self::Flow(e)
    }
}
impl From<FreeSurfaceError3> for FreeSurfaceEnthalpyError {
    fn from(e: FreeSurfaceError3) -> Self {
        Self::Flow(e.into())
    }
}
impl From<MaterialStatePointError> for FreeSurfaceEnthalpyError {
    fn from(e: MaterialStatePointError) -> Self {
        Self::Material(e)
    }
}
impl From<PhaseStateError> for FreeSurfaceEnthalpyError {
    fn from(e: PhaseStateError) -> Self {
        Self::Phase(e)
    }
}
impl From<crate::conjugate::ChtError> for FreeSurfaceEnthalpyError {
    fn from(e: crate::conjugate::ChtError) -> Self {
        Self::Solve(e)
    }
}

type Error = FreeSurfaceEnthalpyError;

/// Native moving liquid, original caloric chart and pinned source law.
///
/// `h` is retained for every wet collision lane, including zero-mass interface
/// ghosts; only `m*h` contributes storage. Unconnected zero-mass ghosts receive
/// a mass-weighted one-ring temperature extension, without moving any energy.
/// A zero-mass carry has no assigned enthalpy unless it participates in a route.
pub struct MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    flow: MaterialFreeSurface3<'s>,
    law: &'m LiquidTemperatureLaw<'m>,
    caloric: &'m EquilibriumEnthalpyPhaseCurve,
    controls: FreeSurfaceEnthalpyControls,
    h: BTreeMap<MassLocation3, f64>,
    states: BTreeMap<MassLocation3, ResolvedLiquidState>,
}

impl<'s> MaterialFreeSurface3<'s> {
    /// Bind liquid enthalpy to this actual VOF state without changing its mass
    /// or populations. The callback supplies specific enthalpy [J/kg] for every
    /// initial wet lane, including any zero-mass interface ghost.
    ///
    /// The caloric chart and liquid law must name the same original material
    /// card as the surface source. The source must explicitly neglect capillarity
    /// and use the same viscosity basis and complete initial query context.
    /// Existing nonzero conversion carry requires a prior thermal history and
    /// cannot be initialized by assigning it an arbitrary temperature here.
    pub fn with_enthalpy<'m>(
        self,
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        initial_h: impl Fn([usize; 3]) -> f64,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceEnthalpy3<'s, 'm>, Error> {
        poll(gate)?;
        if controls.max_cells == 0
            || controls.max_cells > 262_144
            || controls.max_mass_routes == 0
            || controls.max_mass_routes > 16_000_000
            || controls.max_iterations == 0
            || controls.max_iterations > 100_000
            || !controls.relative_tolerance.is_finite()
            || !(1e-14..=1e-6).contains(&controls.relative_tolerance)
        {
            return Err(Error::Invalid(
                "invalid bounded enthalpy transport controls",
            ));
        }
        let source = self.scaling.material();
        if source.capillarity() != LiquidCapillarity::Neglected
            || self.scaling.planar_wetting.is_some()
        {
            return Err(Error::Invalid(
                "thermal free surface needs explicitly neglected capillarity and no wetting law",
            ));
        }
        if source.material().card_identity() != law.reference_state().material().card_identity()
            || caloric.material_card_identity() != source.material().card_identity()
            || source.viscosity_basis() != law.viscosity_basis()
            || source.material().query_point() != law.reference_state().material().query_point()
            || source.material().axis_quantities()
                != law.reference_state().material().axis_quantities()
        {
            return Err(Error::Invalid(
                "surface, caloric and temperature laws must retain one material and query context",
            ));
        }
        let viscosity = match source.viscosity_basis() {
            fs_material::liquid_temperature::LiquidViscosityBasis::Dynamic => "dynamic_viscosity",
            fs_material::liquid_temperature::LiquidViscosityBasis::Kinematic => {
                "kinematic_viscosity"
            }
        };
        for name in ["density", viscosity] {
            if source
                .material()
                .property(name)
                .map(|p| p.answer().receipt.selected)
                != law
                    .reference_state()
                    .material()
                    .property(name)
                    .map(|p| p.answer().receipt.selected)
            {
                return Err(Error::Invalid(
                    "temperature law must retain the surface's original density and viscosity claims",
                ));
            }
        }
        let masses = self.lattice.mass_snapshot();
        admit_masses(&masses, controls.max_cells)?;
        self.ledger_mass_kg()?;
        if masses[&MassLocation3::Carry] != 0.0 {
            return Err(Error::Invalid(
                "initial conversion carry needs an existing enthalpy history",
            ));
        }
        let mut h = BTreeMap::new();
        for &location in masses.keys() {
            poll(gate)?;
            if let Some(coordinates) = FreeSurface3::transport_coordinates(location) {
                h.insert(location, initial_h(coordinates));
            }
        }
        let mut candidate = MaterialFreeSurfaceEnthalpy3 {
            flow: self,
            law,
            caloric,
            controls,
            h,
            states: BTreeMap::new(),
        };
        candidate.states = candidate.resolve_states(&candidate.h, &masses, gate)?.0;
        Ok(candidate)
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Original flow solver and physical scaling, available read-only.
    #[must_use]
    pub const fn flow(&self) -> &MaterialFreeSurface3<'_> {
        &self.flow
    }

    /// Accepted specific enthalpy [J/kg] at a native wet cell.
    #[must_use]
    pub fn specific_enthalpy_at(&self, coordinates: [usize; 3]) -> Option<f64> {
        self.location_at(coordinates)
            .and_then(|key| self.h.get(&key).copied())
    }

    /// Accepted source tuple and original receipts at that cell's caloric T.
    #[must_use]
    pub fn material_at(&self, coordinates: [usize; 3]) -> Option<&ResolvedLiquidState> {
        self.location_at(coordinates)
            .and_then(|key| self.states.get(&key))
    }

    /// Total `mass*h` [J], including any actual native conversion carry.
    pub fn total_enthalpy_j(&self) -> Result<f64, Error> {
        total(
            &self.flow.lattice.mass_snapshot(),
            &self.h,
            self.flow.scaling.mass_unit_kg(),
        )
    }

    /// One native flow step and conservative enthalpy transport with no heat.
    pub fn step(&mut self, gate: &CancelGate) -> Result<FreeSurfaceEnthalpyReport, Error> {
        self.step_with_heat(&[], gate)
    }

    /// Apply supplied cell heat, then advance original momentum and extensive
    /// enthalpy. All source, mass, energy and cancellation gates precede commit.
    ///
    /// Heat is an operator-split beginning-of-step source. The transport is
    /// first-order implicit upwind and may diffuse temperature; its numerical
    /// mixing is not reported as physical dissipation. Conversion uses the
    /// original global redistribution pool, including its numerical nonlocality.
    pub fn step_with_heat(
        &mut self,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<FreeSurfaceEnthalpyReport, Error> {
        poll(gate)?;
        if heat.len() > self.controls.max_cells {
            return Err(Error::Invalid("cell-heat budget exhausted"));
        }
        let masses = self.flow.lattice.mass_snapshot();
        let mass_unit = self.flow.scaling.mass_unit_kg();
        let before = total(&masses, &self.h, mass_unit)?;
        let mut heated = self.h.clone();
        let mut supplied = BTreeSet::new();
        let mut external = 0.0;
        for source in heat {
            poll(gate)?;
            let location = self.location_at(source.coordinates).ok_or(Error::Invalid(
                "cell heat requires an existing occupied cell",
            ))?;
            let mass = masses[&location] * mass_unit;
            if mass <= 0.0
                || !mass.is_finite()
                || !source.energy_j.is_finite()
                || !supplied.insert(location)
            {
                return Err(Error::Invalid(
                    "cell heat must be finite, unique and attached to positive liquid mass",
                ));
            }
            *heated.get_mut(&location).expect("accepted wet cell") += source.energy_j / mass;
            external += source.energy_j;
        }
        if !external.is_finite() {
            return Err(Error::Invalid("total supplied heat overflows"));
        }
        let (_, tau) = self.resolve_states(&heated, &masses, gate)?;
        // Native activation/conversion and the companion solve operate on one
        // private candidate. No accepted populations, carry or source state move
        // until every endpoint query and conservation gate has passed.
        let mut lattice = self.flow.lattice.clone();
        let trace = lattice.step_with_mass_transport(gate, &tau, self.controls.max_mass_routes)?;
        admit_masses(&trace.after, self.controls.max_cells)?;
        let (mut h, outcome) = transport::advect(&trace, &heated, self.controls, gate)?;
        h.retain(|key, _| trace.after.contains_key(key));
        extend_empty_ghosts(&trace.after, &mut h, &heated)?;
        let (states, _) = self.resolve_states(&h, &trace.after, gate)?;
        let after = total(&trace.after, &h, mass_unit)?;
        let defect = after - before - external;
        let scale = absolute_total(&masses, &self.h, mass_unit)?
            + absolute_total(&trace.after, &h, mass_unit)?
            + heat.iter().map(|q| q.energy_j.abs()).sum::<f64>();
        if !defect.is_finite()
            || !scale.is_finite()
            || defect.abs() > self.controls.relative_tolerance * scale
        {
            return Err(Error::Invalid(
                "extensive enthalpy does not close to the declared tolerance",
            ));
        }
        let mass_before_kg = trace.before.values().sum::<f64>() * mass_unit;
        let mass_after_kg = trace.after.values().sum::<f64>() * mass_unit;
        if !mass_before_kg.is_finite() || !mass_after_kg.is_finite() {
            return Err(Error::Invalid("physical liquid mass overflows"));
        }
        let report = FreeSurfaceEnthalpyReport {
            mass_before_kg,
            mass_after_kg,
            enthalpy_before_j: before,
            enthalpy_after_j: after,
            external_heat_j: external,
            balance_defect_j: defect,
            relative_transport_residual: outcome.relative_residual,
            transport_iterations: outcome.iterations,
            mass_routes: trace.routes.len(),
        };
        poll(gate)?;
        self.flow.lattice = lattice;
        self.h = h;
        self.states = states;
        Ok(report)
    }

    fn location_at(&self, coordinates: [usize; 3]) -> Option<MassLocation3> {
        let key = FreeSurface3::transport_location(coordinates)?;
        self.h.contains_key(&key).then_some(key)
    }

    fn resolve_states(
        &self,
        h: &BTreeMap<MassLocation3, f64>,
        masses: &BTreeMap<MassLocation3, f64>,
        gate: &CancelGate,
    ) -> Result<
        (
            BTreeMap<MassLocation3, ResolvedLiquidState>,
            BTreeMap<(u64, usize), f64>,
        ),
        Error,
    > {
        let mut states = BTreeMap::new();
        let mut tau = BTreeMap::new();
        let rho = self.flow.scaling.material().density_kg_m3();
        for (&location, &mass) in masses {
            poll(gate)?;
            if location == MassLocation3::Carry && mass == 0.0 {
                continue;
            }
            let specific = *h
                .get(&location)
                .ok_or(Error::Invalid("wet mass has no enthalpy history"))?;
            let phase = self.caloric.state_at_specific_enthalpy(specific)?;
            if phase.phase() != SolidLiquidPhase::Liquid
                || phase.bulk_density_kg_m3().to_bits() != rho.to_bits()
            {
                return Err(Error::Invalid(
                    "free-surface heat needs pure liquid at the original fixed density",
                ));
            }
            let state = self
                .law
                .at_temperature(Temperature::new(phase.temperature_k()))?;
            if state.density_kg_m3().to_bits() != rho.to_bits() {
                return Err(Error::Invalid(
                    "temperature-dependent density needs a conservative volume-change owner",
                ));
            }
            let (plan, _) = plan_momentum_scaling(
                state.kinematic_viscosity_m2_s(),
                self.flow.scaling.cell_size_m(),
                self.flow.scaling.timestep_s(),
                self.flow.scaling.flow.u_lattice * self.flow.scaling.velocity_scale,
            )?;
            if let MassLocation3::Cell(key, lane) = location {
                tau.insert((key, lane), plan.tau);
            }
            states.insert(location, state);
        }
        Ok((states, tau))
    }
}

fn poll(gate: &CancelGate) -> Result<(), Error> {
    if gate.is_requested() {
        Err(FreeSurfaceError3::Cancelled.into())
    } else {
        Ok(())
    }
}

fn admit_masses(masses: &BTreeMap<MassLocation3, f64>, max_cells: usize) -> Result<(), Error> {
    if masses.len() > max_cells || masses.values().any(|m| !m.is_finite() || *m < 0.0) {
        return Err(Error::Invalid(
            "enthalpy storage requires bounded nonnegative native cell and carry mass",
        ));
    }
    if !masses.values().any(|m| *m > 0.0) {
        return Err(Error::Invalid(
            "enthalpy transport has no positive liquid mass",
        ));
    }
    Ok(())
}

fn total(
    masses: &BTreeMap<MassLocation3, f64>,
    h: &BTreeMap<MassLocation3, f64>,
    unit: f64,
) -> Result<f64, Error> {
    let mut result = 0.0;
    for (&key, &mass) in masses {
        if mass == 0.0 {
            continue;
        }
        result += mass
            * *h.get(&key)
                .ok_or(Error::Invalid("stored mass has no enthalpy"))?
            * unit;
    }
    if result.is_finite() {
        Ok(result)
    } else {
        Err(Error::Invalid("extensive enthalpy overflows"))
    }
}

fn absolute_total(
    masses: &BTreeMap<MassLocation3, f64>,
    h: &BTreeMap<MassLocation3, f64>,
    unit: f64,
) -> Result<f64, Error> {
    let absolute: BTreeMap<_, _> = h.iter().map(|(&key, value)| (key, value.abs())).collect();
    total(masses, &absolute, unit)
}

fn extend_empty_ghosts(
    masses: &BTreeMap<MassLocation3, f64>,
    h: &mut BTreeMap<MassLocation3, f64>,
    old_h: &BTreeMap<MassLocation3, f64>,
) -> Result<(), Error> {
    let donors: BTreeMap<_, _> = masses
        .iter()
        .filter_map(|(&key, &mass)| {
            if mass > 0.0 {
                FreeSurface3::transport_coordinates(key).map(|p| (p, (mass, h[&key])))
            } else {
                None
            }
        })
        .collect();
    for (&key, &mass) in masses {
        if mass != 0.0 || key == MassLocation3::Carry || h.contains_key(&key) {
            continue;
        }
        let p = FreeSurface3::transport_coordinates(key).expect("cell");
        let mut energy = 0.0;
        let mut weight = 0.0;
        for &(dx, dy, dz) in &crate::d3q19::E3[1..] {
            let q = [
                p[0] as i64 + i64::from(dx),
                p[1] as i64 + i64::from(dy),
                p[2] as i64 + i64::from(dz),
            ];
            if q.iter().any(|v| *v < 0) {
                continue;
            }
            if let Some(&(m, specific)) = donors.get(&q.map(|v| v as usize)) {
                energy += m * specific;
                weight += m;
            }
        }
        let value = if weight > 0.0 {
            energy / weight
        } else {
            *old_h.get(&key).ok_or(Error::Invalid(
                "new empty interface has no supported thermal neighbor",
            ))?
        };
        h.insert(key, value);
    }
    Ok(())
}
