//! Moving-liquid extensive enthalpy on the original VOF mass routes.
//!
//! Each step applies declared heat to the accepted liquid, refreshes the SAME
//! source claims at its caloric temperature, then runs native collision with
//! that cell's viscosity. Backward-Euler donor transport carries `mass*h` over
//! the native face and conversion/carry routes. The existing sparse Krylov
//! owner solves that transport; no velocity-derived second mass flux exists.
//!
//! This is the single-liquid, constant-density, stationary-wall rung. A sourced
//! validity-wide constant surface tension retains the native capillary mechanics.
//! Capillary work and surface-energy accuracy remain the native model's scope;
//! no surface heat or total-energy closure is invented here. Optional
//! `thermocapillary` admits source-varying sigma on checked planar, height-graph
//! and detached local surfaces.
//! Optional `wetting` evolves the original planar contact law using
//! local liquid temperature. Phase mobility and pressure work are not supplied. Optional `viscous`
//! converts the molecular stress power into heat on the same VOF mass routes.
//! Optional `conduction` reuses the shared thermal timeline on explicit wetted
//! face geometry. Optional `buoyancy` admits bounded source-density departures
//! while retaining reference mass and volume. In particular, h is enthalpy, not internal or total energy. The
//! incompressible lattice pressure is not a thermodynamic EOS. A negative native
//! accepted interface/carry mass refuses rather than becoming thermal storage.
//! A source-domain, numerical-budget or cancellation refusal publishes nothing.

use std::collections::{BTreeMap, BTreeSet};

use fs_exec::CancelGate;
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, SURFACE_TENSION_DIMS,
    SURFACE_TENSION_PROPERTY,
};
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, PhaseStateError, SolidLiquidPhase};
use fs_material::phase::conductivity::LiquidMassFractionConductivity;
use fs_material::state_point::MaterialStatePointError;
use fs_qty::Temperature;

use super::{MaterialFreeSurface3, MaterialLbmError};
use super::wetting::{FreeSurfaceWettingError, WettingBinding};
use crate::d3q19::freesurface3::transport::MassLocation3;
use crate::d3q19::freesurface3::{FreeSurface3, FreeSurfaceError3};
use crate::material::plan_momentum_scaling;

mod transport;

/// Conservative condensation/evaporation ingress on existing native interfaces.
pub mod mass_exchange;
pub use mass_exchange::{
    FreeSurfaceMassExchange, FreeSurfaceMassExchangeControls, FreeSurfaceMassExchangeReport,
    FreeSurfaceMassExchangeSample,
};

mod phase_conduction;
pub use phase_conduction::FreeSurfaceConductivity;

/// Source-density buoyancy with explicitly bounded reference-mass storage.
pub mod buoyancy;
use buoyancy::{BoussinesqLiquid, resolve_liquid};

/// Molecular stress power deposited into the moving liquid's finite enthalpy.
pub mod viscous;

/// Original T/shear source tables on actual native collision and molecular heat.
pub mod rheology;

/// Local temperature and contact-line history on the original planar boundary.
pub mod wetting;

/// Source-driven tangential stress on an explicitly bounded planar free surface.
pub mod thermocapillary;

/// Source-conductivity exchange on the existing thermal timeline before flow.
#[cfg(feature = "enthalpy-coupling")]
pub mod conduction;

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
    /// Absolute L2 residual allowance [J/kg] for the enthalpies actually
    /// stored after restoring their datum. The Krylov iteration still uses
    /// `relative_tolerance`; zero opts into a strict relative endpoint check.
    /// Must be finite and nonnegative.
    pub absolute_transport_tolerance_j_kg: f64,
}

impl Default for FreeSurfaceEnthalpyControls {
    fn default() -> Self {
        Self {
            max_cells: 16_384,
            max_mass_routes: 400_000,
            max_iterations: 2_000,
            relative_tolerance: 1e-11,
            absolute_transport_tolerance_j_kg: 1e-8,
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
    /// Admitted external cell heat, native Joule heat and any explicitly coupled
    /// environmental reservoir supply [J]. Internal wall and molecular heat stay separate.
    pub external_heat_j: f64,
    /// Constitutive molecular stress work converted into heat [J]. Zero for
    /// steps that do not explicitly select the viscous-heating entry point.
    pub viscous_heat_j: f64,
    /// Signed heat received from explicitly coupled finite solid cells [J].
    /// This is internal exchange when the accounting boundary includes them.
    pub body_heat_j: f64,
    /// `H_after - H_before - external_heat - viscous_heat - body_heat` [J]. No balancing
    /// heat or numerical kinetic-energy residual is added.
    pub balance_defect_j: f64,
    /// True relative residual of the component-centered enthalpy system,
    /// recomputed after restoring the datum to the enthalpies actually stored.
    /// An arbitrary caloric datum is excluded from its normalization; native
    /// mass-closure defects remain. May exceed the relative tolerance only
    /// within the separately declared absolute transport allowance.
    pub relative_transport_residual: f64,
    /// Actual L2 norm of the represented, row-scaled transport residual [J/kg].
    /// Endpoint acceptance requires this to be at most the absolute transport
    /// tolerance plus `relative_tolerance * norm(centered_rhs)`.
    pub absolute_transport_residual_j_kg: f64,
    /// Actual sparse transport iterations.
    pub transport_iterations: usize,
    /// Count of native nonzero mass routes, including conversion redistribution.
    pub mass_routes: usize,
    /// Actual native motion on an imposed-B core, separately observed from
    /// the electrical solver's sampled Lorentz work. No kinetic defect is heat.
    #[cfg(feature = "enthalpy-coupling")]
    pub magnetic_motion: Option<conduction::electrical::FreeSurfaceMagneticMotion>,
}

/// Typed source, native-flow, numerical or physical-domain refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum FreeSurfaceEnthalpyError {
    /// The original material/free-surface scaling or flow owner refused.
    Flow(MaterialLbmError),
    /// Original source query refusal, including pinned-claim support.
    Material(MaterialStatePointError),
    /// The original planar source/geometry/hysteresis owner refused.
    Wetting(FreeSurfaceWettingError),
    /// Original finite-body material/temperature binding refused.
    ThermalMaterial(fs_material::thermal::source::ThermalMaterialError),
    /// Original caloric chart refusal; no temperature extrapolation occurs.
    Phase(PhaseStateError),
    /// The existing sparse transport solver refused.
    Solve(crate::conjugate::ChtError),
    /// The shared molecular viscosity/shear solve exhausted its declared budget.
    RheologyConvergence {
        iterations: usize,
        relaxation_defect: f64,
    },
    /// The existing conservative thermal timeline refused conduction.
    #[cfg(feature = "enthalpy-coupling")]
    Thermal(fs_conduction::enthalpy::DepositionError),
    /// The original finite-atmosphere internal-energy owner refused.
    #[cfg(feature = "enthalpy-coupling")]
    Gas(fs_material::thermochemical_gas::energy::GasEnergyError),
    /// The original scalar electrical source could not adopt the liquid bundle.
    #[cfg(feature = "enthalpy-coupling")]
    Conductor(fs_material::conductor::ConductorError),
    /// The original implicit thermal/nonlinear owners refused surface exchange.
    #[cfg(feature = "implicit-enthalpy-coupling")]
    ImplicitSurface(conduction::surface::implicit::FreeSurfaceImplicitFailure),
    /// Missing physical closure, invalid input or exhausted explicit budget.
    Invalid(&'static str),
}

impl core::fmt::Display for FreeSurfaceEnthalpyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Flow(e) => write!(f, "free-surface enthalpy: {e}"),
            Self::Material(e) => write!(f, "free-surface enthalpy source: {e}"),
            Self::Wetting(e) => write!(f, "free-surface enthalpy wetting: {e}"),
            Self::ThermalMaterial(e) => write!(f, "free-surface thermal body: {e}"),
            Self::Phase(e) => write!(f, "free-surface enthalpy chart: {e}"),
            Self::Solve(e) => write!(f, "free-surface enthalpy solve: {e}"),
            Self::RheologyConvergence {
                iterations,
                relaxation_defect,
            } => write!(
                f,
                "free-surface rheology did not converge in {iterations} queries (tau defect {relaxation_defect})",
            ),
            #[cfg(feature = "enthalpy-coupling")]
            Self::Thermal(e) => write!(f, "free-surface conduction: {e}"),
            #[cfg(feature = "enthalpy-coupling")]
            Self::Gas(e) => write!(f, "free-surface finite atmosphere: {e}"),
            #[cfg(feature = "enthalpy-coupling")]
            Self::Conductor(e) => write!(f, "free-surface electrical source: {e}"),
            #[cfg(feature = "implicit-enthalpy-coupling")]
            Self::ImplicitSurface(e) => write!(f, "free-surface implicit heat: {e}"),
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
impl From<FreeSurfaceWettingError> for FreeSurfaceEnthalpyError {
    fn from(e: FreeSurfaceWettingError) -> Self {
        Self::Wetting(e)
    }
}
impl From<fs_material::thermal::source::ThermalMaterialError> for FreeSurfaceEnthalpyError {
    fn from(e: fs_material::thermal::source::ThermalMaterialError) -> Self {
        Self::ThermalMaterial(e)
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
#[cfg(feature = "enthalpy-coupling")]
impl From<fs_conduction::enthalpy::DepositionError> for FreeSurfaceEnthalpyError {
    fn from(e: fs_conduction::enthalpy::DepositionError) -> Self {
        if e == fs_conduction::enthalpy::DepositionError::Cancelled {
            FreeSurfaceError3::Cancelled.into()
        } else {
            Self::Thermal(e)
        }
    }
}

#[cfg(feature = "enthalpy-coupling")]
impl From<fs_material::thermochemical_gas::energy::GasEnergyError> for FreeSurfaceEnthalpyError {
    fn from(error: fs_material::thermochemical_gas::energy::GasEnergyError) -> Self {
        Self::Gas(error)
    }
}

type Error = FreeSurfaceEnthalpyError;

/// Native moving liquid, original caloric chart and pinned source law.
///
/// `h` is retained for every wet collision lane, including zero-mass interface
/// ghosts; only `m*h` contributes storage. Unconnected zero-mass ghosts receive
/// a mass-weighted one-ring temperature extension, without moving any energy.
/// A zero-mass carry has no assigned enthalpy unless it participates in a route.
#[derive(Clone)]
pub struct MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    flow: MaterialFreeSurface3<'s>,
    law: &'m LiquidTemperatureLaw<'m>,
    caloric: &'m EquilibriumEnthalpyPhaseCurve,
    controls: FreeSurfaceEnthalpyControls,
    time_origin_s: f64,
    accepted_steps: u64,
    h: BTreeMap<MassLocation3, f64>,
    states: BTreeMap<MassLocation3, ResolvedLiquidState>,
    phase_conduction: Option<LiquidMassFractionConductivity>,
    buoyancy: Option<BoussinesqLiquid>,
    rheology: Option<rheology::RheologyBinding>,
    rheology_samples: BTreeMap<MassLocation3, rheology::FreeSurfaceRheologySample>,
    wetting: Option<WettingBinding<'m>>,
    thermocapillary: Option<thermocapillary::ThermocapillarySurface>,
    // Once surface T drives stress, a later step must retain its full boundary.
    // This flag publishes only with the complete native thermal transaction.
    surface_temperature_forces: bool,
    #[cfg(feature = "enthalpy-coupling")]
    electrical: Option<conduction::electrical::ElectricalCore>,
}

type SurfaceTensionResolver<'a> = dyn FnMut(
    &FreeSurface3,
    &BTreeMap<MassLocation3, ResolvedLiquidState>,
    &CancelGate,
) -> Result<BTreeMap<(u64, usize), f64>, Error> + 'a;

struct ResolvedFlowMaterials {
    states: BTreeMap<MassLocation3, ResolvedLiquidState>,
    tau: BTreeMap<(u64, usize), f64>,
    forces: Option<BTreeMap<(u64, usize), [f64; 3]>>,
    surface_tensions: Option<BTreeMap<(u64, usize), f64>>,
    rheology: BTreeMap<MassLocation3, rheology::FreeSurfaceRheologySample>,
}

impl<'s> MaterialFreeSurface3<'s> {
    /// Bind liquid enthalpy to this actual VOF state without changing its mass
    /// or populations. The callback supplies specific enthalpy [J/kg] for every
    /// initial wet lane, including any zero-mass interface ghost.
    ///
    /// The caloric chart and liquid law must name the same original material
    /// card as the surface source and retain its viscosity basis and complete
    /// initial query context. Capillarity must be explicitly neglected or use a
    /// source-declared constant sigma retained in the SAME full thermal bundle.
    /// Equal values sampled from a varying source do not establish constancy.
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
        self.with_enthalpy_model(law, caloric, controls, initial_h, None, None, None, gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn with_enthalpy_model<'m>(
        self,
        law: &'m LiquidTemperatureLaw<'m>,
        caloric: &'m EquilibriumEnthalpyPhaseCurve,
        controls: FreeSurfaceEnthalpyControls,
        initial_h: impl Fn([usize; 3]) -> f64,
        buoyancy: Option<crate::material::transport::boussinesq::BoussinesqConfig<3>>,
        rheology: Option<rheology::RheologyBinding>,
        thermocapillary: Option<thermocapillary::ThermocapillarySurface>,
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
            || !controls.absolute_transport_tolerance_j_kg.is_finite()
            || controls.absolute_transport_tolerance_j_kg < 0.0
        {
            return Err(Error::Invalid(
                "invalid bounded enthalpy transport controls",
            ));
        }
        let source = self.scaling.material();
        if source.rheology_axes().is_some() != rheology.is_some() {
            return Err(Error::Invalid(
                "evolving liquid source requires its runtime rheology owner",
            ));
        }
        if self.scaling.planar_wetting.is_some() {
            return Err(Error::Invalid(
                "thermal free surface must bind wetting through its evolving interface owner",
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
        admit_capillary_source(source, law, thermocapillary.is_some())?;
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
        let buoyancy = buoyancy
            .map(|config| BoussinesqLiquid::try_new(&self, law, config))
            .transpose()?;
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
            time_origin_s: 0.0,
            accepted_steps: 0,
            h,
            states: BTreeMap::new(),
            phase_conduction: None,
            buoyancy,
            rheology,
            rheology_samples: BTreeMap::new(),
            wetting: None,
            thermocapillary,
            surface_temperature_forces: false,
            #[cfg(feature = "enthalpy-coupling")]
            electrical: None,
        };
        let resolved =
            candidate.resolve_states(&candidate.flow.lattice, &candidate.h, &masses, gate)?;
        if let Some(surface_tensions) = resolved.surface_tensions {
            candidate.flow.lattice.rebind_surface_tensions(surface_tensions, gate)?;
        }
        if let Some(forces) = resolved.forces {
            candidate.flow.lattice.rebind_force_field(forces, gate)?;
        }
        candidate.states = resolved.states;
        candidate.rheology_samples = resolved.rheology;
        Ok(candidate)
    }
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    /// Set the physical time at which this thermal/native trajectory was bound.
    /// Populations, mass, material state and enthalpy are unchanged. A trajectory
    /// that has already advanced cannot have its physical clock reset.
    pub fn with_time_origin(mut self, time_origin_s: f64) -> Result<Self, Error> {
        if !time_origin_s.is_finite() || time_origin_s < 0.0 || self.accepted_steps != 0 {
            return Err(Error::Invalid("time origin requires an unadvanced liquid and finite nonnegative time"));
        }
        self.time_origin_s = time_origin_s;
        Ok(self)
    }

    /// Actual accepted VOF/enthalpy intervals since the original binding.
    /// Thermal subcycles and rejected candidates do not increment this count.
    pub const fn accepted_steps(&self) -> u64 {
        self.accepted_steps
    }

    /// Physical time retained by every original heat, gas, wall and electrical API.
    pub fn physical_time_s(&self) -> f64 {
        (self.accepted_steps as f64).mul_add(self.flow.scaling.timestep_s(), self.time_origin_s)
    }

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
        let (heated, external, absolute_heat) = self.prepare_heat(heat, gate)?;
        self.advance_prepared_enthalpy(heated, external, 0.0, 0.0, absolute_heat, None, gate)
    }

    fn prepare_heat(
        &self,
        heat: &[FreeSurfaceCellHeat],
        gate: &CancelGate,
    ) -> Result<(BTreeMap<MassLocation3, f64>, f64, f64), Error> {
        poll(gate)?;
        #[cfg(feature = "enthalpy-coupling")]
        if self.electrical.is_some() {
            return Err(Error::Invalid(
                "an electrical core requires a conduction timeline; its drive cannot be skipped",
            ));
        }
        if heat.len() > self.controls.max_cells {
            return Err(Error::Invalid("cell-heat budget exhausted"));
        }
        let masses = self.flow.lattice.mass_snapshot();
        let mass_unit = self.flow.scaling.mass_unit_kg();
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
        Ok((
            heated,
            external,
            heat.iter().map(|q| q.energy_j.abs()).sum(),
        ))
    }

    // Both external heating and source-conductivity exchange feed the exact
    // same native collision/advection/acceptance owner. No thermal stage can
    // publish before this complete mass and energy candidate succeeds.
    #[allow(clippy::too_many_arguments)]
    fn advance_prepared_enthalpy(
        &mut self,
        heated: BTreeMap<MassLocation3, f64>,
        external: f64,
        viscous: f64,
        body_heat: f64,
        absolute_heat: f64,
        additional_forces: Option<&BTreeMap<MassLocation3, [f64; 3]>>,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceEnthalpyReport, Error> {
        self.advance_prepared_enthalpy_with_surface(
            heated, external, viscous, body_heat, absolute_heat,
            additional_forces, None, None, gate,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_prepared_enthalpy_with_surface(
        &mut self,
        heated: BTreeMap<MassLocation3, f64>,
        external: f64,
        viscous: f64,
        body_heat: f64,
        absolute_heat: f64,
        additional_forces: Option<&BTreeMap<MassLocation3, [f64; 3]>>,
        mut surface: Option<&mut SurfaceTensionResolver<'_>>,
        wetting_outer_gas: Option<&fs_material::thermochemical_gas::ResolvedThermochemicalGasState>,
        gate: &CancelGate,
    ) -> Result<FreeSurfaceEnthalpyReport, Error> {
        poll(gate)?;
        let next_step = self.accepted_steps.checked_add(1)
            .filter(|count| *count <= (1_u64 << 53))
            .ok_or(Error::Invalid("native liquid clock exhausted"))?;
        let next_time = (next_step as f64).mul_add(self.flow.scaling.timestep_s(), self.time_origin_s);
        if !next_time.is_finite() || next_time <= self.physical_time_s() {
            return Err(Error::Invalid("native liquid time must advance representably"));
        }
        if self.surface_temperature_forces && surface.is_none() {
            return Err(Error::Invalid(
                "surface-temperature stress requires its complete solved heat boundary",
            ));
        }
        // Local heat also refreshes geometry for constant sigma. Only a
        // temperature-dependent force field needs that boundary on later steps.
        let uses_surface_temperature =
            surface.is_some() && self.temperature_dependent_thermocapillarity();
        let masses = self.flow.lattice.mass_snapshot();
        let mass_unit = self.flow.scaling.mass_unit_kg();
        let before = total(&masses, &self.h, mass_unit)?;
        let mut resolved = self.resolve_states_with_surface(
            &self.flow.lattice, &heated, &masses, &mut surface, gate,
        )?;
        if let Some(additional) = additional_forces {
            let reference = self.flow.lattice.reference_force();
            let forces = resolved.forces.get_or_insert_with(|| masses.keys().filter_map(|location| {
                if let MassLocation3::Cell(key, lane) = *location { Some(((key, lane), reference)) }
                else { None }
            }).collect());
            for (&location, extra) in additional {
                poll(gate)?;
                let MassLocation3::Cell(key, lane) = location else {
                    return Err(Error::Invalid("unlocated carry cannot receive a body force"));
                };
                let force = forces.get_mut(&(key, lane))
                    .ok_or(Error::Invalid("additional body force requires an accepted wet cell"))?;
                for axis in 0..3 {
                    force[axis] += extra[axis];
                    if !force[axis].is_finite() {
                        return Err(Error::Invalid("combined native body force is not finite"));
                    }
                }
            }
        }
        // Native activation/conversion and the companion solve operate on one
        // private candidate. No accepted populations, carry or source state move
        // until every endpoint query and conservation gate has passed.
        let mut lattice = self.flow.lattice.clone();
        if let Some(surface_tensions) = resolved.surface_tensions {
            lattice.rebind_surface_tensions(surface_tensions, gate)?;
        }
        let prepared_wetting = self.prepare_wetting(
            &mut lattice, &resolved.states, wetting_outer_gas, gate,
        )?;
        let mut trace = lattice.step_with_mass_transport(
            gate,
            &resolved.tau,
            resolved.forces,
            self.controls.max_mass_routes,
        )?;
        if self.thermocapillary.is_some() {
            // Resolve native conversion hysteresis through the SAME mass
            // routes before reconstructing geometry or advecting enthalpy.
            // The candidate stays private until all endpoint gates pass.
            lattice.resolve_interface_fills(&mut trace, gate)?;
        }
        admit_masses(&trace.after, self.controls.max_cells)?;
        let (mut h, outcome) = transport::advect(&trace, &heated, self.controls, gate)?;
        h.retain(|key, _| trace.after.contains_key(key));
        extend_empty_ghosts(&trace.after, &mut h, &heated)?;
        let endpoint = self.resolve_states_with_surface(
            &lattice, &h, &trace.after, &mut surface, gate,
        )?;
        if let Some(surface_tensions) = endpoint.surface_tensions {
            lattice.rebind_surface_tensions(surface_tensions, gate)?;
        }
        #[cfg(feature = "enthalpy-coupling")]
        if let Some(core) = &self.electrical {
            core.admit(&lattice, &endpoint.states, gate)?;
        }
        let endpoint_forces = endpoint.forces.or_else(|| additional_forces.map(|_| {
            let reference = lattice.reference_force();
            trace.after.keys().filter_map(|location| {
                if let MassLocation3::Cell(key, lane) = *location { Some(((key, lane), reference)) }
                else { None }
            }).collect()
        }));
        if let Some(forces) = endpoint_forces {
            // This changes the Guo population frame without adding an impulse.
            // The accepted force now belongs to the accepted advected h state.
            lattice.rebind_force_field(forces, gate)?;
        }
        #[cfg(feature = "enthalpy-coupling")]
        let magnetic_motion = self.electrical.as_ref().map(|core|
            core.admit_motion(&self.flow.lattice, &lattice, &self.flow.scaling, gate)
        ).transpose()?.flatten();
        let accepted_wetting = match prepared_wetting {
            Some(binding) => {
                let temperatures = wetting::temperatures(
                    &endpoint.states, self.law.temperature_axis(), gate,
                )?;
                Some(binding.at_endpoint(
                    &mut lattice,
                    &self.flow.scaling,
                    Some((self.law.temperature_axis(), &temperatures)),
                    gate,
                )?)
            }
            None => None,
        };
        let endpoint_masses = lattice.mass_snapshot();
        admit_masses(&endpoint_masses, self.controls.max_cells)?;
        let after = total(&endpoint_masses, &h, mass_unit)?;
        let defect = after - before - external - viscous - body_heat;
        let scale = absolute_total(&masses, &self.h, mass_unit)?
            + absolute_total(&endpoint_masses, &h, mass_unit)?
            + absolute_heat;
        if !defect.is_finite()
            || !scale.is_finite()
            || defect.abs() > self.controls.relative_tolerance * scale
        {
            return Err(Error::Invalid(
                "extensive enthalpy does not close to the declared tolerance",
            ));
        }
        let mass_before_kg = trace.before.values().sum::<f64>() * mass_unit;
        let mass_after_kg = endpoint_masses.values().sum::<f64>() * mass_unit;
        if !mass_before_kg.is_finite() || !mass_after_kg.is_finite() {
            return Err(Error::Invalid("physical liquid mass overflows"));
        }
        let report = FreeSurfaceEnthalpyReport {
            mass_before_kg,
            mass_after_kg,
            enthalpy_before_j: before,
            enthalpy_after_j: after,
            external_heat_j: external,
            viscous_heat_j: viscous,
            body_heat_j: body_heat,
            balance_defect_j: defect,
            relative_transport_residual: outcome.relative_residual,
            absolute_transport_residual_j_kg: outcome.absolute_residual_j_kg,
            transport_iterations: outcome.iterations,
            mass_routes: trace.routes.len(),
            #[cfg(feature = "enthalpy-coupling")]
            magnetic_motion,
        };
        poll(gate)?;
        self.flow.lattice = lattice;
        self.h = h;
        self.states = endpoint.states;
        self.rheology_samples = endpoint.rheology;
        self.wetting = accepted_wetting;
        self.surface_temperature_forces |= uses_surface_temperature;
        self.accepted_steps = next_step;
        Ok(report)
    }

    fn location_at(&self, coordinates: [usize; 3]) -> Option<MassLocation3> {
        let key = FreeSurface3::transport_location(coordinates)?;
        self.h.contains_key(&key).then_some(key)
    }

    fn resolve_states(
        &self,
        lattice: &FreeSurface3,
        h: &BTreeMap<MassLocation3, f64>,
        masses: &BTreeMap<MassLocation3, f64>,
        gate: &CancelGate,
    ) -> Result<ResolvedFlowMaterials, Error> {
        self.resolve_states_with_surface(lattice, h, masses, &mut None, gate)
    }

    fn resolve_states_with_surface(
        &self,
        lattice: &FreeSurface3,
        h: &BTreeMap<MassLocation3, f64>,
        masses: &BTreeMap<MassLocation3, f64>,
        surface: &mut Option<&mut SurfaceTensionResolver<'_>>,
        gate: &CancelGate,
    ) -> Result<ResolvedFlowMaterials, Error> {
        let mut states = BTreeMap::new();
        let mut tau = BTreeMap::new();
        let mut forces = self.buoyancy.map(|_| BTreeMap::new());
        let rho = self.flow.scaling.material().density_kg_m3();
        let rheology_model = self.rheology_model();
        let mut rheology = BTreeMap::new();
        for (&location, &mass) in masses {
            poll(gate)?;
            if location == MassLocation3::Carry && mass == 0.0 {
                continue;
            }
            let specific = *h
                .get(&location)
                .ok_or(Error::Invalid("wet mass has no enthalpy history"))?;
            let (state, local_tau) = if let Some(model) = &rheology_model {
                if location == MassLocation3::Carry {
                    // Carry has caloric storage but no resolved stress or
                    // collision. Assigning it a constitutive rate would invent
                    // a physical state; its receiving cell resolves the source.
                    buoyancy::liquid_temperature(self.caloric, specific, rho)?;
                    continue;
                }
                let sample =
                    model.sample(specific, self.rheology_stress(lattice, location)?, || {
                        poll(gate)
                    })?;
                rheology.insert(location, sample.sample);
                (sample.material, sample.sample.collision_tau)
            } else {
                let (state, transport) =
                    resolve_liquid(self.law, self.caloric, specific, rho, self.buoyancy)?;
                let (plan, _) = plan_momentum_scaling(
                    transport.kinematic_viscosity_m2_s(),
                    self.flow.scaling.cell_size_m(),
                    self.flow.scaling.timestep_s(),
                    self.flow.scaling.flow.u_lattice * self.flow.scaling.velocity_scale,
                )?;
                (state, plan.tau)
            };
            phase_conduction::resolve(&state, self.phase_conduction)?;
            if let MassLocation3::Cell(key, lane) = location {
                tau.insert((key, lane), local_tau);
                if let (Some(binding), Some(forces)) = (self.buoyancy, &mut forces) {
                    forces.insert((key, lane), binding.force(&self.flow, &state)?);
                }
            }
            states.insert(location, state);
        }
        let solved_tensions = match surface.as_mut() {
            Some(resolve) => Some(resolve(lattice, &states, gate)?),
            None => None,
        };
        let surface_tensions = self.resolve_thermocapillarity(
            lattice, &states, &mut forces, solved_tensions, gate,
        )?;
        Ok(ResolvedFlowMaterials {
            states,
            tau,
            forces,
            surface_tensions,
            rheology,
        })
    }
}

// The original immutable bundle refreshes this same selected claim at every
// thermal/rheological stage. Native scaling already owns the admitted sigma;
// there is no second capillary state or independently queried source epoch.
fn admit_capillary_source(
    source: &ResolvedLiquidSurfaceState,
    law: &LiquidTemperatureLaw<'_>,
    evolving: bool,
) -> Result<(), Error> {
    if source.capillarity() == LiquidCapillarity::Neglected {
        if evolving {
            return Err(Error::Invalid("thermocapillarity requires an explicitly sourced surface tension"));
        }
        return Ok(());
    }
    let original = source.material().property(SURFACE_TENSION_PROPERTY)
        .ok_or(Error::Invalid("missing original capillary source"))?;
    let retained = law.reference_state().material().property(SURFACE_TENSION_PROPERTY)
        .ok_or(Error::Invalid("thermal law must retain the original capillary source"))?;
    for property in [original, retained] {
        let constant = property.source_axes().is_empty()
            && property.answer().receipt.decision == fs_matdb::EvaluationDecision::ConstantWithinValidity;
        let thermal = evolving && !property.source_axes().is_empty()
            && property.source_axes().iter().all(|(axis, dims)|
                axis == law.temperature_axis() && *dims == fs_qty::Temperature::DIMS);
        if !(constant || thermal) || property.requirement().quantity()
                != fs_qty::QuantitySpec::dimensional(SURFACE_TENSION_DIMS)
            || !property.value_si().is_finite()
            || property.value_si() <= 0.0
        {
            return Err(Error::Invalid(
                "thermal sigma requires a source-declared constant or explicit surface temperature evolution",
            ));
        }
    }
    if original.requirement() != retained.requirement()
        || original.answer().receipt.selected != retained.answer().receipt.selected
        || original.value_si().to_bits() != retained.value_si().to_bits()
        || source.surface_tension_n_m().map(f64::to_bits)
            != Some(retained.value_si().to_bits())
    {
        return Err(Error::Invalid(
            "thermal law must retain the surface's exact sigma requirement and source claim",
        ));
    }
    Ok(())
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
