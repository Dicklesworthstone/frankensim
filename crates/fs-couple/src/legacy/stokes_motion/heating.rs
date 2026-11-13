//! Fixed-reference Stokes viscous heating of the original finite thermal cells.
//!
//! The native accepted endpoint supplies `2 mu D:D`; its per-triangle rates
//! are held over the same motion interval while the existing thermal owner
//! subcycles its authored sources, contacts, reservoirs and events. This is a
//! first-order, one-way split with frozen transport. It adds no advection,
//! pressure work, temperature-dependent viscosity or global passivity claim.
//! An explicit wall route also deposits the accepted Navier-slip dissipation.
//! Numerical flow damping and body dashpots are not deposited.
//! Guided and planar motion share cell admission, slip partition and thermal
//! stepping. Each publishes its fluid, body and heat histories together.

use std::collections::{BTreeMap, BTreeSet};

use fs_flux::{MaterialFluxError, MaterialStokesSolution, TriMesh};
use fs_material::fluid_transport::FluidTransport;
use fs_thermal::observers::{CoupledObservers, ObserverCouplingError};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig,
};

use super::planar::{PlanarStokesError, PlanarStokesStep};
use super::{
    GuidedStokesControls, GuidedStokesError, GuidedStokesMotion, GuidedStokesStep,
    MaterialStokesBoundary, nonnegative, poll, positive,
};

/// Explicit correspondence to the original x/y mesh extruded along world +z.
/// Each triangle owns exactly one distinct active thermal cell. Additional
/// thermal cells may carry supports and heat sinks through original links.
#[derive(Clone, Copy, Debug)]
pub struct StokesThermalCells<'a> {
    /// Thermal cell IDs in original `TriMesh::tris` order.
    pub triangle_cells: &'a [CellId],
    /// Extrusion start plane [m]; the end plane is this plus the motion depth.
    pub plane_z_m: f64,
}

/// Thermal ownership of one original Navier-slip boundary edge.
///
/// The fluid receiver is the mapped triangle adjacent to this edge. The wall
/// receiver must be an active cell outside the fluid map; multiple edges may
/// heat the same wall cell. Its geometry, caloric law and relationship to this
/// physical wall remain caller-authored. This route adds no conductive link.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StokesSlipHeatRoute {
    /// Original mesh boundary edge carrying `MaterialStokesBoundary::NavierSlip`.
    pub edge: usize,
    /// Original finite wall cell receiving the complementary heat fraction.
    pub wall_cell: CellId,
    /// Authored fraction of interfacial heat entering the fluid, in [0, 1].
    /// This is physical input, not inferred from the cells' masses or temperatures.
    pub fluid_heat_fraction: f64,
}

/// A refused coupled step changes neither accepted owner.
#[derive(Clone, Debug, PartialEq)]
pub enum GuidedStokesHeatError {
    /// Incomplete, ambiguous, geometrically inconsistent or unsynchronized map.
    Binding(&'static str),
    /// Original flow/body refusal, including cancellation between flow solves.
    Motion(GuidedStokesError),
    /// Original thermal material, budget, energy or cancellation refusal.
    Thermal(DepositionError),
}

impl std::fmt::Display for GuidedStokesHeatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Stokes thermal coupling refusal: {self:?}")
    }
}
impl std::error::Error for GuidedStokesHeatError {}
impl From<GuidedStokesError> for GuidedStokesHeatError {
    fn from(error: GuidedStokesError) -> Self {
        Self::Motion(error)
    }
}
impl From<DepositionError> for GuidedStokesHeatError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// A refused observed step preserves motion, thermal state and observer history.
#[derive(Clone, Debug, PartialEq)]
pub enum GuidedStokesObservedHeatError<E> {
    /// Original map, motion, thermal, heat-accounting or cancellation refusal.
    Heat(GuidedStokesHeatError),
    /// Refusal of the additional owner, retaining its original error type.
    Observer(E),
}

impl<E> From<GuidedStokesHeatError> for GuidedStokesObservedHeatError<E> {
    fn from(error: GuidedStokesHeatError) -> Self { Self::Heat(error) }
}
impl<E> From<GuidedStokesError> for GuidedStokesObservedHeatError<E> {
    fn from(error: GuidedStokesError) -> Self { Self::Heat(error.into()) }
}
impl<E> From<DepositionError> for GuidedStokesObservedHeatError<E> {
    fn from(error: DepositionError) -> Self { Self::Heat(error.into()) }
}
impl<E> From<StokesHeatFailure> for GuidedStokesObservedHeatError<E> {
    fn from(error: StokesHeatFailure) -> Self { Self::Heat(error.into()) }
}
impl<E> From<ObserverCouplingError<StokesHeatFailure, E>> for GuidedStokesObservedHeatError<E> {
    fn from(error: ObserverCouplingError<StokesHeatFailure, E>) -> Self {
        match error {
            ObserverCouplingError::Second(error) => Self::Observer(error),
            ObserverCouplingError::First(error) => Self::Heat(error.into()),
            ObserverCouplingError::Thermal(error) => Self::Heat(error.into()),
            ObserverCouplingError::ConflictingConductance(_) => Self::Heat(
                GuidedStokesHeatError::Binding("conflicting thermal-link conductance owners"),
            ),
            ObserverCouplingError::ConflictingReservoirConductance(_) => Self::Heat(
                GuidedStokesHeatError::Binding("conflicting reservoir-conductance owners"),
            ),
            ObserverCouplingError::ConflictingAmbientConductance(_) => Self::Heat(
                GuidedStokesHeatError::Binding("conflicting ambient-conductance owners"),
            ),
        }
    }
}
impl<E: std::fmt::Display> std::fmt::Display for GuidedStokesObservedHeatError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Heat(error) => std::fmt::Display::fmt(error, f),
            Self::Observer(error) => write!(f, "Stokes thermal observer refusal: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for GuidedStokesObservedHeatError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Heat(error) => Some(error),
            Self::Observer(error) => Some(error),
        }
    }
}

pub(super) fn unobserved_error(
    error: GuidedStokesObservedHeatError<DepositionError>,
) -> GuidedStokesHeatError {
    match error {
        GuidedStokesObservedHeatError::Heat(error) => error,
        GuidedStokesObservedHeatError::Observer(error) => error.into(),
    }
}

/// Receipts from the two original owners, committed at one identical time.
pub struct GuidedStokesHeatStep {
    /// Original body/fluid endpoint and its separate work accounts.
    pub motion: GuidedStokesStep,
    /// Native thermal audit, including viscous heat, authored sources and any
    /// additional observer's heat. Internal heat transfers remain conservative.
    pub thermal: EnergyAudit,
    /// Physical bulk heat actually supplied during native thermal substeps [J].
    /// Already included in `thermal.external_heat_j`; never add it twice.
    pub deposited_viscous_heat_j: f64,
    /// Physical interfacial Navier-slip heat supplied to fluid and wall cells [J].
    /// Zero for `step_with_heating`. Already included in `thermal.external_heat_j`.
    pub deposited_slip_heat_j: f64,
}

/// A refused planar fluid/body/thermal step preserves every accepted owner.
#[derive(Clone, Debug, PartialEq)]
pub enum PlanarStokesHeatError {
    /// Incomplete, ambiguous, geometrically inconsistent or unsynchronized map.
    Binding(&'static str),
    /// Original planar flow/body refusal, including composition cancellation.
    Motion(PlanarStokesError),
    /// Original thermal material, budget, energy or cancellation refusal.
    Thermal(DepositionError),
}

impl std::fmt::Display for PlanarStokesHeatError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "planar Stokes thermal coupling refusal: {self:?}")
    }
}
impl std::error::Error for PlanarStokesHeatError {}
impl From<PlanarStokesError> for PlanarStokesHeatError {
    fn from(error: PlanarStokesError) -> Self {
        Self::Motion(error)
    }
}
impl From<DepositionError> for PlanarStokesHeatError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

/// A refused observed planar step preserves every accepted physical history.
#[derive(Clone, Debug, PartialEq)]
pub enum PlanarStokesObservedHeatError<E> {
    /// Original map, motion, thermal, heat-accounting or cancellation refusal.
    Heat(PlanarStokesHeatError),
    /// Refusal of the additional owner, retaining its original error type.
    Observer(E),
}

impl<E> From<PlanarStokesHeatError> for PlanarStokesObservedHeatError<E> {
    fn from(error: PlanarStokesHeatError) -> Self { Self::Heat(error) }
}
impl<E> From<PlanarStokesError> for PlanarStokesObservedHeatError<E> {
    fn from(error: PlanarStokesError) -> Self { Self::Heat(error.into()) }
}
impl<E> From<DepositionError> for PlanarStokesObservedHeatError<E> {
    fn from(error: DepositionError) -> Self { Self::Heat(error.into()) }
}
impl<E> From<StokesHeatFailure> for PlanarStokesObservedHeatError<E> {
    fn from(error: StokesHeatFailure) -> Self { Self::Heat(error.into()) }
}
impl<E> From<ObserverCouplingError<StokesHeatFailure, E>> for PlanarStokesObservedHeatError<E> {
    fn from(error: ObserverCouplingError<StokesHeatFailure, E>) -> Self {
        match error {
            ObserverCouplingError::Second(error) => Self::Observer(error),
            ObserverCouplingError::First(error) => Self::Heat(error.into()),
            ObserverCouplingError::Thermal(error) => Self::Heat(error.into()),
            ObserverCouplingError::ConflictingConductance(_) => Self::Heat(
                PlanarStokesHeatError::Binding("conflicting thermal-link conductance owners"),
            ),
            ObserverCouplingError::ConflictingReservoirConductance(_) => Self::Heat(
                PlanarStokesHeatError::Binding("conflicting reservoir-conductance owners"),
            ),
            ObserverCouplingError::ConflictingAmbientConductance(_) => Self::Heat(
                PlanarStokesHeatError::Binding("conflicting ambient-conductance owners"),
            ),
        }
    }
}
impl<E: std::fmt::Display> std::fmt::Display for PlanarStokesObservedHeatError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Heat(error) => std::fmt::Display::fmt(error, f),
            Self::Observer(error) => write!(f, "planar Stokes thermal observer refusal: {error}"),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for PlanarStokesObservedHeatError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Heat(error) => Some(error),
            Self::Observer(error) => Some(error),
        }
    }
}

pub(super) fn unobserved_planar_error(
    error: PlanarStokesObservedHeatError<DepositionError>,
) -> PlanarStokesHeatError {
    match error {
        PlanarStokesObservedHeatError::Heat(error) => error,
        PlanarStokesObservedHeatError::Observer(error) => error.into(),
    }
}

/// Original planar momentum and finite thermal receipts at one accepted time.
pub struct PlanarStokesHeatStep {
    /// Original planar flow/body endpoint and its separate work accounts.
    pub motion: PlanarStokesStep,
    /// Native thermal audit, including viscous heat, authored sources and any
    /// additional observer's heat. Internal heat transfers remain conservative.
    pub thermal: EnergyAudit,
    /// Physical bulk heat supplied during native thermal substeps [J].
    /// Already included in `thermal.external_heat_j`; never add it twice.
    pub deposited_viscous_heat_j: f64,
    /// Physical Navier-slip heat supplied to fluid and wall cells [J].
    /// Zero for `step_with_heating`; included in `thermal.external_heat_j`.
    pub deposited_slip_heat_j: f64,
}

// Shared thermal failures retain each motion owner's original public error.
pub(crate) enum StokesHeatFailure {
    Binding(&'static str),
    Flow(MaterialFluxError),
    Thermal(DepositionError),
    Cancelled,
}
impl From<DepositionError> for StokesHeatFailure {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl From<MaterialFluxError> for StokesHeatFailure {
    fn from(error: MaterialFluxError) -> Self {
        Self::Flow(error)
    }
}
impl From<StokesHeatFailure> for GuidedStokesHeatError {
    fn from(error: StokesHeatFailure) -> Self {
        match error {
            StokesHeatFailure::Binding(message) => Self::Binding(message),
            StokesHeatFailure::Flow(error) => Self::Motion(GuidedStokesError::Flow(error)),
            StokesHeatFailure::Thermal(error) => Self::Thermal(error),
            StokesHeatFailure::Cancelled => Self::Motion(GuidedStokesError::Cancelled),
        }
    }
}
impl From<StokesHeatFailure> for PlanarStokesHeatError {
    fn from(error: StokesHeatFailure) -> Self {
        match error {
            StokesHeatFailure::Binding(message) => Self::Binding(message),
            StokesHeatFailure::Flow(error) => Self::Motion(PlanarStokesError::Flow(error)),
            StokesHeatFailure::Thermal(error) => Self::Thermal(error),
            StokesHeatFailure::Cancelled => Self::Motion(PlanarStokesError::Cancelled),
        }
    }
}

/// Frozen geometry, transport and start clock, without a mutable physical owner.
pub(crate) struct StokesHeatBinding<'m> {
    pub(crate) mesh: &'m TriMesh,
    pub(crate) density_kg_m3: f64,
    pub(crate) depth_m: f64,
    pub(crate) time_s: f64,
}

impl<F: FluidTransport + ?Sized> GuidedStokesMotion<'_, '_, '_, F> {
    /// Advance native guided flow/body motion and its finite thermal receivers.
    ///
    /// A receiver's centreline must run from the original triangle centroid at
    /// `plane_z_m` to that centroid plus the original depth along +z. Its volume
    /// and total surface area must describe that same triangular prism; mass
    /// must equal the frozen fluid density times that volume. Admission allows
    /// only floating-point roundoff, not a configurable geometric tolerance.
    /// The caller retains responsibility for choosing the same physical fluid's
    /// caloric law; a transport snapshot alone supplies no caloric identity.
    /// No thermal law, boundary source or link is replaced by this composition.
    ///
    /// Heat uses the original endpoint strain field throughout this interval.
    /// Refine both motion and thermal timesteps for temporal accuracy. Both
    /// clocks must initially agree exactly; the native thermal timeline advances
    /// to the exact proposed motion endpoint, including its scheduled events.
    ///
    /// # Errors
    /// Mapping, flow, mechanics, thermal-domain, energy, budget and cancellation
    /// failures preserve both owners, even after successful thermal substeps.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_heating(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GuidedStokesHeatStep, GuidedStokesHeatError> {
        self.step_heating(
            thermal,
            &mut (),
            cells,
            None,
            timestep_s,
            force_density,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            controls,
            thermal_config,
            cancelled,
        )
        .map_err(unobserved_error)
    }

    /// Advance with both bulk strain heat and physical Navier-slip heating.
    ///
    /// Supply exactly one route for EVERY exterior `NavierSlip` edge, including
    /// perfect-slip edges (whose dissipation is zero). Other boundary kinds and
    /// the no-slip moving body cannot be routed. The declared fluid/wall split
    /// is held over this interval; no flash-temperature or partition model is
    /// inferred. Ordinary thermal contacts still exchange heat independently.
    ///
    /// Heat is `depth * integral(beta * |(u-wall)_t|^2)` from the SAME accepted
    /// endpoint reaction that supplies physical wall effort. Neither material
    /// nor boundary callbacks are queried again. No-slip numerical penalties,
    /// body work and convective transport are not interfacial heat sources.
    /// The fixed-reference, frozen-transport and time-refinement requirements
    /// of [`Self::step_with_heating`] also apply.
    ///
    /// # Errors
    /// Incomplete/duplicate routes, inactive or fluid-owned wall receivers,
    /// invalid partitions, and original solver/thermal/cancellation failures
    /// preserve both owners. All route admission precedes the native flow solves.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_slip_heating(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        routes: &[StokesSlipHeatRoute],
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GuidedStokesHeatStep, GuidedStokesHeatError> {
        self.step_heating(
            thermal,
            &mut (),
            cells,
            Some(routes),
            timestep_s,
            force_density,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            controls,
            thermal_config,
            cancelled,
        )
        .map_err(unobserved_error)
    }

    /// Advance bulk viscous heating with another owner on the same thermal clock.
    ///
    /// Every observer hook is composed through the original thermal protocol,
    /// including source transitions, accepted-interval sources, conductance
    /// overrides, conservative exchanges and timestep restrictions. Nested
    /// `CoupledObservers` may supply several additional owners. Its heat is
    /// counted once in the native thermal audit; Stokes heat receipts continue
    /// to account only for the native viscous loss.
    ///
    /// The observer's `Clone` must isolate mutable physical history as required
    /// by `DepositionObserver`. Motion, thermal state and the cloned observer
    /// publish together only after all acceptance and cancellation checks pass.
    /// A refusal preserves all three inputs and retains the other owner's
    /// original error type. The fixed-reference, frozen-transport and
    /// time-refinement requirements of [`Self::step_with_heating`] still apply;
    /// additional thermal physics does not rebind the fluid's viscosity.
    ///
    /// # Errors
    /// Original map, motion, thermal or observer failures refuse the complete
    /// step. No accepted history changes, including on final cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_heating_observed<O: DepositionObserver>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut O,
        cells: StokesThermalCells<'_>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GuidedStokesHeatStep, GuidedStokesObservedHeatError<O::Error>> {
        self.step_heating(
            thermal, observer, cells, None, timestep_s, force_density, exterior,
            pressure_offsets_pa, external_force_n, controls, thermal_config, cancelled,
        )
    }

    /// The same observed advance with explicit Navier-slip heat routes.
    ///
    /// [`Self::step_with_slip_heating`] owns route admission and physical heat
    /// partition. [`Self::step_with_heating_observed`] owns observer composition
    /// and the shared transaction; both sets of requirements apply.
    ///
    /// # Errors
    /// Invalid routes and all original motion, thermal, observer or cancellation
    /// failures preserve every accepted owner.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_slip_heating_observed<O: DepositionObserver>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut O,
        cells: StokesThermalCells<'_>,
        routes: &[StokesSlipHeatRoute],
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        thermal_config: StepConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<GuidedStokesHeatStep, GuidedStokesObservedHeatError<O::Error>> {
        self.step_heating(
            thermal, observer, cells, Some(routes), timestep_s, force_density, exterior,
            pressure_offsets_pa, external_force_n, controls, thermal_config, cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_heating<O: DepositionObserver>(
        &mut self,
        thermal: &mut DepositionSimulation<'_>,
        observer: &mut O,
        cells: StokesThermalCells<'_>,
        routes: Option<&[StokesSlipHeatRoute]>,
        timestep_s: f64,
        force_density: &dyn Fn([f64; 2]) -> [f64; 2],
        exterior: &[MaterialStokesBoundary<'_>],
        pressure_offsets_pa: &[Option<f64>],
        external_force_n: f64,
        controls: GuidedStokesControls,
        thermal_config: StepConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GuidedStokesHeatStep, GuidedStokesObservedHeatError<O::Error>> {
        poll(&mut cancelled)?;
        let binding = StokesHeatBinding {
            mesh: self.system.mesh,
            density_kg_m3: self.fluid.transport().density_kg_m3(),
            depth_m: self.config.depth_m,
            time_s: self.state.time_s,
        };
        binding.check_thermal_cells(thermal, cells, &mut cancelled)?;
        let routes = match routes {
            Some(routes) => binding.check_slip_heat_routes(
                thermal, cells, &self.edges, exterior, routes, &mut cancelled,
            )?,
            None => BTreeMap::new(),
        };
        let (state, motion) = self.propose_step(
            timestep_s,
            force_density,
            exterior,
            pressure_offsets_pa,
            external_force_n,
            controls,
            &mut cancelled,
        )?;
        let (candidate, accepted_observer, audit, heat_j) = binding.advance_heat_observed(
            thermal,
            observer,
            cells,
            &routes,
            &motion.flow,
            pressure_offsets_pa,
            timestep_s,
            state.time_s,
            motion.bulk_viscous_dissipation_j,
            thermal_config,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        self.state = state;
        *thermal = candidate;
        *observer = accepted_observer;
        Ok(GuidedStokesHeatStep {
            motion,
            thermal: audit,
            deposited_viscous_heat_j: heat_j[0],
            deposited_slip_heat_j: heat_j[1],
        })
    }
}

impl StokesHeatBinding<'_> {
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn advance_heat_observed<'a, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'a>,
        observer: &O,
        cells: StokesThermalCells<'_>,
        routes: &BTreeMap<usize, StokesSlipHeatRoute>,
        flow: &MaterialStokesSolution,
        pressure_offsets_pa: &[Option<f64>],
        timestep_s: f64,
        end_time_s: f64,
        bulk_viscous_heat_j: f64,
        thermal_config: StepConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(DepositionSimulation<'a>, O, EnergyAudit, [f64; 2]),
        ObserverCouplingError<StokesHeatFailure, O::Error>> {
        let mut source = ViscousHeat {
            watts: BTreeMap::new(),
            heat_j: [0.0; 2],
            supplied_heat_j: 0.0,
        };
        let heating = flow.viscous_dissipation().cell_heating_w_m3();
        if heating.len() != cells.triangle_cells.len() {
            return Err(ObserverCouplingError::First(StokesHeatFailure::Binding(
                "incomplete native strain field",
            )));
        }
        for (index, (&id, &rate)) in cells.triangle_cells.iter().zip(heating).enumerate() {
            if index % 256 == 0 {
                heat_poll(cancelled).map_err(ObserverCouplingError::First)?;
            }
            let watts = rate * self.mesh.areas[index] * self.depth_m;
            if !nonnegative(watts) {
                return Err(DepositionError::Arithmetic.into());
            }
            source.watts.insert(id, [watts, 0.0]);
        }
        let mut slip_power_w = 0.0;
        for (index, (&edge, route)) in routes.iter().enumerate() {
            if index % 256 == 0 {
                heat_poll(cancelled).map_err(ObserverCouplingError::First)?;
            }
            let loss = flow
                .boundary_loads(&[edge], [0.0; 2], pressure_offsets_pa)
                .map_err(|error| ObserverCouplingError::First(StokesHeatFailure::Flow(error)))?
                .slip_dissipation_w_per_m;
            let watts = loss * self.depth_m;
            if !nonnegative(watts) || (loss > 0.0 && watts == 0.0) {
                return Err(DepositionError::Arithmetic.into());
            }
            // Compute the smaller share directly, retaining it even when the
            // complementary share rounds to the entire interface power.
            let fraction = route.fluid_heat_fraction;
            let (fluid_w, wall_w) = if fraction <= 0.5 {
                let fluid_w = watts * fraction;
                (fluid_w, watts - fluid_w)
            } else {
                let wall_w = watts * (1.0 - fraction);
                (watts - wall_w, wall_w)
            };
            if watts > 0.0
                && ((fraction > 0.0 && fluid_w == 0.0) || (fraction < 1.0 && wall_w == 0.0))
            {
                return Err(DepositionError::Arithmetic.into());
            }
            let fluid_cell = cells.triangle_cells[self.mesh.edges[edge].tris.0];
            source.add_slip(fluid_cell, fluid_w)?;
            source.add_slip(route.wall_cell, wall_w)?;
            slip_power_w += watts;
        }
        let slip_heat_j = timestep_s * slip_power_w;
        let expected_heat_j = bulk_viscous_heat_j + slip_heat_j;
        if !nonnegative(slip_heat_j)
            || !nonnegative(expected_heat_j)
            || (slip_power_w > 0.0 && slip_heat_j == 0.0)
        {
            return Err(DepositionError::Arithmetic.into());
        }
        let owners = CoupledObservers::new(source, observer.clone());
        let (candidate, owners, audit) =
            thermal.advance_observed(&owners, end_time_s, thermal_config, || {
                if cancelled() {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })?;
        let (source, accepted_observer) = owners.into_parts();
        // This checks conversion and substep integration, not a new fluid
        // energy theorem. Thermal storage has its own original energy gate.
        let residual = [
            source.heat_j[0] - bulk_viscous_heat_j,
            source.heat_j[1] - slip_heat_j,
            source.supplied_heat_j - expected_heat_j,
        ]
        .into_iter()
        .max_by(|a, b| a.abs().total_cmp(&b.abs()))
        .unwrap_or(0.0);
        if candidate.time_s() != end_time_s || residual.abs() > thermal_config.energy_tolerance_j
        {
            return Err(DepositionError::EnergyBalance(residual).into());
        }
        Ok((candidate, accepted_observer, audit, source.heat_j))
    }

    pub(crate) fn check_slip_heat_routes(
        &self,
        thermal: &DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        body_edges: &[usize],
        exterior: &[MaterialStokesBoundary<'_>],
        routes: &[StokesSlipHeatRoute],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<BTreeMap<usize, StokesSlipHeatRoute>, StokesHeatFailure> {
        let mut slip_edges = BTreeSet::new();
        for (index, boundary) in exterior.iter().enumerate() {
            if index % 256 == 0 {
                heat_poll(cancelled)?;
            }
            if let MaterialStokesBoundary::NavierSlip { edge, .. } = boundary
                && !slip_edges.insert(*edge)
            {
                return Err(StokesHeatFailure::Binding(
                    "duplicate Navier-slip boundary",
                ));
            }
        }
        let fluid_cells: BTreeSet<_> = cells.triangle_cells.iter().copied().collect();
        let body_edges: BTreeSet<_> = body_edges.iter().copied().collect();
        let mut admitted = BTreeMap::new();
        for (index, route) in routes.iter().enumerate() {
            if index % 256 == 0 {
                heat_poll(cancelled)?;
            }
            if !slip_edges.remove(&route.edge)
                || route.edge >= self.mesh.edges.len()
                || !self.mesh.is_boundary(route.edge)
                || body_edges.contains(&route.edge)
                || !route.fluid_heat_fraction.is_finite()
                || !(0.0..=1.0).contains(&route.fluid_heat_fraction)
                || !thermal.cells().contains_key(&route.wall_cell)
                || fluid_cells.contains(&route.wall_cell)
            {
                return Err(StokesHeatFailure::Binding(
                    "unique Navier-slip edge, active distinct wall receiver and fraction in [0,1] required",
                ));
            }
            admitted.insert(route.edge, *route);
        }
        if !slip_edges.is_empty() {
            return Err(StokesHeatFailure::Binding(
                "every Navier-slip edge needs a heat route",
            ));
        }
        Ok(admitted)
    }

    pub(crate) fn check_thermal_cells(
        &self,
        thermal: &DepositionSimulation<'_>,
        cells: StokesThermalCells<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), StokesHeatFailure> {
        let mesh = self.mesh;
        let depth = self.depth_m;
        let end_z = cells.plane_z_m + depth;
        if thermal.time_s() != self.time_s
            || cells.triangle_cells.len() != mesh.tris.len()
            || !cells.plane_z_m.is_finite()
            || !end_z.is_finite()
            || !same(end_z - cells.plane_z_m, depth)
        {
            return Err(StokesHeatFailure::Binding(
                "matching clocks, complete triangle map and representable extrusion required",
            ));
        }
        let mut seen = std::collections::BTreeSet::new();
        for (index, &id) in cells.triangle_cells.iter().enumerate() {
            if index % 256 == 0 {
                heat_poll(cancelled)?;
            }
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(StokesHeatFailure::Binding(
                    "inactive thermal triangle receiver",
                ))?;
            if !seen.insert(id) {
                return Err(StokesHeatFailure::Binding(
                    "duplicate thermal triangle receiver",
                ));
            }
            let [x, y] = mesh.centroids[index];
            let start = [x, y, cells.plane_z_m];
            let end = [x, y, end_z];
            let volume = mesh.areas[index] * depth;
            let mass = self.density_kg_m3 * volume;
            let perimeter: f64 = mesh.tri_edges[index]
                .iter()
                .map(|&(edge, _)| mesh.edges[edge].len)
                .sum();
            let surface = 2.0 * mesh.areas[index] + perimeter * depth;
            if !positive(volume)
                || !positive(mass)
                || !positive(surface)
                || !same(cell.geometry.volume_m3, volume)
                || !same(cell.geometry.surface_area_m2, surface)
                || !same(cell.mass_kg, mass)
                || !(0..3).all(|i| {
                    same(cell.geometry.start_m[i], start[i]) && same(cell.geometry.end_m[i], end[i])
                })
            {
                return Err(StokesHeatFailure::Binding(
                    "thermal geometry and mass must match the original fluid triangle prism",
                ));
            }
        }
        Ok(())
    }
}

fn heat_poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), StokesHeatFailure> {
    if cancelled() {
        Err(StokesHeatFailure::Cancelled)
    } else {
        Ok(())
    }
}

fn same(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() <= 64.0 * f64::EPSILON * a.abs().max(b.abs())
}

#[derive(Clone)]
struct ViscousHeat {
    // Separate bulk and interfacial powers, summed only at the thermal source.
    watts: BTreeMap<CellId, [f64; 2]>,
    heat_j: [f64; 2],
    supplied_heat_j: f64,
}

impl ViscousHeat {
    fn add_slip(&mut self, id: CellId, watts: f64) -> Result<(), DepositionError> {
        let entry = self.watts.entry(id).or_default();
        entry[1] += watts;
        if !nonnegative(entry[1]) {
            return Err(DepositionError::Arithmetic);
        }
        Ok(())
    }
}

impl DepositionObserver for ViscousHeat {
    type Error = StokesHeatFailure;

    fn sources(
        &mut self,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let mut sources = BTreeMap::new();
        for (index, (&id, watts)) in self.watts.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let total = watts[0] + watts[1];
            if !nonnegative(total) {
                return Err(DepositionError::Arithmetic.into());
            }
            sources.insert(id, total);
        }
        Ok(sources)
    }

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        // Native events can change boundaries/links and add other cells, but
        // cannot replace an existing receiver's mass, material or geometry.
        Ok(())
    }

    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        for (index, watts) in self.watts.values().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            for (heat, power) in self.heat_j.iter_mut().zip(watts) {
                *heat += dt * power;
            }
            self.supplied_heat_j += dt * (watts[0] + watts[1]);
            if !self.heat_j.iter().all(|heat| heat.is_finite()) || !self.supplied_heat_j.is_finite()
            {
                return Err(DepositionError::Arithmetic.into());
            }
        }
        Ok(())
    }
}
