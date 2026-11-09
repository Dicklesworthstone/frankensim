//! Volumetric shocks with evolving, source-owned chemical composition.
//!
//! Every species density, all three momenta and full referenced total energy
//! share the original mixture Riemann solve. The spatial topology, six-trace
//! limiter, conservative divergence and boundary quadrature are shared with
//! the frozen 3-D field. The default update conserves energy; the explicit
//! contact-equilibrium option accounts for each admitted caloric energy change.
//! [`moving`] advances the same field on translating and piston-driven grids.
//! This Estimated Cartesian lane does not assert a fully discrete
//! entropy/invariant-domain theorem or general contact preservation near shocks.

pub use super::super::mixture::MAX_MIXTURE_WORKSPACE_BYTES;
use super::super::mixture::finite::FiniteStage;
use super::super::{
    cartesian as stepping,
    mixture::{allocated, sum},
};
use super::forcing::{BodyForce3d, ForceChange3d, ForceControl3d};
use super::{
    Boundaries3d, BoundaryCondition, BoundarySide, CartesianGrid, Primitive, cartesian, tangents,
};
pub use super::{MeshMotion, Reconstruction, StepControl};
pub use fs_gas::multicomponent::pressure_equilibrium::{
    ContactEquilibriumControl, ContactEquilibriumStepReport, ContactWorkspace,
};
pub use fs_gas::multicomponent::{Error, MixtureGas};
use fs_gas::multicomponent::{FaceFlux, MAX_SPECIES, SpatialCell, SpatialPrimitive};
use fs_qty::Time;

const WIDTH_LIMIT: usize = MAX_SPECIES + 4;

/// Conservative mechanical forcing within the original accepted-step owner.
pub mod forcing;
/// Reciprocal pressure impulse and work with two finite guided planes.
pub mod inertial_walls;
/// Original species diffusion, Fourier heat and full three-dimensional stress.
pub mod molecular;
/// Conservative original-species transport on moving Cartesian volumes.
pub mod moving;
/// Original finite-rate chemistry coupled atomically to this volumetric field.
pub mod reactive;
/// Native vorticity and production using the evolving composition's calorics.
pub mod vorticity;
pub use forcing::ForcingWorkspace3d;
mod reservoir;
pub use super::super::mixture::{
    FiniteReservoirHeatReport, FiniteReservoirInventory, FiniteReservoirTransfer, WallControl,
    WallThermal,
};
pub use reservoir::{
    FiniteReservoir, FiniteReservoirHeat, FiniteReservoirReport, ReservoirBoundary,
};
/// Prepared physical-face reservoir traces for spatially varying inflow.
pub mod reservoir_profiles;
pub use reservoir_profiles::ReservoirFaceState;
mod reservoir_stages;

/// Borrowed accepted cell in the immutable source's canonical species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cell<'a> {
    /// Species mass densities [kg/m³], including exact zeros.
    pub partial_densities: &'a [f64],
    /// Cartesian x/y/z momentum density [kg/(m² s)].
    pub momentum: [f64; 3],
    /// Full kinetic plus species-referenced thermal energy density [J/m³].
    pub energy: f64,
}
impl<'a> Cell<'a> {
    fn normal(self, axis: usize) -> SpatialCell<'a> {
        let [a, b] = tangents(axis);
        SpatialCell {
            partial_densities: self.partial_densities,
            normal_momentum: self.momentum[axis],
            transverse_momentum: [self.momentum[a], self.momentum[b]],
            energy: self.energy,
        }
    }
    /// Recover the current composition-dependent calorics without a second EOS.
    pub fn primitive(
        self,
        gas: &MixtureGas,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        let p = gas.spatial_primitive(self.normal(0), cancelled)?;
        Ok(Primitive {
            density: p.density,
            velocity: [
                p.velocity,
                self.momentum[1] / p.density,
                self.momentum[2] / p.density,
            ],
            pressure: p.pressure,
            temperature: p.temperature,
            sound_speed: p.sound_speed,
        })
    }
}

/// Construct `[rho_species..., px, py, pz, E_referenced]` through the source.
/// `out` must have exactly `gas.species_count()+4` entries.
pub fn state_into(
    gas: &MixtureGas,
    partial_densities: &[f64],
    velocity: [f64; 3],
    pressure: f64,
    out: &mut [f64],
) -> Result<(), Error> {
    gas.spatial_conserved_into(
        SpatialPrimitive {
            partial_densities,
            normal_velocity: velocity[0],
            transverse_velocity: [velocity[1], velocity[2]],
            pressure,
        },
        out,
    )?;
    Ok(())
}

/// Extensive conserved integrals in the original species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Totals {
    /// Species masses [kg]; unused array entries are exactly zero.
    pub species_masses: [f64; MAX_SPECIES],
    /// Number of original source species.
    pub species_count: usize,
    /// Cartesian momentum [kg m/s].
    pub momentum: [f64; 3],
    /// Species-referenced total energy [J].
    pub energy: f64,
}
impl Totals {
    /// Total mass is derived from the conserved species, never an extra equation.
    pub fn mass(&self) -> f64 {
        sum(self.species_masses[..self.species_count].iter().copied())
    }
    /// Original chemical plus kinetic energy under the source's reference chart.
    pub fn absolute_energy(&self, gas: &MixtureGas) -> Result<f64, Error> {
        if self.species_count != gas.species_count() {
            return Err(Error::Invalid("volumetric total/source species axis"));
        }
        let result = sum(std::iter::once(self.energy).chain(
            self.species_masses[..self.species_count]
                .iter()
                .zip(gas.energy_references())
                .map(|(mass, reference)| mass * reference),
        ));
        if !result.is_finite() {
            return Err(Error::Invalid("volumetric physical energy overflow"));
        }
        Ok(result)
    }
    fn from_row(row: &[f64], n: usize) -> Self {
        let mut species_masses = [0.0; MAX_SPECIES];
        species_masses[..n].copy_from_slice(&row[..n]);
        Self {
            species_masses,
            species_count: n,
            momentum: [row[n], row[n + 1], row[n + 2]],
            energy: row[n + 3],
        }
    }
    fn component(&self, d: usize) -> f64 {
        if d < self.species_count {
            self.species_masses[d]
        } else if d < self.species_count + 3 {
            self.momentum[d - self.species_count]
        } else {
            self.energy
        }
    }
}

/// Outward transfer through one physical face over one accepted RK interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryImpulse<'a> {
    /// Outward species masses [kg].
    pub species_masses: &'a [f64],
    /// Outward momentum [kg m/s], including pressure support impulse.
    pub momentum: [f64; 3],
    /// Outward species-referenced total energy [J].
    pub energy: f64,
}

/// Only a completely admitted step publishes cells, clock and boundary loads.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReport {
    /// Actual accepted duration [s].
    pub elapsed: Time,
    /// Largest combined wave CFL, including prescribed exterior acoustic bounds
    /// in stage-time steps and kicked-cell acoustic bounds in forced steps.
    pub maximum_cfl: f64,
    /// Rejected and halved attempts.
    pub halvings: usize,
    /// Accepted-stage HLL fallback faces.
    pub hll_fallback_faces: usize,
    /// Accepted-stage cells with contracted six-face reconstruction.
    pub limited_cells: usize,
    /// A source refusal switched the complete trial to first order.
    pub first_order_retry: bool,
    /// Initial extensive state.
    pub before: Totals,
    /// Accepted extensive state.
    pub after: Totals,
    /// Net inward boundary transfer integrated with the same RK weights.
    pub boundary_change: Totals,
    /// Applied body impulse and mechanical work, separate from boundary exchange.
    pub source_change: ForceChange3d,
    /// After minus before, boundary change and applied mechanical source.
    pub conservation_residual: Totals,
    /// Residual under the original species-to-physical-energy transform [J].
    pub absolute_energy_residual_j: f64,
}

struct Workspace {
    primitive: Vec<Primitive>,
    slopes: Vec<f64>,
    flux: [Vec<f64>; 3],
    metadata: [Vec<FaceFlux>; 3],
    stage: Vec<f64>,
    candidate: Vec<f64>,
    boundary_trial: Vec<f64>,
}

/// Evolving original species on a physical three-dimensional Cartesian grid.
///
/// Flat prepared buffers avoid heap allocation per cell or accepted step.
/// Density, gas constant, heat capacity and sound speed follow each evolving
/// composition; formation-energy references remain owned by `MixtureGas`.
pub struct MixtureEuler3d {
    gas: MixtureGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: reservoir::Reservoirs,
    states: Vec<f64>,
    time: f64,
    last_step: Option<StepReport>,
    boundary_impulses: Vec<f64>,
    workspace_bytes: usize,
    work: Workspace,
}
impl MixtureEuler3d {
    /// Numerical storage for accepted cells, stages, slopes, faces, cached
    /// primitives and both boundary buffers. Immutable source storage and
    /// bounded stack scratch are separate from this numerical-buffer budget.
    pub fn required_workspace_bytes(
        gas: &MixtureGas,
        grid: &CartesianGrid,
    ) -> Result<usize, Error> {
        let checked = || -> Option<usize> {
            let width = gas.species_count().checked_add(4)?;
            let faces = grid
                .face_len(0)
                .checked_add(grid.face_len(1))?
                .checked_add(grid.face_len(2))?;
            let scalars = grid
                .len()
                .checked_mul(6)?
                .checked_add(faces)?
                .checked_add(grid.boundary_len().checked_mul(2)?)?
                .checked_mul(width)?
                .checked_add(grid.len().checked_mul(3)?)?;
            scalars
                .checked_mul(size_of::<f64>())?
                .checked_add(grid.len().checked_mul(size_of::<Primitive>())?)?
                .checked_add(faces.checked_mul(size_of::<FaceFlux>())?)?
                .checked_add(size_of::<
                    super::super::mixture::finite::FiniteReservoirs<3, 6>,
                >())
        };
        checked().ok_or(Error::Invalid("volumetric species workspace size overflow"))
    }
    /// Admit x-fast flat cell averages and an explicit total numerical budget.
    /// Every row is species, global px/py/pz, then referenced total energy.
    pub fn new(
        gas: MixtureGas,
        grid: CartesianGrid,
        states: Vec<f64>,
        boundaries: Boundaries3d,
        maximum_workspace_bytes: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Self::new_with_reservoirs(
            gas,
            grid,
            states,
            boundaries,
            &[],
            maximum_workspace_bytes,
            cancelled,
        )
    }

    /// Numerical capacity including immutable original-source exterior traces.
    /// Source and side admission remain the constructor's responsibility.
    pub fn required_workspace_bytes_with_reservoirs(
        gas: &MixtureGas,
        grid: &CartesianGrid,
        reservoir_count: usize,
    ) -> Result<usize, Error> {
        Self::required_workspace_bytes(gas, grid)?
            .checked_add(reservoir::Reservoirs::required_bytes(gas, reservoir_count)?)
            .ok_or(Error::Invalid(
                "volumetric reservoir workspace size overflow",
            ))
    }

    /// Supply complete gas states on selected transmissive physical sides.
    ///
    /// Original canonical species densities, Cartesian velocity and absolute
    /// pressure enter the existing spatial Riemann solve. Reservoirs stay fixed
    /// during all stages; their species, momentum and enthalpy/pressure work use
    /// the original conservative boundary account. Unbound sides are unchanged.
    /// This supplies neither finite reservoir storage nor a characteristic-only
    /// boundary. The molecular owner still requires its original closed faces;
    /// an advective trace alone does not define a species/thermal bath.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_reservoirs(
        gas: MixtureGas,
        grid: CartesianGrid,
        states: Vec<f64>,
        boundaries: Boundaries3d,
        reservoirs: &[ReservoirBoundary<'_>],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        for pair in [boundaries.x, boundaries.y, boundaries.z] {
            pair.admit::<Error>()?;
        }
        let bytes = Self::required_workspace_bytes_with_reservoirs(&gas, &grid, reservoirs.len())?
            .checked_add(
                (states.capacity() - states.len())
                    .checked_mul(size_of::<f64>())
                    .ok_or(Error::Invalid("volumetric input capacity overflow"))?,
            )
            .ok_or(Error::Invalid("volumetric species workspace size overflow"))?;
        let n = gas.species_count();
        let width = n + 4;
        let count = grid
            .len()
            .checked_mul(width)
            .ok_or(Error::Invalid("volumetric species state size overflow"))?;
        if bytes > maximum_workspace_bytes
            || bytes > MAX_MIXTURE_WORKSPACE_BYTES
            || states.len() != count
        {
            return Err(Error::Invalid(
                "volumetric species state dimensions or total workspace budget",
            ));
        }
        let mut primitive = allocated(grid.len(), Primitive::default())?;
        for (q, p) in states.chunks_exact(width).zip(&mut primitive) {
            *p = cell(q, n).primitive(&gas, &mut cancelled)?;
        }
        checked_totals(&states, &grid, n, &mut cancelled)?;
        let reservoirs = reservoir::Reservoirs::new(&gas, boundaries, reservoirs, &mut cancelled)?;
        let result = Self {
            gas,
            grid,
            boundaries,
            reservoirs,
            states,
            time: 0.0,
            last_step: None,
            workspace_bytes: bytes,
            boundary_impulses: allocated(grid.boundary_len() * width, 0.0)?,
            work: Workspace {
                primitive,
                // Density and all mass fractions share one contracted slope.
                slopes: allocated(3 * (count + grid.len()), 0.0)?,
                flux: [
                    allocated(grid.face_len(0) * width, 0.0)?,
                    allocated(grid.face_len(1) * width, 0.0)?,
                    allocated(grid.face_len(2) * width, 0.0)?,
                ],
                metadata: [
                    allocated(grid.face_len(0), FaceFlux::default())?,
                    allocated(grid.face_len(1), FaceFlux::default())?,
                    allocated(grid.face_len(2), FaceFlux::default())?,
                ],
                stage: allocated(count, 0.0)?,
                candidate: allocated(count, 0.0)?,
                boundary_trial: allocated(grid.boundary_len() * width, 0.0)?,
            },
        };
        stepping::poll::<Error>(&mut cancelled)?;
        Ok(result)
    }
    /// Select compact, stationary stagnation reservoirs on already bound sides.
    ///
    /// Each selected side must have a uniform original-source reservoir with
    /// exactly zero Cartesian velocity. Its density, pressure and composition
    /// specify the upstream stagnation state. The original one-sided gas law
    /// resolves inflow, outflow, choking and reversal; outward entropy and shear
    /// remain interior-owned. Unselected reservoirs retain their full-state
    /// Riemann boundary. An empty selection restores that original behavior on
    /// held reservoirs; finite tanks always retain their compact pressure law.
    ///
    /// Held reservoirs retain infinite storage. Finite tank bindings preserve
    /// their compact law; neither selection supplies a molecular boundary bath.
    /// Spatial/time-dependent exterior profiles are not selected by this entry.
    /// Selection reuses fixed metadata and allocates no stepping workspace.
    pub fn with_stagnation_reservoirs(
        mut self,
        sides: &[BoundarySide],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        self.reservoirs
            .set_stagnation(&self.gas, sides, &mut cancelled)?;
        Ok(self)
    }

    /// Whether this uniform side uses the compact stagnation-reservoir law.
    pub fn uses_stagnation_reservoir(&self, side: BoundarySide) -> bool {
        self.reservoirs.is_stagnation(side.axis(), side.sign())
    }

    /// Original immutable species, calorics, elements, charge and references.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Physical grid and its original spacings.
    pub const fn grid(&self) -> &CartesianGrid {
        &self.grid
    }
    /// Boundary treatments on all three axes.
    pub const fn boundaries(&self) -> Boundaries3d {
        self.boundaries
    }
    /// Original uniform exterior gas on a selected physical side.
    /// For prepared spatial profiles, use `reservoir_face` instead.
    pub fn reservoir(&self, side: BoundarySide) -> Option<ReservoirBoundary<'_>> {
        self.reservoirs.boundary(side, self.gas.species_count())
    }
    /// Actual admitted numerical storage, including excess input capacity.
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    /// Flat accepted rows in canonical species/global momentum order.
    pub fn states(&self) -> &[f64] {
        &self.states
    }
    /// Number of physical cells.
    pub fn cell_count(&self) -> usize {
        self.grid.len()
    }
    /// Borrow one original-source cell without allocation.
    pub fn cell(&self, index: usize) -> Option<Cell<'_>> {
        let n = self.gas.species_count();
        let width = n + 4;
        (index < self.grid.len()).then(|| cell(&self.states[index * width..(index + 1) * width], n))
    }
    /// Accepted physical time [s].
    pub fn time(&self) -> Time {
        Time::new(self.time)
    }
    /// Latest completely accepted step, absent before the first advance.
    pub const fn last_step(&self) -> Option<&StepReport> {
        self.last_step.as_ref()
    }
    /// Outward face transfers in the grid's original tangential order.
    pub fn boundary_impulses(
        &self,
        side: BoundarySide,
    ) -> impl ExactSizeIterator<Item = BoundaryImpulse<'_>> {
        boundary_view(
            &self.boundary_impulses,
            &self.grid,
            self.gas.species_count(),
            side,
        )
    }
    /// Compensated extensive totals in fixed canonical cell/component order.
    pub fn totals(&self) -> Totals {
        checked_totals(
            &self.states,
            &self.grid,
            self.gas.species_count(),
            &mut || false,
        )
        .expect("admitted finite accepted totals")
    }
    /// Advance at most the requested duration under one combined face CFL.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        self.step_impl(maximum_duration, controls, None, cancelled)
    }

    /// Keep eligible advected material contacts at their admitted pressure
    /// through the original mixture calorics, with explicit numerical energy
    /// accounting. Species, all momenta and physical face fluxes are unchanged.
    ///
    /// Every incident reconstructed pair must match the original cell pressure
    /// and full velocity in both stages. Physical boundary cells and any HLL
    /// fallback remain conservative. Reacting, molecular and forced owners do
    /// not silently inherit this optional correction.
    pub fn step_with_contact_equilibrium(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        contact_control: ContactEquilibriumControl,
        workspace: &mut ContactWorkspace,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ContactEquilibriumStepReport<StepReport>, Error> {
        workspace.admit(self.grid.len(), contact_control)?;
        if self
            .workspace_bytes
            .checked_add(workspace.workspace_bytes())
            .is_none_or(|bytes| bytes > MAX_MIXTURE_WORKSPACE_BYTES)
        {
            return Err(Error::Invalid("volumetric contact workspace budget"));
        }
        let transport = self.step_impl_with_contact(
            maximum_duration,
            controls,
            None,
            Some((contact_control, &mut *workspace)),
            None,
            cancelled,
        )?;
        let correction = workspace.receipt();
        Ok(ContactEquilibriumStepReport {
            corrected_energy_residual_j: transport.absolute_energy_residual_j
                - correction.signed_energy_j,
            transport,
            correction,
        })
    }

    /// Initial pressure-wave bound shared by Euler and complete split owners.
    fn pressure_outlet_rate(&self, cancelled: &mut impl FnMut() -> bool) -> Result<f64, Error> {
        let n = self.gas.species_count();
        let width = n + 4;
        let mut rate: f64 = 0.0;
        let mut outlet_flux = [0.0; WIDTH_LIMIT];
        for side in BoundarySide::ALL {
            let axis = side.axis();
            let sign = side.sign();
            let super::AxisBoundaries::Ends { lower, upper } = self.boundaries.axis(axis) else {
                continue;
            };
            let Some((pressure_pa, reference_velocity)) =
                (if sign < 0.0 { lower } else { upper }).pressure_control()
            else {
                continue;
            };
            let [a, b] = tangents(axis);
            for face in 0..self.grid.face_count(side) {
                stepping::poll::<Error>(cancelled)?;
                let mut coordinate = [0; 3];
                coordinate[axis] = if sign < 0.0 {
                    0
                } else {
                    self.grid.shape[axis] - 1
                };
                coordinate[a] = face % self.grid.shape[a];
                coordinate[b] = face / self.grid.shape[a];
                let index = self.grid.index(coordinate);
                let q = cell(&self.states[index * width..(index + 1) * width], n);
                let flux = if let Some(velocity) = reference_velocity {
                    self.gas.characteristic_outlet_spatial_flux(
                        q.normal(axis),
                        sign,
                        pressure_pa,
                        velocity,
                        &mut outlet_flux[..width],
                        &mut *cancelled,
                    )?
                } else {
                    self.gas.pressure_outlet_spatial_flux(
                        q.normal(axis),
                        sign,
                        pressure_pa,
                        &mut outlet_flux[..width],
                        &mut *cancelled,
                    )?
                };
                let density = sum(q.partial_densities.iter().copied());
                // The characteristic speed includes the sampled face sound
                // speed; shear is carried unchanged through that wave.
                rate = rate.max(
                    flux.signal_speed / self.grid.spacing[axis]
                        + ((q.momentum[a] / density).abs() + flux.signal_speed)
                            / self.grid.spacing[a]
                        + ((q.momentum[b] / density).abs() + flux.signal_speed)
                            / self.grid.spacing[b],
                );
            }
        }
        Ok(rate.max(self.reservoirs.stagnation_rate(
            &self.gas,
            &self.grid,
            &self.states,
            cancelled,
        )?))
    }

    fn step_impl(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        forcing: Option<forcing::Context<'_, '_>>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        self.step_impl_with_contact(maximum_duration, controls, forcing, None, None, cancelled)
    }

    fn step_impl_with_contact(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        mut forcing: Option<forcing::Context<'_, '_>>,
        mut contact: Option<(ContactEquilibriumControl, &mut ContactWorkspace)>,
        mut profiles: Option<&mut reservoir_stages::Context<'_>>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        let n = self.gas.species_count();
        let width = n + 4;
        if let Some(context) = profiles.as_mut() {
            context.prepare_start(
                self.time,
                &self.gas,
                &self.grid,
                &mut self.reservoirs,
                &mut cancelled,
            )?;
        }
        let mut rate = self
            .reservoirs
            .maximum_rate_at(&self.grid, profiles.is_some());
        for (q, p) in self
            .states
            .chunks_exact(width)
            .zip(&mut self.work.primitive)
        {
            *p = cell(q, n).primitive(&self.gas, &mut cancelled)?;
            if let Some((_, control, _)) = forcing.as_ref() {
                control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            }
            rate = rate.max(
                (0..3)
                    .map(|a| (p.velocity[a].abs() + p.sound_speed) / self.grid.spacing[a])
                    .sum(),
            );
        }
        rate = rate.max(self.pressure_outlet_rate(&mut cancelled)?);
        let before = checked_totals(&self.states, &self.grid, n, &mut cancelled)?;
        let mut source_change = ForceChange3d::default();
        let accepted = stepping::advance(
            maximum_duration.value(),
            self.time,
            rate,
            controls,
            &mut cancelled,
            |dt, reconstruction, cancelled| {
                let mut run = || -> Result<stepping::StageStats, reservoir_stages::AttemptError> {
                    if let Some((control, workspace)) = contact.as_mut() {
                        workspace.begin_attempt(*control)?;
                    }
                    if let Some((force, control, workspace)) = forcing.as_mut() {
                        copy_values(&self.states, &mut workspace.kicked, cancelled)?;
                        let mut source = forcing::Accumulator::default();
                        forcing::kick_field(
                            &self.gas,
                            &self.grid,
                            &mut workspace.kicked,
                            *force,
                            *control,
                            0.5 * dt,
                            &mut source,
                            false,
                            cancelled,
                        )?;
                        // The original spatial primitive sweep admits the first
                        // kicked state without a duplicate caloric recovery.
                        let mut stats = if let Some(context) = profiles.as_mut() {
                            context.attempt(
                                self.time,
                                &self.gas,
                                &self.grid,
                                self.boundaries,
                                &mut self.reservoirs,
                                &workspace.kicked,
                                &mut self.work,
                                dt,
                                controls.cfl,
                                reconstruction,
                                Some(*control),
                                cancelled,
                            )?
                        } else {
                            attempt(
                                &self.gas,
                                &self.grid,
                                self.boundaries,
                                &mut self.reservoirs,
                                &workspace.kicked,
                                &mut self.work,
                                dt,
                                controls.cfl,
                                reconstruction,
                                Some(*control),
                                cancelled,
                            )?
                        };
                        let final_rate = forcing::kick_field(
                            &self.gas,
                            &self.grid,
                            &mut self.work.candidate,
                            *force,
                            *control,
                            0.5 * dt,
                            &mut source,
                            true,
                            cancelled,
                        )?;
                        stepping::check_cfl::<Error>(final_rate * dt, controls.cfl)?;
                        stats.rate = stats.rate.max(final_rate);
                        source_change = source.finish()?;
                        Ok(stats)
                    } else if let Some(context) = profiles.as_mut() {
                        context.attempt(
                            self.time,
                            &self.gas,
                            &self.grid,
                            self.boundaries,
                            &mut self.reservoirs,
                            &self.states,
                            &mut self.work,
                            dt,
                            controls.cfl,
                            reconstruction,
                            None,
                            cancelled,
                        )
                    } else {
                        attempt_with_contact(
                            &self.gas,
                            &self.grid,
                            self.boundaries,
                            &mut self.reservoirs,
                            &self.states,
                            &mut self.work,
                            dt,
                            controls.cfl,
                            reconstruction,
                            None,
                            contact.as_mut().map(|(_, workspace)| &mut **workspace),
                            cancelled,
                        )
                        .map_err(reservoir_stages::AttemptError::retryable)
                    }
                };
                run()
            },
        )
        .map_err(|failure| failure.error)?;
        let after = checked_totals(&self.work.candidate, &self.grid, n, &mut cancelled)?;
        let outward = sum_rows(&self.work.boundary_trial, width, 1.0, &mut cancelled)?;
        let boundary_change = Totals::from_row(&outward.map(|v| -v), n);
        let residual = std::array::from_fn::<_, WIDTH_LIMIT, _>(|d| {
            if d < width {
                let delta = after.component(d) - before.component(d) - boundary_change.component(d);
                if forcing.is_some() && d >= n {
                    delta
                        - if d < n + 3 {
                            source_change.impulse[d - n]
                        } else {
                            source_change.work
                        }
                } else {
                    delta
                }
            } else {
                0.0
            }
        });
        if !residual.iter().all(|v| v.is_finite()) {
            return Err(Error::Invalid("volumetric species balance overflow"));
        }
        let conservation_residual = Totals::from_row(&residual[..width], n);
        let report = StepReport {
            elapsed: Time::new(accepted.dt),
            maximum_cfl: accepted.stats.rate * accepted.dt,
            halvings: accepted.halvings,
            hll_fallback_faces: accepted.stats.fallback,
            limited_cells: accepted.stats.limited,
            first_order_retry: accepted.reconstruction != controls.reconstruction,
            before,
            after,
            boundary_change,
            source_change,
            absolute_energy_residual_j: conservation_residual.absolute_energy(&self.gas)?,
            conservation_residual,
        };
        if let Some((_, workspace)) = contact.as_ref() {
            if !(report.absolute_energy_residual_j - workspace.receipt().signed_energy_j)
                .is_finite()
            {
                return Err(Error::Invalid(
                    "volumetric corrected energy balance overflow",
                ));
            }
        }
        stepping::poll::<Error>(&mut cancelled)?;
        if profiles.is_some() {
            self.reservoirs
                .profiles
                .as_mut()
                .expect("prepared stage profiles")
                .commit();
        }
        self.reservoirs.finite.publish();
        std::mem::swap(&mut self.states, &mut self.work.candidate);
        self.boundary_impulses
            .copy_from_slice(&self.work.boundary_trial);
        self.time += accepted.dt;
        self.last_step = Some(report);
        Ok(report)
    }
}

#[allow(clippy::too_many_arguments)]
fn attempt(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &mut reservoir::Reservoirs,
    original: &[f64],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    force_control: Option<ForceControl3d>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    attempt_with_contact(
        gas,
        grid,
        boundaries,
        reservoirs,
        original,
        work,
        dt,
        cfl,
        reconstruction,
        force_control,
        None,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn attempt_with_contact(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &mut reservoir::Reservoirs,
    original: &[f64],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    force_control: Option<ForceControl3d>,
    mut contact: Option<&mut ContactWorkspace>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    stepping::rk2(cancelled, |final_stage, cancelled| {
        rk_stage(
            gas,
            grid,
            boundaries,
            reservoirs,
            false,
            original,
            work,
            dt,
            cfl,
            reconstruction,
            force_control,
            contact.as_deref_mut(),
            final_stage,
            cancelled,
        )
    })
}

/// One original SSP stage; only the borrowed exterior profile buffer varies.
#[allow(clippy::too_many_arguments)]
fn rk_stage(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &mut reservoir::Reservoirs,
    trial_profiles: bool,
    original: &[f64],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    force_control: Option<ForceControl3d>,
    mut contact: Option<&mut ContactWorkspace>,
    final_stage: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    if !final_stage {
        reservoirs.finite.begin_heat(gas, dt, cancelled)?;
    }
    let mut stats = spatial_with_contact(
        gas,
        grid,
        boundaries,
        reservoirs,
        trial_profiles,
        if final_stage {
            FiniteStage::Trial
        } else {
            FiniteStage::TransportStart
        },
        if final_stage { &work.stage } else { original },
        &mut work.primitive,
        &mut work.slopes,
        &mut work.flux,
        &mut work.metadata,
        reconstruction,
        force_control,
        contact
            .as_deref_mut()
            .map(|workspace| (workspace, !final_stage)),
        cancelled,
    )?;
    stepping::check_cfl::<Error>(stats.rate * dt, cfl)?;
    let stage_rate = if final_stage {
        update_with_contact(
            gas,
            grid,
            original,
            &work.stage,
            &work.flux,
            &mut work.candidate,
            dt,
            true,
            force_control,
            contact.as_deref_mut(),
            cancelled,
        )?
    } else {
        update_with_contact(
            gas,
            grid,
            original,
            original,
            &work.flux,
            &mut work.stage,
            dt,
            false,
            force_control,
            contact.as_deref_mut(),
            cancelled,
        )?
    };
    if force_control.is_some() {
        stepping::check_cfl::<Error>(stage_rate * dt, cfl)?;
        stats.rate = stats.rate.max(stage_rate);
    }
    let width = gas.species_count() + 4;
    cartesian::capture::<Error>(
        grid,
        width,
        0.5 * dt,
        final_stage,
        |axis, face, d| work.flux[axis][face * width + d],
        |boundary, d, delta, add| {
            let out = &mut work.boundary_trial[boundary * width + d];
            *out = (if add { *out } else { 0.0 }) + delta;
        },
        cancelled,
    )?;
    reservoirs.finite.stage(
        gas,
        &work.boundary_trial,
        BoundarySide::ALL.map(|side| grid.face_count(side)),
        if final_stage { 1.0 } else { 2.0 },
        cancelled,
    )?;
    if final_stage {
        reservoirs.finite.finish_heat(gas, dt, cancelled)?;
    }
    Ok(stats)
}

#[allow(clippy::too_many_arguments)]
fn spatial(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &reservoir::Reservoirs,
    states: &[f64],
    primitive: &mut [Primitive],
    slopes: &mut [f64],
    fluxes: &mut [Vec<f64>; 3],
    metadata: &mut [Vec<FaceFlux>; 3],
    reconstruction: Reconstruction,
    force_control: Option<ForceControl3d>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    spatial_with_contact(
        gas,
        grid,
        boundaries,
        reservoirs,
        false,
        FiniteStage::Accepted,
        states,
        primitive,
        slopes,
        fluxes,
        metadata,
        reconstruction,
        force_control,
        None,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn spatial_with_contact(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &reservoir::Reservoirs,
    trial_profiles: bool,
    finite_stage: FiniteStage,
    states: &[f64],
    primitive: &mut [Primitive],
    slopes: &mut [f64],
    fluxes: &mut [Vec<f64>; 3],
    metadata: &mut [Vec<FaceFlux>; 3],
    reconstruction: Reconstruction,
    force_control: Option<ForceControl3d>,
    contact: Option<(&mut ContactWorkspace, bool)>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    spatial_with_motion(
        gas,
        grid,
        boundaries,
        reservoirs,
        trial_profiles,
        finite_stage,
        states,
        primitive,
        slopes,
        fluxes,
        metadata,
        reconstruction,
        force_control,
        contact,
        None,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn spatial_with_motion(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    reservoirs: &reservoir::Reservoirs,
    trial_profiles: bool,
    finite_stage: FiniteStage,
    states: &[f64],
    primitive: &mut [Primitive],
    slopes: &mut [f64],
    fluxes: &mut [Vec<f64>; 3],
    metadata: &mut [Vec<FaceFlux>; 3],
    reconstruction: Reconstruction,
    force_control: Option<ForceControl3d>,
    mut contact: Option<(&mut ContactWorkspace, bool)>,
    motion: Option<MeshMotion>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    let n = gas.species_count();
    let width = n + 4;
    let slope_width = width + 1;
    let mut cell_rate: f64 = 0.0;
    for (k, (q, p)) in states
        .chunks_exact(width)
        .zip(primitive.iter_mut())
        .enumerate()
    {
        *p = cell(q, n).primitive(gas, &mut *cancelled)?;
        if let Some((workspace, first_stage)) = contact.as_mut() {
            workspace.set_reference(k, p.pressure, p.velocity, p.sound_speed, *first_stage)?;
        }
        if let Some(control) = force_control {
            control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
            let local = if let Some(motion) = motion {
                forcing::moving_acoustic_rate(*p, grid, k, motion)?
            } else {
                forcing::acoustic_rate(*p, grid)?
            };
            cell_rate = cell_rate.max(local);
        }
    }
    let value = |k: usize, d: usize| {
        if d == 0 {
            primitive[k].density
        } else if d <= n {
            states[k * width + d - 1] / primitive[k].density
        } else if d < n + 4 {
            primitive[k].velocity[d - n - 1]
        } else {
            primitive[k].pressure
        }
    };
    let reconstruct = |k: usize, sign: f64, slope: &[f64], out: &mut [f64]| {
        let density = value(k, 0) + sign * slope[0];
        for d in 0..n {
            out[d] = density * (value(k, d + 1) + sign * slope[d + 1]);
        }
        for d in n..width {
            out[d] = value(k, d + 1) + sign * slope[d + 1];
        }
    };
    let mut trace = [0.0; WIDTH_LIMIT];
    let mut check = [0.0; WIDTH_LIMIT];
    let wall_velocity = motion.map_or([[0.0; 2]; 3], |m| {
        std::array::from_fn(|axis| [m.lower_velocity[axis], m.upper_velocity[axis]])
    });
    let limited = cartesian::slopes_with_wall_velocity::<Error>(
        grid,
        boundaries,
        slope_width,
        [n + 1, n + 2, n + 3],
        reconstruction,
        wall_velocity,
        value,
        |k, lo, hi, slope| {
            slope[0] = stepping::mc(value(k, 0) - value(lo, 0), value(hi, 0) - value(k, 0));
            fs_gas::multicomponent::reconstruction::mass_fraction_slopes(
                &states[lo * width..lo * width + n],
                &states[k * width..k * width + n],
                &states[hi * width..hi * width + n],
                &mut slope[1..=n],
            )?;
            Ok(n + 1)
        },
        slopes,
        |k, axis, sign, slope| {
            reconstruct(k, sign, slope, &mut trace[..width]);
            gas.spatial_conserved_into(
                normal_primitive(&trace[..width], n, axis),
                &mut check[..width],
            )
            .is_ok()
        },
        cancelled,
    )?;
    let mut stats = stepping::StageStats {
        limited,
        ..stepping::StageStats::default()
    };
    let mut finite_rates = [0.0; 6];
    let mut finite_corrections = [0.0; 6];
    let mut left = [0.0; WIDTH_LIMIT];
    let mut right = [0.0; WIDTH_LIMIT];
    let mut result = [0.0; WIDTH_LIMIT];
    cartesian::faces::<Error, _>(grid, boundaries, cancelled, |face, cancelled| {
        let axis = face.axis;
        let flux = &mut fluxes[axis];
        let metadata = &mut metadata[axis];
        if let cartesian::FaceKind::PeriodicCopy { first } = face.kind {
            for d in 0..width {
                flux[face.index * width + d] = flux[first * width + d];
            }
            metadata[face.index] = metadata[first];
            return Ok(());
        }
        let face_trace = |k: usize, sign: f64, out: &mut [f64]| {
            let offset = (3 * k + axis) * slope_width;
            reconstruct(k, sign, &slopes[offset..offset + slope_width], out);
        };
        let q = |k: usize| cell(&states[k * width..(k + 1) * width], n).normal(axis);
        let first_order = reconstruction == Reconstruction::FirstOrder;
        let receipt = if let Some(motion) = motion {
            let speed = motion.face_velocity(grid, axis, face.index);
            // Use the source's primitive constructor for the physical face
            // sound speed. No caloric inversion or new face storage is needed.
            let mut admit_trace = |values: &[f64]| -> Result<(), Error> {
                if let Some(control) = force_control {
                    let p = gas.spatial_conserved_into(
                        normal_primitive(values, n, axis),
                        &mut check[..width],
                    )?;
                    control.admit_velocity::<Error>(
                        [values[n], values[n + 1], values[n + 2]],
                        p.sound_speed,
                    )?;
                }
                Ok(())
            };
            match face.kind {
                cartesian::FaceKind::Interior { left: l, right: r } => {
                    if first_order {
                        gas.moving_hllc_spatial_flux(
                            q(l),
                            q(r),
                            speed,
                            &mut result[..width],
                            &mut *cancelled,
                        )?
                    } else {
                        face_trace(l, 0.5, &mut left[..width]);
                        face_trace(r, -0.5, &mut right[..width]);
                        admit_trace(&left[..width])?;
                        admit_trace(&right[..width])?;
                        gas.moving_hllc_spatial_flux_from_primitive(
                            normal_primitive(&left[..width], n, axis),
                            normal_primitive(&right[..width], n, axis),
                            speed,
                            &mut result[..width],
                            &mut *cancelled,
                        )?
                    }
                }
                cartesian::FaceKind::Boundary {
                    cell: k,
                    sign,
                    condition,
                } => {
                    if !first_order {
                        face_trace(k, 0.5 * sign, &mut left[..width]);
                        admit_trace(&left[..width])?;
                    }
                    match (first_order, condition) {
                        (true, BoundaryCondition::Reflecting) => gas.moving_wall_spatial_flux(
                            q(k),
                            sign,
                            speed,
                            &mut result[..width],
                            &mut *cancelled,
                        )?,
                        (true, BoundaryCondition::Transmissive) => gas.moving_hllc_spatial_flux(
                            q(k),
                            q(k),
                            speed,
                            &mut result[..width],
                            &mut *cancelled,
                        )?,
                        (false, BoundaryCondition::Reflecting) => gas
                            .moving_wall_spatial_flux_from_primitive(
                                normal_primitive(&left[..width], n, axis),
                                sign,
                                speed,
                                &mut result[..width],
                                &mut *cancelled,
                            )?,
                        (false, BoundaryCondition::Transmissive) => {
                            let p = normal_primitive(&left[..width], n, axis);
                            gas.moving_hllc_spatial_flux_from_primitive(
                                p,
                                p,
                                speed,
                                &mut result[..width],
                                &mut *cancelled,
                            )?
                        }
                        _ => unreachable!("moving mixture topology admitted before the trial"),
                    }
                }
                cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
            }
        } else {
            match face.kind {
                cartesian::FaceKind::Interior { left: l, right: r } => {
                    if first_order {
                        gas.hllc_spatial_flux(q(l), q(r), &mut result[..width], &mut *cancelled)?
                    } else {
                        face_trace(l, 0.5, &mut left[..width]);
                        face_trace(r, -0.5, &mut right[..width]);
                        gas.hllc_spatial_flux_from_primitive(
                            normal_primitive(&left[..width], n, axis),
                            normal_primitive(&right[..width], n, axis),
                            &mut result[..width],
                            &mut *cancelled,
                        )?
                    }
                }
                cartesian::FaceKind::Boundary {
                    cell: k,
                    sign,
                    condition,
                } => {
                    let c = grid.coordinates(k);
                    let [a, b] = tangents(axis);
                    let boundary_face = c[a] + grid.shape[a] * c[b];
                    if let Some(reservoir) =
                        reservoirs.stagnation_cell_at(axis, sign, n, finite_stage)
                    {
                        if first_order {
                            gas.reservoir_spatial_flux(
                                q(k),
                                reservoir,
                                sign,
                                &mut result[..width],
                                &mut *cancelled,
                            )?
                        } else {
                            face_trace(k, 0.5 * sign, &mut left[..width]);
                            gas.reservoir_spatial_flux_from_primitive(
                                normal_primitive(&left[..width], n, axis),
                                reservoir,
                                sign,
                                &mut result[..width],
                                &mut *cancelled,
                            )?
                        }
                    } else if let Some(exterior) =
                        reservoirs.normal_face_at(axis, sign, boundary_face, n, trial_profiles)
                    {
                        let interior = if first_order {
                            let [a, b] = tangents(axis);
                            SpatialPrimitive {
                                partial_densities: &states[k * width..k * width + n],
                                normal_velocity: primitive[k].velocity[axis],
                                transverse_velocity: [
                                    primitive[k].velocity[a],
                                    primitive[k].velocity[b],
                                ],
                                pressure: primitive[k].pressure,
                            }
                        } else {
                            face_trace(k, 0.5 * sign, &mut left[..width]);
                            normal_primitive(&left[..width], n, axis)
                        };
                        let (left, right) = if sign < 0.0 {
                            (exterior, interior)
                        } else {
                            (interior, exterior)
                        };
                        gas.hllc_spatial_flux_from_primitive(
                            left,
                            right,
                            &mut result[..width],
                            &mut *cancelled,
                        )?
                    } else {
                        match (first_order, condition) {
                            (true, BoundaryCondition::Reflecting) => gas
                                .stationary_wall_spatial_flux(
                                    q(k),
                                    sign,
                                    &mut result[..width],
                                    &mut *cancelled,
                                )?,
                            (true, BoundaryCondition::Transmissive) => gas.hllc_spatial_flux(
                                q(k),
                                q(k),
                                &mut result[..width],
                                &mut *cancelled,
                            )?,
                            (true, BoundaryCondition::PressureOutlet { pressure_pa }) => gas
                                .pressure_outlet_spatial_flux(
                                    q(k),
                                    sign,
                                    pressure_pa,
                                    &mut result[..width],
                                    &mut *cancelled,
                                )?,
                            (
                                true,
                                BoundaryCondition::CharacteristicOutlet {
                                    reference_pressure_pa,
                                    reference_velocity_m_s,
                                },
                            ) => gas.characteristic_outlet_spatial_flux(
                                q(k),
                                sign,
                                reference_pressure_pa,
                                reference_velocity_m_s,
                                &mut result[..width],
                                &mut *cancelled,
                            )?,
                            (false, BoundaryCondition::Reflecting) => {
                                face_trace(k, 0.5 * sign, &mut left[..width]);
                                gas.stationary_wall_spatial_flux_from_primitive(
                                    normal_primitive(&left[..width], n, axis),
                                    sign,
                                    &mut result[..width],
                                    &mut *cancelled,
                                )?
                            }
                            (false, BoundaryCondition::Transmissive) => {
                                face_trace(k, 0.5 * sign, &mut left[..width]);
                                let p = normal_primitive(&left[..width], n, axis);
                                gas.hllc_spatial_flux_from_primitive(
                                    p,
                                    p,
                                    &mut result[..width],
                                    &mut *cancelled,
                                )?
                            }
                            (false, BoundaryCondition::PressureOutlet { pressure_pa }) => {
                                face_trace(k, 0.5 * sign, &mut left[..width]);
                                gas.pressure_outlet_spatial_flux_from_primitive(
                                    normal_primitive(&left[..width], n, axis),
                                    sign,
                                    pressure_pa,
                                    &mut result[..width],
                                    &mut *cancelled,
                                )?
                            }
                            (
                                false,
                                BoundaryCondition::CharacteristicOutlet {
                                    reference_pressure_pa,
                                    reference_velocity_m_s,
                                },
                            ) => {
                                face_trace(k, 0.5 * sign, &mut left[..width]);
                                gas.characteristic_outlet_spatial_flux_from_primitive(
                                    normal_primitive(&left[..width], n, axis),
                                    sign,
                                    reference_pressure_pa,
                                    reference_velocity_m_s,
                                    &mut result[..width],
                                    &mut *cancelled,
                                )?
                            }
                        }
                    }
                }
                cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
            }
        };
        if let cartesian::FaceKind::Boundary { sign, .. } = face.kind {
            let slot = 2 * axis + usize::from(sign > 0.0);
            if let Some(volume) = reservoirs.finite.volume(slot) {
                let increment =
                    0.5 * grid.face_area(BoundarySide::ALL[slot]) * receipt.signal_speed / volume;
                let y = increment - finite_corrections[slot];
                let next = finite_rates[slot] + y;
                finite_corrections[slot] = (next - finite_rates[slot]) - y;
                finite_rates[slot] = next;
                if !finite_rates[slot].is_finite() {
                    return Err(Error::Invalid("finite volumetric reservoir CFL overflow"));
                }
            }
        }
        if let Some((workspace, _)) = contact.as_mut() {
            match face.kind {
                cartesian::FaceKind::Interior { left: l, right: r } => {
                    let (pl, vl, pr, vr) = if first_order {
                        (
                            primitive[l].pressure,
                            primitive[l].velocity,
                            primitive[r].pressure,
                            primitive[r].velocity,
                        )
                    } else {
                        (
                            left[n + 3],
                            [left[n], left[n + 1], left[n + 2]],
                            right[n + 3],
                            [right[n], right[n + 1], right[n + 2]],
                        )
                    };
                    workspace.observe_pair(l, pl, vl, pr, vr, receipt.used_hll_fallback)?;
                    workspace.observe_pair(r, pl, vl, pr, vr, receipt.used_hll_fallback)?;
                }
                cartesian::FaceKind::Boundary { cell, .. } => workspace.exclude(cell)?,
                cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
            }
        }
        stats.fallback += usize::from(receipt.used_hll_fallback);
        let momentum = [result[n], result[n + 1], result[n + 2]];
        let [a, b] = tangents(axis);
        result[n + axis] = momentum[0];
        result[n + a] = momentum[1];
        result[n + b] = momentum[2];
        flux[face.index * width..(face.index + 1) * width].copy_from_slice(&result[..width]);
        metadata[face.index] = receipt;
        Ok(())
    })?;
    stats.rate = cartesian::rate::<Error>(
        grid,
        |axis, face| metadata[axis][face].signal_speed,
        cancelled,
    )?
    .max(finite_rates.into_iter().fold(0.0, f64::max));
    if force_control.is_some() {
        stats.rate = stats.rate.max(cell_rate);
    }
    Ok(stats)
}

#[allow(clippy::too_many_arguments)]
fn update(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    original: &[f64],
    stage: &[f64],
    fluxes: &[Vec<f64>; 3],
    output: &mut [f64],
    dt: f64,
    final_stage: bool,
    force_control: Option<ForceControl3d>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let n = gas.species_count();
    let width = n + 4;
    let mut rate: f64 = 0.0;
    cartesian::update::<WIDTH_LIMIT, Error, _>(
        grid,
        width,
        dt,
        final_stage,
        |k, d| original[k * width + d],
        |k, d| stage[k * width + d],
        |axis, face, d| fluxes[axis][face * width + d],
        cancelled,
        |k, row, cancelled| {
            let p = cell(row, n).primitive(gas, &mut *cancelled)?;
            if let Some(control) = force_control {
                control.admit_velocity::<Error>(p.velocity, p.sound_speed)?;
                rate = rate.max(forcing::acoustic_rate(p, grid)?);
            }
            output[k * width..(k + 1) * width].copy_from_slice(row);
            Ok(())
        },
    )?;
    Ok(rate)
}

#[allow(clippy::too_many_arguments)]
fn update_with_contact(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    original: &[f64],
    stage: &[f64],
    fluxes: &[Vec<f64>; 3],
    output: &mut [f64],
    dt: f64,
    final_stage: bool,
    force_control: Option<ForceControl3d>,
    contact: Option<&mut ContactWorkspace>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, Error> {
    let Some(workspace) = contact else {
        return update(
            gas,
            grid,
            original,
            stage,
            fluxes,
            output,
            dt,
            final_stage,
            force_control,
            cancelled,
        );
    };
    let n = gas.species_count();
    let width = n + 4;
    cartesian::update_energy_projected::<WIDTH_LIMIT, Error, _>(
        grid,
        width,
        dt,
        final_stage,
        |k, d| original[k * width + d],
        |k, d| stage[k * width + d],
        |axis, face, d| fluxes[axis][face * width + d],
        cancelled,
        |k, row, weight, cancelled| {
            workspace.project_energy(
                gas,
                k,
                &row[..n],
                &row[n..n + 3],
                row[n + 3],
                grid.cell_volume(),
                weight,
                cancelled,
            )
        },
        |k, row, cancelled| {
            cell(row, n).primitive(gas, &mut *cancelled)?;
            output[k * width..(k + 1) * width].copy_from_slice(row);
            Ok(())
        },
    )?;
    Ok(0.0)
}

fn normal_primitive(q: &[f64], n: usize, axis: usize) -> SpatialPrimitive<'_> {
    let [a, b] = tangents(axis);
    SpatialPrimitive {
        partial_densities: &q[..n],
        normal_velocity: q[n + axis],
        transverse_velocity: [q[n + a], q[n + b]],
        pressure: q[n + 3],
    }
}
fn cell(q: &[f64], n: usize) -> Cell<'_> {
    Cell {
        partial_densities: &q[..n],
        momentum: [q[n], q[n + 1], q[n + 2]],
        energy: q[n + 3],
    }
}
fn boundary_view<'a>(
    values: &'a [f64],
    grid: &CartesianGrid,
    species: usize,
    side: BoundarySide,
) -> impl ExactSizeIterator<Item = BoundaryImpulse<'a>> {
    let width = species + 4;
    let range = grid.boundary_range(side);
    values[range.start * width..range.end * width]
        .chunks_exact(width)
        .map(move |q| BoundaryImpulse {
            species_masses: &q[..species],
            momentum: [q[species], q[species + 1], q[species + 2]],
            energy: q[species + 3],
        })
}
fn copy_values(
    source: &[f64],
    target: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    for (source, target) in source.chunks(4096).zip(target.chunks_mut(4096)) {
        stepping::poll::<Error>(cancelled)?;
        target.copy_from_slice(source);
    }
    Ok(())
}
fn sum_rows(
    values: &[f64],
    width: usize,
    scale: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[f64; WIDTH_LIMIT], Error> {
    let mut sums = [0.0; WIDTH_LIMIT];
    let mut correction = [0.0; WIDTH_LIMIT];
    for row in values.chunks_exact(width) {
        stepping::poll::<Error>(cancelled)?;
        for d in 0..width {
            let y = row[d] * scale - correction[d];
            let next = sums[d] + y;
            correction[d] = (next - sums[d]) - y;
            sums[d] = next;
        }
    }
    if !sums.iter().all(|v| v.is_finite()) {
        return Err(Error::Invalid("volumetric species extensive overflow"));
    }
    Ok(sums)
}
fn checked_totals(
    states: &[f64],
    grid: &CartesianGrid,
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Totals, Error> {
    Ok(Totals::from_row(
        &sum_rows(states, n + 4, grid.cell_volume(), cancelled)?,
        n,
    ))
}

#[cfg(test)]
mod tests;
