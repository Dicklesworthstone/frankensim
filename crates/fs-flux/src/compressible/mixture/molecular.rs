//! Molecular species, heat and momentum transport on the original Cartesian field.
//!
//! Spatial topology, conservative divergence and RK boundary integration belong
//! to the existing Cartesian owner. The original source kernel owns all mole/
//! mass conversion, Maxwell--Stefan diffusion, partial enthalpy, Fourier heat,
//! Newtonian stress and constitutive admission. This module composes those
//! owners with prepared flat buffers and atomic requested-duration publication.

use super::{
    Boundaries2d, BoundaryCondition, BoundaryImpulse, BoundarySide, CartesianGrid,
    MAX_MIXTURE_WORKSPACE_BYTES, MixtureEuler2d, Primitive, Totals, WIDTH_LIMIT, allocated,
    boundary_view, cartesian, cell, checked_totals, copy_values, sum_rows, update,
};
use crate::compressible::AxisBoundaries;
pub use crate::compressible::molecular::MolecularWall2d;
use crate::compressible::stress::{self, Viscosity};
pub(in crate::compressible) mod adaptive;
pub(in crate::compressible) mod interval;
pub(in crate::compressible) mod stabilized;
pub use adaptive::{AdaptiveDissipationControl, AdaptiveMolecularReport};
pub use fs_gas::multicomponent::dissipation::{
    DissipationControl, DissipationError, PlanarTransport, TransportLaw, WallThermalCondition,
};
use fs_gas::multicomponent::{MixtureGas, NormalPrimitive};
use fs_qty::Time;
use interval::{Action, Stage, parabolic};
pub use stabilized::{StabilizedDissipationControl, StabilizedMolecularReport};

mod coupled;
pub use coupled::{
    MolecularEuler2d, MolecularStepControl, MolecularStepReport, StabilizedMolecularStepControl,
    StabilizedMolecularStepReport,
};

// The shared face/capture loops require only cancellation and invalid-input
// cases. Remaining trait methods retain the original molecular refusal types.
impl cartesian::Failure for DissipationError {
    fn invalid(reason: &'static str) -> Self {
        Self::Invalid(reason)
    }
    fn cancelled() -> Self {
        Self::Cancelled
    }
    fn cfl(observed: f64, maximum: f64) -> Self {
        Self::ParabolicNumber { observed, maximum }
    }
    fn retries(_attempts: usize, last: Self) -> Self {
        Self::StepFailure {
            elapsed_time_s: 0.0,
            step_s: 0.0,
            last: Box::new(last),
        }
    }
    fn is_cancelled(&self) -> bool {
        matches!(self, Self::Cancelled)
    }
    fn is_cfl(&self) -> bool {
        matches!(self, Self::ParabolicNumber { .. })
    }
}

/// Complete molecular update; the physical clock remains with the composing owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularReport {
    /// Complete represented physical duration; the Euler clock is unchanged.
    pub elapsed: Time,
    /// Accepted SSP RK2 substeps or opt-in RKL2 macrosteps.
    pub accepted_steps: usize,
    /// Failed internal trials before eventual successful substeps.
    pub rejected_steps: usize,
    /// Largest source/stencil parabolic estimate over accepted stages.
    pub maximum_parabolic_number: f64,
    /// Largest original Maxwell--Stefan row residual [1/m].
    pub maximum_gradient_residual_per_m: f64,
    /// Largest absolute net species diffusion mass flux [kg/(m² s)].
    pub maximum_mass_frame_residual_kg_m2_s: f64,
    /// Complete initial extensive state.
    pub before: Totals,
    /// Complete accepted extensive state.
    pub after: Totals,
    /// Inward exchange integrated over all accepted RK stages and substeps.
    pub boundary_change: Totals,
    /// `after-before-boundary_change` in original canonical component order.
    pub conservation_residual: Totals,
    /// Original species-reference transform of the conservative residual [J].
    pub absolute_energy_residual_j: f64,
}

struct Workspace {
    current: Vec<f64>,
    stage: Vec<f64>,
    candidate: Vec<f64>,
    primitive: Vec<Primitive>,
    /// Row-major derivative axes, then velocity component: du_i/dx_j.
    gradient: Vec<[f64; 4]>,
    viscosity: Vec<Viscosity>,
    x_flux: Vec<f64>,
    y_flux: Vec<f64>,
    // Retained lazily by RKL2, never allocated per stage or retry.
    stabilized_x_flux: Vec<f64>,
    stabilized_y_flux: Vec<f64>,
    stabilized_boundary: Vec<f64>,
    adaptive: adaptive::Workspace,
    boundary_step: Vec<f64>,
    boundary_trial: Vec<f64>,
    boundary_correction: Vec<f64>,
}

/// Prepared evolving-species transport with complete planar Newtonian stress.
///
/// Periodic faces share identical flux bits. Reflecting faces default to
/// impermeable, adiabatic free-slip walls. A side may instead prescribe no-slip
/// or finite Navier slip, tangential wall velocity and conductive thermal data.
/// Finite slip explicitly partitions interfacial friction heat between the
/// fluid and wall; all species remain impermeable. Open faces need
/// the field's explicit original-source reservoir used by its HLLC boundary.
/// Species, enthalpy, heat and full viscous work use one two-trace molecular
/// flux at half-cell normal distance. Constitutive values use the shared
/// midpoint approximation; boundary velocity and referenced species enthalpy
/// supply the original bath energy effort. Refine the boundary layer independently.
/// Unbound Transmissive faces still refuse.
///
/// On fully periodic or adiabatic free-slip domains, source viscosities are
/// sampled at cells and combined with the compatible full-tensor stress stencil.
/// Its molecular contribution dissipates kinetic energy semidiscretely even
/// for varying positive shear/bulk viscosity. Species and heat keep the original
/// face law. Any bath, no-slip or finite-slip side retains the original full
/// mechanical operator; no unsupported boundary-work or finite-step claim is made.
pub struct MolecularTransport2d<L: TransportLaw> {
    source: PlanarTransport<L>,
    grid: CartesianGrid,
    boundaries: Boundaries2d,
    walls: [MolecularWall2d; 4],
    maximum_workspace_bytes: usize,
    boundary_impulses: Vec<f64>,
    last_step: Option<MolecularReport>,
    work: Workspace,
}
impl<L: TransportLaw> MolecularTransport2d<L> {
    /// Required combined numerical buffers for the supplied flow and operator.
    /// Includes source-kernel/MS work, every field/stage/face/gradient buffer and
    /// incoming state excess capacity, including sixteen viscosity bytes per cell.
    /// Immutable sources and the provider's own prepared coefficients,
    /// bounded stack scratch and arbitrary TransportLaw storage/allocations are
    /// excluded: this is not a whole-process or provider-allocation budget.
    pub fn required_workspace_bytes(
        flow: &MixtureEuler2d,
        source: &PlanarTransport<L>,
    ) -> Result<usize, DissipationError> {
        let [nx, ny] = flow.grid.shape;
        let width = flow.gas.species_count() + 3;
        let checked = || -> Option<usize> {
            let faces = (nx + 1)
                .checked_mul(ny)?
                .checked_add(nx.checked_mul(ny + 1)?)?;
            let boundary = nx.checked_add(ny)?.checked_mul(2)?;
            let scalars = flow
                .grid
                .len()
                .checked_mul(3)?
                .checked_add(faces)?
                .checked_add(boundary.checked_mul(4)?)?
                .checked_mul(width)?;
            let bytes = scalars
                .checked_mul(size_of::<f64>())?
                .checked_add(flow.grid.len().checked_mul(
                    size_of::<Primitive>() + size_of::<[f64; 4]>() + size_of::<Viscosity>(),
                )?)?
                .checked_add(source.workspace_bytes())?
                .checked_add(flow.workspace_bytes())?;
            Some(bytes)
        };
        checked().ok_or(DissipationError::Invalid(
            "planar molecular workspace size overflow",
        ))
    }

    /// Bind the exact source/domain, Cartesian geometry and admitted boundaries.
    /// The byte limit covers this operator together with the supplied flow.
    pub fn new(
        flow: &MixtureEuler2d,
        source: PlanarTransport<L>,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if source.gas() != flow.gas() {
            return Err(DissipationError::Invalid(
                "planar molecular source or caloric path differs",
            ));
        }
        flow.reservoirs.admit_molecular(flow.boundaries)?;
        let bytes = Self::required_workspace_bytes(flow, &source)?;
        if bytes > maximum_workspace_bytes || bytes > MAX_MIXTURE_WORKSPACE_BYTES {
            return Err(DissipationError::Invalid(
                "planar molecular total workspace budget",
            ));
        }
        let grid = flow.grid;
        let [nx, ny] = grid.shape;
        let width = flow.gas.species_count() + 3;
        let count = grid.len() * width;
        let boundary = 2 * (nx + ny) * width;
        let result = Self {
            source,
            grid,
            boundaries: flow.boundaries,
            walls: [MolecularWall2d::SlipAdiabatic; 4],
            maximum_workspace_bytes,
            boundary_impulses: allocated(boundary, 0.0)?,
            last_step: None,
            work: Workspace {
                current: allocated(count, 0.0)?,
                stage: allocated(count, 0.0)?,
                candidate: allocated(count, 0.0)?,
                primitive: allocated(grid.len(), Primitive::default())?,
                gradient: allocated(grid.len(), [0.0; 4])?,
                viscosity: allocated(grid.len(), Viscosity::default())?,
                x_flux: allocated((nx + 1) * ny * width, 0.0)?,
                y_flux: allocated(nx * (ny + 1) * width, 0.0)?,
                stabilized_x_flux: Vec::new(),
                stabilized_y_flux: Vec::new(),
                stabilized_boundary: Vec::new(),
                adaptive: adaptive::Workspace::default(),
                boundary_step: allocated(boundary, 0.0)?,
                boundary_trial: allocated(boundary, 0.0)?,
                boundary_correction: allocated(boundary, 0.0)?,
            },
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Configure one physical reflecting side using the shared molecular wall
    /// law. Periodic and open sides have no such wall. Normal wall motion needs
    /// a separate conservative moving-geometry owner and is refused here.
    pub fn with_wall(
        mut self,
        side: BoundarySide,
        wall: MolecularWall2d,
    ) -> Result<Self, DissipationError> {
        let AxisBoundaries::Ends { lower, upper } = self.boundaries.axis(side.axis()) else {
            return Err(DissipationError::Invalid(
                "periodic molecular side has no physical wall",
            ));
        };
        let condition = if side.sign() < 0.0 { lower } else { upper };
        if condition != BoundaryCondition::Reflecting {
            return Err(DissipationError::Invalid(
                "molecular wall model requires a reflecting Euler side",
            ));
        }
        wall.validate(
            side.axis(),
            self.gas().domain().minimum_temperature,
            self.gas().domain().maximum_temperature,
        )
        .map_err(DissipationError::Invalid)?;
        self.walls[wall_index(side.axis(), side.sign())] = wall;
        Ok(self)
    }

    /// Configured molecular model, including inert placeholders on periodic sides.
    pub fn wall(&self, side: BoundarySide) -> MolecularWall2d {
        self.walls[wall_index(side.axis(), side.sign())]
    }

    /// Whether every physical side admits compatible molecular stress.
    /// This is a semidiscrete molecular identity, not a time-step certificate.
    pub fn uses_compatible_molecular_stress(&self) -> bool {
        compatible_boundaries(self.boundaries, &self.walls)
    }

    /// Exact immutable source species and caloric path.
    pub fn gas(&self) -> &MixtureGas {
        self.source.gas()
    }
    /// Original transport law and its source/material identities.
    pub fn law(&self) -> &L {
        self.source.law()
    }
    /// Actual owned numerical capacities; excludes the external flow and the
    /// immutable/provider allocations described by required_workspace_bytes.
    pub fn workspace_bytes(&self) -> usize {
        self.source.workspace_bytes()
            + self.work.adaptive.workspace_bytes()
            + (self.boundary_impulses.capacity()
                + self.work.current.capacity()
                + self.work.stage.capacity()
                + self.work.candidate.capacity()
                + self.work.x_flux.capacity()
                + self.work.y_flux.capacity()
                + self.work.stabilized_x_flux.capacity()
                + self.work.stabilized_y_flux.capacity()
                + self.work.stabilized_boundary.capacity()
                + self.work.boundary_step.capacity()
                + self.work.boundary_trial.capacity()
                + self.work.boundary_correction.capacity())
                * size_of::<f64>()
            + self.work.primitive.capacity() * size_of::<Primitive>()
            + self.work.gradient.capacity() * size_of::<[f64; 4]>()
            + self.work.viscosity.capacity() * size_of::<Viscosity>()
    }
    /// Last completely published requested-duration molecular interval.
    pub const fn last_step(&self) -> Option<&MolecularReport> {
        self.last_step.as_ref()
    }
    /// Outward molecular exchange on each face of the last complete interval.
    pub fn boundary_impulses(
        &self,
        side: BoundarySide,
    ) -> impl ExactSizeIterator<Item = BoundaryImpulse<'_>> {
        boundary_view(
            &self.boundary_impulses,
            &self.grid,
            self.gas().species_count(),
            side,
        )
    }

    /// Advance the entire requested molecular duration atomically. Original
    /// Euler time, advection report and advection impulses remain unchanged.
    pub fn advance(
        &mut self,
        flow: &mut MixtureEuler2d,
        duration: Time,
        controls: DissipationControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        let report = self.advance_unpublished(
            flow,
            duration,
            controls,
            self.maximum_workspace_bytes,
            &mut cancelled,
        )?;
        self.publish(report);
        Ok(report)
    }

    pub(super) fn trial_impulses(&self) -> &[f64] {
        &self.work.boundary_trial
    }
    pub(super) fn publish(&mut self, report: MolecularReport) {
        std::mem::swap(&mut self.boundary_impulses, &mut self.work.boundary_trial);
        self.last_step = Some(report);
    }
    pub(super) fn advance_unpublished(
        &mut self,
        flow: &mut MixtureEuler2d,
        duration: Time,
        controls: DissipationControl,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        poll(cancelled)?;
        controls.validate()?;
        let before = self.prepare_interval(flow, duration, maximum_workspace_bytes, cancelled)?;
        let species = self.gas().species_count();
        let interval = interval::advance(
            duration.value(),
            controls,
            cancelled,
            |action, _elapsed, cancelled| {
                if matches!(action, Action::First) {
                    return spatial(
                        &mut self.source,
                        &self.grid,
                        self.boundaries,
                        &self.walls,
                        &flow.reservoirs,
                        &self.work.current,
                        &mut self.work.primitive,
                        &mut self.work.gradient,
                        &mut self.work.viscosity,
                        &mut self.work.x_flux,
                        &mut self.work.y_flux,
                        controls,
                        cancelled,
                    );
                }
                let Action::Trial { dt, first, cached } = action else {
                    unreachable!()
                };
                let first = if cached {
                    first
                } else {
                    spatial(
                        &mut self.source,
                        &self.grid,
                        self.boundaries,
                        &self.walls,
                        &flow.reservoirs,
                        &self.work.current,
                        &mut self.work.primitive,
                        &mut self.work.gradient,
                        &mut self.work.viscosity,
                        &mut self.work.x_flux,
                        &mut self.work.y_flux,
                        controls,
                        cancelled,
                    )?
                };
                parabolic(first.rate * dt, controls.parabolic_cfl)?;
                update(
                    self.source.gas(),
                    &self.grid,
                    &self.work.current,
                    &self.work.current,
                    &self.work.x_flux,
                    &self.work.y_flux,
                    &mut self.work.stage,
                    dt,
                    false,
                    None,
                    cancelled,
                )?;
                capture(
                    &self.grid,
                    species,
                    &self.work.x_flux,
                    &self.work.y_flux,
                    0.5 * dt,
                    false,
                    &mut self.work.boundary_step,
                    cancelled,
                )?;
                let second = spatial(
                    &mut self.source,
                    &self.grid,
                    self.boundaries,
                    &self.walls,
                    &flow.reservoirs,
                    &self.work.stage,
                    &mut self.work.primitive,
                    &mut self.work.gradient,
                    &mut self.work.viscosity,
                    &mut self.work.x_flux,
                    &mut self.work.y_flux,
                    controls,
                    cancelled,
                )?;
                parabolic(second.rate * dt, controls.parabolic_cfl)?;
                update(
                    self.source.gas(),
                    &self.grid,
                    &self.work.current,
                    &self.work.stage,
                    &self.work.x_flux,
                    &self.work.y_flux,
                    &mut self.work.candidate,
                    dt,
                    true,
                    None,
                    cancelled,
                )?;
                capture(
                    &self.grid,
                    species,
                    &self.work.x_flux,
                    &self.work.y_flux,
                    0.5 * dt,
                    true,
                    &mut self.work.boundary_step,
                    cancelled,
                )?;
                for (index, ((total, correction), change)) in self
                    .work
                    .boundary_trial
                    .iter_mut()
                    .zip(&mut self.work.boundary_correction)
                    .zip(&self.work.boundary_step)
                    .enumerate()
                {
                    if index % 4096 == 0 {
                        poll(cancelled)?;
                    }
                    let adjusted = change - *correction;
                    let next = *total + adjusted;
                    *correction = (next - *total) - adjusted;
                    *total = next;
                }
                std::mem::swap(&mut self.work.current, &mut self.work.candidate);
                Ok(first.merge(second))
            },
        )?;
        self.complete_interval(flow, duration, before, interval, cancelled)
    }

    // Both time integrators retain the original source/geometry admission,
    // requested-duration state gathering and reset of unpublished transfers.
    fn prepare_interval(
        &mut self,
        flow: &MixtureEuler2d,
        duration: Time,
        maximum_workspace_bytes: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Totals, DissipationError> {
        if self.gas() != flow.gas() || self.grid != flow.grid || self.boundaries != flow.boundaries
        {
            return Err(DissipationError::Invalid(
                "planar molecular field ownership differs",
            ));
        }
        flow.reservoirs.admit_molecular(flow.boundaries)?;
        let bytes = self
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "planar molecular workspace size overflow",
            ))?;
        if bytes > maximum_workspace_bytes || bytes > MAX_MIXTURE_WORKSPACE_BYTES {
            return Err(DissipationError::Invalid(
                "planar molecular total workspace budget",
            ));
        }
        if !duration.value().is_finite() || duration.value() < 0.0 {
            return Err(DissipationError::Invalid(
                "nonnegative finite molecular duration required",
            ));
        }
        let species = self.gas().species_count();
        let before = checked_totals(&flow.states, &self.grid, species, cancelled)?;
        copy_values(&flow.states, &mut self.work.current, cancelled)?;
        for values in [
            &mut self.work.boundary_trial,
            &mut self.work.boundary_correction,
        ] {
            for tile in values.chunks_mut(4096) {
                poll(cancelled)?;
                tile.fill(0.0);
            }
        }
        Ok(before)
    }

    fn complete_interval(
        &mut self,
        flow: &mut MixtureEuler2d,
        duration: Time,
        before: Totals,
        interval: interval::Interval,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MolecularReport, DissipationError> {
        let species = self.gas().species_count();
        let after = checked_totals(&self.work.current, &self.grid, species, cancelled)?;
        let outward = sum_rows(&self.work.boundary_trial, species + 3, 1.0, cancelled)?;
        let boundary_change = Totals::from_row(&outward.map(|value| -value), species);
        let residual = std::array::from_fn::<_, WIDTH_LIMIT, _>(|d| {
            if d < species + 3 {
                after.component(d) - before.component(d) - boundary_change.component(d)
            } else {
                0.0
            }
        });
        if !residual.iter().all(|value| value.is_finite()) {
            return Err(DissipationError::Invalid(
                "planar molecular balance overflow",
            ));
        }
        let conservation_residual = Totals::from_row(&residual[..species + 3], species);
        let report = MolecularReport {
            elapsed: duration,
            accepted_steps: interval.accepted_steps,
            rejected_steps: interval.rejected_steps,
            maximum_parabolic_number: interval.maximum_parabolic_number,
            maximum_gradient_residual_per_m: interval.maximum_gradient_residual_per_m,
            maximum_mass_frame_residual_kg_m2_s: interval.maximum_mass_frame_residual_kg_m2_s,
            before,
            after,
            boundary_change,
            conservation_residual,
            absolute_energy_residual_j: conservation_residual.absolute_energy(self.gas())?,
        };
        // Publication is an O(1) exchange, so cancellation never exposes a
        // partially copied field. Excess input capacity follows its real owner.
        let bytes = flow.workspace_bytes - flow.states.capacity() * size_of::<f64>()
            + self.work.current.capacity() * size_of::<f64>();
        poll(cancelled)?;
        std::mem::swap(&mut flow.states, &mut self.work.current);
        flow.workspace_bytes = bytes;
        Ok(report)
    }
}

#[allow(clippy::too_many_arguments)]
fn spatial<L: TransportLaw>(
    source: &mut PlanarTransport<L>,
    grid: &CartesianGrid,
    boundaries: Boundaries2d,
    walls: &[MolecularWall2d; 4],
    reservoirs: &super::reservoir::Reservoirs,
    states: &[f64],
    primitive: &mut [Primitive],
    gradient: &mut [[f64; 4]],
    viscosity: &mut [Viscosity],
    x_flux: &mut [f64],
    y_flux: &mut [f64],
    controls: DissipationControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, DissipationError> {
    let species = source.gas().species_count();
    let width = species + 3;
    for (row, value) in states.chunks_exact(width).zip(primitive.iter_mut()) {
        *value = cell(row, species).primitive(source.gas(), &mut *cancelled)?;
    }
    let compatible = compatible_boundaries(boundaries, walls);
    let mut mechanical_scale = 0.0_f64;
    if compatible {
        let mut minimum_density = f64::INFINITY;
        for (index, coefficients) in viscosity.iter_mut().enumerate() {
            let p = primitive[index];
            let sampled = source.coefficients_from_primitive(
                NormalPrimitive {
                    partial_densities: &states[index * width..index * width + species],
                    normal_velocity: p.velocity[0],
                    transverse_velocity: p.velocity[1],
                    pressure: p.pressure,
                },
                &mut *cancelled,
            )?;
            *coefficients = Viscosity {
                shear: sampled.dynamic_viscosity.value(),
                bulk: sampled.bulk_viscosity.value(),
            };
            mechanical_scale = mechanical_scale.max(coefficients.scale::<2>());
            minimum_density = minimum_density.min(p.density);
        }
        // A neighboring high viscosity acts on both cells, including the less
        // dense one. Do not bound only each coefficient divided by its own rho.
        mechanical_scale /= minimum_density;
        if !mechanical_scale.is_finite() || mechanical_scale < 0.0 {
            return Err(DissipationError::Invalid(
                "planar compatible molecular scale overflow",
            ));
        }
    }
    for (index, row) in gradient.iter_mut().enumerate() {
        poll(cancelled)?;
        for axis in 0..2 {
            let low = reservoirs
                .velocity_neighbour(grid, boundaries, walls, primitive, index, axis, false)?;
            let high = reservoirs
                .velocity_neighbour(grid, boundaries, walls, primitive, index, axis, true)?;
            for component in 0..2 {
                row[2 * axis + component] =
                    (high[component] - low[component]) / (2.0 * grid.spacing[axis]);
            }
        }
    }
    let normal = |index: usize, axis: usize| NormalPrimitive {
        partial_densities: &states[index * width..index * width + species],
        normal_velocity: primitive[index].velocity[axis],
        transverse_velocity: primitive[index].velocity[1 - axis],
        pressure: primitive[index].pressure,
    };
    let mut stats = Stage::default();
    cartesian::faces::<DissipationError, _>(grid, boundaries, cancelled, |face, cancelled| {
        let axis = face.axis;
        let tangent = 1 - axis;
        let flux = if axis == 0 {
            &mut *x_flux
        } else {
            &mut *y_flux
        };
        if let cartesian::FaceKind::PeriodicCopy { first } = face.kind {
            flux.copy_within(first * width..(first + 1) * width, face.index * width);
            return Ok(());
        }
        let output = &mut flux[face.index * width..(face.index + 1) * width];
        let mut receipt = match face.kind {
            cartesian::FaceKind::Interior { left, right } => {
                let tangential = std::array::from_fn(|component| {
                    let global = if component == 0 { axis } else { tangent };
                    0.5 * gradient[left][2 * tangent + global]
                        + 0.5 * gradient[right][2 * tangent + global]
                });
                let receipt = source.face_from_primitives(
                    normal(left, axis),
                    normal(right, axis),
                    grid.spacing[axis],
                    tangential,
                    controls.diffusion,
                    output,
                    &mut *cancelled,
                )?;
                let cells = [viscosity[left], viscosity[right]];
                let original = Viscosity {
                    shear: receipt.dynamic_viscosity_pa_s,
                    bulk: receipt.bulk_viscosity_pa_s,
                };
                if compatible && (cells[0] != original || cells[1] != original) {
                    let traction = stress::traction(
                        axis,
                        grid.spacing[axis],
                        [primitive[left].velocity, primitive[right].velocity],
                        cells,
                        |side, component, derivative| {
                            gradient[[left, right][side]][2 * derivative + component]
                        },
                    );
                    // Source output is still in its normal/tangential frame;
                    // the common axis permutation below restores global order.
                    output[species] = -traction[axis];
                    output[species + 1] = -traction[tangent];
                    let mut energy = receipt.nonmechanical_energy_flux_w_m2;
                    for component in 0..2 {
                        let velocity = 0.5 * primitive[left].velocity[component]
                            + 0.5 * primitive[right].velocity[component];
                        energy -= velocity * traction[component];
                    }
                    output[species + 2] = energy;
                    if !output.iter().all(|value| value.is_finite()) {
                        return Err(DissipationError::Invalid(
                            "planar compatible molecular flux overflow",
                        ));
                    }
                }
                receipt
            }
            cartesian::FaceKind::Boundary {
                cell,
                sign,
                condition: BoundaryCondition::Reflecting,
            } => {
                let wall = walls[wall_index(axis, sign)];
                match wall {
                    MolecularWall2d::SlipAdiabatic => source.wall_from_primitive(
                        normal(cell, axis),
                        -sign * 2.0 * primitive[cell].velocity[axis] / grid.spacing[axis],
                        gradient[cell][2 * tangent + tangent],
                        output,
                        &mut *cancelled,
                    )?,
                    MolecularWall2d::NoSlip { velocity, thermal }
                    | MolecularWall2d::NavierSlip {
                        velocity,
                        thermal,
                        slip_length_m: 0.0,
                        ..
                    } => source.no_slip_wall_from_primitive(
                        normal(cell, axis),
                        0.5 * grid.spacing[axis],
                        sign,
                        velocity[tangent],
                        thermal,
                        -sign * 2.0 * primitive[cell].velocity[axis] / grid.spacing[axis],
                        0.0,
                        output,
                        &mut *cancelled,
                    )?,
                    MolecularWall2d::NavierSlip { thermal, .. } => {
                        let distance = 0.5 * grid.spacing[axis];
                        let p = primitive[cell];
                        let ghost = wall
                            .velocity_ghost(axis, distance, p.velocity)
                            .ok_or(DissipationError::Invalid("missing Navier wall fluid trace"))?;
                        let trace = 0.5 * p.velocity[tangent] + 0.5 * ghost[tangent];
                        // The original source owns mole composition, wall/cell
                        // temperature, impermeable species, full normal stress
                        // and conduction. Query it once at the actual trace.
                        let receipt = source.no_slip_wall_from_primitive(
                            normal(cell, axis),
                            distance,
                            sign,
                            trace,
                            thermal,
                            -sign * 2.0 * p.velocity[axis] / grid.spacing[axis],
                            wall.tangential_trace_fraction(distance)
                                * gradient[cell][2 * tangent + tangent],
                            output,
                            &mut *cancelled,
                        )?;
                        // Keep that source's conductive heat exactly. The
                        // shared finite-slip kernel supplies only mechanical
                        // work and the authored interfacial heat partition.
                        let (shear, mechanical_energy) = wall.shear_and_energy(
                            axis,
                            sign,
                            p.velocity,
                            p.temperature,
                            distance,
                            receipt.dynamic_viscosity_pa_s,
                            0.0,
                        );
                        output[species + 1] = -shear[tangent];
                        output[species + 2] = mechanical_energy + receipt.reduced_heat_flux_w_m2;
                        if !output.iter().all(|value| value.is_finite()) {
                            return Err(DissipationError::Invalid(
                                "Navier molecular wall flux overflow",
                            ));
                        }
                        receipt
                    }
                }
            }
            cartesian::FaceKind::Boundary {
                cell,
                sign,
                condition: BoundaryCondition::Transmissive,
            } => {
                let exterior = reservoirs
                    .normal(axis, sign, species)
                    .ok_or(DissipationError::Invalid("unbound molecular reservoir"))?;
                let interior = normal(cell, axis);
                let (left, right) = if sign < 0.0 {
                    (exterior, interior)
                } else {
                    (interior, exterior)
                };
                // Uniform held bath has zero tangential derivatives. The same
                // source kernel uses midpoint constitutive values and exact
                // admitted bath enthalpies/velocity for boundary energy effort.
                let tangential = [
                    0.5 * gradient[cell][2 * tangent + axis],
                    0.5 * gradient[cell][2 * tangent + tangent],
                ];
                source.reservoir_face_from_primitives(
                    left,
                    right,
                    0.5 * grid.spacing[axis],
                    tangential,
                    sign < 0.0,
                    controls.diffusion,
                    output,
                    &mut *cancelled,
                )?
            }
            cartesian::FaceKind::Boundary {
                condition: BoundaryCondition::PressureOutlet { .. },
                ..
            } => {
                return Err(DissipationError::Invalid(
                    "pressure alone does not define a species/thermal molecular reservoir",
                ));
            }
            cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
        };
        if axis == 1 {
            output.swap(species, species + 1);
        }
        if compatible {
            receipt.diffusivity_m2_s = receipt.nonmechanical_diffusivity_m2_s.max(mechanical_scale);
        }
        stats.add(receipt);
        Ok(())
    })?;
    // Covers both axes, full normal/cross-stress stencils and doubled wall
    // gradients. Source admission additionally checks both SSP states. This
    // is an engineering step-size estimate, not an invariant-domain theorem.
    stats.rate *= 4.0 * (grid.spacing[0].recip().powi(2) + grid.spacing[1].recip().powi(2));
    if !stats.rate.is_finite() || stats.rate < 0.0 {
        return Err(DissipationError::Invalid(
            "planar molecular parabolic scale overflow",
        ));
    }
    Ok(stats)
}

#[allow(clippy::too_many_arguments)]
fn capture(
    grid: &CartesianGrid,
    species: usize,
    x_flux: &[f64],
    y_flux: &[f64],
    duration: f64,
    add: bool,
    output: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    let width = species + 3;
    cartesian::capture::<DissipationError>(
        grid,
        width,
        duration,
        add,
        |axis, face, d| {
            if axis == 0 {
                x_flux[face * width + d]
            } else {
                y_flux[face * width + d]
            }
        },
        |face, d, value, add| {
            let out = &mut output[face * width + d];
            *out = if add { *out + value } else { value };
        },
        cancelled,
    )
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipationError> {
    if cancelled() {
        Err(DissipationError::Cancelled)
    } else {
        Ok(())
    }
}

fn wall_index(axis: usize, sign: f64) -> usize {
    2 * axis + usize::from(sign > 0.0)
}

fn compatible_boundaries(boundaries: Boundaries2d, walls: &[MolecularWall2d; 4]) -> bool {
    (0..2).all(|axis| match boundaries.axis(axis) {
        AxisBoundaries::Periodic => true,
        AxisBoundaries::Ends { lower, upper } => {
            [lower, upper]
                .into_iter()
                .enumerate()
                .all(|(side, condition)| {
                    condition == BoundaryCondition::Reflecting
                        && walls[2 * axis + side] == MolecularWall2d::SlipAdiabatic
                })
        }
    })
}
