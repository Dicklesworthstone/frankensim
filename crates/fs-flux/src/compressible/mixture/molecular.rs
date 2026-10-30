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
pub(in crate::compressible) mod interval;
pub use fs_gas::multicomponent::dissipation::{
    DissipationControl, DissipationError, PlanarTransport, TransportLaw,
};
use fs_gas::multicomponent::{MixtureGas, NormalPrimitive};
use fs_qty::Time;
use interval::{Action, Stage, parabolic};

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
    /// Accepted SSP RK2 substeps.
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
    x_flux: Vec<f64>,
    y_flux: Vec<f64>,
    boundary_step: Vec<f64>,
    boundary_trial: Vec<f64>,
    boundary_correction: Vec<f64>,
}

/// Prepared evolving-species transport with complete planar Newtonian stress.
///
/// Periodic faces share identical flux bits. Reflecting faces are impermeable,
/// adiabatic stationary slip walls: normal viscous traction loads their support,
/// while species, tangent traction, heat and wall work vanish. Open faces need
/// the field's explicit original-source reservoir used by its HLLC boundary.
/// Species, enthalpy, heat and full viscous work use one two-trace molecular
/// flux at half-cell normal distance. Constitutive values use the shared
/// midpoint approximation; boundary velocity and referenced species enthalpy
/// supply the original bath energy effort. Refine the boundary layer independently.
/// Unbound Transmissive faces still refuse.
pub struct MolecularTransport2d<L: TransportLaw> {
    source: PlanarTransport<L>,
    grid: CartesianGrid,
    boundaries: Boundaries2d,
    maximum_workspace_bytes: usize,
    boundary_impulses: Vec<f64>,
    last_step: Option<MolecularReport>,
    work: Workspace,
}
impl<L: TransportLaw> MolecularTransport2d<L> {
    /// Required combined numerical buffers for the supplied flow and operator.
    /// Includes source-kernel/MS work, every field/stage/face/gradient buffer and
    /// incoming state excess capacity. Immutable sources, prepared coefficients,
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
                .checked_add(
                    flow.grid
                        .len()
                        .checked_mul(size_of::<Primitive>() + size_of::<[f64; 4]>())?,
                )?
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
            maximum_workspace_bytes,
            boundary_impulses: allocated(boundary, 0.0)?,
            last_step: None,
            work: Workspace {
                current: allocated(count, 0.0)?,
                stage: allocated(count, 0.0)?,
                candidate: allocated(count, 0.0)?,
                primitive: allocated(grid.len(), Primitive::default())?,
                gradient: allocated(grid.len(), [0.0; 4])?,
                x_flux: allocated((nx + 1) * ny * width, 0.0)?,
                y_flux: allocated(nx * (ny + 1) * width, 0.0)?,
                boundary_step: allocated(boundary, 0.0)?,
                boundary_trial: allocated(boundary, 0.0)?,
                boundary_correction: allocated(boundary, 0.0)?,
            },
        };
        poll(&mut cancelled)?;
        Ok(result)
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
            + (self.boundary_impulses.capacity()
                + self.work.current.capacity()
                + self.work.stage.capacity()
                + self.work.candidate.capacity()
                + self.work.x_flux.capacity()
                + self.work.y_flux.capacity()
                + self.work.boundary_step.capacity()
                + self.work.boundary_trial.capacity()
                + self.work.boundary_correction.capacity())
                * size_of::<f64>()
            + self.work.primitive.capacity() * size_of::<Primitive>()
            + self.work.gradient.capacity() * size_of::<[f64; 4]>()
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
        let interval = interval::advance(
            duration.value(),
            controls,
            cancelled,
            |action, cancelled| {
                if matches!(action, Action::First) {
                    return spatial(
                        &mut self.source,
                        &self.grid,
                        self.boundaries,
                        &flow.reservoirs,
                        &self.work.current,
                        &mut self.work.primitive,
                        &mut self.work.gradient,
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
                        &flow.reservoirs,
                        &self.work.current,
                        &mut self.work.primitive,
                        &mut self.work.gradient,
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
                    &flow.reservoirs,
                    &self.work.stage,
                    &mut self.work.primitive,
                    &mut self.work.gradient,
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
    reservoirs: &super::reservoir::Reservoirs,
    states: &[f64],
    primitive: &mut [Primitive],
    gradient: &mut [[f64; 4]],
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
    for (index, row) in gradient.iter_mut().enumerate() {
        poll(cancelled)?;
        for axis in 0..2 {
            for component in 0..2 {
                let low = reservoirs
                    .velocity_neighbour(grid, boundaries, primitive, index, axis, component, false);
                let high = reservoirs
                    .velocity_neighbour(grid, boundaries, primitive, index, axis, component, true);
                row[2 * axis + component] = (high - low) / (2.0 * grid.spacing[axis]);
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
        let receipt = match face.kind {
            cartesian::FaceKind::Interior { left, right } => {
                let tangential = std::array::from_fn(|component| {
                    let global = if component == 0 { axis } else { tangent };
                    0.5 * gradient[left][2 * tangent + global]
                        + 0.5 * gradient[right][2 * tangent + global]
                });
                source.face_from_primitives(
                    normal(left, axis),
                    normal(right, axis),
                    grid.spacing[axis],
                    tangential,
                    controls.diffusion,
                    output,
                    &mut *cancelled,
                )?
            }
            cartesian::FaceKind::Boundary {
                cell,
                sign,
                condition: BoundaryCondition::Reflecting,
            } => source.wall_from_primitive(
                normal(cell, axis),
                -sign * 2.0 * primitive[cell].velocity[axis] / grid.spacing[axis],
                gradient[cell][2 * tangent + tangent],
                output,
                &mut *cancelled,
            )?,
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
            cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
        };
        if axis == 1 {
            output.swap(species, species + 1);
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
