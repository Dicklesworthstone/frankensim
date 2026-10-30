//! Original species diffusion, Fourier heat and complete volumetric stress.
//!
//! The shared three-dimensional grid owns the conservative face divergence
//! and boundary quadrature. The original gas source owns Maxwell--Stefan,
//! partial enthalpies, Fourier heat and Newtonian stress/work. Prepared field,
//! gradient and face buffers compose these owners without a second transport
//! law or an additional viscous-heating source.

use super::{
    Boundaries3d, BoundaryCondition, BoundaryImpulse, BoundarySide, CartesianGrid,
    MAX_MIXTURE_WORKSPACE_BYTES, MixtureEuler3d, Primitive, Totals, WIDTH_LIMIT, allocated,
    boundary_view, cartesian, cell, checked_totals, copy_values, stepping, sum_rows, tangents,
    update,
};
use crate::compressible::AxisBoundaries;
use fs_gas::multicomponent::dissipation::TransportFaceReport;
pub use fs_gas::multicomponent::dissipation::{
    DissipationControl, DissipationError, SpatialTransport, TransportLaw,
};
use fs_gas::multicomponent::{MixtureGas, SpatialPrimitive};
use fs_qty::Time;

/// Complete molecular interval; physical time remains with its composing owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MolecularReport {
    /// Complete represented physical duration; the Euler clock is unchanged.
    pub elapsed: Time,
    /// Accepted SSP RK2 substeps.
    pub accepted_steps: usize,
    /// Rejected internal trials before successful substeps.
    pub rejected_steps: usize,
    /// Largest source/stencil parabolic estimate over accepted stages.
    pub maximum_parabolic_number: f64,
    /// Largest original Maxwell--Stefan row residual [1/m].
    pub maximum_gradient_residual_per_m: f64,
    /// Largest absolute net diffusive species flux [kg/(m² s)].
    pub maximum_mass_frame_residual_kg_m2_s: f64,
    /// Complete initial extensive state.
    pub before: Totals,
    /// Complete accepted extensive state.
    pub after: Totals,
    /// Inward exchange integrated over all accepted stages and substeps.
    pub boundary_change: Totals,
    /// `after-before-boundary_change` in original component order.
    pub conservation_residual: Totals,
    /// Original reference transform of the conservative energy residual [J].
    pub absolute_energy_residual_j: f64,
}

#[derive(Clone, Copy, Default)]
struct Stage {
    rate: f64,
    gradient_residual: f64,
    mass_frame_residual: f64,
}
impl Stage {
    fn add(&mut self, face: TransportFaceReport) {
        self.rate = self.rate.max(face.diffusivity_m2_s);
        self.gradient_residual = self
            .gradient_residual
            .max(face.maximum_gradient_residual_per_m);
        self.mass_frame_residual = self
            .mass_frame_residual
            .max(face.mass_frame_residual_kg_m2_s.abs());
    }
}

struct Workspace {
    current: Vec<f64>,
    stage: Vec<f64>,
    candidate: Vec<f64>,
    primitive: Vec<Primitive>,
    /// Row-major derivative axis, then velocity component: du_i/dx_j.
    gradient: Vec<[[f64; 3]; 3]>,
    flux: [Vec<f64>; 3],
    boundary_step: Vec<f64>,
    boundary_trial: Vec<f64>,
    boundary_correction: Vec<f64>,
}

/// Prepared molecular transport on the original evolving-species 3-D field.
///
/// One stored face flux serves both neighbors, including exact periodic copies.
/// Reflecting faces are stationary adiabatic impermeable slip walls: normal
/// viscous traction loads the support, while species, shear traction, heat and
/// work vanish. Transmissive Euler ends do not define a molecular reservoir
/// and refuse at preparation. All successful field stages reuse flat buffers.
pub struct MolecularTransport3d<L: TransportLaw> {
    source: SpatialTransport<L>,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    maximum_workspace_bytes: usize,
    boundary_impulses: Vec<f64>,
    last_step: Option<MolecularReport>,
    work: Workspace,
}
impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Combined numerical bytes for the supplied flow and prepared operator.
    /// Includes source/MS work, all stages/faces/gradients/boundary buffers and
    /// input state excess capacity. Immutable source/coefficient data, bounded
    /// stack scratch and arbitrary law storage/evaluation allocations are
    /// outside this numerical-buffer budget.
    pub fn required_workspace_bytes(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        let checked = || -> Option<usize> {
            let width = flow.gas.species_count().checked_add(4)?;
            let faces = flow
                .grid
                .face_len(0)
                .checked_add(flow.grid.face_len(1))?
                .checked_add(flow.grid.face_len(2))?;
            let scalars = flow
                .grid
                .len()
                .checked_mul(3)?
                .checked_add(faces)?
                .checked_add(flow.grid.boundary_len().checked_mul(4)?)?
                .checked_mul(width)?;
            scalars
                .checked_mul(size_of::<f64>())?
                .checked_add(
                    flow.grid
                        .len()
                        .checked_mul(size_of::<Primitive>() + size_of::<[[f64; 3]; 3]>())?,
                )?
                .checked_add(source.workspace_bytes())?
                .checked_add(flow.workspace_bytes())
        };
        checked().ok_or(DissipationError::Invalid(
            "volumetric molecular workspace size overflow",
        ))
    }

    /// Bind identical original source/domain, physical grid and boundaries.
    /// The limit covers this operator together with the supplied flow.
    pub fn new(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        if source.gas() != flow.gas() {
            return Err(DissipationError::Invalid(
                "volumetric molecular source or caloric path differs",
            ));
        }
        for boundary in [flow.boundaries.x, flow.boundaries.y, flow.boundaries.z] {
            if matches!(
                boundary,
                AxisBoundaries::Ends {
                    lower: BoundaryCondition::Transmissive,
                    ..
                } | AxisBoundaries::Ends {
                    upper: BoundaryCondition::Transmissive,
                    ..
                }
            ) {
                return Err(DissipationError::Invalid(
                    "transmissive Euler ends do not define a species/thermal molecular reservoir",
                ));
            }
        }
        admit_workspace(
            Self::required_workspace_bytes(flow, &source)?,
            maximum_workspace_bytes,
        )?;
        let grid = flow.grid;
        let width = flow.gas.species_count() + 4;
        let count = grid.len() * width;
        let boundary = grid.boundary_len() * width;
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
                gradient: allocated(grid.len(), [[0.0; 3]; 3])?,
                flux: [
                    allocated(grid.face_len(0) * width, 0.0)?,
                    allocated(grid.face_len(1) * width, 0.0)?,
                    allocated(grid.face_len(2) * width, 0.0)?,
                ],
                boundary_step: allocated(boundary, 0.0)?,
                boundary_trial: allocated(boundary, 0.0)?,
                boundary_correction: allocated(boundary, 0.0)?,
            },
        };
        admit_workspace(
            result.workspace_bytes() + flow.workspace_bytes(),
            maximum_workspace_bytes,
        )?;
        poll(&mut cancelled)?;
        Ok(result)
    }
    /// Exact original caloric source, species and admitted temperature interval.
    pub fn gas(&self) -> &MixtureGas {
        self.source.gas()
    }
    /// Original transport law and its source/material identities.
    pub fn law(&self) -> &L {
        self.source.law()
    }
    /// Actual owned numerical capacities; external flow/provider data excluded.
    pub fn workspace_bytes(&self) -> usize {
        self.source.workspace_bytes()
            + (self.boundary_impulses.capacity()
                + self.work.current.capacity()
                + self.work.stage.capacity()
                + self.work.candidate.capacity()
                + self.work.flux.iter().map(Vec::capacity).sum::<usize>()
                + self.work.boundary_step.capacity()
                + self.work.boundary_trial.capacity()
                + self.work.boundary_correction.capacity())
                * size_of::<f64>()
            + self.work.primitive.capacity() * size_of::<Primitive>()
            + self.work.gradient.capacity() * size_of::<[[f64; 3]; 3]>()
    }
    /// Last completely published requested-duration interval.
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

    /// Advance the complete requested interval atomically. Euler time and its
    /// advection report/loads continue to describe their original interval.
    pub fn advance(
        &mut self,
        flow: &mut MixtureEuler3d,
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
        flow: &mut MixtureEuler3d,
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
                "volumetric molecular field ownership differs",
            ));
        }
        let bytes = self
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "volumetric molecular workspace size overflow",
            ))?;
        admit_workspace(bytes, maximum_workspace_bytes)?;
        if !duration.value().is_finite() || duration.value() < 0.0 {
            return Err(DissipationError::Invalid(
                "nonnegative finite molecular duration required",
            ));
        }
        let species = self.gas().species_count();
        let width = species + 4;
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
        let mut elapsed = 0.0;
        let mut accepted_steps = 0;
        let mut rejected_steps = 0;
        let mut maximum_parabolic_number: f64 = 0.0;
        let mut maximum_gradient_residual_per_m: f64 = 0.0;
        let mut maximum_mass_frame_residual_kg_m2_s: f64 = 0.0;
        while elapsed < duration.value() {
            poll(cancelled)?;
            if accepted_steps + rejected_steps >= controls.maximum_steps {
                return Err(DissipationError::WorkLimit {
                    elapsed_time_s: elapsed,
                    accepted_steps,
                    rejected_steps,
                });
            }
            // This first source/gradient pass both estimates the interval and
            // supplies the first SSP stage. Only retries repeat it.
            let first = spatial(
                &mut self.source,
                &self.grid,
                self.boundaries,
                &self.work.current,
                &mut self.work.primitive,
                &mut self.work.gradient,
                &mut self.work.flux,
                controls,
                cancelled,
            )?;
            let remaining = duration.value() - elapsed;
            let mut dt = remaining.min(controls.maximum_step_s);
            if first.rate > 0.0 {
                dt = dt.min(controls.parabolic_cfl / first.rate);
            }
            for halving in 0..=controls.maximum_halvings {
                poll(cancelled)?;
                if accepted_steps + rejected_steps >= controls.maximum_steps {
                    return Err(DissipationError::WorkLimit {
                        elapsed_time_s: elapsed,
                        accepted_steps,
                        rejected_steps,
                    });
                }
                if !positive(dt)
                    || elapsed + dt <= elapsed
                    || (dt < controls.minimum_step_s && dt != remaining)
                {
                    return Err(DissipationError::Invalid(
                        "unrepresentable volumetric molecular substep",
                    ));
                }
                let attempt = (|| {
                    let first = if halving == 0 {
                        first
                    } else {
                        spatial(
                            &mut self.source,
                            &self.grid,
                            self.boundaries,
                            &self.work.current,
                            &mut self.work.primitive,
                            &mut self.work.gradient,
                            &mut self.work.flux,
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
                        &self.work.flux,
                        &mut self.work.stage,
                        dt,
                        false,
                        cancelled,
                    )?;
                    capture(
                        &self.grid,
                        width,
                        &self.work.flux,
                        0.5 * dt,
                        false,
                        &mut self.work.boundary_step,
                        cancelled,
                    )?;
                    let second = spatial(
                        &mut self.source,
                        &self.grid,
                        self.boundaries,
                        &self.work.stage,
                        &mut self.work.primitive,
                        &mut self.work.gradient,
                        &mut self.work.flux,
                        controls,
                        cancelled,
                    )?;
                    parabolic(second.rate * dt, controls.parabolic_cfl)?;
                    update(
                        self.source.gas(),
                        &self.grid,
                        &self.work.current,
                        &self.work.stage,
                        &self.work.flux,
                        &mut self.work.candidate,
                        dt,
                        true,
                        cancelled,
                    )?;
                    capture(
                        &self.grid,
                        width,
                        &self.work.flux,
                        0.5 * dt,
                        true,
                        &mut self.work.boundary_step,
                        cancelled,
                    )?;
                    Ok::<_, DissipationError>((first, second))
                })();
                match attempt {
                    Ok((first, second)) => {
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
                        maximum_parabolic_number =
                            maximum_parabolic_number.max(dt * first.rate.max(second.rate));
                        maximum_gradient_residual_per_m = maximum_gradient_residual_per_m
                            .max(first.gradient_residual)
                            .max(second.gradient_residual);
                        maximum_mass_frame_residual_kg_m2_s = maximum_mass_frame_residual_kg_m2_s
                            .max(first.mass_frame_residual)
                            .max(second.mass_frame_residual);
                        std::mem::swap(&mut self.work.current, &mut self.work.candidate);
                        elapsed = if dt == remaining {
                            duration.value()
                        } else {
                            elapsed + dt
                        };
                        accepted_steps += 1;
                        break;
                    }
                    Err(DissipationError::Cancelled) => return Err(DissipationError::Cancelled),
                    Err(error) => {
                        rejected_steps += 1;
                        if halving == controls.maximum_halvings
                            || dt * 0.5 < controls.minimum_step_s
                        {
                            return Err(DissipationError::StepFailure {
                                elapsed_time_s: elapsed,
                                step_s: dt,
                                last: Box::new(error),
                            });
                        }
                        dt *= 0.5;
                    }
                }
            }
        }
        let after = checked_totals(&self.work.current, &self.grid, species, cancelled)?;
        let outward = sum_rows(&self.work.boundary_trial, width, 1.0, cancelled)?;
        let boundary_change = Totals::from_row(&outward.map(|value| -value), species);
        let residual = std::array::from_fn::<_, WIDTH_LIMIT, _>(|d| {
            if d < width {
                after.component(d) - before.component(d) - boundary_change.component(d)
            } else {
                0.0
            }
        });
        if !residual.iter().all(|value| value.is_finite()) {
            return Err(DissipationError::Invalid(
                "volumetric molecular balance overflow",
            ));
        }
        let conservation_residual = Totals::from_row(&residual[..width], species);
        let report = MolecularReport {
            elapsed: duration,
            accepted_steps,
            rejected_steps,
            maximum_parabolic_number,
            maximum_gradient_residual_per_m,
            maximum_mass_frame_residual_kg_m2_s,
            before,
            after,
            boundary_change,
            conservation_residual,
            absolute_energy_residual_j: conservation_residual.absolute_energy(self.gas())?,
        };
        // No copied prefix can escape a late cancellation. Excess input
        // capacity follows the actual allocation across this O(1) exchange.
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
    source: &mut SpatialTransport<L>,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    states: &[f64],
    primitive: &mut [Primitive],
    gradient: &mut [[[f64; 3]; 3]],
    fluxes: &mut [Vec<f64>; 3],
    controls: DissipationControl,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Stage, DissipationError> {
    let species = source.gas().species_count();
    let width = species + 4;
    for (row, value) in states.chunks_exact(width).zip(primitive.iter_mut()) {
        *value = cell(row, species).primitive(source.gas(), &mut *cancelled)?;
    }
    // A single nine-derivative cache supplies both tangential directions of
    // every face. Slip ghosts reflect only the wall-normal velocity.
    for (index, row) in gradient.iter_mut().enumerate() {
        poll(cancelled)?;
        for (axis, derivative) in row.iter_mut().enumerate() {
            let (low, reflect_low) = cartesian::neighbour(grid, boundaries, index, axis, false);
            let (high, reflect_high) = cartesian::neighbour(grid, boundaries, index, axis, true);
            for (component, value) in derivative.iter_mut().enumerate() {
                let low = primitive[low].velocity[component]
                    * if reflect_low && component == axis {
                        -1.0
                    } else {
                        1.0
                    };
                let high = primitive[high].velocity[component]
                    * if reflect_high && component == axis {
                        -1.0
                    } else {
                        1.0
                    };
                *value = (high - low) / (2.0 * grid.spacing[axis]);
            }
        }
    }
    let normal = |index: usize, axis: usize| {
        let [a, b] = tangents(axis);
        SpatialPrimitive {
            partial_densities: &states[index * width..index * width + species],
            normal_velocity: primitive[index].velocity[axis],
            transverse_velocity: [primitive[index].velocity[a], primitive[index].velocity[b]],
            pressure: primitive[index].pressure,
        }
    };
    let mut stats = Stage::default();
    cartesian::faces::<DissipationError, _>(grid, boundaries, cancelled, |face, cancelled| {
        let axis = face.axis;
        let [a, b] = tangents(axis);
        let flux = &mut fluxes[axis];
        if let cartesian::FaceKind::PeriodicCopy { first } = face.kind {
            flux.copy_within(first * width..(first + 1) * width, face.index * width);
            return Ok(());
        }
        let output = &mut flux[face.index * width..(face.index + 1) * width];
        let receipt = match face.kind {
            cartesian::FaceKind::Interior { left, right } => {
                let tangential = [a, b].map(|derivative| {
                    [axis, a, b].map(|component| {
                        0.5 * gradient[left][derivative][component]
                            + 0.5 * gradient[right][derivative][component]
                    })
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
                [gradient[cell][a][a], gradient[cell][b][b]],
                output,
                &mut *cancelled,
            )?,
            cartesian::FaceKind::Boundary { .. } => {
                return Err(DissipationError::Invalid("unbound molecular reservoir"));
            }
            cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
        };
        // Restore the field's global momentum order from the normal/tangent
        // frame without changing species or referenced total-energy fluxes.
        let momentum = [output[species], output[species + 1], output[species + 2]];
        output[species + axis] = momentum[0];
        output[species + a] = momentum[1];
        output[species + b] = momentum[2];
        stats.add(receipt);
        Ok(())
    })?;
    // Covers all three axes, full cross-stress stencils and doubled wall
    // gradients. Both SSP states also undergo original source admission.
    // This engineering interval bound is not a positivity/entropy theorem.
    stats.rate *= 4.0 * grid.spacing.iter().map(|h| h.recip().powi(2)).sum::<f64>();
    if !stats.rate.is_finite() || stats.rate < 0.0 {
        return Err(DissipationError::Invalid(
            "volumetric molecular parabolic scale overflow",
        ));
    }
    Ok(stats)
}

#[allow(clippy::too_many_arguments)]
fn capture(
    grid: &CartesianGrid,
    width: usize,
    fluxes: &[Vec<f64>; 3],
    duration: f64,
    add: bool,
    output: &mut [f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DissipationError> {
    cartesian::capture::<DissipationError>(
        grid,
        width,
        duration,
        add,
        |axis, face, d| fluxes[axis][face * width + d],
        |face, d, value, add| {
            let out = &mut output[face * width + d];
            *out = if add { *out + value } else { value };
        },
        cancelled,
    )
}
fn admit_workspace(bytes: usize, maximum: usize) -> Result<(), DissipationError> {
    if bytes > maximum || bytes > MAX_MIXTURE_WORKSPACE_BYTES {
        Err(DissipationError::Invalid(
            "volumetric molecular total workspace budget",
        ))
    } else {
        Ok(())
    }
}
fn parabolic(observed: f64, maximum: f64) -> Result<(), DissipationError> {
    if !observed.is_finite() || observed > maximum * (1.0 + 32.0 * f64::EPSILON) {
        Err(DissipationError::ParabolicNumber { observed, maximum })
    } else {
        Ok(())
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), DissipationError> {
    stepping::poll::<DissipationError>(cancelled)
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
