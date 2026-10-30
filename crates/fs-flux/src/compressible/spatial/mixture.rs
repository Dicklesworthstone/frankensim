//! Volumetric shocks with evolving, source-owned chemical composition.
//!
//! Every species density, all three momenta and full referenced total energy
//! share the original mixture Riemann solve. The spatial topology, six-trace
//! limiter, conservative divergence and boundary quadrature are shared with
//! the frozen 3-D field. This is an Estimated stationary Cartesian lane; it
//! does not assert pressure-equilibrium preservation at unequal-caloric moving
//! contacts or a fully discrete entropy/invariant-domain theorem.

pub use super::super::mixture::MAX_MIXTURE_WORKSPACE_BYTES;
use super::super::{
    cartesian as stepping,
    mixture::{allocated, sum},
};
use super::{
    Boundaries3d, BoundaryCondition, BoundarySide, CartesianGrid, Primitive, cartesian, tangents,
};
pub use super::{Reconstruction, StepControl};
pub use fs_gas::multicomponent::{Error, MixtureGas};
use fs_gas::multicomponent::{FaceFlux, MAX_SPECIES, SpatialCell, SpatialPrimitive};
use fs_qty::Time;

const WIDTH_LIMIT: usize = MAX_SPECIES + 4;

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
    /// Largest combined x/y/z face CFL across both stages.
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
    /// `after - before - boundary_change`, component by component.
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

/// Evolving original species on a fixed physical three-dimensional grid.
///
/// Flat prepared buffers avoid heap allocation per cell or accepted step.
/// Density, gas constant, heat capacity and sound speed follow each evolving
/// composition; formation-energy references remain owned by `MixtureGas`.
pub struct MixtureEuler3d {
    gas: MixtureGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
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
                .checked_mul(width)?;
            scalars
                .checked_mul(size_of::<f64>())?
                .checked_add(grid.len().checked_mul(size_of::<Primitive>())?)?
                .checked_add(faces.checked_mul(size_of::<FaceFlux>())?)
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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        let bytes = Self::required_workspace_bytes(&gas, &grid)?
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
        let result = Self {
            gas,
            grid,
            boundaries,
            states,
            time: 0.0,
            last_step: None,
            workspace_bytes: bytes,
            boundary_impulses: allocated(grid.boundary_len() * width, 0.0)?,
            work: Workspace {
                primitive,
                slopes: allocated(3 * count, 0.0)?,
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
        let n = self.gas.species_count();
        let width = n + 4;
        let range = self.grid.boundary_range(side);
        self.boundary_impulses[range.start * width..range.end * width]
            .chunks_exact(width)
            .map(move |q| BoundaryImpulse {
                species_masses: &q[..n],
                momentum: [q[n], q[n + 1], q[n + 2]],
                energy: q[n + 3],
            })
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
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        stepping::poll::<Error>(&mut cancelled)?;
        let n = self.gas.species_count();
        let width = n + 4;
        let mut rate: f64 = 0.0;
        for (q, p) in self
            .states
            .chunks_exact(width)
            .zip(&mut self.work.primitive)
        {
            *p = cell(q, n).primitive(&self.gas, &mut cancelled)?;
            rate = rate.max(
                (0..3)
                    .map(|a| (p.velocity[a].abs() + p.sound_speed) / self.grid.spacing[a])
                    .sum(),
            );
        }
        let before = checked_totals(&self.states, &self.grid, n, &mut cancelled)?;
        let accepted = stepping::advance(
            maximum_duration.value(),
            self.time,
            rate,
            controls,
            &mut cancelled,
            |dt, reconstruction, cancelled| {
                attempt(
                    &self.gas,
                    &self.grid,
                    self.boundaries,
                    &self.states,
                    &mut self.work,
                    dt,
                    controls.cfl,
                    reconstruction,
                    cancelled,
                )
            },
        )?;
        let after = checked_totals(&self.work.candidate, &self.grid, n, &mut cancelled)?;
        let outward = sum_rows(&self.work.boundary_trial, width, 1.0, &mut cancelled)?;
        let boundary_change = Totals::from_row(&outward.map(|v| -v), n);
        let residual = std::array::from_fn::<_, WIDTH_LIMIT, _>(|d| {
            if d < width {
                after.component(d) - before.component(d) - boundary_change.component(d)
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
            absolute_energy_residual_j: conservation_residual.absolute_energy(&self.gas)?,
            conservation_residual,
        };
        stepping::poll::<Error>(&mut cancelled)?;
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
    original: &[f64],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    stepping::rk2(cancelled, |final_stage, cancelled| {
        let stats = spatial(
            gas,
            grid,
            boundaries,
            if final_stage { &work.stage } else { original },
            &mut work.primitive,
            &mut work.slopes,
            &mut work.flux,
            &mut work.metadata,
            reconstruction,
            cancelled,
        )?;
        stepping::check_cfl::<Error>(stats.rate * dt, cfl)?;
        if final_stage {
            update(
                gas,
                grid,
                original,
                &work.stage,
                &work.flux,
                &mut work.candidate,
                dt,
                true,
                cancelled,
            )?;
        } else {
            update(
                gas,
                grid,
                original,
                original,
                &work.flux,
                &mut work.stage,
                dt,
                false,
                cancelled,
            )?;
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
        Ok(stats)
    })
}

#[allow(clippy::too_many_arguments)]
fn spatial(
    gas: &MixtureGas,
    grid: &CartesianGrid,
    boundaries: Boundaries3d,
    states: &[f64],
    primitive: &mut [Primitive],
    slopes: &mut [f64],
    fluxes: &mut [Vec<f64>; 3],
    metadata: &mut [Vec<FaceFlux>; 3],
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<stepping::StageStats, Error> {
    let n = gas.species_count();
    let width = n + 4;
    for (q, p) in states.chunks_exact(width).zip(primitive.iter_mut()) {
        *p = cell(q, n).primitive(gas, &mut *cancelled)?;
    }
    let value = |k: usize, d: usize| {
        if d < n {
            states[k * width + d]
        } else if d < n + 3 {
            primitive[k].velocity[d - n]
        } else {
            primitive[k].pressure
        }
    };
    let mut trace = [0.0; WIDTH_LIMIT];
    let mut check = [0.0; WIDTH_LIMIT];
    let limited = cartesian::slopes::<Error>(
        grid,
        boundaries,
        width,
        [n, n + 1, n + 2],
        reconstruction,
        value,
        slopes,
        |k, axis, sign, slope| {
            for d in 0..width {
                trace[d] = value(k, d) + sign * slope[d];
            }
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
        let reconstruct = |k: usize, sign: f64, out: &mut [f64]| {
            for d in 0..width {
                out[d] = value(k, d) + sign * slopes[(3 * k + axis) * width + d];
            }
        };
        let q = |k: usize| cell(&states[k * width..(k + 1) * width], n).normal(axis);
        let first_order = reconstruction == Reconstruction::FirstOrder;
        let receipt = match face.kind {
            cartesian::FaceKind::Interior { left: l, right: r } => {
                if first_order {
                    gas.hllc_spatial_flux(q(l), q(r), &mut result[..width], &mut *cancelled)?
                } else {
                    reconstruct(l, 0.5, &mut left[..width]);
                    reconstruct(r, -0.5, &mut right[..width]);
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
            } => match (first_order, condition) {
                (true, BoundaryCondition::Reflecting) => gas.stationary_wall_spatial_flux(
                    q(k),
                    sign,
                    &mut result[..width],
                    &mut *cancelled,
                )?,
                (true, BoundaryCondition::Transmissive) => {
                    gas.hllc_spatial_flux(q(k), q(k), &mut result[..width], &mut *cancelled)?
                }
                (false, BoundaryCondition::Reflecting) => {
                    reconstruct(k, 0.5 * sign, &mut left[..width]);
                    gas.stationary_wall_spatial_flux_from_primitive(
                        normal_primitive(&left[..width], n, axis),
                        sign,
                        &mut result[..width],
                        &mut *cancelled,
                    )?
                }
                (false, BoundaryCondition::Transmissive) => {
                    reconstruct(k, 0.5 * sign, &mut left[..width]);
                    let p = normal_primitive(&left[..width], n, axis);
                    gas.hllc_spatial_flux_from_primitive(
                        p,
                        p,
                        &mut result[..width],
                        &mut *cancelled,
                    )?
                }
            },
            cartesian::FaceKind::PeriodicCopy { .. } => unreachable!(),
        };
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
    )?;
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
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    let n = gas.species_count();
    let width = n + 4;
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
            cell(row, n).primitive(gas, &mut *cancelled)?;
            output[k * width..(k + 1) * width].copy_from_slice(row);
            Ok(())
        },
    )
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
