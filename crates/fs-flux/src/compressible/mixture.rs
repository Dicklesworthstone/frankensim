//! Conservative spatial evolution of the original canonical gas species.
//!
//! Uses the SAME Cartesian MUSCL/SSP RK2, boundary traversal and retry core as
//! the frozen-composition field. `fs-gas::MixtureGas` alone owns calorics,
//! species energy references and planar HLLC/HLL/wall physics. All runtime
//! arrays are flat and prepared once; no cell owns a heap-allocated species Vec.
//! [`reactive`] composes the original finite-rate chemistry with this field.
//! Molecular diffusion is not yet composed with this evolving-species lane.

/// Source-owned finite-rate chemistry on the same conservative Cartesian field.
pub mod reactive;

use super::{
    Boundaries2d, BoundaryCondition, BoundarySide, CartesianGrid, Primitive, Reconstruction,
    StepControl, cartesian,
};
use cartesian::StageStats;
pub use fs_gas::multicomponent::{Error, MixtureGas};
use fs_gas::multicomponent::{FaceFlux, MAX_SPECIES, NormalCell, NormalPrimitive};
use fs_qty::Time;

const WIDTH_LIMIT: usize = MAX_SPECIES + 3;
/// Hard numerical-workspace ceiling, additionally bounded by the caller.
pub const MAX_MIXTURE_WORKSPACE_BYTES: usize = 1 << 30;

/// Borrowed accepted cell, in the original canonical species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cell<'a> {
    /// Species mass densities [kg/m³], including exact-zero species.
    pub partial_densities: &'a [f64],
    /// Cartesian x/y momentum densities [kg/(m² s)].
    pub momentum: [f64; 2],
    /// Total energy including both kinetic terms, less species references [J/m³].
    pub energy: f64,
}
impl<'a> Cell<'a> {
    fn normal(self, axis: usize) -> NormalCell<'a> {
        NormalCell {
            partial_densities: self.partial_densities,
            normal_momentum: self.momentum[axis],
            transverse_momentum: self.momentum[1 - axis],
            energy: self.energy,
        }
    }
    /// Reconstruct both velocities and original composition-dependent thermodynamics.
    pub fn primitive(
        self,
        gas: &MixtureGas,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Primitive, Error> {
        let p = gas.normal_primitive(self.normal(0), cancelled)?;
        Ok(Primitive {
            density: p.density,
            velocity: [p.velocity, self.momentum[1] / p.density],
            pressure: p.pressure,
            temperature: p.temperature,
            sound_speed: p.sound_speed,
        })
    }
}

/// Construct one flat global-frame cell through the original source owner.
/// `out` has `species_count()+3` entries: species, px, py, referenced energy.
pub fn state_into(
    gas: &MixtureGas,
    partial_densities: &[f64],
    velocity: [f64; 2],
    pressure: f64,
    out: &mut [f64],
) -> Result<(), Error> {
    gas.normal_conserved_into(
        NormalPrimitive {
            partial_densities,
            normal_velocity: velocity[0],
            transverse_velocity: velocity[1],
            pressure,
        },
        out,
    )?;
    Ok(())
}

/// Extensive conserved integrals in the original species order.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Totals {
    /// Species masses [kg]; entries outside `species_count` are exactly zero.
    pub species_masses: [f64; MAX_SPECIES],
    /// Number of original species represented by this value.
    pub species_count: usize,
    /// Cartesian momentum [kg m/s].
    pub momentum: [f64; 2],
    /// Species-referenced total energy [J].
    pub energy: f64,
}
impl Totals {
    /// Mixture mass is derived from the conserved species, never another field.
    pub fn mass(&self) -> f64 {
        sum(self.species_masses[..self.species_count].iter().copied())
    }
    /// Recover physical total chemical plus kinetic energy using original references.
    pub fn absolute_energy(&self, gas: &MixtureGas) -> Result<f64, Error> {
        if self.species_count != gas.species_count() {
            return Err(Error::Invalid("Cartesian total/source species axis"));
        }
        let value = sum(std::iter::once(self.energy).chain(
            self.species_masses[..self.species_count]
                .iter()
                .zip(gas.energy_references())
                .map(|(mass, reference)| mass * reference),
        ));
        if !value.is_finite() {
            return Err(Error::Invalid("Cartesian physical energy overflow"));
        }
        Ok(value)
    }
    fn from_row(row: &[f64], n: usize) -> Self {
        let mut species_masses = [0.0; MAX_SPECIES];
        species_masses[..n].copy_from_slice(&row[..n]);
        Self {
            species_masses,
            species_count: n,
            momentum: [row[n], row[n + 1]],
            energy: row[n + 2],
        }
    }
    fn component(&self, d: usize) -> f64 {
        if d < self.species_count {
            self.species_masses[d]
        } else if d < self.species_count + 2 {
            self.momentum[d - self.species_count]
        } else {
            self.energy
        }
    }
}

/// Borrowed outward transfer through one physical face over an accepted step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundaryImpulse<'a> {
    /// Outward species masses [kg].
    pub species_masses: &'a [f64],
    /// Outward Cartesian momentum [kg m/s], including wall support impulse.
    pub momentum: [f64; 2],
    /// Outward species-referenced total energy [J].
    pub energy: f64,
}

/// Only complete accepted steps publish time, cells, and boundary transfers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StepReport {
    /// Actual accepted duration.
    pub elapsed: Time,
    /// Maximum combined x/y face CFL over both SSP stages.
    pub maximum_cfl: f64,
    /// Rejected and halved attempts.
    pub halvings: usize,
    /// Accepted-stage HLL fallback count.
    pub hll_fallback_faces: usize,
    /// Accepted-stage cells whose common reconstruction theta was reduced.
    pub limited_cells: usize,
    /// A source-admissibility refusal switched the complete trial to first order.
    pub first_order_retry: bool,
    /// Initial extensive state.
    pub before: Totals,
    /// Accepted extensive state.
    pub after: Totals,
    /// Net transfer into the domain, integrated over both stages.
    pub boundary_change: Totals,
    /// `after - before - boundary_change`, component by component.
    pub conservation_residual: Totals,
    /// Physical energy residual including original per-species references [J].
    pub absolute_energy_residual_j: f64,
}

struct Workspace {
    primitive: Vec<Primitive>,
    slopes: Vec<f64>,
    x_flux: Vec<f64>,
    y_flux: Vec<f64>,
    x_metadata: Vec<FaceFlux>,
    y_metadata: Vec<FaceFlux>,
    stage: Vec<f64>,
    candidate: Vec<f64>,
    boundary_trial: Vec<f64>,
}

/// Evolving-species, source-owned Cartesian compressible flow.
///
/// Species mass equations determine mixture density; both momenta and full
/// referenced energy are conservative. Source laws remain immutable while
/// composition, pressure, temperature and sound speed evolve independently in
/// every cell. Admission bounds TOTAL component workspace before allocation.
pub struct MixtureEuler2d {
    gas: MixtureGas,
    grid: CartesianGrid,
    boundaries: Boundaries2d,
    states: Vec<f64>,
    time: f64,
    last_step: Option<StepReport>,
    boundary_impulses: Vec<f64>,
    workspace_bytes: usize,
    work: Workspace,
}
impl MixtureEuler2d {
    /// Numerical buffer requirement, including accepted cells, both stages,
    /// slopes, face fluxes/metadata, primitive cache and both boundary buffers.
    /// Source-model storage and bounded stack scratch are not included.
    pub fn required_workspace_bytes(
        gas: &MixtureGas,
        grid: &CartesianGrid,
    ) -> Result<usize, Error> {
        let checked = || -> Option<usize> {
            let [nx, ny] = grid.shape;
            let width = gas.species_count().checked_add(3)?;
            let faces = (nx + 1)
                .checked_mul(ny)?
                .checked_add(nx.checked_mul(ny + 1)?)?;
            let boundary = nx.checked_add(ny)?.checked_mul(2)?;
            let scalars = grid
                .len()
                .checked_mul(5)?
                .checked_add(faces)?
                .checked_add(boundary.checked_mul(2)?)?
                .checked_mul(width)?;
            scalars
                .checked_mul(size_of::<f64>())?
                .checked_add(grid.len().checked_mul(size_of::<Primitive>())?)?
                .checked_add(faces.checked_mul(size_of::<FaceFlux>())?)
        };
        checked().ok_or(Error::Invalid("Cartesian species workspace size overflow"))
    }
    /// Admit flat row-major cell averages and an explicit numerical memory budget.
    /// Layout per cell is original species, global px, global py, referenced E.
    pub fn new(
        gas: MixtureGas,
        grid: CartesianGrid,
        states: Vec<f64>,
        boundaries: Boundaries2d,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        cartesian::poll::<Error>(&mut cancelled)?;
        let bytes = Self::required_workspace_bytes(&gas, &grid)?
            .checked_add(
                (states.capacity() - states.len())
                    .checked_mul(size_of::<f64>())
                    .ok_or(Error::Invalid("Cartesian species input capacity overflow"))?,
            )
            .ok_or(Error::Invalid("Cartesian species workspace size overflow"))?;
        let width = gas.species_count() + 3;
        let count = grid
            .len()
            .checked_mul(width)
            .ok_or(Error::Invalid("Cartesian species state size overflow"))?;
        if bytes > maximum_workspace_bytes
            || bytes > MAX_MIXTURE_WORKSPACE_BYTES
            || states.len() != count
        {
            return Err(Error::Invalid(
                "Cartesian species state dimensions or total workspace budget",
            ));
        }
        let [nx, ny] = grid.shape;
        let mut primitive = allocated(grid.len(), Primitive::default())?;
        for (q, p) in states.chunks_exact(width).zip(&mut primitive) {
            *p = cell(q, width - 3).primitive(&gas, &mut cancelled)?;
        }
        checked_totals(&states, &grid, width - 3, &mut cancelled)?;
        let result = Self {
            gas,
            grid,
            boundaries,
            states,
            time: 0.0,
            last_step: None,
            workspace_bytes: bytes,
            boundary_impulses: allocated(2 * (nx + ny) * width, 0.0)?,
            work: Workspace {
                primitive,
                slopes: allocated(2 * count, 0.0)?,
                x_flux: allocated((nx + 1) * ny * width, 0.0)?,
                y_flux: allocated(nx * (ny + 1) * width, 0.0)?,
                x_metadata: allocated((nx + 1) * ny, FaceFlux::default())?,
                y_metadata: allocated(nx * (ny + 1), FaceFlux::default())?,
                stage: allocated(count, 0.0)?,
                candidate: allocated(count, 0.0)?,
                boundary_trial: allocated(2 * (nx + ny) * width, 0.0)?,
            },
        };
        cartesian::poll::<Error>(&mut cancelled)?;
        Ok(result)
    }
    /// Exact original species, elemental/charge axes, NASA laws and energy references.
    pub const fn gas(&self) -> &MixtureGas {
        &self.gas
    }
    /// Physical grid, including its extrusion depth.
    pub const fn grid(&self) -> &CartesianGrid {
        &self.grid
    }
    /// Admitted physical boundary conditions.
    pub const fn boundaries(&self) -> Boundaries2d {
        self.boundaries
    }
    /// Actual admitted numerical buffer requirement in bytes.
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }
    /// Flat accepted cells in the documented canonical layout.
    pub fn states(&self) -> &[f64] {
        &self.states
    }
    /// Number of physical cells.
    pub fn cell_count(&self) -> usize {
        self.grid.len()
    }
    /// Borrow one accepted physical cell without allocation.
    pub fn cell(&self, index: usize) -> Option<Cell<'_>> {
        let n = self.gas.species_count();
        let width = n + 3;
        (index < self.grid.len()).then(|| cell(&self.states[index * width..(index + 1) * width], n))
    }
    /// Accepted physical time.
    pub fn time(&self) -> Time {
        Time::new(self.time)
    }
    /// Most recent complete accepted step.
    pub const fn last_step(&self) -> Option<&StepReport> {
        self.last_step.as_ref()
    }
    /// Per-face outward transfer in increasing tangential face order.
    pub fn boundary_impulses(
        &self,
        side: BoundarySide,
    ) -> impl ExactSizeIterator<Item = BoundaryImpulse<'_>> {
        let n = self.gas.species_count();
        let width = n + 3;
        let range = self.grid.boundary_range(side);
        self.boundary_impulses[range.start * width..range.end * width]
            .chunks_exact(width)
            .map(move |q| BoundaryImpulse {
                species_masses: &q[..n],
                momentum: [q[n], q[n + 1]],
                energy: q[n + 2],
            })
    }
    /// Compensated extensive sums in fixed canonical cell/component order.
    pub fn totals(&self) -> Totals {
        checked_totals(
            &self.states,
            &self.grid,
            self.gas.species_count(),
            &mut || false,
        )
        .expect("admitted finite accepted totals")
    }

    /// Advance at most the requested duration with one shared x/y flux budget.
    pub fn step(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StepReport, Error> {
        cartesian::poll::<Error>(&mut cancelled)?;
        let n = self.gas.species_count();
        let width = n + 3;
        let mut rate: f64 = 0.0;
        for (q, p) in self
            .states
            .chunks_exact(width)
            .zip(&mut self.work.primitive)
        {
            *p = cell(q, n).primitive(&self.gas, &mut cancelled)?;
            rate = rate.max(
                (p.velocity[0].abs() + p.sound_speed) / self.grid.spacing[0]
                    + (p.velocity[1].abs() + p.sound_speed) / self.grid.spacing[1],
            );
        }
        let before = checked_totals(&self.states, &self.grid, n, &mut cancelled)?;
        let accepted = cartesian::advance(
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
            return Err(Error::Invalid("Cartesian species balance overflow"));
        }
        let conservation_residual = Totals::from_row(&residual[..width], n);
        let absolute_energy_residual_j = conservation_residual.absolute_energy(&self.gas)?;
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
            conservation_residual,
            absolute_energy_residual_j,
        };
        cartesian::poll::<Error>(&mut cancelled)?;
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
    boundaries: Boundaries2d,
    original: &[f64],
    work: &mut Workspace,
    dt: f64,
    cfl: f64,
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    cartesian::rk2(cancelled, |final_stage, cancelled| {
        let stats = spatial(
            gas,
            grid,
            boundaries,
            if final_stage { &work.stage } else { original },
            &mut work.primitive,
            &mut work.slopes,
            &mut work.x_flux,
            &mut work.y_flux,
            &mut work.x_metadata,
            &mut work.y_metadata,
            reconstruction,
            cancelled,
        )?;
        cartesian::check_cfl::<Error>(stats.rate * dt, cfl)?;
        if final_stage {
            update(
                gas,
                grid,
                original,
                &work.stage,
                &work.x_flux,
                &work.y_flux,
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
                &work.x_flux,
                &work.y_flux,
                &mut work.stage,
                dt,
                false,
                cancelled,
            )?;
        }
        let width = gas.species_count() + 3;
        cartesian::capture::<Error>(
            grid,
            width,
            0.5 * dt,
            final_stage,
            |axis, face, d| {
                if axis == 0 {
                    work.x_flux[face * width + d]
                } else {
                    work.y_flux[face * width + d]
                }
            },
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
    boundaries: Boundaries2d,
    states: &[f64],
    primitive: &mut [Primitive],
    slopes: &mut [f64],
    x_flux: &mut [f64],
    y_flux: &mut [f64],
    x_metadata: &mut [FaceFlux],
    y_metadata: &mut [FaceFlux],
    reconstruction: Reconstruction,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<StageStats, Error> {
    let n = gas.species_count();
    let width = n + 3;
    for (q, p) in states.chunks_exact(width).zip(primitive.iter_mut()) {
        *p = cell(q, n).primitive(gas, &mut *cancelled)?;
    }
    let value = |k: usize, d: usize| {
        if d < n {
            states[k * width + d]
        } else if d < n + 2 {
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
        [n, n + 1],
        reconstruction,
        value,
        slopes,
        |k, axis, sign, slope| {
            for d in 0..width {
                trace[d] = value(k, d) + sign * slope[d];
            }
            gas.normal_conserved_into(
                normal_primitive(&trace[..width], n, axis),
                &mut check[..width],
            )
            .is_ok()
        },
        cancelled,
    )?;
    let mut stats = StageStats {
        limited,
        ..StageStats::default()
    };
    let mut left = [0.0; WIDTH_LIMIT];
    let mut right = [0.0; WIDTH_LIMIT];
    let mut result = [0.0; WIDTH_LIMIT];
    cartesian::faces::<Error, _>(grid, boundaries, cancelled, |face, cancelled| {
        let axis = face.axis;
        let (fluxes, metadata) = if axis == 0 {
            (&mut *x_flux, &mut *x_metadata)
        } else {
            (&mut *y_flux, &mut *y_metadata)
        };
        if let cartesian::FaceKind::PeriodicCopy { first } = face.kind {
            for d in 0..width {
                fluxes[face.index * width + d] = fluxes[first * width + d];
            }
            metadata[face.index] = metadata[first];
            return Ok(());
        }
        let reconstruct = |k: usize, sign: f64, out: &mut [f64]| {
            for d in 0..width {
                out[d] = value(k, d) + sign * slopes[(2 * k + axis) * width + d];
            }
        };
        let q = |k: usize| cell(&states[k * width..(k + 1) * width], n).normal(axis);
        let first_order = reconstruction == Reconstruction::FirstOrder;
        let receipt = match face.kind {
            cartesian::FaceKind::Interior { left: l, right: r } => {
                if first_order {
                    gas.hllc_normal_flux(q(l), q(r), &mut result[..width], &mut *cancelled)?
                } else {
                    reconstruct(l, 0.5, &mut left[..width]);
                    reconstruct(r, -0.5, &mut right[..width]);
                    gas.hllc_normal_flux_from_primitive(
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
                (true, BoundaryCondition::Reflecting) => gas.stationary_wall_normal_flux(
                    q(k),
                    sign,
                    &mut result[..width],
                    &mut *cancelled,
                )?,
                (true, BoundaryCondition::Transmissive) => {
                    gas.hllc_normal_flux(q(k), q(k), &mut result[..width], &mut *cancelled)?
                }
                (false, BoundaryCondition::Reflecting) => {
                    reconstruct(k, 0.5 * sign, &mut left[..width]);
                    gas.stationary_wall_normal_flux_from_primitive(
                        normal_primitive(&left[..width], n, axis),
                        sign,
                        &mut result[..width],
                        &mut *cancelled,
                    )?
                }
                (false, BoundaryCondition::Transmissive) => {
                    reconstruct(k, 0.5 * sign, &mut left[..width]);
                    let p = normal_primitive(&left[..width], n, axis);
                    gas.hllc_normal_flux_from_primitive(
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
        if axis == 1 {
            result.swap(n, n + 1);
        }
        fluxes[face.index * width..(face.index + 1) * width].copy_from_slice(&result[..width]);
        metadata[face.index] = receipt;
        Ok(())
    })?;
    stats.rate = cartesian::rate::<Error>(
        grid,
        |axis, index| {
            if axis == 0 {
                x_metadata[index].signal_speed
            } else {
                y_metadata[index].signal_speed
            }
        },
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
    x_flux: &[f64],
    y_flux: &[f64],
    output: &mut [f64],
    dt: f64,
    final_stage: bool,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), Error> {
    let n = gas.species_count();
    let width = n + 3;
    cartesian::update::<WIDTH_LIMIT, Error, _>(
        grid,
        width,
        dt,
        final_stage,
        |k, d| original[k * width + d],
        |k, d| stage[k * width + d],
        |axis, face, d| {
            if axis == 0 {
                x_flux[face * width + d]
            } else {
                y_flux[face * width + d]
            }
        },
        cancelled,
        |k, row, cancelled| {
            cell(row, n).primitive(gas, &mut *cancelled)?;
            output[k * width..(k + 1) * width].copy_from_slice(row);
            Ok(())
        },
    )
}
fn normal_primitive(q: &[f64], n: usize, axis: usize) -> NormalPrimitive<'_> {
    NormalPrimitive {
        partial_densities: &q[..n],
        normal_velocity: q[n + axis],
        transverse_velocity: q[n + 1 - axis],
        pressure: q[n + 2],
    }
}
fn cell(q: &[f64], n: usize) -> Cell<'_> {
    Cell {
        partial_densities: &q[..n],
        momentum: [q[n], q[n + 1]],
        energy: q[n + 2],
    }
}
pub(super) fn sum(values: impl Iterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        let y = value - correction;
        let next = total + y;
        correction = (next - total) - y;
        total = next;
    }
    total
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
        cartesian::poll::<Error>(cancelled)?;
        for d in 0..width {
            let y = row[d] * scale - correction[d];
            let next = sums[d] + y;
            correction[d] = (next - sums[d]) - y;
            sums[d] = next;
        }
    }
    if !sums.iter().all(|v| v.is_finite()) {
        return Err(Error::Invalid("Cartesian species extensive overflow"));
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
        &sum_rows(states, n + 3, grid.cell_volume(), cancelled)?,
        n,
    ))
}
pub(super) fn allocated<T: Clone>(count: usize, value: T) -> Result<Vec<T>, Error> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(count)
        .map_err(|_| Error::Allocation)?;
    result.resize(count, value);
    Ok(result)
}

#[cfg(test)]
pub(crate) mod tests;
