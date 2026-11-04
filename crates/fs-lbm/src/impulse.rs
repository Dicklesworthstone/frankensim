//! Fixed-position impulse density on the native D2Q9 and D3Q19 populations.
//!
//! An impulse changes momentum now, without a collision, transport, or elapsed
//! timestep. For each fully occupied fluid cell, `rho * (u_plus-u_minus) = J`.
//! The external source supplies the signed kinetic work `J dot (u_plus+u_minus)/2`.
//! No heat, pressure impulse, boundary reaction, or free-surface mass is invented.
//! Subsequent native evolution resolves the resulting pressure/velocity field.
//!
//! Stored populations use the Guo frame `f = f_physical - S(u,F)/2`. Updating
//! just their first moment changes the viscous stress; resetting to equilibrium
//! erases it. The shared kernel instead adds the second-order equilibrium
//! difference AND changes the half-source at the new velocity. Consequently
//! density and the force-corrected nonequilibrium second moment are retained
//! to roundoff. This is a hydrodynamic impulse model, not a translation of all
//! kinetic moments or an incompressible pressure projection.

use fs_exec::{CancelGate, reduce::Compensated};

/// Explicit allocation, low-Mach, and independent readback limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImpulseControls {
    /// Maximum original grid cell count, including inactive cells.
    pub max_cells: usize,
    /// Maximum pre/post speed in lattice units, at most `MACH_LIMIT * sqrt(CS2)`.
    pub maximum_speed: f64,
    /// Sum of local L1 momentum defects plus the global L1 defect, lattice units.
    pub maximum_momentum_defect: f64,
    /// Sum of local absolute work defects plus the global defect, lattice units.
    pub maximum_work_defect: f64,
}

impl ImpulseControls {
    pub(crate) fn admit(self, count: usize) -> Result<(), ImpulseError> {
        if count == 0 || count > self.max_cells {
            return Err(ImpulseError::Invalid("impulse exceeds its cell budget"));
        }
        if !self.maximum_speed.is_finite()
            || self.maximum_speed <= 0.0
            || self.maximum_speed > crate::MACH_LIMIT * crate::CS2.sqrt()
            || !self.maximum_momentum_defect.is_finite()
            || self.maximum_momentum_defect < 0.0
            || !self.maximum_work_defect.is_finite()
            || self.maximum_work_defect < 0.0
        {
            return Err(ImpulseError::Invalid(
                "invalid impulse speed or balance budget",
            ));
        }
        Ok(())
    }
}

/// Independently observed impulse and signed work, summed over unit lattice cells.
/// In 2-D these are per unit out-of-plane lattice depth; in 3-D they are totals.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct ImpulseReport<const D: usize> {
    /// Sum of the prescribed impulse densities; no multiplication by dt.
    pub applied_impulse: [f64; D],
    /// Actual change in force-corrected population momentum.
    pub momentum_change: [f64; D],
    /// Actual mass change (roundoff only), never an authored mass source.
    pub mass_change: f64,
    /// Prescribed impulse dotted with the actual mean endpoint velocity.
    pub work: f64,
    /// Kinetic change read back from actual populations, evaluated without
    /// subtracting large endpoint kinetic energies.
    pub kinetic_energy_change: f64,
    /// Local absolute plus global momentum defect, in lattice momentum units.
    pub absolute_momentum_defect: f64,
    /// Local absolute plus global kinetic/work defect, in lattice energy units.
    pub absolute_work_defect: f64,
    /// Number of fluid cells with a nonzero prescribed impulse.
    pub affected_cells: usize,
}

/// A refused impulse never publishes a partial fluid state.
#[derive(Clone, Debug, PartialEq)]
pub enum ImpulseError {
    /// Invalid shape, allocation, control, or unsupported cell type.
    Invalid(&'static str),
    /// Invalid populations, force or impulse, or exceeded low-Mach speed.
    Cell {
        /// Original x-fastest cell index, including inactive cells.
        index: usize,
    },
    /// Actual momentum/work readback exceeded the explicit budget.
    Balance,
    /// Cancellation at a bounded cell/tile boundary or before publication.
    Cancelled,
}
impl core::fmt::Display for ImpulseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "fluid momentum impulse: {self:?}")
    }
}
impl std::error::Error for ImpulseError {}

pub(crate) fn poll(gate: &CancelGate) -> Result<(), ImpulseError> {
    if gate.is_requested() {
        Err(ImpulseError::Cancelled)
    } else {
        Ok(())
    }
}

pub(crate) struct ImpulseSum<const D: usize> {
    applied: [Compensated; D],
    change: [Compensated; D],
    mass: Compensated,
    work: Compensated,
    kinetic: Compensated,
    momentum_defect: Compensated,
    work_defect: Compensated,
    count: usize,
}
impl<const D: usize> ImpulseSum<D> {
    pub(crate) fn new() -> Self {
        Self {
            applied: [Compensated::zero(); D],
            change: [Compensated::zero(); D],
            mass: Compensated::zero(),
            work: Compensated::zero(),
            kinetic: Compensated::zero(),
            momentum_defect: Compensated::zero(),
            work_defect: Compensated::zero(),
            count: 0,
        }
    }

    pub(crate) fn finish(
        self,
        controls: ImpulseControls,
    ) -> Result<ImpulseReport<D>, ImpulseError> {
        let applied_impulse = self.applied.map(Compensated::value);
        let momentum_change = self.change.map(Compensated::value);
        let work = self.work.value();
        let kinetic_energy_change = self.kinetic.value();
        let absolute_momentum_defect = self.momentum_defect.value()
            + (0..D)
                .map(|a| (momentum_change[a] - applied_impulse[a]).abs())
                .sum::<f64>();
        let absolute_work_defect = self.work_defect.value() + (kinetic_energy_change - work).abs();
        let mass_change = self.mass.value();
        if !applied_impulse
            .iter()
            .chain(&momentum_change)
            .all(|x| x.is_finite())
            || ![
                work,
                kinetic_energy_change,
                mass_change,
                absolute_momentum_defect,
                absolute_work_defect,
            ]
            .iter()
            .all(|x| x.is_finite())
            || absolute_momentum_defect > controls.maximum_momentum_defect
            || absolute_work_defect > controls.maximum_work_defect
        {
            return Err(ImpulseError::Balance);
        }
        Ok(ImpulseReport {
            applied_impulse,
            momentum_change,
            mass_change,
            work,
            kinetic_energy_change,
            absolute_momentum_defect,
            absolute_work_defect,
            affected_cells: self.count,
        })
    }

    /// Shared pointwise second-order Hermite update. D2Q9 gravity contributes
    /// rho*g to the force; D3Q19 passes zero acceleration and its held force density.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn kick<const Q: usize>(
        &mut self,
        populations: &mut [f64; Q],
        directions: &[[f64; D]; Q],
        weights: &[f64; Q],
        acceleration: [f64; D],
        external_force: [f64; D],
        impulse: [f64; D],
        controls: ImpulseControls,
        index: usize,
    ) -> Result<(), ImpulseError> {
        let refuse = || ImpulseError::Cell { index };
        if !impulse.iter().all(|x| x.is_finite()) {
            return Err(refuse());
        }
        if impulse == [0.0; D] {
            return Ok(());
        }
        let observe = |f: &[f64; Q]| {
            let rho = f.iter().sum::<f64>();
            let force: [f64; D] =
                core::array::from_fn(|a| acceleration[a].mul_add(rho, external_force[a]));
            let momentum: [f64; D] = core::array::from_fn(|a| {
                f.iter().zip(directions).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a]
            });
            (rho, force, momentum, momentum.map(|p| p / rho))
        };
        let (rho, force, before, velocity) = observe(populations);
        if !rho.is_finite()
            || rho <= 0.0
            || !populations
                .iter()
                .chain(&force)
                .chain(&velocity)
                .all(|x| x.is_finite())
        {
            return Err(refuse());
        }
        let jump = impulse.map(|j| j / rho);
        let target: [f64; D] = core::array::from_fn(|a| velocity[a] + jump[a]);
        let speed_sq = |u: [f64; D]| u.iter().map(|x| x * x).sum::<f64>();
        let limit_sq = controls.maximum_speed * controls.maximum_speed;
        if !target.iter().all(|x| x.is_finite())
            || speed_sq(velocity) > limit_sq
            || speed_sq(target) > limit_sq
        {
            return Err(refuse());
        }
        let midpoint: [f64; D] = core::array::from_fn(|a| velocity[a] + 0.5 * jump[a]);
        let dot = |a: [f64; D], b: [f64; D]| (0..D).map(|d| a[d] * b[d]).sum::<f64>();
        let uj = dot(midpoint, impulse);
        let df = dot(jump, force);
        let inv_cs2 = 1.0 / crate::CS2;
        for ((f, &e), &w) in populations.iter_mut().zip(directions).zip(weights) {
            let ej = dot(e, impulse);
            // Evaluate differences algebraically: no subtraction of nearly
            // equal equilibria for a small impulse or mean-flow boost.
            let equilibrium_change =
                ej * inv_cs2 + dot(e, midpoint) * ej * inv_cs2 * inv_cs2 - uj * inv_cs2;
            let half_source_change =
                0.5 * (df * inv_cs2 - dot(e, jump) * dot(e, force) * inv_cs2 * inv_cs2);
            *f += w * (equilibrium_change + half_source_change);
        }
        // The exact update has zero zeroth moment. The rest direction is zero
        // in both native stencils, so this repairs mass roundoff without a kick.
        populations[0] += rho - populations.iter().sum::<f64>();
        let (rho_after, _, after, velocity_after) = observe(populations);
        if !rho_after.is_finite()
            || rho_after <= 0.0
            || !populations
                .iter()
                .chain(&after)
                .chain(&velocity_after)
                .all(|x| x.is_finite())
            || speed_sq(velocity_after) > limit_sq * (1.0 + 64.0 * f64::EPSILON)
            || (rho_after - rho).abs() > 64.0 * f64::EPSILON * rho
        {
            return Err(refuse());
        }
        let mean_velocity = core::array::from_fn(|a| 0.5 * velocity[a] + 0.5 * velocity_after[a]);
        let change = core::array::from_fn(|a| after[a] - before[a]);
        let mass_change = rho_after - rho;
        let work = dot(impulse, mean_velocity);
        // Exact finite-difference identity for K = |p|^2/(2 rho), including
        // any remaining mass representation roundoff in the independent readback.
        let kinetic =
            dot(change, mean_velocity) - 0.5 * dot(velocity, velocity_after) * mass_change;
        for a in 0..D {
            self.applied[a] = self.applied[a].accumulate(impulse[a]);
            self.change[a] = self.change[a].accumulate(change[a]);
            self.momentum_defect = self
                .momentum_defect
                .accumulate((change[a] - impulse[a]).abs());
        }
        self.mass = self.mass.accumulate(mass_change);
        self.work = self.work.accumulate(work);
        self.kinetic = self.kinetic.accumulate(kinetic);
        self.work_defect = self.work_defect.accumulate((kinetic - work).abs());
        self.count += 1;
        Ok(())
    }
}
