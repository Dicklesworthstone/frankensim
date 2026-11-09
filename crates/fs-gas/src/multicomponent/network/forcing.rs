//! Native body sources inside the common branch/plenum transaction.

use super::{MixtureNetwork1d, NetworkStepReport};
use crate::body_force::{Accumulator, BodyForce, ForceChange, ForceControl};
use crate::multicomponent::evolve::forcing::update;
use crate::multicomponent::{Cell, Error, StepControl, allocated, poll, sum};
use fs_qty::Time;

/// Prescribed axial acceleration or force density in retained branch/cell order.
///
/// Each branch has its own physical axial coordinate; impulses on unrelated
/// axes are never summed as a network momentum. Linear histories use the common
/// accepted clock. Compact plenums retain their stationary perfectly mixed
/// model: this input does not invent a bulk velocity or accelerate a reservoir.
#[derive(Clone, Copy, Debug)]
pub struct NetworkBodyForces<'a> {
    /// Exactly one original body-force law per branch.
    pub branches: &'a [BodyForce<'a>],
    /// Physical Mach admission at source endpoints and interior speed extrema.
    pub controls: ForceControl,
}

/// Small construction-time receipts; full source fields reuse each branch's
/// original lazy body-force buffer. No second population/state owner is made.
#[derive(Clone, Copy, Default)]
pub(super) struct BranchSource {
    first: ForceChange,
    second: ForceChange,
    pub(super) maximum_rate: f64,
}
impl BranchSource {
    pub(super) fn change(self) -> ForceChange {
        ForceChange {
            impulse: sum([self.first.impulse, self.second.impulse]),
            work: sum([self.first.work, self.second.work]),
            energy_residual: sum([self.first.energy_residual, self.second.energy_residual]),
        }
    }
}

impl MixtureNetwork1d {
    /// Advance distributed body forces and all native finite ports atomically.
    ///
    /// Exact source half maps surround the shared SSP RK2 transport. Every
    /// characteristic branch/plenum flux sees the actually accelerated branch
    /// state; its opposite species and enthalpy exchange enters that same
    /// plenum stage. Both source halves use their actual local mixture mass.
    /// All branches, plenums, clocks and receipts retry and publish together.
    ///
    /// One original flat source buffer per branch is allocated on first active
    /// use and reused. All-zero forcing takes the original path bit for bit.
    /// A prescribed force supplies mechanical work, not heat or chemistry.
    /// General source/transport splitting requires temporal refinement; this
    /// method is not an exact hydrostatic-equilibrium discretization.
    pub fn step_with_body_forces(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        forces: NetworkBodyForces<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        self.step_with_body_forces_and_boundaries(
            maximum_duration, controls, forces, None, None, cancelled,
        )
    }

    /// Combine branch acceleration with prescribed end reservoirs and plenum motion.
    ///
    /// Optional reservoir pairs and volume rates retain the ordering, admission
    /// and physical laws of the original network boundary methods. Sources act
    /// only in resolved branch cells. Plenum pressure work, exterior transported
    /// enthalpy and distributed source work remain separate in the first law.
    /// No intermediate branch or boundary state is published.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_body_forces_and_boundaries(
        &mut self,
        maximum_duration: Time,
        controls: StepControl,
        forces: NetworkBodyForces<'_>,
        reservoirs: Option<&[[Option<Cell<'_>>; 2]]>,
        volume_rates_m3_s: Option<&[f64]>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkStepReport, Error> {
        let active = self.admit_body_forces(forces, &mut cancelled)?;
        let reservoirs = match reservoirs {
            Some(values) if self.admit_reservoirs(values, &mut cancelled)? => Some(values),
            _ => None,
        };
        let volume_rates = match volume_rates_m3_s {
            Some(values) if self.admit_volume_rates(values, &mut cancelled)? => Some(values),
            _ => None,
        };
        self.step_forcing_impl(
            maximum_duration, controls, volume_rates, reservoirs,
            active.then_some(forces), None, cancelled,
        )
    }

    pub(super) fn admit_body_forces(
        &mut self,
        forces: NetworkBodyForces<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, Error> {
        poll(cancelled)?;
        forces.controls.admit::<Error>()?;
        if forces.branches.len() != self.branches.len() {
            return Err(Error::Invalid("one body-force law per network branch required"));
        }
        let mut active = false;
        for (branch, force) in self.branches.iter().zip(forces.branches) {
            active |= force.admit::<Error>(branch.cell_count(), cancelled)?;
        }
        if active {
            for branch in &mut self.branches {
                poll(cancelled)?;
                if branch.body_force_stage.is_none() {
                    branch.body_force_stage = Some(allocated(branch.states.len(), 0.0)?);
                }
            }
        }
        Ok(active)
    }

    pub(super) fn force_half(
        &mut self,
        forces: NetworkBodyForces<'_>,
        dt: f64,
        cfl: f64,
        second: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<f64, Error> {
        let mut maximum_rate = 0.0_f64;
        for (index, branch) in self.branches.iter_mut().enumerate() {
            let mut exchange = Accumulator::default();
            let m = branch.gas.species_count() + 2;
            let destination = if second {
                &mut branch.work.candidate
            } else {
                let forced = branch.body_force_stage.as_mut()
                    .expect("active network force buffers were admitted");
                forced.copy_from_slice(&branch.states);
                forced
            };
            let mut branch_rate = 0.0_f64;
            for (i, next) in destination.chunks_exact_mut(m).enumerate() {
                let peak = update(
                    &branch.gas, next, &mut branch.work.primitive[i], i,
                    forces.branches[index],
                    [branch.time, if second { 0.5 * dt } else { 0.0 }, 0.5 * dt],
                    forces.controls, branch.mesh.volumes()[i], &mut exchange, cancelled,
                )?;
                let area = 0.5 * branch.mesh.face_areas()[i]
                    + 0.5 * branch.mesh.face_areas()[i + 1];
                let scale = area / branch.mesh.volumes()[i];
                for state in core::iter::once(branch.work.primitive[i]).chain(peak) {
                    let rate = (state.velocity.abs() + state.sound_speed) * scale;
                    let observed = dt * rate;
                    if !observed.is_finite()
                        || observed > cfl * (1.0 + 64.0 * f64::EPSILON)
                    {
                        return Err(Error::CflExceeded { observed, maximum: cfl });
                    }
                    branch_rate = branch_rate.max(rate);
                }
            }
            let receipt = &mut self.work.branch_sources[index];
            if second { receipt.second = exchange.finish(); }
            else { receipt.first = exchange.finish(); }
            receipt.maximum_rate = receipt.maximum_rate.max(branch_rate);
            maximum_rate = maximum_rate.max(branch_rate);
        }
        Ok(maximum_rate)
    }
}
