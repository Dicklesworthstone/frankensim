//! Exact affine-force trajectories and their force-conjugate velocity moment.
//!
//! Physical state stays [u,p,y_tau...]. Force records are temporary source
//! coordinates of the prepared generator; each interval seeds them explicitly.
//! A one-way displacement integral reuses the native matrix exponential once,
//! giving every velocity moment without n separate Gramian preparations.

use super::*;
use fs_phs::linear_exact::exponential_with_budget;

/// Accepted physical motion plus the moment conjugate to a linear force ramp.
#[derive(Clone, Debug)]
pub struct PreparedTensorAffineStep {
    pub motion: PreparedTensorStep,
    /// Integral_0^1 s*v(h*s) ds [m/s]. For endpoint forces F0,F1, work is
    /// h*(F0 dot (mean_velocity-first_moment) + F1 dot first_moment).
    /// Native F(t) dot v(t) integration independently checks that identity.
    pub first_velocity_moment_mps: Vec<f64>,
}

pub(super) struct AffinePrepared {
    force_scale_n: f64,
    // Maps [D*x0, Delta F] to mean displacement, one row per physical coordinate.
    mean_displacement: Vec<f64>,
}

/// One native exponential for [conditioned physical/source state, held Delta F,
/// mean-displacement observers]. Observers have no feedback into the trajectory,
/// material state, stored energy or heat. The dimension is s+2*ports.
#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_observer(
    ports: usize,
    generator: &[f64],
    input: &[f64],
    scales: &[f64],
    dt_s: f64,
    force_scale_n: f64,
    budget: DenseExponentialBudget,
    gate: &CancelGate,
) -> Result<AffinePrepared> {
    poll(gate)?;
    let states = scales.len();
    let dimension = states + 2 * ports;
    let mut operator = zeros(dimension * dimension)?;
    for i in 0..states {
        poll(gate)?;
        for j in 0..states {
            operator[i * dimension + j] = dt_s * generator[i * states + j];
        }
        for p in 0..ports {
            operator[i * dimension + states + p] = dt_s * input[i * ports + p];
        }
    }
    for p in 0..ports {
        // d(observed D_u*mean u)/ds = D_u*u(hs): exact normalized integration.
        operator[(states + ports + p) * dimension + p] = 1.0;
    }
    poll(gate)?;
    let exponential = exponential_with_budget(&operator, dimension, budget)?;
    poll(gate)?;
    let width = states + ports;
    let mut mean_displacement = zeros(ports * width)?;
    for p in 0..ports {
        poll(gate)?;
        for j in 0..width {
            mean_displacement[p * width + j] =
                exponential[(states + ports + p) * dimension + j] / scales[p];
        }
    }
    if !finite(&mean_displacement) {
        return Err(invalid("affine displacement observer overflow"));
    }
    Ok(AffinePrepared {
        force_scale_n,
        mean_displacement,
    })
}

impl PreparedTensorHistory<'_> {
    /// Advance the actual continuous endpoint-interpolated force. Finite jumps
    /// between intervals are allowed and do not impulsively change any physical
    /// coordinate. This requires explicit affine preparation. The original
    /// physical state, clock and ledgers publish only after source-endpoint,
    /// velocity-moment work, positive heat and original energy checks all pass.
    pub fn advance_affine(
        &mut self,
        force_start_n: &[f64],
        force_end_n: &[f64],
        gate: &CancelGate,
    ) -> Result<PreparedTensorAffineStep> {
        poll(gate)?;
        let body = self.body;
        let source = body
            .affine
            .as_ref()
            .ok_or_else(|| invalid("affine forces require explicit affine preparation"))?;
        let n = body.coordinates();
        let physical = body.state_dimension();
        let propagated = body.propagation_state_dimension();
        if force_start_n.len() != n
            || force_end_n.len() != n
            || !finite(force_start_n)
            || !finite(force_end_n)
        {
            return Err(invalid("complete finite affine force endpoints required"));
        }
        let step = self
            .steps
            .checked_add(1)
            .filter(|s| *s <= body.controls.maximum_steps)
            .ok_or_else(|| invalid("accepted step budget exhausted"))?;
        let time_s = step as f64 * body.dt_s;
        if !time_s.is_finite() || time_s <= self.time_s() {
            return Err(invalid("representable later fixed clock required"));
        }
        self.scratch.source_state[..physical].copy_from_slice(&self.state);
        self.scratch.source_state[physical..].copy_from_slice(force_start_n);
        for i in 0..n {
            self.scratch.source_increment[i] = force_end_n[i] - force_start_n[i];
            if !self.scratch.source_increment[i].is_finite() {
                return Err(invalid("affine force increment is unrepresentable"));
            }
        }
        body.flow.evaluate_into(
            &self.scratch.source_state,
            &self.scratch.source_increment,
            &mut self.scratch.augmented,
            &mut self.scratch.balanced,
            &mut self.scratch.endpoint,
            &mut self.scratch.integrals,
        )?;
        poll(gate)?;
        for i in 0..n {
            let scale = source
                .force_scale_n
                .max(force_start_n[i].abs())
                .max(force_end_n[i].abs());
            let allowance = body.controls.maximum_relative_operator_error * scale;
            let defect = (self.scratch.endpoint[physical + i] - force_end_n[i]).abs();
            if !defect.is_finite() || defect > allowance {
                return Err(PreparedTensorError::Operator { defect, allowance });
            }
        }
        self.scratch
            .physical_candidate
            .copy_from_slice(&self.scratch.endpoint[..physical]);
        body.admit_displacement(&self.scratch.physical_candidate[..n], gate)?;
        let velocity = body
            .inertia
            .velocity(&self.scratch.physical_candidate[n..2 * n]);
        if !finite(&velocity) {
            return Err(invalid("candidate velocity overflow"));
        }
        let after = body.energy(&self.scratch.physical_candidate, &velocity, gate)?;
        let change = body.energy_change(
            &self.state,
            &self.velocity,
            &self.scratch.physical_candidate,
            &velocity,
            gate,
        )?;
        for i in 0..propagated {
            self.scratch.heat_state[i] = body.heat_scales[i] * self.scratch.source_state[i];
        }
        let mut mean_velocity_mps = zeros(n)?;
        let mut first_velocity_moment_mps = zeros(n)?;
        let width = propagated + n;
        let (mut moment_work, mut work_roundoff) = (0.0, 0.0);
        for i in 0..n {
            poll(gate)?;
            let row = &source.mean_displacement[i * width..(i + 1) * width];
            let initial = (0..propagated).fold(0.0, |sum, j| {
                row[j].mul_add(self.scratch.heat_state[j], sum)
            });
            let mean_displacement = (0..n).fold(initial, |sum, j| {
                row[propagated + j].mul_add(self.scratch.source_increment[j], sum)
            });
            let u1 = self.scratch.physical_candidate[i];
            mean_velocity_mps[i] = (u1 - self.state[i]) / body.dt_s;
            first_velocity_moment_mps[i] = (u1 - mean_displacement) / body.dt_s;
            // Endpoint weights express exactly the same source work without
            // treating a single mean velocity as conjugate to a varying force.
            add(
                &mut moment_work,
                &mut work_roundoff,
                body.dt_s
                    * force_start_n[i]
                    * (mean_velocity_mps[i] - first_velocity_moment_mps[i]),
            );
            add(
                &mut moment_work,
                &mut work_roundoff,
                body.dt_s * force_end_n[i] * first_velocity_moment_mps[i],
            );
        }
        let supplied_work_j = self.scratch.integrals[0];
        let moment_defect = moment_work - supplied_work_j;
        if !finite(&mean_velocity_mps)
            || !finite(&first_velocity_moment_mps)
            || !moment_defect.is_finite()
            || moment_defect.abs() > body.controls.maximum_work_error_j
        {
            return Err(PreparedTensorError::Energy {
                defect_j: moment_defect,
                allowance_j: body.controls.maximum_work_error_j,
            });
        }
        let mut heat_by_group_j = zeros(body.heat.len())?;
        let (mut dissipated_j, mut heat_roundoff) = (0.0, 0.0);
        for (out, cache) in heat_by_group_j.iter_mut().zip(&body.heat) {
            poll(gate)?;
            *out = cache.evaluate_into_controlled(
                &self.scratch.heat_state,
                &self.scratch.source_increment,
                &mut self.scratch.heat_coordinates,
                || gate.is_requested(),
            )?;
            if !out.is_finite() || *out < 0.0 {
                return Err(invalid(
                    "original source heat must be finite and nonnegative",
                ));
            }
            add(&mut dissipated_j, &mut heat_roundoff, *out);
        }
        let (mut work_j, mut work_correction_j) = (self.work_j, self.work_correction_j);
        let (mut heat_j, mut heat_correction_j) = (self.heat_j, self.heat_correction_j);
        add(&mut work_j, &mut work_correction_j, supplied_work_j);
        add(&mut heat_j, &mut heat_correction_j, dissipated_j);
        if !finite(&[
            supplied_work_j,
            dissipated_j,
            work_j,
            work_correction_j,
            heat_j,
            heat_correction_j,
        ]) || heat_j < self.heat_j
        {
            return Err(invalid("physical affine work/heat ledger overflow"));
        }
        let energy = PreparedTensorEnergy {
            equilibrium_before_j: self.energy.equilibrium,
            equilibrium_after_j: after.equilibrium,
            kinetic_before_j: self.energy.kinetic,
            kinetic_after_j: after.kinetic,
            memory_before_j: self.energy.memory,
            memory_after_j: after.memory,
            supplied_work_j,
            dissipated_j,
            cumulative_work_j: work_j,
            cumulative_dissipated_j: heat_j,
            balance_error_j: supplied_work_j - dissipated_j - change,
            cumulative_balance_error_j: work_j - heat_j - (after.total() - self.initial_energy_j),
        };
        for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
            if !defect_j.is_finite() || defect_j.abs() > body.controls.maximum_work_error_j {
                return Err(PreparedTensorError::Energy {
                    defect_j,
                    allowance_j: body.controls.maximum_work_error_j,
                });
            }
        }
        poll(gate)?;
        core::mem::swap(&mut self.state, &mut self.scratch.physical_candidate);
        self.velocity = velocity;
        self.energy = after;
        self.steps = step;
        self.work_j = work_j;
        self.work_correction_j = work_correction_j;
        self.heat_j = heat_j;
        self.heat_correction_j = heat_correction_j;
        Ok(PreparedTensorAffineStep {
            motion: PreparedTensorStep {
                step,
                time_s,
                mean_velocity_mps,
                heat_by_group_j,
                energy,
            },
            first_velocity_moment_mps,
        })
    }
}
