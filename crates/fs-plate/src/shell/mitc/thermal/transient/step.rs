//! Endpoint material changes in one permanent mechanical coordinate basis.
//!
//! With h=dt/2, the endpoint-average momentum equation solves
//! (M+h²K1) v_mid = M v0 + h [f_mid-(K0+K1)q0/2]. Both material
//! operators act on the SAME physical coordinates, reference inertia and
//! algebraic drilling constraint. No remap or temperature-rate lift enters.

use super::{
    EnrichedMitcShell, MitcTransientReport, MitcTransientState, MitcTransientStep,
    PreparedMitcMidpoint, action, check, finite, invalid, product, row_check,
};
use crate::PlateError;
use std::sync::Arc;

impl PreparedMitcMidpoint {
    /// Advance with the original accepted stiffness and this factor's endpoint
    /// stiffness. Both owners must share the exact reference coordinate basis,
    /// physical mass and numerical drilling constraint. This is a finite
    /// endpoint-average force scheme; material/thermal histories and their
    /// timestep refinement remain the caller's responsibility.
    ///
    /// The report exposes `q0^T(K1-K0)q1/2` separately from mechanical load
    /// work. A thermodynamic owner must use its original material potential;
    /// this parameter work is not itself a dissipative or reversible heat.
    pub fn advance_from(
        &self,
        before: &EnrichedMitcShell,
        initial: &MitcTransientState,
        load_before: &[f64],
        load_after: &[f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MitcTransientStep, PlateError> {
        check(&mut cancelled)?;
        let after = &self.shell;
        if !Arc::ptr_eq(&before.model, &after.model)
            || !Arc::ptr_eq(&before.mass, &after.mass)
            || !Arc::ptr_eq(&before.gauges, &after.gauges)
            || before.dimension() != after.dimension()
        {
            return Err(invalid(
                "endpoint dynamics requires one exact reference basis, physical inertia and gauge",
            ));
        }
        before.admit_state(initial)?;
        after.admit_state(initial)?;
        let n = after.dimension();
        if load_before.len() != n
            || load_after.len() != n
            || load_before.iter().chain(load_after).any(|v| !v.is_finite())
        {
            return Err(invalid(
                "one finite endpoint force per enriched shell coordinate required",
            ));
        }
        before.admit_gauge_load(load_before)?;
        after.admit_gauge_load(load_after)?;
        let target = finite(initial.time_s + self.dt_s)?;
        if target <= initial.time_s {
            return Err(invalid("midpoint time increment is not representable"));
        }
        let h = 0.5 * self.dt_s;
        let mut rhs = action(&after.mass, &initial.velocity, &mut cancelled)?;
        let k0q0 = action(&before.stiffness, &initial.coordinates, &mut cancelled)?;
        // The fixed-material entry point retains its original SpMV count and
        // exact zero parameter-work report. No temperature-derivative query is
        // needed by either mechanical path.
        let k1q0 = if Arc::ptr_eq(&before.stiffness, &after.stiffness) {
            None
        } else {
            Some(action(
                &after.stiffness,
                &initial.coordinates,
                &mut cancelled,
            )?)
        };
        let force: Vec<_> = load_before
            .iter()
            .zip(load_after)
            .map(|(a, b)| 0.5 * a + 0.5 * b)
            .collect();
        for i in 0..n {
            row_check(i, &mut cancelled)?;
            let stiffness_force = k1q0
                .as_ref()
                .map_or(k0q0[i], |next| 0.5 * k0q0[i] + 0.5 * next[i]);
            rhs[i] = finite(rhs[i] + h * (force[i] - stiffness_force))?;
        }
        let midpoint_velocity = self
            .factor
            .solve_reduced(&self.effective, &rhs, &mut cancelled)?;
        let mut q = vec![0.0; n];
        let mut v = vec![0.0; n];
        for i in 0..n {
            row_check(i, &mut cancelled)?;
            q[i] = finite(initial.coordinates[i] + self.dt_s * midpoint_velocity[i])?;
            v[i] = finite(2.0 * midpoint_velocity[i] - initial.velocity[i])?;
        }
        let state = after.state_reduced(target, q, v, &mut cancelled)?;
        let displacement: Vec<_> = state
            .coordinates
            .iter()
            .zip(&initial.coordinates)
            .map(|(b, a)| b - a)
            .collect();
        let k1q1 = action(&after.stiffness, &state.coordinates, &mut cancelled)?;
        let mean_stiffness_force: Vec<_> = k0q0
            .iter()
            .zip(k1q1)
            .map(|(a, b)| 0.5 * a + 0.5 * b)
            .collect();
        let stiffness_work = product(&displacement, &mean_stiffness_force)?;
        let stiffness_parameter_work_j = if let Some(next) = k1q0 {
            let delta_force: Vec<_> = next.iter().zip(&k0q0).map(|(b, a)| b - a).collect();
            finite(0.5 * product(&state.coordinates, &delta_force)?)?
        } else {
            0.0
        };
        // Polarization avoids subtracting two large quadratic energy datums.
        // The parameter work comes from the operator difference independently;
        // no mechanical residual is relabeled as a source or as material heat.
        let quadratic_energy_change_j = finite(stiffness_work + stiffness_parameter_work_j)?;
        let load_work_j = product(&displacement, &force)?;
        let velocity_change: Vec<_> = state
            .velocity
            .iter()
            .zip(&initial.velocity)
            .map(|(b, a)| b - a)
            .collect();
        let mean_velocity: Vec<_> = state
            .velocity
            .iter()
            .zip(&initial.velocity)
            .map(|(b, a)| 0.5 * a + 0.5 * b)
            .collect();
        let mean_momentum = action(&after.mass, &mean_velocity, &mut cancelled)?;
        let kinetic_energy_change_j = product(&velocity_change, &mean_momentum)?;
        // This is ΔKE+Δ(qKq/2)-W_parameter-W_load, evaluated without an
        // avoidable cancellation of the same parameter work twice.
        let energy_residual_j = finite(kinetic_energy_change_j + stiffness_work - load_work_j)?;
        check(&mut cancelled)?;
        Ok(MitcTransientStep {
            state,
            report: MitcTransientReport {
                kinetic_energy_change_j,
                load_work_j,
                quadratic_energy_change_j,
                stiffness_parameter_work_j,
                energy_residual_j,
            },
        })
    }
}
