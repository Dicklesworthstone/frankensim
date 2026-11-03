//! Prescribed impulses in the original mass-normalized mechanical coordinates.
//!
//! Integrating M v' = F through an ideal impulse gives v+ = v- + J because
//! these coordinates have unit modal mass. Finite spring/dashpot forces carry
//! no impulse. Displacement, Maxwell internal strain and the proper passive
//! radiation model's states remain continuous. A rigid constraint, impact law,
//! inertance feedthrough or nonlinear contact needs its own coupled reset.

use super::{
    CoupledModalSystem, ModalCouplingError, component_energy, extension, finite, invalid, limit,
    poll,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;

/// Accepted external work at one instantaneous momentum transfer.
/// This is separate from the preceding finite-time frame and port window.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalImpulseReceipt {
    /// Number of completed samples, unchanged by the impulse.
    pub sample: u64,
    /// Stable kinetic increment from actual represented velocity changes [J].
    pub kinetic_energy_change_j: f64,
    /// Requested impulse dotted with mean endpoint modal velocity [J].
    pub external_work_j: f64,
    /// Sum of (represented - requested) impulse times mean velocity [J].
    pub work_residual_j: f64,
    /// Largest coordinate momentum defect, scaled only by that impulse.
    pub maximum_relative_momentum_residual: f64,
}

impl CoupledModalSystem {
    /// Apply simultaneous signed generalized impulses [N s / sqrt(kg)] in
    /// component order, then original mode order. Physical port columns map
    /// newton-seconds into these coordinates exactly as they map held forces.
    /// Inputs already contain the sum of every same-time load on each mode.
    ///
    /// This retains all displacements, constitutive histories, acoustic states,
    /// prepared operators, held inputs, clocks and preceding frame/port records.
    /// Only mechanical velocities change. The returned work is external work;
    /// no physical or numerical dissipation is assigned to this zero-duration
    /// event. Its source/contact model remains the caller's responsibility.
    ///
    /// Every original component/whole-system energy, pressure and endpoint
    /// connection-force limit applies. `solve_relative_tolerance` bounds each
    /// represented momentum change relative to the requested/actual impulse;
    /// a small kick swallowed by large background velocity refuses. The work
    /// gate uses the original energy tolerances on the jump, not total storage.
    /// Scratch is reserved with the network. Staging allocates nothing and
    /// cancellation or any refusal leaves the complete accepted state intact.
    pub fn apply_modal_impulse(
        &mut self,
        impulse_n_s_per_sqrt_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<ModalImpulseReceipt, ModalCouplingError> {
        let receipt = self.stage_modal_impulse(impulse_n_s_per_sqrt_kg, Some(gate))?;
        poll(Some(gate))?;
        self.publish_modal_impulse();
        Ok(receipt)
    }

    // The scheduler stages every voice at the same timestamp before publishing
    // any of them. This boundary cannot escape the native rendering owner.
    pub(in crate::render::schedule::force) fn stage_modal_impulse(
        &mut self,
        impulses: &[f64],
        gate: Option<&CancelGate>,
    ) -> Result<ModalImpulseReceipt, ModalCouplingError> {
        poll(gate)?;
        self.require_time_domain()?;
        if impulses.len() != self.mode_count() || impulses.iter().any(|j| !j.is_finite()) {
            return Err(invalid(
                "modal impulse must cover the complete finite mechanical basis",
            ));
        }
        self.staged_material_trajectory = None;
        let (mut kinetic, mut work, mut defect) = (Sum::default(), Sum::default(), Sum::default());
        let mut worst = 0.0_f64;
        let mut pressure = 0.0;
        for component in 0..self.models.len() {
            poll(gate)?;
            let range = self.offsets[component]..self.offsets[component + 1];
            let mut component_pressure = 0.0;
            for ((index, old), mode) in range
                .clone()
                .zip(self.models[component].states())
                .zip(self.models[component].modes())
            {
                if index % 64 == 0 {
                    poll(gate)?;
                }
                let before = old.velocity_m_sqrt_kg_per_s;
                let after = finite(before + impulses[index])?;
                let change = finite(after - before)?;
                let residual = finite(change - impulses[index])?;
                let scale = impulses[index].abs().max(change.abs());
                let relative = if scale == 0.0 {
                    0.0
                } else {
                    residual.abs() / scale
                };
                if relative > self.config.solve_relative_tolerance {
                    return Err(ModalCouplingError::ImpulseResidual {
                        mode: index,
                        relative,
                        tolerance: self.config.solve_relative_tolerance,
                    });
                }
                worst = worst.max(relative);
                // Keep the half-ULP mean and each product's rounding tail:
                // opposite large works can leave only the quadratic kick energy.
                let mean = Dd::from_f64(0.5 * before) + Dd::from_f64(0.5 * after);
                kinetic.add(Dd::from_f64(change) * mean)?;
                work.add(Dd::from_f64(impulses[index]) * mean)?;
                // Subtract before multiplying: large kinetic backgrounds or
                // cancelling modal works cannot hide a lost momentum input.
                defect.add(Dd::from_f64(residual) * mean)?;
                self.impulse_states[index] = crate::modal_acoustic_time::ModalAcousticState {
                    velocity_m_sqrt_kg_per_s: after,
                    ..*old
                };
                component_pressure = finite(
                    component_pressure
                        + mode.pressure_per_modal_velocity.re * after
                        + mode.pressure_per_modal_velocity.im
                            * mode.angular_frequency_rad_s
                            * old.displacement_m_sqrt_kg,
                )?;
            }
            self.candidates[component]
                .restore_states(&self.impulse_states[range])
                .map_err(|source| ModalCouplingError::Component { component, source })?;
            limit(
                "component observer pressure after impulse",
                component_pressure.abs(),
                self.models[component].budget().maximum_abs_pressure_pa,
            )?;
            pressure = finite(pressure + component_pressure)?;
        }
        let mut energy = component_energy(&self.candidates)?;
        for (index, (column, link)) in self.columns.iter().zip(&self.connections).enumerate() {
            poll(gate)?;
            let displacement = extension(&self.candidates, column, link.rest_extension_m)?;
            energy = finite(energy + 0.5 * link.stiffness_n_m * displacement * displacement)?;
            let dashpot = self.connection_damping_column(index).unwrap_or(column);
            let mut rate = 0.0;
            for (mode, (&shape, state)) in dashpot.iter().zip(&self.impulse_states).enumerate() {
                if mode % 64 == 0 {
                    poll(gate)?;
                }
                rate = finite(rate + shape * state.velocity_m_sqrt_kg_per_s)?;
            }
            let reaction = finite(-link.stiffness_n_m * displacement - link.damping_n_s_m * rate)?;
            limit(
                "connection force after impulse",
                reaction.abs(),
                self.config.maximum_abs_connection_force_n,
            )?;
        }
        // These histories depend on unchanged q and z, not the velocity jump.
        energy = finite(energy + self.material_energy_j()? + self.radiation_energy_j()?)?;
        limit(
            "total energy after impulse",
            energy,
            self.config.maximum_total_energy_j,
        )?;
        limit(
            "observer pressure after impulse",
            pressure.abs(),
            self.config.maximum_abs_pressure_pa,
        )?;
        let kinetic_energy_change_j = kinetic.value()?;
        let external_work_j = work.value()?;
        let work_residual_j = defect.value()?;
        let tolerance = finite(
            self.config.energy_absolute_tolerance_j
                + self.config.energy_relative_tolerance
                    * kinetic_energy_change_j.abs().max(external_work_j.abs()),
        )?;
        if work_residual_j.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance {
                residual_j: work_residual_j,
                tolerance_j: tolerance,
            });
        }
        poll(gate)?;
        Ok(ModalImpulseReceipt {
            sample: self.samples_rendered(),
            kinetic_energy_change_j,
            external_work_j,
            work_residual_j,
            maximum_relative_momentum_residual: worst,
        })
    }

    pub(in crate::render::schedule::force) fn publish_modal_impulse(&mut self) {
        core::mem::swap(&mut self.models, &mut self.candidates);
        self.staged_material_trajectory = None;
    }
}

// Keep products and signed component reductions in the same extended value;
// compensation after rounding a midpoint/product cannot recover its lost work.
#[derive(Default)]
struct Sum(Dd);
impl Sum {
    fn add(&mut self, term: Dd) -> Result<(), ModalCouplingError> {
        self.0 = self.0 + term;
        finite(self.0.hi)?;
        finite(self.0.lo)?;
        Ok(())
    }
    fn value(self) -> Result<f64, ModalCouplingError> {
        finite(self.0.to_f64())
    }
}
