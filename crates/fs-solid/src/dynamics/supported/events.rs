//! Located friction events with original prescribed-support reaction work.

use super::sustained::{Interval, Workspace};
use super::*;
use crate::dynamics::spectral::polynomial::RigidForceCorrection;
use crate::dynamics::sustained as contact;
use contact::events::{
    CollinearCoulombContact, ContactForcePath, EventNative, SolidCoulombEventPath,
    SolidCoulombEventStep, SolidCoulombPhase, TraceResponse,
};
use contact::{SolidSustainedControls, SolidSustainedError};
use fs_math::dd::Dd;

fn finite(value: Dd) -> Result<f64> {
    contact::finite(value).map_err(invalid)
}
fn dot(a: &[f64], b: &[f64]) -> Dd {
    a.iter().zip(b).fold(Dd::ZERO, |sum, (&a, &b)| {
        sum + Dd::from_f64(a) * Dd::from_f64(b)
    })
}

impl SupportedSpectralSolid {
    /// Locate Coulomb stopping, static recapture and reversal while every
    /// prescribed support retains its current constant acceleration.
    ///
    /// Contact covectors and external loads use every represented coordinate.
    /// Free contact mass responses must be independent rigid directions of the
    /// original full stiffness. Cross mass, support-side contact loads, actual
    /// reactions and actuator work remain present. Other free elastic modes
    /// evolve through the original polynomial source; no modal cutoff is used.
    ///
    /// The endpoint acceleration must equal the retained acceleration. Changing
    /// it, a nonrigid contact response or a resolved extra mode change refuses
    /// without modifying either native owner. Geometry/feature admission remains
    /// mandatory. This method advances the original clock once for all phases.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_collinear_coulomb_events<P, E>(
        &mut self,
        support_acceleration_end_mps2: &[f64],
        law: &CollinearCoulombContact<'_>,
        initial_force_n: [f64; 3],
        external_force_n: &[f64],
        controls: SolidSustainedControls,
        cx: &Cx<'_>,
        admit_path: impl FnOnce(&SolidCoulombEventPath<'_, '_, Self>) -> core::result::Result<P, E>,
    ) -> core::result::Result<
        SolidCoulombEventStep<P, SupportedSolidStep>,
        SolidSustainedError<E, SupportedSolidError>,
    > {
        contact::events::advance(
            &mut Interval {
                body: self,
                acceleration_end: support_acceleration_end_mps2,
            },
            law,
            initial_force_n,
            external_force_n,
            controls,
            cx,
            admit_path,
        )
    }
}

impl EventNative for Interval<'_> {
    fn trace_response(&self, weights: [&[f64]; 3], cx: &Cx<'_>) -> Result<TraceResponse> {
        if self.acceleration_end != self.body.support_acceleration {
            return Err(invalid(
                "friction events require unchanged constant support acceleration",
            ));
        }
        let body = &self.body;
        let n = body.coordinates.len();
        let free_weights: [Vec<_>; 3] =
            core::array::from_fn(|a| body.free_rows.iter().map(|&row| weights[a][row]).collect());
        let free = body
            .free
            .trace_response([&free_weights[0], &free_weights[1], &free_weights[2]], cx)?;
        let mut columns: [Vec<f64>; 3] = core::array::from_fn(|_| vec![0.0; n]);
        let mut reactions: [Vec<f64>; 3] = core::array::from_fn(|_| vec![0.0; n]);
        let mut acceleration_offset = [0.0; 3];
        for a in 0..3 {
            for (f, &row) in body.free_rows.iter().enumerate() {
                columns[a][row] = free.columns[a][f];
            }
            let mut acceleration = Dd::ZERO;
            for (c, &row) in body.support_rows.iter().enumerate() {
                poll(cx)?;
                // Full equilibrium: M B = J^T + R, with B_c=0.
                reactions[a][row] = finite(
                    dot(&body.mass[row * n..(row + 1) * n], &columns[a])
                        - Dd::from_f64(weights[a][row]),
                )?;
                acceleration = acceleration
                    - Dd::from_f64(reactions[a][row]) * Dd::from_f64(body.support_acceleration[c]);
            }
            acceleration_offset[a] = finite(acceleration)?;
        }
        Ok(TraceResponse {
            columns,
            modal_columns: free.modal_columns,
            support_reaction_columns: Some(reactions),
            acceleration_offset,
        })
    }

    fn add_contact_moments(
        &self,
        scratch: &mut Workspace,
        response: &TraceResponse,
        path: &ContactForcePath<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<f64> {
        let body = &self.body;
        let trial = scratch
            .trial
            .take()
            .ok_or_else(|| invalid("original supported interval required"))?;
        let reactions = response
            .support_reaction_columns
            .as_ref()
            .ok_or_else(|| invalid("original supported contact reaction required"))?;
        let modal_columns = response
            .modal_columns
            .as_ref()
            .ok_or_else(|| invalid("original free modal contact response required"))?;
        let n = body.coordinates.len();
        let (mut free_contact, mut support_correction) = (Dd::ZERO, Dd::ZERO);
        let (mut surface, mut normal, mut friction) = (Dd::ZERO, Dd::ZERO, Dd::ZERO);
        for phase in path.phases {
            poll(cx)?;
            let dt = phase.duration;
            let midpoint = Dd::from_f64(phase.start_offset_s) + Dd::from_f64(0.5) * dt;
            for a in 0..3 {
                let mut supported_trace = Dd::ZERO;
                let mut reaction_work = Dd::ZERO;
                for (c, &row) in body.support_rows.iter().enumerate() {
                    let displacement = (Dd::from_f64(body.velocity[row])
                        + Dd::from_f64(body.support_acceleration[c]) * midpoint)
                        * dt;
                    supported_trace = supported_trace
                        + Dd::from_f64(path.contact.velocity_weights[a][row]) * displacement;
                    reaction_work = reaction_work + Dd::from_f64(reactions[a][row]) * displacement;
                }
                let force = Dd::from_f64(phase.force_n[a]);
                let world_trace = Dd::from_f64(phase.relative_displacement_m[a])
                    + Dd::from_f64(path.contact.prescribed_surface_velocity_mps[a]) * dt;
                free_contact = free_contact + force * (world_trace - supported_trace);
                support_correction = support_correction + force * reaction_work;
            }
            surface = surface
                + dt * dot(
                    &phase.force_n,
                    &path.contact.prescribed_surface_velocity_mps,
                );
            normal = normal
                + Dd::from_f64(phase.force_n[0]) * Dd::from_f64(phase.relative_displacement_m[0]);
            friction = friction + dot(&phase.force_n[1..], &phase.relative_displacement_m[1..]);
        }
        let mut free_source_correction = Dd::ZERO;
        for a in 0..3 {
            // B^T K_fc=0 was admitted against the original full stiffness.
            // Constant support acceleration makes this effective source constant
            // along B even though all spectator sources may remain quadratic.
            let mut effective = dot(&response.columns[a], path.external);
            for (c, &row) in body.support_rows.iter().enumerate() {
                let mass_coupling = dot(&body.mass[row * n..(row + 1) * n], &response.columns[a]);
                effective = effective - mass_coupling * Dd::from_f64(body.support_acceleration[c]);
            }
            free_source_correction =
                free_source_correction + effective * path.displacement_moment[a];
        }
        let free_supplied = finite(
            Dd::from_f64(trial.free.motion().energy.supplied_work_j)
                + free_source_correction
                + free_contact,
        )?;
        let support_work =
            finite(Dd::from_f64(trial.output.energy.support_work_j) + support_correction)?;
        let last = path
            .phases
            .last()
            .ok_or_else(|| invalid("actual friction phase required"))?;
        let free = body.free.correct_rigid_polynomial(
            trial.free,
            RigidForceCorrection {
                modal_columns,
                impulse: path.impulse,
                displacement: path.displacement_moment,
                displacement_moments: displacement_moments(body.dt_s, path.phases)?,
                endpoint_force: last.force_n,
                supplied_work_j: free_supplied,
            },
            cx,
        )?;
        let (mut displacement, mut velocity) = (trial.displacement, trial.velocity);
        let mut acceleration = vec![0.0; n];
        for (f, &row) in body.free_rows.iter().enumerate() {
            displacement[row] = free.displacement_m()[f];
            velocity[row] = free.velocity_mps()[f];
            acceleration[row] = free.acceleration_mps2()[f];
        }
        for (c, &row) in body.support_rows.iter().enumerate() {
            acceleration[row] = body.support_acceleration[c];
        }
        let mut residual = vec![0.0; n];
        let mut external_work = Dd::ZERO;
        for i in 0..n {
            poll(cx)?;
            let mut force = Dd::from_f64(path.external[i]);
            for a in 0..3 {
                force = force
                    + Dd::from_f64(path.contact.velocity_weights[a][i])
                        * Dd::from_f64(last.force_n[a]);
            }
            residual[i] = finite(
                dot(&body.mass[i * n..(i + 1) * n], &acceleration)
                    + dot(&body.stiffness[i * n..(i + 1) * n], &displacement)
                    - force,
            )?;
            external_work = external_work
                + Dd::from_f64(path.external[i])
                    * Dd::from_f64(displacement[i] - body.displacement[i]);
        }
        let free_residual = body
            .free_rows
            .iter()
            .fold(0.0_f64, |m, &row| m.max(residual[row].abs()));
        let reaction = body.support_rows.iter().map(|&row| residual[row]).collect();
        let applied_work = finite(external_work + surface + normal + friction)?;
        scratch.trial = Some(body.trial_force_path(
            free,
            displacement,
            velocity,
            trial.support_acceleration,
            reaction,
            free_residual,
            applied_work,
            support_work,
            cx,
        )?);
        Ok(support_work)
    }

    fn accept_force_path(
        &mut self,
        scratch: Workspace,
        supplied_j: f64,
        viscosity_j: f64,
        cx: &Cx<'_>,
    ) -> Result<SupportedSolidStep> {
        let trial = scratch
            .trial
            .ok_or_else(|| invalid("admitted supported event path required"))?;
        let defect = supplied_j - trial.output.energy.applied_work_j;
        if viscosity_j != 0.0 || !defect.is_finite() {
            return Err(invalid("conservative supported event work required"));
        }
        if defect.abs() > self.body.controls.maximum_work_error_j {
            return Err(SupportedSolidError::Energy {
                defect_j: defect,
                allowance_j: self.body.controls.maximum_work_error_j,
            });
        }
        poll(cx)?;
        Ok(self.body.accept_affine_supports(trial))
    }
}

// Normalized observations I_k=integral_0^1 s^k*q(h*s) ds of the actual
// piecewise-force correction. Integrating each force's remaining-time kernel
// avoids subtracting large endpoint powers for a short late segment.
fn displacement_moments(h: f64, phases: &[SolidCoulombPhase]) -> Result<[[Dd; 3]; 3]> {
    let h = Dd::from_f64(h);
    let mut result = [[Dd::ZERO; 3]; 3];
    for phase in phases {
        let weights = moment_weights(h, Dd::from_f64(phase.start_offset_s), phase.duration);
        for k in 0..3 {
            for a in 0..3 {
                result[k][a] = result[k][a] + Dd::from_f64(phase.force_n[a]) * weights[k];
                finite(result[k][a])?;
            }
        }
    }
    Ok(result)
}

fn moment_weights(h: Dd, start: Dd, duration: Dd) -> [Dd; 3] {
    let coefficients = [[1.0, 0.0, 0.0], [1.0, -1.0, 0.0], [1.0, -2.0, 1.0]];
    let duration = duration / h;
    let start = (h - start) / h;
    let end = start - duration;
    core::array::from_fn(|k| {
        let mut kernel = Dd::ZERO;
        for j in 0..=k {
            let power = j + 3;
            let mut difference = Dd::ZERO;
            for i in 0..power {
                let mut term = duration;
                for _ in 0..power - 1 - i {
                    term = term * start;
                }
                for _ in 0..i {
                    term = term * end;
                }
                difference = difference + term;
            }
            kernel = kernel
                + Dd::from_f64(coefficients[k][j]) * difference
                    / Dd::from_f64(((j + 1) * (j + 2) * (j + 3)) as f64);
        }
        kernel * h * h
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn original_rigid_force_observations_match_piecewise_polynomial_motion() {
        // h=2, acceleration 3 until t=1 and -1 afterwards. Direct original
        // q(t)=1.5*t², then q(t)=-2+4*t-0.5*t² gives these three integrals.
        let first = moment_weights(Dd::from_f64(2.0), Dd::ZERO, Dd::from_f64(1.0));
        let second = moment_weights(Dd::from_f64(2.0), Dd::from_f64(1.0), Dd::from_f64(1.0));
        for (k, expected) in [5.0 / 3.0, 29.0 / 24.0, 113.0 / 120.0].iter().enumerate() {
            let actual = finite(Dd::from_f64(3.0) * first[k] - second[k]).unwrap();
            assert!(
                (actual - expected).abs() < 1e-14,
                "moment {k}: {actual} != {expected}"
            );
        }
        // A very short late phase retains its actual time difference; a force
        // partition must not change any complete continuous observation.
        let start = Dd::from_f64(1.0 - 2.0_f64.powi(-40));
        let whole = moment_weights(Dd::from_f64(1.0), Dd::ZERO, Dd::from_f64(1.0));
        let first = moment_weights(Dd::from_f64(1.0), Dd::ZERO, start);
        let last = moment_weights(Dd::from_f64(1.0), start, Dd::from_f64(1.0) - start);
        for k in 0..3 {
            assert!(finite(first[k] + last[k] - whole[k]).unwrap().abs() < 1e-30);
        }
    }
}
