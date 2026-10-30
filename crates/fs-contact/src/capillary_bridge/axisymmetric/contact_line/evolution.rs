//! Bounded implicit physical-radius evolution and original surface-work audit.
use super::integration::Gradient;
use super::shape::Work;
use super::*;

#[derive(Clone, Copy)]
struct Residual {
    gradient: Gradient,
    value: f64,
    derivative: f64,
    friction: f64,
    pinning: f64,
}

impl SourceAxisymmetricContactLineBridge {
    pub fn advance(
        &self,
        previous: AxisymmetricContactLineState,
        gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<AxisymmetricContactLineStep, ContactLineError> {
        poll(gate)?;
        self.require_state(previous)?;
        let p = previous.evaluation.profile;
        let old = [p.gap_m, p.wetted_radius_m];
        if !dt_s.is_finite()
            || dt_s <= 0.0
            || dt_s > self.config.maximum_step_s
            || !gap_m.is_finite()
            || (gap_m - p.gap_m).abs() > self.config.maximum_gap_increment_m
        {
            return Err(ContactLineError::Domain(
                "finite circular-line timestep and gap increment",
            ));
        }
        self.check_coordinates(gap_m, p.wetted_radius_m)?;
        let mut work = Work::default();
        let initial = self.discrete_gradient(
            1,
            old,
            [gap_m, old[1]],
            previous.evaluation,
            &mut work,
            gate,
        )?;
        let (regime, threshold) = if initial.value < -self.advancing_per_radius_n_m * old[1] {
            (ContactLineRegime::Advancing, self.advancing_per_radius_n_m)
        } else if initial.value > self.receding_per_radius_n_m * old[1] {
            (ContactLineRegime::Receding, -self.receding_per_radius_n_m)
        } else {
            (ContactLineRegime::Pinned, 0.0)
        };
        let mut b = old[1];
        let mut iterations = 0;
        let mut selected = self.residual(initial, b, old[1], dt_s, threshold)?;
        if regime != ContactLineRegime::Pinned {
            let direction = if regime == ContactLineRegime::Advancing {
                1.0
            } else {
                -1.0
            };
            let mut converged = false;
            for iteration in 0..self.config.maximum_iterations {
                poll(gate)?;
                iterations = iteration;
                let tolerance = self.force_tolerance(selected);
                if selected.value.abs() <= tolerance {
                    converged = true;
                    break;
                }
                let increment = -selected.value / selected.derivative;
                let mut fraction = 1.0;
                let mut accepted = None;
                for _ in 0..self.config.maximum_backtracks {
                    let next = b + fraction * increment;
                    let delta = next - old[1];
                    if next > 0.0
                        && delta * direction >= 0.0
                        && delta.abs() <= self.config.maximum_radius_increment_m
                    {
                        match self.discrete_gradient(
                            1,
                            old,
                            [gap_m, next],
                            previous.evaluation,
                            &mut work,
                            gate,
                        ) {
                            Ok(gradient) => {
                                let trial =
                                    self.residual(gradient, next, old[1], dt_s, threshold)?;
                                if trial.value.abs() < selected.value.abs() {
                                    accepted = Some((next, trial));
                                    break;
                                }
                            }
                            Err(ContactLineError::Domain(_)) => {}
                            Err(error) => return Err(error),
                        }
                    }
                    fraction *= 0.5;
                }
                (b, selected) = accepted.ok_or(ContactLineError::Domain(
                    "circular-line local evolution branch",
                ))?;
            }
            if !converged && selected.value.abs() > self.force_tolerance(selected) {
                return Err(ContactLineError::NotConverged);
            }
        }
        let endpoint = self.at(gap_m, b, previous.evaluation, &mut work, gate)?;
        let normal =
            self.discrete_gradient(0, old, [gap_m, b], previous.evaluation, &mut work, gate)?;
        let db = b - old[1];
        let dh = gap_m - old[0];
        let speed = (db / dt_s).abs();
        let flow_speed = speed.max((dh / dt_s).abs());
        let viscosity = self.original.liquid.dynamic_viscosity_pa_s();
        let capillary = viscosity * flow_speed / self.original.sigma;
        let length = p.gap_m.max(gap_m).max(2.0 * p.wetted_radius_m.max(b));
        let weber = self.original.liquid.density_kg_m3() * flow_speed * flow_speed * length
            / self.original.sigma;
        if !speed.is_finite()
            || speed > self.config.maximum_contact_line_speed_m_s
            || !capillary.is_finite()
            || capillary > self.config.maximum_capillary_number
            || !weber.is_finite()
            || weber > self.config.maximum_weber_number
        {
            return Err(ContactLineError::Domain(
                "line-friction dominated low-speed circular wetting",
            ));
        }
        let viscous = selected.friction * db * db;
        let pinning = selected.pinning * db;
        let dissipation = viscous + pinning;
        let energy_change = endpoint.profile.surface_free_energy_j - p.surface_free_energy_j;
        let liquid_work = normal.value * dh;
        let residual = energy_change + dissipation - liquid_work;
        let arithmetic_scale = energy_change.abs() + dissipation.abs() + liquid_work.abs();
        let energy_error = p.estimated_surface_energy_error_j
            + endpoint.profile.estimated_surface_energy_error_j
            + normal.error * dh.abs()
            + selected.gradient.error * db.abs()
            + 64.0 * f64::EPSILON * arithmetic_scale;
        let scale = p
            .surface_free_energy_j
            .abs()
            .max(endpoint.profile.surface_free_energy_j.abs())
            .max(liquid_work.abs())
            .max(dissipation);
        let tolerance =
            self.config.energy_absolute_tolerance_j + self.config.energy_relative_tolerance * scale;
        if !dissipation.is_finite()
            || viscous < 0.0
            || pinning < 0.0
            || !liquid_work.is_finite()
            || liquid_work.abs() > self.original.limits.maximum_abs_step_work_j
            || !energy_error.is_finite()
            || residual.abs() + energy_error > tolerance
        {
            return Err(ContactLineError::Energy {
                residual_j: residual,
                tolerance_j: tolerance,
            });
        }
        let radius_gap_derivative = if regime == ContactLineRegime::Pinned {
            0.0
        } else {
            -selected.gradient.derivative[0] / selected.derivative
        };
        let tangent = -(normal.derivative[0] + normal.derivative[1] * radius_gap_derivative);
        let footprint_error = if regime == ContactLineRegime::Pinned {
            0.0
        } else {
            (selected.value.abs() + selected.gradient.error) / selected.derivative
        };
        let force_error = normal.error + normal.derivative[1].abs() * footprint_error;
        let force_tolerance = self.config.force_absolute_tolerance_n
            + self.config.force_relative_tolerance * normal.value.abs();
        if !force_error.is_finite() || !tangent.is_finite() || force_error > force_tolerance {
            return Err(ContactLineError::Accuracy {
                estimated_error_n: force_error,
                tolerance_n: force_tolerance,
            });
        }
        if normal.value.abs() > self.original.limits.maximum_abs_force_n {
            return Err(ContactLineError::Domain(
                "original mean normal force envelope",
            ));
        }
        poll(gate)?;
        Ok(AxisymmetricContactLineStep {
            state: self.state_from(endpoint),
            profile: endpoint.profile,
            regime,
            mean_normal_force_n: -normal.value,
            mean_force_gap_derivative_n_per_m: tangent,
            surface_energy_change_j: energy_change,
            contact_line_dissipation_j: dissipation,
            viscous_line_dissipation_j: viscous,
            pinning_dissipation_j: pinning,
            work_on_liquid_j: liquid_work,
            energy_balance_residual_j: residual,
            footprint_force_residual_n: if regime == ContactLineRegime::Pinned {
                0.0
            } else {
                selected.value
            },
            estimated_force_error_n: force_error,
            estimated_energy_error_j: energy_error,
            iterations,
            shape_evaluations: work.shapes,
            quadrature_evaluations: work.quadrature,
        })
    }

    fn residual(
        &self,
        gradient: Gradient,
        b: f64,
        old_b: f64,
        dt: f64,
        threshold: f64,
    ) -> Result<Residual, ContactLineError> {
        let mean_radius = 0.5 * (old_b + b);
        // Integral of the ACTUAL two circular perimeters along the linear
        // radius path. Neither line length nor pinning work is frozen at b0.
        let coefficient = 4.0 * PI * self.zeta / dt;
        let friction = coefficient * mean_radius;
        let pinning = threshold * mean_radius;
        let value = gradient.value + friction * (b - old_b) + pinning;
        let derivative =
            gradient.derivative[1] + friction + 0.5 * coefficient * (b - old_b) + 0.5 * threshold;
        if !derivative.is_finite() || derivative <= 0.0 || !value.is_finite() {
            return Err(ContactLineError::Domain(
                "positive circular-line implicit resistance",
            ));
        }
        Ok(Residual {
            gradient,
            value,
            derivative,
            friction,
            pinning,
        })
    }
    fn force_tolerance(&self, residual: Residual) -> f64 {
        self.config.force_absolute_tolerance_n
            + self.config.force_relative_tolerance
                * residual.gradient.value.abs().max(residual.pinning.abs())
    }
}
