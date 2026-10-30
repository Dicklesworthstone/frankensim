//! Coupled pinning/sliding solve and the original complete surface-work audit.
use super::integration::Gradient;
use super::shape::Work;
use super::*;

#[derive(Clone, Copy)]
struct Selected {
    regimes: [ContactLineRegime; 2],
    forces: [f64; 2],
    jacobian: [[f64; 2]; 2],
    friction: [f64; 2],
    pinning: [f64; 2],
}
impl SourceAsymmetricCircularContactLineBridge {
    pub fn advance(
        &self,
        previous: AsymmetricCircularContactLineState,
        gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<AsymmetricCircularContactLineStep, ContactLineError> {
        poll(gate)?;
        let p = self.evaluate(previous)?;
        let c = self.config;
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > c.maximum_step_s {
            return Err(ContactLineError::Domain(
                "positive admitted physical timestep",
            ));
        }
        if !gap_m.is_finite() || (gap_m - p.gap_m).abs() > c.maximum_gap_increment_m {
            return Err(ContactLineError::Domain("physical gap increment"));
        }
        let b0 = p.wetted_radii_m;
        let old = [p.gap_m, b0[0], b0[1]];
        self.check_coordinates([gap_m, b0[0], b0[1]])?;
        let coefficient = self.zeta.map(|z| 2.0 * PI * z / dt_s);
        if coefficient.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(ContactLineError::Domain(
                "finite positive circular-line resistance",
            ));
        }
        let mut work = Work::default();
        let mut u = [0.0; 2];
        let mut gradient = self.discrete_gradient(
            old,
            [gap_m, b0[0], b0[1]],
            previous.evaluation,
            &mut work,
            gate,
        )?;
        let mut selected = self.residual(b0, u, coefficient, gradient)?;
        let mut iterations = 0;
        for iteration in 0..=c.maximum_iterations {
            iterations = iteration;
            let scale = gradient.value[1]
                .abs()
                .max(gradient.value[2].abs())
                .max((selected.friction[0] * u[0]).abs())
                .max((selected.friction[1] * u[1]).abs());
            let tolerance = self.force_tolerance(scale);
            let pinned_exact =
                (0..2).all(|i| selected.regimes[i] != ContactLineRegime::Pinned || u[i] == 0.0);
            let directions = (0..2).all(|i| match selected.regimes[i] {
                ContactLineRegime::Pinned => true,
                ContactLineRegime::Advancing => u[i] >= 0.0,
                ContactLineRegime::Receding => u[i] <= 0.0,
            });
            let error = gradient.error[1].max(gradient.error[2]);
            if maximum(selected.forces) + error <= tolerance && pinned_exact && directions {
                break;
            }
            if iteration == c.maximum_iterations {
                return Err(ContactLineError::NotConverged);
            }
            let correction = solve(selected.jacobian, selected.forces.map(|x| -x))?;
            let merit = maximum(selected.forces);
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..c.maximum_backtracks {
                poll(gate)?;
                let candidate = core::array::from_fn::<_, 2, _>(|i| {
                    if fraction == 1.0 && selected.regimes[i] == ContactLineRegime::Pinned {
                        0.0
                    } else {
                        u[i] + fraction * correction[i]
                    }
                });
                if candidate
                    .iter()
                    .any(|x| !x.is_finite() || x.abs() > c.maximum_radius_increment_m)
                {
                    fraction *= 0.5;
                    continue;
                }
                let new = [gap_m, b0[0] + candidate[0], b0[1] + candidate[1]];
                match self.discrete_gradient(old, new, previous.evaluation, &mut work, gate) {
                    Ok(next) => {
                        let next_selected = self.residual(b0, candidate, coefficient, next)?;
                        let next_merit = maximum(next_selected.forces);
                        if next_merit < merit
                            || next_merit + next.error[1].max(next.error[2]) <= tolerance
                        {
                            accepted = Some((candidate, next, next_selected));
                            break;
                        }
                    }
                    Err(ContactLineError::Domain(_)) => {}
                    Err(error) => return Err(error),
                }
                fraction *= 0.5;
            }
            let (next_u, next_gradient, next_selected) =
                accepted.ok_or(ContactLineError::NotConverged)?;
            if next_u == u {
                return Err(ContactLineError::NotConverged);
            }
            u = next_u;
            gradient = next_gradient;
            selected = next_selected;
        }
        let b = core::array::from_fn::<_, 2, _>(|i| b0[i] + u[i]);
        let actual_u = core::array::from_fn(|i| b[i] - b0[i]);
        if actual_u != u {
            u = actual_u;
            gradient = self.discrete_gradient(
                old,
                [gap_m, b[0], b[1]],
                previous.evaluation,
                &mut work,
                gate,
            )?;
            selected = self.residual(b0, u, coefficient, gradient)?;
            if maximum(selected.forces) + gradient.error[1].max(gradient.error[2])
                > self.force_tolerance(gradient.value[1].abs().max(gradient.value[2].abs()))
            {
                return Err(ContactLineError::NotConverged);
            }
        }
        let endpoint = self.at([gap_m, b[0], b[1]], previous.evaluation, &mut work, gate)?;
        let response = solve(
            selected.jacobian,
            core::array::from_fn(|i| {
                if selected.regimes[i] == ContactLineRegime::Pinned {
                    0.0
                } else {
                    -gradient.derivative[i + 1][0]
                }
            }),
        )?;
        let mean_force = -gradient.value[0];
        let tangent = -(gradient.derivative[0][0]
            + gradient.derivative[0][1] * response[0]
            + gradient.derivative[0][2] * response[1]);
        // Absolute inverse rows propagate both constitutive residual and shape
        // estimates through the actually selected coupled implicit resistance.
        let column0 = solve(selected.jacobian, [1.0, 0.0])?;
        let column1 = solve(selected.jacobian, [0.0, 1.0])?;
        let force_errors =
            core::array::from_fn::<_, 2, _>(|i| selected.forces[i].abs() + gradient.error[i + 1]);
        let radius_errors = core::array::from_fn::<_, 2, _>(|i| {
            column0[i].abs() * force_errors[0] + column1[i].abs() * force_errors[1]
        });
        let force_error = gradient.error[0]
            + gradient.derivative[0][1].abs() * radius_errors[0]
            + gradient.derivative[0][2].abs() * radius_errors[1];
        let force_tolerance = self.force_tolerance(mean_force.abs());
        if !mean_force.is_finite()
            || !tangent.is_finite()
            || !force_error.is_finite()
            || force_error > force_tolerance
        {
            return Err(ContactLineError::Accuracy {
                estimated_error_n: force_error,
                tolerance_n: force_tolerance,
            });
        }
        let dh = gap_m - p.gap_m;
        let speed = u.map(|x| (x / dt_s).abs());
        let flow_speed = speed[0].max(speed[1]).max((dh / dt_s).abs());
        let capillary =
            self.original.liquid.dynamic_viscosity_pa_s() * flow_speed / self.original.sigma;
        let length = p
            .gap_m
            .max(gap_m)
            .max(2.0 * b0[0].max(b0[1]).max(b[0]).max(b[1]));
        let weber = self.original.liquid.density_kg_m3() * flow_speed * flow_speed * length
            / self.original.sigma;
        if speed
            .iter()
            .any(|s| !s.is_finite() || *s > c.maximum_contact_line_speed_m_s)
            || !capillary.is_finite()
            || capillary > c.maximum_capillary_number
            || !weber.is_finite()
            || weber > c.maximum_weber_number
        {
            return Err(ContactLineError::Domain(
                "line-friction dominated low-speed circular wetting",
            ));
        }
        let viscous = core::array::from_fn::<_, 2, _>(|i| selected.friction[i] * u[i] * u[i]);
        let pinning = core::array::from_fn::<_, 2, _>(|i| selected.pinning[i] * u[i]);
        let viscous_total = viscous[0] + viscous[1];
        let pinning_total = pinning[0] + pinning[1];
        let loss = viscous_total + pinning_total;
        let energy_change = endpoint.profile.surface_free_energy_j - p.surface_free_energy_j;
        let liquid_work = gradient.value[0] * dh;
        let residual = energy_change + loss - liquid_work;
        let energy_error = p.estimated_surface_energy_error_j
            + endpoint.profile.estimated_surface_energy_error_j
            + gradient.error[0] * dh.abs()
            + gradient.error[1] * u[0].abs()
            + gradient.error[2] * u[1].abs()
            + 64.0 * f64::EPSILON * (energy_change.abs() + loss.abs() + liquid_work.abs());
        let scale = p
            .surface_free_energy_j
            .abs()
            .max(endpoint.profile.surface_free_energy_j.abs())
            .max(liquid_work.abs())
            .max(loss);
        let tolerance = c.energy_absolute_tolerance_j + c.energy_relative_tolerance * scale;
        if !loss.is_finite()
            || viscous
                .iter()
                .chain(pinning.iter())
                .any(|x| *x < 0.0 || !x.is_finite())
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
        if mean_force.abs() > self.original.limits.maximum_abs_force_n {
            return Err(ContactLineError::Domain(
                "original mean normal force envelope",
            ));
        }
        poll(gate)?;
        Ok(AsymmetricCircularContactLineStep {
            state: self.state_from(endpoint),
            profile: endpoint.profile,
            regimes: selected.regimes,
            mean_normal_force_n: mean_force,
            mean_force_gap_derivative_n_per_m: tangent,
            surface_energy_change_j: energy_change,
            contact_line_dissipation_j: loss,
            viscous_line_dissipation_j: viscous_total,
            pinning_dissipation_j: pinning_total,
            per_line_viscous_dissipation_j: viscous,
            per_line_pinning_dissipation_j: pinning,
            work_on_liquid_j: liquid_work,
            energy_balance_residual_j: residual,
            footprint_force_residuals_n: selected.forces,
            estimated_force_error_n: force_error,
            estimated_energy_error_j: energy_error,
            iterations,
            shape_evaluations: work.shapes,
            quadrature_evaluations: work.quadrature,
        })
    }
    fn residual(
        &self,
        b0: [f64; 2],
        u: [f64; 2],
        coefficient: [f64; 2],
        g: Gradient,
    ) -> Result<Selected, ContactLineError> {
        let mut selected = Selected {
            regimes: [ContactLineRegime::Pinned; 2],
            forces: [0.0; 2],
            jacobian: [[0.0; 2]; 2],
            friction: [0.0; 2],
            pinning: [0.0; 2],
        };
        for i in 0..2 {
            let mean = b0[i] + 0.5 * u[i];
            selected.friction[i] = coefficient[i] * mean;
            let advancing = self.advancing_per_radius_n_m[i] * mean;
            let receding = self.receding_per_radius_n_m[i] * mean;
            let threshold = if g.value[i + 1] < -advancing {
                selected.regimes[i] = ContactLineRegime::Advancing;
                self.advancing_per_radius_n_m[i]
            } else if g.value[i + 1] > receding {
                selected.regimes[i] = ContactLineRegime::Receding;
                -self.receding_per_radius_n_m[i]
            } else {
                0.0
            };
            selected.jacobian[i][i] = coefficient[i] * (b0[i] + u[i]);
            selected.forces[i] = selected.friction[i] * u[i];
            if selected.regimes[i] == ContactLineRegime::Pinned {
                selected.pinning[i] = -g.value[i + 1];
            } else {
                selected.pinning[i] = threshold * mean;
                selected.forces[i] += g.value[i + 1] + selected.pinning[i];
                selected.jacobian[i][i] += 0.5 * threshold;
                for j in 0..2 {
                    selected.jacobian[i][j] += g.derivative[i + 1][j + 1];
                }
            }
            if !mean.is_finite()
                || mean <= 0.0
                || !selected.forces[i].is_finite()
                || !selected.friction[i].is_finite()
            {
                return Err(ContactLineError::Domain(
                    "finite circular-line physical resistance",
                ));
            }
        }
        // Validate the selected local response before it can enter a tangent
        // or error estimate; this is a local branch, not a global uniqueness claim.
        solve(selected.jacobian, [0.0; 2])?;
        Ok(selected)
    }
    fn force_tolerance(&self, scale: f64) -> f64 {
        self.config.force_absolute_tolerance_n + self.config.force_relative_tolerance * scale
    }
}
fn maximum(x: [f64; 2]) -> f64 {
    x[0].abs().max(x[1].abs())
}
fn solve(a: [[f64; 2]; 2], rhs: [f64; 2]) -> Result<[f64; 2], ContactLineError> {
    let schur = a[1][1] - a[1][0] * (a[0][1] / a[0][0]);
    if a.iter().flatten().any(|x| !x.is_finite())
        || a[0][0] <= 0.0
        || !schur.is_finite()
        || schur <= 0.0
    {
        return Err(ContactLineError::Domain(
            "positive coupled circular-line implicit resistance",
        ));
    }
    let y = (rhs[1] - a[1][0] * (rhs[0] / a[0][0])) / schur;
    let x = (rhs[0] - a[0][1] * y) / a[0][0];
    if !x.is_finite() || !y.is_finite() {
        return Err(ContactLineError::NotConverged);
    }
    Ok([x, y])
}
