//! Three-coordinate inverse geometry, retaining both original meridian halves.
use super::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Evaluation {
    pub profile: AsymmetricCircularContactLineProfile,
    pub q: f64,
    pub alpha: [f64; 2],
    pub geometric: [f64; 3],
    /// Rows h,b0,b1; columns q,alpha0,alpha1, at the original fixed volume.
    pub jacobian: [[f64; 3]; 3],
    pub gradient: [f64; 3],
    pub hessian: [[f64; 3]; 3],
    pub gradient_error: [f64; 3],
}
#[derive(Default)]
pub(super) struct Work {
    pub shapes: usize,
    pub quadrature: u64,
}

impl SourceAsymmetricCircularContactLineBridge {
    pub(super) fn at_parameters(
        &self,
        q: f64,
        alpha: [f64; 2],
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Evaluation, ContactLineError> {
        poll(gate)?;
        let theta = alpha.map(|a| FRAC_PI_2 - a);
        if !q.is_finite()
            || q < self.original.branch.minimum_suction_parameter
            || q > self.original.branch.maximum_suction_parameter
            || (0..2).any(|i| {
                !theta[i].is_finite()
                    || theta[i] < self.config.minimum_angle_radians
                    || theta[i] > self.config.maximum_angle_radians
                    || (theta[i] - self.original.wetting()[i].contact_angle_radians()).abs()
                        > self.config.maximum_angle_departure_radians
            })
        {
            return Err(ContactLineError::Domain(
                "declared suction and independent dynamic-angle branch",
            ));
        }
        let count = if alpha[0] == alpha[1] { 1 } else { 2 };
        if work
            .shapes
            .checked_add(count)
            .is_none_or(|n| n > self.config.maximum_shape_evaluations)
        {
            return Err(ContactLineError::Budget);
        }
        work.shapes += count;
        let lower = super::super::numerics::integrate(alpha[0], q, self.config.shape, gate)
            .map_err(original_error)?;
        let (upper_values, upper_error, upper_evaluations) = if alpha[0] == alpha[1] {
            (lower.values, lower.relative_error, 0)
        } else {
            let upper = super::super::numerics::integrate(alpha[1], q, self.config.shape, gate)
                .map_err(original_error)?;
            (upper.values, upper.relative_error, upper.evaluations)
        };
        work.quadrature = work
            .quadrature
            .checked_add(lower.evaluations + upper_evaluations)
            .ok_or(ContactLineError::Budget)?;
        let sums: [f64; 6] = core::array::from_fn(|i| lower.values[i] + upper_values[i]);
        let [h, v, area, hq, vq, _] = sums;
        let relative = lower.relative_error.max(upper_error) + 64.0 * f64::EPSILON;
        let a = det::pow(self.original.volume / (PI * v), 1.0 / 3.0);
        if !a.is_finite() || a < self.config.minimum_neck_radius_m {
            return Err(ContactLineError::Domain("minimum actual axisymmetric neck"));
        }
        let mut radius = [0.0; 2];
        let mut radius_q = [0.0; 2];
        let mut radius_alpha = [0.0; 2];
        let mut height_alpha = [0.0; 2];
        let mut volume_alpha = [0.0; 2];
        let mut sa = [0.0; 2];
        let mut ca = [0.0; 2];
        for i in 0..2 {
            let (x, xq, d) = super::super::numerics::radius(alpha[i], q);
            (sa[i], ca[i]) = det::sin_cos(alpha[i]);
            radius[i] = x;
            radius_q[i] = xq;
            radius_alpha[i] = x * sa[i] / d;
            height_alpha[i] = x * ca[i] / d;
            volume_alpha[i] = x * x * height_alpha[i];
        }
        let scale_derivatives = [
            -a * vq / (3.0 * v),
            -a * volume_alpha[0] / (3.0 * v),
            -a * volume_alpha[1] / (3.0 * v),
        ];
        let gap = a * h;
        let b = radius.map(|r| a * r);
        self.check_coordinates([gap, b[0], b[1]])?;
        let mut jacobian = [[0.0; 3]; 3];
        for k in 0..3 {
            jacobian[0][k] =
                scale_derivatives[k] * h + a * if k == 0 { hq } else { height_alpha[k - 1] };
            for i in 0..2 {
                jacobian[i + 1][k] = scale_derivatives[k] * radius[i]
                    + a * if k == 0 {
                        radius_q[i]
                    } else if k == i + 1 {
                        radius_alpha[i]
                    } else {
                        0.0
                    };
            }
        }
        let inverse = invert(jacobian)?;
        if jacobian[0][0] >= 0.0 {
            return Err(ContactLineError::Domain("monotone suction branch"));
        }
        let gamma = self.original.sigma;
        let gradient = [
            PI * gamma * a * (2.0 + q),
            2.0 * PI * gamma * b[0] * (sa[0] - self.original.cosine[0]),
            2.0 * PI * gamma * b[1] * (sa[1] - self.original.cosine[1]),
        ];
        let mut parameter_derivatives = [[0.0; 3]; 3];
        for k in 0..3 {
            parameter_derivatives[0][k] =
                PI * gamma * ((2.0 + q) * scale_derivatives[k] + if k == 0 { a } else { 0.0 });
            for i in 0..2 {
                parameter_derivatives[i + 1][k] = 2.0
                    * PI
                    * gamma
                    * ((sa[i] - self.original.cosine[i]) * jacobian[i + 1][k]
                        + if k == i + 1 { b[i] * ca[i] } else { 0.0 });
            }
        }
        let mut hessian: [[f64; 3]; 3] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (0..3)
                    .map(|k| parameter_derivatives[i][k] * inverse[k][j])
                    .sum()
            })
        });
        if hessian.iter().flatten().any(|x| !x.is_finite()) {
            return Err(ContactLineError::Domain(
                "finite fixed-volume shape Hessian",
            ));
        }
        for i in 0..3 {
            for j in 0..i {
                let scale = hessian[i][j]
                    .abs()
                    .max(hessian[j][i].abs())
                    .max(f64::EPSILON * hessian[i][i].abs().max(hessian[j][j].abs()));
                if (hessian[i][j] - hessian[j][i]).abs() > (64.0 * relative).max(4e-10) * scale {
                    return Err(ContactLineError::Domain(
                        "reciprocal circular-footprint Hessian",
                    ));
                }
                let mean = 0.5 * (hessian[i][j] + hessian[j][i]);
                hessian[i][j] = mean;
                hessian[j][i] = mean;
            }
        }
        let schur = hessian[2][2] - hessian[2][1] * hessian[1][2] / hessian[1][1];
        if hessian[1][1] <= 0.0 || !schur.is_finite() || schur <= 0.0 {
            return Err(ContactLineError::Domain(
                "positive coupled circular-footprint Hessian",
            ));
        }
        let margin = (FRAC_PI_2 / h).powi(2) * ca[0].min(ca[1]).powi(3) - 1.0;
        if !margin.is_finite()
            || margin
                <= self.original.branch.minimum_stability_margin + 4.0 * relative * (margin + 1.0)
        {
            return Err(ContactLineError::Domain(
                "axisymmetric fixed-footprint stability screen",
            ));
        }
        let free_area = 2.0 * PI * a * a * area;
        let wet_area = b.map(|r| PI * r * r);
        let wet_energy = core::array::from_fn::<_, 2, _>(|i| {
            self.original.laws[i].wetted_surface_energy_difference_j_m2() * wet_area[i]
        });
        let energy = gamma * free_area + wet_energy[0] + wet_energy[1];
        let energy_scale = gamma * free_area + wet_energy[0].abs() + wet_energy[1].abs();
        let pressure = -gamma * q / a;
        let length = gap.max(2.0 * b[0].max(b[1]));
        let bond = self.original.liquid.density_kg_m3() * self.original.limits.gravity_m_s2 / gamma
            * length
            * length;
        if !energy.is_finite()
            || energy.abs() > self.original.limits.maximum_abs_surface_energy_j
            || !pressure.is_finite()
            || -pressure > self.original.limits.maximum_suction_pa
            || gradient
                .iter()
                .any(|x| !x.is_finite() || x.abs() > self.original.limits.maximum_abs_force_n)
            || !bond.is_finite()
            || bond > self.original.limits.maximum_bond_number
        {
            return Err(ContactLineError::Domain(
                "original axisymmetric force/energy/suction/gravity envelope",
            ));
        }
        let geometric = [gap, b[0], b[1]];
        let coordinate_error = geometric.map(|x| 2.0 * relative * x);
        let mut gradient_error = core::array::from_fn::<_, 3, _>(|i| {
            (0..3)
                .map(|j| hessian[i][j].abs() * coordinate_error[j])
                .sum::<f64>()
                + relative * gradient[i].abs()
        });
        for i in 0..2 {
            gradient_error[i + 1] += 64.0
                * f64::EPSILON
                * 2.0
                * PI
                * gamma
                * b[i]
                * (sa[i].abs() + self.original.cosine[i].abs());
        }
        let energy_error = 4.0 * relative * energy_scale;
        if gradient_error.iter().any(|x| !x.is_finite()) || !energy_error.is_finite() {
            return Err(ContactLineError::Domain(
                "finite asymmetric circular numerical uncertainty",
            ));
        }
        Ok(Evaluation {
            profile: AsymmetricCircularContactLineProfile {
                gap_m: gap,
                wetted_radii_m: b,
                contact_angles_radians: theta,
                suction_parameter: q,
                neck_radius_m: a,
                neck_height_m: a * lower.values[0],
                liquid_volume_m3: self.original.volume,
                liquid_gas_area_m2: free_area,
                wetted_areas_m2: wet_area,
                liquid_minus_gas_pressure_pa: pressure,
                surface_free_energy_j: energy,
                normal_force_n: -gradient[0],
                footprint_energy_gradients_n: [gradient[1], gradient[2]],
                energy_hessian_n_per_m: hessian,
                axisymmetric_stability_margin: margin,
                geometric_gap_residual_m: 0.0,
                geometric_radius_residuals_m: [0.0; 2],
                quadrature_relative_error_estimate: lower.relative_error.max(upper_error),
                estimated_volume_error_m3: 4.0 * relative * self.original.volume,
                estimated_surface_energy_error_j: energy_error,
            },
            q,
            alpha,
            geometric,
            jacobian,
            gradient,
            hessian,
            gradient_error,
        })
    }

    pub(super) fn at_gap_angles(
        &self,
        gap: f64,
        theta: [f64; 2],
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Evaluation, ContactLineError> {
        self.original.check_gap(gap).map_err(original_error)?;
        let mut lower = self.original.branch.minimum_suction_parameter;
        let mut upper = self.original.branch.maximum_suction_parameter;
        let mut q = 0.5 * (lower + upper);
        for _ in 0..self.config.shape.maximum_root_iterations {
            let e = self.at_parameters(q, theta.map(|t| FRAC_PI_2 - t), work, gate)?;
            let residual = e.geometric[0] - gap;
            if residual.abs() <= self.coordinate_tolerance(gap) {
                return Ok(e.with_coordinates([gap, e.geometric[1], e.geometric[2]]));
            }
            let next = q - residual / e.jacobian[0][0];
            if residual > 0.0 {
                lower = q;
            } else {
                upper = q;
            }
            let width = upper - lower;
            q = if next > lower + 0.01 * width && next < upper - 0.01 * width {
                next
            } else {
                0.5 * (lower + upper)
            };
        }
        Err(ContactLineError::NotConverged)
    }
    pub(super) fn at(
        &self,
        target: [f64; 3],
        seed: Evaluation,
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Evaluation, ContactLineError> {
        self.check_coordinates(target)?;
        let mut e = seed;
        for _ in 0..self.config.shape.maximum_root_iterations {
            poll(gate)?;
            let residual = core::array::from_fn::<_, 3, _>(|i| e.geometric[i] - target[i]);
            let norm = self.geometry_norm(residual, target);
            if norm <= 1.0 {
                return Ok(e.with_coordinates(target));
            }
            let inverse = invert(e.jacobian)?;
            let delta = core::array::from_fn::<_, 3, _>(|i| {
                (0..3).map(|j| inverse[i][j] * residual[j]).sum::<f64>()
            });
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..self.config.maximum_backtracks {
                let q = e.q - fraction * delta[0];
                let alpha = core::array::from_fn(|i| e.alpha[i] - fraction * delta[i + 1]);
                match self.at_parameters(q, alpha, work, gate) {
                    Ok(next) => {
                        if self.geometry_norm(
                            core::array::from_fn(|i| next.geometric[i] - target[i]),
                            target,
                        ) < norm
                        {
                            accepted = Some(next);
                            break;
                        }
                    }
                    Err(ContactLineError::Domain(_)) => {}
                    Err(error) => return Err(error),
                }
                fraction *= 0.5;
            }
            e = accepted.ok_or(ContactLineError::Domain(
                "joint circular-footprint inverse continuation",
            ))?;
        }
        Err(ContactLineError::NotConverged)
    }
    fn geometry_norm(&self, residual: [f64; 3], target: [f64; 3]) -> f64 {
        (0..3)
            .map(|i| residual[i].abs() / self.coordinate_tolerance(target[i]))
            .fold(0.0, f64::max)
    }
    fn coordinate_tolerance(&self, x: f64) -> f64 {
        self.config.geometry_absolute_tolerance_m
            + self.config.shape.gap_relative_tolerance * x.abs()
    }
    pub(super) fn check_coordinates(&self, x: [f64; 3]) -> Result<(), ContactLineError> {
        if !x[0].is_finite()
            || x[0] < self.original.limits.minimum_gap_m
            || x[0] > self.original.limits.maximum_gap_m
            || x[1..].iter().any(|b| {
                !b.is_finite() || *b <= 0.0 || *b > self.original.limits.maximum_wetted_radius_m
            })
        {
            Err(ContactLineError::Domain(
                "original gap and independent circular footprint envelope",
            ))
        } else {
            Ok(())
        }
    }
}
impl Evaluation {
    fn with_coordinates(mut self, x: [f64; 3]) -> Self {
        let residual = core::array::from_fn::<_, 3, _>(|i| self.geometric[i] - x[i]);
        let old = [
            self.profile.geometric_gap_residual_m,
            self.profile.geometric_radius_residuals_m[0],
            self.profile.geometric_radius_residuals_m[1],
        ];
        let extra =
            core::array::from_fn::<_, 3, _>(|i| (residual[i].abs() - old[i].abs()).max(0.0));
        for i in 0..3 {
            self.gradient_error[i] += (0..3)
                .map(|j| self.hessian[i][j].abs() * extra[j])
                .sum::<f64>();
            self.profile.estimated_surface_energy_error_j += self.gradient[i].abs() * extra[i];
        }
        self.profile.gap_m = x[0];
        self.profile.wetted_radii_m = [x[1], x[2]];
        self.profile.geometric_gap_residual_m = residual[0];
        self.profile.geometric_radius_residuals_m = [residual[1], residual[2]];
        self
    }
}

fn invert(j: [[f64; 3]; 3]) -> Result<[[f64; 3]; 3], ContactLineError> {
    let cofactor = core::array::from_fn::<_, 3, _>(|i| {
        core::array::from_fn::<_, 3, _>(|k| {
            let a = (i + 1) % 3;
            let b = (i + 2) % 3;
            let c = (k + 1) % 3;
            let d = (k + 2) % 3;
            j[a][c] * j[b][d] - j[a][d] * j[b][c]
        })
    });
    let determinant = (0..3).map(|k| j[0][k] * cofactor[0][k]).sum::<f64>();
    let scale = (0..3)
        .map(|k| {
            j[0][k].abs()
                * ((j[1][(k + 1) % 3] * j[2][(k + 2) % 3]).abs()
                    + (j[1][(k + 2) % 3] * j[2][(k + 1) % 3]).abs())
        })
        .sum::<f64>();
    if j.iter().flatten().any(|x| !x.is_finite())
        || !determinant.is_finite()
        || determinant.abs() <= 128.0 * f64::EPSILON * scale
    {
        return Err(ContactLineError::Domain(
            "nonsingular three-coordinate fixed-volume shape",
        ));
    }
    Ok(core::array::from_fn(|i| {
        core::array::from_fn(|k| cofactor[k][i] / determinant)
    }))
}
