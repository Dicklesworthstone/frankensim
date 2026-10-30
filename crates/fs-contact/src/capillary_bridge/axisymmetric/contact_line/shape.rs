//! Inverse (gap, circular radius) geometry using the original first integral.
use super::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Evaluation {
    pub profile: AxisymmetricContactLineProfile,
    pub q: f64,
    pub alpha: f64,
    pub geometric: [f64; 2],
    /// Rows h,b; columns q,alpha. Volume is held exactly in the scale law.
    pub jacobian: [[f64; 2]; 2],
    pub gradient: [f64; 2],
    pub hessian: [[f64; 2]; 2],
    pub gradient_error: [f64; 2],
}

#[derive(Default)]
pub(super) struct Work {
    pub shapes: usize,
    pub quadrature: u64,
}

impl SourceAxisymmetricContactLineBridge {
    pub(super) fn at_parameters(
        &self,
        q: f64,
        alpha: f64,
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Evaluation, ContactLineError> {
        poll(gate)?;
        if work.shapes >= self.config.maximum_shape_evaluations {
            return Err(ContactLineError::Budget);
        }
        work.shapes += 1;
        let theta = FRAC_PI_2 - alpha;
        if !q.is_finite()
            || q < self.original.branch.minimum_suction_parameter
            || q > self.original.branch.maximum_suction_parameter
            || !theta.is_finite()
            || theta < self.config.minimum_angle_radians
            || theta > self.config.maximum_angle_radians
            || (theta - self.original.wetting()[0].contact_angle_radians()).abs()
                > self.config.maximum_angle_departure_radians
        {
            return Err(ContactLineError::Domain(
                "declared suction and dynamic-angle branch",
            ));
        }
        // Exact original constant-mean-curvature quadrature; mirror symmetry
        // needs only ONE half-profile evaluation, not two copied shape owners.
        let integral = super::super::numerics::integrate(alpha, q, self.config.shape, gate)
            .map_err(original_error)?;
        work.quadrature = work
            .quadrature
            .checked_add(integral.evaluations)
            .ok_or(ContactLineError::Budget)?;
        let [h, v, area, hq, vq, _] = integral.values;
        let (x, xq, d) = super::super::numerics::radius(alpha, q);
        let (sa, ca) = det::sin_cos(alpha);
        let ha = x * ca / d;
        let va = x * x * ha;
        let xa = x * sa / d;
        let a = det::pow(self.original.volume / (2.0 * PI * v), 1.0 / 3.0);
        let aq = -a * vq / (3.0 * v);
        let aa = -a * va / (3.0 * v);
        let gap = 2.0 * a * h;
        let b = a * x;
        self.check_coordinates(gap, b)?;
        if !a.is_finite() || a < self.config.minimum_neck_radius_m {
            return Err(ContactLineError::Domain("minimum actual axisymmetric neck"));
        }
        let jacobian = [
            [2.0 * (aq * h + a * hq), 2.0 * (aa * h + a * ha)],
            [aq * x + a * xq, aa * x + a * xa],
        ];
        let inverse = invert(jacobian)?;
        if jacobian[0][0] >= 0.0 {
            return Err(ContactLineError::Domain("monotone suction branch"));
        }
        let gamma = self.original.sigma;
        let eh = PI * gamma * a * (2.0 + q);
        let eb = 4.0 * PI * gamma * b * (sa - self.original.cosine[0]);
        let eq = PI * gamma * (a + (2.0 + q) * aq);
        let ea = PI * gamma * (2.0 + q) * aa;
        let hh = eq * inverse[0][0] + ea * inverse[1][0];
        let hb = eq * inverse[0][1] + ea * inverse[1][1];
        let bh = 4.0 * PI * gamma * b * ca * inverse[1][0];
        let bb = 4.0 * PI * gamma * ((sa - self.original.cosine[0]) + b * ca * inverse[1][1]);
        let relative = integral.relative_error + 64.0 * f64::EPSILON;
        if !bb.is_finite()
            || bb <= 0.0
            || !hh.is_finite()
            || !hb.is_finite()
            || !bh.is_finite()
            || (hb - bh).abs() > (32.0 * relative).max(2e-10) * hb.abs().max(bh.abs())
        {
            return Err(ContactLineError::Domain(
                "positive reciprocal circular-footprint Hessian",
            ));
        }
        let mixed = 0.5 * (hb + bh);
        let hessian = [[hh, mixed], [mixed, bb]];
        // The fixed-footprint perturbation has u=0 on both plates. The same
        // conservative Poincare lower bound used by the original equilibrium
        // owner is sufficient for this narrower perturbation family. Actual
        // dynamic theta replaces Young theta in the radial graph coefficients.
        let margin = (FRAC_PI_2 / (2.0 * h)).powi(2) * ca.powi(3) - 1.0;
        if !margin.is_finite()
            || margin
                <= self.original.branch.minimum_stability_margin + 4.0 * relative * (margin + 1.0)
        {
            return Err(ContactLineError::Domain(
                "axisymmetric fixed-footprint stability screen",
            ));
        }
        let free_area = 4.0 * PI * a * a * area;
        let wet_area = PI * b * b;
        let wet_energy = self
            .original
            .laws
            .map(|law| law.wetted_surface_energy_difference_j_m2() * wet_area);
        let energy = gamma * free_area + wet_energy[0] + wet_energy[1];
        let energy_scale = gamma * free_area + wet_energy[0].abs() + wet_energy[1].abs();
        let pressure = -gamma * q / a;
        let length = gap.max(2.0 * b);
        let bond = self.original.liquid.density_kg_m3() * self.original.limits.gravity_m_s2 / gamma
            * length
            * length;
        if !energy.is_finite()
            || energy.abs() > self.original.limits.maximum_abs_surface_energy_j
            || !pressure.is_finite()
            || -pressure > self.original.limits.maximum_suction_pa
            || !eh.is_finite()
            || eh > self.original.limits.maximum_abs_force_n
            || !eb.is_finite()
            || eb.abs() > self.original.limits.maximum_abs_force_n
            || !bond.is_finite()
            || bond > self.original.limits.maximum_bond_number
        {
            return Err(ContactLineError::Domain(
                "original axisymmetric force/energy/suction/gravity envelope",
            ));
        }
        let coordinate_error = [2.0 * relative * gap, relative * b];
        let gradient_error = [
            hh.abs() * coordinate_error[0]
                + mixed.abs() * coordinate_error[1]
                + relative * eh.abs(),
            mixed.abs() * coordinate_error[0]
                + bb.abs() * coordinate_error[1]
                + 64.0
                    * f64::EPSILON
                    * 4.0
                    * PI
                    * gamma
                    * b
                    * (sa.abs() + self.original.cosine[0].abs()),
        ];
        let energy_error = 4.0 * relative * energy_scale;
        if gradient_error.iter().any(|v| !v.is_finite()) || !energy_error.is_finite() {
            return Err(ContactLineError::Domain(
                "finite axisymmetric numerical uncertainty",
            ));
        }
        Ok(Evaluation {
            profile: AxisymmetricContactLineProfile {
                gap_m: gap,
                wetted_radius_m: b,
                contact_angle_radians: theta,
                suction_parameter: q,
                neck_radius_m: a,
                neck_height_m: 0.5 * gap,
                liquid_volume_m3: self.original.volume,
                liquid_gas_area_m2: free_area,
                wetted_areas_m2: [wet_area; 2],
                liquid_minus_gas_pressure_pa: pressure,
                surface_free_energy_j: energy,
                normal_force_n: -eh,
                footprint_energy_gradient_n: eb,
                normal_force_gap_derivative_n_per_m: -hh,
                normal_force_radius_derivative_n_per_m: -mixed,
                footprint_gradient_radius_derivative_n_per_m: bb,
                axisymmetric_stability_margin: margin,
                geometric_gap_residual_m: 0.0,
                geometric_radius_residual_m: 0.0,
                quadrature_relative_error_estimate: integral.relative_error,
                estimated_volume_error_m3: 4.0 * relative * self.original.volume,
                estimated_surface_energy_error_j: energy_error,
            },
            q,
            alpha,
            geometric: [gap, b],
            jacobian,
            gradient: [eh, eb],
            hessian,
            gradient_error,
        })
    }

    pub(super) fn at_gap_angle(
        &self,
        gap: f64,
        theta: f64,
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Evaluation, ContactLineError> {
        self.original.check_gap(gap).map_err(original_error)?;
        let mut lower = self.original.branch.minimum_suction_parameter;
        let mut upper = self.original.branch.maximum_suction_parameter;
        let mut q = 0.5 * (lower + upper);
        for _ in 0..self.config.shape.maximum_root_iterations {
            let e = self.at_parameters(q, FRAC_PI_2 - theta, work, gate)?;
            let residual = e.geometric[0] - gap;
            if residual.abs() <= self.coordinate_tolerance(gap) {
                return Ok(e.with_coordinates(gap, e.geometric[1]));
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
        gap: f64,
        b: f64,
        seed: Evaluation,
        work: &mut Work,
        gate: &CancelGate,
    ) -> Result<Evaluation, ContactLineError> {
        self.check_coordinates(gap, b)?;
        let target = [gap, b];
        let mut e = seed;
        for _ in 0..self.config.shape.maximum_root_iterations {
            poll(gate)?;
            let residual = [e.geometric[0] - gap, e.geometric[1] - b];
            let norm = self.geometry_norm(residual, target);
            if norm <= 1.0 {
                return Ok(e.with_coordinates(gap, b));
            }
            let inverse = invert(e.jacobian)?;
            let delta = core::array::from_fn::<_, 2, _>(|i| {
                inverse[i][0] * residual[0] + inverse[i][1] * residual[1]
            });
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..self.config.maximum_backtracks {
                let q = e.q - fraction * delta[0];
                let alpha = e.alpha - fraction * delta[1];
                match self.at_parameters(q, alpha, work, gate) {
                    Ok(next) => {
                        let next_norm = self.geometry_norm(
                            [next.geometric[0] - gap, next.geometric[1] - b],
                            target,
                        );
                        if next_norm < norm {
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
                "inverse circular-footprint branch continuation",
            ))?;
        }
        Err(ContactLineError::NotConverged)
    }

    fn geometry_norm(&self, residual: [f64; 2], target: [f64; 2]) -> f64 {
        (residual[0].abs() / self.coordinate_tolerance(target[0]))
            .max(residual[1].abs() / self.coordinate_tolerance(target[1]))
    }
    fn coordinate_tolerance(&self, coordinate: f64) -> f64 {
        self.config.geometry_absolute_tolerance_m
            + self.config.shape.gap_relative_tolerance * coordinate.abs()
    }
    pub(super) fn check_coordinates(&self, h: f64, b: f64) -> Result<(), ContactLineError> {
        if !h.is_finite()
            || h < self.original.limits.minimum_gap_m
            || h > self.original.limits.maximum_gap_m
            || !b.is_finite()
            || b <= 0.0
            || b > self.original.limits.maximum_wetted_radius_m
        {
            return Err(ContactLineError::Domain(
                "original continuum gap and circular footprint envelope",
            ));
        }
        Ok(())
    }
}

impl Evaluation {
    fn with_coordinates(mut self, h: f64, b: f64) -> Self {
        let residual = [self.geometric[0] - h, self.geometric[1] - b];
        // The fixed-volume meridian is numerical. Canonical motion must not
        // inherit its residual as accumulating native gap/footprint drift.
        // Keep residual-induced work/force uncertainty alongside the shape.
        let old_residual = [
            self.profile.geometric_gap_residual_m.abs(),
            self.profile.geometric_radius_residual_m.abs(),
        ];
        for i in 0..2 {
            self.gradient_error[i] += self.hessian[i][0].abs()
                * (residual[0].abs() - old_residual[0]).max(0.0)
                + self.hessian[i][1].abs() * (residual[1].abs() - old_residual[1]).max(0.0);
        }
        self.profile.estimated_surface_energy_error_j += self.gradient[0].abs()
            * (residual[0].abs() - old_residual[0]).max(0.0)
            + self.gradient[1].abs() * (residual[1].abs() - old_residual[1]).max(0.0);
        self.profile.gap_m = h;
        self.profile.wetted_radius_m = b;
        self.profile.geometric_gap_residual_m = residual[0];
        self.profile.geometric_radius_residual_m = residual[1];
        self
    }
}

fn invert(j: [[f64; 2]; 2]) -> Result<[[f64; 2]; 2], ContactLineError> {
    let determinant = j[0][0] * j[1][1] - j[0][1] * j[1][0];
    let scale = (j[0][0] * j[1][1]).abs() + (j[0][1] * j[1][0]).abs();
    if j.iter().flatten().any(|v| !v.is_finite())
        || !determinant.is_finite()
        || determinant.abs() <= 128.0 * f64::EPSILON * scale
    {
        return Err(ContactLineError::Domain(
            "nonsingular fixed-volume shape coordinates",
        ));
    }
    Ok([
        [j[1][1] / determinant, -j[0][1] / determinant],
        [-j[1][0] / determinant, j[0][0] / determinant],
    ])
}
