//! Exact first-integral parameterization; no ODE shooting or profile fitting.
//! The local contact work owner already uses Simpson refinement for estimated
//! work accuracy. Here a fixed six-component dyadic rule integrates height,
//! volume, surface area and their q derivatives together, without allocation.
use super::*;

pub(super) struct Integral {
    pub values: [f64; 6],
    pub relative_error: f64,
    pub evaluations: u64,
}
pub(super) struct Shape {
    pub q: f64,
    pub profile: AxisymmetricProfile,
    pub gap_q_derivative: f64,
    gap_q_derivative_error: f64,
    pub relative_error: f64,
    pub energy_error: f64,
    pub evaluations: u64,
}

/// x=r/a solves q*x²/2+x*cos(beta)=1+q/2 on the one-neck branch.
/// Rationalizing x-1 preserves the thin-gap/large-q geometry and sensitivity.
pub(super) fn radius(beta: f64, q: f64) -> (f64, f64, f64) {
    let c = det::cos(beta);
    let d = (c * c + q * (2.0 + q)).sqrt();
    let s = det::sin(0.5 * beta);
    let excess = 4.0 * s * s / (c + q + d);
    let x = 1.0 + excess;
    let xq = -excess * (x + 1.0) / (2.0 * d);
    (x, xq, d)
}
fn integrand(beta: f64, q: f64) -> Result<[f64; 6], AxisymmetricError> {
    let (x, xq, d) = radius(beta, q);
    let h = x * det::cos(beta) / d;
    let v = x * x * h;
    let area = x * x / d;
    let dx = xq / x;
    let dd = (1.0 + q) / (d * d);
    let values = [
        h,
        v,
        area,
        h * (dx - dd),
        v * (3.0 * dx - dd),
        area * (2.0 * dx - dd),
    ];
    for value in values {
        finite(value, "axisymmetric shape integrand")?;
    }
    Ok(values)
}

pub(super) fn integrate(
    alpha: f64,
    q: f64,
    budget: AxisymmetricNumerics,
    gate: &CancelGate,
) -> Result<Integral, AxisymmetricError> {
    poll(gate)?;
    let left = integrand(0.0, q)?;
    let right = integrand(alpha, q)?;
    let mut trap = std::array::from_fn::<_, 6, _>(|i| alpha * (left[i] + right[i]) * 0.5);
    let mut previous = None;
    let mut panels = 1_u32;
    while panels < budget.maximum_quadrature_panels {
        panels *= 2;
        let step = alpha / f64::from(panels);
        let mut sum = [0.0; 6];
        for index in (1..panels).step_by(2) {
            if index % 64 == 1 {
                poll(gate)?;
            }
            let values = integrand(step * f64::from(index), q)?;
            for i in 0..6 {
                sum[i] += values[i];
            }
        }
        let next = std::array::from_fn::<_, 6, _>(|i| 0.5 * trap[i] + step * sum[i]);
        let simpson = std::array::from_fn::<_, 6, _>(|i| (4.0 * next[i] - trap[i]) / 3.0);
        if let Some(old) = previous {
            let old: [f64; 6] = old;
            let mut relative: f64 = 0.0;
            let corrected = std::array::from_fn::<_, 6, _>(|i| {
                let difference = simpson[i] - old[i];
                relative = relative.max(difference.abs() / (15.0 * simpson[i].abs()));
                simpson[i] + difference / 15.0
            });
            if panels >= 8 && relative <= budget.quadrature_relative_tolerance {
                for value in corrected {
                    finite(value, "axisymmetric shape integral")?;
                }
                poll(gate)?;
                return Ok(Integral {
                    values: corrected,
                    relative_error: relative,
                    evaluations: u64::from(panels) + 1,
                });
            }
        }
        previous = Some(simpson);
        trap = next;
    }
    Err(AxisymmetricError::QuadratureBudget)
}

impl SourceAxisymmetricCapillaryBridge {
    pub(super) fn shape(
        &self,
        q: f64,
        budget: AxisymmetricNumerics,
        gate: &CancelGate,
    ) -> Result<Shape, AxisymmetricError> {
        let lower = integrate(self.alpha[0], q, budget, gate)?;
        let upper = integrate(self.alpha[1], q, budget, gate)?;
        let sums = std::array::from_fn::<_, 6, _>(|i| lower.values[i] + upper.values[i]);
        let [h, v, s, hq, vq, _] = sums;
        positive(h, "dimensionless bridge height")?;
        positive(v, "dimensionless bridge volume")?;
        let a = positive(
            det::pow(
                positive(self.volume / (PI * v), "axisymmetric cubic scale")?,
                1.0 / 3.0,
            ),
            "axisymmetric neck scale",
        )?;
        let gap = positive(a * h, "axisymmetric profile gap")?;
        let aq = -a * vq / (3.0 * v);
        let gapq = finite(a * hq + aq * h, "axisymmetric gap sensitivity")?;
        let margin = finite(
            (FRAC_PI_2 / h).powi(2) * self.minimum_sine.powi(3) - 1.0,
            "axisymmetric stability margin",
        )?;
        let relative_error = lower.relative_error.max(upper.relative_error);
        let gapq_error = finite(
            (a * hq.abs() + aq.abs() * h) * (3.0 * relative_error + 32.0 * f64::EPSILON),
            "axisymmetric gap sensitivity error estimate",
        )?;
        // A positive caller margin is retained above the quadrature estimate;
        // this does not promote an estimate into an interval certificate.
        if margin <= self.branch.minimum_stability_margin + 4.0 * relative_error * (margin + 1.0) {
            return Err(AxisymmetricError::StabilityUnestablished {
                margin,
                required: self.branch.minimum_stability_margin,
            });
        }
        if gapq >= 0.0 {
            return Err(AxisymmetricError::TurningBranch);
        }
        let mut radii = [0.0; 2];
        let mut radius_gap = [0.0; 2];
        for i in 0..2 {
            let (x, xq, _) = radius(self.alpha[i], q);
            radii[i] = positive(a * x, "axisymmetric contact radius")?;
            radius_gap[i] = finite((aq * x + a * xq) / gapq, "contact radius sensitivity")?;
        }
        let areas = radii.map(|r| PI * r * r);
        let area = positive(2.0 * PI * a * a * s, "axisymmetric free surface area")?;
        let force = finite(
            -PI * self.sigma * a * (2.0 + q),
            "axisymmetric normal force",
        )?;
        let derivative = finite(
            -PI * self.sigma * (a + aq * (2.0 + q)) / gapq,
            "axisymmetric force tangent",
        )?;
        let energy = finite(
            self.sigma * area
                + self.laws[0].wetted_surface_energy_difference_j_m2() * areas[0]
                + self.laws[1].wetted_surface_energy_difference_j_m2() * areas[1],
            "axisymmetric free energy",
        )?;
        let energy_scale = finite(
            self.sigma * (area + self.cosine[0] * areas[0] + self.cosine[1] * areas[1]),
            "surface energy scale",
        )?;
        let energy_error = finite(
            (3.0 * relative_error + 32.0 * f64::EPSILON) * energy_scale,
            "surface energy error estimate",
        )?;
        Ok(Shape {
            q,
            gap_q_derivative: gapq,
            gap_q_derivative_error: gapq_error,
            relative_error,
            energy_error,
            evaluations: lower.evaluations + upper.evaluations,
            profile: AxisymmetricProfile {
                gap_m: gap,
                liquid_volume_m3: self.volume,
                neck_radius_m: a,
                neck_height_m: a * lower.values[0],
                contact_radii_m: radii,
                contact_angles_radians: self.alpha.map(|alpha| FRAC_PI_2 - alpha),
                liquid_gas_area_m2: area,
                wetted_areas_m2: areas,
                liquid_minus_gas_pressure_pa: finite(-self.sigma * q / a, "axisymmetric pressure")?,
                normal_force_n: force,
                force_gap_derivative_n_per_m: derivative,
                contact_radius_gap_derivatives: radius_gap,
                surface_free_energy_j: energy,
                axisymmetric_stability_margin: margin,
            },
        })
    }

    /// Integrate F dh without subtracting endpoint energies or nearby forces.
    /// On t in [0,1], q=q0+t(q1-q0), w=-dh/dq is positive. Then
    /// F_mean = integral(F w dt)/integral(w dt). The common q increment
    /// cancels analytically, so neither integral shrinks with the timestep.
    /// With H(t)=integral_0^t w ds, the exact mean-force endpoint derivative is
    /// integral(F_h H w dt)/integral(w dt)^2. This avoids another cancellation
    /// in (F_end-F_mean)/Delta_h while retaining the finite-path tangent.
    pub(super) fn integrate_force_work(
        &self,
        previous: AxisymmetricState,
        next: AxisymmetricState,
        budget: AxisymmetricWorkBudget,
        gate: &CancelGate,
    ) -> Result<(f64, f64, f64, f64), AxisymmetricError> {
        let dh = next.profile.gap_m - previous.profile.gap_m;
        if previous.q == next.q {
            return Err(AxisymmetricError::WorkUnresolved {
                estimated_error_j: previous.report.surface_energy_error_estimate_j
                    + next.report.surface_energy_error_estimate_j,
                tolerance_j: budget.absolute_tolerance_j,
            });
        }
        let mut panels = 1;
        let mut coarse = self.force_work_rule(previous, next, panels, budget.equilibrium, gate)?;
        while panels < budget.equilibrium.maximum_quadrature_panels {
            panels *= 2;
            let fine = self.force_work_rule(previous, next, panels, budget.equilibrium, gate)?;
            let force_error = finite(
                fine.force_error + (fine.force - coarse.force).abs() / 15.0,
                "axisymmetric integrated force error",
            )?;
            let work = finite(-fine.force * dh, "axisymmetric directly integrated work")?;
            let error = finite(
                dh.abs() * force_error + 8.0 * f64::EPSILON * work.abs(),
                "axisymmetric directly integrated work error",
            )?;
            let tolerance = budget.absolute_tolerance_j + budget.relative_tolerance * work.abs();
            let derivative_tolerance = (budget.equilibrium.quadrature_relative_tolerance
                + 64.0 * f64::EPSILON)
                * fine.derivative.abs().max(coarse.derivative.abs());
            if error <= tolerance
                && error < work.abs()
                && (fine.derivative - coarse.derivative).abs() / 15.0 <= derivative_tolerance
            {
                poll(gate)?;
                return Ok((work, error, fine.force, fine.derivative));
            }
            coarse = fine;
        }
        Err(AxisymmetricError::WorkQuadratureBudget)
    }

    fn force_work_rule(
        &self,
        previous: AxisymmetricState,
        next: AxisymmetricState,
        panels: u32,
        budget: AxisymmetricNumerics,
        gate: &CancelGate,
    ) -> Result<ForceWork, AxisymmetricError> {
        let dq = next.q - previous.q;
        let dt = 1.0 / f64::from(panels);
        let left = self.force_work_sample(previous.q, budget, gate)?;
        // Both normalized integrals are invariant under a common weight scale.
        // Remove physical length here so H*w cannot underflow for a very small
        // admitted bridge or overflow for a large one. Never divide by dq.
        let weight_scale = left.weight;
        let mut left = left.normalize_weight(weight_scale)?;
        let initial_gap_error = (left.gap - previous.profile.gap_m).abs();
        let mut weight = 0.0;
        let mut force = 0.0;
        let mut derivative = 0.0;
        let mut weight_error = 0.0;
        let mut force_error = 0.0;
        let mut maximum_derivative = left.derivative.abs();
        for panel in 0..panels {
            poll(gate)?;
            let middle_q = previous.q + (f64::from(panel) + 0.5) * dt * dq;
            let right_q = if panel + 1 == panels {
                next.q
            } else {
                previous.q + f64::from(panel + 1) * dt * dq
            };
            let middle = self
                .force_work_sample(middle_q, budget, gate)?
                .normalize_weight(weight_scale)?;
            let right = self
                .force_work_sample(right_q, budget, gate)?
                .normalize_weight(weight_scale)?;
            // The integral of the same quadratic weight interpolant at the
            // midpoint is dt*(5*w0+8*wm-w1)/24. Carry H cumulatively; never
            // subtract two nearly equal absolute profile heights.
            let mid_weight =
                weight + dt * (5.0 * left.weight + 8.0 * middle.weight - right.weight) / 24.0;
            let end_weight = weight + dt * (left.weight + 4.0 * middle.weight + right.weight) / 6.0;
            force += dt
                * (left.force * left.weight
                    + 4.0 * middle.force * middle.weight
                    + right.force * right.weight)
                / 6.0;
            derivative += dt
                * (left.derivative * weight * left.weight
                    + 4.0 * middle.derivative * mid_weight * middle.weight
                    + right.derivative * end_weight * right.weight)
                / 6.0;
            weight_error +=
                dt * (left.weight_error + 4.0 * middle.weight_error + right.weight_error) / 6.0;
            force_error += dt
                * (left.weighted_force_error()?
                    + 4.0 * middle.weighted_force_error()?
                    + right.weighted_force_error()?)
                / 6.0;
            weight = end_weight;
            maximum_derivative = maximum_derivative
                .max(middle.derivative.abs())
                .max(right.derivative.abs());
            left = right;
        }
        positive(weight - weight_error, "resolved axisymmetric work weight")?;
        let mean = finite(
            force / weight,
            "axisymmetric directly integrated mean force",
        )?;
        let derivative = finite(
            (derivative / weight) / weight,
            "axisymmetric integrated mean-force tangent",
        )?;
        let profile_gap_error = initial_gap_error + (left.gap - next.profile.gap_m).abs();
        let error = finite(
            (force_error + mean.abs() * weight_error) / (weight - weight_error)
                + maximum_derivative * profile_gap_error
                + 32.0 * f64::EPSILON * mean.abs(),
            "axisymmetric weighted force error",
        )?;
        Ok(ForceWork {
            force: mean,
            derivative,
            force_error: error,
        })
    }

    fn force_work_sample(
        &self,
        q: f64,
        budget: AxisymmetricNumerics,
        gate: &CancelGate,
    ) -> Result<ForceWorkSample, AxisymmetricError> {
        let shape = self.shape(q, budget, gate)?;
        let weight = -shape.gap_q_derivative;
        let weight_error = shape.gap_q_derivative_error;
        let force = shape.profile.normal_force_n;
        let force_error = force.abs() * (shape.relative_error / 3.0 + 32.0 * f64::EPSILON);
        let sample = ForceWorkSample {
            gap: shape.profile.gap_m,
            force,
            derivative: shape.profile.force_gap_derivative_n_per_m,
            weight,
            weight_error,
            force_error,
        };
        // Every sampled profile obeys the same original force, gravity and
        // continuum limits as the accepted endpoint equilibria.
        self.accept(shape, 0.0, 0)?;
        Ok(sample)
    }

    pub(super) fn solve_gap(
        &self,
        gap: f64,
        guess: Option<f64>,
        budget: AxisymmetricNumerics,
        gate: &CancelGate,
    ) -> Result<AxisymmetricState, AxisymmetricError> {
        budget.validate()?;
        poll(gate)?;
        self.check_gap(gap)?;
        let mut lower = self.shape(self.branch.minimum_suction_parameter, budget, gate)?;
        let mut upper = self.shape(self.branch.maximum_suction_parameter, budget, gate)?;
        let tolerance = positive(
            budget.gap_relative_tolerance * gap,
            "physical gap tolerance",
        )?;
        let mut count = lower.evaluations + upper.evaluations;
        let minimum = upper.profile.gap_m;
        let maximum = lower.profile.gap_m;
        if gap < minimum - tolerance || gap > maximum + tolerance {
            return Err(AxisymmetricError::GapOutsideBranch {
                minimum_m: minimum,
                maximum_m: maximum,
            });
        }
        if (lower.profile.gap_m - gap).abs() <= tolerance {
            let residual = lower.profile.gap_m - gap;
            lower.evaluations = count;
            return self.accept(lower, residual, 0);
        }
        if (upper.profile.gap_m - gap).abs() <= tolerance {
            let residual = upper.profile.gap_m - gap;
            upper.evaluations = count;
            return self.accept(upper, residual, 0);
        }
        let mut q = guess
            .filter(|q| *q > lower.q && *q < upper.q)
            .unwrap_or(0.5 * lower.q + 0.5 * upper.q);
        for iteration in 1..=budget.maximum_root_iterations {
            let mut shape = self.shape(q, budget, gate)?;
            count += shape.evaluations;
            let residual = shape.profile.gap_m - gap;
            if residual.abs() <= tolerance {
                shape.evaluations = count;
                poll(gate)?;
                return self.accept(shape, residual, iteration);
            }
            let next = q - residual / shape.gap_q_derivative;
            if residual > 0.0 {
                lower = shape;
            } else {
                upper = shape;
            }
            let width = upper.q - lower.q;
            // A safeguarded step strictly contracts the bracket even when the
            // original q interval spans very thin and thick profiles.
            q = if next.is_finite()
                && next > lower.q + 0.05 * width
                && next < upper.q - 0.05 * width
            {
                next
            } else {
                0.5 * lower.q + 0.5 * upper.q
            };
            if q == lower.q || q == upper.q {
                break;
            }
        }
        Err(AxisymmetricError::RootBudget)
    }
}

struct ForceWorkSample {
    gap: f64,
    force: f64,
    derivative: f64,
    weight: f64,
    weight_error: f64,
    force_error: f64,
}
impl ForceWorkSample {
    fn normalize_weight(mut self, scale: f64) -> Result<Self, AxisymmetricError> {
        self.weight = positive(self.weight / scale, "normalized axisymmetric work weight")?;
        self.weight_error = finite(
            self.weight_error / scale,
            "normalized axisymmetric work weight error",
        )?;
        Ok(self)
    }

    fn weighted_force_error(&self) -> Result<f64, AxisymmetricError> {
        Ok(finite(
            self.force_error * (self.weight + self.weight_error)
                + self.force.abs() * self.weight_error,
            "axisymmetric work integrand error",
        )?)
    }
}

struct ForceWork {
    force: f64,
    derivative: f64,
    force_error: f64,
}
