//! Coupled axial and azimuthal flow of one original thermal Carreau law.
//!
//! The invariant rate is hypot(du_z/dr, r*domega/dr). Both stress components
//! therefore share ONE viscosity. Independent nonlinear axial and rotational
//! solutions cannot be superposed. Momentum fixes tau_z=C_z/r-G*r/2 and
//! tau_theta=C_theta/r^2. A bounded two-variable Newton solve enforces both
//! wall-speed differences, using the exact positive inverse constitutive tensor.
//!
//! Log-radius quadrature evaluates velocity, flow, heat and the symmetric
//! wall mobility together, reusing every constitutive query and refinement.
//! The squared stress magnitude is convex in log(radius), so admitting both
//! wall magnitudes admits the entire gap. All work is explicitly bounded;
//! preparation, solves and profile samples allocate no scratch.
//!
//! This is a steady, fully developed, no-slip laminar reduction at one frozen
//! temperature. It supplies no entrance effects, radial pressure/density,
//! centrifugal instability, elastic normal stresses, slip or thermal evolution.

use super::super::pipe::{CarreauPipeControls, RateState, product_ratio};
use super::{
    LaminarCarreauAnnulus, MaterialError, log_ratio_from_gap, parameters, positive, state_error,
};
use fs_math::det;

/// Independent nonlinear, quadrature and total constitutive-query budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HelicalCarreauControls {
    pub local: CarreauPipeControls,
    /// Applied to axial speed and inner-radius times angular-speed difference.
    pub absolute_velocity_tolerance_m_s: f64,
    pub relative_velocity_tolerance: f64,
    pub maximum_newton_iterations: u32,
    pub maximum_backtracks: u32,
    /// Each query has its own bounded original-law constitutive inverse.
    pub maximum_constitutive_queries: u64,
}
impl Default for HelicalCarreauControls {
    fn default() -> Self {
        Self {
            local: CarreauPipeControls::default(),
            absolute_velocity_tolerance_m_s: 1e-12,
            relative_velocity_tolerance: 1e-9,
            maximum_newton_iterations: 64,
            maximum_backtracks: 32,
            maximum_constitutive_queries: 1_000_000,
        }
    }
}
impl HelicalCarreauControls {
    fn admit(self) -> Result<(), MaterialError> {
        self.local.admit()?;
        if !positive(self.absolute_velocity_tolerance_m_s)
            || !positive(self.relative_velocity_tolerance)
            || self.relative_velocity_tolerance >= 1.0
            || !(1..=128).contains(&self.maximum_newton_iterations)
            || self.maximum_backtracks > 64
            || self.maximum_constitutive_queries == 0
        {
            return Err(parameters("invalid helical Carreau work bounds"));
        }
        Ok(())
    }
}

/// Accepted cylindrical momentum and mechanical/thermal ports.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauHelicalFlow {
    pub pressure_drop_pa: f64,
    pub requested_axial_wall_velocities_m_s: [f64; 2],
    pub accepted_axial_wall_velocities_m_s: [f64; 2],
    pub requested_wall_angular_velocities_rad_s: [f64; 2],
    pub accepted_wall_angular_velocities_rad_s: [f64; 2],
    /// Actual differences before adding common translation or rotation.
    pub relative_axial_velocity_m_s: f64,
    pub relative_angular_velocity_rad_s: f64,
    /// Accepted minus requested [axial speed, inner-radius * angular speed].
    pub wall_velocity_residual_m_s: [f64; 2],
    pub volumetric_flow_m3_s: f64,
    /// Relative to the arithmetic mean of the ACCEPTED axial wall velocities.
    pub relative_volume_flow_m3_s: f64,
    /// Fluid-domain outward-normal tractions [axial, azimuthal] at each wall.
    pub wall_tractions_pa: [[f64; 2]; 2],
    pub axial_wall_forces_n: [f64; 2],
    /// Equal/opposite torques exerted on the fluid.
    pub wall_torques_n_m: [f64; 2],
    /// Sum of pressure and requested-wall work, with common-motion cancellation
    /// removed analytically before floating-point evaluation.
    pub mechanical_power_w: f64,
    /// Positive volume integral of the original stress dot strain rate.
    pub dissipated_power_w: f64,
    pub work_residual_w: f64,
    /// At fixed pressure: derivative of [(Fouter-Finner)/2, Touter/a]
    /// with respect to [Uouter-Uinner, a*(Omega_outer-Omega_inner)].
    /// Symmetric and positive definite; a is the original inner radius.
    pub wall_resistance_n_s_m: [[f64; 2]; 2],
    pub maximum_shear_rate_per_s: f64,
    /// Observed quadrature refinement change, not a certified error bound.
    pub quadrature_relative_change: f64,
    pub newton_iterations: u32,
    pub backtracks: u32,
    pub constitutive_queries: u64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauHelicalProfile {
    annulus: LaminarCarreauAnnulus,
    controls: HelicalCarreauControls,
    flow: CarreauHelicalFlow,
    inner_stress_pa: [f64; 2],
    pressure_stress_pa: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CarreauHelicalSample {
    pub axial_velocity_m_s: f64,
    pub angular_velocity_rad_s: f64,
    pub azimuthal_velocity_m_s: f64,
    /// [du_z/dr, r*domega/dr], excluding common rigid motion.
    pub shear_rate_components_per_s: [f64; 2],
    pub shear_rate_per_s: f64,
    pub radial_shear_stress_pa: [f64; 2],
    pub dynamic_viscosity_pa_s: f64,
    pub dissipation_w_m3: f64,
    pub quadrature_relative_change: f64,
}
struct Budget {
    used: u64,
    maximum: u64,
}
impl Budget {
    fn query(
        &mut self,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<(), MaterialError> {
        check()?;
        if self.used == self.maximum {
            return Err(state_error(
                "helical Carreau constitutive-query budget exhausted",
            ));
        }
        self.used += 1;
        Ok(())
    }
}
struct Evaluation {
    // du, a*dOmega, relative flow/area, heat, J00, J01, J11.
    values: [f64; 7],
    maximum_rate: f64,
    change: f64,
}
fn finite(value: f64) -> Result<f64, MaterialError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(state_error("helical Carreau response is unrepresentable"))
    }
}
fn product<const N: usize>(factors: [f64; N]) -> Result<f64, MaterialError> {
    if factors.iter().any(|v| !v.is_finite()) {
        return Err(state_error("helical Carreau factor is nonfinite"));
    }
    if factors.contains(&0.0) {
        return Ok(0.0);
    }
    let negative = factors.iter().filter(|v| v.is_sign_negative()).count() % 2 != 0;
    let value = product_ratio(&factors.map(f64::abs), &[])?;
    Ok(if negative { -value } else { value })
}
fn weighted_ratio(geometry: f64, weight: f64, denominator: f64) -> Result<f64, MaterialError> {
    if weight == 0.0 {
        Ok(0.0)
    } else {
        product_ratio(&[geometry, weight], &[denominator])
    }
}
fn compensated(values: &[f64]) -> Result<f64, MaterialError> {
    let (mut sum, mut correction) = (0.0, 0.0);
    for &value in values {
        let adjusted = finite(value)? - correction;
        let next = finite(sum + adjusted)?;
        correction = (next - sum) - adjusted;
        sum = next;
    }
    Ok(sum)
}

/// Normalized symmetric solve avoids squaring dimensionful matrix entries.
fn solve(m: [f64; 3], rhs: [f64; 2]) -> Result<[f64; 2], MaterialError> {
    if !positive(m[0]) || !positive(m[2]) {
        return Err(state_error("helical Carreau mobility is not positive"));
    }
    let s = [det::sqrt(m[0]), det::sqrt(m[2])];
    let correlation = (m[1] / s[0]) / s[1];
    let determinant = (1.0 - correlation) * (1.0 + correlation);
    if !positive(determinant) {
        return Err(state_error(
            "helical Carreau mobility lost positive definiteness",
        ));
    }
    let r = [rhs[0] / s[0], rhs[1] / s[1]];
    Ok([
        finite((-correlation).mul_add(r[1], r[0]) / determinant / s[0])?,
        finite((-correlation).mul_add(r[0], r[1]) / determinant / s[1])?,
    ])
}

impl LaminarCarreauAnnulus {
    /// Simultaneous pressure-driven flow, axial wall motion and rotation.
    pub fn helical_at_pressure_and_walls(
        &self,
        pressure_drop_pa: f64,
        axial_wall_velocities_m_s: [f64; 2],
        wall_angular_velocities_rad_s: [f64; 2],
        controls: HelicalCarreauControls,
    ) -> Result<CarreauHelicalProfile, MaterialError> {
        self.helical_with_check(
            pressure_drop_pa,
            axial_wall_velocities_m_s,
            wall_angular_velocities_rad_s,
            controls,
            || Ok(()),
        )
    }

    /// Bounded source queries and transactional cancellation; the owner is immutable.
    pub fn helical_with_check(
        &self,
        pressure: f64,
        axial: [f64; 2],
        angular: [f64; 2],
        controls: HelicalCarreauControls,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauHelicalProfile, MaterialError> {
        controls.admit()?;
        check()?;
        if [pressure, axial[0], axial[1], angular[0], angular[1]]
            .into_iter()
            .any(|v| !v.is_finite())
        {
            return Err(parameters(
                "helical Carreau pressure and wall speeds must be finite",
            ));
        }
        let target = [
            finite(axial[1] - axial[0])?,
            product([self.inner_radius_m, finite(angular[1] - angular[0])?])?,
        ];
        let pressure_stress = if pressure == 0.0 {
            0.0
        } else {
            product_ratio(&[pressure.abs(), self.gap_m], &[2.0, self.length_m])?.copysign(pressure)
        };
        let ceiling = self.prepared.ceiling.stress();
        if pressure_stress.abs() > ceiling {
            return Err(state_error(
                "helical pressure exceeds every source-supported stress field",
            ));
        }
        let mean_rate = target[0].hypot(target[1]) / self.gap_m;
        if !mean_rate.is_finite() || mean_rate > self.prepared.ceiling.rate {
            return Err(state_error(
                "helical wall motion exceeds every source-supported rate field",
            ));
        }
        // A balanced-pressure seed always belongs to the convex stress domain.
        let mut stress = [pressure_stress, 0.0];
        if self.prepared.constant() {
            let mu = self.prepared.rest.viscosity;
            let uniform = if target[0] == 0.0 {
                0.0
            } else {
                product_ratio(
                    &[mu, target[0].abs()],
                    &[self.mean_radius_m, self.log_ratio],
                )?
                .copysign(target[0])
            };
            let offset = (-pressure_stress).mul_add(self.langevin, uniform);
            let inner_axial = product([
                self.mean_radius_m / self.inner_radius_m,
                offset + pressure_stress,
            ])?;
            let span = (self.gap_m / self.outer_radius_m)
                * (1.0 + self.inner_radius_m / self.outer_radius_m);
            let inner_rotation = if target[1] == 0.0 {
                0.0
            } else {
                product_ratio(&[2.0, mu, target[1].abs()], &[self.inner_radius_m, span])?
                    .copysign(target[1])
            };
            stress = [inner_axial, inner_rotation];
        }
        if !self.helical_admitted(stress, pressure_stress) {
            return Err(state_error(
                "helical constant-law wall stress leaves source support",
            ));
        }
        let mut budget = Budget {
            used: 0,
            maximum: controls.maximum_constitutive_queries,
        };
        let mut backtracks = 0;
        let mut evaluation =
            self.helical_evaluate(stress, pressure_stress, controls, &mut budget, &mut check)?;
        for iteration in 0..=controls.maximum_newton_iterations {
            let residual = [
                evaluation.values[0] - target[0],
                evaluation.values[1] - target[1],
            ];
            let matrix = [
                evaluation.values[4],
                evaluation.values[5],
                evaluation.values[6],
            ];
            let correction = solve(matrix, residual)?;
            let norm = residual[0].hypot(residual[1]);
            let scale = target[0]
                .hypot(target[1])
                .max(product([self.gap_m, evaluation.maximum_rate])?);
            let tolerance = controls.absolute_velocity_tolerance_m_s
                + controls.relative_velocity_tolerance * scale;
            let stress_scale = stress[0].hypot(stress[1]).max(pressure_stress.abs());
            if norm <= tolerance
                && (norm == 0.0
                    || correction[0].hypot(correction[1])
                        <= controls.relative_velocity_tolerance * stress_scale)
            {
                return self.helical_finish(
                    pressure,
                    axial,
                    angular,
                    target,
                    stress,
                    pressure_stress,
                    evaluation,
                    controls,
                    iteration,
                    backtracks,
                    budget.used,
                    &mut check,
                );
            }
            if iteration == controls.maximum_newton_iterations {
                return Err(state_error("helical Carreau Newton budget exhausted"));
            }
            let mut accepted = None;
            let mut damping = 1.0;
            for attempt in 0..=controls.maximum_backtracks {
                check()?;
                let candidate = [
                    (-damping).mul_add(correction[0], stress[0]),
                    (-damping).mul_add(correction[1], stress[1]),
                ];
                if candidate != stress && self.helical_admitted(candidate, pressure_stress) {
                    let trial = self.helical_evaluate(
                        candidate,
                        pressure_stress,
                        controls,
                        &mut budget,
                        &mut check,
                    )?;
                    let trial_norm =
                        (trial.values[0] - target[0]).hypot(trial.values[1] - target[1]);
                    if trial_norm <= (1.0 - 1e-4 * damping) * norm {
                        accepted = Some((candidate, trial));
                        break;
                    }
                }
                if attempt < controls.maximum_backtracks {
                    backtracks += 1;
                    damping *= 0.5;
                }
            }
            let (next_stress, next_evaluation) = accepted.ok_or_else(|| {
                state_error("helical Carreau line search or source span exhausted")
            })?;
            stress = next_stress;
            evaluation = next_evaluation;
        }
        Err(state_error("helical Carreau solve did not produce a state"))
    }

    fn helical_stress(&self, x: f64, inner: [f64; 2], pressure: f64) -> [f64; 2] {
        if x == 0.0 {
            return inner;
        }
        let (r, y) = self.point(x);
        let ar = self.inner_radius_m / r;
        // Recover physical gap fraction from normalized annular area without
        // subtracting nearly equal radii.
        let am = self.inner_radius_m / self.mean_radius_m;
        let gm = self.gap_m / self.mean_radius_m;
        let fraction = if x == 1.0 {
            1.0
        } else {
            2.0 * y / (am + det::sqrt(am * am + 2.0 * gm * y))
        };
        let p = pressure * fraction;
        [ar.mul_add(inner[0] - p, -p), ar * ar * inner[1]]
    }

    fn helical_admitted(&self, inner: [f64; 2], pressure: f64) -> bool {
        let outer = self.helical_stress(1.0, inner, pressure);
        inner.into_iter().chain(outer).all(f64::is_finite)
            && inner[0].hypot(inner[1]) <= self.prepared.ceiling.stress()
            && outer[0].hypot(outer[1]) <= self.prepared.ceiling.stress()
    }

    fn helical_local(
        &self,
        x: f64,
        stress: [f64; 2],
        pressure: f64,
        controls: HelicalCarreauControls,
        budget: &mut Budget,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<([f64; 2], RateState), MaterialError> {
        budget.query(check)?;
        let tau = self.helical_stress(x, stress, pressure);
        let magnitude = tau[0].hypot(tau[1]);
        let local = self
            .prepared
            .inverse(magnitude, self.prepared.ceiling, controls.local)?
            .0;
        Ok((tau, local))
    }

    fn helical_integral(
        &self,
        stress: [f64; 2],
        pressure: f64,
        bounds: [f64; 2],
        controls: HelicalCarreauControls,
        budget: &mut Budget,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<Evaluation, MaterialError> {
        let [lower, upper] = bounds;
        if lower == upper {
            return Ok(Evaluation {
                values: [0.0; 7],
                maximum_rate: 0.0,
                change: 0.0,
            });
        }
        let (cuts, cut_count) = self.helical_cuts(stress, pressure, bounds, budget, check)?;
        let mut maximum_rate = 0.0_f64;
        let geometry = product([self.inner_radius_m, self.log_ratio])?;
        let mut node = |x: f64| -> Result<[f64; 7], MaterialError> {
            let (tau, local) = self.helical_local(x, stress, pressure, controls, budget, check)?;
            maximum_rate = maximum_rate.max(local.rate);
            let (r, y) = self.point(x);
            let magnitude = tau[0].hypot(tau[1]);
            let direction = if magnitude == 0.0 {
                [1.0, 0.0]
            } else {
                [tau[0] / magnitude, tau[1] / magnitude]
            };
            let gradient = direction.map(|v| v * local.rate);
            let ratio = self.inner_radius_m / r;
            let j00 = weighted_ratio(geometry, direction[1] * direction[1], local.viscosity)?
                + weighted_ratio(geometry, direction[0] * direction[0], local.tangent)?;
            let j11 = weighted_ratio(geometry, (ratio * direction[0]).powi(2), local.viscosity)?
                + weighted_ratio(geometry, (ratio * direction[1]).powi(2), local.tangent)?;
            let contrast = local.viscosity - local.tangent;
            let cross = direction[0] * direction[1];
            let j01 = if contrast == 0.0 || cross == 0.0 {
                0.0
            } else {
                product_ratio(
                    &[geometry, ratio, cross.abs(), contrast.abs()],
                    &[local.viscosity, local.tangent],
                )?
                .copysign(cross * contrast)
            };
            let axial = product([self.log_ratio, r, gradient[0]])?;
            Ok([
                axial,
                product([geometry, gradient[1]])?,
                product([0.5 - y, axial])?,
                product([
                    core::f64::consts::TAU,
                    self.length_m,
                    self.log_ratio,
                    r,
                    r,
                    magnitude,
                    local.rate,
                ])?,
                finite(j00)?,
                j01,
                finite(j11)?,
            ])
        };
        let mut values = [0.0; 7];
        let mut change = 0.0_f64;
        for interval in cuts[..cut_count].windows(2) {
            let (part, estimate) =
                integrate(&mut node, [interval[0], interval[1]], controls.local)?;
            for i in 0..7 {
                values[i] = finite(values[i] + part[i])?;
            }
            change = change.max(estimate);
        }
        Ok(Evaluation {
            values,
            maximum_rate,
            change,
        })
    }

    // Stress-squared is convex in log radius: one minimum, at most two
    // crossings of the original material transition. Resolve these geometric
    // sites before quadrature so a sharp Yasuda crossover cannot hide between
    // early nodes. Fixed bisections evaluate momentum, never a substituted law.
    fn helical_cuts(
        &self,
        stress: [f64; 2],
        pressure: f64,
        bounds: [f64; 2],
        budget: &mut Budget,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<([f64; 6], usize), MaterialError> {
        let derivative = |x: f64| {
            let tau = self.helical_stress(x, stress, pressure);
            let scale = tau[0].abs().max(tau[1].abs());
            if scale == 0.0 {
                return 0.0;
            }
            if pressure == 0.0
                || tau[0] == 0.0
                || pressure.is_sign_negative() == tau[0].is_sign_negative()
            {
                return -1.0;
            }
            let z = tau[0] / scale;
            let t = tau[1] / scale;
            // Sign of -(gap/r)*(tau_z^2+2*tau_theta^2)-2*P*tau_z.
            // Compare positive terms in logs without squaring large stresses.
            let lhs = 2.0 * det::ln(scale) + det::ln(self.gap_m) - det::ln(self.point(x).0)
                + det::ln(z * z + 2.0 * t * t);
            let rhs = core::f64::consts::LN_2 + det::ln(pressure.abs()) + det::ln(tau[0].abs());
            rhs - lhs
        };
        let minimum = if derivative(0.0) >= 0.0 {
            0.0
        } else if derivative(1.0) <= 0.0 {
            1.0
        } else {
            let (mut lo, mut hi) = (0.0, 1.0);
            for _ in 0..64 {
                check()?;
                let mid = 0.5 * lo + 0.5 * hi;
                if mid == lo || mid == hi {
                    break;
                }
                if derivative(mid) < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * lo + 0.5 * hi
        };
        let mut cuts = [bounds[0]; 6];
        cuts[1] = bounds[1];
        let mut count = 2;
        let mut push = |x: f64| {
            if x > bounds[0] && x < bounds[1] && !cuts[..count].contains(&x) {
                cuts[count] = x;
                count += 1;
            }
        };
        push(minimum);
        if pressure != 0.0 {
            let y = 0.5 * (self.inner_radius_m / self.mean_radius_m) * (stress[0] / pressure);
            if y > 0.0 && y < 1.0 {
                push(self.coordinate_at_area(y));
            }
        }
        if let Some(rate) = self.prepared.transition_rate(self.prepared.ceiling.rate) {
            budget.query(check)?;
            let transition = self.prepared.sample_rate(rate)?.stress();
            for [start, end] in [[0.0, minimum], [minimum, 1.0]] {
                let norm = |x: f64| {
                    let t = self.helical_stress(x, stress, pressure);
                    t[0].hypot(t[1])
                };
                let (left, right) = (norm(start), norm(end));
                if start == end || transition <= left.min(right) || transition >= left.max(right) {
                    continue;
                }
                let (mut lo, mut hi) = (start, end);
                for _ in 0..64 {
                    check()?;
                    let mid = 0.5 * lo + 0.5 * hi;
                    if mid == lo || mid == hi {
                        break;
                    }
                    if (norm(mid) < transition) == (left < right) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                push(0.5 * lo + 0.5 * hi);
            }
        }
        cuts[..count].sort_unstable_by(f64::total_cmp);
        Ok((cuts, count))
    }

    fn helical_evaluate(
        &self,
        stress: [f64; 2],
        pressure: f64,
        controls: HelicalCarreauControls,
        budget: &mut Budget,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<Evaluation, MaterialError> {
        self.helical_integral(stress, pressure, [0.0, 1.0], controls, budget, check)
    }

    #[allow(clippy::too_many_arguments)]
    fn helical_finish(
        &self,
        pressure: f64,
        axial: [f64; 2],
        angular: [f64; 2],
        target: [f64; 2],
        stress: [f64; 2],
        pressure_stress: f64,
        evaluation: Evaluation,
        controls: HelicalCarreauControls,
        iterations: u32,
        backtracks: u32,
        queries: u64,
        check: &mut impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauHelicalProfile, MaterialError> {
        let v = evaluation.values;
        let residual = [finite(v[0] - target[0])?, finite(v[1] - target[1])?];
        let relative_flow = product([self.area_m2, v[2]])?;
        let mean = finite(axial[0] + 0.5 * v[0])?;
        let flow = finite(product([self.area_m2, mean])? + relative_flow)?;
        let k = product([core::f64::consts::TAU, self.length_m, self.inner_radius_m])?;
        let inner_force = -product([k, stress[0]])?;
        let pressure_force = product([pressure, self.area_m2])?;
        let outer_force = finite(-inner_force - pressure_force)?;
        let relative_force = finite(-inner_force - 0.5 * pressure_force)?;
        let rotation_force = product([k, stress[1]])?;
        let torque = product([rotation_force, self.inner_radius_m])?;
        let power = compensated(&[
            product([pressure, relative_flow])?,
            product([relative_force, v[0]])?,
            product([rotation_force, v[1]])?,
            -product([outer_force, residual[0]])?,
            -product([rotation_force, residual[1]])?,
        ])?;
        let mobility = [v[4], v[5], v[6]];
        let column0 = solve(mobility, [k, 0.0])?;
        let column1 = solve(mobility, [0.0, k])?;
        let cross = 0.5 * column0[1] + 0.5 * column1[0];
        let resistance = [[column0[0], cross], [cross, column1[1]]];
        let angular_difference = finite(v[1] / self.inner_radius_m)?;
        let outer = self.helical_stress(1.0, stress, pressure_stress);
        let accepted_axial = [axial[0], finite(axial[0] + v[0])?];
        let accepted_angular = [angular[0], finite(angular[0] + angular_difference)?];
        check()?;
        Ok(CarreauHelicalProfile {
            annulus: *self,
            controls,
            inner_stress_pa: stress,
            pressure_stress_pa: pressure_stress,
            flow: CarreauHelicalFlow {
                pressure_drop_pa: pressure,
                requested_axial_wall_velocities_m_s: axial,
                accepted_axial_wall_velocities_m_s: accepted_axial,
                requested_wall_angular_velocities_rad_s: angular,
                accepted_wall_angular_velocities_rad_s: accepted_angular,
                relative_axial_velocity_m_s: v[0],
                relative_angular_velocity_rad_s: angular_difference,
                wall_velocity_residual_m_s: residual,
                volumetric_flow_m3_s: flow,
                relative_volume_flow_m3_s: relative_flow,
                wall_tractions_pa: [[-stress[0], -stress[1]], outer],
                axial_wall_forces_n: [inner_force, outer_force],
                wall_torques_n_m: [-torque, torque],
                mechanical_power_w: power,
                dissipated_power_w: v[3],
                work_residual_w: finite(power - v[3])?,
                wall_resistance_n_s_m: resistance,
                maximum_shear_rate_per_s: evaluation.maximum_rate,
                quadrature_relative_change: evaluation.change,
                newton_iterations: iterations,
                backtracks,
                constitutive_queries: queries,
            },
        })
    }
}

impl CarreauHelicalProfile {
    pub const fn annulus(&self) -> &LaminarCarreauAnnulus {
        &self.annulus
    }
    pub const fn flow(&self) -> &CarreauHelicalFlow {
        &self.flow
    }
    pub fn sample_gap_fraction(
        &self,
        fraction: f64,
    ) -> Result<CarreauHelicalSample, MaterialError> {
        self.sample_gap_fraction_with_check(fraction, || Ok(()))
    }
    pub fn sample_gap_fraction_with_check(
        &self,
        fraction: f64,
        mut check: impl FnMut() -> Result<(), MaterialError>,
    ) -> Result<CarreauHelicalSample, MaterialError> {
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(parameters(
                "helical sample must lie in the physical annular gap",
            ));
        }
        check()?;
        let a = &self.annulus;
        let x = if fraction == 1.0 {
            1.0
        } else {
            (log_ratio_from_gap(a.inner_radius_m, fraction * a.gap_m) / a.log_ratio).clamp(0.0, 1.0)
        };
        let mut budget = Budget {
            used: 0,
            maximum: self.controls.maximum_constitutive_queries,
        };
        let (tau, local) = a.helical_local(
            x,
            self.inner_stress_pa,
            self.pressure_stress_pa,
            self.controls,
            &mut budget,
            &mut check,
        )?;
        let magnitude = tau[0].hypot(tau[1]);
        let gradient = if magnitude == 0.0 {
            [0.0; 2]
        } else {
            tau.map(|v| (v / magnitude) * local.rate)
        };
        let (axial, angular, change) = if fraction == 0.0 {
            (
                self.flow.accepted_axial_wall_velocities_m_s[0],
                self.flow.accepted_wall_angular_velocities_rad_s[0],
                0.0,
            )
        } else if fraction == 1.0 {
            (
                self.flow.accepted_axial_wall_velocities_m_s[1],
                self.flow.accepted_wall_angular_velocities_rad_s[1],
                0.0,
            )
        } else {
            let from_inner = x <= 0.5;
            let bounds = if from_inner { [0.0, x] } else { [x, 1.0] };
            let evaluated = a.helical_integral(
                self.inner_stress_pa,
                self.pressure_stress_pa,
                bounds,
                self.controls,
                &mut budget,
                &mut check,
            )?;
            let sign = if from_inner { 1.0 } else { -1.0 };
            let wall = usize::from(!from_inner);
            (
                finite(
                    self.flow.accepted_axial_wall_velocities_m_s[wall] + sign * evaluated.values[0],
                )?,
                finite(
                    self.flow.accepted_wall_angular_velocities_rad_s[wall]
                        + sign * evaluated.values[1] / a.inner_radius_m,
                )?,
                evaluated.change,
            )
        };
        let radius = a.gap_m.mul_add(fraction, a.inner_radius_m);
        check()?;
        Ok(CarreauHelicalSample {
            axial_velocity_m_s: axial,
            angular_velocity_rad_s: angular,
            azimuthal_velocity_m_s: product([radius, angular])?,
            shear_rate_components_per_s: gradient,
            shear_rate_per_s: local.rate,
            radial_shear_stress_pa: tau,
            dynamic_viscosity_pa_s: local.viscosity,
            dissipation_w_m3: product([magnitude, local.rate])?,
            quadrature_relative_change: change,
        })
    }
}

// Shared-node nested Simpson integration. Signed channels use their positive
// absolute integral as the error scale, so zero flux/cross-mobility is admissible.
// All quadrature weights remain positive, including the heat channel.
fn integrate<const N: usize>(
    node: &mut impl FnMut(f64) -> Result<[f64; N], MaterialError>,
    bounds: [f64; 2],
    controls: CarreauPipeControls,
) -> Result<([f64; N], f64), MaterialError> {
    let [lower, upper] = bounds;
    let width = upper - lower;
    let first = node(lower)?;
    let last = node(upper)?;
    let mut trapezoid = [0.0; N];
    let mut absolute = [0.0; N];
    for i in 0..N {
        trapezoid[i] = 0.5 * finite(first[i])? + 0.5 * finite(last[i])?;
        absolute[i] = 0.5 * first[i].abs() + 0.5 * last[i].abs();
    }
    let mut previous: Option<[f64; N]> = None;
    let mut agreements = 0;
    for level in 1..=controls.maximum_quadrature_refinements {
        let panels = 1_u32 << level;
        let h = 1.0 / f64::from(panels);
        let (mut added, mut magnitudes, mut compensation) = ([0.0; N], [0.0; N], [0.0; N]);
        for index in (1..panels).step_by(2) {
            let value = node(width.mul_add(f64::from(index) * h, lower))?;
            for i in 0..N {
                let adjusted = finite(value[i])? - compensation[i];
                let next = finite(added[i] + adjusted)?;
                compensation[i] = (next - added[i]) - adjusted;
                added[i] = next;
                magnitudes[i] = finite(magnitudes[i] + value[i].abs())?;
            }
        }
        let mut simpson = [0.0; N];
        let mut change = 0.0_f64;
        for i in 0..N {
            simpson[i] = finite((trapezoid[i] + 4.0 * h * added[i]) / 3.0)?;
            let scale = finite((absolute[i] + 4.0 * h * magnitudes[i]) / 3.0)?;
            if let Some(old) = previous {
                let defect = (simpson[i] - old[i]).abs();
                let relative = if defect == 0.0 { 0.0 } else { defect / scale };
                finite(relative)?;
                change = change.max(relative);
            }
        }
        if previous.is_some() && change <= controls.relative_quadrature_tolerance {
            agreements += 1;
            if agreements == 2 {
                for value in &mut simpson {
                    *value = finite(*value * width)?;
                }
                return Ok((simpson, change));
            }
        } else {
            agreements = 0;
        }
        previous = Some(simpson);
        for i in 0..N {
            trapezoid[i] = 0.5 * trapezoid[i] + h * added[i];
            absolute[i] = 0.5 * absolute[i] + h * magnitudes[i];
        }
    }
    Err(state_error("helical Carreau quadrature budget exhausted"))
}
