//! Complete-spectrum cubic forcing and exact continuous displacement moments.
//!
//! Forces are F0+F1*s+F2*s^2+F3*s^3, s=t/h in [0,1]. Every coefficient is
//! in newtons. No temporal quadrature, stiffness inverse, or mode truncation is
//! used. The returned normalized moments integrate s^k*u(h*s) for k=0,1,2;
//! these close work from polynomial sources and accelerating prescribed supports.

use super::*;

/// Accepted motion and original-coordinate continuous trajectory observations.
#[derive(Clone, Debug)]
pub struct PolynomialSolidStep {
    pub motion: ExactLinearSolidStep,
    /// I[k] = integral_0^1 s^k u(h*s) ds [m], k=0,1,2.
    /// Multiply by h for an integral in physical time. I[0] is mean position.
    pub displacement_moments_m: [Vec<f64>; 3],
    /// Endpoint acceleration from the admitted complete elastic pencil [m/s²].
    pub acceleration_mps2: Vec<f64>,
}

/// Read-only candidate consumed only after its enclosing physical owner admits
/// every additional reaction/work check. It cannot be constructed externally.
pub(crate) struct PolynomialTrial {
    displacement: Vec<f64>,
    velocity: Vec<f64>,
    modal_displacement: Vec<f64>,
    modal_velocity: Vec<f64>,
    cumulative_work_j: f64,
    correction_j: f64,
    output: PolynomialSolidStep,
}
impl PolynomialTrial {
    pub(crate) fn displacement_m(&self) -> &[f64] {
        &self.displacement
    }
    pub(crate) fn velocity_mps(&self) -> &[f64] {
        &self.velocity
    }
    pub(crate) fn displacement_moments_m(&self) -> &[Vec<f64>; 3] {
        &self.output.displacement_moments_m
    }
    pub(crate) fn acceleration_mps2(&self) -> &[f64] {
        &self.output.acceleration_mps2
    }
    pub(crate) fn motion(&self) -> &ExactLinearSolidStep {
        &self.output.motion
    }
}

// Columns multiply q0, v0 and the four mass-normalized source coefficients.
pub(super) struct PolynomialOscillator {
    endpoint: [[f64; 6]; 2],
    moments: [[f64; 6]; 3],
    lambda: f64,
}
impl PolynomialOscillator {
    fn prepare(lambda: f64, h: f64) -> Result<Self> {
        let z = lambda.sqrt() * h;
        let z2 = z * z;
        if !z2.is_finite() {
            return Err(invalid("unrepresentable cubic oscillator phase"));
        }
        // E_j(z) = sum_m (-z²)^m / (2m+j)!. The derivatives of
        // s^j E_j(zs) are s^(j-1) E_(j-1)(zs). A broad series region avoids
        // cancellation from repeated inverse-z² recurrences up through E8.
        let mut factorial = [1.0; 9];
        for j in 1..=8 {
            factorial[j] = factorial[j - 1] * j as f64;
        }
        let mut e = [0.0; 9];
        if z <= 4.0 {
            for j in 0..=8 {
                let mut term = 1.0 / factorial[j];
                let (mut sum, mut correction) = (term, 0.0);
                // At |z|=4 the first omitted term is below 3e-51 for E0.
                // This fixed bound is included in the cold product envelope.
                for m in 1..32 {
                    term *= -z2 / ((2 * m + j - 1) * (2 * m + j)) as f64;
                    add(&mut sum, &mut correction, term);
                }
                e[j] = sum;
            }
        } else {
            let (sine, cosine) = det::sin_cos(z);
            e[0] = cosine;
            e[1] = sine / z;
            for j in 0..=6 {
                e[j + 2] = (1.0 / factorial[j] - e[j]) / z2;
            }
        }
        let mut endpoint = [[0.0; 6]; 2];
        endpoint[0][0] = e[0];
        endpoint[0][1] = h * e[1];
        endpoint[1][0] = -(h * lambda) * e[1];
        endpoint[1][1] = e[0];
        for r in 0..4 {
            endpoint[0][r + 2] = h * (h * (factorial[r] * e[r + 2]));
            endpoint[1][r + 2] = h * (factorial[r] * e[r + 1]);
        }
        let mut moments = [[0.0; 6]; 3];
        for (k, row) in moments.iter_mut().enumerate() {
            // Repeated integration by parts, exact for k<=2:
            // integral s^k*s^p E_p(zs) ds = E_(p+1)-k E_(p+2)+k(k-1)E_(p+3).
            let moment = |p: usize| {
                e[p + 1] - k as f64 * e[p + 2] + (k * k.saturating_sub(1)) as f64 * e[p + 3]
            };
            row[0] = moment(0);
            row[1] = h * moment(1);
            for r in 0..4 {
                row[r + 2] = h * (h * (factorial[r] * moment(r + 2)));
            }
        }
        if endpoint
            .iter()
            .flatten()
            .chain(moments.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err(invalid("unrepresentable cubic oscillator coefficients"));
        }
        Ok(Self {
            endpoint,
            moments,
            lambda,
        })
    }
}

impl CompleteSpectralSolid {
    /// Prepare exact cubic forcing with all four source coefficients in the
    /// normalized interval coordinate s=t/h. Basis and original-state admission
    /// is identical to prepare_with_basis; existing affine arithmetic is kept.
    ///
    /// Explicit envelopes are 64*n²+2048*n simultaneous bytes,
    /// 4*n³+32*n²+4096*n cold products, and 24*n²+1024*n interval products.
    /// Extra coefficients are cached once; no large exponential is introduced.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_cubic_forces_with_basis(
        assembly: &TetElasticAssembly,
        eigenvalues_rad2_s2: &[f64],
        basis_columns: &[f64],
        dt_s: f64,
        initial_displacement_m: &[f64],
        initial_velocity_mps: &[f64],
        controls: CompleteSpectralSolidControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        let n = assembly.free_dofs.len();
        let _ = n
            .checked_mul(n)
            .ok_or_else(|| invalid("cubic shape overflow"))?;
        let nn = n as u128;
        bound(n as u128, controls.maximum_free_dofs as u128)?;
        bound(
            64 * nn * nn + 2048 * nn,
            controls.maximum_preparation_bytes as u128,
        )?;
        bound(
            4 * nn * nn * nn + 32 * nn * nn + 4096 * nn,
            controls.maximum_preparation_multiply_terms as u128,
        )?;
        bound(
            24 * nn * nn + 1024 * nn,
            controls.maximum_step_multiply_terms as u128,
        )?;
        let mut solid = Self::prepare_with_basis(
            assembly,
            eigenvalues_rad2_s2,
            basis_columns,
            dt_s,
            initial_displacement_m,
            initial_velocity_mps,
            controls,
            cx,
        )?;
        let mut coefficients = Vec::new();
        coefficients
            .try_reserve_exact(n)
            .map_err(|_| invalid("cubic coefficient allocation refused"))?;
        for &lambda in eigenvalues_rad2_s2 {
            poll(cx)?;
            coefficients.push(PolynomialOscillator::prepare(lambda, dt_s)?);
        }
        poll(cx)?;
        solid.polynomial_oscillators = Some(coefficients);
        Ok(solid)
    }

    /// Advance F(s)=F0+F1*s+F2*s²+F3*s³ in complete free-coordinate order.
    /// Coefficients are forces [N], not time derivatives. Finite force jumps
    /// between accepted intervals introduce no impulse. State and observations
    /// are published together after the independent original-energy work gate.
    pub fn advance_cubic(
        &mut self,
        force_coefficients_n: [&[f64]; 4],
        cx: &Cx<'_>,
    ) -> Result<PolynomialSolidStep> {
        let trial = self.trial_polynomial(force_coefficients_n, cx)?;
        Ok(self.commit_polynomial(trial))
    }

    pub(crate) fn trial_polynomial(
        &self,
        coefficients: [&[f64]; 4],
        cx: &Cx<'_>,
    ) -> Result<PolynomialTrial> {
        poll(cx)?;
        let cache = self
            .polynomial_oscillators
            .as_ref()
            .ok_or_else(|| invalid("cubic forcing requires explicit cubic preparation"))?;
        let n = self.free_dofs.len();
        if coefficients
            .iter()
            .any(|c| c.len() != n || c.iter().any(|v| !v.is_finite()))
        {
            return Err(invalid(
                "four complete finite cubic force coefficients required",
            ));
        }
        let step = self
            .steps
            .checked_add(1)
            .filter(|s| *s <= self.controls.maximum_steps)
            .ok_or_else(|| invalid("accepted step budget exhausted"))?;
        let time_s = step as f64 * self.dt_s;
        if !time_s.is_finite() || time_s <= self.time_s() {
            return Err(invalid("representable later clock required"));
        }
        let forces = [
            project(&self.basis, coefficients[0], cx)?,
            project(&self.basis, coefficients[1], cx)?,
            project(&self.basis, coefficients[2], cx)?,
            project(&self.basis, coefficients[3], cx)?,
        ];
        let (mut q, mut v, mut acceleration) = (zeros(n)?, zeros(n)?, zeros(n)?);
        let mut moments = [zeros(n)?, zeros(n)?, zeros(n)?];
        for i in 0..n {
            poll(cx)?;
            let a = &cache[i];
            let state = [
                self.modal_displacement[i],
                self.modal_velocity[i],
                forces[0][i],
                forces[1][i],
                forces[2][i],
                forces[3][i],
            ];
            let dot = |row: &[f64; 6]| (0..6).fold(0.0, |s, j| row[j].mul_add(state[j], s));
            q[i] = dot(&a.endpoint[0]);
            v[i] = dot(&a.endpoint[1]);
            for k in 0..3 {
                moments[k][i] = dot(&a.moments[k]);
            }
            let end_force = (0..4).fold(0.0, |s, r| s + forces[r][i]);
            acceleration[i] = (-a.lambda).mul_add(q[i], end_force);
        }
        let displacement = reconstruct(&self.basis, &q, cx)?;
        let velocity = reconstruct(&self.basis, &v, cx)?;
        let acceleration_mps2 = reconstruct(&self.basis, &acceleration, cx)?;
        let displacement_moments_m = [
            reconstruct(&self.basis, &moments[0], cx)?,
            reconstruct(&self.basis, &moments[1], cx)?,
            reconstruct(&self.basis, &moments[2], cx)?,
        ];
        let (mut supplied_work_j, mut roundoff) = (0.0, 0.0);
        let mut mean_velocity_mps = zeros(n)?;
        for i in 0..n {
            poll(cx)?;
            let du = displacement[i] - self.displacement[i];
            // Integral F du = F0*Delta u + sum_r Fr*(u1-r*I_(r-1)).
            // Every I is a continuous trajectory moment, independent of H.
            add(&mut supplied_work_j, &mut roundoff, coefficients[0][i] * du);
            for r in 1..4 {
                add(
                    &mut supplied_work_j,
                    &mut roundoff,
                    coefficients[r][i]
                        * (displacement[i] - r as f64 * displacement_moments_m[r - 1][i]),
                );
            }
            mean_velocity_mps[i] = du / self.dt_s;
        }
        if q.iter()
            .chain(&v)
            .chain(&mean_velocity_mps)
            .any(|x| !x.is_finite())
        {
            return Err(invalid("unrepresentable cubic propagated state"));
        }
        let stored_after_j = quadratic(&self.stiffness, &displacement, cx)?;
        let kinetic_after_j = quadratic(&self.mass, &velocity, cx)?;
        let change = quadratic_change(&self.stiffness, &self.displacement, &displacement, cx)?
            + quadratic_change(&self.mass, &self.velocity, &velocity, cx)?;
        let (mut cumulative_work_j, mut correction_j) =
            (self.cumulative_work_j, self.work_correction_j);
        add(&mut cumulative_work_j, &mut correction_j, supplied_work_j);
        let energy = ExactLinearSolidEnergy {
            stored_before_j: self.stored_j,
            stored_after_j,
            kinetic_before_j: self.kinetic_j,
            kinetic_after_j,
            supplied_work_j,
            dissipated_j: 0.0,
            cumulative_dissipated_j: 0.0,
            balance_error_j: supplied_work_j - change,
            cumulative_balance_error_j: cumulative_work_j
                - (stored_after_j + kinetic_after_j - self.initial_energy_j),
        };
        for defect_j in [energy.balance_error_j, energy.cumulative_balance_error_j] {
            if !defect_j.is_finite() || defect_j.abs() > self.controls.maximum_work_error_j {
                return Err(CompleteSpectralSolidError::Energy {
                    defect_j,
                    allowance_j: self.controls.maximum_work_error_j,
                });
            }
        }
        if !cumulative_work_j.is_finite() || !correction_j.is_finite() {
            return Err(invalid("unrepresentable cubic accumulated work"));
        }
        poll(cx)?;
        Ok(PolynomialTrial {
            displacement,
            velocity,
            modal_displacement: q,
            modal_velocity: v,
            cumulative_work_j,
            correction_j,
            output: PolynomialSolidStep {
                motion: ExactLinearSolidStep {
                    step,
                    time_s,
                    mean_velocity_mps,
                    energy,
                },
                displacement_moments_m,
                acceleration_mps2,
            },
        })
    }

    // Only the public immediate advance or the enclosing support owner may
    // consume a trial, with exclusive access and no intervening state change.
    // All fallible work and the final Cx poll precede this publication boundary.
    pub(crate) fn commit_polynomial(&mut self, trial: PolynomialTrial) -> PolynomialSolidStep {
        debug_assert_eq!(self.steps + 1, trial.output.motion.step);
        self.displacement = trial.displacement;
        self.velocity = trial.velocity;
        self.modal_displacement = trial.modal_displacement;
        self.modal_velocity = trial.modal_velocity;
        self.steps = trial.output.motion.step;
        self.stored_j = trial.output.motion.energy.stored_after_j;
        self.kinetic_j = trial.output.motion.energy.kinetic_after_j;
        self.cumulative_work_j = trial.cumulative_work_j;
        self.work_correction_j = trial.correction_j;
        trial.output
    }
}
