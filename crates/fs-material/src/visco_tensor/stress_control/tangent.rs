//! Exact load derivatives, retaining original branch stress reconstruction.

use super::super::WEIGHTS;
use super::{
    Dd, PreparedTensorStressControl, PronyRamp, Result, Tangent6, apply, dd_dot, finite, invalid,
    modal_integrals, poll, zeros,
};

/// Derivative of the physical ramp endpoint with respect to its six controls.
///
/// `strain[i][j]` and `stress[i][j]` map a change in prescribed endpoint
/// component j to resulting physical strain/stress component i. A column's
/// input is stress [Pa] or strain according to the prepared control partition.
/// All input and output shears are tensor shears, never engineering shears.
/// Stress-controlled stress rows and strain-controlled strain rows are the
/// corresponding identity rows, including at zero duration.
///
/// Initial strain, original branch memory, material coefficients and duration
/// are fixed. This matrix is independent of the initial history because the
/// original constitutive law is linear. It is not a history adjoint or a
/// material/time derivative, and does not admit a finite load outside the
/// original strain domain. Native factors and floating-point admission are
/// shared with the primal; no forward-error enclosure is claimed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorControlTangent {
    /// Endpoint strain derivative in physical tensor-shear coordinates.
    pub strain: Tangent6,
    /// Endpoint total stress derivative, including strain-control reactions.
    pub stress: Tangent6,
}

impl PreparedTensorStressControl {
    /// Exact endpoint stress/strain derivative for a prescribed stress ramp.
    ///
    /// The original initial history and coefficients remain fixed. This uses
    /// the same compliance factors and analytical relaxation responses as the
    /// primal; no perturbed history advances or numerical differentiation run.
    /// The small six-column map can be cached by a caller using a fixed clock.
    ///
    /// # Errors
    /// Invalid duration, unrepresentable response, unresolved original-source
    /// reconstruction, or cancellation. No history is advanced.
    #[allow(clippy::too_many_lines)] // One bounded exact source-to-endpoint map.
    pub fn ramp_tangent(
        &self,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorControlTangent> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(invalid("finite nonnegative tangent duration required"));
        }
        let n = self.rates.len();
        let mut integrals = zeros(n)?;
        let mut ramps = Vec::new();
        ramps
            .try_reserve_exact(n)
            .map_err(|_| invalid("load tangent allocation refused"))?;
        for (mode, &rate) in self.rates.iter().enumerate() {
            poll(&mut cancelled)?;
            let r = duration_s * rate;
            if !r.is_finite() || (duration_s > 0.0 && r == 0.0) {
                return Err(invalid("unrepresentable tangent relaxation number"));
            }
            let ramp = PronyRamp::new(r);
            if duration_s > 0.0 {
                let (_, mean, _) = ramp.apply(0.0, 1.0);
                integrals[mode] = modal_integrals(0.0, 1.0, rate, duration_s, mean, &ramp).0;
                if !integrals[mode].is_finite() {
                    return Err(invalid("unrepresentable integrated load tangent"));
                }
            }
            ramps.push(ramp);
        }
        let mut strain = [[0.0; 6]; 6];
        let mut stress = [[0.0; 6]; 6];
        let mut dq = vec![[0.0; 6]; self.model.branches().len()];
        let tolerance = 8192.0 * f64::EPSILON * (n + 6) as f64;
        for col in 0..6 {
            poll(&mut cancelled)?;
            let strain_input = core::array::from_fn(|i| {
                if i == col && !self.stress_components[col] {
                    1.0
                } else {
                    0.0
                }
            });
            let coupled = apply(&self.instantaneous_stiffness, &strain_input);
            let rhs = core::array::from_fn(|i| {
                (if i == col && self.stress_components[col] {
                    1.0
                } else {
                    0.0
                }) - coupled[i]
            });
            let elastic = self.compliance(rhs)?;
            let mut de =
                core::array::from_fn::<_, 6, _>(|i| Dd::from_f64(elastic[i] + strain_input[i]));
            for (mode, &integral) in integrals.iter().enumerate() {
                for (axis, value) in de.iter_mut().enumerate() {
                    if self.stress_components[axis] {
                        *value = *value
                            + Dd::from_f64(self.forcing[mode][axis])
                                * Dd::from_f64(self.forcing[mode][col])
                                * Dd::from_f64(integral);
                    }
                }
            }
            let de = de.map(Dd::to_f64);
            if !finite(&de) {
                return Err(invalid("unrepresentable strain load tangent"));
            }
            for row in 0..n {
                poll(&mut cancelled)?;
                dq[row / 6][row % 6] = if duration_s == 0.0 {
                    let branch = &self.model.branches()[row / 6];
                    let b = branch.resolved_gram_rows(&self.model.factors[row / 6]);
                    dd_dot(&b[row % 6], &de).to_f64()
                } else {
                    let mut value = Dd::ZERO;
                    for (mode, ramp) in ramps.iter().enumerate() {
                        let coefficient =
                            self.modes[row * n + mode] / self.branch_rate_roots[row / 6];
                        value = value
                            + Dd::from_f64(
                                ramp.spring_increment(0.0, coefficient * self.forcing[mode][col]),
                            );
                    }
                    value.to_f64()
                };
            }
            // Reconstruct stress from equilibrium and ORIGINAL memory. In a
            // fully relaxed stiff limit C0 minus a branch correction can lose
            // a weak equilibrium spring; positive source contributions retain it.
            let (ds, magnitude) = self.source_stress(&de, &dq)?;
            for axis in 0..6 {
                let identity = if axis == col { 1.0 } else { 0.0 };
                if self.stress_components[axis]
                    && (ds[axis] - identity).abs() > tolerance * (magnitude[axis] + identity)
                {
                    return Err(invalid("original-source load tangent is unresolved"));
                }
                let scale = WEIGHTS[col] / WEIGHTS[axis];
                strain[axis][col] = scale * de[axis];
                stress[axis][col] = if self.stress_components[axis] {
                    identity
                } else {
                    scale * ds[axis]
                };
                if !strain[axis][col].is_finite() || !stress[axis][col].is_finite() {
                    return Err(invalid("unrepresentable physical load tangent"));
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(TensorControlTangent { strain, stress })
    }
}
