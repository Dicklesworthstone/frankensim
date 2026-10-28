//! Full fixed-state viscous forces inside the existing modal contact/flow solve.
//! The shared linear algebra owner factors M + h C + h² K once, where h=dt/2.
//! Pressure and contact use that SAME response, never a lagged damping force.
use super::{
    AcousticRealizeError, MAX_MODES, ModalAperture, ModalApertureSpec, invalid, nonlinear,
};
use crate::plate_viscosity::ModalPlateViscosity;
use fs_exec::CancelGate;
use fs_la::factor::{Cholesky, cholesky};

pub(super) struct ApertureViscosity {
    matrix: Vec<f64>,
    factor: Option<Cholesky>,
    response: Vec<f64>,
    scale: Vec<f64>,
    pressure: Vec<f64>,
    contact: Vec<f64>,
    passive_columns: Vec<Vec<f64>>,
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() {
        Err(AcousticRealizeError::Cancelled)
    } else {
        Ok(())
    }
}
impl ApertureViscosity {
    fn new(
        spec: &ModalApertureSpec,
        matrix: Vec<f64>,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let n = spec.modes.len();
        if n == 0
            || n > MAX_MODES
            || matrix.len() != n * n
            || spec.modes.iter().any(|m| m.damping_n_s_m != 0.0)
            || matrix.iter().any(|c| !c.is_finite())
        {
            return Err(invalid(
                "full modal viscosity needs one finite square matrix and zero separate diagonal loss",
            ));
        }
        for i in 0..n {
            for j in 0..n {
                if matrix[n * i + i] < 0.0
                    || matrix[n * i + j] != matrix[n * j + i]
                    || (matrix[n * i + i] == 0.0 && matrix[n * i + j] != 0.0)
                {
                    return Err(invalid(
                        "modal viscosity must be symmetric with passive diagonal and explicit zero rows",
                    ));
                }
            }
        }
        // Explicitly undamped coordinates are allowed (e.g. one lossless blade).
        // Coupled nullspaces also occur in localized/interface damping. Retain
        // the original SPD arithmetic; the shared PSD owner admits those laws
        // without a loss floor or an independent modal decomposition.
        let active: Vec<_> = (0..n).filter(|&i| matrix[n * i + i] > 0.0).collect();
        if active.is_empty() {
            return Ok(Self {
                matrix,
                factor: None,
                response: vec![],
                scale: vec![],
                pressure: vec![],
                contact: vec![],
                passive_columns: vec![],
            });
        }
        let mut normalized = Vec::with_capacity(active.len() * active.len());
        for &i in &active {
            for &j in &active {
                normalized
                    .push(matrix[n * i + j] / matrix[n * i + i].sqrt() / matrix[n * j + j].sqrt());
            }
        }
        let passive_columns = if cholesky(&normalized, active.len()).is_ok() {
            vec![]
        } else {
            // The shared owner's C is mass-normalized [1/s]. Transform back
            // before evaluating physical work, and check every original entry:
            // its rank tolerance must not discard a weak physical coordinate.
            let root_mass: Vec<_> = spec.modes.iter().map(|m| m.mass_kg.sqrt()).collect();
            let per_s: Vec<_> = (0..n * n)
                .map(|ij| matrix[ij] / root_mass[ij / n] / root_mass[ij % n])
                .collect();
            let passive = ModalPlateViscosity::from_psd_modal_matrix(per_s, n, 1.0, gate)?;
            let columns: Vec<Vec<f64>> = passive.channel_shapes().iter()
                .map(|b| b.iter().zip(&root_mass).map(|(b, m)| b * m).collect())
                .collect();
            for i in 0..n {
                for j in 0..n {
                    let actual: f64 = columns.iter().map(|b| b[i] * b[j]).sum();
                    let allowance = 1e-10 * matrix[n * i + i].sqrt()
                        * matrix[n * j + j].sqrt();
                    if !actual.is_finite() || !allowance.is_finite()
                        || (actual - matrix[n * i + j]).abs() > allowance
                    {
                        return Err(invalid(
                            "passive channels do not reconstruct the original physical viscosity",
                        ));
                    }
                }
            }
            columns
        };
        let h = 0.5 * spec.time_step_s;
        let scale: Vec<_> = spec
            .modes
            .iter()
            .enumerate()
            .map(|(i, m)| (m.mass_kg + h * matrix[n * i + i] + h * h * m.stiffness_n_m).sqrt())
            .collect();
        if scale.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(invalid("full viscous midpoint diagonal is unrepresentable"));
        }
        let mut response = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                response[n * i + j] = if i == j {
                    1.0
                } else {
                    h * matrix[n * i + j] / scale[i] / scale[j]
                };
            }
        }
        if response.iter().any(|x| !x.is_finite()) {
            return Err(invalid("full viscous midpoint matrix is unrepresentable"));
        }
        let factor = cholesky(&response, n).map_err(nonlinear)?;
        let mut result = Self {
            matrix,
            factor: Some(factor),
            response,
            scale,
            pressure: spec.modes.iter().map(|m| h * m.pressure_area_m2).collect(),
            contact: spec.modes.iter().map(|m| h * m.gap_projection).collect(),
            passive_columns,
        };
        let solve = |rhs: &mut [f64]| {
            solve_scaled(
                result.factor.as_ref().expect("active factor"),
                &result.response,
                &result.scale,
                rhs,
            )
        };
        solve(&mut result.pressure)?;
        solve(&mut result.contact)?;
        poll(gate)?;
        Ok(result)
    }
    pub(super) fn is_active(&self) -> bool {
        self.factor.is_some()
    }
    pub(super) fn pressure_response(&self) -> &[f64] {
        &self.pressure
    }

    /// Reuse the same full-C mechanical factor for changed acoustic faces.
    /// This updates traction and its condensed swept-flow feedback together.
    pub(super) fn rebind_pressure(&mut self, spec: &ModalApertureSpec)
        -> Result<(), AcousticRealizeError>
    {
        if !self.is_active() { return Ok(()); }
        let h = 0.5 * spec.time_step_s;
        let mut pressure: Vec<_> = spec.modes.iter().map(|m| h * m.pressure_area_m2).collect();
        self.solve_free(&mut pressure)?;
        self.pressure = pressure;
        Ok(())
    }
    pub(super) fn contact_response(&self) -> &[f64] {
        &self.contact
    }
    pub(super) fn solve_free(&self, rhs: &mut [f64]) -> Result<(), AcousticRealizeError> {
        solve_scaled(
            self.factor.as_ref().expect("active factor"),
            &self.response,
            &self.scale,
            rhs,
        )
    }
    pub(super) fn force_and_power(
        &self,
        velocity: &[f64],
        force: &mut [f64],
    ) -> Result<f64, AcousticRealizeError> {
        for (row, out) in self
            .matrix
            .chunks_exact(velocity.len())
            .zip(force.iter_mut())
        {
            *out = row.iter().zip(velocity).map(|(c, v)| c * v).sum();
        }
        let power: f64 = if self.passive_columns.is_empty() {
            force.iter().zip(velocity).map(|(f, v)| f * v).sum()
        } else {
            // Sum passive channel powers, so cancellation in a coupled
            // undamped subspace cannot become negative dissipated work.
            self.passive_columns.iter().map(|b| {
                let speed: f64 = b.iter().zip(velocity).map(|(b, v)| b * v).sum();
                speed * speed
            }).sum()
        };
        if force.iter().any(|x| !x.is_finite()) || !power.is_finite() || power < 0.0 {
            return Err(invalid("full viscous force or passive power is unresolved"));
        }
        Ok(power)
    }
}
fn solve_scaled(
    factor: &Cholesky,
    matrix: &[f64],
    scale: &[f64],
    rhs: &mut [f64],
) -> Result<(), AcousticRealizeError> {
    for (v, s) in rhs.iter_mut().zip(scale) {
        *v /= s;
    }
    let mut before = [0.0; MAX_MODES];
    before[..rhs.len()].copy_from_slice(rhs);
    factor.solve(rhs);
    for (i, row) in matrix.chunks_exact(rhs.len()).enumerate() {
        let actual: f64 = row.iter().zip(rhs.iter()).map(|(a, x)| a * x).sum();
        let magnitude = before[i].abs()
            + row
                .iter()
                .zip(rhs.iter())
                .map(|(a, x)| (a * x).abs())
                .sum::<f64>();
        if !magnitude.is_finite()
            || !actual.is_finite()
            || (actual - before[i]).abs()
                > 128.0 * f64::EPSILON * rhs.len() as f64 * magnitude.max(f64::MIN_POSITIVE)
        {
            return Err(invalid(
                "full viscous linear response has an unresolved backward residual",
            ));
        }
    }
    for (v, s) in rhs.iter_mut().zip(scale) {
        *v /= s;
    }
    if rhs.iter().any(|v| !v.is_finite()) {
        return Err(invalid("full viscous modal response overflowed"));
    }
    Ok(())
}
impl ModalAperture {
    pub(super) fn admit_linear_clock(&self) -> Result<(), AcousticRealizeError> {
        // ExactZoh evolves the original complete structural pencil. Pressure
        // and contact remain held sample forces with their own refinement error.
        if self.exact.is_some() {
            return Ok(());
        }
        let modes = &self.spec.modes;
        let matrix = self.viscosity_matrix_n_s_m();
        let area_sum = modes
            .iter()
            .map(|m| m.pressure_area_m2.abs() / fs_math::det::sqrt(m.mass_kg))
            .sum::<f64>();
        for (i, m) in modes.iter().enumerate() {
            let root_mass = fs_math::det::sqrt(m.mass_kg);
            let resistance = if let Some(c) = matrix {
                modes
                    .iter()
                    .enumerate()
                    .map(|(j, other)| {
                        c[i * modes.len() + j].abs() / root_mass / fs_math::det::sqrt(other.mass_kg)
                    })
                    .sum::<f64>()
            } else {
                m.damping_n_s_m / m.mass_kg
            };
            // Full linear-generator norm, including characteristic-load Z A A^T.
            let rate = fs_math::det::sqrt(m.stiffness_n_m / m.mass_kg)
                + resistance
                + self.spec.impedance_pa_s_m3 * m.pressure_area_m2.abs() / root_mass * area_sum;
            fs_phs::linear_exact::admit_midpoint_band(rate, self.spec.time_step_s)
                .map_err(AcousticRealizeError::TimeResolution)?;
        }
        Ok(())
    }

    /// Select a complete constant C [N s/m], row-major in these physical modes.
    /// Midpoint uses -C v_mid and dt*v_mid^T C v_mid. ExactZoh includes C in
    /// the complete cached structural exponential and integrates v(t)^T C v(t)
    /// along that same held-force trajectory. Mass, stiffness and port projections,
    /// initial motion and connected acoustic models are not changed or retuned.
    ///
    /// This raw-matrix API confers no material identity. Material-driven callers
    /// must project C on the SAME eigenbasis, never recompute modes separately.
    /// Positive semidefinite laws include coupled undamped subspaces, such as
    /// a relative-motion dashpot. The shared passive-channel owner resolves
    /// singular laws in mass-normalized coordinates, with reconstruction of
    /// every physical C entry to 1e-10 of its diagonal geometric scale. Midpoint
    /// loss uses channel squares; the original C remains in force and response.
    /// Unresolved laws refuse; no artificial damping fills their nullspaces.
    /// Zero C uses the original diagonal-zero arithmetic. No new per-step matrix
    /// factorization or allocation is added. Memory is bounded by 32² entries
    /// per retained matrix; the shared Cholesky owner performs cold factor work.
    ///
    /// # Errors
    /// Late/repeated selection, separate diagonal loss, nonfinite/asymmetric/
    /// nonpassive C, unrepresentable response or cancellation. Consumes `self`
    /// but never advances physical time or publishes partial mechanical state.
    pub fn with_viscosity(
        mut self,
        matrix_n_s_m: Vec<f64>,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        if self.samples != 0 || self.viscosity.is_some() || self.relaxation.is_some() {
            return Err(invalid(
                "full modal viscosity requires one unadvanced attachment without a separate material spectrum",
            ));
        }
        self.viscosity = Some(ApertureViscosity::new(&self.spec, matrix_n_s_m, gate)?);
        if self.exact.is_some() {
            self.exact = Some(super::exact::ExactModes::new(
                &self.spec,
                self.viscosity_matrix_n_s_m(),
            )?);
            poll(gate)?;
        }
        self.admit_linear_clock()?;
        Ok(self)
    }
    /// Complete admitted C in the source coordinates; None means diagonal loss.
    /// This read does not infer measured-material authority for authored matrices.
    pub fn viscosity_matrix_n_s_m(&self) -> Option<&[f64]> {
        self.viscosity.as_ref().map(|c| c.matrix.as_slice())
    }
}
