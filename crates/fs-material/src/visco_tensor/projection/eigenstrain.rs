//! Integrated free strain in the original quadrature and factor coordinates.
use super::*;

/// Affine expansion coefficient and its derivative in one fixed Mandel frame.
/// This numerical descriptor does not establish material-source authority:
/// the source owner must admit the complete selected paths, frame and interval.
/// Actual integrated eigenstrain remains owned by the original strain model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TensorThermalExpansion {
    reference_temperature_k: f64,
    reference_alpha_mandel_per_k: Voigt,
    alpha_slope_mandel_per_k2: Voigt,
    temperature_k: f64,
    alpha_mandel_per_k: Voigt,
    interval_k: [f64; 2],
    maximum_abs_alpha_mandel_per_k: Voigt,
}

impl TensorThermalExpansion {
    pub fn try_new(
        reference_temperature_k: f64,
        alpha_mandel_per_k: Voigt,
        alpha_slope_mandel_per_k2: Voigt,
        interval_k: [f64; 2],
    ) -> Result<Self> {
        let [low, high] = interval_k;
        if !reference_temperature_k.is_finite()
            || !low.is_finite() || !high.is_finite() || low <= 0.0 || high <= low
            || !(low..=high).contains(&reference_temperature_k)
            || !finite(&alpha_mandel_per_k) || !finite(&alpha_slope_mandel_per_k2)
        {
            return Err(invalid("thermal expansion needs a finite bounded absolute-temperature law"));
        }
        let endpoint = |t| core::array::from_fn::<_, 6, _>(|i|
            alpha_slope_mandel_per_k2[i].mul_add(t - reference_temperature_k, alpha_mandel_per_k[i]));
        let cold = endpoint(low);
        let warm = endpoint(high);
        if !finite(&cold) || !finite(&warm) {
            return Err(invalid("thermal expansion interval derivative overflow"));
        }
        Ok(Self {
            reference_temperature_k,
            reference_alpha_mandel_per_k: alpha_mandel_per_k,
            alpha_slope_mandel_per_k2,
            temperature_k: reference_temperature_k,
            alpha_mandel_per_k,
            interval_k,
            maximum_abs_alpha_mandel_per_k: core::array::from_fn(|i| cold[i].abs().max(warm[i].abs())),
        })
    }

    /// Refresh the derivative snapshot without source queries or allocation.
    pub fn at_temperature(&self, temperature_k: f64) -> Result<Self> {
        if !temperature_k.is_finite() || !(self.interval_k[0]..=self.interval_k[1]).contains(&temperature_k) {
            return Err(invalid("thermal expansion derivative leaves its admitted interval"));
        }
        let alpha = core::array::from_fn(|i| self.alpha_slope_mandel_per_k2[i]
            .mul_add(temperature_k - self.reference_temperature_k, self.reference_alpha_mandel_per_k[i]));
        if !finite(&alpha) { return Err(invalid("thermal expansion derivative overflow")); }
        Ok(Self { temperature_k, alpha_mandel_per_k: alpha, ..*self })
    }

    pub const fn temperature_k(&self) -> f64 { self.temperature_k }
    pub const fn alpha_mandel_per_k(&self) -> Voigt { self.alpha_mandel_per_k }
    pub const fn alpha_slope_mandel_per_k2(&self) -> Voigt { self.alpha_slope_mandel_per_k2 }
    pub const fn temperature_interval_k(&self) -> [f64; 2] { self.interval_k }

    pub(super) fn admit(
        &self, derivative: &source::thermal::TensorTemperatureDerivative, temperature_k: f64,
    ) -> Result<()> {
        if self.temperature_k.to_bits() != temperature_k.to_bits()
            || self.interval_k != derivative.interval_k
        {
            return Err(invalid("thermal expansion needs the matching derivative temperature and source interval"));
        }
        Ok(())
    }
}

impl ProjectedTensorProny {
    /// Bind an explicitly admitted total stress-free strain in this point's
    /// Mandel frame. Geometry, reference mass and every original strain column
    /// remain unchanged. The caller owns its source and temperature history.
    /// A changing Kelvin eigenstrain needs a resolved thermal strain rate and
    /// is deliberately outside the endpoint Maxwell continuation model.
    pub fn with_eigenstrain_mandel(&self, free: Voigt) -> Result<Self> {
        self.model
            .validate_mandel_strain(&free)
            .map_err(|_| invalid("free strain exceeds original small-strain domain"))?;
        if free.iter().any(|x| *x != 0.0) && self.model.kv_pa_s().is_some() {
            return Err(invalid(
                "Kelvin thermal expansion requires a resolved strain-rate coupling",
            ));
        }
        let mut next = self.clone();
        next.free_strain_mandel = free;
        Ok(next)
    }

    /// Total stress-free strain in the same Mandel frame as the material.
    pub const fn free_strain_mandel(&self) -> &Voigt {
        &self.free_strain_mandel
    }

    /// Equilibrium affine load and constant free energy: psi_eq =
    /// .5 u^T K u - u^T load + constant. These are internal material storage,
    /// not an external force port or a source of dissipative heat.
    pub fn equilibrium_eigenstrain_terms(&self) -> Result<(Vec<f64>, f64)> {
        let stress = apply(&self.model.c_eq, &self.free_strain_mandel);
        let load: Vec<_> = self
            .columns
            .iter()
            .map(|column| self.volume_m3 * dot(&mandel(column), &stress))
            .collect();
        let constant = 0.5 * self.volume_m3 * dot(&self.free_strain_mandel, &stress);
        if load.iter().any(|v| !v.is_finite()) || !constant.is_finite() || constant < 0.0 {
            return Err(invalid("equilibrium eigenstrain storage overflow"));
        }
        Ok((load, constant))
    }

    /// Original volume-scaled branch factors applied to the total free strain.
    /// Subtract these from A*u only when explicitly initializing instantaneous
    /// memory. Existing histories must use transfer_factor_memory instead.
    pub fn eigenstrain_factor_coordinates(&self) -> Result<Vec<f64>> {
        let root = det::sqrt(self.volume_m3);
        let mut values = Vec::with_capacity(6 * self.model.branches.len());
        for (branch, factor) in self.model.branches.iter().zip(&self.model.factors) {
            values.extend(
                apply(&branch.resolved_gram_rows(factor), &self.free_strain_mandel)
                    .map(|v| root * v),
            );
        }
        if values.iter().any(|v| !v.is_finite()) {
            return Err(invalid("eigenstrain factor coordinates overflow"));
        }
        Ok(values)
    }

    /// Positive magnitude of -dU/dT / T at fixed total and viscous strains,
    /// for CONSTANT stored tensors and constant free-strain derivative alpha.
    /// The complete equilibrium and Maxwell tensors contribute, including
    /// constrained directions absent from a reduced displacement basis.
    pub fn expansion_capacity_coefficient_j_k2(&self, alpha_mandel: Voigt) -> Result<f64> {
        if !finite(&alpha_mandel) {
            return Err(invalid("nonfinite expansion derivative"));
        }
        let mut value = dot(&alpha_mandel, &apply(&self.model.c_eq, &alpha_mandel));
        for branch in &self.model.branches {
            value += dot(&alpha_mandel, &apply(&branch.c_pa, &alpha_mandel));
        }
        value *= self.volume_m3;
        if !value.is_finite() || value < 0.0 {
            return Err(invalid("thermal capacity contribution overflow"));
        }
        Ok(value)
    }

    /// Conservative upper bound on the mechanical psi_TT contribution [J/K^2]
    /// for affine stored tensors and constant stress-free derivative alpha.
    ///
    /// The caller owns the original source identity/frame and must supply the
    /// temperature at which THIS point's tensors were resolved. At fixed total
    /// and viscous strains,
    /// psi_TT/V = alpha:C_inst:alpha - 2 alpha:C_eq':e
    ///            - 2 sum_j alpha:C_j':q_j.
    /// Native strain admission bounds ||e|| by L. Every varying branch used
    /// with expansion is separately admitted with ||q_j|| <= 2L by the shared
    /// thermal-state evaluator. Absolute matrix arithmetic bounds these terms
    /// over the COMPLETE source interval, without sampling candidate strains.
    pub fn affine_expansion_capacity_bound_j_k2(
        &self,
        derivative: &source::thermal::TensorTemperatureDerivative,
        temperature_k: f64,
        alpha_mandel: Voigt,
    ) -> Result<f64> {
        if derivative.resolved_temperature_k.is_some() {
            return Err(invalid(
                "non-affine storage needs the general thermal capacity bound",
            ));
        }
        self.thermal_capacity_bound_j_k2(derivative, temperature_k, alpha_mandel)
    }

    /// Upper bound on psi_TT over the complete admitted source interval and
    /// native strain envelope, including zero-expansion nonlinear storage.
    /// Source-wide symmetric operator-norm bounds give
    /// .5 L^2 ||C_TT|| + 2 ||alpha|| L ||C_T|| + ||alpha||^2 ||C||.
    /// Affine tensors retain their original entrywise bound and constant-law
    /// limit. The receiving caloric owner must admit C_caloric - T*bound > 0;
    /// this value is not an extra caloric energy or dissipative heat source.
    pub fn thermal_capacity_bound_j_k2(
        &self,
        derivative: &source::thermal::TensorTemperatureDerivative,
        temperature_k: f64,
        alpha_mandel: Voigt,
    ) -> Result<f64> {
        self.capacity_bound_j_k2(derivative, temperature_k, alpha_mandel, None)
    }

    /// Bound the complete source interval for an affine expansion coefficient.
    /// Includes -alpha':C:e and every Maxwell -alpha':Cj:qj contribution.
    /// Constant expansion retains the original exact admission and arithmetic.
    pub fn thermal_capacity_bound_with_expansion_j_k2(
        &self,
        derivative: &source::thermal::TensorTemperatureDerivative,
        temperature_k: f64,
        expansion: &TensorThermalExpansion,
    ) -> Result<f64> {
        expansion.admit(derivative, temperature_k)?;
        let beta = expansion.alpha_slope_mandel_per_k2;
        if beta.iter().all(|v| *v == 0.0) {
            return self.thermal_capacity_bound_j_k2(derivative, temperature_k, expansion.alpha_mandel_per_k);
        }
        self.capacity_bound_j_k2(
            derivative, temperature_k, expansion.maximum_abs_alpha_mandel_per_k, Some(beta),
        )
    }

    fn capacity_bound_j_k2(
        &self,
        derivative: &source::thermal::TensorTemperatureDerivative,
        temperature_k: f64,
        alpha_mandel: Voigt,
        expansion_curvature: Option<Voigt>,
    ) -> Result<f64> {
        if !finite(&alpha_mandel)
            || !temperature_k.is_finite()
            || !(derivative.interval_k[0]..=derivative.interval_k[1]).contains(&temperature_k)
            || derivative
                .resolved_temperature_k
                .is_some_and(|t| t.to_bits() != temperature_k.to_bits())
            || derivative.branches_pa_k.len() != self.model.branches.len()
        {
            return Err(invalid(
                "thermal capacity needs the matching bounded source state",
            ));
        }
        let varying = derivative
            .equilibrium_pa_k
            .iter()
            .flatten()
            .chain(derivative.branches_pa_k.iter().flatten().flatten())
            .chain(derivative.equilibrium_pa_k2.iter().flatten())
            .any(|v| *v != 0.0)
            || derivative
                .equilibrium_bounds
                .as_ref()
                .is_some_and(|b| b.slope_norm_pa_k != 0.0 || b.curvature_norm_pa_k2 != 0.0);
        if !varying && expansion_curvature.is_none() {
            return self.expansion_capacity_coefficient_j_k2(alpha_mandel);
        }
        let [low, high] = derivative.interval_k;
        let alpha = alpha_mandel.map(f64::abs);
        let beta = expansion_curvature.unwrap_or([0.0; 6]).map(f64::abs);
        let strain_limit = self.model.maximum_strain_norm();
        let bound = |tensor: &Tangent6, slope: &Tangent6, strain_bound: f64| -> Result<f64> {
            let mut storage = 0.0;
            let mut cross = 0.0;
            let mut expansion_curve = 0.0;
            for i in 0..6 {
                for j in 0..6 {
                    let a = (low - temperature_k).mul_add(slope[i][j], tensor[i][j]);
                    let b = (high - temperature_k).mul_add(slope[i][j], tensor[i][j]);
                    if !a.is_finite() || !b.is_finite() || !slope[i][j].is_finite() {
                        return Err(invalid("affine expansion endpoint tensor overflows"));
                    }
                    storage += alpha[i] * alpha[j] * a.abs().max(b.abs());
                    // Each strain component is bounded by its native Mandel norm.
                    cross += alpha[i] * slope[i][j].abs();
                    expansion_curve += beta[i] * a.abs().max(b.abs());
                }
            }
            let value = storage + 2.0 * strain_bound * cross + strain_bound * expansion_curve;
            if !value.is_finite() || value < 0.0 {
                return Err(invalid("affine expansion capacity bound overflows"));
            }
            Ok(value)
        };
        let mut value = if let Some(bounds) = &derivative.equilibrium_bounds {
            if [
                bounds.stiffness_norm_pa,
                bounds.slope_norm_pa_k,
                bounds.curvature_norm_pa_k2,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
            {
                return Err(invalid("thermal tensor interval bounds are inadmissible"));
            }
            let alpha_squared = dot(&alpha, &alpha);
            let alpha_norm = det::sqrt(alpha_squared);
            0.5 * strain_limit * strain_limit * bounds.curvature_norm_pa_k2
                + 2.0 * alpha_norm * strain_limit * bounds.slope_norm_pa_k
                + alpha_squared * bounds.stiffness_norm_pa
                + strain_limit * det::sqrt(dot(&beta, &beta)) * bounds.stiffness_norm_pa
        } else {
            if derivative
                .equilibrium_pa_k2
                .iter()
                .flatten()
                .any(|v| *v != 0.0)
            {
                return Err(invalid(
                    "nonlinear tensor curvature needs complete source interval bounds",
                ));
            }
            bound(&self.model.c_eq, &derivative.equilibrium_pa_k, strain_limit)?
        };
        for (branch, slope) in self
            .model
            .branches
            .iter()
            .zip(derivative.branches_pa_k.iter())
        {
            value += bound(&branch.c_pa, slope, 2.0 * strain_limit)?;
        }
        value *= self.volume_m3;
        if expansion_curvature.is_some() {
            // Same allowance as the observable Maxwell memory admission.
            value *= THERMAL_EXPANSION_ENVELOPE_GUARD;
        }
        if !value.is_finite() || value < 0.0 {
            return Err(invalid("volume-scaled thermal capacity bound overflows"));
        }
        Ok(value)
    }
}
