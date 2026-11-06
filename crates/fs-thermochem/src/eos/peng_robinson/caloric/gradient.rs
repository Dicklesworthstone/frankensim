//! Original PR derivatives and sufficient homogeneous rectangle admission.
//!
//! These scalar operations retain the same prepared source and fixed mixture.
//! They do not choose a phase or replace a vapor-liquid stability/flash solve.

use super::{
    PengRobinsonCaloricReferenceV1, PengRobinsonDensityCaloricEvaluationV1,
    PengRobinsonErrorV1, R, Sum, Temperature, det, unrepresentable,
};

// The unique minimum on x>1 of
// f(x)=(x*x+2*x-1)^2 / (2*(x+1)*(x-1)^2) occurs at
// x^3-3*x^2-3*x-3=0. This is the same PR spinodal geometry used in saturation.
const CRITICAL_VOLUME_OVER_B: f64 = 3.951_373_035_591_441;

impl PengRobinsonDensityCaloricEvaluationV1 {
    /// Original mixture covolume b [m3/mol].
    #[must_use]
    pub const fn covolume_m3_per_mol(&self) -> f64 {
        self.caloric.parameters.b_mix
    }

    /// (d2p/dT2)_v [Pa/K2] and d2p/(dT dv) [Pa mol/(K m3)].
    /// These are analytic derivatives of this SAME density-based source state;
    /// the volume coordinate is molar volume, not mass-specific or total volume.
    #[must_use]
    pub fn pressure_temperature_second_derivatives(&self) -> (f64, f64) {
        let v = self.molar_volume_m3_per_mol();
        let p = self.caloric.parameters;
        let b = p.b_mix;
        let d = v * v + 2.0 * b * v - b * b;
        (
            -p.d2a_mix_dt2 / d,
            -R / ((v - b) * (v - b))
                + p.da_mix_dt * (2.0 * v + 2.0 * b) / (d * d),
        )
    }

    /// Molar (du/dv)_T = T*(dp/dT)_v-p [Pa]. Ideal-reference u(T) adds zero.
    /// Evaluate the departure directly to avoid cancelling the ideal pressure,
    /// including in the dilute limit where the residual is much smaller.
    #[must_use]
    pub fn internal_energy_volume_derivative_pa(&self) -> f64 {
        let v = self.molar_volume_m3_per_mol();
        let p = self.caloric.parameters;
        let b = p.b_mix;
        let d = v * v + 2.0 * b * v - b * b;
        (p.a_mix - self.temperature().value() * p.da_mix_dt) / d
    }
}

impl PengRobinsonCaloricReferenceV1<'_> {
    /// Admit a sufficient continuous mechanically stable, positive-pressure
    /// rectangle between the retained densities and the supplied T bounds.
    /// Return a conservative lower bound for residual Cv [J/(mol K)] there;
    /// the material owner must add and admit its ORIGINAL ideal-reference Cv.
    ///
    /// Bounds may be reversed or equal but must include the reference T.
    /// Source alpha zeros/sign changes are refused, because the branch-local
    /// temperature derivatives are not smooth through those points.
    ///
    /// Each component's |C/sqrt(T)-D| is bounded by its two endpoint values.
    /// Nonnegative original mixing weights then bound a(T)/T. The exact PR
    /// spinodal threshold has one interior minimum in v/b, so checking that
    /// minimum admits the whole volume interval rather than only its corners.
    /// Positive pressure at the largest volume, with this same attraction
    /// bound, then implies positive pressure throughout the rectangle.
    ///
    /// This sufficient bound can refuse a locally admissible rectangle: narrow
    /// the temperature/volume range. Successful admission does not establish
    /// global phase equilibrium; metastable homogeneous states remain possible.
    /// Work is allocation-free and bounded by three original-size pair loops
    /// and four scalar corners, with no iterative phase solve or sampling grid.
    pub fn admit_stable_rectangle(
        &self,
        temperature_bounds: [Temperature; 2],
    ) -> Result<f64, PengRobinsonErrorV1> {
        let [ta, tb] = temperature_bounds.map(|t| t.value());
        for t in [ta, tb] {
            if !t.is_finite() || t <= 0.0 {
                return Err(PengRobinsonErrorV1::InvalidTemperature { bits: t.to_bits() });
            }
        }
        let lo = ta.min(tb);
        let hi = ta.max(tb);
        let initial_t = self.initial.temperature().value();
        if initial_t < lo || initial_t > hi {
            return Err(unrepresentable("rectangle omits the original reference temperature"));
        }
        let root_lo = det::sqrt(lo);
        let root_hi = det::sqrt(hi);
        let n = self.mixture.components.len();
        let mut maxima = [0.0; super::MAX_PENG_ROBINSON_COMPONENTS_V1];
        for (i, maximum) in maxima[..n].iter_mut().enumerate() {
            let c = self.components[i];
            let at_lo = c.constant / root_lo - c.root_coefficient;
            let at_hi = c.constant / root_hi - c.root_coefficient;
            // Absolute padding matters near a component alpha zero; a relative
            // bound on the nearly cancelled result would not protect its sign.
            let scale = (c.constant / root_lo).abs()
                + (c.constant / root_hi).abs() + 2.0 * c.root_coefficient.abs();
            let padding = 128.0 * f64::EPSILON * scale;
            if ![at_lo, at_hi, padding].iter().all(|x| x.is_finite())
                || at_lo.abs() <= padding || at_hi.abs() <= padding
                || at_lo.is_sign_negative() != at_hi.is_sign_negative()
            {
                return Err(unrepresentable("rectangle crosses or approaches a component alpha zero"));
            }
            *maximum = at_lo.abs().max(at_hi.abs()) + padding;
        }
        let mut attraction = Sum::default();
        let mut half_b_magnitude = Sum::default();
        for i in 0..n {
            for j in 0..n {
                let weight = self.mixture.mole_fractions[i]
                    * self.mixture.mole_fractions[j]
                    * (1.0 - self.mixture.binary_interactions[i * n + j]);
                attraction.add(weight * maxima[i] * maxima[j]);
                let left = self.components[i];
                let right = self.components[j];
                half_b_magnitude.add(weight * (
                    (0.5 * left.constant * right.root_coefficient).abs()
                        + (0.5 * right.constant * left.root_coefficient).abs()
                ));
            }
        }
        let b = self.initial.covolume_m3_per_mol();
        let v0 = self.initial.molar_volume_m3_per_mol();
        let v1 = self.mixture.molar_mass.value() / self.final_density.value();
        let eta_min = b / v0.max(v1);
        let eta_max = b / v0.min(v1);
        let roundoff = (128.0 + 4.0 * (n * n) as f64) * f64::EPSILON;
        let alpha = attraction.value() / (R * b) * (1.0 + roundoff);
        if !alpha.is_finite() || alpha <= 0.0
            || !eta_min.is_finite() || eta_min <= 0.0
            || !eta_max.is_finite() || eta_max >= 1.0
        {
            return Err(unrepresentable("rectangle dimensionless attraction or covolume"));
        }
        // eta=b/v avoids squaring an enormous dilute v/b ratio.
        let eta = (1.0 / CRITICAL_VOLUME_OVER_B).clamp(eta_min, eta_max);
        let d = 1.0 + 2.0 * eta - eta * eta;
        let free = 1.0 - eta;
        let destabilizing = alpha * (2.0 * eta * (1.0 + eta) * free * free);
        if !destabilizing.is_finite() || destabilizing >= d * d * (1.0 - roundoff) {
            return Err(PengRobinsonErrorV1::MechanicallyUnstable);
        }
        let d_end = 1.0 + 2.0 * eta_min - eta_min * eta_min;
        let attractive_pressure = alpha * eta_min * (1.0 - eta_min);
        if !attractive_pressure.is_finite()
            || attractive_pressure >= d_end * (1.0 - roundoff)
        {
            return Err(unrepresentable("rectangle cannot guarantee positive pressure"));
        }

        // With no alpha sign change, a''(T)=B/(4*T^(3/2)) on this exact
        // branch, hence Cv^R=B*I(v)/(4*sqrt(T)). I(v) is positive and
        // decreasing. Its extrema are at the four original corners even if
        // B is negative. Two source mixing loops serve both densities.
        let mut minimum = f64::INFINITY;
        let mut cv_scale = 0.0_f64;
        for t in [lo, hi] {
            let mixed = self.mixture.mixing_caloric(t, |_| {})?;
            for density in [self.initial.density(), self.final_density] {
                let corner = mixed.density_caloric_evaluation(
                    Temperature::new(t), density, self.mixture.molar_mass,
                )?;
                let cv = corner.residual_cv_j_mol_k();
                let magnitude = (half_b_magnitude.value() / (2.0 * det::sqrt(t)))
                    * corner.caloric.residual_geometry.attraction_integral;
                if !cv.is_finite() || !magnitude.is_finite() {
                    return Err(unrepresentable("rectangle residual heat-capacity bound"));
                }
                minimum = minimum.min(cv);
                cv_scale = cv_scale.max(cv.abs()).max(magnitude);
            }
        }
        let bound = minimum - roundoff * cv_scale;
        if !bound.is_finite() {
            return Err(unrepresentable("rectangle residual heat-capacity bound"));
        }
        Ok(bound)
    }
}

#[cfg(test)]
mod tests;
