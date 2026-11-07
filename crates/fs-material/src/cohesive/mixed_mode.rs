//! Mode-dependent bilinear fracture in the original cohesive point/solver path.
//!
//! Intact traction has separate Kn and Kt. Damage initiates on
//! `(Kn <gn>/Tn)^2 + (Kt |gt|/Tt)^2 = 1`. At the FIRST outward crossing of
//! that surface, beta = Kt |gt|^2/(Kn <gn>^2 + Kt |gt|^2) selects
//! the explicitly selected Benzeggagh--Kenane or power-law fracture envelope.
//! BK uses `Gc = GIc + (GIIc-GIc) beta^eta`; power law uses
//! `(GI/GIc)^alpha + (GII/GIIc)^alpha = 1`, with GI=(1-beta)Gc.
//!
//! At fixed coefficients this implementation FREEZES beta and its envelope at
//! initiation. The physical active-separation direction is retained as well:
//! an explicit `temperature` transfer reevaluates the envelope on that original
//! direction, never on the current possibly rotated separation.
//! Subsequent damage is driven by the maximum energy norm
//! `r = sqrt(Kn <gn>^2 + Kt |gt|^2)`. This makes arbitrary subsequent straight
//! jump paths passive, with exact analytic work, while reproducing the chosen
//! pure-mode and proportional mixed-mode fracture energies. It is a declared
//! onset-mix extension, not a claim that coupon data calibrate arbitrary
//! nonproportional histories. Changing direction after onset does not reselect
//! toughness or restore failed bonds. Both shear directions are isotropic.
//!
//! Quadratic initiation and BK interaction follow Camanho and Davila (2002),
//! NASA/TM-2002-211737, https://ntrs.nasa.gov/citations/20020053651.
//! The frozen-onset energy norm and its full initiation-sensitive algorithmic
//! tangent are explicit choices here. Softening tangents can be nonsymmetric;
//! the existing native Newton/FGMRES owner supports them.

use super::{
    BilinearCohesive, CohesiveEnergy, CohesiveError, CohesivePoint, CohesiveStep, active, det, dot,
    compression_work, elastic_work, finite, fracture_increment, loading_storage_change, norm,
    outward_crossing, positive,
};

mod fracture;
pub use fracture::{MixedModeCoefficients, MixedModeFractureCriterion};
pub(super) use fracture::MixedModeLaw;

/// Explicit BK ordered-interface coefficients. All values are coherent SI.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedModeParameters {
    /// Undamaged opening stiffness, Pa/m.
    pub normal_stiffness_pa_m: f64,
    /// Undamaged stiffness in each shear direction, Pa/m.
    pub tangential_stiffness_pa_m: f64,
    /// Pure opening peak traction, Pa.
    pub normal_strength_pa: f64,
    /// Pure shear peak traction (rotationally isotropic in-plane), Pa.
    pub tangential_strength_pa: f64,
    /// Pure opening fracture work, J/m2.
    pub mode_i_fracture_energy_j_m2: f64,
    /// Pure shear fracture work, J/m2.
    pub mode_ii_fracture_energy_j_m2: f64,
    /// BK exponent, at least one so initiation sensitivities remain finite.
    pub bk_exponent: f64,
}

/// Retained physical initiation mode and the envelope at current coefficients.
/// Ordinary mechanical updates never change this record. An explicit thermal
/// transfer preserves the first direction and energy fraction, while rebuilding
/// the current envelope using that direction and the new source coefficients.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedModeHistory {
    /// Original shear share of undamaged energy at first initiation, in [0,1].
    /// This historical observation does not change with temperature.
    pub energy_shear_fraction: f64,
    /// Unit active physical separation at first initiation: [opening, shear1,
    /// shear2]. Compression contributes zero; shear signs are retained.
    pub onset_direction: [f64; 3],
    /// Shear energy fraction on the ORIGINAL direction at current coefficients.
    /// Different temperature dependence of Kn and Kt changes this weighting,
    /// without changing the retained physical initiation mode.
    pub envelope_shear_fraction: f64,
    /// Current energy-norm initiation coordinate, sqrt(J/m2), not a distance.
    pub onset_sqrt_energy: f64,
    /// Complete-failure energy-norm coordinate, sqrt(J/m2).
    pub failure_sqrt_energy: f64,
    /// Selected fracture work on the retained mode at current coefficients, J/m2.
    pub fracture_energy_j_m2: f64,
}

impl BilinearCohesive {
    /// Select anisotropic stiffness/strength and mode-dependent toughness.
    /// Existing scalar getters describe the pure-opening envelope. The complete
    /// selected mixed envelope is available on each point after initiation.
    /// Compression retains the existing independent undamaged penalty.
    ///
    /// Both pure-mode failure distances must fit the separation limit. The
    /// complete BK mode interval must retain a positive softening branch at
    /// its own quadratic initiation energy. Coefficients from different modes
    /// are never combined into a fictitious worst-case adhesive. The range
    /// limit still applies if a nonproportional path rotates after initiation.
    pub fn new_mixed_mode(
        parameters: MixedModeParameters,
        compression_stiffness_pa_m: f64,
        maximum_separation_m: f64,
    ) -> Result<Self, CohesiveError> {
        Self::new_mixed_law(parameters.into(), compression_stiffness_pa_m, maximum_separation_m)
    }

    fn new_mixed_law(
        p: MixedModeLaw,
        compression_stiffness_pa_m: f64,
        maximum_separation_m: f64,
    ) -> Result<Self, CohesiveError> {
        for value in [
            p.normal_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
            p.normal_strength_pa,
            p.tangential_strength_pa,
            p.mode_i_fracture_energy_j_m2,
            p.mode_ii_fracture_energy_j_m2,
        ] {
            positive(value)?;
        }
        p.admit_criterion()?;
        let mut law = Self::new(
            p.normal_stiffness_pa_m,
            p.normal_strength_pa,
            p.mode_i_fracture_energy_j_m2,
            compression_stiffness_pa_m,
            maximum_separation_m,
        )?;
        let shear_onset = positive(p.tangential_strength_pa / p.tangential_stiffness_pa_m)?;
        let shear_failure =
            positive(2.0 * (p.mode_ii_fracture_energy_j_m2 / p.tangential_strength_pa))?;
        if shear_failure <= shear_onset || shear_failure > maximum_separation_m {
            return Err(CohesiveError::Invalid(
                "shear onset < failure <= separation limit required",
            ));
        }
        let a = positive(p.normal_strength_pa / det::sqrt(p.normal_stiffness_pa_m))?;
        let b = positive(p.tangential_strength_pa / det::sqrt(p.tangential_stiffness_pa_m))?;
        admit_softening_interval(p, positive((0.5 * a) * a)?, positive((0.5 * b) * b)?)?;
        p.admit_power_range(a.min(b), maximum_separation_m)?;
        positive((p.tangential_stiffness_pa_m * maximum_separation_m) * maximum_separation_m)?;
        law.mixed_mode = Some(p);
        Ok(law)
    }

    /// Original BK inputs, absent for common-mode and power-law models.
    /// Use `mixed_mode_coefficients` and `mixed_mode_criterion` for either family.
    #[must_use]
    pub const fn mixed_mode_parameters(&self) -> Option<MixedModeParameters> {
        match self.mixed_mode {
            Some(MixedModeLaw { coefficients: p,
                criterion: MixedModeFractureCriterion::BenzeggaghKenane { exponent } }) =>
                Some(MixedModeParameters {
                    normal_stiffness_pa_m: p.normal_stiffness_pa_m,
                    tangential_stiffness_pa_m: p.tangential_stiffness_pa_m,
                    normal_strength_pa: p.normal_strength_pa,
                    tangential_strength_pa: p.tangential_strength_pa,
                    mode_i_fracture_energy_j_m2: p.mode_i_fracture_energy_j_m2,
                    mode_ii_fracture_energy_j_m2: p.mode_ii_fracture_energy_j_m2,
                    bk_exponent: exponent,
                }),
            _ => None,
        }
    }
}

// At shear-energy fraction beta the quadratic strength criterion gives
// E0(beta) = 1 / ((1-beta)/EI0 + beta/EII0), where EI0 = Tn^2/(2Kn).
// A positive bilinear branch requires Gc(beta)/E0(beta) > 1 throughout [0,1].
// Normalize both positive factors before multiplication to retain SI rescaling
// and avoid large reciprocal energies. Ordinary mechanical steps do no work here.
fn admit_softening_interval(
    p: MixedModeLaw,
    normal_onset_energy: f64,
    shear_onset_energy: f64,
) -> Result<(), CohesiveError> {
    let bk_exponent = match p.criterion {
        MixedModeFractureCriterion::BenzeggaghKenane { exponent } => exponent,
        MixedModeFractureCriterion::PowerLaw { exponent } =>
            return p.admit_power_softening(normal_onset_energy, shear_onset_energy, exponent),
    };
    // Preserve the inexpensive sufficient screen for ordinary well-separated
    // envelopes, including repeated source-bound thermal coefficient refreshes.
    let relative_guard = 1.0 + 64.0 * f64::EPSILON;
    if p.mode_i_fracture_energy_j_m2
        .min(p.mode_ii_fracture_energy_j_m2)
        / normal_onset_energy.max(shear_onset_energy)
        > relative_guard
    {
        return Ok(());
    }
    let minimum_onset = normal_onset_energy.min(shear_onset_energy);
    let maximum_toughness = p
        .mode_i_fracture_energy_j_m2
        .max(p.mode_ii_fracture_energy_j_m2);
    let a = positive(minimum_onset / normal_onset_energy)?;
    let b = positive(minimum_onset / shear_onset_energy)?;
    let gi = positive(p.mode_i_fracture_energy_j_m2 / maximum_toughness)?;
    let gii = positive(p.mode_ii_fracture_energy_j_m2 / maximum_toughness)?;
    let threshold = positive(minimum_onset / maximum_toughness)?;
    let mut minimum = (gi * a).min(gii * b);
    // Except when toughness AND initiation energy increase towards shear with
    // eta>1, the product is monotone or concave, so the endpoints suffice.
    // In the remaining case its second derivative changes sign at most once;
    // the convex prefix contains the only possible interior minimum.
    if gii > gi && a > b && bk_exponent > 1.0 {
        let inv = 1.0 / (bk_exponent + 1.0);
        let convex_end = (((bk_exponent - 1.0) * inv) * (a / (a - b))).min(1.0);
        let derivative = |beta: f64| {
            let power = det::pow(beta, bk_exponent);
            let toughness = (1.0 - power) * gi + power * gii;
            let compliance = (1.0 - beta) * a + beta * b;
            (gii - gi) * (1.0 - inv) * det::pow(beta, bk_exponent - 1.0) * compliance
                - (a - b) * inv * toughness
        };
        if derivative(convex_end) > 0.0 {
            let mut low = 0.0;
            let mut high = convex_end;
            // Locate the sole minimum on the convex interval. This is a
            // construction-time bounded scalar solve, independent of meshes.
            for _ in 0..64 {
                let mid = low + 0.5 * (high - low);
                if derivative(mid) > 0.0 {
                    high = mid;
                } else {
                    low = mid;
                }
            }
            // Gc increases and the reciprocal onset energy decreases over
            // this bracket. Their separate lower endpoints bound the entire
            // unresolved interval, even if eta concentrates the minimum closer
            // to beta=1 than one floating-point spacing can represent.
            let power = det::pow(low, bk_exponent);
            let lower_toughness = (1.0 - power) * gi + power * gii;
            let lower_compliance = (1.0 - high) * a + high * b;
            minimum = minimum.min(lower_toughness * lower_compliance);
        }
    }
    // An unresolved near-zero softening interval is refused rather than being
    // admitted on a last-bit comparison of two independently rounded products.
    if !minimum.is_finite() || minimum <= threshold * relative_guard {
        return Err(CohesiveError::Invalid(
            "mixed-mode envelope has no resolved positive softening interval",
        ));
    }
    Ok(())
}

impl MixedModeHistory {
    pub(super) fn intact(self, maximum: f64) -> f64 {
        if maximum <= self.onset_sqrt_energy {
            1.0
        } else if maximum >= self.failure_sqrt_energy {
            0.0
        } else {
            (self.onset_sqrt_energy / maximum)
                * ((self.failure_sqrt_energy - maximum)
                    / (self.failure_sqrt_energy - self.onset_sqrt_energy))
        }
    }
    pub(super) fn fracture(self, maximum: f64) -> f64 {
        if maximum <= self.onset_sqrt_energy {
            0.0
        } else if maximum >= self.failure_sqrt_energy {
            self.fracture_energy_j_m2
        } else {
            self.fracture_energy_j_m2
                * ((maximum - self.onset_sqrt_energy)
                    / (self.failure_sqrt_energy - self.onset_sqrt_energy))
        }
    }
    /// Intersect the new quadratic strength surface with the SAME physical
    /// initiation ray. Its energy weighting changes with Kn/Kt, so temperature
    /// dependent fracture data use that weighting rather than a rotated current gap.
    pub(super) fn at_parameters(self, p: MixedModeLaw) -> Result<Self, CohesiveError> {
        let radius = positive(coordinate(p, self.onset_direction)?)?;
        let onset = positive(radius / positive(strength_coordinate(p, self.onset_direction)?)?)?;
        let shear = det::sqrt(p.tangential_stiffness_pa_m)
            * norm([0.0, self.onset_direction[1], self.onset_direction[2]])?;
        let beta = finite((shear / radius) * (shear / radius))?.clamp(0.0, 1.0);
        let fracture = p.toughness(self.onset_direction, radius, beta)?;
        let failure = positive(2.0 * (fracture / onset))?;
        if failure <= onset {
            return Err(CohesiveError::Invalid(
                "thermal mixed envelope leaves positive-softening domain",
            ));
        }
        Ok(Self {
            envelope_shear_fraction: beta,
            onset_sqrt_energy: onset,
            failure_sqrt_energy: failure,
            fracture_energy_j_m2: fracture,
            ..self
        })
    }
    fn envelope_work(self, start: f64, end: f64) -> f64 {
        let a = start.min(self.onset_sqrt_energy);
        let b = end.min(self.onset_sqrt_energy);
        let elastic = 0.5 * (b - a) * (b + a);
        let a = start
            .max(self.onset_sqrt_energy)
            .min(self.failure_sqrt_energy);
        let b = end
            .max(self.onset_sqrt_energy)
            .min(self.failure_sqrt_energy);
        elastic
            + self.onset_sqrt_energy
                * ((b - a) / (self.failure_sqrt_energy - self.onset_sqrt_energy))
                * (0.5 * (self.failure_sqrt_energy - a) + 0.5 * (self.failure_sqrt_energy - b))
    }
}

impl<'a> CohesivePoint<'a> {
    /// Envelope selected at initiation, absent before initiation and for the
    /// original mode-independent law. Returned observations cannot alter history.
    #[must_use]
    pub const fn mixed_mode_history(&self) -> Option<MixedModeHistory> {
        self.mixed_history
    }

    pub(super) fn mixed_integrity(&self) -> f64 {
        self.mixed_history.map_or(1.0, |h| h.intact(self.maximum))
    }
    pub(super) fn mixed_stored_energy(&self) -> f64 {
        let p = self.law.mixed_mode.expect("mixed-mode owner");
        let g = active(self.jump);
        let recoverable = 0.5 * (p.normal_stiffness_pa_m * g[0]) * g[0]
            + 0.5 * (p.tangential_stiffness_pa_m * g[1]) * g[1]
            + 0.5 * (p.tangential_stiffness_pa_m * g[2]) * g[2];
        let compression = self.jump[0].min(0.0);
        self.mixed_integrity() * recoverable
            + (0.5 * self.law.compression_stiffness * compression) * compression
    }
    pub(super) fn mixed_traction(&self) -> [f64; 3] {
        let p = self.law.mixed_mode.expect("mixed-mode owner");
        let intact = self.mixed_integrity();
        [
            if self.jump[0] < 0.0 {
                self.law.compression_stiffness * self.jump[0]
            } else {
                intact * p.normal_stiffness_pa_m * self.jump[0]
            },
            intact * p.tangential_stiffness_pa_m * self.jump[1],
            intact * p.tangential_stiffness_pa_m * self.jump[2],
        ]
    }

    pub(super) fn trial_mixed(
        &self,
        jump: [f64; 3],
        tolerance: f64,
    ) -> Result<CohesiveStep<'a>, CohesiveError> {
        let p = self.law.mixed_mode.expect("mixed-mode owner");
        let radius = coordinate(p, jump)?;
        let old_radius = coordinate(p, self.jump)?;
        let initiation = if self.mixed_history.is_none() && strength_coordinate(p, jump)? >= 1.0 {
            Some(initiate(p, self.jump, jump)?)
        } else {
            None
        };
        let history = self.mixed_history.or(initiation.map(|x| x.0));
        let old_maximum = self.mixed_history.map_or_else(
            || history.map_or(0.0, |h| h.onset_sqrt_energy),
            |_| self.maximum,
        );
        let maximum = history.map_or(0.0, |h| old_maximum.max(radius).min(h.failure_sqrt_energy));
        let point = Self {
            law: self.law,
            jump,
            maximum,
            mixed_history: history,
            ..*self
        };
        let old_intact = self.mixed_integrity();
        let stiffness = [
            p.normal_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
        ];
        let (work, stored_change) = match history {
            Some(h) if radius > old_maximum && radius > h.onset_sqrt_energy => {
                let retained =
                    0.5 * old_intact * (old_maximum - old_radius) * (old_maximum + old_radius);
                (retained + h.envelope_work(old_maximum, radius),
                    retained + loading_storage_change(1.0, h.onset_sqrt_energy,
                        h.failure_sqrt_energy, old_maximum, radius))
            }
            _ => {
                let elastic = elastic_work(stiffness.map(|k| old_intact*k), self.jump, jump)?;
                (elastic, elastic)
            }
        };
        let compression = compression_work(self.law.compression_stiffness, self.jump[0], jump[0])?;
        let work = finite(work + compression)?;
        let stored_change = finite(stored_change + compression)?;
        let dissipated = finite(history.map_or(0.0, |h|
            fracture_increment(h.fracture_energy_j_m2, h.onset_sqrt_energy,
                h.failure_sqrt_energy, old_maximum, maximum)))?;
        let defect = finite(work - stored_change - dissipated)?;
        if defect.abs() > tolerance {
            return Err(CohesiveError::EnergyDefect(defect));
        }
        let traction_pa = point.traction_pa();
        for value in traction_pa {
            finite(value)?;
        }
        let active_jump = active(jump);
        let intact = point.mixed_integrity();
        let mut tangent = [[0.0; 3]; 3];
        for i in 0..3 {
            tangent[i][i] = intact * stiffness[i];
        }
        if jump[0] < 0.0 {
            tangent[0][0] = self.law.compression_stiffness;
        }
        let mut classical = jump[0] != 0.0;
        if let Some(h) = history {
            classical &= !(radius == h.onset_sqrt_energy
                || radius == h.failure_sqrt_energy
                || (radius == old_maximum && radius > h.onset_sqrt_energy && !point.failed()));
            if radius > old_maximum
                && radius > h.onset_sqrt_energy
                && radius < h.failure_sqrt_energy
            {
                let a = h.onset_sqrt_energy;
                let b = 2.0 * h.fracture_energy_j_m2;
                let denominator = b - a * a;
                let derivative_r = -(a / radius) * (b / denominator) / radius;
                let mut derivative = core::array::from_fn::<_, 3, _>(|i| {
                    derivative_r * stiffness[i] * active_jump[i] / radius
                });
                if let Some((_, lambda)) = initiation {
                    classical &= p.smooth_initiation(h.onset_direction);
                    let (da, dgc) = initiation_gradient(p, self.jump, jump, lambda, h)?;
                    let derivative_a =
                        (b / denominator) * ((b + a * a) / radius - 2.0 * a) / denominator;
                    let derivative_b = (a * a / denominator) * (1.0 - a / radius) / denominator;
                    for i in 0..3 {
                        derivative[i] += derivative_a * da[i] + 2.0 * derivative_b * dgc[i];
                    }
                }
                for i in 0..3 {
                    for j in 0..3 {
                        tangent[i][j] += stiffness[i] * active_jump[i] * derivative[j];
                    }
                }
            }
        }
        for row in tangent {
            for value in row {
                finite(value)?;
            }
        }
        Ok(CohesiveStep {
            point,
            traction_pa,
            tangent_pa_m: classical.then_some(tangent),
            branch_tangent_pa_m: tangent,
            energy: CohesiveEnergy {
                work,
                stored_change,
                dissipated,
                defect,
            },
        })
    }
}

pub(super) fn coordinate(p: MixedModeLaw, g: [f64; 3]) -> Result<f64, CohesiveError> {
    norm([
        det::sqrt(p.normal_stiffness_pa_m) * g[0].max(0.0),
        det::sqrt(p.tangential_stiffness_pa_m) * g[1],
        det::sqrt(p.tangential_stiffness_pa_m) * g[2],
    ])
}
pub(super) fn strength_coordinate(
    p: MixedModeLaw,
    g: [f64; 3],
) -> Result<f64, CohesiveError> {
    norm([
        g[0].max(0.0) / (p.normal_strength_pa / p.normal_stiffness_pa_m),
        g[1] / (p.tangential_strength_pa / p.tangential_stiffness_pa_m),
        g[2] / (p.tangential_strength_pa / p.tangential_stiffness_pa_m),
    ])
}
pub(super) fn initiate(
    p: MixedModeLaw,
    old: [f64; 3],
    end: [f64; 3],
) -> Result<(MixedModeHistory, f64), CohesiveError> {
    let scales = [
        p.normal_strength_pa / p.normal_stiffness_pa_m,
        p.tangential_strength_pa / p.tangential_stiffness_pa_m,
        p.tangential_strength_pa / p.tangential_stiffness_pa_m,
    ];
    let lambda = outward_crossing(
        core::array::from_fn(|i| old[i] / scales[i]),
        core::array::from_fn(|i| end[i] / scales[i]),
        1.0,
        0.0,
    )?;
    let g = core::array::from_fn(|i| (1.0 - lambda) * old[i] + lambda * end[i]);
    let onset = positive(coordinate(p, g)?)?;
    let shear = det::sqrt(p.tangential_stiffness_pa_m) * norm([0.0, g[1], g[2]])?;
    let beta = finite((shear / onset) * (shear / onset))?.clamp(0.0, 1.0);
    let fracture = p.toughness(g, onset, beta)?;
    let failure = positive(2.0 * (fracture / onset))?;
    if failure <= onset {
        return Err(CohesiveError::Invalid(
            "mixed initiation leaves positive-softening domain",
        ));
    }
    let active_g = active(g);
    let length = positive(norm(active_g)?)?;
    Ok((
        MixedModeHistory {
            energy_shear_fraction: beta,
            onset_direction: active_g.map(|component| component / length),
            envelope_shear_fraction: beta,
            onset_sqrt_energy: onset,
            failure_sqrt_energy: failure,
            fracture_energy_j_m2: fracture,
        },
        lambda,
    ))
}
pub(super) fn initiation_gradient(
    p: MixedModeLaw,
    old: [f64; 3],
    end: [f64; 3],
    lambda: f64,
    h: MixedModeHistory,
) -> Result<([f64; 3], [f64; 3]), CohesiveError> {
    let delta = core::array::from_fn::<_, 3, _>(|i| end[i] - old[i]);
    let g = active(core::array::from_fn(|i| old[i] + lambda * delta[i]));
    let onset = [
        p.normal_strength_pa / p.normal_stiffness_pa_m,
        p.tangential_strength_pa / p.tangential_stiffness_pa_m,
        p.tangential_strength_pa / p.tangential_stiffness_pa_m,
    ];
    let coefficient = core::array::from_fn::<_, 3, _>(|i| {
        if i == 0 && g[0] == 0.0 {
            0.0
        } else {
            2.0 * (g[i] / onset[i]) / onset[i]
        }
    });
    let denominator = positive(dot(coefficient, delta))?;
    let stiffness = [
        p.normal_stiffness_pa_m,
        p.tangential_stiffness_pa_m,
        p.tangential_stiffness_pa_m,
    ];
    let radius = h.onset_sqrt_energy;
    let beta = h.energy_shear_fraction;
    let gradient_radius = core::array::from_fn::<_, 3, _>(|i| stiffness[i] * g[i] / radius);
    let gradient_gc = p.toughness_gradient(g, radius, beta, h.fracture_energy_j_m2)?;
    let mut da = [0.0; 3];
    let mut dgc = [0.0; 3];
    for j in 0..3 {
        let dlambda = -lambda * coefficient[j] / denominator;
        for i in 0..3 {
            let dg = if i == j { lambda } else { 0.0 } + delta[i] * dlambda;
            da[j] += gradient_radius[i] * dg;
            dgc[j] += gradient_gc[i] * dg;
        }
        finite(da[j])?;
        finite(dgc[j])?;
    }
    Ok((da, dgc))
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod power_law_tests;
