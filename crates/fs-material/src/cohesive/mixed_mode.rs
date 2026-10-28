//! Mode-dependent bilinear fracture in the original cohesive point/solver path.
//!
//! Intact traction has separate Kn and Kt. Damage initiates on
//! `(Kn <gn>/Tn)^2 + (Kt |gt|/Tt)^2 = 1`. At the FIRST outward crossing of
//! that surface, beta = Kt |gt|^2/(Kn <gn>^2 + Kt |gt|^2) selects
//! `Gc = GIc + (GIIc-GIc) beta^eta` (Benzeggagh--Kenane).
//!
//! This implementation explicitly FREEZES beta and its envelope at initiation.
//! Subsequent damage is driven by the maximum energy norm
//! `r = sqrt(Kn <gn>^2 + Kt |gt|^2)`. This makes arbitrary subsequent straight
//! jump paths passive, with exact analytic work, while reproducing the chosen
//! pure-mode and proportional mixed-mode fracture energies. It is a declared
//! onset-mix extension, not a claim that BK coupon data calibrate arbitrary
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
    finite, norm, outward_crossing, positive,
};

/// Explicit ordered-interface coefficients. All values are coherent SI.
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

/// Immutable envelope chosen at the first strength-surface crossing.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixedModeHistory {
    /// Shear share of undamaged recoverable energy at initiation, in [0,1].
    pub energy_shear_fraction: f64,
    /// Energy-norm initiation coordinate, sqrt(J/m2), not a distance.
    pub onset_sqrt_energy: f64,
    /// Complete-failure energy-norm coordinate, sqrt(J/m2).
    pub failure_sqrt_energy: f64,
    /// Frozen BK fracture work for this point, J/m2.
    pub fracture_energy_j_m2: f64,
}

impl BilinearCohesive {
    /// Select anisotropic stiffness/strength and mode-dependent toughness.
    /// Existing scalar getters describe the pure-opening envelope. The complete
    /// selected mixed envelope is available on each point after initiation.
    /// Compression retains the existing independent undamaged penalty.
    ///
    /// Both pure-mode failure distances must fit the separation limit. A
    /// conservative all-mode screen requires twice the smaller toughness to
    /// exceed the larger pure-mode initiation energy norm squared. The range
    /// limit still applies if a nonproportional path rotates after initiation.
    pub fn new_mixed_mode(
        parameters: MixedModeParameters,
        compression_stiffness_pa_m: f64,
        maximum_separation_m: f64,
    ) -> Result<Self, CohesiveError> {
        let p = parameters;
        for value in [
            p.normal_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
            p.normal_strength_pa,
            p.tangential_strength_pa,
            p.mode_i_fracture_energy_j_m2,
            p.mode_ii_fracture_energy_j_m2,
            p.bk_exponent,
        ] {
            positive(value)?;
        }
        if p.bk_exponent < 1.0 {
            return Err(CohesiveError::Invalid("BK exponent must be at least one"));
        }
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
        let maximum_onset_energy = positive(a.max(b) * a.max(b))?;
        if 0.5 * maximum_onset_energy
            >= p.mode_i_fracture_energy_j_m2
                .min(p.mode_ii_fracture_energy_j_m2)
        {
            return Err(CohesiveError::Invalid(
                "mixed-mode envelope has no positive softening interval",
            ));
        }
        positive((p.tangential_stiffness_pa_m * maximum_separation_m) * maximum_separation_m)?;
        law.mixed_mode = Some(p);
        Ok(law)
    }

    /// Explicit mode-dependent data, absent for the original isotropic law.
    #[must_use]
    pub const fn mixed_mode_parameters(&self) -> Option<MixedModeParameters> {
        self.mixed_mode
    }
}

impl MixedModeHistory {
    fn intact(self, maximum: f64) -> f64 {
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
        };
        let old_intact = self.mixed_integrity();
        let work = match history {
            Some(h) if radius > old_maximum => {
                0.5 * old_intact * (old_maximum - old_radius) * (old_maximum + old_radius)
                    + h.envelope_work(old_maximum, radius)
            }
            _ => 0.5 * old_intact * (radius - old_radius) * (radius + old_radius),
        };
        let a = self.jump[0].min(0.0);
        let b = jump[0].min(0.0);
        let work = finite(work + 0.5 * self.law.compression_stiffness * (b - a) * (b + a))?;
        let stored_change = finite(point.stored_energy_j_m2() - self.stored_energy_j_m2())?;
        let dissipated = finite(point.dissipated_j_m2() - self.dissipated_j_m2())?;
        let defect = finite(work - stored_change - dissipated)?;
        if defect.abs() > tolerance {
            return Err(CohesiveError::EnergyDefect(defect));
        }
        let traction_pa = point.traction_pa();
        for value in traction_pa {
            finite(value)?;
        }
        let stiffness = [
            p.normal_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
        ];
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

fn coordinate(p: MixedModeParameters, g: [f64; 3]) -> Result<f64, CohesiveError> {
    norm([
        det::sqrt(p.normal_stiffness_pa_m) * g[0].max(0.0),
        det::sqrt(p.tangential_stiffness_pa_m) * g[1],
        det::sqrt(p.tangential_stiffness_pa_m) * g[2],
    ])
}
fn strength_coordinate(p: MixedModeParameters, g: [f64; 3]) -> Result<f64, CohesiveError> {
    norm([
        g[0].max(0.0) / (p.normal_strength_pa / p.normal_stiffness_pa_m),
        g[1] / (p.tangential_strength_pa / p.tangential_stiffness_pa_m),
        g[2] / (p.tangential_strength_pa / p.tangential_stiffness_pa_m),
    ])
}
pub(super) fn initiate(
    p: MixedModeParameters,
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
    let fracture = positive(
        p.mode_i_fracture_energy_j_m2
            + (p.mode_ii_fracture_energy_j_m2 - p.mode_i_fracture_energy_j_m2)
                * det::pow(beta, p.bk_exponent),
    )?;
    let failure = positive(2.0 * (fracture / onset))?;
    if failure <= onset {
        return Err(CohesiveError::Invalid(
            "mixed initiation leaves positive-softening domain",
        ));
    }
    Ok((
        MixedModeHistory {
            energy_shear_fraction: beta,
            onset_sqrt_energy: onset,
            failure_sqrt_energy: failure,
            fracture_energy_j_m2: fracture,
        },
        lambda,
    ))
}
fn initiation_gradient(
    p: MixedModeParameters,
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
    let gradient_beta = [
        -2.0 * beta * stiffness[0] * g[0] / (radius * radius),
        2.0 * (1.0 - beta) * stiffness[1] * g[1] / (radius * radius),
        2.0 * (1.0 - beta) * stiffness[2] * g[2] / (radius * radius),
    ];
    let slope = (p.mode_ii_fracture_energy_j_m2 - p.mode_i_fracture_energy_j_m2)
        * p.bk_exponent
        * if beta == 0.0 {
            if p.bk_exponent == 1.0 { 1.0 } else { 0.0 }
        } else {
            det::pow(beta, p.bk_exponent - 1.0)
        };
    let mut da = [0.0; 3];
    let mut dgc = [0.0; 3];
    for j in 0..3 {
        let dlambda = -lambda * coefficient[j] / denominator;
        for i in 0..3 {
            let dg = if i == j { lambda } else { 0.0 } + delta[i] * dlambda;
            da[j] += gradient_radius[i] * dg;
            dgc[j] += slope * gradient_beta[i] * dg;
        }
        finite(da[j])?;
        finite(dgc[j])?;
    }
    Ok((da, dgc))
}

#[cfg(test)]
mod tests;
