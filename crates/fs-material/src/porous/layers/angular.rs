//! Angular power on the existing planar, source-owned layer responses.
//!
//! Different directions are mutually incoherent; their POWERS add, not their
//! pressure reflection coefficients. Isotropic diffuse incidence supplies
//! `dP proportional to cos(theta) sin(theta) dtheta`. With `mu=cos(theta)`,
//! the normalized hemisphere average is `2 integral_0^1 mu A(mu) dmu`.
//! A declared cone `mu >= mu_min` divides by `(1-mu_min^2)/2` instead.
//! See Sgard, Atalla and Robin (2024), Eq. (9),
//! <https://doi.org/10.3389/facou.2024.1414356>.
//!
//! This averages the original extended-reaction solutions, including their
//! internal refraction. It is not a locally reacting impedance fit or a
//! reverberation-room/finite-panel prediction. Polar resolution is explicit;
//! refinement is required near resonances and critical refraction. No
//! quadrature error bound or experimentally validated diffuse field is implied.

use super::{RigidBackedLayers, TransmittingLayers};
use crate::MaterialError;
use fs_math::dd::Dd;

/// Hard ceiling on directional responses in one angular-power query.
pub const MAX_ANGULAR_POWER_SAMPLES: usize = 65_536;

/// Caller-owned work allowance, including zero-power declared directions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AngularPowerBudget {
    /// At most this many directional responses; must lie in `1..=65536`.
    pub max_samples: usize,
}

/// An incident direction and its power crossing the specimen's physical face.
/// Power already contains projected area, solid angle and any quadrature
/// weight. No further cosine, area or pressure-amplitude factor is applied.
/// Isotropic layers need only the polar cosine; different azimuths may be
/// combined by adding their incident powers. Directions are incoherent.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncidentPower {
    /// Positive incidence cosine from the exposed-face normal, in `(0,1]`.
    pub cosine: f64,
    /// Finite nonnegative incident normal power [W].
    pub power_w: f64,
}

/// Isotropic intensity over the stated incident cone at one frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiffuseIncidence {
    /// `0` includes the complete hemisphere; values in `(0,1)` explicitly
    /// exclude more grazing incidence. `1` has zero solid-angle measure and
    /// refuses; use a single `IncidentPower` for normal incidence instead.
    pub minimum_cosine: f64,
    /// Composite midpoint panels in cosine, with positive exact panel weights
    /// for the incident measure `mu dmu`. This is a convergence axis, not a
    /// tolerance or certificate. No exact grazing singularity is sampled.
    pub polar_samples: usize,
}

/// Separately integrated powers and their incident-power-normalized fractions.
/// Nothing is clipped or rescaled to force `R+T+A=1`; the balance residual
/// retains the original solvers' floating-point error. No complex pressure
/// coefficient exists for this incoherent angular mixture.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularPower {
    pub incident_power_w: f64,
    pub reflected_power_w: f64,
    pub transmitted_power_w: f64,
    pub absorbed_power_w: f64,
    pub reflected_fraction: f64,
    pub transmitted_fraction: f64,
    pub absorbed_fraction: f64,
    pub balance_residual_w: f64,
    /// Nonzero-power directions for which the original response was evaluated.
    pub evaluated_directions: usize,
}

/// Original material failure, explicit input/work refusal, or cancellation.
#[derive(Clone, Debug, PartialEq)]
pub enum AngularPowerError {
    Material(MaterialError),
    Invalid(&'static str),
    Cancelled,
}
impl From<MaterialError> for AngularPowerError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl core::fmt::Display for AngularPowerError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "layer angular power refused: {self:?}")
    }
}
impl std::error::Error for AngularPowerError {}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AngularPowerError> {
    if cancelled() {
        Err(AngularPowerError::Cancelled)
    } else {
        Ok(())
    }
}

fn admit_samples(samples: usize, budget: AngularPowerBudget) -> Result<(), AngularPowerError> {
    if budget.max_samples == 0
        || budget.max_samples > MAX_ANGULAR_POWER_SAMPLES
        || samples == 0
        || samples > budget.max_samples
    {
        Err(AngularPowerError::Invalid(
            "angular power needs nonempty samples within the explicit and hard sample budgets",
        ))
    } else {
        Ok(())
    }
}

// One allocation-free summation owner serves both terminations and both
// angular measures. For diffuse quadrature the weights are dimensionless;
// otherwise they are the caller's already projected normal powers [W].
fn integrate(
    samples: impl ExactSizeIterator<Item = IncidentPower>,
    budget: AngularPowerBudget,
    declared_power_w: Option<f64>,
    cancelled: &mut impl FnMut() -> bool,
    mut response: impl FnMut(f64) -> Result<[f64; 3], MaterialError>,
) -> Result<AngularPower, AngularPowerError> {
    poll(cancelled)?;
    admit_samples(samples.len(), budget)?;
    let mut sums = [Dd::ZERO; 4];
    let mut evaluated = 0;
    for sample in samples {
        poll(cancelled)?;
        if !sample.cosine.is_finite()
            || sample.cosine <= 0.0
            || sample.cosine > 1.0
            || !sample.power_w.is_finite()
            || sample.power_w < 0.0
        {
            return Err(AngularPowerError::Invalid(
                "incident angular power needs cosine in (0,1] and finite nonnegative watts",
            ));
        }
        if sample.power_w == 0.0 {
            continue;
        }
        let fractions = response(sample.cosine)?;
        sums[0] = sums[0] + Dd::from_f64(sample.power_w);
        for (sum, fraction) in sums[1..].iter_mut().zip(fractions) {
            let power = sample.power_w * fraction;
            if !power.is_finite() {
                return Err(AngularPowerError::Invalid(
                    "angular power contribution is unrepresentable",
                ));
            }
            *sum = *sum + Dd::from_f64(power);
        }
        evaluated += 1;
    }
    let weights = sums.map(Dd::to_f64);
    if !weights.into_iter().all(f64::is_finite) || weights[0] <= 0.0 {
        return Err(AngularPowerError::Invalid(
            "angular power needs a positive representable total incident weight",
        ));
    }
    let fractions = [
        weights[1] / weights[0],
        weights[2] / weights[0],
        weights[3] / weights[0],
    ];
    let powers = if let Some(incident) = declared_power_w {
        [
            incident,
            incident * fractions[0],
            incident * fractions[1],
            incident * fractions[2],
        ]
    } else {
        weights
    };
    let residual = (Dd::from_f64(powers[0])
        - Dd::from_f64(powers[1])
        - Dd::from_f64(powers[2])
        - Dd::from_f64(powers[3]))
    .to_f64();
    if !powers.into_iter().all(f64::is_finite) || !residual.is_finite() {
        return Err(AngularPowerError::Invalid(
            "integrated angular power is unrepresentable",
        ));
    }
    poll(cancelled)?;
    Ok(AngularPower {
        incident_power_w: powers[0],
        reflected_power_w: powers[1],
        transmitted_power_w: powers[2],
        absorbed_power_w: powers[3],
        reflected_fraction: fractions[0],
        transmitted_fraction: fractions[1],
        absorbed_fraction: fractions[2],
        balance_residual_w: residual,
        evaluated_directions: evaluated,
    })
}

fn diffuse(
    incident_power_w: f64,
    incidence: DiffuseIncidence,
    budget: AngularPowerBudget,
    cancelled: &mut impl FnMut() -> bool,
    response: impl FnMut(f64) -> Result<[f64; 3], MaterialError>,
) -> Result<AngularPower, AngularPowerError> {
    poll(cancelled)?;
    admit_samples(incidence.polar_samples, budget)?;
    if !incident_power_w.is_finite()
        || incident_power_w < 0.0
        || !incidence.minimum_cosine.is_finite()
        || incidence.minimum_cosine < 0.0
        || incidence.minimum_cosine >= 1.0
    {
        return Err(AngularPowerError::Invalid(
            "diffuse incidence needs finite nonnegative watts and minimum cosine in [0,1)",
        ));
    }
    let width = (1.0 - incidence.minimum_cosine) / incidence.polar_samples as f64;
    integrate(
        (0..incidence.polar_samples).map(|i| {
            let cosine = incidence.minimum_cosine + (i as f64 + 0.5) * width;
            // The common panel width and azimuthal 2*pi cancel in the
            // normalized incident-power average. Keeping mu itself avoids
            // dividing by the tiny solid angle of a narrow declared cone.
            IncidentPower {
                cosine,
                power_w: cosine,
            }
        }),
        budget,
        Some(incident_power_w),
        cancelled,
        response,
    )
}

impl RigidBackedLayers<'_> {
    /// Incoherently sum the authored normal incident powers. The returned
    /// reflected and absorbed watts use the original angle-resolved solution;
    /// transmission is exactly zero. All-zero input refuses because it defines
    /// no angular power distribution. No material is re-evaluated or allocated.
    pub fn angular_power(
        &self,
        samples: &[IncidentPower],
        budget: AngularPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularPower, AngularPowerError> {
        integrate(
            samples.iter().copied(),
            budget,
            None,
            &mut cancelled,
            |cosine| {
                let value = self.reflection(cosine)?;
                Ok([
                    value.pressure_reflection.norm_sq(),
                    0.0,
                    value.absorbed_fraction,
                ])
            },
        )
    }

    /// Integrate isotropic intensity over the explicitly declared hemisphere
    /// or cone. `incident_power_w` is total power through the face in that cone;
    /// the cosine projection is included once in the quadrature. Zero incident
    /// power gives zero watts while retaining the distribution's coefficients.
    pub fn diffuse_power(
        &self,
        incident_power_w: f64,
        incidence: DiffuseIncidence,
        budget: AngularPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularPower, AngularPowerError> {
        diffuse(
            incident_power_w,
            incidence,
            budget,
            &mut cancelled,
            |cosine| {
                let value = self.reflection(cosine)?;
                Ok([
                    value.pressure_reflection.norm_sq(),
                    0.0,
                    value.absorbed_fraction,
                ])
            },
        )
    }
}

impl TransmittingLayers<'_> {
    /// Incoherently sum authored normal powers through the SAME original
    /// incident/layer/receiving media. Evanescent receiving waves contribute
    /// zero transmitted normal power; they are never relabeled as absorption.
    pub fn angular_power(
        &self,
        samples: &[IncidentPower],
        budget: AngularPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularPower, AngularPowerError> {
        integrate(
            samples.iter().copied(),
            budget,
            None,
            &mut cancelled,
            |cosine| {
                let value = self.scattering(cosine)?;
                Ok([
                    value.reflected_fraction,
                    value.transmitted_fraction,
                    value.absorbed_fraction,
                ])
            },
        )
    }

    /// The same positive diffuse incident-power quadrature as the rigid owner,
    /// retaining reflected, transmitted and absorbed powers independently.
    /// Critical refraction is handled by the original scattering solve; resolve
    /// its angular neighborhood through the caller's polar refinement ladder.
    pub fn diffuse_power(
        &self,
        incident_power_w: f64,
        incidence: DiffuseIncidence,
        budget: AngularPowerBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularPower, AngularPowerError> {
        diffuse(
            incident_power_w,
            incidence,
            budget,
            &mut cancelled,
            |cosine| {
                let value = self.scattering(cosine)?;
                Ok([
                    value.reflected_fraction,
                    value.transmitted_fraction,
                    value.absorbed_fraction,
                ])
            },
        )
    }
}
