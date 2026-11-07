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
//! internal refraction. Absorption is the positive local viscous/thermal work,
//! not a difference of nearly equal incident and outgoing powers. Each layer's
//! work can be retained for heating under the same incoherent angular measure.
//! It is not a locally reacting impedance fit or a
//! reverberation-room/finite-panel prediction. Polar resolution is explicit;
//! refinement is required near resonances and critical refraction. No
//! quadrature error bound or experimentally validated diffuse field is implied.

use super::{
    LayerAbsorbedPower, LayerAbsorption, LayerAbsorptionBudget, LayerAbsorptionError, MAX_LAYERS,
    RigidBackedLayers, TransmittingLayers,
};
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

/// Independent reflected/transmitted powers and physical heat in every layer.
/// Layer fractions are normalized by the TOTAL incident angular power, not by
/// an individual direction's power. Viscous and thermal mechanisms are those
/// of the original equivalent fluid; they do not assign heat to separate pore
/// fluid/frame temperatures or evolve a thermal state.
#[derive(Clone, Debug, PartialEq)]
pub struct AngularLayerAbsorption {
    pub power: AngularPower,
    layers: [LayerAbsorbedPower; MAX_LAYERS],
    layer_count: usize,
}
impl AngularLayerAbsorption {
    /// Original exposed-to-rear order, including zero-thickness/lossless rows.
    pub fn layers(&self) -> &[LayerAbsorbedPower] {
        &self.layers[..self.layer_count]
    }
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
impl From<LayerAbsorptionError> for AngularPowerError {
    fn from(error: LayerAbsorptionError) -> Self {
        match error {
            LayerAbsorptionError::Material(error) => Self::Material(error),
            LayerAbsorptionError::Invalid(what) => Self::Invalid(what),
            LayerAbsorptionError::Cancelled => Self::Cancelled,
        }
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
fn integrate<C: FnMut() -> bool>(
    samples: impl ExactSizeIterator<Item = IncidentPower>,
    budget: AngularPowerBudget,
    layer_count: usize,
    layer_budget: LayerAbsorptionBudget,
    declared_power_w: Option<f64>,
    cancelled: &mut C,
    mut response: impl FnMut(f64, &mut C) -> Result<LayerAbsorption, LayerAbsorptionError>,
) -> Result<AngularLayerAbsorption, AngularPowerError> {
    poll(cancelled)?;
    admit_samples(samples.len(), budget)?;
    if layer_budget.max_layers > MAX_LAYERS || layer_count > layer_budget.max_layers {
        return Err(AngularPowerError::Invalid(
            "angular absorption layer work budget exceeded",
        ));
    }
    let mut sums = [Dd::ZERO; 3];
    let mut layer_sums = [[Dd::ZERO; 2]; MAX_LAYERS];
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
        let value = response(sample.cosine, cancelled)?;
        sums[0] = sums[0] + Dd::from_f64(sample.power_w);
        for (sum, fraction) in sums[1..]
            .iter_mut()
            .zip([value.reflected_fraction, value.transmitted_fraction])
        {
            let power = sample.power_w * fraction;
            if !power.is_finite() {
                return Err(AngularPowerError::Invalid(
                    "angular power contribution is unrepresentable",
                ));
            }
            *sum = *sum + Dd::from_f64(power);
        }
        for (sum, layer) in layer_sums[..layer_count].iter_mut().zip(value.layers()) {
            for (mechanism, fraction) in sum.iter_mut().zip([
                layer.viscous_absorbed_fraction,
                layer.thermal_absorbed_fraction,
            ]) {
                let power = sample.power_w * fraction;
                if !power.is_finite() {
                    return Err(AngularPowerError::Invalid(
                        "angular layer absorption is unrepresentable",
                    ));
                }
                *mechanism = *mechanism + Dd::from_f64(power);
            }
        }
        evaluated += 1;
    }
    let absorbed = layer_sums[..layer_count]
        .iter()
        .fold(Dd::ZERO, |total, row| total + row[0] + row[1]);
    let weights = [sums[0], sums[1], sums[2], absorbed].map(Dd::to_f64);
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
    let mut layers = [LayerAbsorbedPower::default(); MAX_LAYERS];
    for (row, sum) in layers[..layer_count].iter_mut().zip(layer_sums) {
        let viscous = sum[0].to_f64() / weights[0];
        let thermal = sum[1].to_f64() / weights[0];
        let power = if let Some(incident) = declared_power_w {
            incident * (viscous + thermal)
        } else {
            (sum[0] + sum[1]).to_f64()
        };
        if !power.is_finite() {
            return Err(AngularPowerError::Invalid(
                "integrated angular layer power is unrepresentable",
            ));
        }
        *row = LayerAbsorbedPower {
            viscous_absorbed_fraction: viscous,
            thermal_absorbed_fraction: thermal,
            absorbed_power_w: power,
        };
    }
    poll(cancelled)?;
    Ok(AngularLayerAbsorption {
        power: AngularPower {
            incident_power_w: powers[0],
            reflected_power_w: powers[1],
            transmitted_power_w: powers[2],
            absorbed_power_w: powers[3],
            reflected_fraction: fractions[0],
            transmitted_fraction: fractions[1],
            absorbed_fraction: fractions[2],
            balance_residual_w: residual,
            evaluated_directions: evaluated,
        },
        layers,
        layer_count,
    })
}

fn diffuse<C: FnMut() -> bool>(
    incident_power_w: f64,
    incidence: DiffuseIncidence,
    budget: AngularPowerBudget,
    layer_count: usize,
    layer_budget: LayerAbsorptionBudget,
    cancelled: &mut C,
    response: impl FnMut(f64, &mut C) -> Result<LayerAbsorption, LayerAbsorptionError>,
) -> Result<AngularLayerAbsorption, AngularPowerError> {
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
        layer_count,
        layer_budget,
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<AngularPower, AngularPowerError> {
        self.angular_layer_absorption(
            samples,
            budget,
            LayerAbsorptionBudget {
                max_layers: self.layer_count(),
            },
            cancelled,
        )
        .map(|report| report.power)
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<AngularPower, AngularPowerError> {
        self.diffuse_layer_absorption(
            incident_power_w,
            incidence,
            budget,
            LayerAbsorptionBudget {
                max_layers: self.layer_count(),
            },
            cancelled,
        )
        .map(|report| report.power)
    }

    /// Integrate physical viscous/thermal work separately in every original
    /// layer. Uses one bounded propagation and one analytic loss pass per
    /// nonzero direction, with fixed storage and cancellation in both passes.
    /// Zero-power directions count against the budget but add no heat.
    pub fn angular_layer_absorption(
        &self,
        samples: &[IncidentPower],
        budget: AngularPowerBudget,
        layer_budget: LayerAbsorptionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularLayerAbsorption, AngularPowerError> {
        integrate(
            samples.iter().copied(),
            budget,
            self.layer_count(),
            layer_budget,
            None,
            &mut cancelled,
            |cosine, cancelled| self.layer_absorption(cosine, 1.0, layer_budget, cancelled),
        )
    }

    /// The same cosine-projected diffuse measure as [`Self::diffuse_power`],
    /// retaining each layer's positive work instead of inferring it from R/T.
    /// Zero incident watts retain the angularly averaged loss coefficients.
    pub fn diffuse_layer_absorption(
        &self,
        incident_power_w: f64,
        incidence: DiffuseIncidence,
        budget: AngularPowerBudget,
        layer_budget: LayerAbsorptionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularLayerAbsorption, AngularPowerError> {
        diffuse(
            incident_power_w,
            incidence,
            budget,
            self.layer_count(),
            layer_budget,
            &mut cancelled,
            |cosine, cancelled| self.layer_absorption(cosine, 1.0, layer_budget, cancelled),
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<AngularPower, AngularPowerError> {
        self.angular_layer_absorption(
            samples,
            budget,
            LayerAbsorptionBudget {
                max_layers: self.prepared_layers().layer_count(),
            },
            cancelled,
        )
        .map(|report| report.power)
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<AngularPower, AngularPowerError> {
        self.diffuse_layer_absorption(
            incident_power_w,
            incidence,
            budget,
            LayerAbsorptionBudget {
                max_layers: self.prepared_layers().layer_count(),
            },
            cancelled,
        )
        .map(|report| report.power)
    }

    /// Independent reflected/transmitted power and positive viscous/thermal
    /// heat per layer for authored incoherent directions. The original fluid
    /// snapshots are reused; no material query or allocation occurs per angle.
    /// The two native per-layer passes retain bounded cancellation latency.
    pub fn angular_layer_absorption(
        &self,
        samples: &[IncidentPower],
        budget: AngularPowerBudget,
        layer_budget: LayerAbsorptionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularLayerAbsorption, AngularPowerError> {
        integrate(
            samples.iter().copied(),
            budget,
            self.prepared_layers().layer_count(),
            layer_budget,
            None,
            &mut cancelled,
            |cosine, cancelled| self.layer_absorption(cosine, 1.0, layer_budget, cancelled),
        )
    }

    /// Resolve each layer's heat under isotropic intensity in the declared
    /// cone. Absorption below the precision of the boundary-power difference
    /// survives through the original analytic local-loss integrals. This adds
    /// no finite-panel, frame-motion, or transient thermal approximation.
    pub fn diffuse_layer_absorption(
        &self,
        incident_power_w: f64,
        incidence: DiffuseIncidence,
        budget: AngularPowerBudget,
        layer_budget: LayerAbsorptionBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AngularLayerAbsorption, AngularPowerError> {
        diffuse(
            incident_power_w,
            incidence,
            budget,
            self.prepared_layers().layer_count(),
            layer_budget,
            &mut cancelled,
            |cosine, cancelled| self.layer_absorption(cosine, 1.0, layer_budget, cancelled),
        )
    }
}
