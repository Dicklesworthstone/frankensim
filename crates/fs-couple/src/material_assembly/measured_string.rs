//! Measured whole-string sections in the original prestressed-beam owner.
//!
//! A commercial wound string may have source-supplied line mass, flexural
//! rigidity, diameter and operating tension while its construction remains
//! undisclosed. Those coefficients define the transverse continuum. Converting
//! them to an equivalent bulk density or Young modulus would invent axial and
//! torsional behavior, so this adapter never does so. Source selectors, state,
//! uncertainty and measurement receipts remain attached without promotion.
//!
//! The admitted law is a frozen, uniform scalar Euler--Bernoulli string with
//! measured fixed tension and pinned ends. It does not infer intrinsic loss
//! from a mean Q, reconstruct a winding, or supply shear/axial/thermal response.
//! Intrinsic loss exists only as source rows selected by name, through
//! [`ResolvedMeasuredStringSpecimen::with_loss_factors`].
//! Gas drag and reciprocal supports remain independent physical consumers.
//! Adding independent fluid inertia requires source-qualified dry structural
//! coefficients: effective mass or stiffness inferred from in-air modes is
//! not de-embedded by this adapter. The original source qualification remains
//! authoritative; retaining coefficients does not promote how they were obtained.

use crate::acoustic_realize::AcousticRealizeError;
use fs_matdb::EvaluationDecision;
use fs_material::state_point::{
    MaterialStatePointError, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement,
};
use fs_qty::{Dims, Force, Length, QuantitySpec};
use fs_scenario::PrestressedString;

/// Linear mass density [kg/m], distinct from volumetric material density.
pub const LINEAR_DENSITY_DIMS: Dims = Dims([-1, 1, 0, 0, 0, 0]);
/// Whole-section flexural rigidity EI [N m²], distinct from Young modulus.
pub const BENDING_STIFFNESS_DIMS: Dims = Dims([3, 1, -2, 0, 0, 0]);

/// Explicit property roles in one original measured specimen's source bundle.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MeasuredStringProperties {
    /// Loaded test span [m], or the speaking length bound by
    /// [`ResolvedMeasuredStringSpecimen::on_span_from_speaking_length`].
    pub length: String,
    /// Finished outside diameter [m].
    pub diameter: String,
    /// Source operating force [N].
    pub tension: String,
    /// Whole-string line density [kg/m].
    pub linear_density: String,
    /// Whole-string flexural rigidity [N m²].
    pub bending_stiffness: String,
}

impl MeasuredStringProperties {
    fn roles(&self) -> [(&str, Dims); 5] {
        [
            (&self.length, Length::DIMS),
            (&self.diameter, Length::DIMS),
            (&self.tension, Force::DIMS),
            (&self.linear_density, LINEAR_DENSITY_DIMS),
            (&self.bending_stiffness, BENDING_STIFFNESS_DIMS),
        ]
    }

    /// The exact source query required by this transverse-only consumer.
    /// No bulk modulus, shear law, viscosity or surface identity is requested.
    pub fn requirements(&self) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
        let roles = self.roles();
        let mut requirements = Vec::with_capacity(roles.len());
        for (index, &(name, dims)) in roles.iter().enumerate() {
            if roles[..index].iter().any(|(prior, _)| *prior == name) {
                return Err(MaterialStatePointError::InvalidDerived {
                    quantity: "measured string roles must select distinct physical properties",
                });
            }
            requirements.push(ScalarPropertyRequirement::try_new(
                name,
                dims,
                ScalarAdmissibility::StrictlyPositive,
            )?);
        }
        Ok(requirements)
    }
}

/// Source rows of one string's two constant loss factors.
///
/// Woodhouse, Plucked guitar transients, Acta Acustica united with Acustica 90
/// (2004) 945-965, eq. (8), fits a string's loss factor as
/// `[T (eta_F + eta_A / omega) + B eta_B k²] / [T + B k²]`. These selectors name
/// the dimensionless friction and bending factors. Its `eta_A` has no role: the
/// renderer computes air drag from the diameter and the declared gas, so adding
/// a fitted air term would count air twice, and the paper calls that term's form
/// and origin open to question. The computed drag does not reproduce it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MeasuredStringLossProperties {
    /// Dimensionless tension loss factor `eta_F`.
    pub friction: String,
    /// Dimensionless bending loss factor `eta_B`.
    pub bending: String,
}

impl MeasuredStringLossProperties {
    /// The exact source query for the two factors, friction then bending.
    pub fn requirements(&self) -> Result<Vec<ScalarPropertyRequirement>, MaterialStatePointError> {
        if self.friction == self.bending {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "measured string loss roles must select distinct physical properties",
            });
        }
        [&self.friction, &self.bending]
            .into_iter()
            .map(|name| {
                ScalarPropertyRequirement::try_new(
                    name,
                    Dims::NONE,
                    ScalarAdmissibility::StrictlyPositive,
                )
            })
            .collect()
    }
}

/// One measured uniform string at its actual source conditions and tension.
/// Missing internal loss is an explicit model limitation, never a zero-valued
/// material claim, until [`Self::with_loss_factors`] selects source rows.
/// The owning source variant refuses unsupported extra laws.
#[derive(Clone, Debug)]
pub struct ResolvedMeasuredStringSpecimen {
    measurements: ResolvedMaterialStatePoint,
    properties: MeasuredStringProperties,
    string: PrestressedString,
    omega_band_rad_s: Option<(f64, f64)>,
    speaking_length_m: Option<f64>,
    loss_properties: Option<MeasuredStringLossProperties>,
}

impl ResolvedMeasuredStringSpecimen {
    /// Admit source scalar coefficients without reconstructing hidden layers.
    ///
    /// Every role is a positive dimension-only SI scalar that is constant in
    /// its source validity. Sampled frequency/state curves are refused. Any
    /// declared frequency bands are intersected and checked on every retained
    /// mode; absence of a band remains unknown rather than becoming evidence
    /// for an unlimited physical range. Numerical mode/clock limits still apply.
    pub fn try_new(
        n_modes: usize,
        measurements: &ResolvedMaterialStatePoint,
        properties: MeasuredStringProperties,
    ) -> Result<Self, AcousticRealizeError> {
        properties
            .requirements()
            .map_err(|_| invalid("invalid measured string property selectors"))?;
        if n_modes == 0 || n_modes > fs_modal::time::MAX_TIME_DOMAIN_ACOUSTIC_MODES {
            return Err(invalid(
                "measured string requires the original bounded nonempty modal basis",
            ));
        }
        let mut values = [0.0; 5];
        let mut band: Option<(f64, f64)> = None;
        for (index, (name, dims)) in properties.roles().into_iter().enumerate() {
            values[index] = source_scalar(measurements, name, dims, &mut band)?;
        }
        admit_band(band)?;
        let [
            length_m,
            width_m,
            tension_n,
            lin_density_kg_m,
            bending_stiffness_n_m2,
        ] = values;
        let string = PrestressedString {
            length_m,
            width_m,
            tension_n,
            lin_density_kg_m,
            bending_stiffness_n_m2,
            n_modes,
            // No measured EA is available. Zero disables that channel in the
            // linear descriptor; the source owner refuses requests for it.
            axial_stiffness_n: 0.0,
            damping_ratio: 0.0,
            rayleigh: None,
            kelvin_voigt_bending: None,
            relaxation_bending: None,
            polarization_detune: 0.0,
            moving_end: false,
            tension_loss_factor: None,
            bending_loss_factor: None,
        };
        if !positive(0.5 * width_m) || !positive(lin_density_kg_m * length_m) {
            return Err(invalid(
                "measured string geometry and mass are not representable",
            ));
        }
        let result = Self {
            measurements: measurements.clone(),
            properties,
            string,
            omega_band_rad_s: band,
            speaking_length_m: None,
            loss_properties: None,
        };
        result.admit_retained_modes()?;
        Ok(result)
    }

    /// Bind a speaking-length measurement to the performance's full pinned span.
    ///
    /// Some sources measure over the speaking length, from the nut to an interior support (a guitar's scale
    /// length), not between the string's two pinned ends. The same uniform string continues past the support, so
    /// its measured tension, line mass and flexural rigidity hold on the whole `span_m`. The support's station
    /// times the span must reproduce the measured length to rounding: a span or station that disagrees refuses
    /// and is never rescaled. Every retained mode of the longer span is re-admitted against the source band.
    pub fn on_span_from_speaking_length(
        mut self,
        span_m: f64,
        support_fraction: f64,
    ) -> Result<Self, AcousticRealizeError> {
        if self.speaking_length_m.is_some() {
            return Err(invalid("measured speaking length is already bound to a span"));
        }
        if !positive(span_m) || !(support_fraction > 0.0 && support_fraction < 1.0) {
            return Err(invalid(
                "speaking length needs a positive span and an interior support station",
            ));
        }
        let measured = self.string.length_m;
        let speaking = span_m * support_fraction;
        if (speaking - measured).abs() > 8.0 * f64::EPSILON * speaking.max(measured) {
            return Err(invalid(
                "measured speaking length disagrees with the support station times the span",
            ));
        }
        self.string.length_m = span_m;
        self.speaking_length_m = Some(measured);
        self.admit_retained_modes()?;
        Ok(self)
    }

    /// Bind the string's own loss factors from rows of this measured source.
    ///
    /// The rows must already be resolved in this specimen's state point, so they carry the same receipts, and meet
    /// the same admission as every other role. Each factor's stated frequency band narrows the specimen's, and
    /// every retained mode is re-admitted against it.
    pub fn with_loss_factors(
        mut self,
        loss: MeasuredStringLossProperties,
    ) -> Result<Self, AcousticRealizeError> {
        if self.loss_properties.is_some() {
            return Err(invalid("measured string loss factors are already bound"));
        }
        loss.requirements()
            .map_err(|_| invalid("invalid measured string loss selectors"))?;
        let friction =
            source_scalar(&self.measurements, &loss.friction, Dims::NONE, &mut self.omega_band_rad_s)?;
        let bending =
            source_scalar(&self.measurements, &loss.bending, Dims::NONE, &mut self.omega_band_rad_s)?;
        admit_band(self.omega_band_rad_s)?;
        self.string.tension_loss_factor = Some(friction);
        self.string.bending_loss_factor = Some(bending);
        self.loss_properties = Some(loss);
        self.admit_retained_modes()?;
        Ok(self)
    }

    /// The source's measured speaking length once bound to a longer span; `None` for a test-span measurement.
    #[must_use]
    pub const fn speaking_length_m(&self) -> Option<f64> {
        self.speaking_length_m
    }

    /// The selected loss rows; `None` is the explicit no-intrinsic-loss limitation.
    #[must_use]
    pub const fn loss_properties(&self) -> Option<&MeasuredStringLossProperties> {
        self.loss_properties.as_ref()
    }

    /// Original uniform continuum, with no manufactured constitutive coefficients.
    #[must_use]
    pub fn string(&self) -> PrestressedString {
        self.string.clone()
    }

    /// Complete immutable material-card query, uncertainty and usage receipts.
    #[must_use]
    pub const fn measurements(&self) -> &ResolvedMaterialStatePoint {
        &self.measurements
    }

    /// Exact source-role selection; property names never imply a material law.
    #[must_use]
    pub const fn properties(&self) -> &MeasuredStringProperties {
        &self.properties
    }

    /// Actual measured outside radius used by the shared fluid boundary problem.
    #[must_use]
    pub fn radius_m(&self) -> f64 {
        0.5 * self.string.width_m
    }

    /// Intersected source applicability when stated. `None` means unqualified.
    #[must_use]
    pub const fn omega_band_rad_s(&self) -> Option<(f64, f64)> {
        self.omega_band_rad_s
    }

    /// Shared numerical/source-band admission, including after fluid loading.
    pub fn admit_angular_frequency(&self, omega: f64) -> Result<(), AcousticRealizeError> {
        if !positive(omega)
            || self
                .omega_band_rad_s
                .is_some_and(|(low, high)| omega < low || omega > high)
        {
            return Err(invalid(
                "measured string mode is unrepresentable or outside its source frequency band",
            ));
        }
        Ok(())
    }

    fn admit_retained_modes(&self) -> Result<(), AcousticRealizeError> {
        for mode in 1..=self.string.n_modes {
            self.admit_angular_frequency(fs_nlmodal::prestressed_beam_omega(
                self.string.length_m,
                self.string.tension_n,
                self.string.lin_density_kg_m,
                self.string.bending_stiffness_n_m2,
                mode,
            ))?;
        }
        Ok(())
    }
}

/// One source role's value, narrowing `band` by the frequency domain it states.
fn source_scalar(
    measurements: &ResolvedMaterialStatePoint,
    name: &str,
    dims: Dims,
    band: &mut Option<(f64, f64)>,
) -> Result<f64, AcousticRealizeError> {
    let property = measurements
        .property(name)
        .ok_or_else(|| invalid("measured string is missing an explicitly selected source property"))?;
    let value = property.value_si();
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !value.is_finite()
        || value <= 0.0
        || !property.source_axes().is_empty()
        || property.answer().receipt.decision != EvaluationDecision::ConstantWithinValidity
    {
        return Err(invalid(
            "measured string needs positive dimension-only SI scalars with validity-wide constant support",
        ));
    }
    if let Some(next) = property
        .angular_frequency_band_rad_s()
        .map_err(|_| invalid("measured string has invalid frequency applicability semantics"))?
    {
        *band = Some(band.map_or(next, |old| (old.0.max(next.0), old.1.min(next.1))));
    }
    Ok(value)
}

fn admit_band(band: Option<(f64, f64)>) -> Result<(), AcousticRealizeError> {
    if band.is_some_and(|(low, high)| {
        !low.is_finite() || !high.is_finite() || low < 0.0 || high <= low
    }) {
        return Err(invalid(
            "measured string source frequency domains have no nonempty intersection",
        ));
    }
    Ok(())
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

#[cfg(test)]
mod tests;
