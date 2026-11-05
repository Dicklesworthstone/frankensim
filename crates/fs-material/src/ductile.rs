//! Irreversible ductile-damage **initiation** from native plastic history.
//!
//! The accumulated criterion is `omega = integral d(alpha)/epsilon_f`, where
//! `alpha` is the native accumulated equivalent PLASTIC strain. Each proposed
//! increment holds its endpoint stress triaxiality and temperature and its
//! interval-averaged plastic strain rate fixed: `delta_omega=delta_alpha/epsilon_f`.
//! This is endpoint quadrature, not exact integration of an unresolved stress
//! path. Newton trials borrow the same accepted history; a consumer commits
//! mechanics and this history together only after all its acceptance gates.
//!
//! The Johnson--Cook convention here is
//! `(d1+d2*exp(-d3*eta))*(1+d4*ln(rate/reference_rate))*(1+d5*theta)`,
//! with tensile mean stress positive and `eta=mean(sigma)/von_mises(sigma)`.
//! This is the convention in SIMULIA, *Damage Initiation for Ductile Metals*,
//! <https://docs.software.vt.edu/abaqusv2025/English/SIMACAEMATRefMap/simamat-c-damageinitductile.htm>.
//! Its d3 sign is opposite to the original Johnson--Cook (1985) convention.
//! Parameters must be converted explicitly when importing that convention.
//!
//! Initiation is reached at `omega >= 1`; overshoot is retained. These laws
//! supply no post-initiation stress degradation, fracture energy, length-scale
//! regularization, crack topology, healing, necking solution or fragment state.
//! They omit Lode-angle dependence and cannot predict every shear failure.
//! Coupon elongation at break and cyclic fatigue-ductility coefficients are
//! not local equivalent plastic strain at damage initiation. Native plasticity
//! and material calibration retain their own applicability and evidence limits.

use crate::Voigt;
use fs_math::det;

/// Original material-card initiation curves and their retained source support.
pub mod source;

/// Progressive length-scaled ductile degradation on the native J2 owner.
pub mod evolution;

/// Complete rectangular calibration bounds. Every interval is closed; a
/// degenerate interval represents a fixed calibrated condition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileDomain {
    /// Tensile-positive mean stress divided by Mises equivalent stress.
    pub triaxiality: [f64; 2],
    /// Equivalent plastic strain rate [1/s], with strictly positive endpoints.
    pub strain_rate_s_inv: [f64; 2],
    /// Absolute temperature [K], with strictly positive endpoints.
    pub temperature_k: [f64; 2],
}

impl DuctileDomain {
    fn validate(&self) -> Result<(), DuctileError> {
        if !ordered(self.triaxiality)
            || !ordered(self.strain_rate_s_inv)
            || self.strain_rate_s_inv[0] <= 0.0
            || !ordered(self.temperature_k)
            || self.temperature_k[0] <= 0.0
        {
            return Err(DuctileError::Invalid(
                "calibration intervals must be finite and ordered; rate and temperature must be positive",
            ));
        }
        Ok(())
    }

    fn check_temperature(&self, temperature_k: f64) -> Result<(), DuctileError> {
        in_domain("temperature_k", temperature_k, self.temperature_k)
    }

    fn check(&self, eta: f64, rate: f64, temperature_k: f64) -> Result<(), DuctileError> {
        in_domain("triaxiality", eta, self.triaxiality)?;
        in_domain("plastic_strain_rate_s_inv", rate, self.strain_rate_s_inv)?;
        self.check_temperature(temperature_k)
    }
}

/// A local equivalent-plastic-strain initiation locus, independent of the
/// constitutive owner that produces plastic strain and physical stress.
pub trait DuctileInitiationLaw {
    /// The original admitted bounds, including temperature during elastic holds.
    fn domain(&self) -> &DuctileDomain;

    /// Return positive, finite initiation strain or refuse unsupported conditions.
    /// No implementation may silently extrapolate outside its original support.
    fn fracture_strain(
        &self,
        triaxiality: f64,
        strain_rate_s_inv: f64,
        temperature_k: f64,
    ) -> Result<f64, DuctileError>;
}

impl<L: DuctileInitiationLaw + ?Sized> DuctileInitiationLaw for &L {
    fn domain(&self) -> &DuctileDomain {
        (**self).domain()
    }

    fn fracture_strain(&self, eta: f64, rate: f64, temperature: f64) -> Result<f64, DuctileError> {
        (**self).fracture_strain(eta, rate, temperature)
    }
}

/// Explicit coefficients of the `exp(-d3*eta)` Johnson--Cook initiation law.
/// Signed coefficients are allowed only when all three factors stay positive
/// and finite throughout the supplied calibration domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JohnsonCookParameters {
    pub d1: f64,
    pub d2: f64,
    pub d3: f64,
    pub d4: f64,
    pub d5: f64,
    /// Reference equivalent plastic strain rate [1/s].
    pub reference_strain_rate_s_inv: f64,
    /// At or below this temperature, the normalized temperature is zero [K].
    pub transition_temperature_k: f64,
    /// At or above this temperature, the normalized temperature is one [K].
    /// This algebraic clamp does not supply a solid/liquid phase model.
    pub melting_temperature_k: f64,
}

/// Immutable phenomenological initiation law with no implicit material defaults.
#[derive(Clone, Debug, PartialEq)]
pub struct JohnsonCookInitiation {
    parameters: JohnsonCookParameters,
    domain: DuctileDomain,
}

impl JohnsonCookInitiation {
    /// Admit the coefficients over the COMPLETE requested rectangular domain.
    /// Each separate factor is monotone (possibly constant), so both endpoints
    /// suffice for factor positivity; all eight corners check representability
    /// of their product. No coefficient fit or physical validation is inferred.
    pub fn new(
        parameters: JohnsonCookParameters,
        domain: DuctileDomain,
    ) -> Result<Self, DuctileError> {
        domain.validate()?;
        let p = parameters;
        if [p.d1, p.d2, p.d3, p.d4, p.d5]
            .iter()
            .any(|x| !x.is_finite())
            || !positive(p.reference_strain_rate_s_inv)
            || !positive(p.transition_temperature_k)
            || !positive(p.melting_temperature_k)
            || p.melting_temperature_k <= p.transition_temperature_k
        {
            return Err(DuctileError::Invalid(
                "Johnson-Cook coefficients must be finite, reference rate positive, and 0 < transition < melting temperature",
            ));
        }
        let law = Self { parameters, domain };
        for eta in domain.triaxiality {
            for rate in domain.strain_rate_s_inv {
                for temperature in domain.temperature_k {
                    law.evaluate(eta, rate, temperature)?;
                }
            }
        }
        Ok(law)
    }

    /// Original coefficients, including the stated d3 sign convention.
    pub const fn parameters(&self) -> &JohnsonCookParameters {
        &self.parameters
    }

    /// Original calibration domain; construction does not expand it.
    pub const fn domain(&self) -> &DuctileDomain {
        &self.domain
    }

    fn evaluate(&self, eta: f64, rate: f64, temperature: f64) -> Result<f64, DuctileError> {
        let p = self.parameters;
        let stress_factor = if p.d2 == 0.0 {
            p.d1
        } else {
            let exponent = -p.d3 * eta;
            if !exponent.is_finite() {
                return Err(DuctileError::Numerical("Johnson-Cook exponent is not representable"));
            }
            p.d1 + p.d2 * det::exp(exponent)
        };
        // Taking the ratio first preserves nearby rates. The difference of
        // logs handles widely separated finite rates without ratio overflow.
        let rate_factor = if p.d4 == 0.0 {
            1.0
        } else {
            let ratio = rate / p.reference_strain_rate_s_inv;
            let logarithm = if positive(ratio) {
                det::ln(ratio)
            } else {
                det::ln(rate) - det::ln(p.reference_strain_rate_s_inv)
            };
            1.0 + p.d4 * logarithm
        };
        let theta = if temperature <= p.transition_temperature_k {
            0.0
        } else if temperature >= p.melting_temperature_k {
            1.0
        } else {
            (temperature - p.transition_temperature_k)
                / (p.melting_temperature_k - p.transition_temperature_k)
        };
        let temperature_factor = 1.0 + p.d5 * theta;
        if !positive(stress_factor) || !positive(rate_factor) || !positive(temperature_factor) {
            return Err(DuctileError::Invalid(
                "every Johnson-Cook factor must remain finite and strictly positive over its calibration domain",
            ));
        }
        let strain = (stress_factor * rate_factor) * temperature_factor;
        if !positive(strain) {
            return Err(DuctileError::Numerical("initiation strain is not representable"));
        }
        Ok(strain)
    }
}

impl DuctileInitiationLaw for JohnsonCookInitiation {
    fn domain(&self) -> &DuctileDomain {
        &self.domain
    }

    fn fracture_strain(&self, eta: f64, rate: f64, temperature: f64) -> Result<f64, DuctileError> {
        self.domain.check(eta, rate, temperature)?;
        self.evaluate(eta, rate, temperature)
    }
}

/// Piecewise-linear initiation strain versus triaxiality, with explicitly
/// bounded rate/temperature applicability. It makes no rate or temperature
/// dependence beyond the original supplied support and no extrapolation.
#[derive(Clone, Debug, PartialEq)]
pub struct TabulatedDuctileLocus {
    knots: Vec<(f64, f64)>,
    domain: DuctileDomain,
}

impl TabulatedDuctileLocus {
    /// Knots are `(triaxiality, local equivalent plastic strain at initiation)`.
    /// There must be at least two strictly increasing, finite triaxialities and
    /// positive finite strains; the domain must lie within the knot hull.
    pub fn new(knots: Vec<(f64, f64)>, domain: DuctileDomain) -> Result<Self, DuctileError> {
        domain.validate()?;
        if knots.len() < 2
            || knots.iter().any(|(eta, strain)| !eta.is_finite() || !positive(*strain))
            || knots.windows(2).any(|pair| {
                let width = pair[1].0 - pair[0].0;
                !positive(width)
            })
        {
            return Err(DuctileError::Invalid(
                "ductile locus needs at least two strictly ordered finite knots and positive finite initiation strains",
            ));
        }
        if domain.triaxiality[0] < knots[0].0
            || domain.triaxiality[1] > knots[knots.len() - 1].0
        {
            return Err(DuctileError::Invalid("triaxiality domain exceeds the original curve hull"));
        }
        Ok(Self { knots, domain })
    }

    /// Unmodified supplied curve, in its original ordered knot sequence.
    pub fn knots(&self) -> &[(f64, f64)] {
        &self.knots
    }

    /// Explicit calibration support, never extended by endpoint clamping.
    pub const fn domain(&self) -> &DuctileDomain {
        &self.domain
    }
}

impl DuctileInitiationLaw for TabulatedDuctileLocus {
    fn domain(&self) -> &DuctileDomain {
        &self.domain
    }

    fn fracture_strain(&self, eta: f64, rate: f64, temperature: f64) -> Result<f64, DuctileError> {
        self.domain.check(eta, rate, temperature)?;
        let upper = self.knots.partition_point(|(x, _)| *x < eta);
        if self.knots[upper].0 == eta {
            return Ok(self.knots[upper].1);
        }
        // Domain admission and the original knot hull guarantee 0 < upper < n.
        let (a, fa) = self.knots[upper - 1];
        let (b, fb) = self.knots[upper];
        let fraction = (eta - a) / (b - a);
        let strain = (1.0 - fraction) * fa + fraction * fb;
        if !positive(strain) {
            return Err(DuctileError::Numerical("interpolated initiation strain is not representable"));
        }
        Ok(strain)
    }
}

/// Accepted local initiation history. Only evaluated trials can produce a
/// non-virgin value; a constitutive alpha cannot reconstruct lost ductility
/// history without its original stress, rate, temperature and calibration.
/// Consumers retain the SAME immutable law/source alongside this state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileState {
    accumulated_plastic_strain: f64,
    utilization: f64,
}

impl Default for DuctileState {
    fn default() -> Self {
        Self::virgin()
    }
}

/// Immutable local proposal. A consumer must inspect initiation before
/// publishing a mechanical endpoint that assumes undamaged material response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DuctileStep {
    pub state: DuctileState,
    /// Original native plastic increment; creep strain is never substituted.
    pub plastic_increment: f64,
    /// Interval-average equivalent plastic strain rate, zero on elastic steps.
    pub strain_rate_s_inv: f64,
    pub increment_utilization: f64,
    /// Undefined and unnecessary for zero plastic flow.
    pub triaxiality: Option<f64>,
    /// Not evaluated for zero plastic flow.
    pub fracture_strain: Option<f64>,
    /// Fraction of this PLASTIC increment that would reach omega=1 if its
    /// endpoint conditions were held fixed. This is not a resolved event time
    /// or a substep equilibrium solve. Already-initiated history returns zero.
    pub initiation_fraction: Option<f64>,
}

impl DuctileState {
    /// A material point with no prior plastic flow or consumed ductility.
    pub const fn virgin() -> Self {
        Self { accumulated_plastic_strain: 0.0, utilization: 0.0 }
    }

    pub const fn accumulated_plastic_strain(&self) -> f64 {
        self.accumulated_plastic_strain
    }

    /// Accumulated dimensionless initiation criterion, without clipping at one.
    pub const fn utilization(&self) -> f64 {
        self.utilization
    }

    /// Damage initiation, not complete fracture or zero carrying capacity.
    pub fn initiated(&self) -> bool {
        self.utilization >= 1.0
    }

    /// Evaluate endpoint quadrature from an accepted native plastic history.
    ///
    /// `endpoint_alpha` must be the accumulated equivalent PLASTIC strain, not
    /// signed plastic strain, total strain, creep or a cycle amplitude. `dt_s`
    /// is the actual positive physical duration; no pseudo-time is inferred.
    /// Stress is the actual symmetric Cauchy or Kirchhoff tensor, NOT effective
    /// stress minus backstress or first Piola stress. Their common positive
    /// Jacobian scaling leaves triaxiality unchanged. Rate/eta support is checked
    /// only when plastic flow occurs; temperature support applies on every step.
    /// All refusals leave this accepted state unchanged.
    pub fn trial<L: DuctileInitiationLaw + ?Sized>(
        &self,
        law: &L,
        endpoint_alpha: f64,
        dt_s: f64,
        temperature_k: f64,
        stress_pa: &Voigt,
    ) -> Result<DuctileStep, DuctileError> {
        if !endpoint_alpha.is_finite() || endpoint_alpha < 0.0
            || !positive(dt_s) || !positive(temperature_k)
            || stress_pa.iter().any(|x| !x.is_finite())
        {
            return Err(DuctileError::Invalid(
                "endpoint alpha must be finite and nonnegative, stress finite, and duration/temperature positive",
            ));
        }
        if endpoint_alpha < self.accumulated_plastic_strain {
            return Err(DuctileError::NonMonotonePlasticStrain {
                previous: self.accumulated_plastic_strain, proposed: endpoint_alpha,
            });
        }
        law.domain().check_temperature(temperature_k)?;
        let plastic_increment = endpoint_alpha - self.accumulated_plastic_strain;
        if plastic_increment == 0.0 {
            return Ok(DuctileStep {
                state: *self, plastic_increment, strain_rate_s_inv: 0.0,
                increment_utilization: 0.0, triaxiality: None, fracture_strain: None,
                initiation_fraction: if self.initiated() { Some(0.0) } else { None },
            });
        }
        let rate = plastic_increment / dt_s;
        if !positive(rate) {
            return Err(DuctileError::Numerical("plastic strain rate is not representable"));
        }
        let eta = stress_triaxiality(stress_pa)?.ok_or(DuctileError::HydrostaticPlasticFlow)?;
        law.domain().check(eta, rate, temperature_k)?;
        let strain = law.fracture_strain(eta, rate, temperature_k)?;
        if !positive(strain) {
            return Err(DuctileError::Invalid("initiation law returned nonpositive or nonfinite strain"));
        }
        let increment_utilization = plastic_increment / strain;
        let utilization = self.utilization + increment_utilization;
        if !positive(increment_utilization) || !utilization.is_finite() {
            return Err(DuctileError::Numerical("ductility accumulation is not representable"));
        }
        let initiation_fraction = if self.initiated() {
            Some(0.0)
        } else if utilization >= 1.0 {
            Some(((1.0 - self.utilization) / increment_utilization).min(1.0))
        } else {
            None
        };
        Ok(DuctileStep {
            state: Self { accumulated_plastic_strain: endpoint_alpha, utilization },
            plastic_increment, strain_rate_s_inv: rate, increment_utilization,
            triaxiality: Some(eta), fracture_strain: Some(strain), initiation_fraction,
        })
    }
}

/// Tensile-positive stress triaxiality from tensor-Voigt Cauchy/Kirchhoff stress.
/// Returns `None` for exactly hydrostatic stress (including zero). Two scaling
/// stages avoid overflow of stress squares and underflow of a small deviator.
/// No arbitrary equivalent-stress threshold or compressive cutoff is applied.
pub fn stress_triaxiality(stress_pa: &Voigt) -> Result<Option<f64>, DuctileError> {
    if stress_pa.iter().any(|x| !x.is_finite()) {
        return Err(DuctileError::Invalid("stress tensor must be finite"));
    }
    if stress_pa[0] == stress_pa[1] && stress_pa[1] == stress_pa[2]
        && stress_pa[3..].iter().all(|x| *x == 0.0)
    {
        return Ok(None);
    }
    let scale = stress_pa.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    let s = stress_pa.map(|x| x / scale);
    let differences = [s[0] - s[1], s[1] - s[2], s[2] - s[0]];
    let deviator_scale = differences.iter().chain(&s[3..])
        .fold(0.0_f64, |m, x| m.max(x.abs()));
    if deviator_scale == 0.0 {
        return Err(DuctileError::Numerical("non-hydrostatic stress ratio is not representable"));
    }
    let diagonal = differences.map(|x| x / deviator_scale);
    let shear = [s[3], s[4], s[5]].map(|x| x / deviator_scale);
    let q = deviator_scale * det::sqrt(
        0.5 * diagonal.iter().map(|x| x * x).sum::<f64>()
            + 3.0 * shear.iter().map(|x| x * x).sum::<f64>(),
    );
    let eta = ((s[0] + s[1] + s[2]) / 3.0) / q;
    if !eta.is_finite() {
        return Err(DuctileError::Numerical("stress triaxiality is not representable"));
    }
    Ok(Some(eta))
}

/// Typed refusal without mutation of accepted material or mechanics history.
#[derive(Clone, Debug, PartialEq)]
pub enum DuctileError {
    Invalid(&'static str),
    Domain { quantity: &'static str, value: f64, bounds: [f64; 2] },
    NonMonotonePlasticStrain { previous: f64, proposed: f64 },
    HydrostaticPlasticFlow,
    Numerical(&'static str),
}

impl core::fmt::Display for DuctileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "ductile initiation refused: {self:?}")
    }
}
impl std::error::Error for DuctileError {}

fn positive(value: f64) -> bool { value.is_finite() && value > 0.0 }
fn ordered(bounds: [f64; 2]) -> bool {
    bounds.iter().all(|x| x.is_finite()) && bounds[0] <= bounds[1]
}
fn in_domain(quantity: &'static str, value: f64, bounds: [f64; 2]) -> Result<(), DuctileError> {
    if !value.is_finite() || value < bounds[0] || value > bounds[1] {
        return Err(DuctileError::Domain { quantity, value, bounds });
    }
    Ok(())
}

#[cfg(test)]
mod tests;
