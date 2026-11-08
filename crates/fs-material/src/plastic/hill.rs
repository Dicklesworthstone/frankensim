//! Hill's quadratic orthotropic metal yield surface with associated plastic flow.
//!
//! In the fixed material axes,
//! `q² = F(syy-szz)² + G(szz-sxx)² + H(sxx-syy)²`
//! `     + 2L syz² + 2M szx² + 2N sxy²`.
//! The yield radius is `yield_stress + hardening * alpha`, with an optional
//! energetic Voce increment `Q*(1-exp(-b*alpha))`. Here `alpha` is
//! the accumulated work-conjugate multiplier: `(sigma-X):d(eps_p) = q d(alpha)`
//! with backstress X=0 unless linear kinematic hardening is explicitly selected.
//! It is NOT generally the von Mises equivalent plastic strain. Individual
//! F/G/H coefficients may be negative provided the quadratic form is positive
//! definite on deviatoric stresses. Plastic flow preserves volume; hydrostatic
//! compression cannot by itself plastically deform this dense-metal model.
//!
//! The five deviatoric Hill modes are prepared once. Isotropic elasticity makes
//! the backward-Euler return a monotone scalar solve with analytical bounds,
//! followed by the exact algorithmic tangent. All point evaluations use fixed
//! arrays, without allocation, tensor Newton iteration or numerical derivatives.
//! The same `SmallStrainLaw` serves the native plane-stress, solid and dynamic
//! consumers. Material axes must rotate with the specimen under a change of
//! coordinates; this is a small-strain, small-rotation law, not texture evolution
//! or finite-strain forming. It supplies neither fracture nor forming limits.
//!
//! Optional energetic linear kinematic hardening translates this same surface
//! by `X = C * eps_p`: yielding uses `q(sigma-X) = yield_stress + R(alpha)`.
//! C multiplies TENSOR plastic strain, not engineering shear or a uniaxial
//! tangent modulus. Its recoverable storage is `C eps_p:eps_p / 2`.
//! Existing constructors select C=0; opt in before retaining any history.
//! Saturation uses the same scalar return and the endpoint hardening slope in
//! its exact tangent. Its recoverable storage is the integral of the Voce
//! increment; neither that storage nor integration loss becomes plastic heat.
//!
//! References: Hill (1948), doi:10.1098/rspa.1948.0045; Sandia LAME,
//! <https://www.sandia.gov/files/sierra/SM_LAME_5_28/main/models/hill.html>.
//! The Prager coefficient convention follows LAME's elastic-plastic model,
//! equation 4.27: <https://www.sandia.gov/files/sierra/SM_LAME_5_28/main/models/elastic_plastic.html>.

mod calibration;
/// Fixed material orientation without changing the local history convention.
pub mod oriented;
/// Explicit-time associative linear overstress on the same anisotropic return.
pub mod rate;
/// Original material-card support for the directional constitutive parameters.
pub mod source;
/// Source-bound thermal softening and reciprocal heat for fixed material axes.
pub mod thermal;

use super::energy::{EnergyTolerances, IncrementEnergy, StoredEnergy};
use super::voce::hardening;
use crate::elastic::{IsotropicElastic, bulk_shear_stiffness};
use crate::tensor::{contract, frobenius_norm, mandel_stiffness_to_voigt, trace};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use core::fmt;
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

const INV_SQRT_6: f64 = 0.4082482904638631;
const MAX_RETURN_ITERATIONS: usize = 128;

/// Dimensionless coefficients in the documented material-axis convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hill48Coefficients {
    /// Multiplier of `(sigma_yy-sigma_zz)^2`.
    pub f: f64,
    /// Multiplier of `(sigma_zz-sigma_xx)^2`.
    pub g: f64,
    /// Multiplier of `(sigma_xx-sigma_yy)^2`.
    pub h: f64,
    /// Half the multiplier of `sigma_yz^2`.
    pub l: f64,
    /// Half the multiplier of `sigma_zx^2`.
    pub m: f64,
    /// Half the multiplier of `sigma_xy^2`.
    pub n: f64,
}

impl Hill48Coefficients {
    /// The exact coefficient convention that recovers von Mises J2 yielding.
    pub const fn isotropic() -> Self {
        Self {
            f: 0.5,
            g: 0.5,
            h: 0.5,
            l: 1.5,
            m: 1.5,
            n: 1.5,
        }
    }

    /// Ratios in native `[xx, yy, zz, xy, yz, zx]` order. Normal yield stresses
    /// are relative to the reference yield stress; SHEAR yield stresses are
    /// relative to `reference_yield / sqrt(3)`. Thus six unit ratios recover J2.
    /// Directional yield data do not automatically identify plastic r-values.
    pub fn from_yield_ratios(ratios: [f64; 6]) -> Result<Self, HillError> {
        if ratios.iter().any(|x| !positive(*x)) {
            return Err(HillError::Parameters(
                "yield ratios must be finite and positive",
            ));
        }
        let inverse = ratios.map(|x| (1.0 / x).powi(2));
        let value = Self {
            f: 0.5 * (inverse[1] + inverse[2] - inverse[0]),
            g: 0.5 * (inverse[2] + inverse[0] - inverse[1]),
            h: 0.5 * (inverse[0] + inverse[1] - inverse[2]),
            l: 1.5 * inverse[4],
            m: 1.5 * inverse[5],
            n: 1.5 * inverse[3],
        };
        Modes::new(value)?;
        Ok(value)
    }
}

/// Accepted specimen history, in the original material axes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct HillState {
    /// Isochoric plastic strain, native tensor-Voigt `[xx,yy,zz,xy,yz,zx]`.
    pub plastic_strain: Voigt,
    /// Accumulated Hill work-conjugate multiplier, dimensionless.
    pub alpha: f64,
}

/// Immutable coefficients and once-prepared scalar return data.
#[derive(Clone, Debug, PartialEq)]
pub struct Hill48 {
    elastic: IsotropicElastic,
    yield_stress: f64,
    hardening: f64,
    kinematic: f64,
    saturation_pa: f64,
    saturation_rate_per_strain: f64,
    coefficients: Hill48Coefficients,
    bulk: f64,
    twice_shear: f64,
    return_twice_shear: f64,
    modes: Modes,
    return_scale: f64,
    return_ratios: [f64; 5],
}

/// One checked candidate; evaluating it does not mutate the accepted history.
#[derive(Clone, Debug, PartialEq)]
pub struct HillResponse {
    /// Candidate Cauchy stress in the fixed material axes [Pa].
    pub stress_pa: Voigt,
    /// Exact derivative with respect to the native TENSOR shear coordinates.
    pub tangent_pa: Tangent6,
    /// History to commit only after the surrounding physical step succeeds.
    pub state: HillState,
    /// Increment of the Hill work-conjugate accumulated strain.
    pub plastic_increment: f64,
    /// Signed endpoint `q(sigma-X) - yield_stress - R(alpha)`, in Pa.
    pub yield_residual_pa: f64,
}

/// Endpoint work, recoverable storage and separately identified physical heat.
#[derive(Clone, Debug, PartialEq)]
pub struct HillEnergyIncrement {
    /// Total candidate tensor strain in the fixed material axes.
    pub strain: Voigt,
    /// Original constitutive candidate; no second return is performed.
    pub response: HillResponse,
    /// Work and storage per reference volume, in J/m3.
    pub energy: IncrementEnergy,
}

// Stress-only line searches and accepted-history updates use the same checked
// return without constructing its tangent. Retain the scalar solve's data so
// a requested tangent still differentiates that exact candidate.
struct HillTrial {
    stress_pa: Voigt,
    state: HillState,
    plastic_increment: f64,
    yield_residual_pa: f64,
    tangent: Option<PlasticTangent>,
}

#[derive(Clone, Debug, PartialEq)]
struct LinearRateReturn {
    modulus_pa: f64,
    scale: f64,
    ratios: [f64; 5],
}

struct PlasticTangent {
    return_scale: f64,
    radius_ratio: f64,
    denominators: [f64; 5],
    normalized: [f64; 5],
    ratios: [f64; 5],
    derivative_ratios: [f64; 5],
}

/// Parameter, history, numerical or explicit work-account refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum HillError {
    /// Invalid or unrepresentable fixed material coefficients.
    Parameters(&'static str),
    /// Invalid retained history or unsupported deformation.
    State(&'static str),
    /// An intermediate or result cannot be represented by finite arithmetic.
    Numerical(&'static str),
    /// The safeguarded scalar return exhausted its bounded work allowance.
    NonConvergence {
        /// Number of scalar root iterations attempted.
        iterations: usize,
        /// Last dimensionless scalar yield residual.
        residual: f64,
    },
    /// An independently supplied yield gate was not satisfied.
    YieldResidual {
        /// Signed mismatch in Pa.
        residual_pa: f64,
        /// Allowed absolute mismatch in Pa.
        tolerance_pa: f64,
    },
    /// An independently supplied endpoint-work gate was not satisfied.
    EnergyResidual {
        /// Signed first-law mismatch in J/m3.
        residual_j_m3: f64,
        /// Allowed absolute mismatch in J/m3.
        tolerance_j_m3: f64,
    },
}

impl fmt::Display for HillError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parameters(s) => write!(f, "Hill parameters: {s}"),
            Self::State(s) => write!(f, "Hill state: {s}"),
            Self::Numerical(s) => write!(f, "Hill arithmetic: {s}"),
            Self::NonConvergence {
                iterations,
                residual,
            } => write!(
                f,
                "Hill return failed after {iterations} iterations; relative residual {residual}"
            ),
            Self::YieldResidual {
                residual_pa,
                tolerance_pa,
            } => write!(
                f,
                "Hill yield residual {residual_pa} Pa exceeds {tolerance_pa} Pa"
            ),
            Self::EnergyResidual {
                residual_j_m3,
                tolerance_j_m3,
            } => write!(
                f,
                "Hill work residual {residual_j_m3} J/m3 exceeds {tolerance_j_m3} J/m3"
            ),
        }
    }
}
impl std::error::Error for HillError {}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

fn norm<const N: usize>(v: &[f64; N]) -> f64 {
    if v.iter().any(|x| x.is_nan()) {
        return f64::NAN;
    }
    let scale = v.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    if scale == 0.0 || !scale.is_finite() {
        return scale;
    }
    scale * det::sqrt(v.iter().map(|x| (x / scale).powi(2)).sum())
}

/// Orthonormal deviatoric Mandel modes. Normal modes need only one 2x2 Jacobi
/// rotation; shear modes are already diagonal in the specified material frame.
#[derive(Clone, Debug, PartialEq)]
struct Modes {
    cosine: f64,
    sine: f64,
    eigenvalues: [f64; 5],
    roots: [f64; 5],
    vectors: [Voigt; 5],
}

impl Modes {
    fn new(coefficients: Hill48Coefficients) -> Result<Self, HillError> {
        let Hill48Coefficients { f, g, h, l, m, n } = coefficients;
        if [f, g, h, l, m, n].iter().any(|x| !x.is_finite())
            || !positive(l)
            || !positive(m)
            || !positive(n)
        {
            return Err(HillError::Parameters(
                "Hill coefficients must be finite with positive shear coefficients",
            ));
        }
        let scale = f.abs().max(g.abs()).max(h.abs());
        if !positive(scale) {
            return Err(HillError::Parameters(
                "normal yield surface must be positive definite on deviatoric stresses",
            ));
        }
        let (f, g, h) = (f / scale, g / scale, h / scale);
        let a = 0.5 * (f + g) + 2.0 * h;
        let d = 1.5 * (f + g);
        let b = 0.5 * det::sqrt(3.0) * (g - f);
        let determinant = 3.0 * (f * g + f * h + g * h);
        if !positive(a) || !positive(d) || !positive(determinant) {
            return Err(HillError::Parameters(
                "Hill normal quadratic form must be positive definite on deviatoric stresses",
            ));
        }
        let (cosine, sine, eigen) = if b == 0.0 {
            (1.0, 0.0, [a, d])
        } else {
            let delta = 0.5 * (d - a);
            let radius = norm(&[delta, b]);
            let t = b / if delta >= 0.0 {
                delta + radius
            } else {
                delta - radius
            };
            let cosine = 1.0 / det::sqrt(1.0 + t * t);
            let sine = t * cosine;
            // Recover the small eigenvalue by product/high, avoiding subtraction
            // of nearly equal diagonal values near an anisotropic degeneracy.
            let eigen = if delta >= 0.0 {
                let high = d + t * b;
                [determinant / high, high]
            } else {
                let high = a - t * b;
                [high, determinant / high]
            };
            (cosine, sine, eigen)
        };
        let eigenvalues = [scale * eigen[0], scale * eigen[1], n, l, m];
        if eigenvalues.iter().any(|x| !positive(*x)) {
            return Err(HillError::Parameters(
                "Hill modal coefficients are not representable",
            ));
        }
        let mut value = Self {
            cosine,
            sine,
            eigenvalues,
            roots: eigenvalues.map(det::sqrt),
            vectors: [[0.0; 6]; 5],
        };
        for mode in 0..5 {
            let mut unit = [0.0; 5];
            unit[mode] = 1.0;
            value.vectors[mode] = crate::tensor::to_mandel(&value.tensor(&unit));
        }
        Ok(value)
    }

    fn coordinates(&self, v: &Voigt) -> [f64; 5] {
        // Differences remove a hydrostatic offset before any projection.
        let a = (v[0] - v[1]) * core::f64::consts::FRAC_1_SQRT_2;
        let b = ((v[0] - v[2]) + (v[1] - v[2])) * INV_SQRT_6;
        [
            self.cosine * a - self.sine * b,
            self.sine * a + self.cosine * b,
            core::f64::consts::SQRT_2 * v[3],
            core::f64::consts::SQRT_2 * v[4],
            core::f64::consts::SQRT_2 * v[5],
        ]
    }

    fn tensor(&self, v: &[f64; 5]) -> Voigt {
        let a = (self.cosine * v[0] + self.sine * v[1]) * core::f64::consts::FRAC_1_SQRT_2;
        let b = (-self.sine * v[0] + self.cosine * v[1]) * INV_SQRT_6;
        let (xx, yy) = (a + b, -a + b);
        [
            xx,
            yy,
            -(xx + yy),
            v[2] * core::f64::consts::FRAC_1_SQRT_2,
            v[3] * core::f64::consts::FRAC_1_SQRT_2,
            v[4] * core::f64::consts::FRAC_1_SQRT_2,
        ]
    }

    fn equivalent(&self, v: &[f64; 5]) -> f64 {
        norm(&core::array::from_fn::<_, 5, _>(|i| self.roots[i] * v[i]))
    }
}

impl Hill48 {
    /// Admit fixed elastic/yield/hardening coefficients and a strictly convex
    /// deviatoric Hill surface. No scalar malleability, hardness conversion or
    /// tensile coupon elongation is substituted for these independent data.
    pub fn new(
        elastic: IsotropicElastic,
        yield_stress_pa: f64,
        hardening_modulus_pa: f64,
        coefficients: Hill48Coefficients,
    ) -> Result<Self, HillError> {
        if !positive(elastic.youngs)
            || !positive(elastic.strain_limit)
            || !(elastic.poisson.is_finite() && elastic.poisson > -1.0 && elastic.poisson < 0.5)
            || !positive(yield_stress_pa)
            || !(hardening_modulus_pa.is_finite() && hardening_modulus_pa >= 0.0)
        {
            return Err(HillError::Parameters(
                "elasticity, yield stress, hardening or strain domain is inadmissible",
            ));
        }
        let bulk = elastic.youngs / (3.0 * (1.0 - 2.0 * elastic.poisson));
        let twice_shear = elastic.youngs / (1.0 + elastic.poisson);
        let modes = Modes::new(coefficients)?;
        let returns = modes
            .eigenvalues
            .map(|v| hardening_modulus_pa + twice_shear * v);
        if !positive(bulk)
            || !positive(twice_shear)
            || returns.iter().any(|v| !positive(*v))
            || modes.eigenvalues.iter().any(|v| !positive(twice_shear * v))
            || bulk_shear_stiffness(bulk, 0.5 * twice_shear)
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(HillError::Parameters(
                "elastic or return moduli are not representable",
            ));
        }
        let return_scale = returns.iter().copied().fold(0.0_f64, f64::max);
        let return_ratios = returns.map(|v| v / return_scale);
        if return_ratios.iter().any(|v| !positive(*v)) {
            return Err(HillError::Parameters(
                "return modulus ratios are not representable",
            ));
        }
        Ok(Self {
            elastic,
            yield_stress: yield_stress_pa,
            hardening: hardening_modulus_pa,
            kinematic: 0.0,
            saturation_pa: 0.0,
            saturation_rate_per_strain: 1.0,
            coefficients,
            bulk,
            twice_shear,
            return_twice_shear: twice_shear,
            modes,
            return_scale,
            return_ratios,
        })
    }

    /// Select energetic linear Prager hardening before retaining any history.
    /// The backstress is `X=C*eps_p` in tensor coordinates and its storage is
    /// `C*eps_p:eps_p/2`. C is finite and nonnegative [Pa]; it is NOT the J2
    /// uniaxial hardening slope (that convention uses C=2*H_k/3). No kinematic
    /// saturation, dynamic recovery, pre-existing backstress or history rebinding
    /// is supplied.
    /// The original anisotropic metric and scalar return solver are retained.
    pub fn with_linear_kinematic_hardening(mut self, modulus_pa: f64) -> Result<Self, HillError> {
        if !modulus_pa.is_finite() || modulus_pa < 0.0 {
            return Err(HillError::Parameters(
                "kinematic modulus C must be finite and nonnegative",
            ));
        }
        let effective = self.twice_shear + modulus_pa;
        let returns = self
            .modes
            .eigenvalues
            .map(|value| self.hardening + effective * value);
        if !positive(effective)
            || returns.iter().any(|value| !positive(*value))
            || self
                .modes
                .eigenvalues
                .iter()
                .any(|value| !positive(effective * value))
        {
            return Err(HillError::Parameters(
                "kinematic return moduli are not representable",
            ));
        }
        let scale = returns.iter().copied().fold(0.0_f64, f64::max)
            + self.saturation_pa * self.saturation_rate_per_strain;
        let ratios = returns.map(|value| value / scale);
        if !positive(scale) || ratios.iter().any(|value| !positive(*value)) {
            return Err(HillError::Parameters(
                "kinematic return modulus ratios are not representable",
            ));
        }
        self.kinematic = modulus_pa;
        self.return_twice_shear = effective;
        self.return_scale = scale;
        self.return_ratios = ratios;
        Ok(self)
    }

    /// Select energetic Voce isotropic hardening before retaining any history.
    /// Q >= 0 [Pa] adds `Q*(1-exp(-b*alpha))` to the existing yield radius;
    /// b > 0 is per Hill work-conjugate strain, not generally J2 equivalent
    /// strain. Linear isotropic and Prager terms remain independently active.
    /// Q=0 preserves the original return, tangent and energy arithmetic.
    /// Explicit coefficients do not establish a measured-metal calibration.
    pub fn with_voce_isotropic_hardening(
        mut self,
        saturation_pa: f64,
        rate_per_strain: f64,
    ) -> Result<Self, HillError> {
        let slope = saturation_pa * rate_per_strain;
        if !saturation_pa.is_finite()
            || saturation_pa < 0.0
            || !positive(rate_per_strain)
            || !slope.is_finite()
            || (saturation_pa > 0.0 && slope == 0.0)
            || !(self.yield_stress + saturation_pa).is_finite()
        {
            return Err(HillError::Parameters(
                "Voce Q must be nonnegative and b positive with representable strength and slope",
            ));
        }
        let returns = self
            .modes
            .eigenvalues
            .map(|value| self.hardening + self.return_twice_shear * value);
        // Include the maximum saturation slope in the prepared scale. Both
        // the initial and endpoint normalized return derivatives stay <= 1.
        let scale = returns.iter().copied().fold(0.0_f64, f64::max) + slope;
        let ratios = returns.map(|value| value / scale);
        if !positive(scale) || ratios.iter().any(|value| !positive(*value)) {
            return Err(HillError::Parameters(
                "Voce return modulus ratios are not representable",
            ));
        }
        self.saturation_pa = saturation_pa;
        self.saturation_rate_per_strain = rate_per_strain;
        self.return_scale = scale;
        self.return_ratios = ratios;
        Ok(self)
    }

    /// Original fixed isotropic elastic coefficients and strain domain.
    pub const fn elastic(&self) -> &IsotropicElastic {
        &self.elastic
    }
    /// Reference Hill yield radius, not an inferred directional tensile yield.
    pub const fn yield_stress_pa(&self) -> f64 {
        self.yield_stress
    }
    /// Linear energetic modulus against the Hill work-conjugate accumulation.
    pub const fn hardening_modulus_pa(&self) -> f64 {
        self.hardening
    }
    /// C [Pa] in `X=C*eps_p`; zero retains the original isotropic-hardening law.
    pub const fn kinematic_modulus_pa(&self) -> f64 {
        self.kinematic
    }
    /// Additional isotropic saturation strength Q [Pa]; zero disables Voce.
    pub const fn voce_saturation_pa(&self) -> f64 {
        self.saturation_pa
    }
    /// Approach rate b per Hill work-conjugate strain; unused when Q=0.
    pub const fn voce_rate_per_strain(&self) -> f64 {
        self.saturation_rate_per_strain
    }
    /// Energetic isotropic force H*alpha + Q*(1-exp(-b*alpha)) [Pa].
    pub fn isotropic_hardening_pa(&self, alpha: f64) -> f64 {
        let linear = self.hardening * alpha;
        if self.saturation_pa == 0.0 {
            linear
        } else {
            linear + hardening::strength(self.saturation_pa, self.saturation_rate_per_strain, alpha)
        }
    }
    /// Endpoint isotropic slope H + Q*b*exp(-b*alpha) [Pa].
    /// The independent kinematic modulus does not expand the yield radius.
    pub fn isotropic_hardening_derivative_pa(&self, alpha: f64) -> f64 {
        if self.saturation_pa == 0.0 {
            self.hardening
        } else {
            self.hardening
                + hardening::slope(self.saturation_pa, self.saturation_rate_per_strain, alpha)
        }
    }
    /// Original dimensionless quadratic coefficients in the material axes.
    pub const fn coefficients(&self) -> &Hill48Coefficients {
        &self.coefficients
    }

    /// Virgin material history. It has no inferred prestrain or annealing state.
    pub fn initial_state(&self) -> HillState {
        HillState::default()
    }

    /// Equivalent stress of the supplied tensor in fixed material axes.
    /// With kinematic hardening, yielding uses the relative tensor sigma-X;
    /// this state-free observation does not subtract a backstress implicitly.
    pub fn equivalent_stress_pa(&self, stress: &Voigt) -> Result<f64, HillError> {
        if stress.iter().any(|v| !v.is_finite()) {
            return Err(HillError::State("stress must be finite"));
        }
        // Subtract nearby normal stresses before rescaling: these differences
        // can be exact even when their common hydrostatic part is enormous.
        // Normalizing first would discard a substantial part of that signal.
        let direct = self.modes.equivalent(&self.modes.coordinates(stress));
        if direct.is_finite() {
            return Ok(direct);
        }
        // A finite effective stress can still have overflowing coordinate or
        // Mandel intermediates. Normalize only when that direct path fails.
        let scale = stress.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
        if scale == 0.0 {
            return Ok(0.0);
        }
        let normalized = stress.map(|value| value / scale);
        let value = scale * self.modes.equivalent(&self.modes.coordinates(&normalized));
        if !value.is_finite() {
            return Err(HillError::Numerical(
                "equivalent stress is not representable",
            ));
        }
        Ok(value)
    }

    fn radius(&self, alpha: f64) -> Result<f64, HillError> {
        let radius = self.yield_stress + self.isotropic_hardening_pa(alpha);
        if !positive(radius) {
            return Err(HillError::Numerical("yield radius is not representable"));
        }
        Ok(radius)
    }

    fn admit(&self, strain: &Voigt, state: &HillState) -> Result<Voigt, HillError> {
        let limit = self.elastic.strain_limit.min(0.05);
        if strain.iter().any(|v| !v.is_finite())
            || frobenius_norm(strain) > limit
            || state.plastic_strain.iter().any(|v| !v.is_finite())
            || !(state.alpha.is_finite() && state.alpha >= 0.0)
            || (state.alpha == 0.0 && state.plastic_strain.iter().any(|v| *v != 0.0))
        {
            return Err(HillError::State(
                "strain or history exceeds the finite small-strain domain",
            ));
        }
        let magnitude = frobenius_norm(&state.plastic_strain);
        if trace(&state.plastic_strain).abs() > 128.0 * f64::EPSILON * magnitude {
            return Err(HillError::State(
                "Hill plastic history must preserve volume",
            ));
        }
        // Associated increments have unit A^-1 norm. The triangle inequality
        // therefore bounds the retained offset by the accumulated multiplier,
        // including reversals and closed plastic paths.
        let modes = self.modes.coordinates(&state.plastic_strain);
        let required = norm(&core::array::from_fn::<_, 5, _>(|i| {
            modes[i] / self.modes.roots[i]
        }));
        if !required.is_finite()
            || required - state.alpha > 256.0 * f64::EPSILON * required.max(state.alpha)
        {
            return Err(HillError::State(
                "plastic offset exceeds accumulated Hill work-conjugate strain",
            ));
        }
        self.radius(state.alpha)?;
        let elastic = core::array::from_fn(|i| strain[i] - state.plastic_strain[i]);
        if elastic.iter().any(|v| !v.is_finite()) || frobenius_norm(&elastic) > limit {
            return Err(HillError::State(
                "elastic predictor exceeds the declared small-strain bound",
            ));
        }
        Ok(elastic)
    }

    fn stress_at_fixed_history(&self, elastic: &Voigt) -> Voigt {
        let modes = self
            .modes
            .coordinates(elastic)
            .map(|v| self.twice_shear * v);
        let mut stress = self.modes.tensor(&modes);
        let pressure = self.bulk * trace(elastic);
        for component in &mut stress[..3] {
            *component += pressure;
        }
        stress
    }

    fn relative_modes(&self, elastic: &Voigt, state: &HillState) -> [f64; 5] {
        let mut modes = self
            .modes
            .coordinates(elastic)
            .map(|value| self.twice_shear * value);
        if self.kinematic > 0.0 {
            let plastic = self.modes.coordinates(&state.plastic_strain);
            for i in 0..5 {
                modes[i] -= self.kinematic * plastic[i];
            }
        }
        modes
    }

    /// Check one accepted state without projecting it or advancing its history.
    pub fn validate_committed_history(
        &self,
        strain: &Voigt,
        state: &HillState,
        yield_tolerance_pa: f64,
    ) -> Result<(), HillError> {
        if !positive(yield_tolerance_pa) {
            return Err(HillError::Parameters(
                "yield allowance must be finite and positive",
            ));
        }
        let elastic = self.admit(strain, state)?;
        let equivalent = if self.kinematic == 0.0 {
            self.equivalent_stress_pa(&self.stress_at_fixed_history(&elastic))?
        } else {
            self.modes.equivalent(&self.relative_modes(&elastic, state))
        };
        if !equivalent.is_finite() {
            return Err(HillError::Numerical("relative stress is not representable"));
        }
        let residual = equivalent - self.radius(state.alpha)?;
        if residual > yield_tolerance_pa {
            return Err(HillError::YieldResidual {
                residual_pa: residual,
                tolerance_pa: yield_tolerance_pa,
            });
        }
        Ok(())
    }

    fn returned(&self, strain: &Voigt, state: &HillState) -> Result<HillTrial, HillError> {
        self.returned_with_rate(strain, state, None)
    }

    fn returned_with_rate(
        &self,
        strain: &Voigt,
        state: &HillState,
        rate: Option<&LinearRateReturn>,
    ) -> Result<HillTrial, HillError> {
        let (return_scale, return_ratios, incremental_hardening) = match rate {
            Some(rate) => (rate.scale, rate.ratios, self.hardening + rate.modulus_pa),
            None => (self.return_scale, self.return_ratios, self.hardening),
        };
        let elastic = self.admit(strain, state)?;
        let trial = self.relative_modes(&elastic, state);
        let q_trial = self.modes.equivalent(&trial);
        let y0 = self.radius(state.alpha)?;
        let pressure = self.bulk * trace(&elastic);
        if !q_trial.is_finite() || !pressure.is_finite() || trial.iter().any(|v| !v.is_finite()) {
            return Err(HillError::Numerical("trial stress is not representable"));
        }
        let largest_root = self.modes.roots.iter().copied().fold(0.0_f64, f64::max);
        let subtraction_scale = if self.kinematic == 0.0 {
            (self.twice_shear * largest_root)
                * (frobenius_norm(strain) + frobenius_norm(&state.plastic_strain))
        } else {
            (self.twice_shear * largest_root) * frobenius_norm(strain)
                + (self.return_twice_shear * largest_root) * frobenius_norm(&state.plastic_strain)
        };
        let yield_roundoff = 64.0 * f64::EPSILON * q_trial.max(y0).max(subtraction_scale);
        if q_trial - y0 <= yield_roundoff {
            let stress = self.stress_at_fixed_history(&elastic);
            if stress.iter().any(|v| !v.is_finite()) {
                return Err(HillError::Numerical("elastic stress is not representable"));
            }
            return Ok(HillTrial {
                stress_pa: stress,
                state: state.clone(),
                plastic_increment: 0.0,
                yield_residual_pa: q_trial - y0,
                tangent: None,
            });
        }

        // In mode i, xi_i = y*xi_trial_i / (y + (2mu+C)*lambda_i*dl),
        // where y=y0+(H+eta/dt)*dl+Delta_R_Voce for selected linear kinetics.
        // Its denominator derivative H+eta/dt+(2mu+C)*lambda_i+R'_Voce
        // is always strictly positive. eta/dt is zero in the native law.
        // Thus sum_i [sqrt(lambda_i)*xi_trial_i/denominator_i]^2 = 1.
        // Normalize stresses and moduli independently before solving for
        // z = return_scale*dl/q_trial. This retains stress-unit covariance.
        let unit: [f64; 5] = core::array::from_fn(|i| self.modes.roots[i] * trial[i] / q_trial);
        let radius_ratio = y0 / q_trial;
        if !positive(radius_ratio) {
            return Err(HillError::Numerical(
                "yield-to-trial ratio is not representable",
            ));
        }
        // return_scale bounds every denominator derivative from above;
        // the minimum linear modal ratio bounds it from below, also for Voce.
        let mut lower = (q_trial - y0) / q_trial;
        let minimum = return_ratios.iter().copied().fold(f64::INFINITY, f64::min);
        let mut upper = lower / minimum;
        if !positive(upper) {
            return Err(HillError::Numerical(
                "scalar return bracket is not representable",
            ));
        }
        let mut z = lower;
        let mut normalized = [0.0; 5];
        let mut denominators = [0.0; 5];
        let mut derivative_ratios = return_ratios;
        let (remaining, initial_slope) = if self.saturation_pa == 0.0 {
            (0.0, 0.0)
        } else {
            let decay = det::exp(-self.saturation_rate_per_strain * state.alpha);
            (
                self.saturation_pa * decay,
                (self.saturation_pa * self.saturation_rate_per_strain) * decay / return_scale,
            )
        };
        let mut saturation_increment_ratio = 0.0;
        let mut solved = false;
        let mut residual = f64::INFINITY;
        for _ in 0..MAX_RETURN_ITERATIONS {
            if self.saturation_pa > 0.0 {
                let increment = (q_trial / return_scale) * z;
                saturation_increment_ratio =
                    hardening::strength(remaining, self.saturation_rate_per_strain, increment)
                        / q_trial;
                let slope_ratio =
                    initial_slope * det::exp(-self.saturation_rate_per_strain * increment);
                derivative_ratios = return_ratios.map(|ratio| ratio + slope_ratio);
            }
            for i in 0..5 {
                denominators[i] = radius_ratio + return_ratios[i] * z + saturation_increment_ratio;
                normalized[i] = unit[i] / denominators[i];
            }
            let value = norm(&normalized);
            residual = value - 1.0;
            if !value.is_finite() {
                return Err(HillError::Numerical(
                    "scalar return norm is not representable",
                ));
            }
            if residual.abs() <= 32.0 * f64::EPSILON {
                solved = true;
                break;
            }
            if residual > 0.0 {
                lower = z;
            } else {
                upper = z;
            }
            let derivative = -value
                * (0..5)
                    .map(|i| {
                        (normalized[i] / value).powi(2) * (derivative_ratios[i] / denominators[i])
                    })
                    .sum::<f64>();
            let newton = z - residual / derivative;
            z = if newton.is_finite() && newton > lower && newton < upper {
                newton
            } else {
                0.5 * lower + 0.5 * upper
            };
        }
        if !solved {
            return Err(HillError::NonConvergence {
                iterations: MAX_RETURN_ITERATIONS,
                residual,
            });
        }
        let increment = (q_trial / return_scale) * z;
        let alpha = state.alpha + increment;
        if !positive(increment) || !alpha.is_finite() || alpha == state.alpha {
            return Err(HillError::Numerical(
                "plastic increment or accumulation is not representable",
            ));
        }
        let yield_ratio =
            radius_ratio + (incremental_hardening / return_scale) * z + saturation_increment_ratio;
        let ratios: [f64; 5] = core::array::from_fn(|i| yield_ratio / denominators[i]);
        let returned: [f64; 5] = core::array::from_fn(|i| trial[i] * ratios[i]);
        let offset = self.modes.tensor(&core::array::from_fn(|i| {
            increment * self.modes.roots[i] * normalized[i]
        }));
        let candidate = HillState {
            plastic_strain: core::array::from_fn(|i| state.plastic_strain[i] + offset[i]),
            alpha,
        };
        self.admit(strain, &candidate)?;
        let mut stress = self.modes.tensor(&returned);
        if self.kinematic > 0.0 {
            for (value, plastic) in stress.iter_mut().zip(candidate.plastic_strain) {
                *value += self.kinematic * plastic;
            }
        }
        for component in &mut stress[..3] {
            *component += pressure;
        }

        let yield_residual = self.modes.equivalent(&returned) - self.radius(alpha)?;
        if stress.iter().any(|v| !v.is_finite()) || !yield_residual.is_finite() {
            return Err(HillError::Numerical(
                "returned stress or yield residual is not representable",
            ));
        }
        Ok(HillTrial {
            stress_pa: stress,
            state: candidate,
            plastic_increment: increment,
            yield_residual_pa: yield_residual,
            tangent: Some(PlasticTangent {
                return_scale,
                // Differentiation uses y-dl*y', not the old yield radius or
                // a secant slope. The concave Voce increment adds a positive
                // remainder; its shared series avoids cancellation near yield.
                radius_ratio: radius_ratio
                    + if self.saturation_pa == 0.0 {
                        0.0
                    } else {
                        hardening::slope_remainder(
                            remaining,
                            self.saturation_rate_per_strain,
                            increment,
                        ) / q_trial
                    },
                denominators,
                normalized,
                ratios,
                derivative_ratios,
            }),
        })
    }

    fn returned_tangent(&self, trial: &HillTrial) -> Result<Tangent6, HillError> {
        let Some(data) = &trial.tangent else {
            return Ok(bulk_shear_stiffness(self.bulk, 0.5 * self.twice_shear));
        };
        let PlasticTangent {
            return_scale,
            radius_ratio,
            denominators,
            normalized,
            ratios,
            derivative_ratios,
        } = data;

        // Analytic rank-one correction in the prepared orthonormal modes.
        // Dividing by return_scale/q_trial before forming products avoids
        // squared stress moduli and preserves the symmetric Mandel tangent.
        let b = (0..5)
            .map(|i| normalized[i].powi(2) * derivative_ratios[i] / denominators[i])
            .sum::<f64>();
        if !positive(b) {
            return Err(HillError::Numerical(
                "return derivative is not representable",
            ));
        }
        let w: [f64; 5] = core::array::from_fn(|i| {
            normalized[i]
                * det::sqrt(
                    ((self.twice_shear * self.modes.eigenvalues[i] / *return_scale)
                        / denominators[i]
                        / b)
                        * (*radius_ratio / denominators[i]),
                )
        });
        let direction: Voigt =
            core::array::from_fn(|i| (0..5).map(|k| w[k] * self.modes.vectors[k][i]).sum());
        // xi=2mu*eps_dev-(2mu+C)*eps_p. Reconstructing physical stress
        // adds C*eps_p; its derivative modifies only the modal diagonal.
        // The rank term already retains the physical driving modulus 2mu.
        let diagonal = if self.kinematic == 0.0 {
            *ratios
        } else {
            ratios.map(|ratio| {
                self.kinematic / self.return_twice_shear
                    + (self.twice_shear / self.return_twice_shear) * ratio
            })
        };
        let mandel: Tangent6 = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                let dev = (0..5)
                    .map(|k| diagonal[k] * self.modes.vectors[k][i] * self.modes.vectors[k][j])
                    .sum::<f64>()
                    - direction[i] * direction[j];
                (if i < 3 && j < 3 { self.bulk } else { 0.0 }) + self.twice_shear * dev
            })
        });
        let tangent = mandel_stiffness_to_voigt(&mandel);
        if tangent.iter().flatten().any(|v| !v.is_finite()) {
            return Err(HillError::Numerical(
                "returned tangent is not representable",
            ));
        }
        Ok(tangent)
    }

    /// Checked backward-Euler stress, algorithmic tangent and candidate history.
    /// The elastic-side tangent is selected at first yield and within its
    /// floating-point stress envelope. The envelope includes cancellation in
    /// subtracting retained plastic strain, so inspecting an accepted endpoint
    /// cannot invent a new plastic step from roundoff. This is not a physical
    /// transition band. The scalar return has a fixed iteration cap.
    pub fn response(&self, strain: &Voigt, state: &HillState) -> Result<HillResponse, HillError> {
        self.response_with_rate(strain, state, None)
    }

    fn response_with_rate(
        &self,
        strain: &Voigt,
        state: &HillState,
        rate: Option<&LinearRateReturn>,
    ) -> Result<HillResponse, HillError> {
        let trial = self.returned_with_rate(strain, state, rate)?;
        let tangent = self.returned_tangent(&trial)?;
        Ok(HillResponse {
            stress_pa: trial.stress_pa,
            tangent_pa: tangent,
            state: trial.state,
            plastic_increment: trial.plastic_increment,
            yield_residual_pa: trial.yield_residual_pa,
        })
    }

    /// Recoverable elastic and explicitly ENERGETIC isotropic/kinematic storage.
    /// The Prager contribution is C*eps_p:eps_p/2; it is never plastic heat.
    /// This thermodynamic interpretation does not calibrate a real heat fraction.
    pub fn stored_energy(
        &self,
        strain: &Voigt,
        state: &HillState,
    ) -> Result<StoredEnergy, HillError> {
        let elastic = self.admit(strain, state)?;
        let modes = self.modes.coordinates(&elastic);
        let mut isotropic_storage = 0.5 * self.hardening * state.alpha * state.alpha;
        if self.saturation_pa > 0.0 {
            isotropic_storage += hardening::storage(
                self.saturation_pa,
                self.saturation_rate_per_strain,
                state.alpha,
            );
        }
        let value = StoredEnergy {
            elastic_j_m3: 0.5 * self.bulk * trace(&elastic).powi(2)
                + 0.5 * self.twice_shear * modes.iter().map(|v| v * v).sum::<f64>(),
            hardening_j_m3: isotropic_storage,
            kinematic_hardening_j_m3: if self.kinematic == 0.0 {
                0.0
            } else {
                0.5 * self.kinematic * contract(&state.plastic_strain, &state.plastic_strain)
            },
        };
        if !value.total_j_m3().is_finite() {
            return Err(HillError::Numerical("stored energy is not representable"));
        }
        Ok(value)
    }

    /// Audit endpoint work as storage change + `yield_stress*delta_alpha`
    /// physical dissipation + backward-Euler numerical loss. Hardening storage
    /// and numerical loss are never silently converted into physical heat.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &HillState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<HillEnergyIncrement, HillError> {
        self.advance_with_energy_from(self, old_strain, old_state, strain, tolerances)
    }

    /// Native endpoint-frozen work account for a checked coefficient change.
    /// The thermal owner adds the frozen-history coefficient storage change
    /// and mechanical entropy. A softened endpoint must not reject an old
    /// state that was admissible at its actual previous temperature.
    pub(crate) fn advance_with_energy_from(
        &self,
        previous: &Self,
        old_strain: &Voigt,
        old_state: &HillState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<HillEnergyIncrement, HillError> {
        self.advance_with_rate_energy_from(
            previous, old_strain, old_state, strain, tolerances, None,
        )
    }

    fn advance_with_rate_energy_from(
        &self,
        previous: &Self,
        old_strain: &Voigt,
        old_state: &HillState,
        strain: &Voigt,
        tolerances: EnergyTolerances,
        rate: Option<&LinearRateReturn>,
    ) -> Result<HillEnergyIncrement, HillError> {
        if !positive(tolerances.energy_j_m3) || !positive(tolerances.yield_pa) {
            return Err(HillError::Parameters(
                "energy and yield allowances must be finite and positive",
            ));
        }
        if self.coefficients != previous.coefficients {
            return Err(HillError::Parameters(
                "thermal Hill history requires the same fixed anisotropic metric",
            ));
        }
        if rate.is_some() {
            // A valid viscous endpoint can lie outside the static yield surface.
            // Check its actual old coefficient/history domain without projecting
            // it or inventing another physical interval.
            previous.admit(old_strain, old_state)?;
        } else {
            previous.validate_committed_history(old_strain, old_state, tolerances.yield_pa)?;
        }
        let response = self.response_with_rate(strain, old_state, rate)?;
        let overstress = rate.map_or(0.0, |value| value.modulus_pa * response.plastic_increment);
        let consistency_residual = response.yield_residual_pa - overstress;
        if !overstress.is_finite() || !consistency_residual.is_finite() {
            return Err(HillError::Numerical(
                "viscous consistency is not representable",
            ));
        }
        if consistency_residual > tolerances.yield_pa
            || (response.plastic_increment > 0.0
                && consistency_residual.abs() > tolerances.yield_pa)
        {
            return Err(HillError::YieldResidual {
                residual_pa: consistency_residual,
                tolerance_pa: tolerances.yield_pa,
            });
        }
        let old_elastic = self.admit(old_strain, old_state)?;
        let new_elastic = self.admit(strain, &response.state)?;
        let reconstruction = self.stress_at_fixed_history(&new_elastic);
        if frobenius_norm(&core::array::from_fn(|i| {
            reconstruction[i] - response.stress_pa[i]
        })) > tolerances.yield_pa
        {
            return Err(HillError::State(
                "returned stress disagrees with committed elastic history",
            ));
        }
        let de: Voigt = core::array::from_fn(|i| new_elastic[i] - old_elastic[i]);
        let total: Voigt = core::array::from_fn(|i| strain[i] - old_strain[i]);
        let modes = self.modes.coordinates(&de);
        let da = response.plastic_increment;
        let kinematic_loss = if self.kinematic == 0.0 {
            0.0
        } else {
            let dp = core::array::from_fn(|i| {
                response.state.plastic_strain[i] - old_state.plastic_strain[i]
            });
            0.5 * self.kinematic * contract(&dp, &dp)
        };
        let mut energy = IncrementEnergy {
            initial: self.stored_energy(old_strain, old_state)?,
            final_state: self.stored_energy(strain, &response.state)?,
            mechanical_work_j_m3: contract(&response.stress_pa, &total),
            physical_dissipation_j_m3: self.yield_stress * da + overstress * da,
            numerical_dissipation_j_m3: 0.5 * self.bulk * trace(&de).powi(2)
                + 0.5 * self.twice_shear * modes.iter().map(|v| v * v).sum::<f64>()
                + 0.5 * self.hardening * da * da
                + kinematic_loss,
            residual_j_m3: 0.0,
        };
        if self.saturation_pa > 0.0 {
            energy.numerical_dissipation_j_m3 += hardening::numerical_loss(
                self.saturation_pa,
                self.saturation_rate_per_strain,
                old_state.alpha,
                da,
            );
        }
        energy.residual_j_m3 = energy.mechanical_work_j_m3
            - energy.stored_change_j_m3()
            - energy.physical_dissipation_j_m3
            - energy.numerical_dissipation_j_m3;
        if !energy.residual_j_m3.is_finite() || energy.residual_j_m3.abs() > tolerances.energy_j_m3
        {
            return Err(HillError::EnergyResidual {
                residual_j_m3: energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            });
        }
        Ok(HillEnergyIncrement {
            strain: *strain,
            response,
            energy,
        })
    }
}

impl SmallStrainLaw for Hill48 {
    type State = HillState;
    fn initial_state(&self) -> HillState {
        Self::initial_state(self)
    }
    fn stress(&self, strain: &Voigt, state: &HillState) -> Voigt {
        self.returned(strain, state)
            .map_or([f64::NAN; 6], |r| r.stress_pa)
    }
    fn tangent(&self, strain: &Voigt, state: &HillState) -> Tangent6 {
        self.response(strain, state)
            .map_or([[f64::NAN; 6]; 6], |r| r.tangent_pa)
    }
    fn update_state(&self, strain: &Voigt, state: &HillState) -> HillState {
        self.returned(strain, state)
            .map_or_else(|_| state.clone(), |r| r.state)
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "Hill48 quadratic orthotropic yield; no fracture or forming-limit prediction",
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.hill48-plasticity",
            "0.1.0",
            Ambition::Solid,
            vec![
                "small strain and rotation; fixed material axes".into(),
                "associative orthotropic pressure-independent yielding".into(),
                if self.saturation_pa == 0.0 {
                    "isotropic elasticity and explicitly selected energetic linear isotropic/Prager hardening"
                } else {
                    "isotropic elasticity and explicitly selected energetic linear/Voce isotropic and Prager hardening"
                }
                .into(),
            ],
            ValidityDomain::unconstrained().with(
                "strain-magnitude",
                0.0,
                self.elastic.strain_limit.min(0.05),
            ),
            vec![
                if self.saturation_pa == 0.0 {
                    "no evolving texture, nonlinear cyclic hardening, rate or thermal dependence"
                } else {
                    "no evolving texture, nonlinear backstress recovery, rate or thermal dependence"
                }
                .into(),
                "no calibrated material accuracy, damage, necking or finite-strain forming".into(),
            ],
            0.05,
        )
    }
}
