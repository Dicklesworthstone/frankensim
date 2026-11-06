//! Hyperelastic energy potentials (plan §8.2): compressible Neo-Hookean
//! and Mooney–Rivlin, written ONCE as generic-scalar stored energies.
//! Stress is the exact gradient (fs-ad duals), and the consistent tangent
//! is the exact Hessian (nested duals) — derivative consistency is by
//! construction, then independently FD-gated in conformance. Ogden uses
//! coalescence-safe spectral derivatives; HGO uses tension-only structural
//! strains. Gent retains finite chain extensibility. Every law enforces its
//! declared principal-stretch domain and its own constitutive limits.

use crate::{MaterialAdmissibility, MaterialError};
use fs_ad::{Dual, Dual64, Real};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

/// Finite-deformation acoustic tensors and signed directional wave modes.
pub mod acoustic;
/// Dimensioned constitutive-model cards and graph registration.
pub mod cards;
/// Objective two-state elastic work and consistent discrete-gradient tangents.
pub mod discrete_gradient;
/// Exact incompressibility with traction-free thickness and transverse shears.
pub mod incompressible_plane_stress;
mod increment;
/// Original decoupled energies in a mixed displacement/pressure potential.
pub mod mixed;
/// Traction-free finite-strain thickness and transverse-shear relaxation.
pub mod plane_stress;
/// Prescribed stress-free reference deformation and original-volume pullbacks.
pub mod reference;
mod spectral;
/// Exact natural-strain thermodynamics and finite isentropic acoustic tangents.
pub mod thermoelastic;

/// Maximum number of independently supplied Ogden terms.
pub const MAX_OGDEN_TERMS: usize = 8;
/// Maximum number of independently supplied fibre families.
pub const MAX_FIBRE_FAMILIES: usize = 4;

/// An Ogden coefficient in the original `mu/alpha` convention (Pa).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OgdenTerm {
    /// Stress coefficient, which may be negative when alpha is negative.
    pub mu_pa: f64,
    /// Dimensionless exponent; `mu_pa * alpha` must be positive.
    pub alpha: f64,
}

/// Decoupled volumetric stored energy per reference volume (Pa).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Volumetric {
    /// `U(J) = kappa/2 (J-1)^2`.
    Quadratic { kappa: f64 },
    /// `U(J) = kappa/2 (ln J)^2`.
    LogSquared { kappa: f64 },
    /// `U(J) = kappa/4 (J^2 - 1 - 2 ln J)` (Simo--Taylor).
    /// Its curvature is positive for every `J > 0`, and compression energy
    /// diverges as `J` tends to zero. The pressure inverse is single-valued.
    /// See Hartmann & Neff (2003), doi:10.1016/S0020-7683(03)00086-6,
    /// Table 4. This is an explicit material choice, not a replacement for
    /// either existing bulk law or an experimental calibration.
    SimoTaylor { kappa: f64 },
}

impl Volumetric {
    /// Infinitesimal bulk modulus (Pa).
    #[must_use]
    pub fn kappa(self) -> f64 {
        match self {
            Self::Quadratic { kappa } | Self::LogSquared { kappa } | Self::SimoTaylor { kappa } => {
                kappa
            }
        }
    }

    /// Stored energy, with `J > 0` required by the caller.
    pub fn energy<T: Real>(self, j: T) -> T {
        if let Self::SimoTaylor { kappa } = self {
            let increment = j - T::one();
            return if increment.value().abs() <= 0.125 {
                simo_taylor_small_energy(kappa, increment)
            } else {
                T::from_f64(kappa)
                    * (T::from_f64(0.25) * increment * (j + T::one()) - T::from_f64(0.5) * j.ln())
            };
        }
        let strain = match self {
            Self::Quadratic { .. } => j - T::one(),
            Self::LogSquared { .. } => j.ln(),
            Self::SimoTaylor { .. } => unreachable!(),
        };
        T::from_f64(0.5 * self.kappa()) * strain * strain
    }
}

// U/kappa = x^2/2 [1 - x/3 + x^2/4 - ...], x=J-1.
// Retain weak compression without subtracting O(x) terms. On |x|<=1/8
// the omitted second derivative is below roundoff; nested duals differentiate
// this same expression, including the exact reference stress and bulk modulus.
fn simo_taylor_small_energy<T: Real>(kappa: f64, x: T) -> T {
    let mut tail = T::from_f64(1.0 / 22.0);
    for n in (1..20).rev() {
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        tail = tail * x + T::from_f64(sign / (n + 2) as f64);
    }
    T::from_f64(0.5 * kappa) * x * x * (T::one() + x * tail)
}

/// Structural invariants used by the fibre contribution. The matrix always
/// uses the isochoric invariant. Full invariants avoid suppressing fibres
/// under hydrostatic expansion (Nolan et al., 2014, doi:10.1016/j.jmbbm.2014.06.016).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FibreInvariants {
    /// Use `Cbar = J^(-2/3) F^T F` (HGO 2000).
    Isochoric,
    /// Use `C = F^T F` for the structural strain (compressible extension).
    Full,
}

/// The implemented energy potentials.
#[derive(Debug, Clone, PartialEq)]
pub enum HyperelasticModel {
    /// Compressible Neo-Hookean:
    /// `W = μ/2 (I₁ − 3) − μ ln J + λ/2 (ln J)²`.
    NeoHookean {
        /// Shear modulus μ.
        mu: f64,
        /// Lamé λ. It may be negative for an auxetic reference response;
        /// the infinitesimal bulk modulus `λ + 2μ/3` must remain positive.
        lambda: f64,
    },
    /// Mooney–Rivlin (decoupled):
    /// `W = c10 (Ī₁ − 3) + c01 (Ī₂ − 3) + κ/2 (J − 1)²`.
    MooneyRivlin {
        /// Ī₁ coefficient.
        c10: f64,
        /// Ī₂ coefficient.
        c01: f64,
        /// Volumetric stiffness κ.
        kappa: f64,
    },
    /// Gent finite-chain elasticity with an explicitly selected bulk energy:
    /// `W = -mu*Jm/2 * ln(1 - (I1bar-3)/Jm) + U(J)`.
    ///
    /// `Jm` bounds the isochoric first-invariant increment, not a principal
    /// stretch or a rupture strain. The strict domain is `I1bar-3 < Jm`.
    /// Gent (1996), doi:10.5254/1.3538357; the incompressible convention is
    /// reproduced in Horgan & Murphy (2022), doi:10.1007/s10659-022-09908-1,
    /// Eq. (6.11). Replacing I1 by I1bar and adding U is the declared
    /// compressible extension, not a bulk law inferred from that paper.
    Gent {
        /// Infinitesimal shear modulus [Pa].
        mu: f64,
        /// Positive dimensionless limiting invariant increment `Jm`.
        locking_invariant: f64,
        /// Original decoupled volumetric storage.
        volumetric: Volumetric,
    },
    /// `sum mu_p/alpha_p (sum lambda_bar_i^alpha_p - 3) + U(J)`.
    Ogden {
        terms: Vec<OgdenTerm>,
        volumetric: Volumetric,
    },
    /// HGO/GOH: `c10 (I1bar-3) + sum k1/(2 k2) (exp(k2 <E>_+^2)-1)`
    /// plus `U(J)`, where `E = dispersion*(I1star-3) +
    /// (1-3*dispersion)*(I4star-1)`. At zero dispersion the fibres are
    /// exactly inactive for `I4star <= 1`. At nonzero dispersion tension
    /// refers to the distributed structural strain E, not the mean direction.
    /// Positive `k2` approaching zero retains the quadratic fibre limit
    /// `k1 <E>_+^2 / 2`, including its stress and consistent tangent.
    HolzapfelGasserOgden {
        c10: f64,
        k1: f64,
        k2: f64,
        dispersion: f64,
        fibres: Vec<[f64; 3]>,
        invariants: FibreInvariants,
        volumetric: Volumetric,
    },
}

/// A hyperelastic law over the deformation gradient F (row-major 3×3).
#[derive(Debug, Clone, PartialEq)]
pub struct Hyperelastic {
    /// The energy model.
    pub model: HyperelasticModel,
    /// Calibrated stretch bound (validity domain on principal stretch).
    pub stretch_limit: f64,
}

fn check_admissible(f: &[f64; 9]) -> Result<(), MaterialError> {
    let j = det3(f);
    if !f.iter().all(|x| x.is_finite()) || !(j.is_finite() && j > 0.0) {
        return Err(MaterialError::State {
            what: format!("det F = {j} (orientation lost or degenerate)"),
        });
    }
    Ok(())
}

fn det3<T: Real>(f: &[T; 9]) -> T {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

/// tr(FᵀF) and tr((FᵀF)²) — the invariants I₁ and the I₂ ingredient.
fn cauchy_green_invariants<T: Real>(f: &[T; 9]) -> (T, T) {
    // C = FᵀF, c_ij = Σ_k f_ki f_kj (row-major: f[3k + i]).
    let mut c = [T::zero(); 9];
    for i in 0..3 {
        for j in 0..3 {
            let mut s = T::zero();
            for k in 0..3 {
                s = s + f[3 * k + i] * f[3 * k + j];
            }
            c[3 * i + j] = s;
        }
    }
    let i1 = c[0] + c[4] + c[8];
    let mut tr_c2 = T::zero();
    for i in 0..3 {
        for j in 0..3 {
            tr_c2 = tr_c2 + c[3 * i + j] * c[3 * j + i];
        }
    }
    (i1, tr_c2)
}

// C-I near any rigid rotation, shared by the coupled and decoupled laws.
// Retain the quadratic squares even when a column's leading square is one.
pub(crate) fn small_objective_strain<T: Real>(f: &[T; 9], j: T) -> Option<[T; 9]> {
    if !(j.value() > 0.0) {
        return None;
    }
    let mut a = [T::zero(); 9];
    for i in 0..3 {
        let pivot = (1..3).fold(0, |p, k| {
            if f[3 * k + i].value().abs() > f[3 * p + i].value().abs() {
                k
            } else {
                p
            }
        });
        let lead = f[3 * pivot + i];
        let mut diagonal = lead.mul_add(lead, -T::one());
        for k in 0..3 {
            if k != pivot {
                diagonal = f[3 * k + i].mul_add(f[3 * k + i], diagonal);
            }
        }
        a[3 * i + i] = diagonal;
        for j in 0..i {
            let value = (0..3).fold(T::zero(), |s, k| f[3 * k + i].mul_add(f[3 * k + j], s));
            a[3 * i + j] = value;
            a[3 * j + i] = value;
        }
    }
    a.iter().all(|v| v.value().abs() <= 1.0 / 64.0).then_some(a)
}

// ln(1+d) and d-ln(1+d), without cancelling their first-order terms.
// The caller bounds |d/(2+d)| < 0.073. The omitted atanh tail and its
// first two derivatives are <1e-32, below f64 roundoff in this domain.
pub(crate) fn small_log1p_and_remainder<T: Real>(d: T) -> (T, T) {
    let s = d / (T::from_f64(2.0) + d);
    let s2 = s * s;
    let mut power = s;
    let mut tail = T::zero();
    for k in 1..16 {
        power = power * s2;
        tail = tail + power / T::from_f64(f64::from(2 * k + 1));
    }
    tail = T::from_f64(2.0) * tail;
    (
        T::from_f64(2.0) * s + tail,
        d * d / (T::from_f64(2.0) + d) - tail,
    )
}

// Real has no expm1 operation. This scalar-generic series retains the same
// potential for AD; at |x| <= 1/8 its omitted term and first two derivatives
// are below 3e-26. Outside that interval keep the original exponential.
fn exp_minus_one<T: Real>(x: T) -> T {
    if x.value().abs() > 0.125 {
        return x.exp() - T::one();
    }
    x * small_exp_relative(x)
}

// The analytic continuation of (exp(x)-1)/x, including at zero. Keep the
// constant term separate so a vanishing exponent still retains finite
// quadratic fibre storage and the derivatives of that same potential.
// Only called for |x| <= 1/8, with the same series budget as exp_minus_one.
fn small_exp_relative<T: Real>(x: T) -> T {
    let mut factorial = 20_922_789_888_000.0; // 16!, exactly representable.
    let mut tail = T::from_f64(1.0 / factorial);
    for degree in (2..16).rev() {
        factorial /= f64::from(degree + 1);
        tail = tail.mul_add(x, T::from_f64(1.0 / factorial));
    }
    T::one() + x * tail
}

struct IsochoricStrain<T> {
    determinant_root: T,
    full_trace: T,
    full_quadratic: T,
    full_cubic: T,
    i1_minus_three: T,
    i2_minus_three: T,
}

// Write C=m(I+B), with m=tr(C)/3 and tr(B)=0, before normalizing det(C).
// Then det(I+B)=1+Q(B)+det(B), whose increment starts at quadratic order.
// This avoids cancellation between volumetric and distortional strains,
// including a tiny shear superposed on a much larger uniform dilation.
fn small_isochoric_strain<T: Real>(a: &[T; 9]) -> IsochoricStrain<T> {
    let three = T::from_f64(3.0);
    let full_trace = a[0] + a[4] + a[8];
    let mean = T::one() + full_trace / three;
    let mut b = a.map(|entry| entry / mean);
    let denominator = three * mean;
    b[0] = ((a[0] - a[4]) + (a[0] - a[8])) / denominator;
    b[4] = ((a[4] - a[0]) + (a[4] - a[8])) / denominator;
    b[8] = -b[0] - b[4];
    let quadratic = -T::from_f64(0.5)
        * (b[0] * b[0]
            + b[4] * b[4]
            + b[8] * b[8]
            + T::from_f64(2.0) * (b[1] * b[1] + b[2] * b[2] + b[5] * b[5]));
    let (log_det, _) = small_log1p_and_remainder(quadratic + det3(&b));
    let scale_increment = exp_minus_one(-log_det / three);
    let square_increment = exp_minus_one(-T::from_f64(2.0) * log_det / three);
    IsochoricStrain {
        determinant_root: mean / (T::one() + scale_increment),
        full_trace,
        full_quadratic: a[0] * a[4] + a[0] * a[8] + a[4] * a[8]
            - a[1] * a[3]
            - a[2] * a[6]
            - a[5] * a[7],
        full_cubic: det3(a),
        i1_minus_three: three * scale_increment,
        i2_minus_three: three * square_increment + quadratic * (T::one() + square_increment),
    }
}

// Gent's domain and energy must use the same invariant, including the
// cancellation-resistant expression near any rigid rotation.
fn gent_trace_shift<T: Real>(f: &[T; 9]) -> T {
    let j = det3(f);
    if let Some(a) = small_objective_strain(f, j) {
        small_isochoric_strain(&a).i1_minus_three
    } else {
        let i1 = f
            .iter()
            .fold(T::zero(), |sum, value| value.mul_add(*value, sum));
        (T::from_f64(-2.0 / 3.0) * j.ln()).exp() * i1 - T::from_f64(3.0)
    }
}

fn gent_energy<T: Real>(shift: T, mu: f64, locking_invariant: f64) -> T {
    let x = shift / T::from_f64(locking_invariant);
    // Factor out the invariant increment. This avoids both subtracting a
    // sub-ULP x from one and forming mu*Jm, which can overflow even when
    // energy, stress and tangent are representable in the neo-Hookean limit.
    let logarithm_ratio = if x.value().abs() <= 0.125 {
        // -ln(1-x)/x = 2/(2-x) * sum (x/(2-x))^(2k)/(2k+1).
        // At |x|<=1/8, the omitted tail and its first two derivatives are
        // below 1e-32. The expression also differentiates correctly at x=0.
        let denominator = T::from_f64(2.0) - x;
        let s = x / denominator;
        let square = s * s;
        let mut series = T::from_f64(1.0 / 31.0);
        for k in (0..15).rev() {
            series = series.mul_add(square, T::from_f64(1.0 / f64::from(2 * k + 1)));
        }
        T::from_f64(2.0) * series / denominator
    } else {
        -(T::one() - x).ln() / x
    };
    T::from_f64(0.5 * mu) * shift * logarithm_ratio
}

// Rationalize I4/det(C)^(1/3)-1 using I4^3-det(C). Group equal orders
// in C-I so axial fibres under pure dilation or unchanged by J=1 simple
// shear retain zero strain and the declared inactive-side Hessian at zero.
fn isochoric_structural_shift<T: Real>(shift: T, bar: &IsochoricStrain<T>) -> T {
    let square = shift * shift;
    let numerator = (T::from_f64(3.0) * shift - bar.full_trace)
        + (T::from_f64(3.0) * square - bar.full_quadratic)
        + (shift * square - bar.full_cubic);
    let i4 = T::one() + shift;
    let root = bar.determinant_root;
    numerator / (root * (i4 * i4 + i4 * root + root * root))
}

// a^T(Cstar-I)a + (|a|^2-1). Keep the actual supplied direction norm:
// the constructor's tolerance is an admission rule, not a normalization.
fn structural_invariant_shift<T: Real>(strain: &[T; 9], direction: &[f64; 3]) -> T {
    let norm_shift = direction.iter().fold(0.0, |sum, x| sum + x * x) - 1.0;
    let mut shift = T::from_f64(norm_shift);
    for i in 0..3 {
        for j in 0..3 {
            shift = T::from_f64(direction[i] * direction[j]).mul_add(strain[3 * i + j], shift);
        }
    }
    shift
}

fn hgo_fibre_energy<T: Real>(strain: T, k1: f64, k2: f64) -> T {
    // Absent fibres store exactly zero even if their exponent would overflow.
    if k1 != 0.0 && strain.value() > 0.0 {
        let exponent = T::from_f64(k2) * strain * strain;
        if exponent.value() <= 0.125 {
            // W = k1 E²/2 * exprel(k2 E²). Dividing by k2 first can
            // overflow or erase energy when the exponent underflows, although
            // W, its stress and its tangent have a finite quadratic limit.
            T::from_f64(0.5 * k1) * strain * strain * small_exp_relative(exponent)
        } else {
            T::from_f64(k1 / (2.0 * k2)) * exp_minus_one(exponent)
        }
    } else {
        T::zero()
    }
}

// The original coupled NH potential, evaluated in objective strain coordinates
// near ANY rigid rotation. Directly subtracting I1-3 and ln(J) loses quadratic
// storage at acoustic strain amplitudes even while the AD stress stays nonzero.
fn neo_hookean_energy<T: Real>(f: &[T; 9], mu: f64, lambda: f64) -> T {
    let j = det3(f);
    if let Some(a) = small_objective_strain(f, j) {
        let quadratic =
            a[0] * a[4] + a[0] * a[8] + a[4] * a[8] - a[1] * a[3] - a[2] * a[6] - a[5] * a[7];
        let cubic = det3(&a);
        let d = a[0] + a[4] + a[8] + quadratic + cubic; // det(C)-1.
        let (log_c, remainder) = small_log1p_and_remainder(d);
        // ln det(C)=2 ln J and tr(A)-ln det(C)=-Q(A)-det(A)
        // +[d-ln(1+d)]. AD differentiates this same stable expression.
        T::from_f64(0.5 * mu) * (-quadratic - cubic + remainder)
            + T::from_f64(0.125 * lambda) * log_c * log_c
    } else {
        let (i1, _) = cauchy_green_invariants(f);
        let ln_j = j.ln();
        T::from_f64(0.5 * mu) * (i1 - T::from_f64(3.0)) - T::from_f64(mu) * ln_j
            + T::from_f64(0.5 * lambda) * ln_j * ln_j
    }
}

impl Hyperelastic {
    /// Construct with parameter checks.
    /// Neo-Hookean requires positive reference shear and bulk moduli,
    /// `mu > 0` and `lambda + 2*mu/3 > 0`; a negative Lamé `lambda`
    /// represents an auxetic reference state. This local condition does not
    /// certify stability over the entire calibrated finite-stretch domain.
    ///
    /// # Errors
    /// [`MaterialError::Parameters`] on inadmissible moduli or a
    /// principal-stretch limit that is nonfinite or smaller than one.
    pub fn new(model: HyperelasticModel, stretch_limit: f64) -> Result<Self, MaterialError> {
        if !stretch_limit.is_finite() || stretch_limit < 1.0 {
            return Err(MaterialError::Parameters {
                what: "principal-stretch limit must be finite and at least one".into(),
            });
        }
        let ok = match &model {
            HyperelasticModel::NeoHookean { mu, lambda } => {
                let bulk = *lambda + (2.0 / 3.0) * *mu;
                mu.is_finite()
                    && lambda.is_finite()
                    && *mu > 0.0
                    && bulk.is_finite()
                    && bulk > 0.0
                    && (2.0 * *mu + *lambda).is_finite()
            }
            HyperelasticModel::MooneyRivlin { c10, c01, kappa } => {
                let shear = 2.0 * (*c10 + *c01);
                [*c10, *c01, *kappa, shear]
                    .iter()
                    .all(|value| value.is_finite())
                    && *c10 >= 0.0
                    && *c01 >= 0.0
                    && *kappa > 0.0
                    && shear > 0.0
                    && (*kappa + (4.0 / 3.0) * shear).is_finite()
            }
            HyperelasticModel::Gent {
                mu,
                locking_invariant,
                volumetric,
            } => {
                [*mu, *locking_invariant, volumetric.kappa()]
                    .iter()
                    .all(|value| value.is_finite() && *value > 0.0)
                    && (volumetric.kappa() + (4.0 / 3.0) * *mu).is_finite()
            }
            HyperelasticModel::Ogden { terms, volumetric } => {
                !terms.is_empty()
                    && terms.len() <= MAX_OGDEN_TERMS
                    && terms.iter().all(|p| {
                        p.mu_pa.is_finite()
                            && p.alpha.is_finite()
                            && (p.mu_pa / p.alpha).is_finite()
                            && (p.mu_pa * p.alpha).is_finite()
                            && p.mu_pa * p.alpha > 0.0
                    })
                    && terms
                        .iter()
                        .map(|p| 0.5 * p.mu_pa * p.alpha)
                        .sum::<f64>()
                        .is_finite()
                    && volumetric.kappa().is_finite()
                    && volumetric.kappa() > 0.0
            }
            HyperelasticModel::HolzapfelGasserOgden {
                c10,
                k1,
                k2,
                dispersion,
                fibres,
                volumetric,
                ..
            } => {
                [*c10, *k1, *k2, *dispersion, volumetric.kappa()]
                    .iter()
                    .all(|v| v.is_finite())
                    && *c10 > 0.0
                    && *k1 >= 0.0
                    && *k2 > 0.0
                    && (0.0..=1.0 / 3.0).contains(dispersion)
                    && volumetric.kappa() > 0.0
                    && fibres.len() <= MAX_FIBRE_FAMILIES
                    && fibres.iter().all(|a| {
                        a.iter().all(|v| v.is_finite())
                            && (fs_math::det::sqrt(a.iter().map(|x| x * x).sum()) - 1.0).abs()
                                <= 1e-12
                    })
            }
        };
        if !ok {
            return Err(MaterialError::Parameters {
                what: format!("{model:?}: inadmissible or unrepresentable moduli"),
            });
        }
        Ok(Hyperelastic {
            model,
            stretch_limit,
        })
    }

    /// The stored energy at a generic-scalar deformation gradient.
    ///
    /// Requires an admissible F. Ogden's spectral lift carries the exact
    /// value, gradient and Hessian at the primal point, including repeated
    /// stretches; derivatives above second order are not claimed. All native
    /// laws use equivalent objective strain expressions near rigid rotations
    /// so small stored energies survive unit-invariant subtraction.
    pub fn energy<T: Real>(&self, f: &[T; 9]) -> T {
        match &self.model {
            HyperelasticModel::NeoHookean { mu, lambda } => neo_hookean_energy(f, *mu, *lambda),
            HyperelasticModel::MooneyRivlin { .. }
            | HyperelasticModel::Gent { .. }
            | HyperelasticModel::Ogden { .. }
            | HyperelasticModel::HolzapfelGasserOgden { .. } => {
                self.isochoric_energy(f) + self.volumetric_energy(det3(f))
            }
        }
    }

    /// Original stored-energy change per reference volume, `W(after)-W(before)`.
    ///
    /// The finite-pair evaluation retains small elastic work independently of
    /// large unchanged storage. It uses the same native potential, including
    /// active fibres and original bulk terms; it is not a stress quadrature or
    /// a tangent approximation. Both endpoint domains are checked, including
    /// coincident endpoints. Invalid states or unrepresentable storage return
    /// NaN, which energy-only consumers must refuse.
    pub fn energy_increment(&self, before: &[f64; 9], after: &[f64; 9]) -> f64 {
        if self.check_state(before).is_err() || self.check_state(after).is_err() {
            return f64::NAN;
        }
        increment::energy_change(self, before, after, true)
    }

    /// Change of the same matrix/fibre remainder returned by
    /// [`Self::isochoric_energy`], without evaluating an unused bulk penalty.
    /// Full-invariant fibres retain their actual volume dependence. The
    /// endpoint-domain and nonfinite conventions match [`Self::energy_increment`].
    pub fn isochoric_energy_increment(&self, before: &[f64; 9], after: &[f64; 9]) -> f64 {
        if self.check_state(before).is_err() || self.check_state(after).is_err() {
            return f64::NAN;
        }
        increment::energy_change(self, before, after, false)
    }

    /// Distortional/matrix-plus-fibre energy for a mixed displacement/pressure
    /// consumer. For full-invariant HGO the fibre part intentionally depends
    /// on J. This accessor does not add a bulk penalty.
    pub fn isochoric_energy<T: Real>(&self, f: &[T; 9]) -> T {
        match &self.model {
            HyperelasticModel::NeoHookean { mu, .. } => neo_hookean_energy(f, *mu, 0.0),
            HyperelasticModel::MooneyRivlin { c10, c01, .. } => {
                let j = det3(f);
                if let Some(a) = small_objective_strain(f, j) {
                    let bar = small_isochoric_strain(&a);
                    return T::from_f64(*c10) * bar.i1_minus_three
                        + T::from_f64(*c01) * bar.i2_minus_three;
                }
                let ln_j = j.ln();
                let jm23 = (T::from_f64(-2.0 / 3.0) * ln_j).exp();
                let jm43 = (T::from_f64(-4.0 / 3.0) * ln_j).exp();
                let (i1, tr_c2) = cauchy_green_invariants(f);
                T::from_f64(*c10) * (jm23 * i1 - T::from_f64(3.0))
                    + T::from_f64(0.5 * c01) * (jm43 * (i1 * i1 - tr_c2) - T::from_f64(6.0))
            }
            HyperelasticModel::Gent {
                mu,
                locking_invariant,
                ..
            } => gent_energy(gent_trace_shift(f), *mu, *locking_invariant),
            HyperelasticModel::Ogden { terms, .. } => {
                let j = det3(f);
                if let Some(a) = small_objective_strain(f, j)
                    && let Some(energy) = spectral::small_ogden_energy(&a, terms)
                {
                    return energy;
                }
                let jm23 = (T::from_f64(-2.0 / 3.0) * j.ln()).exp();
                let cbar = core::array::from_fn(|ij| {
                    let (i, k) = (ij / 3, ij % 3);
                    jm23 * (0..3).fold(T::zero(), |s, r| s + f[3 * r + i] * f[3 * r + k])
                });
                spectral::ogden_energy(&cbar, terms)
            }
            HyperelasticModel::HolzapfelGasserOgden {
                c10,
                k1,
                k2,
                dispersion,
                fibres,
                invariants,
                ..
            } => {
                let j = det3(f);
                if let Some(a) = small_objective_strain(f, j) {
                    let bar = small_isochoric_strain(&a);
                    let mut w = T::from_f64(*c10) * bar.i1_minus_three;
                    if *k1 == 0.0 {
                        return w;
                    }
                    let trace = match invariants {
                        FibreInvariants::Isochoric => bar.i1_minus_three,
                        FibreInvariants::Full => bar.full_trace,
                    };
                    for direction in fibres {
                        let shift = structural_invariant_shift(&a, direction);
                        let shift = match invariants {
                            FibreInvariants::Isochoric => isochoric_structural_shift(shift, &bar),
                            FibreInvariants::Full => shift,
                        };
                        let e = T::from_f64(*dispersion) * trace
                            + T::from_f64(1.0 - 3.0 * dispersion) * shift;
                        w = w + hgo_fibre_energy(e, *k1, *k2);
                    }
                    return w;
                }
                let jm23 = (T::from_f64(-2.0 / 3.0) * j.ln()).exp();
                let (i1, _) = cauchy_green_invariants(f);
                let scale = match invariants {
                    FibreInvariants::Isochoric => jm23,
                    FibreInvariants::Full => T::one(),
                };
                let mut w = T::from_f64(*c10) * (jm23 * i1 - T::from_f64(3.0));
                if *k1 == 0.0 {
                    return w;
                }
                for a in fibres {
                    let mut i4 = T::zero();
                    for r in 0..3 {
                        let fa =
                            (0..3).fold(T::zero(), |s, c| s + f[3 * r + c] * T::from_f64(a[c]));
                        i4 = i4 + fa * fa;
                    }
                    let e = T::from_f64(*dispersion) * (scale * i1 - T::from_f64(3.0))
                        + T::from_f64(1.0 - 3.0 * dispersion) * (scale * i4 - T::one());
                    w = w + hgo_fibre_energy(e, *k1, *k2);
                }
                w
            }
        }
    }

    /// Original scalar volume term. NH retains its coupled energy split:
    /// its `lambda/2 * ln(J)^2` can be negative for an auxetic material and
    /// is not an independent bulk storage law. MR and the new variants
    /// expose their decoupled bulk penalty directly.
    pub fn volumetric_energy<T: Real>(&self, j: T) -> T {
        match &self.model {
            HyperelasticModel::NeoHookean { lambda, .. } => {
                T::from_f64(0.5 * lambda) * j.ln() * j.ln()
            }
            HyperelasticModel::MooneyRivlin { kappa, .. } => {
                Volumetric::Quadratic { kappa: *kappa }.energy(j)
            }
            HyperelasticModel::Gent { volumetric, .. }
            | HyperelasticModel::Ogden { volumetric, .. }
            | HyperelasticModel::HolzapfelGasserOgden { volumetric, .. } => volumetric.energy(j),
        }
    }

    /// Sorted singular values of F, using the shared symmetric eigensolver.
    /// No eigenvector differentiation or perturbation of repeated roots.
    pub fn principal_stretches(&self, f: &[f64; 9]) -> Result<[f64; 3], MaterialError> {
        check_admissible(f)?;
        spectral::principal_stretches(f)
    }

    /// Admit orientation, the principal-stretch domain and constitutive limits.
    ///
    /// This performs the same deformation checks as stress and tangent
    /// evaluation without computing their derivatives. Energy-only consumers
    /// can retain the original material domain during backtracking.
    ///
    /// # Errors
    /// Invalid public stretch limits return [`MaterialError::Parameters`].
    /// Nonfinite, inverted or out-of-domain deformation returns
    /// [`MaterialError::State`].
    pub fn check_state(&self, f: &[f64; 9]) -> Result<(), MaterialError> {
        // These fields are public, so a caller may have changed the limit
        // since construction. An invalid domain must not disable admission.
        if !self.stretch_limit.is_finite() || self.stretch_limit < 1.0 {
            return Err(MaterialError::Parameters {
                what: "principal-stretch limit must be finite and at least one".into(),
            });
        }
        let stretches = self.principal_stretches(f)?;
        // Eight ULPs accommodate eigensolver/rotation roundoff at the
        // inclusive boundary; 1e-12 beyond the limit still refuses.
        let roundoff = 8.0 * f64::EPSILON;
        if stretches[0] < (1.0 - roundoff) / self.stretch_limit
            || stretches[2] > self.stretch_limit * (1.0 + roundoff)
        {
            return Err(MaterialError::State {
                what: format!(
                    "principal stretches {stretches:?} outside [{}, {}]",
                    1.0 / self.stretch_limit,
                    self.stretch_limit
                ),
            });
        }
        if let HyperelasticModel::Gent {
            locking_invariant, ..
        } = &self.model
        {
            // Public model fields can be changed after construction. An
            // invalid chain limit must never disable its strict domain gate.
            if !(locking_invariant.is_finite() && *locking_invariant > 0.0) {
                return Err(MaterialError::Parameters {
                    what: "Gent locking invariant must be finite and positive".into(),
                });
            }
            let shift = gent_trace_shift(f);
            if !(shift.is_finite() && shift < *locking_invariant) {
                return Err(MaterialError::State {
                    what: format!(
                        "Gent chain limit: isochoric I1-3 = {shift} must be below {locking_invariant}"
                    ),
                });
            }
        }
        Ok(())
    }

    /// First Piola–Kirchhoff stress `P = ∂W/∂F` (row-major 9-vector).
    ///
    /// # Errors
    /// [`MaterialError::State`] when det F ≤ 0, F is nonfinite, or a
    /// principal stretch exceeds the declared domain. An invalid public
    /// `stretch_limit` returns [`MaterialError::Parameters`].
    pub fn piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        self.check_state(f)?;
        let (_, grad) = fs_ad::gradient::<9>(*f, |x| self.energy(&x));
        if grad.iter().all(|x| x.is_finite()) {
            Ok(grad)
        } else {
            Err(MaterialError::State {
                what: "nonfinite Piola stress".to_string(),
            })
        }
    }

    /// The consistent material tangent `A = ∂²W/∂F∂F` (9×9, symmetric),
    /// exact via nested duals — one inner-gradient sweep per column.
    ///
    /// # Errors
    /// [`MaterialError::State`] when det F ≤ 0, F is nonfinite, or a
    /// principal stretch exceeds the declared domain. An invalid public
    /// `stretch_limit` returns [`MaterialError::Parameters`].
    pub fn tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        self.check_state(f)?;
        let mut a = [[0.0f64; 9]; 9];
        for col in 0..9 {
            // Inner duals carry the full gradient; the outer single lane
            // differentiates it along e_col.
            let mut vars = [Dual::<Dual64<9>, 1>::constant(Dual64::constant(0.0)); 9];
            for (i, v) in vars.iter_mut().enumerate() {
                *v = Dual {
                    re: Dual64::variable(f[i], i),
                    eps: [if i == col {
                        Dual64::constant(1.0)
                    } else {
                        Dual64::constant(0.0)
                    }],
                };
            }
            let w = self.energy(&vars);
            for (row, out) in a.iter_mut().enumerate() {
                out[col] = w.eps[0].eps[row];
            }
        }
        if a.iter().flatten().all(|x| x.is_finite()) {
            Ok(a)
        } else {
            Err(MaterialError::State {
                what: "nonfinite material tangent".to_string(),
            })
        }
    }

    /// Thermodynamic declarations.
    #[must_use]
    pub fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true, // elastic: zero dissipation
            polyconvex: match &self.model {
                // Sufficient convex-in-(F,cof F,J) conditions over the
                // declared stretch domain, never a near-reference spot test.
                HyperelasticModel::NeoHookean { mu, lambda } => {
                    // W=mu/2 |F|²+g(J), with
                    // g''(J)=[mu+lambda*(1-ln J)]/J² and
                    // limit^-3 <= J <= limit³. Positive lambda has its
                    // worst numerator at maximum J; auxetic negative lambda
                    // has it at minimum J. Positive reference bulk alone
                    // does not certify this whole finite-strain domain.
                    let log_volume_span = 3.0 * fs_math::det::ln(self.stretch_limit);
                    let convex = if *lambda == 0.0 {
                        true
                    } else if *lambda > 0.0 {
                        log_volume_span <= 1.0 + mu / lambda
                    } else {
                        log_volume_span <= mu / (-lambda) - 1.0
                    };
                    convex.then_some(true)
                }
                HyperelasticModel::MooneyRivlin { c01, .. } => {
                    // |F|²/J^(2/3) is a convex power perspective. The
                    // |cof F|²/J^(4/3) term lacks that sufficient condition;
                    // absence of this certificate is not non-polyconvexity.
                    (*c01 == 0.0).then_some(true)
                }
                HyperelasticModel::Gent { volumetric, .. } => {
                    // -ln(1-q/Jm) is convex and increasing on q<Jm.
                    // Compose it with the convex power perspective
                    // q=|F|²/J^(2/3)-3; the chain domain is its strict sublevel.
                    self.convex_bulk(*volumetric).then_some(true)
                }
                HyperelasticModel::Ogden { terms, volumetric } => {
                    // Schatten powers times J^-b are convex for p>=1,
                    // 0<=b<=p-1 (power perspective); positive exponents use
                    // F, negative exponents cof F. Ball (1976), polyconvexity.
                    let powers = terms.iter().all(|p| p.alpha >= 1.5 || p.alpha <= -3.0);
                    (powers && self.convex_bulk(*volumetric)).then_some(true)
                }
                HyperelasticModel::HolzapfelGasserOgden {
                    invariants,
                    volumetric,
                    k1,
                    fibres,
                    ..
                } => {
                    // Full structural E is a convex quadratic in F; the
                    // exponential of <E>_+² is convex and nondecreasing.
                    ((*invariants == FibreInvariants::Full || *k1 == 0.0 || fibres.is_empty())
                        && self.convex_bulk(*volumetric))
                    .then_some(true)
                }
            },
            tangent_symmetric: true,
            failure_envelope: if matches!(self.model, HyperelasticModel::Gent { .. }) {
                "no damage law; strict Gent chain limit and principal-stretch domain"
            } else {
                "none declared (validity ends at stretch_limit)"
            },
        }
    }

    fn convex_bulk(&self, volumetric: Volumetric) -> bool {
        match volumetric {
            Volumetric::Quadratic { .. } | Volumetric::SimoTaylor { .. } => true,
            // U'' = kappa (1-ln J)/J². On the declared stretch box,
            // J <= limit³, so this is a sufficient domain-local condition.
            Volumetric::LogSquared { .. } => 3.0 * fs_math::det::ln(self.stretch_limit) <= 1.0,
        }
    }

    /// The model card.
    #[must_use]
    pub fn card(&self) -> ModelCard {
        let (name, disc) = match &self.model {
            HyperelasticModel::NeoHookean { .. } => ("material.neo-hookean", 0.03),
            HyperelasticModel::MooneyRivlin { .. } => ("material.mooney-rivlin", 0.03),
            // No empirical discrepancy band is inferred from a model family.
            HyperelasticModel::Gent { .. } => ("material.gent-decoupled", f64::INFINITY),
            HyperelasticModel::Ogden { .. } => ("material.ogden-decoupled", 0.03),
            HyperelasticModel::HolzapfelGasserOgden { .. } => ("material.hgo-dispersed", 0.03),
        };
        let symmetry = match &self.model {
            HyperelasticModel::HolzapfelGasserOgden { invariants, .. } => format!(
                "reference-frame fibre families; {invariants:?} fibre invariants; tension-only structural strain"
            ),
            _ => "isotropic".to_string(),
        };
        let mut assumptions = vec![
            "isothermal".to_string(),
            symmetry,
            "rate independence".to_string(),
        ];
        let mut failures = vec!["no damage/Mullins effect".to_string()];
        if let HyperelasticModel::Gent {
            locking_invariant, ..
        } = &self.model
        {
            assumptions.push(format!("strict isochoric I1-3 < {locking_invariant}"));
            failures.push("no empirical model-discrepancy bound supplied".to_string());
        }
        ModelCard::new(
            name,
            "0.1.0",
            Ambition::Solid,
            assumptions,
            ValidityDomain::unconstrained().with(
                "stretch",
                1.0 / self.stretch_limit,
                self.stretch_limit,
            ),
            failures,
            disc,
        )
    }
}
