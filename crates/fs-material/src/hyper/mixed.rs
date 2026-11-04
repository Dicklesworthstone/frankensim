//! Mixed displacement/pressure potentials from the original decoupled laws.
//!
//! Pressure is tensile-positive: `p=U'(J)` and
//! `Pi(F,p)=W_remainder(F)+p*(J-1)-U_star(p)`. The pressure residual is
//! `J-J(p)`, and both cross derivatives equal `cof(F)`. The native remainder
//! retains full-invariant HGO fibres, so its pressure need not equal mean
//! Cauchy stress. Coupled Neo-Hookean is deliberately not split here.
//!
//! Quadratic bulk energy has its exact quadratic conjugate. Log-squared bulk
//! energy uses ONLY the locally convex branch `ln(J)<1`: its positive-pressure
//! inverse is multivalued outside that branch, and a global finite Legendre
//! transform is not claimed. No bulk modulus is replaced. Exact incompressibility
//! is an explicit alternative, with zero pressure compliance and no bulk energy.
//! A spatial mixed discretization must separately establish inf-sup stability.

use super::{Hyperelastic, HyperelasticModel, Volumetric, det3};
use crate::MaterialError;
use fs_ad::{Dual, Dual64};
use fs_math::det;

/// The original finite-compliance material or an explicitly imposed constraint.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VolumetricTreatment {
    /// `J=1`, with pressure as its independent Lagrange multiplier.
    ExactIncompressible,
    /// Preserve the original `U(J)` on its invertible, positive-curvature branch.
    OriginalMaterial,
}

/// Exact local pressure conjugate of the selected original volumetric energy.
#[derive(Clone, Copy, Debug)]
pub struct PressurePotential {
    /// `p*(J(p)-1)-U(J(p))`, in Pa.
    pub complementary_pa: f64,
    /// Original physical bulk storage `U(J(p))`, in Pa.
    pub bulk_stored_pa: f64,
    /// Volume reconstructed from pressure; independent of the current trial F.
    pub volume_ratio: f64,
    /// `J(p)-1`, evaluated without cancellation near the reference state.
    pub volume_increment: f64,
    /// `dJ(p)/dp=1/U''(J(p))`, in Pa^-1.
    pub compliance_per_pa: f64,
}

/// Mixed saddle potential and residual at one reference-volume material point.
#[derive(Clone, Copy, Debug)]
pub struct MixedPointState {
    /// Complete mixed potential per reference volume, in Pa.
    pub potential_pa: f64,
    /// Elastic remainder, or INCREMENTAL remainder for a frozen viscous step.
    /// The latter is not physical stored energy or heat.
    pub remainder_potential_pa: f64,
    /// Physical bulk storage at the pressure-reconstructed volume, in Pa.
    pub bulk_stored_pa: f64,
    /// Pressure-reconstructed volume `J(p)`; one for an exact constraint.
    pub bulk_volume_ratio: f64,
    /// Pressure complementary potential, in Pa; zero for an exact constraint.
    pub bulk_complementary_pa: f64,
    /// First Piola stress, row-major, in Pa.
    pub first_piola_pa: [f64; 9],
    /// `dPi/dp = J(F)-J(p)`, dimensionless.
    pub pressure_residual: f64,
    /// Actual trial volume `det(F)`; weak pressure equations need not make
    /// this equal to `bulk_volume_ratio` at each quadrature point.
    pub volume_ratio: f64,
}

/// All four blocks of the symmetric mixed constitutive Jacobian.
#[derive(Clone, Copy, Debug)]
pub struct MixedPointResponse {
    /// Potential and first derivatives at the same point.
    pub state: MixedPointState,
    /// `dP/dF` at FIXED pressure, in Pa; includes `p*d²J/dF²`.
    pub tangent_pa: [[f64; 9]; 9],
    /// Both `dP/dp` and `d(pressure_residual)/dF`, in row-major order.
    pub pressure_coupling: [f64; 9],
    /// `d(pressure_residual)/dp=-dJ(p)/dp`, in Pa^-1.
    pub pressure_tangent_per_pa: f64,
}

/// Original bulk-free remainder and its exact Hessian at fixed history/time.
#[derive(Clone, Copy, Debug)]
pub struct RemainderResponse {
    /// Stored or incremental remainder potential, in Pa.
    pub potential_pa: f64,
    /// Its first Piola gradient, row-major, in Pa.
    pub first_piola_pa: [f64; 9],
    /// Its consistent Hessian, in Pa.
    pub tangent_pa: [[f64; 9]; 9],
}

fn state_error(what: &str) -> MaterialError {
    MaterialError::State {
        what: format!("mixed finite-strain potential: {what}"),
    }
}

pub(crate) fn original_volume(law: &Hyperelastic) -> Result<Volumetric, MaterialError> {
    let volume = match &law.model {
        HyperelasticModel::NeoHookean { .. } => {
            return Err(MaterialError::Parameters {
                what: "coupled Neo-Hookean has no admitted mixed isochoric/volumetric split".into(),
            });
        }
        HyperelasticModel::MooneyRivlin { kappa, .. } => Volumetric::Quadratic { kappa: *kappa },
        HyperelasticModel::Ogden { volumetric, .. }
        | HyperelasticModel::HolzapfelGasserOgden { volumetric, .. } => *volumetric,
    };
    if !(volume.kappa().is_finite() && volume.kappa() > 0.0) {
        return Err(MaterialError::Parameters {
            what: "mixed bulk modulus must be finite and positive".into(),
        });
    }
    Ok(volume)
}

/// Fixed scalar-inversion work cap; it does not alter any constitutive coefficient.
pub const MAX_PRESSURE_INVERSE_ITERATIONS: usize = 96;

// Solve y=x*exp(-x), x<1. Close to the fold, t=1-x and
// 1-(1-t)*exp(t)=sum_{n>=2}(n-1)t^n/n! avoid a cancellation of order t².
fn inverse_log_pressure(y: f64) -> Result<(f64, f64), MaterialError> {
    let ceiling = det::exp(-1.0);
    if !y.is_finite() || y >= ceiling {
        return Err(state_error(
            "log-squared pressure requires p/kappa < exp(-1)",
        ));
    }
    if y == 0.0 {
        return Ok((0.0, 1.0));
    }
    if y > 0.125 {
        let target = (ceiling - y) / ceiling;
        let (mut lo, mut hi) = (0.0, 1.0);
        let mut t = det::sqrt(2.0 * target).min(1.0);
        for _ in 0..MAX_PRESSURE_INVERSE_ITERATIONS {
            let value = if t <= 0.5 {
                let mut term = 0.5 * t * t;
                let mut sum = term;
                for n in 3..=20 {
                    term *= t * (n - 1) as f64 / (n * (n - 2)) as f64;
                    sum += term;
                }
                sum
            } else {
                1.0 + (t - 1.0) * det::exp(t)
            };
            let residual = value - target;
            let correction = residual / (t * det::exp(t));
            if correction.abs() <= 16.0 * f64::EPSILON * t {
                return Ok((1.0 - t, t));
            }
            if residual < 0.0 {
                lo = t;
            } else {
                hi = t;
            }
            let next = t - correction;
            t = if next > lo && next < hi {
                next
            } else {
                0.5 * (lo + hi)
            };
        }
    } else {
        let (mut lo, mut hi) = if y > 0.0 {
            (0.0, 1.0)
        } else if y >= -1.0 {
            (y, 0.0)
        } else {
            (-det::ln(1.0 - y), 0.0)
        };
        let mut x = if y < -1.0 {
            lo + det::ln(1.0 - lo)
        } else {
            y / (1.0 - y)
        };
        for _ in 0..MAX_PRESSURE_INVERSE_ITERATIONS {
            let term = y * det::exp(x);
            let residual = x - term;
            let correction = residual / (1.0 - term);
            if correction.abs() <= 16.0 * f64::EPSILON * x.abs() {
                return Ok((x, 1.0 - x));
            }
            if residual < 0.0 {
                lo = x;
            } else {
                hi = x;
            }
            let next = x - correction;
            x = if next > lo && next < hi {
                next
            } else {
                0.5 * (lo + hi)
            };
        }
    }
    Err(state_error(
        "original pressure inverse exhausted its scalar iteration budget",
    ))
}

/// Original local volume/pressure dual. Positive pressure is tensile.
/// Log-squared volume uses its locally convex `ln(J)<1` branch only.
pub fn pressure_potential(
    volume: Volumetric,
    pressure_pa: f64,
) -> Result<PressurePotential, MaterialError> {
    let kappa = volume.kappa();
    if !(kappa.is_finite() && kappa > 0.0) {
        return Err(MaterialError::Parameters {
            what: "mixed bulk modulus must be finite and positive".into(),
        });
    }
    if !pressure_pa.is_finite() {
        return Err(state_error("pressure must be finite"));
    }
    let result = match volume {
        Volumetric::Quadratic { .. } => {
            let increment = pressure_pa / kappa;
            let stored = 0.5 * pressure_pa * increment;
            PressurePotential {
                complementary_pa: stored,
                bulk_stored_pa: stored,
                volume_ratio: 1.0 + increment,
                volume_increment: increment,
                compliance_per_pa: 1.0 / kappa,
            }
        }
        Volumetric::LogSquared { .. } => {
            let (x, distance) = inverse_log_pressure(pressure_pa / kappa)?;
            let j = det::exp(x);
            let increment = det::expm1(x);
            let stored = 0.5 * kappa * x * x;
            let mut compliance = (j * j) / (kappa * distance);
            if !(compliance.is_finite() && compliance > 0.0) {
                compliance = det::exp(2.0 * x - det::ln(kappa) - det::ln(distance));
            }
            PressurePotential {
                complementary_pa: pressure_pa * increment - stored,
                bulk_stored_pa: stored,
                volume_ratio: j,
                volume_increment: increment,
                compliance_per_pa: compliance,
            }
        }
    };
    if ![
        result.complementary_pa,
        result.bulk_stored_pa,
        result.volume_ratio,
        result.volume_increment,
        result.compliance_per_pa,
    ]
    .iter()
    .all(|v| v.is_finite())
        || result.volume_ratio <= 0.0
        || result.compliance_per_pa <= 0.0
    {
        return Err(state_error(
            "original pressure dual is outside its positive-volume branch or unrepresentable",
        ));
    }
    Ok(result)
}

fn pressure_data(
    law: &Hyperelastic,
    f: &[f64; 9],
    p: f64,
    treatment: VolumetricTreatment,
) -> Result<PressurePotential, MaterialError> {
    let volume = original_volume(law)?;
    if !p.is_finite() {
        return Err(state_error("pressure must be finite"));
    }
    match treatment {
        VolumetricTreatment::ExactIncompressible => Ok(PressurePotential {
            complementary_pa: 0.0,
            bulk_stored_pa: 0.0,
            volume_ratio: 1.0,
            volume_increment: 0.0,
            compliance_per_pa: 0.0,
        }),
        VolumetricTreatment::OriginalMaterial => {
            if matches!(volume, Volumetric::LogSquared { .. }) && det::ln(det3(f)) >= 1.0 {
                return Err(state_error(
                    "actual log-squared material volume requires ln(det F)<1 on the selected branch",
                ));
            }
            pressure_potential(volume, p)
        }
    }
}

pub(crate) fn cofactor(f: &[f64; 9]) -> [f64; 9] {
    [
        f[4] * f[8] - f[5] * f[7],
        f[5] * f[6] - f[3] * f[8],
        f[3] * f[7] - f[4] * f[6],
        f[2] * f[7] - f[1] * f[8],
        f[0] * f[8] - f[2] * f[6],
        f[1] * f[6] - f[0] * f[7],
        f[1] * f[5] - f[2] * f[4],
        f[2] * f[3] - f[0] * f[5],
        f[0] * f[4] - f[1] * f[3],
    ]
}

fn combine_state(
    f: &[f64; 9],
    p: f64,
    remainder: ([f64; 9], f64),
    bulk: PressurePotential,
) -> Result<MixedPointState, MaterialError> {
    let cof = cofactor(f);
    let j = det3(f);
    let potential = remainder.1 + p * (j - 1.0) - bulk.complementary_pa;
    let stress = core::array::from_fn(|i| remainder.0[i] + p * cof[i]);
    if !potential.is_finite() || stress.iter().any(|v| !v.is_finite()) {
        return Err(state_error("nonfinite mixed potential or stress"));
    }
    Ok(MixedPointState {
        potential_pa: potential,
        remainder_potential_pa: remainder.1,
        bulk_stored_pa: bulk.bulk_stored_pa,
        bulk_volume_ratio: bulk.volume_ratio,
        bulk_complementary_pa: bulk.complementary_pa,
        first_piola_pa: stress,
        pressure_residual: (j - 1.0) - bulk.volume_increment,
        volume_ratio: j,
    })
}

/// Bulk-free finite-strain energy owner, including a frozen viscous value
/// function when selected. `validate_mixed` admits parameters/history once at
/// spatial preparation; each point still checks original deformation domains.
pub trait MixedFiniteStrainPotential {
    /// Original immutable equilibrium law, owning the selected bulk energy.
    fn equilibrium_law(&self) -> &Hyperelastic;
    /// Validate the original public material parameters and mixed split.
    fn validate_mixed(&self) -> Result<(), MaterialError> {
        let law = self.equilibrium_law();
        Hyperelastic::new(law.model.clone(), law.stretch_limit)?;
        original_volume(law).map(|_| ())
    }
    /// Bulk-free remainder gradient and value, without its Hessian.
    fn remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError>;
    /// Bulk-free remainder gradient, value and Hessian from one prepared point.
    fn remainder_response(&self, f: &[f64; 9]) -> Result<RemainderResponse, MaterialError>;
    /// Original pressure/volume dual without evaluating the remainder Hessian.
    /// A prescribed-displacement pressure solve needs only this scalar law.
    /// Original deformation and bulk-branch restrictions still apply; callers
    /// separately admit the remainder/history through the native residual.
    fn mixed_pressure_potential(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<PressurePotential, MaterialError> {
        check_split_state(self.equilibrium_law(), f)?;
        pressure_data(self.equilibrium_law(), f, pressure_pa, treatment)
    }
    /// Mixed potential and residuals for a spatial residual line search.
    fn mixed_state(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedPointState, MaterialError> {
        let remainder = self.remainder_stress_and_potential(f)?;
        let bulk = pressure_data(self.equilibrium_law(), f, pressure_pa, treatment)?;
        combine_state(f, pressure_pa, remainder, bulk)
    }
    /// All consistent mixed blocks at fixed accepted history and interval.
    fn mixed_response(
        &self,
        f: &[f64; 9],
        pressure_pa: f64,
        treatment: VolumetricTreatment,
    ) -> Result<MixedPointResponse, MaterialError> {
        let remainder = self.remainder_response(f)?;
        let bulk = pressure_data(self.equilibrium_law(), f, pressure_pa, treatment)?;
        let state = combine_state(
            f,
            pressure_pa,
            (remainder.first_piola_pa, remainder.potential_pa),
            bulk,
        )?;
        let mut tangent = remainder.tangent_pa;
        // d²J/dF_iJ dF_kL = eps_ikm eps_JLN F_mN. Only distinct
        // rows and columns survive; no determinant Hessian AD is necessary.
        for (a, row) in tangent.iter_mut().enumerate() {
            for (b, value) in row.iter_mut().enumerate() {
                let (i, j, k, l) = (a / 3, a % 3, b / 3, b % 3);
                if i != k && j != l {
                    let sign_i = if (i + 1) % 3 == k { 1.0 } else { -1.0 };
                    let sign_j = if (j + 1) % 3 == l { 1.0 } else { -1.0 };
                    *value += pressure_pa * sign_i * sign_j * f[3 * (3 - i - k) + (3 - j - l)];
                }
            }
        }
        if tangent.iter().flatten().any(|v| !v.is_finite()) {
            return Err(state_error("nonfinite mixed tangent"));
        }
        Ok(MixedPointResponse {
            state,
            tangent_pa: tangent,
            pressure_coupling: cofactor(f),
            pressure_tangent_per_pa: -bulk.compliance_per_pa,
        })
    }
}

pub(crate) fn check_split_state(law: &Hyperelastic, f: &[f64; 9]) -> Result<(), MaterialError> {
    original_volume(law)?;
    law.check_state(f)
}

pub(crate) fn checked_stress(
    value: f64,
    stress: [f64; 9],
) -> Result<([f64; 9], f64), MaterialError> {
    if !value.is_finite() || stress.iter().any(|v| !v.is_finite()) {
        return Err(state_error("nonfinite original remainder energy or stress"));
    }
    Ok((stress, value))
}

// Reuse the original generic-scalar energy. Value and gradient come from the
// first Hessian sweep, so a viscous caller prepares the branch return only once.
pub(crate) fn response_from_energy(
    f: &[f64; 9],
    energy: impl Fn(&[Dual<Dual64<9>, 1>; 9]) -> Dual<Dual64<9>, 1>,
) -> Result<RemainderResponse, MaterialError> {
    let mut result = RemainderResponse {
        potential_pa: 0.0,
        first_piola_pa: [0.0; 9],
        tangent_pa: [[0.0; 9]; 9],
    };
    for col in 0..9 {
        let vars = core::array::from_fn(|row| Dual {
            re: Dual64::variable(f[row], row),
            eps: [Dual64::constant(if row == col { 1.0 } else { 0.0 })],
        });
        let w = energy(&vars);
        if col == 0 {
            result.potential_pa = w.re.re;
            result.first_piola_pa = w.re.eps;
        }
        for row in 0..9 {
            result.tangent_pa[row][col] = w.eps[0].eps[row];
        }
    }
    checked_stress(result.potential_pa, result.first_piola_pa)?;
    if result.tangent_pa.iter().flatten().any(|v| !v.is_finite()) {
        return Err(state_error("nonfinite original remainder tangent"));
    }
    Ok(result)
}

impl MixedFiniteStrainPotential for Hyperelastic {
    fn equilibrium_law(&self) -> &Hyperelastic {
        self
    }
    fn remainder_stress_and_potential(
        &self,
        f: &[f64; 9],
    ) -> Result<([f64; 9], f64), MaterialError> {
        check_split_state(self, f)?;
        let (value, stress) = fs_ad::gradient::<9>(*f, |x| self.isochoric_energy(&x));
        checked_stress(value, stress)
    }
    fn remainder_response(&self, f: &[f64; 9]) -> Result<RemainderResponse, MaterialError> {
        check_split_state(self, f)?;
        response_from_energy(f, |x| self.isochoric_energy(x))
    }
}
