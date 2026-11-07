//! Multiplicative pressure-sensitive plasticity with retained plastic volume.
//!
//! F=Fe Fp, C=Cp^-1 and b_trial=F C F^T. The shared porous return acts on
//! e_trial=.5 log(b_trial), with Kirchhoff stress and energy per ORIGINAL
//! volume. The corrected metric is C_new=F^-1 exp(2e_new) F^-T: its determinant
//! is never normalized. Hydrostatic compression therefore changes Jp.
//!
//! The exponential/logarithmic construction follows the isotropic return-map
//! family of Simo (1992), doi:10.1016/0045-7825(92)90123-2. The normalized
//! Deshpande--Fleck ellipse and energetic hardening are shared with
//! plastic::porous::return_map. This is an explicitly authored finite law:
//! small-strain or Cauchy-yield calibration is not silently promoted to it.
//!
//! Large rotation and bounded permanent volume change are admitted. The model
//! does not evolve porosity, cell-wall densification, damage, pore-fluid
//! pressure, temperature or rate dependence. The yield surface is symmetric
//! in mean Kirchhoff stress. All physical applicability bounds are explicit.

/// Explicit source binding for the native finite pressure-sensitive law.
pub mod source;

/// Thermal expansion, coefficient evolution and heat through the same compaction history.
pub mod thermal;

mod voce;

use super::voce::FiniteVoceHardening;

use super::{
    FiniteJ2Error, IDENTITY6, ROUND_OFF, Spectrum, determinant, deviator, dot,
    elastic_metric_spectrum, expand, inverse, left_metric_increment, multiply, norm, pack,
    transpose,
};
use crate::hyper::cards::FiniteStrainPotential;
use crate::plastic::porous::return_map::{PressureReturn, ReturnParameters};
use crate::{MaterialError, Voigt};
use core::fmt;
use fs_math::det;

/// Fixed coefficients. Strengths and hardening are Kirchhoff/reference-volume quantities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinitePorousParameters {
    pub shear_modulus_pa: f64,
    pub bulk_modulus_pa: f64,
    pub yield_stress_pa: f64,
    pub hardening_modulus_pa: f64,
    /// Dimensionless normalized Deshpande--Fleck ellipse parameter; zero is J2.
    pub pressure_sensitivity: f64,
}

/// Authored bounds; no material-dependent values are inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FinitePorousLimits {
    /// Strict lower bound on det(F), between zero and one.
    pub minimum_jacobian: f64,
    /// Maximum Frobenius norm of the returned elastic logarithmic strain.
    pub maximum_elastic_log_strain: f64,
    /// Maximum Frobenius norm of .5 log(Cp^-1), including its volume part.
    pub maximum_plastic_log_strain: f64,
    /// Maximum accumulated work-conjugate plastic strain, including reversals.
    pub maximum_accumulated_plastic_strain: f64,
}

/// Complete caller-owned history. The metric is SPD and may change volume.
#[derive(Clone, Debug, PartialEq)]
pub struct FinitePorousState {
    /// Tensor-Voigt [xx,yy,zz,xy,yz,zx], NOT an isochoric metric when a>0.
    pub inverse_plastic_cauchy_green: Voigt,
    pub alpha: f64,
}

/// Immutable native law. One common stress-space return serves both kinematics.
#[derive(Clone, Debug, PartialEq)]
pub struct FinitePorous {
    parameters: FinitePorousParameters,
    limits: FinitePorousLimits,
    kernel: ReturnParameters,
    voce: Option<FiniteVoceHardening>,
}

/// Borrow the same native law from an authored or source-resolved material.
pub trait FinitePorousMaterial {
    /// Immutable coefficients and native compression/history limits.
    fn finite_porous(&self) -> &FinitePorous;
}

impl FinitePorousMaterial for FinitePorous {
    fn finite_porous(&self) -> &FinitePorous {
        self
    }
}

impl<T: FinitePorousMaterial + ?Sized> FinitePorousMaterial for &T {
    fn finite_porous(&self) -> &FinitePorous {
        (*self).finite_porous()
    }
}

/// Caller commits the candidate only after its complete coupled step succeeds.
#[derive(Clone, Debug, PartialEq)]
pub struct FinitePorousResponse {
    pub state: FinitePorousState,
    /// Work conjugate to F, row-major, per reference area.
    pub first_piola_pa: [f64; 9],
    /// tau=J sigma, row-major.
    pub kirchhoff_pa: [f64; 9],
    /// Exact dP/dF at frozen accepted history, including geometric terms.
    pub tangent_pa: [[f64; 9]; 9],
    /// Stored endpoint energy plus yield_stress * plastic_increment.
    /// Its gradient is P; this is not the endpoint stored energy.
    pub incremental_potential_j_m3: f64,
    pub stored_energy_j_m3: f64,
    /// Physical yield heat per original volume, excluding numerical work defect.
    pub plastic_dissipation_j_m3: f64,
    pub plastic_increment: f64,
    pub equivalent_kirchhoff_pa: f64,
    pub yield_residual_pa: f64,
    pub jacobian: f64,
    /// det(Fp), retained independently of the total current volume.
    pub plastic_jacobian: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FinitePorousError {
    Parameters(&'static str),
    State(&'static str),
    Domain {
        quantity: &'static str,
        value: f64,
        limit: f64,
    },
    YieldResidual {
        residual_pa: f64,
        tolerance_pa: f64,
    },
    /// The shared native SPD logarithm, inverse or metric operation refused.
    Metric(FiniteJ2Error),
    /// The shared pressure-sensitive scalar return refused.
    Return(MaterialError),
    Numerical(&'static str),
}
impl fmt::Display for FinitePorousError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parameters(x) => write!(f, "finite porous parameters: {x}"),
            Self::State(x) => write!(f, "finite porous history: {x}"),
            Self::Domain {
                quantity,
                value,
                limit,
            } => write!(f, "finite porous {quantity}={value}, bound={limit}"),
            Self::YieldResidual {
                residual_pa,
                tolerance_pa,
            } => write!(
                f,
                "finite porous yield residual {residual_pa} Pa exceeds {tolerance_pa} Pa"
            ),
            Self::Metric(x) => write!(f, "finite porous metric operation: {x}"),
            Self::Return(x) => write!(f, "finite porous constitutive return: {x}"),
            Self::Numerical(x) => write!(f, "finite porous numerical refusal: {x}"),
        }
    }
}
impl std::error::Error for FinitePorousError {}
impl From<FiniteJ2Error> for FinitePorousError {
    fn from(e: FiniteJ2Error) -> Self {
        Self::Metric(e)
    }
}
impl From<MaterialError> for FinitePorousError {
    fn from(e: MaterialError) -> Self {
        Self::Return(e)
    }
}

struct Endpoint {
    jacobian: f64,
    log_plastic_jacobian: f64,
    spectrum: Spectrum,
    elastic: Voigt,
}
struct Returned {
    trial: Endpoint,
    pressure: PressureReturn,
    state: FinitePorousState,
    inverse_transpose: [f64; 9],
    tau: [f64; 9],
    piola: [f64; 9],
    stored: f64,
    heat: f64,
    potential: f64,
    plastic_jacobian: f64,
}

fn bound(quantity: &'static str, value: f64, limit: f64) -> Result<(), FinitePorousError> {
    if !value.is_finite() || value > limit {
        Err(FinitePorousError::Domain {
            quantity,
            value,
            limit,
        })
    } else {
        Ok(())
    }
}
fn trace(a: &[f64; 9]) -> f64 {
    a[0] + a[4] + a[8]
}
fn volume(log_jp: f64) -> Result<f64, FinitePorousError> {
    let jp = det::exp(log_jp);
    if !(jp.is_finite() && jp > 0.0) {
        Err(FinitePorousError::Numerical(
            "plastic volume ratio is unrepresentable",
        ))
    } else {
        Ok(jp)
    }
}

impl FinitePorous {
    pub fn new(
        parameters: FinitePorousParameters,
        limits: FinitePorousLimits,
    ) -> Result<Self, FinitePorousError> {
        if !(limits.minimum_jacobian.is_finite()
            && limits.minimum_jacobian > 0.0
            && limits.minimum_jacobian < 1.0)
            || [
                limits.maximum_elastic_log_strain,
                limits.maximum_plastic_log_strain,
                limits.maximum_accumulated_plastic_strain,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(FinitePorousError::Parameters(
                "positive finite strain bounds and 0<minimum_jacobian<1 required",
            ));
        }
        let p = &parameters;
        let kernel = ReturnParameters::new(
            p.bulk_modulus_pa,
            p.shear_modulus_pa,
            p.yield_stress_pa,
            p.hardening_modulus_pa,
            p.pressure_sensitivity,
        )?;
        if !(p.yield_stress_pa + p.hardening_modulus_pa * limits.maximum_accumulated_plastic_strain)
            .is_finite()
        {
            return Err(FinitePorousError::Parameters(
                "yield strength must remain representable throughout the history domain",
            ));
        }
        Ok(Self {
            parameters,
            limits,
            kernel,
            voce: None,
        })
    }

    pub fn parameters(&self) -> &FinitePorousParameters {
        &self.parameters
    }
    pub fn limits(&self) -> &FinitePorousLimits {
        &self.limits
    }
    pub fn initial_state(&self) -> FinitePorousState {
        FinitePorousState {
            inverse_plastic_cauchy_green: IDENTITY6,
            alpha: 0.0,
        }
    }

    /// Necessary reachability bound in the weighted congruence-invariant SPD metric.
    /// A closed plastic path may restore the identity with positive accumulation.
    pub fn validate_history(&self, state: &FinitePorousState) -> Result<(), FinitePorousError> {
        self.history_log_volume(state).map(|_| ())
    }

    fn history_log_volume(&self, state: &FinitePorousState) -> Result<f64, FinitePorousError> {
        if !(state.alpha.is_finite() && state.alpha >= 0.0)
            || state
                .inverse_plastic_cauchy_green
                .iter()
                .any(|x| !x.is_finite())
        {
            return Err(FinitePorousError::State(
                "finite metric and nonnegative accumulation required",
            ));
        }
        bound(
            "accumulated plastic strain",
            state.alpha,
            self.limits.maximum_accumulated_plastic_strain,
        )?;
        if state.inverse_plastic_cauchy_green == IDENTITY6 {
            return Ok(0.0);
        }
        if state.alpha == 0.0 {
            return Err(FinitePorousError::State(
                "zero accumulation requires the virgin plastic metric",
            ));
        }
        let spectrum = Spectrum::spd(&expand(&state.inverse_plastic_cauchy_green))?;
        let h = spectrum.logarithm.map(|x| 0.5 * x);
        let magnitude = norm(&h);
        bound(
            "plastic logarithmic strain",
            magnitude,
            self.limits.maximum_plastic_log_strain,
        )?;
        let tr = trace(&h);
        let trace_roundoff = ROUND_OFF * (1.0 + magnitude);
        let a = self.parameters.pressure_sensitivity;
        let volumetric = if a == 0.0 {
            if tr.abs() > trace_roundoff {
                return Err(FinitePorousError::State(
                    "zero pressure sensitivity requires an isochoric plastic metric",
                ));
            }
            0.0
        } else {
            (tr.abs() - trace_roundoff).max(0.0) / a
        };
        let dev = det::sqrt(2.0 / 3.0) * norm(&deviator(&h));
        let required = det::sqrt(self.kernel.c) * norm(&[dev, volumetric]);
        if !required.is_finite()
            || required > state.alpha + ROUND_OFF * (1.0 + required + state.alpha)
        {
            return Err(FinitePorousError::State(
                "plastic metric exceeds accumulated work-conjugate flow",
            ));
        }
        volume(-tr)?;
        Ok(-tr)
    }

    pub fn plastic_jacobian(&self, state: &FinitePorousState) -> Result<f64, FinitePorousError> {
        volume(self.history_log_volume(state)?)
    }

    fn endpoint(
        &self,
        f: &[f64; 9],
        state: &FinitePorousState,
    ) -> Result<Endpoint, FinitePorousError> {
        let log_plastic_jacobian = self.history_log_volume(state)?;
        if f.iter().any(|x| !x.is_finite()) {
            return Err(FinitePorousError::Numerical("F must be finite"));
        }
        let jacobian = determinant(f);
        if !jacobian.is_finite() || jacobian <= self.limits.minimum_jacobian {
            return Err(FinitePorousError::Domain {
                quantity: "det(F), strict lower bound",
                value: jacobian,
                limit: self.limits.minimum_jacobian,
            });
        }
        let cp_inv = expand(&state.inverse_plastic_cauchy_green);
        let spectrum = elastic_metric_spectrum(f, &cp_inv, &left_metric_increment(f))?;
        let elastic = pack(&spectrum.logarithm.map(|x| 0.5 * x));
        Ok(Endpoint {
            jacobian,
            log_plastic_jacobian,
            spectrum,
            elastic,
        })
    }

    fn elastic_bound(&self, strain: &Voigt) -> Result<(), FinitePorousError> {
        bound(
            "elastic logarithmic strain",
            norm(&expand(strain)),
            self.limits.maximum_elastic_log_strain,
        )
    }

    fn energy(&self, elastic: &Voigt, alpha: f64) -> Result<f64, FinitePorousError> {
        let e = expand(elastic);
        let dev = deviator(&e);
        let p = &self.parameters;
        let result = p.shear_modulus_pa * dot(&dev, &dev)
            + 0.5 * p.bulk_modulus_pa * trace(&e).powi(2)
            + 0.5 * p.hardening_modulus_pa * alpha * alpha
            + self.voce.map_or(0.0, |v| v.stored_energy(alpha));
        if !result.is_finite() {
            Err(FinitePorousError::Numerical(
                "stored energy is unrepresentable",
            ))
        } else {
            Ok(result)
        }
    }

    fn equivalent(&self, stress: &Voigt) -> f64 {
        let tau = expand(stress);
        let root_c = det::sqrt(self.kernel.c);
        norm(&[
            det::sqrt(1.5) * norm(&deviator(&tau)) / root_c,
            (self.parameters.pressure_sensitivity / root_c) * (trace(&tau) / 3.0),
        ])
    }

    /// Actual storage at the supplied endpoint; never repairs history by another return.
    pub fn stored_energy(
        &self,
        f: &[f64; 9],
        state: &FinitePorousState,
    ) -> Result<f64, FinitePorousError> {
        let trial = self.endpoint(f, state)?;
        self.elastic_bound(&trial.elastic)?;
        self.energy(&trial.elastic, state.alpha)
    }

    pub fn validate_committed(
        &self,
        f: &[f64; 9],
        state: &FinitePorousState,
        yield_tolerance_pa: f64,
    ) -> Result<(), FinitePorousError> {
        if !(yield_tolerance_pa.is_finite() && yield_tolerance_pa >= 0.0) {
            return Err(FinitePorousError::Parameters(
                "finite nonnegative yield allowance required",
            ));
        }
        let trial = self.endpoint(f, state)?;
        self.elastic_bound(&trial.elastic)?;
        self.energy(&trial.elastic, state.alpha)?;
        let e = expand(&trial.elastic);
        let dev = deviator(&e);
        let stress = pack(&core::array::from_fn(|i| {
            2.0 * self.parameters.shear_modulus_pa * dev[i]
                + if i % 4 == 0 {
                    self.parameters.bulk_modulus_pa * trace(&e)
                } else {
                    0.0
                }
        }));
        let residual_pa = self.equivalent(&stress) - self.yield_radius(state.alpha);
        if !residual_pa.is_finite() {
            return Err(FinitePorousError::Numerical(
                "accepted yield residual is unrepresentable",
            ));
        }
        if residual_pa > yield_tolerance_pa {
            Err(FinitePorousError::YieldResidual {
                residual_pa,
                tolerance_pa: yield_tolerance_pa,
            })
        } else {
            Ok(())
        }
    }

    fn return_map(
        &self,
        f: &[f64; 9],
        old: &FinitePorousState,
    ) -> Result<Returned, FinitePorousError> {
        let trial = self.endpoint(f, old)?;
        let p = &self.parameters;
        // Logarithm/reconstruction roundoff at yield must not create flow on
        // exact holds or pure rotations. This is not a physical overstress law.
        let zero = 64.0
            * f64::EPSILON
            * (2.0 * p.shear_modulus_pa)
                .max(p.bulk_modulus_pa)
                .max(self.yield_radius(old.alpha));
        let pressure = match self.voce {
            Some(ref saturation) => {
                self.kernel
                    .response_with_saturation(&trial.elastic, old.alpha, zero, saturation)?
            }
            None => self.kernel.response(&trial.elastic, old.alpha, zero)?,
        };
        self.elastic_bound(&pressure.elastic_strain)?;
        bound(
            "accumulated plastic strain",
            pressure.alpha,
            self.limits.maximum_accumulated_plastic_strain,
        )?;
        let inverse = inverse(f, trial.jacobian)?;
        let inverse_transpose = transpose(&inverse);
        let state = if pressure.multiplier == 0.0 {
            old.clone()
        } else {
            // Isotropic associated return commutes with the trial elastic metric.
            // Reuse its eigensystem; no extra eigenproblem or determinant repair.
            let elastic = expand(&pressure.elastic_strain);
            let rotated = multiply(
                &multiply(&transpose(&trial.spectrum.vectors), &elastic),
                &trial.spectrum.vectors,
            );
            let roots = [0, 4, 8].map(|i| det::exp(2.0 * rotated[i]));
            if roots.iter().any(|x| !x.is_finite() || *x <= 0.0) {
                return Err(FinitePorousError::Numerical(
                    "returned elastic metric is unrepresentable",
                ));
            }
            let be = trial.spectrum.reconstruct(&roots);
            let cp_inv = multiply(&multiply(&inverse, &be), &inverse_transpose);
            FinitePorousState {
                inverse_plastic_cauchy_green: pack(&cp_inv),
                alpha: pressure.alpha,
            }
        };
        let log_jp = self.history_log_volume(&state)?;
        let volume_flow = pressure.plastic_increment[..3].iter().sum::<f64>();
        let mismatch = log_jp - trial.log_plastic_jacobian - volume_flow;
        if !mismatch.is_finite()
            || mismatch.abs()
                > ROUND_OFF
                    * (1.0 + log_jp.abs() + trial.log_plastic_jacobian.abs() + volume_flow.abs())
        {
            return Err(FinitePorousError::Numerical(
                "plastic metric volume does not match the associated flow",
            ));
        }
        let plastic_jacobian = volume(log_jp)?;
        let tau = expand(&pressure.stress);
        let piola = multiply(&tau, &inverse_transpose);
        let stored = self.energy(&pressure.elastic_strain, pressure.alpha)?;
        let heat = p.yield_stress_pa * pressure.multiplier;
        let potential = stored + heat;
        if tau.iter().chain(piola.iter()).any(|x| !x.is_finite())
            || !potential.is_finite()
            || !(heat.is_finite() && heat >= 0.0)
            || (pressure.multiplier > 0.0 && heat == 0.0)
        {
            return Err(FinitePorousError::Numerical(
                "stress, potential or physical heat is unrepresentable",
            ));
        }
        Ok(Returned {
            trial,
            pressure,
            state,
            inverse_transpose,
            tau,
            piola,
            stored,
            heat,
            potential,
            plastic_jacobian,
        })
    }

    /// Exact discrete tangent through the shared Frechet logarithm and scalar return.
    pub fn response(
        &self,
        f: &[f64; 9],
        old: &FinitePorousState,
    ) -> Result<FinitePorousResponse, FinitePorousError> {
        let r = self.return_map(f, old)?;
        let f_cp = multiply(f, &expand(&old.inverse_plastic_cauchy_green));
        let differences = r.trial.spectrum.divided_differences();
        let mut tangent = [[0.0; 9]; 9];
        for column in 0..9 {
            let (row, col) = (column / 3, column % 3);
            let db = core::array::from_fn(|ij| {
                let (i, j) = (ij / 3, ij % 3);
                (if i == row { f_cp[3 * j + col] } else { 0.0 })
                    + (if j == row { f_cp[3 * i + col] } else { 0.0 })
            });
            let de = pack(&r.trial.spectrum.frechet(&db, &differences).map(|x| 0.5 * x));
            let dtau = expand(&core::array::from_fn(|i| {
                r.pressure.tangent[i]
                    .iter()
                    .zip(&de)
                    .map(|(c, e)| c * e)
                    .sum()
            }));
            let dp = multiply(&dtau, &r.inverse_transpose);
            for ij in 0..9 {
                tangent[ij][column] =
                    dp[ij] - r.piola[3 * (ij / 3) + col] * r.inverse_transpose[3 * row + ij % 3];
            }
        }
        if tangent.iter().flatten().any(|x| !x.is_finite()) {
            return Err(FinitePorousError::Numerical(
                "consistent Piola tangent is unrepresentable",
            ));
        }
        let equivalent_kirchhoff_pa = self.equivalent(&r.pressure.stress);
        if !equivalent_kirchhoff_pa.is_finite() {
            return Err(FinitePorousError::Numerical(
                "returned equivalent stress is unrepresentable",
            ));
        }
        Ok(FinitePorousResponse {
            state: r.state,
            first_piola_pa: r.piola,
            kirchhoff_pa: r.tau,
            tangent_pa: tangent,
            incremental_potential_j_m3: r.potential,
            stored_energy_j_m3: r.stored,
            plastic_dissipation_j_m3: r.heat,
            plastic_increment: r.pressure.multiplier,
            equivalent_kirchhoff_pa,
            yield_residual_pa: r.pressure.yield_residual,
            jacobian: r.trial.jacobian,
            plastic_jacobian: r.plastic_jacobian,
        })
    }

    /// Same return/history admission, without assembling the 9-by-9 Piola tangent.
    pub fn stress_and_potential(
        &self,
        f: &[f64; 9],
        old: &FinitePorousState,
    ) -> Result<([f64; 9], f64), FinitePorousError> {
        let r = self.return_map(f, old)?;
        Ok((r.piola, r.potential))
    }

    /// Borrow an immutable accepted history for the native finite-element solver.
    /// Every Newton trial sees this same history. Commit the checked candidate
    /// only after all elements and the surrounding physical step are accepted.
    pub fn frozen<'a>(
        &'a self,
        state: &'a FinitePorousState,
    ) -> Result<FrozenPorousTrial<'a>, FinitePorousError> {
        self.validate_history(state)?;
        Ok(FrozenPorousTrial {
            law: self,
            accepted: state,
        })
    }
}

/// Native assembly adapter. Its potential is incremental, not stored energy.
/// No state update occurs during residual, tangent or line-search evaluations.
#[derive(Clone, Copy, Debug)]
pub struct FrozenPorousTrial<'a> {
    law: &'a FinitePorous,
    accepted: &'a FinitePorousState,
}
impl FiniteStrainPotential for FrozenPorousTrial<'_> {
    fn energy_density(&self, f: &[f64; 9]) -> f64 {
        self.law
            .stress_and_potential(f, self.accepted)
            .map_or(f64::NAN, |(_, w)| w)
    }
    fn first_piola(&self, f: &[f64; 9]) -> Result<[f64; 9], MaterialError> {
        self.law
            .stress_and_potential(f, self.accepted)
            .map(|(p, _)| p)
            .map_err(|e| MaterialError::State {
                what: e.to_string(),
            })
    }
    fn consistent_tangent(&self, f: &[f64; 9]) -> Result<[[f64; 9]; 9], MaterialError> {
        self.law
            .response(f, self.accepted)
            .map(|r| r.tangent_pa)
            .map_err(|e| MaterialError::State {
                what: e.to_string(),
            })
    }
    fn stress_and_potential(&self, f: &[f64; 9]) -> Result<([f64; 9], f64), MaterialError> {
        self.law
            .stress_and_potential(f, self.accepted)
            .map_err(|e| MaterialError::State {
                what: e.to_string(),
            })
    }
}
