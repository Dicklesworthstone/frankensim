//! Permanent, sub-yield, isotropic Norton creep with an Arrhenius temperature
//! shift. This is the secondary-creep law, independent of plastic yield:
//! `alpha_dot = rate_ref (q/stress_ref)^n exp[Q/R (1/T_ref - 1/T)]`,
//! `epsilon_c_dot = (3/2) dev(sigma)/q alpha_dot`.
//!
//! Backward Euler reduces to the strictly monotone scalar equation
//! `q + 3 G dt alpha_dot(q) = q_trial`. Its consistent tangent uses the
//! derivative of this discrete solve. The native Voigt convention stores tensor
//! shear; contraction doubles its shear entries. Trial evaluations never advance
//! accepted history. Use [`NortonCreep::advance_with_energy`] to admit a candidate
//! before committing it. The infallible legacy trait returns NaNs on refusal.
//!
//! Reference formulations: MOOSE's PowerLawCreepStressUpdate and COMSOL's
//! Creep and Viscoplasticity theory (Norton, Arrhenius and backward Euler):
//! <https://mooseframework.inl.gov/source/materials/PowerLawCreepStressUpdate.html>
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.sme/sme_ug_theory.06.039.html>.
//! Parameters are authored model inputs until separately bound to source data.

/// Simultaneous secondary creep and energetic isotropic J2 plasticity.
pub mod plastic;
/// Immutable material-card binding for the same permanent creep law.
pub mod source;
/// Source-bound uniform physical-heat feedback with explicit zero expansion.
pub mod thermal;
/// Free-sided axial restriction for the existing fiber-section material contract.
pub mod uniaxial;

use crate::elastic::IsotropicElastic;
use crate::tensor::{contract, trace};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use core::fmt;
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

/// Molar gas constant in coherent SI units, J/(mol K).
pub const GAS_CONSTANT_J_MOL_K: f64 = 8.314_462_618_153_24;

/// Coefficients of an explicitly normalized secondary-creep model.
#[derive(Clone, Debug, PartialEq)]
pub struct NortonCreepParameters {
    /// Positive equivalent reference stress [Pa].
    pub reference_stress_pa: f64,
    /// Positive equivalent creep rate at reference stress and temperature [1/s].
    pub reference_rate_s: f64,
    /// Dimensionless stress exponent, at least one.
    pub stress_exponent: f64,
    /// Nonnegative molar activation energy [J/mol].
    pub activation_energy_j_mol: f64,
    /// Absolute temperature at which the reference rate is specified [K].
    pub reference_temperature_k: f64,
    /// Closed admitted absolute-temperature interval [K].
    pub temperature_range_k: [f64; 2],
}

/// Per-material-point permanent history, separate from plastic hardening.
#[derive(Clone, Debug, PartialEq)]
pub struct CreepState {
    /// Signed, traceless creep strain in tensor-Voigt coordinates.
    pub creep_strain: Voigt,
    /// Nondecreasing path length, integral of the equivalent creep rate.
    pub equivalent_creep_strain: f64,
}

/// Absolute numerical admission tolerances, not material-accuracy claims.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepTolerances {
    /// Endpoint flow and committed stress consistency [Pa].
    pub stress_pa: f64,
    /// Incremental work-balance residual [J/m3].
    pub energy_j_m3: f64,
}

/// Per-reference-volume work account [J/m3]. Elastic/hardening free energy,
/// reversible thermal storage/transfer and irreversible losses remain distinct.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CreepEnergy {
    pub initial_stored_j_m3: f64,
    pub final_stored_j_m3: f64,
    /// Endpoint stress contracted with the actual total-strain increment.
    pub mechanical_work_j_m3: f64,
    /// Sum of the separate creep and plastic physical dissipation accounts.
    pub physical_dissipation_j_m3: f64,
    /// Endpoint `q delta_alpha_c`, available as physical creep heat.
    pub creep_dissipation_j_m3: f64,
    /// Yield-offset work `sigma_y delta_alpha_p`, zero for creep alone.
    pub plastic_dissipation_j_m3: f64,
    /// Backward-Euler elastic and, when selected, hardening loss; never heat.
    pub numerical_dissipation_j_m3: f64,
    /// Change in T times mechanical entropy; zero for isothermal storage laws.
    pub thermoelastic_storage_change_j_m3: f64,
    /// Signed reversible transfer into the caloric chart; never dissipation.
    pub reversible_heat_j_m3: f64,
    /// Work minus free/internal storage changes, heat transfer and numerical loss.
    pub residual_j_m3: f64,
}

impl CreepEnergy {
    /// Physical creep/plastic heat plus signed reciprocal thermal exchange.
    pub fn caloric_heat_j_m3(&self) -> f64 {
        self.physical_dissipation_j_m3 + self.reversible_heat_j_m3
    }
}

/// Admitted candidate; the caller atomically commits state after its solve.
#[derive(Clone, Debug, PartialEq)]
pub struct CreepIncrement<S = CreepState> {
    pub strain: Voigt,
    pub stress_pa: Voigt,
    pub tangent_pa: Tangent6,
    pub state: S,
    pub energy: CreepEnergy,
}

/// A physical creep interval on the original constitutive history owner.
/// Consumers freeze this history through all trial evaluations, then audit one
/// increment before accepting it. A combined law retains its plastic history.
pub trait CreepLaw: SmallStrainLaw {
    fn kinetics(&self) -> &NortonCreep;
    fn validate_creep_history(&self, strain: &Voigt, state: &Self::State) -> Result<()>;
    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &Self::State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<Self::State>>;
}

/// Parameter/domain, representability or independently measured residual failure.
#[derive(Clone, Debug, PartialEq)]
pub enum CreepError {
    Invalid(&'static str),
    /// Original source refusal while resolving a thermal history endpoint.
    Source(crate::state_point::MaterialStatePointError),
    /// Original elastic-potential coefficient or state refusal.
    Elastic(crate::MaterialError),
    /// Native energetic J2 coefficient/history/storage refusal.
    Plastic(crate::plastic::energy::EnergyError),
    StressResidual {
        residual_pa: f64,
        tolerance_pa: f64,
    },
    EnergyResidual {
        residual_j_m3: f64,
        tolerance_j_m3: f64,
    },
}
impl fmt::Display for CreepError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(what) => write!(f, "Norton creep: {what}"),
            Self::Plastic(error) => error.fmt(f),
            Self::Source(error) => error.fmt(f),
            Self::Elastic(error) => error.fmt(f),
            Self::StressResidual {
                residual_pa,
                tolerance_pa,
            } => write!(
                f,
                "Norton creep stress residual {residual_pa} Pa exceeds {tolerance_pa} Pa"
            ),
            Self::EnergyResidual {
                residual_j_m3,
                tolerance_j_m3,
            } => write!(
                f,
                "Norton creep work residual {residual_j_m3} J/m3 exceeds {tolerance_j_m3} J/m3"
            ),
        }
    }
}
impl std::error::Error for CreepError {}
impl From<crate::plastic::energy::EnergyError> for CreepError {
    fn from(error: crate::plastic::energy::EnergyError) -> Self {
        Self::Plastic(error)
    }
}
type Result<T> = core::result::Result<T, CreepError>;
fn require(ok: bool, what: &'static str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(CreepError::Invalid(what))
    }
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn difference(a: &Voigt, b: &Voigt) -> Voigt {
    core::array::from_fn(|i| a[i] - b[i])
}

// Scale first: squaring representable large stresses must not overflow early.
fn norm(v: &Voigt) -> f64 {
    let scale = v.iter().fold(0.0_f64, |m, x| m.max(x.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    let unit = core::array::from_fn(|i| v[i] / scale);
    scale * det::sqrt(contract(&unit, &unit))
}
fn logistic_inverse(log_x: f64) -> f64 {
    if log_x >= 0.0 {
        let e = det::exp(-log_x);
        e / (1.0 + e)
    } else {
        1.0 / (1.0 + det::exp(log_x))
    }
}

/// One explicit physical step at fixed temperature and elastic coefficients.
/// Rebuild with the same parameters for each actual duration/temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct NortonCreep {
    elastic: IsotropicElastic,
    parameters: NortonCreepParameters,
    temperature_k: f64,
    duration_s: f64,
    log_rate_at_reference_stress: f64,
    log_reference_stress: f64,
    log_duration: f64,
    log_step_scale: f64,
    shear_pa: f64,
    bulk_pa: f64,
}
struct Update {
    stress: Voigt,
    tangent: Tangent6,
    state: CreepState,
    q_trial: f64,
    q: f64,
    delta_alpha: f64,
}

struct Trial {
    pressure: f64,
    deviator: Voigt,
    q: f64,
}

// One scalar kinetic owner for tensor relaxation (3G) and its free-sided axial
// restriction (E). The dimensional stiffness is supplied through its prepared
// logarithmic step scale, never through a second constitutive flow rule.
struct Relaxation {
    q: f64,
    delta_alpha: f64,
    ratio: f64,
    radial: f64,
}

impl NortonCreep {
    pub fn new(
        elastic: IsotropicElastic,
        parameters: NortonCreepParameters,
        temperature_k: f64,
        duration_s: f64,
    ) -> Result<Self> {
        let p = &parameters;
        let [lo, hi] = p.temperature_range_k;
        require(
            positive(p.reference_stress_pa)
                && positive(p.reference_rate_s)
                && p.stress_exponent.is_finite()
                && p.stress_exponent >= 1.0
                && p.activation_energy_j_mol.is_finite()
                && p.activation_energy_j_mol >= 0.0,
            "reference stress/rate must be positive, exponent >= 1 and activation energy >= 0",
        )?;
        require(
            positive(lo)
                && hi.is_finite()
                && hi >= lo
                && p.reference_temperature_k.is_finite()
                && p.reference_temperature_k >= lo
                && p.reference_temperature_k <= hi
                && temperature_k.is_finite()
                && temperature_k >= lo
                && temperature_k <= hi
                && positive(duration_s),
            "temperature/reference temperature must lie in the positive domain; dt must be positive",
        )?;
        require(
            positive(elastic.youngs)
                && elastic.poisson.is_finite()
                && elastic.poisson > -1.0
                && elastic.poisson < 0.5
                && positive(elastic.strain_limit)
                && elastic.strain_limit <= 0.05,
            "elasticity requires positive E, -1 < nu < 0.5 and a strain ceiling <= 0.05",
        )?;
        let (lambda, shear_pa) = elastic.lame();
        let bulk_pa = lambda + 2.0 * shear_pa / 3.0;
        require(
            positive(shear_pa)
                && positive(3.0 * shear_pa)
                && positive(bulk_pa)
                && elastic.stiffness().iter().flatten().all(|x| x.is_finite()),
            "elastic stiffness is not representable",
        )?;
        let shift = (p.activation_energy_j_mol / GAS_CONSTANT_J_MOL_K)
            * ((temperature_k - p.reference_temperature_k)
                / p.reference_temperature_k
                / temperature_k);
        let log_rate_at_reference_stress = det::ln(p.reference_rate_s) + shift;
        require(
            log_rate_at_reference_stress.is_finite(),
            "Arrhenius shift is not representable",
        )?;
        let log_reference_stress = det::ln(p.reference_stress_pa);
        let log_duration = det::ln(duration_s);
        let log_step_scale = det::ln(3.0 * shear_pa) + log_duration + log_rate_at_reference_stress
            - log_reference_stress;
        Ok(Self {
            elastic,
            parameters,
            temperature_k,
            duration_s,
            log_rate_at_reference_stress,
            log_reference_stress,
            log_duration,
            log_step_scale,
            shear_pa,
            bulk_pa,
        })
    }
    pub const fn elastic(&self) -> &IsotropicElastic {
        &self.elastic
    }
    pub const fn parameters(&self) -> &NortonCreepParameters {
        &self.parameters
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    fn log_rate(&self, q: f64) -> f64 {
        self.log_rate_at_reference_stress
            + self.parameters.stress_exponent * (det::ln(q) - self.log_reference_stress)
    }
    /// Evaluate equivalent rate at a nonnegative stress inside the selected law.
    pub fn equivalent_creep_rate_s(&self, q: f64) -> Result<f64> {
        require(
            q.is_finite() && q >= 0.0,
            "equivalent stress must be finite and nonnegative",
        )?;
        if q == 0.0 {
            return Ok(0.0);
        }
        let rate = det::exp(self.log_rate(q));
        require(positive(rate), "creep rate is not representable")?;
        Ok(rate)
    }

    /// Admit native tensor history. Total, creep and elastic Frobenius strain
    /// norms stay within the same small-strain ceiling; accumulated path length
    /// is finite and at least the equivalent norm of the signed creep strain.
    pub fn validate_history(&self, strain: &Voigt, state: &CreepState) -> Result<()> {
        self.validate_strain(strain)?;
        self.validate_creep_state(state)?;
        self.validate_strain(&difference(strain, &state.creep_strain))
    }

    fn validate_strain(&self, strain: &Voigt) -> Result<()> {
        let magnitude = norm(strain);
        require(
            strain.iter().all(|x| x.is_finite())
                && magnitude.is_finite()
                && magnitude <= self.elastic.strain_limit,
            "total, elastic or inelastic strain exceeds the small-strain domain",
        )
    }

    fn validate_creep_state(&self, state: &CreepState) -> Result<()> {
        require(
            state.creep_strain.iter().all(|x| x.is_finite())
                && state.equivalent_creep_strain.is_finite()
                && state.equivalent_creep_strain >= 0.0,
            "strain and creep history must be finite and accumulated creep nonnegative",
        )?;
        let creep_norm = norm(&state.creep_strain);
        let allowance = 256.0 * f64::EPSILON * creep_norm;
        require(
            trace(&state.creep_strain).abs() <= allowance
                && det::sqrt(2.0 / 3.0) * creep_norm <= state.equivalent_creep_strain + allowance,
            "creep history must be traceless and consistent with its accumulated path length",
        )?;
        self.validate_strain(&state.creep_strain)
    }

    // q/q_trial = exp(y), exp(y)+exp(log_k+n*y)=1. Log-space evaluation avoids
    // ever forming q^n or a potentially overflowing kinetic prefactor. The
    // bracket has width ln(2); safeguarded Newton usually takes a few steps.
    fn stress_ratio(&self, q_trial: f64, log_step_scale: f64) -> Result<f64> {
        let n = self.parameters.stress_exponent;
        let log_k = log_step_scale + (n - 1.0) * (det::ln(q_trial) - self.log_reference_stress);
        require(
            log_k.is_finite(),
            "dimensionless creep step is not representable",
        )?;
        if n == 1.0 {
            return Ok(logistic_inverse(log_k));
        }
        let mut upper = if log_k > 0.0 { -log_k / n } else { 0.0 };
        let mut lower = upper - core::f64::consts::LN_2;
        let mut y = upper;
        for _ in 0..64 {
            let second = log_k + n * y;
            let shift = y.max(second);
            let first_weight = det::exp(y - shift);
            let second_weight = det::exp(second - shift);
            let sum = first_weight + second_weight;
            let residual = shift + det::ln(sum);
            let roundoff = 8.0 * f64::EPSILON * (1.0 + log_k.abs() + (n * y).abs());
            if residual.abs() <= roundoff {
                let ratio = det::exp(y);
                require(
                    positive(ratio) && ratio <= 1.0,
                    "relaxed stress is not representable",
                )?;
                return Ok(ratio);
            }
            if residual > 0.0 {
                upper = y;
            } else {
                lower = y;
            }
            let slope = (first_weight + n * second_weight) / sum;
            let candidate = y - residual / slope;
            y = if candidate.is_finite() && candidate > lower && candidate < upper {
                candidate
            } else {
                0.5 * lower + 0.5 * upper
            };
        }
        Err(CreepError::Invalid("bounded creep root did not converge"))
    }

    fn relaxation(&self, q_trial: f64, log_step_scale: f64) -> Result<Relaxation> {
        let n = self.parameters.stress_exponent;
        let ratio = if q_trial > 0.0 {
            self.stress_ratio(q_trial, log_step_scale)?
        } else if n == 1.0 {
            logistic_inverse(log_step_scale)
        } else {
            1.0
        };
        let q = q_trial * ratio;
        let delta_alpha = if q_trial == 0.0 {
            0.0
        } else {
            det::exp(self.log_duration + self.log_rate(q))
        };
        require(
            q_trial == 0.0 || positive(q) && delta_alpha.is_finite() && delta_alpha >= 0.0,
            "creep increment or relaxed stress is not representable",
        )?;
        let radial = if q_trial == 0.0 {
            ratio
        } else {
            ratio / (ratio + n * (1.0 - ratio))
        };
        Ok(Relaxation {
            q,
            delta_alpha,
            ratio,
            radial,
        })
    }

    fn predictor(&self, elastic_strain: &Voigt) -> Result<Trial> {
        let pressure = self.bulk_pa * trace(elastic_strain);
        let two_g = 2.0 * self.shear_pa;
        // Form the isotropic split from strain differences. A generic elastic
        // matvec can round equal normal stresses differently; subtracting those
        // rounded stresses would turn exactly hydrostatic strain into false creep.
        let sx = two_g
            * ((elastic_strain[0] - elastic_strain[1]) / 3.0
                + (elastic_strain[0] - elastic_strain[2]) / 3.0);
        let sy = two_g
            * ((elastic_strain[1] - elastic_strain[0]) / 3.0
                + (elastic_strain[1] - elastic_strain[2]) / 3.0);
        let s = [
            sx,
            sy,
            -(sx + sy),
            two_g * elastic_strain[3],
            two_g * elastic_strain[4],
            two_g * elastic_strain[5],
        ];
        require(
            pressure.is_finite() && s.iter().all(|x| x.is_finite()),
            "trial stress is not representable",
        )?;
        let s_norm = norm(&s);
        let q_trial = det::sqrt(1.5) * s_norm;
        require(
            q_trial.is_finite(),
            "equivalent trial stress is not representable",
        )?;
        Ok(Trial {
            pressure,
            deviator: s,
            q: q_trial,
        })
    }

    fn radial_response(&self, trial: &Trial, ratio: f64, radial: f64) -> Result<(Voigt, Tangent6)> {
        let s = trial.deviator;
        let stress =
            core::array::from_fn(|i| ratio * s[i] + if i < 3 { trial.pressure } else { 0.0 });
        let s_norm = norm(&s);
        let direction: Voigt = if s_norm > 0.0 {
            core::array::from_fn(|i| s[i] / s_norm)
        } else {
            [0.0; 6]
        };
        let tangent = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                let identity = if i == j { 1.0 } else { 0.0 };
                let volumetric = if i < 3 && j < 3 { 1.0 } else { 0.0 };
                self.bulk_pa * volumetric
                    + 2.0
                        * self.shear_pa
                        * (ratio * (identity - volumetric / 3.0)
                            + (radial - ratio)
                                * direction[i]
                                * direction[j]
                                * if j < 3 { 1.0 } else { 2.0 })
            })
        });
        require(
            stress
                .iter()
                .chain(tangent.iter().flatten())
                .all(|x| x.is_finite()),
            "endpoint stress or tangent is not representable",
        )?;
        Ok((stress, tangent))
    }

    fn update(&self, strain: &Voigt, old: &CreepState) -> Result<Update> {
        self.validate_history(strain, old)?;
        let trial = self.predictor(&difference(strain, &old.creep_strain))?;
        let q_trial = trial.q;
        let Relaxation {
            q,
            delta_alpha,
            ratio,
            radial,
        } = self.relaxation(q_trial, self.log_step_scale)?;
        let mut state = old.clone();
        if delta_alpha > 0.0 {
            for i in 0..6 {
                state.creep_strain[i] += (trial.deviator[i] / q_trial) * (1.5 * delta_alpha);
            }
            // Keep the dependent deviatoric component through nearly cancelled
            // load reversals, whose remaining signed history can approach zero.
            state.creep_strain[2] = -(state.creep_strain[0] + state.creep_strain[1]);
            state.equivalent_creep_strain += delta_alpha;
            // Near zero-load equilibrium a positive increment may round away
            // against accepted history. Retain its analytic heat; the checked
            // stress/work budgets audit this rounding. Rejecting every such
            // trial would prevent Newton from completing stress-free unloading.
        }
        self.validate_history(strain, &state)?;
        let (stress, tangent) = self.radial_response(&trial, ratio, radial)?;
        Ok(Update {
            stress,
            tangent,
            state,
            q_trial,
            q,
            delta_alpha,
        })
    }

    fn stored_energy(&self, e: &Voigt) -> Result<f64> {
        let value = 0.5 * contract(&self.elastic.stress(e, &()), e);
        require(
            value.is_finite() && value >= 0.0,
            "stored elastic energy is not representable",
        )?;
        Ok(value)
    }

    /// Evaluate one candidate without changing accepted history. Elasticity is
    /// fixed across this work account; temperature controls the creep kinetics.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old_state: &CreepState,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement> {
        require(
            positive(tolerances.stress_pa) && positive(tolerances.energy_j_m3),
            "stress and energy tolerances must be positive and finite",
        )?;
        self.validate_history(old_strain, old_state)?;
        let next = self.update(strain, old_state)?;
        let old_e = difference(old_strain, &old_state.creep_strain);
        let new_e = difference(strain, &next.state.creep_strain);
        let stored_stress = self.elastic.stress(&new_e, &());
        let mut residual = (next.q + 3.0 * self.shear_pa * next.delta_alpha - next.q_trial).abs();
        for i in 0..6 {
            residual = residual.max((stored_stress[i] - next.stress[i]).abs());
        }
        if !residual.is_finite() || residual > tolerances.stress_pa {
            return Err(CreepError::StressResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.stress_pa,
            });
        }
        let mut energy = CreepEnergy {
            initial_stored_j_m3: self.stored_energy(&old_e)?,
            final_stored_j_m3: self.stored_energy(&new_e)?,
            mechanical_work_j_m3: contract(&next.stress, &difference(strain, old_strain)),
            physical_dissipation_j_m3: next.q * next.delta_alpha,
            creep_dissipation_j_m3: next.q * next.delta_alpha,
            plastic_dissipation_j_m3: 0.0,
            numerical_dissipation_j_m3: self.stored_energy(&difference(&new_e, &old_e))?,
            thermoelastic_storage_change_j_m3: 0.0,
            reversible_heat_j_m3: 0.0,
            residual_j_m3: 0.0,
        };
        energy.residual_j_m3 = energy.mechanical_work_j_m3
            - (energy.final_stored_j_m3 - energy.initial_stored_j_m3)
            - energy.physical_dissipation_j_m3
            - energy.numerical_dissipation_j_m3;
        require(
            energy.mechanical_work_j_m3.is_finite()
                && energy.physical_dissipation_j_m3.is_finite()
                && energy.physical_dissipation_j_m3 >= 0.0,
            "creep work or heat is not representable",
        )?;
        if !energy.residual_j_m3.is_finite() || energy.residual_j_m3.abs() > tolerances.energy_j_m3
        {
            return Err(CreepError::EnergyResidual {
                residual_j_m3: energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            });
        }
        Ok(CreepIncrement {
            strain: *strain,
            stress_pa: next.stress,
            tangent_pa: next.tangent,
            state: next.state,
            energy,
        })
    }
}

impl CreepLaw for NortonCreep {
    fn kinetics(&self) -> &NortonCreep {
        self
    }
    fn validate_creep_history(&self, strain: &Voigt, state: &CreepState) -> Result<()> {
        self.validate_history(strain, state)
    }
    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &CreepState,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement> {
        self.advance_with_energy(old_strain, old, strain, tolerances)
    }
}

impl SmallStrainLaw for NortonCreep {
    type State = CreepState;
    fn initial_state(&self) -> CreepState {
        CreepState {
            creep_strain: [0.0; 6],
            equivalent_creep_strain: 0.0,
        }
    }
    fn stress(&self, strain: &Voigt, state: &CreepState) -> Voigt {
        self.update(strain, state)
            .map(|v| v.stress)
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, strain: &Voigt, state: &CreepState) -> Tangent6 {
        self.update(strain, state)
            .map(|v| v.tangent)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, strain: &Voigt, state: &CreepState) -> CreepState {
        self.update(strain, state)
            .map(|v| v.state)
            .unwrap_or(CreepState {
                creep_strain: [f64::NAN; 6],
                equivalent_creep_strain: f64::NAN,
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "bounded small-strain isotropic secondary creep; no rupture or damage law",
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.norton-arrhenius-creep",
            "0.1.0",
            Ambition::Solid,
            vec![
                "isotropic small-strain J2 secondary creep without a yield threshold".into(),
                "Norton exponent >= 1; Arrhenius shift at absolute temperature".into(),
                "backward Euler with separate physical heat and numerical loss".into(),
            ],
            ValidityDomain::unconstrained()
                .with("strain-magnitude", 0.0, self.elastic.strain_limit)
                .with(
                    "temperature-k",
                    self.parameters.temperature_range_k[0],
                    self.parameters.temperature_range_k[1],
                ),
            vec![
                "no primary/tertiary creep, finite deformation, plasticity or rupture".into(),
                "no measured material qualification without source binding".into(),
            ],
            0.05,
        )
    }
}
