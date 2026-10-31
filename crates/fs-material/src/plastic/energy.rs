//! Incremental work of the existing J2 return map, per reference volume.
//!
//! This is the explicitly ENERGETIC linear-hardening interpretation:
//! psi = 1/2 e:C:e + 1/2 H_iso alpha^2 + H_k dev(eps_p):dev(eps_p)/3.
//! The separately selected Prager term is recoverable backstress storage.
//! The yield offset dissipates
//! sigma_y delta_alpha; selected linear overstress adds eta delta_alpha^2/dt
//! to physical dissipation. Hardening energy is stored. An empirical hardening
//! curve alone does not establish this heat partition for a real specimen.
//!
//! Endpoint stress work closes as delta_psi + physical_dissipation +
//! backward_Euler_loss. The last term is numerical loss, never a heat source.
//! Elastic coefficients and both hardening moduli must be identical across this increment.
//! The yield offset may change without changing the stored-energy function.
//! No radial return, flow rule or tangent is duplicated here.

use core::fmt;
use super::{J2Plasticity, J2State};
use crate::{SmallStrainLaw, Voigt};
use crate::tensor::{contract, trace};
use fs_math::det;

/// Independent absolute numerical gates; neither certifies material accuracy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnergyTolerances {
    /// Maximum absolute work-account residual [J/m3].
    pub energy_j_m3: f64,
    /// Maximum absolute endpoint flow residual [Pa]. For linear overstress this
    /// is `q - sigma_y - H alpha - eta delta_alpha/dt`, not zero overstress.
    pub yield_pa: f64,
}

/// Stored energy of the declared energetic-hardening model [J/m3].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StoredEnergy {
    /// Recoverable quadratic elastic energy density.
    pub elastic_j_m3: f64,
    /// Stored energy density of the explicitly energetic hardening law.
    pub hardening_j_m3: f64,
    /// Recoverable Prager storage H_k dev(eps_p):dev(eps_p)/3, never plastic heat.
    pub kinematic_hardening_j_m3: f64,
}
impl StoredEnergy {
    /// Sum of recoverable elastic and hardening energy density.
    #[must_use]
    pub fn total_j_m3(self) -> f64 {
        let original = self.elastic_j_m3 + self.hardening_j_m3;
        if self.kinematic_hardening_j_m3 == 0.0 { original }
        else { original + self.kinematic_hardening_j_m3 }
    }
}

/// Actual endpoint work and separately identified losses [J/m3].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IncrementEnergy {
    /// Stored energy at the old committed strain and history.
    pub initial: StoredEnergy,
    /// Stored energy at the candidate strain and history.
    pub final_state: StoredEnergy,
    /// Endpoint stress contracted with the total strain increment.
    pub mechanical_work_j_m3: f64,
    /// Yield-offset and, when selected, linear-overstress dissipation.
    pub physical_dissipation_j_m3: f64,
    /// 1/2 delta_e:C:delta_e + 1/2 H_iso (delta_alpha)^2
    /// + H_k dev(delta_eps_p):dev(delta_eps_p)/3.
    pub numerical_dissipation_j_m3: f64,
    /// work - delta_stored - physical_dissipation - numerical_dissipation.
    pub residual_j_m3: f64,
}
impl IncrementEnergy {
    /// Stored-energy change, subtracting the two storage components separately.
    #[must_use]
    pub fn stored_change_j_m3(self) -> f64 {
        let original = (self.final_state.elastic_j_m3 - self.initial.elastic_j_m3)
            + (self.final_state.hardening_j_m3 - self.initial.hardening_j_m3);
        let kinematic = self.final_state.kinematic_hardening_j_m3
            - self.initial.kinematic_hardening_j_m3;
        if kinematic == 0.0 { original } else { original + kinematic }
    }
}

/// Native stress and committed internal variables, with their work account.
#[derive(Clone, Debug, PartialEq)]
pub struct EnergyIncrement {
    /// Candidate total strain in the native tensor-Voigt convention.
    pub strain: Voigt,
    /// Native endpoint stress [Pa].
    pub stress_pa: Voigt,
    /// Native committed plastic variables for this candidate.
    pub state: J2State,
    /// Incremental mechanical work and storage account per reference volume.
    pub energy: IncrementEnergy,
}

/// Explicit model, arithmetic or endpoint consistency refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum EnergyError {
    /// Invalid coefficients, strain, history or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The native endpoint missed the independently supplied yield allowance.
    YieldResidual {
        /// Signed yield-function residual [Pa].
        residual_pa: f64,
        /// Caller-supplied absolute yield allowance [Pa].
        tolerance_pa: f64,
    },
    /// The incremental work identity missed its independently supplied allowance.
    EnergyResidual {
        /// Signed work-account residual [J/m3].
        residual_j_m3: f64,
        /// Caller-supplied absolute work allowance [J/m3].
        tolerance_j_m3: f64,
    },
}
impl fmt::Display for EnergyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "J2 energy: {s}"),
            Self::YieldResidual { residual_pa, tolerance_pa } =>
                write!(f, "J2 endpoint yield residual {residual_pa} Pa exceeds {tolerance_pa} Pa"),
            Self::EnergyResidual { residual_j_m3, tolerance_j_m3 } =>
                write!(f, "J2 work residual {residual_j_m3} J/m3 exceeds {tolerance_j_m3} J/m3"),
        }
    }
}
impl std::error::Error for EnergyError {}

impl J2Plasticity {
    /// Admit finite small-strain history without imposing instantaneous yield.
    /// A viscous plastic solid may retain positive overstress between steps.
    pub fn validate_overstress_history(&self, strain: &Voigt, state: &J2State)
        -> Result<(), EnergyError> {
        self.validate_energy_law()?;
        let elastic = self.admit_energy_history(strain, state)?;
        if self.elastic.stress(&elastic, &()).iter().any(|x| !x.is_finite()) {
            return Err(EnergyError::Invalid("committed overstress is unrepresentable"));
        }
        Ok(())
    }
    /// Admit a committed history at THIS law's temperature and yield offset.
    /// Initial histories need this check; a thermal trial may lower the yield
    /// offset and intentionally make the previously committed state overstress.
    pub fn validate_committed_history(
        &self, strain: &Voigt, state: &J2State, yield_tolerance_pa: f64,
    ) -> Result<(), EnergyError> {
        self.validate_energy_law()?;
        if !positive(yield_tolerance_pa) {
            return Err(EnergyError::Invalid("yield allowance must be finite and positive"));
        }
        let e = self.admit_energy_history(strain, state)?;
        let stress = self.elastic.stress(&e, &());
        let residual = finite(self.yield_function_with_state(&stress, state))?;
        if residual > yield_tolerance_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: residual, tolerance_pa: yield_tolerance_pa,
            });
        }
        Ok(())
    }

    /// Propose a native return-map increment and audit its endpoint work.
    ///
    /// The old history must have been committed with the SAME elastic and
    /// hardening energies. It need not be admissible at this endpoint's yield
    /// offset: temperature-induced yielding at fixed total strain is supported.
    /// No accepted strain or history is mutated on success or refusal.
    pub fn advance_with_energy(
        &self, old_strain: &Voigt, old_state: &J2State, strain: &Voigt,
        tolerances: EnergyTolerances,
    ) -> Result<EnergyIncrement, EnergyError> {
        self.advance_with_overstress_energy(old_strain, old_state, strain, tolerances, 0.0)
    }

    pub(super) fn advance_with_overstress_energy(
        &self, old_strain: &Voigt, old_state: &J2State, strain: &Voigt,
        tolerances: EnergyTolerances, overstress_modulus: f64,
    ) -> Result<EnergyIncrement, EnergyError> {
        self.validate_energy_law()?;
        if !positive(tolerances.energy_j_m3) || !positive(tolerances.yield_pa) {
            return Err(EnergyError::Invalid("energy and yield allowances must be finite and positive"));
        }
        let old_e = self.admit_energy_history(old_strain, old_state)?;
        self.admit_energy_strain(strain)?;
        // The constitutive owner alone performs the predictor and radial return.
        let stress_pa = self.radial_return(strain, old_state, overstress_modulus).stress;
        let state = self.update_state_with_overstress(strain, old_state, overstress_modulus);
        let new_e = self.admit_energy_history(strain, &state)?;
        if stress_pa.iter().any(|v| !v.is_finite()) {
            return Err(EnergyError::Invalid("native stress is not representable"));
        }
        let delta_alpha = finite(state.alpha - old_state.alpha)?;
        if delta_alpha < 0.0
            || delta_alpha == 0.0 && state.plastic_strain != old_state.plastic_strain {
            return Err(EnergyError::Invalid("plastic history increment is not representable or monotone"));
        }
        let yield_residual = finite(self.yield_function_with_state(&stress_pa, &state)
            - product(overstress_modulus, delta_alpha)?)?;
        if yield_residual > tolerances.yield_pa
            || delta_alpha > 0.0 && yield_residual.abs() > tolerances.yield_pa {
            return Err(EnergyError::YieldResidual {
                residual_pa: yield_residual, tolerance_pa: tolerances.yield_pa,
            });
        }
        // Also admit the committed elastic stress, not just the corrected
        // predictor: storing an inconsistent internal state must refuse.
        if overstress_modulus == 0.0 {
            self.validate_committed_history(strain, &state, tolerances.yield_pa)?;
        } else {
            self.validate_overstress_history(strain, &state)?;
        }
        let initial = self.energy_of(&old_e, old_state.alpha, &old_state.plastic_strain)?;
        let final_state = self.energy_of(&new_e, state.alpha, &state.plastic_strain)?;
        let delta_e = subtract(&new_e, &old_e);
        let delta_strain = subtract(strain, old_strain);
        let delta_plastic = subtract(&state.plastic_strain, &old_state.plastic_strain);
        let numerical = self.energy_of(&delta_e, delta_alpha, &delta_plastic)?;
        let mut energy = IncrementEnergy {
            initial, final_state,
            mechanical_work_j_m3: finite(contract(&stress_pa, &delta_strain))?,
            physical_dissipation_j_m3: finite(product(self.yield_stress, delta_alpha)?
                + product(overstress_modulus, product(delta_alpha, delta_alpha)?)?)?,
            numerical_dissipation_j_m3: finite(numerical.total_j_m3())?,
            residual_j_m3: 0.0,
        };
        energy.residual_j_m3 = finite(energy.mechanical_work_j_m3
            - energy.stored_change_j_m3() - energy.physical_dissipation_j_m3
            - energy.numerical_dissipation_j_m3)?;
        if energy.residual_j_m3.abs() > tolerances.energy_j_m3 {
            return Err(EnergyError::EnergyResidual {
                residual_j_m3: energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            });
        }
        Ok(EnergyIncrement { strain: *strain, stress_pa, state, energy })
    }

    fn validate_energy_law(&self) -> Result<(), EnergyError> {
        let (lambda, mu) = self.elastic.lame();
        let bulk = lambda + 2.0 * mu / 3.0;
        if !positive(self.elastic.youngs)
            || !(self.elastic.poisson.is_finite()
                && self.elastic.poisson > -1.0 && self.elastic.poisson < 0.5)
            || !positive(self.elastic.strain_limit)
            || !positive(self.yield_stress)
            || !(self.hardening.is_finite() && self.hardening >= 0.0)
            || !(self.kinematic_hardening.is_finite() && self.kinematic_hardening >= 0.0)
            || !positive(mu) || !positive(bulk)
            || !positive(3.0 * mu + self.combined_hardening())
            || self.elastic.stiffness().iter().flatten().any(|x| !x.is_finite()) {
            return Err(EnergyError::Invalid("energetic coefficients are not admissible or representable"));
        }
        Ok(())
    }

    pub(super) fn admit_energy_strain(&self, strain: &Voigt) -> Result<(), EnergyError> {
        if strain.iter().any(|x| !x.is_finite())
            || norm(strain)? > self.elastic.strain_limit.min(0.05) {
            return Err(EnergyError::Invalid("tensor strain exceeds the declared small-strain bound"));
        }
        Ok(())
    }

    fn admit_energy_history(&self, strain: &Voigt, state: &J2State)
        -> Result<Voigt, EnergyError> {
        self.admit_energy_strain(strain)?;
        self.validate_plastic_history(state)?;
        let e = subtract(strain, &state.plastic_strain);
        self.admit_energy_strain(&e)?;
        Ok(e)
    }

    pub(super) fn validate_plastic_history(&self, state: &J2State)
        -> Result<(), EnergyError> {
        if !state.alpha.is_finite() || state.alpha < 0.0
            || state.plastic_strain.iter().any(|x| !x.is_finite()) {
            return Err(EnergyError::Invalid("plastic history must be finite with nonnegative accumulated strain"));
        }
        let plastic_norm = norm(&state.plastic_strain)?;
        let scale = state.alpha.max(plastic_norm).max(f64::MIN_POSITIVE);
        let roundoff = 64.0 * f64::EPSILON * scale;
        if trace(&state.plastic_strain).abs() > roundoff
            || det::sqrt(2.0 / 3.0) * plastic_norm > state.alpha + roundoff {
            return Err(EnergyError::Invalid("plastic strain violates the deviatoric or accumulated-history bound"));
        }
        Ok(())
    }

    pub(super) fn energy_of(&self, elastic_strain: &Voigt, alpha: f64, plastic_strain: &Voigt)
        -> Result<StoredEnergy, EnergyError> {
        let stress = self.elastic.stress(elastic_strain, &());
        let elastic_j_m3 = product(0.5, finite(contract(elastic_strain, &stress))?)?;
        let hardening_j_m3 = product(0.5, product(product(self.hardening, alpha)?, alpha)?)?;
        let kinematic_hardening_j_m3 = if self.kinematic_hardening == 0.0 { 0.0 } else {
            let plastic = crate::tensor::deviator(plastic_strain);
            product(1.0 / 3.0, product(self.kinematic_hardening, finite(contract(&plastic, &plastic))?)?)?
        };
        if elastic_j_m3 < 0.0 || hardening_j_m3 < 0.0
            || kinematic_hardening_j_m3 < 0.0
            || !(elastic_j_m3 + hardening_j_m3 + kinematic_hardening_j_m3).is_finite() {
            return Err(EnergyError::Invalid("stored elastic or hardening energy is not admissible"));
        }
        Ok(StoredEnergy { elastic_j_m3, hardening_j_m3, kinematic_hardening_j_m3 })
    }
}

fn subtract(a: &Voigt, b: &Voigt) -> Voigt {
    core::array::from_fn(|i| a[i] - b[i])
}
fn norm(a: &Voigt) -> Result<f64, EnergyError> {
    finite(det::sqrt(finite(contract(a, a))?))
}
fn positive(x: f64) -> bool { x.is_finite() && x > 0.0 }
fn finite(x: f64) -> Result<f64, EnergyError> {
    if x.is_finite() { Ok(x) } else { Err(EnergyError::Invalid("nonfinite energy arithmetic")) }
}
fn product(a: f64, b: f64) -> Result<f64, EnergyError> {
    let x = finite(a * b)?;
    if a != 0.0 && b != 0.0 && x == 0.0 {
        return Err(EnergyError::Invalid("nonzero energy underflowed"));
    }
    Ok(x)
}
