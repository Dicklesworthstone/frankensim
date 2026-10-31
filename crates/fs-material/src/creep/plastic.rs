//! Simultaneous Norton creep and energetic, linear-isotropic J2 plasticity.
//!
//! The elastic strain is `e = epsilon - epsilon_c - epsilon_p`. Both flows are
//! coaxial with the endpoint deviatoric stress, so one radial update satisfies
//! `q + 3G (delta_c + delta_p) = q_trial`, `delta_c = dt rate(q)` and the J2
//! complementarity conditions. The plastic branch reduces exactly to the same
//! bounded logarithmic Norton root with effective stiffness `3G H/(3G+H)`.
//! This is a simultaneous solve, not sequential creep/plastic return maps.
//!
//! Storage and loss use the original native energetic J2 interpretation:
//! `psi = e:C:e/2 + H alpha_p^2/2`; physical heat is
//! `q delta_c + sigma_y delta_p`. The backward-Euler elastic/hardening loss
//! stays separate. Elastic coefficients and H are fixed across the accounting
//! interval; the endpoint yield offset may change, including thermal yielding
//! at held strain. Every trial starts from the same immutable accepted history.
//!
//! Small-strain isotropy only. Prager hardening is refused: its shifted plastic
//! flow is generally not coaxial with Norton flow. No finite-strain, damage,
//! tertiary-creep, rupture or service-life claim follows from this model.

/// Original material-card and claim binding for the combined law.
pub mod source;
/// Exact free-sided axial restriction for the native fiber material contract.
pub mod uniaxial;

use super::{
    CreepEnergy, CreepError, CreepIncrement, CreepLaw, CreepState, CreepTolerances, NortonCreep,
    Relaxation, Result, difference, positive, require,
};
use crate::plastic::{J2Plasticity, J2State};
use crate::tensor::contract;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_math::det;

/// Distinct native histories; neither mechanism resets or absorbs the other.
#[derive(Clone, Debug, PartialEq)]
pub struct NortonJ2State {
    pub creep: CreepState,
    pub plastic: J2State,
}

/// One uncommitted response, with the exact derivative of this coupled update.
#[derive(Clone, Debug, PartialEq)]
pub struct NortonJ2Response {
    pub strain: Voigt,
    pub stress_pa: Voigt,
    pub tangent_pa: Tangent6,
    pub state: NortonJ2State,
}

/// One physical duration, original Norton kinetics and original isotropic J2.
#[derive(Clone, Debug, PartialEq)]
pub struct NortonJ2 {
    creep: NortonCreep,
    plastic: J2Plasticity,
    stiffness: CoupledStiffness,
}

// Prepared dimensional stiffness of the same coupled radial solve: 3G in
// three dimensions, E after exact elimination of free transverse stress.
#[derive(Clone, Debug, PartialEq)]
struct CoupledStiffness {
    modulus_pa: f64,
    log_modulus: f64,
    log_step_scale: f64,
    hardening_weight: f64,
    log_effective_step_scale: f64,
}

impl CoupledStiffness {
    fn new(creep: &NortonCreep, hardening: f64, modulus_pa: f64) -> Self {
        let log_modulus = det::ln(modulus_pa);
        let log_step_scale = log_modulus + creep.log_duration + creep.log_rate_at_reference_stress
            - creep.log_reference_stress;
        let denominator = modulus_pa + hardening;
        // Prepare in logarithms: h*H/(h+H) can under/overflow in direct arithmetic.
        let log_effective_step_scale = if hardening > 0.0 {
            log_step_scale + det::ln(hardening) - det::ln(denominator)
        } else {
            f64::NEG_INFINITY
        };
        Self {
            modulus_pa,
            log_modulus,
            log_step_scale,
            hardening_weight: hardening / denominator,
            log_effective_step_scale,
        }
    }
}

struct CoupledRelaxation {
    q: f64,
    delta_c: f64,
    delta_p: f64,
    ratio: f64,
    radial: f64,
}

struct Update {
    response: NortonJ2Response,
    q_trial: f64,
    relaxed: CoupledRelaxation,
}

impl NortonJ2 {
    /// Admit one common elastic backbone and the explicitly energetic,
    /// isotropic hardening model. The stored native laws remain immutable.
    pub fn new(creep: NortonCreep, plastic: J2Plasticity) -> Result<Self> {
        plastic.validate_energy_law()?;
        require(
            creep.elastic() == &plastic.elastic,
            "creep and plasticity must use the identical elastic backbone",
        )?;
        require(
            plastic.kinematic_hardening == 0.0,
            "simultaneous radial creep/plasticity does not admit kinematic hardening",
        )?;
        let stiffness = CoupledStiffness::new(&creep, plastic.hardening, 3.0 * creep.shear_pa);
        Ok(Self {
            creep,
            plastic,
            stiffness,
        })
    }

    pub const fn creep(&self) -> &NortonCreep {
        &self.creep
    }
    pub const fn plastic(&self) -> &J2Plasticity {
        &self.plastic
    }

    fn elastic_strain(strain: &Voigt, state: &NortonJ2State) -> Voigt {
        core::array::from_fn(|i| {
            strain[i] - state.creep.creep_strain[i] - state.plastic.plastic_strain[i]
        })
    }

    /// Validate each native history separately and their shared elastic strain.
    /// Yield admissibility is checked on the candidate, not this old history:
    /// a changed endpoint yield offset may legitimately cause new plastic flow.
    pub fn validate_history(&self, strain: &Voigt, state: &NortonJ2State) -> Result<()> {
        self.creep.validate_strain(strain)?;
        self.validate_state(state)?;
        self.creep
            .validate_strain(&Self::elastic_strain(strain, state))
    }

    fn validate_state(&self, state: &NortonJ2State) -> Result<()> {
        self.creep.validate_creep_state(&state.creep)?;
        self.plastic.validate_plastic_history(&state.plastic)?;
        self.creep.validate_strain(&state.plastic.plastic_strain)?;
        require(
            positive(self.plastic.yield_stress + self.plastic.hardening * state.plastic.alpha),
            "hardened yield stress is not representable",
        )
    }

    fn relaxation(
        &self,
        q_trial: f64,
        alpha: f64,
        stiffness: &CoupledStiffness,
    ) -> Result<CoupledRelaxation> {
        let h = stiffness.modulus_pa;
        let yield_stress = self.plastic.yield_stress + self.plastic.hardening * alpha;
        // Monotonicity selects the branch with one rate evaluation, avoiding a
        // preliminary creep-only root. Logs avoid an overflowing rate at yield.
        let inactive = q_trial <= yield_stress || {
            let log_rate = self.creep.log_rate(yield_stress);
            let log_flow = stiffness.log_modulus + self.creep.log_duration + log_rate;
            let log_gap = det::ln(q_trial - yield_stress);
            // At sustained yield, log rounding can spuriously activate plasticity
            // and produce a negative roundoff-sized increment. Select the
            // creep-only side of this floating-point tie; the checked candidate
            // still audits actual stress balance, yield and work. No flow is clipped.
            let tie = 8.0
                * f64::EPSILON
                * (stiffness.log_modulus.abs()
                    + self.creep.log_duration.abs()
                    + log_rate.abs()
                    + log_gap.abs());
            log_flow >= log_gap || log_flow.is_finite() && log_gap - log_flow <= tie
        };
        if inactive {
            let Relaxation {
                q,
                delta_alpha,
                ratio,
                radial,
            } = self.creep.relaxation(q_trial, stiffness.log_step_scale)?;
            return Ok(CoupledRelaxation {
                q,
                delta_c: delta_alpha,
                delta_p: 0.0,
                ratio,
                radial,
            });
        }
        let (q, delta_c, radial) = if self.plastic.hardening == 0.0 {
            (
                yield_stress,
                det::exp(self.creep.log_duration + self.creep.log_rate(yield_stress)),
                0.0,
            )
        } else {
            let effective_trial =
                yield_stress + stiffness.hardening_weight * (q_trial - yield_stress);
            let relaxed = self
                .creep
                .relaxation(effective_trial, stiffness.log_effective_step_scale)?;
            (
                relaxed.q,
                relaxed.delta_alpha,
                stiffness.hardening_weight * relaxed.radial,
            )
        };
        let delta_p = ((q_trial - yield_stress) - h * delta_c) / (h + self.plastic.hardening);
        require(
            positive(q)
                && delta_c.is_finite()
                && delta_c >= 0.0
                && delta_p.is_finite()
                && delta_p >= 0.0,
            "coupled creep/plastic increments are not representable and nonnegative",
        )?;
        Ok(CoupledRelaxation {
            q,
            delta_c,
            delta_p,
            ratio: q / q_trial,
            radial,
        })
    }

    fn update(&self, strain: &Voigt, old: &NortonJ2State) -> Result<Update> {
        self.validate_history(strain, old)?;
        let trial = self.creep.predictor(&Self::elastic_strain(strain, old))?;
        let relaxed = self.relaxation(trial.q, old.plastic.alpha, &self.stiffness)?;
        let mut state = old.clone();
        if relaxed.delta_c > 0.0 {
            for i in 0..6 {
                state.creep.creep_strain[i] +=
                    (trial.deviator[i] / trial.q) * (1.5 * relaxed.delta_c);
            }
            state.creep.creep_strain[2] =
                -(state.creep.creep_strain[0] + state.creep.creep_strain[1]);
            state.creep.equivalent_creep_strain += relaxed.delta_c;
        }
        if relaxed.delta_p > 0.0 {
            for i in 0..6 {
                state.plastic.plastic_strain[i] +=
                    (trial.deviator[i] / trial.q) * (1.5 * relaxed.delta_p);
            }
            state.plastic.plastic_strain[2] =
                -(state.plastic.plastic_strain[0] + state.plastic.plastic_strain[1]);
            state.plastic.alpha += relaxed.delta_p;
        }
        self.validate_history(strain, &state)?;
        let (stress_pa, tangent_pa) =
            self.creep
                .radial_response(&trial, relaxed.ratio, relaxed.radial)?;
        Ok(Update {
            response: NortonJ2Response {
                strain: *strain,
                stress_pa,
                tangent_pa,
                state,
            },
            q_trial: trial.q,
            relaxed,
        })
    }

    pub fn response(&self, strain: &Voigt, old: &NortonJ2State) -> Result<NortonJ2Response> {
        self.update(strain, old).map(|next| next.response)
    }

    /// Propose and audit one simultaneous interval without changing old history.
    /// The total stored account includes elastic and energetic hardening storage.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old: &NortonJ2State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<NortonJ2State>> {
        require(
            positive(tolerances.stress_pa) && positive(tolerances.energy_j_m3),
            "stress and energy tolerances must be positive and finite",
        )?;
        self.validate_history(old_strain, old)?;
        let next = self.update(strain, old)?;
        let old_e = Self::elastic_strain(old_strain, old);
        let new_e = Self::elastic_strain(strain, &next.response.state);
        let stored_stress = self.creep.elastic.stress(&new_e, &());
        let q = next.relaxed.q;
        let dc = next.relaxed.delta_c;
        let dp = next.relaxed.delta_p;
        let mut residual = (q + 3.0 * self.creep.shear_pa * (dc + dp) - next.q_trial).abs();
        let yield_residual = q
            - self.plastic.yield_stress
            - self.plastic.hardening * next.response.state.plastic.alpha;
        residual = residual.max(if dp > 0.0 {
            yield_residual.abs()
        } else {
            yield_residual.max(0.0)
        });
        for (stored, returned) in stored_stress.iter().zip(&next.response.stress_pa) {
            residual = residual.max((stored - returned).abs());
        }
        if !residual.is_finite() || residual > tolerances.stress_pa {
            return Err(CreepError::StressResidual {
                residual_pa: residual,
                tolerance_pa: tolerances.stress_pa,
            });
        }
        let initial =
            self.plastic
                .energy_of(&old_e, old.plastic.alpha, &old.plastic.plastic_strain)?;
        let final_state = self.plastic.energy_of(
            &new_e,
            next.response.state.plastic.alpha,
            &next.response.state.plastic.plastic_strain,
        )?;
        let numerical = self.plastic.energy_of(
            &difference(&new_e, &old_e),
            dp,
            &difference(
                &next.response.state.plastic.plastic_strain,
                &old.plastic.plastic_strain,
            ),
        )?;
        let creep_heat = q * dc;
        let plastic_heat = self.plastic.yield_stress * dp;
        let mut energy = CreepEnergy {
            initial_stored_j_m3: initial.total_j_m3(),
            final_stored_j_m3: final_state.total_j_m3(),
            mechanical_work_j_m3: contract(
                &next.response.stress_pa,
                &difference(strain, old_strain),
            ),
            physical_dissipation_j_m3: creep_heat + plastic_heat,
            creep_dissipation_j_m3: creep_heat,
            plastic_dissipation_j_m3: plastic_heat,
            thermoelastic_storage_change_j_m3: 0.0,
            reversible_heat_j_m3: 0.0,
            numerical_dissipation_j_m3: numerical.total_j_m3(),
            residual_j_m3: 0.0,
        };
        require(
            energy.mechanical_work_j_m3.is_finite()
                && energy.physical_dissipation_j_m3.is_finite()
                && creep_heat >= 0.0
                && plastic_heat >= 0.0,
            "coupled mechanical work or physical heat is not representable",
        )?;
        energy.residual_j_m3 = energy.mechanical_work_j_m3
            - (energy.final_stored_j_m3 - energy.initial_stored_j_m3)
            - energy.physical_dissipation_j_m3
            - energy.numerical_dissipation_j_m3;
        if !energy.residual_j_m3.is_finite() || energy.residual_j_m3.abs() > tolerances.energy_j_m3
        {
            return Err(CreepError::EnergyResidual {
                residual_j_m3: energy.residual_j_m3,
                tolerance_j_m3: tolerances.energy_j_m3,
            });
        }
        Ok(CreepIncrement {
            strain: *strain,
            stress_pa: next.response.stress_pa,
            tangent_pa: next.response.tangent_pa,
            state: next.response.state,
            energy,
        })
    }
}

impl CreepLaw for NortonJ2 {
    fn kinetics(&self) -> &NortonCreep {
        &self.creep
    }
    fn validate_creep_history(&self, strain: &Voigt, state: &NortonJ2State) -> Result<()> {
        self.validate_history(strain, state)
    }
    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &NortonJ2State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<NortonJ2State>> {
        self.advance_with_energy(old_strain, old, strain, tolerances)
    }
}

impl SmallStrainLaw for NortonJ2 {
    type State = NortonJ2State;
    fn initial_state(&self) -> NortonJ2State {
        NortonJ2State {
            creep: self.creep.initial_state(),
            plastic: self.plastic.initial_state(),
        }
    }
    fn stress(&self, strain: &Voigt, old: &NortonJ2State) -> Voigt {
        self.update(strain, old)
            .map(|v| v.response.stress_pa)
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, strain: &Voigt, old: &NortonJ2State) -> Tangent6 {
        self.update(strain, old)
            .map(|v| v.response.tangent_pa)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, strain: &Voigt, old: &NortonJ2State) -> NortonJ2State {
        self.update(strain, old)
            .map(|v| v.response.state)
            .unwrap_or(NortonJ2State {
                creep: CreepState {
                    creep_strain: [f64::NAN; 6],
                    equivalent_creep_strain: f64::NAN,
                },
                plastic: J2State {
                    plastic_strain: [f64::NAN; 6],
                    alpha: f64::NAN,
                },
            })
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "bounded small-strain secondary creep and energetic isotropic J2 plasticity",
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "material.norton-isotropic-j2",
            "0.1.0",
            Ambition::Solid,
            vec![
                "simultaneous Norton/Arrhenius creep and associative isotropic J2 plasticity"
                    .into(),
                "separate permanent histories and energetic linear hardening".into(),
                "backward Euler at explicit physical duration; numerical loss is not heat".into(),
            ],
            ValidityDomain::unconstrained()
                .with("strain-magnitude", 0.0, self.creep.elastic.strain_limit)
                .with(
                    "temperature-k",
                    self.creep.parameters.temperature_range_k[0],
                    self.creep.parameters.temperature_range_k[1],
                ),
            vec![
                "no kinematic hardening, finite deformation, damage, tertiary creep or rupture"
                    .into(),
                "no measured material qualification without original source binding".into(),
            ],
            0.05,
        )
    }
}
