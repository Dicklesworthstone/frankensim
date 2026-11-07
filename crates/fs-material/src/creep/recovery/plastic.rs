//! Simultaneous primary-resistance recovery and energetic isotropic J2 flow.
//!
//! The original joint owner solves creep, plastic complementarity and stress
//! balance at one endpoint. Only its primary kinetic memory changes: with the
//! actual solved creep increment `dc`, `r1 = (r0+dc)/(1+k(T) dt)`. The prepared
//! primary scale is `a*(1+k(T) dt)`, so the original yield-crossing test,
//! fixed-stress perfect-plastic branch and hardening branch remain authoritative.
//!
//! Permanent creep alpha and plastic alpha never become resistance memory.
//! Native elastic/isotropic-hardening storage and the physical heat split
//! `q*dc + yield_stress*dp` are retained, with no extra resistance energy or heat.
//! Recovery parameters require their own calibration. This is small-strain
//! phenomenological behavior, not a tertiary-damage or rupture model.

/// Original source binding for all joint creep, recovery and plastic claims.
pub mod source;
/// Exact free-sided axial restriction retaining both permanent histories.
pub mod uniaxial;

use super::{CreepRecoveryParameters, RecoveryStep, validate_resistance};
use crate::creep::plastic::{NortonJ2, NortonJ2State};
use crate::creep::thermal::expansion::CreepHistoryStrain;
use crate::creep::{CreepIncrement, CreepLaw, CreepState, CreepTolerances, NortonCreep, Result};
use crate::plastic::J2State;
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::{Ambition, ModelCard, ValidityDomain};

/// Original permanent creep/plastic histories and separate kinetic resistance.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveringNortonJ2State {
    pub native: NortonJ2State,
    /// Nonnegative primary memory, bounded by the permanent creep path alone.
    pub resistance_strain: f64,
}

impl CreepHistoryStrain for RecoveringNortonJ2State {
    const STRAIN_FIELDS: usize = 2;

    fn inelastic_strain(&self) -> Voigt {
        self.native.inelastic_strain()
    }
}

/// One uncommitted response of the simultaneous recovering creep/J2 update.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveringNortonJ2Response {
    pub strain: Voigt,
    pub stress_pa: Voigt,
    pub tangent_pa: Tangent6,
    pub state: RecoveringNortonJ2State,
}

/// One endpoint temperature and duration, retaining the authored native laws.
#[derive(Clone, Debug, PartialEq)]
pub struct RecoveringNortonJ2 {
    native: NortonJ2,
    prepared: NortonJ2,
    recovery: RecoveryStep,
}

impl RecoveringNortonJ2 {
    pub fn new(native: NortonJ2, recovery: CreepRecoveryParameters) -> Result<Self> {
        let (creep, recovery) = RecoveryStep::prepare(native.creep(), recovery)?;
        let prepared = NortonJ2::new(creep, native.plastic().clone())?;
        Ok(Self {
            native,
            prepared,
            recovery,
        })
    }

    /// Original joint law. Its coefficients never include the timestep-scaled
    /// primary strain used privately by the recovery return.
    pub const fn norton_j2(&self) -> &NortonJ2 {
        &self.native
    }

    pub const fn recovery_parameters(&self) -> &CreepRecoveryParameters {
        &self.recovery.parameters
    }

    pub const fn recovery_rate_s(&self) -> f64 {
        self.recovery.rate_s
    }

    pub(super) const fn prepared_norton_j2(&self) -> &NortonJ2 {
        &self.prepared
    }

    pub(super) fn next_resistance(&self, old: f64, delta: f64, alpha: f64) -> Result<f64> {
        self.recovery.next_resistance(old, delta, alpha)
    }

    /// Admit both original physical histories and their separate creep resistance.
    pub fn validate_history(&self, strain: &Voigt, state: &RecoveringNortonJ2State) -> Result<()> {
        self.native.validate_history(strain, &state.native)?;
        validate_resistance(
            state.resistance_strain,
            state.native.creep.equivalent_creep_strain,
        )
    }

    /// Evaluate one coupled trial without mutating either accepted history.
    pub fn response(
        &self,
        strain: &Voigt,
        old: &RecoveringNortonJ2State,
    ) -> Result<RecoveringNortonJ2Response> {
        validate_resistance(
            old.resistance_strain,
            old.native.creep.equivalent_creep_strain,
        )?;
        let (next, delta) =
            self.prepared
                .response_with_memory(strain, &old.native, old.resistance_strain)?;
        let resistance_strain = self.next_resistance(
            old.resistance_strain,
            delta,
            next.state.creep.equivalent_creep_strain,
        )?;
        Ok(RecoveringNortonJ2Response {
            strain: next.strain,
            stress_pa: next.stress_pa,
            tangent_pa: next.tangent_pa,
            state: RecoveringNortonJ2State {
                native: next.state,
                resistance_strain,
            },
        })
    }

    /// Admit the original joint work/heat account and evolve resistance from
    /// the same solved creep increment. No second constitutive solve is used.
    pub fn advance_with_energy(
        &self,
        old_strain: &Voigt,
        old: &RecoveringNortonJ2State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<RecoveringNortonJ2State>> {
        validate_resistance(
            old.resistance_strain,
            old.native.creep.equivalent_creep_strain,
        )?;
        let (next, delta) = self.prepared.advance_with_memory(
            old_strain,
            &old.native,
            strain,
            tolerances,
            old.resistance_strain,
        )?;
        let resistance_strain = self.next_resistance(
            old.resistance_strain,
            delta,
            next.state.creep.equivalent_creep_strain,
        )?;
        Ok(CreepIncrement {
            strain: next.strain,
            stress_pa: next.stress_pa,
            tangent_pa: next.tangent_pa,
            state: RecoveringNortonJ2State {
                native: next.state,
                resistance_strain,
            },
            energy: next.energy,
        })
    }
}

impl CreepLaw for RecoveringNortonJ2 {
    fn kinetics(&self) -> &NortonCreep {
        self.native.creep()
    }

    fn validate_creep_history(&self, strain: &Voigt, state: &Self::State) -> Result<()> {
        self.validate_history(strain, state)
    }

    fn creep_increment(
        &self,
        old_strain: &Voigt,
        old: &Self::State,
        strain: &Voigt,
        tolerances: CreepTolerances,
    ) -> Result<CreepIncrement<Self::State>> {
        self.advance_with_energy(old_strain, old, strain, tolerances)
    }
}

impl SmallStrainLaw for RecoveringNortonJ2 {
    type State = RecoveringNortonJ2State;

    fn initial_state(&self) -> Self::State {
        RecoveringNortonJ2State {
            native: self.native.initial_state(),
            resistance_strain: 0.0,
        }
    }

    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.response(strain, state)
            .map(|v| v.stress_pa)
            .unwrap_or([f64::NAN; 6])
    }

    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.response(strain, state)
            .map(|v| v.tangent_pa)
            .unwrap_or([[f64::NAN; 6]; 6])
    }

    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.response(strain, state)
            .map(|v| v.state)
            .unwrap_or(RecoveringNortonJ2State {
                native: NortonJ2State {
                    creep: CreepState {
                        creep_strain: [f64::NAN; 6],
                        equivalent_creep_strain: f64::NAN,
                    },
                    plastic: J2State {
                        plastic_strain: [f64::NAN; 6],
                        alpha: f64::NAN,
                    },
                },
                resistance_strain: f64::NAN,
            })
    }

    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "bounded small-strain recovering primary creep and energetic isotropic J2 plasticity; no damage or rupture law",
        }
    }

    fn card(&self) -> ModelCard {
        let creep = self.native.creep();
        ModelCard::new(
            "material.recovering-primary-norton-isotropic-j2",
            "0.1.0",
            Ambition::Solid,
            vec![
                "simultaneous Norton creep and energetic isotropic J2 plasticity at one endpoint stress".into(),
                "independent primary resistance with first-order Arrhenius recovery".into(),
                "separate permanent creep/plastic paths; backward Euler at explicit duration".into(),
                "native elastic/hardening storage and creep/plastic heat; resistance stores no energy".into(),
            ],
            ValidityDomain::unconstrained()
                .with("strain-magnitude", 0.0, creep.elastic().strain_limit)
                .with(
                    "temperature-k",
                    creep.parameters().temperature_range_k[0],
                    creep.parameters().temperature_range_k[1],
                ),
            vec![
                "no kinematic hardening, finite deformation, viscoelastic strain recovery, tertiary damage or rupture".into(),
                "independent recovery calibration required; no measured qualification without original source binding".into(),
            ],
            0.05,
        )
    }
}
