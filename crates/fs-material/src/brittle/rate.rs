//! Physical fracture-plane viscosity on the original brittle crack band.
//!
//! The native cohesive law supplies `t_v = (1-d) eta * jump_rate`, its complete
//! straight-path work integral and consistent endpoint/path derivatives. Eta is
//! explicitly authored in Pa s/m; no bulk viscosity or numerical relaxation
//! time is converted into a fracture-plane coefficient. Intrinsic fracture
//! work and physical viscous loss remain separate. The fixed width and modal
//! projection still map area work to volume work by division by h.
//!
//! This is a declared rate-sensitive traction law, not a rate-dependent source
//! strength, evolving crack direction or proof of arbitrary-width objectivity.
//! In particular, the damage integrity depends on the original M/h stiffness;
//! changing h need not preserve viscous work. The intrinsic Gc already excludes
//! the additional dashpot loss. Previously failed bonds have no dashpots and
//! retain the native independent compression law without an adhesive rate cap.
//!
//! Timed endpoint traction uses the interval secant strain rate. An accepted
//! midpoint trajectory can have a different instantaneous endpoint velocity;
//! `stress_at_rate` observes that actual rate without a fictitious time step.

use super::{
    BrittleCrackBand, BrittleError, BrittleState, BrittleStep, add_outer, contract, finite,
    finite_tangent, finite_voigt, positive,
};
use crate::cohesive::{CohesiveRateStep, CohesiveViscosity};
use crate::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_evidence::ModelCard;

/// Original intrinsic response plus the total endpoint response at the
/// interval's constant secant strain rate. No second damage history is stored.
#[derive(Clone, Debug, PartialEq)]
pub struct BrittleRateStep<'a> {
    /// Original equilibrium stress, tangents, history and intrinsic fracture work.
    pub equilibrium: BrittleStep<'a>,
    /// Total endpoint stress at `(new_strain-old_strain)/dt`, Pa; this need not
    /// be the stress at a dynamics owner's instantaneous endpoint velocity.
    pub stress_pa: Voigt,
    /// Classical endpoint derivative at fixed old state and fixed positive dt.
    pub tangent_pa: Option<Tangent6>,
    /// Native semismooth endpoint branch, including integrity/rate sensitivity.
    pub branch_tangent_pa: Tangent6,
    /// Integrated nonnegative physical dashpot loss, J/m3; not fracture or heat.
    pub viscous_dissipated_j_m3: f64,
}

impl BrittleRateStep<'_> {
    /// Intrinsic path work plus the separately represented viscous loss.
    #[must_use]
    pub fn total_work_j_m3(&self) -> f64 {
        self.equilibrium.energy.work + self.viscous_dissipated_j_m3
    }
}

/// Timed endpoint history and the distinct work-conjugate mean total stress.
#[derive(Clone, Debug, PartialEq)]
pub struct BrittleRatePathStep<'a> {
    pub endpoint: BrittleRateStep<'a>,
    /// Analytic mean total stress on the declared constant-rate strain path, Pa.
    pub mean_stress_pa: Voigt,
    /// `mean_stress : (new_strain-old_strain)`, J/m3.
    pub path_work_j_m3: f64,
    /// Mean-stress work minus intrinsic work and viscous loss, J/m3.
    pub path_work_defect_j_m3: f64,
}

/// Endpoint-force adapter with an explicit fixed physical interval.
#[derive(Clone, Copy, Debug)]
pub struct BrittleRateSmallStrain<'a> {
    owner: &'a BrittleCrackBand,
    dt_s: f64,
}

/// Analytic path-mean force adapter at fixed accepted history and physical dt.
#[derive(Clone, Copy, Debug)]
pub struct BrittleRateMeanSmallStrain<'a> {
    owner: &'a BrittleCrackBand,
    dt_s: f64,
}

impl BrittleCrackBand {
    /// Attach explicit fracture-plane viscosities before creating any history.
    /// Coefficients retain Pa s/m units; they are not divided by band width.
    /// Intrinsic Gc must exclude the additional measured viscous dissipation.
    #[must_use]
    pub fn with_viscosity(mut self, viscosity: CohesiveViscosity) -> Self {
        self.cohesive = self.cohesive.with_viscosity(viscosity);
        self
    }

    /// Original optional physical fracture-plane viscosity and rate bound.
    #[must_use]
    pub const fn viscosity(&self) -> Option<CohesiveViscosity> {
        self.cohesive.viscosity()
    }

    /// Borrow a timed endpoint-force law; all trial iterates use the same dt.
    pub fn small_strain_law_over_time(
        &self,
        dt_s: f64,
    ) -> Result<BrittleRateSmallStrain<'_>, BrittleError> {
        positive(dt_s)?;
        Ok(BrittleRateSmallStrain { owner: self, dt_s })
    }

    /// Borrow the work-conjugate mean-force law for the explicit interval.
    /// A missing classical path derivative remains an explicit refusal through
    /// the existing SmallStrainLaw nonfinite-tangent convention.
    pub fn mean_small_strain_law_over_time(
        &self,
        dt_s: f64,
    ) -> Result<BrittleRateMeanSmallStrain<'_>, BrittleError> {
        positive(dt_s)?;
        Ok(BrittleRateMeanSmallStrain { owner: self, dt_s })
    }

    /// Total endpoint stress on one straight constant-rate strain increment.
    /// The original owner, strain limits, damage, rate and work gates apply.
    pub fn trial_over_time<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
        dt_s: f64,
    ) -> Result<BrittleRateStep<'a>, BrittleError> {
        let jump = self.jump_for_trial(strain, state)?;
        let step = state
            .point
            .trial_over_time(jump, dt_s, self.area_energy_tolerance)?;
        self.finish_rate_trial(strain, state, step)
    }

    /// The same timed endpoint plus the native analytic mean TOTAL stress.
    /// Its complete tensor work includes intrinsic fracture and viscous loss;
    /// multiplying endpoint stress by the increment is not this integral.
    pub fn trial_over_time_with_mean_stress<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
        dt_s: f64,
    ) -> Result<BrittleRatePathStep<'a>, BrittleError> {
        let jump = self.jump_for_trial(strain, state)?;
        let (step, mean) = state.point.trial_over_time_with_mean_traction(
            jump,
            dt_s,
            self.area_energy_tolerance,
        )?;
        self.finish_rate_path(strain, state, step, mean)
    }

    /// Timed endpoint and mean, plus d(mean stress)/d(new tensor strain) at
    /// fixed old history and dt. Endpoint and path derivatives stay distinct.
    /// Native None at contact/damage switches is never replaced by an endpoint
    /// tangent or a zero Jacobian; the nonlinear consumer owns its fallback.
    pub fn trial_over_time_with_mean_stress_tangent<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
        dt_s: f64,
    ) -> Result<(BrittleRatePathStep<'a>, Option<Tangent6>), BrittleError> {
        let jump = self.jump_for_trial(strain, state)?;
        let (step, mean, local_tangent) = state.point.trial_over_time_with_mean_traction_tangent(
            jump,
            dt_s,
            self.area_energy_tolerance,
        )?;
        let path = self.finish_rate_path(strain, state, step, mean)?;
        let tangent = local_tangent
            .map(|local| self.rate_mean_tangent_from_modes(local))
            .transpose()?;
        Ok((path, tangent))
    }

    /// Observe total stress at an accepted state and an actual instantaneous
    /// tensor strain rate [1/s], without advancing damage or charging work.
    /// The original jump-rate bound applies only while physical bonds survive.
    pub fn stress_at_rate(
        &self,
        state: &BrittleState<'_>,
        strain_rate: &Voigt,
    ) -> Result<Voigt, BrittleError> {
        self.jump_for_trial(state.strain(), state)?;
        finite_voigt(strain_rate)?;
        // Inactive dashpots must not overflow an otherwise irrelevant rate
        // projection after separation or for a purely intrinsic material.
        let jump_rate = if self.viscosity().is_some() && !state.failed() {
            let rate = self
                .modes
                .map(|mode| self.settings.band_width_m * contract(&mode, strain_rate));
            for value in rate {
                finite(value)?;
            }
            rate
        } else {
            [0.0; 3]
        };
        let traction = state.point.traction_at_rate(jump_rate)?;
        self.stress_from_modes(state.strain(), traction)
    }

    fn finish_rate_trial<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
        step: CohesiveRateStep<'a>,
    ) -> Result<BrittleRateStep<'a>, BrittleError> {
        let stress_pa = self.stress_from_modes(strain, step.traction_pa)?;
        let branch_tangent_pa = self.tangent_from_modes(step.branch_tangent_pa_m)?;
        let tangent_pa = step.tangent_pa_m.map(|_| branch_tangent_pa);
        let viscous_dissipated_j_m3 =
            finite(step.viscous_dissipated_j_m2 / self.settings.band_width_m)?;
        if step.viscous_dissipated_j_m2 > 0.0 && viscous_dissipated_j_m3 == 0.0 {
            return Err(BrittleError::Invalid(
                "viscous fracture-plane work per volume is unrepresentable",
            ));
        }
        let equilibrium = self.finish_trial(strain, state, step.equilibrium)?;
        let total_work = finite(equilibrium.energy.work + viscous_dissipated_j_m3)?;
        let defect = finite(
            ((total_work - equilibrium.energy.stored_change) - equilibrium.energy.dissipated)
                - viscous_dissipated_j_m3,
        )?;
        if defect.abs() > self.settings.energy_tolerance_j_m3 {
            return Err(BrittleError::EnergyDefect(defect));
        }
        Ok(BrittleRateStep {
            equilibrium,
            stress_pa,
            tangent_pa,
            branch_tangent_pa,
            viscous_dissipated_j_m3,
        })
    }

    fn finish_rate_path<'a>(
        &'a self,
        strain: &Voigt,
        state: &BrittleState<'a>,
        step: CohesiveRateStep<'a>,
        mean_traction: [f64; 3],
    ) -> Result<BrittleRatePathStep<'a>, BrittleError> {
        let midpoint = core::array::from_fn(|i| 0.5 * state.strain[i] + 0.5 * strain[i]);
        let mean_stress_pa = self.stress_from_modes(&midpoint, mean_traction)?;
        let endpoint = self.finish_rate_trial(strain, state, step)?;
        let change = core::array::from_fn(|i| strain[i] - state.strain[i]);
        let path_work_j_m3 = finite(contract(&mean_stress_pa, &change))?;
        let path_work_defect_j_m3 = finite(path_work_j_m3 - endpoint.total_work_j_m3())?;
        if path_work_defect_j_m3.abs() > self.settings.energy_tolerance_j_m3 {
            return Err(BrittleError::EnergyDefect(path_work_defect_j_m3));
        }
        Ok(BrittleRatePathStep {
            endpoint,
            mean_stress_pa,
            path_work_j_m3,
            path_work_defect_j_m3,
        })
    }

    fn rate_mean_tangent_from_modes(&self, local: [[f64; 3]; 3]) -> Result<Tangent6, BrittleError> {
        let mut tangent = self
            .parallel_tangent
            .map(|row| row.map(|value| 0.5 * value));
        for (i, row) in local.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                add_outer(
                    &mut tangent,
                    &self.modes[i],
                    &self.modes[j],
                    finite(self.settings.band_width_m * value)?,
                );
            }
        }
        finite_tangent(&tangent)?;
        Ok(tangent)
    }
}

impl<'a> BrittleRateSmallStrain<'a> {
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }

    /// Checked endpoint response without publishing the returned history.
    pub fn trial(
        &self,
        strain: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<BrittleRateStep<'a>, BrittleError> {
        self.owner.trial_over_time(strain, state, self.dt_s)
    }
}

impl<'a> BrittleRateMeanSmallStrain<'a> {
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.dt_s
    }

    /// Checked physical trial and optional classical path-mean derivative.
    pub fn trial(
        &self,
        strain: &Voigt,
        state: &BrittleState<'a>,
    ) -> Result<(BrittleRatePathStep<'a>, Option<Tangent6>), BrittleError> {
        self.owner
            .trial_over_time_with_mean_stress_tangent(strain, state, self.dt_s)
    }
}

impl<'a> SmallStrainLaw for BrittleRateSmallStrain<'a> {
    type State = BrittleState<'a>;
    fn initial_state(&self) -> Self::State {
        self.owner.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.trial(strain, state)
            .map(|step| step.stress_pa)
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.trial(strain, state)
            .map(|step| step.branch_tangent_pa)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.trial(strain, state)
            .map(|step| step.equilibrium.state)
            .unwrap_or_else(|_| invalid_state(self.owner))
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.small_strain_law().admissibility()
    }
    fn card(&self) -> ModelCard {
        rate_card(self.owner, false)
    }
}

impl<'a> SmallStrainLaw for BrittleRateMeanSmallStrain<'a> {
    type State = BrittleState<'a>;
    fn initial_state(&self) -> Self::State {
        self.owner.initial_state()
    }
    fn stress(&self, strain: &Voigt, state: &Self::State) -> Voigt {
        self.owner
            .trial_over_time_with_mean_stress(strain, state, self.dt_s)
            .map(|path| path.mean_stress_pa)
            .unwrap_or([f64::NAN; 6])
    }
    fn tangent(&self, strain: &Voigt, state: &Self::State) -> Tangent6 {
        self.trial(strain, state)
            .ok()
            .and_then(|(_, tangent)| tangent)
            .unwrap_or([[f64::NAN; 6]; 6])
    }
    fn update_state(&self, strain: &Voigt, state: &Self::State) -> Self::State {
        self.owner
            .trial_over_time(strain, state, self.dt_s)
            .map(|step| step.equilibrium.state)
            .unwrap_or_else(|_| invalid_state(self.owner))
    }
    fn admissibility(&self) -> MaterialAdmissibility {
        self.owner.small_strain_law().admissibility()
    }
    fn card(&self) -> ModelCard {
        rate_card(self.owner, true)
    }
}

fn invalid_state(owner: &BrittleCrackBand) -> BrittleState<'_> {
    BrittleState {
        owner,
        strain: [f64::NAN; 6],
        point: owner.cohesive.bonded_point(),
    }
}

fn rate_card(owner: &BrittleCrackBand, mean: bool) -> ModelCard {
    let mut card = owner.small_strain_law().card();
    card.name = if mean {
        "material.brittle-fixed-plane-rate-path-mean"
    } else {
        "material.brittle-fixed-plane-rate-endpoint"
    }
    .into();
    card.assumptions.push("one straight tensor-strain path over an explicit positive physical dt; native surviving-bond viscosity in Pa s/m and intrinsic fracture work retained separately".into());
    card.known_failures.push("interval secant-rate endpoint stress is not an instantaneous endpoint-velocity observation; no arbitrary-width objectivity of additional viscous work is inferred".into());
    if mean {
        card.known_failures.push("path-mean stress is not endpoint stress; the native classical path derivative can be absent at contact or damage switches, requiring the consumer's admitted fallback".into());
    }
    card.assumptions.sort_unstable();
    card.known_failures.sort_unstable();
    card
}
