//! Finite-rate mineral exchange with fast ORIGINAL aqueous equilibria.
//!
//! The caller declares a fixed-surface formula-unit rate K_A [mol/s]. The
//! selected original mineral supplies Omega=IAP/K and the rate is
//! du/dt=K_A(1-Omega); every trial re-solves the original aqueous charge and
//! acid/base equations. No apparent neutral solubility, frozen pH, empirical
//! coefficient, nucleation mechanism or excess caloric model is inferred.
//!
//! This is the affinity form described in the USGS PHREEQC theory manual,
//! "Rate Expressions":
//! <https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-20.html>.
//! A numerical coefficient alone does not establish empirical applicability.
//! The original fixed-temperature, ideal-activity and fixed-solvent source
//! boundaries continue to apply. Time integration is first order, controlled
//! by bounded full-step/two-half-step comparison; its estimate is not a
//! rigorous global error bound. Only one selected solid evolves.

use super::{HydrolysisPrecipitation, PrecipitationOptions, salt::SaltPrecipitation};
use super::super::{
    AcidBaseError, AcidBaseSpeciation, ConservedReaction, IdealAcidBaseSystem, Result,
    SpeciesId, SpeciationOptions, SpectatorIon, StandardReactionGibbs, Sum, charge,
    solve::{AcidBaseTrial, PreparedAcidBase},
};
use fs_math::det;
use fs_qty::{Amount, Dimensionless, Time};
use std::sync::Arc;

mod solve;

/// Growth requires an existing crystal or an explicitly persistent seed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MineralGrowthPolicyV1 {
    ExistingSolidOnly,
    SeededSurface,
}

/// Explicit fixed-area affinity coefficient [mol of original solid/s].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MineralSurfaceRateV1 {
    rate_mol_s: f64,
    growth_policy: MineralGrowthPolicyV1,
}
impl MineralSurfaceRateV1 {
    pub fn new(rate_mol_s: f64, growth_policy: MineralGrowthPolicyV1) -> Result<Self> {
        nonnegative(rate_mol_s, "finite nonnegative mineral surface rate required")?;
        Ok(Self { rate_mol_s, growth_policy })
    }
    pub const fn rate_mol_s(&self) -> f64 { self.rate_mol_s }
    pub const fn growth_policy(&self) -> MineralGrowthPolicyV1 { self.growth_policy }
}

/// Original conserved equation, standard-state sources and aqueous owner.
#[derive(Clone, Debug)]
pub enum MineralKineticSourceV1 {
    Salt(SaltPrecipitation),
    Hydrolysis(HydrolysisPrecipitation),
}
impl MineralKineticSourceV1 {
    pub fn aqueous_model(&self) -> &IdealAcidBaseSystem {
        match self { Self::Salt(s) => s.aqueous_model(), Self::Hydrolysis(s) => s.aqueous_model() }
    }
    pub fn standard(&self) -> &StandardReactionGibbs {
        match self { Self::Salt(s) => s.standard(), Self::Hydrolysis(s) => s.standard() }
    }
    pub fn reaction(&self) -> &ConservedReaction {
        match self { Self::Salt(s) => s.reaction(), Self::Hydrolysis(s) => s.reaction() }
    }
    pub fn solid(&self) -> &SpeciesId {
        match self { Self::Salt(s) => s.solid(), Self::Hydrolysis(s) => s.solid() }
    }
    pub fn cation(&self) -> &SpeciesId {
        match self { Self::Salt(s) => s.cation(), Self::Hydrolysis(s) => s.cation() }
    }
    pub fn family_id(&self) -> Option<&str> {
        match self { Self::Salt(s) => Some(s.family_id()), Self::Hydrolysis(_) => None }
    }
    pub fn cation_stoichiometry(&self) -> u32 {
        match self { Self::Salt(s) => s.cation_stoichiometry(), Self::Hydrolysis(_) => 1 }
    }
    pub fn anion_stoichiometry(&self) -> u32 {
        match self { Self::Salt(s) => s.anion_stoichiometry(), Self::Hydrolysis(_) => 0 }
    }
}

#[derive(Debug)]
struct Owner {
    source: MineralKineticSourceV1,
    rate: MineralSurfaceRateV1,
    cation_charge: i32,
    family_index: Option<usize>,
}

/// One immutable source/rate binding. Clones preserve exact owner identity.
#[derive(Clone, Debug)]
pub struct MineralKineticsV1 { owner: Arc<Owner> }
impl MineralKineticsV1 {
    pub fn from_salt(source: SaltPrecipitation, rate: MineralSurfaceRateV1) -> Result<Self> {
        Self::new(MineralKineticSourceV1::Salt(source), rate)
    }
    pub fn from_hydrolysis(source: HydrolysisPrecipitation, rate: MineralSurfaceRateV1) -> Result<Self> {
        Self::new(MineralKineticSourceV1::Hydrolysis(source), rate)
    }
    fn new(source: MineralKineticSourceV1, rate: MineralSurfaceRateV1) -> Result<Self> {
        let cation_charge = charge(source.reaction(), source.cation())?;
        let family_index = source.family_id().map(|id| {
            source.aqueous_model().families().iter().position(|family| family.id() == id)
                .ok_or(AcidBaseError::Invalid("missing original mineral family"))
        }).transpose()?;
        Ok(Self { owner: Arc::new(Owner { source, rate, cation_charge, family_index }) })
    }
    pub fn source(&self) -> &MineralKineticSourceV1 { &self.owner.source }
    pub fn rate(&self) -> &MineralSurfaceRateV1 { &self.owner.rate }

    /// Declare actual dissolved analytical pools and an independent solid stock.
    /// Complete totals are formed ONCE and retained thereafter. Dissolved
    /// inputs exclude the solid's stoichiometric contribution. The selected
    /// cation must not also appear among caller spectators or family members.
    /// Aqueous equilibration does not move either mineral phase.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare(
        &self,
        dissolved_cation: Amount,
        solid: Amount,
        dissolved_family_totals: Vec<(String, Amount)>,
        spectators: Vec<SpectatorIon>,
        amount_at_unit_activity: Amount,
        options: SpeciationOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MineralKineticStateV1> {
        poll(&mut cancelled)?;
        nonnegative(dissolved_cation.value(), "invalid dissolved mineral cation")?;
        nonnegative(solid.value(), "invalid actual mineral solid")?;
        if spectators.iter().any(|ion| ion.species == *self.source().cation() || ion.species == *self.source().solid()) {
            return Err(AcidBaseError::Invalid("duplicated owned mineral participant"));
        }
        // The original aqueous owner validates the complete identities and
        // canonical family/spectator layout before any kinetic state exists.
        let mut ions = spectators.clone();
        ions.push(SpectatorIon { species: self.source().cation().clone(), charge: self.owner.cation_charge, amount: dissolved_cation });
        self.source().aqueous_model().prepare_speciation(
            &dissolved_family_totals, &ions, amount_at_unit_activity, options, &mut cancelled,
        )?;
        let mut families: Vec<_> = self.source().aqueous_model().families().iter().map(|family| {
            dissolved_family_totals.iter().find(|(id, _)| id == family.id()).cloned()
                .ok_or(AcidBaseError::Invalid("missing original analytical family"))
        }).collect::<Result<_>>()?;
        let selected = self.owner.family_index.map(|i| families[i].1.value());
        let (u, excess, cation_limits) = match &self.owner.source {
            MineralKineticSourceV1::Salt(source) => super::salt::partition_basis(
                dissolved_cation.value(), selected.expect("admitted salt family"),
                source.cation_stoichiometry(), source.anion_stoichiometry(),
            )?,
            MineralKineticSourceV1::Hydrolysis(_) => (dissolved_cation.value(), 0.0, true),
        };
        let p = f64::from(self.source().cation_stoichiometry());
        let q = f64::from(self.source().anion_stoichiometry());
        let component_total = finite(dissolved_cation.value() + p * solid.value(), "complete mineral cation total")?;
        if let Some(index) = self.owner.family_index {
            families[index].1 = Amount::new(finite(families[index].1.value() + q * solid.value(), "complete mineral family total")?);
        }
        let mut spectators = spectators;
        spectators.sort_by(|a, b| a.species.cmp(&b.species));
        let stocks = Arc::new(Stocks {
            component_total: Amount::new(component_total), family_totals: families.into_boxed_slice(),
            spectators: spectators.into_boxed_slice(),
            total_formula: finite(u + solid.value(), "complete mineral formula total")?,
            anchor_u: u, anchor_metal: dissolved_cation.value(), anchor_family: selected,
            excess, cation_limits,
        });
        let phase = Phase { dissolved: u, solid: solid.value() };
        let workspace = Workspace::new(self, &stocks, amount_at_unit_activity, options, &mut cancelled)?;
        workspace.finish(self, Arc::clone(&stocks), phase, amount_at_unit_activity, options, &mut cancelled).map(|(state, _)| state)
    }

    /// Change the original solvent/activity basis while freezing both mineral
    /// stocks. Only fast aqueous acid/base and charge equilibrium is recomputed.
    pub fn respeciate_at_basis(
        &self, initial: &MineralKineticStateV1, amount_at_unit_activity: Amount,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MineralKineticStateV1> {
        poll(&mut cancelled)?;
        self.admit_owner(initial)?;
        let workspace = Workspace::new(self, &initial.stocks, amount_at_unit_activity, initial.options, &mut cancelled)?;
        workspace.finish(self, Arc::clone(&initial.stocks), initial.phase, amount_at_unit_activity, initial.options, &mut cancelled).map(|(state, _)| state)
    }
    fn admit_owner(&self, state: &MineralKineticStateV1) -> Result<()> {
        if !Arc::ptr_eq(&self.owner, &state.model.owner) {
            Err(AcidBaseError::Invalid("mineral state requires its exact original source/rate owner"))
        } else { Ok(()) }
    }
}

#[derive(Debug)]
struct Stocks {
    component_total: Amount,
    family_totals: Box<[(String, Amount)]>,
    spectators: Box<[SpectatorIon]>,
    total_formula: f64,
    anchor_u: f64,
    anchor_metal: f64,
    anchor_family: Option<f64>,
    excess: f64,
    cation_limits: bool,
}
#[derive(Clone, Copy, Debug)]
struct Phase { dissolved: f64, solid: f64 }

/// Actual finite mineral stocks and fast aqueous equilibrium, with no mineral
/// equilibrium-status claim. Positive IAP/K may coexist with zero or positive
/// solid according to the explicit growth policy and elapsed finite time.
#[derive(Clone, Debug)]
pub struct MineralKineticStateV1 {
    model: MineralKineticsV1,
    stocks: Arc<Stocks>,
    phase: Phase,
    basis: Amount,
    options: SpeciationOptions,
    aqueous: AcidBaseSpeciation,
    metal: Amount,
    family: Option<Amount>,
    log_omega: Option<Dimensionless>,
    spread: Option<Dimensionless>,
    cation_residual: Amount,
    family_residual: Option<Amount>,
}
impl MineralKineticStateV1 {
    pub const fn model(&self) -> &MineralKineticsV1 { &self.model }
    pub const fn aqueous(&self) -> &AcidBaseSpeciation { &self.aqueous }
    pub fn component_total(&self) -> Amount { self.stocks.component_total }
    pub fn family_totals(&self) -> &[(String, Amount)] { &self.stocks.family_totals }
    pub fn spectators(&self) -> &[SpectatorIon] { &self.stocks.spectators }
    pub const fn dissolved_cation(&self) -> Amount { self.metal }
    pub const fn dissolved_family(&self) -> Option<Amount> { self.family }
    pub fn solid_amount(&self) -> Amount { Amount::new(self.phase.solid) }
    pub const fn amount_at_unit_activity(&self) -> Amount { self.basis }
    pub const fn options(&self) -> SpeciationOptions { self.options }
    /// Natural log of dissolution IAP/K. None is the exact IAP=0 limit.
    pub const fn log_saturation_ratio(&self) -> Option<Dimensionless> { self.log_omega }
    pub const fn log_saturation_numerical_spread(&self) -> Option<Dimensionless> { self.spread }
    pub const fn cation_balance_residual(&self) -> Amount { self.cation_residual }
    pub const fn family_balance_residual(&self) -> Option<Amount> { self.family_residual }
}

/// Bounded aqueous work, implicit extent resolution and temporal accuracy.
#[derive(Clone, Copy, Debug)]
pub struct MineralKineticControlsV1 {
    pub aqueous: SpeciationOptions,
    pub initial_step: Time,
    pub minimum_step: Time,
    pub maximum_step: Time,
    pub relative_tolerance: Dimensionless,
    pub amount_absolute_tolerance: Amount,
    /// Numerical implicit-root error allowance, INCLUDING inner pH uncertainty.
    /// A residual or a monotone remaining-stock displacement estimate qualifies
    /// the root; growth also limits relative donor error and honors the phase
    /// amount tolerances. This is not an interval-arithmetic certificate and
    /// does not floor a donor.
    pub extent_absolute_tolerance: Amount,
    pub maximum_stage_iterations: usize,
    pub maximum_attempts: usize,
}
impl Default for MineralKineticControlsV1 {
    fn default() -> Self {
        Self {
            aqueous: PrecipitationOptions::default().aqueous,
            initial_step: Time::new(0.1), minimum_step: Time::new(1e-8), maximum_step: Time::new(1.0),
            relative_tolerance: Dimensionless::new(1e-5), amount_absolute_tolerance: Amount::new(1e-12),
            extent_absolute_tolerance: Amount::new(1e-14), maximum_stage_iterations: 256, maximum_attempts: 16384,
        }
    }
}

/// Accepted native mineral evolution and bounded work receipts.
#[derive(Clone, Debug)]
pub struct MineralKineticOutcomeV1 {
    final_state: MineralKineticStateV1,
    duration: Time,
    transfer: Amount,
    accepted: usize,
    rejected: usize,
    counts: Counts,
    maximum_error: f64,
    blocked: bool,
    exhausted: bool,
}
impl MineralKineticOutcomeV1 {
    pub const fn final_state(&self) -> &MineralKineticStateV1 { &self.final_state }
    pub fn into_state(self) -> MineralKineticStateV1 { self.final_state }
    pub const fn duration(&self) -> Time { self.duration }
    /// Signed transferred original solid formula units; dissolution is positive.
    pub const fn transferred_to_solution(&self) -> Amount { self.transfer }
    pub const fn accepted_steps(&self) -> usize { self.accepted }
    pub const fn rejected_steps(&self) -> usize { self.rejected }
    pub const fn aqueous_solves(&self) -> usize { self.counts.solves }
    pub const fn aqueous_iterations(&self) -> usize { self.counts.iterations }
    pub const fn maximum_accepted_error_norm(&self) -> f64 { self.maximum_error }
    pub const fn precipitation_blocked(&self) -> bool { self.blocked }
    /// A finite solid was exhausted during accepted stages. This does not
    /// claim an exact continuous exhaustion time from an implicit trajectory.
    pub const fn solid_exhausted(&self) -> bool { self.exhausted }
}

#[derive(Clone, Copy, Debug, Default)]
struct Counts { solves: usize, iterations: usize }
struct Workspace<'a> {
    model: &'a MineralKineticsV1,
    stocks: &'a Stocks,
    prepared: PreparedAcidBase<'a>,
    cation_index: usize,
    basis: f64,
    counts: Counts,
}
struct Trial {
    phase: Phase,
    metal: f64,
    family: Option<f64>,
    aqueous: AcidBaseTrial,
    log_omega: Option<f64>,
    spread: f64,
    log_roundoff: f64,
}
impl<'a> Workspace<'a> {
    fn new(model: &'a MineralKineticsV1, stocks: &'a Stocks, basis: Amount, options: SpeciationOptions, cancelled: &mut impl FnMut() -> bool) -> Result<Self> {
        let mut families = stocks.family_totals.to_vec();
        if let Some(index) = model.owner.family_index { families[index].1 = Amount::new(stocks.anchor_family.expect("original family")); }
        let mut ions = stocks.spectators.to_vec();
        ions.push(SpectatorIon { species: model.source().cation().clone(), charge: model.owner.cation_charge, amount: Amount::new(stocks.anchor_metal) });
        let prepared = model.source().aqueous_model().prepare_speciation(&families, &ions, basis, options, cancelled)?;
        let cation_index = prepared.spectator_index(model.source().cation()).ok_or(AcidBaseError::Invalid("missing owned kinetic cation"))?;
        Ok(Self { model, stocks, prepared, cation_index, basis: basis.value(), counts: Counts::default() })
    }
    fn evaluate(&mut self, phase: Phase, cancelled: &mut impl FnMut() -> bool) -> Result<Trial> {
        poll(cancelled)?;
        balance(phase.dissolved, phase.solid, self.stocks.total_formula)?;
        let p = f64::from(self.model.source().cation_stoichiometry());
        let q = f64::from(self.model.source().anion_stoichiometry());
        let (metal, family) = if phase.dissolved == self.stocks.anchor_u {
            (self.stocks.anchor_metal, self.stocks.anchor_family)
        } else if self.model.owner.family_index.is_none() {
            (phase.dissolved, None)
        } else if self.stocks.cation_limits {
            (p * phase.dissolved, Some(self.stocks.excess + q * phase.dissolved))
        } else {
            (self.stocks.excess + p * phase.dissolved, Some(q * phase.dissolved))
        };
        nonnegative(metal, "kinetic cation loses range")?;
        if let Some(family) = family {
            nonnegative(family, "kinetic family loses range")?;
            self.prepared.set_family_total(self.model.owner.family_index.expect("salt family"), Amount::new(family))?;
        }
        if phase.dissolved > 0.0 && (metal == 0.0 || family == Some(0.0)) {
            return Err(AcidBaseError::Arithmetic("positive dissolved mineral component underflows"));
        }
        self.prepared.set_spectator_amount(self.cation_index, Amount::new(metal))?;
        let aqueous = self.prepared.solve(cancelled)?;
        self.counts.solves += 1;
        self.counts.iterations += aqueous.iterations();
        let width = aqueous.ph_bracket().iter().map(|ph| (ph - aqueous.ph()).abs()).fold(0.0, f64::max);
        let (log_omega, spread, log_roundoff) = match &self.model.owner.source {
            MineralKineticSourceV1::Salt(source) => {
                let index = self.model.owner.family_index.expect("salt family");
                let anion = self.prepared.terminal_activity(index, &aqueous)?;
                let (logarithm, rounding) = if metal == 0.0 || anion == 0.0 { (None, 0.0) } else {
                    let terms = [
                        p * det::ln(metal), -p * det::ln(self.basis),
                        q * det::ln(anion), source.standard().log_equilibrium_constant().value(),
                    ];
                    let mut sum = Sum::default();
                    for value in terms { sum.add(value)?; }
                    (Some(sum.finish()?), logarithmic_roundoff(&terms))
                };
                (logarithm, q * source.aqueous_model().families()[index].steps().len() as f64 * core::f64::consts::LN_10 * width, rounding)
            }
            MineralKineticSourceV1::Hydrolysis(source) => {
                let z = f64::from(self.model.owner.cation_charge);
                let (logarithm, rounding) = if metal == 0.0 { (None, 0.0) } else {
                    let log_metal = det::ln(metal);
                    let log_basis = det::ln(self.basis);
                    let proton = -z * det::ln(aqueous.proton_activity());
                    let constant = source.standard().log_equilibrium_constant().value();
                    let mut sum = Sum::default();
                    sum.add(log_metal - log_basis)?;
                    sum.add(proton)?;
                    sum.add(constant)?;
                    (Some(sum.finish()?), logarithmic_roundoff(&[log_metal, log_basis, proton, constant]))
                };
                (logarithm, z * core::f64::consts::LN_10 * width, rounding)
            }
        };
        let spread = finite(spread, "kinetic saturation numerical spread")?;
        let log_roundoff = finite(log_roundoff, "kinetic logarithmic roundoff scale")?;
        Ok(Trial { phase, metal, family, aqueous, log_omega, spread, log_roundoff })
    }
    fn finish(mut self, model: &MineralKineticsV1, stocks: Arc<Stocks>, phase: Phase, basis: Amount, options: SpeciationOptions, cancelled: &mut impl FnMut() -> bool) -> Result<(MineralKineticStateV1, Counts)> {
        let trial = self.evaluate(phase, cancelled)?;
        let p = f64::from(model.source().cation_stoichiometry());
        let q = f64::from(model.source().anion_stoichiometry());
        let cation_residual = balance(trial.metal, p * phase.solid, stocks.component_total.value())?;
        let family_residual = trial.family.map(|family| {
            balance(family, q * phase.solid, stocks.family_totals[model.owner.family_index.expect("salt family")].1.value()).map(Amount::new)
        }).transpose()?;
        let aqueous = self.prepared.finish(trial.aqueous, cancelled)?;
        poll(cancelled)?;
        Ok((MineralKineticStateV1 {
            model: model.clone(), stocks, phase, basis, options, aqueous,
            metal: Amount::new(trial.metal), family: trial.family.map(Amount::new),
            log_omega: trial.log_omega.map(Dimensionless::new),
            spread: trial.log_omega.map(|_| Dimensionless::new(trial.spread)),
            cation_residual: Amount::new(cation_residual), family_residual,
        }, self.counts))
    }
}

// Numerical allowance for elementary logarithms and cancellation in a sum.
// Scale BEFORE adding magnitudes to avoid an unnecessary large-source overflow.
// This is a floating-point estimate, not an interval-arithmetic certificate.
fn logarithmic_roundoff(terms: &[f64]) -> f64 {
    let factor = 32.0 * f64::EPSILON;
    terms.iter().fold(factor, |sum, term| sum + factor * term.abs())
}

fn finite(value: f64, message: &'static str) -> Result<f64> {
    if value.is_finite() { Ok(value) } else { Err(AcidBaseError::Arithmetic(message)) }
}
fn nonnegative(value: f64, message: &'static str) -> Result<()> {
    if value.is_finite() && value >= 0.0 { Ok(()) } else { Err(AcidBaseError::Invalid(message)) }
}
fn balance(first: f64, second: f64, total: f64) -> Result<f64> {
    nonnegative(first, "negative mineral dissolved stock")?;
    nonnegative(second, "negative mineral solid stock")?;
    let mut sum = Sum::default();
    sum.add(first)?; sum.add(second)?; sum.add(-total)?;
    let residual = sum.finish()?;
    if residual.abs() > 64.0 * f64::EPSILON * total {
        return Err(AcidBaseError::Arithmetic("kinetic analytical component loses closure"));
    }
    Ok(residual)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() { Err(AcidBaseError::Cancelled) } else { Ok(()) }
}
