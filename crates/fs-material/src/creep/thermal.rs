//! Uniform heat feedback for the explicitly selected zero-expansion reduction.
//!
//! One immutable card supplies Norton kinetics, constant elastic coefficients,
//! constant reference density and the original single-solid caloric chart.
//! Thermal expansion is explicitly omitted by this model selection; missing
//! expansion data never becomes a sourced zero. Uniform temperature/strain and
//! fixed reference volume are assumptions, with boundary lumping admitted by
//! the transport owner. No phase change, pressure work or finite deformation.
//!
//! Every within-step trial starts from the same accepted creep strain, total
//! strain, enthalpy and clock. Only physical creep dissipation heats the chart.
//! The final mechanics are recomputed at the actual returned temperature, and
//! a failed/cancelled/budget-exhausted step publishes no partially aged state.

use super::source::{CreepSourceError, NortonCreepSource};
use super::{CreepError, CreepIncrement, CreepState, CreepTolerances};
use crate::phase::{
    EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve, UniformEnthalpyStep,
    UniformEnthalpyStepInput,
};
use crate::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use crate::{SmallStrainLaw, Voigt};
use core::{convert::Infallible, fmt};
use fs_matdb::MaterialCard;
use fs_qty::{Dims, QuantitySpec};

pub const MAX_FEEDBACK_ITERATIONS: usize = 256;

/// Accepted native mechanical history and fixed reference specimen geometry.
#[derive(Clone, Debug)]
pub struct InitialState {
    pub strain: Voigt,
    pub history: CreepState,
    pub specific_enthalpy_j_kg: f64,
    pub mass_kg: f64,
    pub surface_area_m2: f64,
    pub time_s: f64,
}

/// Independent numerical gates; this is a bounded fixed point, not a stability proof.
#[derive(Clone, Copy, Debug)]
pub struct FeedbackConfig {
    pub max_iterations: usize,
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    pub stress_tolerance_pa: f64,
    /// Fixed-point damping in (0,1].
    pub relaxation: f64,
}

#[derive(Debug)]
pub enum FeedbackError<E = Infallible> {
    Source(CreepSourceError),
    Mechanical(CreepError),
    Caloric(PhaseStateError),
    Transport(E),
    Invalid(&'static str),
    Cancelled,
    NotConverged {
        iterations: usize,
        temperature_defect_k: f64,
        heat_defect_j: f64,
        thermal_residual_j: f64,
        combined_residual_j: f64,
    },
}
impl<E: fmt::Display> fmt::Display for FeedbackError<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Source(e) => write!(f, "thermocreep source: {e}"),
            Self::Mechanical(e) => write!(f, "thermocreep mechanics: {e}"),
            Self::Caloric(e) => write!(f, "thermocreep caloric state: {e}"),
            Self::Transport(e) => write!(f, "thermocreep transport: {e}"),
            Self::Invalid(what) => write!(f, "thermocreep: {what}"),
            Self::Cancelled => f.write_str("thermocreep: cancelled"),
            Self::NotConverged {
                iterations,
                temperature_defect_k,
                heat_defect_j,
                thermal_residual_j,
                combined_residual_j,
            } => write!(
                f,
                "thermocreep did not converge in {iterations} proposals: dT={temperature_defect_k} K, dQ={heat_defect_j} J, thermal={thermal_residual_j} J, combined={combined_residual_j} J"
            ),
        }
    }
}
impl<E: std::error::Error + 'static> std::error::Error for FeedbackError<E> {}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn require<E>(ok: bool, reason: &'static str) -> Result<(), FeedbackError<E>> {
    if ok {
        Ok(())
    } else {
        Err(FeedbackError::Invalid(reason))
    }
}
fn poll<E>(cancelled: &mut impl FnMut() -> bool) -> Result<(), FeedbackError<E>> {
    if cancelled() {
        Err(FeedbackError::Cancelled)
    } else {
        Ok(())
    }
}

/// The type name explicitly selects zero thermal expansion, constant reference
/// density and the source family's constant elastic storage. This selection is
/// a reduction of the physical model, not an inferred material-property value.
#[derive(Clone, Debug)]
pub struct ZeroExpansionCaloricPoint<'a> {
    source: NortonCreepSource<'a>,
    caloric: &'a ResolvedHeatCapacityCurve,
    strain: Voigt,
    stress_pa: Voigt,
    history: CreepState,
    state: EquilibriumPhaseState,
    mass_kg: f64,
    volume_m3: f64,
    surface_area_m2: f64,
    time_s: f64,
}

/// Proposed specimen and actual transport result must be accepted together.
#[derive(Debug)]
pub struct FeedbackStep<'a, R> {
    pub next: ZeroExpansionCaloricPoint<'a>,
    /// Fixed-temperature native tangent; not a coupled thermomechanical tangent.
    pub mechanical: CreepIncrement,
    pub transport: UniformEnthalpyStep<R>,
    pub iterations: usize,
    pub physical_heat_j: f64,
    pub numerical_loss_j: f64,
    pub applied_heat_j: f64,
    pub temperature_defect_k: f64,
    pub heat_defect_j: f64,
    pub thermal_residual_j: f64,
    pub combined_residual_j: f64,
}

impl<'a> ZeroExpansionCaloricPoint<'a> {
    /// Bind source/caloric identities and all source conditions before creating
    /// history. The caloric temperature interval must be a subset of the creep
    /// interval; every original density knot must have the same positive value.
    pub fn new(
        card: &'a MaterialCard,
        source: NortonCreepSource<'a>,
        caloric: &'a ResolvedHeatCapacityCurve,
        initial: InitialState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FeedbackError> {
        poll(&mut cancelled)?;
        let density = admit_zero_expansion_caloric(card, &source, caloric, &mut cancelled)?;
        require(
            positive(initial.mass_kg)
                && positive(initial.surface_area_m2)
                && initial.time_s.is_finite()
                && initial.time_s >= 0.0,
            "mass/area must be positive and the initial clock nonnegative",
        )?;
        let volume_m3 = initial.mass_kg / density;
        require(positive(volume_m3), "reference volume is not representable")?;
        let state = caloric
            .curve()
            .state_at_specific_enthalpy(initial.specific_enthalpy_j_kg)
            .map_err(FeedbackError::Caloric)?;
        let bound = source
            .at_temperature(state.temperature_k(), 1.0)
            .map_err(FeedbackError::Source)?;
        bound
            .law()
            .validate_history(&initial.strain, &initial.history)
            .map_err(FeedbackError::Mechanical)?;
        // A birth has no physical duration. Reconstruct accepted elastic stress,
        // never call a duration-bearing creep trial and silently age the state.
        let elastic_strain =
            core::array::from_fn(|i| initial.strain[i] - initial.history.creep_strain[i]);
        let stress_pa = bound.law().elastic().stress(&elastic_strain, &());
        require(
            stress_pa.iter().all(|s| s.is_finite()),
            "initial elastic stress is not representable",
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            source,
            caloric,
            strain: initial.strain,
            stress_pa,
            history: initial.history,
            state,
            mass_kg: initial.mass_kg,
            volume_m3,
            surface_area_m2: initial.surface_area_m2,
            time_s: initial.time_s,
        })
    }
    pub fn source(&self) -> &NortonCreepSource<'a> {
        &self.source
    }
    pub fn history(&self) -> &CreepState {
        &self.history
    }
    pub const fn strain(&self) -> &Voigt {
        &self.strain
    }
    pub const fn stress_pa(&self) -> &Voigt {
        &self.stress_pa
    }
    pub const fn caloric_state(&self) -> &EquilibriumPhaseState {
        &self.state
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }

    /// Reuse the caloric owner's signed-energy operation for an insulated body.
    pub fn advance_adiabatic(
        &self,
        strain: Voigt,
        duration_s: f64,
        config: FeedbackConfig,
        cancelled: impl FnMut() -> bool,
    ) -> Result<FeedbackStep<'a, ()>, FeedbackError<PhaseStateError>> {
        self.advance(
            strain,
            duration_s,
            config,
            |input| {
                let state = input.curve.advance_specific_energy(
                    input.initial,
                    input.internal_heat_j / input.mass_kg,
                )?;
                Ok(UniformEnthalpyStep {
                    state,
                    external_heat_j: 0.0,
                    energy_residual_tolerance_j: config.energy_tolerance_j,
                    report: (),
                })
            },
            cancelled,
        )
    }

    /// Advance prescribed native total strain and physical time. A transport
    /// callback proposes heat exchange; it must not debit a finite reservoir
    /// until this whole transaction succeeds. History is immutable on all paths.
    pub fn advance<R, E>(
        &self,
        strain: Voigt,
        duration_s: f64,
        config: FeedbackConfig,
        mut transport: impl FnMut(UniformEnthalpyStepInput<'_>) -> Result<UniformEnthalpyStep<R>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FeedbackStep<'a, R>, FeedbackError<E>> {
        poll(&mut cancelled)?;
        let final_time = self.time_s + duration_s;
        let material_tolerance = config.energy_tolerance_j / self.volume_m3;
        require(
            config.max_iterations > 0
                && config.max_iterations <= MAX_FEEDBACK_ITERATIONS
                && positive(config.temperature_tolerance_k)
                && positive(config.energy_tolerance_j)
                && positive(config.stress_tolerance_pa)
                && positive(config.relaxation)
                && config.relaxation <= 1.0
                && positive(duration_s)
                && positive(material_tolerance)
                && final_time.is_finite()
                && final_time > self.time_s,
            "invalid feedback controls or physical interval",
        )?;
        let tolerances = CreepTolerances {
            stress_pa: config.stress_tolerance_pa,
            energy_j_m3: material_tolerance,
        };
        let propose = |temperature| {
            let bound = self
                .source
                .at_temperature(temperature, duration_s)
                .map_err(FeedbackError::Source)?;
            bound
                .law()
                .advance_with_energy(&self.strain, &self.history, &strain, tolerances)
                .map_err(FeedbackError::Mechanical)
        };
        let mut temperature = self.state.temperature_k();
        let mut defects = [f64::INFINITY; 4];
        for iterations in 1..=config.max_iterations {
            poll(&mut cancelled)?;
            let trial = propose(temperature)?;
            let applied_heat_j = self.volume_m3 * trial.energy.physical_dissipation_j_m3;
            require(
                applied_heat_j.is_finite(),
                "physical heat is not representable",
            )?;
            let proposal = transport(UniformEnthalpyStepInput {
                curve: self.caloric.curve(),
                initial: self.state,
                mass_kg: self.mass_kg,
                volume_m3: self.volume_m3,
                surface_area_m2: self.surface_area_m2,
                internal_heat_j: applied_heat_j,
                duration_s,
            })
            .map_err(FeedbackError::Transport)?;
            poll(&mut cancelled)?;
            let checked = self
                .caloric
                .curve()
                .state_at_specific_enthalpy(proposal.state.specific_enthalpy_j_kg())
                .map_err(FeedbackError::Caloric)?;
            require(
                checked == proposal.state
                    && proposal.external_heat_j.is_finite()
                    && proposal.energy_residual_tolerance_j.is_finite()
                    && proposal.energy_residual_tolerance_j >= 0.0
                    && proposal.energy_residual_tolerance_j <= config.energy_tolerance_j,
                "transport returned a foreign state or invalid energy allowance",
            )?;
            let thermal_change = self.mass_kg
                * (proposal.state.specific_enthalpy_j_kg() - self.state.specific_enthalpy_j_kg());
            let transport_residual = thermal_change - applied_heat_j - proposal.external_heat_j;
            require(
                transport_residual.is_finite()
                    && transport_residual.abs() <= proposal.energy_residual_tolerance_j,
                "transport violates its declared heat balance",
            )?;
            let final_temperature = proposal.state.temperature_k();
            let mechanical = propose(final_temperature)?;
            let physical_heat_j = self.volume_m3 * mechanical.energy.physical_dissipation_j_m3;
            let numerical_loss_j = self.volume_m3 * mechanical.energy.numerical_dissipation_j_m3;
            let thermal_residual = thermal_change - physical_heat_j - proposal.external_heat_j;
            let combined_residual = self.volume_m3
                * (mechanical.energy.final_stored_j_m3
                    - mechanical.energy.initial_stored_j_m3
                    - mechanical.energy.mechanical_work_j_m3)
                + thermal_change
                + numerical_loss_j
                - proposal.external_heat_j;
            defects = [
                (final_temperature - temperature).abs(),
                (physical_heat_j - applied_heat_j).abs(),
                thermal_residual,
                combined_residual,
            ];
            require(
                defects.iter().all(|x| x.is_finite()) && numerical_loss_j.is_finite(),
                "combined heat/work account is not representable",
            )?;
            if defects[0] <= config.temperature_tolerance_k
                && defects[1] <= config.energy_tolerance_j
                && defects[2].abs() <= config.energy_tolerance_j
                && defects[3].abs() <= config.energy_tolerance_j
            {
                poll(&mut cancelled)?;
                let next = Self {
                    source: self.source.clone(),
                    caloric: self.caloric,
                    strain,
                    stress_pa: mechanical.stress_pa,
                    history: mechanical.state.clone(),
                    state: proposal.state,
                    mass_kg: self.mass_kg,
                    volume_m3: self.volume_m3,
                    surface_area_m2: self.surface_area_m2,
                    time_s: final_time,
                };
                return Ok(FeedbackStep {
                    next,
                    mechanical,
                    transport: proposal,
                    iterations,
                    physical_heat_j,
                    numerical_loss_j,
                    applied_heat_j,
                    temperature_defect_k: defects[0],
                    heat_defect_j: defects[1],
                    thermal_residual_j: defects[2],
                    combined_residual_j: defects[3],
                });
            }
            temperature =
                (1.0 - config.relaxation) * temperature + config.relaxation * final_temperature;
        }
        poll(&mut cancelled)?;
        Err(FeedbackError::NotConverged {
            iterations: config.max_iterations,
            temperature_defect_k: defects[0],
            heat_defect_j: defects[1],
            thermal_residual_j: defects[2],
            combined_residual_j: defects[3],
        })
    }
}

/// Admit the original caloric and creep sources for the explicit zero-expansion
/// reduction, returning their constant reference density [kg/m³]. Uniform and
/// spatial consumers share this same-card, fixed-context, pinned Cp/density and
/// fully solid temperature-domain admission; absent expansion data is not zero.
pub fn admit_zero_expansion_caloric(
    card: &MaterialCard,
    source: &NortonCreepSource<'_>,
    caloric: &ResolvedHeatCapacityCurve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, FeedbackError> {
    let curve = caloric.curve();
    let knots = curve.knots();
    let range = source.parameters().temperature_range_k;
    require(
        curve.material_card_identity() == card.content_hash()
            && source.reference().card_identity() == card.content_hash()
            && knots.len() >= 2
            && knots.len() <= crate::phase::MAX_HEAT_CAPACITY_ENTHALPY_KNOTS
            && knots[0].temperature_k >= range[0]
            && knots[knots.len() - 1].temperature_k <= range[1],
        "caloric and creep laws require one card and a supported bounded thermal interval",
    )?;
    let density = knots[0].bulk_density_kg_m3;
    require(
        positive(density)
            && knots
                .iter()
                .all(|k| k.liquid_mass_fraction == 0.0 && k.bulk_density_kg_m3 == density)
            && knots
                .windows(2)
                .all(|p| p[0].temperature_k < p[1].temperature_k),
        "zero-expansion feedback requires a strictly increasing fixed-density solid chart",
    )?;
    let receipts = caloric.receipts();
    require(
        receipts.len() >= 4 && receipts.len() % 2 == 0,
        "caloric chart lacks paired source Cp/density receipts",
    )?;
    let context = source
        .at_temperature(knots[0].temperature_k, 1.0)
        .map_err(FeedbackError::Source)?;
    let point = context.material().query();
    let axis = source.temperature_axis();
    for (index, receipt) in receipts.iter().enumerate() {
        if index % 128 == 0 {
            poll(cancelled)?;
        }
        require(
            receipt.axis_quantities == *context.material().axis_quantities()
                && receipt
                    .query_point
                    .iter()
                    .filter(|(name, _)| name != axis)
                    .eq(context
                        .material()
                        .query_point()
                        .iter()
                        .filter(|(name, _)| name != axis))
                && receipt.selected == receipts[index % 2].selected,
            "caloric and mechanical sources must share all fixed conditions and pinned claims",
        )?;
    }
    for (index, dims) in [Dims([2, 0, -2, -1, 0, 0]), Dims([-3, 1, 0, 0, 0, 0])]
        .into_iter()
        .enumerate()
    {
        poll(cancelled)?;
        let claim = card
            .claims()
            .claim(receipts[index].selected)
            .ok_or(FeedbackError::Invalid(
                "caloric source claim is absent from the original card",
            ))?;
        require(
            claim.key.quantity() == QuantitySpec::dimensional(dims),
            "caloric source has incorrect units",
        )?;
        let requirement = ScalarPropertyRequirement::try_with_key(
            &claim.key,
            ScalarAdmissibility::StrictlyPositive,
        )
        .map_err(|e| FeedbackError::Source(e.into()))?;
        let resolved = resolve_material_state_point(
            card,
            &point,
            &[requirement],
            MaterialPropertySelection::PinnedByProperty(vec![(
                claim.key.name().to_owned(),
                receipts[index].selected,
            )]),
        )
        .map_err(|e| FeedbackError::Source(e.into()))?;
        require(
            resolved.properties()[0].has_only_fixed_external_conditions(Some(axis), false),
            "caloric coefficients cannot freeze an unowned evolving material coordinate",
        )?;
    }
    Ok(density)
}
