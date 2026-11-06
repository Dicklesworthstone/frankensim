//! Finite caloric feedback for an original hardness-limited sliding contact.
//!
//! The existing massless-interface heat residual solves one interface temperature
//! and the original finite-cell energies. Every nonlinear trial refreshes the
//! SAME accepted material/contact history to those trial temperatures, advances
//! Coulomb sliding and Archard removal from that refreshed pressure, and then
//! re-equilibrates the worn normal contact. Thermal coefficients are implicit
//! endpoint values; wear geometry uses the native beginning-pressure split and
//! requires physical step refinement. This is not a fully implicit wear solve.
//!
//! Two uniform bulk temperatures and one interface temperature cover the native
//! active cells. The default requires a fixed active footprint. The evolving-area
//! constructor admits recruitment and opening on the original grid: final worn
//! endpoint area supplies actual half-cell conductance and ambient coverage,
//! while friction/wear still use the native refreshed beginning pressure. This
//! first-order split requires time and grid refinement; discontinuous cell
//! activation can prevent nonlinear convergence. Complete detachment refuses.
//! The thermal masses and geometry stay fixed. A cumulative removal fraction bounds
//! this small-removal approximation; removed material carries no modeled mass,
//! enthalpy, or debris transport. Zero permitted removal gives a no-wear rung.
//!
//! Native friction dissipation heats the interface. An optional explicit
//! plastic-work partition deposits heat directly into the two original bodies;
//! the unconverted work remains a separate nonthermal account. No universal
//! Taylor--Quinney fraction or assignment of relative indentation to either
//! material is inferred. Signed temperature/geometry parameter work remains
//! separate, with a bound on the SUM of absolute parameter-work supplies.
//! There is no reciprocal thermoelastic free-energy model.
//! Caller-owned tangent compliance supplies no body storage here: local contact
//! identities and the fixed-inventory thermal balance do not establish a full
//! body-energy, entropy, or passivity theorem.

use core::convert::Infallible;
use std::cell::Cell;

use fs_contact::normal_patch::plastic_source::wear::sliding::temperature::{
    SourcePlasticSlidingThermalReceipt, SourcePlasticSlidingThermalRefresh,
};
use fs_contact::normal_patch::plastic_source::wear::sliding::{
    PlasticCellSliding, PlasticSlidingInput, SourcePlasticSlidingError, SourcePlasticSlidingWear,
    SourcePlasticSlidingWearAdvance, SourcePlasticSlidingWearState,
};
use fs_contact::normal_patch::{
    FiniteGapPlasticContinuationReceipt, FiniteGapPlasticDrive, FiniteGapPlasticDriveReceipt,
    FiniteGapPlasticLoadControls, FiniteGapPlasticLoadReceipt,
};
use fs_exec::Cx;
use fs_matdb::{MaterialCard, PropertyKey};
use fs_material::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use fs_material::state_point::{DENSITY_PROPERTY, ResolvedMaterialStatePoint};
use fs_qty::Density;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};
use fs_tribo::TangentialSolveConfig;

use super::{
    FrictionHeatInput, FrictionHeatReport, FrictionHeatRoute, FrictionHeatStep,
    FrictionHeatTolerances,
};

/// Physical limits of the fixed-inventory caloric approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlasticSlidingHeatLimits {
    /// Cumulative removed volume divided by each ORIGINAL receiving cell volume.
    /// Each bound is in [0, 1); zero admits only exactly zero cumulative removal.
    /// This bounds an omitted effect, not its thermal error or physical validity.
    pub maximum_removed_volume_fraction: [f64; 2],
    /// Sum of absolute temperature and recession parameter work per step [J].
    /// These are explicit external work accounts, never friction heat sources.
    pub maximum_absolute_parameter_work_j: f64,
}

/// Declared fractions of incremental normal plastic work heating bodies A/B.
///
/// Fractions are nonnegative and sum to at most one. Their remainder represents
/// nonthermal plastic work, not heat or numerical damping. This is an explicit
/// constitutive/partition assumption, not an inference from hardness, a material
/// calibration, or a resolved subsurface plastic-zone temperature model.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PlasticWorkHeatPartition {
    body_fractions: [f64; 2],
}

impl PlasticWorkHeatPartition {
    /// Declare both receiving-body fractions. The default converts no work.
    pub fn new(body_fractions: [f64; 2]) -> Result<Self, PlasticSlidingHeatError> {
        if body_fractions
            .iter()
            .any(|fraction| !fraction.is_finite() || !(0.0..=1.0).contains(fraction))
            || body_fractions[0] + body_fractions[1] > 1.0
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "plastic heat fractions must be finite, nonnegative and sum to at most one",
            ));
        }
        Ok(Self { body_fractions })
    }

    /// Fractions in the original ordered material/body roles.
    #[must_use]
    pub const fn body_fractions(self) -> [f64; 2] {
        self.body_fractions
    }

    fn split<E>(self, work_j: f64) -> Result<([f64; 2], f64), PlasticSlidingHeatError<E>> {
        if !work_j.is_finite() || work_j < 0.0 {
            return Err(PlasticSlidingHeatError::Invalid(
                "only finite nonnegative incremental plastic work can supply heat",
            ));
        }
        let heat_j = self.body_fractions.map(|fraction| fraction * work_j);
        let nonthermal_j = (1.0 - (self.body_fractions[0] + self.body_fractions[1])) * work_j;
        Ok((heat_j, nonthermal_j))
    }
}

/// Normal boundary retained during both temperature refresh and sliding/removal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlasticSlidingHeatBoundary {
    HeldApproach,
    MaintainedLoad {
        force_n: f64,
        controls: FiniteGapPlasticLoadControls,
    },
    FixedDrive {
        drive: FiniteGapPlasticDrive,
        controls: FiniteGapPlasticLoadControls,
    },
}

/// Spatial scope of the two uniform finite caloric receivers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlasticSlidingHeatFootprint {
    /// Every originally loaded cell stays loaded, and no other cell loads.
    FixedActiveCells,
    /// Loaded cells may change within the same grid on both original faces.
    /// Heat uses final endpoint area, with no swept residence or lateral field.
    EvolvingActiveCells,
}

/// Independent numerical work and contact/thermal residual limits.
#[derive(Clone, Copy, Debug)]
pub struct PlasticSlidingHeatControls {
    pub friction: TangentialSolveConfig,
    pub maximum_cells: usize,
    pub maximum_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Numerical finite-difference perturbation in energy coordinates [J].
    pub difference_step_j: f64,
    /// Scaling of the interface-temperature unknown [J/K], not a heat capacity.
    pub temperature_scale_j_per_k: f64,
    pub heat: FrictionHeatTolerances,
    /// Bound for each native normal identity and actual tangential work [J].
    pub contact_work_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

/// Endpoint conductivity and OTHER heat/conduction from the physical source owner.
pub struct PlasticSlidingHeatSource<P> {
    /// Ordered original bulk conductivities at their own endpoint temperatures.
    /// The callback must retain its actual query/selection evidence in `payload`;
    /// this adapter neither infers conductivity from hardness nor selects claims.
    pub conductivity_w_m_k: [f64; 2],
    /// Excludes friction heat, the configured normal plastic heat, and the
    /// original route's direct conductance. Each physical source is supplied once.
    pub transport: StepCoefficients<P>,
}

#[derive(Debug)]
pub enum PlasticSlidingHeatError<E = Infallible> {
    Invalid(&'static str),
    Contact(SourcePlasticSlidingError),
    Phase(PhaseStateError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    TrialBudget,
}
impl<E> From<SourcePlasticSlidingError> for PlasticSlidingHeatError<E> {
    fn from(error: SourcePlasticSlidingError) -> Self {
        Self::Contact(error)
    }
}
impl<E> From<PhaseStateError> for PlasticSlidingHeatError<E> {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl<E> From<DepositionError> for PlasticSlidingHeatError<E> {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}
impl<E> From<NewtonError> for PlasticSlidingHeatError<E> {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for PlasticSlidingHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hardness-limited sliding/finite heat: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for PlasticSlidingHeatError<E> {}

/// The two ORIGINAL native operations; their histories and work remain inspectable.
#[derive(Clone, Debug, PartialEq)]
pub struct PlasticSlidingHeatResponse<B = ()> {
    pub temperature: SourcePlasticSlidingThermalReceipt<B>,
    pub sliding: SourcePlasticSlidingWearAdvance<B>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum PlasticSlidingHeatMechanics {
    HeldApproach(PlasticSlidingHeatResponse),
    MaintainedLoad(PlasticSlidingHeatResponse<FiniteGapPlasticLoadReceipt>),
    FixedDrive(PlasticSlidingHeatResponse<FiniteGapPlasticDriveReceipt>),
}

impl PlasticSlidingHeatMechanics {
    pub fn temperature_normal(&self) -> &FiniteGapPlasticContinuationReceipt {
        match self {
            Self::HeldApproach(r) => &r.temperature.normal,
            Self::MaintainedLoad(r) => &r.temperature.normal,
            Self::FixedDrive(r) => &r.temperature.normal,
        }
    }
    pub fn wear_normal(&self) -> &FiniteGapPlasticContinuationReceipt {
        match self {
            Self::HeldApproach(r) => &r.sliding.normal,
            Self::MaintainedLoad(r) => &r.sliding.normal,
            Self::FixedDrive(r) => &r.sliding.normal,
        }
    }
    pub fn candidate(&self) -> &SourcePlasticSlidingWearState {
        match self {
            Self::HeldApproach(r) => &r.sliding.candidate,
            Self::MaintainedLoad(r) => &r.sliding.candidate,
            Self::FixedDrive(r) => &r.sliding.candidate,
        }
    }
    pub fn frictional_work_j(&self) -> f64 {
        match self {
            Self::HeldApproach(r) => r.sliding.frictional_work_j,
            Self::MaintainedLoad(r) => r.sliding.frictional_work_j,
            Self::FixedDrive(r) => r.sliding.frictional_work_j,
        }
    }
    pub fn cells(&self) -> &[PlasticCellSliding] {
        match self {
            Self::HeldApproach(r) => &r.sliding.cells,
            Self::MaintainedLoad(r) => &r.sliding.cells,
            Self::FixedDrive(r) => &r.sliding.cells,
        }
    }
}

/// Accepted original normal/shear/wear history and its immutable caloric binding.
#[derive(Clone)]
pub struct PlasticSlidingHeatContact<'a> {
    source: SourcePlasticSlidingWear,
    state: SourcePlasticSlidingWearState,
    material_cards: [&'a MaterialCard; 2],
    caloric: [&'a ResolvedHeatCapacityCurve; 2],
    route: FrictionHeatRoute,
    original_body_ids: [u64; 2],
    original_volume_m3: [f64; 2],
    active_cells: Vec<bool>,
    footprint: PlasticSlidingHeatFootprint,
    limits: PlasticSlidingHeatLimits,
    plastic_work_heat: PlasticWorkHeatPartition,
    time_s: f64,
}

pub struct PlasticSlidingHeatReport<P> {
    pub mechanics: PlasticSlidingHeatMechanics,
    pub heat: FrictionHeatReport,
    /// Exact original conductivity/transport callback payload at the accepted root.
    pub source: P,
    /// Exact original material conductivity supplied at the accepted thermal root.
    pub conductivity_w_m_k: [f64; 2],
    pub previous_contact_area_m2: f64,
    /// Native temperature-refreshed pressure area before this interval's wear.
    pub refreshed_contact_area_m2: f64,
    /// Final native worn pressure area; the evolving route uses it for heat.
    pub endpoint_contact_area_m2: f64,
    pub nonlinear: NewtonReport,
    pub trials: usize,
    /// Largest absolute local normal or actual tangential work residual [J].
    /// This does not audit the bodies behind caller-supplied tangent compliance.
    pub contact_work_residual_j: f64,
    pub frictional_work_j: f64,
    pub temperature_parameter_work_j: f64,
    pub wear_parameter_work_j: f64,
    /// Incremental normal plastic work, before its explicit heat partition [J].
    pub normal_plastic_work_j: f64,
    /// Incremental plastic heat deposited directly into original bodies A/B [J].
    /// Friction remains at the interface; these body sources are not added to
    /// `heat.dissipated_j` or `heat.body_heat_j` a second time.
    pub normal_plastic_heat_j: [f64; 2],
    /// Plastic work left outside the caloric model by the declared partition [J].
    pub nonthermal_plastic_work_j: f64,
    /// Authored/environment heat excluding BOTH friction and normal plastic heat.
    /// The lower-level `heat.external_heat_j` includes the direct plastic sources
    /// because they are external to that thermal-only accounting boundary.
    pub other_heat_j: f64,
    /// Normal boundary work into contact, including both native operations [J].
    pub normal_work_j: f64,
    /// Cumulative original A/B removed volumes, with no thermal mass export [m^3].
    pub removed_volume_m3: [f64; 2],
    /// The same cumulative volumes divided by the ORIGINAL receiving cell volumes.
    pub removed_volume_fraction: [f64; 2],
}

pub struct PlasticSlidingHeatCandidate<'a, 't, P> {
    contact: PlasticSlidingHeatContact<'a>,
    thermal: DepositionSimulation<'t>,
    report: PlasticSlidingHeatReport<P>,
}
impl<'a, 't, P> PlasticSlidingHeatCandidate<'a, 't, P> {
    pub const fn contact(&self) -> &PlasticSlidingHeatContact<'a> {
        &self.contact
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn report(&self) -> &PlasticSlidingHeatReport<P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        PlasticSlidingHeatContact<'a>,
        DepositionSimulation<'t>,
        PlasticSlidingHeatReport<P>,
    ) {
        (self.contact, self.thermal, self.report)
    }
}

struct PreparedMechanics {
    source: SourcePlasticSlidingWear,
    mechanics: PlasticSlidingHeatMechanics,
    contact_work_residual_j: f64,
    removed_volume_fraction: [f64; 2],
    refreshed_contact_area_m2: f64,
    endpoint_contact_area_m2: f64,
    normal_plastic_heat_j: [f64; 2],
    nonthermal_plastic_work_j: f64,
}

struct Staged<P> {
    prepared: PreparedMechanics,
    payload: P,
    conductivity_w_m_k: [f64; 2],
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), DepositionError> {
    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
}
fn finite<E>(value: f64) -> Result<f64, PlasticSlidingHeatError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic.into())
    }
}
fn numerical<E>(
    error: FiniteDifferenceError<PlasticSlidingHeatError<E>>,
) -> PlasticSlidingHeatError<E> {
    match error {
        FiniteDifferenceError::Source(error) => error,
        FiniteDifferenceError::Invalid(message) => PlasticSlidingHeatError::Invalid(message),
    }
}

impl<'a> PlasticSlidingHeatContact<'a> {
    /// Bind already accepted contact and heat histories to their original sources.
    ///
    /// The exact Cp/density charts used by the thermal cells must belong to the
    /// original ordered cards, fixed conditions, and selected density claims.
    /// The caller declares this fixed-footprint, small-removal approximation;
    /// a valid binding does not establish its experimental applicability.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        source: SourcePlasticSlidingWear,
        state: SourcePlasticSlidingWearState,
        material_cards: [&'a MaterialCard; 2],
        caloric: [&'a ResolvedHeatCapacityCurve; 2],
        route: &FrictionHeatRoute,
        limits: PlasticSlidingHeatLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, PlasticSlidingHeatError> {
        Self::new_with_footprint(
            thermal,
            source,
            state,
            material_cards,
            caloric,
            route,
            limits,
            PlasticSlidingHeatFootprint::FixedActiveCells,
            cx,
        )
    }

    /// Bind changing loaded cells on the same projected finite-gap grid.
    ///
    /// The entire original grid must fit both uniform prismatic thermal faces;
    /// every current positive-pressure cell contributes its actual grid area.
    /// No pressure/hardness area surrogate or effective conductivity is used.
    /// The shared uniform-body model resolves total area, not lateral temperature
    /// gradients, contact residence, constriction resistance, or cell-edge events.
    #[allow(clippy::too_many_arguments)]
    pub fn new_evolving_area(
        thermal: &DepositionSimulation<'_>,
        source: SourcePlasticSlidingWear,
        state: SourcePlasticSlidingWearState,
        material_cards: [&'a MaterialCard; 2],
        caloric: [&'a ResolvedHeatCapacityCurve; 2],
        route: &FrictionHeatRoute,
        limits: PlasticSlidingHeatLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, PlasticSlidingHeatError> {
        Self::new_with_footprint(
            thermal,
            source,
            state,
            material_cards,
            caloric,
            route,
            limits,
            PlasticSlidingHeatFootprint::EvolvingActiveCells,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_with_footprint(
        thermal: &DepositionSimulation<'_>,
        source: SourcePlasticSlidingWear,
        state: SourcePlasticSlidingWearState,
        material_cards: [&'a MaterialCard; 2],
        caloric: [&'a ResolvedHeatCapacityCurve; 2],
        route: &FrictionHeatRoute,
        limits: PlasticSlidingHeatLimits,
        footprint: PlasticSlidingHeatFootprint,
        cx: &Cx<'_>,
    ) -> Result<Self, PlasticSlidingHeatError> {
        checkpoint(cx)?;
        if limits
            .maximum_removed_volume_fraction
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0 || *v >= 1.0)
            || !limits.maximum_absolute_parameter_work_j.is_finite()
            || limits.maximum_absolute_parameter_work_j < 0.0
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "finite nonnegative wear fractions below one and parameter-work budget required",
            ));
        }
        route.validate(thermal)?;
        let cells = route.cells();
        if thermal.cells()[&cells[0]].body_id == thermal.cells()[&cells[1]].body_id {
            return Err(PlasticSlidingHeatError::Invalid(
                "thermal destinations must retain two distinct original bodies",
            ));
        }
        let normal = source.wear_source().contact().normal_source().elastic();
        let bulk = [normal.surface_a_state(), normal.surface_b_state()];
        for i in 0..2 {
            let resolved = bulk[i].resolved();
            let density_name = bulk[i]
                .elastic_state()
                .derivation()
                .map_or(DENSITY_PROPERTY, |d| d.density_requirement().name());
            let density = resolved
                .property_by_key(&PropertyKey::new(density_name, Density::DIMS))
                .ok_or(PlasticSlidingHeatError::Invalid(
                    "original normal density source is missing",
                ))?;
            if material_cards[i].content_hash() != resolved.card_identity()
                || caloric[i].curve().material_card_identity() != resolved.card_identity()
                || caloric[i].receipts().is_empty()
            {
                return Err(PlasticSlidingHeatError::Invalid(
                    "heat and hardness must retain the original ordered material cards and sourced calorics",
                ));
            }
            for (index, receipt) in caloric[i].receipts().iter().enumerate() {
                if index % 64 == 0 {
                    checkpoint(cx)?;
                }
                if receipt.axis_quantities != *resolved.axis_quantities()
                    || receipt
                        .query_point
                        .iter()
                        .filter(|(axis, _)| axis != "T")
                        .ne(resolved
                            .query_point()
                            .iter()
                            .filter(|(axis, _)| axis != "T"))
                    || (index % 2 == 1
                        && (receipt.property != density_name
                            || receipt.selected != density.answer().receipt.selected))
                {
                    return Err(PlasticSlidingHeatError::Invalid(
                        "caloric and normal sources require identical fixed conditions, quantity axes and density claims",
                    ));
                }
            }
        }
        let temperature = temperatures(thermal, caloric, route)?;
        if (0..2).any(|i| bulk[i].resolved().state_coordinate("T") != Some(temperature[i])) {
            return Err(PlasticSlidingHeatError::Invalid(
                "accepted bulk source temperatures must equal the original thermal checkpoint",
            ));
        }
        // This original-source no-op performs the native opaque-state ownership
        // admission without reconstructing or relabeling a physical history.
        let admitted = source.refresh_at_temperatures(
            &state,
            material_cards,
            temperature,
            source.interface_temperature_k(),
            cx,
        )?;
        if admitted.source() != &source || admitted.candidate() != &state {
            return Err(PlasticSlidingHeatError::Invalid(
                "thermal binding must not change accepted source or physical history",
            ));
        }
        let pressure = state
            .wear_state()
            .contact_state()
            .contact_state()
            .pressure_pa();
        let mut active_cells = Vec::with_capacity(pressure.len());
        let mut count = 0usize;
        for (i, p) in pressure.iter().enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            let active = *p > 0.0;
            active_cells.push(active);
            count += usize::from(active);
        }
        let grid = source.wear_source().contact().model().config().grid;
        let area = finite(count as f64 * (grid.cell_width_m * grid.cell_depth_m))?;
        let thermal_area = route.contact_area_m2();
        if count == 0 || (area - thermal_area).abs() > 64.0 * f64::EPSILON * area.max(thermal_area)
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "actual loaded grid cells must exactly cover the original thermal face area",
            ));
        }
        let bound = Self {
            source,
            state,
            material_cards,
            caloric,
            route: route.clone(),
            original_body_ids: cells.map(|id| thermal.cells()[&id].body_id),
            original_volume_m3: cells.map(|id| thermal.cells()[&id].geometry.volume_m3),
            active_cells,
            footprint,
            limits,
            plastic_work_heat: PlasticWorkHeatPartition::default(),
            time_s: thermal.time_s(),
        };
        bound.admit_footprint_support()?;
        bound.admit_removal(&bound.state)?;
        checkpoint(cx)?;
        Ok(bound)
    }

    pub const fn source(&self) -> &SourcePlasticSlidingWear {
        &self.source
    }
    pub const fn state(&self) -> &SourcePlasticSlidingWearState {
        &self.state
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn limits(&self) -> PlasticSlidingHeatLimits {
        self.limits
    }
    pub const fn route(&self) -> &FrictionHeatRoute {
        &self.route
    }
    pub const fn footprint(&self) -> PlasticSlidingHeatFootprint {
        self.footprint
    }

    /// Configure conversion of FUTURE incremental indentation work into heat.
    ///
    /// Accepted cumulative work is never reheated. The same fractions apply to
    /// every thermal trial and survive accepted steps. No plastic depth is
    /// assigned to a body: the caller explicitly owns this energy partition.
    #[must_use]
    pub fn with_plastic_work_heat(mut self, partition: PlasticWorkHeatPartition) -> Self {
        self.plastic_work_heat = partition;
        self
    }

    /// Configured conversion/partition of future normal plastic work.
    pub const fn plastic_work_heat(&self) -> PlasticWorkHeatPartition {
        self.plastic_work_heat
    }

    /// Propose one heat/contact interval from the exact accepted histories.
    ///
    /// The existing heat owner supplies cell-energy and massless-interface
    /// residuals; the existing Newton--Krylov solver supplies iteration. Every
    /// residual, Jacobian probe, and final evaluation starts from `self`, never
    /// from a previous trial. Endpoint temperatures refresh original sources
    /// before the native beginning-pressure sliding/wear split. The callback
    /// supplies only original conductivity and OTHER heat/conduction.
    ///
    /// Return/adopt the contact and thermal candidates together. A source,
    /// footprint, cumulative-removal, parameter-work, numerical, or cancellation
    /// refusal returns no partial physical step. Cut intervals at thermal events.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_candidate<'t, P, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        motion: &[PlasticSlidingInput],
        boundary: PlasticSlidingHeatBoundary,
        controls: PlasticSlidingHeatControls,
        energy_guess_j: &[f64],
        interface_temperature_guess_k: f64,
        resolve: impl Fn(&DepositionSimulation<'t>, f64) -> Result<PlasticSlidingHeatSource<P>, E>,
        cx: &Cx<'_>,
    ) -> Result<PlasticSlidingHeatCandidate<'a, 't, P>, PlasticSlidingHeatError<E>> {
        checkpoint(cx)?;
        self.admit_thermal(thermal)?;
        let previous_contact_area_m2 = self.admit_footprint(&self.state, cx)?;
        let bound_area = self.route.contact_area_m2();
        if (previous_contact_area_m2 - bound_area).abs()
            > 64.0 * f64::EPSILON * previous_contact_area_m2.max(bound_area)
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "accepted native pressure area differs from its thermal contact",
            ));
        }
        controls.newton.validate()?;
        let cells = thermal.cells().len();
        let dimension = cells
            .checked_add(1)
            .ok_or(PlasticSlidingHeatError::Invalid(
                "thermal dimension overflow",
            ))?;
        if controls.maximum_cells == 0
            || cells > controls.maximum_cells
            || controls.maximum_trials == 0
            || controls.maximum_newton_iterations == 0
            || !controls.difference_step_j.is_finite()
            || controls.difference_step_j <= 0.0
            || !controls.temperature_scale_j_per_k.is_finite()
            || controls.temperature_scale_j_per_k <= 0.0
            || !interface_temperature_guess_k.is_finite()
            || interface_temperature_guess_k <= 0.0
            || energy_guess_j.len() != cells
            || energy_guess_j.iter().any(|v| !v.is_finite())
            || [
                controls.heat.cell_j,
                controls.heat.thermal_j,
                controls.heat.interface_j,
                controls.heat.total_j,
                controls.contact_work_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "finite bounded thermal unknowns, iteration controls and independent energy tolerances required",
            ));
        }
        let check = || checkpoint(cx);
        let step = FrictionHeatStep::new(thermal, target_s, &self.route, check)?;
        let duration_s = step.duration_s();
        let trials = Cell::new(0usize);
        let evaluate =
            |x: &[f64]| {
                check()?;
                if trials.get() >= controls.maximum_trials {
                    return Err(PlasticSlidingHeatError::TrialBudget);
                }
                trials.set(trials.get() + 1);
                if x.len() != dimension || x.iter().any(|v| !v.is_finite()) {
                    return Err(PlasticSlidingHeatError::Invalid(
                        "complete finite energy/interface-temperature iterate required",
                    ));
                }
                let interface_temperature = finite(
                    interface_temperature_guess_k + x[cells] / controls.temperature_scale_j_per_k,
                )?;
                let prepare = |endpoint: &DepositionSimulation<'t>, interface_temperature|
                -> Result<(f64, PreparedMechanics), PlasticSlidingHeatError<E>> {
                let bulk_temperature = temperatures(endpoint, self.caloric, &self.route)?;
                let (source, mechanics, contact_areas) = self.advance_mechanics(
                    bulk_temperature, interface_temperature, motion, duration_s,
                    boundary, controls.friction, cx,
                )?;
                let contact_work_residual_j =
                    mechanics.admit_work(controls.contact_work_tolerance_j, cx)?;
                let removed_volume_fraction = self.admit_removal(mechanics.candidate())?;
                let [refreshed_contact_area_m2, endpoint_contact_area_m2] = contact_areas;
                let plastic_work_j = finite(
                    mechanics.temperature_normal().response.plastic_work_j
                        + mechanics.wear_normal().response.plastic_work_j,
                )?;
                let (normal_plastic_heat_j, nonthermal_plastic_work_j) =
                    self.plastic_work_heat.split(plastic_work_j)?;
                Ok((endpoint_contact_area_m2, PreparedMechanics {
                    source, mechanics, contact_work_residual_j, removed_volume_fraction,
                    refreshed_contact_area_m2, endpoint_contact_area_m2,
                    normal_plastic_heat_j, nonthermal_plastic_work_j,
                }))
            };
                let heat_input = |endpoint: &DepositionSimulation<'t>,
                                  interface_temperature,
                                  prepared: PreparedMechanics|
                 -> Result<
                    FrictionHeatInput<Staged<P>>,
                    PlasticSlidingHeatError<E>,
                > {
                    let PlasticSlidingHeatSource {
                        conductivity_w_m_k,
                        transport,
                    } = resolve(endpoint, interface_temperature)
                        .map_err(PlasticSlidingHeatError::Source)?;
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = transport;
                    for (side, cell) in self.route.cells().into_iter().enumerate() {
                        let heat_j = prepared.normal_plastic_heat_j[side];
                        if heat_j > 0.0 {
                            let power_w = finite(heat_j / duration_s)?;
                            if power_w == 0.0 {
                                return Err(DepositionError::Arithmetic.into());
                            }
                            let source_w = sources_w.entry(cell).or_insert(0.0);
                            *source_w = finite(*source_w + power_w)?;
                        }
                    }
                    Ok(FrictionHeatInput {
                        // Plastic heat enters its declared bulk receiver above;
                        // the interface receives friction once. Parameter work
                        // is never silently included in either heat source.
                        dissipated_j: prepared.mechanics.frictional_work_j(),
                        conductivity_w_m_k,
                        transport: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: Staged {
                                prepared,
                                payload,
                                conductivity_w_m_k,
                            },
                        },
                    })
                };
                if self.footprint == PlasticSlidingHeatFootprint::EvolvingActiveCells {
                    step.evaluate_with_area(
                        &x[..cells],
                        interface_temperature,
                        prepare,
                        heat_input,
                        check,
                    )
                    .map(|area_trial| {
                        let (trial, route) = area_trial.into_trial_and_route();
                        (trial, Some(route))
                    })
                } else {
                    step.evaluate(
                        &x[..cells],
                        interface_temperature,
                        |endpoint, temperature| {
                            let (_, prepared) = prepare(endpoint, temperature)?;
                            heat_input(endpoint, temperature, prepared)
                        },
                        check,
                    )
                    .map(|trial| (trial, None))
                }
            };
        let problem =
            FallibleFiniteDifference::new(dimension, controls.difference_step_j, |x: &[f64]| {
                evaluate(x).map(|(trial, _)| trial.residual_j().to_vec())
            })
            .map_err(numerical)?;
        let mut initial = energy_guess_j.to_vec();
        initial.push(0.0);
        let construction = NewtonKrylovState::new(&problem, initial, controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                cx.checkpoint().is_err() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        check()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(PlasticSlidingHeatError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(PlasticSlidingHeatError::NotConverged(result.report));
        }
        let (trial, successor_route) = evaluate(&solver.x)?;
        let (thermal, staged, heat) = trial.accept(controls.heat, check)?;
        let Staged {
            prepared,
            payload,
            conductivity_w_m_k,
        } = staged;
        let PreparedMechanics {
            source,
            mechanics,
            contact_work_residual_j,
            removed_volume_fraction,
            refreshed_contact_area_m2,
            endpoint_contact_area_m2,
            normal_plastic_heat_j,
            nonthermal_plastic_work_j,
        } = prepared;
        let t = mechanics.temperature_normal();
        let w = mechanics.wear_normal();
        let temperature_parameter_work_j = t.parameter_work_j;
        let wear_parameter_work_j = w.parameter_work_j;
        let normal_plastic_work_j = finite(t.response.plastic_work_j + w.response.plastic_work_j)?;
        let normal_work_j = finite(t.normal_work_j + w.normal_work_j)?;
        let removed_volume_m3 = mechanics.candidate().wear_state().volumes_m3();
        let frictional_work_j = mechanics.frictional_work_j();
        let other_heat_j =
            finite(heat.external_heat_j - normal_plastic_heat_j[0] - normal_plastic_heat_j[1])?;
        let mut active_cells = Vec::with_capacity(self.active_cells.len());
        for (i, p) in mechanics
            .candidate()
            .wear_state()
            .contact_state()
            .contact_state()
            .pressure_pa()
            .iter()
            .enumerate()
        {
            if i % 64 == 0 {
                check()?;
            }
            active_cells.push(*p > 0.0);
        }
        let contact = Self {
            source,
            // Clone the accepted candidate once for continuing the sealed
            // owner; all Newton trials retained only their native candidate.
            state: mechanics.candidate().clone(),
            material_cards: self.material_cards,
            caloric: self.caloric,
            route: successor_route.unwrap_or_else(|| self.route.clone()),
            original_body_ids: self.original_body_ids,
            original_volume_m3: self.original_volume_m3,
            active_cells,
            footprint: self.footprint,
            limits: self.limits,
            plastic_work_heat: self.plastic_work_heat,
            time_s: target_s,
        };
        check()?;
        Ok(PlasticSlidingHeatCandidate {
            contact,
            thermal,
            report: PlasticSlidingHeatReport {
                mechanics,
                heat,
                source: payload,
                conductivity_w_m_k,
                previous_contact_area_m2,
                refreshed_contact_area_m2,
                endpoint_contact_area_m2,
                nonlinear: result.report,
                trials: trials.get(),
                contact_work_residual_j,
                frictional_work_j,
                temperature_parameter_work_j,
                wear_parameter_work_j,
                normal_plastic_work_j,
                normal_plastic_heat_j,
                nonthermal_plastic_work_j,
                other_heat_j,
                normal_work_j,
                removed_volume_m3,
                removed_volume_fraction,
            },
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_mechanics<E>(
        &self,
        bulk_temperature: [f64; 2],
        interface_temperature: f64,
        motion: &[PlasticSlidingInput],
        duration_s: f64,
        boundary: PlasticSlidingHeatBoundary,
        controls: TangentialSolveConfig,
        cx: &Cx<'_>,
    ) -> Result<
        (
            SourcePlasticSlidingWear,
            PlasticSlidingHeatMechanics,
            [f64; 2],
        ),
        PlasticSlidingHeatError<E>,
    > {
        match boundary {
            PlasticSlidingHeatBoundary::HeldApproach => {
                let refresh = self.source.refresh_at_temperatures(
                    &self.state,
                    self.material_cards,
                    bulk_temperature,
                    interface_temperature,
                    cx,
                )?;
                self.complete_motion(
                    refresh,
                    cx,
                    |source, state| source.advance(state, motion, duration_s, controls, cx),
                    PlasticSlidingHeatMechanics::HeldApproach,
                )
            }
            PlasticSlidingHeatBoundary::MaintainedLoad {
                force_n,
                controls: normal,
            } => {
                let refresh = self.source.refresh_at_temperatures_at_load(
                    &self.state,
                    self.material_cards,
                    bulk_temperature,
                    interface_temperature,
                    force_n,
                    normal,
                    cx,
                )?;
                self.complete_motion(
                    refresh,
                    cx,
                    |source, state| {
                        source.advance_at_load(
                            state, motion, duration_s, controls, force_n, normal, cx,
                        )
                    },
                    PlasticSlidingHeatMechanics::MaintainedLoad,
                )
            }
            PlasticSlidingHeatBoundary::FixedDrive {
                drive,
                controls: normal,
            } => {
                let refresh = self.source.refresh_at_temperatures_with_drive(
                    &self.state,
                    self.material_cards,
                    bulk_temperature,
                    interface_temperature,
                    drive,
                    normal,
                    cx,
                )?;
                self.complete_motion(
                    refresh,
                    cx,
                    |source, state| {
                        source.advance_with_drive(
                            state, motion, duration_s, controls, drive, normal, cx,
                        )
                    },
                    PlasticSlidingHeatMechanics::FixedDrive,
                )
            }
        }
    }

    fn complete_motion<B, E>(
        &self,
        refresh: SourcePlasticSlidingThermalRefresh<B>,
        cx: &Cx<'_>,
        advance: impl FnOnce(
            &SourcePlasticSlidingWear,
            &SourcePlasticSlidingWearState,
        )
            -> Result<SourcePlasticSlidingWearAdvance<B>, SourcePlasticSlidingError>,
        wrap: impl FnOnce(PlasticSlidingHeatResponse<B>) -> PlasticSlidingHeatMechanics,
    ) -> Result<
        (
            SourcePlasticSlidingWear,
            PlasticSlidingHeatMechanics,
            [f64; 2],
        ),
        PlasticSlidingHeatError<E>,
    > {
        let (source, before, temperature) = refresh.into_parts();
        let refreshed_area = self.admit_footprint(&before, cx)?;
        if temperature.normal.parameter_work_j.abs() > self.limits.maximum_absolute_parameter_work_j
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "temperature parameter work exceeds the explicit external-work budget",
            ));
        }
        let sliding = advance(&source, &before)?;
        let endpoint_area = self.admit_footprint(&sliding.candidate, cx)?;
        self.admit_removal(&sliding.candidate)?;
        let supplied = finite(
            temperature.normal.parameter_work_j.abs() + sliding.normal.parameter_work_j.abs(),
        )?;
        if supplied > self.limits.maximum_absolute_parameter_work_j {
            return Err(PlasticSlidingHeatError::Invalid(
                "temperature and wear parameter work exceed their absolute external-work budget",
            ));
        }
        checkpoint(cx)?;
        Ok((
            source,
            wrap(PlasticSlidingHeatResponse {
                temperature,
                sliding,
            }),
            [refreshed_area, endpoint_area],
        ))
    }

    fn admit_thermal<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<[f64; 2], PlasticSlidingHeatError<E>> {
        self.route.validate(thermal)?;
        self.admit_footprint_support()?;
        if self.time_s != thermal.time_s()
            || self.route.cells().map(|id| thermal.cells()[&id].body_id) != self.original_body_ids
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "contact and thermal clocks or original body identities differ",
            ));
        }
        let temperature = temperatures(thermal, self.caloric, &self.route)?;
        let bulk = original_bulk(&self.source);
        if (0..2).any(|i| bulk[i].state_coordinate("T") != Some(temperature[i])) {
            return Err(PlasticSlidingHeatError::Invalid(
                "accepted bulk source temperatures differ from the thermal checkpoint",
            ));
        }
        self.admit_removal(&self.state)?;
        Ok(temperature)
    }

    fn admit_footprint<E>(
        &self,
        state: &SourcePlasticSlidingWearState,
        cx: &Cx<'_>,
    ) -> Result<f64, PlasticSlidingHeatError<E>> {
        let pressure = state
            .wear_state()
            .contact_state()
            .contact_state()
            .pressure_pa();
        if pressure.len() != self.active_cells.len() {
            return Err(PlasticSlidingHeatError::Invalid(
                "contact grid changed within fixed thermal registration",
            ));
        }
        let mut count = 0usize;
        for (i, (&p, &active)) in pressure.iter().zip(&self.active_cells).enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            if !p.is_finite() || p < 0.0 {
                return Err(PlasticSlidingHeatError::Invalid(
                    "native thermal footprint requires finite nonnegative pressure",
                ));
            }
            count += usize::from(p > 0.0);
            if self.footprint == PlasticSlidingHeatFootprint::FixedActiveCells
                && (p > 0.0) != active
            {
                return Err(PlasticSlidingHeatError::Invalid(
                    "temperature or wear changed the admitted active thermal footprint",
                ));
            }
        }
        if count == 0 {
            return Err(PlasticSlidingHeatError::Invalid(
                "detached normal contact has no massless thermal interface temperature",
            ));
        }
        let grid = self.source.wear_source().contact().model().config().grid;
        finite(count as f64 * (grid.cell_width_m * grid.cell_depth_m))
    }

    fn admit_footprint_support<E>(&self) -> Result<(), PlasticSlidingHeatError<E>> {
        if self.footprint == PlasticSlidingHeatFootprint::EvolvingActiveCells {
            let grid = self.source.wear_source().contact().model().config().grid;
            let support =
                finite(self.active_cells.len() as f64 * (grid.cell_width_m * grid.cell_depth_m))?;
            if support <= 0.0 {
                return Err(PlasticSlidingHeatError::Invalid(
                    "finite positive original grid area required",
                ));
            }
            for i in 0..2 {
                let face =
                    finite(self.route.geometry[i].volume_m3 / (2.0 * self.route.distances_m[i]))?;
                if support - face > 64.0 * f64::EPSILON * support.max(face) {
                    return Err(PlasticSlidingHeatError::Invalid(
                        "evolving contact grid exceeds an original uniform thermal face",
                    ));
                }
            }
        }
        Ok(())
    }

    fn admit_removal<E>(
        &self,
        state: &SourcePlasticSlidingWearState,
    ) -> Result<[f64; 2], PlasticSlidingHeatError<E>> {
        let volumes = state.wear_state().volumes_m3();
        let mut fractions = [0.0; 2];
        for i in 0..2 {
            fractions[i] = finite(volumes[i] / self.original_volume_m3[i])?;
            if volumes[i] < 0.0
                || (volumes[i] > 0.0 && fractions[i] == 0.0)
                || fractions[i] > self.limits.maximum_removed_volume_fraction[i]
            {
                return Err(PlasticSlidingHeatError::Invalid(
                    "cumulative removal exceeds original fixed-inventory thermal envelope",
                ));
            }
        }
        Ok(fractions)
    }
}

impl PlasticSlidingHeatMechanics {
    fn admit_work<E>(
        &self,
        tolerance: f64,
        cx: &Cx<'_>,
    ) -> Result<f64, PlasticSlidingHeatError<E>> {
        let mut maximum = 0.0_f64;
        for normal in [self.temperature_normal(), self.wear_normal()] {
            let residual = finite(
                normal.elastic_energy_change_j + normal.response.plastic_work_j
                    - normal.normal_work_j
                    - normal.parameter_work_j,
            )?;
            maximum = maximum.max(residual.abs());
        }
        let mut work = 0.0;
        let mut correction = 0.0;
        for (i, cell) in self.cells().iter().enumerate() {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            let response = &cell.response;
            let supplied = finite(
                response.traction_n[0] * response.slip_increment_m[0]
                    + response.traction_n[1] * response.slip_increment_m[1],
            )?;
            let residual = finite(response.stored_change_j + response.dissipated_j - supplied)?;
            maximum = maximum.max(residual.abs());
            let next = work + response.dissipated_j;
            correction += if work.abs() >= response.dissipated_j.abs() {
                (work - next) + response.dissipated_j
            } else {
                (response.dissipated_j - next) + work
            };
            work = next;
        }
        maximum = maximum.max(finite(work + correction - self.frictional_work_j())?.abs());
        if maximum > tolerance {
            return Err(PlasticSlidingHeatError::Invalid(
                "native contact work exceeds its independent energy tolerance",
            ));
        }
        Ok(maximum)
    }
}

fn original_bulk(source: &SourcePlasticSlidingWear) -> [&ResolvedMaterialStatePoint; 2] {
    let normal = source.wear_source().contact().normal_source().elastic();
    [
        normal.surface_a_state().resolved(),
        normal.surface_b_state().resolved(),
    ]
}

fn temperatures<E>(
    thermal: &DepositionSimulation<'_>,
    caloric: [&ResolvedHeatCapacityCurve; 2],
    route: &FrictionHeatRoute,
) -> Result<[f64; 2], PlasticSlidingHeatError<E>> {
    let mut result = [0.0; 2];
    for (i, id) in route.cells().into_iter().enumerate() {
        let cell = thermal
            .cells()
            .get(&id)
            .ok_or(PlasticSlidingHeatError::Invalid(
                "original receiving thermal cell is inactive",
            ))?;
        let curve = caloric[i].curve();
        if !thermal.uses_material_law(cell.material_id, curve) {
            return Err(PlasticSlidingHeatError::Invalid(
                "thermal timeline substituted an original caloric object",
            ));
        }
        let phase = curve.values_at_specific_enthalpy(cell.specific_h_j_kg)?;
        if phase.liquid_mass_fraction != 0.0 {
            return Err(PlasticSlidingHeatError::Invalid(
                "dry hardness-limited sliding refuses melted material",
            ));
        }
        result[i] = phase.temperature_k;
    }
    Ok(result)
}
