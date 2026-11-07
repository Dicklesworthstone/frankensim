//! Finite frictional heating with original temperature/depth hardness and wear.
//!
//! Each native caloric trial starts from the same accepted depth/shear history,
//! refreshes original E(T), nu(T), H(T,z), k(T) and mu(T,p), then removes material
//! using actual Coulomb slip and the original depth-hardness integral. The same
//! finite-cell and massless-interface residual used by ordinary plastic wear
//! solves both body temperatures and the interface temperature together.
//!
//! Temperature refresh, recession and exposure of a new hardness cap each retain
//! their own native normal-work receipt. All three contribute to the independent
//! work audit, absolute parameter-work budget and optional explicit plastic-heat
//! partition. Signed parameter work never supplies friction heat.
//!
//! This retains the native beginning-pressure wear split, uniform body thermal
//! receivers and either explicit fixed or original source-resolved contact
//! resistance. It does not resolve lateral thermal fields or certify a
//! reciprocal thermoelastic free-energy model.
//! Geometry and thermal mass stay fixed during an interval; cumulative removal
//! is constrained against the original receiving volumes.

/// Conservative export of newly worn mass and its current signed enthalpy.
pub mod removal;

/// Original resistance at actual final depth-contact pressure and temperature.
pub mod resistance;

/// Adaptive native heat/wear steps with per-body temperature and hardness control.
pub mod adaptive;

use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::temperature::DepthSlidingThermalRefresh;
use fs_contact::normal_patch::plastic_source::wear::sliding::hardness::{
    HardnessSlidingAdvance, HardnessSlidingContinuation, HardnessSlidingControl,
    SourceDepthSlidingWear,
};

use super::*;

/// The three original native normal continuations and actual tangential interval.
#[derive(Clone, Debug, PartialEq)]
pub struct DepthSlidingHeatMechanics {
    pub temperature: DepthSlidingThermalRefresh,
    pub sliding: HardnessSlidingAdvance<SourceDepthSlidingWear>,
}
impl DepthSlidingHeatMechanics {
    pub fn temperature_normal(&self) -> &FiniteGapPlasticContinuationReceipt {
        &self.temperature.normal.normal
    }
    /// Recession followed by exposure of the original new depth-hardness cap.
    pub fn wear_normals(&self) -> &[HardnessSlidingContinuation] {
        &self.sliding.normal
    }
    pub fn candidate(&self) -> &SourceDepthSlidingWear {
        &self.sliding.candidate
    }
    pub fn frictional_work_j(&self) -> f64 {
        self.sliding.sliding.frictional_work_j
    }
    pub fn cells(&self) -> &[PlasticCellSliding] {
        &self.sliding.sliding.cells
    }
}

/// Shared heat accounts with the complete original depth mechanics payload.
/// `wear_parameter_work_j` sums recession and cap exposure; their individual
/// signed accounts remain in `mechanics.wear_normals()`.
pub type DepthSlidingHeatReport<P> = PlasticSlidingHeatReport<P, DepthSlidingHeatMechanics>;

/// Accepted depth/shear/wear history bound to two original finite caloric cells.
#[derive(Clone)]
pub struct DepthSlidingHeatContact<'a> {
    source: SourceDepthSlidingWear,
    material_cards: [&'a MaterialCard; 2],
    caloric: [&'a ResolvedHeatCapacityCurve; 2],
    route: FrictionHeatRoute,
    original_body_ids: [u64; 2],
    original_volume_m3: [f64; 2],
    thermal_wear_checkpoint_m3: [f64; 2],
    active_cells: Vec<bool>,
    footprint: PlasticSlidingHeatFootprint,
    limits: PlasticSlidingHeatLimits,
    plastic_work_heat: PlasticWorkHeatPartition,
    resistance: Option<PlasticThermalResistance>,
    time_s: f64,
}

/// Accept the contact and thermal successors together after all original gates.
pub struct DepthSlidingHeatCandidate<'a, 't, P> {
    contact: DepthSlidingHeatContact<'a>,
    thermal: DepositionSimulation<'t>,
    report: DepthSlidingHeatReport<P>,
}
impl<'a, 't, P> DepthSlidingHeatCandidate<'a, 't, P> {
    pub const fn contact(&self) -> &DepthSlidingHeatContact<'a> {
        &self.contact
    }
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn report(&self) -> &DepthSlidingHeatReport<P> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepthSlidingHeatContact<'a>,
        DepositionSimulation<'t>,
        DepthSlidingHeatReport<P>,
    ) {
        (self.contact, self.thermal, self.report)
    }
}

impl<'a> DepthSlidingHeatContact<'a> {
    /// Bind an accepted original depth/sliding owner to exact caloric objects,
    /// ordered material cards and a fixed active thermal footprint.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        source: SourceDepthSlidingWear,
        material_cards: [&'a MaterialCard; 2],
        caloric: [&'a ResolvedHeatCapacityCurve; 2],
        route: &FrictionHeatRoute,
        limits: PlasticSlidingHeatLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, PlasticSlidingHeatError> {
        Self::new_with_footprint(
            thermal,
            source,
            material_cards,
            caloric,
            route,
            limits,
            PlasticSlidingHeatFootprint::FixedActiveCells,
            cx,
        )
    }

    /// Admit recruitment/opening on the original grid inside both uniform faces.
    /// Final native pressure area drives conductance and ambient coverage; no
    /// swept residence or lateral heat field is inferred. Detachment refuses.
    #[allow(clippy::too_many_arguments)]
    pub fn new_evolving_area(
        thermal: &DepositionSimulation<'_>,
        source: SourceDepthSlidingWear,
        material_cards: [&'a MaterialCard; 2],
        caloric: [&'a ResolvedHeatCapacityCurve; 2],
        route: &FrictionHeatRoute,
        limits: PlasticSlidingHeatLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, PlasticSlidingHeatError> {
        Self::new_with_footprint(
            thermal,
            source,
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
        source: SourceDepthSlidingWear,
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
        let normal = source.wear_state().normal_source();
        binding::admit_caloric_sources(normal, material_cards, caloric, cx)?;
        let temperature = temperatures(thermal, caloric, route)?;
        let bulk = [
            normal.elastic().surface_a_state(),
            normal.elastic().surface_b_state(),
        ];
        if (0..2).any(|i| bulk[i].resolved().state_coordinate("T") != Some(temperature[i])) {
            return Err(PlasticSlidingHeatError::Invalid(
                "accepted bulk source temperatures must equal the original thermal checkpoint",
            ));
        }
        let admitted = source.refresh_at_temperatures(
            material_cards,
            temperature,
            source.interface_temperature_k(),
            HardnessSlidingControl::Approach(source.wear_state().contact_state().approach_m()),
            cx,
        )?;
        if admitted.candidate != source {
            return Err(PlasticSlidingHeatError::Invalid(
                "thermal binding must not change accepted source or physical history",
            ));
        }
        let pressure = source.wear_state().contact_state().pressure_pa();
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
        let grid = source.wear_state().model().config().grid;
        let area = finite(count as f64 * (grid.cell_width_m * grid.cell_depth_m))?;
        let thermal_area = route.contact_area_m2();
        if count == 0 || (area - thermal_area).abs() > 64.0 * f64::EPSILON * area.max(thermal_area)
        {
            return Err(PlasticSlidingHeatError::Invalid(
                "actual loaded grid cells must exactly cover the original thermal face area",
            ));
        }
        let thermal_wear_checkpoint_m3 = source.wear_state().volumes_m3();
        let bound = Self {
            source,
            material_cards,
            caloric,
            route: route.clone(),
            original_body_ids: cells.map(|id| thermal.cells()[&id].body_id),
            original_volume_m3: cells.map(|id| thermal.cells()[&id].geometry.volume_m3),
            thermal_wear_checkpoint_m3,
            active_cells,
            footprint,
            limits,
            plastic_work_heat: PlasticWorkHeatPartition::default(),
            resistance: None,
            time_s: thermal.time_s(),
        };
        bound.admit_footprint_support()?;
        bound.admit_removal(&bound.source)?;
        checkpoint(cx)?;
        Ok(bound)
    }

    pub const fn source(&self) -> &SourceDepthSlidingWear {
        &self.source
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn route(&self) -> &FrictionHeatRoute {
        &self.route
    }
    pub const fn limits(&self) -> PlasticSlidingHeatLimits {
        self.limits
    }
    pub const fn footprint(&self) -> PlasticSlidingHeatFootprint {
        self.footprint
    }

    /// Declare conversion of future incremental native plastic work only.
    #[must_use]
    pub fn with_plastic_work_heat(mut self, partition: PlasticWorkHeatPartition) -> Self {
        self.plastic_work_heat = partition;
        self
    }
    pub const fn plastic_work_heat(&self) -> PlasticWorkHeatPartition {
        self.plastic_work_heat
    }

    /// Solve finite caloric feedback from the same accepted histories on every
    /// trial. Original full-depth support and actual loaded friction queries
    /// must accept. The callback supplies conductivity and OTHER heat once.
    /// Native shear, wear, normal and caloric candidates are accepted together.
    #[allow(clippy::too_many_arguments)]
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
    ) -> Result<DepthSlidingHeatCandidate<'a, 't, P>, PlasticSlidingHeatError<E>> {
        checkpoint(cx)?;
        self.admit_thermal(thermal)?;
        let previous_area = self.admit_footprint(&self.source, cx)?;
        let control = match boundary {
            PlasticSlidingHeatBoundary::HeldApproach => HardnessSlidingControl::Approach(
                self.source.wear_state().contact_state().approach_m(),
            ),
            PlasticSlidingHeatBoundary::MaintainedLoad { force_n, controls } => {
                HardnessSlidingControl::Load { force_n, controls }
            }
            PlasticSlidingHeatBoundary::FixedDrive { drive, controls } => {
                HardnessSlidingControl::Drive { drive, controls }
            }
        };
        let solved = engine::solve(
            thermal,
            target_s,
            &self.route,
            previous_area,
            self.footprint == PlasticSlidingHeatFootprint::EvolvingActiveCells
                || self.resistance.is_some(),
            controls,
            energy_guess_j,
            interface_temperature_guess_k,
            |endpoint, interface_temperature, duration_s| {
                let body_temperatures = temperatures(endpoint, self.caloric, &self.route)?;
                let temperature = self.source.refresh_at_temperatures(
                    self.material_cards,
                    body_temperatures,
                    interface_temperature,
                    control,
                    cx,
                )?;
                let refreshed_contact_area_m2 = self.admit_footprint(&temperature.candidate, cx)?;
                let mut supplied = temperature.normal.normal.parameter_work_j.abs();
                if supplied > self.limits.maximum_absolute_parameter_work_j {
                    return Err(PlasticSlidingHeatError::Invalid(
                        "temperature parameter work exceeds the explicit external-work budget",
                    ));
                }
                let sliding = temperature.candidate.advance(
                    motion,
                    duration_s,
                    controls.friction,
                    control,
                    cx,
                )?;
                if sliding.normal.len() != 2 {
                    return Err(PlasticSlidingHeatError::Invalid(
                        "depth wear requires original recession and exposed-hardness continuations",
                    ));
                }
                for normal in &sliding.normal {
                    supplied = finite(supplied + normal.normal.parameter_work_j.abs())?;
                }
                if supplied > self.limits.maximum_absolute_parameter_work_j {
                    return Err(PlasticSlidingHeatError::Invalid(
                        "temperature, recession and depth-hardness parameter work exceed their absolute external-work budget",
                    ));
                }
                let endpoint_contact_area_m2 = self.admit_footprint(&sliding.candidate, cx)?;
                let removed_volume_fraction = self.admit_removal(&sliding.candidate)?;
                let removed_volume_m3 = sliding.candidate.wear_state().volumes_m3();
                let mechanics = DepthSlidingHeatMechanics {
                    temperature,
                    sliding,
                };
                let work = engine::mechanical_work(
                    mechanics.temperature_normal(),
                    &[
                        &mechanics.wear_normals()[0].normal,
                        &mechanics.wear_normals()[1].normal,
                    ],
                    mechanics.cells(),
                    mechanics.frictional_work_j(),
                    controls.contact_work_tolerance_j,
                    cx,
                )?;
                let (normal_plastic_heat_j, nonthermal_plastic_work_j) =
                    self.plastic_work_heat.split(work.normal_plastic_work_j)?;
                Ok(engine::PreparedMechanics {
                    source: (),
                    mechanics,
                    work,
                    removed_volume_m3,
                    removed_volume_fraction,
                    refreshed_contact_area_m2,
                    endpoint_contact_area_m2,
                    normal_plastic_heat_j,
                    nonthermal_plastic_work_j,
                })
            },
            |prepared, interface_temperature| {
                self.resistance
                    .as_ref()
                    .map(|resistance| {
                        resistance.sample_depth(
                            prepared.mechanics.candidate(),
                            interface_temperature,
                            cx,
                        )
                    })
                    .transpose()
            },
            resolve,
            cx,
        )?;
        let engine::Solved {
            thermal,
            successor_route,
            report,
            ..
        } = solved;
        let source = report.mechanics.candidate().clone();
        let mut active_cells = Vec::with_capacity(self.active_cells.len());
        for (i, p) in source
            .wear_state()
            .contact_state()
            .pressure_pa()
            .iter()
            .enumerate()
        {
            if i % 64 == 0 {
                checkpoint(cx)?;
            }
            active_cells.push(*p > 0.0);
        }
        let contact = Self {
            source,
            material_cards: self.material_cards,
            caloric: self.caloric,
            route: successor_route.unwrap_or_else(|| self.route.clone()),
            original_body_ids: self.original_body_ids,
            original_volume_m3: self.original_volume_m3,
            thermal_wear_checkpoint_m3: self.thermal_wear_checkpoint_m3,
            active_cells,
            footprint: self.footprint,
            limits: self.limits,
            plastic_work_heat: self.plastic_work_heat,
            resistance: self.resistance.clone(),
            time_s: target_s,
        };
        checkpoint(cx)?;
        Ok(DepthSlidingHeatCandidate {
            contact,
            thermal,
            report,
        })
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
        let normal = self.source.wear_state().normal_source().elastic();
        let bulk = [
            normal.surface_a_state().resolved(),
            normal.surface_b_state().resolved(),
        ];
        if (0..2).any(|i| bulk[i].state_coordinate("T") != Some(temperature[i])) {
            return Err(PlasticSlidingHeatError::Invalid(
                "accepted bulk source temperatures differ from the thermal checkpoint",
            ));
        }
        self.admit_removal(&self.source)?;
        Ok(temperature)
    }

    fn admit_footprint<E>(
        &self,
        source: &SourceDepthSlidingWear,
        cx: &Cx<'_>,
    ) -> Result<f64, PlasticSlidingHeatError<E>> {
        let grid = source.wear_state().model().config().grid;
        native_footprint(
            source.wear_state().contact_state().pressure_pa(),
            &self.active_cells,
            self.footprint,
            grid.cell_width_m * grid.cell_depth_m,
            cx,
        )
    }

    fn admit_footprint_support<E>(&self) -> Result<(), PlasticSlidingHeatError<E>> {
        let grid = self.source.wear_state().model().config().grid;
        footprint_support(
            &self.route,
            self.footprint,
            self.active_cells.len(),
            grid.cell_width_m * grid.cell_depth_m,
        )
    }

    fn admit_removal<E>(
        &self,
        source: &SourceDepthSlidingWear,
    ) -> Result<[f64; 2], PlasticSlidingHeatError<E>> {
        removal_fractions(
            source.wear_state().volumes_m3(),
            self.original_volume_m3,
            self.limits,
        )
    }
}
