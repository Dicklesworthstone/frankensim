//! Permanent creep and spatial heat on the original tetrahedral/enthalpy owners.
//!
//! Each reference tetrahedron binds one actual caloric cell and one original
//! Norton material source. The explicitly selected zero-expansion reduction
//! holds E, nu, density and reference geometry fixed; local temperature changes
//! the source-pinned Arrhenius kinetics. Conduction, reservoir exchange, thermal
//! events and observer histories remain in the existing thermal owner.
//!
//! Every mechanical and thermal proposal starts from the accepted state. The
//! final mechanics are solved at the actual returned cell temperatures. Only
//! physical creep heat enters caloric storage; backward-Euler loss is separate.
//! All local heat, temperature and complete work gates precede publication.
//! This is small-strain quasistatic secondary creep, not finite deformation,
//! phase change, variable elastic storage, thermal expansion or a life prediction.

use std::collections::{BTreeMap, BTreeSet};
use std::convert::Infallible;
use std::sync::Arc;

use fs_exec::Cx;
use fs_matdb::MaterialCard;
use fs_material::SmallStrainLaw;
use fs_material::creep::source::{CreepSourceError, NortonCreepSource, ResolvedNortonCreep};
use fs_material::creep::thermal::{
    FeedbackError, MAX_FEEDBACK_ITERATIONS, admit_zero_expansion_caloric,
};
use fs_material::creep::{CreepEnergy, CreepState, CreepTolerances};
use fs_material::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::{LoadError3, TetLoadGeometry3};
use fs_solid::plastic3::creep::{CreepTetEnergy, CreepTetError, CreepTetStep, CreepTetTolerances};
use fs_solid::plastic3::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetProblem, SmallStrainTetSolution,
};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    Event, HeatTransfer, StepConfig, ThermalLaw,
};

#[derive(Debug)]
pub enum CreepCouplingError<E = Infallible> {
    Invalid(&'static str),
    Source(CreepSourceError),
    Material(FeedbackError),
    Caloric(PhaseStateError),
    Geometry(LoadError3),
    Mechanical(CreepTetError),
    Thermal(DepositionError),
    Observer(E),
    Cancelled,
    NotConverged {
        iterations: usize,
        temperature_defect_k: f64,
        heat_defect_j: f64,
        thermal_residual_j: f64,
        combined_residual_j: f64,
    },
}
impl<E> From<DepositionError> for CreepCouplingError<E> {
    fn from(error: DepositionError) -> Self {
        if matches!(error, DepositionError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}
impl<E: std::fmt::Debug> std::fmt::Display for CreepCouplingError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "spatial creep/enthalpy coupling refused: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CreepCouplingError<E> {}

/// Source and thermal cell in tetrahedron connectivity order. The caloric curve
/// must be the actual object registered with the existing thermal simulation.
pub struct CreepCell<'a> {
    pub cell: CellId,
    pub card: &'a MaterialCard,
    pub source: NortonCreepSource<'a>,
    pub caloric: &'a ResolvedHeatCapacityCurve,
}

#[derive(Clone, Copy, Debug)]
pub struct CreepConfig {
    pub maximum_iterations: usize,
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    pub stress_tolerance_pa: f64,
    /// Relaxation of trial cell enthalpy in (0,1].
    pub enthalpy_relaxation: f64,
}

#[derive(Clone, Debug)]
struct BoundCell<'a> {
    source: NortonCreepSource<'a>,
    caloric: &'a ResolvedHeatCapacityCurve,
    accepted: DepositedCell,
}

/// Accepted native creep history with an explicit zero-expansion thermal model.
#[derive(Clone, Debug)]
pub struct ZeroExpansionCreepSolid<'a> {
    nodes_m: &'a [[f64; 3]],
    tetrahedra: &'a [[usize; 4]],
    geometry: Arc<TetLoadGeometry3<'a>>,
    cells: Vec<BoundCell<'a>>,
    materials: Vec<ResolvedNortonCreep>,
    solution: SmallStrainTetSolution<CreepState>,
    budget: TetAssemblyBudget,
    volume_m3: f64,
    time_s: f64,
}

impl<'a> ZeroExpansionCreepSolid<'a> {
    /// Bind a virgin, stress-free solid to the actual active thermal inventory.
    /// Original source admission is shared with the uniform creep/heat owner.
    pub fn try_new(
        thermal: &DepositionSimulation<'_>,
        nodes_m: &'a [[f64; 3]],
        tetrahedra: &'a [[usize; 4]],
        sources: Vec<CreepCell<'a>>,
        budget: TetAssemblyBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, CreepCouplingError> {
        poll(cx)?;
        if sources.len() != tetrahedra.len() || sources.is_empty() {
            return Err(CreepCouplingError::Invalid(
                "one source cell per nonempty tetrahedral field required",
            ));
        }
        let face_entries = tetrahedra
            .len()
            .checked_mul(4)
            .ok_or(CreepCouplingError::Invalid("face count overflow"))?;
        let geometry = Arc::new(
            TetLoadGeometry3::new(nodes_m, tetrahedra, budget, face_entries, cx)
                .map_err(CreepCouplingError::Geometry)?,
        );
        let mut used = BTreeSet::new();
        let mut cells = Vec::with_capacity(sources.len());
        let mut materials = Vec::with_capacity(sources.len());
        let mut volume_m3 = 0.0;
        for (element, entry) in sources.into_iter().enumerate() {
            poll(cx)?;
            let cell = thermal
                .cells()
                .get(&entry.cell)
                .ok_or(CreepCouplingError::Invalid("inactive structural cell"))?;
            if !used.insert(entry.cell)
                || !thermal
                    .uses_material_law(cell.material_id, entry.caloric.curve() as &dyn ThermalLaw)
            {
                return Err(CreepCouplingError::Invalid(
                    "unique cell mapping and original caloric object required",
                ));
            }
            let density =
                admit_zero_expansion_caloric(entry.card, &entry.source, entry.caloric, &mut || {
                    cx.checkpoint().is_err()
                })
                .map_err(CreepCouplingError::Material)?;
            let volume = geometry.element_volumes_m3()[element];
            if !super::same_measure(cell.geometry.volume_m3, volume)
                || !super::same_measure(cell.mass_kg, density * volume)
            {
                return Err(CreepCouplingError::Invalid(
                    "cell volume and mass must match native geometry and source density",
                ));
            }
            let state = entry
                .caloric
                .curve()
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(CreepCouplingError::Caloric)?;
            materials.push(
                entry
                    .source
                    .at_temperature(state.temperature_k(), 1.0)
                    .map_err(CreepCouplingError::Source)?,
            );
            volume_m3 = finite(volume_m3 + volume)?;
            cells.push(BoundCell {
                source: entry.source,
                caloric: entry.caloric,
                accepted: cell.clone(),
            });
        }
        if !positive(volume_m3) {
            return Err(CreepCouplingError::Invalid(
                "positive finite structural volume required",
            ));
        }
        let solution = SmallStrainTetSolution {
            displacement_m: vec![[0.0; 3]; nodes_m.len()],
            element_strain: vec![[0.0; 6]; tetrahedra.len()],
            element_stress: vec![[0.0; 6]; tetrahedra.len()],
            element_states: materials
                .iter()
                .map(SmallStrainLaw::initial_state)
                .collect(),
            increments: Vec::new(),
        };
        poll(cx)?;
        Ok(Self {
            nodes_m,
            tetrahedra,
            geometry,
            cells,
            materials,
            solution,
            budget,
            volume_m3,
            time_s: thermal.time_s(),
        })
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn geometry(&self) -> &TetLoadGeometry3<'a> {
        &self.geometry
    }
    pub fn solution(&self) -> &SmallStrainTetSolution<CreepState> {
        &self.solution
    }
    pub fn materials(&self) -> &[ResolvedNortonCreep] {
        &self.materials
    }

    fn check_thermal<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<(), CreepCouplingError<E>> {
        for bound in &self.cells {
            poll(cx)?;
            let old = &bound.accepted;
            let cell = thermal
                .cells()
                .get(&old.id)
                .ok_or(CreepCouplingError::Invalid(
                    "structural cell became inactive",
                ))?;
            if cell.id != old.id
                || cell.segment_id != old.segment_id
                || cell.body_id != old.body_id
                || cell.material_id != old.material_id
                || cell.geometry != old.geometry
                || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
                || accepted && cell.specific_h_j_kg.to_bits() != old.specific_h_j_kg.to_bits()
                || !thermal
                    .uses_material_law(cell.material_id, bound.caloric.curve() as &dyn ThermalLaw)
            {
                return Err(CreepCouplingError::Invalid(
                    "structural identity, geometry, mass or accepted enthalpy changed",
                ));
            }
            bound
                .caloric
                .curve()
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(CreepCouplingError::Caloric)?;
        }
        Ok(())
    }

    fn enthalpies<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, CreepCouplingError<E>> {
        self.check_thermal(thermal, accepted, cx)?;
        Ok(self
            .cells
            .iter()
            .map(|cell| thermal.cells()[&cell.accepted.id].specific_h_j_kg)
            .collect())
    }

    #[allow(clippy::too_many_arguments)]
    fn mechanics<E>(
        &self,
        enthalpies: &[f64],
        duration_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        config: CreepConfig,
        cx: &Cx<'_>,
    ) -> Result<MechanicalProposal, CreepCouplingError<E>> {
        let mut materials = Vec::with_capacity(self.cells.len());
        let mut temperatures_k = Vec::with_capacity(self.cells.len());
        for (cell, &h) in self.cells.iter().zip(enthalpies) {
            poll(cx)?;
            let temperature = cell
                .caloric
                .curve()
                .state_at_specific_enthalpy(h)
                .map_err(CreepCouplingError::Caloric)?
                .temperature_k();
            materials.push(
                cell.source
                    .at_temperature(temperature, duration_s)
                    .map_err(CreepCouplingError::Source)?,
            );
            temperatures_k.push(temperature);
        }
        let problem = SmallStrainTetProblem {
            nodes_m: self.nodes_m,
            tetrahedra: self.tetrahedra,
            laws: SmallStrainLawField::PerElement(&materials),
            prescribed_m,
            nodal_forces_n,
            load_path: &[1.0],
            budget: self.budget,
            settings,
        };
        let step = problem
            .solve_creep_from(
                self.solution.displacement_m.as_flattened(),
                &self.solution.element_states,
                CreepTetTolerances {
                    material: CreepTolerances {
                        stress_pa: config.stress_tolerance_pa,
                        energy_j_m3: config.energy_tolerance_j / self.volume_m3,
                    },
                    energy_j: config.energy_tolerance_j,
                },
                cx,
            )
            .map_err(CreepCouplingError::Mechanical)?;
        let heat_j = step
            .element_energy
            .iter()
            .zip(&step.element_volume_m3)
            .map(|(e, v)| finite(e.physical_dissipation_j_m3 * v))
            .collect::<Result<_, _>>()?;
        Ok(MechanicalProposal {
            materials,
            temperatures_k,
            step,
            heat_j,
        })
    }

    /// One atomic loaded-solid/thermal proposal. Loads and supports are endpoint
    /// values; native load geometry can prepare pressure, traction or body forces.
    /// Refine physical intervals and thermal substeps separately. Feedback
    /// convergence bounds the discrete account, not continuum/model error.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        thermal_config: StepConfig,
        config: CreepConfig,
        cx: &Cx<'_>,
    ) -> Result<CreepCandidate<'a, 't, O>, CreepCouplingError<O::Error>> {
        poll(cx)?;
        let duration_s = target_s - self.time_s;
        if thermal.time_s() != self.time_s
            || !target_s.is_finite()
            || !positive(duration_s)
            || config.maximum_iterations == 0
            || config.maximum_iterations > MAX_FEEDBACK_ITERATIONS
            || !positive(config.temperature_tolerance_k)
            || !positive(config.energy_tolerance_j)
            || !positive(config.stress_tolerance_pa)
            || !positive(config.enthalpy_relaxation)
            || config.enthalpy_relaxation > 1.0
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(CreepCouplingError::Invalid(
                "matching clocks, positive physical interval and bounded tolerances required",
            ));
        }
        let mut enthalpies = self.enthalpies(thermal, true, cx)?;
        let mut defects = [f64::INFINITY; 4];
        for iteration in 1..=config.maximum_iterations {
            poll(cx)?;
            let trial = self.mechanics(
                &enthalpies,
                duration_s,
                prescribed_m,
                nodal_forces_n,
                settings,
                config,
                cx,
            )?;
            let heat = CreepHeat {
                solid: self,
                inner: observer.clone(),
                cx,
                watts: self
                    .cells
                    .iter()
                    .zip(&trial.heat_j)
                    .map(|(cell, q)| Ok((cell.accepted.id, finite(*q / duration_s)?)))
                    .collect::<Result<_, CreepCouplingError<O::Error>>>()?,
            };
            let (next_thermal, heat, thermal_audit) =
                thermal.advance_observed(&heat, target_s, thermal_config, || {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                })?;
            let final_h = self.enthalpies(&next_thermal, false, cx)?;
            let mut temperature_defect: f64 = 0.0;
            for ((cell, &h), &trial_temperature) in
                self.cells.iter().zip(&final_h).zip(&trial.temperatures_k)
            {
                poll(cx)?;
                let actual_temperature = cell
                    .caloric
                    .curve()
                    .state_at_specific_enthalpy(h)
                    .map_err(CreepCouplingError::Caloric)?
                    .temperature_k();
                temperature_defect =
                    temperature_defect.max(finite((actual_temperature - trial_temperature).abs())?);
            }
            // A failed temperature gate cannot publish. Avoid a second global
            // Newton solve until a pair can be accepted, or the final iteration
            // needs complete refusal diagnostics. The thermal fixed-point map
            // and all final-temperature stress/heat/work gates stay unchanged.
            if temperature_defect > config.temperature_tolerance_k
                && iteration < config.maximum_iterations
            {
                for (h, end_h) in enthalpies.iter_mut().zip(final_h) {
                    poll(cx)?;
                    *h = finite(
                        (1.0 - config.enthalpy_relaxation) * *h
                            + config.enthalpy_relaxation * end_h,
                    )?;
                }
                continue;
            }
            let actual = self.mechanics(
                &final_h,
                duration_s,
                prescribed_m,
                nodal_forces_n,
                settings,
                config,
                cx,
            )?;
            let applied_heat_j = finite(trial.heat_j.iter().sum())?;
            let physical_heat_j = actual.step.energy.physical_dissipation_j;
            let external_heat_j = finite(thermal_audit.external_heat_j - applied_heat_j)?;
            let thermal_residual_j = finite(
                thermal_audit.stored_energy_change_j
                    - thermal_audit.incoming_energy_j
                    - external_heat_j
                    - physical_heat_j,
            )?;
            let energy = actual.step.energy;
            let combined_residual_j = finite(
                energy.final_stored_j - energy.initial_stored_j
                    + thermal_audit.stored_energy_change_j
                    + energy.numerical_dissipation_j
                    - energy.applied_work_j
                    - energy.support_work_j
                    - external_heat_j
                    - thermal_audit.incoming_energy_j,
            )?;
            defects = [
                temperature_defect,
                finite(
                    trial
                        .heat_j
                        .iter()
                        .zip(&actual.heat_j)
                        .map(|(a, b)| (a - b).abs())
                        .sum(),
                )?,
                thermal_residual_j,
                combined_residual_j,
            ];
            if defects[0] <= config.temperature_tolerance_k
                && defects[1] <= config.energy_tolerance_j
                && defects[2].abs() <= config.energy_tolerance_j
                && defects[3].abs() <= config.energy_tolerance_j
            {
                let mut solid = self.clone();
                solid.time_s = target_s;
                solid.materials = actual.materials;
                solid.solution = actual.step.solution;
                for bound in &mut solid.cells {
                    bound.accepted = next_thermal.cells()[&bound.accepted.id].clone();
                }
                let report = CreepReport {
                    iterations: iteration,
                    thermal_audit,
                    mechanical_energy: energy,
                    element_energy: actual.step.element_energy,
                    applied_heat_j,
                    external_heat_j,
                    temperature_defect_k: defects[0],
                    heat_defect_j: defects[1],
                    thermal_residual_j,
                    combined_residual_j,
                };
                poll(cx)?;
                return Ok(CreepCandidate {
                    thermal: next_thermal,
                    solid,
                    observer: heat.inner,
                    report,
                });
            }
            for (h, end_h) in enthalpies.iter_mut().zip(final_h) {
                *h = finite(
                    (1.0 - config.enthalpy_relaxation) * *h + config.enthalpy_relaxation * end_h,
                )?;
            }
        }
        poll(cx)?;
        Err(CreepCouplingError::NotConverged {
            iterations: config.maximum_iterations,
            temperature_defect_k: defects[0],
            heat_defect_j: defects[1],
            thermal_residual_j: defects[2],
            combined_residual_j: defects[3],
        })
    }
}

struct MechanicalProposal {
    materials: Vec<ResolvedNortonCreep>,
    temperatures_k: Vec<f64>,
    step: CreepTetStep,
    heat_j: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct CreepReport {
    pub iterations: usize,
    /// Thermal input includes physical creep heat; external_heat_j removes it.
    pub thermal_audit: EnergyAudit,
    pub mechanical_energy: CreepTetEnergy,
    /// Native material energy densities [J/m³], in connectivity order.
    pub element_energy: Vec<CreepEnergy>,
    pub applied_heat_j: f64,
    pub external_heat_j: f64,
    pub temperature_defect_k: f64,
    /// Sum of absolute per-cell generated-heat defects; cancellation cannot hide them.
    pub heat_defect_j: f64,
    pub thermal_residual_j: f64,
    pub combined_residual_j: f64,
}

/// Accept the actual solid, thermal history and observer state together.
pub struct CreepCandidate<'m, 't, O> {
    thermal: DepositionSimulation<'t>,
    solid: ZeroExpansionCreepSolid<'m>,
    observer: O,
    report: CreepReport,
}
impl<'m, 't, O> CreepCandidate<'m, 't, O> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn solid(&self) -> &ZeroExpansionCreepSolid<'m> {
        &self.solid
    }
    pub fn observer(&self) -> &O {
        &self.observer
    }
    pub fn report(&self) -> &CreepReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ZeroExpansionCreepSolid<'m>,
        O,
        CreepReport,
    ) {
        (self.thermal, self.solid, self.observer, self.report)
    }
}

#[derive(Clone)]
struct CreepHeat<'s, 'm, 'cx, O> {
    solid: &'s ZeroExpansionCreepSolid<'m>,
    inner: O,
    cx: &'s Cx<'cx>,
    watts: BTreeMap<CellId, f64>,
}
impl<O: DepositionObserver> DepositionObserver for CreepHeat<'_, '_, '_, O> {
    type Error = CreepCouplingError<O::Error>;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner.begin().map_err(CreepCouplingError::Observer)
    }
    fn requires_step_interval(&self) -> bool {
        self.inner.requires_step_interval()
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .step_sources(thermal, start, dt, end, check)
            .map_err(CreepCouplingError::Observer)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.solid.check_thermal(thermal, false, self.cx)?;
        let mut rates = self
            .inner
            .sources(thermal, check)
            .map_err(CreepCouplingError::Observer)?;
        for (&cell, &watts) in &self.watts {
            check()?;
            let rate = rates.entry(cell).or_insert(0.0);
            *rate = finite(*rate + watts)?;
        }
        Ok(rates)
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner
            .conductances(thermal, check)
            .map_err(CreepCouplingError::Observer)
    }
    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .reservoir_conductances(thermal, check)
            .map_err(CreepCouplingError::Observer)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner
            .heat_transfers(thermal, check)
            .map_err(CreepCouplingError::Observer)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.inner
            .limit_timestep_s(thermal, proposed, check)
            .map_err(CreepCouplingError::Observer)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.inner
            .event(event, thermal)
            .map_err(CreepCouplingError::Observer)?;
        self.solid.check_thermal(thermal, false, self.cx)
    }
    fn step(
        &mut self,
        start: f64,
        dt: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner
            .step(start, dt, old, thermal, check)
            .map_err(CreepCouplingError::Observer)?;
        self.solid.check_thermal(thermal, false, self.cx)
    }
}

fn poll<E>(cx: &Cx<'_>) -> Result<(), CreepCouplingError<E>> {
    cx.checkpoint().map_err(|_| CreepCouplingError::Cancelled)
}
fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn finite<E>(x: f64) -> Result<f64, CreepCouplingError<E>> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(CreepCouplingError::Invalid(
            "unrepresentable creep work or heat arithmetic",
        ))
    }
}
