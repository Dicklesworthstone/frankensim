//! Source-bound thermal plasticity on the existing frame and caloric owners.
//!
//! The original force-based fibers retain permanent strain, Prager memory,
//! every cyclic backstress branch and temperature. Each finite caloric cell prepares one actual-temperature law
//! per coupling trial; all fibers in that cell share its immutable interval.
//! The native frame solver and `DepositionSimulation` still own mechanics and
//! heat transport. Rejected trials advance neither accepted history nor clocks.
//!
//! Physical material dissipation and signed reversible exchange enter each
//! original caloric inventory once. Hardening storage, mechanical entropy and
//! material/inertial numerical loss remain separate. The endpoint coupling is
//! first order: refine physical intervals and thermal substeps independently.
//! The native small-strain, linear-geometry and solid-source domains remain in
//! force; cells represent the caller's explicit thermal spatial approximation.

mod observer;
mod material;

/// Reciprocal finite heating and native source-bound laminate equilibrium.
pub mod laminate;

pub use material::FrameThermalSource;

use fs_exec::Cx;
use fs_exec::reduce::Compensated;
use fs_material::phase::{PhaseStateError, ResolvedHeatCapacityCurve};
use fs_material::plastic::thermal::{FeedbackError, MAX_FEEDBACK_ITERATIONS};
use fs_solid::frame2d::Frame2d;
use fs_solid::frame2d::creep::{
    FrameCreepDynamicEnergy, FrameCreepEnergy, FrameCreepError, FrameCreepTolerances,
    FrameFiberCreepEnergy,
};
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit,
    StepConfig, ThermalLaw,
};
use std::collections::{BTreeMap, BTreeSet};
use std::convert::Infallible;

/// One original thermal inventory and the material family that owns its data.
/// Many original frame quadrature fibers may share this cell. The supplied
/// caloric object must be the actual one registered with the thermal owner.
#[derive(Clone, Debug)]
pub struct FrameThermalCell<'a> {
    pub cell: CellId,
    pub source: FrameThermalSource<'a>,
}

/// Explicit correspondence from an original frame fiber to a caloric cell.
/// Entry order is immaterial; the stored order is the native frame's order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameFiberCell {
    pub element: usize,
    pub section: usize,
    pub fiber: usize,
    pub cell: CellId,
}

/// Use the original force-equilibrium or inertial frame owner for this interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameThermalMechanics {
    Quasistatic,
    BackwardEuler,
}

/// Bounded coupled-temperature and complete-energy gates.
#[derive(Clone, Copy, Debug)]
pub struct FrameEnthalpyConfig {
    pub maximum_iterations: usize,
    pub temperature_tolerance_k: f64,
    pub energy_tolerance_j: f64,
    /// Relaxation of specific caloric enthalpy in (0,1].
    pub enthalpy_relaxation: f64,
}

#[derive(Debug)]
pub enum FrameEnthalpyError<E = Infallible> {
    Invalid(&'static str),
    Material(FeedbackError),
    Caloric(PhaseStateError),
    Mechanical(FrameCreepError),
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

impl<E> From<DepositionError> for FrameEnthalpyError<E> {
    fn from(error: DepositionError) -> Self {
        if matches!(error, DepositionError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}

impl<E: std::fmt::Debug> std::fmt::Display for FrameEnthalpyError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "frame/enthalpy coupling refused: {self:?}")
    }
}

impl<E: std::error::Error + 'static> std::error::Error for FrameEnthalpyError<E> {}

#[derive(Clone, Debug)]
struct BoundCell<'a> {
    source: FrameThermalSource<'a>,
    caloric: &'a ResolvedHeatCapacityCurve,
    accepted: DepositedCell,
}

#[derive(Clone, Copy, Debug)]
struct BoundFiber {
    location: (usize, usize, usize),
    cell: usize,
    volume_m3: f64,
}

/// Accepted nonlinear frame and the original finite thermal correspondence.
/// The frame is read-only outside this owner; `advance` returns a private joint
/// candidate that can be accepted by taking all its parts together.
#[derive(Clone, Debug)]
pub struct ThermoplasticFrame<'a> {
    frame: Frame2d,
    cells: Vec<BoundCell<'a>>,
    fibers: Vec<BoundFiber>,
    maximum_fiber_points: usize,
    time_s: f64,
}

impl<'a> ThermoplasticFrame<'a> {
    /// Bind a complete original thermal-fiber frame, including nonvirgin
    /// histories, to matching actual cell temperatures and reference masses.
    /// Every structural fiber and every supplied cell must be used exactly once
    /// in the correspondence; sharing a caloric cell aggregates its volume.
    pub fn try_new(
        frame: Frame2d,
        thermal: &DepositionSimulation<'_>,
        sources: Vec<FrameThermalCell<'a>>,
        mapping: Vec<FrameFiberCell>,
        maximum_fiber_points: usize,
        cx: &Cx<'_>,
    ) -> Result<Self, FrameEnthalpyError> {
        poll(cx)?;
        if sources.is_empty()
            || maximum_fiber_points == 0
            || sources.len() > maximum_fiber_points
            || mapping.len() > maximum_fiber_points
        {
            return Err(FrameEnthalpyError::Invalid(
                "nonempty bounded thermal sources and complete fiber mapping required",
            ));
        }
        let points = frame
            .thermal_inelastic_fibers(maximum_fiber_points, || cx.checkpoint().is_err())
            .map_err(FrameEnthalpyError::Mechanical)?;
        if mapping.len() != points.len() {
            return Err(FrameEnthalpyError::Invalid(
                "one thermal correspondence per original frame fiber required",
            ));
        }
        let mut positions = BTreeMap::new();
        for item in mapping {
            poll(cx)?;
            if positions
                .insert((item.element, item.section, item.fiber), item.cell)
                .is_some()
            {
                return Err(FrameEnthalpyError::Invalid(
                    "duplicate fiber correspondence",
                ));
            }
        }
        let mut cells = Vec::with_capacity(sources.len());
        let mut cell_index = BTreeMap::new();
        let mut prepared = Vec::with_capacity(sources.len());
        let mut temperatures = Vec::with_capacity(sources.len());
        for item in sources {
            poll(cx)?;
            if cell_index.insert(item.cell, cells.len()).is_some() {
                return Err(FrameEnthalpyError::Invalid("duplicate thermal source cell"));
            }
            let accepted = thermal
                .cells()
                .get(&item.cell)
                .ok_or(FrameEnthalpyError::Invalid(
                    "inactive structural thermal cell",
                ))?;
            let caloric = item.source.base().caloric_source();
            if !thermal.uses_material_law(accepted.material_id, caloric.curve() as &dyn ThermalLaw)
            {
                return Err(FrameEnthalpyError::Invalid(
                    "the original registered caloric object is required",
                ));
            }
            let temperature = caloric
                .curve()
                .state_at_specific_enthalpy(accepted.specific_h_j_kg)
                .map_err(FrameEnthalpyError::Caloric)?
                .temperature_k();
            prepared.push(
                item.source
                    .prepare(temperature, temperature, 1.0)
                    .map_err(FrameEnthalpyError::Material)?,
            );
            temperatures.push(temperature);
            cells.push(BoundCell {
                source: item.source,
                caloric,
                accepted: accepted.clone(),
            });
        }
        let mut volumes = vec![Sum::default(); cells.len()];
        let mut masses = vec![Sum::default(); cells.len()];
        let mut member_masses = BTreeMap::<usize, (f64, Sum)>::new();
        let mut used = BTreeSet::new();
        let mut fibers = Vec::with_capacity(points.len());
        for point in points {
            poll(cx)?;
            let location = (point.element, point.section, point.fiber);
            let id = positions
                .remove(&location)
                .ok_or(FrameEnthalpyError::Invalid(
                    "thermal mapping does not name the original frame topology",
                ))?;
            let index = *cell_index.get(&id).ok_or(FrameEnthalpyError::Invalid(
                "fiber correspondence names an unbound caloric cell",
            ))?;
            let law = &prepared[index];
            point
                .law
                .admit_successor(law)
                .map_err(FrameEnthalpyError::Mechanical)?;
            let temperature = temperatures[index];
            if point.law.temperature_k() != temperature
                || point.law.law_temperature_k() != temperature
                || !positive(point.reference_volume_m3)
                || !positive(point.member_reference_mass_kg)
            {
                return Err(FrameEnthalpyError::Invalid(
                    "accepted fiber temperature and positive reference geometry must match the caloric cell",
                ));
            }
            let mass = finite(point.reference_volume_m3 * point.law.density_kg_m3())?;
            volumes[index].add(point.reference_volume_m3)?;
            masses[index].add(mass)?;
            let member = member_masses
                .entry(point.element)
                .or_insert((point.member_reference_mass_kg, Sum::default()));
            if member.0 != point.member_reference_mass_kg {
                return Err(FrameEnthalpyError::Invalid(
                    "inconsistent native member mass",
                ));
            }
            member.1.add(mass)?;
            used.insert(index);
            fibers.push(BoundFiber {
                location,
                cell: index,
                volume_m3: point.reference_volume_m3,
            });
        }
        if !positions.is_empty() || used.len() != cells.len() {
            return Err(FrameEnthalpyError::Invalid(
                "every original fiber and supplied caloric cell must participate",
            ));
        }
        for (index, bound) in cells.iter().enumerate() {
            poll(cx)?;
            if !same_measure(bound.accepted.geometry.volume_m3, volumes[index].value()?)
                || !same_measure(bound.accepted.mass_kg, masses[index].value()?)
            {
                return Err(FrameEnthalpyError::Invalid(
                    "caloric reference volume and mass must equal their complete original fiber inventory",
                ));
            }
        }
        for (_, (original, sourced)) in member_masses {
            poll(cx)?;
            if !same_measure(original, sourced.value()?) {
                return Err(FrameEnthalpyError::Invalid(
                    "original member inertia must use the same source density and reference volume",
                ));
            }
        }
        poll(cx)?;
        Ok(Self {
            frame,
            cells,
            fibers,
            maximum_fiber_points,
            time_s: thermal.time_s(),
        })
    }

    pub fn frame(&self) -> &Frame2d {
        &self.frame
    }

    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    fn check_thermal<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<(), FrameEnthalpyError<E>> {
        for bound in &self.cells {
            poll(cx)?;
            let old = &bound.accepted;
            let cell = thermal
                .cells()
                .get(&old.id)
                .ok_or(FrameEnthalpyError::Invalid(
                    "structural thermal cell became inactive",
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
                return Err(FrameEnthalpyError::Invalid(
                    "original structural cell identity, geometry, mass or accepted enthalpy changed",
                ));
            }
            bound
                .caloric
                .curve()
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(FrameEnthalpyError::Caloric)?;
        }
        Ok(())
    }

    fn enthalpies<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FrameEnthalpyError<E>> {
        self.check_thermal(thermal, accepted, cx)?;
        Ok(self
            .cells
            .iter()
            .map(|cell| thermal.cells()[&cell.accepted.id].specific_h_j_kg)
            .collect())
    }

    fn mechanics<E>(
        &self,
        enthalpies: &[f64],
        duration_s: f64,
        load: &[f64],
        mode: FrameThermalMechanics,
        tolerances: FrameCreepTolerances,
        cx: &Cx<'_>,
    ) -> Result<MechanicalProposal, FrameEnthalpyError<E>> {
        if enthalpies.len() != self.cells.len() {
            return Err(FrameEnthalpyError::Invalid(
                "one trial enthalpy per original cell required",
            ));
        }
        let mut temperatures_k = Vec::with_capacity(self.cells.len());
        let mut prepared = Vec::with_capacity(self.cells.len());
        for (cell, &h) in self.cells.iter().zip(enthalpies) {
            poll(cx)?;
            let temperature = cell
                .caloric
                .curve()
                .state_at_specific_enthalpy(h)
                .map_err(FrameEnthalpyError::Caloric)?
                .temperature_k();
            let old_temperature = cell
                .caloric
                .curve()
                .state_at_specific_enthalpy(cell.accepted.specific_h_j_kg)
                .map_err(FrameEnthalpyError::Caloric)?
                .temperature_k();
            prepared.push(
                cell.source
                    .prepare(old_temperature, temperature, duration_s)
                    .map_err(FrameEnthalpyError::Material)?,
            );
            temperatures_k.push(temperature);
        }
        let mut laws = Vec::with_capacity(self.fibers.len());
        for fiber in &self.fibers {
            poll(cx)?;
            laws.push(prepared[fiber.cell].clone());
        }
        // Thermal laws retain their own original coordinate names. The legacy
        // Norton-axis argument has no role in their successor admission.
        let frame = self
            .frame
            .with_creep_laws(&laws, "", || cx.checkpoint().is_err())
            .map_err(FrameEnthalpyError::Mechanical)?;
        let (frame, energy, dynamics, fibers, iterations) = match mode {
            FrameThermalMechanics::Quasistatic => {
                let step = frame
                    .solve_creep_step(load, tolerances, || cx.checkpoint().is_err())
                    .map_err(FrameEnthalpyError::Mechanical)?;
                (
                    step.frame,
                    step.energy,
                    None,
                    step.fibers,
                    step.equilibrium.iterations,
                )
            }
            FrameThermalMechanics::BackwardEuler => {
                let step = frame
                    .advance_creep_backward_euler(load, duration_s, tolerances, || {
                        cx.checkpoint().is_err()
                    })
                    .map_err(FrameEnthalpyError::Mechanical)?;
                (
                    step.frame,
                    step.energy.material,
                    Some(step.energy),
                    step.fibers,
                    step.iterations,
                )
            }
        };
        if fibers.len() != self.fibers.len() {
            return Err(FrameEnthalpyError::Invalid(
                "native material audit omitted an original fiber",
            ));
        }
        let mut cell_heat = vec![Sum::default(); self.cells.len()];
        for (bound, audit) in self.fibers.iter().zip(&fibers) {
            poll(cx)?;
            if bound.location != (audit.element, audit.section, audit.fiber)
                || !same_measure(bound.volume_m3, audit.reference_volume_m3)
            {
                return Err(FrameEnthalpyError::Invalid(
                    "native material heat left its original fiber inventory",
                ));
            }
            cell_heat[bound.cell]
                .add(finite(bound.volume_m3 * audit.energy.caloric_heat_j_m3())?)?;
        }
        let heat_j = cell_heat
            .into_iter()
            .map(Sum::value)
            .collect::<Result<Vec<_>, _>>()?;
        Ok(MechanicalProposal {
            frame,
            energy,
            dynamics,
            fibers,
            iterations,
            temperatures_k,
            heat_j,
        })
    }

    /// Advance the actual thermal inventory and retained nonlinear frame together.
    /// Loads are absolute endpoint forces/moments. Both mechanics modes reuse
    /// the same accepted history and physical interval throughout feedback.
    /// Only the final-temperature frame and matching transport history can be
    /// published. Native cell transport and events retain their own substeps.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        target_s: f64,
        load: &[f64],
        mode: FrameThermalMechanics,
        mechanical: FrameCreepTolerances,
        thermal_config: StepConfig,
        config: FrameEnthalpyConfig,
        cx: &Cx<'_>,
    ) -> Result<FrameEnthalpyCandidate<'a, 't, O>, FrameEnthalpyError<O::Error>> {
        poll(cx)?;
        let duration_s = target_s - self.time_s;
        if thermal.time_s() != self.time_s
            || !target_s.is_finite()
            || !positive(duration_s)
            || config.maximum_iterations == 0
            || config.maximum_iterations > MAX_FEEDBACK_ITERATIONS
            || !positive(config.temperature_tolerance_k)
            || !positive(config.energy_tolerance_j)
            || !positive(config.enthalpy_relaxation)
            || config.enthalpy_relaxation > 1.0
            || !positive(mechanical.energy_j)
            || mechanical.energy_j > config.energy_tolerance_j
            || mechanical.maximum_fiber_points > self.maximum_fiber_points
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(FrameEnthalpyError::Invalid(
                "matching clocks, positive physical interval and bounded component/coupling controls required",
            ));
        }
        let mut enthalpies = self.enthalpies(thermal, true, cx)?;
        let mut defects = [f64::INFINITY; 4];
        for iteration in 1..=config.maximum_iterations {
            poll(cx)?;
            let trial = self.mechanics(&enthalpies, duration_s, load, mode, mechanical, cx)?;
            let mut watts = BTreeMap::new();
            for (cell, &heat) in self.cells.iter().zip(&trial.heat_j) {
                poll(cx)?;
                watts.insert(cell.accepted.id, finite(heat / duration_s)?);
            }
            let heat = observer::FrameHeat {
                frame: self,
                inner: observer.clone(),
                cx,
                watts,
            };
            let (next_thermal, heat, thermal_audit) =
                thermal.advance_observed(&heat, target_s, thermal_config, || {
                    cx.checkpoint().map_err(|_| DepositionError::Cancelled)
                })?;
            let actual_h = self.enthalpies(&next_thermal, false, cx)?;
            let mut temperature_defect: f64 = 0.0;
            for ((cell, &h), &trial_temperature) in
                self.cells.iter().zip(&actual_h).zip(&trial.temperatures_k)
            {
                poll(cx)?;
                let actual_temperature = cell
                    .caloric
                    .curve()
                    .state_at_specific_enthalpy(h)
                    .map_err(FrameEnthalpyError::Caloric)?
                    .temperature_k();
                temperature_defect =
                    temperature_defect.max(finite((actual_temperature - trial_temperature).abs())?);
            }
            // A failed temperature gate cannot publish. Delay the second native
            // equilibrium solve until it can close the work account.
            if temperature_defect > config.temperature_tolerance_k
                && iteration < config.maximum_iterations
            {
                update_enthalpies(&mut enthalpies, &actual_h, config.enthalpy_relaxation, cx)?;
                continue;
            }
            let actual = self.mechanics(&actual_h, duration_s, load, mode, mechanical, cx)?;
            let applied_heat_j = sum(&trial.heat_j)?;
            let caloric_heat_j = finite(actual.energy.caloric_heat_j())?;
            let external_heat_j = finite(thermal_audit.external_heat_j - applied_heat_j)?;
            let thermal_residual_j = finite(
                thermal_audit.stored_energy_change_j
                    - thermal_audit.incoming_energy_j
                    - external_heat_j
                    - caloric_heat_j,
            )?;
            let kinetic_change = actual
                .dynamics
                .map_or(0.0, |energy| energy.kinetic_change_j);
            let inertial_numerical_loss = actual
                .dynamics
                .map_or(0.0, |energy| energy.inertial_numerical_dissipation_j);
            let combined_residual_j = finite(
                actual.energy.final_stored_j - actual.energy.initial_stored_j
                    + actual.energy.thermoelastic_storage_change_j
                    + kinetic_change
                    + thermal_audit.stored_energy_change_j
                    + actual.energy.numerical_dissipation_j
                    + inertial_numerical_loss
                    - actual.energy.applied_work_j
                    - actual.energy.support_work_j
                    - external_heat_j
                    - thermal_audit.incoming_energy_j,
            )?;
            let mut heat_defect = Sum::default();
            for (&applied, &actual_heat) in trial.heat_j.iter().zip(&actual.heat_j) {
                poll(cx)?;
                heat_defect.add(finite((applied - actual_heat).abs())?)?;
            }
            // Also close quadrature aggregation against the native whole-frame
            // account. Signed per-cell heat cannot cancel this independent gate.
            let aggregation_defect = finite((sum(&actual.heat_j)? - caloric_heat_j).abs())?;
            defects = [
                temperature_defect,
                heat_defect.value()?.max(aggregation_defect),
                thermal_residual_j,
                combined_residual_j,
            ];
            if defects[0] <= config.temperature_tolerance_k
                && defects[1] <= config.energy_tolerance_j
                && defects[2].abs() <= config.energy_tolerance_j
                && defects[3].abs() <= config.energy_tolerance_j
            {
                let report = FrameEnthalpyReport {
                    iterations: iteration,
                    mechanical_iterations: actual.iterations,
                    thermal_audit,
                    mechanical_energy: actual.energy,
                    dynamic_energy: actual.dynamics,
                    fiber_energy: actual.fibers,
                    cell_heat_j: self
                        .cells
                        .iter()
                        .zip(actual.heat_j)
                        .map(|(cell, heat)| (cell.accepted.id, heat))
                        .collect(),
                    applied_heat_j,
                    external_heat_j,
                    temperature_defect_k: defects[0],
                    heat_defect_j: defects[1],
                    thermal_residual_j,
                    combined_residual_j,
                };
                let mut next = Self {
                    frame: actual.frame,
                    cells: self.cells.clone(),
                    fibers: self.fibers.clone(),
                    maximum_fiber_points: self.maximum_fiber_points,
                    time_s: target_s,
                };
                for bound in &mut next.cells {
                    poll(cx)?;
                    bound.accepted = next_thermal.cells()[&bound.accepted.id].clone();
                }
                poll(cx)?;
                return Ok(FrameEnthalpyCandidate {
                    thermal: next_thermal,
                    frame: next,
                    observer: heat.inner,
                    report,
                });
            }
            if iteration < config.maximum_iterations {
                update_enthalpies(&mut enthalpies, &actual_h, config.enthalpy_relaxation, cx)?;
            }
        }
        poll(cx)?;
        Err(FrameEnthalpyError::NotConverged {
            iterations: config.maximum_iterations,
            temperature_defect_k: defects[0],
            heat_defect_j: defects[1],
            thermal_residual_j: defects[2],
            combined_residual_j: defects[3],
        })
    }
}

struct MechanicalProposal {
    frame: Frame2d,
    energy: FrameCreepEnergy,
    dynamics: Option<FrameCreepDynamicEnergy>,
    fibers: Vec<FrameFiberCreepEnergy>,
    iterations: usize,
    temperatures_k: Vec<f64>,
    heat_j: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct FrameEnthalpyReport {
    pub iterations: usize,
    pub mechanical_iterations: usize,
    /// Native transport input includes generated material heat. The separate
    /// external heat field below removes that internal transfer exactly once.
    pub thermal_audit: EnergyAudit,
    pub mechanical_energy: FrameCreepEnergy,
    pub dynamic_energy: Option<FrameCreepDynamicEnergy>,
    pub fiber_energy: Vec<FrameFiberCreepEnergy>,
    pub cell_heat_j: BTreeMap<CellId, f64>,
    pub applied_heat_j: f64,
    pub external_heat_j: f64,
    pub temperature_defect_k: f64,
    pub heat_defect_j: f64,
    pub thermal_residual_j: f64,
    pub combined_residual_j: f64,
}

/// Joint acceptance transfers the frame, actual thermal history and observer.
pub struct FrameEnthalpyCandidate<'m, 't, O> {
    thermal: DepositionSimulation<'t>,
    frame: ThermoplasticFrame<'m>,
    observer: O,
    report: FrameEnthalpyReport,
}

impl<'m, 't, O> FrameEnthalpyCandidate<'m, 't, O> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn frame(&self) -> &ThermoplasticFrame<'m> {
        &self.frame
    }
    pub fn observer(&self) -> &O {
        &self.observer
    }
    pub fn report(&self) -> &FrameEnthalpyReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticFrame<'m>,
        O,
        FrameEnthalpyReport,
    ) {
        (self.thermal, self.frame, self.observer, self.report)
    }
}

fn update_enthalpies<E>(
    current: &mut [f64],
    image: &[f64],
    relaxation: f64,
    cx: &Cx<'_>,
) -> Result<(), FrameEnthalpyError<E>> {
    for (h, &value) in current.iter_mut().zip(image) {
        poll(cx)?;
        *h = finite((1.0 - relaxation) * *h + relaxation * value)?;
    }
    poll(cx)
}

// Compensate grouped quadrature and signed heat without allocating histories
// or allowing large opposite exchanges to conceal an unrepresentable account.
#[derive(Clone, Copy, Debug)]
struct Sum(Compensated);
impl Default for Sum {
    fn default() -> Self {
        Self(Compensated::zero())
    }
}
impl Sum {
    fn add<E>(&mut self, value: f64) -> Result<(), FrameEnthalpyError<E>> {
        self.0 = self.0.accumulate(value);
        finite::<E>(self.0.sum)?;
        finite::<E>(self.0.comp)?;
        Ok(())
    }
    fn value<E>(self) -> Result<f64, FrameEnthalpyError<E>> {
        finite(self.0.value())
    }
}
fn sum<E>(values: &[f64]) -> Result<f64, FrameEnthalpyError<E>> {
    let mut result = Sum::default();
    for &value in values {
        result.add(value)?;
    }
    result.value()
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn same_measure(a: f64, b: f64) -> bool {
    positive(a) && positive(b) && (a - b).abs() <= 512.0 * f64::EPSILON * a.max(b)
}
fn finite<E>(value: f64) -> Result<f64, FrameEnthalpyError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FrameEnthalpyError::Invalid(
            "unrepresentable frame heat/work arithmetic",
        ))
    }
}
fn poll<E>(cx: &Cx<'_>) -> Result<(), FrameEnthalpyError<E>> {
    cx.checkpoint().map_err(|_| FrameEnthalpyError::Cancelled)
}
