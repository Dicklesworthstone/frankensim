//! Implicit pressure/temperature-dependent heat exchange at prescribed supports.
//!
//! One complete exterior reference face contacts a finite thermal support. The
//! face displacement is prescribed; its reactions supply the actual support
//! work. Native small-strain stress supplies compression `p = -n . sigma . n`
//! using the geometry owner's outward unit normal. This is a fixed-footprint,
//! compressive-contact reduction, not an opening/sliding contact solver.
//! Both body temperatures, stress, source resistance and plastic/reversible heat
//! are evaluated at the SAME endpoint. Every trial starts from the same accepted
//! histories, and only a balanced thermal/solid pair can be returned.

use super::*;
use fs_conduction::ConductionError;
use fs_conduction::interface::{InterfaceResistance, SeriesThermalResistance};
use fs_conduction::{ResistanceUncertainty, ThermalResistanceTerm};
use fs_matdb::ClaimId;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermal::source::ThermalMaterialError;
use fs_qty::Pressure;
use fs_qty::semantic::{QuantityKind, ValueForm};
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use fs_thermal::enthalpy::backward_euler::{BackwardEulerStep, StepCoefficients};
use std::cell::Cell;

/// Source-bound whole-face contact. The link's A cell is the owning tetrahedron;
/// B is a separate finite thermal support outside this solid's mechanical mesh.
/// Both caloric objects must be the original live thermal laws. A source claim
/// is explicit and remains pinned across calls, temperatures and load increments.
pub struct ThermalSupportContact<'a> {
    pub link: u64,
    pub face_nodes: [usize; 3],
    pub point: ThermalInterfacePoint<'a>,
    pub caloric: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    pub claim: ClaimId,
    /// Existing, distinct pressure coordinate, expressed in Pa. Other declared
    /// conditions (finish, medium, history) are retained unchanged.
    pub pressure_axis: &'a str,
}

#[derive(Clone, Copy, Debug)]
pub struct ContactSolve {
    pub maximum_cells: usize,
    pub maximum_mechanical_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    /// Independent caloric and complete mechanical/thermal balance gate.
    pub energy_tolerance_j: f64,
    pub yield_tolerance_pa: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum ContactError<E> {
    Invalid(&'static str),
    Mechanical(ThermoplasticError),
    State(ThermalMaterialError),
    Phase(PhaseStateError),
    Conduction(ConductionError),
    Thermal(DepositionError),
    Source(E),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
    MechanicalBudget,
    Cancelled,
}
impl<E> From<DepositionError> for ContactError<E> {
    fn from(e: DepositionError) -> Self {
        if matches!(e, DepositionError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(e)
        }
    }
}
impl<E> From<ThermoplasticError> for ContactError<E> {
    fn from(e: ThermoplasticError) -> Self {
        if matches!(e, ThermoplasticError::Cancelled) {
            Self::Cancelled
        } else {
            Self::Mechanical(e)
        }
    }
}
impl<E: std::fmt::Debug> core::fmt::Display for ContactError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "thermoplastic contact: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for ContactError<E> {}

/// Actual endpoint constitutive inputs and signed heat from A into B.
#[derive(Clone, Debug)]
pub struct ContactSample {
    pub link: u64,
    pub phases: [EquilibriumPhaseState; 2],
    pub compression_pa: f64,
    pub outward_normal: [f64; 3],
    pub area_m2: f64,
    pub resistance: InterfaceResistance,
    /// Original bulk-only coefficient when finite layers participate in series.
    /// None denotes the original interface-only support reduction.
    pub bulk_conductance_w_k: Option<f64>,
    pub conductance_w_k: f64,
    pub heat_a_to_b_j: f64,
}

pub struct ContactReport<S> {
    pub balance: ThermoplasticReport,
    pub contacts: Vec<ContactSample>,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub mechanical_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Work of prescribed support reactions, included in balance.mechanical_work_j.
    pub support_work_j: f64,
}
pub struct ContactCandidate<'m, 't, S> {
    thermal: DepositionSimulation<'t>,
    solid: ThermoplasticSolid<'m>,
    report: ContactReport<S>,
}
impl<'m, 't, S> ContactCandidate<'m, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn solid(&self) -> &ThermoplasticSolid<'m> {
        &self.solid
    }
    pub const fn report(&self) -> &ContactReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermoplasticSolid<'m>,
        ContactReport<S>,
    ) {
        (self.thermal, self.solid, self.report)
    }
}

fn number<E>(value: f64) -> Result<f64, ContactError<E>> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ContactError::Invalid(
            "nonfinite contact energy or transport",
        ))
    }
}

impl ThermalSupportContact<'_> {
    pub(super) fn sample<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        stress: &[f64; 6],
        normal: [f64; 3],
        duration_s: f64,
    ) -> Result<ContactSample, ContactError<E>> {
        let link = &thermal.links()[&self.link];
        let states = [link.a, link.b]
            .into_iter()
            .zip(self.caloric)
            .map(|(id, curve)| {
                let cell = &thermal.cells()[&id];
                if !thermal.uses_material_law(cell.material_id, curve) {
                    return Err(ContactError::Invalid(
                        "contact must retain both original caloric objects",
                    ));
                }
                curve
                    .state_at_specific_enthalpy(cell.specific_h_j_kg)
                    .map_err(ContactError::Phase)
            })
            .collect::<Result<Vec<_>, _>>()?;
        let phases = [states[0], states[1]];
        let mut point = self.point.at(phases).map_err(ContactError::State)?;
        let [nx, ny, nz] = normal;
        // Native Voigt order is xx, yy, zz, xy, yz, zx (tensor shear, not Mandel).
        let pressure = number(
            -(stress[0] * nx * nx
                + stress[1] * ny * ny
                + stress[2] * nz * nz
                + 2.0 * (stress[3] * nx * ny + stress[4] * ny * nz + stress[5] * nz * nx)),
        )?;
        if pressure < 0.0 {
            return Err(ContactError::Invalid(
                "tensile normal traction requires an opening/contact model",
            ));
        }
        point = point
            .with(self.pressure_axis, pressure)
            .map_err(|_| ContactError::Invalid("contact pressure coordinate refused"))?;
        let resistance = InterfaceResistance::from_card_pinned(
            "thermal support",
            self.point.interface,
            &point,
            self.claim,
        )
        .map_err(ContactError::Conduction)?;
        let series = SeriesThermalResistance::new(vec![
            resistance
                .term_for_area("thermal support", link.area_m2)
                .map_err(ContactError::Conduction)?,
        ])
        .map_err(ContactError::Conduction)?;
        let conductance = number(1.0 / series.budget().value_k_per_w)?;
        if conductance <= 0.0 {
            return Err(ContactError::Invalid("contact conductance underflow"));
        }
        Ok(ContactSample {
            link: self.link,
            phases,
            compression_pa: pressure,
            outward_normal: normal,
            area_m2: link.area_m2,
            resistance,
            bulk_conductance_w_k: None,
            conductance_w_k: conductance,
            heat_a_to_b_j: number(
                duration_s * conductance * (phases[0].temperature_k() - phases[1].temperature_k()),
            )?,
        })
    }
}

impl ContactSample {
    // A zero bulk coefficient explicitly insulates; missing bulk data is not
    // zero. The original resistance owner composes the two physical terms.
    pub(super) fn with_bulk<E>(mut self, bulk: f64, duration_s: f64)
        -> Result<Self, ContactError<E>> {
        if !bulk.is_finite() || bulk < 0.0 {
            return Err(ContactError::Invalid("finite nonnegative bulk conductance required"));
        }
        let conductance = if bulk == 0.0 { 0.0 } else {
            let term = ThermalResistanceTerm::declared("bulk", 1.0 / bulk,
                ResistanceUncertainty::Unstated,
                "original endpoint bulk conductance; source report retained by caller")
                .map_err(ContactError::Conduction)?;
            let series = SeriesThermalResistance::new(vec![term,
                self.resistance.term_for_area("thermal support", self.area_m2)
                    .map_err(ContactError::Conduction)?])
                .map_err(ContactError::Conduction)?;
            let g = number(1.0 / series.budget().value_k_per_w)?;
            if g <= 0.0 { return Err(ContactError::Invalid("series conductance underflow")); }
            g
        };
        self.bulk_conductance_w_k = Some(bulk);
        self.conductance_w_k = conductance;
        self.heat_a_to_b_j = number(duration_s * conductance
            * (self.phases[0].temperature_k() - self.phases[1].temperature_k()))?;
        Ok(self)
    }
}

impl<'m> ThermoplasticSolid<'m> {
    /// One backward-Euler thermal interval coupled to native quasistatic mechanics.
    /// The callback composes other endpoint material transport and external heat;
    /// it must not supply these contact links or the mechanical heat added here.
    /// Authored boundaries and unrelated links remain under the original thermal
    /// owner. Its event semantics apply: split at topology/source events first.
    /// Newton/finite-difference proposals never accumulate strain history or heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_thermal_contacts<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        contacts: &[ThermalSupportContact<'_>],
        target_s: f64,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: ContactSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cx: &Cx<'_>,
    ) -> Result<ContactCandidate<'m, 't, S>, ContactError<E>> {
        self.advance_contacts_candidate(thermal, contacts, target_s, prescribed_m,
            nodal_forces_n, settings, controls, guess_j, resolve, false, None, cx)
            .map(|(candidate, _)| candidate)
    }

    /// Combine solved support pressure with finite bulk conduction in series.
    /// The callback MUST provide a fresh bulk-only coefficient for every support
    /// link; it must not already include that interface resistance. Both source
    /// records remain available: bulk in `report.source`, interface and actual
    /// combined conductance/heat in `report.contacts`. Explicit zero bulk
    /// conductivity insulates without skipping source or geometry admission.
    ///
    /// This shares the original contact/solid/enthalpy solve and all energy gates.
    /// It does not change the fixed-footprint, compressive, prescribed-support
    /// assumptions or infer opening, friction, contact area or electrical forces.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_series_thermal_contacts<'t, S, E>(
        &self, thermal: &DepositionSimulation<'t>, contacts: &[ThermalSupportContact<'_>],
        target_s: f64, prescribed_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings, controls: ContactSolve, guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cx: &Cx<'_>,
    ) -> Result<ContactCandidate<'m, 't, S>, ContactError<E>> {
        self.advance_contacts_candidate(thermal, contacts, target_s, prescribed_m,
            nodal_forces_n, settings, controls, guess_j, resolve, true, None, cx)
            .map(|(candidate, _)| candidate)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_contacts_candidate<'t, S, E>(
        &self, thermal: &DepositionSimulation<'t>, contacts: &[ThermalSupportContact<'_>],
        target_s: f64, prescribed_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings, controls: ContactSolve, guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        series_bulk: bool, velocity_mps: Option<&[[f64; 3]]>, cx: &Cx<'_>,
    ) -> Result<(ContactCandidate<'m, 't, S>, Option<dynamics::Motion>), ContactError<E>> {
        let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        checkpoint()?;
        controls.newton.validate().map_err(ContactError::Newton)?;
        if self.time_s != thermal.time_s()
            || contacts.is_empty()
            || contacts.len() > 64
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || guess_j.len() != thermal.cells().len()
            || controls.maximum_mechanical_trials == 0
            || controls.maximum_newton_iterations == 0
            || !positive(controls.energy_tolerance_j)
            || !positive(controls.yield_tolerance_pa)
            || !controls.cell_energy_tolerance_j.is_finite()
            || controls.cell_energy_tolerance_j < 0.0
            || !positive(controls.difference_step_j)
        {
            return Err(ContactError::Invalid(
                "common clocks and bounded cells/trials/energy controls required",
            ));
        }
        self.enthalpies::<Infallible>(thermal, true, cx)?;
        let routes = self.support_contact_routes::<E>(thermal, contacts, prescribed_m, cx)?;
        self.solve_contacts_candidate(thermal, contacts, target_s, prescribed_m,
            nodal_forces_n, settings, controls, guess_j, resolve, series_bulk, &routes,
            velocity_mps, cx)
    }

    // Both powered and unpowered solids use these exact face/source checks.
    // Returned normals and indices belong to the native reference mesh, never
    // to a second caller-authored pressure map.
    pub(super) fn support_contact_routes<E>(
        &self, thermal: &DepositionSimulation<'_>, contacts: &[ThermalSupportContact<'_>],
        prescribed_m: &[(usize, f64)], cx: &Cx<'_>,
    ) -> Result<Vec<(usize, [f64; 3])>, ContactError<E>> {
        let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        checkpoint()?;
        if contacts.len() > 64 {
            return Err(ContactError::Invalid("bounded support contact count required"));
        }
        let mut links = BTreeSet::new();
        let mut faces = BTreeSet::new();
        let mut routes = Vec::with_capacity(contacts.len());
        for binding in contacts {
            checkpoint()?;
            let mut nodes = binding.face_nodes;
            nodes.sort_unstable();
            let face = self
                .geometry
                .faces()
                .iter()
                .find(|face| face.nodes == nodes)
                .ok_or(ContactError::Invalid(
                    "thermal contact must name a complete exterior face",
                ))?;
            let element = self
                .tetrahedra
                .iter()
                .position(|tet| nodes.iter().all(|n| tet.contains(n)))
                .ok_or(ContactError::Invalid(
                    "contact face has no structural owner",
                ))?;
            let link = thermal
                .links()
                .get(&binding.link)
                .ok_or(ContactError::Invalid("inactive contact link"))?;
            if !links.insert(binding.link)
                || !faces.insert(nodes)
                || link.a != self.cells[element].accepted.id
                || self.cells.iter().any(|cell| cell.accepted.id == link.b)
                || !super::super::same_measure(link.area_m2, face.area_m2)
                || nodes
                    .iter()
                    .any(|n| (0..3).any(|c| !prescribed_m.iter().any(|(d, _)| *d == 3 * n + c)))
            {
                return Err(ContactError::Invalid(
                    "unique whole-face link to finite prescribed support required",
                ));
            }
            let axis = binding.pressure_axis;
            if axis.trim().is_empty()
                || !binding.point.context.axes().contains_key(axis)
                || binding
                    .point
                    .context
                    .axis_quantities()
                    .get(axis)
                    .is_some_and(|q| q.dims() != Pressure::DIMS)
                || binding.point.temperature_axes.contains(&axis)
                || binding.point.liquid_fraction_axes.contains(&Some(axis))
                || binding.point.mean_temperature_axis == Some(axis)
            {
                return Err(ContactError::Invalid(
                    "existing distinct pressure coordinate in pascals required",
                ));
            }
            let claim = binding
                .point
                .interface
                .claims()
                .claim(binding.claim)
                .ok_or(ContactError::Invalid(
                    "pinned thermal interface claim absent",
                ))?;
            for (axis, value) in binding.point.context.axes() {
                let semantic = binding
                    .point
                    .context
                    .axis_quantities()
                    .get(axis)
                    .and_then(|q| q.semantic_type());
                if axis == "omega"
                    || semantic.is_some_and(|s| matches!(s.kind(), QuantityKind::Frequency(_)))
                {
                    if *value != 0.0
                        || semantic.is_some_and(|s| s.form() != ValueForm::Static)
                        || axis == binding.pressure_axis
                    {
                        return Err(ContactError::Invalid(
                            "thermal support resistance requires a static source",
                        ));
                    }
                }
            }
            if !claim.value.axes().all(|(name, dims)| {
                (name == axis && dims == Pressure::DIMS)
                    || binding.point.is_thermal_axis(name, dims)
            }) {
                return Err(ContactError::Invalid(
                    "contact resistance must follow declared static pressure/thermal coordinates",
                ));
            }
            // Admission uses actual accepted stress/temperatures, never a stale
            // pressure copied from the card context. Trials retain the same pin.
            binding.sample::<E>(
                thermal,
                &self.solution.element_stress[element],
                face.outward_normal,
                0.0,
            )?;
            routes.push((element, face.outward_normal));
        }
        Ok(routes)
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_contacts_candidate<'t, S, E>(
        &self, thermal: &DepositionSimulation<'t>, contacts: &[ThermalSupportContact<'_>],
        target_s: f64, prescribed_m: &[(usize, f64)], nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings, controls: ContactSolve, guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        series_bulk: bool, routes: &[(usize, [f64; 3])],
        velocity_mps: Option<&[[f64; 3]]>, cx: &Cx<'_>,
    ) -> Result<(ContactCandidate<'m, 't, S>, Option<dynamics::Motion>), ContactError<E>> {
        let checkpoint = || cx.checkpoint().map_err(|_| DepositionError::Cancelled);
        let mechanical_config = ThermoplasticConfig {
            maximum_iterations: 1,
            temperature_tolerance_k: f64::EPSILON,
            energy_tolerance_j: controls.energy_tolerance_j,
            yield_tolerance_pa: controls.yield_tolerance_pa,
            enthalpy_relaxation: 1.0,
        };
        let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
        let dt = step.duration_s();
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            step.evaluate(
                x,
                |endpoint| {
                    checkpoint()?;
                    if trials.get() >= controls.maximum_mechanical_trials {
                        return Err(ContactError::MechanicalBudget);
                    }
                    trials.set(trials.get() + 1);
                    let h = self.enthalpies::<Infallible>(endpoint, false, cx)?;
                    let mechanical = self.mechanics_with_motion::<Infallible>(
                        &h,
                        prescribed_m,
                        nodal_forces_n,
                        settings,
                        mechanical_config,
                        dt,
                        velocity_mps,
                        cx,
                    )?;
                    let mut coefficients = resolve(endpoint).map_err(ContactError::Source)?;
                    let mut samples = Vec::with_capacity(contacts.len());
                    for (binding, &(element, normal)) in contacts.iter().zip(routes) {
                        checkpoint()?;
                        let mut sample = binding.sample(
                            endpoint,
                            &mechanical.solution.element_stress[element],
                            normal,
                            dt,
                        )?;
                        if series_bulk {
                            let bulk = *coefficients.conductances_w_k.get(&binding.link)
                                .ok_or(ContactError::Invalid("support lacks explicit bulk-only conductance"))?;
                            sample = sample.with_bulk(bulk, dt)?;
                        }
                        if coefficients
                            .conductances_w_k
                            .insert(binding.link, sample.conductance_w_k)
                            .is_some()
                            && !series_bulk
                        {
                            return Err(ContactError::Invalid(
                                "contact conductance supplied twice",
                            ));
                        }
                        samples.push(sample);
                    }
                    let external_heat = number(
                        dt * endpoint.external_step_power_w(
                            &coefficients.sources_w,
                            &coefficients.conductances_w_k,
                            checkpoint,
                        )?,
                    )?;
                    let mut applied_heat = 0.0;
                    for (cell, heat) in self.cells.iter().zip(&mechanical.heat_j) {
                        checkpoint()?;
                        let watts = number(*heat / dt)?;
                        if *heat != 0.0 && watts == 0.0 {
                            return Err(ContactError::Invalid("mechanical heat underflow"));
                        }
                        applied_heat = number(applied_heat + number(dt * watts)?)?;
                        let source = coefficients
                            .sources_w
                            .entry(cell.accepted.id)
                            .or_insert(0.0);
                        *source = number(*source + watts)?;
                    }
                    Ok::<_, ContactError<E>>(StepCoefficients {
                        sources_w: coefficients.sources_w,
                        conductances_w_k: coefficients.conductances_w_k,
                        payload: (
                            mechanical,
                            samples,
                            coefficients.payload,
                            applied_heat,
                            external_heat,
                        ),
                    })
                },
                checkpoint,
            )
        };
        let numerical = |e| match e {
            FiniteDifferenceError::Source(e) => e,
            FiniteDifferenceError::Invalid(s) => ContactError::Invalid(s),
        };
        let problem = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|trial| trial.residual_j().to_vec()),
        )
        .map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem, guess_j.to_vec(), controls.newton);
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        let mut solver = construction.map_err(ContactError::Newton)?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                checkpoint().is_err() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(ContactError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(ContactError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (
            thermal,
            (mechanical, contacts, source, applied_heat_j, external_heat_j),
            thermal_audit,
        ) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.energy_tolerance_j,
            checkpoint,
        )?;
        let caloric_heat = number(mechanical.heat_j.iter().sum())?;
        let thermal_residual_j =
            number(thermal_audit.stored_energy_change_j - external_heat_j - caloric_heat)?;
        let combined_residual_j = number(
            mechanical.stored_change_j
                + mechanical.thermoelastic_storage_change_j
                + thermal_audit.stored_energy_change_j
                + mechanical.numerical_loss_j
                + mechanical.motion.as_ref().map_or(0.0, |motion|
                    motion.energy.kinetic_change_j + motion.energy.numerical_loss_j)
                - mechanical.external_work_j
                - external_heat_j,
        )?;
        if thermal_residual_j.abs() > controls.energy_tolerance_j
            || combined_residual_j.abs() > controls.energy_tolerance_j
            || number(applied_heat_j - caloric_heat)?.abs() > controls.energy_tolerance_j
        {
            return Err(ContactError::Invalid(
                "thermal, support work and complete material energy accounts do not close",
            ));
        }
        let old_u = self.solution.displacement_m.as_flattened();
        let new_u = mechanical.solution.displacement_m.as_flattened();
        let support_work_j = number(
            prescribed_m
                .iter()
                .zip(&mechanical.solution.increments[0].reactions_n)
                .map(|(&(dof, _), &reaction)| reaction * (new_u[dof] - old_u[dof]))
                .sum(),
        )?;
        let balance = ThermoplasticReport {
            iterations: result.report.iterations,
            thermal_audit,
            element_energy: mechanical.energy,
            physical_heat_j: mechanical.physical_heat_j,
            reversible_heat_j: mechanical.reversible_heat_j,
            thermoelastic_storage_change_j: mechanical.thermoelastic_storage_change_j,
            coefficient_storage_change_j: mechanical.coefficient_storage_change_j,
            applied_heat_j,
            external_heat_j,
            stored_mechanical_change_j: mechanical.stored_change_j,
            numerical_loss_j: mechanical.numerical_loss_j,
            mechanical_work_j: mechanical.external_work_j,
            mechanical_residual_j: mechanical.mechanical_residual_j,
            thermal_residual_j,
            combined_residual_j,
            // No staggered endpoint exists: both owners use identical unknowns.
            temperature_defect_k: 0.0,
            heat_defect_j: number(applied_heat_j - caloric_heat)?.abs(),
        };
        let mut solid = self.clone();
        solid.time_s = target_s;
        solid.materials = mechanical.materials;
        solid.solution = mechanical.solution;
        for cell in &mut solid.cells {
            cell.accepted = thermal.cells()[&cell.accepted.id].clone();
        }
        checkpoint()?;
        Ok((ContactCandidate {
            thermal,
            solid,
            report: ContactReport {
                balance,
                contacts,
                source,
                nonlinear: result.report,
                mechanical_trials: trials.get(),
                maximum_cell_residual_j,
                support_work_j,
            },
        }, mechanical.motion))
    }
}
