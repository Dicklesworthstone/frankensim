//! Reciprocal finite-cell heat on the original inertial tetrahedral trajectory.
//!
//! Native average-acceleration dynamics supplies the accepted motion and endpoint
//! stress work. Each cell receives Q_rev = W_stress - delta(Psi + T S), using
//! the same source-derived potential as the quasistatic owner. Kinetic energy
//! remains mechanical storage. Rayleigh loss is explicitly exported; no
//! unsupported spatial deposition of that system-level damping law is invented.
//!
//! Trials restart both owners from their accepted history. The final detached
//! native checkpoint retains matrices, motion, loads and step count without a
//! second solve or a self-referential material lifetime. The scheme uses the
//! existing endpoint coupling and requires time refinement; an energy identity
//! alone is not a continuum or reversible-entropy accuracy certificate.
use super::*;
use fs_solid::loads3::transient::{
    TransientConfig3, response::TransientResponse3, checkpoint::TetTransientCheckpoint3,
};

/// Original source laws plus the actual loaded inertial and caloric history.
pub struct DynamicThermoelasticSolid<'s, 'm> {
    pub(super) source: &'s SolidEnthalpy<'m>,
    laws: Vec<AffineIsotropicThermoelastic>,
    state: SolidEnthalpyState<'m>,
    motion: TetTransientCheckpoint3<'s, 'm>,
    response: TransientResponse3,
    potentials: Vec<ThermoelasticState>,
    fixed_dofs: Vec<usize>,
}

impl<'s, 'm> ThermoelasticSolid<'s, 'm> {
    /// Launch actual inertia from this admitted loaded equilibrium and supplied
    /// nodal velocity. The load must be the original equilibrium load. Its
    /// subsequent changes use native endpoint work, not equilibrium re-solving.
    /// This transition retains the same caloric chart, source intervals, strain
    /// envelope, reference density, geometry and stationary support components.
    pub fn begin_transient(
        self, loads: &TetLoadCase3<'s, 'm>, velocity_m_s: &[[f64; 3]],
        config: TransientConfig3, cx: &Cx<'_>,
    ) -> Result<DynamicThermoelasticSolid<'s, 'm>, ThermoelasticFeedbackError> {
        check(cx)?;
        if loads.nodal_forces_n() != self.nodal_forces_n {
            return Err(ThermoelasticFeedbackError::Invalid("initial transient load must match the admitted loaded equilibrium"));
        }
        let problem = self.state.problem(&self.fixed_dofs);
        let run = loads.begin_transient(&problem, Some(self.state.thermal_strains()),
            self.state.time_s(), self.equilibrium.solution.displacement.displacement_m(),
            velocity_m_s, config, cx).map_err(SolidEnthalpyError::from)?;
        let response = run.recover_response(self.state.mesh_frame,
            self.state.maximum_absolute_mandel_strain, cx).map_err(SolidEnthalpyError::from)?;
        let motion = run.into_checkpoint(cx).map_err(SolidEnthalpyError::from)?;
        check(cx)?;
        Ok(DynamicThermoelasticSolid {
            source: self.source, laws: self.laws, state: self.state, motion, response,
            potentials: self.potentials, fixed_dofs: self.fixed_dofs,
        })
    }
}

impl<'s, 'm> DynamicThermoelasticSolid<'s, 'm> {
    pub fn state(&self) -> &SolidEnthalpyState<'m> { &self.state }
    pub fn response(&self) -> &TransientResponse3 { &self.response }
    pub fn potentials(&self) -> &[ThermoelasticState] { &self.potentials }
    pub fn time_s(&self) -> f64 { self.state.time_s() }
    pub fn step(&self) -> usize { self.motion.state().step }
    pub fn kinetic_energy_j(&self) -> f64 { self.motion.state().kinetic_energy_j }

    /// Local entropy-constrained tangent at the actual accepted inertial state.
    /// This does not advance motion or add heat feedback a second time.
    pub fn adiabatic_state(
        &self, cx: &Cx<'_>,
    ) -> Result<super::acoustic::AdiabaticSolidState<'m>, ThermoelasticFeedbackError> {
        super::acoustic::resolve(self.source, &self.laws, &self.state, &self.response.stress, cx)
    }

    /// Adiabatic constitutive tangent plus the ACTUAL inertial endpoint stress.
    /// No equilibrium solve replaces the motion, no force/stress is resupplied,
    /// and local reversible heat is not applied a second time.
    ///
    /// The returned signed pencil is a frozen conservative snapshot of a
    /// possibly evolving trajectory, not a claim of stationary normal modes.
    pub fn adiabatic_prestressed_state(
        &self, cx: &Cx<'_>,
    ) -> Result<super::acoustic::AdiabaticPrestressedSolidState<'m>, ThermoelasticFeedbackError> {
        let material_state = self.adiabatic_state(cx)?;
        let original = self.state.problem(&self.fixed_dofs);
        let run = self.motion.bind(&original, cx).map_err(SolidEnthalpyError::from)?;
        let tangent = run.assemble_prestressed_constitutive_tangent(
            TetMaterialField::PerElement(material_state.materials()),
            self.state.mesh_frame, self.state.maximum_absolute_mandel_strain, cx,
        ).map_err(SolidEnthalpyError::from)?;
        super::acoustic::with_prestress(material_state, tangent)
    }

    /// Full accepted displacement, velocity and acceleration in native node order.
    pub fn nodal_motion(
        &self, cx: &Cx<'_>,
    ) -> Result<[Vec<[f64; 3]>; 3], ThermoelasticFeedbackError> {
        let problem = self.state.problem(&self.fixed_dofs);
        let run = self.motion.bind(&problem, cx).map_err(SolidEnthalpyError::from)?;
        Ok(run.nodal_motion(cx).map_err(SolidEnthalpyError::from)?)
    }

    fn admit<E>(
        &self, thermal: &DepositionSimulation<'_>, accepted: bool, cx: &Cx<'_>,
    ) -> Result<SolidEnthalpyState<'m>, ThermoelasticFeedbackError<E>> {
        check(cx)?;
        let state = self.source.resolve(thermal, cx)?;
        if accepted && (state.time_s.to_bits() != self.state.time_s.to_bits()
            || state.materials != self.state.materials || state.strains != self.state.strains
            || state.samples.iter().zip(&self.state.samples).any(|(a, b)| a.phase != b.phase)) {
            return Err(ThermoelasticFeedbackError::Invalid("accepted inertial and caloric histories must remain exactly paired"));
        }
        for (law, sample) in self.laws.iter().zip(&state.samples) {
            check(cx)?;
            law.at_materials(&sample.elastic, sample.expansion.current())?;
        }
        Ok(state)
    }

    /// Advance real Newmark inertia and original spatial heat as one candidate.
    /// Every heat iteration starts from the same accepted mechanical checkpoint
    /// and thermal history. Cancellation, range exits, strain limits, nonlinear
    /// refusal or either native energy gate publish neither owner.
    ///
    /// The initial interval is seeded from a loaded static equilibrium, but
    /// later states are fully dynamic. Supplied Rayleigh damping is retained as
    /// explicitly exported energy, not silently converted into local heating.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t, O: DepositionObserver>(
        &self, thermal: &DepositionSimulation<'t>, observer: &O, target_s: f64,
        loads: &TetLoadCase3<'_, '_>, thermal_config: StepConfig,
        config: ThermoelasticSolve, cx: &Cx<'_>,
    ) -> Result<DynamicThermoelasticCandidate<'s, 'm, 't, O>, ThermoelasticFeedbackError<O::Error>> {
        check(cx)?;
        controls(config)?;
        let dt = target_s - self.time_s();
        if !dt.is_finite() || dt <= 0.0 || !target_s.is_finite()
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j < 0.0
            || thermal_config.energy_tolerance_j > config.energy_tolerance_j {
            return Err(ThermoelasticFeedbackError::Invalid("one positive physical interval and compatible thermal energy tolerance required"));
        }
        self.admit(thermal, true, cx)?;
        let accepted_problem = self.state.problem(&self.fixed_dofs);
        let accepted = self.motion.bind(&accepted_problem, cx).map_err(SolidEnthalpyError::from)?;
        let mut applied = vec![0.0; self.source.cells.len()];
        let mut defects = [f64::INFINITY; 2];
        for iteration in 1..=config.maximum_iterations {
            check(cx)?;
            let mut watts = BTreeMap::new();
            for (bound, &heat) in self.source.cells.iter().zip(&applied) {
                let power = finite(heat / dt)?;
                if heat != 0.0 && power == 0.0 {
                    return Err(ThermoelasticFeedbackError::Invalid("reversible heat underflows its power conversion"));
                }
                if power != 0.0 { watts.insert(bound.source.cell, power); }
            }
            let coupled = DynamicHeat { solid: self, inner: observer.clone(), cx, watts };
            let (thermal_next, coupled, audit) = thermal.advance_observed(
                &coupled, target_s, thermal_config,
                || cx.checkpoint().map_err(|_| DepositionError::Cancelled))?;
            let state = self.admit(&thermal_next, false, cx)?;
            let problem = state.problem(&self.fixed_dofs);
            let (run, native) = state.advance_transient_candidate(
                &self.state, &accepted, &problem, loads, cx)?;
            let response = run.recover_response(state.mesh_frame,
                state.maximum_absolute_mandel_strain, cx).map_err(SolidEnthalpyError::from)?;
            let potentials = dynamic_potentials(&state, &response, &self.laws, config, cx)?;
            let mut cells = Vec::with_capacity(self.source.cells.len());
            let mut stress_work = 0.0;
            let mut internal_change = 0.0;
            let mut reversible_heat = 0.0;
            for (i, bound) in self.source.cells.iter().enumerate() {
                check(cx)?;
                let old = &self.response.stress.elements()[i];
                let next = &response.stress.elements()[i];
                let a = self.potentials[i];
                let b = potentials[i];
                let free0 = self.state.strains[i].free_strain_mandel();
                let free1 = state.strains[i].free_strain_mandel();
                let sigma0 = old.stress_mandel_pa();
                let sigma1 = next.stress_mandel_pa();
                let e0 = old.elastic_strain_mandel();
                let e1 = next.elastic_strain_mandel();
                let volume = bound.reference_volume_m3;
                let mechanical_work_j = finite(0.5 * volume * (0..6).map(|k|
                    (sigma0[k] + sigma1[k]) * (e1[k] + free1[k] - e0[k] - free0[k])).sum::<f64>())?;
                let expansion_work_j = finite(-0.5 * volume * (0..6).map(|k|
                    (sigma0[k] + sigma1[k]) * (free1[k] - free0[k])).sum::<f64>())?;
                let c0 = self.state.materials[i].stiffness_mandel_pa();
                let c1 = state.materials[i].stiffness_mandel_pa();
                let stiffness_work_j = finite(0.5 * volume * (0..6).map(|row|
                    e0[row] * (0..6).map(|col|
                        (c1[row][col] - c0[row][col]) * e1[col]).sum::<f64>()).sum::<f64>())?;
                let free_energy_residual_j = finite(b.free_energy_j - a.free_energy_j
                    - mechanical_work_j - expansion_work_j - stiffness_work_j)?;
                let change = finite(b.internal_energy_j - a.internal_energy_j)?;
                let reversible_heat_j = finite(mechanical_work_j - change)?;
                let temperature_k = [
                    self.state.samples[i].phase.temperature_k(), state.samples[i].phase.temperature_k()];
                let potential_residual = finite(reversible_heat_j + expansion_work_j + stiffness_work_j
                    + temperature_k[1] * b.entropy_j_k - temperature_k[0] * a.entropy_j_k)?;
                if free_energy_residual_j.abs() > config.energy_tolerance_j
                    || potential_residual.abs() > config.energy_tolerance_j {
                    return Err(ThermoelasticFeedbackError::Invalid("dynamic source potential and native stress/material work disagree"));
                }
                stress_work = finite(stress_work + mechanical_work_j)?;
                internal_change = finite(internal_change + change)?;
                reversible_heat = finite(reversible_heat + reversible_heat_j)?;
                cells.push(ThermoelasticCellStep {
                    cell: bound.source.cell, temperature_k, potentials: [a, b],
                    mechanical_work_j, expansion_work_j, stiffness_work_j,
                    reversible_heat_j, free_energy_residual_j,
                });
            }
            let kinetic_energy_change_j = finite(run.state().kinetic_energy_j - self.kinetic_energy_j())?;
            let mechanical_work_j = finite(native.native.mechanical_work_j + native.native.support_work_j)?;
            let damping_export_j = native.native.damping_loss_j;
            let mechanical_residual_j = finite(kinetic_energy_change_j + stress_work
                + damping_export_j - mechanical_work_j)?;
            if mechanical_residual_j.abs() > config.energy_tolerance_j {
                return Err(ThermoelasticFeedbackError::Invalid("native stress work, kinetic energy, damping and applied work disagree"));
            }
            let applied_heat_j = finite(applied.iter().sum())?;
            let heat_defect_j = finite(applied.iter().zip(&cells)
                .map(|(used, cell)| (used - cell.reversible_heat_j).abs()).sum())?;
            let external_heat_j = finite(audit.external_heat_j - applied_heat_j)?;
            let thermal_residual_j = finite(audit.stored_energy_change_j - audit.incoming_energy_j
                - external_heat_j - reversible_heat)?;
            let combined_residual_j = finite(kinetic_energy_change_j + internal_change
                + audit.stored_energy_change_j + damping_export_j - audit.incoming_energy_j
                - mechanical_work_j - external_heat_j)?;
            defects = [heat_defect_j, combined_residual_j];
            if heat_defect_j <= config.heat_tolerance_j
                && thermal_residual_j.abs() <= config.energy_tolerance_j
                && combined_residual_j.abs() <= config.energy_tolerance_j {
                let motion = run.into_checkpoint(cx).map_err(SolidEnthalpyError::from)?;
                let report = DynamicThermoelasticReport {
                    iterations: iteration, thermal_audit: audit, native, cells,
                    kinetic_energy_change_j, internal_energy_change_j: internal_change,
                    mechanical_work_j, damping_export_j, reversible_heat_j: reversible_heat,
                    applied_heat_j, external_heat_j, mechanical_residual_j,
                    thermal_residual_j, combined_residual_j, heat_defect_j,
                };
                let solid = Self {
                    source: self.source, laws: self.laws.clone(), state, motion, response,
                    potentials, fixed_dofs: self.fixed_dofs.clone(),
                };
                check(cx)?;
                return Ok(DynamicThermoelasticCandidate {
                    thermal: thermal_next, solid, observer: coupled.inner, report,
                });
            }
            for (value, cell) in applied.iter_mut().zip(&cells) {
                *value = finite(config.heat_relaxation.mul_add(cell.reversible_heat_j - *value, *value))?;
            }
        }
        check(cx)?;
        Err(ThermoelasticFeedbackError::NotConverged {
            iterations: config.maximum_iterations, heat_defect_j: defects[0],
            combined_residual_j: defects[1],
        })
    }
}

fn dynamic_potentials<E>(
    state: &SolidEnthalpyState<'_>, response: &TransientResponse3,
    laws: &[AffineIsotropicThermoelastic], config: ThermoelasticSolve, cx: &Cx<'_>,
) -> Result<Vec<ThermoelasticState>, ThermoelasticFeedbackError<E>> {
    let mut potentials = Vec::with_capacity(laws.len());
    for (i, (law, sample)) in laws.iter().zip(&state.samples).enumerate() {
        check(cx)?;
        let source = law.at_materials(&sample.elastic, sample.expansion.current())?;
        let native = &response.stress.elements()[i];
        let free = state.strains[i].free_strain_mandel();
        let expected = source.free_strain()?;
        for k in 0..6 {
            let target = if k < 3 { expected } else { 0.0 };
            if (free[k] - target).abs() > 4096.0 * f64::EPSILON * law.maximum_absolute_mandel_strain() {
                return Err(ThermoelasticFeedbackError::Invalid("native integrated expansion differs from the dynamic source potential"));
            }
        }
        let strain = core::array::from_fn(|k| native.elastic_strain_mandel()[k] + free[k]);
        let potential = source.state(strain, sample.reference_volume_m3)?;
        if (potential.free_energy_j - native.elastic_energy_j()).abs() > config.energy_tolerance_j
            || potential.stress_mandel_pa.iter().zip(native.stress_mandel_pa())
                .any(|(a, b)| (a - b).abs() > config.stress_tolerance_pa) {
            return Err(ThermoelasticFeedbackError::Invalid("native dynamic stress/energy differs from its source potential"));
        }
        potentials.push(potential);
    }
    Ok(potentials)
}

#[derive(Debug)]
pub struct DynamicThermoelasticReport {
    pub iterations: usize,
    /// Original thermal audit includes the reversible source exactly once.
    pub thermal_audit: EnergyAudit,
    pub native: ThermalTransientStep3,
    pub cells: Vec<ThermoelasticCellStep>,
    pub kinetic_energy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub mechanical_work_j: f64,
    /// Explicitly exported Rayleigh energy; NOT deposited into any thermal cell.
    pub damping_export_j: f64,
    pub reversible_heat_j: f64,
    pub applied_heat_j: f64,
    pub external_heat_j: f64,
    pub mechanical_residual_j: f64,
    pub thermal_residual_j: f64,
    pub combined_residual_j: f64,
    pub heat_defect_j: f64,
}

pub struct DynamicThermoelasticCandidate<'s, 'm, 't, O> {
    pub thermal: DepositionSimulation<'t>,
    pub solid: DynamicThermoelasticSolid<'s, 'm>,
    pub observer: O,
    pub report: DynamicThermoelasticReport,
}
impl<'s, 'm, 't, O> DynamicThermoelasticCandidate<'s, 'm, 't, O> {
    pub fn into_parts(self) -> (
        DepositionSimulation<'t>, DynamicThermoelasticSolid<'s, 'm>, O, DynamicThermoelasticReport,
    ) { (self.thermal, self.solid, self.observer, self.report) }
}

#[derive(Clone)]
struct DynamicHeat<'a, 's, 'm, 'cx, O> {
    solid: &'a DynamicThermoelasticSolid<'s, 'm>,
    inner: O,
    cx: &'a Cx<'cx>,
    watts: BTreeMap<CellId, f64>,
}
impl<O: DepositionObserver> DepositionObserver for DynamicHeat<'_, '_, '_, '_, O> {
    type Error = ThermoelasticFeedbackError<O::Error>;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner.begin().map_err(ThermoelasticFeedbackError::Observer)
    }
    fn requires_step_interval(&self) -> bool { self.inner.requires_step_interval() }
    fn step_sources(&mut self, thermal: &DepositionSimulation<'_>, start: f64, dt: f64, end: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner.step_sources(thermal, start, dt, end, check).map_err(ThermoelasticFeedbackError::Observer)
    }
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.solid.admit(thermal, false, self.cx)?;
        let mut rates = self.inner.sources(thermal, checkpoint).map_err(ThermoelasticFeedbackError::Observer)?;
        for (&cell, &power) in &self.watts {
            checkpoint()?;
            let rate = rates.entry(cell).or_insert(0.0);
            *rate = finite(*rate + power)?;
        }
        Ok(rates)
    }
    fn conductances(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner.conductances(thermal, check).map_err(ThermoelasticFeedbackError::Observer)
    }
    fn heat_transfers(&mut self, thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner.heat_transfers(thermal, check).map_err(ThermoelasticFeedbackError::Observer)
    }
    fn limit_timestep_s(&mut self, thermal: &DepositionSimulation<'_>, proposed: f64,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<f64, Self::Error> {
        self.inner.limit_timestep_s(thermal, proposed, check).map_err(ThermoelasticFeedbackError::Observer)
    }
    fn event(&mut self, event: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.inner.event(event, thermal).map_err(ThermoelasticFeedbackError::Observer)?;
        self.solid.admit(thermal, false, self.cx).map(|_| ())
    }
    fn step(&mut self, start: f64, dt: f64, old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Self::Error> {
        self.inner.step(start, dt, old, thermal, check).map_err(ThermoelasticFeedbackError::Observer)?;
        self.solid.admit(thermal, false, self.cx).map(|_| ())
    }
}
