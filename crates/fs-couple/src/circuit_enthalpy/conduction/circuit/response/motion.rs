//! Original nonlinear moving current on the native charge/flux timeline.
//!
//! A frozen conductor has I=V/R+I_motion. Its exact homogeneous resistance and
//! zero-voltage motional current enter the original descriptor through an
//! internal Norton port. Only actual-voltage native recovery supplies the
//! material-frame current, Joule heat and signed Lorentz work.

use super::*;
use fs_em::conduction::motion::{ConductionMotion, MovingConductionField};
use fs_math::dd::Dd;

pub(super) const TRIAL_WORK: usize = 3;

fn reserve(completed: &mut usize, maximum: usize) -> Result<(), Error> {
    *completed = completed
        .checked_add(TRIAL_WORK)
        .filter(|&next| next <= maximum)
        .ok_or(CircuitEnthalpyError::CircuitBudget)?;
    Ok(())
}

impl ConductionResponseCircuit {
    /// Initialize original RLC storage and the actual moving constitutive state.
    /// Authored charges and fluxes remain exact; no stationary operating-point
    /// admission is imposed before the moving source closure.
    ///
    /// Each constitutive trial costs two native sparse field solves and one
    /// descriptor consistency solve. All three consume maximum_work_trials
    /// before that trial starts. Motion, B, geometry and other graph values are
    /// held snapshots. Contacts retain the native stationary-adjacency rule.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_motion(
        problem: PlanarConductionResponse,
        graph: CircuitGraph,
        terminal: ConductionCircuitTerminal,
        time_s: f64,
        fluxes_wb: Vec<f64>,
        charges_c: Vec<f64>,
        ports: &[f64],
        motion: ConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        maximum_work_trials: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        admit_graph(&graph)?;
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(invalid("finite nonnegative initial circuit time required"));
        }
        let mut completed = 0usize;
        let (candidate, graph) = problem
            .solve_coupled_with_motion_and_materials_and_contacts_candidate(
                problem.conductors().to_vec(),
                problem.contacts().to_vec(),
                assembly,
                field,
                response,
                &mut || (*cancel.borrow_mut())(),
                |native, _| {
                    reserve(&mut completed, maximum_work_trials)?;
                    let affine = native
                        .moving_resistance(
                            &terminal.positive,
                            &terminal.negative,
                            motion,
                            field,
                            &mut || (*cancel.borrow_mut())(),
                        )
                        .map_err(field_error)?;
                    let graph = stamp(&graph, &terminal, affine.resistance_ohm())?;
                    let (a, b, _) = graph.branches[terminal.branch];
                    let norton = super::super::affine::NortonPorts::new(
                        &graph,
                        &[a, b],
                        &[affine.motion_current_a(), -affine.motion_current_a()],
                        ports,
                        field.current_tolerance_a,
                        &cancel,
                    )?;
                    admit_graph(&norton.graph)?;
                    let dae = assemble_circuit(&norton.graph).map_err(electrical)?;
                    let x = consistent(&dae, &fluxes_wb, &charges_c, &norton.values)?;
                    let voltage = |node| {
                        if node == 0 {
                            0.0
                        } else {
                            x[dae.node_potential_index[node - 1]]
                        }
                    };
                    let field = affine
                        .field_at_voltage(finite(voltage(a) - voltage(b))?, &mut || {
                            (*cancel.borrow_mut())()
                        })
                        .map_err(field_error)?;
                    Ok::<_, Error>((field, graph))
                },
            )?;
        let (problem, _) = candidate.into_parts();
        let (problem, sources) = problem.into_frozen_parts();
        check(&cancel)?;
        Ok(Self {
            circuit: ConductionCircuit {
                problem,
                graph,
                terminal,
                fluxes_wb,
                charges_c,
                time_s,
                material_resistors: Vec::new(),
                storage: Vec::new(),
            },
            sources,
        })
    }

    /// Propose a native midpoint RLC interval at prescribed conductor motion.
    /// All material trials start from the same original charge/flux history.
    /// Actual E'=E+u cross B and J close the original bulk/contact laws; only
    /// local irreversible heat is reported for the eventual thermal receivers.
    ///
    /// Signed Lorentz work is positive out to prescribed motion. External
    /// authored circuit work, circuit storage, graph heat, local field heat and
    /// Lorentz work are independently retained and checked. Internal Norton
    /// work is an exact lowering device, not a second physical source.
    /// maximum_work_trials bounds both preparation field solves and each
    /// descriptor step, including discarded constitutive/refinement trials.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_motion_candidate(
        &self,
        target_s: f64,
        ports: &[f64],
        motion: ConductionMotion<'_>,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        circuit: MaterialCircuitSolve,
        maximum_work_trials: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ConductionResponseCircuitCandidate, Error> {
        let cancel = RefCell::new(cancelled);
        let problem = self.response(&cancel)?;
        self.trial(
            &self.circuit,
            &problem,
            problem.conductors().to_vec(),
            problem.contacts().to_vec(),
            self.circuit.material_resistors.clone(),
            target_s,
            ports,
            Some(motion),
            maximum_work_trials,
            assembly,
            field,
            response,
            circuit,
            None,
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn trial_with_motion(
        &self,
        driver: &ConductionCircuit,
        problem: &PlanarConductionResponse,
        sources: Vec<ResolvedConductorState>,
        contact_sources: Vec<fs_em::conduction::ConductionContact>,
        material_resistors: Vec<MaterialResistor>,
        target_s: f64,
        ports: &[f64],
        motion: ConductionMotion<'_>,
        maximum_work_trials: usize,
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        response: ConductionResponseControls,
        controls: MaterialCircuitSolve,
        contact_heat: Option<&contacts::PreparedContacts<'_, '_>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<ConductionResponseCircuitCandidate, Error> {
        let mut completed = 0usize;
        let (candidate, (trial, mechanical_work_j)) = problem
            .solve_coupled_with_motion_and_materials_and_contacts_candidate(
                sources,
                contact_sources,
                assembly,
                field,
                response,
                &mut || (*cancel.borrow_mut())(),
                |native, _| {
                    reserve(&mut completed, maximum_work_trials)?;
                    driver.trial_moving_field(
                        native,
                        material_resistors.clone(),
                        target_s,
                        ports,
                        motion,
                        field,
                        controls,
                        contact_heat,
                        cancel,
                    )
                },
            )?;
        let (problem, report) = candidate.into_parts();
        let (problem, sources) = problem.into_frozen_parts();
        let native = trial.finish(driver, problem);
        check(cancel)?;
        Ok(ConductionResponseCircuitCandidate {
            circuit: Self {
                circuit: native.circuit,
                sources,
            },
            step: ConductionResponseCircuitStep {
                circuit: native.step,
                motion: Some(report.motion),
                mechanical_work_j,
                evaluated: report.field.evaluated,
                evaluated_contacts: report.field.evaluated_contacts,
                field_solves: completed,
                source_evaluations: report.field.source_evaluations,
            },
        })
    }
}

impl ConductionCircuit {
    #[allow(clippy::too_many_arguments)]
    fn trial_moving_field<'p>(
        &self,
        problem: &'p PlanarConduction,
        mut material_resistors: Vec<MaterialResistor>,
        target_s: f64,
        ports: &[f64],
        motion: ConductionMotion<'_>,
        field_controls: ConductionControls,
        controls: MaterialCircuitSolve,
        contacts: Option<&contacts::PreparedContacts<'_, '_>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(MovingConductionField<'p>, (CircuitTrialState, f64)), Error> {
        check(cancel)?;
        if contacts.is_none() && !problem.contacts().is_empty() {
            return Err(invalid(
                "spatial contact fields require their prepared source and heat owner",
            ));
        }
        let dt = finite(target_s - self.time_s)?;
        let affine = problem
            .moving_resistance(
                &self.terminal.positive,
                &self.terminal.negative,
                motion,
                field_controls,
                &mut || (*cancel.borrow_mut())(),
            )
            .map_err(field_error)?;
        let ohms = affine.resistance_ohm();
        let mut graph = stamp(&self.graph, &self.terminal, ohms)?;
        if !self.material_resistors.is_empty() {
            graph = stamp_material_resistors(&graph, &self.terminal, &mut material_resistors)?;
        }
        let (a, b, _) = graph.branches[self.terminal.branch];
        let norton = super::super::affine::NortonPorts::new(
            &graph,
            &[a, b],
            &[affine.motion_current_a(), -affine.motion_current_a()],
            ports,
            field_controls.current_tolerance_a,
            cancel,
        )?;
        admit_graph(&norton.graph)?;
        let dae = assemble_circuit(&norton.graph).map_err(electrical)?;
        let x0 = consistent(&dae, &self.fluxes_wb, &self.charges_c, &norton.values)?;
        check(cancel)?;
        let (native, _) = dae
            .step_audited(&x0, &norton.values, dt)
            .map_err(electrical)?;
        check(cancel)?;
        let stored_before_j = dae.system.hamiltonian(&x0);
        let stored_after_j = dae.system.hamiltonian(&native.x);
        if native.x.len() != x0.len()
            || native.y.len() != norton.values.len()
            || native.x.iter().chain(&native.y).any(|x| !x.is_finite())
            || [
                stored_before_j,
                stored_after_j,
                native.delta_h,
                native.dissipated,
                native.supplied,
                native.solver_residual,
            ]
            .iter()
            .any(|x| !x.is_finite())
            || stored_before_j < 0.0
            || stored_after_j < 0.0
            || native.dissipated < 0.0
        {
            return Err(invalid(
                "native moving circuit storage/loss is unrepresentable",
            ));
        }
        let voltage = |node| {
            if node == 0 {
                0.0
            } else {
                let i = dae.node_potential_index[node - 1];
                f64::midpoint(x0[i], native.x[i])
            }
        };
        let spatial_voltage_v = finite(voltage(a) - voltage(b))?;
        let spatial_current_a = finite(
            (Dd::from_f64(spatial_voltage_v / ohms) + Dd::from_f64(affine.motion_current_a()))
                .to_f64(),
        )?;
        let field = affine
            .field_at_voltage(spatial_voltage_v, &mut || (*cancel.borrow_mut())())
            .map_err(field_error)?;
        let power = field.joule_power_w();
        let terminal_power_w = finite(spatial_voltage_v * spatial_current_a)?;
        if finite(field.terminal_power_w() - terminal_power_w)?.abs() > controls.power_tolerance_w {
            return Err(invalid(
                "native moving field and circuit terminal work disagree",
            ));
        }
        let spatial_heat = finite(dt * power)?;
        if power > 0.0 && spatial_heat == 0.0 {
            return Err(invalid("moving spatial circuit heat underflows"));
        }
        let mechanical_work_j = finite(dt * field.mechanical_power_w())?;
        let spatial = ConductionThermalStep {
            start_s: self.time_s,
            end_s: target_s,
            voltage_v: field.voltage_v().to_vec(),
            electric_field_v_m: field.electric_field_v_m().to_vec(),
            current_density_a_m2: field.current_density_a_m2().to_vec(),
            cell_power_w: field.cell_joule_power_w().to_vec(),
            contacts: contacts
                .map(|contacts| contacts.heat(field.contact_response()))
                .transpose()?
                .unwrap_or_default(),
            nodal_boundary_current_a: field.nodal_boundary_current_a().to_vec(),
            nodal_applied_current_a: field.nodal_applied_current_a().to_vec(),
            boundary_current_a: field.boundary_current_a(),
            maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
            power_balance_residual_w: field.power_balance_residual_w(),
            heat_j: spatial_heat,
            supplied_work_j: finite(dt * field.terminal_power_w())?,
            linear: field.solve_report().cloned(),
        };
        let (graph_resistors, graph_heat_j) = graph_resistor_heat(
            &graph,
            self.terminal.branch,
            &material_resistors,
            voltage,
            dt,
            controls,
            cancel,
        )?;
        let heat_j = finite(spatial_heat + graph_heat_j)?;
        let stored_change_j = crate::material_circuit::capacitance::storage_change(
            &graph,
            &[],
            &dae,
            &x0,
            &native.x,
        )?;
        let work = norton.work(&native.y, dt)?;
        let supplied_work_j = work.authored_j;
        // This exact physical remainder enters the existing finite-heat driver.
        finite(supplied_work_j - mechanical_work_j)?;
        let energy_residual_j = finite(
            (Dd::from_f64(stored_change_j)
                + Dd::from_f64(heat_j)
                + Dd::from_f64(mechanical_work_j)
                - Dd::from_f64(supplied_work_j))
            .to_f64(),
        )?;
        // The native homogeneous loss plus its internal Norton port equals
        // actual terminal work. Neither quantity is itself material heat.
        let homogeneous_heat_j = finite(dt * (spatial_voltage_v / ohms) * spatial_voltage_v)?;
        let native_heat_j = finite(homogeneous_heat_j + graph_heat_j)?;
        let affine_work_residual_j = finite(
            (Dd::from_f64(spatial.supplied_work_j) - Dd::from_f64(homogeneous_heat_j)
                + Dd::from_f64(work.offset_j))
            .to_f64(),
        )?;
        if finite(native_heat_j - native.dissipated)?.abs() > controls.energy_tolerance_j
            || finite(stored_change_j - native.delta_h)?.abs() > controls.energy_tolerance_j
            || finite(work.total_j - native.supplied)?.abs() > controls.energy_tolerance_j
            || affine_work_residual_j.abs() > controls.energy_tolerance_j
            || energy_residual_j.abs() > controls.energy_tolerance_j
        {
            return Err(invalid(
                "native RLC storage, authored work, heat and Lorentz transfer do not close",
            ));
        }
        let fluxes_wb = dae.flux_index.iter().map(|&i| native.x[i]).collect();
        let charges_c = dae.charge_index.iter().map(|&i| native.x[i]).collect();
        let step = ConductionCircuitStep {
            native,
            stored_before_j,
            stored_after_j,
            stored_change_j,
            supplied_work_j,
            spatial,
            spatial_voltage_v,
            spatial_current_a,
            graph_resistors,
            heat_j,
            energy_residual_j,
        };
        check(cancel)?;
        Ok((
            field,
            (
                CircuitTrialState {
                    graph,
                    fluxes_wb,
                    charges_c,
                    material_resistors,
                    step,
                },
                mechanical_work_j,
            ),
        ))
    }
}
