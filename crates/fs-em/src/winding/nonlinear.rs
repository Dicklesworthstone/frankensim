//! Native circuit terminals coupled to the actual nonlinear magnetic field.
//!
//! Each stranded winding owns a signed spatial J/current profile and material
//! resistance. Its original inductor branch becomes a PRIVATE current-source
//! port in the original Kirchhoff DAE. The input is -I (current leaving the
//! first terminal), and its native conjugate output is the winding voltage V.
//! The shared Newton owner solves V = delta(lambda)/dt + R*I, where lambda is
//! the exact P1 field/profile pairing plus declared series end inductance.
//!
//! Every trial restarts the accepted circuit and nonlinear-eddy histories.
//! Native graph work -dt*I*V cancels terminal work into the field/winding. Core
//! eddy heat, stranded Joule heat, graph heat and backward-Euler numerical loss
//! remain separate. Strand cells cannot also own resolved eddy conductivity.
//! This is a fixed-geometry, fixed-H(B), reversible isotropic MQS lane. Only
//! supported same-source conductivity temperature refresh is allowed. No core
//! hysteresis, moving conductor, magnetic temperature work or fitted L is used.
//! The allocating native DAE must admit the held private-current trials; ideal
//! current-only constraint topologies can refuse this bounded partitioned lane.

use super::ohmic::{
    OhmicWindingError, WindingCellHeat, WindingConductor, WindingGraphResistorHeat,
};
use super::{MAX_WINDINGS, WindingError, admit_profiles, ohmic, pairing_on_mesh};
use crate::planar::{
    MagneticError,
    eddy::{
        EddyConductor,
        nonlinear::{
            NonlinearEddyControls, NonlinearEddyError, NonlinearEddyStep, PlanarNonlinearEddy,
            PlanarNonlinearEddyCandidate,
        },
    },
    nonlinear::IsotropicMagneticMaterial,
};
use fs_material::field_storage::admit_fixed_magnetic_refresh;
use fs_phs::{
    StepRecord,
    circuit::{Branch, CircuitDae, CircuitError, CircuitGraph, assemble_circuit},
};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
    finite_difference::{FallibleFiniteDifference, FiniteDifferenceError},
};
use std::{
    cell::{Cell, RefCell},
    rc::Rc,
};

/// Fixed spatial winding geometry, in terminal-current order.
#[derive(Clone, Debug)]
pub struct NonlinearWindings {
    /// Distinct original graph inductor branches replaced by these fields.
    /// Other inductors retain their original native lumped storage.
    pub branch_indices: Vec<usize>,
    /// Signed J_z / terminal current [1/m2], including turns and return paths.
    pub profiles_per_m2: Vec<Vec<f64>>,
    /// Explicit unrepresented end/leakage self-inductance [H], including zero.
    pub end_inductance_h: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct NonlinearWindingConfig {
    pub max_windings: usize,
    pub max_profile_values: usize,
    pub current_balance_tolerance_turns: f64,
}

/// Original native storage coordinates and physical field-current history.
#[derive(Clone, Debug)]
pub struct NonlinearWindingInitial {
    pub currents_a: Vec<f64>,
    /// Unbound original inductor branches, in original branch order [Wb].
    pub remaining_fluxes_wb: Vec<f64>,
    /// Original capacitor branches, in original branch order [C].
    pub charges_c: Vec<f64>,
    /// Original source ports only, including any gaps in their index order.
    pub external_inputs: Vec<f64>,
}

/// Explicit physical acceptance and nested work budgets.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearWindingControls {
    pub field: NonlinearEddyControls,
    /// Each trial performs one native circuit and one nonlinear field step.
    pub maximum_field_trials: usize,
    pub maximum_newton_iterations: usize,
    /// Maximum component perturbation for the shared numerical Jacobian [A].
    pub difference_step_a: f64,
    pub terminal_voltage_tolerance_v: f64,
    pub power_tolerance_w: f64,
    pub energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum NonlinearWindingError {
    Invalid(&'static str),
    Winding(WindingError),
    Conductor(OhmicWindingError),
    Field(NonlinearEddyError),
    Circuit(CircuitError),
    Newton(NewtonError),
    TrialBudget,
    Interrupted(SolverRunProgress),
    NotConverged(NewtonReport),
    Balance(&'static str),
    Cancelled,
}
impl core::fmt::Display for NonlinearWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear winding circuit: {self:?}")
    }
}
impl std::error::Error for NonlinearWindingError {}
impl From<WindingError> for NonlinearWindingError {
    fn from(e: WindingError) -> Self {
        match e {
            WindingError::Magnetic(MagneticError::Cancelled) => Self::Cancelled,
            other => Self::Winding(other),
        }
    }
}
impl From<OhmicWindingError> for NonlinearWindingError {
    fn from(e: OhmicWindingError) -> Self {
        match e {
            OhmicWindingError::Cancelled => Self::Cancelled,
            other => Self::Conductor(other),
        }
    }
}
impl From<NonlinearEddyError> for NonlinearWindingError {
    fn from(e: NonlinearEddyError) -> Self {
        match e {
            NonlinearEddyError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<CircuitError> for NonlinearWindingError {
    fn from(e: CircuitError) -> Self {
        Self::Circuit(e)
    }
}
impl From<NewtonError> for NonlinearWindingError {
    fn from(e: NewtonError) -> Self {
        Self::Newton(e)
    }
}
type Result<T> = core::result::Result<T, NonlinearWindingError>;

struct Layout {
    graph: CircuitGraph,
    dae: CircuitDae,
    windings: NonlinearWindings,
    external_ports: usize,
}

/// Accepted original field plus the original circuit's remaining fluxes,
/// capacitor charges and algebraic solution. Immutable topology is shared.
#[derive(Clone)]
pub struct NonlinearWindingCircuit<'m> {
    layout: Rc<Layout>,
    field: PlanarNonlinearEddy<'m>,
    conductors: Vec<WindingConductor>,
    conductor_windings: Vec<usize>,
    resistance_ohm: Vec<f64>,
    state: Vec<f64>,
    currents_a: Vec<f64>,
    flux_linkages_wb: Vec<f64>,
}

#[derive(Clone, Debug)]
pub struct NonlinearWindingStep {
    /// Native circuit ledger includes PRIVATE interface ports. Their work is
    /// internal; use supplied_work_j for the combined accounting boundary.
    pub native: StepRecord,
    pub field: NonlinearEddyStep,
    /// Spatial stranded heat at the backward-Euler endpoint current.
    pub cells: Vec<WindingCellHeat>,
    pub graph_resistors: Vec<WindingGraphResistorHeat>,
    pub currents_a: Vec<f64>,
    /// Native step-average voltage, in original first-to-second orientation.
    pub terminal_voltage_v: Vec<f64>,
    pub flux_linkages_wb: [Vec<f64>; 2],
    pub terminal_residual_v: Vec<f64>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Core eddy + stranded conductor + native graph physical heat, once each.
    pub heat_j: f64,
    /// Original external source work ONLY; private port work cancels internally.
    pub supplied_work_j: f64,
    /// Magnetic and declared end-inductance BE loss, never physical heating.
    pub numerical_dissipation_j: f64,
    pub energy_residual_j: f64,
    pub field_trials: usize,
    pub report: NewtonReport,
}
pub struct NonlinearWindingCandidate<'m> {
    circuit: NonlinearWindingCircuit<'m>,
    step: NonlinearWindingStep,
}
impl<'m> NonlinearWindingCandidate<'m> {
    pub fn circuit(&self) -> &NonlinearWindingCircuit<'m> {
        &self.circuit
    }
    pub fn step(&self) -> &NonlinearWindingStep {
        &self.step
    }
    pub fn into_parts(self) -> (NonlinearWindingCircuit<'m>, NonlinearWindingStep) {
        (self.circuit, self.step)
    }
}

fn finite(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(NonlinearWindingError::Invalid(
            "unrepresentable winding arithmetic",
        ))
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(NonlinearWindingError::Cancelled)
    } else {
        Ok(())
    }
}
fn numerical(e: FiniteDifferenceError<NonlinearWindingError>) -> NonlinearWindingError {
    match e {
        FiniteDifferenceError::Source(e) => e,
        FiniteDifferenceError::Invalid(s) => NonlinearWindingError::Invalid(s),
    }
}

fn impressed(
    profiles: &[Vec<f64>],
    currents: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    if currents.len() != profiles.len() || currents.iter().any(|v| !v.is_finite()) {
        return Err(NonlinearWindingError::Invalid(
            "complete finite terminal currents required",
        ));
    }
    let mut j = vec![0.0; profiles[0].len()];
    for (profile, &current) in profiles.iter().zip(currents) {
        for (value, &density) in j.iter_mut().zip(profile) {
            poll(cancelled)?;
            *value = finite(density.mul_add(current, *value))?;
        }
    }
    Ok(j)
}

fn linkages(
    field: &PlanarNonlinearEddy<'_>,
    windings: &NonlinearWindings,
    currents: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    windings
        .profiles_per_m2
        .iter()
        .enumerate()
        .map(|(i, profile)| {
            let active = pairing_on_mesh(
                field.problem().mesh(),
                field.potential_wb_m(),
                profile,
                cancelled,
            )?;
            finite(windings.end_inductance_h[i].mul_add(currents[i], active))
        })
        .collect()
}

fn bind(
    field: &PlanarNonlinearEddy<'_>,
    profiles: &[Vec<f64>],
    sources: Vec<WindingConductor>,
    previous: Option<&[WindingConductor]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<WindingConductor>, Vec<usize>, Vec<f64>)> {
    // A stranded current has already paid rho*J^2 in the terminal model. It
    // cannot simultaneously be a resolved floating eddy conductor.
    for source in &sources {
        poll(cancelled)?;
        if !field.insulating_triangles().contains(&source.triangle) {
            return Err(NonlinearWindingError::Invalid(
                "stranded winding and resolved eddy conductor cells must be disjoint",
            ));
        }
    }
    Ok(ohmic::bind_sources(
        field.problem().mesh(),
        profiles,
        sources,
        previous,
        |source| {
            let material = source.conductor.material();
            match &field.problem().materials()[source.triangle] {
                IsotropicMagneticMaterial::Linear(original) => {
                    admit_fixed_magnetic_refresh(original, material)?
                }
                IsotropicMagneticMaterial::Reversible(law) => {
                    law.admit_fixed_material_refresh(material)?
                }
            }
            Ok(())
        },
        cancelled,
    )?)
}

fn ports(layout: &Layout, external: &[f64], currents: &[f64]) -> Result<Vec<f64>> {
    if external.len() != layout.external_ports || external.iter().any(|v| !v.is_finite()) {
        return Err(NonlinearWindingError::Invalid(
            "complete finite original source inputs required",
        ));
    }
    let mut input = external.to_vec();
    input.extend(currents.iter().map(|i| -*i));
    Ok(input)
}

fn native_seed(dae: &CircuitDae, accepted: &[f64], input: &[f64]) -> Result<Vec<f64>> {
    let nd = dae.system.differential_dim();
    let mut seed = vec![0.0; dae.system.state_dim()];
    seed[..nd].copy_from_slice(&accepted[..nd]);
    let initialized = dae.consistent_initial_state(&seed, input)?;
    if initialized[..nd]
        .iter()
        .zip(accepted)
        .any(|(a, b)| !a.is_finite() || (a - b).abs() > 64.0 * f64::EPSILON * b.abs())
    {
        return Err(NonlinearWindingError::Invalid(
            "native constraints changed conserved circuit flux or charge",
        ));
    }
    Ok(initialized)
}

fn end_energy(ends: &[f64], currents: &[f64]) -> Result<f64> {
    ends.iter().zip(currents).try_fold(0.0, |total, (&l, &i)| {
        finite((0.5 * l * i).mul_add(i, total))
    })
}

impl<'m> NonlinearWindingCircuit<'m> {
    /// Bind actual accepted field history to original graph terminals. Bound
    /// inductor values are placeholders replaced by the field; unbound L and C
    /// retain their native values. Initial J must equal the authored profiles
    /// times the supplied physical currents; no field reset is hidden here.
    pub fn new(
        field: PlanarNonlinearEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        sources: Vec<WindingConductor>,
        initial: NonlinearWindingInitial,
        config: NonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        use NonlinearWindingError::Invalid;
        poll(cancelled)?;
        let count = windings.profiles_per_m2.len();
        let values = count.checked_mul(field.problem().mesh().faces().len());
        if count == 0
            || count > config.max_windings
            || config.max_windings > MAX_WINDINGS
            || values.is_none_or(|n| n > config.max_profile_values)
            || windings.branch_indices.len() != count
            || windings.end_inductance_h.len() != count
            || windings
                .end_inductance_h
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
            || !config.current_balance_tolerance_turns.is_finite()
            || config.current_balance_tolerance_turns < 0.0
            || graph.node_count > 32
            || graph.branches.len() > 64
            || graph.transformers.len() > 16
            || graph.branches.iter().any(|&(a, b, _)| a == b)
        {
            return Err(Invalid(
                "invalid winding geometry, end model, graph capacity or limits",
            ));
        }
        let mesh = field.problem().mesh();
        let thickness = match mesh.metric() {
            fs_rep_mesh::Metric2::Planar { thickness } => thickness,
            _ => {
                return Err(Invalid(
                    "nonlinear windings require a planar magnetic section",
                ));
            }
        };
        admit_profiles(
            mesh,
            &windings.profiles_per_m2,
            thickness,
            config.current_balance_tolerance_turns,
            cancelled,
        )?;
        let expected_j = impressed(&windings.profiles_per_m2, &initial.currents_a, cancelled)?;
        if expected_j
            .iter()
            .zip(field.impressed_current_density_a_m2())
            .any(|(a, b)| (a - b).abs() > 128.0 * f64::EPSILON * a.abs().max(b.abs()))
        {
            return Err(Invalid(
                "initial winding currents disagree with accepted impressed field current",
            ));
        }
        let mut external_ports = 0;
        for &(_, _, branch) in &graph.branches {
            if let Branch::VoltageSource { port } | Branch::CurrentSource { port } = branch {
                if port >= 32 {
                    return Err(Invalid("external circuit port exceeds capacity"));
                }
                external_ports = external_ports.max(port + 1);
            }
        }
        if external_ports + count > 32 {
            return Err(Invalid(
                "combined external and winding port capacity exceeded",
            ));
        }
        let mut lowered = graph.clone();
        let mut bound = vec![false; graph.branches.len()];
        for (winding, &branch) in windings.branch_indices.iter().enumerate() {
            poll(cancelled)?;
            if branch >= graph.branches.len()
                || bound[branch]
                || !matches!(graph.branches[branch].2, Branch::Inductor { .. })
            {
                return Err(Invalid(
                    "each winding needs a distinct original inductor branch",
                ));
            }
            bound[branch] = true;
            lowered.branches[branch].2 = Branch::CurrentSource {
                port: external_ports + winding,
            };
        }
        let dae = assemble_circuit(&lowered)?;
        if dae.system.state_dim() > 128 {
            return Err(Invalid("extended circuit state exceeds capacity"));
        }
        if initial.remaining_fluxes_wb.len() != dae.flux_index.len()
            || initial.charges_c.len() != dae.charge_index.len()
            || initial
                .remaining_fluxes_wb
                .iter()
                .chain(&initial.charges_c)
                .any(|v| !v.is_finite())
        {
            return Err(Invalid(
                "complete finite remaining inductor fluxes and capacitor charges required",
            ));
        }
        let (conductors, conductor_windings, resistance_ohm) =
            bind(&field, &windings.profiles_per_m2, sources, None, cancelled)?;
        let flux_linkages_wb = linkages(&field, &windings, &initial.currents_a, cancelled)?;
        let layout = Rc::new(Layout {
            graph,
            dae,
            windings,
            external_ports,
        });
        let input = ports(&layout, &initial.external_inputs, &initial.currents_a)?;
        let mut seed = vec![0.0; layout.dae.system.state_dim()];
        for (&index, &value) in layout
            .dae
            .flux_index
            .iter()
            .zip(&initial.remaining_fluxes_wb)
        {
            seed[index] = value;
        }
        for (&index, &value) in layout.dae.charge_index.iter().zip(&initial.charges_c) {
            seed[index] = value;
        }
        let state = native_seed(&layout.dae, &seed, &input)?;
        let total = finite(
            field.energy_j()
                + layout.dae.system.hamiltonian(&state)
                + end_energy(&layout.windings.end_inductance_h, &initial.currents_a)?,
        )?;
        if total < 0.0 {
            return Err(Invalid("negative initial field/circuit storage"));
        }
        poll(cancelled)?;
        Ok(Self {
            layout,
            field,
            conductors,
            conductor_windings,
            resistance_ohm,
            state,
            currents_a: initial.currents_a,
            flux_linkages_wb,
        })
    }

    pub fn field(&self) -> &PlanarNonlinearEddy<'m> {
        &self.field
    }
    pub fn graph(&self) -> &CircuitGraph {
        &self.layout.graph
    }
    pub fn dae(&self) -> &CircuitDae {
        &self.layout.dae
    }
    pub fn branch_indices(&self) -> &[usize] {
        &self.layout.windings.branch_indices
    }
    pub fn profiles_per_m2(&self) -> &[Vec<f64>] {
        &self.layout.windings.profiles_per_m2
    }
    pub fn end_inductance_h(&self) -> &[f64] {
        &self.layout.windings.end_inductance_h
    }
    pub fn conductors(&self) -> &[WindingConductor] {
        &self.conductors
    }
    pub fn conductor_windings(&self) -> &[usize] {
        &self.conductor_windings
    }
    pub fn resistance_ohm(&self) -> &[f64] {
        &self.resistance_ohm
    }
    pub fn state(&self) -> &[f64] {
        &self.state
    }
    pub fn currents_a(&self) -> &[f64] {
        &self.currents_a
    }
    pub fn flux_linkages_wb(&self) -> &[f64] {
        &self.flux_linkages_wb
    }
    pub fn time_s(&self) -> f64 {
        self.field.time_s()
    }

    /// Simultaneously solve native graph voltage/current and nonlinear field.
    /// Endpoint winding currents are held private circuit inputs over this BE
    /// step; graph storage uses its native discrete gradient. The combined
    /// scheme is first order. Final physical residual and energy checks are
    /// independent of the interface solver's numerical convergence flag.
    pub fn step_to_candidate(
        &self,
        target_s: f64,
        core_conductors: Vec<EddyConductor>,
        winding_conductors: Vec<WindingConductor>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearWindingCandidate<'m>> {
        use NonlinearWindingError::{Balance, Invalid};
        poll(cancelled)?;
        controls.newton.validate()?;
        let dt = target_s - self.time_s();
        if !target_s.is_finite()
            || !positive(dt)
            || controls.maximum_field_trials == 0
            || controls.maximum_newton_iterations == 0
            || !positive(controls.difference_step_a)
            || !positive(controls.terminal_voltage_tolerance_v)
            || !positive(controls.power_tolerance_w)
            || !positive(controls.energy_tolerance_j)
        {
            return Err(Invalid(
                "positive interval, work budgets and physical tolerances required",
            ));
        }
        ports(&self.layout, external_inputs, &self.currents_a)?;
        let (conductors, conductor_windings, resistance_ohm) = bind(
            &self.field,
            self.profiles_per_m2(),
            winding_conductors,
            Some(&self.conductors),
            cancelled,
        )?;
        let cancel = RefCell::new(cancelled);
        let trials = Cell::new(0);
        let evaluate = |currents: &[f64]| -> Result<Trial<'m>> {
            let mut check = || (*cancel.borrow_mut())();
            poll(&mut check)?;
            if trials.get() >= controls.maximum_field_trials {
                return Err(NonlinearWindingError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            let j = impressed(self.profiles_per_m2(), currents, &mut check)?;
            let input = ports(&self.layout, external_inputs, currents)?;
            let x0 = native_seed(&self.layout.dae, &self.state, &input)?;
            poll(&mut check)?;
            let (native, _) = self.layout.dae.step_audited(&x0, &input, dt)?;
            poll(&mut check)?;
            let field = self.field.step_to_candidate(
                target_s,
                core_conductors.clone(),
                &j,
                controls.field,
                &mut check,
            )?;
            let linkage = linkages(field.state(), &self.layout.windings, currents, &mut check)?;
            let voltage = native.y[self.layout.external_ports..].to_vec();
            let residual = currents
                .iter()
                .enumerate()
                .map(|(i, &current)| {
                    finite(
                        voltage[i]
                            - (linkage[i] - self.flux_linkages_wb[i]) / dt
                            - resistance_ohm[i] * current,
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(Trial {
                native,
                field,
                x0,
                linkage,
                voltage,
                residual,
            })
        };
        let problem = FallibleFiniteDifference::new(
            self.currents_a.len(),
            controls.difference_step_a,
            |x: &[f64]| evaluate(x).map(|trial| trial.residual),
        )
        .map_err(numerical)?;
        let construction =
            NewtonKrylovState::new(&problem, self.currents_a.clone(), controls.newton);
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(e) = problem.take_failure() {
            return Err(numerical(e));
        }
        poll(&mut || (*cancel.borrow_mut())())?;
        if result.progress == SolverRunProgress::Paused {
            return Err(NonlinearWindingError::Cancelled);
        }
        if result.progress != SolverRunProgress::Complete {
            return Err(NonlinearWindingError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(NonlinearWindingError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        if trial
            .residual
            .iter()
            .any(|v| v.abs() > controls.terminal_voltage_tolerance_v)
        {
            return Err(Balance(
                "native terminal voltage and actual field linkage do not close",
            ));
        }
        let mut check = || (*cancel.borrow_mut())();
        let dae = &self.layout.dae;
        let midpoint: Vec<_> = trial
            .x0
            .iter()
            .zip(&trial.native.x)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let (cells, winding_heat) = ohmic::spatial_heat(
            self.field.problem().mesh(),
            self.profiles_per_m2(),
            &conductors,
            &conductor_windings,
            &resistance_ohm,
            &solver.x,
            dt,
            controls.power_tolerance_w,
            &mut check,
        )?;
        let (graph_resistors, graph_heat) = ohmic::graph_heat(
            self.graph(),
            dae,
            &midpoint,
            dt,
            controls.power_tolerance_w,
            &mut check,
        )?;
        let external_work = finite(
            dt * external_inputs
                .iter()
                .zip(&trial.native.y)
                .map(|(u, y)| u * y)
                .sum::<f64>(),
        )?;
        let terminal_work = finite(
            dt * solver
                .x
                .iter()
                .zip(&trial.voltage)
                .map(|(i, v)| i * v)
                .sum::<f64>(),
        )?;
        let end_before = end_energy(self.end_inductance_h(), &self.currents_a)?;
        let end_after = end_energy(self.end_inductance_h(), &solver.x)?;
        let increments: Vec<_> = solver
            .x
            .iter()
            .zip(&self.currents_a)
            .map(|(a, b)| a - b)
            .collect();
        let end_numerical = end_energy(self.end_inductance_h(), &increments)?;
        let end_work = finite(
            self.end_inductance_h()
                .iter()
                .zip(&increments)
                .zip(&solver.x)
                .map(|((&l, &di), &i)| l * di * i)
                .sum(),
        )?;
        let field_work = finite(
            solver
                .x
                .iter()
                .enumerate()
                .map(|(i, &current)| {
                    current
                        * (trial.linkage[i]
                            - self.flux_linkages_wb[i]
                            - self.end_inductance_h()[i] * increments[i])
                })
                .sum(),
        )?;
        let field_step = trial.field.step();
        let circuit_before = finite(dae.system.hamiltonian(&trial.x0))?;
        let circuit_after = finite(dae.system.hamiltonian(&trial.native.x))?;
        let stored_before_j = finite(circuit_before + field_step.stored_before_j + end_before)?;
        let stored_after_j = finite(circuit_after + field_step.stored_after_j + end_after)?;
        let heat_j = finite(winding_heat + graph_heat + field_step.heat_j)?;
        let numerical_dissipation_j = finite(field_step.numerical_dissipation_j + end_numerical)?;
        let energy_residual_j = finite(
            stored_after_j - stored_before_j + heat_j + numerical_dissipation_j - external_work,
        )?;
        let defects = [
            trial.native.supplied - (external_work - terminal_work),
            trial.native.dissipated - graph_heat,
            trial.native.delta_h - (circuit_after - circuit_before),
            trial.native.balance_residual(),
            field_work - field_step.supplied_work_j,
            end_after - end_before + end_numerical - end_work,
            terminal_work - field_work - end_work - winding_heat,
            energy_residual_j,
        ];
        if circuit_before < 0.0
            || circuit_after < 0.0
            || heat_j < 0.0
            || numerical_dissipation_j < 0.0
            || trial
                .native
                .x
                .iter()
                .chain(&trial.native.y)
                .any(|v| !v.is_finite())
            || defects
                .iter()
                .any(|v| !v.is_finite() || v.abs() > controls.energy_tolerance_j)
        {
            return Err(Balance(
                "external work, private terminal transfer, field/circuit storage and physical/numerical losses do not close",
            ));
        }
        let (field, field_step) = trial.field.into_parts();
        let step = NonlinearWindingStep {
            native: trial.native.clone(),
            field: field_step,
            cells,
            graph_resistors,
            currents_a: solver.x.clone(),
            terminal_voltage_v: trial.voltage,
            flux_linkages_wb: [self.flux_linkages_wb.clone(), trial.linkage.clone()],
            terminal_residual_v: trial.residual,
            stored_before_j,
            stored_after_j,
            heat_j,
            supplied_work_j: external_work,
            numerical_dissipation_j,
            energy_residual_j,
            field_trials: trials.get(),
            report: result.report,
        };
        let circuit = Self {
            layout: self.layout.clone(),
            field,
            conductors,
            conductor_windings,
            resistance_ohm,
            state: trial.native.x,
            currents_a: solver.x,
            flux_linkages_wb: trial.linkage,
        };
        poll(&mut check)?;
        Ok(NonlinearWindingCandidate { circuit, step })
    }
}

struct Trial<'m> {
    native: StepRecord,
    field: PlanarNonlinearEddyCandidate<'m>,
    x0: Vec<f64>,
    linkage: Vec<f64>,
    voltage: Vec<f64>,
    residual: Vec<f64>,
}
