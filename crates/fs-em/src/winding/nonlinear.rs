//! Native circuit terminals coupled to the actual nonlinear magnetic field.
//!
//! Each winding owns a signed spatial J/current profile and material
//! resistance. Explicit solid-region bindings normalize that profile from the
//! current conductivity and retain resolved diffusion/skin/proximity history.
//! Its original inductor branch becomes a PRIVATE power port in
//! the original Kirchhoff DAE. By default its input is -I (current leaving the
//! first terminal) and its native conjugate output is the winding voltage V.
//! An explicitly voltage-input port instead imposes the actual field voltage
//! and reads the native output -I. Both use the same shared Newton owner and
//! V = F_new^T delta(A)/dt + L_end delta(I)/dt + R_dc*I + V_motional.
//! The moving solid owner derives V_motional from the actual Lorentz field;
//! it is zero for stationary windings. Fixed profiles permit
//! the equivalent delta(lambda)/dt form. A conductivity refresh changes the
//! endpoint linkage diagnostic but cannot invent a delta(F)^T A voltage.
//! Voltage-input ports admit current-driven
//! winding topologies without introducing artificial resistance or storage.
//!
//! Every trial restarts the accepted circuit and nonlinear-eddy histories.
//! Native graph work -dt*I*V cancels terminal work into the field/winding. Core
//! and solid-conductor Joule heat share one actual total-field receipt. Stranded
//! Joule heat, graph heat and backward-Euler numerical loss remain separate.
//! Strand cells cannot also own resolved eddy conductivity.
//! This is a fixed-geometry reversible isotropic MQS lane. A pinned linear mu(T)
//! refresh carries explicit magnetic parameter work; a thermal owner must close
//! it with the field entropy. Nonlinear H(B) stays fixed. No core hysteresis or
//! fitted inductance is used. The moving-conductor step keeps the stranded
//! winding geometry fixed and resolves its motional feedback through this same
//! terminal Newton solve, with signed Lorentz work in the combined account.
//! Caller-declared shared end-flux modes add passive self/mutual storage outside
//! the resolved field, including the solid and moving-conductor histories.
//! The allocating native DAE must admit the declared terminal causality and
//! initial current/voltage receipt; singular ideal-source topologies can refuse.

/// Moving bulk conductors coupled to stationary, spatially resolved windings.
pub mod moving;
/// Source-driven ordinary graph resistance on this same native circuit owner.
pub mod resistors;
/// Source-resolved solid winding diffusion, skin/proximity current and heat.
pub mod solid;
use resistors::GraphResistorMaterial;
use solid::SolidWindingRegion;

use super::end::{EndWindingError, EndWindingInductance, EndWindingMode, EndWindingStep};
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
            PlanarNonlinearEddyCandidate, admit_material_refresh,
        },
    },
    nonlinear::IsotropicMagneticMaterial,
};
use fs_feec::weighted::WeightedAssemblyLimits;
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
    /// Independent unrepresented end/leakage self-inductance [H], including zero.
    /// Shared end-flux modes are supplied separately to new_with_end_modes.
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

/// Causality of a private native port, in winding order. This selects which
/// terminal variable the circuit receives; it does not change the field law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NonlinearWindingPort {
    /// Impose -I and read V. The default for voltage-driven windings; the native
    /// graph must admit independently perturbed winding-current trials.
    CurrentInput,
    /// Impose V = delta(lambda)/dt + R*I and read -I. This supports windings
    /// driven by ideal current sources and mixed current/voltage-driven graphs.
    /// The native graph must admit independently imposed terminal voltages.
    VoltageInput {
        /// Initial first-to-second voltage [V], used only at construction.
        /// It must admit the original native storage and supplied field current.
        /// Field linkage at one instant alone does not determine this voltage.
        initial_voltage_v: f64,
        /// Positive numerical scale [ohm] converting the native-minus-field
        /// current residual to volts for the shared interface Newton solver.
        /// This adds no physical impedance, heat, storage or terminal voltage.
        residual_scale_ohm: f64,
        /// Independent physical current allowance [A], at initialization and
        /// after each finite step; not a scaled Newton convergence threshold.
        current_tolerance_a: f64,
    },
}

impl NonlinearWindingPort {
    fn branch(self, port: usize) -> Branch {
        match self {
            Self::CurrentInput => Branch::CurrentSource { port },
            Self::VoltageInput { .. } => Branch::VoltageSource { port },
        }
    }
}

/// Explicit physical acceptance and nested work budgets.
#[derive(Clone, Copy, Debug)]
pub struct NonlinearWindingControls<F = NonlinearEddyControls> {
    pub field: F,
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
    End(EndWindingError),
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
impl From<EndWindingError> for NonlinearWindingError {
    fn from(error: EndWindingError) -> Self {
        Self::End(error)
    }
}
pub(super) type Result<T> = core::result::Result<T, NonlinearWindingError>;

pub(super) struct Layout {
    pub(super) graph: CircuitGraph,
    pub(super) dae: CircuitDae,
    pub(super) windings: NonlinearWindings,
    pub(super) terminal_ports: Vec<NonlinearWindingPort>,
    pub(super) external_ports: usize,
}

/// Accepted original field plus the original circuit's remaining fluxes,
/// capacitor charges and algebraic solution. Immutable topology is shared.
#[derive(Clone)]
pub struct NonlinearWindingCircuit<'m> {
    layout: Rc<Layout>,
    end: Rc<EndWindingInductance>,
    graph_materials: Vec<GraphResistorMaterial>,
    field: PlanarNonlinearEddy<'m>,
    conductors: Vec<WindingConductor>,
    conductor_windings: Vec<usize>,
    resistance_ohm: Vec<f64>,
    state: Vec<f64>,
    currents_a: Vec<f64>,
    flux_linkages_wb: Vec<f64>,
    solid_regions: Rc<[SolidWindingRegion]>,
    solid_profiles_per_m2: Option<Rc<Vec<Vec<f64>>>>,
    profile_balance_tolerance_turns: f64,
}

#[derive(Clone, Debug)]
pub struct NonlinearWindingStep<F = NonlinearEddyStep> {
    /// Native circuit ledger includes PRIVATE interface ports. Their work is
    /// internal; use supplied_work_j for the combined accounting boundary.
    pub native: StepRecord,
    pub field: F,
    /// Spatial stranded heat at the backward-Euler endpoint current.
    pub cells: Vec<WindingCellHeat>,
    pub graph_resistors: Vec<WindingGraphResistorHeat>,
    /// Complete source bindings used for ordinary graph loss, in branch order.
    /// Empty for the original explicitly authored numeric-resistor path.
    pub graph_materials: Vec<GraphResistorMaterial>,
    pub currents_a: Vec<f64>,
    /// Native step-average voltage, in original first-to-second orientation.
    pub terminal_voltage_v: Vec<f64>,
    /// Endpoint linkage diagnostics. Solid conductivity refresh changes their
    /// weighting; their difference is NOT the induced voltage in that case.
    pub flux_linkages_wb: [Vec<f64>; 2],
    /// Actual F_new^T delta(A)/dt plus explicit end-inductance voltage [V].
    /// Changing conductivity weights in a stationary field creates no EMF.
    pub inductive_voltage_v: Vec<f64>,
    /// Solid-winding back-EMF from -integral(J_unit * (v cross B)_z dV),
    /// including zero-current operation. Zero for the stationary step.
    pub motional_voltage_v: Vec<f64>,
    /// Actual native voltage minus the field/copper-required voltage [V].
    pub terminal_residual_v: Vec<f64>,
    /// Actual native current minus field current [A]. Exactly zero for
    /// current-input ports; voltage-input ports use their declared allowance.
    pub terminal_current_residual_a: Vec<f64>,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Native graph/field change plus the directly evaluated end-mode change.
    /// End energy is never differenced across a large magnetic bias.
    pub stored_change_j: f64,
    /// Actual solid/floating field + stranded + native graph heat, once each.
    pub heat_j: f64,
    /// Original external source work ONLY; private port work cancels internally.
    pub supplied_work_j: f64,
    /// Magnetic and declared end-inductance BE loss, never physical heating.
    pub numerical_dissipation_j: f64,
    /// Signed source-resolved magnetic work, separate from external source
    /// ports and Joule heat. A caloric/material owner must close this account.
    pub parameter_work_j: f64,
    /// Signed work done by the field on moving bulk conductors. Zero for the
    /// stationary operation; negative values are mechanical input to the field.
    pub mechanical_work_j: f64,
    pub energy_residual_j: f64,
    pub field_trials: usize,
    pub report: NewtonReport,
}
pub struct NonlinearWindingCandidate<'m, F = NonlinearEddyStep> {
    circuit: NonlinearWindingCircuit<'m>,
    step: NonlinearWindingStep<F>,
}
impl<'m, F> NonlinearWindingCandidate<'m, F> {
    pub fn circuit(&self) -> &NonlinearWindingCircuit<'m> {
        &self.circuit
    }
    pub fn step(&self) -> &NonlinearWindingStep<F> {
        &self.step
    }
    pub fn into_parts(self) -> (NonlinearWindingCircuit<'m>, NonlinearWindingStep<F>) {
        (self.circuit, self.step)
    }
}

pub(super) fn finite(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(NonlinearWindingError::Invalid(
            "unrepresentable winding arithmetic",
        ))
    }
}
pub(super) fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
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

pub(super) fn impressed(
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

pub(super) trait PlanarInductionHistory {
    fn mesh(&self) -> &fs_rep_mesh::TriComplex2;
    fn potential_wb_m(&self) -> &[f64];
    fn time_s(&self) -> f64;
}

impl PlanarInductionHistory for PlanarNonlinearEddy<'_> {
    fn mesh(&self) -> &fs_rep_mesh::TriComplex2 {
        self.problem().mesh()
    }
    fn potential_wb_m(&self) -> &[f64] {
        self.potential_wb_m()
    }
    fn time_s(&self) -> f64 {
        self.time_s()
    }
}

pub(super) fn linkages(
    field: &impl PlanarInductionHistory,
    profiles: &[Vec<f64>],
    end: &EndWindingInductance,
    currents: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>> {
    let end_linkage = end.flux_linkages_wb(currents)?;
    profiles
        .iter()
        .enumerate()
        .map(|(i, profile)| {
            let active = pairing_on_mesh(field.mesh(), field.potential_wb_m(), profile, cancelled)?;
            finite(end_linkage[i] + active)
        })
        .collect()
}

fn bind(
    field: &PlanarNonlinearEddy<'_>,
    materials: Option<&[IsotropicMagneticMaterial]>,
    profiles: &[Vec<f64>],
    sources: Vec<WindingConductor>,
    previous: Option<&[WindingConductor]>,
    solid: bool,
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
    if !solid
        && field
            .conductors()
            .iter()
            .any(|s| profiles.iter().any(|p| p[s.triangle] != 0.0))
    {
        return Err(NonlinearWindingError::Invalid(
            "active resolved conductors require explicit solid winding bindings",
        ));
    }
    let mut retained = Vec::new();
    let mut stranded_profiles = Vec::new();
    if solid {
        let mut resolved = vec![false; field.problem().mesh().faces().len()];
        for source in field.conductors() {
            resolved[source.triangle] = true;
        }
        for (winding, profile) in profiles.iter().enumerate() {
            poll(cancelled)?;
            if profile
                .iter()
                .enumerate()
                .any(|(cell, &j)| j != 0.0 && !resolved[cell])
            {
                retained.push(winding);
                stranded_profiles.push(
                    profile
                        .iter()
                        .enumerate()
                        .map(|(cell, &j)| if resolved[cell] { 0.0 } else { j })
                        .collect(),
                );
            }
        }
    }
    let (sources, mut winding_indices, resistance) = ohmic::bind_sources(
        field.problem().mesh(),
        if solid { &stranded_profiles } else { profiles },
        sources,
        previous,
        |source| {
            let material = source.conductor.material();
            match &materials.unwrap_or_else(|| field.problem().materials())[source.triangle] {
                IsotropicMagneticMaterial::Linear(original) => {
                    admit_fixed_magnetic_refresh(original, material)?
                }
                IsotropicMagneticMaterial::Reversible(law) => {
                    law.admit_fixed_material_refresh(material)?
                }
                IsotropicMagneticMaterial::Recoil(law) => {
                    law.material().admit_fixed_material_refresh(material)?
                }
            }
            Ok(())
        },
        cancelled,
    )?;
    if solid {
        let mut complete = vec![0.0; profiles.len()];
        for (winding, value) in retained.iter().zip(resistance) {
            complete[*winding] = value;
        }
        for winding in &mut winding_indices {
            *winding = retained[*winding];
        }
        Ok((sources, winding_indices, complete))
    } else {
        Ok((sources, winding_indices, resistance))
    }
}

fn admit_external_inputs(layout: &Layout, external: &[f64]) -> Result<()> {
    if external.len() != layout.external_ports || external.iter().any(|v| !v.is_finite()) {
        return Err(NonlinearWindingError::Invalid(
            "complete finite original source inputs required",
        ));
    }
    Ok(())
}

fn ports(layout: &Layout, external: &[f64], currents: &[f64], voltage: &[f64]) -> Result<Vec<f64>> {
    admit_external_inputs(layout, external)?;
    if currents.len() != layout.terminal_ports.len()
        || voltage.len() != layout.terminal_ports.len()
        || currents.iter().chain(voltage).any(|v| !v.is_finite())
    {
        return Err(NonlinearWindingError::Invalid(
            "complete finite winding terminal currents and voltages required",
        ));
    }
    let mut input = external.to_vec();
    input.extend(
        layout
            .terminal_ports
            .iter()
            .enumerate()
            .map(|(i, port)| match port {
                NonlinearWindingPort::CurrentInput => -currents[i],
                NonlinearWindingPort::VoltageInput { .. } => voltage[i],
            }),
    );
    Ok(input)
}

pub(crate) fn native_seed(dae: &CircuitDae, accepted: &[f64], input: &[f64]) -> Result<Vec<f64>> {
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

impl Layout {
    /// Lower the original winding branches into reciprocal native ports while
    /// retaining every unbound storage coordinate and algebraic constraint.
    pub(super) fn new(
        graph: CircuitGraph,
        windings: NonlinearWindings,
        terminal_ports: Vec<NonlinearWindingPort>,
        initial: &NonlinearWindingInitial,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Rc<Self>, Vec<f64>)> {
        use NonlinearWindingError::Invalid;
        poll(cancelled)?;
        let count = windings.profiles_per_m2.len();
        if count == 0
            || count > MAX_WINDINGS
            || windings.branch_indices.len() != count
            || terminal_ports.len() != count
            || initial.currents_a.len() != count
            || initial.currents_a.iter().any(|v| !v.is_finite())
            || terminal_ports.iter().any(|port| {
                matches!(port, NonlinearWindingPort::VoltageInput {
                    initial_voltage_v, residual_scale_ohm, current_tolerance_a
                } if !initial_voltage_v.is_finite()
                    || !positive(*residual_scale_ohm)
                    || !positive(*current_tolerance_a))
            })
            || graph.node_count > 32
            || graph.branches.len() > 64
            || graph.transformers.len() > 16
            || graph.branches.iter().any(|&(a, b, _)| a == b)
        {
            return Err(Invalid(
                "invalid native winding graph, ports or initial currents",
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
            lowered.branches[branch].2 = terminal_ports[winding].branch(external_ports + winding);
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
        let layout = Rc::new(Self {
            graph,
            dae,
            windings,
            terminal_ports,
            external_ports,
        });
        let initial_voltage: Vec<_> = layout
            .terminal_ports
            .iter()
            .map(|port| match port {
                NonlinearWindingPort::CurrentInput => 0.0,
                NonlinearWindingPort::VoltageInput {
                    initial_voltage_v, ..
                } => *initial_voltage_v,
            })
            .collect();
        let input = ports(
            &layout,
            &initial.external_inputs,
            &initial.currents_a,
            &initial_voltage,
        )?;
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
        // The initializer's admitted effort is midpoint(0, state). Read its
        // original receipt without another solve choosing a different split.
        let (_, _, g) = layout.dae.system.structure();
        for (winding, port) in layout.terminal_ports.iter().enumerate() {
            if let NonlinearWindingPort::VoltageInput {
                current_tolerance_a,
                ..
            } = port
            {
                let p = layout.external_ports + winding;
                let output = (layout.dae.system.differential_dim()..state.len())
                    .try_fold(0.0, |sum, row| {
                        finite(g[row * layout.dae.ports + p].mul_add(0.5 * state[row], sum))
                    })?;
                if finite(-output - initial.currents_a[winding])?.abs() > *current_tolerance_a {
                    return Err(Invalid(
                        "native initial winding current disagrees with accepted field history",
                    ));
                }
            }
        }
        poll(cancelled)?;
        Ok((layout, state))
    }
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
        poll(cancelled)?;
        let count = windings.profiles_per_m2.len();
        if count == 0 || count > config.max_windings || config.max_windings > MAX_WINDINGS {
            return Err(NonlinearWindingError::Invalid(
                "invalid winding count or winding capacity",
            ));
        }
        let terminal_ports = vec![NonlinearWindingPort::CurrentInput; count];
        Self::new_with_ports(
            field,
            graph,
            windings,
            terminal_ports,
            sources,
            initial,
            config,
            cancelled,
        )
    }

    /// Bind explicit per-winding native port causality to the same physical
    /// field and circuit histories. Voltage-input ports support ideal current
    /// drive; current-input ports retain the original voltage-driven operation.
    ///
    /// Each voltage-input port declares an initial voltage. The existing native
    /// initializer must accept that voltage at the supplied original fluxes and
    /// charges, and its admitted current must match the supplied field current.
    /// Ambiguous initial current splits may refuse; no stored coordinate or
    /// field history is projected to conceal incompatible initial conditions.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_ports(
        field: PlanarNonlinearEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        terminal_ports: Vec<NonlinearWindingPort>,
        sources: Vec<WindingConductor>,
        initial: NonlinearWindingInitial,
        config: NonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_with_end_modes(
            field,
            graph,
            windings,
            terminal_ports,
            Vec::new(),
            sources,
            initial,
            config,
            cancelled,
        )
    }

    /// Add explicit shared end-flux modes to the original stranded history.
    /// Each mode adds ell*n*n^T outside the resolved field, with signed winding
    /// participation. Dependent modes require no invented leakage or inverse.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_end_modes(
        field: PlanarNonlinearEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        terminal_ports: Vec<NonlinearWindingPort>,
        end_modes: Vec<EndWindingMode>,
        sources: Vec<WindingConductor>,
        initial: NonlinearWindingInitial,
        config: NonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_complete(
            field,
            graph,
            windings,
            terminal_ports,
            end_modes,
            sources,
            vec![],
            initial,
            config,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn new_complete(
        field: PlanarNonlinearEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        terminal_ports: Vec<NonlinearWindingPort>,
        end_modes: Vec<EndWindingMode>,
        sources: Vec<WindingConductor>,
        solid_regions: Vec<SolidWindingRegion>,
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
            || terminal_ports.len() != count
            || terminal_ports.iter().any(|port| {
                matches!(port, NonlinearWindingPort::VoltageInput {
                    initial_voltage_v, residual_scale_ohm, current_tolerance_a
                } if !initial_voltage_v.is_finite()
                    || !positive(*residual_scale_ohm)
                    || !positive(*current_tolerance_a))
            })
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
        let end = Rc::new(EndWindingInductance::new(
            windings.end_inductance_h.clone(),
            end_modes,
        )?);
        let mesh = field.problem().mesh();
        let thickness = match mesh.metric() {
            fs_rep_mesh::Metric2::Planar { thickness } => thickness,
            _ => {
                return Err(Invalid(
                    "nonlinear windings require a planar magnetic section",
                ));
            }
        };
        let solid = if solid_regions.is_empty() {
            None
        } else {
            Some(solid::prepare(
                &field,
                field.conductors(),
                &windings.profiles_per_m2,
                &solid_regions,
                true,
                cancelled,
            )?)
        };
        let profiles = solid
            .as_ref()
            .map(|s| s.profiles.as_slice())
            .unwrap_or(&windings.profiles_per_m2);
        admit_profiles(
            mesh,
            profiles,
            thickness,
            config.current_balance_tolerance_turns,
            cancelled,
        )?;
        let expected_j = impressed(profiles, &initial.currents_a, cancelled)?;
        if expected_j
            .iter()
            .zip(field.impressed_current_density_a_m2())
            .any(|(a, b)| (a - b).abs() > 128.0 * f64::EPSILON * a.abs().max(b.abs()))
        {
            return Err(Invalid(
                "initial winding currents disagree with accepted impressed field current",
            ));
        }
        let (conductors, conductor_windings, mut resistance_ohm) = bind(
            &field,
            None,
            profiles,
            sources,
            None,
            !solid_regions.is_empty(),
            cancelled,
        )?;
        if let Some(solid) = &solid {
            for (total, &dc) in resistance_ohm.iter_mut().zip(&solid.layout.resistance_ohm) {
                *total = finite(*total + dc)?;
            }
        }
        if resistance_ohm.iter().any(|&r| !positive(r)) {
            return Err(Invalid(
                "every winding requires represented positive DC resistance",
            ));
        }
        let flux_linkages_wb = linkages(&field, profiles, &end, &initial.currents_a, cancelled)?;
        let (layout, state) = Layout::new(graph, windings, terminal_ports, &initial, cancelled)?;
        let total = finite(
            field.energy_j()
                + layout.dae.system.hamiltonian(&state)
                + end.stored_energy_j(&initial.currents_a)?,
        )?;
        if total < 0.0 {
            return Err(Invalid("negative initial field/circuit storage"));
        }
        poll(cancelled)?;
        Ok(Self {
            layout,
            end,
            graph_materials: Vec::new(),
            field,
            conductors,
            conductor_windings,
            resistance_ohm,
            state,
            currents_a: initial.currents_a,
            flux_linkages_wb,
            solid_regions: Rc::from(solid_regions),
            solid_profiles_per_m2: solid.map(|s| s.profiles),
            profile_balance_tolerance_turns: config.current_balance_tolerance_turns,
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
        self.solid_profiles_per_m2
            .as_deref()
            .map(Vec::as_slice)
            .unwrap_or(&self.layout.windings.profiles_per_m2)
    }
    /// Independent diagonal contribution, before shared-mode self/mutual terms.
    pub fn end_inductance_h(&self) -> &[f64] {
        &self.layout.windings.end_inductance_h
    }
    /// Original end-flux model, retained through source and material refresh.
    pub fn end_inductance(&self) -> &EndWindingInductance {
        &self.end
    }
    /// Full row-major end self/mutual matrix in original winding order [H].
    pub fn end_inductance_matrix_h(&self) -> &[f64] {
        self.end.matrix_h()
    }
    pub fn terminal_ports(&self) -> &[NonlinearWindingPort] {
        &self.layout.terminal_ports
    }
    pub fn conductors(&self) -> &[WindingConductor] {
        &self.conductors
    }
    pub fn conductor_windings(&self) -> &[usize] {
        &self.conductor_windings
    }
    /// Represented solid DC plus stranded resistance. Dynamic redistribution
    /// and magnetic storage remain in the retained spatial field history.
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
    /// Endpoint field currents are paired with native step-average voltages.
    /// A current-input port holds -I; a voltage-input port holds the actual
    /// field/copper-required V and solves the returned -I. Original external
    /// inputs are held over the finite interval; graph storage uses its native
    /// discrete gradient. This first-order BE coupling reports finite-step work,
    /// not an instantaneous ideal-source impulse or a reset of stored flux.
    /// Physical voltage/current residuals and energy checks are independent of
    /// the interface solver's numerical convergence flag.
    pub fn step_to_candidate(
        &self,
        target_s: f64,
        core_conductors: Vec<EddyConductor>,
        winding_conductors: Vec<WindingConductor>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearWindingCandidate<'m>> {
        self.step_with_materials(
            target_s,
            None,
            core_conductors,
            winding_conductors,
            external_inputs,
            controls,
            cancelled,
        )
    }

    /// Continue the original field/circuit history with source-pinned linear
    /// magnetic temperature refresh. Complete fixed nonlinear H(B) entries are
    /// retained. The native field's signed parameter work is reported separately
    /// from external source work and physical heat. Rejected trials publish no
    /// field, current, linkage or circuit state. Assembly is explicitly bounded
    /// for every field trial; all trials restart from this accepted history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_materials_candidate(
        &self,
        target_s: f64,
        materials: Vec<IsotropicMagneticMaterial>,
        core_conductors: Vec<EddyConductor>,
        winding_conductors: Vec<WindingConductor>,
        external_inputs: &[f64],
        assembly: WeightedAssemblyLimits,
        controls: NonlinearWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearWindingCandidate<'m>> {
        admit_material_refresh(self.field.problem().materials(), &materials, cancelled)?;
        self.step_with_materials(
            target_s,
            Some((materials, assembly)),
            core_conductors,
            winding_conductors,
            external_inputs,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_materials(
        &self,
        target_s: f64,
        materials: Option<(Vec<IsotropicMagneticMaterial>, WeightedAssemblyLimits)>,
        core_conductors: Vec<EddyConductor>,
        winding_conductors: Vec<WindingConductor>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearWindingCandidate<'m>> {
        let solid = if self.solid_regions.is_empty() {
            None
        } else {
            let solid = solid::prepare(
                &self.field,
                &core_conductors,
                self.profiles_per_m2(),
                &self.solid_regions,
                false,
                cancelled,
            )?;
            admit_profiles(
                self.field.problem().mesh(),
                &solid.profiles,
                self.field.problem().thickness_m(),
                self.profile_balance_tolerance_turns,
                cancelled,
            )?;
            Some(solid)
        };
        let evaluate_field = |currents: &[f64], j: &[f64], cancelled: &mut dyn FnMut() -> bool| {
            let mut check = || cancelled();
            let candidate = if let Some(solid) = &solid {
                let (drives, targets) = solid
                    .layout
                    .drives(currents, &mut check)
                    .map_err(NonlinearEddyError::from)?;
                self.field.step_with_region_drives_candidate(
                    target_s,
                    materials.clone(),
                    core_conductors.clone(),
                    j,
                    &drives,
                    &targets,
                    controls.field,
                    &mut check,
                )?
            } else if let Some((materials, assembly)) = &materials {
                self.field.step_with_materials_candidate(
                    target_s,
                    materials.clone(),
                    core_conductors.clone(),
                    j,
                    *assembly,
                    controls.field,
                    &mut check,
                )?
            } else {
                self.field.step_to_candidate(
                    target_s,
                    core_conductors.clone(),
                    j,
                    controls.field,
                    &mut check,
                )?
            };
            Ok(FieldTrial::stationary(candidate))
        };
        self.step_coupled(
            target_s,
            materials
                .as_ref()
                .map(|(materials, _)| materials.as_slice()),
            solid.as_ref().map(|solid| solid.profiles.clone()),
            solid
                .as_ref()
                .map(|solid| solid.layout.resistance_ohm.as_slice()),
            winding_conductors,
            external_inputs,
            controls,
            evaluate_field,
            cancelled,
        )
    }

    /// One terminal solve and one reciprocity/accounting boundary for the
    /// stationary and moving field owners. The complete physical field step
    /// remains in the result; only its scalar transfers share this adapter.
    #[allow(clippy::too_many_arguments)]
    fn step_coupled<F, C>(
        &self,
        target_s: f64,
        materials: Option<&[IsotropicMagneticMaterial]>,
        refreshed_profiles: Option<Rc<Vec<Vec<f64>>>>,
        solid_resistance_ohm: Option<&[f64]>,
        winding_conductors: Vec<WindingConductor>,
        external_inputs: &[f64],
        controls: NonlinearWindingControls<C>,
        evaluate_field: impl Fn(
            &[f64],
            &[f64],
            &mut dyn FnMut() -> bool,
        ) -> Result<FieldTrial<PlanarNonlinearEddy<'m>, F>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NonlinearWindingCandidate<'m, F>> {
        use NonlinearWindingError::Invalid;
        poll(cancelled)?;
        let profiles = refreshed_profiles
            .as_deref()
            .map(Vec::as_slice)
            .unwrap_or_else(|| self.profiles_per_m2());
        let (conductors, conductor_windings, stranded_resistance_ohm) = bind(
            &self.field,
            materials,
            profiles,
            winding_conductors,
            Some(&self.conductors),
            !self.solid_regions.is_empty(),
            cancelled,
        )?;
        let mut resistance_ohm = stranded_resistance_ohm.clone();
        if let Some(solid) = solid_resistance_ohm {
            if solid.len() != resistance_ohm.len() {
                return Err(Invalid("complete solid resistance coordinates required"));
            }
            for (resistance, &dc) in resistance_ohm.iter_mut().zip(solid) {
                if !dc.is_finite() || dc < 0.0 {
                    return Err(Invalid("nonnegative finite solid resistance required"));
                }
                *resistance = finite(*resistance + dc)?;
            }
        }
        if resistance_ohm.iter().any(|&r| !positive(r)) {
            return Err(Invalid(
                "every winding requires represented positive DC resistance",
            ));
        }
        let accepted = TerminalContext {
            field: &self.field,
            layout: &self.layout,
            end: &self.end,
            graph_materials: &self.graph_materials,
            state: &self.state,
            currents_a: &self.currents_a,
            flux_linkages_wb: &self.flux_linkages_wb,
        }
        .step(
            target_s,
            profiles,
            &conductors,
            &conductor_windings,
            &stranded_resistance_ohm,
            &resistance_ohm,
            external_inputs,
            controls,
            evaluate_field,
            cancelled,
        )?;
        let circuit = Self {
            layout: self.layout.clone(),
            end: Rc::clone(&self.end),
            graph_materials: self.graph_materials.clone(),
            field: accepted.field,
            conductors,
            conductor_windings,
            resistance_ohm,
            state: accepted.state,
            currents_a: accepted.currents_a,
            flux_linkages_wb: accepted.flux_linkages_wb,
            solid_regions: self.solid_regions.clone(),
            solid_profiles_per_m2: refreshed_profiles
                .clone()
                .or_else(|| self.solid_profiles_per_m2.clone()),
            profile_balance_tolerance_turns: self.profile_balance_tolerance_turns,
        };
        Ok(NonlinearWindingCandidate {
            circuit,
            step: accepted.step,
        })
    }
}

/// Borrow the accepted physical histories for one shared native terminal solve.
/// The field callback owns its constitutive memory; no trial publishes it.
pub(super) struct TerminalContext<'a, S> {
    pub(super) field: &'a S,
    pub(super) layout: &'a Layout,
    pub(super) end: &'a EndWindingInductance,
    pub(super) graph_materials: &'a [GraphResistorMaterial],
    pub(super) state: &'a [f64],
    pub(super) currents_a: &'a [f64],
    pub(super) flux_linkages_wb: &'a [f64],
}

pub(super) struct TerminalCandidate<S, F> {
    pub(super) field: S,
    pub(super) state: Vec<f64>,
    pub(super) currents_a: Vec<f64>,
    pub(super) flux_linkages_wb: Vec<f64>,
    pub(super) step: NonlinearWindingStep<F>,
}

impl<S: PlanarInductionHistory> TerminalContext<'_, S> {
    fn terminal_ports(&self) -> &[NonlinearWindingPort] {
        &self.layout.terminal_ports
    }
    fn graph(&self) -> &CircuitGraph {
        &self.layout.graph
    }
    fn branch_indices(&self) -> &[usize] {
        &self.layout.windings.branch_indices
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn step<F, C>(
        &self,
        target_s: f64,
        profiles: &[Vec<f64>],
        conductors: &[WindingConductor],
        conductor_windings: &[usize],
        stranded_resistance_ohm: &[f64],
        resistance_ohm: &[f64],
        external_inputs: &[f64],
        controls: NonlinearWindingControls<C>,
        evaluate_field: impl Fn(&[f64], &[f64], &mut dyn FnMut() -> bool) -> Result<FieldTrial<S, F>>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TerminalCandidate<S, F>> {
        use NonlinearWindingError::{Balance, Invalid};
        poll(cancelled)?;
        controls.newton.validate()?;
        let dt = target_s - self.field.time_s();
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
        admit_external_inputs(self.layout, external_inputs)?;
        // The trial source weighting is fixed throughout this terminal solve.
        // Pair its accepted-field reference once, avoiding a second complete
        // mesh/profile traversal in every finite-difference field trial.
        let reference_linkages =
            linkages(self.field, profiles, self.end, self.currents_a, cancelled)?;
        let cancel = RefCell::new(cancelled);
        let trials = Cell::new(0);
        let voltage_input = self
            .terminal_ports()
            .iter()
            .any(|port| matches!(port, NonlinearWindingPort::VoltageInput { .. }));
        let evaluate = |currents: &[f64]| -> Result<Trial<S, F>> {
            let mut check = || (*cancel.borrow_mut())();
            poll(&mut check)?;
            if trials.get() >= controls.maximum_field_trials {
                return Err(NonlinearWindingError::TrialBudget);
            }
            trials.set(trials.get() + 1);
            let j = impressed(profiles, currents, &mut check)?;
            let field = evaluate_field(currents, &j, &mut check)?;
            let motional_voltage_v = field
                .motional_voltage_v
                .clone()
                .unwrap_or_else(|| vec![0.0; currents.len()]);
            if motional_voltage_v.len() != currents.len()
                || motional_voltage_v.iter().any(|v| !v.is_finite())
            {
                return Err(Invalid(
                    "complete finite winding motional voltages required",
                ));
            }
            let potential_increment = field
                .state
                .potential_wb_m()
                .iter()
                .zip(self.field.potential_wb_m())
                .map(|(&new, &old)| finite(new - old))
                .collect::<Result<Vec<_>>>()?;
            let field_linkage_increments_wb = profiles
                .iter()
                .map(|profile| {
                    pairing_on_mesh(self.field.mesh(), &potential_increment, profile, &mut check)
                        .map_err(NonlinearWindingError::from)
                })
                .collect::<Result<Vec<_>>>()?;
            let end = self.end.step(self.currents_a, currents)?;
            let linkage = field_linkage_increments_wb
                .iter()
                .enumerate()
                .map(|(i, &increment)| {
                    finite(reference_linkages[i] + increment + end.flux_increment_wb[i])
                })
                .collect::<Result<Vec<_>>>()?;
            // A conductivity refresh changes the source functional itself.
            // Pair the new functional with delta(A); delta(F^T A) would create
            // voltage in a stationary biased field merely by changing sigma.
            let inductive_voltage_v = field_linkage_increments_wb
                .iter()
                .enumerate()
                .map(|(i, &increment)| finite((increment + end.flux_increment_wb[i]) / dt))
                .collect::<Result<Vec<_>>>()?;
            let required_voltage = inductive_voltage_v
                .iter()
                .zip(&motional_voltage_v)
                .zip(currents)
                .zip(resistance_ohm)
                .map(|(((&inductive, &motional), &current), &r)| {
                    finite(r.mul_add(current, finite(inductive + motional)?))
                })
                .collect::<Result<Vec<_>>>()?;
            let input = ports(self.layout, external_inputs, currents, &required_voltage)?;
            // A voltage trial is a finite-step effort, not an initial voltage
            // clamp. Reinitializing against it could falsely reject a capacitor
            // or project its accepted charge. Keep the entire accepted native
            // state and let the original descriptor step solve its multipliers.
            let x0 = if voltage_input {
                self.state.to_vec()
            } else {
                native_seed(&self.layout.dae, self.state, &input)?
            };
            poll(&mut check)?;
            let (native, _) = self.layout.dae.step_audited(&x0, &input, dt)?;
            poll(&mut check)?;
            let mut voltage = Vec::with_capacity(currents.len());
            let mut voltage_residual = Vec::with_capacity(currents.len());
            let mut current_residual = Vec::with_capacity(currents.len());
            let mut residual = Vec::with_capacity(currents.len());
            for (winding, port) in self.terminal_ports().iter().enumerate() {
                let output = native.y[self.layout.external_ports + winding];
                let v = match port {
                    NonlinearWindingPort::CurrentInput => finite(output)?,
                    NonlinearWindingPort::VoltageInput { .. } => {
                        let (a, b, _) = self.graph().branches[self.branch_indices()[winding]];
                        let potential = |node: usize| {
                            if node == 0 {
                                0.0
                            } else {
                                let index = self.layout.dae.node_potential_index[node - 1];
                                f64::midpoint(x0[index], native.x[index])
                            }
                        };
                        finite(potential(a) - potential(b))?
                    }
                };
                let dv = finite(v - required_voltage[winding])?;
                let di = match port {
                    NonlinearWindingPort::CurrentInput => 0.0,
                    NonlinearWindingPort::VoltageInput { .. } => {
                        finite(-output - currents[winding])?
                    }
                };
                let scaled = match port {
                    NonlinearWindingPort::CurrentInput => dv,
                    NonlinearWindingPort::VoltageInput {
                        residual_scale_ohm, ..
                    } => finite(residual_scale_ohm * di)?,
                };
                voltage.push(v);
                voltage_residual.push(dv);
                current_residual.push(di);
                residual.push(scaled);
            }
            Ok(Trial {
                native,
                field,
                end,
                x0,
                linkage,
                field_linkage_increments_wb,
                inductive_voltage_v,
                motional_voltage_v,
                voltage,
                voltage_residual,
                current_residual,
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
            NewtonKrylovState::new(&problem, self.currents_a.to_vec(), controls.newton);
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
            .voltage_residual
            .iter()
            .any(|v| v.abs() > controls.terminal_voltage_tolerance_v)
        {
            return Err(Balance(
                "native terminal voltage and actual field linkage do not close",
            ));
        }
        if self
            .terminal_ports()
            .iter()
            .zip(&trial.current_residual)
            .any(|(port, di)| {
                matches!(port, NonlinearWindingPort::VoltageInput { current_tolerance_a, .. }
                if di.abs() > *current_tolerance_a)
            })
        {
            return Err(Balance(
                "native terminal current and actual field current do not close",
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
            self.field.mesh(),
            profiles,
            conductors,
            conductor_windings,
            stranded_resistance_ohm,
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
        let end_before = trial.end.stored_before_j;
        let end_after = trial.end.stored_after_j;
        let end_numerical = trial.end.numerical_dissipation_j;
        let end_work = trial.end.endpoint_work_j;
        let end_change = trial.end.stored_change_j;
        let field_work = finite(
            solver
                .x
                .iter()
                .zip(&trial.field_linkage_increments_wb)
                .map(|(&current, &increment)| current * increment)
                .sum(),
        )?;
        let field_step = &trial.field.account;
        let motional_work = finite(
            dt * solver
                .x
                .iter()
                .zip(&trial.motional_voltage_v)
                .map(|(&current, &voltage)| current * voltage)
                .sum::<f64>(),
        )?;
        let circuit_before = finite(dae.system.hamiltonian(&trial.x0))?;
        let circuit_after = finite(dae.system.hamiltonian(&trial.native.x))?;
        let stored_before_j = finite(circuit_before + field_step.stored_before_j + end_before)?;
        let stored_after_j = finite(circuit_after + field_step.stored_after_j + end_after)?;
        let field_change = field_step.stored_change_j;
        let stored_change_j = finite(trial.native.delta_h + field_change + end_change)?;
        let heat_j = finite(winding_heat + graph_heat + field_step.heat_j)?;
        let numerical_dissipation_j = finite(field_step.numerical_dissipation_j + end_numerical)?;
        let parameter_work_j = field_step.parameter_work_j;
        let mechanical_work_j = field_step.mechanical_work_j;
        let energy_residual_j = finite(
            stored_change_j + heat_j + numerical_dissipation_j + mechanical_work_j
                - external_work
                - parameter_work_j,
        )?;
        let defects = [
            trial.native.supplied - (external_work - terminal_work),
            trial.native.dissipated - graph_heat,
            trial.native.delta_h - (circuit_after - circuit_before),
            trial.native.balance_residual(),
            field_work + field_step.uniform_drive_work_j + field_step.motional_drive_work_j
                - field_step.supplied_work_j,
            motional_work - field_step.motional_drive_work_j,
            end_change + end_numerical - end_work,
            terminal_work
                - field_work
                - end_work
                - winding_heat
                - field_step.uniform_drive_work_j
                - motional_work,
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
                "external work, private terminal transfer, field/circuit storage, mechanical work and physical/numerical losses do not close",
            ));
        }
        let FieldTrial {
            state: field,
            step: field_step,
            ..
        } = trial.field;
        let step = NonlinearWindingStep {
            native: trial.native.clone(),
            field: field_step,
            cells,
            graph_resistors,
            graph_materials: self.graph_materials.to_vec(),
            currents_a: solver.x.clone(),
            terminal_voltage_v: trial.voltage,
            flux_linkages_wb: [self.flux_linkages_wb.to_vec(), trial.linkage.clone()],
            inductive_voltage_v: trial.inductive_voltage_v,
            motional_voltage_v: trial.motional_voltage_v,
            terminal_residual_v: trial.voltage_residual,
            terminal_current_residual_a: trial.current_residual,
            stored_before_j,
            stored_after_j,
            stored_change_j,
            heat_j,
            supplied_work_j: external_work,
            numerical_dissipation_j,
            parameter_work_j,
            mechanical_work_j,
            energy_residual_j,
            field_trials: trials.get(),
            report: result.report,
        };
        let accepted = TerminalCandidate {
            field,
            state: trial.native.x,
            currents_a: solver.x,
            flux_linkages_wb: trial.linkage,
            step,
        };
        poll(&mut check)?;
        Ok(accepted)
    }
}

struct Trial<S, F> {
    native: StepRecord,
    field: FieldTrial<S, F>,
    end: EndWindingStep,
    x0: Vec<f64>,
    linkage: Vec<f64>,
    field_linkage_increments_wb: Vec<f64>,
    inductive_voltage_v: Vec<f64>,
    motional_voltage_v: Vec<f64>,
    voltage: Vec<f64>,
    voltage_residual: Vec<f64>,
    current_residual: Vec<f64>,
    residual: Vec<f64>,
}

pub(super) struct FieldAccount {
    pub(super) stored_before_j: f64,
    pub(super) stored_after_j: f64,
    pub(super) stored_change_j: f64,
    pub(super) heat_j: f64,
    pub(super) supplied_work_j: f64,
    pub(super) uniform_drive_work_j: f64,
    pub(super) motional_drive_work_j: f64,
    pub(super) numerical_dissipation_j: f64,
    pub(super) parameter_work_j: f64,
    pub(super) mechanical_work_j: f64,
}

pub(super) struct FieldTrial<S, F> {
    pub(super) state: S,
    pub(super) step: F,
    pub(super) account: FieldAccount,
    pub(super) motional_voltage_v: Option<Vec<f64>>,
}
impl<'m> FieldTrial<PlanarNonlinearEddy<'m>, NonlinearEddyStep> {
    fn stationary(candidate: PlanarNonlinearEddyCandidate<'m>) -> Self {
        let (state, step) = candidate.into_parts();
        let account = FieldAccount {
            stored_before_j: step.stored_before_j,
            stored_after_j: step.stored_after_j,
            stored_change_j: step.stored_after_j - step.stored_before_j,
            heat_j: step.heat_j,
            supplied_work_j: step.supplied_work_j,
            uniform_drive_work_j: step.uniform_drive_work_j,
            motional_drive_work_j: 0.0,
            numerical_dissipation_j: step.numerical_dissipation_j,
            parameter_work_j: step.parameter_work_j,
            mechanical_work_j: 0.0,
        };
        Self {
            state,
            step,
            account,
            motional_voltage_v: None,
        }
    }
}
