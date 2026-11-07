//! Native winding circuits retain the original spatial vector hysteresis.
//!
//! Each terminal trial starts from the SAME accepted A and every branch
//! polarization. The existing nonlinear terminal solver pairs endpoint current
//! with `F^T delta(A)/dt + L_end delta(I)/dt + R I`; the original field owner
//! supplies remanence, saturation, pinning loss, source-selected intrinsic
//! switching rates and resolved core eddy diffusion.
//! No effective inductance, scalar-memory replacement or empirical loss map is
//! substituted. Private native-port work cancels once against actual field work.
//!
//! The fixed planar geometry contains stationary, disjoint stranded windings.
//! Magnetic fits remain constant on their source-admitted temperature support;
//! source-bound copper/core conductivity and graph resistance may refresh.
//! Physical pinning, viscous and Joule heat stay separate from backward-Euler
//! numerical loss. Motion, solid current injection, changing magnetic fits and physical
//! validation remain outside this estimated, default-off rung.

use super::{
    MAX_WINDINGS, admit_profiles,
    end::{EndWindingInductance, EndWindingMode},
    nonlinear::{
        FieldAccount, FieldTrial, Layout, NonlinearWindingConfig, NonlinearWindingControls,
        NonlinearWindingInitial, NonlinearWindingPort, NonlinearWindingStep, NonlinearWindings,
        PlanarInductionHistory, Result, TerminalContext, finite, impressed, linkages, poll,
        resistors::{self, GraphResistorMaterial},
    },
    ohmic::{self, OhmicWindingError, WindingConductor},
};
use crate::planar::eddy::{
    EddyConductor,
    hysteretic::{HystereticEddyControls, HystereticEddyStep, PlanarHystereticEddy},
};
use fs_material::field_storage::admit_source_temperature_refresh;
use fs_phs::circuit::{CircuitDae, CircuitGraph};
use fs_rep_mesh::{Metric2, TriComplex2};
use std::rc::Rc;

/// Original field, source, native circuit and bounded terminal-solver refusals.
pub use super::nonlinear::NonlinearWindingError as HystereticWindingError;
/// The existing terminal controls with the actual vector-history field limits.
pub type HystereticWindingControls = NonlinearWindingControls<HystereticEddyControls>;
/// All original native graph, winding, vector-history and core-eddy accounts.
pub type HystereticWindingStep = NonlinearWindingStep<HystereticEddyStep>;

/// One accepted circuit and spatial magnetic memory, published atomically.
#[derive(Clone)]
pub struct HystereticWindingCircuit<'m> {
    layout: Rc<Layout>,
    end: Rc<EndWindingInductance>,
    graph_materials: Vec<GraphResistorMaterial>,
    field: PlanarHystereticEddy<'m>,
    conductors: Vec<WindingConductor>,
    conductor_windings: Vec<usize>,
    resistance_ohm: Vec<f64>,
    state: Vec<f64>,
    currents_a: Vec<f64>,
    flux_linkages_wb: Vec<f64>,
}

/// A proposal has no effect on its original field, vector memories or circuit.
pub struct HystereticWindingCandidate<'m> {
    circuit: HystereticWindingCircuit<'m>,
    step: HystereticWindingStep,
}
impl<'m> HystereticWindingCandidate<'m> {
    pub fn circuit(&self) -> &HystereticWindingCircuit<'m> {
        &self.circuit
    }
    pub fn step(&self) -> &HystereticWindingStep {
        &self.step
    }
    pub fn into_parts(self) -> (HystereticWindingCircuit<'m>, HystereticWindingStep) {
        (self.circuit, self.step)
    }
}

impl PlanarInductionHistory for PlanarHystereticEddy<'_> {
    fn mesh(&self) -> &TriComplex2 {
        self.problem().mesh()
    }
    fn potential_wb_m(&self) -> &[f64] {
        self.potential_wb_m()
    }
    fn time_s(&self) -> f64 {
        self.time_s()
    }
}

fn bind(
    field: &PlanarHystereticEddy<'_>,
    profiles: &[Vec<f64>],
    sources: Vec<WindingConductor>,
    previous: Option<&[WindingConductor]>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<WindingConductor>, Vec<usize>, Vec<f64>)> {
    // Native field binding sorts core conductors and admits a complete cell
    // partition. Avoid rescanning every insulating cell for each copper cell.
    for source in &sources {
        poll(cancelled)?;
        if source.triangle >= field.problem().mesh().faces().len()
            || field
                .conductors()
                .binary_search_by_key(&source.triangle, |cell| cell.triangle)
                .is_ok()
        {
            return Err(HystereticWindingError::Invalid(
                "stranded winding and resolved eddy conductor cells must be disjoint",
            ));
        }
    }
    if field.conductors().iter().any(|source| {
        profiles
            .iter()
            .any(|profile| profile[source.triangle] != 0.0)
    }) {
        return Err(HystereticWindingError::Invalid(
            "hysteretic stranded profiles cannot inject current into resolved core conductors",
        ));
    }
    Ok(ohmic::bind_sources(
        field.problem().mesh(),
        profiles,
        sources,
        previous,
        |source| {
            let original = &field.problem().materials()[source.triangle];
            let electrical = source.conductor.material();
            if original.material().card_identity() != electrical.card_identity()
                || original.material().material() != electrical.material()
            {
                return Err(OhmicWindingError::Invalid(
                    "magnetic fit and winding conductor require the same original material card",
                ));
            }
            let checked = original.at_query_point(&electrical.query())?;
            admit_source_temperature_refresh(original.material(), checked.material())?;
            for ((axis, before), (_, after)) in original
                .material()
                .query_point()
                .iter()
                .zip(checked.material().query_point())
            {
                if before.to_bits() != after.to_bits() {
                    original.admit_temperature_interval(
                        axis,
                        [before.min(*after), before.max(*after)],
                    )?;
                }
            }
            Ok(())
        },
        cancelled,
    )?)
}

impl<'m> HystereticWindingCircuit<'m> {
    /// Bind physical field history to original inductor placeholders. No
    /// magnetic state is equilibrated, projected or reset during construction.
    pub fn new(
        field: PlanarHystereticEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        sources: Vec<WindingConductor>,
        initial: NonlinearWindingInitial,
        config: NonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let count = windings.profiles_per_m2.len();
        if count == 0 || count > config.max_windings || config.max_windings > MAX_WINDINGS {
            return Err(HystereticWindingError::Invalid(
                "invalid winding count or capacity",
            ));
        }
        Self::new_with_ports(
            field,
            graph,
            windings,
            vec![NonlinearWindingPort::CurrentInput; count],
            sources,
            initial,
            config,
            cancelled,
        )
    }

    /// Explicit native causality supports both voltage and current excitation.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_ports(
        field: PlanarHystereticEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        ports: Vec<NonlinearWindingPort>,
        sources: Vec<WindingConductor>,
        initial: NonlinearWindingInitial,
        config: NonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        Self::new_with_end_modes(
            field,
            graph,
            windings,
            ports,
            vec![],
            sources,
            initial,
            config,
            cancelled,
        )
    }

    /// Explicit shared end-flux modes add only unresolved self/mutual storage.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_end_modes(
        field: PlanarHystereticEddy<'m>,
        graph: CircuitGraph,
        windings: NonlinearWindings,
        ports: Vec<NonlinearWindingPort>,
        end_modes: Vec<EndWindingMode>,
        sources: Vec<WindingConductor>,
        initial: NonlinearWindingInitial,
        config: NonlinearWindingConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let count = windings.profiles_per_m2.len();
        let mesh = field.problem().mesh();
        if count == 0
            || count > config.max_windings
            || config.max_windings > MAX_WINDINGS
            || count
                .checked_mul(mesh.faces().len())
                .is_none_or(|n| n > config.max_profile_values)
            || windings.branch_indices.len() != count
            || windings.end_inductance_h.len() != count
            || !config.current_balance_tolerance_turns.is_finite()
            || config.current_balance_tolerance_turns < 0.0
        {
            return Err(HystereticWindingError::Invalid(
                "invalid bounded winding geometry or limits",
            ));
        }
        let Metric2::Planar { thickness } = mesh.metric() else {
            return Err(HystereticWindingError::Invalid(
                "hysteretic windings require the planar magnetic metric",
            ));
        };
        admit_profiles(
            mesh,
            &windings.profiles_per_m2,
            thickness,
            config.current_balance_tolerance_turns,
            cancelled,
        )?;
        let expected = impressed(&windings.profiles_per_m2, &initial.currents_a, cancelled)?;
        if expected
            .iter()
            .zip(field.impressed_current_density_a_m2())
            .any(|(a, b)| (a - b).abs() > 128.0 * f64::EPSILON * a.abs().max(b.abs()))
        {
            return Err(HystereticWindingError::Invalid(
                "initial currents disagree with the accepted impressed field history",
            ));
        }
        let end = Rc::new(EndWindingInductance::new(
            windings.end_inductance_h.clone(),
            end_modes,
        )?);
        let (conductors, conductor_windings, resistance_ohm) =
            bind(&field, &windings.profiles_per_m2, sources, None, cancelled)?;
        let flux_linkages_wb = linkages(
            &field,
            &windings.profiles_per_m2,
            &end,
            &initial.currents_a,
            cancelled,
        )?;
        let (layout, state) = Layout::new(graph, windings, ports, &initial, cancelled)?;
        let total = finite(
            field.energy_j()
                + layout.dae.system.hamiltonian(&state)
                + end.stored_energy_j(&initial.currents_a)?,
        )?;
        if total < 0.0 {
            return Err(HystereticWindingError::Invalid(
                "negative initial field/circuit storage",
            ));
        }
        poll(cancelled)?;
        Ok(Self {
            layout,
            end,
            graph_materials: vec![],
            field,
            conductors,
            conductor_windings,
            resistance_ohm,
            state,
            currents_a: initial.currents_a,
            flux_linkages_wb,
        })
    }

    pub fn field(&self) -> &PlanarHystereticEddy<'m> {
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
    pub fn terminal_ports(&self) -> &[NonlinearWindingPort] {
        &self.layout.terminal_ports
    }
    pub fn end_inductance(&self) -> &EndWindingInductance {
        &self.end
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
    pub fn graph_resistor_materials(&self) -> &[GraphResistorMaterial] {
        &self.graph_materials
    }

    /// Original source-bound graph resistance changes without moving stored
    /// flux, charge, vector polarization, winding current or accepted time.
    pub fn with_graph_resistors(
        &self,
        sources: Vec<GraphResistorMaterial>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let (layout, sources) = resistors::rebind(&self.layout, &self.state, sources, cancelled)?;
        let mut candidate = self.clone();
        candidate.layout = layout;
        candidate.graph_materials = sources;
        Ok(candidate)
    }

    /// Solve the actual circuit, vector-history field and core diffusion from
    /// one accepted time. Endpoint current excites the field; source voltages
    /// are held over the interval. Source refresh preserves the same fitted
    /// magnetic law and does not erase pinning memory between terminal trials.
    ///
    /// This first-order composition returns all spatial copper and constitutive
    /// heat once. End/field numerical loss is never physical heating. A thermal
    /// receiver must accept its own corresponding caloric candidate atomically.
    pub fn step_to_candidate(
        &self,
        target_s: f64,
        core_conductors: Vec<EddyConductor>,
        winding_conductors: Vec<WindingConductor>,
        external_inputs: &[f64],
        controls: HystereticWindingControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HystereticWindingCandidate<'m>> {
        poll(cancelled)?;
        let (conductors, conductor_windings, resistance_ohm) = bind(
            &self.field,
            self.profiles_per_m2(),
            winding_conductors,
            Some(&self.conductors),
            cancelled,
        )?;
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
            self.profiles_per_m2(),
            &conductors,
            &conductor_windings,
            &resistance_ohm,
            &resistance_ohm,
            external_inputs,
            controls,
            |_, impressed, cancelled| {
                let (state, step) = self
                    .field
                    .step_to_candidate(
                        target_s,
                        core_conductors.clone(),
                        impressed,
                        controls.field,
                        &mut || cancelled(),
                    )?
                    .into_parts();
                let account = FieldAccount {
                    stored_before_j: step.stored_before_j,
                    stored_after_j: step.stored_after_j,
                    stored_change_j: step.magnetic.stored_energy_change_j,
                    heat_j: step.heat_j,
                    supplied_work_j: step.supplied_work_j,
                    uniform_drive_work_j: 0.0,
                    motional_drive_work_j: 0.0,
                    numerical_dissipation_j: step.numerical_dissipation_j,
                    parameter_work_j: 0.0,
                    mechanical_work_j: 0.0,
                };
                Ok(FieldTrial {
                    state,
                    step,
                    account,
                    motional_voltage_v: None,
                })
            },
            cancelled,
        )?;
        let circuit = Self {
            layout: self.layout.clone(),
            end: self.end.clone(),
            graph_materials: self.graph_materials.clone(),
            field: accepted.field,
            conductors,
            conductor_windings,
            resistance_ohm,
            state: accepted.state,
            currents_a: accepted.currents_a,
            flux_linkages_wb: accepted.flux_linkages_wb,
        };
        Ok(HystereticWindingCandidate {
            circuit,
            step: accepted.step,
        })
    }
}
