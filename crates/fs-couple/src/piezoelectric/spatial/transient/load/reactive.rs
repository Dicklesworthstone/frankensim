//! Source-bound reactive shunts on the original spatial midpoint equations.
//!
//! Each winding retains its flux linkage and each external capacitor its charge.
//! A series R-L branch obeys `lambda_dot = V - R I`, `I = lambda/L`;
//! eliminating its midpoint state gives
//! `I_mid = (lambda_0 + dt V_mid/2)/(L + dt R/2)`.
//! A capacitor obeys `q_dot = I`, `V = q/C`; its midpoint elimination is
//! `dt I_mid/2 = C V_mid - q_0`. These exact companion stamps enter the SAME
//! sparse reciprocal solid/dielectric solve as the existing Ohmic loads.
//! No separate circuit trajectory or artificial mechanical damping is used.
//!
//! Material sources and branch topology stay fixed. Inductors/capacitors use the original
//! `fs-material` uniform-field reductions, with their source and geometry
//! receipts. A resistor is optional on each inductive branch; its Joule heat
//! belongs to that conductor. Mutual inductance, internal circuit nodes,
//! nonlinear switching, thermal refresh and parasitic losses require other
//! models. Stored reactive energy and charge are never reset at step boundaries.

use super::super::acoustic::{TransientAcousticPort, TransientAcousticPortResponse};
use super::*;
use fs_material::field_storage::{ResolvedCapacitor, ResolvedInductor};

/// Fixed material/geometry branch between two original spatial electrodes.
/// Positive current and voltage run from the first electrode to the second.
#[derive(Clone, Debug)]
pub enum ReactiveElectrodeLoad {
    /// Material winding with optional series material resistance. `None` is
    /// an explicitly lossless winding; no resistivity is inferred from its core.
    SeriesRl {
        /// Original electrode indices in current/voltage orientation.
        electrodes: [usize; 2],
        /// Original magnetic source and actual winding geometry.
        inductor: Arc<ResolvedInductor>,
        /// Original conductor source and geometry, if a series resistor exists.
        conductor: Option<Arc<ResolvedConductor>>,
    },
    /// Additional material capacitor in parallel with the electrode pair.
    /// This is separate physical storage, not another copy of ceramic storage.
    Capacitor {
        /// Original electrode indices; charge is on the first terminal.
        electrodes: [usize; 2],
        /// Original dielectric source and actual gap geometry.
        capacitor: Arc<ResolvedCapacitor>,
    },
}

impl ReactiveElectrodeLoad {
    /// Original oriented terminal pair.
    pub const fn electrodes(&self) -> [usize; 2] {
        match self {
            Self::SeriesRl { electrodes, .. } | Self::Capacitor { electrodes, .. } => *electrodes,
        }
    }

    fn energy(&self, state: ReactiveElectrodeState) -> Result<f64> {
        let value = match (self, state) {
            (Self::SeriesRl { inductor, .. }, ReactiveElectrodeState::FluxLinkage(flux)) => {
                inductor.stored_energy_j(flux)
            }
            (Self::Capacitor { capacitor, .. }, ReactiveElectrodeState::Charge(charge)) => {
                capacitor.stored_energy_j(charge)
            }
            _ => {
                return Err(invalid(
                    "reactive storage coordinate must match its original load kind",
                ));
            }
        };
        value.map_err(|_| invalid("reactive load field energy is not representable"))
    }
}

/// Conserved electrical coordinate of one original reactive load.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ReactiveElectrodeState {
    /// Signed winding flux LINKAGE [Wb], not single-turn magnetic flux.
    FluxLinkage(f64),
    /// Signed charge on the first capacitor terminal [C].
    Charge(f64),
}

/// Actual midpoint flow and field-energy account of one named physical load.
#[derive(Clone, Debug)]
pub struct ReactiveElectrodeStep {
    /// Original immutable material sources, geometry and electrode orientation.
    pub load: ReactiveElectrodeLoad,
    /// Original and accepted physical storage coordinates.
    pub states: [ReactiveElectrodeState; 2],
    /// Oriented voltage across the complete branch [V].
    pub midpoint_voltage_v: f64,
    /// Oriented current through the complete branch [A].
    pub midpoint_current_a: f64,
    /// Signed charge transferred from first to second electrode [C].
    pub transferred_charge_c: f64,
    /// Original material field storage at the two endpoints [J].
    pub stored_energy_j: [f64; 2],
    /// Change of the original capacitor or winding storage [J].
    pub stored_energy_change_j: f64,
    /// Positive heat in the optional series conductor; zero for lossless storage.
    pub joule_heat_j: f64,
    /// Storage change plus heat minus actual voltage/charge work [J].
    pub energy_residual_j: f64,
    stored_change: Dd,
}

/// One independent account enclosing body, reactive storage and all conductors.
#[derive(Clone, Debug)]
pub struct ReactivePiezoelectricStep {
    /// Original spatial mechanical/electrical account, including net body work.
    pub body: TetPiezoelectricTransientStep,
    /// Interval-mean force of prescribed supports on the body [N], including
    /// inertia and intrinsic memory. Free mechanical coordinates are zero.
    pub midpoint_support_force_n: Vec<[f64; 3]>,
    /// External prescribed-current/voltage and volume-charge work [J].
    pub electrical_source_work_j: f64,
    /// Sum of load capacitor/inductor energy at the original/accepted endpoints.
    pub reactive_energy_j: [f64; 2],
    /// Sum of all memoryless and series-conductor losses [J].
    pub joule_heat_j: f64,
    /// Delta(body + reactive storage) + heat + signed acoustic outflow
    /// - mechanical/source work [J]. Outgoing waves are never conductor heat.
    pub energy_residual_j: f64,
    /// Original guide waves and physical work, in input port order.
    pub acoustic_ports: Vec<TransientAcousticPortResponse>,
    /// Original memoryless conductor currents and heat, in input order.
    pub resistors: Vec<ElectrodeLoadHeat>,
    /// Original reactive-branch currents, storage and heat, in input order.
    pub reactive: Vec<ReactiveElectrodeStep>,
}

/// Atomic accepted history of the spatial ceramic and its external storage.
/// Immutable source bundles are shared; only physical state is copied per step.
#[derive(Clone)]
pub struct ReactivePiezoelectricTransient<'model, 'mesh> {
    body: TetPiezoelectricTransient<'model, 'mesh>,
    loads: Arc<[ReactiveElectrodeLoad]>,
    states: Vec<ReactiveElectrodeState>,
    resistors: Arc<[ResistiveElectrodeLoad]>,
}

/// Joint candidate. The original body and all electrical histories are untouched.
pub struct ReactivePiezoelectricCandidate<'model, 'mesh> {
    model: ReactivePiezoelectricTransient<'model, 'mesh>,
    step: ReactivePiezoelectricStep,
}

impl<'model, 'mesh> ReactivePiezoelectricCandidate<'model, 'mesh> {
    /// Complete candidate body/load history after all physical audits pass.
    pub const fn model(&self) -> &ReactivePiezoelectricTransient<'model, 'mesh> {
        &self.model
    }
    /// Independently reconstructed storage, current, source and heat account.
    pub const fn step(&self) -> &ReactivePiezoelectricStep {
        &self.step
    }
    /// Publish the complete accepted joint state and its account.
    pub fn into_parts(
        self,
    ) -> (
        ReactivePiezoelectricTransient<'model, 'mesh>,
        ReactivePiezoelectricStep,
    ) {
        (self.model, self.step)
    }
}

impl<'model, 'mesh> ReactivePiezoelectricTransient<'model, 'mesh> {
    /// Bind fixed branches to an existing spatial state without projecting any
    /// charge, flux, displacement or velocity. A capacitor's initial charge must
    /// already agree with its terminal voltage within `charge_tolerance_c`.
    /// Terminals must share the original body's electrical reference component.
    pub fn new(
        body: TetPiezoelectricTransient<'model, 'mesh>,
        loads: Vec<ReactiveElectrodeLoad>,
        states: Vec<ReactiveElectrodeState>,
        resistors: Vec<ResistiveElectrodeLoad>,
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<Self> {
        poll(cx)?;
        validate_controls(controls)?;
        admit(&body, &loads, &states, controls, cx)?;
        for resistor in &resistors {
            poll(cx)?;
            terminal_nodes(&body, resistor.electrodes)?;
            if !resistor
                .conductor
                .source_property()
                .is_static_at_fixed_external_conditions()
            {
                return Err(invalid(
                    "reactive shunt resistors require fixed scalar Ohmic sources",
                ));
            }
        }
        Ok(Self {
            body,
            loads: loads.into(),
            states,
            resistors: resistors.into(),
        })
    }

    /// Release selected ideal voltage sources without advancing time or
    /// changing any physical storage. The native body retains actual conductor
    /// charges, mechanical motion and intrinsic memories; original capacitor
    /// charges and winding fluxes remain exactly unchanged.
    ///
    /// Material branches and their electrode topology remain attached. Later
    /// steps supply Current drives at released electrodes; zero removes the
    /// external source while the original shunts still exchange charge.
    /// The body's reference/representability admission applies transactionally.
    pub fn release_electrodes_candidate(&self, electrodes: &[usize], cx: &Cx<'_>) -> Result<Self> {
        let body = self.body.release_electrodes_candidate(electrodes, cx)?;
        let candidate = Self {
            body,
            loads: Arc::clone(&self.loads),
            states: self.states.clone(),
            resistors: Arc::clone(&self.resistors),
        };
        poll(cx)?;
        Ok(candidate)
    }

    /// Original spatial state. Advance this joint owner to retain load memory.
    pub const fn body(&self) -> &TetPiezoelectricTransient<'model, 'mesh> {
        &self.body
    }
    /// Original immutable source/geometry branches in physical state order.
    pub fn loads(&self) -> &[ReactiveElectrodeLoad] {
        &self.loads
    }
    /// Accepted signed flux linkages and charges in original branch order.
    pub fn states(&self) -> &[ReactiveElectrodeState] {
        &self.states
    }
    /// Original memoryless conductors, separate from series-R load conductors.
    pub fn resistors(&self) -> &[ResistiveElectrodeLoad] {
        &self.resistors
    }

    /// Advance motion, electric field, capacitor charge and winding flux in ONE
    /// midpoint solve. Each source type and all branch topology stay fixed.
    /// Attached intrinsic material memories remain in the native body solve;
    /// their heat is accounted separately from conductor losses. No thermal
    /// update or empirical damping is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(
        &self,
        dt_s: f64,
        drives: &[ElectrodeDrive],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<ReactivePiezoelectricCandidate<'model, 'mesh>> {
        self.advance_candidate(
            dt_s,
            None,
            &[],
            drives,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
    }

    /// Advance the same body/R/L/C system with prescribed support velocities.
    /// Each original constrained mechanical DOF must appear exactly once; its
    /// displacement is integrated by the native midpoint kinematics. Complete
    /// dynamic support work participates in the enclosing storage/heat audit.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_support_motion_candidate(
        &self,
        dt_s: f64,
        endpoint_support_velocities_m_s: &[(usize, f64)],
        drives: &[ElectrodeDrive],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<ReactivePiezoelectricCandidate<'model, 'mesh>> {
        self.advance_candidate(
            dt_s,
            Some(endpoint_support_velocities_m_s),
            &[],
            drives,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
    }

    /// Advance acoustic sensing/radiation and the original R/L/C histories in
    /// the same reciprocal midpoint solve. The fluid's traction acts back on
    /// motion while winding flux and capacitor charge remain conserved state.
    /// Guide work, intrinsic ceramic heat and every conductor's heat retain
    /// separate accounts. No acoustic or electrical history is reset.
    ///
    /// Supply `None` for stationary supports, or `Some` containing the endpoint
    /// velocity of every original support DOF. Moving supports pay their full
    /// dynamic reaction work, including outgoing sound.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_acoustic_ports_candidate(
        &self,
        dt_s: f64,
        acoustic_ports: &[TransientAcousticPort<'_>],
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        drives: &[ElectrodeDrive],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<ReactivePiezoelectricCandidate<'model, 'mesh>> {
        self.advance_candidate(
            dt_s,
            endpoint_support_velocities_m_s,
            acoustic_ports,
            drives,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_candidate(
        &self,
        dt_s: f64,
        endpoint_support_velocities_m_s: Option<&[(usize, f64)]>,
        acoustic_ports: &[TransientAcousticPort<'_>],
        drives: &[ElectrodeDrive],
        midpoint_nodal_forces_n: &[[f64; 3]],
        endpoint_nodal_free_charge_c: &[f64],
        controls: TetPiezoelectricControls,
        cx: &Cx<'_>,
    ) -> Result<ReactivePiezoelectricCandidate<'model, 'mesh>> {
        poll(cx)?;
        let endpoints = endpoint_conditions(&self.body, dt_s, drives, cx)?;
        let (candidate, resistors, reactive) = self.body.advance_with_acoustic_loads(
            dt_s,
            endpoint_support_velocities_m_s,
            &endpoints,
            &self.resistors,
            &self.loads,
            &self.states,
            acoustic_ports,
            midpoint_nodal_forces_n,
            endpoint_nodal_free_charge_c,
            controls,
            cx,
        )?;
        let mut outgoing = vec![Dd::ZERO; drives.len()];
        let mut joule_heat = Dd::ZERO;
        for flow in &resistors {
            poll(cx)?;
            outgoing[flow.load.electrodes[0]] =
                outgoing[flow.load.electrodes[0]] + Dd::from_f64(flow.transferred_charge_c);
            outgoing[flow.load.electrodes[1]] =
                outgoing[flow.load.electrodes[1]] - Dd::from_f64(flow.transferred_charge_c);
            joule_heat = joule_heat + Dd::from_f64(flow.joule_heat_j);
        }
        let mut reactive_energy_j = [0.0; 2];
        let mut reactive_change = Dd::ZERO;
        for flow in &reactive.steps {
            poll(cx)?;
            let [a, b] = flow.load.electrodes();
            outgoing[a] = outgoing[a] + Dd::from_f64(flow.transferred_charge_c);
            outgoing[b] = outgoing[b] - Dd::from_f64(flow.transferred_charge_c);
            joule_heat = joule_heat + Dd::from_f64(flow.joule_heat_j);
            reactive_change = reactive_change + flow.stored_change;
            for (sum, value) in reactive_energy_j.iter_mut().zip(flow.stored_energy_j) {
                *sum += value;
            }
        }
        let source_work = source_work(
            &self.body,
            &candidate.body.state,
            dt_s,
            drives,
            endpoint_nodal_free_charge_c,
            &outgoing,
            cx,
        )?;
        let electrical_source_work_j = source_work.to_f64();
        let joule_heat_j = joule_heat.to_f64();
        let energy_residual_j = (candidate.step.retained_energy_change_j()
            + reactive_change
            + joule_heat
            + Dd::from_f64(candidate.step.material_heat_j)
            + Dd::from_f64(candidate.step.acoustic_outward_work_j)
            - candidate.step.retained_mechanical_work_j()
            - source_work)
            .to_f64();
        if reactive_energy_j
            .iter()
            .chain([&joule_heat_j, &energy_residual_j])
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "joint reactive piezoelectric energy is not representable",
            ));
        }
        if energy_residual_j.abs() > controls.energy_tolerance_j {
            return Err(TetPiezoelectricError::Balance(
                "body/reactive storage, heat, acoustic outflow and external source work disagree",
            ));
        }
        let TetPiezoelectricTransientCandidate {
            body,
            step: body_step,
            midpoint_support_force_n,
            acoustic_ports,
        } = candidate;
        // Repeat endpoint admission before publishing: a charge/voltage defect
        // cannot be concealed by the total energy account.
        admit(&body, &self.loads, &reactive.states, controls, cx)?;
        poll(cx)?;
        Ok(ReactivePiezoelectricCandidate {
            model: Self {
                body,
                loads: Arc::clone(&self.loads),
                states: reactive.states,
                resistors: Arc::clone(&self.resistors),
            },
            step: ReactivePiezoelectricStep {
                body: body_step,
                midpoint_support_force_n,
                electrical_source_work_j,
                reactive_energy_j,
                joule_heat_j,
                energy_residual_j,
                acoustic_ports,
                resistors,
                reactive: reactive.steps,
            },
        })
    }
}

fn terminal_nodes(
    body: &TetPiezoelectricTransient<'_, '_>,
    electrodes: [usize; 2],
) -> Result<[usize; 2]> {
    let [a, b] = electrodes;
    if a == b
        || a >= body.boundary.electrode_nodes.len()
        || b >= body.boundary.electrode_nodes.len()
    {
        return Err(invalid(
            "reactive loads require distinct original in-range electrodes",
        ));
    }
    let nodes = [
        body.boundary.electrode_nodes[a][0],
        body.boundary.electrode_nodes[b][0],
    ];
    if body.boundary.map.voltage_reference_electrode[nodes[0]]
        != body.boundary.map.voltage_reference_electrode[nodes[1]]
    {
        return Err(invalid(
            "reactive terminals must share the original dielectric reference component",
        ));
    }
    Ok(nodes)
}

fn admit(
    body: &TetPiezoelectricTransient<'_, '_>,
    loads: &[ReactiveElectrodeLoad],
    states: &[ReactiveElectrodeState],
    controls: TetPiezoelectricControls,
    cx: &Cx<'_>,
) -> Result<()> {
    if loads.len() != states.len() {
        return Err(invalid(
            "complete original reactive storage coordinates required",
        ));
    }
    for (load, &state) in loads.iter().zip(states) {
        poll(cx)?;
        let [a, b] = terminal_nodes(body, load.electrodes())?;
        load.energy(state)?;
        match (load, state) {
            (
                ReactiveElectrodeLoad::SeriesRl {
                    inductor,
                    conductor,
                    ..
                },
                ReactiveElectrodeState::FluxLinkage(flux),
            ) => {
                if !(flux / inductor.inductance_h()).is_finite()
                    || conductor.as_ref().is_some_and(|r| {
                        !r.source_property().is_static_at_fixed_external_conditions()
                    })
                {
                    return Err(invalid(
                        "reactive winding current and fixed Ohmic source must be representable",
                    ));
                }
            }
            (
                ReactiveElectrodeLoad::Capacitor { capacitor, .. },
                ReactiveElectrodeState::Charge(charge),
            ) => {
                let voltage = body.state.relative_voltage_v[a] - body.state.relative_voltage_v[b];
                let residual = capacitor.capacitance_f().mul_add(voltage, -charge);
                if !residual.is_finite() || residual.abs() > controls.charge_tolerance_c {
                    return Err(invalid(
                        "load capacitor charge is incompatible with the spatial terminal voltage",
                    ));
                }
            }
            _ => unreachable!("state kind checked by original storage energy"),
        }
    }
    Ok(())
}

struct Companion {
    coefficient_f: f64,
    history_c: f64,
}

fn companion(
    load: &ReactiveElectrodeLoad,
    state: ReactiveElectrodeState,
    dt: f64,
) -> Result<Companion> {
    let value = match (load, state) {
        (
            ReactiveElectrodeLoad::SeriesRl {
                inductor,
                conductor,
                ..
            },
            ReactiveElectrodeState::FluxLinkage(flux),
        ) => {
            let half_dt = 0.5 * dt;
            let resistance = conductor.as_ref().map_or(0.0, |r| r.resistance_ohm());
            let denominator = half_dt.mul_add(resistance, inductor.inductance_h());
            if !positive(denominator) {
                return Err(invalid(
                    "reactive winding midpoint impedance is not representable",
                ));
            }
            Companion {
                coefficient_f: half_dt * (half_dt / denominator),
                history_c: half_dt * (flux / denominator),
            }
        }
        (
            ReactiveElectrodeLoad::Capacitor { capacitor, .. },
            ReactiveElectrodeState::Charge(charge),
        ) => Companion {
            coefficient_f: capacitor.capacitance_f(),
            history_c: -charge,
        },
        _ => {
            return Err(invalid(
                "reactive storage coordinate must match its original load kind",
            ));
        }
    };
    if !positive(value.coefficient_f) || !value.history_c.is_finite() {
        return Err(invalid("reactive midpoint companion is not representable"));
    }
    Ok(value)
}

pub(in crate::piezoelectric::spatial::transient) fn stamp(
    body: &TetPiezoelectricTransient<'_, '_>,
    loads: &[ReactiveElectrodeLoad],
    states: &[ReactiveElectrodeState],
    dt: f64,
    boundary: &mut BoundaryMap,
    controls: TetPiezoelectricControls,
    cx: &Cx<'_>,
) -> Result<Vec<(usize, usize, f64)>> {
    admit(body, loads, states, controls, cx)?;
    let mut entries = Vec::with_capacity(
        loads
            .len()
            .checked_mul(4)
            .ok_or_else(|| invalid("reactive load count overflows"))?,
    );
    let nu = 3 * body.model.nodes_m.len();
    for (load, &state) in loads.iter().zip(states) {
        poll(cx)?;
        let [a, b] = terminal_nodes(body, load.electrodes())?;
        let value = companion(load, state, dt)?;
        let initial_voltage = body.state.relative_voltage_v[a] - body.state.relative_voltage_v[b];
        let charge = value
            .coefficient_f
            .mul_add(initial_voltage, value.history_c);
        if !charge.is_finite() {
            return Err(invalid("reactive companion charge is not representable"));
        }
        for (node, sign) in [(a, 1.0), (b, -1.0)] {
            if let Some(row) = boundary.index[nu + node] {
                boundary.rhs[row] += sign * charge;
            }
        }
        for (row, col, sign) in [(a, a, -1.0), (a, b, 1.0), (b, a, 1.0), (b, b, -1.0)] {
            entries.push((nu + row, nu + col, sign * value.coefficient_f));
        }
    }
    Ok(entries)
}

pub(in crate::piezoelectric::spatial::transient) struct ReactiveUpdate {
    states: Vec<ReactiveElectrodeState>,
    steps: Vec<ReactiveElectrodeStep>,
}

#[allow(clippy::too_many_arguments)]
pub(in crate::piezoelectric::spatial::transient) fn recover(
    body: &TetPiezoelectricTransient<'_, '_>,
    loads: &[ReactiveElectrodeLoad],
    states: &[ReactiveElectrodeState],
    dt: f64,
    midpoint_voltage: &[f64],
    endpoint: &mut [ElectrodeCondition],
    normalized_endpoint: &mut [ElectrodeCondition],
    midpoint: &mut [ElectrodeCondition],
    controls: TetPiezoelectricControls,
    cx: &Cx<'_>,
) -> Result<ReactiveUpdate> {
    let mut updated = ReactiveUpdate {
        states: Vec::with_capacity(loads.len()),
        steps: Vec::with_capacity(loads.len()),
    };
    if loads.is_empty() {
        return Ok(updated);
    }
    for (load, &state) in loads.iter().zip(states) {
        poll(cx)?;
        let [a, b] = terminal_nodes(body, load.electrodes())?;
        let voltage = midpoint_voltage[a] - midpoint_voltage[b];
        let value = companion(load, state, dt)?;
        let charge = 2.0 * value.coefficient_f.mul_add(voltage, value.history_c);
        let current = charge / dt;
        let (next, heat) = match (load, state) {
            (
                ReactiveElectrodeLoad::SeriesRl {
                    inductor,
                    conductor,
                    ..
                },
                ReactiveElectrodeState::FluxLinkage(flux),
            ) => {
                let next = ReactiveElectrodeState::FluxLinkage(
                    (2.0 * inductor.inductance_h()).mul_add(current, -flux),
                );
                let heat = match conductor {
                    Some(r) => {
                        dt * r.joule_power_w(current).map_err(|_| {
                            invalid("reactive series-conductor Joule heat is not representable")
                        })?
                    }
                    None => 0.0,
                };
                (next, heat)
            }
            (ReactiveElectrodeLoad::Capacitor { .. }, ReactiveElectrodeState::Charge(old)) => {
                (ReactiveElectrodeState::Charge(old + charge), 0.0)
            }
            _ => unreachable!("reactive sources admitted before assembly"),
        };
        let energy = [load.energy(state)?, load.energy(next)?];
        let (old_coordinate, new_coordinate, storage) = match (load, state, next) {
            (
                ReactiveElectrodeLoad::SeriesRl { inductor, .. },
                ReactiveElectrodeState::FluxLinkage(old),
                ReactiveElectrodeState::FluxLinkage(new),
            ) => (old, new, inductor.inductance_h()),
            (
                ReactiveElectrodeLoad::Capacitor { capacitor, .. },
                ReactiveElectrodeState::Charge(old),
                ReactiveElectrodeState::Charge(new),
            ) => (old, new, capacitor.capacitance_f()),
            _ => unreachable!("reactive storage kind admitted above"),
        };
        let mut products = Products::default();
        products.add(
            (Dd::from_f64(0.5 * old_coordinate) + Dd::from_f64(0.5 * new_coordinate))
                / Dd::from_f64(storage),
            Dd::from_f64(new_coordinate) - Dd::from_f64(old_coordinate),
        );
        let stored_change = products.value();
        let change = stored_change.to_f64();
        let residual = (stored_change + Dd::from_f64(heat)
            - Dd::from_f64(voltage) * Dd::from_f64(charge))
        .to_f64();
        if [voltage, charge, current, heat, change, residual]
            .iter()
            .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "reactive current, storage or work is not representable",
            ));
        }
        if residual.abs() > controls.energy_tolerance_j {
            return Err(TetPiezoelectricError::Balance(
                "reactive field energy and series heat disagree with terminal work",
            ));
        }
        for (id, flow) in [
            (load.electrodes()[0], charge),
            (load.electrodes()[1], -charge),
        ] {
            if let ElectrodeCondition::Charge(q) = &mut endpoint[id] {
                *q -= flow;
            }
        }
        updated.states.push(next);
        updated.steps.push(ReactiveElectrodeStep {
            load: load.clone(),
            states: [state, next],
            midpoint_voltage_v: voltage,
            midpoint_current_a: current,
            transferred_charge_c: charge,
            stored_energy_j: energy,
            stored_energy_change_j: change,
            joule_heat_j: heat,
            energy_residual_j: residual,
            stored_change,
        });
    }
    update_charge_conditions(body, endpoint, normalized_endpoint, midpoint, cx)?;
    Ok(updated)
}
