//! The battery-motor drive's fixtures as declared graphs on the composition
//! kernel.
//!
//! `battery_motor` and its `capacitance` and `charge` lanes couple fs-phs's
//! commutated drive to the enthalpy timeline by hand. A step runs the gate
//! program once, each resistance and switch at its receiving cell's card at
//! the END temperature, with Newton on the cell energies. The graph says the
//! same with no code for the pair:
//! - `fs_system::commutated_drive::CommutatedBatteryMotorDrive` in the
//!   fixture's lane ([`Lane`]), its five resistances following
//!   [`CardResistance`] and its four switches following [`CardSwitch`];
//! - two `elements::Force` sources holding the command's shaft torques;
//! - the seven cells, and the fixture's coolant if it has one ([`Coolant`]),
//!   as one backward-Euler timeline whose cell ports return end-of-step
//!   temperatures. The coolant has no port;
//! - one heat route per heat port into the module's receiving cell;
//! - where the fixture links two cells, the same link in the timeline at the
//!   module's fixed conductance ([`FixedCoolingStandIn`]).
//!
//! The two sides solve the same equations with different solvers, so they
//! agree to within the solves' tolerances.

use std::fmt::Debug;
use std::sync::Arc;

use fs_material::state_point::{
    ResolvedScalarProperty, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_phs::circuit::converter::commutation::battery_motor::PreparedCommutatedBatteryMotorCycle;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    PreparedCapacitiveBatteryMotorCycle, PreparedChargeBatteryMotorCycle,
};
use fs_ports::PortKind;
use fs_system::circuit::ResistanceLaw;
use fs_system::commutated_drive::{
    CommutatedBatteryMotorDrive, DRIVER_HEAT_PORTS, DeviceStorage, DriveLane, DriveParameters,
    DriveResistor, GateProgram, HEAT_PORTS, LOAD_PORT, NoDeviceStorage, RAIL_PORTS, SwitchLaw,
    SwitchPosition,
};
use fs_system::elements::{Force, SignalClock, SignalSource, SourceAccount, Waveform};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::link::{ImplicitLinkLaw, ImplicitLinkState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductiveLink, ConductivityLaw};

use super::*;

/// The cell each heat port heats, in `HEAT_PORTS` order: the module's
/// channels (`drive_resistors`), device packages (`Fixture::devices`) and
/// bearing cell, the grouping of `per_body`.
const RECEIVERS: [CellId; 10] = [
    IDS[0], IDS[1], IDS[2], IDS[3], IDS[4], IDS[5], IDS[6], IDS[5], IDS[6], IDS[0],
];

/// One fixture of a lane, for both sides.
pub(super) struct Case<'s, L: Lane = DeviceStorage> {
    pub(super) fixture: &'static Fixture,
    /// The query context of the module's cells.
    pub(super) context: &'static QueryPoint,
    /// Each package's storage claim, in switch order; none without switch
    /// storage.
    pub(super) storage_keys: &'static [PropertyKey],
    pub(super) storage: L,
    pub(super) steady: &'s SteadyDrive,
    pub(super) start: L::State,
    pub(super) dt_s: f64,
    pub(super) substeps: usize,
    /// The cell that receives the auxiliary load's work; none leaves it
    /// outside.
    pub(super) load_cell: Option<CellId>,
    /// The gates of each step in turn; empty runs the steady gates in every
    /// step.
    pub(super) schedule: Vec<DeviceGateState>,
    /// Whether the cells carry the module's heat-assembly link ([`LINK`]).
    pub(super) linked: bool,
    /// The module's finite liquid coolant, linked to one of the cells.
    pub(super) coolant: Option<&'static Coolant>,
    /// The gate lane's drivers; none in the other lanes.
    pub(super) drivers: Option<&'static Drivers>,
}

/// The module's finite liquid coolant: one more cell with a caloric law of
/// its own, joined by the link [`LINK`] names to the cell `partner` names,
/// at the conductance it gives \[W/K\].
pub(super) struct Coolant {
    pub(super) card: MaterialCard,
    pub(super) caloric: ResolvedHeatCapacityCurve,
    pub(super) context: QueryPoint,
    pub(super) cell: DepositedCell,
    pub(super) partner: (CellId, f64),
}

/// The gate lane's drivers as the module's gate fixture authors them: a
/// cell per gate for its driver's heat, the card keys of the on and off
/// resistances, and the rails' voltages \[V\] in [`RAIL_PORTS`] order.
/// `start` holds every cell at the start, the seven and the four.
pub(super) struct Drivers {
    pub(super) cells: [CellId; 4],
    pub(super) keys: [PropertyKey; 2],
    pub(super) rails: [f64; 8],
    pub(super) start: DepositionSimulation<'static>,
}

/// The case's link: [`LINK`], or with a coolant the same link from its
/// partner's cell to the coolant's.
fn link<L: Lane>(case: &Case<'_, L>) -> Option<(u64, [CellId; 2], f64)> {
    match case.coolant {
        Some(coolant) => {
            let (cell, conductance_w_k) = coolant.partner;
            Some((LINK.0, [cell, coolant.cell.id], conductance_w_k))
        }
        None => case.linked.then_some(LINK),
    }
}

/// The cells with timeline ports: the seven, then the drivers' if any.
fn port_cells<L: Lane>(case: &Case<'_, L>) -> Vec<CellId> {
    IDS.into_iter()
        .chain(case.drivers.into_iter().flat_map(|drivers| drivers.cells))
        .collect()
}

/// The timeline's cells: those with ports, then the coolant if any.
pub(super) fn cell_ids<L: Lane>(case: &Case<'_, L>) -> Vec<CellId> {
    let mut ids = port_cells(case);
    ids.extend(case.coolant.map(|coolant| coolant.cell.id));
    ids
}

/// A drive lane on the graph.
pub(super) trait Lane: DriveLane {
    /// `drive` with its gate drivers' resistances following the case's
    /// driver claims. Only the gate lane has drivers; the others add nothing.
    fn add_drivers(
        drive: CommutatedBatteryMotorDrive<Self>,
        _: &Case<'_, Self>,
    ) -> CommutatedBatteryMotorDrive<Self> {
        drive
    }
}

/// A lane whose owner the harness steps by hand from its authored storage
/// (`owner_step`). The gate lane is not one: its channels and drivers take
/// their cells' temperatures, so `gate/graph.rs` builds its owner.
pub(super) trait HandOwner: Lane {
    /// One owner cycle of `intervals` from `old` under `command`: the new
    /// state, the fourteen heats and the load's work \[J\].
    fn cycle(
        &self,
        original: BatteryMotorDrive,
        devices: [SwitchDevice; 4],
        intervals: &[CommutationInterval],
        substeps: usize,
        command: BatteryMotorSwitchCommand,
        old: &Self::State,
    ) -> (Self::State, [f64; 14], f64);
}

impl Lane for NoDeviceStorage {}

impl HandOwner for NoDeviceStorage {
    fn cycle(
        &self,
        original: BatteryMotorDrive,
        devices: [SwitchDevice; 4],
        intervals: &[CommutationInterval],
        substeps: usize,
        command: BatteryMotorSwitchCommand,
        old: &[f64; 9],
    ) -> ([f64; 9], [f64; 14], f64) {
        let (e, c) = (power_controls(), cycle_controls());
        let mut next = [0.0; 9];
        let step = PreparedCommutatedBatteryMotorCycle::new_controlled(
            original,
            devices,
            intervals,
            substeps,
            e,
            c,
            || false,
        )
        .unwrap()
        .step_interval_commands_controlled(old, &[command], &mut next, || false)
        .unwrap();
        (next, step.heat_joules(), step.load_work_j)
    }
}

impl Lane for DeviceStorage {}

impl HandOwner for DeviceStorage {
    fn cycle(
        &self,
        original: BatteryMotorDrive,
        devices: [SwitchDevice; 4],
        intervals: &[CommutationInterval],
        n: usize,
        command: BatteryMotorSwitchCommand,
        old: &[f64; 11],
    ) -> ([f64; 11], [f64; 14], f64) {
        let (commands, e, c) = ([command], power_controls(), cycle_controls());
        let mut next = [0.0; 11];
        let step = match self {
            Self::Capacitance(storage) => PreparedCapacitiveBatteryMotorCycle::new_controlled(
                original,
                devices,
                *storage,
                intervals,
                n,
                e,
                c,
                || false,
            )
            .unwrap()
            .step_interval_commands_controlled(old, &commands, &mut next, || false),
            Self::Charge {
                curves,
                maximum_cells_per_chord: chords,
            } => PreparedChargeBatteryMotorCycle::new_charge_controlled(
                original,
                devices,
                curves.clone(),
                *chords,
                intervals,
                n,
                e,
                c,
                || false,
            )
            .unwrap()
            .step_interval_commands_controlled(old, &commands, &mut next, || false),
        }
        .unwrap();
        (next, step.heat_joules(), step.load_work_j)
    }
}

/// The link of the module's heat-assembly fixtures: its id, its two cells
/// and the conductance `FixedCooling` gives it \[W/K\].
pub(super) const LINK: (u64, [CellId; 2], f64) = (77, [IDS[0], IDS[6]], 0.25);

impl<L: Lane> Case<'_, L> {
    /// The gates of step `k`.
    pub(super) fn gates(&self, k: usize) -> DeviceGateState {
        self.schedule.get(k).copied().unwrap_or(self.steady.gates)
    }
}

/// One accepted step: the drive's coordinates and the cells' temperatures
/// \[K\], the seven in `IDS` order, then the coolant if any.
#[derive(Debug, Clone, Copy)]
pub(super) struct Row<S = [f64; 11], K = [f64; 7]> {
    pub(super) state: S,
    pub(super) kelvin: K,
}

pub(super) fn leak<T>(value: T) -> &'static T {
    Box::leak(Box::new(value))
}

pub(super) fn text(error: impl Debug) -> String {
    format!("{error:?}")
}

/// The specific enthalpy \[J/kg\] at `kelvin` on the fixture's caloric
/// curve, by the curve's own linear interpolation.
fn enthalpy_at(curve: &ResolvedHeatCapacityCurve, kelvin: f64) -> Result<f64, String> {
    for pair in curve.curve().knots().windows(2) {
        let (low, high) = (&pair[0], &pair[1]);
        if (low.temperature_k..=high.temperature_k).contains(&kelvin)
            && high.temperature_k > low.temperature_k
        {
            let fraction = (kelvin - low.temperature_k) / (high.temperature_k - low.temperature_k);
            return Ok(low.specific_enthalpy_j_kg
                + fraction * (high.specific_enthalpy_j_kg - low.specific_enthalpy_j_kg));
        }
    }
    Err(format!("{kelvin} K is outside the caloric curve"))
}

/// The fixture's card at the phase its caloric curve gives a temperature,
/// queried at the module's cell context: what the module reads for a cell
/// at its end-of-step state.
#[derive(Debug, Clone, Copy)]
pub(super) struct Card {
    pub(super) material: &'static MaterialCard,
    pub(super) curve: &'static ResolvedHeatCapacityCurve,
    pub(super) context: &'static QueryPoint,
}

impl Card {
    /// The module's query for a cell at `kelvin`.
    pub(super) fn query(&self, kelvin: f64) -> Result<QueryPoint, String> {
        let h = enthalpy_at(self.curve, kelvin)?;
        let phase = self
            .curve
            .curve()
            .state_at_specific_enthalpy(h)
            .map_err(text)?;
        let point = ThermalMaterialPoint {
            material: self.material,
            context: self.context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        };
        point.at(phase).map_err(text)
    }

    pub(super) fn values(
        &self,
        kelvin: f64,
        wanted: &[(&PropertyKey, ScalarAdmissibility)],
    ) -> Result<Vec<f64>, String> {
        let query = self.query(kelvin)?;
        let requirements = wanted
            .iter()
            .map(|(key, rule)| ScalarPropertyRequirement::try_with_key(key, *rule))
            .collect::<Result<Vec<_>, _>>()
            .map_err(text)?;
        let snapshot = resolve_material_state_point(
            self.material,
            &query,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .map_err(text)?;
        wanted
            .iter()
            .map(|(key, _)| {
                snapshot
                    .property_by_key(key)
                    .map(ResolvedScalarProperty::value_si)
                    .ok_or_else(|| format!("the card has no {key:?}"))
            })
            .collect()
    }
}

/// A resistance from its receiving cell's card, as the module's channels
/// read it.
#[derive(Debug)]
pub(super) struct CardResistance {
    pub(super) card: Card,
    pub(super) key: &'static PropertyKey,
    /// Names the law in the drive's identity.
    pub(super) tag: u64,
}

impl ResistanceLaw for CardResistance {
    fn ohms(&self, kelvin: f64) -> Result<f64, String> {
        let wanted = [(self.key, ScalarAdmissibility::StrictlyPositive)];
        Ok(self.card.values(kelvin, &wanted)?[0])
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.tag);
    }
}

/// A switch from its package's card, with the module's admissibility (the
/// diode threshold may be zero). The package's storage claim, if any, must
/// admit the temperature too; its value is the drive's fixed storage.
#[derive(Debug)]
struct CardSwitch {
    card: Card,
    keys: &'static [PropertyKey; 7],
    storage: Option<&'static PropertyKey>,
    tag: u64,
}

impl SwitchLaw for CardSwitch {
    fn switch(&self, kelvin: f64) -> Result<SwitchDevice, String> {
        use ScalarAdmissibility::{Finite, NonNegative, StrictlyPositive};
        let mut wanted: Vec<_> = self
            .keys
            .iter()
            .enumerate()
            .map(|(role, key)| {
                (
                    key,
                    if role == 2 {
                        NonNegative
                    } else {
                        StrictlyPositive
                    },
                )
            })
            .collect();
        if let Some(storage) = self.storage {
            wanted.push((storage, Finite));
        }
        let v = self.card.values(kelvin, &wanted)?;
        Ok(SwitchDevice {
            on_resistance_ohm: v[0],
            off_resistance_ohm: v[1],
            body_diode: BodyDiode {
                forward_voltage_v: v[2],
                differential_resistance_ohm: v[3],
            },
            bounds: DeviceOperatingBounds {
                maximum_reverse_voltage_v: v[4],
                maximum_channel_current_a: v[5],
                maximum_diode_current_a: v[6],
            },
        })
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.tag);
    }
}

/// A caloric curve as the implicit timeline's law. No conductive link joins
/// the cells, so the conductivity is never used.
struct Caloric(&'static ResolvedHeatCapacityCurve);

impl ConductivityLaw for Caloric {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(1.0)
    }
}

impl ImplicitConductionLaw for Caloric {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.0.curve()
    }

    fn conductivity_law(&self) -> &dyn ConductivityLaw {
        self
    }

    fn enthalpy_bounds_j_kg(&self) -> [f64; 2] {
        let knots = self.0.curve().knots();
        [
            knots[0].specific_enthalpy_j_kg,
            knots[knots.len() - 1].specific_enthalpy_j_kg,
        ]
    }

    fn state_at_specific_enthalpy(
        &self,
        h: f64,
    ) -> Result<ImplicitConductionState, DepositionError> {
        let knots = self.0.curve().knots();
        let (first, last) = (&knots[0], &knots[knots.len() - 1]);
        // The fixture's specific heat is constant: one slope dT/dh.
        Ok(ImplicitConductionState {
            temperature_k: self.0.curve().temperature_k(h)?,
            temperature_slope_kg_k_j: (last.temperature_k - first.temperature_k)
                / (last.specific_enthalpy_j_kg - first.specific_enthalpy_j_kg),
            conductivity_w_m_k: 1.0,
            conductivity_slope_kg_w_m_k_j: 0.0,
        })
    }
}

/// The graph's stand-in for the module's `FixedCooling`, local to these
/// tests and not a material law: the link at the fixed conductance the
/// heat-assembly fixtures give it, between cells of the two caloric laws.
/// The timeline takes `Q = G (T_a - T_b)`; the law gives its derivatives in
/// both cells' enthalpies.
struct FixedCoolingStandIn([&'static Caloric; 2], f64);

impl ImplicitLinkLaw for FixedCoolingStandIn {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.0.map(|law| law.0.curve() as &dyn ThermalLaw)
    }

    fn state_at_specific_enthalpies(
        &self,
        enthalpies_j_kg: [f64; 2],
        _: ConductiveLink,
        _: ContactPatch,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitLinkState, DepositionError> {
        checkpoint()?;
        let [a, b] = enthalpies_j_kg;
        let a = ImplicitConductionLaw::state_at_specific_enthalpy(self.0[0], a)?;
        let b = ImplicitConductionLaw::state_at_specific_enthalpy(self.0[1], b)?;
        let g = self.1;
        Ok(ImplicitLinkState {
            conductance_w_k: g,
            power_enthalpy_derivatives_w_kg_j: [
                g * a.temperature_slope_kg_k_j,
                -g * b.temperature_slope_kg_k_j,
            ],
        })
    }
}

/// The cells at the start: the seven at 300 K, or the drivers' `start`.
/// With a link they are born again as the module's heat-assembly fixtures
/// bear them, with the coolant if any and one contact.
pub(super) fn cells_start<L: Lane>(case: &Case<'_, L>) -> DepositionSimulation<'static> {
    let f = case.fixture;
    let born = case
        .drivers
        .map_or_else(|| f.thermal(), |drivers| drivers.start.clone());
    let Some((id, [a, b], _)) = link(case) else {
        return born;
    };
    let birth = |cell: &DepositedCell| TimedEvent {
        time_s: 0.0,
        event: Event::Birth(cell.clone()),
    };
    let mut events: Vec<_> = born.cells().values().map(birth).collect();
    events.extend(case.coolant.map(|coolant| birth(&coolant.cell)));
    let cells = events.len();
    let link = ThermalLink {
        id,
        a,
        b,
        area_m2: 1.0,
        conductance_w_k: 0.0,
    };
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(link),
    });
    let budget = DepositionBudget {
        max_cells: cells,
        max_links: 1,
        max_events: cells + 1,
        max_substeps: 1000,
    };
    let config = StepConfig {
        max_dt_s: 0.1,
        energy_tolerance_j: 1e-9,
    };
    let mut laws = vec![(MaterialId(0), f.caloric.curve() as &dyn ThermalLaw)];
    laws.extend(
        case.coolant
            .map(|coolant| (MaterialId(1), coolant.caloric.curve() as &dyn ThermalLaw)),
    );
    let cells = DepositionSimulation::new(events, &laws, budget).unwrap();
    cells.advance_to(0.0, config, || Ok(())).unwrap().0
}

/// Runs `body` with a context and the module's heat assembly around the
/// case's cells, the drivers' and the coolant's included, the link at the
/// conductance its heat-assembly fixtures give it (`FixedCooling`).
pub(super) fn heat_assembly<L: Lane, R>(
    case: &Case<'_, L>,
    start: &DepositionSimulation<'_>,
    body: impl FnOnce(&Cx<'_>, &FixedCooling<'_, '_, '_>) -> R,
) -> R {
    let f = case.fixture;
    let mut cells: Vec<_> = f.cells(case.context).into_iter().collect();
    let drivers = case.drivers.into_iter().flat_map(|drivers| drivers.cells);
    cells.extend(drivers.map(|cell| MaterialCell {
        cell,
        material_id: MaterialId(0),
        curve: f.caloric.curve(),
        point: ThermalMaterialPoint {
            material: &f.card,
            context: case.context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }));
    cells.extend(case.coolant.map(|coolant| MaterialCell {
        cell: coolant.cell.id,
        material_id: MaterialId(1),
        curve: coolant.caloric.curve(),
        point: ThermalMaterialPoint {
            material: &coolant.card,
            context: &coolant.context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }));
    let bulk = MaterialNetwork {
        cells: &cells,
        links: &[],
        resistors: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    with_cx(|cx| {
        let unset = BTreeMap::from([(LINK.0, 0.0)]);
        let prepared = heat
            .prepare_with_link_conductances(cx, start, &unset)
            .unwrap();
        let transport = FixedCooling {
            prepared: &prepared,
            conductance_w_k: link(case).map_or(LINK.2, |(_, _, g)| g),
        };
        body(cx, &transport)
    })
}

/// The module's cells as one backward-Euler timeline whose ports, one for
/// each of `port_cells`, return end-of-step temperatures, and its start. Its
/// longest interval is the cells' own (`Fixture::thermal_at`), so it does
/// not cap the step.
fn timeline<L: Lane>(case: &Case<'_, L>) -> (ThermalTimeline<'static>, TimelineState<'static>) {
    let f = case.fixture;
    let start = cells_start(case);
    let law: &'static Caloric = leak(Caloric(&f.caloric));
    let fluid = case
        .coolant
        .map(|coolant| -> &'static Caloric { leak(Caloric(&coolant.caloric)) });
    let mut conductivity = vec![(MaterialId(0), law as &dyn ConductivityLaw)];
    let mut caloric = vec![(MaterialId(0), law as &dyn ImplicitConductionLaw)];
    if let Some(fluid) = fluid {
        conductivity.push((MaterialId(1), fluid));
        caloric.push((MaterialId(1), fluid));
    }
    let linked = link(case).is_some();
    // The link's geometry is bound as the observer requires; the stand-in
    // replaces the conduction it would give.
    let geometry = [ConductiveLink {
        link_id: LINK.0,
        half_distances_m: [0.5; 2],
        contact_resistance_m2_k_w: 0.0,
    }];
    let observer = ConductionObserver::new(
        &start,
        &conductivity,
        &geometry[..usize::from(linked)],
        &[],
        || Ok(()),
    )
    .unwrap();
    let conduction = ImplicitConductionGraph::new(
        &start,
        &observer,
        &caloric,
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let conduction = if linked {
        let laws = [law, fluid.unwrap_or(law)];
        let g = link(case).map_or(LINK.2, |(_, _, g)| g);
        let link: &'static FixedCoolingStandIn = leak(FixedCoolingStandIn(laws, g));
        let links = [(LINK.0, link as &dyn ImplicitLinkLaw)];
        conduction
            .with_link_laws(&start, &links, || Ok(()))
            .unwrap()
    } else {
        conduction
    };
    let controls = ImplicitConductionControls {
        maximum_dt_s: 1.0,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-9,
        ..Default::default()
    };
    let timeline = ThermalTimeline::implicit(
        port_cells(case),
        conduction,
        controls,
        ImplicitScheme::BackwardEuler,
    )
    .unwrap()
    .returning_end_temperatures()
    .unwrap();
    let start = timeline.start(&start).unwrap();
    (timeline, start)
}

/// The native drive of `steady_drive` running the case's gate program,
/// its resistances and switches following the fixture's card.
fn drive<L: Lane>(case: &Case<'_, L>) -> CommutatedBatteryMotorDrive<L> {
    let f = case.fixture;
    let program = |k| GateProgram {
        intervals: vec![CommutationInterval {
            gates: case.gates(k),
            duration_s: case.dt_s,
        }],
        auxiliary_loads: vec![case.steady.command.auxiliary_load],
        substeps_per_interval: case.substeps,
    };
    let drive_parameters = DriveParameters {
        capacity_c: 100.0,
        soc_voltage: vec![(0.0, 12.0), (1.0, 12.0)],
        battery: BatteryParameters {
            series_resistance_ohm: 0.1,
            polarization_resistance_ohm: 0.1,
            polarization_capacitance_f: 2.0,
            bus_capacitance_f: 0.5,
        },
        lead_resistance_ohm: 0.1,
        converter: parameters(),
        motor: motor(),
    };
    let mut drive = CommutatedBatteryMotorDrive::new(
        drive_parameters,
        case.storage.clone(),
        program(0),
        power_controls(),
        cycle_controls(),
        1e-9,
    )
    .unwrap();
    // A schedule runs step k's program as the drive's program k.
    if case.schedule.len() > 1 {
        let programs = (1..case.schedule.len()).map(program).collect();
        drive = drive.scheduled(programs).unwrap();
    }
    let card = Card {
        material: &f.card,
        curve: &f.caloric,
        context: case.context,
    };
    let resistors = [
        DriveResistor::ConverterWinding,
        DriveResistor::BatterySeries,
        DriveResistor::BatteryPolarization,
        DriveResistor::Lead,
        DriveResistor::MotorWinding,
    ];
    for (tag, resistor) in (0..).zip(resistors) {
        // `drive_resistors` reads every channel from the first key.
        let law = CardResistance {
            card,
            key: &f.keys[0],
            tag,
        };
        drive = drive.following(resistor, Arc::new(law)).unwrap();
    }
    let switches = [
        SwitchPosition::InputHigh,
        SwitchPosition::InputLow,
        SwitchPosition::OutputHigh,
        SwitchPosition::OutputLow,
    ];
    let mut storage = case.storage_keys.iter();
    for (tag, position) in (0..).zip(switches) {
        let law = CardSwitch {
            card,
            keys: &f.keys,
            storage: storage.next(),
            tag,
        };
        drive = drive.switch_following(position, Arc::new(law)).unwrap();
    }
    L::add_drivers(drive, case)
}

/// The case as a graph. Returns the system, the drive and the cells.
pub(super) fn declare<L: Lane>(case: &Case<'_, L>) -> (System, NodeId, NodeId) {
    let mut builder = SystemBuilder::new();
    let node = builder.add("drive", drive(case), case.start).unwrap();
    let torques = case.steady.command.shaft_torque_nm;
    let shafts = ["rotor", "stator"].into_iter().zip(torques);
    let shafts =
        shafts.map(|(name, torque)| (name, torque, PortKind::RotationalTorqueAngularVelocity));
    let rails = case.drivers.map_or(&[][..], |drivers| &drivers.rails[..]);
    let rails = RAIL_PORTS.into_iter().zip(rails.iter().copied());
    let rails = rails.map(|(name, volts)| (name, volts, PortKind::ElectricalVoltageCurrent));
    for (port, effort, kind) in shafts.chain(rails) {
        let source = Force::in_domain(kind, Waveform::constant(effort).unwrap());
        let source = builder.add(port, source, SourceAccount::default()).unwrap();
        let (a, b) = (
            builder.port(source, "p").unwrap(),
            builder.port(node, port).unwrap(),
        );
        builder.bond(a, b).unwrap();
    }
    if case.schedule.len() > 1 {
        // Program k from the start of step k, held over the step.
        let steps = u32::try_from(case.schedule.len()).unwrap();
        let table = (0..steps)
            .map(|k| (case.dt_s * f64::from(k), f64::from(k)))
            .collect();
        let source = SignalSource::new(vec![Waveform::held_table(table).unwrap()]).unwrap();
        let source = builder
            .add_with_tangent("schedule", source, SignalClock::default())
            .unwrap();
        let (a, b) = (
            builder.port(source, "output").unwrap(),
            builder.port(node, "program").unwrap(),
        );
        builder.signal(a, b).unwrap();
    }
    let (timeline, start) = timeline(case);
    let cells = builder.add("cells", timeline, start).unwrap();
    let mut routes: Vec<_> = HEAT_PORTS.into_iter().zip(RECEIVERS.map(Some)).collect();
    // The drive has a load port only when its program carries a load.
    let load = case.steady.command.auxiliary_load;
    routes.extend(load.map(|_| (LOAD_PORT, case.load_cell)));
    let drivers = case.drivers.into_iter().flat_map(|drivers| drivers.cells);
    routes.extend(DRIVER_HEAT_PORTS.into_iter().zip(drivers.map(Some)));
    for &(name, receiver) in &routes {
        let port = builder.port(node, name).unwrap();
        match receiver {
            Some(cell) => {
                let cell = builder.port(cells, &format!("cell{}", cell.0)).unwrap();
                builder.heat_route(port, &[(cell, 1.0)], None, RouteTemperature::SinkEquivalent)
            }
            None => builder.heat_route(
                port,
                &[],
                Some(heat_sink("outside", 300.0).unwrap()),
                RouteTemperature::Fixed(300.0),
            ),
        }
        .unwrap();
    }
    let environment = builder.port(cells, "environment").unwrap();
    builder
        .heat_route(
            environment,
            &[],
            Some(heat_sink("environment", 300.0).unwrap()),
            RouteTemperature::Fixed(300.0),
        )
        .unwrap();
    let mut system = builder.seal().unwrap();
    for &(name, _) in &routes {
        let port = system.port(node, name).unwrap();
        assert!(system.guess(port, &[300.0]));
    }
    (system, node, cells)
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let key = fs_exec::StreamKey {
            seed: 37,
            kernel_id: 66,
            tile: 0,
            iteration: 0,
        };
        let budget = fs_exec::Budget::INFINITE;
        f(&Cx::new(
            &gate,
            arena,
            key,
            budget,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

/// The module's own check of its heat-assembly fixture, on the graph's first
/// step: backward Euler through the link between cells of 1 J/K from 300 K,
/// `X = h G (Q0 - Q6) / (1 + 2 h G)` with the steady heats. The largest
/// difference over the cells \[K\].
pub(super) fn closed_form_step(case: &Case<'_>) -> f64 {
    let (h, g) = (case.dt_s, LINK.2);
    let q = per_body(case.steady.heat_w.map(|p| h * p), 0.0);
    let x = h * g * (q[0] - q[6]) / (1.0 + 2.0 * h * g);
    let exchange = [-x, 0.0, 0.0, 0.0, 0.0, 0.0, x];
    let kelvin = kernel_run::<DeviceStorage, [f64; 7]>(case, 1)[0].kelvin;
    let error = |i: usize| (kelvin[i] - (300.0 + q[i] + exchange[i])).abs();
    (0..7).map(error).fold(0.0, f64::max)
}

/// The owner's step built by hand for step `k` of the case, from `old` at the
/// cells' temperatures `kelvin`: the card in closed form, not a second
/// engine. As `fixture` authors the card, every resistance is 0.1 ohm, a
/// channel is off at 100 ohm, each body diode conducts from
/// 0.5 V + `slope_v_k` (T - 300 K) at its package's temperature through
/// 0.2 ohm, and every limit is 100 (`parameters()`). The battery's flat
/// 12 V curve and `motor()` are `steady_drive`'s. The new state and the
/// heat each cell receives \[J\].
fn owner_step<L: HandOwner>(
    case: &Case<'_, L>,
    k: usize,
    (old, kelvin): (&L::State, &[f64]),
    slope_v_k: f64,
) -> (L::State, [f64; 7]) {
    let p = parameters();
    let devices = core::array::from_fn(|i| {
        let mut device = p.devices[i];
        device.body_diode.forward_voltage_v = 0.5 + slope_v_k * (kelvin[i + 1] - 300.0);
        device
    });
    let intervals = [CommutationInterval {
        gates: case.gates(k),
        duration_s: case.dt_s,
    }];
    let (next, heat, load_work_j) = case.storage.cycle(
        owner_drive(),
        devices,
        &intervals,
        case.substeps,
        case.steady.command,
        old,
    );
    let routed = if case.load_cell.is_some() {
        load_work_j
    } else {
        0.0
    };
    (next, per_body(heat, routed))
}

/// The drive of `owner_step`'s owner: the converter's, battery's and lead's
/// resistances of 0.1 ohm, the battery's flat 12 V curve and `motor()`.
pub(super) fn owner_drive() -> BatteryMotorDrive {
    let p = parameters();
    let converter = BuckBoostParameters {
        inductance_h: p.inductance_h,
        bus_capacitance_f: p.bus_capacitance_f,
        winding_resistance_ohm: 0.1,
        input_on_resistance_ohm: [0.1; 2],
        output_on_resistance_ohm: [0.1; 2],
    };
    let bridge = BridgeModulation {
        input_high_fraction: 0.0,
        output_high_fraction: 0.0,
    };
    let battery = BatteryParameters {
        series_resistance_ohm: 0.1,
        polarization_resistance_ohm: 0.1,
        polarization_capacitance_f: 2.0,
        bus_capacitance_f: 0.5,
    };
    let curve = BatteryOcvCurve::new(100.0, &[(0.0, 12.0), (1.0, 12.0)]).unwrap();
    AveragedBuckBoost::new(converter, bridge)
        .and_then(|c| c.with_battery_source(BatterySource::new(curve, battery)?, 0.1))
        .unwrap()
        .with_dc_motor(DcMotor::new(motor()).unwrap())
        .unwrap()
}

/// Holds the largest differences of a reference check within `bounds`.
fn within(what: &str, worst: [f64; 2], bounds: [f64; 2]) {
    println!(
        "{what}: coordinates {:.1e} relative, cells {:.1e}",
        worst[0], worst[1]
    );
    assert!(
        worst[0] <= bounds[0] && worst[1] <= bounds[1],
        "{what}: {worst:?} beyond {bounds:?}"
    );
}

/// The graph's steps against the owner's step built by hand (`owner_step`);
/// see [`held_to_owner`].
pub(super) fn owner_check<L: HandOwner, K: AsRef<[f64]>>(
    case: &Case<'_, L>,
    rows: &[Row<L::State, K>],
    slope_v_k: f64,
    bounds: [f64; 2],
) {
    let owner = |k, old: &L::State, kelvin: &[f64]| {
        let (next, heat) = owner_step(case, k, (old, kelvin), slope_v_k);
        (next, heat.to_vec())
    };
    held_to_owner(case, rows, owner, bounds);
}

/// The graph's steps against `owner`, the owner's step k built by hand from
/// the graph's old state at its end temperatures (the new state and the
/// heat of each port cell), and each cell's backward-Euler balance: its
/// energy change is the heat routed to it and, where linked, the link's
/// `G h (T_a - T_b)` at the end temperatures. The coolant's energy is its
/// mass times its caloric curve's enthalpy. Bounds on the coordinates
/// (relative) and the cells' energies \[J\].
pub(super) fn held_to_owner<L: Lane, K: AsRef<[f64]>>(
    case: &Case<'_, L>,
    rows: &[Row<L::State, K>],
    owner: impl Fn(usize, &L::State, &[f64]) -> (L::State, Vec<f64>),
    bounds: [f64; 2],
) {
    let ids = cell_ids(case);
    let at = |cell| ids.iter().position(|&id| id == cell).unwrap();
    let mut before = match case.drivers {
        Some(drivers) => port_cells(case)
            .into_iter()
            .map(|id| drivers.start.temperature_k(id).unwrap())
            .collect(),
        None => vec![300.0; 7],
    };
    before.extend(case.coolant.map(|coolant| {
        let h = coolant.cell.specific_h_j_kg;
        coolant.caloric.curve().temperature_k(h).unwrap()
    }));
    let (mut old, mut worst) = (case.start, [0.0_f64; 2]);
    for (k, row) in rows.iter().enumerate() {
        let kelvin = row.kelvin.as_ref();
        let (next, mut heat) = owner(k, &old, kelvin);
        heat.resize(ids.len(), 0.0);
        if let Some((_, [a, b], g)) = link(case) {
            let flux = g * case.dt_s * (kelvin[at(a)] - kelvin[at(b)]);
            (heat[at(a)], heat[at(b)]) = (heat[at(a)] - flux, heat[at(b)] + flux);
        }
        for (a, b) in row.state.as_ref().iter().zip(next.as_ref()) {
            worst[0] = worst[0].max((a - b).abs() / a.abs().max(1.0));
        }
        for (i, q) in heat.into_iter().enumerate() {
            let change = match case.coolant.filter(|_| i + 1 == ids.len()) {
                Some(c) => {
                    let h = |kelvin| enthalpy_at(&c.caloric, kelvin).unwrap();
                    c.cell.mass_kg * (h(kelvin[i]) - h(before[i]))
                }
                None => case.fixture.cp * (kelvin[i] - before[i]),
            };
            worst[1] = worst[1].max((change - q).abs());
        }
        (old, before) = (row.state, kelvin.to_vec());
    }
    within("graph against the owner by hand", worst, bounds);
}

/// The steady fixtures' closed form (`steady_drive`), as the module checks
/// it: after k periods each cell of 1 J/K holds k times its steady heat
/// (and the routed load), and the drive is at its start but for the
/// battery's charge, k dt `cell_current_a` lower. Bounds on the coordinates
/// (relative) and the temperatures \[K\].
pub(super) fn steady_check<L: Lane, K: AsRef<[f64]>>(
    case: &Case<'_, L>,
    rows: &[Row<L::State, K>],
    bounds: [f64; 2],
) {
    let (s, dt) = (case.steady, case.dt_s);
    let load = s
        .command
        .auxiliary_load
        .map_or(0.0, |l| l.requested_power_w());
    let routed = if case.load_cell.is_some() {
        load * dt
    } else {
        0.0
    };
    let q = per_body(s.heat_w.map(|p| p * dt), routed);
    let mut worst = [0.0_f64; 2];
    for (k, row) in (1..).map(f64::from).zip(rows) {
        let mut state = case.start;
        state.as_mut()[0] -= k * dt * s.cell_current_a;
        for (a, b) in row.state.as_ref().iter().zip(state.as_ref()) {
            worst[0] = worst[0].max((a - b).abs() / a.abs().max(1.0));
        }
        for (kelvin, q) in row.kelvin.as_ref().iter().zip(q) {
            worst[1] = worst[1].max((kelvin - 300.0 - k * q).abs());
        }
    }
    within("graph against the closed form", worst, bounds);
}

/// `steps` steps of the case's graph. The kernel's ledger closes.
pub(super) fn kernel_run<L: Lane, K: Default + AsMut<[f64]>>(
    case: &Case<'_, L>,
    steps: usize,
) -> Vec<Row<L::State, K>> {
    let ids = cell_ids(case);
    let (mut system, drive, cells) = declare(case);
    let rows = with_cx(|cx| {
        (0..steps)
            .map(|_| {
                let report = system.step(case.dt_s, cx).unwrap();
                assert_eq!(report.dt_s.to_bits(), case.dt_s.to_bits());
                let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
                let mut kelvin = K::default();
                assert_eq!(kelvin.as_mut().len(), ids.len());
                for (t, &id) in kelvin.as_mut().iter_mut().zip(&ids) {
                    *t = timeline.simulation().temperature_k(id).unwrap();
                }
                Row {
                    state: *system.state::<L::State>(drive).unwrap(),
                    kelvin,
                }
            })
            .collect()
    });
    let residual = system.ledger().energy_residual_j();
    assert!(residual.abs() <= 1e-9, "ledger residual {residual}");
    rows
}

/// One step of `dt_s` on the case's graph, which must be refused with the
/// drive and the cells as they were. Returns the refusal's text.
pub(super) fn kernel_refusal<L: Lane>(case: &Case<'_, L>, dt_s: f64) -> String {
    let (mut system, drive, cells) = declare(case);
    let bits = |system: &System| {
        let state = system.state::<L::State>(drive).unwrap();
        state
            .as_ref()
            .iter()
            .map(|q| q.to_bits())
            .collect::<Vec<_>>()
    };
    let state = bits(&system);
    let timeline = system
        .state::<TimelineState<'static>>(cells)
        .unwrap()
        .clone();
    let error = with_cx(|cx| system.step(dt_s, cx)).unwrap_err();
    assert_eq!(bits(&system), state);
    assert_eq!(
        system.state::<TimelineState<'static>>(cells).unwrap(),
        &timeline
    );
    text(error)
}

/// Holds every step of the kernel's run to the module's: each coordinate
/// within `bounds[0]` of its size (at least 1), each temperature within
/// `bounds[1]` K. Prints the largest differences.
pub(super) fn compare<S: AsRef<[f64]>, K: AsRef<[f64]>>(
    module: &[Row<S, K>],
    kernel: &[Row<S, K>],
    bounds: [f64; 2],
) {
    assert_eq!(module.len(), kernel.len());
    let mut worst = [0.0_f64; 2];
    for (k, (theirs, ours)) in module.iter().zip(kernel).enumerate() {
        let (a, b) = (ours.state.as_ref(), theirs.state.as_ref());
        for (i, (a, b)) in a.iter().zip(b).enumerate() {
            let relative = (a - b).abs() / a.abs().max(1.0);
            worst[0] = worst[0].max(relative);
            assert!(
                relative <= bounds[0],
                "step {k}, coordinate {i}: {a} vs {b}"
            );
        }
        let (a, b) = (ours.kelvin.as_ref(), theirs.kelvin.as_ref());
        for (i, (a, b)) in a.iter().zip(b).enumerate() {
            worst[1] = worst[1].max((a - b).abs());
            assert!(
                (a - b).abs() <= bounds[1],
                "step {k}, cell {i}: {a} K vs {b} K"
            );
        }
    }
    println!(
        "graph against module: coordinates {:.1e} relative, temperatures {:.1e} K",
        worst[0], worst[1]
    );
}
