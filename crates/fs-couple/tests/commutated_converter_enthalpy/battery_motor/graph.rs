//! The battery-motor drive's storage fixtures as declared graphs on the
//! composition kernel.
//!
//! `battery_motor::{capacitance, charge}` couple fs-phs's commutated drive to
//! the enthalpy timeline by hand. A step runs the gate program once, each
//! resistance and switch at its receiving cell's card at the END
//! temperature, with Newton on the cell energies. The graph says the same
//! with no code for the pair:
//! - `fs_system::commutated_drive::CommutatedBatteryMotorDrive` with the
//!   fixture's storage, its five resistances following [`CardResistance`]
//!   and its four switches following [`CardSwitch`];
//! - two `elements::Force` sources holding the command's shaft torques;
//! - the seven cells as one backward-Euler timeline whose cell ports return
//!   end-of-step temperatures;
//! - one heat route per heat port into the module's receiving cell.
//!
//! The two sides solve the same equations with different solvers, so they
//! agree to within the solves' tolerances.

use std::fmt::Debug;
use std::sync::Arc;

use fs_material::state_point::{
    ResolvedScalarProperty, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_ports::PortKind;
use fs_system::circuit::ResistanceLaw;
use fs_system::commutated_drive::{
    CommutatedBatteryMotorDrive, DeviceStorage, DriveParameters, DriveResistor, GateProgram,
    HEAT_PORTS, LOAD_PORT, SwitchLaw, SwitchPosition,
};
use fs_system::elements::{Force, SourceAccount, Waveform};
use fs_system::{NodeId, RouteTemperature, System, SystemBuilder, heat_sink};
use fs_thermal::component::{ImplicitScheme, ThermalTimeline, TimelineState};
use fs_thermal::enthalpy::conduction::implicit::{
    ImplicitConductionControls, ImplicitConductionGraph, ImplicitConductionLaw,
    ImplicitConductionLinearControls, ImplicitConductionState,
};
use fs_thermal::enthalpy::conduction::{ConductionObserver, ConductivityLaw};

use super::*;

/// The cell each heat port heats, in `HEAT_PORTS` order: the module's
/// channels (`drive_resistors`), device packages (`Fixture::devices`) and
/// bearing cell, the grouping of `per_body`.
const RECEIVERS: [CellId; 10] = [
    IDS[0], IDS[1], IDS[2], IDS[3], IDS[4], IDS[5], IDS[6], IDS[5], IDS[6], IDS[0],
];

/// One fixture of a storage file, for both sides.
pub(super) struct Case<'s> {
    pub(super) fixture: &'static Fixture,
    /// The query context of the module's cells.
    pub(super) context: &'static QueryPoint,
    /// Each package's storage claim, in switch order.
    pub(super) storage_keys: &'static [PropertyKey; 4],
    pub(super) storage: DeviceStorage,
    pub(super) steady: &'s SteadyDrive,
    pub(super) start: [f64; 11],
    pub(super) dt_s: f64,
    pub(super) substeps: usize,
    /// The cell that receives the auxiliary load's work; none leaves it
    /// outside.
    pub(super) load_cell: Option<CellId>,
}

/// One accepted step: the drive's eleven coordinates and the seven cells'
/// temperatures \[K\].
#[derive(Debug, Clone, Copy)]
pub(super) struct Row {
    pub(super) state: [f64; 11],
    pub(super) kelvin: [f64; 7],
}

pub(super) fn leak<T>(value: T) -> &'static T {
    Box::leak(Box::new(value))
}

fn text(error: impl Debug) -> String {
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
struct Card {
    material: &'static MaterialCard,
    curve: &'static ResolvedHeatCapacityCurve,
    context: &'static QueryPoint,
}

impl Card {
    fn values(
        &self,
        kelvin: f64,
        wanted: &[(&PropertyKey, ScalarAdmissibility)],
    ) -> Result<Vec<f64>, String> {
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
        let query = point.at(phase).map_err(text)?;
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
struct CardResistance {
    card: Card,
    key: &'static PropertyKey,
    /// Names the law in the drive's identity.
    tag: u64,
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
/// diode threshold may be zero). The package's storage claim must admit the
/// temperature too; its value is the drive's fixed storage.
#[derive(Debug)]
struct CardSwitch {
    card: Card,
    keys: &'static [PropertyKey; 7],
    storage: &'static PropertyKey,
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
        wanted.push((self.storage, Finite));
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

/// The fixture's caloric curve as the implicit timeline's law. No link
/// joins the cells, so the conductivity is never used.
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

/// The module's seven cells as one backward-Euler timeline whose ports
/// return end-of-step temperatures, and its start. Its longest interval is
/// the cells' own (`Fixture::thermal_at`), so it does not cap the step.
fn timeline(f: &'static Fixture) -> (ThermalTimeline<'static>, TimelineState<'static>) {
    let start = f.thermal();
    let law: &'static Caloric = leak(Caloric(&f.caloric));
    let observer = ConductionObserver::new(
        &start,
        &[(MaterialId(0), law as &dyn ConductivityLaw)],
        &[],
        &[],
        || Ok(()),
    )
    .unwrap();
    let conduction = ImplicitConductionGraph::new(
        &start,
        &observer,
        &[(MaterialId(0), law as &dyn ImplicitConductionLaw)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let controls = ImplicitConductionControls {
        maximum_dt_s: 1.0,
        cell_tolerance_j: 1e-10,
        total_tolerance_j: 1e-9,
        ..Default::default()
    };
    let timeline = ThermalTimeline::implicit(
        IDS.to_vec(),
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
fn drive(case: &Case<'_>) -> CommutatedBatteryMotorDrive {
    let f = case.fixture;
    let program = GateProgram {
        intervals: vec![CommutationInterval {
            gates: case.steady.gates,
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
        program,
        power_controls(),
        cycle_controls(),
        1e-9,
    )
    .unwrap();
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
    for ((tag, position), storage) in (0..).zip(switches).zip(case.storage_keys) {
        let law = CardSwitch {
            card,
            keys: &f.keys,
            storage,
            tag,
        };
        drive = drive.switch_following(position, Arc::new(law)).unwrap();
    }

    drive
}

/// The case as a graph. Returns the system, the drive and the cells.
pub(super) fn declare(case: &Case<'_>) -> (System, NodeId, NodeId) {
    let f = case.fixture;
    let mut builder = SystemBuilder::new();
    let node = builder.add("drive", drive(case), case.start).unwrap();
    let torques = case.steady.command.shaft_torque_nm;
    for (shaft, torque) in ["rotor", "stator"].into_iter().zip(torques) {
        let kind = PortKind::RotationalTorqueAngularVelocity;
        let source = Force::in_domain(kind, Waveform::constant(torque).unwrap());
        let source = builder
            .add(shaft, source, SourceAccount::default())
            .unwrap();
        let (a, b) = (
            builder.port(source, "p").unwrap(),
            builder.port(node, shaft).unwrap(),
        );
        builder.bond(a, b).unwrap();
    }
    let (timeline, start) = timeline(f);
    let cells = builder.add("cells", timeline, start).unwrap();
    let mut routes: Vec<_> = HEAT_PORTS.into_iter().zip(RECEIVERS.map(Some)).collect();
    routes.push((LOAD_PORT, case.load_cell));
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

/// `steps` steps of the case's graph. The kernel's ledger closes.
pub(super) fn kernel_run(case: &Case<'_>, steps: usize) -> Vec<Row> {
    let (mut system, drive, cells) = declare(case);
    let rows = with_cx(|cx| {
        (0..steps)
            .map(|_| {
                let report = system.step(case.dt_s, cx).unwrap();
                assert_eq!(report.dt_s.to_bits(), case.dt_s.to_bits());
                let timeline = system.state::<TimelineState<'static>>(cells).unwrap();
                Row {
                    state: *system.state::<[f64; 11]>(drive).unwrap(),
                    kelvin: IDS.map(|id| timeline.simulation().temperature_k(id).unwrap()),
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
pub(super) fn kernel_refusal(case: &Case<'_>, dt_s: f64) -> String {
    let (mut system, drive, cells) = declare(case);
    let bits = |system: &System| system.state::<[f64; 11]>(drive).unwrap().map(f64::to_bits);
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
pub(super) fn compare(module: &[Row], kernel: &[Row], bounds: [f64; 2]) {
    assert_eq!(module.len(), kernel.len());
    let mut worst = [0.0_f64; 2];
    for (k, (theirs, ours)) in module.iter().zip(kernel).enumerate() {
        for (i, (a, b)) in ours.state.iter().zip(theirs.state).enumerate() {
            let relative = (a - b).abs() / a.abs().max(1.0);
            worst[0] = worst[0].max(relative);
            assert!(
                relative <= bounds[0],
                "step {k}, coordinate {i}: {a} vs {b}"
            );
        }
        for (i, (a, b)) in ours.kelvin.iter().zip(theirs.kelvin).enumerate() {
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
