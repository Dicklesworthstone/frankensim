//! Complete motor heat transport, independent balances and atomic refusals.
//! Synthetic fixtures establish discrete integration, not measured motor validity.
use super::*;
use fs_conduction::transient::deposition::contact::{ContactCell, MaterialContactLink};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
};
use fs_couple::circuit_enthalpy::{
    CircuitThermalCoefficients,
    contacts::{
        ContactThermalCoefficients,
        heat::{HeatAssembly, HeatNetworkError, HeatNetworkSample, MaterialHeatNetwork, PreparedHeatNetwork},
        material::ElectricalContactNetwork,
    },
    material::MaterialNetwork,
    motor::heat::{MotorHeatError, MotorHeatNetwork},
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim,
    PropertyValue, Provenance, SelectionPolicy, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_mbd::rotordynamics::transient::drive::voltage::{
    RotorVoltageControl, RotorWindingMotor, RotorWindingState,
};
use fs_qty::{QuantitySpec, Temperature};

const CP: [f64; 2] = [1.0, 0.5];
const INDUCTANCE: f64 = 0.2;

fn voltage_control() -> RotorVoltageControl {
    RotorVoltageControl {
        source: power(1000.0),
        absolute_current_tolerance_a: 1e-10,
    }
}

fn with_case<R>(
    resistance: f64,
    extra_contact: bool,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &MaterialHeatNetwork<'_>, DepositionSimulation<'_>) -> R,
) -> R {
    let (cards, curves) = sources(CP, ALPHA);
    let context = common::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, Some(&key));
    let links = [ConductionLink { link: 7, distance_a_m: 0.5, distance_b_m: 0.5 }];
    let bulk = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk, contacts: &[], electrical_bases: &[],
    };
    let temperature = || QuantitySpec::dimensional(Temperature::DIMS);
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(RP, RD),
        value: PropertyValue::Curve {
            abscissa: "aT".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, resistance), (2000.0, resistance)],
            dims: RD,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("aT", temperature(), 200.0, 2000.0)
            .with_quantity("bT", temperature(), 200.0, 2000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic winding/bridge thermal contact".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }).unwrap();
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec { material: cards[0].id().clone(), texture_frame: "winding".into() },
        SurfaceSpec { material: cards[1].id().clone(), texture_frame: "bridge".into() },
        SystemContext {
            medium: "fixed interface".into(), third_body: None,
            environment: "authored".into(), history: "fixed".into(),
        },
        claims, vec![],
    ).unwrap();
    let contact_context = QueryPoint::new()
        .with_quantity("aT", temperature(), 300.0).unwrap()
        .with_quantity("bT", temperature(), 300.0).unwrap();
    let contacts = [MaterialContactLink {
        link: 7,
        sides: [0, 1].map(|i| ContactCell {
            cell: [WINDING, BRIDGE][i], material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
        }),
        point: ThermalInterfacePoint {
            interface: &interface, materials: [&cards[0], &cards[1]],
            context: &contact_context, temperature_axes: ["aT", "bT"],
            liquid_fraction_axes: [None, None], mean_temperature_axis: None,
        },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: if extra_contact { &contacts } else { &[] },
        radiation: None,
    };
    let thermal = body(&curves, CP, true, 0.0);
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            fs_exec::StreamKey { seed: 81, kernel_id: 53, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &assembly, thermal)
    })
}

// Solve the two caloric equations after eliminating the native free-inertia
// shaft work. The contact is in series with both half-cell Fourier resistances.
fn expected_current(resistance: f64) -> [f64; 2] {
    let g = DT / (1.0 + resistance);
    let q0 = DT * 4.0 * R0;
    let q1 = q0 * ALPHA;
    let bridge_ratio = 1.0 / drive().motoring_efficiency - 1.0;
    let b0 = bridge_ratio * (q0 + DT * DT);
    let b1 = bridge_ratio * q1;
    let a = CP[0] + g - q1;
    let d = CP[1] + g;
    let determinant = a * d - g * (g + b1);
    [(q0 * d + g * b0) / determinant, (a * b0 + (g + b1) * q0) / determinant]
}

// Independent scalar midpoint-current/endpoint-temperature reduction, with
// magnetic storage included in terminal work before bridge efficiency loss.
fn expected_voltage(resistance: f64) -> ([f64; 2], f64) {
    let g = DT / (1.0 + resistance);
    let rates = |winding_rise: f64| {
        let r = R0 * (1.0 + ALPHA * winding_rise);
        let current = DT * 12.0 / (2.0 * INDUCTANCE + DT * r + DT * DT / 4.0);
        let winding = DT * r * current * current;
        let bridge = DT * 12.0 * current * (1.0 / drive().motoring_efficiency - 1.0);
        let bridge_rise = (bridge + g * winding_rise) / (CP[1] + g);
        (current, bridge_rise, (CP[0] + g) * winding_rise - g * bridge_rise - winding)
    };
    let (mut low, mut high) = (0.0, 10.0);
    for _ in 0..80 {
        let middle = 0.5 * (low + high);
        if rates(middle).2 > 0.0 { high = middle; } else { low = middle; }
    }
    let rise = 0.5 * (low + high);
    let (current, bridge, _) = rates(rise);
    ([rise, bridge], current)
}

#[test]
fn contact_reaches_current_and_voltage_winding_feedback_and_both_heat_receivers() {
    let mut endpoints = Vec::new();
    for resistance in [0.5, 8.0] {
        endpoints.push(with_case(resistance, true, |cx, _, assembly, thermal| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let network = MotorHeatNetwork { assembly: &prepared, winding: winding(), bridge_cell: BRIDGE };
            let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
            let native = rotor.at_speed(0.0, || false).unwrap()
                .midpoint_stepper(DT, 1e-12, || false).unwrap();
            let motor = network.bind(cx, &thermal, &native, state(0.0), drive()).unwrap();
            let next = network.step_candidate(
                cx, &thermal, DT, &motor, &[0.0], 2.0, power(1000.0), controls(), &[0.0; 2],
            ).unwrap();
            let expected = expected_current(resistance);
            for (id, rise) in [WINDING, BRIDGE].into_iter().zip(expected) {
                near(next.thermal().temperature_k(id).unwrap(), 300.0 + rise);
            }
            near(next.motor().parameters().winding_resistance_ohm, R0 * (1.0 + ALPHA * expected[0]));
            near(next.motor().state().velocity[0], DT);
            near(next.report().source.network.thermal_contacts[0].conductance_w_k, 1.0 / (1.0 + resistance));
            let step = &next.report().electrical.drive;
            near(next.thermal().energy_j() - thermal.energy_j()
                + step.motion.kinetic_energy_change_j, step.electrical_draw_j);
            near(next.report().total_energy_residual_j, 0.0);

            let voltage_native = native.voltage_motor(
                RotorWindingMotor { drive: drive(), inductance_h: INDUCTANCE }, || false,
            ).unwrap();
            let voltage_motor = network.bind_voltage(cx, &thermal, &voltage_native,
                RotorWindingState { rotor: state(0.0), current_a: 0.0 }).unwrap();
            let voltage_next = network.step_voltage_candidate(
                cx, &thermal, DT, &voltage_motor, &[0.0], 12.0,
                voltage_control(), controls(), &[0.0; 2],
            ).unwrap();
            let (expected, current) = expected_voltage(resistance);
            for (id, rise) in [WINDING, BRIDGE].into_iter().zip(expected) {
                near(voltage_next.thermal().temperature_k(id).unwrap(), 300.0 + rise);
            }
            near(voltage_next.motor().state().current_a, 2.0 * current);
            near(voltage_next.motor().state().rotor.velocity[0], DT * current / 2.0);
            near(voltage_next.motor().parameters().drive.winding_resistance_ohm, R0 * (1.0 + ALPHA * expected[0]));
            near(voltage_next.report().source.network.thermal_contacts[0].conductance_w_k, 1.0 / (1.0 + resistance));
            let voltage_step = &voltage_next.report().electrical.drive;
            near(voltage_next.thermal().energy_j() - thermal.energy_j()
                + voltage_step.motion.kinetic_energy_change_j
                + voltage_step.magnetic_energy_change_j, voltage_step.electrical_draw_j);
            near(voltage_next.report().total_energy_residual_j, 0.0);
            [next.thermal().temperature_k(WINDING).unwrap(), voltage_next.motor().state().current_a]
        }));
    }
    for (low, high) in endpoints[0].iter().zip(&endpoints[1]) {
        assert!((low - high).abs() > 1e-6, "contact resistance must reach the motor endpoint");
    }
}

#[test]
fn empty_extra_paths_recover_both_original_motor_entry_points() {
    with_case(0.5, false, |cx, _, assembly, thermal| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let network = MotorHeatNetwork { assembly: &prepared, winding: winding(), bridge_cell: BRIDGE };
        let bare = MotorThermalNetwork {
            cells: assembly.electrical.bulk.cells, links: assembly.electrical.bulk.links,
            winding: winding(), bridge_cell: BRIDGE,
        };
        let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
        let native = rotor.at_speed(0.0, || false).unwrap()
            .midpoint_stepper(DT, 1e-12, || false).unwrap();
        let motor = network.bind(cx, &thermal, &native, state(0.0), drive()).unwrap();
        let complete = network.step_candidate(cx, &thermal, DT, &motor, &[0.0], 2.0,
            power(1000.0), controls(), &[0.0; 2]).unwrap();
        let original = bare.step_candidate(&thermal, DT, &motor, &[0.0], 2.0,
            power(1000.0), controls(), &[0.0; 2], || false).unwrap();
        assert_eq!(complete.thermal().cells(), original.thermal().cells());
        assert_eq!(complete.motor().state(), original.motor().state());
        assert_eq!(complete.report().electrical, original.report().electrical);

        let voltage_native = native.voltage_motor(
            RotorWindingMotor { drive: drive(), inductance_h: INDUCTANCE }, || false,
        ).unwrap();
        let motor = network.bind_voltage(cx, &thermal, &voltage_native,
            RotorWindingState { rotor: state(0.0), current_a: 0.0 }).unwrap();
        let complete = network.step_voltage_candidate(cx, &thermal, DT, &motor, &[0.0], 12.0,
            voltage_control(), controls(), &[0.0; 2]).unwrap();
        let original = bare.step_voltage_candidate(&thermal, DT, &motor, &[0.0], 12.0,
            voltage_control(), controls(), &[0.0; 2], || false).unwrap();
        assert_eq!(complete.thermal().cells(), original.thermal().cells());
        assert_eq!(complete.motor().state(), original.motor().state());
        assert_eq!(complete.report().electrical, original.report().electrical);
    });
}

struct TaggedSample<S> { heat: S, resolution: usize }
impl<S: AsRef<HeatNetworkSample>> AsRef<HeatNetworkSample> for TaggedSample<S> {
    fn as_ref(&self) -> &HeatNetworkSample { self.heat.as_ref() }
}
#[derive(Debug)]
enum TaggedError<E> { Heat(E), Refused }
impl<E: From<HeatNetworkError>> From<HeatNetworkError> for TaggedError<E> {
    fn from(error: HeatNetworkError) -> Self { Self::Heat(error.into()) }
}
struct Tagged<'a, A: HeatAssembly> {
    inner: &'a A,
    gate: &'a CancelGate,
    calls: Cell<usize>,
    fail_at: Cell<Option<usize>>,
    cancel_at: Cell<Option<usize>>,
}
impl<A: HeatAssembly> HeatAssembly for Tagged<'_, A> {
    type Sample = TaggedSample<A::Sample>;
    type Error = TaggedError<A::Error>;
    fn prepared_heat_network(&self) -> &PreparedHeatNetwork<'_, '_> {
        self.inner.prepared_heat_network()
    }
    fn resolve_heat(&self, cx: &Cx<'_>, thermal: &DepositionSimulation<'_>)
        -> Result<ContactThermalCoefficients<Self::Sample>, CircuitEnthalpyError<Self::Error>> {
        let resolution = self.calls.get() + 1;
        self.calls.set(resolution);
        let result = self.inner.resolve_heat(cx, thermal)
            .map_err(|error| error.map_source(TaggedError::Heat))?;
        if self.fail_at.get() == Some(resolution) {
            return Err(CircuitEnthalpyError::Source(TaggedError::Refused));
        }
        if self.cancel_at.get() == Some(resolution) { self.gate.request(); }
        Ok(ContactThermalCoefficients {
            contacts: result.contacts,
            bulk: CircuitThermalCoefficients {
                resistors: result.bulk.resistors,
                sources_w: result.bulk.sources_w,
                conductances_w_k: result.bulk.conductances_w_k,
                source: TaggedSample { heat: result.bulk.source, resolution },
            },
        })
    }
}

#[test]
fn complete_payload_typed_refusal_and_late_cancel_preserve_both_motor_histories() {
    for voltage in [false, true] {
        with_case(0.5, true, |cx, gate, assembly, thermal| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let tagged = Tagged {
                inner: &prepared, gate, calls: Cell::new(0),
                fail_at: Cell::new(None), cancel_at: Cell::new(None),
            };
            let network = MotorHeatNetwork { assembly: &tagged, winding: winding(), bridge_cell: BRIDGE };
            let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
            let native = rotor.at_speed(0.0, || false).unwrap()
                .midpoint_stepper(DT, 1e-12, || false).unwrap();
            let voltage_native = native.voltage_motor(
                RotorWindingMotor { drive: drive(), inductance_h: INDUCTANCE }, || false,
            ).unwrap();
            let motor = network.bind(cx, &thermal, &native, state(0.0), drive()).unwrap();
            let voltage_motor = network.bind_voltage(cx, &thermal, &voltage_native,
                RotorWindingState { rotor: state(0.0), current_a: 0.0 }).unwrap();
            let before_cells = thermal.cells().clone();
            let before_current = motor.state().clone();
            let before_voltage = voltage_motor.state().clone();
            let attempt = || {
                if voltage {
                    network.step_voltage_candidate(cx, &thermal, DT, &voltage_motor, &[0.0], 12.0,
                        voltage_control(), controls(), &[0.0; 2])
                        .map(|next| next.report().source.network.resolution)
                } else {
                    network.step_candidate(cx, &thermal, DT, &motor, &[0.0], 2.0,
                        power(1000.0), controls(), &[0.0; 2])
                        .map(|next| next.report().source.network.resolution)
                }
            };
            tagged.calls.set(0);
            let last = attempt().unwrap();
            assert_eq!(last, tagged.calls.get());
            assert!(last > 3, "exercise a late source callback after nonlinear work");
            tagged.calls.set(0);
            tagged.fail_at.set(Some(last));
            assert!(matches!(attempt(), Err(CircuitEnthalpyError::Source(
                MotorHeatError::Heat(TaggedError::Refused)
            ))));
            assert_eq!(thermal.cells(), &before_cells);
            assert_eq!(motor.state(), &before_current);
            assert_eq!(voltage_motor.state(), &before_voltage);

            tagged.calls.set(0);
            tagged.fail_at.set(None);
            assert_eq!(attempt().unwrap(), last);
            tagged.calls.set(0);
            tagged.cancel_at.set(Some(last));
            assert!(matches!(attempt(), Err(CircuitEnthalpyError::Cancelled)));
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(thermal.cells(), &before_cells);
            assert_eq!(motor.time_s(), 0.0);
            assert_eq!(motor.state(), &before_current);
            assert_eq!(voltage_motor.time_s(), 0.0);
            assert_eq!(voltage_motor.state(), &before_voltage);
        });
    }
}
