//! Source-bound contact transport reaches both induction windings before acceptance.
//! Synthetic material curves verify the discrete coupling, not a measured machine.
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
    motor::induction::heat::{InductionHeatError, InductionHeatNetwork},
};
use fs_evidence::ValidityDomain;
use fs_exec::{CancelGate, Cx};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterfaceSystemCard, InterpolationPolicy, PropertyClaim,
    PropertyValue, Provenance, SelectionPolicy, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::{QuantitySpec, Temperature};

const DT: f64 = 0.04;
const INPUTS: [f64; 4] = [8.0, 0.0, 0.0, 0.0];

fn with_case<R>(
    contact_resistance: Option<f64>,
    run: impl FnOnce(&Cx<'_>, &CancelGate, &MaterialHeatNetwork<'_>, DepositionSimulation<'_>) -> R,
) -> R {
    let (cards, curves) = sources();
    let context = common::point(300.0);
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let cells = bindings(&cards, &curves, &context, Some(&key));
    let links = [ConductionLink { link: 7, distance_a_m: 0.5, distance_b_m: 0.5 }];
    let bulk = MaterialNetwork { cells: &cells, resistors: &[], links: &links };
    let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[], electrical_bases: &[] };
    let temperature = || QuantitySpec::dimensional(Temperature::DIMS);
    let resistance = contact_resistance.unwrap_or(1.0);
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(RP, RD),
        value: PropertyValue::Curve {
            abscissa: "aT".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, resistance), (2000.0, resistance)], dims: RD,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("aT", temperature(), 200.0, 2000.0)
            .with_quantity("bT", temperature(), 200.0, 2000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic stator/rotor thermal contact".into(),
            license: "CC0-1.0".into(), artifact: None,
        },
    }).unwrap();
    let interface = InterfaceSystemCard::assemble(
        SurfaceSpec { material: cards[0].id().clone(), texture_frame: "stator".into() },
        SurfaceSpec { material: cards[1].id().clone(), texture_frame: "rotor".into() },
        SystemContext {
            medium: "fixed thermal link".into(), third_body: None,
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
            cell: [STATOR, ROTOR][i], material_id: MaterialId(i as u64), curve: curves[i].curve(),
        }),
        point: ThermalInterfacePoint {
            interface: &interface, materials: [&cards[0], &cards[1]], context: &contact_context,
            temperature_axes: ["aT", "bT"], liquid_fraction_axes: [None, None], mean_temperature_axis: None,
        },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }];
    let assembly = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: if contact_resistance.is_some() { &contacts } else { &[] },
        radiation: None,
    };
    let thermal = body(&curves, true, 0.0);
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate, arena,
            fs_exec::StreamKey { seed: 89, kernel_id: 67, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic,
        );
        run(&cx, &gate, &assembly, thermal)
    })
}

#[test]
fn thermal_contact_changes_both_currents_and_closes_induction_energy_across_steps() {
    let mut endpoints = Vec::new();
    for resistance in [0.5, 8.0] {
        endpoints.push(with_case(Some(resistance), |cx, _, assembly, mut thermal| {
            let prepared = assembly.prepare(cx, &thermal).unwrap();
            let network = InductionHeatNetwork {
                assembly: &prepared, stator: winding(0), rotor: winding(1), bearing_cell: STATOR,
            };
            let mut motor = network.bind(cx, &thermal, [0.0; 7], parameters()).unwrap();
            let (mut current, mut rise) = ([0.0; 2], [0.0; 2]);
            // Independent Fourier half-cell resistances are both 0.5 K/W.
            let conductance = 1.0 / (1.0 + resistance);
            for index in 1..=2 {
                let old_state = *motor.state();
                let old_cells = thermal.cells().clone();
                let (expected_current, expected_rise, heat) = reference(
                    current, rise, INPUTS[0], DT, conductance, 0.0,
                );
                let next = network.step_candidate(
                    cx, &thermal, f64::from(index) * DT, &motor, &INPUTS, controls(), &[0.0; 2],
                ).unwrap();
                let point = next.motor().operating_point();
                near(point.stator_current_alpha_beta_a[0], expected_current[0]);
                near(point.rotor_current_dq_a[0], expected_current[1]);
                for i in 0..2 {
                    near(next.thermal().temperature_k([STATOR, ROTOR][i]).unwrap(), 300.0 + expected_rise[i]);
                    near(next.report().source.resistance_ohm[i], R0[i] * (1.0 + ALPHA[i] * expected_rise[i]));
                }
                near(next.report().source.network.thermal_contacts[0].conductance_w_k, conductance);
                let step = &next.report().electrical;
                near(step.stator_copper_heat_j, heat[0]);
                near(step.rotor_copper_heat_j, heat[1]);
                let supplied: f64 = step.electrical_work_j.iter().chain(&step.shaft_work_j).sum();
                near(energy(parameters(), next.motor().state()) - energy(parameters(), &old_state)
                    + next.thermal().energy_j() - thermal.energy_j(), supplied);
                near(next.report().total_energy_residual_j, 0.0);
                near(next.report().non_electrical_heat_j, 0.0);
                assert_eq!(motor.state(), &old_state);
                assert_eq!(thermal.cells(), &old_cells);
                assert_eq!(motor.stator_winding().material().card_identity(),
                    next.motor().stator_winding().material().card_identity());
                assert_eq!(motor.rotor_winding().material().card_identity(),
                    next.motor().rotor_winding().material().card_identity());
                (thermal, motor, _) = next.into_parts();
                (current, rise) = (expected_current, expected_rise);
            }
            [thermal.temperature_k(STATOR).unwrap(), thermal.temperature_k(ROTOR).unwrap(),
                motor.operating_point().stator_current_alpha_beta_a[0], motor.operating_point().rotor_current_dq_a[0]]
        }));
    }
    for (low, high) in endpoints[0].iter().zip(&endpoints[1]) {
        assert!((low - high).abs() > 1e-6, "contact resistance must reach both winding currents and temperatures");
    }
}

#[test]
fn empty_extra_heat_paths_recover_the_original_induction_driver() {
    with_case(None, |cx, _, assembly, thermal| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let network = InductionHeatNetwork {
            assembly: &prepared, stator: winding(0), rotor: winding(1), bearing_cell: STATOR,
        };
        let bare = InductionThermalNetwork {
            cells: assembly.electrical.bulk.cells, links: assembly.electrical.bulk.links,
            stator: winding(0), rotor: winding(1), bearing_cell: STATOR,
        };
        let motor = network.bind(cx, &thermal, [0.0; 7], parameters()).unwrap();
        let expected = bare.step_candidate(&thermal, DT, &motor, &INPUTS, controls(), &[0.0; 2], || false).unwrap();
        let actual = network.step_candidate(cx, &thermal, DT, &motor, &INPUTS, controls(), &[0.0; 2]).unwrap();
        assert_eq!(actual.motor().state(), expected.motor().state());
        assert_eq!(actual.motor().parameters(), expected.motor().parameters());
        assert_eq!(actual.thermal().cells(), expected.thermal().cells());
        assert_eq!(actual.report().electrical, expected.report().electrical);
    })
}

#[test]
fn saturation_survives_complete_heat_assembly_binding_and_repeated_steps() {
    with_case(None, |cx, _, assembly, mut thermal| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let network = InductionHeatNetwork {
            assembly: &prepared, stator: winding(0), rotor: winding(1), bearing_cell: STATOR,
        };
        let bare = InductionThermalNetwork {
            cells: assembly.electrical.bulk.cells, links: assembly.electrical.bulk.links,
            stator: winding(0), rotor: winding(1), bearing_cell: STATOR,
        };
        let curve = magnetizing_curve();
        let mut motor = network.bind_with_magnetizing_curve(
            cx, &thermal, [0.0; 7], parameters(), curve.clone(),
        ).unwrap();
        let bare_motor = bare.bind_with_magnetizing_curve(
            &thermal, [0.0; 7], parameters(), curve.clone(), || false,
        ).unwrap();
        assert_eq!(motor.state(), bare_motor.state());
        assert_eq!(motor.magnetizing_curve(), bare_motor.magnetizing_curve());
        for interval in 1..=3 {
            let target = interval as f64 * DT;
            let expected = bare.step_candidate(
                &thermal, target, &motor, &INPUTS, controls(), &[0.0; 2], || false,
            ).unwrap();
            let actual = network.step_candidate(
                cx, &thermal, target, &motor, &INPUTS, controls(), &[0.0; 2],
            ).unwrap();
            assert_eq!(actual.motor().state(), expected.motor().state());
            assert_eq!(actual.motor().magnetizing_curve(), Some(&curve));
            assert_eq!(actual.thermal().cells(), expected.thermal().cells());
            assert_eq!(actual.report().electrical, expected.report().electrical);
            (thermal, motor, _) = actual.into_parts();
        }
        let point = motor.operating_point();
        assert!((point.stator_current_alpha_beta_a[0] + point.rotor_current_dq_a[0]).abs() > 2.0);
    })
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
fn late_assembly_refusal_cancellation_and_budget_preserve_all_induction_history() {
    with_case(Some(0.5), |cx, gate, assembly, thermal| {
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let tagged = Tagged {
            inner: &prepared, gate, calls: Cell::new(0),
            fail_at: Cell::new(None), cancel_at: Cell::new(None),
        };
        let network = InductionHeatNetwork {
            assembly: &tagged, stator: winding(0), rotor: winding(1), bearing_cell: STATOR,
        };
        let motor = network.bind(cx, &thermal, [0.0; 7], parameters()).unwrap();
        let before_cells = thermal.cells().clone();
        let before_state = *motor.state();
        let before_parameters = motor.parameters();
        let attempt = |solve| network.step_candidate(cx, &thermal, DT, &motor, &INPUTS, solve, &[0.0; 2]);
        tagged.calls.set(0);
        let expected = attempt(controls()).unwrap();
        let last = expected.report().source.network.resolution;
        assert_eq!(last, tagged.calls.get());
        assert!(last > 3, "exercise final source refusal after nonlinear induction trials");
        tagged.calls.set(0);
        tagged.fail_at.set(Some(last));
        assert!(matches!(attempt(controls()), Err(CircuitEnthalpyError::Source(
            InductionHeatError::Heat(TaggedError::Refused)
        ))));
        assert_eq!(thermal.cells(), &before_cells);
        assert_eq!(motor.state(), &before_state);
        assert_eq!(motor.parameters(), before_parameters);

        tagged.calls.set(0);
        tagged.fail_at.set(None);
        let mut limited = controls();
        limited.thermal.maximum_circuit_trials = 1;
        assert!(matches!(attempt(limited), Err(CircuitEnthalpyError::CircuitBudget)));
        tagged.calls.set(0);
        let replay = attempt(controls()).unwrap();
        assert_eq!(replay.motor().state(), expected.motor().state());
        assert_eq!(replay.thermal().cells(), expected.thermal().cells());
        assert_eq!(replay.report().electrical, expected.report().electrical);

        tagged.calls.set(0);
        tagged.cancel_at.set(Some(last));
        assert!(matches!(attempt(controls()), Err(CircuitEnthalpyError::Cancelled)));
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells(), &before_cells);
        assert_eq!(motor.time_s(), 0.0);
        assert_eq!(motor.state(), &before_state);
        assert_eq!(motor.parameters(), before_parameters);
    })
}
