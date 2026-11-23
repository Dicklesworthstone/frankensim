//! The constant-capacitance fixtures as declared graphs (`battery_motor/graph.rs`).

use fs_phs::circuit::converter::commutation::capacitance::CapacitiveCommutationParameters;
use fs_qty::Capacitance;
use fs_system::commutated_drive::DeviceStorage;

use super::graph::{
    Case, closed_form_step, kernel_refusal, kernel_run, leak, owner_check, steady_check,
};
use super::*;

const COSS: [f64; 4] = [0.025, 0.04, 0.035, 0.05];

struct CapacitiveFixture {
    base: Fixture,
    keys: [PropertyKey; 4],
}
fn capacitive_fixture(cp: f64, threshold_slope: f64, cap_upper_k: f64) -> CapacitiveFixture {
    let f = fixture(cp, threshold_slope, 350.0, 100.0);
    let keys = core::array::from_fn(|i| {
        PropertyKey::new(format!("device_output_capacitance_{i}"), Capacitance::DIMS)
    });
    let mut claims = ClaimSet::new();
    for (_, claim) in f.card.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    for (i, key) in keys.iter().enumerate() {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Scalar {
                    value: COSS[i],
                    dims: Capacitance::DIMS,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    299.0,
                    cap_upper_k,
                ),
                uncertainty: UncertaintyModel::Unstated,
                // The Scalar is a constant law: LinearInside has no knot span on a scalar.
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source:
                        "synthetic constant device capacitor fixture, not qualified device data"
                            .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-device".into(),
            phase: "solid".into(),
            process: "authored fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &common::point(299.0),
            upper: &common::point(350.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    CapacitiveFixture {
        base: Fixture {
            card,
            caloric,
            keys: f.keys,
            cp,
        },
        keys,
    }
}
fn initial(steady: &SteadyDrive, node_offset_v: [f64; 2]) -> [f64; 11] {
    let p = parameters();
    let current = steady.state[3] / p.inductance_h;
    let input = steady.state[4] / p.bus_capacitance_f[0];
    let output = steady.state[5] / p.bus_capacitance_f[1];
    // Independent DC branch KCL from the steady fixture; only this explicitly
    // authored initialization is steady. Subsequent nodes are dynamic charges.
    let left = if steady.gates.on[0] {
        (10.0 * input - current) / 10.01
    } else {
        (5.01 * input + 2.5 - current) / 5.02
    };
    let right = left - 0.1 * current;
    let c = CapacitiveCommutationParameters {
        converter: p,
        output_capacitance_f: COSS,
    }
    .state(
        current,
        [input, output],
        [left + node_offset_v[0], right + node_offset_v[1]],
    )
    .unwrap();
    let mut state = [0.0; 11];
    state[..3].copy_from_slice(&steady.state[..3]);
    state[3..8].copy_from_slice(&c);
    state[8..].copy_from_slice(&steady.state[6..]);
    state
}

fn capacitive_case<'s>(
    f: &'static CapacitiveFixture,
    steady: &'s SteadyDrive,
    start: [f64; 11],
    dt_s: f64,
    substeps: usize,
    load_cell: Option<CellId>,
) -> Case<'s> {
    Case {
        fixture: &f.base,
        context: leak(common::point(300.0)),
        storage_keys: &f.keys,
        storage: DeviceStorage::Capacitance(COSS),
        steady,
        start,
        dt_s,
        substeps,
        load_cell,
        schedule: Vec::new(),
        linked: false,
        coolant: None,
        drivers: None,
    }
}

#[test]
fn constant_capacitors_as_a_graph_keep_the_steady_motoring_and_regenerative_drive() {
    for regeneration in [false, true] {
        let f = leak(capacitive_fixture(1.0, 0.0, 350.0));
        let steady = steady_drive(regeneration);
        let start = initial(&steady, [0.0; 2]);
        let case = capacitive_case(f, &steady, start, 1.0 / 64.0, 4, Some(IDS[3]));
        // On 9b41f08404 within 1.1e-15 relative and 2.7e-14 K of the closed
        // form.
        steady_check(&case, &kernel_run(&case, 3), [2e-14, 1e-12]);
    }
}

#[test]
fn constant_capacitors_as_a_graph_follow_the_hot_diode_threshold() {
    // The dynamic fixture: cells of 0.01 J/K, a diode threshold
    // rising 20 mV/K, the nodes started off their steady values.
    let f = leak(capacitive_fixture(0.01, 0.02, 350.0));
    let steady = steady_drive(false);
    let start = initial(&steady, [-0.3, 0.4]);
    let case = capacitive_case(f, &steady, start, 1.0 / 1024.0, 1, None);
    let graph = kernel_run(&case, 3);
    // On 9b41f08404 the owner by hand bit for bit, each cell's balance
    // within 5.5e-16 J.
    owner_check(&case, &graph, 0.02, [1e-14, 1e-14]);
    // The feedback is real: the output high package, whose diode conducts,
    // warms past 300.8 K, so its threshold rises by more than 16 mV.
    assert!(graph[2].kelvin[3] > 300.8, "{:?}", graph[2].kelvin);
}

#[test]
fn constant_capacitors_as_a_graph_follow_a_gate_program_that_changes() {
    // The dynamic fixture continues with only the output high
    // switch on; a third step returns to the first gates.
    let f = leak(capacitive_fixture(0.01, 0.02, 350.0));
    let steady = steady_drive(false);
    let start = initial(&steady, [-0.3, 0.4]);
    let other = DeviceGateState {
        on: [false, false, true, false],
    };
    let dt = 1.0 / 1024.0;
    let case = Case {
        schedule: vec![steady.gates, other, steady.gates],
        ..capacitive_case(f, &steady, start, dt, 1, None)
    };
    let graph = kernel_run(&case, 3);
    // On 9b41f08404 the owner by hand bit for bit, each cell's balance
    // within 5.5e-16 J.
    owner_check(&case, &graph, 0.02, [1e-14, 1e-14]);
    // The schedule matters: with the first gates kept, the second step ends
    // elsewhere.
    let kept = kernel_run(&capacitive_case(f, &steady, start, dt, 1, None), 2);
    assert_ne!(
        kept[1].state.map(f64::to_bits),
        graph[1].state.map(f64::to_bits)
    );
}

#[test]
fn constant_capacitors_as_a_graph_carry_the_link_between_cells() {
    // The heat-assembly fixture: steady motoring with a link of
    // 0.25 W/K between the inductor's cell and the motor's. Then the hot
    // dynamic fixture with the same link.
    let steady = steady_drive(false);
    let fixtures = [
        ((1.0, 0.0), [0.0; 2], 1.0 / 64.0, 4),
        ((0.01, 0.02), [-0.3, 0.4], 1.0 / 1024.0, 1),
    ];
    for ((cp, slope), offset, dt, substeps) in fixtures {
        let f = leak(capacitive_fixture(cp, slope, 350.0));
        let start = initial(&steady, offset);
        let case = Case {
            linked: true,
            ..capacitive_case(f, &steady, start, dt, substeps, None)
        };
        // On 9b41f08404 the owner by hand bit for bit, each cell's balance
        // within 4.4e-14 J.
        owner_check(&case, &kernel_run(&case, 3), slope, [1e-14, 1e-12]);
        // The steady fixture's first step also meets the closed form through
        // the link: within 5.7e-14 K.
        assert!(substeps != 4 || closed_form_step(&case) < 1e-12);
    }
}

#[test]
fn constant_capacitors_as_a_graph_refuse_a_claim_out_of_its_domain_and_a_partial_period() {
    let steady = steady_drive(false);
    let start = initial(&steady, [0.0; 2]);
    let dt = 1.0 / 64.0;
    // Capacitors qualified only up to 300.000001 K, which the step passes.
    let f = leak(capacitive_fixture(1.0, 0.0, 300.000001));
    let case = capacitive_case(f, &steady, start, dt, 4, None);
    let graph = kernel_refusal(&case, dt);
    assert!(graph.contains("NoClaimInDomain"), "{graph}");
    assert!(graph.contains("device_output_capacitance_0"), "{graph}");
    assert!(
        graph.contains(r#"Law { port: "heat_input_high""#),
        "{graph}"
    );
    // A step that is not one gate program.
    let f = leak(capacitive_fixture(1.0, 0.0, 350.0));
    let case = capacitive_case(f, &steady, start, dt, 4, None);
    let graph = kernel_refusal(&case, 2.0 * dt);
    assert!(graph.contains("one period of its program"), "{graph}");
}
