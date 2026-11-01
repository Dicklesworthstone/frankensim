//! Coupled capacitor charge histories through real material/contact/thermal
//! owners. Synthetic source coefficients are numerical fixtures, not fits.
#[path = "electrical_contacts/source.rs"]
mod contact_source;
#[path = "circuit_enthalpy/common.rs"]
mod thermal_source;

use fs_couple::material_circuit::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
};
use fs_matdb::MaterialCard;
use fs_material::conductor::resolve_uniform_conductor;
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::capacitance::{CapacitanceBlock, CapacitorChargeCurve, CapacitorCurve};
use fs_phs::circuit::{Branch, CircuitGraph};

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn resistor(card: &MaterialCard, branch: usize, temperature: f64) -> MaterialResistor {
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(
            card,
            &thermal_source::point(temperature),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}

fn graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 17.0 }),
            (1, 0, Branch::Resistor { ohms: 999.0 }),
            (2, 0, Branch::Resistor { ohms: 999.0 }),
            (2, 0, Branch::Capacitor { farads: 23.0 }),
        ],
        transformers: vec![],
    }
}

fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 1e-10,
        power_tolerance_w: 1e-10,
    }
}

fn charges(c: [f64; 4], voltage: [f64; 2]) -> Vec<f64> {
    vec![
        c[0] * voltage[0] + c[1] * voltage[1],
        c[2] * voltage[0] + c[3] * voltage[1],
    ]
}

fn energy(c: [f64; 4], voltage: [f64; 2]) -> f64 {
    0.5 * (c[0] * voltage[0].powi(2)
        + 2.0 * c[1] * voltage[0] * voltage[1]
        + c[3] * voltage[1].powi(2))
}

fn rc_update(c: [f64; 4], resistance: [f64; 2], voltage: [f64; 2], dt: f64) -> [f64; 2] {
    let a = [
        c[0] + dt / (2.0 * resistance[0]),
        c[1],
        c[2],
        c[3] + dt / (2.0 * resistance[1]),
    ];
    let rhs = [
        (c[0] - dt / (2.0 * resistance[0])) * voltage[0] + c[1] * voltage[1],
        c[2] * voltage[0] + (c[3] - dt / (2.0 * resistance[1])) * voltage[1],
    ];
    let determinant = a[0] * a[3] - a[1] * a[2];
    [
        (a[3] * rhs[0] - a[1] * rhs[1]) / determinant,
        (a[0] * rhs[1] - a[2] * rhs[0]) / determinant,
    ]
}

#[test]
fn picofarad_block_charge_and_material_heat_survive_a_large_held_energy_bias() {
    let c = 1e-12;
    let (card, _) = thermal_source::material(1e9, 1.0, 0.0);
    let block = CapacitanceBlock::try_new(vec![0, 3], vec![2.0 * c, -c, -c, 2.0 * c]).unwrap();
    let dt = 2e-4;
    let solve = MaterialCircuitSolve {
        energy_tolerance_j: 2e-22,
        power_tolerance_w: 1e-15,
    };
    for bias in [0.0, 1e8] {
        let mut graph = graph();
        graph.node_count = 4;
        graph
            .branches
            .push((3, 0, Branch::Capacitor { farads: 1.0 }));
        graph
            .branches
            .push((3, 0, Branch::VoltageSource { port: 0 }));
        let mut circuit = MaterialCircuit::new_with_capacitance_blocks(
            graph,
            vec![resistor(&card, 1, 300.0), resistor(&card, 2, 300.0)],
            vec![],
            vec![],
            vec![block.clone()],
            vec![],
            vec![],
            vec![11.0 * c, -7.0 * c, bias],
            &[bias],
        )
        .unwrap();
        let original = circuit.clone();
        let (mut common, mut differential) = (2.0_f64, 3.0_f64);
        let mut total_heat = 0.0;
        for _ in 0..12 {
            let before = [common + differential, common - differential];
            let energy_before = c * (common.powi(2) + 3.0 * differential.powi(2));
            common *= 9.0 / 11.0;
            differential *= 29.0 / 31.0;
            let after = [common + differential, common - differential];
            let energy_after = c * (common.powi(2) + 3.0 * differential.powi(2));
            let candidate = circuit
                .step_candidate(circuit.resistors().to_vec(), &[bias], dt, solve, || false)
                .unwrap();
            assert_eq!(candidate.circuit().capacitance_blocks(), &[block.clone()]);
            assert!(std::ptr::eq(
                candidate.circuit().capacitance_blocks()[0].matrix_f(),
                block.matrix_f()
            ));
            near(
                candidate.circuit().charges_c()[0],
                c * (common + 3.0 * differential),
                2e-22,
            );
            near(
                candidate.circuit().charges_c()[1],
                c * (common - 3.0 * differential),
                2e-22,
            );
            assert_eq!(candidate.circuit().charges_c()[2].to_bits(), bias.to_bits());
            for i in 0..2 {
                let mean = 0.5 * (before[i] + after[i]);
                near(candidate.step().resistors[i].current_a, mean / 1e9, 2e-19);
                near(
                    candidate.step().resistors[i].heat_j,
                    dt * mean * mean / 1e9,
                    2e-22,
                );
            }
            near(
                candidate.step().stored_energy_change_j(),
                energy_after - energy_before,
                2e-22,
            );
            near(candidate.step().heat_j, energy_before - energy_after, 2e-22);
            near(candidate.step().native.supplied, 0.0, 2e-22);
            if bias != 0.0 {
                assert_eq!(
                    candidate.step().stored_before_j,
                    candidate.step().stored_after_j
                );
            }
            total_heat += candidate.step().heat_j;
            near(total_heat + energy_after, 31.0 * c, 2e-21);
            circuit = candidate.into_parts().0;
        }
        assert_eq!(original.charges_c(), &[11.0 * c, -7.0 * c, bias]);
    }
}

#[test]
fn contact_and_material_resistance_rebuilds_keep_full_capacitance_and_separate_heat() {
    let (a, _) = thermal_source::material(1.0, 1.0, 0.005);
    let (b, _) = thermal_source::material(3.0, 2.0, 0.0);
    let c = [2.0, 0.5, 0.5, 3.0];
    let block = CapacitanceBlock::try_new(vec![0, 3], c.to_vec()).unwrap();
    let mut voltage = [2.0, -0.5];
    let mut circuit = MaterialCircuit::new_with_capacitance_blocks(
        graph(),
        vec![resistor(&a, 1, 300.0), resistor(&b, 2, 300.0)],
        vec![],
        vec![],
        vec![block.clone()],
        vec![],
        vec![],
        charges(c, voltage),
        &[],
    )
    .unwrap();
    let initial_energy = energy(c, voltage);
    let mut heat = 0.0;
    let dt = 0.1;
    for k in 0..8 {
        let contact_r = if k < 4 { 1.0 } else { 2.0 };
        let contact_card = contact_source::card(&a, &a, contact_r, 0.0);
        let before_charge = circuit.charges_c().to_vec();
        circuit = circuit
            .rebind_contacts_candidate(vec![contact_source::contact(9, 1, &contact_card, 1.0)], &[])
            .unwrap();
        assert_eq!(circuit.charges_c(), before_charge);
        let temperature = if k < 4 { 300.0 } else { 400.0 };
        let sources = vec![resistor(&a, 1, temperature), resistor(&b, 2, 300.0)];
        let bulk_r = [1.0 + 0.005 * (temperature - 300.0), 3.0];
        let total_r = [bulk_r[0] + contact_r, bulk_r[1]];
        let next_voltage = rc_update(c, total_r, voltage, dt);
        let candidate = circuit
            .step_with_contacts_candidate(
                sources.clone(),
                circuit.contacts().to_vec(),
                &[],
                dt,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(candidate.circuit().capacitance_blocks(), &[block.clone()]);
        assert!(std::ptr::eq(
            candidate.circuit().capacitance_blocks()[0].matrix_f(),
            block.matrix_f()
        ));
        for (actual, expected) in candidate
            .circuit()
            .charges_c()
            .iter()
            .zip(charges(c, next_voltage))
        {
            near(*actual, expected, 2e-11);
        }
        for i in 0..2 {
            let current = 0.5 * (voltage[i] + next_voltage[i]) / total_r[i];
            near(candidate.step().resistors[i].current_a, current, 2e-12);
            near(
                candidate.step().resistors[i].heat_j,
                dt * bulk_r[i] * current * current,
                2e-12,
            );
            assert_eq!(candidate.step().resistors[i].material, sources[i]);
            if i == 0 {
                near(
                    candidate.step().contacts[0].heat_j,
                    dt * contact_r * current * current,
                    2e-12,
                );
            }
        }
        heat += candidate.step().heat_j;
        near(heat + energy(c, next_voltage), initial_energy, 5e-11);
        circuit = candidate.into_parts().0;
        voltage = next_voltage;
    }
}

#[test]
fn disjoint_scalar_curve_parameter_work_retains_block_and_overlap_refuses_atomically() {
    let (card, _) = thermal_source::material(1.0, 1.0, 0.0);
    let c = [2.0, 0.4, 0.4, 3.0];
    let block = CapacitanceBlock::try_new(vec![0, 3], c.to_vec()).unwrap();
    let mut graph = graph();
    graph.node_count = 4;
    graph
        .branches
        .push((3, 0, Branch::Capacitor { farads: 1.0 }));
    graph
        .branches
        .push((3, 0, Branch::Resistor { ohms: 999.0 }));
    let voltage = [1.0, -0.5];
    let mut q = charges(c, voltage);
    q.push(2.0);
    let original = MaterialCircuit::new_with_capacitance_blocks(
        graph,
        vec![
            resistor(&card, 1, 300.0),
            resistor(&card, 2, 300.0),
            resistor(&card, 5, 300.0),
        ],
        vec![],
        vec![],
        vec![block.clone()],
        vec![],
        vec![],
        q.clone(),
        &[],
    )
    .unwrap();
    let curve = CapacitorChargeCurve::try_new(vec![
        (-4.0, -13.0),
        (-1.0, -1.0),
        (0.0, 0.0),
        (1.0, 1.0),
        (4.0, 13.0),
    ])
    .unwrap();
    let binding = CapacitorCurve { branch: 4, curve };
    let rebind = original
        .rebind_capacitor_curves_candidate(vec![binding.clone()], &[], || false)
        .unwrap();
    near(rebind.parameter_work_j(), 1.5, 1e-14);
    assert_eq!(rebind.work().len(), 1);
    assert_eq!(rebind.work()[0].branch, 4);
    assert_eq!(rebind.circuit().charges_c(), q);
    assert_eq!(rebind.circuit().capacitance_blocks(), &[block.clone()]);
    assert!(std::ptr::eq(
        rebind.circuit().capacitance_blocks()[0].matrix_f(),
        block.matrix_f()
    ));
    let saved = rebind.circuit().clone();
    let overlapping = CapacitorCurve {
        branch: 0,
        curve: binding.curve.clone(),
    };
    assert!(
        rebind
            .circuit()
            .rebind_capacitor_curves_candidate(vec![binding.clone(), overlapping], &[], || false)
            .is_err()
    );
    assert!(matches!(
        rebind
            .circuit()
            .rebind_capacitor_curves_candidate(vec![binding], &[], || true),
        Err(MaterialCircuitError::Cancelled)
    ));
    assert_eq!(rebind.circuit(), &saved);
    assert_eq!(original.charges_c(), q);
    let dt = 0.05;
    let step = rebind
        .circuit()
        .step_candidate(
            rebind.circuit().resistors().to_vec(),
            &[],
            dt,
            controls(),
            || false,
        )
        .unwrap();
    let next_voltage = rc_update(c, [1.0, 1.0], voltage, dt);
    let scalar_q = 0.75 + 1.25 * (1.0 - 2.0 * dt) / (1.0 + 2.0 * dt);
    for (actual, expected) in step
        .circuit()
        .charges_c()
        .iter()
        .zip(charges(c, next_voltage).into_iter().chain([scalar_q]))
    {
        near(*actual, expected, 3e-11);
    }
    let scalar_h = 0.5 + 2.0 * scalar_q * scalar_q - 3.0 * scalar_q + 1.0;
    near(
        step.step().heat_j + energy(c, next_voltage) + scalar_h,
        energy(c, voltage) + 3.5,
        5e-11,
    );
    assert_eq!(step.circuit().capacitance_blocks(), &[block]);
}

#[test]
fn finite_resistor_temperature_feedback_preserves_block_and_per_body_caloric_balance() {
    use fs_couple::circuit_enthalpy::step_candidate;
    use fs_thermal::enthalpy::CellId;
    let (a, ca) = thermal_source::material(1.0, 1.0, 0.005);
    let (b, cb) = thermal_source::material(3.0, 2.0, 0.002);
    let cards = [a, b];
    let curves = [ca, cb];
    let mut thermal = thermal_source::initial(&curves, 1.0, None);
    let c = [2.0, 0.5, 0.5, 3.0];
    let block = CapacitanceBlock::try_new(vec![0, 3], c.to_vec()).unwrap();
    let mut voltage = [3.0, -0.5];
    let mut circuit = MaterialCircuit::new_with_capacitance_blocks(
        graph(),
        thermal_source::coefficients(&thermal, &cards, &curves)
            .unwrap()
            .resistors,
        vec![],
        vec![],
        vec![block.clone()],
        vec![],
        vec![],
        charges(c, voltage),
        &[],
    )
    .unwrap();
    let initial_energy = thermal.energy_j() + energy(c, voltage);
    let dt = 0.25;
    for k in 1..=4 {
        let next = step_candidate(
            &thermal,
            f64::from(k) * dt,
            thermal_source::input(&circuit),
            thermal_source::config(),
            &[0.0; 2],
            |state| thermal_source::coefficients(state, &cards, &curves),
            || false,
        )
        .unwrap();
        assert_eq!(next.circuit().capacitance_blocks(), &[block.clone()]);
        assert!(std::ptr::eq(
            next.circuit().capacitance_blocks()[0].matrix_f(),
            block.matrix_f()
        ));
        let temperature = [
            next.thermal().temperature_k(CellId(30)).unwrap(),
            next.thermal().temperature_k(CellId(10)).unwrap(),
        ];
        let r = [
            1.0 + 0.005 * (temperature[0] - 300.0),
            3.0 * (1.0 + 0.002 * (temperature[1] - 300.0)),
        ];
        let next_voltage = rc_update(c, r, voltage, dt);
        for (actual, expected) in next
            .circuit()
            .charges_c()
            .iter()
            .zip(charges(c, next_voltage))
        {
            near(*actual, expected, 2e-8);
        }
        for i in 0..2 {
            let cell = [CellId(30), CellId(10)][i];
            let mean = 0.5 * (voltage[i] + next_voltage[i]);
            let heat = dt * mean * mean / r[i];
            near(next.report().electrical.resistors[i].heat_j, heat, 2e-8);
            near(
                next.report().electrical.resistors[i]
                    .material
                    .conductor
                    .resistance_ohm(),
                r[i],
                2e-10,
            );
            near(
                [1.0, 2.0][i] * (temperature[i] - thermal.temperature_k(cell).unwrap()),
                heat + dt * (temperature[1 - i] - temperature[i]),
                3e-6,
            );
        }
        near(
            next.thermal().energy_j() + energy(c, next_voltage),
            initial_energy,
            5e-6,
        );
        let (after_thermal, after_circuit, _) = next.into_parts();
        thermal = after_thermal;
        circuit = after_circuit;
        voltage = next_voltage;
    }
    assert!(thermal.temperature_k(CellId(30)).unwrap() > 300.0);
    assert!(circuit.resistors()[0].conductor.resistance_ohm() > 1.0);
}
