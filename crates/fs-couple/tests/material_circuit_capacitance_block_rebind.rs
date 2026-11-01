//! Fixed-charge changes of full reciprocal capacitance. Independent 2x2
//! inverse and midpoint-RC oracles retain cross energy and signed parameter work.
#[path = "electrical_contacts/source.rs"]
mod contact_source;

use fs_couple::material_circuit::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
    storage::{MaterialStorage, MaterialStorageElement},
};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_material::{
    conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor},
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, ResolvedCapacitor,
    },
    state_point::MaterialPropertySelection,
};
use fs_phs::circuit::{
    Branch, CircuitGraph, MutualInductance,
    capacitance::{CapacitanceBlock, CapacitorChargeCurve, CapacitorCurve},
};

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "actual={actual:e}, expected={expected:e}, allowance={allowance:e}"
    );
}

fn block(branches: [usize; 2], c: [f64; 4]) -> CapacitanceBlock {
    CapacitanceBlock::try_new(branches.to_vec(), c.to_vec()).unwrap()
}

fn inverse(c: [f64; 4]) -> [f64; 4] {
    let determinant = c[0] * c[3] - c[1] * c[2];
    [
        c[3] / determinant,
        -c[1] / determinant,
        -c[2] / determinant,
        c[0] / determinant,
    ]
}

fn voltage(c: [f64; 4], q: [f64; 2]) -> [f64; 2] {
    let k = inverse(c);
    [k[0] * q[0] + k[1] * q[1], k[2] * q[0] + k[3] * q[1]]
}

fn energy(c: [f64; 4], q: [f64; 2]) -> f64 {
    (c[3] * q[0].powi(2) - 2.0 * c[1] * q[0] * q[1] + c[0] * q[1].powi(2))
        / (2.0 * (c[0] * c[3] - c[1] * c[2]))
}

fn resistor(card: &MaterialCard, branch: usize) -> MaterialResistor {
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(
            card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}

// (C + dt G/2) v1 = (C - dt G/2) v0, G = diag(1/R).
fn rc_step(c: [f64; 4], r: [f64; 2], v: [f64; 2], dt: f64) -> [f64; 2] {
    let a = [
        c[0] + dt / (2.0 * r[0]),
        c[1],
        c[2],
        c[3] + dt / (2.0 * r[1]),
    ];
    let rhs = [
        (c[0] - dt / (2.0 * r[0])) * v[0] + c[1] * v[1],
        c[2] * v[0] + (c[3] - dt / (2.0 * r[1])) * v[1],
    ];
    let det = a[0] * a[3] - a[1] * a[2];
    [
        (a[3] * rhs[0] - a[1] * rhs[1]) / det,
        (a[0] * rhs[1] - a[2] * rhs[0]) / det,
    ]
}

fn capacitors(count: usize) -> CircuitGraph {
    CircuitGraph {
        node_count: count + 1,
        branches: (1..=count)
            .map(|node| (node, 0, Branch::Capacitor { farads: 1.0 }))
            .collect(),
        transformers: vec![],
    }
}

fn simple(
    graph: CircuitGraph,
    blocks: Vec<CapacitanceBlock>,
    charges: Vec<f64>,
    inputs: &[f64],
) -> MaterialCircuit {
    MaterialCircuit::new_with_capacitance_blocks(
        graph,
        vec![],
        vec![],
        vec![],
        blocks,
        vec![],
        vec![],
        charges,
        inputs,
    )
    .unwrap()
}

#[test]
fn full_mutual_parameter_work_preserves_charge_flux_sources_and_subsequent_contact_heat() {
    let old_c = [2.0, 0.5, 0.5, 3.0];
    let new_c = [1.25, -0.25, -0.25, 2.5];
    let a = contact_source::scalar_card("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, 2.0);
    let b = contact_source::scalar_card("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, 5.0);
    let initial_contact = contact_source::card(&a, &a, 3.0, 0.0);
    let next_contact = contact_source::card(&a, &a, 4.0, 0.0);
    let graph = CircuitGraph {
        node_count: 6,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 17.0 }),
            (3, 0, Branch::Inductor { henries: 2.0 }),
            (1, 0, Branch::Resistor { ohms: 999.0 }),
            (2, 0, Branch::Capacitor { farads: 23.0 }),
            (4, 0, Branch::Inductor { henries: 3.0 }),
            (2, 0, Branch::Resistor { ohms: 999.0 }),
            (3, 0, Branch::VoltageSource { port: 0 }),
            (4, 0, Branch::VoltageSource { port: 1 }),
            (5, 0, Branch::Capacitor { farads: 19.0 }),
            (5, 0, Branch::VoltageSource { port: 2 }),
        ],
        transformers: vec![],
    };
    let inputs = [0.0, 0.0, 0.5];
    let flux = vec![0.7, -0.4];
    let original_block = block([0, 3], old_c);
    let proposed_block = block([0, 3], new_c);
    for sign in [-1.0, 1.0] {
        let q = [sign * 2.0, sign * 0.25];
        let original = MaterialCircuit::new_with_capacitance_blocks(
            graph.clone(),
            vec![resistor(&a, 2), resistor(&b, 5)],
            vec![contact_source::contact(9, 2, &initial_contact, 1.0)],
            vec![MutualInductance {
                branches: (1, 4),
                henries: 0.25,
            }],
            vec![original_block.clone()],
            vec![CapacitorCurve {
                branch: 8,
                curve: CapacitorChargeCurve::try_new(vec![(-2.0, -4.0), (0.0, 0.0), (2.0, 4.0)])
                    .unwrap(),
            }],
            flux.clone(),
            vec![q[0], q[1], 0.25],
            &inputs,
        )
        .unwrap();
        let saved = original.clone();
        let update = original
            .rebind_capacitance_blocks_candidate(vec![proposed_block.clone()], &inputs, || false)
            .unwrap();
        let updated = update.circuit();
        assert_eq!(updated.charges_c(), original.charges_c());
        assert_eq!(updated.fluxes_wb(), original.fluxes_wb());
        assert_eq!(updated.graph(), original.graph());
        assert_eq!(updated.mutual_inductance(), original.mutual_inductance());
        assert_eq!(updated.resistors(), original.resistors());
        assert_eq!(updated.contacts(), original.contacts());
        assert_eq!(updated.capacitor_curves(), original.capacitor_curves());
        let expected_work = energy(new_c, q) - energy(old_c, q);
        assert!(expected_work > 0.0);
        near(update.parameter_work_j(), expected_work, 2e-14);
        assert_eq!(
            update.work().iter().map(|w| w.branch).collect::<Vec<_>>(),
            vec![0, 3]
        );
        let (v0, v1) = (voltage(old_c, q), voltage(new_c, q));
        // A coupled coordinate's signed contribution can be negative while the
        // complete SPD energy remains positive. Rejecting this loses valid C.
        assert!(0.5 * q[1] * v0[1] < 0.0);
        for (i, row) in update.work().iter().enumerate() {
            near(row.stored_before_j, 0.5 * q[i] * v0[i], 2e-14);
            near(row.stored_after_j, 0.5 * q[i] * v1[i], 2e-14);
            near(row.supplied_work_j, 0.5 * q[i] * (v1[i] - v0[i]), 2e-14);
        }
        near(
            update.work().iter().map(|w| w.supplied_work_j).sum(),
            expected_work,
            2e-14,
        );
        let reverse = updated
            .rebind_capacitance_blocks_candidate(vec![original_block.clone()], &inputs, || false)
            .unwrap();
        near(reverse.parameter_work_j(), -expected_work, 2e-14);
        assert_eq!(reverse.circuit(), &saved);
        assert_eq!(original, saved);

        // Continue through an actual material/contact rebuild. Parameter work
        // is paid before this step; only resistive loss is deposited as heat.
        let accepted = update
            .into_parts()
            .0
            .rebind_contacts_candidate(
                vec![contact_source::contact(9, 2, &next_contact, 1.0)],
                &inputs,
            )
            .unwrap();
        assert_eq!(accepted.capacitance_blocks(), &[proposed_block.clone()]);
        let dt = 0.1;
        let resistance = [6.0, 5.0];
        let next_v = rc_step(new_c, resistance, v1, dt);
        let next_q = [
            new_c[0] * next_v[0] + new_c[1] * next_v[1],
            new_c[2] * next_v[0] + new_c[3] * next_v[1],
        ];
        let candidate = accepted
            .step_with_contacts_candidate(
                accepted.resistors().to_vec(),
                accepted.contacts().to_vec(),
                &inputs,
                dt,
                MaterialCircuitSolve {
                    energy_tolerance_j: 1e-10,
                    power_tolerance_w: 1e-10,
                },
                || false,
            )
            .unwrap();
        near(candidate.circuit().charges_c()[0], next_q[0], 2e-12);
        near(candidate.circuit().charges_c()[1], next_q[1], 2e-12);
        assert_eq!(
            candidate.circuit().charges_c()[2].to_bits(),
            0.25_f64.to_bits()
        );
        assert_eq!(candidate.circuit().fluxes_wb(), flux);
        let current = [
            0.5 * (v1[0] + next_v[0]) / resistance[0],
            0.5 * (v1[1] + next_v[1]) / resistance[1],
        ];
        let heat = dt * (6.0 * current[0].powi(2) + 5.0 * current[1].powi(2));
        near(
            candidate.step().contacts[0].heat_j,
            dt * 4.0 * current[0].powi(2),
            2e-13,
        );
        near(
            candidate.step().resistors[0].heat_j,
            dt * 2.0 * current[0].powi(2),
            2e-13,
        );
        near(candidate.step().heat_j, heat, 2e-13);
        near(
            candidate.step().stored_energy_change_j(),
            energy(new_c, next_q) - energy(new_c, q),
            2e-13,
        );
        near(
            energy(new_c, next_q) + candidate.step().heat_j,
            energy(old_c, q) + expected_work,
            2e-13,
        );
    }
}

#[test]
fn equivalent_reordered_blocks_reuse_original_preparation_and_require_exactly_zero_work() {
    let a = block([0, 2], [2.0, 0.5, 0.5, 3.0]);
    let b = block([1, 3], [4.0, -0.25, -0.25, 2.0]);
    let original = simple(
        capacitors(4),
        vec![a.clone(), b.clone()],
        vec![1.0, -2.0, 3.0, 0.5],
        &[],
    );
    // Reorder both the block list and each block's rows/columns. Re-preparation
    // may round an inverse differently; physical C equality must reuse the
    // original admitted preparation, not create artificial parameter work.
    let reversed_a = block([2, 0], [3.0, 0.5, 0.5, 2.0]);
    let reversed_b = block([3, 1], [2.0, -0.25, -0.25, 4.0]);
    let update = original
        .rebind_capacitance_blocks_candidate(vec![reversed_b, reversed_a], &[], || false)
        .unwrap();
    assert_eq!(update.parameter_work_j(), 0.0);
    assert!(update.work().iter().all(|w| w.supplied_work_j == 0.0));
    assert_eq!(update.circuit().charges_c(), original.charges_c());
    for before in [a, b] {
        let after = update
            .circuit()
            .capacitance_blocks()
            .iter()
            .find(|candidate| candidate.branches() == before.branches())
            .unwrap();
        assert!(std::ptr::eq(before.matrix_f(), after.matrix_f()));
        assert!(std::ptr::eq(
            before.inverse_matrix_per_f(),
            after.inverse_matrix_per_f()
        ));
    }
}

#[test]
fn tiny_mutual_inverse_change_survives_equal_biased_energies_and_rounded_endpoint_voltage() {
    let tiny = 2.0_f64.powi(-80);
    let large = 2.0_f64.powi(40);
    let before = block([0, 1], [1.0, 0.0, 0.0, 1.0]);
    let after = block([0, 1], [1.0, tiny, tiny, 1.0]);
    let charge = [large, 1.0];
    let original = simple(capacitors(2), vec![before.clone()], charge.to_vec(), &[]);
    let update = original
        .rebind_capacitance_blocks_candidate(vec![after.clone()], &[], || false)
        .unwrap();
    let k0 = before.inverse_matrix_per_f();
    let k1 = after.inverse_matrix_per_f();
    assert_eq!(k0[0], k1[0]);
    assert_eq!(k0[3], k1[3]);
    near(k1[1], -tiny, 2e-15 * tiny);
    // Oracle uses the native represented inverse coefficients, before any
    // multiplication by huge q; mathematically unrepresented O(tiny^2) changes
    // are not part of that admitted native Hamiltonian.
    let expected =
        0.5 * (k1[0] - k0[0]) * large * large + (k1[1] - k0[1]) * large + 0.5 * (k1[3] - k0[3]);
    near(expected, -2.0_f64.powi(-40), 1e-24);
    near(update.parameter_work_j(), expected, 1e-24);
    assert_eq!(
        before.stored_energy_j(&charge).unwrap(),
        after.stored_energy_j(&charge).unwrap()
    );
    assert_eq!(
        before.voltages_v(&charge).unwrap()[0],
        after.voltages_v(&charge).unwrap()[0]
    );
    assert_eq!(
        update.work()[0].stored_before_j,
        update.work()[0].stored_after_j
    );
    for row in update.work() {
        near(row.supplied_work_j, 0.5 * expected, 1e-24);
    }
    assert_eq!(update.circuit().charges_c(), charge);
    let reverse = update
        .circuit()
        .rebind_capacitance_blocks_candidate(vec![before], &[], || false)
        .unwrap();
    near(reverse.parameter_work_j(), -expected, 1e-24);
}

#[test]
fn cancellation_membership_source_authority_and_held_voltage_refuse_without_history_changes() {
    let before = block([0, 1], [2.0, 0.5, 0.5, 3.0]);
    let after = block([0, 1], [1.0, -0.25, -0.25, 2.0]);
    let original = simple(
        capacitors(3),
        vec![before.clone()],
        vec![2.0, -1.0, 0.25],
        &[],
    );
    let saved = original.clone();
    let mut visits = 0;
    let first = original
        .rebind_capacitance_blocks_candidate(vec![after.clone()], &[], || {
            visits += 1;
            false
        })
        .unwrap();
    assert!(visits > 1);
    let mut cancel_visit = 0;
    assert!(matches!(
        original.rebind_capacitance_blocks_candidate(vec![after.clone()], &[], || {
            cancel_visit += 1;
            cancel_visit == visits
        }),
        Err(MaterialCircuitError::Cancelled)
    ));
    assert_eq!(original, saved);
    let replay = original
        .rebind_capacitance_blocks_candidate(vec![after], &[], || false)
        .unwrap();
    assert_eq!(replay.circuit(), first.circuit());
    assert_eq!(replay.work(), first.work());
    assert_eq!(replay.parameter_work_j(), first.parameter_work_j());
    for replacement in [
        vec![],
        vec![block([0, 2], [2.0, 0.5, 0.5, 3.0])],
        vec![before.clone(), before.clone()],
    ] {
        assert!(
            original
                .rebind_capacitance_blocks_candidate(replacement, &[], || false)
                .is_err()
        );
        assert_eq!(original, saved);
    }

    let mut clamped_graph = capacitors(2);
    clamped_graph
        .branches
        .push((1, 0, Branch::VoltageSource { port: 0 }));
    let clamped = simple(
        clamped_graph,
        vec![block([0, 1], [1.0, 0.0, 0.0, 1.0])],
        vec![1.0, 0.0],
        &[1.0],
    );
    let held = clamped.clone();
    // At unchanged q this would jump the first voltage from 1 to 0.5 V. The
    // ideal source has supplied no impulse, so the native constraint refuses.
    assert!(
        clamped
            .rebind_capacitance_blocks_candidate(
                vec![block([0, 1], [2.0, 0.0, 0.0, 1.0])],
                &[1.0],
                || false
            )
            .is_err()
    );
    assert_eq!(clamped, held);

    let card = contact_source::scalar_card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        1.0,
    );
    let capacitor = ResolvedCapacitor::resolve(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
        1.0,
    )
    .unwrap();
    let source_bound = MaterialCircuit::new_with_storage(
        capacitors(2),
        vec![],
        (0..2)
            .map(|branch| MaterialStorage {
                branch,
                element: MaterialStorageElement::Capacitor(capacitor.clone()),
            })
            .collect(),
        vec![],
        vec![1.0, 0.0],
        &[],
    )
    .unwrap();
    let source_saved = source_bound.clone();
    assert!(
        source_bound
            .rebind_capacitance_blocks_candidate(
                vec![block([0, 1], [1.0, 0.25, 0.25, 1.0])],
                &[],
                || false
            )
            .is_err()
    );
    assert_eq!(source_bound, source_saved);

    let curved = MaterialCircuit::new_with_capacitance_blocks(
        capacitors(3),
        vec![],
        vec![],
        vec![],
        vec![before],
        vec![CapacitorCurve {
            branch: 2,
            curve: CapacitorChargeCurve::try_new(vec![(-2.0, -2.0), (0.0, 0.0), (2.0, 2.0)])
                .unwrap(),
        }],
        vec![],
        vec![2.0, -1.0, 0.25],
        &[],
    )
    .unwrap();
    let curved_saved = curved.clone();
    assert!(
        curved
            .rebind_capacitance_blocks_candidate(
                vec![block([0, 2], [2.0, 0.5, 0.5, 3.0])],
                &[],
                || false
            )
            .is_err()
    );
    assert_eq!(curved, curved_saved);
}
