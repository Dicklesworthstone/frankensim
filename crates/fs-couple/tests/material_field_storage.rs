//! G1/G3 source-derived field storage on the real Kirchhoff circuit owner.
use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
use fs_couple::material_circuit::{MaterialCircuit, MaterialCircuitSolve};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::field_storage::*;
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};
use fs_qty::Dims;

fn card(property: &str, dims: Dims, value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(property, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic linear circuit material".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic field storage".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn capacitor(branch: usize, epsilon: f64) -> MaterialStorage {
    let source = card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        epsilon,
    );
    MaterialStorage {
        branch,
        element: MaterialStorageElement::Capacitor(
            ResolvedCapacitor::resolve(
                &source,
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                1.0,
                1.0,
            )
            .unwrap(),
        ),
    }
}
fn inductor(branch: usize, mu: f64) -> MaterialStorage {
    let source = card(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        mu,
    );
    MaterialStorage {
        branch,
        element: MaterialStorageElement::Inductor(
            ResolvedInductor::resolve(
                &source,
                &QueryPoint::new(),
                MaterialPropertySelection::SingleClaimOnly,
                1.0,
                1.0,
                1,
            )
            .unwrap(),
        ),
    }
}
fn graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: f64::NAN }),
            (1, 0, Branch::Inductor { henries: f64::NAN }),
        ],
        transformers: vec![],
    }
}
fn lc(c: f64, l: f64, flux: f64, charge: f64) -> MaterialCircuit {
    // Reverse source order deliberately; source/coordinate order is canonical.
    MaterialCircuit::new_with_storage(
        graph(),
        vec![],
        vec![inductor(1, l), capacitor(0, c)],
        vec![flux],
        vec![charge],
        &[],
    )
    .unwrap()
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 1e-10,
        power_tolerance_w: 1e-10,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

#[test]
fn data_only_material_swaps_change_the_actual_lc_oscillation() {
    let (c, l, q) = (0.002_f64, 0.05_f64, 0.02_f64);
    let quarter = 0.5 * core::f64::consts::PI * (c * l).sqrt();
    let dt = quarter / 512.0;
    for (cap_factor, core_factor) in [(1.0_f64, 1.0_f64), (2.0, 1.0), (1.0, 2.0)] {
        let mut circuit = lc(c * cap_factor, l * core_factor, 0.0, q);
        let sources = circuit.storage().to_vec();
        let h0 = 0.5 * q * q / (c * cap_factor);
        for _ in 0..512 {
            let candidate = circuit
                .step_candidate(vec![], &[], dt, controls(), || false)
                .unwrap();
            near(candidate.step().stored_after_j, h0, 1e-9);
            assert_eq!(candidate.step().heat_j, 0.0);
            assert_eq!(candidate.circuit().storage(), &sources);
            circuit = candidate.into_parts().0;
        }
        // Independent continuous LC solution, with midpoint's O(dt^2) error.
        let angle = quarter / (c * cap_factor * l * core_factor).sqrt();
        near(circuit.charges_c()[0], q * angle.cos(), 3e-8);
    }
}

#[test]
fn live_storage_update_preserves_coordinates_and_accounts_parameter_work() {
    let original = lc(0.002, 0.05, 0.01, 0.02);
    let before = original.clone();
    let change = original
        .rebind_storage_candidate(vec![capacitor(0, 0.004), inductor(1, 0.1)], &[], || false)
        .unwrap();
    assert_eq!(change.circuit().charges_c(), original.charges_c());
    assert_eq!(change.circuit().fluxes_wb(), original.fluxes_wb());
    assert_eq!(original, before);
    let h0 = 0.5 * 0.02_f64.powi(2) / 0.002 + 0.5 * 0.01_f64.powi(2) / 0.05;
    near(change.parameter_work_j(), -0.5 * h0, 1e-14);
    for work in change.work() {
        near(work.stored_after_j, 0.5 * work.stored_before_j, 1e-14);
        near(
            work.supplied_work_j,
            work.stored_after_j - work.stored_before_j,
            1e-14,
        );
    }
    let (circuit, _, work) = change.into_parts();
    let next = circuit
        .step_candidate(vec![], &[], 0.001, controls(), || false)
        .unwrap();
    near(next.step().stored_before_j, h0 + work, 1e-14);
    near(next.step().stored_after_j, h0 + work, 1e-10);
    let reverse = circuit
        .rebind_storage_candidate(original.storage().to_vec(), &[], || false)
        .unwrap();
    near(reverse.parameter_work_j(), -work, 1e-14);
    assert_eq!(reverse.circuit(), &original);
}

#[test]
fn shared_magnetic_material_update_uses_full_inverse_inductance_parameter_work() {
    let graph = CircuitGraph { node_count:3, branches:vec![
        (1,0,Branch::Inductor {henries:f64::NAN}),
        (2,0,Branch::Inductor {henries:f64::NAN}),
        (1,0,Branch::VoltageSource {port:0}),
        (2,0,Branch::VoltageSource {port:1}),
    ], transformers:vec![] };
    let links = vec![MutualInductance {branches:(0,1),henries:1.0}];
    let original = MaterialCircuit::new_with_storage_and_mutual_inductance(graph,vec![],
        vec![inductor(0,2.0),inductor(1,3.0)],links.clone(),vec![1.0,4.0],vec![],&[0.0;2]).unwrap();
    let saved = original.clone();
    let change = original.rebind_storage_candidate(vec![inductor(0,4.0),inductor(1,5.0)],
        &[0.0;2],||false).unwrap();
    // L^-1=(L1 L2-M^2)^-1 [[L2,-M],[-M,L1]], at fixed total lambda=(1,4).
    let energy = |l1:f64,l2:f64| 0.5*(l2-8.0+16.0*l1)/(l1*l2-1.0);
    let h0 = energy(2.0,3.0);
    let h1 = energy(4.0,5.0);
    near(change.parameter_work_j(),h1-h0,1e-13);
    near(change.work().iter().map(|w|w.stored_before_j).sum(),h0,1e-13);
    near(change.work().iter().map(|w|w.stored_after_j).sum(),h1,1e-13);
    // Individual shared-field contributions need not be positive.
    near(change.work()[0].stored_before_j,-0.1,1e-13);
    assert_eq!(change.circuit().mutual_inductance(),links);
    assert_eq!(change.circuit().fluxes_wb(),original.fluxes_wb());
    let next = change.circuit().step_candidate(vec![],&[0.0;2],0.01,controls(),||false).unwrap();
    near(next.step().stored_before_j,h1,1e-13);
    near(next.step().stored_after_j,h1,1e-13);
    let reverse = change.circuit().rebind_storage_candidate(original.storage().to_vec(),
        &[0.0;2],||false).unwrap();
    near(reverse.parameter_work_j(),h0-h1,1e-13);
    assert_eq!(reverse.circuit(),&original);
    // A self-L refresh may destroy positive definiteness at held M; refuse it.
    assert!(original.rebind_storage_candidate(vec![inductor(0,0.1),inductor(1,0.1)],
        &[0.0;2],||false).is_err());
    assert_eq!(original,saved);
}

#[test]
fn inconsistent_ideal_source_and_invalid_binding_refuse_without_resetting_state() {
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 99.0 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
        ],
        transformers: vec![],
    };
    let clamped = MaterialCircuit::new_with_storage(
        graph,
        vec![],
        vec![capacitor(0, 0.002)],
        vec![],
        vec![0.02],
        &[10.0],
    )
    .unwrap();
    let before = clamped.clone();
    assert!(
        clamped
            .rebind_storage_candidate(vec![capacitor(0, 0.004)], &[10.0], || false)
            .is_err()
    );
    assert_eq!(clamped, before);

    let circuit = lc(0.002, 0.05, 0.01, 0.02);
    for bindings in [
        vec![],
        vec![capacitor(0, 0.002)],
        vec![capacitor(0, 0.002), capacitor(0, 0.002)],
        vec![inductor(0, 0.05), capacitor(1, 0.002)],
        vec![capacitor(0, 0.002), inductor(99, 0.05)],
    ] {
        assert!(
            circuit
                .rebind_storage_candidate(bindings, &[], || false)
                .is_err()
        );
    }
    let mut checks = 0;
    circuit
        .rebind_storage_candidate(circuit.storage().to_vec(), &[], || {
            checks += 1;
            false
        })
        .unwrap();
    let mut now = 0;
    assert!(
        circuit
            .rebind_storage_candidate(circuit.storage().to_vec(), &[], || {
                now += 1;
                now == checks
            })
            .is_err()
    );
    assert_eq!(circuit, lc(0.002, 0.05, 0.01, 0.02));
}
