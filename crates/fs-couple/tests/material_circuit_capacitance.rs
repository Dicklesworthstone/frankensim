//! Authored nonlinear storage with genuine material, contact and native DAE
//! owners. Synthetic coefficients verify physics/numerics, not measured fits.
#[path = "electrical_contacts/source.rs"]
mod contact_source;

use fs_couple::material_circuit::conductivity::{ConductivitySolve, ConductorResponse};
use fs_couple::material_circuit::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor};
use fs_material::state_point::MaterialPropertySelection;
use fs_phs::circuit::capacitance::{CapacitorChargeCurve, CapacitorCurve};
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};
use fs_qty::{QuantitySpec, Temperature};

fn card() -> MaterialCard {
    let dims = ELECTRICAL_RESISTIVITY_DIMS;
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("electrical_resistivity", dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 1.0), (400.0, 2.0)],
                dims,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                300.0,
                400.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic rho(T)".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic capacitor load".into(),
            phase: "solid".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resistor(card: &MaterialCard, branch: usize, t: f64) -> MaterialResistor {
    let point = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap();
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(
            card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}
fn curve(scale: f64) -> CapacitorChargeCurve {
    CapacitorChargeCurve::try_new(vec![
        (-4.0 * scale, -13.0),
        (-scale, -1.0),
        (0.0, 0.0),
        (scale, 1.0),
        (4.0 * scale, 13.0),
    ])
    .unwrap()
}
fn rc_graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 0, Branch::Resistor { ohms: 999.0 }),
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
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn nonlinear_rc_curve_survives_temperature_resistance_and_contact_rebuilds() {
    let card = card();
    let curves = vec![CapacitorCurve {
        branch: 0,
        curve: curve(1.0),
    }];
    let mut circuit = MaterialCircuit::new_with_capacitor_curves(
        rc_graph(),
        vec![resistor(&card, 1, 300.0)],
        vec![],
        vec![],
        curves.clone(),
        vec![],
        vec![2.0],
        &[],
    )
    .unwrap();
    let mut heat = 0.0;
    let mut exposure = 0.0_f64;
    let dt = 0.005;
    for k in 0..200 {
        circuit = circuit.rebind_contacts_candidate(vec![], &[]).unwrap();
        let t = if k < 100 { 300.0 } else { 400.0 };
        let source = resistor(&card, 1, t);
        exposure += dt / source.conductor.resistance_ohm();
        let proposal = circuit
            .step_candidate(vec![source.clone()], &[], dt, controls(), || false)
            .unwrap();
        assert_eq!(proposal.circuit().capacitor_curves(), curves);
        assert_eq!(proposal.step().resistors[0].material, source);
        let crossing = 5.0_f64.ln() / 4.0;
        let exact = if exposure <= crossing {
            0.75 + 1.25 * (-4.0 * exposure).exp()
        } else {
            (-(exposure - crossing)).exp()
        };
        near(proposal.circuit().charges_c()[0], exact, 4e-5);
        heat += proposal.step().heat_j;
        near(
            heat + curves[0]
                .curve
                .stored_energy_j(proposal.circuit().charges_c()[0])
                .unwrap(),
            3.5,
            2e-10,
        );
        circuit = proposal.into_parts().0;
    }
    assert!(
        circuit.charges_c()[0] < 1.0,
        "trajectory must cross the nonlinear knot"
    );
}

#[test]
fn nonlinear_capacitor_retains_separate_bulk_and_refreshed_contact_heat() {
    let card = card();
    let interface = contact_source::card(&card, &card, 2.0, 0.0);
    let contact = contact_source::contact(9, 1, &interface, 1.0);
    let curves = vec![CapacitorCurve {
        branch: 0,
        curve: curve(1.0),
    }];
    let circuit = MaterialCircuit::new_with_capacitor_curves(
        rc_graph(),
        vec![resistor(&card, 1, 300.0)],
        vec![contact.clone()],
        vec![],
        curves.clone(),
        vec![],
        vec![2.0],
        &[],
    )
    .unwrap();
    let step = circuit
        .step_with_contacts_candidate(
            circuit.resistors().to_vec(),
            vec![contact],
            &[],
            0.1,
            controls(),
            || false,
        )
        .unwrap();
    near(
        step.circuit().charges_c()[0],
        0.75 + 1.25 * (1.0 - 0.2 / 3.0) / (1.0 + 0.2 / 3.0),
        1e-12,
    );
    near(
        step.step().contacts[0].heat_j,
        2.0 * step.step().resistors[0].heat_j,
        1e-12,
    );
    let updated_interface = contact_source::card(&card, &card, 3.0, 0.0);
    let rebound = step
        .circuit()
        .rebind_contacts_candidate(
            vec![contact_source::contact(9, 1, &updated_interface, 1.0)],
            &[],
        )
        .unwrap();
    assert_eq!(rebound.charges_c(), step.circuit().charges_c());
    assert_eq!(rebound.capacitor_curves(), curves);
    let next = rebound
        .step_with_contacts_candidate(
            rebound.resistors().to_vec(),
            rebound.contacts().to_vec(),
            &[],
            0.1,
            controls(),
            || false,
        )
        .unwrap();
    near(
        next.step().contacts[0].heat_j,
        3.0 * next.step().resistors[0].heat_j,
        1e-12,
    );
    assert_eq!(next.circuit().capacitor_curves(), curves);
}

#[test]
fn refreshed_material_conductivity_trials_keep_the_same_nonlinear_storage() {
    let card = card();
    let curves = vec![CapacitorCurve {
        branch: 0,
        curve: curve(1.0),
    }];
    let circuit = MaterialCircuit::new_with_capacitor_curves(
        rc_graph(),
        vec![resistor(&card, 1, 300.0)],
        vec![],
        vec![],
        curves.clone(),
        vec![],
        vec![2.0],
        &[],
    )
    .unwrap();
    let accepted = circuit.clone();
    let candidate = circuit
        .conductivity_with_materials_candidate(
            vec![resistor(&card, 1, 400.0)],
            &[ConductorResponse {
                branch: 1,
                card: &card,
                axes: &[],
            }],
            &[],
            0.1,
            controls(),
            ConductivitySolve {
                maximum_iterations: 8,
                absolute_resistance_tolerance_ohm: 1e-12,
                relative_resistance_tolerance: 1e-12,
            },
            || false,
        )
        .unwrap();
    near(
        candidate.circuit().charges_c()[0],
        0.75 + 1.25 * 0.9 / 1.1,
        1e-12,
    );
    assert_eq!(candidate.circuit().capacitor_curves(), curves);
    assert_eq!(circuit, accepted);
}

#[test]
fn nonlinear_electric_energy_and_mutual_magnetic_energy_survive_material_trials() {
    let card = card();
    let curves = vec![CapacitorCurve {
        branch: 0,
        curve: curve(1.0),
    }];
    let mutual = vec![MutualInductance {
        branches: (1, 3),
        henries: 0.5,
    }];
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 7.0 }),
            (1, 0, Branch::Inductor { henries: 2.0 }),
            (1, 0, Branch::Resistor { ohms: 999.0 }),
            (2, 0, Branch::Inductor { henries: 3.0 }),
            (2, 0, Branch::Resistor { ohms: 999.0 }),
        ],
        transformers: vec![],
    };
    let energy = |c: &MaterialCircuit| {
        let x = c.fluxes_wb();
        curves[0].curve.stored_energy_j(c.charges_c()[0]).unwrap()
            + 0.5 * (3.0 * x[0] * x[0] - x[0] * x[1] + 2.0 * x[1] * x[1]) / 5.75
    };
    let mut circuit = MaterialCircuit::new_with_capacitor_curves(
        graph,
        vec![resistor(&card, 2, 300.0), resistor(&card, 4, 400.0)],
        vec![],
        mutual.clone(),
        curves.clone(),
        vec![1.0, 0.5],
        vec![2.0],
        &[],
    )
    .unwrap();
    let initial = energy(&circuit);
    let mut heat = 0.0;
    for _ in 0..40 {
        circuit = circuit.rebind_contacts_candidate(vec![], &[]).unwrap();
        let candidate = circuit
            .step_candidate(circuit.resistors().to_vec(), &[], 0.01, controls(), || {
                false
            })
            .unwrap();
        heat += candidate.step().heat_j;
        near(energy(candidate.circuit()) + heat, initial, 2e-10);
        assert_eq!(candidate.circuit().mutual_inductance(), mutual);
        assert_eq!(candidate.circuit().capacitor_curves(), curves);
        circuit = candidate.into_parts().0;
    }
}

#[test]
fn fixed_charge_law_replacement_reports_signed_integral_work_and_is_atomic() {
    for charge in [-2.0, 2.0] {
        let graph = CircuitGraph {
            node_count: 2,
            branches: vec![(1, 0, Branch::Capacitor { farads: 1.0 })],
            transformers: vec![],
        };
        let circuit = MaterialCircuit::new(graph, vec![], vec![], vec![charge], &[]).unwrap();
        let original = circuit.clone();
        let laws = vec![CapacitorCurve {
            branch: 0,
            curve: curve(1.0),
        }];
        let update = circuit
            .rebind_capacitor_curves_candidate(laws.clone(), &[], || false)
            .unwrap();
        assert_eq!(update.circuit().charges_c(), &[charge]);
        near(update.parameter_work_j(), 1.5, 1e-14);
        near(update.work()[0].stored_before_j, 2.0, 1e-14);
        near(update.work()[0].stored_after_j, 3.5, 1e-14);
        let reverse = update
            .circuit()
            .rebind_capacitor_curves_candidate(vec![], &[], || false)
            .unwrap();
        near(reverse.parameter_work_j(), -1.5, 1e-14);
        assert_eq!(reverse.circuit(), &original);
        let mut calls = 0;
        circuit
            .rebind_capacitor_curves_candidate(laws.clone(), &[], || {
                calls += 1;
                false
            })
            .unwrap();
        let mut seen = 0;
        assert!(matches!(
            circuit.rebind_capacitor_curves_candidate(laws.clone(), &[], || {
                seen += 1;
                seen == calls
            }),
            Err(MaterialCircuitError::Cancelled)
        ));
        assert_eq!(circuit, original);
        assert!(
            circuit
                .rebind_capacitor_curves_candidate(
                    vec![laws[0].clone(), laws[0].clone()],
                    &[],
                    || false
                )
                .is_err()
        );
        assert!(
            circuit
                .rebind_capacitor_curves_candidate(
                    vec![CapacitorCurve {
                        branch: 0,
                        curve: curve(0.1)
                    }],
                    &[],
                    || false
                )
                .is_err()
        );
    }
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
        ],
        transformers: vec![],
    };
    let clamped = MaterialCircuit::new(graph, vec![], vec![], vec![2.0], &[2.0]).unwrap();
    assert!(
        clamped
            .rebind_capacitor_curves_candidate(
                vec![CapacitorCurve {
                    branch: 0,
                    curve: curve(1.0)
                }],
                &[2.0],
                || false
            )
            .is_err()
    );
}

#[test]
fn parameter_work_is_resolved_below_large_endpoint_energy_spacing() {
    let q = 1e8;
    let before = CapacitorChargeCurve::try_new(vec![(-2e8, -2e8), (0.0, 0.0), (2e8, 2e8)]).unwrap();
    let after = CapacitorChargeCurve::try_new(vec![
        (-2e8, -2e8),
        (0.0, 0.0),
        (q - 1.0, q - 1.0),
        (q, q + 0.25),
        (2e8, 2e8 + 0.25),
    ])
    .unwrap();
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![(1, 0, Branch::Capacitor { farads: 1.0 })],
        transformers: vec![],
    };
    let circuit = MaterialCircuit::new_with_capacitor_curves(
        graph,
        vec![],
        vec![],
        vec![],
        vec![CapacitorCurve {
            branch: 0,
            curve: before,
        }],
        vec![],
        vec![q],
        &[],
    )
    .unwrap();
    let update = circuit
        .rebind_capacitor_curves_candidate(
            vec![CapacitorCurve {
                branch: 0,
                curve: after,
            }],
            &[],
            || false,
        )
        .unwrap();
    assert_eq!(
        update.work()[0].stored_before_j,
        update.work()[0].stored_after_j
    );
    near(update.parameter_work_j(), 0.125, 1e-14);
}

#[test]
fn retained_linear_material_authority_cannot_silently_describe_a_different_curve() {
    use fs_couple::material_circuit::storage::{MaterialStorage, MaterialStorageElement};
    use fs_material::field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, ResolvedCapacitor,
    };
    let dielectric = contact_source::scalar_card(
        ELECTRIC_PERMITTIVITY_PROPERTY,
        ELECTRIC_PERMITTIVITY_DIMS,
        1.0,
    );
    let source = ResolvedCapacitor::resolve(
        &dielectric,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
        1.0,
    )
    .unwrap();
    let graph = CircuitGraph {
        node_count: 2,
        branches: vec![(1, 0, Branch::Capacitor { farads: 1.0 })],
        transformers: vec![],
    };
    let circuit = MaterialCircuit::new_with_storage(
        graph,
        vec![],
        vec![MaterialStorage {
            branch: 0,
            element: MaterialStorageElement::Capacitor(source),
        }],
        vec![],
        vec![2.0],
        &[],
    )
    .unwrap();
    let saved = circuit.clone();
    assert!(
        circuit
            .rebind_capacitor_curves_candidate(
                vec![CapacitorCurve {
                    branch: 0,
                    curve: curve(1.0)
                }],
                &[],
                || false
            )
            .is_err()
    );
    assert_eq!(circuit, saved);
}

#[test]
fn picofarad_knot_crossing_survives_a_large_independent_stored_bias() {
    let card =
        contact_source::scalar_card("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, 1e6);
    let resistor = MaterialResistor {
        branch: 3,
        conductor: resolve_uniform_conductor(
            &card,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    };
    let law = curve(1e-12);
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
            (2, 0, Branch::Capacitor { farads: 1.0 }),
            (2, 0, Branch::Resistor { ohms: 999.0 }),
        ],
        transformers: vec![],
    };
    let circuit = MaterialCircuit::new_with_capacitor_curves(
        graph,
        vec![resistor],
        vec![],
        vec![],
        vec![CapacitorCurve {
            branch: 2,
            curve: law.clone(),
        }],
        vec![],
        vec![1e8, 2e-12],
        &[1e8],
    )
    .unwrap();
    let candidate = circuit
        .step_candidate(
            circuit.resistors().to_vec(),
            &[1e8],
            0.5e-6,
            MaterialCircuitSolve {
                energy_tolerance_j: 1e-22,
                power_tolerance_w: 1e-16,
            },
            || false,
        )
        .unwrap();
    let charge = candidate.circuit().charges_c()[1];
    assert!(
        charge > 0.0 && charge < 1e-12,
        "physical step must cross the curve knot"
    );
    assert_eq!(
        candidate.step().stored_before_j,
        candidate.step().stored_after_j
    );
    let change = law.energy_change_j(2e-12, charge).unwrap();
    near(candidate.step().stored_energy_change_j(), change, 1e-23);
    near(change + candidate.step().heat_j, 0.0, 1e-22);
    assert!(candidate.step().heat_j > 0.0);
}
