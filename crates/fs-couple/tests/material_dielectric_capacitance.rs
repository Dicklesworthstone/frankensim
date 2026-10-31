//! One original dielectric source supplies native stationary capacitor history.
//! Authored source values verify the bridge, not measured device accuracy.
use fs_couple::material_circuit::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve, MaterialResistor,
    capacitance::dielectric::MaterialDielectricCapacitor,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor},
    field_storage::dielectric::reversible::{
        ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS, REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
        ReversibleDielectricMaterial,
    },
    state_point::MaterialPropertySelection,
};
use fs_phs::circuit::capacitance::{CapacitorChargeCurve, CapacitorCurve};
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};
use fs_qty::{
    QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
fn card_with_value(value: PropertyValue) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, value, validity| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored stationary nonlinear dielectric fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    let validity = ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 400.0);
    add(
        REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
        ELECTRIC_DISPLACEMENT_DIMS,
        value,
        validity.clone().with_quantity(
            "E",
            QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
            0.0,
            13.0,
        ),
    );
    add(
        "electrical_resistivity",
        ELECTRICAL_RESISTIVITY_DIMS,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 1.0), (400.0, 2.0)],
            dims: ELECTRICAL_RESISTIVITY_DIMS,
        },
        validity,
    );
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic nonlinear dielectric".into(),
            phase: "solid".into(),
            process: "stationary source fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn card(scale: f64) -> MaterialCard {
    card_with_value(PropertyValue::Curve {
        abscissa: "E".into(),
        abscissa_dims: ELECTRIC_FIELD_DIMS,
        knots: vec![(0.0, 0.0), (1.0, scale), (13.0, 4.0 * scale)],
        dims: ELECTRIC_DISPLACEMENT_DIMS,
    })
}
fn material(card: &MaterialCard) -> ReversibleDielectricMaterial {
    ReversibleDielectricMaterial::resolve(
        card,
        &point(300.0)
            .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
            .unwrap(),
        "E",
        13.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn capacitor(
    card: &MaterialCard,
    branch: usize,
    area: f64,
    gap: f64,
) -> MaterialDielectricCapacitor {
    MaterialDielectricCapacitor::new(branch, material(card), area, gap).unwrap()
}
fn resistor(card: &MaterialCard, branch: usize, t: f64) -> MaterialResistor {
    MaterialResistor {
        branch,
        conductor: resolve_uniform_conductor(
            card,
            &point(t),
            MaterialPropertySelection::SingleClaimOnly,
            1.0,
            1.0,
        )
        .unwrap(),
    }
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 1e-9,
        power_tolerance_w: 1e-9,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}
fn rc(source: &MaterialCard) -> MaterialCircuit {
    MaterialCircuit::new_with_dielectric_capacitors(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::Resistor { ohms: f64::NAN }),
            ],
            transformers: vec![],
        },
        vec![resistor(source, 1, 300.0)],
        vec![],
        vec![],
        vec![],
        vec![capacitor(source, 0, 1.0, 1.0)],
        vec![],
        vec![2.0],
        &[],
        || false,
    )
    .unwrap()
}

#[test]
fn exact_geometry_mapping_retains_source_and_both_native_charge_polarities() {
    let source = card(1.0);
    let capacitor = capacitor(&source, 0, 0.2, 0.05);
    assert_eq!(
        capacitor.material().material().card_identity(),
        source.content_hash()
    );
    assert_eq!(capacitor.area_m2(), 0.2);
    assert_eq!(capacitor.gap_m(), 0.05);
    assert_eq!(capacitor.maximum_charge_c(), 0.8);
    assert_eq!(capacitor.maximum_voltage_v(), 0.65);
    for sign in [-1.0, 1.0] {
        near(capacitor.voltage_v(sign * 0.4).unwrap(), sign * 0.25, 1e-15);
        near(capacitor.charge_c(sign * 0.25).unwrap(), sign * 0.4, 1e-15);
        // w(D=2)=1/2+integral_1^2(4D-3)dD=3.5; volume=0.01.
        near(capacitor.stored_energy_j(sign * 0.4).unwrap(), 0.035, 1e-15);
    }
    assert!((capacitor.stored_energy_j(0.4).unwrap() - 0.5 * 0.4 * 0.25).abs() > 0.01);
    let copy = capacitor.clone();
    assert!(std::ptr::eq(
        copy.curve().knots(),
        capacitor.curve().knots()
    ));
    assert!(capacitor.voltage_v(0.81).is_err());
}

#[test]
fn native_rc_discharge_crosses_source_knot_and_preserves_source_through_rebuilds() {
    let source = card(1.0);
    let mut circuit = rc(&source);
    let original = circuit.dielectric_capacitors().to_vec();
    let mut exposure = 0.0;
    let mut heat = 0.0;
    for i in 0..120 {
        circuit = circuit.rebind_contacts_candidate(vec![], &[]).unwrap();
        let t = if i < 60 { 300.0 } else { 400.0 };
        let resistor = resistor(&source, 1, t);
        exposure += 0.01 / resistor.conductor.resistance_ohm();
        let candidate = circuit
            .step_candidate(vec![resistor], &[], 0.01, controls(), || false)
            .unwrap();
        assert_eq!(candidate.circuit().dielectric_capacitors(), original);
        let crossing = 5.0_f64.ln() / 4.0;
        let exact = if exposure <= crossing {
            0.75 + 1.25 * (-4.0 * exposure).exp()
        } else {
            (-(exposure - crossing)).exp()
        };
        near(candidate.circuit().charges_c()[0], exact, 2e-4);
        heat += candidate.step().heat_j;
        near(
            heat + original[0]
                .stored_energy_j(candidate.circuit().charges_c()[0])
                .unwrap(),
            3.5,
            1e-9,
        );
        circuit = candidate.into_parts().0;
    }
    assert!(circuit.charges_c()[0] < 1.0);
}

#[test]
fn source_rebind_retains_native_histories_other_curves_and_mutual_storage_with_signed_work() {
    let source = card(1.0);
    let replacement = card(2.0);
    let authored = CapacitorCurve {
        branch: 3,
        curve: CapacitorChargeCurve::try_new(vec![(-2.0, -2.0), (0.0, 0.0), (2.0, 2.0)]).unwrap(),
    };
    let mutual = vec![MutualInductance {
        branches: (1, 2),
        henries: 0.4,
    }];
    let circuit = MaterialCircuit::new_with_dielectric_capacitors(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (2, 0, Branch::Inductor { henries: 2.0 }),
                (2, 0, Branch::Inductor { henries: 3.0 }),
                (2, 0, Branch::Capacitor { farads: 10.0 }),
                (1, 0, Branch::Resistor { ohms: 1.0 }),
                (2, 0, Branch::Resistor { ohms: 1.0 }),
            ],
            transformers: vec![],
        },
        vec![resistor(&source, 4, 300.0), resistor(&source, 5, 300.0)],
        vec![],
        mutual.clone(),
        vec![authored.clone()],
        vec![capacitor(&source, 0, 1.0, 1.0)],
        vec![0.2, -0.1],
        vec![2.0, 0.25],
        &[],
        || false,
    )
    .unwrap();
    let accepted = circuit.clone();
    let refreshed = capacitor(&replacement, 0, 1.0, 1.0);
    let candidate = circuit
        .rebind_dielectric_capacitors_candidate(vec![refreshed.clone()], &[], || false)
        .unwrap();
    assert_eq!(candidate.circuit().charges_c(), circuit.charges_c());
    assert_eq!(candidate.circuit().fluxes_wb(), circuit.fluxes_wb());
    assert_eq!(candidate.circuit().mutual_inductance(), mutual);
    assert_eq!(candidate.circuit().resistors(), circuit.resistors());
    assert_eq!(candidate.circuit().capacitor_curves()[1], authored);
    assert_eq!(candidate.circuit().dielectric_capacitors(), &[refreshed]);
    // At retained q=2 C: original U=3.5 J; the doubled D(E) source gives U=1 J.
    near(candidate.parameter_work_j(), -2.5, 1e-14);
    near(candidate.work()[0].stored_before_j, 3.5, 1e-14);
    near(candidate.work()[0].stored_after_j, 1.0, 1e-14);
    assert_eq!(candidate.work()[1].supplied_work_j, 0.0);
    let reverse = candidate
        .circuit()
        .rebind_dielectric_capacitors_candidate(
            circuit.dielectric_capacitors().to_vec(),
            &[],
            || false,
        )
        .unwrap();
    near(reverse.parameter_work_j(), 2.5, 1e-14);
    assert_eq!(reverse.circuit().charges_c(), circuit.charges_c());
    let unchanged = circuit
        .rebind_capacitor_curves_candidate(circuit.capacitor_curves().to_vec(), &[], || false)
        .unwrap();
    assert_eq!(
        unchanged.circuit().dielectric_capacitors(),
        circuit.dielectric_capacitors()
    );
    assert_eq!(unchanged.parameter_work_j(), 0.0);
    assert!(
        circuit
            .rebind_capacitor_curves_candidate(vec![authored], &[], || false)
            .is_err()
    );
    assert!(
        circuit
            .rebind_dielectric_capacitors_candidate(vec![], &[], || false)
            .is_err()
    );
    assert!(
        circuit
            .rebind_dielectric_capacitors_candidate(
                vec![capacitor(&replacement, 0, 2.0, 1.0)],
                &[],
                || false
            )
            .is_err()
    );
    assert_eq!(circuit, accepted);
}

#[test]
fn clamp_and_cancellation_refuse_atomically_and_unowned_thermal_dependence_is_not_frozen() {
    let source = card(1.0);
    let replacement = card(2.0);
    let clamped = MaterialCircuit::new_with_dielectric_capacitors(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        vec![],
        vec![],
        vec![],
        vec![],
        vec![capacitor(&source, 0, 1.0, 1.0)],
        vec![],
        vec![2.0],
        &[5.0],
        || false,
    )
    .unwrap();
    let accepted_clamp = clamped.clone();
    assert!(
        clamped
            .rebind_dielectric_capacitors_candidate(
                vec![capacitor(&replacement, 0, 1.0, 1.0)],
                &[5.0],
                || false
            )
            .is_err()
    );
    assert_eq!(clamped, accepted_clamp);
    let circuit = rc(&source);
    let accepted = circuit.clone();
    let next = capacitor(&replacement, 0, 1.0, 1.0);
    let mut calls = 0;
    circuit
        .rebind_dielectric_capacitors_candidate(vec![next.clone()], &[], || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls] {
        let mut observed = 0;
        assert!(matches!(
            circuit.rebind_dielectric_capacitors_candidate(vec![next.clone()], &[], || {
                observed += 1;
                observed >= stop
            }),
            Err(MaterialCircuitError::Cancelled)
        ));
        assert_eq!(circuit, accepted);
    }
    let dependent = card_with_value(PropertyValue::Surface2D {
        axes: ["E".into(), "T".into()],
        axis_dims: [ELECTRIC_FIELD_DIMS, Temperature::DIMS],
        samples: vec![
            ([0.0, 300.0], 0.0),
            ([13.0, 300.0], 4.0),
            ([0.0, 400.0], 0.0),
            ([13.0, 400.0], 8.0),
        ],
        triangles: vec![[0, 1, 3], [0, 3, 2]],
        dims: ELECTRIC_DISPLACEMENT_DIMS,
    });
    assert!(MaterialDielectricCapacitor::new(0, material(&dependent), 1.0, 1.0).is_err());
    assert!(MaterialDielectricCapacitor::new(0, material(&source), 0.0, 1.0).is_err());
}
