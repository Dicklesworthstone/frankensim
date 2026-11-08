//! Independent branch-work checks for source-bound cycle-mean AC heat.
use super::{controls, material};
use fs_couple::material_circuit::{
    MaterialCircuitError, MaterialCircuitSolve, harmonic::account_harmonic_resistor_heat,
};
use fs_em::planar::eddy::harmonic::circuit::HarmonicCircuitBranch;
use fs_math::c64::C64;
use fs_phs::circuit::{Branch, CircuitGraph};

fn graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 2,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 0, Branch::Resistor { ohms: 2.0 }),
            (1, 0, Branch::Inductor { henries: 1.0 }),
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (1, 0, Branch::Resistor { ohms: 4.0 }),
        ],
        transformers: vec![],
    }
}

// Parallel RLC at omega=1 under exp(-i omega t). The storage currents cancel,
// while their individual reactive work remains nonzero and is not heat. These
// independently specified exact efforts are deliberately in mixed graph order.
fn branches(voltage: C64) -> Vec<HarmonicCircuitBranch> {
    [
        (4, voltage.scale(0.25)),
        (2, C64::new(-voltage.im, voltage.re)),
        (0, voltage.scale(-0.75)),
        (3, C64::new(voltage.im, -voltage.re)),
        (1, voltage.scale(0.5)),
    ]
    .into_iter()
    .map(|(branch, current_a)| HarmonicCircuitBranch {
        branch,
        voltage_v: voltage,
        current_a,
        absorbed_power_va: (voltage * current_a.conj()).scale(0.5),
    })
    .collect()
}

#[test]
fn peak_phasors_preserve_source_identity_phase_and_original_resistor_order() {
    let graph = graph();
    let sources = [material(4, 4.0), material(1, 2.0)];
    for voltage in [C64::new(6.0, 8.0), C64::new(-8.0, 6.0)] {
        let branches = branches(voltage);
        let (heat, total) = account_harmonic_resistor_heat(
            &graph,
            &sources,
            &branches,
            37.5,
            0.4,
            controls(),
            1e-12,
            &mut || false,
        )
        .unwrap();
        assert_eq!(heat.len(), 2);
        assert_eq!(heat[0].material, sources[1]);
        assert_eq!(heat[1].material, sources[0]);
        assert_eq!(heat[0].voltage_v, voltage);
        assert_eq!(heat[0].current_a, voltage.scale(0.5));
        assert_eq!(heat[1].current_a, voltage.scale(0.25));
        assert_eq!(heat[0].mean_power_w, 25.0);
        assert_eq!(heat[1].mean_power_w, 12.5);
        assert_eq!(heat[0].heat_j, 10.0);
        assert_eq!(heat[1].heat_j, 5.0);
        assert_eq!(total, 15.0);
        assert_eq!(branches[1].absorbed_power_va.im, -50.0);
        assert_eq!(branches[3].absorbed_power_va.im, 50.0);
        assert_eq!(-0.4 * branches[2].absorbed_power_va.re, total);
    }
}

#[test]
fn resistor_heat_does_not_reinterpret_native_field_storage_placeholders() {
    let mut graph = graph();
    graph.branches[2].2 = Branch::Inductor { henries: f64::NAN };
    let sources = [material(1, 2.0), material(4, 4.0)];
    let (heat, total) = account_harmonic_resistor_heat(
        &graph,
        &sources,
        &branches(C64::new(6.0, 8.0)),
        37.5,
        0.4,
        controls(),
        1e-12,
        &mut || false,
    )
    .unwrap();
    assert_eq!(heat.len(), 2);
    assert_eq!(total, 15.0);
}

#[test]
fn instantaneous_field_or_current_response_cannot_be_frozen_as_ac_resistance() {
    use super::{
        ClaimSet, ELECTRICAL_RESISTIVITY_DIMS, InterpolationPolicy, MaterialCard,
        MaterialPropertySelection, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, QueryPoint, UncertaintyModel, ValidityDomain, resolve_uniform_conductor,
    };
    use fs_couple::material_circuit::MaterialResistor;
    use fs_qty::{Dims, QuantitySpec};

    for axis_dims in [Dims([1, 1, -3, 0, -1, 0]), Dims([-2, 0, 0, 0, 1, 0])] {
        for response in [false, true] {
            let dims = ELECTRICAL_RESISTIVITY_DIMS;
            let mut claims = ClaimSet::new();
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new("electrical_resistivity", dims),
                    value: if response {
                        PropertyValue::Curve {
                            abscissa: "drive".into(),
                            abscissa_dims: axis_dims,
                            knots: vec![(0.0, 2.0), (10.0, 3.0)],
                            dims,
                        }
                    } else {
                        PropertyValue::Scalar { value: 2.0, dims }
                    },
                    validity: ValidityDomain::unconstrained().with_quantity(
                        "drive",
                        QuantitySpec::dimensional(axis_dims),
                        0.0,
                        10.0,
                    ),
                    interpolation: if response {
                        InterpolationPolicy::LinearInside
                    } else {
                        InterpolationPolicy::ConstantWithinValidity
                    },
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic operating-dependent AC resistor".into(),
                        license: "CC0-1.0".into(),
                        artifact: None,
                    },
                })
                .unwrap();
            let card = MaterialCard::assemble(
                MaterialStateId {
                    chemistry: "synthetic nonlinear resistor".into(),
                    phase: "solid".into(),
                    process: "numerical fixture".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap();
            let point = QueryPoint::new()
                .with_quantity("drive", QuantitySpec::dimensional(axis_dims), 0.0)
                .unwrap();
            let conductor = resolve_uniform_conductor(
                &card,
                &point,
                MaterialPropertySelection::SingleClaimOnly,
                1.0,
                1.0,
            )
            .unwrap();
            assert!(matches!(
                account_harmonic_resistor_heat(
                    &graph(),
                    &[
                        MaterialResistor {
                            branch: 1,
                            conductor
                        },
                        material(4, 4.0)
                    ],
                    &branches(C64::new(6.0, 8.0)),
                    37.5,
                    0.4,
                    controls(),
                    1e-12,
                    &mut || false,
                ),
                Err(MaterialCircuitError::Invalid(
                    "AC resistor heat requires an instantaneous scalar source"
                ))
            ));
        }
    }
}

#[test]
fn changed_material_current_complex_work_or_native_loss_cannot_be_hidden() {
    let graph = graph();
    let sources = [material(1, 2.0), material(4, 4.0)];
    let original = branches(C64::new(6.0, 8.0));
    let account = |graph: &CircuitGraph, records: &[HarmonicCircuitBranch], power| {
        account_harmonic_resistor_heat(
            graph,
            &sources,
            records,
            power,
            0.4,
            controls(),
            1e-12,
            &mut || false,
        )
    };
    let mut wrong_current = original.clone();
    // Same magnitude and source-owned Joule heat, but the wrong phase/current.
    wrong_current[4].current_a.im = -wrong_current[4].current_a.im;
    assert!(matches!(
        account(&graph, &wrong_current, 37.5),
        Err(MaterialCircuitError::Balance(_))
    ));
    for addition in [C64::from_re(0.25), C64::new(0.0, 0.25)] {
        let mut wrong_power = original.clone();
        wrong_power[4].absorbed_power_va = wrong_power[4].absorbed_power_va + addition;
        assert!(matches!(
            account(&graph, &wrong_power, 37.5),
            Err(MaterialCircuitError::Balance(_))
        ));
    }
    assert!(matches!(
        account(&graph, &original, 38.5),
        Err(MaterialCircuitError::Balance(_))
    ));
    let mut wrong_graph = graph.clone();
    wrong_graph.branches[1].2 = Branch::Resistor { ohms: 3.0 };
    assert!(matches!(
        account(&wrong_graph, &original, 37.5),
        Err(MaterialCircuitError::Invalid(_))
    ));
    let mut repeated_branch = original.clone();
    repeated_branch[0] = repeated_branch[4].clone();
    assert!(matches!(
        account(&graph, &repeated_branch, 37.5),
        Err(MaterialCircuitError::Invalid(_))
    ));
    assert!(matches!(
        account(&graph, &original[..4], 37.5),
        Err(MaterialCircuitError::Invalid(_))
    ));
    assert!(
        account_harmonic_resistor_heat(
            &graph,
            &[sources[0].clone(), sources[0].clone()],
            &original,
            37.5,
            0.4,
            controls(),
            1e-12,
            &mut || false,
        )
        .is_err()
    );
}

#[test]
fn interval_controls_and_late_cancellation_refuse_without_consuming_sources() {
    let graph = graph();
    let sources = [material(1, 2.0), material(4, 4.0)];
    let records = branches(C64::new(6.0, 8.0));
    let account = |dt, controls, cancel: &mut dyn FnMut() -> bool| {
        account_harmonic_resistor_heat(
            &graph,
            &sources,
            &records,
            37.5,
            dt,
            controls,
            1e-12,
            &mut || cancel(),
        )
    };
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            account(dt, controls(), &mut || false),
            Err(MaterialCircuitError::Invalid(_))
        ));
    }
    assert!(matches!(
        account(
            0.4,
            MaterialCircuitSolve {
                power_tolerance_w: f64::NAN,
                ..controls()
            },
            &mut || false,
        ),
        Err(MaterialCircuitError::Invalid(_))
    ));
    let mut count = 0;
    account(0.4, controls(), &mut || {
        count += 1;
        false
    })
    .unwrap();
    for stop in [1, count] {
        let mut calls = 0;
        assert!(matches!(
            account(0.4, controls(), &mut || {
                calls += 1;
                calls == stop
            }),
            Err(MaterialCircuitError::Cancelled)
        ));
    }
    let (heat, total) = account(0.4, controls(), &mut || false).unwrap();
    assert_eq!(total, 15.0);
    assert_eq!(heat[0].material, sources[0]);
    assert_eq!(heat[1].material, sources[1]);
}
