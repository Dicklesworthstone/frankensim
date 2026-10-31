//! Explicit source work preserves nonlinear instantaneous and Debye histories.
//! All source data below are authored numerical fixtures, not measurements.
use fs_couple::material_circuit::{
    MaterialCircuitError, MaterialResistor, capacitance::dielectric::MaterialDielectricCapacitor,
    dielectric::DebyeCircuit,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ELECTRICAL_RESISTIVITY_DIMS, resolve_uniform_conductor},
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::{
            DebyePoleProperties, ResolvedDebyeCapacitor,
            reversible::{
                ELECTRIC_DISPLACEMENT_DIMS, ELECTRIC_FIELD_DIMS,
                REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY, ReversibleDielectricMaterial,
            },
        },
    },
    state_point::{MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement},
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{
    QuantitySpec, Temperature, Time,
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
fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 400.0)
}
fn names() -> Vec<DebyePoleProperties> {
    (0..2)
        .map(|i| DebyePoleProperties {
            permittivity_increment: format!("chi_{i}"),
            relaxation_time: format!("tau_{i}"),
        })
        .collect()
}
fn card(scale: f64, increments: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, value, interpolation, validity| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored nonlinear source-work fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(
        "epsilon_inf",
        ELECTRIC_PERMITTIVITY_DIMS,
        PropertyValue::Scalar {
            value: scale,
            dims: ELECTRIC_PERMITTIVITY_DIMS,
        },
        InterpolationPolicy::ConstantWithinValidity,
        validity(),
    );
    add(
        REVERSIBLE_ELECTRIC_DISPLACEMENT_PROPERTY,
        ELECTRIC_DISPLACEMENT_DIMS,
        PropertyValue::Curve {
            abscissa: "E".into(),
            abscissa_dims: ELECTRIC_FIELD_DIMS,
            knots: vec![
                (0.0, 0.0),
                (1.0, scale),
                (4.0, 2.0 * scale),
                (8.0, 4.0 * scale),
            ],
            dims: ELECTRIC_DISPLACEMENT_DIMS,
        },
        InterpolationPolicy::LinearInside,
        validity().with_quantity(
            "E",
            QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS),
            0.0,
            8.0,
        ),
    );
    for (i, increment) in increments.into_iter().enumerate() {
        add(
            &format!("chi_{i}"),
            ELECTRIC_PERMITTIVITY_DIMS,
            PropertyValue::Scalar {
                value: increment,
                dims: ELECTRIC_PERMITTIVITY_DIMS,
            },
            InterpolationPolicy::ConstantWithinValidity,
            validity(),
        );
        let tau = [0.2, 0.7][i];
        add(
            &format!("tau_{i}"),
            Time::DIMS,
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, tau), (400.0, 2.0 * tau)],
                dims: Time::DIMS,
            },
            InterpolationPolicy::LinearInside,
            validity(),
        );
    }
    add(
        "electrical_resistivity",
        ELECTRICAL_RESISTIVITY_DIMS,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 1.0), (400.0, 3.0)],
            dims: ELECTRICAL_RESISTIVITY_DIMS,
        },
        InterpolationPolicy::LinearInside,
        validity(),
    );
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: format!("explicit instantaneous scale {scale}"),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(card: &MaterialCard, t: f64) -> ResolvedDebyeCapacitor {
    let source = ResolvedDebyeCapacitor::resolve(
        card,
        &point(t),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names(),
        1.0,
        1.0,
    )
    .unwrap();
    // Carry an additional source claim to catch bundle-truncating rebuilds.
    let material = source
        .material()
        .with_additional_requirements(
            card,
            &[ScalarPropertyRequirement::try_new(
                "electrical_resistivity",
                ELECTRICAL_RESISTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    ResolvedDebyeCapacitor::try_from_material(material, "epsilon_inf", &names(), 1.0, 1.0).unwrap()
}
fn instantaneous(card: &MaterialCard, t: f64, maximum_field: f64) -> MaterialDielectricCapacitor {
    let point = point(t)
        .with_quantity("E", QuantitySpec::dimensional(ELECTRIC_FIELD_DIMS), 0.0)
        .unwrap();
    MaterialDielectricCapacitor::new(
        0,
        ReversibleDielectricMaterial::resolve(
            card,
            &point,
            "E",
            maximum_field,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        1.0,
        1.0,
    )
    .unwrap()
}
fn resistors(card: &MaterialCard, t: f64) -> Vec<MaterialResistor> {
    [2, 4]
        .into_iter()
        .map(|branch| MaterialResistor {
            branch,
            conductor: resolve_uniform_conductor(
                card,
                &point(t),
                MaterialPropertySelection::SingleClaimOnly,
                5.0,
                1.0,
            )
            .unwrap(),
        })
        .collect()
}
fn circuit(card: &MaterialCard, charge: f64, unrelated_charge: f64) -> DebyeCircuit {
    DebyeCircuit::new_with_reversible_instantaneous(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (2, 0, Branch::Inductor { henries: 4.0 }),
                (2, 0, Branch::Resistor { ohms: 5.0 }),
                (2, 0, Branch::Capacitor { farads: 0.5 }),
                (1, 0, Branch::Resistor { ohms: 10.0 }),
            ],
            transformers: vec![],
        },
        source(card, 300.0),
        instantaneous(card, 300.0, 8.0),
        vec![0.2],
        vec![charge, unrelated_charge],
        vec![0.6, -0.2],
        &[],
        || false,
    )
    .unwrap()
}
fn bits(circuit: &DebyeCircuit) -> Vec<u64> {
    circuit.state().iter().map(|v| v.to_bits()).collect()
}
fn same_history(before: &DebyeCircuit, after: &DebyeCircuit) {
    assert_eq!(before.dae().charge_index, after.dae().charge_index);
    assert_eq!(before.dae().flux_index, after.dae().flux_index);
    for &i in before
        .dae()
        .charge_index
        .iter()
        .chain(&before.dae().flux_index)
    {
        assert_eq!(before.state()[i].to_bits(), after.state()[i].to_bits());
    }
    assert_eq!(before.graph().node_count, after.graph().node_count);
    assert_eq!(before.graph().branches.len(), after.graph().branches.len());
    assert_eq!(before.graph().transformers, after.graph().transformers);
    for (a, b) in before.graph().branches.iter().zip(&after.graph().branches) {
        assert_eq!((a.0, a.1), (b.0, b.1));
    }
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn explicit_source_pair_changes_storage_with_signed_work_and_exact_native_history() {
    let old_card = card(1.0, [2.0, 0.5]);
    let new_card = card(2.0, [1.0, 1.0]);
    for (charge, unrelated) in [(2.5, -0.25), (-2.5, 1e12)] {
        let initial = circuit(&old_card, charge, unrelated);
        let (before, zero) = initial
            .rebind_reversible_sources_and_resistors_candidate(
                source(&old_card, 300.0),
                instantaneous(&old_card, 300.0, 8.0),
                resistors(&old_card, 300.0),
                &[],
                || false,
            )
            .unwrap()
            .into_parts();
        assert_eq!(zero, 0.0);
        let accepted = bits(&before);
        let next_source = source(&new_card, 300.0);
        let next_law = instantaneous(&new_card, 300.0, 8.0);
        let candidate = before
            .rebind_reversible_sources_candidate(next_source.clone(), next_law.clone(), &[], || {
                false
            })
            .unwrap();
        let after = candidate.circuit();
        same_history(&before, after);
        assert_eq!(after.source(), &next_source);
        assert_eq!(after.reversible_instantaneous(), Some(&next_law));
        assert_eq!(after.material_resistors(), before.material_resistors());
        assert_eq!(after.source().material().properties().len(), 6);
        // Original V(q) segments integrate to 5.25 at |q|=2.5.
        // Doubling D(E) gives 2*[1/2 + 1/4 + 3/2*(1/4)^2].
        // Pole work independently equals .6²/2*(1-1/2)+.2²/2*(1-2).
        let instantaneous_work = 1.6875 - 5.25;
        let polarization_work = 0.6_f64.powi(2) * 0.25 - 0.2_f64.powi(2) * 0.5;
        near(
            candidate.parameter_work_j(),
            instantaneous_work + polarization_work,
        );
        near(before.dielectric_energy_j().unwrap(), 5.38);
        near(after.dielectric_energy_j().unwrap(), 1.8875);
        near(after.terminal_voltage_v(), charge.signum() * 1.75);
        let reverse = after
            .rebind_reversible_sources_candidate(
                source(&old_card, 300.0),
                instantaneous(&old_card, 300.0, 8.0),
                &[],
                || false,
            )
            .unwrap();
        near(reverse.parameter_work_j(), -candidate.parameter_work_j());
        same_history(&before, reverse.circuit());
        assert_eq!(bits(&before), accepted);
    }
}

#[test]
fn temperature_refresh_requeries_complete_sources_but_tau_changes_do_no_storage_work() {
    let card = card(1.0, [2.0, 0.5]);
    let initial = circuit(&card, 2.5, -0.25);
    let (before, _) = initial
        .rebind_reversible_sources_and_resistors_candidate(
            source(&card, 300.0),
            instantaneous(&card, 300.0, 8.0),
            resistors(&card, 300.0),
            &[],
            || false,
        )
        .unwrap()
        .into_parts();
    let old_law = before.reversible_instantaneous().unwrap();
    let query = old_law
        .material()
        .material()
        .query()
        .with("T", 350.0)
        .unwrap();
    let refreshed = old_law.at_temperature_query_point(&query).unwrap();
    assert_eq!(refreshed.curve(), old_law.curve());
    assert!(std::ptr::eq(
        refreshed.curve().knots().as_ptr(),
        old_law.curve().knots().as_ptr()
    ));
    assert_eq!(
        refreshed.material().material().state_coordinate("T"),
        Some(350.0)
    );
    let candidate = before
        .rebind_reversible_sources_and_resistors_candidate(
            before
                .source()
                .at_query_point(&card, &point(350.0))
                .unwrap(),
            refreshed,
            resistors(&card, 350.0),
            &[],
            || false,
        )
        .unwrap();
    let after = candidate.circuit();
    same_history(&before, after);
    assert_eq!(candidate.parameter_work_j(), 0.0);
    assert_eq!(
        before.dielectric_energy_j().unwrap(),
        after.dielectric_energy_j().unwrap()
    );
    for (old, new) in before.source().poles().iter().zip(after.source().poles()) {
        near(new.relaxation_time_s, 1.5 * old.relaxation_time_s);
        near(new.resistance_ohm, 1.5 * old.resistance_ohm);
    }
    for r in after.material_resistors() {
        near(r.conductor.resistance_ohm(), 10.0);
        assert_eq!(r.conductor.material().state_coordinate("T"), Some(350.0));
    }
    assert!(
        old_law
            .at_temperature_query_point(&query.clone().with("E", 0.1).unwrap())
            .is_err()
    );
    assert!(
        old_law
            .at_temperature_query_point(&query.with("T", 450.0).unwrap())
            .is_err()
    );
    assert!(
        before
            .rebind_source_candidate(source(&card, 350.0), &[], || false)
            .is_err()
    );
}

#[test]
fn incompatible_source_pair_geometry_support_clamp_and_cancellation_are_atomic() {
    let old_card = card(1.0, [2.0, 0.5]);
    let new_card = card(2.0, [1.0, 1.0]);
    let before = circuit(&old_card, 2.5, -0.25);
    let accepted = bits(&before);
    let old_source = before.source().clone();
    let old_law = before.reversible_instantaneous().unwrap().clone();
    let proposal = || {
        (
            source(&new_card, 300.0),
            instantaneous(&new_card, 300.0, 8.0),
        )
    };
    let (next_source, law) = proposal();
    let mut polls = 0;
    before
        .rebind_reversible_sources_candidate(next_source, law, &[], || {
            polls += 1;
            false
        })
        .unwrap();
    for cancel_at in [1, polls] {
        let (source, law) = proposal();
        let mut calls = 0;
        assert!(matches!(
            before.rebind_reversible_sources_candidate(source, law, &[], || {
                calls += 1;
                calls == cancel_at
            }),
            Err(MaterialCircuitError::Cancelled)
        ));
    }
    assert!(
        before
            .rebind_reversible_sources_candidate(
                source(&new_card, 300.0),
                old_law.clone(),
                &[],
                || false
            )
            .is_err()
    );
    let moved = MaterialDielectricCapacitor::new(0, old_law.material().clone(), 2.0, 1.0).unwrap();
    assert!(
        before
            .rebind_reversible_sources_candidate(old_source.clone(), moved, &[], || false)
            .is_err()
    );
    assert!(
        before
            .rebind_reversible_sources_candidate(
                source(&new_card, 300.0),
                instantaneous(&new_card, 300.0, 1.0),
                &[],
                || false
            )
            .is_err()
    );
    let reversed = names().into_iter().rev().collect::<Vec<_>>();
    let swapped_poles = ResolvedDebyeCapacitor::try_from_material(
        old_source.material().clone(),
        "epsilon_inf",
        &reversed,
        1.0,
        1.0,
    )
    .unwrap();
    assert!(
        before
            .rebind_reversible_sources_candidate(swapped_poles, old_law.clone(), &[], || false)
            .is_err()
    );
    assert!(
        before
            .rebind_reversible_sources_and_resistors_candidate(
                old_source.clone(),
                old_law.clone(),
                vec![resistors(&old_card, 300.0).remove(0)],
                &[],
                || false
            )
            .is_err()
    );
    let clamped = DebyeCircuit::new_with_reversible_instantaneous(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        old_source,
        old_law,
        vec![],
        vec![2.5],
        vec![0.6, -0.2],
        &[5.0],
        || false,
    )
    .unwrap();
    let clamp_state = bits(&clamped);
    let (source, law) = proposal();
    assert!(
        clamped
            .rebind_reversible_sources_candidate(source, law, &[5.0], || false)
            .is_err()
    );
    assert_eq!(bits(&clamped), clamp_state);
    assert_eq!(bits(&before), accepted);
}
