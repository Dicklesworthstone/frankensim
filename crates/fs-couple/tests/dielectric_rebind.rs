//! Temperature refresh preserves native storage history and complete sources.
//! Authored numerical fixtures, not measured material properties.
use fs_couple::material_circuit::{
    MaterialCircuitError, MaterialResistor, dielectric::DebyeCircuit,
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
        dielectric::{DebyePoleProperties, ResolvedDebyeCapacitor},
    },
    state_point::MaterialPropertySelection,
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_qty::{
    QuantitySpec, Temperature, Time,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
        .with_quantity("p", quantity(QuantityKind::Pressure), 1e5)
        .unwrap()
}
fn fixture(variant: u32) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut add = |name: &str, dims, knots| {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        quantity(QuantityKind::AbsoluteTemperature),
                        250.0,
                        500.0,
                    )
                    .with_quantity("p", quantity(QuantityKind::Pressure), 0.5e5, 2e5),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored dielectric rebind fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    add(
        "epsilon_inf",
        ELECTRIC_PERMITTIVITY_DIMS,
        vec![(250.0, 1.0), (500.0, 6.0)],
    );
    add(
        "delta",
        ELECTRIC_PERMITTIVITY_DIMS,
        vec![(250.0, 0.5), (500.0, 3.0)],
    );
    add("tau", Time::DIMS, vec![(250.0, 0.2), (500.0, 0.7)]);
    add(
        "electrical_resistivity",
        ELECTRICAL_RESISTIVITY_DIMS,
        vec![(250.0, 1.0), (500.0, 6.0)],
    );
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic dielectric".into(),
            phase: "solid".into(),
            process: format!("numerical fixture {variant}"),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(card: &MaterialCard, t: f64) -> ResolvedDebyeCapacitor {
    ResolvedDebyeCapacitor::resolve(
        card,
        &point(t),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &[DebyePoleProperties {
            permittivity_increment: "delta".into(),
            relaxation_time: "tau".into(),
        }],
        1.0,
        1.0,
    )
    .unwrap()
}
fn resistor(card: &MaterialCard, context: &QueryPoint, length: f64) -> MaterialResistor {
    MaterialResistor {
        branch: 2,
        conductor: resolve_uniform_conductor(
            card,
            context,
            MaterialPropertySelection::SingleClaimOnly,
            length,
            0.1,
        )
        .unwrap(),
    }
}
fn circuit(card: &MaterialCard) -> DebyeCircuit {
    DebyeCircuit::new(
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (2, 0, Branch::Inductor { henries: 4.0 }),
                (2, 0, Branch::Resistor { ohms: 5.0 }),
                (2, 0, Branch::Capacitor { farads: 0.5 }),
            ],
            transformers: vec![],
        },
        0,
        source(card, 300.0),
        vec![0.2],
        vec![2.0, -0.25],
        vec![0.25],
        &[],
        || false,
    )
    .unwrap()
}
fn bits(circuit: &DebyeCircuit) -> Vec<u64> {
    circuit.state().iter().map(|v| v.to_bits()).collect()
}
fn same_history(before: &DebyeCircuit, after: &DebyeCircuit) {
    assert_eq!(before.dae().flux_index, after.dae().flux_index);
    assert_eq!(before.dae().charge_index, after.dae().charge_index);
    for &index in before
        .dae()
        .flux_index
        .iter()
        .chain(&before.dae().charge_index)
    {
        assert_eq!(
            before.state()[index].to_bits(),
            after.state()[index].to_bits()
        );
    }
}
fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-12, "{actual} != {expected}");
}

#[test]
fn refresh_preserves_every_native_charge_and_flux_with_signed_parameter_work() {
    let card = fixture(0);
    let before = circuit(&card);
    let accepted = bits(&before);
    let candidate = before
        .rebind_source_candidate(source(&card, 350.0), &[], || false)
        .unwrap();
    let after = candidate.circuit();
    assert_eq!(after.original_branch_count(), 4);
    assert_eq!(after.dielectric_branch(), 0);
    same_history(&before, after);
    assert_eq!(after.instantaneous_charge_c(), 2.0);
    assert_eq!(after.polarization_charges_c(), vec![0.25]);
    near(after.terminal_voltage_v(), 2.0 / 3.0);
    // Independent fixed-charge oracle: C_inf 2 -> 3 F, C_p 1 -> 1.5 F.
    // The other capacitor and the inductor retain their energy exactly.
    let work = 0.5 * 2.0_f64.powi(2) * (1.0 / 3.0 - 1.0 / 2.0)
        + 0.5 * 0.25_f64.powi(2) * (1.0 / 1.5 - 1.0);
    near(candidate.parameter_work_j(), work);
    near(
        after.dae().system.hamiltonian(after.state())
            - before.dae().system.hamiltonian(before.state()),
        work,
    );
    assert!(work < 0.0);
    let reverse = after
        .rebind_source_candidate(before.source().clone(), &[], || false)
        .unwrap();
    near(reverse.parameter_work_j(), -work);
    same_history(&before, reverse.circuit());
    assert_eq!(bits(&before), accepted);
    assert_eq!(before.source(), &source(&card, 300.0));
}

#[test]
fn incompatible_ideal_voltage_clamp_refuses_without_projecting_charge() {
    let card = fixture(0);
    let before = DebyeCircuit::new(
        CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Capacitor { farads: f64::NAN }),
                (1, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        0,
        source(&card, 300.0),
        vec![],
        vec![2.0],
        vec![0.25],
        &[1.0],
        || false,
    )
    .unwrap();
    let accepted = bits(&before);
    // Conserved 2 C on the refreshed 3 F capacitance cannot remain at 1 V.
    assert!(
        before
            .rebind_source_candidate(source(&card, 350.0), &[1.0], || false)
            .is_err()
    );
    assert_eq!(bits(&before), accepted);
    assert_eq!(before.source(), &source(&card, 300.0));
    assert_eq!(before.polarization_charges_c(), vec![0.25]);
    near(before.terminal_voltage_v(), 1.0);
}

#[test]
fn complete_conductor_sources_refresh_only_temperature_and_cancel_atomically() {
    let card = fixture(0);
    let original = circuit(&card);
    let initial_source = source(&card, 300.0);
    let initial_resistor = resistor(&card, &point(300.0), 0.5);
    for bindings in [
        vec![],
        vec![initial_resistor.clone(), initial_resistor.clone()],
    ] {
        assert!(
            original
                .rebind_source_and_resistors_candidate(
                    initial_source.clone(),
                    bindings,
                    &[],
                    || false
                )
                .is_err()
        );
    }
    for branch in [1, original.original_branch_count()] {
        let mut invalid = initial_resistor.clone();
        invalid.branch = branch;
        assert!(
            original
                .rebind_source_and_resistors_candidate(
                    initial_source.clone(),
                    vec![invalid],
                    &[],
                    || false
                )
                .is_err()
        );
    }
    let (bound, work) = original
        .rebind_source_and_resistors_candidate(
            initial_source,
            vec![initial_resistor.clone()],
            &[],
            || false,
        )
        .unwrap()
        .into_parts();
    assert_eq!(work, 0.0);
    assert_eq!(bound.material_resistors(), &[initial_resistor]);
    assert_eq!(bound.graph().branches[2].2, Branch::Resistor { ohms: 10.0 });
    same_history(&original, &bound);
    let accepted = bits(&bound);
    let accepted_records = bound.material_resistors().to_vec();
    let refreshed = source(&card, 350.0);
    let conductor = resistor(&card, &point(350.0), 0.5);
    let mut polls = 0;
    let candidate = bound
        .rebind_source_and_resistors_candidate(
            refreshed.clone(),
            vec![conductor.clone()],
            &[],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls >= 3);
    assert_eq!(
        candidate.circuit().material_resistors(),
        &[conductor.clone()]
    );
    assert_eq!(
        candidate.circuit().graph().branches[2].2,
        Branch::Resistor { ohms: 15.0 }
    );
    same_history(&bound, candidate.circuit());
    near(candidate.parameter_work_j(), -0.34375);
    let retained = bound
        .rebind_source_candidate(refreshed.clone(), &[], || false)
        .unwrap();
    assert_eq!(retained.circuit().material_resistors(), accepted_records);
    assert_eq!(
        retained.circuit().graph().branches[2].2,
        Branch::Resistor { ohms: 10.0 }
    );
    let substituted_card = fixture(1);
    for bindings in [
        vec![],
        vec![resistor(&card, &point(350.0), 1.0)],
        vec![resistor(
            &card,
            &point(350.0).with("p", 1.2e5).unwrap(),
            0.5,
        )],
        vec![resistor(&substituted_card, &point(350.0), 0.5)],
    ] {
        assert!(
            bound
                .rebind_source_and_resistors_candidate(refreshed.clone(), bindings, &[], || false)
                .is_err()
        );
    }
    // Exercise cancellation at entry, during preparation, and immediately before publication.
    for stop in [1, 3, polls] {
        let mut observed = 0;
        assert!(matches!(
            bound.rebind_source_and_resistors_candidate(
                refreshed.clone(),
                vec![conductor.clone()],
                &[],
                || {
                    observed += 1;
                    observed >= stop
                }
            ),
            Err(MaterialCircuitError::Cancelled)
        ));
        assert_eq!(observed, stop);
    }
    assert_eq!(bits(&bound), accepted);
    assert_eq!(bound.material_resistors(), accepted_records);
    assert_eq!(bound.source(), &source(&card, 300.0));
}
