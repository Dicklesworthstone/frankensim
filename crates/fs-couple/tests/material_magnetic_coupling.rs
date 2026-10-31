//! Synthetic source/energy checks for material-derived reciprocal winding links.
use fs_couple::material_circuit::{
    MaterialCircuit, MaterialCircuitError, MaterialCircuitSolve,
    conductivity::ConductivitySolve,
    storage::{MaterialMagneticCoupling, MaterialStorage, MaterialStorageElement},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY, ResolvedCapacitor, ResolvedInductor,
    },
    state_point::{MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement},
};
use fs_phs::circuit::{Branch, CircuitGraph, MutualInductance};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

fn claim(name: &str, dims: Dims, value: PropertyValue) -> PropertyClaim {
    let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else {
        InterpolationPolicy::LinearInside
    };
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value,
        validity: ValidityDomain::unconstrained(),
        interpolation,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored shared-permeability fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card(claims: ClaimSet, revision: u32) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic shared field".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn magnetic_claim() -> PropertyClaim {
    claim(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 1.0), (500.0, 4.0)],
            dims: MAGNETIC_PERMEABILITY_DIMS,
        },
    )
}
fn magnetic_card(revision: u32) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(magnetic_claim()).unwrap();
    let density = Dims([-3, 1, 0, 0, 0, 0]);
    claims
        .insert_claim(claim(
            "density",
            density,
            PropertyValue::Scalar {
                value: 1000.0,
                dims: density,
            },
        ))
        .unwrap();
    card(claims, revision)
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
}
fn inductor(source: &MaterialCard, at: &QueryPoint, turns: u32) -> ResolvedInductor {
    ResolvedInductor::resolve(
        source,
        at,
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
        1.0,
        turns,
    )
    .unwrap()
}
fn storage(source: &MaterialCard, temperature: f64, turns: &[u32]) -> Vec<MaterialStorage> {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(
            ELECTRIC_PERMITTIVITY_PROPERTY,
            ELECTRIC_PERMITTIVITY_DIMS,
            PropertyValue::Scalar {
                value: 1.0,
                dims: ELECTRIC_PERMITTIVITY_DIMS,
            },
        ))
        .unwrap();
    let dielectric = card(claims, 0);
    let capacitor = ResolvedCapacitor::resolve(
        &dielectric,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        1.0,
        1.0,
    )
    .unwrap();
    let mut storage = Vec::new();
    for (index, &turns) in turns.iter().enumerate() {
        storage.push(MaterialStorage {
            branch: 2 * index,
            element: MaterialStorageElement::Inductor(inductor(source, &point(temperature), turns)),
        });
        storage.push(MaterialStorage {
            branch: 2 * index + 1,
            element: MaterialStorageElement::Capacitor(capacitor.clone()),
        });
    }
    storage
}
fn coupling(a: usize, b: usize, geometry_m: f64, constant_h: f64) -> MaterialMagneticCoupling {
    MaterialMagneticCoupling {
        branches: (a, b),
        geometry_m,
        constant_h,
    }
}
fn circuit(
    storage: Vec<MaterialStorage>,
    held: Vec<MutualInductance>,
    couplings: Vec<MaterialMagneticCoupling>,
    fluxes: Vec<f64>,
) -> Result<MaterialCircuit, MaterialCircuitError> {
    let windings = fluxes.len();
    let mut branches = Vec::new();
    for node in 1..=windings {
        branches.push((node, 0, Branch::Inductor { henries: f64::NAN }));
        branches.push((node, 0, Branch::Capacitor { farads: f64::NAN }));
    }
    MaterialCircuit::new_with_magnetic_coupling(
        CircuitGraph {
            node_count: windings + 1,
            branches,
            transformers: vec![],
        },
        vec![],
        storage,
        held,
        couplings,
        fluxes,
        vec![0.0; windings],
        &[],
    )
}
fn controls() -> MaterialCircuitSolve {
    MaterialCircuitSolve {
        energy_tolerance_j: 1e-11,
        power_tolerance_w: 1e-11,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 5e-12 * b.abs().max(1.0), "{a:e} != {b:e}");
}
fn partition(mu: f64, m: f64, flux: [f64; 2]) -> [f64; 2] {
    let determinant = 4.0 * mu * mu - m * m;
    let current = [
        (4.0 * mu * flux[0] - m * flux[1]) / determinant,
        (mu * flux[1] - m * flux[0]) / determinant,
    ];
    [0.5 * flux[0] * current[0], 0.5 * flux[1] * current[1]]
}

#[test]
fn shared_source_refresh_changes_full_matrix_at_retained_total_linkage() {
    let source = magnetic_card(0);
    let flux = [1.2, -0.4];
    let original = circuit(
        storage(&source, 300.0, &[1, 2]),
        vec![],
        vec![coupling(2, 0, -0.5, 0.125)],
        flux.to_vec(),
    )
    .unwrap();
    assert_eq!(
        original.magnetic_couplings(),
        &[coupling(0, 2, -0.5, 0.125)]
    );
    near(original.mutual_inductance()[0].henries, -0.875);
    let accepted = original.clone();
    let candidate = original
        .rebind_storage_candidate(storage(&source, 400.0, &[1, 2]), &[], || false)
        .unwrap();
    assert_eq!(original, accepted);
    assert_eq!(candidate.circuit().fluxes_wb(), flux);
    assert_eq!(
        candidate.circuit().magnetic_couplings(),
        original.magnetic_couplings()
    );
    near(candidate.circuit().mutual_inductance()[0].henries, -1.375);
    let before = partition(2.0, -0.875, flux);
    let after = partition(3.0, -1.375, flux);
    for (branch, slot) in [(0, 0), (2, 1)] {
        let work = candidate
            .work()
            .iter()
            .find(|work| work.branch == branch)
            .unwrap();
        near(work.stored_before_j, before[slot]);
        near(work.stored_after_j, after[slot]);
        near(work.supplied_work_j, after[slot] - before[slot]);
    }
    near(
        candidate.parameter_work_j(),
        after.iter().sum::<f64>() - before.iter().sum::<f64>(),
    );
    let reverse = candidate
        .circuit()
        .rebind_storage_candidate(original.storage().to_vec(), &[], || false)
        .unwrap();
    near(reverse.parameter_work_j(), -candidate.parameter_work_j());
    assert_eq!(
        reverse.circuit().mutual_inductance(),
        original.mutual_inductance()
    );

    // An explicit, paid material reconfiguration may replace both common sources.
    let replacement = magnetic_card(1);
    let replaced = original
        .rebind_storage_candidate(storage(&replacement, 400.0, &[1, 2]), &[], || false)
        .unwrap();
    near(replaced.parameter_work_j(), candidate.parameter_work_j());
    let mut mismatched = storage(&source, 400.0, &[1, 2]);
    mismatched[2].element = MaterialStorageElement::Inductor(inductor(&source, &point(300.0), 2));
    assert!(
        original
            .rebind_storage_candidate(mismatched, &[], || false)
            .is_err()
    );
    assert!(
        original
            .rebind_storage_candidate(original.storage().to_vec(), &[], || true)
            .is_err()
    );
    assert_eq!(original, accepted);
}

#[test]
fn held_links_and_source_metadata_survive_steps_contacts_and_conductivity() {
    let source = magnetic_card(0);
    let held = MutualInductance {
        branches: (4, 0),
        henries: 0.2,
    };
    let original = circuit(
        storage(&source, 300.0, &[1, 2, 3]),
        vec![held],
        vec![coupling(0, 2, 0.5, 0.0)],
        vec![0.4, 0.2, -0.3],
    )
    .unwrap();
    let rebound = original
        .rebind_storage_candidate(storage(&source, 400.0, &[1, 2, 3]), &[], || false)
        .unwrap();
    assert_eq!(rebound.circuit().mutual_inductance()[0], held);
    near(rebound.circuit().mutual_inductance()[1].henries, 1.5);
    let contacts = rebound
        .circuit()
        .rebind_contacts_candidate(vec![], &[])
        .unwrap();
    let step = contacts
        .step_candidate(vec![], &[], 0.01, controls(), || false)
        .unwrap();
    assert_eq!(
        step.circuit().magnetic_couplings(),
        original.magnetic_couplings()
    );
    assert_eq!(
        step.circuit().mutual_inductance(),
        contacts.mutual_inductance()
    );
    near(step.step().energy_residual_j, 0.0);
    let conductivity = step
        .circuit()
        .conductivity_candidate(
            &[],
            &[],
            0.01,
            controls(),
            ConductivitySolve {
                maximum_iterations: 1,
                absolute_resistance_tolerance_ohm: 1e-12,
                relative_resistance_tolerance: 0.0,
            },
            || false,
        )
        .unwrap();
    assert_eq!(
        conductivity.circuit().magnetic_couplings(),
        original.magnetic_couplings()
    );
    assert_eq!(
        conductivity.circuit().mutual_inductance(),
        contacts.mutual_inductance()
    );
    let final_refresh = conductivity
        .circuit()
        .rebind_storage_candidate(storage(&source, 350.0, &[1, 2, 3]), &[], || false)
        .unwrap();
    near(final_refresh.circuit().mutual_inductance()[1].henries, 1.25);
}

#[test]
fn same_actual_claim_allows_pinned_policy_and_unrelated_bundle_properties() {
    let source = magnetic_card(0);
    let mut owners = storage(&source, 300.0, &[1, 2]);
    let base = inductor(&source, &point(300.0), 2);
    let material = base
        .material()
        .at_query_point(&source, &point(300.0))
        .unwrap()
        .with_additional_requirements(
            &source,
            &[ScalarPropertyRequirement::try_new(
                "density",
                Dims([-3, 1, 0, 0, 0, 0]),
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    owners[2].element = MaterialStorageElement::Inductor(
        ResolvedInductor::try_from_material(material, 1.0, 1.0, 2).unwrap(),
    );
    assert!(
        circuit(
            owners,
            vec![],
            vec![coupling(0, 2, 0.5, 0.0)],
            vec![0.1, 0.2]
        )
        .is_ok()
    );
}

#[test]
fn equal_numeric_coefficients_cannot_replace_shared_source_identity() {
    let source = magnetic_card(0);
    let mut owners = storage(&source, 300.0, &[1, 1]);
    owners[2].element =
        MaterialStorageElement::Inductor(inductor(&magnetic_card(1), &point(300.0), 1));
    assert!(
        circuit(
            owners,
            vec![],
            vec![coupling(0, 2, 0.5, 0.0)],
            vec![0.0, 0.0]
        )
        .is_err()
    );

    let mut claims = ClaimSet::new();
    let first = claims.insert_claim(magnetic_claim()).unwrap();
    let mut alternative = magnetic_claim();
    alternative.provenance.source = "independent equal-valued claim".into();
    let second = claims.insert_claim(alternative).unwrap();
    let ambiguous = card(claims, 0);
    let mut owners = storage(&source, 300.0, &[1, 1]);
    for (slot, selected) in [(0, first), (2, second)] {
        owners[slot].element = MaterialStorageElement::Inductor(
            ResolvedInductor::resolve(
                &ambiguous,
                &point(300.0),
                MaterialPropertySelection::PinnedByProperty(vec![(
                    MAGNETIC_PERMEABILITY_PROPERTY.into(),
                    selected,
                )]),
                1.0,
                1.0,
                1,
            )
            .unwrap(),
        );
    }
    assert!(
        circuit(
            owners,
            vec![],
            vec![coupling(0, 2, 0.5, 0.0)],
            vec![0.0, 0.0]
        )
        .is_err()
    );

    // Identical constant mu still cannot share two physical contexts or schemas.
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            PropertyValue::Scalar {
                value: 2.0,
                dims: MAGNETIC_PERMEABILITY_DIMS,
            },
        ))
        .unwrap();
    let constant = card(claims, 0);
    for query in [point(301.0), QueryPoint::new().with("T", 300.0).unwrap()] {
        let mut owners = storage(&constant, 300.0, &[1, 1]);
        owners[2].element = MaterialStorageElement::Inductor(inductor(&constant, &query, 1));
        assert!(
            circuit(
                owners,
                vec![],
                vec![coupling(0, 2, 0.5, 0.0)],
                vec![0.0, 0.0]
            )
            .is_err()
        );
    }
}

#[test]
fn full_matrix_admission_refuses_nonpassive_or_ambiguous_couplings() {
    let source = magnetic_card(0);
    // Every pair satisfies |M|<sqrt(L1*L2), but the three-way sum mode has L<0.
    let links = vec![
        coupling(0, 2, -0.75, 0.0),
        coupling(0, 4, -0.75, 0.0),
        coupling(2, 4, -0.75, 0.0),
    ];
    assert!(
        circuit(
            storage(&source, 300.0, &[1, 1, 1]),
            vec![],
            links,
            vec![0.0; 3]
        )
        .is_err()
    );
    for links in [
        vec![coupling(0, 2, 1.0, 0.0)],
        vec![coupling(0, 0, 0.5, 0.0)],
        vec![coupling(0, 1, 0.5, 0.0)],
        vec![coupling(0, 2, f64::NAN, 0.0)],
        vec![coupling(0, 2, 0.5, f64::INFINITY)],
        vec![coupling(0, 2, f64::MAX, 0.0)],
        vec![coupling(0, 2, 0.5, 0.0), coupling(2, 0, 0.5, 0.0)],
    ] {
        assert!(
            circuit(
                storage(&source, 300.0, &[1, 1]),
                vec![],
                links,
                vec![0.0; 2]
            )
            .is_err()
        );
    }
    assert!(
        circuit(
            storage(&source, 300.0, &[1, 1]),
            vec![MutualInductance {
                branches: (2, 0),
                henries: 0.0
            }],
            vec![coupling(0, 2, 0.5, 0.0)],
            vec![0.0; 2]
        )
        .is_err()
    );
}

#[test]
fn real_zero_crossing_is_retained_but_nonzero_underflow_is_refused() {
    let source = magnetic_card(0);
    let original = circuit(
        storage(&source, 300.0, &[1, 1]),
        vec![],
        vec![coupling(0, 2, 0.5, -1.0)],
        vec![0.1, 0.2],
    )
    .unwrap();
    assert_eq!(original.mutual_inductance()[0].henries, 0.0);
    assert_eq!(original.magnetic_couplings().len(), 1);
    let warmed = original
        .rebind_storage_candidate(storage(&source, 400.0, &[1, 1]), &[], || false)
        .unwrap();
    near(warmed.circuit().mutual_inductance()[0].henries, 0.5);

    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            PropertyValue::Scalar {
                value: 1e-100,
                dims: MAGNETIC_PERMEABILITY_DIMS,
            },
        ))
        .unwrap();
    let tiny = card(claims, 0);
    assert!(
        circuit(
            storage(&tiny, 300.0, &[1, 1]),
            vec![],
            vec![coupling(0, 2, f64::from_bits(1), 0.0)],
            vec![0.0; 2],
        )
        .is_err()
    );
}
