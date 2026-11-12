//! Endpoint viscosity refreshes original claims while reference inventory stays fixed.
use super::*;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::sheet::MembraneLaw;
use fs_material::state_point::resolve_mean_thermal_expansion_strain;
use fs_qty::{QuantitySpec, Temperature};

const VALUES: [f64; 6] = [4.0, 3.0, 2.0, 1.0, 0.4, 0.2];

struct Fixture {
    card: MaterialCard,
    reference: ResolvedSheet,
    thermal: ResolvedSheet,
    viscosity: PlaneStressViscosity,
}

fn fixture(
    original_max_temperature: f64,
    alternative: bool,
    active_endpoint: bool,
    endpoint_temperature: f64,
) -> Fixture {
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let provenance = Provenance {
        source: "synthetic thermal viscosity endpoint; no specimen authority".into(),
        license: "test".into(),
        artifact: None,
    };
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [
        (AREAL_DENSITY, 6.0, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, 0.002, Dims([1, 0, 0, 0, 0, 0])),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                // This reference inventory intentionally has NO hot-endpoint support.
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature,
                    299.0,
                    301.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    for (i, (name, value)) in KEYS.into_iter().zip(VALUES).enumerate() {
        let at = |t: f64| {
            if active_endpoint && i == 3 {
                value + (t - 300.0) * 99.0 / 30.0
            } else {
                value * (1.0 + 2.0 * (t - 300.0) / 30.0)
            }
        };
        let mut knots = vec![(300.0, value)];
        if original_max_temperature > 330.0 {
            knots.push((330.0, at(330.0)));
        }
        knots.push((original_max_temperature, at(original_max_temperature)));
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, DynViscosity::DIMS),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots,
                    dims: DynViscosity::DIMS,
                },
                interpolation: InterpolationPolicy::LinearInside,
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature, 300.0, original_max_temperature)
                    .with("omega", 0.0, 1000.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
        if alternative {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, DynViscosity::DIMS),
                    value: PropertyValue::Scalar {
                        value: 99.0 * value,
                        dims: DynViscosity::DIMS,
                    },
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    validity: ValidityDomain::unconstrained()
                        .with_quantity("T", temperature, 325.0, 340.0)
                        .with("omega", 0.0, 1000.0),
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: provenance.clone(),
                })
                .unwrap();
        }
        requirements.push(
            ScalarPropertyRequirement::try_new(
                name,
                DynViscosity::DIMS,
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let expansion_names = [
        "mean_warp_expansion",
        "mean_weft_expansion",
        "mean_normal_expansion",
    ];
    for (name, alpha) in expansion_names.into_iter().zip([0.001, 0.002, 0.0005]) {
        let dims = Dims([0, 0, 0, -1, 0, 0]);
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value: alpha, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T_ref", temperature, 299.0, 301.0)
                    .with_quantity("T_hot", temperature, 300.0, 340.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let model = ConstitutiveModelCard {
        law: LawId("neo-hookean-log".into()),
        law_version: 1,
        parameters: [
            ("mu", 3.0e5, pa),
            ("lambda", 8.0e5, pa),
            ("stretch_limit", 5.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with_quantity("T", temperature, 290.0, 340.0),
        sources: vec![],
        provenance,
    };
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic hyperelastic sheet".into(),
            phase: "solid".into(),
            process: "pinned thermal loss verification".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with_quantity("T", temperature, 300.0)
            .unwrap()
            .with("omega", 0.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let viscosity =
        PlaneStressViscosity::from_material(material.clone(), MaterialPlane::Axes12, KEYS).unwrap();
    let reference = ResolvedSheet::try_from_material(material)
        .unwrap()
        .with_hyperelastic_membrane(
            &card,
            0,
            PlaneStressSettings {
                relative_tolerance: 2e-13,
                ..PlaneStressSettings::default()
            },
        )
        .unwrap();

    // Carry the exact reference query's other coordinates into the existing
    // thermal reference owner. Frequency remains only a held loss-applicability
    // coordinate; none of the expansion sources depends on it.
    let reference_point = reference.material().query();
    let mut expansion_point = QueryPoint::new()
        .with_quantity("T_ref", temperature, 300.0)
        .unwrap()
        .with_quantity("T_hot", temperature, endpoint_temperature)
        .unwrap();
    for (axis, value) in reference_point.axes() {
        if axis == "T" {
            continue;
        }
        expansion_point = match reference_point.axis_quantities().get(axis) {
            Some(quantity) => expansion_point.with_quantity(axis, *quantity, *value),
            None => expansion_point.with(axis, *value),
        }
        .unwrap();
    }
    let expansion = expansion_names.map(|name| {
        resolve_mean_thermal_expansion_strain(
            &card,
            &PropertyKey::new(name, Dims([0, 0, 0, -1, 0, 0])),
            &expansion_point,
            "T_ref",
            "T_hot",
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    });
    let thermal = reference
        .clone()
        .with_hyperelastic_mean_thermal_expansion(&card, expansion, "T_ref", "T_hot")
        .unwrap();
    Fixture {
        card,
        reference,
        thermal,
        viscosity,
    }
}

fn patch_from_sheet(sheet: ResolvedSheet) -> MembranePatch {
    MembranePatch::try_new(
        Mesh2 {
            nodes: vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
            elems: vec![vec![0, 1, 2]],
            patches: vec![],
        },
        vec![sheet],
        vec![0.37],
        MembraneDomain {
            minimum_area_ratio: 0.01,
            maximum_abs_strain: 0.5,
        },
        MembraneBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        || false,
    )
    .unwrap()
}

fn thermal_state() -> Vec<[f64; 3]> {
    vec![[0.0; 3], [0.20, 0.08, 0.19], [-0.03, 0.15, 0.25]]
}

fn endpoint(sheet: &ResolvedSheet) -> QueryPoint {
    let MembraneLaw::Hyperelastic(elastic) = sheet.membrane_law().unwrap() else {
        panic!("fixture must retain its hyperelastic owner");
    };
    elastic.model().context().clone()
}

fn stretches(sheet: &ResolvedSheet) -> [f64; 3] {
    let MembraneLaw::Hyperelastic(elastic) = sheet.membrane_law().unwrap() else {
        panic!("fixture must retain its hyperelastic owner");
    };
    elastic.reference_stretches()
}

#[test]
fn endpoint_tensor_changes_force_heat_and_full_modal_loss_without_refreshing_inventory() {
    // Another hot claim is available, but the original curve remains pinned.
    let fixture = fixture(340.0, true, false, 330.0);
    let p = patch_from_sheet(fixture.thermal.clone());
    let expected_context = endpoint(&fixture.thermal);
    let original_identity = fixture.viscosity.source().identity();
    assert!(
        fixture
            .viscosity
            .source()
            .at_query_point(&fixture.card, &expected_context)
            .is_err()
    );
    assert!(
        MembraneViscosity::new(&p, vec![fixture.viscosity.clone()], (0.0, 1000.0), || false)
            .is_err()
    );
    let actual = MembraneViscosity::new_with_thermal_endpoints(
        &p,
        vec![fixture.viscosity.clone()],
        &[&fixture.card],
        (0.0, 1000.0),
        || false,
    )
    .unwrap();
    assert!(core::ptr::eq(actual.patch(), &p));
    assert_eq!(
        actual.reference_laws()[0].source().identity(),
        original_identity
    );
    assert_eq!(
        actual.reference_laws()[0].source(),
        p.sheets()[0].material()
    );
    assert_eq!(actual.laws()[0].source().query(), expected_context);
    assert_eq!(
        actual.laws()[0].source().axis_quantities(),
        fixture.viscosity.source().axis_quantities()
    );
    assert_eq!(
        actual.laws()[0].source().card_identity(),
        fixture.card.content_hash()
    );
    assert_eq!(actual.laws()[0].source().properties().len(), 6);
    assert!(actual.laws()[0].source().property(AREAL_DENSITY).is_none());
    assert_eq!(p.sheets()[0].material().state_coordinate("T"), Some(300.0));
    assert_eq!(actual.laws()[0].source().state_coordinate("T"), Some(330.0));
    assert_eq!(
        actual.laws()[0].source().state_coordinate("omega"),
        Some(0.0)
    );
    for name in KEYS {
        let reference = fixture.viscosity.source().property(name).unwrap();
        let current = actual.laws()[0].source().property(name).unwrap();
        assert_eq!(
            current.answer().receipt.selected,
            reference.answer().receipt.selected
        );
        assert_eq!(current.requirement(), reference.requirement());
        assert_eq!(current.source_axes(), reference.source_axes());
        close(current.value_si(), 3.0 * reference.value_si());
    }

    // Hold the SAME natural G and mechanics, but use the original-temperature
    // viscosity as a causal control. Only the source viscosity changes.
    let cold = patch_from_sheet(
        fixture
            .reference
            .clone()
            .with_hyperelastic_reference_stretches(stretches(&fixture.thermal))
            .unwrap(),
    );
    let control = MembraneViscosity::new(
        &cold,
        vec![fixture.viscosity.clone()],
        (0.0, 1000.0),
        || false,
    )
    .unwrap();
    assert_eq!(p.reference_mass_kg(), cold.reference_mass_kg());
    assert_eq!(p.sheets()[0].thickness_m(), cold.sheets()[0].thickness_m());
    let u = thermal_state();
    close(
        p.evaluate(&u, || false).unwrap().strain_energy_j,
        cold.evaluate(&u, || false).unwrap().strain_energy_j,
    );
    let velocity = [[0.0; 3], [0.0, 0.0, 1.3], [0.0, 0.0, -0.7]];
    let hot = actual.evaluate(&u, &velocity, || false).unwrap();
    let old = control.evaluate(&u, &velocity, || false).unwrap();
    close(hot.power_w, 3.0 * old.power_w);
    close(hot.element_power_w[0], 3.0 * old.element_power_w[0]);
    for (current, original) in hot
        .resisting_force_n
        .iter()
        .flatten()
        .zip(old.resisting_force_n.iter().flatten())
    {
        close(*current, 3.0 * original);
    }
    let work: f64 = hot
        .resisting_force_n
        .iter()
        .flatten()
        .zip(velocity.iter().flatten())
        .map(|(f, v)| f * v)
        .sum();
    close(work, hot.power_w);
    let hot_modes = modes(&p, &u, &[5, 8]);
    let cold_modes = modes(&cold, &u, &[5, 8]);
    let hot_matrix = actual.project_modes(&hot_modes, 100_000, || false).unwrap();
    let cold_matrix = control
        .project_modes(&cold_modes, 100_000, || false)
        .unwrap();
    assert!(hot_matrix.matrix_per_s()[1].abs() > 1e-8);
    for (current, original) in hot_matrix
        .matrix_per_s()
        .iter()
        .zip(cold_matrix.matrix_per_s())
    {
        close(*current, 3.0 * original);
    }
    let rate = [0.7, -0.3];
    let modal_velocity = hot_modes.reconstruct(&rate).unwrap();
    let response = actual.evaluate(&u, &modal_velocity, || false).unwrap();
    let gram_power: f64 = hot_matrix
        .dissipation_rows()
        .chunks_exact(2)
        .map(|row| (row[0] * rate[0] + row[1] * rate[1]).powi(2))
        .sum();
    close(gram_power, response.power_w);
    assert_eq!(fixture.viscosity.source().identity(), original_identity);
}

#[test]
fn expired_original_claims_cannot_switch_to_available_hot_claims() {
    let fixture = fixture(320.0, true, false, 330.0);
    let p = patch_from_sheet(fixture.thermal.clone());
    let requirements: Vec<_> = KEYS
        .map(|name| {
            fixture
                .viscosity
                .source()
                .property(name)
                .unwrap()
                .requirement()
                .clone()
        })
        .into_iter()
        .collect();
    // A loose endpoint query could succeed using different claims.
    let alternative = resolve_material_state_point(
        &fixture.card,
        &endpoint(&fixture.thermal),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    for name in KEYS {
        assert_ne!(
            alternative
                .property(name)
                .unwrap()
                .answer()
                .receipt
                .selected,
            fixture
                .viscosity
                .source()
                .property(name)
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    let result = MembraneViscosity::new_with_thermal_endpoints(
        &p,
        vec![fixture.viscosity.clone()],
        &[&fixture.card],
        (0.0, 1000.0),
        || false,
    );
    assert!(matches!(
        result,
        Err(MembraneViscosityError::Source { element: 0, .. })
    ));
    assert_eq!(fixture.viscosity.source(), p.sheets()[0].material());
    assert_eq!(
        fixture.viscosity.source().state_coordinate("T"),
        Some(300.0)
    );
}

#[test]
fn source_supported_but_nonpassive_endpoint_refuses_without_repair() {
    let fixture = fixture(340.0, false, true, 330.0);
    let p = patch_from_sheet(fixture.thermal);
    let result = MembraneViscosity::new_with_thermal_endpoints(
        &p,
        vec![fixture.viscosity.clone()],
        &[&fixture.card],
        (0.0, 1000.0),
        || false,
    );
    assert!(matches!(
        result,
        Err(MembraneViscosityError::Material { element: 0, .. })
    ));
    assert_eq!(
        fixture.viscosity.matrix_pa_s(),
        &[4.0, 1.0, 0.4, 1.0, 3.0, 0.2, 0.4, 0.2, 2.0]
    );
}

#[test]
fn zero_temperature_change_preserves_original_force_heat_and_modal_arithmetic() {
    let fixture = fixture(340.0, false, false, 300.0);
    let ordinary = patch_from_sheet(fixture.reference);
    let selected = patch_from_sheet(fixture.thermal);
    let a = MembraneViscosity::new(
        &ordinary,
        vec![fixture.viscosity.clone()],
        (0.0, 1000.0),
        || false,
    )
    .unwrap();
    let b = MembraneViscosity::new_with_thermal_endpoints(
        &selected,
        vec![fixture.viscosity],
        &[&fixture.card],
        (0.0, 1000.0),
        || false,
    )
    .unwrap();
    let u = thermal_state();
    let velocity = [[0.0; 3], [0.0, 0.0, 1.3], [0.0, 0.0, -0.7]];
    let left = a.evaluate(&u, &velocity, || false).unwrap();
    let right = b.evaluate(&u, &velocity, || false).unwrap();
    assert_eq!(left.resisting_force_n, right.resisting_force_n);
    assert_eq!(left.power_w, right.power_w);
    assert_eq!(left.element_power_w, right.element_power_w);
    let left_modes = modes(&ordinary, &u, &[5, 8]);
    let right_modes = modes(&selected, &u, &[5, 8]);
    let left = a.project_modes(&left_modes, 100_000, || false).unwrap();
    let right = b.project_modes(&right_modes, 100_000, || false).unwrap();
    assert_eq!(left.matrix_per_s(), right.matrix_per_s());
    assert_eq!(left.dissipation_rows(), right.dissipation_rows());
}

#[test]
fn endpoint_binding_keeps_original_card_axes_context_band_and_atomic_cancellation() {
    let fixture = fixture(340.0, false, false, 330.0);
    let p = patch_from_sheet(fixture.thermal);
    let foreign = MaterialCard::supersede(
        &fixture.card,
        fixture.card.claims().clone(),
        fixture.card.models().to_vec(),
    )
    .unwrap();
    assert!(
        MembraneViscosity::new_with_thermal_endpoints(
            &p,
            vec![fixture.viscosity.clone()],
            &[&foreign],
            (0.0, 1000.0),
            || false,
        )
        .is_err()
    );
    assert!(
        MembraneViscosity::new_with_thermal_endpoints(
            &p,
            vec![fixture.viscosity.clone()],
            &[],
            (0.0, 1000.0),
            || false,
        )
        .is_err()
    );
    let axes = PlaneStressViscosity::from_material(
        fixture.viscosity.source().clone(),
        MaterialPlane::Axes21,
        KEYS,
    )
    .unwrap();
    assert!(
        MembraneViscosity::new_with_thermal_endpoints(
            &p,
            vec![axes],
            &[&fixture.card],
            (0.0, 1000.0),
            || false,
        )
        .is_err()
    );
    let moved = fixture
        .viscosity
        .source()
        .at_query_point(
            &fixture.card,
            &fixture
                .viscosity
                .source()
                .query()
                .with("omega", 10.0)
                .unwrap(),
        )
        .unwrap();
    let wrong_context =
        PlaneStressViscosity::from_material(moved, MaterialPlane::Axes12, KEYS).unwrap();
    assert!(
        MembraneViscosity::new_with_thermal_endpoints(
            &p,
            vec![wrong_context],
            &[&fixture.card],
            (0.0, 1000.0),
            || false,
        )
        .is_err()
    );
    assert!(
        MembraneViscosity::new_with_thermal_endpoints(
            &p,
            vec![fixture.viscosity.clone()],
            &[&fixture.card],
            (0.0, 1001.0),
            || false,
        )
        .is_err()
    );
    let identity = fixture.viscosity.source().identity();
    let mut polls = 0;
    let cancelled = MembraneViscosity::new_with_thermal_endpoints(
        &p,
        vec![fixture.viscosity.clone()],
        &[&fixture.card],
        (0.0, 1000.0),
        || {
            polls += 1;
            polls == 6
        },
    );
    assert!(matches!(
        cancelled,
        Err(MembraneViscosityError::Membrane(MembraneError::Cancelled))
    ));
    assert_eq!(fixture.viscosity.source().identity(), identity);
    assert_eq!(fixture.viscosity.source(), p.sheets()[0].material());
    MembraneViscosity::new_with_thermal_endpoints(
        &p,
        vec![fixture.viscosity],
        &[&fixture.card],
        (0.0, 1000.0),
        || false,
    )
    .unwrap();
}
