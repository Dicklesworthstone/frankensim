//! Synthetic source/stiffness/mass oracles, not measured material validation.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::laminate::thermal::source::{
    IntegratedLaminaExpansion, LaminaElasticSource, LaminaExpansionSource, LaminaSourceInput,
    LaminateThermalSource, LaminateThermalSourceError,
};
use fs_material::laminate::{LaminateError, LayerLaw, MAX_LAMINA_LAYERS, ReferenceSurface};
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, IsotropicElasticParameters,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as CTE, MaterialPropertySelection,
    ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES as CTE3, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use std::sync::Arc;

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn at(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
}
fn domain(lo: f64, hi: f64) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity("T", temperature_quantity(), lo, hi)
}
fn claim(name: &str, dims: Dims, knots: &[(f64, f64)]) -> PropertyClaim {
    let constant = knots.len() == 1;
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if constant {
            PropertyValue::Scalar {
                value: knots[0].1,
                dims,
            }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: knots.to_vec(),
                dims,
            }
        },
        interpolation: if constant {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        validity: domain(300.0, 400.0),
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored laminate thermal source oracle".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn iso_claims(reference_density_only: bool) -> Vec<PropertyClaim> {
    let mut density = claim(DENSITY_PROPERTY, Density::DIMS, &[(300.0, 7800.0)]);
    if reference_density_only {
        density.validity = domain(300.0, 300.0);
    }
    vec![
        density,
        claim(
            "young_modulus",
            Pressure::DIMS,
            &[(300.0, 120e9), (400.0, 60e9)],
        ),
        claim("poisson_ratio", Dims::NONE, &[(300.0, 0.25)]),
        claim(
            CTE,
            INVERSE_TEMPERATURE_DIMS,
            &[(300.0, 1e-5), (350.0, 4e-5), (400.0, 2e-5)],
        ),
    ]
}
fn card(claims: Vec<PropertyClaim>) -> Arc<MaterialCard> {
    let mut set = ClaimSet::new();
    for c in claims {
        set.insert_claim(c).unwrap();
    }
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic thermal laminate".into(),
                phase: "solid".into(),
                process: "numerical source test".into(),
                revision: 0,
            },
            set,
            vec![],
        )
        .unwrap(),
    )
}
fn req(name: &str, dims: Dims, admissibility: ScalarAdmissibility) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, admissibility).unwrap()
}
fn input(card: Arc<MaterialCard>, h: f64) -> LaminaSourceInput {
    LaminaSourceInput {
        card,
        reference: at(300.0),
        thickness_m: h,
        orientation_rad: 0.0,
        elasticity: LaminaElasticSource::Isotropic(IsotropicElasticParameters::YoungPoisson {
            young_modulus: req(
                "young_modulus",
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            ),
            poisson_ratio: req("poisson_ratio", Dims::NONE, ScalarAdmissibility::Finite),
        }),
        expansion: LaminaExpansionSource::Isotropic,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-20),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_current_layer_stiffness_and_integrated_cte_keep_reference_mass() {
    let original = card(iso_claims(true));
    let source = LaminateThermalSource::try_new(
        vec![
            input(Arc::clone(&original), 0.001),
            input(Arc::clone(&original), 0.003),
        ],
        ReferenceSurface::Offset(0.0017),
        || false,
    )
    .unwrap();
    let current = source.at_temperatures(&[400.0, 350.0], || false).unwrap();
    let prepared = current.prepared();
    // Independent reciprocal plane-stress Q and exact trapezoidal CTE integrals.
    let mut force = 0.0;
    let mut moment = 0.0;
    for (i, (e, h, z, free)) in [
        (60e9, 0.001, -0.0012, 0.00275),
        (90e9, 0.003, 0.0008, 0.00125),
    ]
    .into_iter()
    .enumerate()
    {
        let n = e / (1.0 - 0.25) * free * h;
        force += n;
        moment += n * z;
        near(
            prepared.section().layer_integrals()[i].q_bar_pa[0],
            e / (1.0 - 0.25 * 0.25),
            2e-15,
        );
        near(prepared.profiles()[i].faces()[0][0], free, 2e-15);
        assert!(Arc::ptr_eq(&source.layers()[i].input().card, &original));
        assert_eq!(
            current.layers()[i].elasticity().card_identity(),
            original.content_hash()
        );
        assert!(
            current.layers()[i]
                .elasticity()
                .property(DENSITY_PROPERTY)
                .is_none()
        );
        assert_eq!(source.layers()[i].reference_density().value_si(), 7800.0);
        assert_eq!(
            current.source().layers()[i].reference_material(),
            source.layers()[i].reference_material()
        );
    }
    near(prepared.membrane_force_n_m()[0], force, 2e-15);
    near(prepared.bending_moment_n()[0], moment, 2e-15);
    let reference = source.reference_section();
    assert_eq!(
        prepared.section().i0_kg_m2().to_bits(),
        reference.i0_kg_m2().to_bits()
    );
    assert_eq!(
        prepared.section().i1_kg_m().to_bits(),
        reference.i1_kg_m().to_bits()
    );
    assert_eq!(
        prepared.section().i2_kg().to_bits(),
        reference.i2_kg().to_bits()
    );
    assert_ne!(prepared.section().a_pa_m(), reference.a_pa_m());
    // The original density cannot even be queried at 400 K; the update still
    // succeeds because a closed reference inventory is not an endpoint density.
    assert!(
        resolve_material_state_point(
            &original,
            &at(400.0),
            &[req(
                DENSITY_PROPERTY,
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive
            ),],
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn g1_heating_and_cooling_use_each_reference_and_preserve_typed_context() {
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    let mut claims = iso_claims(false);
    for claim in &mut claims {
        claim.validity =
            claim
                .validity
                .clone()
                .with_quantity("ambient_pressure", pressure, 1e5, 3e5);
    }
    let original = card(claims);
    let mut heating = input(Arc::clone(&original), 0.002);
    heating.reference = at(300.0)
        .with_quantity("ambient_pressure", pressure, 2e5)
        .unwrap();
    let mut cooling = heating.clone();
    cooling.reference = cooling.reference.with("T", 400.0).unwrap();
    let source = LaminateThermalSource::try_new(
        vec![heating, cooling],
        ReferenceSurface::Midsurface,
        || false,
    )
    .unwrap();
    let state = source.at_temperatures(&[400.0, 300.0], || false).unwrap();
    for (i, sign) in [1.0, -1.0].into_iter().enumerate() {
        near(
            state.prepared().profiles()[i].faces()[0][0],
            sign * 0.00275,
            2e-15,
        );
        for resolved in [
            state.layers()[i].elasticity(),
            state.layers()[i].expansion().current(),
        ] {
            assert_eq!(resolved.state_coordinate("ambient_pressure"), Some(2e5));
            assert_eq!(
                resolved.axis_quantities().get("ambient_pressure"),
                Some(&pressure)
            );
            assert_eq!(
                resolved.axis_quantities().get("T"),
                Some(&temperature_quantity())
            );
        }
    }
}

#[test]
fn g0_claims_remain_pinned_when_alternatives_enter_the_temperature_domain() {
    let mut claims = iso_claims(false);
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, 1e9),
        (CTE, INVERSE_TEMPERATURE_DIMS, 7e-5),
    ] {
        let mut alternate = claim(name, dims, &[(350.0, value), (450.0, value)]);
        alternate.validity = domain(350.0, 450.0);
        claims.push(alternate);
    }
    let original = card(claims);
    let source = LaminateThermalSource::try_new(
        vec![input(Arc::clone(&original), 0.001)],
        ReferenceSurface::Midsurface,
        || false,
    )
    .unwrap();
    let hot = source.at_temperatures(&[400.0], || false).unwrap();
    for reference in source.layers()[0].reference_elasticity().properties() {
        assert_eq!(
            reference.answer().receipt.selected,
            hot.layers()[0]
                .elasticity()
                .property(reference.requirement().name())
                .unwrap()
                .answer()
                .receipt
                .selected
        );
    }
    let expansion = hot.layers()[0].expansion();
    assert_eq!(
        expansion
            .reference()
            .property(CTE)
            .unwrap()
            .answer()
            .receipt
            .selected,
        expansion
            .current()
            .property(CTE)
            .unwrap()
            .answer()
            .receipt
            .selected
    );
    near(hot.prepared().profiles()[0].faces()[0][0], 0.00275, 2e-15);
    assert!(source.at_temperatures(&[425.0], || false).is_err());
    assert!(
        resolve_material_state_point(
            &original,
            &at(400.0),
            &[req(
                "young_modulus",
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive
            ),],
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    // The caller can also pin the original complete set explicitly. The owner
    // partitions that admitted plan without losing or reselecting any claim.
    let mut pinned = input(original, 0.001);
    pinned.selection = MaterialPropertySelection::PinnedByProperty(
        source.layers()[0]
            .reference_material()
            .properties()
            .iter()
            .map(|p| {
                (
                    p.requirement().name().to_owned(),
                    p.answer().receipt.selected,
                )
            })
            .collect(),
    );
    let replay =
        LaminateThermalSource::try_new(vec![pinned], ReferenceSurface::Midsurface, || false)
            .unwrap();
    assert_eq!(
        replay
            .at_temperatures(&[400.0], || false)
            .unwrap()
            .prepared(),
        hot.prepared()
    );
}

#[test]
fn g1_rotated_orthotropic_plane_needs_no_invented_third_axis_elastic_constants() {
    let mut claims = vec![claim(DENSITY_PROPERTY, Density::DIMS, &[(300.0, 1600.0)])];
    for (name, lo, hi) in [
        ("young_modulus_1", 100e9, 80e9),
        ("young_modulus_2", 40e9, 20e9),
        ("shear_modulus_12", 12e9, 12e9),
        ("shear_modulus_31", 6e9, 6e9),
        ("shear_modulus_23", 4e9, 4e9),
    ] {
        claims.push(claim(name, Pressure::DIMS, &[(300.0, lo), (400.0, hi)]));
    }
    claims.push(claim("poisson_ratio_12", Dims::NONE, &[(300.0, 0.2)]));
    for (name, value) in CTE3.into_iter().zip([2e-5, -1e-5, 3e-5]) {
        claims.push(claim(name, INVERSE_TEMPERATURE_DIMS, &[(300.0, value)]));
    }
    claims.push(claim(CTE, INVERSE_TEMPERATURE_DIMS, &[(300.0, 1e-5)]));
    let original = card(claims);
    assert!(original.claims_for("young_modulus_3").is_empty());
    assert!(original.claims_for("poisson_ratio_13").is_empty());
    let mut ply = input(original, 0.002);
    ply.elasticity = LaminaElasticSource::OrthotropicPlane12;
    ply.expansion = LaminaExpansionSource::Orthotropic;
    ply.orientation_rad = core::f64::consts::FRAC_PI_4;
    let source =
        LaminateThermalSource::try_new(vec![ply.clone()], ReferenceSurface::Midsurface, || false)
            .unwrap();
    let snapshot = source.at_temperatures(&[400.0], || false).unwrap();
    let prepared = snapshot.prepared();
    let free = prepared.layer_free_strain_in_section(0).unwrap()[0];
    for (actual, expected) in free.into_iter().zip([0.0005, 0.0005, 0.003]) {
        near(actual, expected, 2e-15);
    }
    let denominator = 1.0 - 0.2 * 0.2 * 20e9 / 80e9;
    let s1 = (80e9 * 0.002 + 0.2 * 20e9 * (-0.001)) / denominator;
    let s2 = (0.2 * 20e9 * 0.002 + 20e9 * (-0.001)) / denominator;
    let forces = [0.001 * (s1 + s2), 0.001 * (s1 + s2), 0.001 * (s1 - s2)];
    for (actual, expected) in prepared.membrane_force_n_m().into_iter().zip(forces) {
        near(actual, expected, 2e-15);
    }
    near(
        prepared.restrained_energy_j_m2(),
        0.001 * (0.002 * s1 - 0.001 * s2),
        3e-15,
    );
    match &prepared.section().layers()[0].law {
        LayerLaw::PlaneStress {
            transverse_shear_pa,
            ..
        } => assert_eq!(*transverse_shear_pa, [6e9, 4e9]),
        other => panic!("unexpected law {other:?}"),
    }
    assert!(matches!(
        snapshot.layers()[0].expansion(),
        IntegratedLaminaExpansion::Orthotropic(_)
    ));
    // Elastic orthotropy does not silently require orthotropic expansion.
    ply.expansion = LaminaExpansionSource::Isotropic;
    let scalar = LaminateThermalSource::try_new(vec![ply], ReferenceSurface::Midsurface, || false)
        .unwrap()
        .at_temperatures(&[400.0], || false)
        .unwrap();
    assert!(matches!(
        scalar.layers()[0].expansion(),
        IntegratedLaminaExpansion::Isotropic(_)
    ));
}

#[test]
fn g1_bulk_shear_source_roles_and_directional_expansion_are_independent() {
    let mut claims = vec![
        claim(DENSITY_PROPERTY, Density::DIMS, &[(300.0, 1000.0)]),
        claim("bulk", Pressure::DIMS, &[(300.0, 120e9), (400.0, 80e9)]),
        claim("shear", Pressure::DIMS, &[(300.0, 40e9), (400.0, 20e9)]),
    ];
    for (name, value) in CTE3.into_iter().zip([1e-5, 2e-5, -1e-5]) {
        claims.push(claim(name, INVERSE_TEMPERATURE_DIMS, &[(300.0, value)]));
    }
    let mut ply = input(card(claims), 0.001);
    ply.elasticity = LaminaElasticSource::Isotropic(IsotropicElasticParameters::BulkShear {
        bulk_modulus: req(
            "bulk",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
        shear_modulus: req(
            "shear",
            Pressure::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        ),
    });
    ply.expansion = LaminaExpansionSource::Orthotropic;
    let source =
        LaminateThermalSource::try_new(vec![ply], ReferenceSurface::Midsurface, || false).unwrap();
    let current = source.at_temperatures(&[400.0], || false).unwrap();
    let e = 9.0 * 80e9 * 20e9 / (3.0 * 80e9 + 20e9);
    let nu = (3.0 * 80e9 - 2.0 * 20e9) / (2.0 * (3.0 * 80e9 + 20e9));
    near(
        current.prepared().section().layer_integrals()[0].q_bar_pa[0],
        e / (1.0 - nu * nu),
        2e-15,
    );
    assert_eq!(current.layers()[0].elasticity().properties().len(), 2);
    assert!(
        current.layers()[0]
            .elasticity()
            .property("young_modulus")
            .is_none()
    );
    near(current.prepared().profiles()[0].faces()[0][1], 0.002, 2e-15);
}

#[test]
fn g0_source_domain_context_resolution_and_cancellation_refuse_atomically() {
    let ply = input(card(iso_claims(false)), 0.001);
    let source =
        LaminateThermalSource::try_new(vec![ply.clone()], ReferenceSurface::Midsurface, || false)
            .unwrap();
    let good = source.at_temperatures(&[350.0], || false).unwrap();
    let before = good.prepared().clone();
    assert!(matches!(
        source.at_temperatures(&[], || false),
        Err(LaminateThermalSourceError::TemperatureCount { .. })
    ));
    for t in [0.0, f64::NAN, 401.0] {
        assert!(source.at_temperatures(&[t], || false).is_err());
    }
    let mut polls = 0;
    assert!(matches!(
        source.at_temperatures(&[360.0], || {
            polls += 1;
            polls >= 3
        }),
        Err(LaminateThermalSourceError::Cancelled)
    ));
    assert_eq!(good.prepared(), &before);
    let mut evolving = ply.clone();
    evolving.reference = evolving.reference.with("strain", 0.0).unwrap();
    assert!(matches!(
        LaminateThermalSource::try_new(vec![evolving], ReferenceSurface::Midsurface, || false),
        Err(LaminateThermalSourceError::Invalid { .. })
    ));
    assert!(matches!(
        LaminateThermalSource::try_new(
            vec![ply.clone(); MAX_LAMINA_LAYERS + 1],
            ReferenceSurface::Midsurface,
            || false
        ),
        Err(LaminateThermalSourceError::Section(
            LaminateError::LayerCount { .. }
        ))
    ));
    assert!(matches!(
        LaminateThermalSource::try_new(vec![ply], ReferenceSurface::Midsurface, || true),
        Err(LaminateThermalSourceError::Cancelled)
    ));
}
