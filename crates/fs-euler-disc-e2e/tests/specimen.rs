//! G0/G3 checks for one-source-of-truth Euler-disc specimen resolution.

use fs_blake3::hash_domain;
use fs_conduction::lumped::{
    BiotGate, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig, LumpedThermalTransport,
    solve_lumped_enthalpy,
};
use fs_conduction::material::CONDUCTIVITY_DIMS;
use fs_conduction::radiation::SURFACE_EMISSIVITY_PROPERTY;
use fs_euler_disc_e2e::specimen::{
    DiscPhaseGeometryRegime, DiscProfileError, DiscProfileSpec, DiscThermalCouplingError,
    PhaseDiscBindingError, SolidGeometryEvolutionError, UniformIsotropicFreeExpansionLaw,
};
use fs_evidence::ValidityDomain;
use fs_exec::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, integrate_isotropic_thermal_expansion,
    resolve_isotropic_elastic_state_point, resolve_isotropic_solid_state_point,
};
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_rep_frep::SquatDiscEdgeTreatment;
use fs_solid::TetThermalStrainState;

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x4555_4c45_525f_5052,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn assert_close(actual: f64, expected: f64) {
    let tolerance = 3.0e-10 * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.17e}"
    );
}

fn isotropic_card(chemistry: &str, young_modulus_pa: f64) -> MaterialCard {
    isotropic_card_with_density(chemistry, young_modulus_pa, 8_000.0)
}

fn isotropic_card_with_density(
    chemistry: &str,
    young_modulus_pa: f64,
    density_kg_m3: f64,
) -> MaterialCard {
    isotropic_card_with_density_and_thermal(chemistry, young_modulus_pa, density_kg_m3, None)
}

fn isotropic_card_with_density_and_thermal(
    chemistry: &str,
    young_modulus_pa: f64,
    density_kg_m3: f64,
    thermal: Option<(f64, f64)>,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, density_kg_m3),
        ("young_modulus", Pressure::DIMS, young_modulus_pa),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("yield_stress", Pressure::DIMS, 200.0e6),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("synthetic {chemistry} {name}"),
                    license: "CC0-1.0".to_owned(),
                    artifact: None,
                },
            })
            .expect("synthetic material claim");
    }
    if let Some((conductivity_w_per_m_k, emissivity)) = thermal {
        for (name, dims, value) in [
            (
                "thermal_conductivity",
                CONDUCTIVITY_DIMS,
                conductivity_w_per_m_k,
            ),
            (SURFACE_EMISSIVITY_PROPERTY, Dims::NONE, emissivity),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("T", 290.0, 1_000.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    observations: Vec::new(),
                    provenance: Provenance {
                        source: format!("synthetic {chemistry} {name}"),
                        license: "CC0-1.0".to_owned(),
                        artifact: None,
                    },
                })
                .expect("synthetic thermal material claim");
        }
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: chemistry.to_owned(),
            phase: "solid".to_owned(),
            process: "synthetic".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("synthetic material card")
}

fn isotropic_card_with_density_curve(
    chemistry: &str,
    density_knots: Vec<(f64, f64)>,
) -> MaterialCard {
    let lower_temperature_k = density_knots.first().unwrap().0;
    let upper_temperature_k = density_knots.last().unwrap().0;
    let validity =
        ValidityDomain::unconstrained().with("T", lower_temperature_k, upper_temperature_k);
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("density", Density::DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".to_owned(),
                abscissa_dims: Temperature::DIMS,
                knots: density_knots,
                dims: Density::DIMS,
            },
            validity: validity.clone(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("synthetic {chemistry} density curve"),
                license: "CC0-1.0".to_owned(),
                artifact: None,
            },
        })
        .expect("synthetic density curve");
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, 193.0e9),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("yield_stress", Pressure::DIMS, 200.0e6),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: format!("synthetic {chemistry} {name}"),
                    license: "CC0-1.0".to_owned(),
                    artifact: None,
                },
            })
            .expect("synthetic elastic claim");
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
            ),
            value: PropertyValue::Curve {
                abscissa: "T".to_owned(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (lower_temperature_k, 10.0e-6),
                    (upper_temperature_k, 20.0e-6),
                ],
                dims: INVERSE_TEMPERATURE_DIMS,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("synthetic {chemistry} expansion curve"),
                license: "CC0-1.0".to_owned(),
                artifact: None,
            },
        })
        .expect("synthetic expansion curve");
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: chemistry.to_owned(),
            phase: "solid".to_owned(),
            process: "synthetic-density-curve".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("synthetic temperature-dependent material card")
}

#[test]
fn g0_solid_cylinder_resolution_matches_closed_form_mass_and_inertia() {
    let radius = 0.038;
    let thickness = 0.006;
    let density = 7_800.0;
    let resolved = with_cx(|cx| {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: radius,
            thickness_m: thickness,
            edge_treatment: SquatDiscEdgeTreatment::Sharp,
        }
        .resolve(density, cx)
        .expect("admitted sharp cylinder")
    });
    let volume = core::f64::consts::PI * radius.powi(2) * thickness;
    let mass = density * volume;
    assert_close(resolved.mass_properties.volume, volume);
    assert_close(resolved.mass_properties.mass, mass);
    assert_close(resolved.mass_properties.center_of_mass.z, 0.0);
    assert_close(
        resolved.mass_properties.principal_inertia.axial,
        0.5 * mass * radius.powi(2),
    );
    assert_close(
        resolved.mass_properties.principal_inertia.transverse,
        mass * (3.0 * radius.powi(2) + thickness.powi(2)) / 12.0,
    );
}

#[test]
fn g0_lumped_thermal_geometry_uses_the_same_profile_area_and_volume() {
    let radius = 0.038;
    let thickness = 0.006;
    let spec = DiscProfileSpec::SolidCylinder {
        outer_radius_m: radius,
        thickness_m: thickness,
        edge_treatment: SquatDiscEdgeTreatment::Sharp,
    };
    let (steel, tungsten, thermal, tungsten_thermal) = with_cx(|cx| {
        let steel = spec.resolve(7_800.0, cx).expect("steel profile");
        let tungsten = spec.resolve(19_250.0, cx).expect("tungsten profile");
        let thermal = steel.thermal_geometry(cx).expect("thermal geometry");
        let tungsten_thermal = tungsten
            .thermal_geometry(cx)
            .expect("same-shape thermal geometry");
        (steel, tungsten, thermal, tungsten_thermal)
    });

    let expected_volume = core::f64::consts::PI * radius * radius * thickness;
    let expected_area = 2.0 * core::f64::consts::PI * radius * (radius + thickness);
    assert_close(thermal.volume_m3, expected_volume);
    assert_close(thermal.surface_area_m2, expected_area);
    assert_close(
        thermal.characteristic_length_m,
        expected_volume / expected_area,
    );
    assert_eq!(thermal, tungsten_thermal);
    assert_ne!(
        steel.content_identities().profile,
        tungsten.content_identities().profile
    );
}

#[test]
fn g0_material_specimen_uses_card_density_without_aliasing_equal_density_materials() {
    let point = QueryPoint::new().with("T", 293.15).expect("state point");
    let copper = resolve_isotropic_solid_state_point(
        &isotropic_card("copper-c110", 117.0e9),
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .expect("copper state");
    let steel = resolve_isotropic_solid_state_point(
        &isotropic_card("stainless-316l", 193.0e9),
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .expect("steel state");
    let spec = DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.038,
        thickness_m: 0.006,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    };
    let (copper, steel) = with_cx(|cx| {
        (
            spec.resolve_with_material_state(&copper, cx)
                .expect("copper specimen"),
            spec.resolve_with_material_state(&steel, cx)
                .expect("steel specimen"),
        )
    });
    assert_eq!(
        copper.profile.content_identities(),
        steel.profile.content_identities(),
        "equal geometry and density retain equal mass geometry"
    );
    assert_ne!(
        copper.identity, steel.identity,
        "different complete material states cannot alias merely because density matches"
    );
    assert_eq!(copper.material.young_modulus_pa(), 117.0e9);
    assert_eq!(steel.material.young_modulus_pa(), 193.0e9);
}

#[test]
fn g0_phase_state_invalidates_fixed_solid_mechanics_at_first_liquid_fraction() {
    let card = isotropic_card("phase-bound-test", 193.0e9);
    let point = QueryPoint::new().with("T", 293.15).expect("state point");
    let mechanical = resolve_isotropic_solid_state_point(
        &card,
        &point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .expect("solid mechanics state");
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 293.15,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 8_000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 100_000.0,
                temperature_k: 293.15,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 7_500.0,
            },
        ],
    )
    .expect("phase curve");
    let spec = DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.038,
        thickness_m: 0.006,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    };
    let solid_phase = curve.state_at_specific_enthalpy(0.0).unwrap();
    let solid = with_cx(|cx| spec.resolve_with_phase_state(solid_phase, cx).unwrap());
    let bound = solid.try_bind_fixed_solid(&mechanical).unwrap();
    assert_eq!(
        bound.profile.content_identities(),
        solid.profile.content_identities()
    );

    let partially_liquid = curve.state_at_specific_enthalpy(1.0).unwrap();
    let evolving = with_cx(|cx| spec.resolve_with_phase_state(partially_liquid, cx).unwrap());
    assert!(matches!(
        evolving.try_bind_fixed_solid(&mechanical),
        Err(PhaseDiscBindingError::EvolvingPhaseRequired {
            liquid_mass_fraction
        }) if liquid_mass_fraction > 0.0
    ));
}

#[test]
fn g0_phase_updates_preserve_mass_and_demand_the_correct_geometry_rung() {
    let card = isotropic_card("mass-conserving-phase-test", 193.0e9);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 293.15,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 8_000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 50_000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 7_900.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 150_000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 7_400.0,
            },
        ],
    )
    .expect("phase curve");
    let reference_state = curve.state_at_specific_enthalpy(0.0).unwrap();
    let specimen = with_cx(|cx| {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        }
        .resolve_with_phase_state(reference_state, cx)
        .expect("reference specimen")
    });

    let reference = specimen.mass_conserving_state(reference_state).unwrap();
    assert_eq!(
        reference.geometry_regime,
        DiscPhaseGeometryRegime::ReferenceGeometry
    );
    assert_eq!(
        reference.required_volume_m3.to_bits(),
        specimen.profile.mass_properties.volume.to_bits()
    );

    let warmer_solid = specimen
        .mass_conserving_state(curve.state_at_specific_enthalpy(25_000.0).unwrap())
        .unwrap();
    assert_eq!(
        warmer_solid.geometry_regime,
        DiscPhaseGeometryRegime::SolidThermomechanicalUpdateRequired
    );
    assert!(warmer_solid.volume_ratio > 1.0);
    assert_close(
        warmer_solid.phase_state.bulk_density_kg_m3() * warmer_solid.required_volume_m3,
        warmer_solid.invariant_mass_kg,
    );
    let expansion_law = UniformIsotropicFreeExpansionLaw::try_new(
        0.05,
        hash_domain(
            "fs-euler-disc-e2e/test/free-isotropic-expansion-authority/v1",
            b"synthetic-isotropic-solid",
        ),
    )
    .expect("bounded isotropic free-expansion law");
    let mut tampered_solid = warmer_solid;
    tampered_solid.volume_ratio *= 1.01;
    assert!(matches!(
        with_cx(|cx| specimen.resolve_uniform_isotropic_free_expansion(
            tampered_solid,
            expansion_law,
            cx,
        )),
        Err(SolidGeometryEvolutionError::PhaseStateMismatch)
    ));
    let evolved = with_cx(|cx| {
        specimen
            .resolve_uniform_isotropic_free_expansion(warmer_solid, expansion_law, cx)
            .expect("warmer solid geometry")
    });
    let expected_scale = warmer_solid.volume_ratio.cbrt();
    assert_close(evolved.linear_scale(), expected_scale);
    assert_close(
        evolved.profile().dimensions.outer_radius_m,
        expected_scale * specimen.profile.dimensions.outer_radius_m,
    );
    assert_close(
        evolved.profile().dimensions.thickness_m,
        expected_scale * specimen.profile.dimensions.thickness_m,
    );
    assert_close(
        evolved.profile().mass_properties.mass,
        specimen.profile.mass_properties.mass,
    );
    assert_close(
        evolved.profile().mass_properties.principal_inertia.axial,
        expected_scale.powi(2) * specimen.profile.mass_properties.principal_inertia.axial,
    );
    let too_narrow_law = UniformIsotropicFreeExpansionLaw::try_new(
        1.0e-6,
        hash_domain(
            "fs-euler-disc-e2e/test/free-isotropic-expansion-authority/v1",
            b"too-narrow",
        ),
    )
    .unwrap();
    assert!(matches!(
        with_cx(|cx| specimen.resolve_uniform_isotropic_free_expansion(
            warmer_solid,
            too_narrow_law,
            cx,
        )),
        Err(SolidGeometryEvolutionError::LawValidityExceeded { .. })
    ));

    let partially_liquid = specimen
        .mass_conserving_state(curve.state_at_specific_enthalpy(100_000.0).unwrap())
        .unwrap();
    assert_eq!(
        partially_liquid.geometry_regime,
        DiscPhaseGeometryRegime::EvolvingFreeSurfaceRequired
    );
    assert!(partially_liquid.phase_state.liquid_mass_fraction() > 0.0);
    assert_close(
        partially_liquid.phase_state.bulk_density_kg_m3() * partially_liquid.required_volume_m3,
        partially_liquid.invariant_mass_kg,
    );
    assert!(matches!(
        with_cx(|cx| specimen.resolve_uniform_isotropic_free_expansion(
            partially_liquid,
            expansion_law,
            cx,
        )),
        Err(SolidGeometryEvolutionError::EvolvingFreeSurfaceRequired { .. })
    ));
}

#[test]
fn g0_evolved_solid_rebinds_same_card_state_for_contact_modes_and_sound() {
    let card = isotropic_card_with_density_curve(
        "temperature-dependent-solid",
        vec![(293.15, 8_000.0), (500.0, 7_900.0)],
    );
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 293.15,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 8_000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 50_000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 7_900.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 100_000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 7_800.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 150_000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 7_500.0,
            },
        ],
    )
    .expect("solid phase curve");
    let reference = curve.state_at_specific_enthalpy(0.0).unwrap();
    let warmer = curve.state_at_specific_enthalpy(50_000.0).unwrap();
    let specimen = with_cx(|cx| {
        DiscProfileSpec::OuterFilletedAnnularCylinder {
            outer_radius_m: 0.038,
            inner_radius_m: 0.012,
            thickness_m: 0.006,
            outer_fillet_radius_m: 0.001,
        }
        .resolve_with_phase_state(reference, cx)
        .expect("reference specimen")
    });
    let evolved_state = specimen.mass_conserving_state(warmer).unwrap();
    let law = UniformIsotropicFreeExpansionLaw::try_new(
        0.05,
        hash_domain(
            "fs-euler-disc-e2e/test/free-isotropic-expansion-authority/v1",
            b"temperature-dependent-solid",
        ),
    )
    .unwrap();
    let evolved = with_cx(|cx| {
        specimen
            .resolve_uniform_isotropic_free_expansion(evolved_state, law, cx)
            .expect("evolved solid")
    });
    let warm_point = QueryPoint::new().with("T", 500.0).unwrap();
    let solid = resolve_isotropic_solid_state_point(
        &card,
        &warm_point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .expect("warm solid state");
    let elastic = resolve_isotropic_elastic_state_point(
        &card,
        &warm_point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .expect("warm elastic state");
    let cold_point = QueryPoint::new().with("T", 293.15).unwrap();
    let expansion = integrate_isotropic_thermal_expansion(
        &card,
        &cold_point,
        &warm_point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .expect("one evidence-bearing expansion curve spans cold to warm");
    let thermal_strain =
        TetThermalStrainState::try_from_isotropic_expansion(&elastic, &expansion, 0.05)
            .expect("the integrated path and tangent elasticity share one current state");
    assert_eq!(
        thermal_strain.free_strain_mandel(),
        [
            expansion.free_linear_strain(),
            expansion.free_linear_strain(),
            expansion.free_linear_strain(),
            0.0,
            0.0,
            0.0,
        ]
    );
    let contact_specimen = evolved
        .try_bind_isotropic_solid(&solid)
        .expect("contact specimen");
    let acoustic_specimen = evolved
        .try_bind_isotropic_elastic(&elastic)
        .expect("structural-acoustic specimen");
    assert_eq!(
        contact_specimen.profile.content_identities(),
        evolved.profile().content_identities()
    );
    assert_eq!(
        acoustic_specimen.profile.content_identities(),
        evolved.profile().content_identities()
    );
    assert_eq!(
        acoustic_specimen.material_card_identity,
        card.content_hash()
    );

    let stale_elastic = resolve_isotropic_elastic_state_point(
        &card,
        &cold_point,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(matches!(
        evolved.try_bind_isotropic_elastic(&stale_elastic),
        Err(PhaseDiscBindingError::TemperatureMismatch { .. })
    ));
}

#[test]
fn g1_hot_ambient_phase_march_is_bound_to_exact_disc_geometry_and_escalates_at_melting() {
    let card = isotropic_card_with_density_and_thermal(
        "lead-like-thermal-coupling",
        16.0e9,
        11_340.0,
        Some((35.0, 0.8)),
    );
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 293.15,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 11_340.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 53_000.0,
                temperature_k: 600.6,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 10_900.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 76_000.0,
                temperature_k: 600.6,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 10_650.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 150_000.0,
                temperature_k: 1_000.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 10_200.0,
            },
        ],
    )
    .expect("lead-like phase curve");
    let reference_state = curve.state_at_specific_enthalpy(0.0).unwrap();
    let specimen = with_cx(|cx| {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.010,
            thickness_m: 0.001,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.0002 },
        }
        .resolve_with_phase_state(reference_state, cx)
        .expect("lead-like reference specimen")
    });
    let thermal = with_cx(|cx| specimen.profile.thermal_geometry(cx).unwrap());
    let transport = LumpedThermalTransport::from_material_card(
        &card,
        "thermal_conductivity",
        &[293.15, 600.6, 1_000.0],
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("same-card thermal transport");
    let body = LumpedEnthalpyBody::try_new_with_transport(
        "same resolved lead-like disc",
        specimen.profile.mass_properties.mass,
        thermal.surface_area_m2,
        20.0,
        thermal.characteristic_length_m,
        transport,
        &curve,
    )
    .expect("Biot-gated body inputs");
    let config = LumpedEnthalpyMarchConfig {
        initial_specific_enthalpy_j_kg: 0.0,
        ambient_temperature_k: 1_200.0,
        internal_power_w: 0.0,
        duration_s: 4.0,
        maximum_step_s: 0.01,
        maximum_steps: 400,
        enthalpy_tolerance_j_kg: 1.0e-8,
    };
    let march = with_cx(|cx| {
        solve_lumped_enthalpy(cx, &body, BiotGate::corpus_default(), config)
            .expect("hot-ambient enthalpy march")
    });
    let coupled = with_cx(|cx| {
        specimen
            .bind_lumped_enthalpy_march(&body, &march, cx)
            .expect("exact specimen/thermal binding")
    });

    let initial = coupled.samples.first().unwrap().mass_conserving_state;
    let final_state = coupled.samples.last().unwrap().mass_conserving_state;
    assert_eq!(
        initial.geometry_regime,
        DiscPhaseGeometryRegime::ReferenceGeometry
    );
    assert_eq!(
        final_state.geometry_regime,
        DiscPhaseGeometryRegime::EvolvingFreeSurfaceRequired
    );
    assert!(final_state.phase_state.liquid_mass_fraction() > 0.0);
    assert!(final_state.volume_ratio > 1.0);
    assert!(coupled.maximum_biot <= 0.1);
    assert!(coupled.samples.iter().all(|sample| {
        let state = sample.mass_conserving_state;
        (state.phase_state.bulk_density_kg_m3() * state.required_volume_m3
            - state.invariant_mass_kg)
            .abs()
            <= 1.0e-12
    }));

    let surrogate = LumpedEnthalpyBody::try_new_with_transport(
        "wrong-area surrogate",
        specimen.profile.mass_properties.mass,
        0.5 * thermal.surface_area_m2,
        20.0,
        thermal.characteristic_length_m,
        body.transport().clone(),
        &curve,
    )
    .unwrap();
    let surrogate_march = with_cx(|cx| {
        solve_lumped_enthalpy(
            cx,
            &surrogate,
            BiotGate::corpus_default(),
            LumpedEnthalpyMarchConfig {
                duration_s: 0.0,
                maximum_steps: 1,
                ..config
            },
        )
        .unwrap()
    });
    assert!(matches!(
        with_cx(|cx| specimen.bind_lumped_enthalpy_march(&surrogate, &surrogate_march, cx)),
        Err(DiscThermalCouplingError::SpecimenQuantityMismatch {
            field: "surface_area_m2",
            ..
        })
    ));
}

#[test]
fn g0_annular_cylinder_resolution_uses_bore_for_mass_not_only_inertia() {
    let outer = 0.038;
    let inner = 0.021;
    let thickness = 0.006;
    let density = 2_680.0;
    let resolved = with_cx(|cx| {
        DiscProfileSpec::AnnularCylinder {
            outer_radius_m: outer,
            inner_radius_m: inner,
            thickness_m: thickness,
        }
        .resolve(density, cx)
        .expect("admitted annular cylinder")
    });
    let volume = core::f64::consts::PI * (outer.powi(2) - inner.powi(2)) * thickness;
    let mass = density * volume;
    assert!(
        !resolved.chart.construction_certificate().touches_axis,
        "an annular profile must retain its bore rather than gain an axis closure"
    );
    assert_close(resolved.mass_properties.volume, volume);
    assert_close(resolved.mass_properties.mass, mass);
    assert_close(
        resolved.mass_properties.principal_inertia.axial,
        0.5 * mass * (outer.powi(2) + inner.powi(2)),
    );
    assert_close(
        resolved.mass_properties.principal_inertia.transverse,
        mass * (3.0 * (outer.powi(2) + inner.powi(2)) + thickness.powi(2)) / 12.0,
    );
}

#[test]
fn g0_outer_filleted_annulus_uses_one_chart_for_equal_mass_geometry_and_inertia() {
    let outer = 0.038;
    let inner = 0.021;
    let thickness = 0.006;
    let fillet = 0.001;
    let solid_density = 7_800.0;
    let solid = DiscProfileSpec::SolidCylinder {
        outer_radius_m: outer,
        thickness_m: thickness,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: fillet },
    };
    let annulus = DiscProfileSpec::OuterFilletedAnnularCylinder {
        outer_radius_m: outer,
        inner_radius_m: inner,
        thickness_m: thickness,
        outer_fillet_radius_m: fillet,
    };
    let (solid, unit_density_annulus) = with_cx(|cx| {
        (
            solid
                .resolve(solid_density, cx)
                .expect("filleted solid baseline"),
            annulus
                .resolve(1.0, cx)
                .expect("unit-density filleted annulus"),
        )
    });
    let equal_mass_density =
        solid.mass_properties.mass / unit_density_annulus.mass_properties.volume;
    let ring = with_cx(|cx| {
        annulus
            .resolve(equal_mass_density, cx)
            .expect("equal-mass filleted annulus")
    });

    assert!(
        !ring.chart.construction_certificate().touches_axis,
        "the equal-mass control must retain a real through bore"
    );
    assert!(
        ring.chart
            .segments()
            .iter()
            .filter(|segment| matches!(segment, fs_rep_frep::MeridianSegment::Arc { .. }))
            .count()
            == 2,
        "the ring's outer contact geometry must be circular rather than sharp"
    );
    assert!(
        equal_mass_density > solid_density,
        "equal mass is achieved by the declared density control, not an inertia multiplier"
    );
    assert_close(ring.mass_properties.mass, solid.mass_properties.mass);
    assert_close(ring.mass_properties.center_of_mass.z, 0.0);
    assert!(
        ring.mass_properties.principal_inertia.axial
            > solid.mass_properties.principal_inertia.axial
    );
    assert!(
        ring.mass_properties
            .principal_inertia
            .transverse
            .is_finite()
    );
}

#[test]
fn g0_true_symmetric_bicone_matches_closed_forms_and_is_centered() {
    let radius = 0.038;
    let thickness = 0.006;
    let density = 2_680.0;
    let resolved = with_cx(|cx| {
        DiscProfileSpec::SymmetricTapered {
            outer_radius_m: radius,
            face_radius_m: 0.0,
            thickness_m: thickness,
        }
        .resolve(density, cx)
        .expect("admitted symmetric bicone")
    });
    // The profile is two equal right cones, each height thickness / 2.
    let volume = core::f64::consts::PI * radius.powi(2) * thickness / 3.0;
    let mass = density * volume;
    assert_close(resolved.mass_properties.volume, volume);
    assert_close(resolved.mass_properties.mass, mass);
    assert_close(resolved.mass_properties.center_of_mass.z, 0.0);
    assert_close(
        resolved.mass_properties.principal_inertia.axial,
        0.3 * mass * radius.powi(2),
    );
    assert_close(
        resolved.mass_properties.principal_inertia.transverse,
        mass * (6.0 * radius.powi(2) + thickness.powi(2)) / 40.0,
    );
    assert!(
        resolved.chart.construction_certificate().touches_axis,
        "the zero-face-radius bicone must retain its two axis tips"
    );
}

#[test]
fn g0_full_height_chamfers_omit_the_collapsed_cylindrical_band() {
    let resolved = with_cx(|cx| {
        DiscProfileSpec::ChamferedCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            chamfer_radial_m: 0.001,
            chamfer_axial_m: 0.003,
        }
        .resolve(2_680.0, cx)
        .expect("full-height chamfers form a valid double-frustum rather than a zero line")
    });
    let certificate = resolved.chart.construction_certificate();
    assert_eq!(certificate.input_feature_count, 5);
    assert_eq!(certificate.surfaced_feature_count, 4);
    assert!(resolved.mass_properties.mass > 0.0);
}

#[test]
fn g0_zero_fillet_is_the_same_resolved_sharp_profile() {
    let sharp = with_cx(|cx| {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::Sharp,
        }
        .resolve(7_800.0, cx)
        .expect("sharp")
    });
    let zero_fillet = with_cx(|cx| {
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.0 },
        }
        .resolve(7_800.0, cx)
        .expect("zero fillet")
    });
    assert_eq!(sharp.identity, zero_fillet.identity);
    assert_eq!(sharp.chart.segments(), zero_fillet.chart.segments());
    assert_eq!(
        sharp.mass_properties.mass.to_bits(),
        zero_fillet.mass_properties.mass.to_bits()
    );
    assert_eq!(
        sharp.mass_properties.principal_inertia,
        zero_fillet.mass_properties.principal_inertia
    );
    assert_eq!(
        sharp.content_identities(),
        zero_fillet.content_identities(),
        "semantically identical resolved charts and densities need one durable identity"
    );
}

#[test]
fn g0_strong_profile_identities_bind_complete_chart_density_and_mass_inputs() {
    let spec = DiscProfileSpec::SolidCylinder {
        outer_radius_m: 0.038,
        thickness_m: 0.006,
        edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
    };
    let (steel, replay, tungsten) = with_cx(|cx| {
        (
            spec.resolve(7_800.0, cx).expect("steel profile"),
            spec.resolve(7_800.0, cx).expect("steel replay"),
            spec.resolve(19_250.0, cx).expect("tungsten profile"),
        )
    });
    let steel_ids = steel.content_identities();
    let replay_ids = replay.content_identities();
    let tungsten_ids = tungsten.content_identities();

    assert_eq!(steel_ids, replay_ids, "re-resolution must be deterministic");
    assert_eq!(steel_ids.chart, tungsten_ids.chart);
    assert_ne!(steel_ids.profile, tungsten_ids.profile);
    assert_ne!(steel_ids.mass_properties, tungsten_ids.mass_properties);
    for identity in [
        steel_ids.chart,
        steel_ids.profile,
        steel_ids.mass_properties,
    ] {
        assert!(identity.as_bytes().iter().any(|byte| *byte != 0));
    }

    let chamfer = with_cx(|cx| {
        DiscProfileSpec::ChamferedCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            chamfer_radial_m: 0.001,
            chamfer_axial_m: 0.001,
        }
        .resolve(7_800.0, cx)
        .expect("chamfered profile")
    });
    assert_ne!(steel_ids.chart, chamfer.content_identities().chart);
}

#[test]
fn g3_uniform_scaling_preserves_dimensional_mass_laws_for_every_profile_family() {
    let profiles = [
        DiscProfileSpec::SolidCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            edge_treatment: SquatDiscEdgeTreatment::CircularFillet { radius: 0.001 },
        },
        DiscProfileSpec::AnnularCylinder {
            outer_radius_m: 0.038,
            inner_radius_m: 0.012,
            thickness_m: 0.006,
        },
        DiscProfileSpec::OuterFilletedAnnularCylinder {
            outer_radius_m: 0.038,
            inner_radius_m: 0.012,
            thickness_m: 0.006,
            outer_fillet_radius_m: 0.001,
        },
        DiscProfileSpec::SymmetricTapered {
            outer_radius_m: 0.038,
            face_radius_m: 0.020,
            thickness_m: 0.006,
        },
        DiscProfileSpec::ChamferedCylinder {
            outer_radius_m: 0.038,
            thickness_m: 0.006,
            chamfer_radial_m: 0.001,
            chamfer_axial_m: 0.0015,
        },
    ];
    let scale = 2.5;
    for spec in profiles {
        let scaled_spec = spec.uniformly_scaled(scale).expect("scaled profile");
        let (base, scaled) = with_cx(|cx| {
            (
                spec.resolve(2_680.0, cx).expect("base profile"),
                scaled_spec.resolve(2_680.0, cx).expect("scaled profile"),
            )
        });
        assert_close(
            scaled.mass_properties.volume,
            scale.powi(3) * base.mass_properties.volume,
        );
        assert_close(
            scaled.mass_properties.mass,
            scale.powi(3) * base.mass_properties.mass,
        );
        assert_close(
            scaled.mass_properties.principal_inertia.axial,
            scale.powi(5) * base.mass_properties.principal_inertia.axial,
        );
        assert_close(
            scaled.mass_properties.principal_inertia.transverse,
            scale.powi(5) * base.mass_properties.principal_inertia.transverse,
        );
    }
}

#[test]
fn g0_invalid_profile_relationships_refuse_before_mass_properties_publish() {
    with_cx(|cx| {
        assert!(matches!(
            DiscProfileSpec::SolidCylinder {
                outer_radius_m: 0.038,
                thickness_m: 0.006,
                edge_treatment: SquatDiscEdgeTreatment::Sharp,
            }
            .uniformly_scaled(f64::INFINITY),
            Err(DiscProfileError::InvalidParameter {
                field: "linear_scale",
                ..
            })
        ));
        assert!(matches!(
            DiscProfileSpec::AnnularCylinder {
                outer_radius_m: 0.038,
                inner_radius_m: 0.038,
                thickness_m: 0.006,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidRelationship { .. })
        ));
        assert!(matches!(
            DiscProfileSpec::SymmetricTapered {
                outer_radius_m: 0.038,
                face_radius_m: 0.038,
                thickness_m: 0.006,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidRelationship { .. })
        ));
        assert!(matches!(
            DiscProfileSpec::ChamferedCylinder {
                outer_radius_m: 0.038,
                thickness_m: 0.006,
                chamfer_radial_m: 0.001,
                chamfer_axial_m: 0.0031,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidRelationship { .. })
        ));
        assert!(matches!(
            DiscProfileSpec::OuterFilletedAnnularCylinder {
                outer_radius_m: 0.038,
                inner_radius_m: 0.021,
                thickness_m: 0.006,
                outer_fillet_radius_m: 0.0,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidParameter {
                field: "outer_fillet_radius_m",
                ..
            })
        ));
        assert!(matches!(
            DiscProfileSpec::OuterFilletedAnnularCylinder {
                outer_radius_m: 0.038,
                inner_radius_m: 0.021,
                thickness_m: 0.006,
                outer_fillet_radius_m: 0.0031,
            }
            .resolve(2_680.0, cx),
            Err(DiscProfileError::InvalidRelationship { .. })
        ));
        assert!(matches!(
            DiscProfileSpec::SolidCylinder {
                outer_radius_m: 0.038,
                thickness_m: 0.006,
                edge_treatment: SquatDiscEdgeTreatment::Sharp,
            }
            .resolve(0.0, cx),
            Err(DiscProfileError::InvalidParameter {
                field: "density_kg_per_m3",
                ..
            })
        ));
    });
}
