//! G0/G3 checks for one-source-of-truth Euler-disc specimen resolution.

use fs_conduction::lumped::{
    BiotGate, LumpedEnthalpyBody, LumpedEnthalpyMarchConfig, LumpedThermalTransport,
    solve_lumped_enthalpy,
};
use fs_conduction::material::CONDUCTIVITY_DIMS;
use fs_conduction::radiation::SURFACE_EMISSIVITY_PROPERTY;
use fs_euler_disc_e2e::specimen::{
    DiscPhaseGeometryRegime, DiscProfileError, DiscProfileSpec, DiscThermalCouplingError,
    PhaseDiscBindingError,
};
use fs_evidence::ValidityDomain;
use fs_exec::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::{MaterialPropertySelection, resolve_isotropic_solid_state_point};
use fs_qty::{Density, Dims, Pressure};
use fs_rep_frep::SquatDiscEdgeTreatment;

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
fn g3_uniform_scaling_preserves_exact_dimensional_mass_laws_for_chamfers() {
    let base = DiscProfileSpec::ChamferedCylinder {
        outer_radius_m: 0.038,
        thickness_m: 0.006,
        chamfer_radial_m: 0.001,
        chamfer_axial_m: 0.0015,
    };
    let scale = 2.5;
    let scaled = DiscProfileSpec::ChamferedCylinder {
        outer_radius_m: 0.038 * scale,
        thickness_m: 0.006 * scale,
        chamfer_radial_m: 0.001 * scale,
        chamfer_axial_m: 0.0015 * scale,
    };
    let (base, scaled) = with_cx(|cx| {
        (
            base.resolve(2_680.0, cx).expect("base chamfer"),
            scaled.resolve(2_680.0, cx).expect("scaled chamfer"),
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

#[test]
fn g0_invalid_profile_relationships_refuse_before_mass_properties_publish() {
    with_cx(|cx| {
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
