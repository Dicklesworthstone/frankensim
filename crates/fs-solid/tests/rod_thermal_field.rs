//! G1/G3 sourced temperature -> local constitutive response -> shared equilibrium.
//! Synthetic coefficients exercise wiring and equations, not material calibration.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use fs_qty::{Density, Pressure, Temperature};
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};
use fs_solid::rod::{TipLoad, WorldTipLoad};

fn at(temperature: f64, expansion_scales: &[f64]) -> (ConcentricShearSection, ConcentricThermoelasticSection) {
    let point = QueryPoint::new().with("T", temperature).unwrap();
    let reference = QueryPoint::new().with("T", 300.0).unwrap();
    let mut layers = Vec::new();
    let mut paths = Vec::new();
    for (i, scale) in expansion_scales.iter().enumerate() {
        let mut claims = ClaimSet::new();
        for (name, lo, hi, dims) in [
            ("density", 1000.0, 1000.0, Density::DIMS),
            ("young_modulus", 1e6, 0.5e6, Pressure::DIMS),
            ("shear_modulus", 4e5, 2e5, Pressure::DIMS),
            ("linear_thermal_expansion_coefficient", scale * 1e-5, scale * 3e-5, INVERSE_TEMPERATURE_DIMS),
        ] {
            let (value, interpolation) = if name == "density" {
                (PropertyValue::Scalar { value: lo, dims }, InterpolationPolicy::ConstantWithinValidity)
            } else {
                (PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, lo), (400.0, hi)], dims }, InterpolationPolicy::LinearInside)
            };
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims), value, interpolation,
                validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                uncertainty: UncertaintyModel::Unstated, observations: vec![],
                provenance: Provenance { source: "synthetic thermoelastic rod".into(), license: "test-only".into(), artifact: None },
            }).unwrap();
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: format!("synthetic layer {i}"),
            phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
        let requirements = [("density", Density::DIMS), ("young_modulus", Pressure::DIMS),
            ("shear_modulus", Pressure::DIMS)].map(|(name, dims)|
                ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
        let material = resolve_material_state_point(&card, &point, &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        paths.push(integrate_isotropic_thermal_expansion(&card, &reference, &point,
            MaterialPropertySelection::SingleClaimOnly).unwrap());
        layers.push(ConcentricMaterialLayer {
            outer_radius_m: 0.05 * (i + 1) as f64 / expansion_scales.len() as f64, material,
        });
    }
    let elastic = ConcentricElasticSection::try_new(0.0, layers).unwrap();
    let thermal = ConcentricThermoelasticSection::try_new(elastic.clone(), paths).unwrap();
    (elastic.with_shear_modulus("shear_modulus", 0.8).unwrap(), thermal)
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!((actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "actual={actual}, expected={expected}");
}

#[test]
fn local_heating_changes_expansion_and_compliance_in_the_same_solve() {
    let (sa, ta) = at(350.0, &[1.0]);
    let (sb, tb) = at(400.0, &[1.0]);
    let plain = SpatialMaterialRod::try_straight(1.0, vec![sa, sb]).unwrap();
    let mut heated = plain.with_thermal_sections(vec![ta, tb], 0.02, &CancelGate::new()).unwrap();
    heated.solve_static(&TipLoad { force: [0.5, 0.0, 0.0], moment: [0.0; 3] },
        1, 1e-7, &CancelGate::new()).unwrap();
    let area = core::f64::consts::PI * 0.05_f64.powi(2);
    let mut extension = 0.0;
    for (i, dt) in [50.0, 100.0].into_iter().enumerate() {
        // Independently integrate the DECLARED linear alpha(T), not alpha(T_end)*dT.
        let eigenstrain = 1e-5 * dt + 0.5 * 2e-7 * dt * dt;
        let young = 1e6 - 5000.0 * dt;
        let expected = eigenstrain + 0.5 / (young * area);
        near(heated.segment_strains(i).unwrap().0[0], expected, 3e-5);
        near(heated.segment_resultants(i).unwrap().0[0], 0.5, 3e-5);
        near(heated.thermal_axial_response(i).unwrap().unwrap().layers()[0].thermal_strain(), eigenstrain, 1e-12);
        extension += 0.5 * expected;
    }
    near(heated.positions_m()[2][0] - 1.0, extension, 3e-5);
    assert_eq!(plain.positions_m()[2], [1.0, 0.0, 0.0]);
    assert!(plain.thermal_sections().is_none());
}

#[test]
fn world_held_torsion_uses_each_heated_segments_actual_shear_modulus() {
    let gate = CancelGate::new();
    let (sa, ta) = at(350.0, &[1.0]);
    let (sb, tb) = at(400.0, &[1.0]);
    let sources = vec![sa, sb];
    let mut rod = SpatialMaterialRod::try_straight(1.0, sources.clone())
        .unwrap()
        .with_thermal_sections(vec![ta, tb], 0.02, &gate)
        .unwrap();
    let load = WorldTipLoad {
        force: [0.0; 3],
        moment: [0.03, 0.0, 0.0],
    };
    rod.solve_static_world(&load, 1, 1e-9, &gate).unwrap();
    let polar = core::f64::consts::PI * 0.05_f64.powi(4) / 2.0;
    let mut expected_twist = 0.0;
    for (i, shear) in [3e5, 2e5].into_iter().enumerate() {
        let twist_rate = load.moment[0] / (shear * polar);
        near(rod.segment_strains(i).unwrap().1[0], twist_rate, 2e-7);
        near(
            rod.segment_resultants(i).unwrap().1[0],
            load.moment[0],
            2e-7,
        );
        expected_twist += 0.5 * twist_rate;
    }
    near(
        rod.frames()[2].log().unwrap().angular.x,
        expected_twist,
        2e-7,
    );
    assert_eq!(rod.sections(), sources);
    assert!(
        rod.support_reactions_world(&load, &gate)
            .unwrap()
            .free_residual_norm
            < 1e-9
    );
}

#[test]
fn free_layered_member_retains_opposing_stresses_and_nonzero_stored_energy() {
    let (shear, thermal) = at(400.0, &[1.0, 2.0]);
    let mut rod = SpatialMaterialRod::try_straight(0.5, vec![shear]).unwrap()
        .with_thermal_sections(vec![thermal], 0.02, &CancelGate::new()).unwrap();
    rod.solve_static(&TipLoad::default(), 1, 1e-7, &CancelGate::new()).unwrap();
    let response = rod.thermal_axial_response(0).unwrap().unwrap();
    let expected_strain = 0.25 * 0.002 + 0.75 * 0.004;
    near(rod.positions_m()[1][0] - 0.5, 0.5 * expected_strain, 3e-5);
    assert!(response.force_n().abs() < 1e-6);
    near(response.layers()[0].stress_pa(), 750.0, 3e-5);
    near(response.layers()[1].stress_pa(), -250.0, 3e-5);
    assert!(rod.energy_j().unwrap() > 0.0, "zero resultant is not zero layer energy");
    near(rod.energy_j().unwrap(), 0.5 * response.strain_energy_per_length_j_m(), 1e-12);
}

#[test]
fn the_same_temperature_change_updates_torsion_not_only_free_expansion() {
    let mut twists = Vec::new();
    for t in [300.0, 400.0] {
        let (shear, thermal) = at(t, &[1.0]);
        let mut rod = SpatialMaterialRod::try_straight(0.5, vec![shear]).unwrap()
            .with_thermal_sections(vec![thermal], 0.02, &CancelGate::new()).unwrap();
        rod.solve_static(&TipLoad { force: [0.0; 3], moment: [0.01, 0.0, 0.0] },
            1, 1e-7, &CancelGate::new()).unwrap();
        twists.push(rod.segment_strains(0).unwrap().1[0]);
        near(rod.segment_resultants(0).unwrap().1[0], 0.01, 3e-5);
        assert_eq!(rod.thermal_sections().unwrap()[0].section(), rod.sections()[0].elastic());
    }
    near(twists[1] / twists[0], 2.0, 3e-5);
}

#[test]
fn summary_and_layer_outputs_share_all_constitutive_arithmetic_and_refusals() {
    let (_, source) = at(400.0, &[1.0, 2.0]);
    for strain in [-0.005, 0.0, 0.0035, 0.008] {
        let full = source.response_at_axial_strain(strain, 0.02).unwrap();
        let summary = source.axial_force_and_energy_at_strain(strain, 0.02).unwrap();
        assert_eq!(summary.0.to_bits(), full.force_n().to_bits());
        assert_eq!(summary.1.to_bits(), full.strain_energy_per_length_j_m().to_bits());
    }
    for (strain, limit) in [(f64::NAN, 0.02), (0.0, 0.0), (0.03, 0.02), (0.0, 0.001)] {
        assert_eq!(source.response_at_axial_strain(strain, limit).unwrap_err(),
            source.axial_force_and_energy_at_strain(strain, limit).unwrap_err());
    }
}

#[test]
fn stale_sources_excessive_strain_and_cancellation_publish_no_changed_pose() {
    let (shear, thermal) = at(400.0, &[1.0]);
    let (_, stale) = at(350.0, &[1.0]);
    let plain = SpatialMaterialRod::try_straight(0.5, vec![shear]).unwrap();
    let gate = CancelGate::new();
    assert!(plain.with_thermal_sections(vec![stale], 0.02, &gate).is_err());
    assert!(plain.with_thermal_sections(vec![], 0.02, &gate).is_err());
    assert!(plain.with_thermal_sections(vec![thermal.clone()], 0.001, &gate).is_err());
    let mut rod = plain.with_thermal_sections(vec![thermal], 0.02, &gate).unwrap();
    rod.solve_static(&TipLoad::default(), 1, 1e-7, &gate).unwrap();
    let before = rod.positions_m().to_vec();
    let energy = rod.energy_j().unwrap();
    assert!(matches!(rod.solve_static(&TipLoad { force: [1e4, 0.0, 0.0], moment: [0.0; 3] },
        1, 1e-7, &gate), Err(SpatialRodError::Material(_))));
    assert_eq!(rod.positions_m(), before);
    assert_eq!(rod.energy_j().unwrap(), energy);
    gate.request();
    assert_eq!(rod.solve_static(&TipLoad::default(), 1, 1e-7, &gate), Err(SpatialRodError::Cancelled));
    assert_eq!(rod.positions_m(), before);
}

#[test]
fn restraining_nonuniform_heating_produces_compatible_compression_and_support_force() {
    use fs_ga::So3;
    use fs_solid::rod::material_field::supports::RodSupport;
    let (sa, ta) = at(350.0, &[1.0]);
    let (sb, tb) = at(400.0, &[1.0]);
    let gate = CancelGate::new();
    let free = SpatialMaterialRod::try_straight(1.0, vec![sa, sb]).unwrap()
        .with_thermal_sections(vec![ta, tb], 0.02, &gate).unwrap();
    let mut restrained = free.with_supports(&[RodSupport { node: 2,
        position_m: [Some(1.0), Some(0.0), Some(0.0)], frame: Some(So3::identity()) }], &gate).unwrap();
    let load = TipLoad::default();
    restrained.solve_static(&load, 1, 1e-7, &gate).unwrap();
    let area = core::f64::consts::PI * 0.05_f64.powi(2);
    let force = -(0.5 * 0.00075 + 0.5 * 0.002)
        / (0.5 / (0.75e6 * area) + 0.5 / (0.5e6 * area));
    for (i, (young, thermal)) in [(0.75e6, 0.00075), (0.5e6, 0.002)].into_iter().enumerate() {
        near(restrained.segment_strains(i).unwrap().0[0], thermal + force / (young * area), 3e-5);
        near(restrained.segment_resultants(i).unwrap().0[0], force, 3e-5);
    }
    assert_eq!(restrained.positions_m()[2], [1.0, 0.0, 0.0]);
    assert_eq!(restrained.frames()[2], So3::identity());
    let reactions = restrained.support_reactions(&load, &gate).unwrap();
    near(reactions.nodes[0].force_n[0].unwrap(), force, 3e-5);
    assert!(reactions.free_residual_norm < 1e-7);
    assert_eq!(free.positions_m()[1], [0.5, 0.0, 0.0]);
    assert!(free.supports().is_empty());
}

#[test]
fn fully_prescribed_pose_keeps_real_reactions_when_there_are_no_free_dofs() {
    use fs_ga::{So3, So3Tangent, Vec3};
    use fs_solid::rod::material_field::supports::RodSupport;
    let (shear, _) = at(300.0, &[1.0]);
    let axial = shear.elastic().axial_stiffness_n();
    let torsion = shear.torsional_stiffness_n_m2();
    let gate = CancelGate::new();
    let pose = So3::identity().body_plus(So3Tangent::new(Vec3::new(0.005, 0.0, 0.0))).unwrap();
    let mut rod = SpatialMaterialRod::try_straight(0.5, vec![shear]).unwrap()
        .with_supports(&[RodSupport { node: 1, position_m: [Some(0.501), Some(0.0), Some(0.0)],
            frame: Some(pose) }], &gate).unwrap();
    let load = TipLoad { force: [0.5, 0.0, 0.0], moment: [0.01, 0.0, 0.0] };
    assert_eq!(rod.solve_static(&load, 2, 1e-7, &gate).unwrap(), vec![0.0, 0.0]);
    assert_eq!(rod.positions_m()[1], [0.501, 0.0, 0.0]);
    assert_eq!(rod.frames()[1], pose);
    let response = rod.support_reactions(&load, &gate).unwrap();
    assert_eq!(response.free_residual_norm, 0.0);
    near(response.nodes[0].force_n[0].unwrap(), axial * 0.002 - 0.5, 3e-5);
    near(response.nodes[0].moment_body_nm.unwrap()[0], torsion * 0.01 - 0.01, 3e-5);
    assert!(rod.energy_j().unwrap() > 0.0);
}

#[test]
fn rotation_clamp_leaves_translation_free_and_retains_the_applied_moment_reaction() {
    use fs_ga::So3;
    use fs_solid::rod::material_field::supports::RodSupport;
    let (shear, _) = at(300.0, &[1.0]);
    let axial = shear.elastic().axial_stiffness_n();
    let gate = CancelGate::new();
    let mut rod = SpatialMaterialRod::try_straight(0.5, vec![shear]).unwrap()
        .with_supports(&[RodSupport { node: 1, position_m: [None; 3],
            frame: Some(So3::identity()) }], &gate).unwrap();
    let load = TipLoad { force: [0.5, 0.0, 0.0], moment: [0.01, 0.0, 0.0] };
    rod.solve_static(&load, 1, 1e-7, &gate).unwrap();
    near(rod.positions_m()[1][0] - 0.5, 0.5 * 0.5 / axial, 3e-5);
    assert_eq!(rod.frames()[1], So3::identity());
    let response = rod.support_reactions(&load, &gate).unwrap();
    assert_eq!(response.nodes[0].force_n, [None; 3]);
    near(response.nodes[0].moment_body_nm.unwrap()[0], -0.01, 3e-5);
    assert!(response.free_residual_norm < 1e-7);
}

#[test]
fn malformed_or_out_of_regime_supports_cannot_mutate_the_material_specimen() {
    use fs_solid::rod::material_field::supports::RodSupport;
    let (shear, thermal) = at(400.0, &[1.0]);
    let gate = CancelGate::new();
    let rod = SpatialMaterialRod::try_straight(0.5, vec![shear]).unwrap()
        .with_thermal_sections(vec![thermal], 0.02, &gate).unwrap();
    let support = RodSupport { node: 1, position_m: [Some(0.5), None, None], frame: None };
    for bad in [RodSupport { node: 0, ..support }, RodSupport { node: 2, ..support },
        RodSupport { position_m: [None; 3], ..support },
        RodSupport { position_m: [Some(f64::NAN), None, None], ..support },
        RodSupport { position_m: [Some(0.6), None, None], ..support }]
    {
        assert!(rod.with_supports(&[bad], &gate).is_err());
    }
    assert!(rod.with_supports(&[support, support], &gate).is_err());
    let bound = rod.with_supports(&[support], &gate).unwrap();
    assert_eq!(bound.supports(), &[support]);
    assert!(bound.with_supports(&[], &gate).unwrap().supports().is_empty());
    gate.request();
    assert!(matches!(rod.with_supports(&[support], &gate), Err(SpatialRodError::Cancelled)));
    assert!(rod.supports().is_empty());
    assert_eq!(rod.positions_m()[1], [0.5, 0.0, 0.0]);
}
