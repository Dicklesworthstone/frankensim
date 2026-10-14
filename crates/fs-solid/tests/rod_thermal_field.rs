//! G1/G3 sourced temperature -> local constitutive response -> shared equilibrium.
//! Synthetic coefficients exercise wiring and equations, not material calibration.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, integrate_isotropic_thermal_expansion,
    resolve_material_state_point};
use fs_qty::{Density, Pressure, Temperature};
use fs_solid::rod::TipLoad;
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};

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
