//! G1/G3 material-field statics against series compliance and beam equations.
//! All property values below are synthetic, not measured material qualifications.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Pressure};
use fs_solid::rod::{Rod, RodSection, TipLoad};
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};

#[path = "rod_material_field/reference.rs"]
mod reference;

fn section(scale: f64) -> ConcentricShearSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [("density", 1000.0, Density::DIMS),
        ("young_modulus", 1e6 * scale, Pressure::DIMS),
        ("shear_modulus", 4e5 * scale, Pressure::DIMS)]
    {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(), uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "synthetic spatial rod".into(), license: "test-only".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic rod".into(),
        phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    let material = resolve_material_state_point(&card, &QueryPoint::new(), &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    ConcentricElasticSection::try_new(0.0,
        vec![ConcentricMaterialLayer { outer_radius_m: 0.05, material }])
        .unwrap().with_shear_modulus("shear_modulus", 0.8).unwrap()
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!((actual - expected).abs() <= relative * expected.abs().max(1e-12),
        "actual={actual}, expected={expected}");
}

#[test]
fn uniform_field_uses_the_identical_equilibrium_path() {
    let source = section(1.0);
    let coefficients = RodSection::try_from_concentric_material(&source).unwrap();
    let mut uniform = Rod::straight(1.0, 2, coefficients);
    let mut spatial = SpatialMaterialRod::try_straight(1.0, vec![source.clone(); 2]).unwrap();
    let load = TipLoad { force: [0.5, 0.001, 0.0], moment: [0.002, 0.0, 0.0005] };
    let history = uniform.solve_static(&load, 2, 1e-7).unwrap();
    let actual = spatial.solve_static(&load, 2, 1e-7, &CancelGate::new()).unwrap();
    assert_eq!(actual, history);
    assert_eq!(spatial.positions_m(), uniform.positions.as_slice());
    assert_eq!(spatial.frames(), uniform.frames.as_slice());
    assert_eq!(spatial.energy_j().unwrap().to_bits(), uniform.energy().unwrap().to_bits());
    assert_eq!(spatial.sections(), &[source.clone(), source]);
}

#[test]
fn each_segment_contributes_its_actual_axial_and_torsional_compliance() {
    let sections = vec![section(1.0), section(2.0), section(0.5)];
    let mut rod = SpatialMaterialRod::try_straight(1.2, sections.clone()).unwrap();
    let load = TipLoad { force: [0.5, 0.0, 0.0], moment: [0.01, 0.0, 0.0] };
    rod.solve_static(&load, 1, 1e-7, &CancelGate::new()).unwrap();
    let mut extension = 0.0;
    let mut angle = 0.0;
    let mut expected_energy = 0.0;
    for (i, source) in sections.iter().enumerate() {
        let (gamma, kappa) = rod.segment_strains(i).unwrap();
        let ea = source.elastic().axial_stiffness_n();
        let gj = source.torsional_stiffness_n_m2();
        near(gamma[0], 0.5 / ea, 3e-5);
        near(kappa[0], 0.01 / gj, 3e-5);
        let (force, moment) = rod.segment_resultants(i).unwrap();
        near(force[0], 0.5, 3e-5);
        near(moment[0], 0.01, 3e-5);
        extension += 0.4 * 0.5 / ea;
        angle += 0.4 * kappa[0];
        expected_energy += 0.5 * 0.4 * (0.5_f64.powi(2) / ea + 0.01_f64.powi(2) / gj);
    }
    near(rod.positions_m()[3][0] - 1.2, extension, 3e-5);
    let exact_angle: f64 = sections.iter().map(|s| 0.4 * 0.01 / s.torsional_stiffness_n_m2()).sum();
    near(angle, exact_angle, 3e-5);
    near(rod.energy_j().unwrap(), expected_energy, 3e-5);
    assert_eq!(rod.sections(), sections);
}

#[test]
fn moving_the_soft_material_changes_bending_without_changing_average_properties() {
    let soft = section(0.5);
    let stiff = section(2.0);
    let distributions = [vec![soft.clone(), stiff.clone(), stiff.clone(), stiff.clone()],
        vec![stiff.clone(), stiff.clone(), stiff, soft]];
    let mut tips = Vec::new();
    for sources in distributions {
        // The midpoint Cosserat discretization's small-deflection limit:
        // delta_y = P sum L_i ((L-x_mid_i)^2/EI_i + 1/GA_i).
        let expected: f64 = sources.iter().enumerate().map(|(i, s)| {
            let lever = 1.0 - (i as f64 + 0.5) * 0.25;
            0.01 * 0.25 * (lever * lever / s.elastic().bending_stiffness_n_m2()
                + 1.0 / s.shear_stiffness_n())
        }).sum();
        let mut rod = SpatialMaterialRod::try_straight(1.0, sources).unwrap();
        rod.solve_static(&TipLoad { force: [0.0, 0.01, 0.0], moment: [0.0; 3] },
            1, 1e-7, &CancelGate::new()).unwrap();
        let tip = rod.positions_m()[4][1];
        near(tip, expected, 0.005);
        tips.push(tip);
    }
    assert!(tips[0] > 1.5 * tips[1], "a modulus average would lose this spatial effect");
}

#[test]
fn invalid_inputs_and_cancellation_preserve_the_last_accepted_pose() {
    let source = section(1.0);
    assert!(SpatialMaterialRod::try_straight(1.0, vec![]).is_err());
    assert!(SpatialMaterialRod::try_straight(1.0, vec![source.clone(); 65]).is_err());
    assert!(SpatialMaterialRod::try_straight(f64::NAN, vec![source.clone()]).is_err());
    let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source.clone(); 2]).unwrap();
    rod.solve_static(&TipLoad { force: [0.5, 0.0, 0.0], moment: [0.0; 3] },
        1, 1e-7, &CancelGate::new()).unwrap();
    let before = rod.positions_m().to_vec();
    let energy = rod.energy_j().unwrap();
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(rod.solve_static(&TipLoad::default(), 1, 1e-7, &gate), Err(SpatialRodError::Cancelled));
    for (steps, tolerance) in [(0, 1e-7), (257, 1e-7), (1, 0.0), (1, f64::NAN)] {
        assert!(rod.solve_static(&TipLoad::default(), steps, tolerance, &CancelGate::new()).is_err());
    }
    assert!(rod.segment_strains(usize::MAX).is_err());
    assert!(rod.segment_resultants(2).is_err());
    assert_eq!(rod.positions_m(), before);
    assert_eq!(rod.energy_j().unwrap(), energy);
    assert_eq!(rod.sections(), &[source.clone(), source]);
}
