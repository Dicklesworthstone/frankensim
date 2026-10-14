//! G1/G3 reference-mass loading through the existing material/Newton owners.
//! Synthetic data only: not measured gravity sag or material qualification.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point, INVERSE_TEMPERATURE_DIMS};
use fs_qty::{Density, Pressure, Temperature};
use fs_solid::rod::TipLoad;
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};
use fs_solid::rod::material_field::supports::RodSupport;
use fs_solid::rod::material_field::thermal::source::{RodLayerSource, RodSegmentSource, RodMaterialSource};

fn card(rho: f64, hot_density_scale: f64, hot_stiffness_scale: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, cold, hot, dims) in [
        ("density", rho, rho*hot_density_scale, Density::DIMS),
        ("young_modulus", 1e7, 1e7*hot_stiffness_scale, Pressure::DIMS),
        ("shear_modulus", 4e6, 4e6*hot_stiffness_scale, Pressure::DIMS),
        ("linear_thermal_expansion_coefficient", 0.0, 0.0, INVERSE_TEMPERATURE_DIMS),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, cold), (400.0, hot)], dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside, observations: vec![],
            provenance: Provenance { source: "synthetic reference-mass gravity".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic rod".into(), phase: "solid".into(),
        process: "test".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn point(t: f64) -> QueryPoint { QueryPoint::new().with("T", t).unwrap() }
fn section(card: &MaterialCard) -> ConcentricShearSection {
    let requirements = [("density", Density::DIMS), ("young_modulus", Pressure::DIMS),
        ("shear_modulus", Pressure::DIMS)].into_iter().map(|(key, dims)|
            ScalarPropertyRequirement::try_new(key, dims, ScalarAdmissibility::StrictlyPositive).unwrap())
        .collect::<Vec<_>>();
    let material = resolve_material_state_point(card, &point(300.0), &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: 0.05, material,
    }]).unwrap().with_shear_modulus("shear_modulus", 0.8).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual-expected).abs() <= tolerance*expected.abs().max(1e-9),
        "actual={actual}, expected={expected}");
}
fn source(card: &MaterialCard, n: usize) -> RodMaterialSource<'_> {
    let segments = (0..n).map(|_| RodSegmentSource { inner_radius_m: 0.0, shear_correction: 0.8,
        layers: vec![RodLayerSource { outer_radius_m: 0.05, card, reference_point: point(300.0),
            mechanical_selection: MaterialPropertySelection::SingleClaimOnly,
            expansion_selection: MaterialPropertySelection::SingleClaimOnly }] }).collect();
    RodMaterialSource::try_new(1.0, segments, "shear_modulus", 0.02, &CancelGate::new()).unwrap()
}

#[test]
fn g1_spatial_mass_and_first_moment_include_the_clamped_root() {
    let a = section(&card(1000.0, 1.0, 1.0));
    let b = section(&card(2000.0, 1.0, 1.0));
    let rod = SpatialMaterialRod::try_straight(1.0, vec![a.clone(), b.clone()]).unwrap();
    let mass = rod.reference_mass();
    let m0 = 0.5*a.elastic().line_density_kg_m();
    let m1 = 0.5*b.elastic().line_density_kg_m();
    assert_eq!(mass.segment_mass_kg(), &[m0, m1]);
    assert_eq!(mass.nodal_mass_kg(), &[0.5*m0, 0.5*(m0+m1), 0.5*m1]);
    assert_eq!(mass.sections(), &[a.clone(), b]);
    near(mass.total_mass_kg(), m0+m1, 2e-15);
    let nodal_first_moment: f64 = mass.nodal_mass_kg().iter().enumerate().map(|(i,m)| 0.5*i as f64*m).sum();
    near(nodal_first_moment, 0.25*m0+0.75*m1, 2e-15);
    let split = SpatialMaterialRod::try_straight(1.0, vec![a.clone(); 4]).unwrap();
    let whole = SpatialMaterialRod::try_straight(1.0, vec![a]).unwrap();
    near(split.reference_mass().total_mass_kg(), whole.reference_mass().total_mass_kg(), 2e-15);
}

#[test]
fn g1_axial_body_force_drives_series_equilibrium_and_root_reaction() {
    let gate = CancelGate::new();
    let mut rod = SpatialMaterialRod::try_straight(1.0, vec![
        section(&card(1000.0, 1.0, 1.0)), section(&card(2000.0, 1.0, 1.0)),
    ]).unwrap().with_gravity([0.1, 0.0, 0.0], &gate).unwrap();
    let weights = rod.gravity_nodal_forces_n().unwrap().to_vec();
    let tip = TipLoad { force: [0.2, 0.0, 0.0], moment: [0.0; 3] };
    rod.solve_static(&tip, 2, 1e-7, &gate).unwrap();
    let mut extension = 0.0;
    for segment in 0..2 {
        let expected_force = weights[segment+1..].iter().map(|f| f[0]).sum::<f64>() + tip.force[0];
        near(rod.segment_resultants(segment).unwrap().0[0], expected_force, 2e-5);
        extension += 0.5*expected_force/rod.sections()[segment].elastic().axial_stiffness_n();
    }
    near(rod.positions_m()[2][0]-1.0, extension, 2e-5);
    let root = rod.root_support_reaction(&gate).unwrap();
    near(root.force_n[0].unwrap(), -rod.reference_mass().total_mass_kg()*0.1-tip.force[0], 2e-5);
    assert!(rod.support_reactions(&tip, &gate).unwrap().free_residual_norm < 1e-7);
    assert!(rod.energy_j().unwrap() > 0.0);
    assert!(rod.gravity_potential_j().unwrap() < 0.0);
}

#[test]
fn g1_g3_density_only_swap_changes_actual_cantilever_sag() {
    let gate = CancelGate::new();
    let mut displacements = Vec::new();
    for rho in [1000.0, 2000.0] {
        let s = section(&card(rho, 1.0, 1.0));
        let q = 0.1*s.elastic().line_density_kg_m();
        let expected = q*(1.0/(8.0*s.elastic().bending_stiffness_n_m2()) + 1.0/(2.0*s.shear_stiffness_n()));
        let mut rod = SpatialMaterialRod::try_straight(1.0, vec![s; 4]).unwrap()
            .with_gravity([0.0, 0.0, -0.1], &gate).unwrap();
        rod.solve_static(&TipLoad::default(), 1, 1e-7, &gate).unwrap();
        let sag = -rod.positions_m()[4][2];
        near(sag, expected, 0.02); // small-load Timoshenko limit; not a stability claim
        displacements.push(sag);
        let root = rod.root_support_reaction(&gate).unwrap();
        near(root.force_n[2].unwrap(), q, 2e-5);
        // Moment balances weight at the accepted deformed nodal positions.
        let moment = rod.positions_m().iter().zip(rod.gravity_nodal_forces_n().unwrap())
            .map(|(x,f)| x[2]*f[0]-x[0]*f[2]).sum::<f64>();
        near(root.moment_body_nm.unwrap()[1], -moment, 2e-5);
    }
    near(displacements[1]/displacements[0], 2.0, 0.002);
}

#[test]
fn g3_heating_changes_stiffness_but_not_reference_weight_or_receipts() {
    let gate = CancelGate::new();
    // Deliberately strong density variation catches accidental rho(T)*V0 mass.
    let card = card(1000.0, 0.5, 0.5);
    let source = source(&card, 2);
    let mut cold = source.resolve_at(&[300.0; 2], &gate).unwrap()
        .with_gravity([0.1, 0.0, 0.0], &gate).unwrap();
    let mut hot = source.resolve_at(&[400.0; 2], &gate).unwrap()
        .with_gravity([0.1, 0.0, 0.0], &gate).unwrap();
    assert_eq!(cold.reference_mass(), hot.reference_mass());
    assert!(std::ptr::eq(cold.reference_mass(), hot.reference_mass()));
    assert_eq!(cold.gravity_nodal_forces_n(), hot.gravity_nodal_forces_n());
    assert_ne!(cold.sections(), hot.sections());
    near(hot.sections()[0].elastic().line_density_kg_m(),
        0.5*cold.sections()[0].elastic().line_density_kg_m(), 2e-15);
    for rod in [&mut cold, &mut hot] {
        rod.solve_static(&TipLoad::default(), 1, 1e-7, &gate).unwrap();
    }
    near(hot.positions_m()[2][0]-1.0, 2.0*(cold.positions_m()[2][0]-1.0), 2e-5);
    near(hot.reference_mass().total_mass_kg(), source.reference_mass_kg(), 2e-15);
}

#[test]
fn g1_fully_constrained_gravity_reactions_include_every_mass_and_no_fake_energy() {
    let gate = CancelGate::new();
    let rod = SpatialMaterialRod::try_straight(1.0, vec![section(&card(1000.0, 1.0, 1.0)); 2]).unwrap();
    let supports = (1..=2).map(|node| RodSupport { node,
        position_m: [Some(0.5*node as f64), Some(0.0), Some(0.0)], frame: Some(So3::identity()),
    }).collect::<Vec<_>>();
    let mut rod = rod.with_supports(&supports, &gate).unwrap().with_gravity([0.0, 0.0, -1.0], &gate).unwrap();
    assert_eq!(rod.solve_static(&TipLoad::default(), 1, 1e-7, &gate).unwrap(), vec![0.0]);
    let mut reactions = rod.support_reactions(&TipLoad::default(), &gate).unwrap().nodes;
    reactions.insert(0, rod.root_support_reaction(&gate).unwrap());
    for (reaction, force) in reactions.iter().zip(rod.gravity_nodal_forces_n().unwrap()) {
        near(reaction.force_n[2].unwrap(), -force[2], 2e-12);
    }
    assert_eq!(rod.energy_j().unwrap(), 0.0);
    assert_eq!(rod.gravity_potential_j().unwrap(), 0.0);
}

#[test]
fn g0_g4_zero_gravity_legacy_equivalence_and_failures_leave_pose_intact() {
    let gate = CancelGate::new();
    let mut original = SpatialMaterialRod::try_straight(1.0, vec![section(&card(1000.0, 1.0, 1.0))]).unwrap();
    let mut zero = original.with_gravity([0.0; 3], &gate).unwrap();
    let load = TipLoad { force: [0.1, 0.0, 0.0], moment: [0.0; 3] };
    original.solve_static(&load, 1, 1e-7, &gate).unwrap();
    zero.solve_static(&load, 1, 1e-7, &gate).unwrap();
    assert_eq!(zero.positions_m(), original.positions_m());
    assert_eq!(zero.frames(), original.frames());
    assert_eq!(zero.energy_j().unwrap().to_bits(), original.energy_j().unwrap().to_bits());
    let pose = original.positions_m().to_vec();
    for g in [[f64::NAN, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0], [f64::MAX, 0.0, 0.0]] {
        assert!(original.with_gravity(g, &gate).is_err());
    }
    let weighted = original.with_gravity([0.0, 0.0, -1.0], &gate).unwrap();
    let removed = weighted.with_gravity([0.0; 3], &gate).unwrap();
    assert!(removed.gravity_nodal_forces_n().is_none());
    assert_eq!(removed.reference_mass(), original.reference_mass());
    gate.request();
    assert!(matches!(original.with_gravity([0.0, 0.0, -1.0], &gate), Err(SpatialRodError::Cancelled)));
    assert!(matches!(original.solve_static(&load, 1, 1e-7, &gate), Err(SpatialRodError::Cancelled)));
    assert_eq!(original.positions_m(), pose);
}

#[test]
fn g1_small_thermal_member_probes_do_not_exceed_the_material_strain_limit() {
    let card = card(1000.0, 0.5, 0.5);
    let gate = CancelGate::new();
    let original = source(&card, 4);
    let short = RodMaterialSource::try_new(0.04, original.segments().to_vec(),
        "shear_modulus", 0.03, &gate).unwrap();
    let mut rod = short.resolve_at(&[350.0; 4], &gate).unwrap()
        .with_gravity([10.0, 0.0, 0.0], &gate).unwrap();
    // The former 0.0003 m tangent plus 0.000003 m residual probe created
    // 0.0303 strain on these 0.01 m segments and refused before Newton.
    rod.solve_static(&TipLoad::default(), 1, 1e-7, &gate).unwrap();
    let expected = 1000.0*10.0*0.04_f64.powi(2)/(2.0*7.5e6);
    near(rod.positions_m()[4][0]-0.04, expected, 2e-5);
    assert!(rod.support_reactions(&TipLoad::default(), &gate).unwrap().free_residual_norm < 1e-7);
    near(rod.root_support_reaction(&gate).unwrap().force_n[0].unwrap(),
        -10.0*rod.reference_mass().total_mass_kg(), 2e-5);
}
