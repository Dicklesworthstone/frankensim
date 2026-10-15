//! Synthetic material coefficients test section physics, not measured alloys.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, DynViscosity, Pressure};

fn material(rho: f64, shear: f64, eta: f64, band: (f64, f64)) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new(); let mut requirements = Vec::new();
    for (name, dims, value) in [("density", Density::DIMS, rho),
        ("young_modulus", Pressure::DIMS, 200e9), ("G", Pressure::DIMS, shear),
        ("eta_s", DynViscosity::DIMS, eta)] {
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", band.0, band.1),
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source: "synthetic torsion".into(), license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "test".into(), phase: "solid".into(),
        process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap();
    resolve_material_state_point(&card, &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn section(a: ResolvedMaterialStatePoint, b: ResolvedMaterialStatePoint) -> ConcentricElasticSection {
    ConcentricElasticSection::try_new(0.0, vec![
        ConcentricMaterialLayer { outer_radius_m: 0.0002, material: a },
        ConcentricMaterialLayer { outer_radius_m: 0.0003, material: b },
    ]).unwrap()
}
fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-14 * b.abs().max(1e-30), "{a:e} != {b:e}"); }

#[test]
fn original_layer_integrals_control_rotary_mass_stiffness_and_distinct_shear_loss() {
    let a = material(7800.0, 80e9, 11.0, (0.0, 50_000.0));
    let b = material(15_600.0, 40e9, 29.0, (100.0, 40_000.0));
    let elastic = section(a, b);
    let j: Vec<_> = elastic.layer_integrals().iter().map(|g| 2.0*g.second_moment_m4()).collect();
    let shear = elastic.clone().with_shear_modulus("G", 0.7).unwrap();
    let gj = shear.torsional_stiffness_n_m2();
    let torsion = shear.with_torsional_viscosity("eta_s").unwrap();
    assert_eq!(torsion.shear().elastic(), &elastic);
    assert_eq!(torsion.stiffness_n_m2(), gj);
    near(torsion.rotary_inertia_kg_m(), 7800.0*j[0] + 15_600.0*j[1]);
    near(torsion.viscous_stiffness_n_m2_s(), 11.0*j[0] + 29.0*j[1]);
    assert_eq!(torsion.omega_band_rad_s(), (100.0, 40_000.0));
    assert_eq!(torsion.viscosity_property(), "eta_s");
}

#[test]
fn artificial_layer_splitting_and_shear_correction_do_not_change_torsion() {
    let m = material(7800.0, 80e9, 10.0, (0.0, 50_000.0));
    let whole = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: 0.0003, material: m.clone(),
    }]).unwrap().with_shear_modulus("G", 1.0).unwrap().with_torsional_viscosity("eta_s").unwrap();
    let split = section(m.clone(), m).with_shear_modulus("G", 0.5).unwrap()
        .with_torsional_viscosity("eta_s").unwrap();
    near(whole.rotary_inertia_kg_m(), split.rotary_inertia_kg_m());
    near(whole.stiffness_n_m2(), split.stiffness_n_m2());
    near(whole.viscous_stiffness_n_m2_s(), split.viscous_stiffness_n_m2_s());
    near(whole.stiffness_n_m2()/whole.rotary_inertia_kg_m(), 80e9/7800.0);
}

#[test]
fn density_and_shear_only_swaps_do_not_reuse_old_rotary_coefficients() {
    let make = |rho, shear| { let m=material(rho,shear,10.0,(0.0,50_000.0));
        section(m.clone(),m).with_shear_modulus("G",1.0).unwrap().with_torsional_viscosity("eta_s").unwrap() };
    let a=make(7800.0,80e9); let dense=make(15_600.0,80e9); let stiff=make(7800.0,160e9);
    assert_eq!(a.stiffness_n_m2(),dense.stiffness_n_m2());
    near(dense.rotary_inertia_kg_m(),2.0*a.rotary_inertia_kg_m());
    near(stiff.stiffness_n_m2(),2.0*a.stiffness_n_m2());
    assert_eq!(stiff.rotary_inertia_kg_m(),a.rotary_inertia_kg_m());
    assert_eq!(stiff.viscous_stiffness_n_m2_s(),a.viscous_stiffness_n_m2_s());
}

#[test]
fn zero_shear_loss_is_explicit_and_missing_wrong_or_negative_loss_refuses() {
    for eta in [0.0, -1.0] {
        let m=material(7800.0,80e9,eta,(0.0,50_000.0));
        let s=section(m.clone(),m).with_shear_modulus("G",1.0).unwrap();
        assert!(s.clone().with_torsional_viscosity("missing").is_err());
        assert!(s.clone().with_torsional_viscosity("G").is_err());
        let result=s.with_torsional_viscosity("eta_s");
        if eta==0.0 { assert_eq!(result.unwrap().viscous_stiffness_n_m2_s(),0.0); }
        else { assert!(result.is_err()); }
    }
}
