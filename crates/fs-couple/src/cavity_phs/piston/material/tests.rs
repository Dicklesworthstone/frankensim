//! Density/geometry source changes exercise the actual coupled piston solver.
use super::*;
use crate::cavity_phs::piston::tests::{controls, drive, gas};
use fs_evidence::ValidityDomain;
use fs_ga::Vec3;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd_material::{RegionBinding, VolumeMoments, resolve_body};
use fs_qty::{Density, Temperature};

fn card(density: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("density", Density::DIMS),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0,density),(400.0,density)], dims: Density::DIMS },
        validity: ValidityDomain::unconstrained().with("T",200.0,400.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic specimen density".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic specimen".into(),
        phase: "solid".into(), process: "authored".into(), revision: 0 }, claims, vec![]).unwrap()
}
pub(in crate::cavity_phs::piston) fn body(density: f64) -> ResolvedBody {
    let card = card(density);
    let point = QueryPoint::new().with("T",300.0).unwrap();
    resolve_body(&[RegionBinding { key: "wall", material: &card, point: &point,
        geometry: VolumeMoments::cuboid(Vec3::new(0.1,0.1,0.01)).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly }],1,||false).unwrap()
}
fn near(a: f64, b: f64) { assert!((a-b).abs() <= 2e-12*a.abs().max(b.abs()).max(1.0)); }

#[test]
fn density_card_swap_changes_inertia_and_the_coupled_thermal_trajectory() {
    let run = |density| {
        let mut wall = MaterialGasPiston::new(gas(0.75,0.001),body(density),0.01,0.0).unwrap();
        let report = wall.step(drive(20.0),controls(),||false).unwrap(); (wall,report)
    };
    let (light,lr) = run(1000.0); let (dense,dr) = run(8000.0);
    near(light.piston().mass_kg(),0.1); near(dense.piston().mass_kg(),0.8);
    assert!(light.piston().velocity_m_s() > dense.piston().velocity_m_s());
    assert!(light.gas().volume_m3() > dense.gas().volume_m3());
    // More boundary work leaves less energy in the same gas after the same heat.
    assert!(lr.gas.work_by_gas_j > dr.gas.work_by_gas_j);
    assert!(light.gas().state().state().temperature < dense.gas().state().state().temperature);
    assert!(light.gas().state().state().sound_speed < dense.gas().state().state().sound_speed);
    for (wall,report) in [(&light,&lr),(&dense,&dr)] {
        assert_eq!(wall.piston().mass_kg(),wall.body().inertia().mass());
        assert_eq!(wall.body().regions().len(),1);
        assert_eq!(wall.gas().mass_kg(),gas(0.75,0.001).mass_kg());
        assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
    }
}

#[test]
fn composite_regions_use_the_shared_inertia_and_identical_momentum_solver() {
    let a=card(1000.0); let b=card(8000.0);
    let point=QueryPoint::new().with("T",300.0).unwrap();
    // Adjacent non-overlapping tetrahedral regions in one constrained body.
    let geometry = |z0| VolumeMoments::tetrahedron([
        Vec3::new(0.0,0.0,z0), Vec3::new(0.1,0.0,z0),
        Vec3::new(0.0,0.1,z0), Vec3::new(0.0,0.0,z0+0.03),
    ]).unwrap();
    let regions=[
        RegionBinding { key: "a", material: &a, point: &point, geometry: geometry(0.0),
            selection: MaterialPropertySelection::SingleClaimOnly },
        RegionBinding { key: "b", material: &b, point: &point, geometry: geometry(0.03),
            selection: MaterialPropertySelection::SingleClaimOnly },
    ];
    let body=resolve_body(&regions,2,||false).unwrap();
    let reversed=resolve_body(&[regions[1].clone(),regions[0].clone()],2,||false).unwrap();
    assert_eq!(body,reversed);
    near(body.inertia().mass(),0.00005*9000.0);
    let mut numeric=GasPiston::new(gas(0.75,0.001),body.inertia().mass(),0.01,-0.1).unwrap();
    let mut sourced=MaterialGasPiston::new(gas(0.75,0.001),body.clone(),0.01,-0.1).unwrap();
    let expected=numeric.step(drive(20.0),controls(),||false).unwrap();
    let actual=sourced.step(drive(20.0),controls(),||false).unwrap();
    assert_eq!(actual,expected);
    assert_eq!(sourced.piston(),&numeric);
    assert_eq!(sourced.body(),&body);
    for (region,density) in sourced.body().regions().iter().zip([1000.0,8000.0]) {
        near(region.material().property("density").unwrap().value_si(),density);
    }
}

#[test]
fn material_bound_state_and_sources_survive_failures_and_cancellation() {
    let mut wall=MaterialGasPiston::new(gas(0.75,0.001),body(8000.0),0.01,0.0).unwrap();
    let before=wall.clone();
    assert!(matches!(wall.step(drive(20.0),controls(),||true),Err(PistonError::Cancelled)));
    assert_eq!(wall,before);
    assert!(wall.step(drive(1e9),controls(),||false).is_err());
    assert_eq!(wall,before);
    let mut polls=0; let mut probe=wall.clone();
    probe.step(drive(20.0),controls(),||{polls+=1;false}).unwrap();
    let mut count=0;
    assert!(matches!(wall.step(drive(20.0),controls(),||{count+=1;count==polls}),Err(PistonError::Cancelled)));
    assert_eq!(wall,before);
    wall.step(drive(20.0),controls(),||false).unwrap();
    // Gas heating may not silently re-evaluate or heat the rigid wall source.
    assert_eq!(wall.body(),before.body());
    assert_ne!(wall.gas().state(),before.gas().state());
}

mod double_acting {
    //! Reuse the original source fixture and inertia owner.
    use super::{body, controls, drive, gas, near};
    use crate::cavity_phs::piston::double_acting::{
        DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve,
    };
    use crate::cavity_phs::piston::material::{MaterialDoubleActingGasPiston, MaterialGasPiston};
    use crate::cavity_phs::piston::{GasPiston, PistonError};

    fn input() -> DoublePistonInput {
        DoublePistonInput { duration_s: drive(0.0).duration_s, heat_j: [20.0, 0.0], external_force_n: 0.0 }
    }
    fn solve() -> DoublePistonSolve {
        DoublePistonSolve { piston: controls(), opposing_gas: controls().gas, volume_tolerance_m3: 1e-18 }
    }

    #[test]
    fn solid_density_data_changes_one_inertia_and_both_material_gases() {
        let run = |density| {
            let mut wall = MaterialDoubleActingGasPiston::new(
                gas(0.75, 0.001), gas(0.25, 0.001), body(density), [0.01, 0.01], 0.0,
            ).unwrap();
            let report = wall.step(input(), solve(), || false).unwrap();
            (wall, report)
        };
        let (light, lr) = run(1000.0);
        let (dense, dr) = run(8000.0);
        near(light.piston().primary().mass_kg(), 0.1);
        near(dense.piston().primary().mass_kg(), 0.8);
        assert!(light.piston().primary().velocity_m_s() > dense.piston().primary().velocity_m_s());
        assert!(light.gases()[0].state().state().temperature < dense.gases()[0].state().state().temperature);
        assert!(light.gases()[1].state().state().temperature > dense.gases()[1].state().state().temperature);
        assert!(light.gases()[1].state().state().sound_speed > dense.gases()[1].state().state().sound_speed);
        for (wall, report) in [(&light, &lr), (&dense, &dr)] {
            assert_eq!(wall.piston().primary().mass_kg(), wall.body().inertia().mass());
            assert_eq!(wall.body().regions().len(), 1);
            assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
            for material in report.midpoint.iter() { assert_eq!(material.components().len(), 2); }
        }
    }

    #[test]
    fn material_binding_matches_the_original_numeric_dynamics_exactly() {
        let specimen = body(8000.0);
        let primary = gas(0.75, 0.001);
        let opposing = gas(0.25, 0.001);
        let mut numeric = DoubleActingGasPiston::new(
            GasPiston::new(primary.clone(), specimen.inertia().mass(), 0.01, -0.1).unwrap(),
            opposing.clone(), 0.006,
        ).unwrap();
        let mut sourced = MaterialDoubleActingGasPiston::new(
            primary, opposing, specimen.clone(), [0.01, 0.006], -0.1,
        ).unwrap();
        let expected = numeric.step(input(), solve(), || false).unwrap();
        let actual = sourced.step(input(), solve(), || false).unwrap();
        assert_eq!(actual, expected);
        assert_eq!(sourced.piston(), &numeric);
        assert_eq!(sourced.body(), &specimen);
    }

    #[test]
    fn adding_opposing_gas_preserves_the_existing_source_and_momentum_history() {
        let mut single = MaterialGasPiston::new(gas(0.75, 0.001), body(8000.0), 0.01, 0.0).unwrap();
        single.step(drive(5.0), controls(), || false).unwrap();
        let specimen = single.body().clone();
        let accepted = single.piston().clone();
        let opposing = gas(0.25, 0.001);
        let bound = single.with_opposing_gas(opposing.clone(), 0.01).unwrap();
        assert_eq!(bound.body(), &specimen);
        assert_eq!(bound.piston().primary(), &accepted);
        assert_eq!(bound.piston().opposing(), &opposing);
    }

    #[test]
    fn source_failure_and_final_cancellation_cannot_publish_half_a_material_system() {
        let mut wall = MaterialDoubleActingGasPiston::new(
            gas(0.75, 0.001), gas(0.25, 0.001), body(8000.0), [0.01, 0.01], 0.0,
        ).unwrap();
        let before = wall.clone();
        assert!(wall.step(DoublePistonInput { heat_j: [0.0, 1e9], ..input() }, solve(), || false).is_err());
        assert_eq!(wall, before);
        let mut probe = wall.clone();
        let mut polls = 0;
        probe.step(input(), solve(), || { polls += 1; false }).unwrap();
        let mut actual = 0;
        assert!(matches!(wall.step(input(), solve(), || {
            actual += 1; actual == polls
        }), Err(PistonError::Cancelled)));
        assert_eq!(wall, before);
        wall.step(input(), solve(), || false).unwrap();
        assert_eq!(wall.body(), before.body());
        assert_ne!(wall.gases()[0].state(), before.gases()[0].state());
        assert_ne!(wall.gases()[1].state(), before.gases()[1].state());
    }
}
