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
fn sourced_wall_density_changes_real_gas_motion_without_changing_either_inventory() {
    use crate::cavity_phs::piston::tests::real_gas::{real_controls, real_parcel};

    let initial = real_parcel(300.0, 5e6);
    let controls = real_controls();
    let input = PistonInput {
        // Resolve the dense gas spring even for the 0.1 kg specimen.
        duration_s: 5e-6,
        heat_j: 200.0,
        exterior_pressure_pa: initial.state().pressure_pa(),
        external_force_n: 0.0,
    };
    let run = |density| {
        let specimen = body(density);
        let mut wall = MaterialGasPiston::new(
            initial.clone(), specimen.clone(), 0.01, 0.0,
        ).unwrap();
        let report = wall.step(input, controls, || false).unwrap();
        assert_eq!(wall.body(), &specimen);
        (wall, report)
    };
    let (light, light_report) = run(1000.0);
    let (dense, dense_report) = run(8000.0);
    near(light.piston().mass_kg(), 0.1);
    near(dense.piston().mass_kg(), 0.8);
    assert!(light.piston().velocity_m_s() > dense.piston().velocity_m_s());
    assert!(dense.piston().velocity_m_s() > 0.0);
    assert!(light.gas().volume_m3() > dense.gas().volume_m3());
    assert!(light_report.gas.work_by_gas_j > dense_report.gas.work_by_gas_j);
    for (wall, report) in [(&light, &light_report), (&dense, &dense_report)] {
        assert_eq!(wall.piston().mass_kg(), wall.body().inertia().mass());
        assert_eq!(wall.gas().mass_kg(), initial.mass_kg());
        near(wall.gas().state().density_kg_m3(),
            initial.mass_kg() / wall.gas().volume_m3());
        assert_ne!(wall.gas().state().density_kg_m3(),
            wall.gas().state().source().state().density);
        let stored = initial.mass_kg() * initial.model()
            .specific_internal_energy_change_between_densities(
                initial.density_kg_m3(), initial.state().temperature_k(),
                wall.gas().density_kg_m3(), wall.gas().state().temperature_k(),
            ).unwrap();
        let independent_balance = stored
            + 0.5 * wall.body().inertia().mass() * wall.piston().velocity_m_s().powi(2)
            + input.exterior_pressure_pa * (wall.gas().volume_m3() - initial.volume_m3())
            - input.heat_j;
        assert!(independent_balance.abs() <= controls.energy_tolerance_j);
        assert!(report.total_energy_residual_j.abs() <= controls.energy_tolerance_j);
        assert!(report.momentum_residual_n_s.abs() <= controls.momentum_tolerance_n_s);
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

mod nonideal {
    //! Source-bound rigid inertia on the existing manufactured PR gas fixture.
    //! The native parcel fixture is synthetic; no empirical validation is claimed.
    use super::body;
    use crate::cavity_phs::piston::material::MaterialGasPiston;
    use crate::cavity_phs::piston::nonideal_tests::{controls, drive as native_drive, parcel};
    use crate::cavity_phs::piston::thermal::layers::MaterialConductionLayers;
    use crate::cavity_phs::piston::thermal::layers::tests::{
        material as thermal_material, point as thermal_point,
    };
    use crate::cavity_phs::piston::thermal::{HeatExchangeSolve, ThermalPistonError};
    use crate::cavity_phs::piston::{
        GasPiston, PistonError, PistonInput, PistonSolve, PistonSupport,
    };
    use fs_material::phase::UniformEnthalpyStepInput;
    use fs_material::state_point::MaterialPropertySelection;
    use fs_ports::{ModuleWindowAudit, PortedModule};

    fn drive(heat_j: f64) -> PistonInput {
        // Resolve the faster acoustic/mechanical clock of the lighter specimen.
        PistonInput {
            duration_s: 2e-6,
            ..native_drive(heat_j)
        }
    }

    fn same_ports(left: &impl PortedModule, right: &impl PortedModule) {
        assert_eq!(left.port_schemas(), right.port_schemas());
        assert_eq!(
            left.port_coordinate_charts(),
            right.port_coordinate_charts()
        );
        assert_eq!(left.junction_schemas(), right.junction_schemas());
        assert_eq!(
            left.last_window_duration_s(),
            right.last_window_duration_s()
        );
        assert_eq!(left.last_window_transfers(), right.last_window_transfers());
        assert_eq!(
            left.last_window_transfer_error(),
            right.last_window_transfer_error()
        );
        assert_eq!(
            left.audit_last_window().unwrap(),
            right.audit_last_window().unwrap()
        );
    }

    #[test]
    fn real_gas_binding_reuses_the_exact_native_owner_and_canonical_inertia() {
        let gas = parcel();
        let specimen = body(8000.0);
        let mut native =
            GasPiston::new(gas.clone(), specimen.inertia().mass(), 0.01, -0.03).unwrap();
        let mut sourced = MaterialGasPiston::new(gas, specimen.clone(), 0.01, -0.03).unwrap();
        same_ports(&sourced, &native);
        for heat_j in [25.0, 0.0, -5.0] {
            let expected = native.step(drive(heat_j), controls(), || false).unwrap();
            let actual = sourced.step(drive(heat_j), controls(), || false).unwrap();
            assert_eq!(actual, expected);
            assert_eq!(sourced.piston(), &native);
            assert_eq!(sourced.body(), &specimen);
            assert_eq!(sourced.piston().mass_kg(), specimen.inertia().mass());
            assert_eq!(
                sourced.gas().point().axes().get("specimen-lot"),
                Some(&42.0)
            );
            same_ports(&sourced, &native);
        }
        assert!(matches!(
            sourced.audit_last_window().unwrap(),
            ModuleWindowAudit::Available(_)
        ));
    }

    #[test]
    fn real_gas_fixed_support_and_failed_updates_preserve_the_specimen_and_ports() {
        let specimen = body(8000.0);
        let gas = parcel();
        let mut direct = gas.clone();
        let mut wall = MaterialGasPiston::new(gas.clone(), specimen.clone(), 0.01, 0.0)
            .unwrap()
            .with_fixed_support()
            .unwrap();
        let input = PistonInput {
            duration_s: 0.5,
            ..drive(200.0)
        };
        let fixed = PistonSolve {
            volume_bracket_m3: [gas.volume_m3(); 2],
            ..controls()
        };
        let step = wall.step(input, fixed, || false).unwrap();
        direct.add_heat(input.heat_j, fixed.gas, || false).unwrap();
        assert_eq!(wall.body(), &specimen);
        assert_eq!(wall.piston().support(), PistonSupport::Fixed);
        assert_eq!(wall.piston().velocity_m_s(), 0.0);
        assert_eq!(wall.gas().volume_m3(), gas.volume_m3());
        assert_eq!(step.support_work_j, 0.0);
        assert_eq!(step.gas.work_by_gas_j, 0.0);
        assert!((wall.gas().state().temperature_k() - direct.state().temperature_k()).abs() < 1e-8);
        same_ports(&wall, wall.piston());
        let before = wall.clone();
        assert!(matches!(
            wall.step(
                PistonInput {
                    heat_j: 1e12,
                    ..input
                },
                fixed,
                || false
            ),
            Err(PistonError::Gas(_))
        ));
        assert_eq!(wall, before);
        same_ports(&wall, &before);
        let mut probe = before.clone();
        let mut calls = 0;
        let expected = probe
            .step(input, fixed, || {
                calls += 1;
                false
            })
            .unwrap();
        let mut observed = 0;
        assert_eq!(
            wall.step(input, fixed, || {
                observed += 1;
                observed == calls
            }),
            Err(PistonError::Cancelled)
        );
        assert_eq!(wall, before);
        same_ports(&wall, &before);
        assert_eq!(wall.step(input, fixed, || false).unwrap(), expected);
        assert_eq!(wall, probe);
        same_ports(&wall, &probe);
    }

    #[test]
    fn real_gas_thermal_candidate_preserves_body_identity_and_complete_native_history() {
        let mut wall = MaterialGasPiston::new(parcel(), body(8000.0), 0.01, 0.03).unwrap();
        wall.step(drive(10.0), controls(), || false).unwrap();
        let before = wall.clone();
        let (card, key, caloric) = thermal_material(Some(0.2), true, false);
        let context = thermal_point(300.0);
        let curve = caloric.curve();
        let layers = || MaterialConductionLayers {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
            conductivity: &key,
            selection: MaterialPropertySelection::SingleClaimOnly,
            area_m2: 0.02,
            gas_distance_m: 1e-5,
            body_distance_m: 1e-4,
        };
        let input = drive(5.0);
        let thermal = UniformEnthalpyStepInput {
            curve,
            initial: curve.state_at_specific_enthalpy(150000.0).unwrap(),
            mass_kg: 0.001,
            volume_m3: 0.001 / 4000.0,
            surface_area_m2: 0.02,
            internal_heat_j: 0.0,
            duration_s: input.duration_s,
        };
        let exchange = HeatExchangeSolve {
            heat_bracket_j: [-1.0, 1.0],
            source_tolerance_j: 1e-9,
            energy_tolerance_j: 1e-6,
            maximum_iterations: 128,
        };
        let native = wall
            .piston()
            .thermal_candidate(
                input,
                thermal,
                controls(),
                exchange,
                |gas, body| layers().resolve(gas, body),
                || false,
            )
            .unwrap();
        let mut calls = 0;
        let candidate = wall
            .thermal_layers_candidate(input, thermal, controls(), exchange, layers(), || {
                calls += 1;
                false
            })
            .unwrap();
        assert_eq!(candidate.piston().piston(), native.piston());
        assert_eq!(candidate.piston().body(), before.body());
        assert_eq!(candidate.thermal(), native.thermal());
        assert_eq!(candidate.step(), native.step());
        assert!(candidate.step().exchanged_heat_j > 0.0);
        assert!(candidate.step().total_energy_residual_j.abs() <= exchange.energy_tolerance_j);
        same_ports(candidate.piston(), native.piston());
        assert_eq!(wall, before);
        same_ports(&wall, &before);
        let mut observed = 0;
        assert!(matches!(
            wall.thermal_layers_candidate(input, thermal, controls(), exchange, layers(), || {
                observed += 1;
                observed == calls
            }),
            Err(ThermalPistonError::Cancelled)
        ));
        assert_eq!(wall, before);
        same_ports(&wall, &before);
        let (accepted, accepted_thermal, report) = candidate.into_parts();
        assert_eq!(accepted.body(), before.body());
        assert_eq!(accepted_thermal.report.card_identity(), card.content_hash());
        assert_eq!(
            accepted_thermal.report.query_point().get("T"),
            Some(&report.thermal_midpoint.temperature_k())
        );
        assert_eq!(accepted_thermal.report.query_point().get("lot"), Some(&7.0));
        assert_eq!(accepted_thermal.external_heat_j, -report.exchanged_heat_j);
        same_ports(&accepted, native.piston());
    }
}
