//! G1/G3/G4: one finite gas couples source-backed solids and moving liquid.
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceBodyGasFace;

#[path = "solid/radiation.rs"]
mod radiation;

fn solid_face() -> FreeSurfaceBodyGasFace {
    FreeSurfaceBodyGasFace {
        wall: [-1, 0, 3],
        gas: [0, 0, 3],
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }
}

#[test]
fn g1_one_gas_closes_solid_liquid_and_both_film_exchanges() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let mut gas = parcel(350.0, 1e-7, 100.0);
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let mut bodies =
            [body::wall(&solid, &chart, &context, 450.0).with_gas_faces(vec![solid_face()])];
        let before_liquid = thermal.total_enthalpy_j().unwrap();
        let before_solid = bodies[0].enthalpy_j();
        let before_gas = gas.state().specific_internal_energy_j_kg();
        let report = thermal
            .step_with_finite_gas_environment(
                &[],
                &[],
                controls(),
                &mut bodies,
                &[body::contact()],
                &[],
                &mut gas,
                &[exposed],
                gas_controls(),
                false,
                gate,
            )
            .unwrap()
            .conduction;

        // Independent explicit three-store calculation: two liquid cells,
        // one finite wall, and one NASA gas. All three rates use the SAME old
        // endpoint, before any of these four temperatures changes.
        let (mut twet, mut tliquid, mut ts, mut tg) = (400.0, 400.0, 450.0, 350.0);
        let (mut gas_liquid, mut gas_solid, mut solid_liquid) = (0.0, 0.0, 0.0);
        let r = 8.314_462_618_153_24 / 0.028;
        let mg = 100000.0 * 1e-7 / (r * 350.0);
        for _ in 0..4 {
            let kg = 2.6 * (0.5 + (100000.0 * tg / 350.0 - 1000.0) / 249000.0);
            let ks = 0.5 * (ts - 200.0);
            let gsl = 1e-6 / (0.0005 / (twet - 200.0) + 0.0005 / ks);
            let ggl = 1e-6 / (0.0005 / (tliquid - 200.0) + 0.0005 / kg);
            let ggs = 1e-6 / (0.0005 / ks + 0.0005 / kg);
            let qsl = 0.00025 * gsl * (ts - twet);
            let qgl = 0.00025 * ggl * (tg - tliquid);
            let qgs = 0.00025 * ggs * (tg - ts);
            twet += qsl / 0.002;
            tliquid += qgl / 0.002;
            ts += (qgs - qsl) / 0.001;
            tg = (((2.0_f64 + 0.001 * tg).powi(2) - 0.002 * (qgl + qgs) / (mg * r)).sqrt() - 2.0)
                / 0.001;
            gas_liquid += qgl;
            gas_solid += qgs;
            solid_liquid += qsl;
        }
        let films = report.body_environment.as_ref().unwrap();
        let walls = report.bodies.as_ref().unwrap();
        let environment = report.environment.as_ref().unwrap();
        let finite = report.finite_gas.as_ref().unwrap();
        assert_eq!(report.thermal.substeps, 4);
        assert_eq!(finite.substeps, 4);
        near(gas.state().state().temperature, tg, 2e-10);
        near(films.thermal_endpoint[0].solid_temperature_k, ts, 2e-11);
        near(
            environment.thermal_endpoint[0].liquid_temperature_k,
            tliquid,
            2e-11,
        );
        near(films.heat_to_solids_j, gas_solid, 2e-11);
        near(environment.heat_to_liquid_j, gas_liquid, 2e-11);
        near(walls.heat_to_liquid_j, solid_liquid, 2e-11);
        near(walls.heat_from_gas_j, gas_solid, 2e-11);
        near(report.flow.external_heat_j, gas_liquid, 2e-11);
        near(report.flow.body_heat_j, solid_liquid, 2e-11);
        near(
            bodies[0].enthalpy_j() - before_solid,
            gas_solid - solid_liquid,
            2e-11,
        );
        near(
            thermal.total_enthalpy_j().unwrap() - before_liquid,
            gas_liquid + solid_liquid,
            2e-10,
        );
        let gas_change = gas.mass_kg() * (gas.state().specific_internal_energy_j_kg() - before_gas);
        near(gas_change + gas_liquid + gas_solid, 0.0, 2e-10);
        near(
            thermal.total_enthalpy_j().unwrap() - before_liquid + bodies[0].enthalpy_j()
                - before_solid
                + gas_change,
            0.0,
            2e-9,
        );
        assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
        assert_eq!(
            films.thermal_endpoint[0].solid.card_identity(),
            solid.content_hash()
        );
        near(
            films.thermal_endpoint[0].gas.state().pressure,
            gas.state().state().pressure,
            2e-13,
        );
        assert_eq!(thermal.accepted_steps(), 1);
        near(thermal.physical_time_s(), 0.001, 1e-13);
    });
}

#[test]
fn g3_solid_material_and_finite_gas_capacity_reach_native_liquid_motion() {
    let run = |scale: f64, volume: f64| {
        let liquid = conducting_card(1.0);
        let solid = body::wall_card(scale, 500.0);
        let chart = body::wall_chart(&solid);
        let context = point(450.0);
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let mut gas = parcel(350.0, volume, 100.0);
            let exposed = gas_contact(gas.model(), gas.point()).face;
            let mut bodies =
                [body::wall(&solid, &chart, &context, 450.0).with_gas_faces(vec![solid_face()])];
            // No direct wall/liquid contact: a solid-data change reaches the
            // liquid only through the common finite atmosphere.
            let report = thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[],
                    &[],
                    &mut gas,
                    &[exposed],
                    gas_controls(),
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal
                    .material_at(EXPOSED)
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                gas.state().state().temperature,
                report.thermal.substeps,
                report.finite_gas.unwrap().combined_thermal_balance_defect_j,
            )
        })
    };
    let insulating = run(0.0001, 1e-7);
    let conducting = run(1.0, 1e-7);
    assert_ne!(insulating.0, conducting.0);
    assert_ne!(insulating.1.to_bits(), conducting.1.to_bits());
    assert!(conducting.2 > insulating.2);
    let small = run(1.0, 1e-9);
    assert!(
        small.3 > 4,
        "aggregate solid+liquid gas conductance must tighten CFL"
    );
    assert!(small.2 >= 350.0 && small.2 <= 450.0);
    assert!(small.4.abs() < 1e-9);
}

#[test]
fn g4_solid_film_domains_geometry_and_late_refusals_preserve_all_histories() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let mut gas = parcel(350.0, 1e-7, 100.0);
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let original = body::wall(&solid, &chart, &context, 450.0);
        let before = thermal.total_enthalpy_j().unwrap();
        let bits = thermal.flow().lattice().grid().state_bits();
        let gas_before = gas.state().clone();
        let wall_before = original.specific_enthalpy_j_kg();
        let wall_source = original.material().clone();
        let mut oversized = body::contact();
        oversized.area_m2 = 1e-6;
        // Six distinct exposed native wall faces total the full declared
        // body area. Adding a wetted contact must refuse shared coverage.
        let covered: Vec<_> = (1..=6)
            .map(|z| FreeSurfaceBodyGasFace {
                wall: [-1, 0, z],
                gas: [0, 0, z as usize],
                ..solid_face()
            })
            .collect();
        let mut wet = solid_face();
        wet.wall = [-1, 3, 3];
        wet.gas = [0, 3, 3];
        for (faces, contacts) in [
            (vec![solid_face(), solid_face()], vec![]),
            (vec![wet], vec![]),
            (covered, vec![oversized]),
        ] {
            let mut bodies = [original.clone().with_gas_faces(faces)];
            assert!(
                thermal
                    .step_with_finite_gas_environment(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &contacts,
                        &[],
                        &mut gas,
                        &[exposed],
                        gas_controls(),
                        false,
                        gate,
                    )
                    .is_err()
            );
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
            assert_eq!(bodies[0].material(), &wall_source);
            assert_eq!(gas.state(), &gas_before);
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        }
        let mut bodies = [original.with_gas_faces(vec![solid_face()])];
        assert!(
            thermal
                .step_with_body_conduction(&[], &[], controls(), &mut bodies, &[], gate,)
                .is_err(),
            "a body cannot silently lose its finite atmosphere"
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[],
                    &[],
                    &mut gas,
                    &[exposed],
                    gas_controls(),
                    false,
                    &cancelled,
                )
                .is_err()
        );
        // Restrict gas energy inversion after admitted initial film queries;
        // heat arrives from the hot solid and crosses this bound during a step.
        let tight = GasEnergySolve {
            maximum_temperature_k: 350.0001,
            ..gas_controls()
        };
        assert!(
            thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[],
                    &[],
                    &mut gas,
                    &[exposed],
                    tight,
                    false,
                    gate,
                )
                .is_err()
        );
        assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.accepted_steps(), 0);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
        assert_eq!(bodies[0].material(), &wall_source);
        assert_eq!(gas.state(), &gas_before);
        thermal
            .step_with_finite_gas_environment(
                &[],
                &[],
                controls(),
                &mut bodies,
                &[],
                &[],
                &mut gas,
                &[],
                gas_controls(),
                false,
                gate,
            )
            .unwrap();
        assert_eq!(thermal.accepted_steps(), 1);
        assert!(gas.state().state().temperature > 350.0);
        assert!(bodies[0].specific_enthalpy_j_kg() < wall_before);
    });
}
