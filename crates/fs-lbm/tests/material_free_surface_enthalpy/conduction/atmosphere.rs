//! Independent finite gas caloric reference, native data feedback and rollback.
use super::*;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, IsochoricGasParcel};

#[path = "atmosphere/solid.rs"]
mod solid;

fn gas_controls() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 690.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

fn parcel(temperature: f64, volume: f64, conductivity: f64) -> IsochoricGasParcel {
    // Nonzero formation reference must remain stored; it is never heat input.
    IsochoricGasParcel::new(
        gas_source::model(100000.0, 250000.0, conductivity),
        gas_source::point(temperature, 100000.0),
        volume,
    )
    .unwrap()
}

#[test]
fn g1_finite_gas_and_moving_liquid_match_independent_nasa_energy_and_pressure() {
    let liquid = conducting_card(1.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut gas = parcel(450.0, 1e-8, 100.0);
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let before = thermal.total_enthalpy_j().unwrap();
        let packs = gas.state().source_packs().to_vec();
        let report = thermal
            .step_with_finite_gas_environment(
                &[],
                &[],
                controls(),
                &mut [],
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
        let environment = report.environment.as_ref().unwrap();
        let finite = report.finite_gas.as_ref().unwrap();
        // Authored NASA9: cv/R=2+0.001T and p/rho=R*T. Independent inversion
        // of the integrated quadratic u(T) keeps the formation datum intact.
        let r = 8.314_462_618_153_24 / 0.028;
        let mass = 100000.0 * 1e-8 / (r * 450.0);
        let (mut tl, mut tg) = (400.0, 450.0);
        let mut total_heat = 0.0;
        let mut applied_g = 0.0;
        for _ in 0..4 {
            let pressure = 100000.0 * tg / 450.0;
            let kg = 2.6 * (0.5 + (pressure - 1000.0) / 249000.0);
            applied_g = 1e-6 / (0.0005 / (tl - 200.0) + 0.0005 / kg);
            let heat = 0.00025 * applied_g * (tg - tl);
            tl += heat / 0.002;
            tg =
                (((2.0_f64 + 0.001 * tg).powi(2) - 0.002 * heat / (mass * r)).sqrt() - 2.0) / 0.001;
            total_heat += heat;
        }
        assert_eq!(report.thermal.substeps, 4);
        assert_eq!(finite.substeps, 4);
        near(
            environment.last_applied[0].conductance_w_k,
            applied_g,
            2e-12,
        );
        near(
            environment.thermal_endpoint[0].liquid_temperature_k,
            tl,
            2e-11,
        );
        near(gas.state().state().temperature, tg, 2e-11);
        near(gas.state().state().pressure, 100000.0 * tg / 450.0, 2e-11);
        near(
            environment.thermal_endpoint[0].gas.state().pressure,
            gas.state().state().pressure,
            2e-13,
        );
        near(
            environment.thermal_endpoint[0].gas.state().temperature,
            0.5 * (tl + tg),
            2e-11,
        );
        near(environment.heat_to_liquid_j, total_heat, 2e-12);
        assert!(
            (finite.internal_energy_change_j + total_heat).abs()
                <= finite.gas_energy_tolerance_j + 1e-14
        );
        near(
            thermal.total_enthalpy_j().unwrap(),
            before + total_heat,
            2e-12,
        );
        near(report.flow.external_heat_j, total_heat, 2e-12);
        assert!(finite.gas_balance_defect_j.abs() <= finite.gas_energy_tolerance_j);
        assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
        near(gas.mass_kg(), mass, 2e-13);
        assert_eq!(gas.volume_m3(), 1e-8);
        assert_eq!(gas.state().source_packs(), packs.as_slice());
    });
}

#[test]
fn g3_gas_formation_reference_cannot_reject_or_change_liquid_heat_and_motion() {
    let run = |offset, temperature| {
        let liquid = conducting_card(1.0);
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let mut gas = IsochoricGasParcel::new(
                gas_source::model(offset, 250000.0, 100.0),
                gas_source::point(temperature, 100000.0), 1e-8,
            ).unwrap();
            let exposed = gas_contact(gas.model(), gas.point()).face;
            let report = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut [], &[], &[], &mut gas,
                &[exposed], gas_controls(), false, gate,
            ).unwrap().conduction;
            let finite = report.finite_gas.unwrap();
            let endpoint = gas.state().state().temperature;
            // Independent integral of the fixture's cv/R, with no absolute
            // formation reference in the expected transferred thermal energy.
            let expected = gas.mass_kg() * (8.314_462_618_153_24 / 0.028)
                * (endpoint - temperature) * (2.0 + 0.0005 * (endpoint + temperature));
            near(finite.internal_energy_change_j, expected, 2e-12);
            assert!(finite.gas_balance_defect_j.abs() <= finite.gas_energy_tolerance_j);
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal.total_enthalpy_j().unwrap(),
                endpoint,
                gas.state().state().pressure,
                finite,
            )
        })
    };
    for temperature in [350.0, 450.0] {
        let reference = run(0.0, temperature);
        assert_ne!(reference.2, temperature, "the finite gas must actually exchange heat");
        for offset in [-1.0e12, 1.0e12] {
            assert_eq!(run(offset, temperature), reference);
        }
    }
}

#[test]
fn g3_finite_atmosphere_storage_and_source_data_change_actual_liquid_motion() {
    let run = |volume, conductivity| {
        let liquid = conducting_card(1.0);
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let mut gas = parcel(450.0, volume, conductivity);
            let exposed = gas_contact(gas.model(), gas.point()).face;
            let report = thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut [],
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
                report.environment.unwrap().heat_to_liquid_j,
                gas.state().state().temperature,
            )
        })
    };
    let slow = run(1e-8, 10.0);
    let fast = run(1e-8, 1000.0);
    let larger = run(1e-6, 1000.0);
    assert!(fast.2 > slow.2 && larger.2 > fast.2);
    assert!(fast.1 < slow.1 && larger.1 < fast.1);
    assert!(fast.3 < slow.3 && larger.3 > fast.3);
    assert_ne!(
        slow.0, fast.0,
        "gas material must reach actual native populations"
    );
    assert_ne!(
        fast.0, larger.0,
        "finite gas capacity must reach actual native populations"
    );
}

#[test]
fn g1_finite_gas_wall_interface_molecular_and_authored_heat_share_one_interval() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    let interface_source = body::interface_card(&liquid, &solid, 1e-5, 500.0, false);
    let interface_point = body::interface_point();
    let interface = body::interface(&liquid, &solid, &interface_source, &interface_point);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let mut gas = parcel(350.0, 1e-8, 100.0);
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let gas_before = gas.state().specific_internal_energy_j_kg();
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let before = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j();
        let heat = FreeSurfaceCellHeat {
            coordinates: EXPOSED,
            energy_j: 0.001,
        };
        let face = FreeSurfaceConductionFace {
            cells: [[2, 3, 3], EXPOSED],
            area_m2: 1e-6,
            distances_m: [0.0005; 2],
        };
        let report = thermal
            .step_with_finite_gas_environment(
                &[face],
                &[heat],
                controls(),
                &mut bodies,
                &[body::contact()],
                std::slice::from_ref(&interface),
                &mut gas,
                &[exposed],
                gas_controls(),
                true,
                gate,
            )
            .unwrap()
            .conduction;
        let environment = report.environment.as_ref().unwrap();
        let wall = report.bodies.as_ref().unwrap();
        let finite = report.finite_gas.as_ref().unwrap();
        let gas_change = gas.mass_kg() * (gas.state().specific_internal_energy_j_kg() - gas_before);
        assert!(environment.heat_to_liquid_j < 0.0);
        assert!(gas_change > 0.0 && wall.heat_to_liquid_j > 0.0);
        assert!(report.flow.viscous_heat_j > 0.0);
        near(
            thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j() - before + gas_change,
            heat.energy_j + report.flow.viscous_heat_j,
            2e-9,
        );
        near(
            report.flow.external_heat_j,
            heat.energy_j + environment.heat_to_liquid_j,
            2e-12,
        );
        assert!((gas_change + environment.heat_to_liquid_j).abs() <= finite.gas_energy_tolerance_j);
        assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
        assert_eq!(
            wall.last_applied[0]
                .interface
                .as_ref()
                .unwrap()
                .card_identity(),
            interface_source.content_hash()
        );
        assert_eq!(report.thermal.substeps, finite.substeps);
    });
}

#[test]
fn g1_all_exposed_patches_share_one_small_gas_capacity_and_monotone_step_bound() {
    let liquid = conducting_card(1.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut gas = parcel(450.0, 1e-9, 100.0);
        let a = gas_contact(gas.model(), gas.point()).face;
        let b = FreeSurfaceGasFace {
            liquid: [3, 2, 3],
            gas: [4, 2, 3],
            ..a
        };
        let maximum_heat = gas.mass_kg()
            * (8.314_462_618_153_24 / 0.028)
            * (2.0 * 50.0 + 0.0005 * (450.0 * 450.0 - 400.0 * 400.0));
        let report = thermal
            .step_with_finite_gas_environment(
                &[],
                &[],
                controls(),
                &mut [],
                &[],
                &[],
                &mut gas,
                &[a, b],
                gas_controls(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let environment = report.environment.as_ref().unwrap();
        assert!(
            report.thermal.substeps > 4,
            "finite gas capacity must tighten the original thermal clock"
        );
        assert!(gas.state().state().temperature >= 400.0 - 1e-9);
        assert!(gas.state().state().temperature < 450.0);
        assert!(environment.heat_to_liquid_j <= maximum_heat + 1e-12);
        for sample in &environment.thermal_endpoint {
            assert!(sample.liquid_temperature_k >= 400.0 && sample.liquid_temperature_k <= 450.0);
            near(
                sample.bath_temperature_k,
                gas.state().state().temperature,
                2e-13,
            );
        }
        assert!(
            report
                .finite_gas
                .unwrap()
                .combined_thermal_balance_defect_j
                .abs()
                < 1e-9
        );
    });
}

#[test]
fn g4_late_gas_film_refusal_and_cancellation_preserve_every_native_participant() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_exposed_liquid(&liquid, true, |thermal, gate| {
        let mut gas = IsochoricGasParcel::new(
            gas_source::model_with_temperature_limit(100000.0, 250000.0, 100.0, 450.0),
            gas_source::point(450.0, 100000.0),
            1e-7,
        )
        .unwrap();
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let bits = thermal.flow().lattice().grid().state_bits();
        let before = thermal.total_enthalpy_j().unwrap();
        let wall_before = bodies[0].specific_enthalpy_j_kg();
        let wall_source = bodies[0].material().clone();
        let gas_before = gas.state().clone();
        let gas_point = gas.point().clone();
        let controls_gas = GasEnergySolve {
            maximum_temperature_k: 450.0,
            ..gas_controls()
        };
        // Several admitted substeps precede an unsupported film above 450 K.
        // The original parcel and solid cannot publish those partial updates.
        let heat = FreeSurfaceCellHeat {
            coordinates: EXPOSED,
            energy_j: 0.18,
        };
        let error = thermal
            .step_with_finite_gas_environment(
                &[],
                &[heat],
                controls(),
                &mut bodies,
                &[body::contact()],
                &[],
                &mut gas,
                &[exposed],
                controls_gas,
                false,
                gate,
            )
            .unwrap_err();
        assert!(matches!(
            error,
            FreeSurfaceEnthalpyError::Flow(fs_lbm::material::MaterialLbmError::Transport(_))
        ));
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), before);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
        assert_eq!(bodies[0].material(), &wall_source);
        assert_eq!(gas.state(), &gas_before);
        assert_eq!(gas.point(), &gas_point);
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[body::contact()],
                    &[],
                    &mut gas,
                    &[exposed],
                    controls_gas,
                    false,
                    &cancelled,
                )
                .is_err()
        );
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
        assert_eq!(gas.state(), &gas_before);
        thermal
            .step_with_finite_gas_environment(
                &[],
                &[],
                controls(),
                &mut bodies,
                &[body::contact()],
                &[],
                &mut gas,
                &[exposed],
                controls_gas,
                false,
                gate,
            )
            .unwrap();
        assert!(gas.state().state().temperature < gas_before.state().temperature);
    });
}
