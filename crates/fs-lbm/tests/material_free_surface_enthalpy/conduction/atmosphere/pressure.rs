//! Original finite gas heat/EOS pressure drives the SAME native liquid step.
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceInterfaceFace, FreeSurfaceInterfaceSolve,
};

// The 8^3 voxel vessel contains 4*6*6 fully filled liquid cells initially.
// The physical gas volume is (512-144)*dx^3, fixed at constant liquid mass/rho.
// It is never obtained from the evolving numerical fill or lattice density.
const GAS_VOLUME: f64 = 3.68e-7;

fn with_pressure<R>(
    card: &MaterialCard,
    reference: Option<f64>,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    with_exposed_law_and_pressure(
        card, false, &caloric(card), 200000.0, &law, reference, run,
    )
}

// Independent four-step caloric/film calculation. The fixture has
// cv/R=2+0.001*T, k_liquid=T-200 and pressure-dependent gas conductivity.
fn thermal_reference(initial_gas_k: f64, volume: f64, conductivity: f64) -> (f64, f64) {
    let r = 8.314_462_618_153_24 / 0.028;
    let mass = 100000.0 * volume / (r * initial_gas_k);
    let (mut liquid, mut gas) = (400.0, initial_gas_k);
    let mut supplied = 0.0;
    for _ in 0..4 {
        let pressure = 100000.0 * gas / initial_gas_k;
        let kg = 0.026 * conductivity * (0.5 + (pressure - 1000.0) / 249000.0);
        let conductance = 1e-6 / (0.0005 / (liquid - 200.0) + 0.0005 / kg);
        let heat = 0.00025 * conductance * (gas - liquid);
        liquid += heat / 0.002;
        gas = (((2.0 + 0.001 * gas).powi(2) - 0.002 * heat / (mass * r)).sqrt() - 2.0)
            / 0.001;
        supplied += heat;
    }
    (gas, supplied)
}

fn front_velocity(thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> f64 {
    thermal.flow().sample_at_si(3, 3, 3).unwrap().unwrap().velocity_m_s[0]
}

#[test]
fn g1_finite_gas_heating_and_cooling_drive_the_post_thermal_pressure_traction() {
    let liquid = conducting_card(1.0);
    for initial in [350.0, 450.0] {
        with_pressure(&liquid, Some(100000.0), |thermal, gate| {
            let mut unloaded = thermal.clone();
            let mut gas = parcel(initial, GAS_VOLUME, 0.01);
            let mut unpressurized_gas = gas.clone();
            let original_mass = gas.mass_kg();
            let exposed = gas_contact(gas.model(), gas.point()).face;
            let without = unloaded.step_with_finite_gas_environment(
                &[], &[], controls(), &mut [], &[], &[], &mut unpressurized_gas,
                &[exposed], gas_controls(), false, gate,
            ).unwrap().conduction;
            let actual = thermal.step_with_pressurized_finite_gas_environment(
                &[], &[], controls(), &mut [], &[], &[], &mut gas,
                &[exposed], gas_controls(), 0.05, false, gate,
            ).unwrap().conduction;
            let (expected_t, expected_heat) = thermal_reference(initial, GAS_VOLUME, 0.01);
            let expected_pressure = 100000.0 * expected_t / initial;
            near(gas.state().state().temperature, expected_t, 2e-11);
            near(gas.state().state().pressure, expected_pressure, 2e-11);
            near(actual.environment.as_ref().unwrap().heat_to_liquid_j, expected_heat, 2e-12);
            assert_eq!(actual.thermal.substeps, 4);
            assert_eq!(gas.state(), unpressurized_gas.state());
            assert_eq!(actual.thermal.stored_energy_change_j, without.thermal.stored_energy_change_j);
            assert_eq!(actual.thermal.external_heat_j, without.thermal.external_heat_j);
            assert_eq!(actual.thermal.substeps, without.thermal.substeps);
            let finite = actual.finite_gas.as_ref().unwrap();
            let pressure = finite.pressure.unwrap();
            assert_eq!(pressure.initial_pressure_pa, 100000.0);
            assert_eq!(pressure.applied_pressure_pa, gas.state().state().pressure);
            assert_eq!(pressure.volume_m3, GAS_VOLUME);
            assert_eq!(pressure.maximum_density_excursion, 0.05);
            assert_eq!(without.finite_gas.as_ref().unwrap().pressure, None);
            assert_eq!(gas.mass_kg(), original_mass);
            assert_eq!(gas.volume_m3(), GAS_VOLUME);
            // At this planar front the five incoming -x gas weights sum to
            // 1/6. With rho_b=1+3*delta_p the independently summed populations
            // have rho=1+delta_p and momentum_x=-delta_p. Here dx/dt=1 m/s.
            let perturbation = (expected_pressure - 100000.0) / 1000.0;
            assert!((front_velocity(thermal) + perturbation / (1.0 + perturbation)).abs() < 2e-9);
            assert!(front_velocity(&unloaded).abs() < 2e-13);
            assert!((initial - 400.0) * front_velocity(thermal) > 0.0);
            assert_ne!(
                thermal.flow().lattice().grid().state_bits(),
                unloaded.flow().lattice().grid().state_bits(),
            );
            assert!(finite.gas_balance_defect_j.abs() <= finite.gas_energy_tolerance_j);
            assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
            assert_eq!(thermal.physical_time_s(), 0.001);
        });
    }
}

#[test]
fn g1_radiant_surface_heat_pressurizes_the_same_finite_atmosphere_and_liquid() {
    let liquid = conducting_card(1.0);
    with_pressure(&liquid, Some(100000.0), |thermal, gate| {
        let mut unloaded = thermal.clone();
        let mut gas = parcel(400.0, GAS_VOLUME, 1.0);
        let mut unpressurized_gas = gas.clone();
        let surface = FreeSurfaceInterfaceFace {
            face: gas_contact(gas.model(), gas.point()).face,
            emissivity: 0.8.into(),
            surroundings_k: 600.0,
        };
        let without = unloaded.step_with_finite_gas_surface_exchange(
            &[], &[], controls(), &mut [], &[], &[], &mut unpressurized_gas,
            &[surface], gas_controls(), FreeSurfaceInterfaceSolve::default(), false, gate,
        ).unwrap().conduction;
        let actual = thermal.step_with_pressurized_finite_gas_surface_exchange(
            &[], &[], controls(), &mut [], &[], &[], &mut gas, &[surface], gas_controls(),
            FreeSurfaceInterfaceSolve::default(), 0.05, false, gate,
        ).unwrap().conduction;
        assert!(gas.state().state().temperature > 400.0);
        assert!(gas.state().state().pressure > 100000.0);
        assert_eq!(gas.state(), unpressurized_gas.state());
        assert_eq!(actual.thermal.stored_energy_change_j, without.thermal.stored_energy_change_j);
        assert_eq!(actual.thermal.external_heat_j, without.thermal.external_heat_j);
        assert_eq!(actual.thermal.substeps, without.thermal.substeps);
        let exchange = actual.surface.as_ref().unwrap();
        assert!(exchange.heat_from_radiation_j > 0.0);
        assert!(exchange.heat_from_gas_j < 0.0);
        let finite = actual.finite_gas.as_ref().unwrap();
        assert!(finite.internal_energy_change_j > 0.0);
        assert!((finite.internal_energy_change_j + exchange.heat_from_gas_j).abs()
            <= finite.gas_energy_tolerance_j + 1e-14);
        assert_eq!(
            finite.pressure.unwrap().applied_pressure_pa,
            gas.state().state().pressure,
        );
        let perturbation = (gas.state().state().pressure - 100000.0) / 1000.0;
        assert!((front_velocity(thermal) + perturbation / (1.0 + perturbation)).abs() < 2e-9);
        assert!(front_velocity(thermal) < 0.0);
        assert!(front_velocity(&unloaded).abs() < 2e-13);
        assert!(finite.combined_thermal_balance_defect_j.abs() < 1e-9);
    });
}

#[test]
fn g3_finite_pressure_feedback_retains_energy_gauges_and_storage_size_causality() {
    let run = |reference, volume, repeats| {
        let liquid = conducting_card(1.0);
        with_pressure(&liquid, Some(100000.0), |thermal, gate| {
            let mut gas = IsochoricGasParcel::new(
                gas_source::model(reference, 250000.0, 0.001),
                gas_source::point(450.0, 100000.0),
                volume,
            ).unwrap();
            let exposed = gas_contact(gas.model(), gas.point()).face;
            for _ in 0..repeats {
                thermal.step_with_pressurized_finite_gas_environment(
                    &[], &[], controls(), &mut [], &[], &[], &mut gas,
                    &[exposed], gas_controls(), 0.05, false, gate,
                ).unwrap();
            }
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal.total_enthalpy_j().unwrap(),
                gas.state().state().temperature,
                gas.state().state().pressure,
                front_velocity(thermal),
                thermal.physical_time_s(),
            )
        })
    };
    let baseline = run(0.0, GAS_VOLUME, 3);
    for reference in [-1e12, 1e12] {
        assert_eq!(run(reference, GAS_VOLUME, 3), baseline);
    }
    let small = run(0.0, GAS_VOLUME, 1);
    // The larger case declares an additional connected rigid plenum, whose
    // volume is external physical input rather than inferred liquid expansion.
    let large = run(0.0, 2.0 * GAS_VOLUME, 1);
    assert!(small.2 < large.2 && small.3 < large.3);
    assert!(small.4 > large.4 && large.4 > 0.0);
    assert_ne!(small.0, large.0);
}

#[test]
fn g4_late_pressure_refusal_preserves_gas_liquid_solids_and_exact_continuation() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_pressure(&liquid, Some(100000.0), |thermal, gate| {
        let mut gas = parcel(450.0, GAS_VOLUME, 0.01);
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let mut walls = [body::wall(&solid, &chart, &context, 450.0)];
        let mut expected = thermal.clone();
        let mut expected_gas = gas.clone();
        let mut expected_walls = walls.clone();
        let bits = thermal.flow().lattice().grid().state_bits();
        let h = thermal.total_enthalpy_j().unwrap();
        let gas_before = gas.clone();
        let wall_h = walls[0].specific_enthalpy_j_kg();
        let wall_source = walls[0].material().clone();
        let time = thermal.physical_time_s();
        let error = thermal.step_with_pressurized_finite_gas_environment(
            &[], &[], controls(), &mut walls, &[body::contact()], &[], &mut gas,
            &[exposed], gas_controls(), 1e-8, false, gate,
        ).unwrap_err();
        assert!(format!("{error:?}").contains("InvalidGasPressure"));
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(thermal.total_enthalpy_j().unwrap(), h);
        assert_eq!(thermal.physical_time_s(), time);
        assert_eq!(walls[0].specific_enthalpy_j_kg(), wall_h);
        assert_eq!(walls[0].material(), &wall_source);
        assert_eq!(gas.state(), gas_before.state());
        assert_eq!(gas.point(), gas_before.point());
        assert_eq!(gas.mass_kg(), gas_before.mass_kg());
        assert_eq!(gas.volume_m3(), gas_before.volume_m3());

        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(thermal.step_with_pressurized_finite_gas_environment(
            &[], &[], controls(), &mut walls, &[body::contact()], &[], &mut gas,
            &[exposed], gas_controls(), 0.05, false, &cancelled,
        ).is_err());
        assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
        assert_eq!(walls[0].specific_enthalpy_j_kg(), wall_h);
        assert_eq!(gas.state(), gas_before.state());
        for (flow, parcel, solids) in [
            (&mut *thermal, &mut gas, &mut walls),
            (&mut expected, &mut expected_gas, &mut expected_walls),
        ] {
            flow.step_with_pressurized_finite_gas_environment(
                &[], &[], controls(), solids, &[body::contact()], &[], parcel,
                &[exposed], gas_controls(), 0.05, false, gate,
            ).unwrap();
        }
        assert_eq!(
            thermal.flow().lattice().grid().state_bits(),
            expected.flow().lattice().grid().state_bits(),
        );
        assert_eq!(thermal.total_enthalpy_j().unwrap(), expected.total_enthalpy_j().unwrap());
        assert_eq!(thermal.physical_time_s(), expected.physical_time_s());
        assert_eq!(gas.state(), expected_gas.state());
        assert_eq!(walls[0].specific_enthalpy_j_kg(), expected_walls[0].specific_enthalpy_j_kg());
        assert_eq!(walls[0].material(), expected_walls[0].material());
    });
    with_pressure(&liquid, None, |thermal, gate| {
        let mut gas = parcel(450.0, GAS_VOLUME, 0.01);
        let exposed = gas_contact(gas.model(), gas.point()).face;
        let before = gas.state().clone();
        let error = thermal.step_with_pressurized_finite_gas_environment(
            &[], &[], controls(), &mut [], &[], &[], &mut gas,
            &[exposed], gas_controls(), 0.05, false, gate,
        ).unwrap_err();
        assert!(format!("{error:?}").contains("pressure reference"));
        assert_eq!(gas.state(), &before);
        assert_eq!(thermal.physical_time_s(), 0.0);
    });
}
