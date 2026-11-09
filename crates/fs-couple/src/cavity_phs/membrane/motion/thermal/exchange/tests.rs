//! Source parity, physical stream/work accounting and joint publication.
use super::*;
use crate::cavity_phs::membrane::tests::motion::{motion_controls, real_gas};
use crate::cavity_phs::membrane::tests::{near, pyramid, source_gas};
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_material::thermochemical_gas::nonideal::energy::parcel::transfer::RealGasMassTransferSolve;

fn controls() -> GasMembraneThermalControls {
    let old = motion_controls();
    GasMembraneThermalControls {
        gas: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        exterior_gas: None,
        time: old.time,
        maximum_evaluations: 2000,
        maximum_steps: old.maximum_steps,
        maximum_krylov_restart: old.maximum_krylov_restart,
        energy_tolerance_j: 1e-5,
        pressure_work_tolerance_j: 1e-8,
    }
}
fn ideal_exchange() -> GasNozzleStepSolve {
    GasNozzleStepSolve {
        nozzle: controls().gas,
        transfer: GasMassTransferSolve {
            donor: controls().gas,
            receiver: controls().gas,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        maximum_donor_fraction: 0.05,
    }
}
fn real_exchange() -> RealGasNozzleStepSolve {
    RealGasNozzleStepSolve {
        entropy: motion_controls().gas,
        transfer: RealGasMassTransferSolve {
            donor: controls().gas,
            receiver: controls().gas,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-7,
        },
        minimum_density_fraction: 0.5,
        pressure_tolerance_pa: 1e-4,
        sonic_tolerance_m2_s2: 1e-4,
        maximum_iterations: 128,
        maximum_donor_fraction: 0.05,
    }
}
fn nozzle() -> GasNozzle {
    GasNozzle {
        area_m2: 1e-6,
        discharge_coefficient: 0.8,
    }
}
fn reservoir(volume: f64) -> MembraneGas {
    let mut gas = source_gas(0.25, 0.002);
    gas.set_volume_isentropic(volume, motion_controls().gas, || false)
        .unwrap();
    gas.into()
}
fn stored(gas: &MembraneGas) -> f64 {
    gas.mass_kg()
        * match gas {
            MembraneGas::Ideal(gas) => gas.state().specific_internal_energy_j_kg(),
            MembraneGas::PengRobinson(gas) => gas.state().specific_internal_energy_j_kg(),
        }
}

#[test]
fn g1_inflow_changes_actual_mixture_pressure_motion_and_continuing_history() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        100000.0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let loads = [[0.0; 3]; 5];
    let mut motion = GasMembraneMotion::new(
        &model,
        &reference,
        &[[0.0; 3]; 5],
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    motion
        .step_discrete_gradient(1e-5, &loads, 0.0, None, controls(), || false)
        .unwrap();
    let old = motion.accepted().clone();
    let old_state = motion.state().clone();
    let old_kinetic = motion.kinetic_energy_j().unwrap();
    let mut tank = reservoir(0.0005);
    let old_tank = tank.clone();
    let (MembraneGas::Ideal(a), MembraneGas::Ideal(b)) = (&old.gas, &old_tank) else {
        unreachable!()
    };
    let (native, flow) = b
        .through_nozzle_to(a, nozzle(), 1e-5, ideal_exchange(), || false)
        .unwrap();
    let (expected_tank, post_transfer, transfer) = native.into_parts();
    let report = motion
        .step_discrete_gradient_with_nozzle(
            &mut tank,
            nozzle(),
            1e-5,
            &loads,
            0.01,
            None,
            GasMembraneNozzleSolve::Ideal(ideal_exchange()),
            controls(),
            1e-5,
            || false,
        )
        .unwrap();
    assert_eq!(
        report.direction,
        Some(GasMembraneNozzleDirection::IntoChamber)
    );
    assert!(flow.choked);
    assert_eq!(report.flow, GasMembraneNozzleFlow::Ideal(flow));
    assert_eq!(report.transfer, transfer);
    assert_eq!(tank, MembraneGas::Ideal(expected_tank));
    assert!(motion.accepted().gas.mass_kg() > old.gas.mass_kg());
    assert!(motion.state().v[0] > old_state.v[0]);
    assert!(motion.accepted().pressure.pressure_pa > old.pressure.pressure_pa);
    assert_eq!(motion.state().steps, old_state.steps + 1);
    assert_eq!(
        &motion.state().history[..old_state.history.len()],
        &old_state.history
    );
    assert_eq!(report.motion.gas.supplied_heat_j, 0.01);

    let native = post_transfer
        .discrete_gradient_heat_volume(
            motion.accepted().gas.volume_m3(),
            0.01,
            controls().gas,
            || false,
        )
        .unwrap();
    assert_eq!(
        report.motion.work_pressure_pa,
        native.step().work_pressure_pa
    );
    let (expected_chamber, _) = native.into_parts();
    assert_eq!(motion.accepted().gas, MembraneGas::Ideal(expected_chamber));
    let MembraneGas::Ideal(after) = &motion.accepted().gas else {
        unreachable!()
    };
    assert_ne!(
        after.state().state().specific_gas_constant,
        a.state().state().specific_gas_constant
    );
    assert_eq!(after.point().axes().get("specimen-lot"), Some(&42.0));
    near(
        report.reservoir_energy_change_j,
        -report.transfer.transferred_enthalpy_j,
        1e-8,
        0.0,
    );
    near(
        report.chamber_energy_change_j + report.motion.gas.work_by_gas_j,
        report.transfer.transferred_enthalpy_j + 0.01,
        1e-8,
        0.0,
    );
    let independent = stored(&motion.accepted().gas) - stored(&old.gas) + stored(&tank)
        - stored(&old_tank)
        + motion.kinetic_energy_j().unwrap()
        - old_kinetic
        + motion.accepted().strain_energy_j
        - old.strain_energy_j
        + 100000.0 * report.motion.swept_volume_change_m3
        - 0.01;
    near(report.total_energy_residual_j, independent, 2e-9, 0.0);
    assert!(report.total_mass_residual_kg.abs() <= 1e-14);
    assert!(report.motion.pressure_work_defect_j.abs() <= controls().pressure_work_tolerance_j);

    let composition = after
        .state()
        .thermodynamics()
        .receipt()
        .mole_fraction_bits()
        .to_vec();
    let mass = after.mass_kg();
    motion
        .step_discrete_gradient(1e-5, &loads, 0.0, None, controls(), || false)
        .unwrap();
    let MembraneGas::Ideal(after) = &motion.accepted().gas else {
        unreachable!()
    };
    assert_eq!(after.mass_kg(), mass);
    assert_eq!(
        after
            .state()
            .thermodynamics()
            .receipt()
            .mole_fraction_bits(),
        composition.as_slice()
    );
    assert_eq!(motion.state().history.len(), old_state.history.len() + 2);
}

#[test]
fn g1_outflow_and_heated_opposing_gas_keep_three_inventory_work_balances() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::with_sealed_exterior(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        [0.0007, 0.0015],
        SealedMembraneExterior {
            gas: source_gas(0.75, 0.002).into(),
            volume_interval_m3: [0.0007, 0.0015],
            entropy: motion_controls().gas,
        },
        || false,
    )
    .unwrap();
    let settings = GasMembraneThermalControls {
        exterior_gas: Some(controls().gas),
        ..controls()
    };
    let loads = [[0.0; 3]; 5];
    let mut motion = GasMembraneMotion::new(
        &model,
        &reference,
        &[[0.0; 3]; 5],
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    let old = motion.accepted().clone();
    let mut tank = reservoir(0.002);
    let old_tank = tank.clone();
    let report = motion
        .step_discrete_gradient_with_nozzle(
            &mut tank,
            nozzle(),
            1e-5,
            &loads,
            0.01,
            Some(0.02),
            GasMembraneNozzleSolve::Ideal(ideal_exchange()),
            settings,
            1e-5,
            || false,
        )
        .unwrap();
    assert_eq!(
        report.direction,
        Some(GasMembraneNozzleDirection::OutOfChamber)
    );
    assert!(motion.accepted().gas.mass_kg() < old.gas.mass_kg());
    assert!(tank.mass_kg() > old_tank.mass_kg());
    let before_other = old.exterior_gas.as_ref().unwrap();
    let after_other = motion.accepted().exterior_gas.as_ref().unwrap();
    assert_eq!(before_other.mass_kg(), after_other.mass_kg());
    near(
        after_other.volume_m3() - before_other.volume_m3(),
        -report.motion.swept_volume_change_m3,
        2e-18,
        1e-8,
    );
    assert_eq!(report.motion.exterior_pressure_work_j, 0.0);
    assert_eq!(report.motion.gas.supplied_heat_j, 0.01);
    assert_eq!(report.motion.exterior_gas.unwrap().supplied_heat_j, 0.02);
    near(
        report.reservoir_energy_change_j,
        report.transfer.transferred_enthalpy_j,
        1e-8,
        0.0,
    );
    near(
        report.chamber_energy_change_j + report.motion.gas.work_by_gas_j,
        0.01 - report.transfer.transferred_enthalpy_j,
        1e-8,
        0.0,
    );
    let independent = stored(&motion.accepted().gas) - stored(&old.gas) + stored(&tank)
        - stored(&old_tank)
        + stored(after_other)
        - stored(before_other)
        + motion.kinetic_energy_j().unwrap()
        + motion.accepted().strain_energy_j
        - old.strain_energy_j
        - 0.01
        - 0.02;
    near(report.total_energy_residual_j, independent, 2e-9, 0.0);
    let pressure = motion
        .accepted()
        .gas
        .pressure_response()
        .unwrap()
        .pressure_pa;
    let other_pressure = after_other.pressure_response().unwrap().pressure_pa;
    assert_eq!(motion.accepted().pressure.pressure_pa, pressure);
    assert_eq!(motion.accepted().exterior_pressure_pa, other_pressure);
    let internal = patch
        .evaluate(&motion.accepted().displacement_m, || false)
        .unwrap();
    let pressure_force = patch
        .pressure_load_n(
            &motion.accepted().displacement_m,
            pressure - other_pressure,
            || false,
        )
        .unwrap();
    // Includes prescribed support rows, not only the free coordinate.
    for (i, row) in motion.accepted().residual_n.iter().enumerate() {
        for axis in 0..3 {
            near(
                row[axis],
                internal.internal_gradient_n[i][axis] - pressure_force[i][axis],
                1e-10,
                1e-12,
            );
        }
    }
}

#[test]
fn g3_pr_stream_and_gibbs_work_reuse_both_original_nonideal_owners() {
    let chamber = real_gas();
    let exterior = chamber.state().pressure_pa();
    let mut tank = real_gas();
    tank.add_heat(tank.mass_kg() * 30000.0, controls().gas, || false)
        .unwrap();
    let old_tank = tank.clone();
    let dt = 1e-6;
    let opening = GasNozzle {
        area_m2: 2e-7,
        ..nozzle()
    };
    let (native, flow) = old_tank
        .clone()
        .into_isochoric()
        .through_nozzle_to(
            &chamber.clone().into_isochoric(),
            opening,
            dt,
            real_exchange(),
            || false,
        )
        .unwrap();
    let (expected_tank, post_transfer, transfer) = native.into_parts();
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        chamber.clone(),
        &reference,
        &fixed,
        exterior,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let loads = [[0.0; 3]; 5];
    let mut motion = GasMembraneMotion::new(
        &model,
        &reference,
        &[[0.0; 3]; 5],
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    let old = motion.accepted().clone();
    let mut tank = MembraneGas::PengRobinson(tank);
    let report = motion
        .step_discrete_gradient_with_nozzle(
            &mut tank,
            opening,
            dt,
            &loads,
            0.2,
            None,
            GasMembraneNozzleSolve::PengRobinson(real_exchange()),
            controls(),
            1e-5,
            || false,
        )
        .unwrap();
    assert_eq!(
        report.direction,
        Some(GasMembraneNozzleDirection::IntoChamber)
    );
    assert_eq!(report.flow, GasMembraneNozzleFlow::PengRobinson(flow));
    assert_eq!(report.transfer, transfer);
    assert_eq!(
        tank,
        MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(expected_tank))
    );
    let source = post_transfer
        .discrete_gradient_heat_volume(
            motion.accepted().gas.volume_m3(),
            0.2,
            controls().gas,
            || false,
        )
        .unwrap();
    assert_eq!(
        report.motion.work_pressure_pa,
        source.step().work_pressure_pa
    );
    assert_eq!(
        report.motion.gas_sampling_pressure_pa,
        source.thermodynamic_mean().pressure_pa()
    );
    assert_eq!(
        report.motion.gas_sampling_temperature_k,
        source.thermodynamic_mean().temperature_k()
    );
    let (endpoint, _) = source.into_parts();
    assert_eq!(
        motion.accepted().gas,
        MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(endpoint)),
    );
    assert_eq!(
        motion.accepted().pressure.pressure_pa,
        motion
            .accepted()
            .gas
            .pressure_response()
            .unwrap()
            .pressure_pa
    );
    let independent = stored(&motion.accepted().gas) - stored(&old.gas) + stored(&tank)
        - stored(&MembraneGas::PengRobinson(old_tank))
        + motion.kinetic_energy_j().unwrap()
        + motion.accepted().strain_energy_j
        - old.strain_energy_j
        + exterior * report.motion.swept_volume_change_m3
        - 0.2;
    near(report.total_energy_residual_j, independent, 2e-8, 0.0);
    assert!(report.total_mass_residual_kg.abs() <= 1e-14);
}

#[test]
fn g3_closed_nozzle_preserves_exact_native_motion_and_accepted_history() {
    for source in [
        MembraneGas::Ideal(source_gas(0.75, 0.002)),
        MembraneGas::PengRobinson(real_gas()),
    ] {
        let solve = match &source {
            MembraneGas::Ideal(_) => GasMembraneNozzleSolve::Ideal(ideal_exchange()),
            MembraneGas::PengRobinson(_) => GasMembraneNozzleSolve::PengRobinson(real_exchange()),
        };
        let (patch, reference, fixed) = pyramid();
        let model = GasMembrane::new(
            &patch,
            source.clone(),
            &reference,
            &fixed,
            source.pressure_response().unwrap().pressure_pa,
            [0.0007, 0.0015],
            || false,
        )
        .unwrap();
        let loads = [[0.0; 3]; 5];
        let mut prefix = GasMembraneMotion::new(
            &model,
            &reference,
            &[[0.0; 3]; 5],
            &loads,
            motion_controls(),
            || false,
        )
        .unwrap()
        .into_thermal();
        prefix
            .step_discrete_gradient(1e-6, &loads, 0.01, None, controls(), || false)
            .unwrap();
        for opening in [
            GasNozzle {
                area_m2: 0.0,
                ..nozzle()
            },
            GasNozzle {
                discharge_coefficient: 0.0,
                ..nozzle()
            },
        ] {
            let mut expected = prefix.clone();
            let direct = expected
                .step_discrete_gradient(1e-6, &loads, 0.02, None, controls(), || false)
                .unwrap();
            let mut actual = prefix.clone();
            let mut tank = source.clone();
            let report = actual
                .step_discrete_gradient_with_nozzle(
                    &mut tank,
                    opening,
                    1e-6,
                    &loads,
                    0.02,
                    None,
                    solve,
                    controls(),
                    1e-5,
                    || false,
                )
                .unwrap();
            assert_eq!(report.direction, None);
            assert_eq!(report.transfer.transferred_mass_kg, 0.0);
            assert_eq!(tank, source);
            assert_eq!(actual.accepted(), expected.accepted());
            assert_eq!(actual.state(), expected.state());
            assert_eq!(report.motion.work_pressure_pa, direct.work_pressure_pa);
            assert_eq!(report.total_energy_residual_j, direct.integration_defect_j);
        }
    }
}

#[test]
fn g4_source_native_solve_and_last_poll_refusals_publish_neither_gas_nor_history() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        100000.0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let loads = [[0.0; 3]; 5];
    let mut original = GasMembraneMotion::new(
        &model,
        &reference,
        &[[0.0; 3]; 5],
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap()
    .into_thermal();
    original
        .step_discrete_gradient(1e-5, &loads, 0.01, None, controls(), || false)
        .unwrap();
    let original_tank = reservoir(0.0005);
    let mut expected = original.clone();
    let mut expected_tank = original_tank.clone();
    let mut polls = 0;
    let reference_step = expected
        .step_discrete_gradient_with_nozzle(
            &mut expected_tank,
            nozzle(),
            1e-5,
            &loads,
            0.01,
            None,
            GasMembraneNozzleSolve::Ideal(ideal_exchange()),
            controls(),
            1e-5,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 2);
    let mut actual = original.clone();
    let mut tank = original_tank.clone();
    let mut calls = 0;
    assert!(matches!(
        actual.step_discrete_gradient_with_nozzle(
            &mut tank,
            nozzle(),
            1e-5,
            &loads,
            0.01,
            None,
            GasMembraneNozzleSolve::Ideal(ideal_exchange()),
            controls(),
            1e-5,
            || {
                calls += 1;
                calls == polls
            },
        ),
        Err(GasMembraneMotionError::Cancelled),
    ));
    assert_eq!(actual.state(), original.state());
    assert_eq!(actual.accepted(), original.accepted());
    assert_eq!(tank, original_tank);
    for (solve, settings) in [
        (
            GasMembraneNozzleSolve::Ideal(GasNozzleStepSolve {
                maximum_donor_fraction: 1e-16,
                ..ideal_exchange()
            }),
            controls(),
        ),
        (
            GasMembraneNozzleSolve::Ideal(ideal_exchange()),
            GasMembraneThermalControls {
                maximum_evaluations: 1,
                ..controls()
            },
        ),
        (
            GasMembraneNozzleSolve::PengRobinson(real_exchange()),
            controls(),
        ),
    ] {
        assert!(
            actual
                .step_discrete_gradient_with_nozzle(
                    &mut tank,
                    nozzle(),
                    1e-5,
                    &loads,
                    0.01,
                    None,
                    solve,
                    settings,
                    1e-5,
                    || false,
                )
                .is_err()
        );
        assert_eq!(actual.state(), original.state());
        assert_eq!(actual.accepted(), original.accepted());
        assert_eq!(tank, original_tank);
    }
    if reference_step.total_energy_residual_j != 0.0 {
        assert!(matches!(
            actual.step_discrete_gradient_with_nozzle(
                &mut tank,
                nozzle(),
                1e-5,
                &loads,
                0.01,
                None,
                GasMembraneNozzleSolve::Ideal(ideal_exchange()),
                controls(),
                0.5 * reference_step.total_energy_residual_j.abs(),
                || false,
            ),
            Err(GasMembraneMotionError::Energy { .. }),
        ));
        assert_eq!(actual.state(), original.state());
        assert_eq!(actual.accepted(), original.accepted());
        assert_eq!(tank, original_tank);
    }
    let retry = actual
        .step_discrete_gradient_with_nozzle(
            &mut tank,
            nozzle(),
            1e-5,
            &loads,
            0.01,
            None,
            GasMembraneNozzleSolve::Ideal(ideal_exchange()),
            controls(),
            1e-5,
            || false,
        )
        .unwrap();
    assert_eq!(actual.accepted(), expected.accepted());
    assert_eq!(actual.state(), expected.state());
    assert_eq!(tank, expected_tank);
    assert_eq!(retry.transfer, reference_step.transfer);
    assert_eq!(
        retry.total_energy_residual_j,
        reference_step.total_energy_residual_j
    );
}
