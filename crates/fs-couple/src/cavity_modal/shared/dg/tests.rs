use super::*;
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh, MAX_DG_PROPAGATOR_DOFS, PortId, PortInput};
use fs_duct::{Duct, Segment};
use fs_material::gas::GasSpec;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}
fn config(pressure: f64) -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 32,
        max_connections: 32,
        max_setup_terms: 262144,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 0.01,
        maximum_abs_pressure_pa: pressure,
        maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-14,
        energy_relative_tolerance: 1e-8,
    }
}
fn cavity(rate: u32, memory: bool, pressure: f64) -> SharedCavity {
    let source = HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: 1e-6,
            neck_radius_m: 0.001,
            neck_length_m: 0.02,
            inner_end: fs_scenario::NeckEnd::Flanged,
            outer_end: fs_scenario::NeckEnd::Flanged,
        },
        gas(),
        CavityNeckLoss::ExteriorMouth,
    )
    .unwrap();
    let surfaces = vec![CavitySurface {
        model: ModalAcousticTimeModel::try_new(
            rate,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 1400.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap(),
        outward_volume_shapes: vec![1e-4],
    }];
    let initial = ModalAcousticState {
        displacement_m_sqrt_kg: 3e-8,
        velocity_m_sqrt_kg_per_s: 1e-4,
    };
    let gate = CancelGate::new_clock_free();
    if memory {
        SharedCavity::new_with_cylindrical_thermal_wall(
            source,
            surfaces,
            rate,
            initial,
            ModalAcousticTimeBudget::audible_reference(),
            config(pressure),
            Some(ViscousNeckOptions {
                poles: 8,
                band_rad_s: [100.0, 10000.0],
                maximum_relative_error: 0.1,
            }),
            CylindricalCavityThermalOptions {
                radius_m: 0.003,
                poles: 8,
                band_rad_s: [100.0, 10000.0],
                maximum_relative_error: 0.1,
            },
            Some(vec![0.0; 8]),
            &gate,
        )
        .unwrap()
    } else {
        SharedCavity::new(
            source,
            surfaces,
            rate,
            initial,
            ModalAcousticTimeBudget::audible_reference(),
            config(pressure),
            &gate,
        )
        .unwrap()
    }
}
fn load(rate: u32, walls: bool, radius: f64) -> DgPortLoad {
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius,
                length: 0.05,
            }],
        },
        &gas(),
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::IdealOpen,
    )
    .unwrap();
    let duct = if walls {
        duct.with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap()
    } else {
        duct
    };
    DgPortLoad::new(duct, 1.0 / f64::from(rate), 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap()
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + 2e-10 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn downstream_duct_retains_both_original_propagators_and_real_wall_history() {
    let gate = CancelGate::new_clock_free();
    let mut original_cavity = cavity(48000, true, 1e6);
    let mut coupled = cavity(48000, true, 1e6)
        .with_downstream_duct(load(48000, true, 0.001), &gate)
        .unwrap();
    let propagator = coupled.load().propagator().clone();
    let mut original_duct = coupled.load().state().clone();
    let initial = coupled.stored_energy_j().unwrap();
    let (mut loss, mut work) = (0.0, 0.0);
    for tick in 0..32 {
        let force = if tick < 8 { 0.001 } else { 0.0 };
        let frame = coupled.step(&[force], 0.0, &gate).unwrap().clone();
        original_cavity
            .step(&[force], frame.mean_pressure_pa, &gate)
            .unwrap();
        let expected = propagator
            .step(
                &mut original_duct,
                &[PortInput {
                    id: PortId(7),
                    effort: frame.mean_flow_m3_s,
                }],
            )
            .unwrap();
        propagator.sync_nodal(&mut original_duct).unwrap();
        assert_eq!(frame.duct, expected);
        assert_eq!(coupled.load().state(), &original_duct);
        assert_eq!(
            coupled.cavity().network().material_memory(),
            original_cavity.network().material_memory()
        );
        for (a, b) in coupled
            .cavity()
            .network()
            .components()
            .iter()
            .zip(original_cavity.network().components())
        {
            assert_eq!(a.states(), b.states());
        }
        let internal = frame.mean_pressure_pa * frame.mean_flow_m3_s / 48000.0;
        close(
            frame.modal.external_work_j,
            frame.authored_external_work_j - internal,
            1e-22,
        );
        close(frame.duct.ledger.port_work_step_j, internal, 1e-22);
        assert!(frame.port_work_residual_j.abs() <= frame.balance_tolerance_j);
        assert!(frame.balance_residual_j.abs() <= frame.balance_tolerance_j);
        assert!(!frame.duct.ledger.balance_checked);
        assert_eq!(frame.terminal_work_j, 0.0);
        assert_eq!(frame.modal.sample, coupled.load().state().step_count());
        close(
            frame.stored_energy_j,
            coupled.stored_energy_j().unwrap(),
            1e-22,
        );
        loss += frame.inferred_total_dissipation_j;
        work += frame.authored_external_work_j;
    }
    assert!(coupled.load().state().wall_stored_energy_j() > 0.0);
    assert!(
        coupled
            .cavity()
            .network()
            .material_memory()
            .iter()
            .any(|x| *x != 0.0)
    );
    assert!(
        coupled
            .cavity()
            .network()
            .last_frame()
            .unwrap()
            .material_energy_j
            > 0.0
    );
    assert!(loss > 0.0);
    close(
        coupled.stored_energy_j().unwrap() + loss,
        initial + work,
        1e-17,
    );
}

fn trajectory(rate: u32) -> Vec<f64> {
    let gate = CancelGate::new_clock_free();
    let mut coupled = cavity(rate, false, 1e6)
        .with_downstream_duct(load(rate, false, 0.001), &gate)
        .unwrap();
    for _ in 0..(rate / 375) {
        coupled.step(&[0.0], 0.0, &gate).unwrap();
    }
    let mut result = Vec::new();
    for state in coupled
        .cavity()
        .network()
        .components()
        .iter()
        .flat_map(|m| m.states())
    {
        result.extend([
            1400.0 * state.displacement_m_sqrt_kg,
            state.velocity_m_sqrt_kg_per_s,
        ]);
    }
    // Fixed physical energy units for pressure/flow, independent of the clock.
    let volume = core::f64::consts::PI * 0.001_f64.powi(2) * 0.05;
    let pressure_scale = det::sqrt(volume / (gas().density * gas().sound_speed.powi(2)));
    result.extend(
        coupled
            .load()
            .state()
            .pressure()
            .iter()
            .map(|p| p * pressure_scale),
    );
    result
}

#[test]
fn three_coupling_clock_rungs_converge_with_the_same_spatial_duct() {
    // Only the coupling clock changes. This original DG owner retains angular
    // frequencies approximately [10769, 37824, 47649] rad/s. The impulsively
    // established neck flow excites all of them. Coarse 24 -> 48 kHz refinement
    // was nonmonotone: endpoint norm errors 2.293e-4 -> 2.513e-4 against 384 kHz.
    // Resolve the highest spatial mode with omega*dt <= 0.50, 0.25, 0.13 here;
    // the 1536 kHz numerical reference has omega*dt < 0.032. This is convergence
    // of this fixture, not a general clock certificate or continuum solution.
    let reference = trajectory(1536000);
    let errors: Vec<_> = [96000, 192000, 384000]
        .into_iter()
        .map(|rate| {
            let values = trajectory(rate);
            det::sqrt(
                values
                    .iter()
                    .zip(&reference)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum(),
            )
        })
        .collect();
    assert!(errors[0] > 0.0);
    assert!(errors[1] < 0.4 * errors[0], "clock errors: {errors:?}");
    assert!(errors[2] < 0.4 * errors[1], "clock errors: {errors:?}");
}

#[test]
fn late_pressure_refusal_cancellation_and_invalid_input_publish_neither_owner() {
    let gate = CancelGate::new_clock_free();
    let mut coupled = cavity(48000, false, 1e-12)
        .with_downstream_duct(load(48000, false, 0.001), &gate)
        .unwrap();
    let original = coupled.cavity().observe().unwrap();
    let duct = coupled.load().state().clone();
    // Observer pressure is zero. This physical inlet-pressure refusal happens
    // AFTER both original candidates have propagated, before either publishes.
    assert!(matches!(
        coupled.step(&[0.0], 0.0, &gate),
        Err(SharedCavityDgError::Modal(ModalCouplingError::Budget {
            what: "DG mean inlet pressure",
            ..
        }))
    ));
    assert_eq!(coupled.cavity().observe().unwrap(), original);
    assert_eq!(coupled.load().state(), &duct);
    assert!(coupled.last_frame().is_none());
    let mut coupled = cavity(48000, true, 1e6)
        .with_downstream_duct(load(48000, true, 0.001), &gate)
        .unwrap();
    coupled.step(&[0.001], 0.1, &gate).unwrap();
    let frame = coupled.last_frame().unwrap().clone();
    let states: Vec<_> = coupled
        .cavity()
        .network()
        .components()
        .iter()
        .map(|m| m.states().to_vec())
        .collect();
    let memory = coupled.cavity().network().material_memory().to_vec();
    let duct = coupled.load().state().clone();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        coupled.step(&[0.0], 0.0, &cancelled),
        Err(SharedCavityDgError::Modal(ModalCouplingError::Cancelled))
    ));
    assert!(coupled.step(&[0.0], f64::NAN, &gate).is_err());
    assert!(coupled.step(&[0.0], 1e30, &gate).is_err());
    assert_eq!(coupled.last_frame().unwrap(), &frame);
    assert_eq!(coupled.load().state(), &duct);
    assert_eq!(coupled.cavity().network().material_memory(), memory);
    for (model, expected) in coupled.cavity().network().components().iter().zip(states) {
        assert_eq!(model.states(), expected);
    }
}

#[test]
fn connection_rejects_unmodeled_area_change_clock_change_and_duplicate_exterior() {
    let gate = CancelGate::new_clock_free();
    assert!(
        cavity(48000, false, 1e6)
            .with_downstream_duct(load(48000, false, 0.002), &gate)
            .is_err()
    );
    assert!(
        cavity(48000, false, 1e6)
            .with_downstream_duct(load(96000, false, 0.001), &gate)
            .is_err()
    );
    let radiating = cavity(48000, false, 1e6)
        .with_exterior_radiation(
            &RadiationModel {
                ports: 1,
                poles: vec![crate::render::plate::impact::radiation::Pole {
                    omega: 700.0,
                    zeta: 0.1,
                    coupling: vec![0.0],
                }],
            },
            &[1],
            1,
            &gate,
        )
        .unwrap();
    assert!(
        radiating
            .with_downstream_duct(load(48000, false, 0.001), &gate)
            .is_err()
    );
    let mut advanced = cavity(48000, false, 1e6);
    advanced.step(&[0.0], 0.0, &gate).unwrap();
    assert!(
        advanced
            .with_downstream_duct(load(48000, false, 0.001), &gate)
            .is_err()
    );
}
