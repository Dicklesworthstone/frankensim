//! Omitted string storage is in series with the actual material film.
//! All source coefficients here are manufactured analytical fixtures.
use super::*;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionModalSystem, ModalFrictionPort,
    multiple::{MultiFrictionModalSystem, surfaces::geometry::FrictionGeometry},
};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment};
use fs_math::c64::C64;

const TAIL: [f64; 4] = [2e-5, 5e-6, 5e-6, 4e-5];

fn quadratic(force: [f64; 2], g: &[f64]) -> f64 {
    0.5 * (g[0] * force[0].powi(2) + 2.0 * g[1] * force[0] * force[1] + g[3] * force[1].powi(2))
}

fn analytical_system(
    sources: Vec<ResolvedShearLayer>,
    residual_first: bool,
    gate: &CancelGate,
) -> MultiFrictionModalSystem {
    let models = [180.0, 260.0]
        .into_iter()
        .map(|hz| {
            ModalAcousticTimeModel::try_new(
                48_000,
                vec![ModalAcousticMode {
                    angular_frequency_rad_s: core::f64::consts::TAU * hz,
                    damping_ratio: 0.0,
                    pressure_per_modal_velocity: C64::ZERO,
                }],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap()
        })
        .collect();
    let network = CoupledModalSystem::new(models, vec![], options().network, gate).unwrap();
    let port = |component| ModalFrictionPort {
        left: ModalAttachment {
            component,
            shapes: vec![100.0],
        },
        right: None,
        law: law(),
    };
    let maps = (0..2)
        .map(|i| FrictionGeometry {
            left: port(i).left,
            right: None,
        })
        .collect();
    let mut system = FrictionModalSystem::new(network, port(0), options().friction, gate)
        .unwrap()
        .with_additional_ports(vec![(port(1), options().friction)], joint_limits(), gate)
        .unwrap()
        .with_series_compliant_drives(vec![(0, 4e4, 0.4), (1, 7e4, 0.8)], gate)
        .unwrap();
    let bindings = sources.into_iter().enumerate().collect();
    if residual_first {
        system = system
            .with_friction_residual_ports(maps, TAIL.to_vec(), gate)
            .unwrap()
            .with_material_shear_layers(bindings, gate)
            .unwrap();
    } else {
        system = system
            .with_material_shear_layers(bindings, gate)
            .unwrap()
            .with_friction_residual_ports(maps, TAIL.to_vec(), gate)
            .unwrap();
    }
    system
}

#[test]
fn residual_film_and_series_drive_match_an_independent_two_port_solution() {
    let gate = CancelGate::new_clock_free();
    for plastic in [false, true] {
        for memory in [false, true] {
            for residual_first in [false, true] {
                let sources = [(8e6, 10.0), (16e6, 40.0)]
                    .into_iter()
                    .map(|(g, eta)| {
                        layer_with_memory(g, if plastic { 0.0 } else { 2e6 }, eta, 1e6, memory)
                    })
                    .collect();
                let mut system = analytical_system(sources, residual_first, &gate);
                let dt = system.sample_period_s();
                let k = [8e5, 16e5]; // G A/h from the original source geometry.
                let resistance = [1.0, 4.0]; // eta A/h, not a friction coefficient.
                let hair_k = [4e4, 7e4];
                let hair_c = [0.4, 0.8];
                let mut a = TAIL;
                let mut material_compliance = [0.0; 2];
                let mut drive_compliance = [0.0; 2];
                for i in 0..2 {
                    // Exact undamped oscillator displacement under a held force,
                    // derived independently of the network's response compiler.
                    let omega = core::f64::consts::TAU * [180.0, 260.0][i];
                    let retained = 10000.0 * (1.0 - (omega * dt).cos()) / omega.powi(2);
                    let midpoint_k = 0.5 * k[i]
                        + if memory {
                            0.5 * k[i] / (2.0 + dt / 0.0002)
                        } else {
                            0.0
                        };
                    material_compliance[i] = 1.0 / midpoint_k;
                    drive_compliance[i] = 1.0 / (0.5 * hair_k[i] + hair_c[i] / dt);
                    a[3 * i] += retained
                        + material_compliance[i]
                        + drive_compliance[i]
                        + if plastic { dt / resistance[i] } else { 0.0 };
                }
                let distance = 0.02 * dt;
                let determinant = a[0] * a[3] - a[1] * a[2];
                let force = [
                    -a[1] * distance / determinant,
                    a[0] * distance / determinant,
                ];
                let frame = system
                    .step(&[0.0, 0.0], &[drive(0.0, 1.0), drive(0.02, 1.0)])
                    .unwrap()
                    .clone();
                let (mut film_energy, mut film_heat, mut drive_energy, mut drive_heat) =
                    (0.0, 0.0, 0.0, 0.0);
                for i in 0..2 {
                    let port = frame.ports[i];
                    let film = port.shear_layer.unwrap();
                    let elastic = material_compliance[i] * force[i];
                    let plastic_slip = if plastic {
                        dt * force[i] / resistance[i]
                    } else {
                        0.0
                    };
                    let arm_slip = if memory {
                        elastic / (1.0 + dt / (2.0 * 0.0002))
                    } else {
                        0.0
                    };
                    let stored = 0.5 * k[i] * elastic.powi(2) + 0.25 * k[i] * arm_slip.powi(2);
                    let relaxation_heat = dt / 0.0002 * 0.5 * k[i] * (0.5 * arm_slip).powi(2);
                    let heat = resistance[i] * plastic_slip.powi(2) / dt + relaxation_heat;
                    let stretch = drive_compliance[i] * force[i];
                    assert!((port.traction_n - force[i]).abs() < 2e-10);
                    assert!((film.elastic_slip_m - elastic).abs() < 2e-13);
                    assert!((film.plastic_slip_m - plastic_slip).abs() < 2e-13);
                    assert!((port.slip_distance_m - elastic - plastic_slip).abs() < 2e-13);
                    assert!((film.stored_energy_j - stored).abs() < 1e-14);
                    assert!((film.dissipated_j - heat).abs() < 1e-14);
                    assert!((film.relaxation_dissipated_j - relaxation_heat).abs() < 1e-14);
                    assert_eq!(film.plastic_flow, plastic);
                    film_energy += stored;
                    film_heat += heat;
                    drive_energy += 0.5 * hair_k[i] * stretch.powi(2);
                    drive_heat += hair_c[i] * stretch.powi(2) / dt;
                }
                let tail = quadratic(force, &TAIL);
                assert!((frame.friction_residual_energy_j - tail).abs() < 1e-14);
                assert_eq!(
                    frame.friction_residual_loss_j,
                    frame.friction_residual_energy_j
                );
                assert!((frame.friction_dissipation_j - film_heat - drive_heat).abs() < 1e-14);
                assert!(
                    (system.total_energy_j().unwrap()
                        - frame.network_energy_j
                        - film_energy
                        - drive_energy
                        - tail)
                        .abs()
                        < 1e-14
                );
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
        }
    }
}

fn instrument(
    cx: &Cx<'_>,
    gate: &CancelGate,
    residual_first: bool,
    g: f64,
    viscosity: f64,
) -> MultiBowStringBank {
    let bank = base(cx, gate, 7800.0)
        .with_residual_flexibility(gate)
        .unwrap()
        .with_additional_bows(vec![extra()], joint_limits(), gate)
        .unwrap();
    let bindings = vec![(0, layer_with_memory(g, 2000.0, viscosity, 1e6, true))];
    if residual_first {
        bank.with_bow_residual_flexibility(gate)
            .unwrap()
            .with_bow_shear_layers(bindings, gate)
            .unwrap()
    } else {
        bank.with_bow_shear_layers(bindings, gate)
            .unwrap()
            .with_bow_residual_flexibility(gate)
            .unwrap()
    }
}

#[test]
fn residual_film_history_reaches_original_pressure_and_survives_rejected_opening() {
    with_cx(|cx, gate| {
        let mut a = instrument(cx, gate, true, 8e6, 10.0);
        let mut b = instrument(cx, gate, false, 8e6, 10.0);
        let mut changed = instrument(cx, gate, true, 16e6, 40.0);
        let zero = vec![0.0; a.system().mode_count()];
        let g = a.system().friction_residual_compliance().unwrap().to_vec();
        let (mut work, mut loss, mut numerical, mut pressure_difference) = (0.0, 0.0, 0.0, 0.0);
        let mut old = [0.0; 2];
        for sample in 0..128 {
            let drivers = [
                drive(if sample < 64 { 0.001 } else { -0.001 }, 0.1),
                drive(0.0, 0.0),
            ];
            let frame = a.step(&zero, &drivers).unwrap().clone();
            assert_eq!(&frame, b.step(&zero, &drivers).unwrap());
            let other = changed.step(&zero, &drivers).unwrap();
            pressure_difference += (frame.observer_pressure_pa - other.observer_pressure_pa).abs();
            let force = [frame.ports[0].traction_n, frame.ports[1].traction_n];
            let increment = [force[0] - old[0], force[1] - old[1]];
            assert!((frame.friction_residual_energy_j - quadratic(force, &g)).abs() < 1e-14);
            assert!((frame.friction_residual_loss_j - quadratic(increment, &g)).abs() < 1e-14);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            work += frame.drive_work_j;
            loss += frame.network_dissipation_j
                + frame.friction_dissipation_j
                + frame.friction_residual_loss_j;
            numerical += frame.friction_residual_loss_j;
            assert!((a.system().total_energy_j().unwrap() + loss - work).abs() < 2e-10);
            old = force;
        }
        assert!(numerical > 0.0 && pressure_difference > 1e-15);
        let accepted = a.system().last_frame().unwrap().clone();
        let history = a.system().material_shear_state(0).unwrap().clone();
        let energy = a.system().total_energy_j().unwrap();
        let original = states(&a);
        assert!(a.step(&zero, &[drive(0.0, 0.0); 2]).is_err());
        assert_eq!(a.system().last_frame(), Some(&accepted));
        assert_eq!(a.system().material_shear_state(0), Some(&history));
        assert_eq!(a.system().total_energy_j().unwrap(), energy);
        assert_eq!(states(&a), original);
        let retry = [drive(0.0, 0.1), drive(0.0, 0.0)];
        assert_eq!(
            a.step(&zero, &retry).unwrap(),
            b.step(&zero, &retry).unwrap()
        );
    });
}
