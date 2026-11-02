//! G1: a condensed support and both physical surfaces share one virtual work.
//! Frequencies and layer coefficients are manufactured analytical fixtures.
use super::*;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionModalSystem, ModalFrictionPort, multiple::surfaces::geometry::FrictionGeometry,
};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment};
use fs_math::c64::C64;

#[test]
fn receiving_residual_refuses_incomplete_or_foreign_retained_corrections() {
    let gate = CancelGate::new_clock_free();
    for case in 0..3 {
        let models = (0..3)
            .map(|_| {
                ModalAcousticTimeModel::try_new(
                    48_000,
                    vec![ModalAcousticMode {
                        angular_frequency_rad_s: 1000.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    }],
                    ModalAcousticTimeBudget::audible_reference(),
                )
                .unwrap()
            })
            .collect();
        let network = CoupledModalSystem::new(models, vec![], options().network, &gate).unwrap();
        let port = ModalFrictionPort {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0],
            },
            right: Some(ModalAttachment {
                component: 1,
                shapes: vec![1.0],
            }),
            law: law(),
        };
        let geometry = vec![FrictionGeometry {
            left: port.left.clone(),
            right: port.right.clone(),
        }];
        let mut corrections = vec![ModalAttachment {
            component: 2,
            shapes: vec![0.25],
        }];
        match case {
            0 => corrections.clear(),
            1 => corrections[0].component = 3,
            _ => corrections[0].shapes.push(0.25),
        }
        let system = FrictionModalSystem::new(network, port, options().friction, &gate)
            .unwrap()
            .with_additional_ports(vec![], joint_limits(), &gate)
            .unwrap();
        assert!(
            system
                .with_friction_residual_ports_and_corrections(
                    geometry,
                    corrections,
                    vec![1e-5],
                    &gate
                )
                .is_err()
        );
    }
}

#[test]
fn receiving_residual_film_matches_independent_three_participant_reactions_and_work() {
    let gate = CancelGate::new_clock_free();
    // Two string coordinates, one shared bridge body, two receiving fibers.
    // The shared body contributes a signed cross-port compliance; omitting
    // its force or replacing the fiber with it changes this analytic solution.
    let hz = [180.0, 260.0, 90.0, 330.0, 410.0];
    let columns = [
        [100.0, 0.0, 35.0, -60.0, 0.0],
        [0.0, 120.0, -20.0, 0.0, -75.0],
    ];
    let tail = [2e-5, 5e-6, 5e-6, 4e-5];
    for film_first in [false, true] {
        let models = hz
            .iter()
            .map(|frequency| {
                ModalAcousticTimeModel::try_new(
                    48_000,
                    vec![ModalAcousticMode {
                        angular_frequency_rad_s: core::f64::consts::TAU * frequency,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    }],
                    ModalAcousticTimeBudget::audible_reference(),
                )
                .unwrap()
            })
            .collect();
        let network = CoupledModalSystem::new(models, vec![], options().network, &gate).unwrap();
        let port = |i: usize| ModalFrictionPort {
            left: ModalAttachment {
                component: i,
                shapes: vec![columns[i][i]],
            },
            right: Some(ModalAttachment {
                component: 3 + i,
                shapes: vec![-columns[i][3 + i]],
            }),
            law: law(),
        };
        let mut system = FrictionModalSystem::new(network, port(0), options().friction, &gate)
            .unwrap()
            .with_additional_ports(vec![(port(1), options().friction)], joint_limits(), &gate)
            .unwrap();
        let maps = (0..2)
            .map(|i| FrictionGeometry {
                left: port(i).left,
                right: port(i).right,
            })
            .collect();
        let corrections = (0..2)
            .map(|i| ModalAttachment {
                component: 2,
                shapes: vec![columns[i][2]],
            })
            .collect();
        let layers = vec![
            (0, layer(8e6, 0.0, 10.0, 1e6)),
            (1, layer(16e6, 0.0, 40.0, 1e6)),
        ];
        if film_first {
            system = system
                .with_material_shear_layers(layers.clone(), &gate)
                .unwrap();
        }
        system = system
            .with_friction_residual_ports_and_corrections(maps, corrections, tail.to_vec(), &gate)
            .unwrap();
        if !film_first {
            system = system.with_material_shear_layers(layers, &gate).unwrap();
        }
        let dt = system.sample_period_s();
        let modal_response: Vec<_> = hz
            .iter()
            .map(|frequency| {
                let omega = core::f64::consts::TAU * frequency;
                (1.0 - (omega * dt).cos()) / omega.powi(2)
            })
            .collect();
        let stiffness = [8e5, 16e5]; // Source G*A/h.
        let resistance = [1.0, 4.0]; // Source eta*A/h.
        let mut a = tail;
        for i in 0..2 {
            for j in 0..2 {
                a[2 * i + j] += (0..5)
                    .map(|mode| columns[i][mode] * columns[j][mode] * modal_response[mode])
                    .sum::<f64>();
            }
            a[3 * i] += 2.0 / stiffness[i] + dt / resistance[i];
        }
        let speeds = [0.01, 0.02];
        let determinant = a[0] * a[3] - a[1] * a[2];
        let force = [
            dt * (a[3] * speeds[0] - a[1] * speeds[1]) / determinant,
            dt * (a[0] * speeds[1] - a[2] * speeds[0]) / determinant,
        ];
        let frame = system
            .step(&[0.0; 5], &[drive(speeds[0], 1.0), drive(speeds[1], 1.0)])
            .unwrap()
            .clone();
        let mut mechanical_energy = 0.0;
        for mode in 0..5 {
            let load = columns[0][mode] * force[0] + columns[1][mode] * force[1];
            let omega = core::f64::consts::TAU * hz[mode];
            let q = load * modal_response[mode];
            let v = load * (omega * dt).sin() / omega;
            let actual = system.components()[mode].states()[0];
            assert!(
                (actual.displacement_m_sqrt_kg - q).abs() < 2e-14,
                "mode {mode}: reciprocal displacement"
            );
            assert!(
                (actual.velocity_m_sqrt_kg_per_s - v).abs() < 2e-10,
                "mode {mode}: reciprocal force"
            );
            mechanical_energy += 0.5 * (v * v + omega * omega * q * q);
        }
        let (mut film_energy, mut physical_heat, mut supplied) = (0.0, 0.0, 0.0);
        for i in 0..2 {
            assert_eq!(
                system.port(i).unwrap().right,
                port(i).right,
                "the physical fiber receiver was replaced by the condensed body"
            );
            assert!((frame.ports[i].traction_n - force[i]).abs() < 2e-10);
            let elastic = 2.0 * force[i] / stiffness[i];
            let plastic = dt * force[i] / resistance[i];
            let film = frame.ports[i].shear_layer.unwrap();
            assert!((film.elastic_slip_m - elastic).abs() < 2e-13);
            assert!((film.plastic_slip_m - plastic).abs() < 2e-13);
            film_energy += 0.5 * stiffness[i] * elastic * elastic;
            physical_heat += resistance[i] * plastic * plastic / dt;
            supplied += force[i] * speeds[i] * dt;
        }
        let tail_energy = 0.5
            * (tail[0] * force[0] * force[0]
                + 2.0 * tail[1] * force[0] * force[1]
                + tail[3] * force[1] * force[1]);
        assert!((frame.network_energy_j - mechanical_energy).abs() < 1e-14);
        assert!((frame.friction_dissipation_j - physical_heat).abs() < 1e-14);
        assert!((frame.friction_residual_energy_j - tail_energy).abs() < 1e-14);
        assert_eq!(
            frame.friction_residual_loss_j,
            frame.friction_residual_energy_j
        );
        assert!(
            (mechanical_energy + film_energy + physical_heat + 2.0 * tail_energy - supplied).abs()
                < 1e-14
        );
        assert!(
            (system.total_energy_j().unwrap() - mechanical_energy - film_energy - tail_energy)
                .abs()
                < 1e-14
        );
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
}
