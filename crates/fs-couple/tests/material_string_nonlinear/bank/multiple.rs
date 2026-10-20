//! Joint prescribed bows, original material stretching and full network work.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    MaterialBowContact, MultiBowStringBank,
};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::contact_load::BowContactLoad;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::stops::{MaterialStringStop, StringStopSide};
use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirection;
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::{
    FrictionObstacleConfig, motion::ObstacleMotionLimits,
};
use fs_dcontact::Obstacle;

mod stations;

fn joint() -> MultiFrictionConfig {
    MultiFrictionConfig { max_ports: 2, max_sweeps: 96, max_setup_terms: 1_048_576 }
}
fn add_second(bank: ReciprocalStringBank, gate: &CancelGate) -> MultiBowStringBank {
    let friction = bank.source().options().friction;
    bank.with_additional_bows(vec![MaterialBowContact {
        location: StringBowLocation { string: 1, station_fraction: 0.31 },
        law: FrictionLaw::Coulomb { static_mu: 0.7, kinetic_mu: 0.3 },
        friction,
    }], joint(), gate).unwrap()
}
fn second() -> FrictionDrive {
    FrictionDrive { speed_m_s: -0.025, normal_force_n: 0.12 }
}

#[test]
fn both_bows_share_source_strain_body_torsion_fluid_and_actual_work() {
    with_cx(|cx, gate| {
        let original = physical_bank(cx, gate, 200e9, gas(), 2)
            .with_geometric_stretching(0.2, nonlinear_config(), gate).unwrap();
        let modes: Vec<_> = original.system().components().iter()
            .map(|c| c.modes().to_vec()).collect();
        let weights = original.system().acceleration_pressure_weights().unwrap().to_vec();
        let mut bank = add_second(original, gate);
        assert_eq!(bank.system().components().iter().map(|c| c.modes().to_vec()).collect::<Vec<_>>(), modes);
        assert_eq!(bank.system().acceleration_pressure_weights(), Some(weights.as_slice()));
        assert_eq!(bank.source().bank().geometric_stretching().unwrap().potentials().len(), 2);
        let mut external = vec![0.0; bank.system().mode_count()];
        let (mut work, mut loss, mut maximum_stretch) = (0.0, 0.0, 0.0_f64);
        for k in 0..320 {
            external[0] = if k < 96 { 0.025 } else { 0.0 };
            let drives = [if k < 192 { active() } else { released() },
                if (32..224).contains(&k) { second() } else { released() }];
            let frame = bank.step_under_gate(&external, &drives, gate).unwrap().clone();
            let nonlinear = frame.nonlinear.unwrap();
            assert!(nonlinear.reaction_residual <= nonlinear.reaction_tolerance);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert_eq!(frame.ports.len(), 2);
            for (port, drive) in frame.ports.iter().zip(drives) {
                let supply = port.traction_n * (drive.speed_m_s / 48_000.0);
                assert!((port.drive_work_j - supply).abs() < 1e-18);
                assert!(port.dissipation_j >= 0.0);
                if drive.normal_force_n == 0.0 { assert_eq!(port.traction_n, 0.0); }
            }
            assert_eq!(frame.drive_work_j, frame.ports.iter().map(|p| p.drive_work_j).sum::<f64>());
            work += frame.external_work_j + frame.drive_work_j;
            loss += frame.network_dissipation_j + frame.friction_dissipation_j;
            assert!((bank.system().total_energy_j().unwrap() + loss - work).abs() < 1e-9);
            let mut geometric_energy = 0.0;
            for string in 0..2 {
                let card = &bank.cards()[string];
                let q = bank.system().components()[string].states();
                // Independent continuous integral of both physical sine fields.
                // Fluid-added mass, but not torsional polar inertia, normalizes q.
                let mass = card.linear_density_kg_m * card.length_m * 0.5;
                let strain = 0.25 * (core::f64::consts::PI / card.length_m).powi(2)
                    * (q[0].displacement_m_sqrt_kg.powi(2) + q[1].displacement_m_sqrt_kg.powi(2)) / mass;
                let observation = bank.stretching_observation(string).unwrap();
                let ea = 200e9 * core::f64::consts::PI * 0.0003_f64.powi(2);
                assert!((observation.additional_strain - strain).abs() < 1e-17);
                assert!((observation.tension_n - card.tension_n - ea * strain).abs() < 1e-11);
                geometric_energy += 0.5 * ea * card.length_m * strain * strain;
            }
            assert!((nonlinear.nonlinear_energy_j - geometric_energy).abs() < 1e-15);
            maximum_stretch = maximum_stretch.max(geometric_energy);
        }
        assert!(maximum_stretch > 0.0);
        assert!(loss > 0.0);
        for string in 0..2 {
            assert_ne!(bank.system().components()[string].states()[2].velocity_m_sqrt_kg_per_s, 0.0,
                "a material torsional coordinate was disconnected from its bow");
        }
        assert!(bank.system().components()[2].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        assert!(bank.system().components()[3..].iter().flat_map(|c| c.states())
            .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0), "causal gas memory was dropped");
    });
}

#[test]
fn empty_additions_preserve_the_exact_nonlinear_single_port_history() {
    with_cx(|cx, gate| {
        let build = || physical_bank(cx, gate, 200e9, gas(), 1)
            .with_geometric_stretching(0.2, nonlinear_config(), gate).unwrap();
        let mut scalar = build();
        let mut original = build();
        let zero = vec![0.0; scalar.system().mode_count()];
        for _ in 0..12 {
            assert_eq!(scalar.step(&zero, active()).unwrap(), original.step(&zero, active()).unwrap());
        }
        let mut forwarded = original.with_additional_bows(vec![], joint(), gate).unwrap();
        assert_eq!(forwarded.system().samples_rendered(), 12);
        for k in 0..128 {
            let drive = if k < 64 { active() } else { released() };
            let expected = *scalar.step(&zero, drive).unwrap();
            let actual = forwarded.step(&zero, &[drive]).unwrap().clone();
            let diagnostic = scalar.system().nonlinear_frame().unwrap();
            assert_eq!(actual.sample, expected.sample);
            assert_eq!(actual.observer_pressure_pa.to_bits(), expected.observer_pressure_pa.to_bits());
            assert_eq!(actual.network_energy_j.to_bits(), expected.network_energy_j.to_bits());
            assert_eq!(actual.drive_work_j.to_bits(), expected.drive_work_j.to_bits());
            assert_eq!(actual.nonlinear.unwrap().nonlinear_energy_j, diagnostic.nonlinear_energy_j);
            assert_eq!(actual.nonlinear.unwrap().nonlinear_iterations, diagnostic.nonlinear_iterations);
            for (a, b) in forwarded.system().components().iter().zip(scalar.system().components()) {
                assert_eq!(a.states(), b.states());
            }
        }
    });
}

#[test]
fn joint_nonlinear_refusal_restores_potentials_coordinates_reports_and_clock() {
    with_cx(|cx, gate| {
        for (slope, config) in [
            (0.2, NonlinearModalConfig { max_iterations: 1,
                force_absolute_tolerance: 1e-30, force_relative_tolerance: 1e-15,
                ..nonlinear_config() }),
            (1e-15, nonlinear_config()),
        ] {
            let mut bank = add_second(physical_bank(cx, gate, 200e9, gas(), 2)
                .with_geometric_stretching(slope, config, gate).unwrap(), gate);
            let zero = vec![0.0; bank.system().mode_count()];
            // Establish a complete accepted report before the refused active trial.
            let accepted = bank.step(&zero, &[released(), released()]).unwrap().clone();
            let states: Vec<_> = bank.system().components().iter().map(|c| c.states().to_vec()).collect();
            assert!(bank.step_under_gate(&zero, &[active(), second()], gate).is_err());
            assert_eq!(bank.system().samples_rendered(), 1);
            assert_eq!(bank.system().last_frame(), Some(&accepted));
            assert_eq!(bank.system().components().iter().map(|c| c.states().to_vec()).collect::<Vec<_>>(), states);
            // A valid retry still has the original constitutive participant.
            let retry = bank.step(&zero, &[released(), released()]).unwrap();
            assert_eq!(retry.sample, 2);
            assert!(retry.nonlinear.is_some());
            assert_eq!(retry.network_energy_j, 0.0);
        }
    });
}

#[test]
fn material_stretching_normal_bow_reaction_and_translating_support_share_one_bank() {
    with_cx(|cx, gate| {
        let build = |young| {
            let original = physical_bank(cx, gate, young, gas(), 1)
                .with_geometric_stretching(0.2, nonlinear_config(), gate).unwrap();
            let modes: Vec<_> = original.system().components().iter().map(|c| c.modes().to_vec()).collect();
            let bank = original.with_additional_bows(vec![], joint(), gate).unwrap()
                .with_string_stops(vec![MaterialStringStop {
                    location: StringBowLocation { string: 0, station_fraction: 0.2 },
                    // Perpendicular to t=(.6,.8), with n x t=+1 for the retained twist.
                    side: StringStopSide::Along(TransverseDirection::try_new([0.8, -0.6]).unwrap()),
                    body: None,
                    law: Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], 5000.0, 1.0,
                        "synthetic prescribed bow surface, not measured hair compliance".into()).unwrap()
                        .with_internal_loss(0.2).unwrap(),
                    limits: ModalContactConfig { max_iterations: 96, maximum_force_n: 10.0,
                        maximum_penetration_m: 0.001, force_absolute_tolerance_n: 1e-11,
                        force_relative_tolerance: 1e-10 },
                }], FrictionObstacleConfig { max_contacts: 1, max_setup_terms: 1_048_576 }, gate).unwrap()
                .with_bow_contact_loads(vec![BowContactLoad { bow: 0, stop: 0 }], gate).unwrap()
                .with_stop_motion(ObstacleMotionLimits { maximum_abs_offset_m: 0.001,
                    maximum_speed_m_s: 0.1 }, gate).unwrap();
            assert_eq!(bank.system().components().iter().map(|c| c.modes().to_vec()).collect::<Vec<_>>(), modes);
            assert_eq!(bank.source().bank().cylinder_gas_loading().unwrap().strings().len(), 1);
            bank
        };
        let mut runs = [build(200e9), build(50e9)];
        let zero = vec![0.0; runs[0].system().mode_count()];
        let drive = [FrictionDrive { speed_m_s: 0.04, normal_force_n: 0.0 }];
        let mut work = [0.0; 2]; let mut loss = [0.0; 2];
        let mut previous = 0.0;
        let (mut active, mut separated, mut recontact) = (0, 0, false);
        let mut max_stretch = 0.0_f64; let mut reaction_difference = 0.0_f64;
        for sample in 0..320 {
            let offset = if sample < 64 { -3e-5 * (sample+1) as f64 / 64.0 }
                else if sample < 128 { -3e-5 }
                else if sample < 192 { -3e-5 + 6e-5 * (sample-127) as f64 / 64.0 }
                else if sample < 256 { 3e-5 }
                else { 3e-5 - 6e-5 * (sample-255) as f64 / 64.0 };
            let mut normal = [0.0; 2];
            for (index, bank) in runs.iter_mut().enumerate() {
                let frame = bank.step_with_stop_offsets_under_gate(&zero, &drive, &[offset], gate).unwrap().clone();
                let contact = bank.system().last_obstacle_frames().unwrap()[0];
                normal[index] = contact.normal_force_n;
                assert_eq!(frame.ports[0].normal_force_n, contact.normal_force_n);
                assert_eq!(bank.system().obstacle_offsets_m(), Some([offset].as_slice()));
                assert_eq!(contact.actuator_work_j, -contact.normal_force_n * (offset-previous));
                assert_eq!(frame.obstacle_work_j, contact.actuator_work_j);
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
                work[index] += frame.external_work_j + frame.drive_work_j + frame.obstacle_work_j;
                loss[index] += frame.network_dissipation_j + frame.friction_dissipation_j + frame.obstacle_dissipation_j;
                assert!((bank.system().total_energy_j().unwrap() + loss[index] - work[index]).abs() < 1e-9);
                assert!((bank.system().total_energy_j().unwrap() - frame.network_energy_j - frame.obstacle_energy_j).abs() < 1e-13);
                let card = &bank.cards()[0]; let q = bank.system().components()[0].states();
                let mass = card.linear_density_kg_m * card.length_m * 0.5;
                let strain = 0.25 * (core::f64::consts::PI / card.length_m).powi(2)
                    * (q[0].displacement_m_sqrt_kg.powi(2) + q[1].displacement_m_sqrt_kg.powi(2)) / mass;
                let ea = [200e9, 50e9][index] * core::f64::consts::PI * 0.0003_f64.powi(2);
                let observed = bank.stretching_observation(0).unwrap();
                assert!((observed.additional_strain - strain).abs() < 1e-17);
                assert!((observed.tension_n - card.tension_n - ea*strain).abs() < 1e-11);
                assert!((frame.nonlinear.unwrap().nonlinear_energy_j - 0.5*ea*card.length_m*strain*strain).abs() < 1e-15);
                max_stretch = max_stretch.max(frame.nonlinear.unwrap().nonlinear_energy_j);
                if index == 0 {
                    if contact.normal_force_n == 0.0 {
                        separated += 1; assert_eq!(frame.ports[0].traction_n, 0.0);
                    } else { active += 1; recontact |= sample >= 256; }
                }
            }
            reaction_difference += (normal[0]-normal[1]).abs();
            previous = offset;
        }
        assert!(active > 0 && separated > 0 && recontact && max_stretch > 0.0);
        assert!(reaction_difference > 1e-10, "material change did not reach the actual contact load");
        for bank in &runs {
            assert_ne!(bank.system().components()[0].states()[2].velocity_m_sqrt_kg_per_s, 0.0);
            assert!(bank.system().components()[1].states().iter().any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
            assert!(bank.system().components()[2..].iter().flat_map(|c| c.states())
                .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0));
        }
    });
}
