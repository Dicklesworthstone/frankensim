//! Joint prescribed bows, original material stretching and full network work.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    MaterialBowContact, MultiBowStringBank,
};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig;

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
