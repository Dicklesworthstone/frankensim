//! G1: moving-phase ownership precedes physical endpoint arrest.

use core::f64::consts::PI;
use fs_tribo::sliding_patch::contact::multiple::rolling::*;
use fs_tribo::sliding_patch::contact::multiple::*;
use fs_tribo::sliding_patch::contact::*;

fn law(static_mu: f64, kinetic_mu: f64) -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 1.0,
        contact_radius_m: 0.2,
        static_mu,
        kinetic_mu,
    })
    .unwrap()
}
fn config() -> HertzWrenchConfig {
    HertzWrenchConfig {
        maximum_static_iterations: 96,
        maximum_iterations: 256,
        maximum_evaluations: 1024,
        force_tolerance_n: 2e-10,
        moment_tolerance_nm: 2e-11,
        displacement_tolerance_m: 2e-13,
        energy_tolerance_j: 2e-11,
        relative_matrix_tolerance: 1e-13,
        maximum_force_n: 20.0,
        maximum_spin_moment_nm: 4.0,
    }
}
fn multiple() -> MultipleHertzConfig {
    MultipleHertzConfig {
        maximum_contacts: 2,
        maximum_matrix_terms: 36,
        maximum_sweeps: 128,
        maximum_evaluations: 16384,
        local: config(),
        total_energy_tolerance_j: 4e-11,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tol {tolerance:e}"
    );
}

#[test]
fn high_static_capacity_does_not_capture_moving_translation_before_or_at_arrest() {
    let law = law(5.0, 1.0);
    let original = law.parameters();
    // v0=W=Fk=1 gives the exact physical arrest time 1s. Larger static
    // capacity could incorrectly choose x_mid=0 already at0.4s.
    for dt in [0.75, 1.0] {
        let c = 0.5 * dt * dt;
        let input = HertzWrenchInput {
            free_displacement_m: [dt, 0.0],
            free_spin_angle_rad: 0.0,
            scaled_compliance_m_per_n: [[c, 0.0, 0.0], [0.0, c, 0.0], [0.0, 0.0, c]],
        };
        let reflected = law.advance(&input, config(), || false).unwrap();
        assert_eq!(reflected.mode, HertzWrenchMode::Sticking);
        let physical = law
            .advance_phase(&input, HertzWrenchPhase::Kinetic, config(), || false)
            .unwrap();
        assert_eq!(physical.mode, HertzWrenchMode::Sliding);
        close(physical.force_n[0], -1.0, 2e-15);
        close(physical.slip_displacement_m[0], dt - 0.5 * dt * dt, 2e-15);
        close(
            2.0 * physical.slip_displacement_m[0] / dt - 1.0,
            1.0 - dt,
            2e-15,
        );
        close(physical.dissipated_energy_j, dt - 0.5 * dt * dt, 2e-15);
    }
    assert_eq!(law.parameters(), original);
}

#[test]
fn high_static_spin_and_zero_kinetic_limits_keep_the_original_moving_surface() {
    let law = law(5.0, 1.0);
    let capacity = 3.0 * PI / 16.0;
    let dt = 1.0 / capacity;
    let c = 0.5 * dt * dt;
    let input = HertzWrenchInput {
        free_displacement_m: [0.0; 2],
        free_spin_angle_rad: dt / 0.2,
        scaled_compliance_m_per_n: [[c, 0.0, 0.0], [0.0, c, 0.0], [0.0, 0.0, c]],
    };
    let spin = law
        .advance_phase(&input, HertzWrenchPhase::Kinetic, config(), || false)
        .unwrap();
    assert_eq!(spin.mode, HertzWrenchMode::Sliding);
    close(spin.spin_moment_nm, -0.2 * capacity, 2e-15);
    close(2.0 * 0.2 * spin.spin_angle_rad / dt - 1.0, 0.0, 3e-15);
    close(spin.dissipated_energy_j, 0.5, 3e-15);
    let dry = HertzWrenchContact::new(HertzWrenchParameters {
        kinetic_mu: 0.0,
        ..law.parameters()
    })
    .unwrap();
    let free = dry
        .advance_phase(&input, HertzWrenchPhase::Kinetic, config(), || false)
        .unwrap();
    assert_eq!(free.mode, HertzWrenchMode::Frictionless);
    assert_eq!(free.spin_moment_nm, 0.0);
    assert_eq!(free.dissipated_energy_j, 0.0);
    close(free.spin_angle_rad, input.free_spin_angle_rad, 2e-15);
    assert_eq!(dry.parameters().static_mu, 5.0);
}

#[test]
fn shared_moving_and_stationary_patches_keep_distinct_original_capacities() {
    let laws = [law(10.0, 1.0), law(1.0, 0.1)];
    let mut c = vec![0.0; 36];
    for i in 0..6 {
        c[i * 6 + i] = 0.25;
    }
    c[3] = 0.125;
    c[18] = 0.125;
    let prepared = PreparedHertzContacts::new(2, c, multiple(), || false).unwrap();
    let sites = [
        HertzContactInput {
            law: &laws[0],
            free_displacement_m: [0.5, 0.0],
            free_spin_angle_rad: 0.0,
        },
        HertzContactInput {
            law: &laws[1],
            free_displacement_m: [0.0; 2],
            free_spin_angle_rad: 0.0,
        },
    ];
    let phases = [
        HertzWrenchPhase::Kinetic,
        HertzWrenchPhase::StaticAdmissible,
    ];
    let step = prepared.advance_phases(&sites, &phases, || false).unwrap();
    assert_eq!(step.contacts[0].mode, HertzWrenchMode::Sliding);
    assert_eq!(step.contacts[1].mode, HertzWrenchMode::Sticking);
    close(step.contacts[0].force_n[0], -1.0, 2e-14);
    close(step.contacts[1].force_n[0], 0.5, 2e-14);
    assert!(step.contacts[1].force_n[0] > laws[1].parameters().kinetic_mu);
    close(step.contacts[0].slip_displacement_m[0], 0.3125, 2e-14);
    assert_eq!(step.contacts[1].dissipated_energy_j, 0.0);
    assert_eq!(
        prepared.advance_phases(&sites, &phases, || false).unwrap(),
        step
    );
    assert!(
        prepared
            .advance_phases(&sites, &phases[..1], || false)
            .is_err()
    );
    assert!(matches!(
        prepared.advance_phases(&sites, &phases, || true),
        Err(MultipleHertzError::Cancelled)
    ));
}

#[test]
fn joint_rolling_reduction_preserves_moving_hertz_phase_and_source_parameters() {
    let original = law(5.0, 1.0);
    let mut c = vec![0.0; 25];
    for i in 0..5 {
        c[i * 5 + i] = 0.5;
    }
    let prepared = PreparedHertzRollingContacts::new(
        1,
        c,
        vec![2.0],
        MultipleHertzRollingConfig {
            hertz: multiple(),
            maximum_matrix_terms: 100,
            rolling_force_tolerance_n: 2e-10,
            total_energy_tolerance_j: 4e-11,
        },
        || false,
    )
    .unwrap();
    let sites = [HertzRollingContactInput {
        law: &original,
        free_coordinates_m: [1.0, 0.0, 0.0, 0.4, 0.0],
    }];
    let step = prepared
        .advance_phases(&sites, &[HertzWrenchPhase::Kinetic], || false)
        .unwrap();
    assert_eq!(step.hertz.contacts[0].mode, HertzWrenchMode::Sliding);
    close(step.hertz.contacts[0].force_n[0], -1.0, 2e-15);
    close(step.surface_heat_j, 0.5, 2e-15);
    close(step.rolling[0].scaled_displacement_m[0], 0.2, 2e-15);
    close(step.bulk_heat_j, 0.08, 2e-15);
    assert_eq!(original.parameters().static_mu, 5.0);
    assert!(prepared.advance_phases(&sites, &[], || false).is_err());
}
