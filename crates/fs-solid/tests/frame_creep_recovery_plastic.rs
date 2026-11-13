//! Native member creep, plasticity and resistance recovery share one interval.
//! Independent scalar references use manufactured coefficients, not service data.

use fs_material::IsotropicElastic;
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::plastic::uniaxial::UniaxialNortonJ2;
use fs_material::creep::recovery::CreepRecoveryParameters;
use fs_material::creep::recovery::plastic::RecoveringNortonJ2;
use fs_material::creep::recovery::plastic::uniaxial::UniaxialRecoveringNortonJ2;
use fs_material::creep::{
    CreepTolerances, NortonCreep, NortonCreepParameters, StrainHardeningCreep,
};
use fs_material::plastic::J2Plasticity;
use fs_solid::fiber::material::FiberMaterial;
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::frame2d::creep::{FrameCreepError, FrameCreepTolerances};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel, SolidError};
use std::cell::Cell;

const E: f64 = 1e9;
const STRESS: f64 = 1e6;
const YIELD: f64 = 0.8e6;
const HARDENING: f64 = 2e8;
const RATE: f64 = 1e-5;
const PRIMARY: f64 = 4e-4;
const RECOVERY: f64 = 0.2;
const DT: f64 = 10.0;
const AREA: f64 = 5e-4;
const ARM: f64 = 0.05;
const FORCE: f64 = 2.0 * AREA * STRESS;

fn law(dt: f64, recovery: f64, yield_stress: f64, hardening: f64) -> UniaxialRecoveringNortonJ2 {
    let elastic = IsotropicElastic::new(E, 0.25, 0.03).unwrap();
    let creep = NortonCreep::new(
        elastic.clone(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: 1.0,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 300.0,
            temperature_range_k: [300.0, 500.0],
        },
        300.0,
        dt,
    )
    .unwrap()
    .with_strain_hardening(StrainHardeningCreep {
        reference_strain: PRIMARY,
        exponent: 1.0,
        residual_rate_fraction: 0.0,
    })
    .unwrap();
    let plastic = J2Plasticity::new(elastic, yield_stress, hardening).unwrap();
    UniaxialRecoveringNortonJ2::new(
        RecoveringNortonJ2::new(
            NortonJ2::new(creep, plastic).unwrap(),
            CreepRecoveryParameters {
                reference_rate_s: recovery,
                activation_energy_j_mol: 0.0,
                reference_temperature_k: 300.0,
            },
        )
        .unwrap(),
    )
}

fn frame(
    law: &UniaxialRecoveringNortonJ2,
    mass: f64,
    formulation: FrameMassFormulation,
) -> Frame2d {
    let material = FiberMaterial::authored(
        FiberLaw::RecoveringNortonJ2(law.clone()),
        mass / (2.0 * AREA),
    )
    .unwrap();
    assert!(material.source().is_none());
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    frame.fix(1, [false, true, true]);
    frame
        .add_element(0, 1, 3, mass, &|| {
            SectionModel::Fiber(Section::new(
                [-ARM, ARM]
                    .into_iter()
                    .map(|y| Fiber {
                        y,
                        area: AREA,
                        law: material.law().clone(),
                    })
                    .collect(),
            ))
        })
        .unwrap();
    frame
}

fn controls() -> FrameCreepTolerances {
    FrameCreepTolerances {
        material: CreepTolerances {
            stress_pa: 1e-3,
            energy_j_m3: 1e-5,
        },
        energy_j: 1e-7,
        force_n: 1e-5,
        moment_n_m: 1e-6,
        maximum_fiber_points: 6,
    }
}

fn load(force: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, force, 0.0, 0.0]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

fn creep_increment(stress: f64, resistance: f64, dt: f64) -> f64 {
    let a = PRIMARY * (1.0 + RECOVERY * dt);
    let b = a + resistance;
    let c = dt * RATE * stress.abs() / STRESS * a;
    2.0 * c / (b + (b * b + 4.0 * c).sqrt())
}

#[test]
fn loaded_joint_member_retains_plastic_hardening_while_dwell_renews_creep() {
    let mut accepted = frame(
        &law(DT, RECOVERY, YIELD, HARDENING),
        1.0,
        FrameMassFormulation::Lumped,
    );
    let plastic = (STRESS - YIELD) / HARDENING;
    let (mut creep, mut resistance) = (0.0, 0.0);
    for index in 0..3 {
        let before = accepted
            .creep_recovery_plastic_fiber(0, 0, 0)
            .unwrap()
            .1
            .clone();
        let dc = creep_increment(STRESS, resistance, DT);
        let step = accepted
            .solve_creep_step(&load(FORCE), controls(), || false)
            .unwrap();
        assert_eq!(
            accepted.creep_recovery_plastic_fiber(0, 0, 0).unwrap().1,
            &before
        );
        creep += dc;
        resistance = (resistance + dc) / (1.0 + RECOVERY * DT);
        close(
            step.frame.displacements()[3],
            STRESS / E + creep + plastic,
            2e-11,
        );
        close(step.frame.element_forces(0)[0], FORCE, 1e-6);
        close(step.energy.creep_dissipation_j, FORCE * dc, 1e-8);
        close(
            step.energy.plastic_dissipation_j,
            if index == 0 {
                2.0 * AREA * YIELD * plastic
            } else {
                0.0
            },
            1e-8,
        );
        close(
            step.energy.final_stored_j,
            AREA * (STRESS * STRESS / E + HARDENING * plastic * plastic),
            1e-8,
        );
        close(step.energy.balance_residual_j, 0.0, controls().energy_j);
        for section in 0..3 {
            for fiber in 0..2 {
                let state = step
                    .frame
                    .creep_recovery_plastic_fiber(0, section, fiber)
                    .unwrap()
                    .1;
                close(state.creep_strain(), creep, 2e-11);
                close(state.equivalent_creep_strain(), creep, 2e-11);
                close(state.plastic_strain(), plastic, 2e-11);
                close(state.equivalent_plastic_strain(), plastic, 2e-11);
                close(state.resistance_strain(), resistance, 2e-12);
            }
        }
        accepted = step.frame;
    }
    let mut held = accepted.clone();
    held.fix(1, [true; 3]);
    let relaxed = held
        .solve_creep_step(&[0.0; 6], controls(), || false)
        .unwrap();
    assert_eq!(relaxed.frame.displacements(), held.displacements());
    assert!(relaxed.frame.element_forces(0)[0] < FORCE);
    assert!(relaxed.frame.element_forces(0)[0] > 0.0);
    close(
        relaxed
            .frame
            .creep_recovery_plastic_fiber(0, 0, 0)
            .unwrap()
            .1
            .plastic_strain(),
        plastic,
        2e-11,
    );
    assert!(relaxed.energy.creep_dissipation_j > 0.0);
    close(relaxed.energy.plastic_dissipation_j, 0.0, 1e-8);
    close(relaxed.energy.balance_residual_j, 0.0, controls().energy_j);

    let released = accepted
        .solve_creep_step(&load(0.0), controls(), || false)
        .unwrap();
    close(released.frame.displacements()[3], creep + plastic, 2e-11);
    close(released.energy.physical_dissipation_j, 0.0, 1e-8);
    let immediate = released
        .frame
        .solve_creep_step(&load(FORCE), controls(), || false)
        .unwrap();
    let original = released
        .frame
        .creep_recovery_plastic_fiber(0, 0, 0)
        .unwrap()
        .1
        .clone();
    let mut rested = released.frame;
    for _ in 0..4 {
        let before = rested
            .creep_recovery_plastic_fiber(0, 0, 0)
            .unwrap()
            .1
            .clone();
        let step = rested
            .solve_creep_step(&load(0.0), controls(), || false)
            .unwrap();
        let after = step.frame.creep_recovery_plastic_fiber(0, 0, 0).unwrap().1;
        close(after.creep_strain(), original.creep_strain(), 2e-11);
        close(
            after.equivalent_creep_strain(),
            original.equivalent_creep_strain(),
            2e-11,
        );
        close(after.plastic_strain(), original.plastic_strain(), 2e-11);
        close(
            after.equivalent_plastic_strain(),
            original.equivalent_plastic_strain(),
            2e-11,
        );
        close(
            after.resistance_strain(),
            before.resistance_strain() / 3.0,
            2e-12,
        );
        close(step.energy.physical_dissipation_j, 0.0, 1e-8);
        rested = step.frame;
    }
    let reload = rested
        .solve_creep_step(&load(FORCE), controls(), || false)
        .unwrap();
    let after = reload
        .frame
        .creep_recovery_plastic_fiber(0, 0, 0)
        .unwrap()
        .1;
    let without_dwell = immediate
        .frame
        .creep_recovery_plastic_fiber(0, 0, 0)
        .unwrap()
        .1;
    assert!(after.creep_strain() > without_dwell.creep_strain());
    assert!(after.equivalent_creep_strain() > without_dwell.equivalent_creep_strain());
    close(after.plastic_strain(), original.plastic_strain(), 2e-11);
    close(
        after.equivalent_plastic_strain(),
        original.equivalent_plastic_strain(),
        2e-11,
    );
    assert!(after.resistance_strain() < after.equivalent_creep_strain());
    close(reload.energy.plastic_dissipation_j, 0.0, 1e-8);
    close(reload.energy.balance_residual_j, 0.0, controls().energy_j);
}

#[test]
fn joint_motion_shares_the_recovery_clock_and_preserves_all_histories_on_retiming() {
    let yield_stress = 0.3e6;
    let member_mass = 1e4;
    // The checks below hold the displacement to 2e-11 m. On this path the free
    // residual's stiffness is 6.2e5 N/m at the consistent-mass steps (member
    // 2.9e5 N/m plus m/dt^2 3.3e5 N/m), so the shared 1e-5 N allowance admits
    // 1.6e-11 m in one step, and the next step carries that error on. The solve
    // asks for 1e-8 N here, which puts its own stopping error well below the checks.
    let tolerances = FrameCreepTolerances {
        force_n: 1e-8,
        ..controls()
    };
    for (formulation, mass) in [
        (FrameMassFormulation::Lumped, member_mass / 2.0),
        (FrameMassFormulation::ConsistentReference, member_mass / 3.0),
    ] {
        let mut accepted = frame(
            &law(0.1, RECOVERY, yield_stress, HARDENING),
            member_mass,
            formulation,
        );
        let (mut u, mut v, mut creep, mut alpha_c, mut plastic, mut alpha_p, mut resistance) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        for (dt, force) in [
            (0.1, FORCE),
            (0.1, 1.2 * FORCE),
            (0.05, FORCE),
            (0.1, 0.2 * FORCE),
            (0.1, 0.0),
        ] {
            let before = accepted
                .creep_recovery_plastic_fiber(0, 0, 0)
                .unwrap()
                .1
                .clone();
            let prepared = accepted
                .with_creep_laws(
                    &vec![law(dt, RECOVERY, yield_stress, HARDENING); 6],
                    "T",
                    || false,
                )
                .unwrap();
            assert_eq!(
                prepared.creep_recovery_plastic_fiber(0, 0, 0).unwrap().1,
                &before
            );
            assert_eq!(prepared.velocities(), accepted.velocities());
            let increments = |stress: f64| {
                (
                    creep_increment(stress, resistance, dt),
                    ((stress.abs() - yield_stress) / HARDENING - alpha_p).max(0.0),
                )
            };
            let strain = |stress: f64| {
                let (dc, dp) = increments(stress);
                stress / E + creep + plastic + stress.signum() * (dc + dp)
            };
            let residual = |stress: f64| {
                2.0 * AREA * stress + mass * (strain(stress) - u - dt * v) / (dt * dt) - force
            };
            let (mut lo, mut hi) = (-4.0 * STRESS, 4.0 * STRESS);
            assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
            for _ in 0..100 {
                let mid = 0.5 * (lo + hi);
                if residual(mid) > 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                }
            }
            let stress = 0.5 * (lo + hi);
            let next_u = strain(stress);
            let acceleration = (force - 2.0 * AREA * stress) / mass;
            let next_v = v + dt * acceleration;
            let (dc, dp) = increments(stress);
            let step = prepared
                .advance_creep_backward_euler(&load(force), dt, tolerances, || false)
                .unwrap();
            assert_eq!(
                prepared.creep_recovery_plastic_fiber(0, 0, 0).unwrap().1,
                &before
            );
            close(step.frame.displacements()[3], next_u, 2e-11);
            close(step.frame.velocities()[3], next_v, 2e-10);
            close(step.frame.accelerations()[3], acceleration, 2e-9);
            close(
                step.energy.material.creep_dissipation_j,
                2.0 * AREA * stress.abs() * dc,
                1e-7,
            );
            close(
                step.energy.material.plastic_dissipation_j,
                2.0 * AREA * yield_stress * dp,
                1e-7,
            );
            close(
                step.energy.inertial_numerical_dissipation_j,
                0.5 * mass * (next_v - v).powi(2),
                1e-7,
            );
            close(
                step.energy.caloric_heat_j(),
                step.energy.material.physical_dissipation_j,
                0.0,
            );
            close(step.energy.balance_residual_j, 0.0, tolerances.energy_j);
            creep += stress.signum() * dc;
            alpha_c += dc;
            plastic += stress.signum() * dp;
            alpha_p += dp;
            resistance = (resistance + dc) / (1.0 + RECOVERY * dt);
            for section in 0..3 {
                for fiber in 0..2 {
                    let (law, state) = step
                        .frame
                        .creep_recovery_plastic_fiber(0, section, fiber)
                        .unwrap();
                    close(state.creep_strain(), creep, 2e-11);
                    close(state.equivalent_creep_strain(), alpha_c, 2e-11);
                    close(state.plastic_strain(), plastic, 2e-11);
                    close(state.equivalent_plastic_strain(), alpha_p, 2e-11);
                    close(state.resistance_strain(), resistance, 2e-12);
                    assert_eq!(law.kinetics().duration_s(), dt);
                }
            }
            (u, v) = (next_u, next_v);
            accepted = step.frame;
        }
        assert!(alpha_p > 0.0 && alpha_c > 0.0 && resistance < alpha_c);
    }
}

#[test]
fn joint_rebinding_cannot_replace_storage_recovery_parameters_or_the_physical_clock() {
    let original_law = law(DT, RECOVERY, YIELD, HARDENING);
    let accepted = frame(&original_law, 1.0, FrameMassFormulation::Lumped)
        .solve_creep_step(&load(FORCE), controls(), || false)
        .unwrap()
        .frame;
    let before = accepted
        .creep_recovery_plastic_fiber(0, 0, 0)
        .unwrap()
        .1
        .clone();
    for changed in [
        law(DT, 2.0 * RECOVERY, YIELD, HARDENING),
        law(DT, RECOVERY, 1.1 * YIELD, HARDENING),
        law(DT, RECOVERY, YIELD, 2.0 * HARDENING),
    ] {
        assert!(
            accepted
                .with_creep_laws(&vec![changed; 6], "T", || false)
                .is_err()
        );
    }
    let no_recovery = UniaxialNortonJ2::new(original_law.backbone().norton_j2().clone());
    assert!(
        accepted
            .with_creep_laws(&vec![no_recovery; 6], "T", || false)
            .is_err()
    );
    let prepared = accepted
        .with_creep_laws(
            &vec![law(DT / 2.0, RECOVERY, YIELD, HARDENING); 6],
            "T",
            || false,
        )
        .unwrap();
    assert_eq!(
        prepared.creep_recovery_plastic_fiber(0, 0, 0).unwrap().1,
        &before
    );
    assert!(matches!(
        prepared.advance_creep_backward_euler(&load(FORCE), DT, controls(), || false),
        Err(FrameCreepError::Invalid(_))
    ));
    let mut mixed = vec![law(DT, RECOVERY, YIELD, HARDENING); 6];
    mixed[1] = law(DT / 2.0, RECOVERY, YIELD, HARDENING);
    let unequal = accepted.with_creep_laws(&mixed, "T", || false).unwrap();
    assert!(matches!(
        unequal.solve_creep_step(&load(FORCE), controls(), || false),
        Err(FrameCreepError::Invalid(_))
    ));
    let polls = Cell::new(0);
    assert!(matches!(
        accepted.solve_creep_step(&load(FORCE), controls(), || {
            polls.set(polls.get() + 1);
            polls.get() >= 12
        }),
        Err(FrameCreepError::Cancelled)
    ));
    let mut untimed = accepted.clone();
    assert!(matches!(
        untimed.newmark(&[0.0, 0.0], DT, (1.0, 0.0), Rayleigh { a0: 0.0, a1: 0.0 }),
        Err(SolidError::MaterialRefused { .. })
    ));
    assert_eq!(
        accepted.creep_recovery_plastic_fiber(0, 0, 0).unwrap().1,
        &before
    );
    assert!(
        accepted
            .solve_creep_step(&load(FORCE), controls(), || false)
            .is_ok()
    );
}
