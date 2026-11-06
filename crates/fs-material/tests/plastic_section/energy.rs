use super::*;
use fs_material::plastic::J2State;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{PowerLawOverstressJ2, VoceOverstressJ2};
use fs_material::plastic::section::energy::{NativeSectionEnergyLaw, SectionEnergyError};

fn gates() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-9,
        yield_pa: 1e-8,
    }
}

fn yield_heat<L: NativeSectionEnergyLaw<State = J2State>>(law: L) {
    let section = section(law, 4);
    let old = section.initial_state();
    let next = section
        .evaluate([0.002, 0.0, 0.001, 0.08, 0.0, 0.01], &old)
        .unwrap();
    let account = section
        .audit_energy(&old, &next.state, gates(), 1e-10)
        .unwrap();
    let expected_heat: f64 = section
        .points()
        .iter()
        .zip(&next.state.points)
        .map(|(point, state)| point.weight_m * 3.0 * state.history.alpha)
        .sum();
    assert!(expected_heat > 0.0);
    near(account.physical_dissipation_j_m2, expected_heat, 1e-12);
    assert!(account.final_state.hardening_j_m2 > 0.0);
    assert!(account.final_state.kinematic_hardening_j_m2 > 0.0);
    assert!(account.numerical_dissipation_j_m2 > 0.0);
    let work: f64 = next
        .resultants
        .iter()
        .zip(next.state.generalized_strain)
        .map(|(force, strain)| force * strain)
        .sum();
    near(account.generalized_work_j_m2, work, 1e-11);
    near(
        work + account.transverse_residual_work_j_m2,
        account.stored_change_j_m2() + expected_heat + account.numerical_dissipation_j_m2,
        1e-10,
    );
    assert_eq!(old, section.initial_state());
    let mut substituted = next.state.clone();
    substituted.points[0].history.alpha += 1e-4;
    assert!(matches!(
        section.audit_energy(&old, &substituted, gates(), 1e-10),
        Err(SectionEnergyError::State { .. })
    ));
}

#[test]
fn native_j2_and_voce_heat_excludes_recoverable_hardening_and_numerical_loss() {
    yield_heat(plastic());
    yield_heat(VoceJ2::new(plastic(), 4.0, 30.0).unwrap());
    let law = ChabocheJ2::new(
        VoceJ2::new(plastic(), 4.0, 30.0).unwrap(),
        [BackstressBranch::new(80.0, 20.0).unwrap()],
    )
    .unwrap();
    let section = section(law.clone(), 4);
    let old = section.initial_state();
    let next = section
        .evaluate([0.01, -0.002, 0.004, 0.02, 0.0, 0.01], &old)
        .unwrap();
    let account = section
        .audit_energy(&old, &next.state, gates(), 1e-10)
        .unwrap();
    let mut physical = 0.0;
    let mut recovery = 0.0;
    for (i, point) in section.points().iter().enumerate() {
        let native = law
            .advance_with_energy(
                &old.points[i].strain,
                &old.points[i].history,
                &next.state.points[i].strain,
                gates(),
            )
            .unwrap();
        assert_eq!(native.state, next.state.points[i].history);
        physical += point.weight_m * native.energy.physical_dissipation_j_m3;
        recovery += point.weight_m * native.recovery_dissipation_j_m3;
    }
    assert!(recovery > 0.0);
    near(account.physical_dissipation_j_m2, physical, 1e-12);
    assert!((account.physical_dissipation_j_m2 - physical - recovery).abs() > 1e-12);
}

fn timed<L: NativeSectionEnergyLaw<State = J2State>>(law: L) {
    let section = section(law, 2);
    let g = [0.015, 0.0, 0.002, 0.0, 0.0, 0.0];
    let first = section.evaluate(g, &section.initial_state()).unwrap();
    let saved = first.state.clone();
    let held = section.evaluate(g, &first.state).unwrap();
    let account = section
        .audit_energy(&first.state, &held.state, gates(), 1e-10)
        .unwrap();
    assert!(held.state.points[0].history.alpha > first.state.points[0].history.alpha);
    assert!(account.physical_dissipation_j_m2 > 0.0);
    near(account.generalized_work_j_m2, 0.0, 0.0);
    near(account.balance_residual_j_m2, 0.0, 1e-10);
    // A hold is a real interval for an overstressed solid. Replacing its
    // returned history with the old one must fail instead of losing its heat.
    assert!(
        section
            .audit_energy(&first.state, &first.state, gates(), 1e-10)
            .is_err()
    );
    assert_eq!(first.state, saved);
}

#[test]
fn kinetic_accounts_replay_exactly_one_interval_and_keep_valid_old_overstress() {
    timed(LinearOverstressJ2::new(plastic(), 20.0, 0.1).unwrap());
    timed(PowerLawOverstressJ2::new(plastic(), 20.0, 5.0, 2.0, 0.1).unwrap());
    timed(
        VoceOverstressJ2::new_power_law(
            VoceJ2::new(plastic(), 4.0, 30.0).unwrap(),
            20.0,
            5.0,
            2.0,
            0.1,
        )
        .unwrap(),
    );
}

#[test]
fn finite_transverse_residual_work_is_signed_and_not_an_extra_heat_source() {
    let law = plastic();
    let section = LayeredPlaneStressSection::new(
        vec![PlaneStressLayer {
            bottom_z_m: -0.1,
            top_z_m: 0.1,
            subdivisions: 2,
            law: law.clone(),
        }],
        ThicknessQuadrature::Gauss3,
        ReductionConfig {
            stress_tolerance_pa: 1e-2,
            ..controls()
        },
        budget(),
    )
    .unwrap();
    let old = section.initial_state();
    let tight = super::section(law.clone(), 2);
    let mut next = tight
        .evaluate([0.01, 0.0, 0.002, 0.0, 0.0, 0.0], &old)
        .unwrap()
        .state;
    // Construct an endpoint within the explicitly loose transverse allowance,
    // retaining exactly the native update at its complete tensor strain.
    for (after, before) in next.points.iter_mut().zip(&old.points) {
        after.strain[2] += 1e-6;
        after.history = law.update_state(&after.strain, &before.history);
    }
    let account = section.audit_energy(&old, &next, gates(), 1e-10).unwrap();
    assert!(account.transverse_residual_work_j_m2.abs() > 1e-10);
    near(
        account.material_work_j_m2 - account.generalized_work_j_m2,
        account.transverse_residual_work_j_m2,
        1e-11,
    );
    near(account.balance_residual_j_m2, 0.0, 1e-10);
    let expected: f64 = section
        .points()
        .iter()
        .zip(&next.points)
        .map(|(p, s)| p.weight_m * law.yield_stress * s.history.alpha)
        .sum();
    near(account.physical_dissipation_j_m2, expected, 1e-12);
}
