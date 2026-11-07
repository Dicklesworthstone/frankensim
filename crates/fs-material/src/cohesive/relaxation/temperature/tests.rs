use super::*;
use crate::cohesive::relaxation::ViscoelasticBondLimits;
use crate::visco::GeneralizedMaxwell;

fn law(
    normal: [(f64, f64); 2],
    shear: [(f64, f64); 1],
    equilibrium: [f64; 2],
) -> IntactViscoelasticBond {
    IntactViscoelasticBond::new(
        GeneralizedMaxwell::new(equilibrium[0], normal.to_vec()).unwrap(),
        GeneralizedMaxwell::new(equilibrium[1], shear.to_vec()).unwrap(),
        0.001,
        ViscoelasticBondLimits {
            maximum_jump_m: 1e-4,
            maximum_rate_m_s: 1.0,
            maximum_elapsed_s: 10.0,
        },
    )
    .unwrap()
}
fn original() -> IntactViscoelasticBond {
    law([(2e6, 0.1), (3e5, 0.8)], [(8e5, 0.25)], [1e6, 4e5])
}
fn loaded(law: &IntactViscoelasticBond) -> ViscoelasticBondPoint<'_> {
    let mut point = law.bonded_point();
    for (dt, gap) in [
        (0.02, [2e-5, -1e-5, 0.5e-5]),
        (0.07, [-1e-5, 2e-5, -0.4e-5]),
        (0.13, [0.5e-5, -0.7e-5, 1e-5]),
    ] {
        point = point
            .trial(&law.prepare(dt).unwrap(), gap, 1e-10)
            .unwrap()
            .point;
    }
    point
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}
fn same_memory(a: &ViscoelasticBondPoint<'_>, b: &ViscoelasticBondPoint<'_>) {
    assert_eq!(a.elapsed_s(), b.elapsed_s());
    assert_eq!(a.jump_m(), b.jump_m());
    for i in 0..3 {
        assert_eq!(a.states[i].q, b.states[i].q);
        assert_eq!(a.states[i].strain, b.states[i].strain);
        assert_eq!(a.states[i].work, b.states[i].work);
        assert_eq!(a.states[i].dissipated, b.states[i].dissipated);
    }
}

#[test]
fn changed_coefficients_preserve_physical_memory_and_close_signed_spring_work() {
    let cold = original();
    let hot = law([(1.1e6, 0.005), (7e5, 0.03)], [(3e5, 0.04)], [6e5, 2e5]);
    let before = loaded(&cold);
    let step = before.transfer_temperature(&hot, 1e-12).unwrap();
    same_memory(&before, &step.point);
    assert!(core::ptr::eq(step.point.law(), &hot));
    let mut expected_work = 0.0;
    for i in 0..3 {
        let mode = usize::from(i > 0);
        let a = &cold.spectra[mode];
        let b = &hot.spectra[mode];
        let s = &before.states[i];
        expected_work += 0.5
            * cold.thickness
            * ((b.e_inf - a.e_inf) * s.strain * s.strain
                + (0..s.q.len())
                    .map(|j| (b.terms[j].0 - a.terms[j].0) * s.q[j] * s.q[j])
                    .sum::<f64>());
        let expected_stress =
            b.e_inf * s.strain + (0..s.q.len()).map(|j| b.terms[j].0 * s.q[j]).sum::<f64>();
        close(step.point.traction_pa()[i], expected_stress, 1e-14);
    }
    assert!(expected_work < 0.0);
    close(step.parameter_work_j_m2, expected_work, 1e-14);
    close(
        step.stored_change_j_m2,
        step.point.stored_energy_j_m2() - before.stored_energy_j_m2(),
        1e-14,
    );
    close(step.defect_j_m2, 0.0, 1e-14);
    let back = step.point.transfer_temperature(&cold, 1e-12).unwrap();
    same_memory(&before, &back.point);
    close(
        step.parameter_work_j_m2 + back.parameter_work_j_m2,
        0.0,
        1e-14,
    );
    assert_eq!(before.traction_pa(), back.point.traction_pa());
    assert_eq!(before.stored_energy_j_m2(), back.point.stored_energy_j_m2());
}

#[test]
fn times_may_cross_without_reordering_branches_or_resetting_creep() {
    let cold = original();
    // The same branch 0 is now slower than branch 1. No sorting/remapping.
    let hot = law([(2e6, 0.7), (3e5, 0.003)], [(8e5, 0.01)], [1e6, 4e5]);
    let before = loaded(&cold);
    assert_ne!(before.states[0].q[0], before.states[0].q[1]);
    let transfer = before.transfer_temperature(&hot, 1e-15).unwrap();
    same_memory(&before, &transfer.point);
    assert_eq!(transfer.point.traction_pa(), before.traction_pa());
    assert_eq!(transfer.parameter_work_j_m2, 0.0);
    assert_eq!(transfer.stored_change_j_m2, 0.0);
    assert_eq!(transfer.defect_j_m2, 0.0);
    let dt = 0.2;
    let step = transfer
        .point
        .trial(&hot.prepare(dt).unwrap(), before.jump_m(), 1e-10)
        .unwrap();
    assert_eq!(step.work_j_m2, 0.0);
    assert!(step.viscous_dissipated_j_m2 > 0.0);
    for i in 0..3 {
        let spectrum = &hot.spectra[usize::from(i > 0)];
        let state = &before.states[i];
        let expected = spectrum.e_inf * state.strain
            + spectrum
                .terms
                .iter()
                .zip(&state.q)
                .map(|(&(e, tau), &q)| e * q * (-dt / tau).exp())
                .sum::<f64>();
        close(step.traction_pa[i], expected, 5e-13);
    }
    assert_eq!(step.point.elapsed_s(), before.elapsed_s() + dt);
    close(
        step.stored_change_j_m2,
        -step.viscous_dissipated_j_m2,
        1e-12,
    );
}

#[test]
fn tiny_changing_shear_work_survives_large_unchanged_normal_storage() {
    let cold = law([(2e6, 0.1), (3e5, 0.8)], [(8e5, 0.25)], [1e24, 4e5]);
    let hot = law([(2e6, 0.1), (3e5, 0.8)], [(8e5, 0.25)], [1e24, 4e5 + 1e-4]);
    // Both points follow actual history; normal energy dwarfs the shear change.
    let before = cold
        .bonded_point()
        .trial(&cold.prepare(0.1).unwrap(), [1e-5, 2e-5, 0.0], 100.0)
        .unwrap()
        .point;
    let transfer = before.transfer_temperature(&hot, 1e-20).unwrap();
    let epsilon = before.states[1].strain;
    let expected =
        0.5 * cold.thickness * (hot.spectra[1].e_inf - cold.spectra[1].e_inf) * epsilon * epsilon;
    assert!(expected > 0.0);
    assert_eq!(
        transfer.point.stored_energy_j_m2() - before.stored_energy_j_m2(),
        0.0
    );
    assert!((transfer.parameter_work_j_m2 - expected).abs() < expected * 1e-14);
    assert!((transfer.stored_change_j_m2 - expected).abs() < expected * 1e-14);
    same_memory(&before, &transfer.point);
}

#[test]
fn held_refresh_refuses_geometry_branch_or_limit_changes_without_mutation() {
    let cold = original();
    let before = loaded(&cold);
    let mut geometry = cold.clone();
    geometry.thickness *= 1.1;
    let mut count = cold.clone();
    count.spectra[0].terms.pop();
    let mut limits = cold.clone();
    limits.limits.maximum_elapsed_s *= 0.5;
    for invalid in [&geometry, &count, &limits] {
        assert!(before.transfer_temperature(invalid, 1e-12).is_err());
    }
    for tolerance in [0.0, -1.0, f64::NAN] {
        assert!(before.transfer_temperature(&cold, tolerance).is_err());
    }
    let accepted = before.transfer_temperature(&cold, 1e-12).unwrap();
    same_memory(&before, &accepted.point);
    assert_eq!(before.traction_pa(), accepted.point.traction_pa());
    assert_eq!(accepted.parameter_work_j_m2, 0.0);
    let retry = before
        .trial(&cold.prepare(0.04).unwrap(), [1e-5, 0.0, 0.0], 1e-10)
        .unwrap();
    let control = loaded(&cold)
        .trial(&cold.prepare(0.04).unwrap(), [1e-5, 0.0, 0.0], 1e-10)
        .unwrap();
    assert_eq!(retry.traction_pa, control.traction_pa);
    assert_eq!(retry.work_j_m2, control.work_j_m2);
    assert_eq!(
        retry.viscous_dissipated_j_m2,
        control.viscous_dissipated_j_m2
    );
}
