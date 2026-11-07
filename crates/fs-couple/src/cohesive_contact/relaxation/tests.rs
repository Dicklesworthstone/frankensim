use super::*;
use crate::cohesive_contact::CohesiveContactConfig;
use fs_material::cohesive::relaxation::{IntactViscoelasticBond, ViscoelasticBondLimits};
use fs_material::visco::GeneralizedMaxwell;

fn config() -> CohesiveContactConfig {
    CohesiveContactConfig {
        max_iterations: 8,
        max_material_evaluations: 32,
        max_backtracks: 8,
        jump_absolute_tolerance_m: 1e-14,
        jump_relative_tolerance: 1e-11,
        maximum_jump_increment_m: 1e-4,
        maximum_force_n: 100.0,
        material_energy_tolerance_j_m2: 1e-9,
    }
}
fn law() -> IntactViscoelasticBond {
    IntactViscoelasticBond::new(
        GeneralizedMaxwell::new(1e6, vec![(2e6, 0.01)]).unwrap(),
        GeneralizedMaxwell::new(4e5, vec![(8e5, 0.025)]).unwrap(),
        0.001,
        ViscoelasticBondLimits {
            maximum_jump_m: 1e-4,
            maximum_rate_m_s: 1.0,
            maximum_elapsed_s: 5.0,
        },
    )
    .unwrap()
}
fn close(a: f64, b: f64, t: f64) {
    assert!(
        (a - b).abs() <= t * (1.0 + a.abs().max(b.abs())),
        "{a:e} != {b:e}"
    );
}

#[test]
fn full_mechanical_cross_compliance_recovers_the_original_affine_memory_force() {
    let law = law();
    let prepared = law.prepare(0.03).unwrap();
    let point = law
        .bonded_point()
        .trial(&prepared, [1e-5, 2e-5, -1e-5], 1e-9)
        .unwrap()
        .point;
    let desired = [2e-5, -1e-5, 0.5e-5];
    let area = 2e-4;
    let material = point.trial(&prepared, desired, 1e-9).unwrap();
    let compliance = [[2e-6, 3e-7, -2e-7], [3e-7, 1e-6, 1e-7], [-2e-7, 1e-7, 3e-6]];
    let free = core::array::from_fn(|i| {
        desired[i]
            + (0..3)
                .map(|j| compliance[i][j] * area * material.mean_traction_pa[j])
                .sum::<f64>()
    });
    let mut solver = CohesiveLoadLine::new(area, compliance, config()).unwrap();
    let result = solver
        .trial_relaxation(&point, &prepared, free, || true)
        .unwrap();
    for i in 0..3 {
        close(result.material.point.jump_m()[i], desired[i], 1e-14);
        close(
            result.mean_force_n[i],
            area * material.mean_traction_pa[i],
            1e-12,
        );
    }
    assert_ne!(
        result.material.traction_pa,
        result.material.mean_traction_pa
    );
    assert!(
        solver
            .trial_relaxation(&point, &prepared, free, || false)
            .is_err()
    );
    let mut polls = 0;
    assert!(
        solver
            .trial_relaxation(&point, &prepared, free, || {
                polls += 1;
                polls < 3
            })
            .is_err()
    );
    assert_eq!(point.elapsed_s(), 0.03);
}

#[test]
fn original_load_line_advances_bonded_mass_and_balances_work_storage_and_heat() {
    let law = law();
    let dt = 0.001;
    let mass = 1.0;
    let area = 1e-4;
    let prepared = law.prepare(dt).unwrap();
    let compliance = 0.5 * dt * dt / mass;
    let mut solver = CohesiveLoadLine::new(
        area,
        [
            [compliance, 0.0, 0.0],
            [0.0, compliance, 0.0],
            [0.0, 0.0, compliance],
        ],
        config(),
    )
    .unwrap();
    let mut point = law.bonded_point();
    let mut velocity = [0.0; 3];
    let force = [0.5, -0.2, 0.1];
    let mut external_work = 0.0;
    let mut heat = 0.0;
    for _ in 0..2000 {
        let old = point.jump_m();
        let free = core::array::from_fn(|i| old[i] + dt * velocity[i] + compliance * force[i]);
        let step = solver
            .trial_relaxation(&point, &prepared, free, || true)
            .unwrap();
        let new = step.material.point.jump_m();
        let next_velocity =
            core::array::from_fn::<_, 3, _>(|i| 2.0 * (new[i] - old[i]) / dt - velocity[i]);
        let work = (0..3).map(|i| force[i] * (new[i] - old[i])).sum::<f64>();
        let kinetic = (0..3)
            .map(|i| {
                0.5 * mass * (next_velocity[i] - velocity[i]) * (next_velocity[i] + velocity[i])
            })
            .sum::<f64>();
        close(
            work,
            kinetic
                + area * (step.material.stored_change_j_m2 + step.material.viscous_dissipated_j_m2),
            2e-12,
        );
        external_work += work;
        heat += area * step.material.viscous_dissipated_j_m2;
        velocity = next_velocity;
        point = step.material.point;
    }
    let kinetic = 0.5 * mass * velocity.iter().map(|v| v * v).sum::<f64>();
    close(
        external_work,
        kinetic + area * point.stored_energy_j_m2() + heat,
        2e-11,
    );
    assert!(heat > 0.0);
    for i in 0..3 {
        let relaxed = if i == 0 { 1e6 } else { 4e5 };
        close(
            point.jump_m()[i],
            force[i] * law.thickness_m() / (area * relaxed),
            2e-10,
        );
        assert!(velocity[i].abs() < 1e-9);
    }
}
