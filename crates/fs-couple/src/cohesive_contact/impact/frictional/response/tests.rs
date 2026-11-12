//! Independent mass, Coulomb-law and work identities through the original native solvers.
use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_solver::LinearOp;
use fs_sparse::precond::IdentityPrecond;
use std::cell::Cell;

/// Analytic original mass operators; the production checked CG performs every inverse action.
struct Mass3([f64; 9]);
impl LinearOp for Mass3 {
    fn n(&self) -> usize {
        3
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for (row, value) in y.iter_mut().enumerate() {
            *value = (0..3).map(|col| self.0[3 * row + col] * x[col]).sum();
        }
    }
    // Every instance in this file (isotropic, and the 0.4-coupled mass) is a symmetric matrix, so its transpose is itself.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}
fn isotropic(mass: f64) -> Mass3 {
    Mass3([mass, 0.0, 0.0, 0.0, mass, 0.0, 0.0, 0.0, mass])
}
fn kinetic(mass: &Mass3, velocity: [f64; 3]) -> f64 {
    let mut twice_energy = 0.0;
    for row in 0..3 {
        for col in 0..3 {
            twice_energy += velocity[row] * mass.0[3 * row + col] * velocity[col];
        }
    }
    0.5 * twice_energy
}
fn controls() -> FrictionalResponseControls {
    FrictionalResponseControls {
        mass_solve: CheckedCgConfig {
            tolerance: 1e-13,
            max_iterations: 64,
            max_corrections: 2,
        },
        contact_solve: FrictionalContactConfig::default(),
        maximum_bodies: 4,
        maximum_trace_terms: 128,
        maximum_workspace_entries: 2_000_000,
        velocity_tolerance_m_s: 1e-9,
        energy_tolerance_j: 1e-9,
    }
}
fn trace(body: usize, coefficient: f64) -> [[TraceTerm; 1]; 3] {
    core::array::from_fn(|dof| {
        [TraceTerm {
            body,
            dof,
            coefficient,
        }]
    })
}
fn contact<'a>(id: u64, trace: [&'a [TraceTerm]; 3], mu: f64) -> FrictionalContact<'a> {
    FrictionalContact {
        id,
        gap_m: 0.0,
        activation_tolerance_m: 0.0,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
        prescribed_velocity_m_s: [0.0; 3],
        trace,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-8,
        "{actual:e} != {expected:e}"
    );
}
fn close3(actual: &[f64], expected: [f64; 3]) {
    assert_eq!(actual.len(), 3);
    for axis in 0..3 {
        close(actual[axis], expected[axis]);
    }
}

#[test]
fn sliding_and_sticking_use_final_slip_and_one_simultaneous_impulse() {
    let mass = isotropic(1.0);
    let rows = trace(0, 1.0);
    let c = contact(5, [&rows[0], &rows[1], &rows[2]], 0.5);
    for (before, expected_impulse, after, loss) in [
        ([-2.0, 3.0, 0.0], [2.0, -1.0, 0.0], [0.0, 2.0, 0.0], 4.5),
        ([-2.0, 0.4, -0.3], [2.0, -0.4, 0.3], [0.0; 3], 2.125),
    ] {
        let bodies = [Body {
            mass: &mass,
            velocity: &before,
            preconditioner: &IdentityPrecond,
        }];
        let plan = prepare(&bodies, &[c], controls(), &mut || false).unwrap();
        close3(&plan.body_impulses[0], expected_impulse);
        let report = finish(
            plan,
            &[c],
            &[&after],
            &[kinetic(&mass, after) - kinetic(&mass, before)],
            &[0.0],
            controls(),
            &mut || false,
        )
        .unwrap();
        close3(&report.contact.impulses[0], expected_impulse);
        close3(&report.relative_velocities_after_m_s[0], after);
        close(report.contact.dissipated_energy, loss);
        close(report.kinetic_change_j, -loss);
        close(report.work_energy_residual_j, 0.0);
        close(report.dissipation_residual_j, 0.0);
    }
}

#[test]
fn nondiagonal_mass_couples_normal_and_tangent_without_normal_first_clipping() {
    let mass = Mass3([2.0, 0.4, 0.0, 0.4, 1.0, 0.0, 0.0, 0.0, 1.0]);
    let before = [-1.0, 1.0, 0.0];
    let bodies = [Body {
        mass: &mass,
        velocity: &before,
        preconditioner: &IdentityPrecond,
    }];
    let rows = trace(0, 1.0);
    let c = contact(9, [&rows[0], &rows[1], &rows[2]], 1.0);
    let plan = prepare(&bodies, &[c], controls(), &mut || false).unwrap();
    // The sticking endpoint is zero, so the exact simultaneous impulse is -M*v.
    // A normal-first update would instead produce the incorrect normal impulse 1.84.
    close3(&plan.body_impulses[0], [1.6, -0.6, 0.0]);
    let report = finish(
        plan,
        &[c],
        &[&[0.0; 3]],
        &[-1.1],
        &[0.0],
        controls(),
        &mut || false,
    )
    .unwrap();
    let inverse = [
        1.0 / 1.84,
        -0.4 / 1.84,
        0.0,
        -0.4 / 1.84,
        2.0 / 1.84,
        0.0,
        0.0,
        0.0,
        1.0,
    ];
    for (&actual, expected) in report.contact_mobility_per_kg.iter().zip(inverse) {
        close(actual, expected);
    }
    close(report.contact.dissipated_energy, 1.1);
    close(report.work_energy_residual_j, 0.0);
}

#[test]
fn moving_surface_accounts_for_normal_and_tangential_external_work() {
    let mass = isotropic(1.0);
    let bodies = [Body {
        mass: &mass,
        velocity: &[0.0; 3],
        preconditioner: &IdentityPrecond,
    }];
    let rows = trace(0, 1.0);
    let mut c = contact(6, [&rows[0], &rows[1], &rows[2]], 0.5);
    c.prescribed_velocity_m_s = [-2.0, -3.0, 0.0];
    let plan = prepare(&bodies, &[c], controls(), &mut || false).unwrap();
    close3(&plan.body_impulses[0], [2.0, 1.0, 0.0]);
    let report = finish(
        plan,
        &[c],
        &[&[2.0, 1.0, 0.0]],
        &[2.5],
        &[0.0],
        controls(),
        &mut || false,
    )
    .unwrap();
    close3(&report.relative_velocities_after_m_s[0], [0.0, -2.0, 0.0]);
    close(report.prescribed_work_j, 7.0);
    close(report.relative_impulse_work_j, -4.5);
    close(report.kinetic_change_j, 2.5);
    close(report.contact.dissipated_energy, 4.5);
    close(report.work_energy_residual_j, 0.0);
}

#[test]
fn unequal_bodies_recoil_with_exact_restitution_and_friction_momentum_balance() {
    let masses = [isotropic(2.0), isotropic(3.0)];
    let before = [[1.0, 2.0, 0.0], [-1.0, -1.0, 0.0]];
    let bodies: [Body<'_>; 2] = core::array::from_fn(|body| Body {
        mass: &masses[body],
        velocity: &before[body],
        preconditioner: &IdentityPrecond,
    });
    let rows: [[TraceTerm; 2]; 3] = core::array::from_fn(|dof| {
        [
            TraceTerm {
                body: 0,
                dof,
                coefficient: -1.0,
            },
            TraceTerm {
                body: 1,
                dof,
                coefficient: 1.0,
            },
        ]
    });
    let mut c = contact(8, [&rows[0], &rows[1], &rows[2]], 0.25);
    c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 };
    let plan = prepare(&bodies, &[c], controls(), &mut || false).unwrap();
    close3(&plan.body_impulses[0], [-3.6, -0.9, 0.0]);
    close3(&plan.body_impulses[1], [3.6, 0.9, 0.0]);
    let after = [[-0.8, 1.55, 0.0], [0.2, -0.7, 0.0]];
    let changes: [f64; 2] = core::array::from_fn(|body| {
        kinetic(&masses[body], after[body]) - kinetic(&masses[body], before[body])
    });
    let report = finish(
        plan,
        &[c],
        &[&after[0], &after[1]],
        &changes,
        &[0.0; 2],
        controls(),
        &mut || false,
    )
    .unwrap();
    close3(&report.contact.impulses[0], [3.6, 0.9, 0.0]);
    close3(&report.relative_velocities_after_m_s[0], [1.0, -2.25, 0.0]);
    for axis in 0..3 {
        close(
            2.0 * after[0][axis] + 3.0 * after[1][axis],
            2.0 * before[0][axis] + 3.0 * before[1][axis],
        );
    }
    close(report.contact.dissipated_energy, 4.1625);
    close(report.kinetic_change_j, -4.1625);
    close(report.work_energy_residual_j, 0.0);
}

#[test]
fn shared_body_cross_contact_terms_preserve_momentum_during_two_sticking_impacts() {
    let masses = [isotropic(1.0), isotropic(2.0), isotropic(3.0)];
    let before = [[3.0, 2.0, 0.0], [0.0; 3], [-1.0, -1.0, 0.0]];
    let bodies: [Body<'_>; 3] = core::array::from_fn(|body| Body {
        mass: &masses[body],
        velocity: &before[body],
        preconditioner: &IdentityPrecond,
    });
    let rows: [[[TraceTerm; 2]; 3]; 2] = core::array::from_fn(|pair| {
        core::array::from_fn(|dof| {
            [
                TraceTerm {
                    body: pair,
                    dof,
                    coefficient: -1.0,
                },
                TraceTerm {
                    body: pair + 1,
                    dof,
                    coefficient: 1.0,
                },
            ]
        })
    });
    let contacts: [FrictionalContact<'_>; 2] = core::array::from_fn(|pair| {
        contact(
            10 + pair as u64,
            [&rows[pair][0], &rows[pair][1], &rows[pair][2]],
            1.0,
        )
    });
    let plan = prepare(&bodies, &contacts, controls(), &mut || false).unwrap();
    close3(&plan.body_impulses[0], [-3.0, -13.0 / 6.0, 0.0]);
    close3(&plan.body_impulses[1], [0.0, -1.0 / 3.0, 0.0]);
    close3(&plan.body_impulses[2], [3.0, 2.5, 0.0]);
    let after = [[0.0, -1.0 / 6.0, 0.0]; 3];
    let changes: [f64; 3] = core::array::from_fn(|body| {
        kinetic(&masses[body], after[body]) - kinetic(&masses[body], before[body])
    });
    let report = finish(
        plan,
        &contacts,
        &[&after[0], &after[1], &after[2]],
        &changes,
        &[0.0; 3],
        controls(),
        &mut || false,
    )
    .unwrap();
    // Every same-axis pair has W=[[3/2,-1/2],[-1/2,5/6]].
    for axis in 0..3 {
        close(report.contact_mobility_per_kg[axis * 6 + axis], 1.5);
        close(report.contact_mobility_per_kg[axis * 6 + axis + 3], -0.5);
        close(report.contact_mobility_per_kg[(axis + 3) * 6 + axis], -0.5);
        close(
            report.contact_mobility_per_kg[(axis + 3) * 6 + axis + 3],
            5.0 / 6.0,
        );
        close(
            after[0][axis] + 2.0 * after[1][axis] + 3.0 * after[2][axis],
            before[0][axis] + 2.0 * before[1][axis] + 3.0 * before[2][axis],
        );
    }
    close3(&report.contact.impulses[0], [3.0, 13.0 / 6.0, 0.0]);
    close3(&report.contact.impulses[1], [3.0, 2.5, 0.0]);
    close(report.contact.dissipated_energy, 113.0 / 12.0);
    close(report.work_energy_residual_j, 0.0);
}

struct CountedMass(Cell<usize>);
impl LinearOp for CountedMass {
    fn n(&self) -> usize {
        3
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.set(self.0.get() + 1);
        for axis in 0..3 {
            y[axis] = 2.0 * x[axis];
        }
    }
    // A counting diagonal, so its transpose is itself; through apply, the count is unchanged.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

#[test]
fn zero_friction_has_no_tangent_impulse_and_admission_or_endpoint_failure_is_atomic() {
    let mass = CountedMass(Cell::new(0));
    let before = [-2.0, 3.0, -4.0];
    let bodies = [Body {
        mass: &mass,
        velocity: &before,
        preconditioner: &IdentityPrecond,
    }];
    let rows = trace(0, 1.0);
    let c = contact(2, [&rows[0], &rows[1], &rows[2]], 0.0);
    let mut limited = controls();
    limited.maximum_workspace_entries = 0;
    assert!(matches!(
        prepare(&bodies, &[c], limited, &mut || false),
        Err(ResponseError::Budget(_))
    ));
    assert!(matches!(
        prepare(&bodies, &[c], controls(), &mut || true),
        Err(ResponseError::Cancelled)
    ));
    assert!(matches!(
        prepare(&bodies, &[c, c], controls(), &mut || false),
        Err(ResponseError::Invalid(_))
    ));
    assert_eq!(
        mass.0.get(),
        0,
        "admission must precede original mass actions"
    );
    let plan = prepare(&bodies, &[c], controls(), &mut || false).unwrap();
    close3(&plan.body_impulses[0], [4.0, 0.0, 0.0]);
    let report = finish(
        plan,
        &[c],
        &[&[0.0, 3.0, -4.0]],
        &[-4.0],
        &[0.0],
        controls(),
        &mut || false,
    )
    .unwrap();
    close3(&report.contact.impulses[0], [4.0, 0.0, 0.0]);
    close(report.contact.dissipated_energy, 4.0);
    let plan = prepare(&bodies, &[c], controls(), &mut || false).unwrap();
    assert!(matches!(
        finish(
            plan,
            &[c],
            &[&[0.0, 3.1, -4.0]],
            &[-3.39],
            &[0.0],
            controls(),
            &mut || false,
        ),
        Err(ResponseError::Audit {
            quantity: "native endpoint velocity [m/s]",
            ..
        })
    ));
    let plan = prepare(&bodies, &[c], controls(), &mut || false).unwrap();
    assert!(matches!(
        finish(
            plan,
            &[c],
            &[&[0.0, 3.0, -4.0]],
            &[-4.0],
            &[1e-4],
            controls(),
            &mut || false,
        ),
        Err(ResponseError::Audit {
            quantity: "absolute native body work [J]",
            ..
        })
    ));

    // A numerical owner's permitted positive work is still a real physical
    // energy increase. The assembly must not subtract that allowance from its
    // actual-state dissipation gate or publish a fictitious conservative event.
    let target_mass = isotropic(1.0);
    let target_bodies = [Body {
        mass: &target_mass,
        velocity: &[-1.0, 0.0, 0.0],
        preconditioner: &IdentityPrecond,
    }];
    let mut target = c;
    target.normal_law = NormalVelocityLaw::Target { speed: 2.0 };
    let mut permissive_native = controls();
    permissive_native.contact_solve.energy_tolerance = 2.0;
    let plan = prepare(&target_bodies, &[target], permissive_native, &mut || false).unwrap();
    close3(&plan.body_impulses[0], [3.0, 0.0, 0.0]);
    assert!(matches!(
        finish(
            plan,
            &[target],
            &[&[2.0, 0.0, 0.0]],
            &[1.5],
            &[0.0],
            permissive_native,
            &mut || false,
        ),
        Err(ResponseError::Audit {
            quantity: "native contact dissipation [J]",
            residual,
            ..
        }) if (residual - 1.5).abs() <= 1e-10
    ));
}
