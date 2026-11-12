//! G0/G1: independent mass/impact identities through the original native solvers.
use super::*;
use fs_contact::multicontact::MultiImpactDeterminacy;
use fs_solver::LinearOp;
use fs_sparse::precond::IdentityPrecond;
use std::cell::Cell;

struct Scalar(f64);
impl LinearOp for Scalar {
    fn n(&self) -> usize {
        1
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y[0] = self.0 * x[0];
    }
    // A 1x1 scalar, so its transpose is itself.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

fn controls() -> ResponseControls {
    ResponseControls {
        mass_solve: CheckedCgConfig {
            tolerance: 1e-13,
            max_iterations: 64,
            max_corrections: 2,
        },
        contact_solve: MultiImpactConfig::default(),
        maximum_bodies: 4,
        maximum_trace_terms: 32,
        maximum_workspace_entries: 20_000,
        velocity_tolerance_m_s: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}
fn row(id: u64, trace: &[TraceTerm]) -> NormalContact<'_> {
    NormalContact {
        id,
        gap_m: 0.0,
        activation_tolerance_m: 0.0,
        restitution: 0.0,
        prescribed_normal_velocity_m_s: 0.0,
        trace,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-10,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn unequal_masses_exchange_one_opposed_impulse_with_actual_restitution_loss() {
    let masses = [Scalar(2.0), Scalar(3.0)];
    let velocities = [[-2.0], [1.0]];
    let bodies = [
        Body {
            mass: &masses[0],
            velocity: &velocities[0],
            preconditioner: &IdentityPrecond,
        },
        Body {
            mass: &masses[1],
            velocity: &velocities[1],
            preconditioner: &IdentityPrecond,
        },
    ];
    let trace = [
        TraceTerm {
            body: 0,
            dof: 0,
            coefficient: 1.0,
        },
        TraceTerm {
            body: 1,
            dof: 0,
            coefficient: -1.0,
        },
    ];
    let mut contact = row(17, &trace);
    contact.restitution = 0.5;
    let plan = prepare(&bodies, &[contact], controls(), &mut || false).unwrap();
    close(plan.body_impulses[0][0], 5.4);
    close(plan.body_impulses[1][0], -5.4);
    let after = [[0.7], [-0.8]];
    let changes = [
        0.5 * 2.0 * (0.7_f64.powi(2) - 4.0),
        0.5 * 3.0 * ((-0.8_f64).powi(2) - 1.0),
    ];
    let report = finish(
        plan,
        &[contact],
        &[&after[0], &after[1]],
        &changes,
        &[0.0; 2],
        controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.mass_solves.len(), 2);
    close(report.contact_mobility_per_kg[0], 5.0 / 6.0);
    close(report.relative_velocities_after_m_s[0], 1.5);
    close(report.kinetic_change_j, -4.05);
    close(report.contact.dissipated_energy, 4.05);
    close(report.work_energy_residual_j, 0.0);
    close(report.dissipation_residual_j, 0.0);
    close(2.0 * after[0][0] + 3.0 * after[1][0], -1.0);
}

struct Coupled;
impl LinearOp for Coupled {
    fn n(&self) -> usize {
        2
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y[0] = 2.0 * x[0] + 0.5 * x[1];
        y[1] = 0.5 * x[0] + x[1];
    }
    // A fixed symmetric matrix (off-diagonal 0.5 both ways), so its transpose is itself.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

#[test]
fn cross_contact_response_uses_the_full_nondiagonal_original_mass() {
    let bodies = [Body {
        mass: &Coupled,
        velocity: &[-1.0, -2.0],
        preconditioner: &IdentityPrecond,
    }];
    let first = [TraceTerm {
        body: 0,
        dof: 0,
        coefficient: 1.0,
    }];
    let second = [TraceTerm {
        body: 0,
        dof: 1,
        coefficient: 1.0,
    }];
    let contacts = [row(4, &first), row(9, &second)];
    let plan = prepare(&bodies, &contacts, controls(), &mut || false).unwrap();
    // Direct 2x2 inverse, and -M*v for the simultaneous zero-velocity endpoint.
    for (&actual, expected) in
        plan.mobility
            .matrix
            .iter()
            .zip([4.0 / 7.0, -2.0 / 7.0, -2.0 / 7.0, 8.0 / 7.0])
    {
        close(actual, expected);
    }
    close(plan.body_impulses[0][0], 3.0);
    close(plan.body_impulses[0][1], 2.5);
    let report = finish(
        plan,
        &contacts,
        &[&[0.0, 0.0]],
        &[-4.0],
        &[0.0],
        controls(),
        &mut || false,
    )
    .unwrap();
    assert_eq!(report.mass_solves.len(), 2);
    close(report.contact.dissipated_energy, 4.0);
    close(report.work_energy_residual_j, 0.0);
}

#[test]
fn moving_support_supplies_work_even_when_impact_dissipates_relative_energy() {
    let bodies = [Body {
        mass: &Scalar(2.0),
        velocity: &[0.0],
        preconditioner: &IdentityPrecond,
    }];
    let trace = [TraceTerm {
        body: 0,
        dof: 0,
        coefficient: 1.0,
    }];
    let mut contact = row(2, &trace);
    contact.prescribed_normal_velocity_m_s = -3.0;
    let plan = prepare(&bodies, &[contact], controls(), &mut || false).unwrap();
    close(plan.body_impulses[0][0], 6.0);
    let report = finish(
        plan,
        &[contact],
        &[&[3.0]],
        &[9.0],
        &[0.0],
        controls(),
        &mut || false,
    )
    .unwrap();
    close(report.relative_velocities_after_m_s[0], 0.0);
    close(report.prescribed_work_j, 18.0);
    close(report.relative_impulse_work_j, -9.0);
    close(report.kinetic_change_j, 9.0);
    close(report.contact.dissipated_energy, 9.0);
    close(report.work_energy_residual_j, 0.0);
}

#[test]
fn redundant_contact_rows_keep_native_velocity_determinacy_without_load_authority() {
    let bodies = [Body {
        mass: &Scalar(2.0),
        velocity: &[-3.0],
        preconditioner: &IdentityPrecond,
    }];
    let trace = [TraceTerm {
        body: 0,
        dof: 0,
        coefficient: 1.0,
    }];
    let contacts = [row(1, &trace), row(2, &trace)];
    let plan = prepare(&bodies, &contacts, controls(), &mut || false).unwrap();
    close(plan.body_impulses[0][0], 6.0);
    let report = finish(
        plan,
        &contacts,
        &[&[0.0]],
        &[-9.0],
        &[0.0],
        controls(),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        report.contact.determinacy,
        MultiImpactDeterminacy::UniqueVelocity(_)
    ));
    close(report.contact.normal_impulses.iter().sum(), 6.0);
    close(report.contact.dissipated_energy, 9.0);
}

struct Counted(Cell<usize>);
impl LinearOp for Counted {
    fn n(&self) -> usize {
        1
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.set(self.0.get() + 1);
        y[0] = 2.0 * x[0];
    }
    // A counting 1x1 scalar, so its transpose is itself; through apply, the count is unchanged.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

#[test]
fn admission_and_actual_endpoint_failures_return_no_accepted_response() {
    let mass = Counted(Cell::new(0));
    let bodies = [Body {
        mass: &mass,
        velocity: &[-1.0],
        preconditioner: &IdentityPrecond,
    }];
    let trace = [TraceTerm {
        body: 0,
        dof: 0,
        coefficient: 1.0,
    }];
    let contact = row(3, &trace);
    let mut limited = controls();
    limited.maximum_workspace_entries = 0;
    assert!(matches!(
        prepare(&bodies, &[contact], limited, &mut || false),
        Err(ResponseError::Budget(_))
    ));
    assert!(matches!(
        prepare(&bodies, &[contact], controls(), &mut || true),
        Err(ResponseError::Cancelled)
    ));
    assert!(matches!(
        prepare(&bodies, &[contact, contact], controls(), &mut || false),
        Err(ResponseError::Invalid(_))
    ));
    let duplicate = [trace[0], trace[0]];
    assert!(matches!(
        prepare(&bodies, &[row(4, &duplicate)], controls(), &mut || false),
        Err(ResponseError::Invalid(_))
    ));
    assert_eq!(
        mass.0.get(),
        0,
        "admission must precede original operator work"
    );
    let plan = prepare(&bodies, &[contact], controls(), &mut || false).unwrap();
    assert!(matches!(
        finish(
            plan,
            &[contact],
            &[&[0.1]],
            &[-0.99],
            &[0.0],
            controls(),
            &mut || false
        ),
        Err(ResponseError::Audit {
            quantity: "native endpoint velocity [m/s]",
            ..
        })
    ));
    let plan = prepare(&bodies, &[contact], controls(), &mut || false).unwrap();
    assert!(matches!(
        finish(
            plan,
            &[contact],
            &[&[0.0]],
            &[-1.0],
            &[1e-4],
            controls(),
            &mut || false
        ),
        Err(ResponseError::Audit {
            quantity: "absolute native body work [J]",
            ..
        })
    ));
}

#[test]
fn three_body_chain_resolves_shared_cross_terms_for_elastic_and_inelastic_impacts() {
    let masses = [Scalar(1.0), Scalar(2.0), Scalar(3.0)];
    let velocities = [[3.0], [0.0], [-1.0]];
    let bodies: [Body<'_>; 3] = core::array::from_fn(|body| Body {
        mass: &masses[body],
        velocity: &velocities[body],
        preconditioner: &IdentityPrecond,
    });
    let first = [
        TraceTerm {
            body: 0,
            dof: 0,
            coefficient: -1.0,
        },
        TraceTerm {
            body: 1,
            dof: 0,
            coefficient: 1.0,
        },
    ];
    let second = [
        TraceTerm {
            body: 1,
            dof: 0,
            coefficient: -1.0,
        },
        TraceTerm {
            body: 2,
            dof: 0,
            coefficient: 1.0,
        },
    ];
    for restitution in [0.0, 1.0] {
        let mut contacts = [row(10, &first), row(20, &second)];
        for contact in &mut contacts {
            contact.restitution = restitution;
        }
        let plan = prepare(&bodies, &contacts, controls(), &mut || false).unwrap();
        // W=[[3/2,-1/2],[-1/2,5/6]], det(W)=1. Hence both
        // simultaneous impulses are 3*(1+e); separate pairwise solves differ.
        let magnitude = 3.0 * (1.0 + restitution);
        close(plan.body_impulses[0][0], -magnitude);
        close(plan.body_impulses[1][0], 0.0);
        close(plan.body_impulses[2][0], magnitude);
        let after: [[f64; 1]; 3] = core::array::from_fn(|body| {
            [velocities[body][0] + plan.body_impulses[body][0] / masses[body].0]
        });
        close(after[0][0], -3.0 * restitution);
        close(after[1][0], 0.0);
        close(after[2][0], restitution);
        let changes: [f64; 3] = core::array::from_fn(|body| {
            0.5 * masses[body].0 * (after[body][0].powi(2) - velocities[body][0].powi(2))
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
        for (&actual, expected) in
            report
                .contact_mobility_per_kg
                .iter()
                .zip([1.5, -0.5, -0.5, 5.0 / 6.0])
        {
            close(actual, expected);
        }
        close(report.contact.normal_impulses[0], magnitude);
        close(report.contact.normal_impulses[1], magnitude);
        close(
            report.kinetic_change_j,
            -6.0 * (1.0 - restitution * restitution),
        );
        close(
            report.contact.dissipated_energy,
            6.0 * (1.0 - restitution * restitution),
        );
        close(report.work_energy_residual_j, 0.0);
        close(report.dissipation_residual_j, 0.0);
        close(after[0][0] + 2.0 * after[1][0] + 3.0 * after[2][0], 0.0);
    }
}
