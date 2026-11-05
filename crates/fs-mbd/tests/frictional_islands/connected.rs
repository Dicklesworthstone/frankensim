//! Exact rigid-body balances for one connected contact graph beyond 21 rows.
use super::*;
use fs_contact::frictional_multicontact::FrictionalContactDeterminacy;

const COUNT: usize = 32;
const RADIUS: f64 = 0.5;
const INERTIA: f64 = 0.1;
const CLOSING_SPEED: f64 = 0.5;

fn stack(
    mu: f64,
) -> (
    Vec<RigidBodyState>,
    Vec<MassProperties>,
    Vec<RigidFrictionalContact>,
) {
    let property =
        MassProperties::new(1.0, Vec3::ZERO, Vec3::new(INERTIA, INERTIA, INERTIA)).unwrap();
    let bodies = (0..COUNT)
        .map(|i| {
            body(
                Vec3::new(0.0, 0.0, RADIUS + i as f64),
                Vec3::new(1.0, 0.0, -0.125 - CLOSING_SPEED),
            )
        })
        .collect();
    let contacts = (0..COUNT)
        .map(|i| RigidFrictionalContact {
            id: 100 + i as u64,
            point: RigidImpactPoint {
                body: i,
                arm_body: Vec3::new(0.0, 0.0, -RADIUS),
            },
            other: if i == 0 {
                RigidImpactCounterpart::PrescribedPlane {
                    point_world: Vec3::ZERO,
                    velocity_world: Vec3::new(0.25, 0.0, -0.125),
                }
            } else {
                RigidImpactCounterpart::Body(RigidImpactPoint {
                    body: i - 1,
                    arm_body: Vec3::new(0.0, 0.0, RADIUS),
                })
            },
            friction_coefficient: mu,
            ..plane(i, 100 + i as u64)
        })
        .collect();
    (bodies, vec![property; COUNT], contacts)
}

fn large_controls() -> RigidImpactControls {
    let mut controls = controls();
    controls.contact.maximum_contacts = COUNT;
    controls.contact.matrix.max_variables = 3 * COUNT;
    // Budget the larger coupled Krylov solve without changing its accuracy.
    controls.contact.newton.linear_restart = 64;
    controls.contact.newton.max_linear_cycles = 16;
    controls
}

fn close_vector(actual: Vec3, expected: Vec3) {
    close(actual.x, expected.x);
    close(actual.y, expected.y);
    close(actual.z, expected.z);
}

#[test]
fn g3_large_connected_stack_transfers_impulse_spin_and_moving_support_work() {
    for mu in [0.0, 0.25] {
        let (initial, properties, contacts) = stack(mu);
        let mut states = initial.clone();
        let report = apply_rigid_contact_islands(
            &mut states,
            &properties,
            &contacts,
            large_controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(report.contact.islands.len(), 1);
        let island = &report.contact.islands[0];
        assert_eq!(island.body_indices, (0..COUNT).collect::<Vec<_>>());
        assert_eq!(island.contact_indices, (0..COUNT).collect::<Vec<_>>());
        assert!(matches!(
            &island.impact.contact.determinacy,
            FrictionalContactDeterminacy::UniqueImpulses
        ));

        // Momentum balance above each interface fixes p_n=(N-i)*m*closing.
        // For sticking tangential motion, eliminate actual body translation
        // and spin. The independent recurrence is b*p_(i-1)+2*a*p_i+
        // b*p_(i+1)=0, with p_N=0 and a*p_0+b*p_1=-(v_x-U_x).
        // a=1/m+R²/I, b=-1/m+R²/I. Its decaying root gives every traction,
        // including alternating internal signs; no production solver is used.
        let a = 1.0 + RADIUS * RADIUS / INERTIA;
        let b = -1.0 + RADIUS * RADIUS / INERTIA;
        let root = (-a + (a * a - b * b).sqrt()) / b;
        let incoming_slip = 0.75;
        let amplitude = if mu == 0.0 {
            0.0
        } else {
            -incoming_slip
                / (a * (1.0 - root.powi(2 * COUNT as i32))
                    + b * (root - root.powi(2 * COUNT as i32 - 1)))
        };
        let tangent: Vec<_> = (0..=COUNT)
            .map(|i| amplitude * (root.powi(i as i32) - root.powi((2 * COUNT - i) as i32)))
            .collect();
        if mu > 0.0 {
            assert!(tangent[0] < -0.2 && tangent[1] > 0.05);
        }
        for i in 0..COUNT {
            let normal = (COUNT - i) as f64 * CLOSING_SPEED;
            let impulse = report.contact.impulses[i];
            close(impulse[0], normal);
            close(impulse[1], tangent[i]);
            close(impulse[2], 0.0);
            assert!(tangent[i].abs() <= mu * normal);
            let momentum = states[i].linear_momentum_world();
            close(momentum.x, 1.0 + tangent[i] - tangent[i + 1]);
            close(momentum.y, 0.0);
            close(momentum.z, -0.125);
            let spin = states[i].angular_momentum_body();
            close(spin.x, 0.0);
            close(spin.y, -RADIUS * (tangent[i] + tangent[i + 1]));
            close(spin.z, 0.0);
            assert_eq!(states[i].pose(), initial[i].pose());
        }
        let normal_total = COUNT as f64 * CLOSING_SPEED;
        let loss =
            0.5 * COUNT as f64 * CLOSING_SPEED * CLOSING_SPEED - 0.5 * incoming_slip * tangent[0];
        let support_work = 0.25 * tangent[0] - 0.125 * normal_total;
        close_vector(
            report.support_impulse_ns,
            Vec3::new(tangent[0], 0.0, normal_total),
        );
        close(report.support_work_j, support_work);
        close(report.contact.dissipated_energy, loss);
        close(report.kinetic_change_j, support_work - loss);
        close(report.total_energy_residual_j, 0.0);
        close(report.work_energy_residual_j, 0.0);
        close_vector(report.linear_momentum_residual_ns, Vec3::ZERO);
        close_vector(report.angular_momentum_residual_nms, Vec3::ZERO);
        close_vector(report.support_angular_impulse_nms, Vec3::ZERO);

        let mut reversed = contacts.clone();
        reversed.reverse();
        let mut replay = initial;
        let permuted = apply_rigid_contact_islands(
            &mut replay,
            &properties,
            &reversed,
            large_controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(replay, states);
        assert_eq!(permuted.contact_ids, report.contact_ids);
        assert_eq!(permuted.contact.impulses, report.contact.impulses);
    }
}

#[test]
fn g0_large_connected_contacts_require_explicit_budgets_without_partial_state() {
    let (initial, properties, contacts) = stack(0.0);
    let mut states = initial.clone();
    assert!(matches!(
        apply_rigid_contact_islands(
            &mut states,
            &properties,
            &contacts,
            controls(),
            &CancelGate::new(),
        ),
        Err(RigidFrictionalImpactError::Contact(
            IslandContactError::ContactLimit {
                contacts: COUNT,
                maximum: 21,
                ..
            }
        ))
    ));
    assert_eq!(states, initial);
    let mut scalar_limit = large_controls();
    scalar_limit.contact.matrix.max_variables = 3 * COUNT - 1;
    assert!(
        apply_rigid_contact_islands(
            &mut states,
            &properties,
            &contacts,
            scalar_limit,
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(states, initial);
}
