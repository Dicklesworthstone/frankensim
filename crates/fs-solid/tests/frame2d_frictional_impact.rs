//! G0/G1/G3: actual planar Coulomb impulses, original inertia and continuation.
use fs_contact::frictional_multicontact::{
    FrictionalContactConfig, FrictionalContactDeterminacy, FrictionalContactError,
    NormalVelocityLaw,
};
use fs_contact::multicontact::MultiImpactConfig;
use fs_solid::frame2d::{
    FrameFrictionalContact, FrameFrictionalImpact, FrameFrictionalImpactError,
    FrameMassFormulation, FrameNormalContact,
};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

const UNDAMPED: Rayleigh = Rayleigh { a0: 0.0, a1: 0.0 };
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance={tolerance:.3e}"
    );
}
fn vector_close(a: &[f64], b: &[f64], tolerance: f64) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        close(a, b, tolerance);
    }
}
fn contact<'a>(id: u64, n: &'a [f64], t: &'a [f64], mu: f64) -> FrameFrictionalContact<'a> {
    FrameFrictionalContact {
        id,
        normal_covector: n,
        tangent_covector: t,
        gap: 0.0,
        activation_tolerance: 0.0,
        obstacle_velocity: [0.0; 2],
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}
fn event(
    frame: &mut Frame2d,
    contacts: &[FrameFrictionalContact<'_>],
) -> Result<FrameFrictionalImpact, FrameFrictionalImpactError> {
    let load = vec![0.0; frame.displacements().len()];
    frame.apply_frictional_impacts(
        contacts,
        &load,
        UNDAMPED,
        0.01,
        FrictionalContactConfig::default(),
        || false,
    )
}
fn seed(frame: &mut Frame2d, impulse: &[f64]) {
    let load = vec![0.0; impulse.len()];
    frame.apply_impulse(impulse, &load, UNDAMPED, 0.01).unwrap();
}
fn point_mass(mass: f64, velocity: [f64; 2]) -> Frame2d {
    let mut f = Frame2d::new(Geometry::Linear);
    f.add_node(0.0, 0.0);
    f.add_mass(0, mass);
    f.fix(0, [false, false, true]);
    if velocity != [0.0; 2] {
        seed(&mut f, &[mass * velocity[0], mass * velocity[1], 0.0]);
    }
    f
}
fn elastic() -> SectionModel {
    SectionModel::Elastic {
        ea: 4000.0,
        ei: 20.0,
    }
}
fn rotary_frame() -> Frame2d {
    let mut f = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::Lumped);
    f.add_node(0.0, 0.0);
    f.add_node(1.0, 0.0);
    f.fix(0, [true; 3]);
    // Native lumped endpoint inventory: m=2 in both translations, I=1.
    f.add_element_with_rotary_inertia(0, 1, 5, 4.0, 2.0, &elastic)
        .unwrap();
    f
}
fn coupled_frame() -> Frame2d {
    let mut f =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    f.add_node(0.0, 0.0);
    f.add_node(1.0, 0.0);
    f.fix(0, [false, false, true]);
    f.fix(1, [false, false, true]);
    f.add_element(0, 1, 5, 6.0, &elastic).unwrap();
    // M_x=[[2,1],[1,2]], M_y=[[78,27],[27,78]]/35.
    // Independently seed v_x=[-1,1/4], v_y=[3,-2].
    seed(&mut f, &[-1.75, 36.0 / 7.0, 0.0, -0.5, -15.0 / 7.0, 0.0]);
    f
}
fn planar_receipt(result: &FrameFrictionalImpact) {
    assert!(matches!(
        result.contact.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    assert!(result.contact.impulse_work_witness.is_some());
    for (p, v) in result
        .contact
        .impulses
        .iter()
        .zip(&result.contact.relative_velocities_after)
    {
        assert_eq!(p[2], 0.0);
        assert_eq!(v[2], 0.0);
    }
}

#[test]
fn g1_sliding_point_mass_matches_the_independent_coulomb_impulse_and_energy() {
    let mut f = point_mass(2.0, [-1.0, 3.0]);
    let n = [1.0, 0.0, 0.0];
    let t = [0.0, 1.0, 0.0];
    let result = event(&mut f, &[contact(1, &n, &t, 0.5)]).unwrap();
    vector_close(&result.contact.impulses[0], &[2.0, -1.0, 0.0], 2e-9);
    vector_close(f.velocities(), &[0.0, 2.5, 0.0], 2e-9);
    close(result.impulse.kinetic_energy_change_j, -3.75, 3e-9);
    close(result.dissipated_energy_j, 3.75, 3e-9);
    close(result.obstacle_work_j, 0.0, 0.0);
    close(result.residuals.sliding_impulse_residual, 0.0, 2e-9);
    planar_receipt(&result);
}

#[test]
fn g1_off_center_sliding_couples_normal_impulse_friction_and_nodal_moment() {
    let mut f = rotary_frame();
    seed(&mut f, &[0.0, 0.0, 0.0, 6.0, -2.0, 0.0]);
    // Contact offset r=(1,1), normal=+y, tangent=+x:
    // v_n=v_y+omega, v_t=v_x-omega. W=[[3/2,-1],[-1,3/2]].
    let n = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0];
    let t = [0.0, 0.0, 0.0, 1.0, 0.0, -1.0];
    let result = event(&mut f, &[contact(1, &n, &t, 0.25)]).unwrap();
    // Simultaneous sliding: p_t=-p_n/4 and (3/2+1/4)*p_n=1.
    vector_close(
        &result.contact.impulses[0],
        &[4.0 / 7.0, -1.0 / 7.0, 0.0],
        2e-9,
    );
    vector_close(
        f.velocities(),
        &[0.0, 0.0, 0.0, 41.0 / 14.0, -5.0 / 7.0, 5.0 / 7.0],
        2e-9,
    );
    close(result.relative_velocities_after[0][0], 0.0, 2e-9);
    close(result.relative_velocities_after[0][1], 31.0 / 14.0, 2e-9);
    close(
        result.impulse.generalized_momentum_change[5],
        5.0 / 7.0,
        2e-9,
    );
    close(result.dissipated_energy_j, 129.0 / 196.0, 3e-9);
    // A normal-only 2/3 kick followed by friction would violate the final
    // normal complementarity equation; this is not sequential clipping.
    assert!((result.contact.impulses[0][0] - 2.0 / 3.0).abs() > 0.05);
    planar_receipt(&result);
}

#[test]
fn g1_coupled_sticking_can_reverse_the_sign_selected_from_initial_slip_alone() {
    let mut f = rotary_frame();
    seed(&mut f, &[0.0, 0.0, 0.0, 0.5, -2.0, 0.0]);
    let n = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0];
    let t = [0.0, 0.0, 0.0, 1.0, 0.0, -1.0];
    let result = event(&mut f, &[contact(1, &n, &t, 0.6)]).unwrap();
    // W^-1*[1,-1/4] gives p=[1,1/2], strictly inside |p_t|<0.6*p_n.
    // The normal kick would otherwise drive the initially positive slip
    // backwards. Contact sticks while the original node still rotates.
    vector_close(&result.contact.impulses[0], &[1.0, 0.5, 0.0], 2e-8);
    vector_close(&result.relative_velocities_after[0], &[0.0, 0.0], 2e-8);
    vector_close(f.velocities(), &[0.0, 0.0, 0.0, 0.5, -0.5, 0.5], 2e-8);
    close(result.dissipated_energy_j, 7.0 / 16.0, 3e-8);
    assert!(result.contact.impulses[0][1].abs() < 0.6 * result.contact.impulses[0][0]);
    planar_receipt(&result);
}

#[test]
fn g1_simultaneous_consistent_contacts_keep_cross_mass_support_moments_and_id_replay() {
    let original = coupled_frame();
    let n0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let t0 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let n1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let t1 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let c0 = contact(20, &n0, &t0, 0.4);
    let c1 = contact(10, &n1, &t1, 0.4);
    let mut f = original.clone();
    let result = event(&mut f, &[c0, c1]).unwrap();
    vector_close(&result.contact.impulses[0], &[1.75, -0.7, 0.0], 2e-8);
    vector_close(&result.contact.impulses[1], &[0.5, 0.2, 0.0], 2e-8);
    vector_close(
        f.velocities(),
        &[0.0, 133.0 / 51.0, 0.0, 0.0, -181.0 / 102.0, 0.0],
        2e-8,
    );
    close(result.impulse.support_impulse[2], -83.0 / 1020.0, 2e-8);
    close(result.impulse.support_impulse[5], 1.0 / 510.0, 2e-8);
    close(result.dissipated_energy_j, 12863.0 / 4080.0, 4e-8);
    assert_eq!(result.contact.active_indices, vec![1, 0]);
    let mut reordered = original;
    let swapped = event(&mut reordered, &[c1, c0]).unwrap();
    assert_eq!(f.velocities(), reordered.velocities());
    assert_eq!(f.accelerations(), reordered.accelerations());
    assert_eq!(result.impulse, swapped.impulse);
    assert_eq!(result.relative_work_j, swapped.relative_work_j);
    assert_eq!(result.contact.impulses[0], swapped.contact.impulses[1]);
    assert_eq!(result.contact.impulses[1], swapped.contact.impulses[0]);
    planar_receipt(&result);
}

#[test]
fn g1_inactive_friction_contact_keeps_its_full_coupled_endpoint_without_receiving_load() {
    let mut f = coupled_frame();
    let n0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let t0 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let n1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let t1 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let mut distant = contact(2, &n1, &t1, 0.4);
    distant.gap = 1.0;
    let result = event(&mut f, &[contact(1, &n0, &t0, 0.4), distant]).unwrap();
    vector_close(&result.contact.impulses[0], &[1.5, -0.6, 0.0], 2e-9);
    assert_eq!(result.contact.impulses[1], [0.0; 3]);
    vector_close(
        &result.relative_velocities_after[1],
        &[-0.25, -161.0 / 85.0],
        2e-9,
    );
    assert_eq!(result.contact.active_indices, vec![0]);
}

#[test]
fn g3_moving_obstacle_work_and_relative_loss_transform_under_a_common_velocity_boost() {
    let n = [1.0, 0.0, 0.0];
    let t = [0.0, 1.0, 0.0];
    let run = |boost: [f64; 2]| {
        let mut f = point_mass(2.0, boost);
        let mut c = contact(1, &n, &t, 0.5);
        c.obstacle_velocity = [1.0 + boost[0], 3.0 + boost[1]];
        let result = event(&mut f, &[c]).unwrap();
        vector_close(&result.contact.impulses[0], &[2.0, 1.0, 0.0], 2e-9);
        vector_close(f.velocities(), &[1.0 + boost[0], 0.5 + boost[1], 0.0], 2e-9);
        let transport = 2.0 * boost[0] + boost[1];
        close(result.obstacle_work_j, 5.0 + transport, 2e-8);
        close(
            result.impulse.kinetic_energy_change_j,
            1.25 + transport,
            2e-8,
        );
        close(result.dissipated_energy_j, 3.75, 2e-8);
        close(
            result.impulse.kinetic_energy_change_j,
            result.obstacle_work_j + result.relative_work_j,
            2e-8,
        );
        result
    };
    let a = run([0.0, 0.0]);
    let b = run([7.0, -5.0]);
    vector_close(&a.contact.impulses[0], &b.contact.impulses[0], 2e-9);
    vector_close(&a.impulse.velocity_change, &b.impulse.velocity_change, 2e-9);
}

#[test]
fn g1_internal_offset_contact_preserves_both_linear_and_reference_angular_momentum() {
    let mut f = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::Lumped);
    f.add_node(0.0, 0.0);
    f.add_node(1.0, 0.0);
    f.add_element_with_rotary_inertia(0, 1, 5, 4.0, 2.0, &elastic)
        .unwrap();
    seed(&mut f, &[2.0, 4.0, 0.0, -2.0, -2.0, 0.0]);
    // Both point traces meet at (1/2,0): r_A=(1/2,0), r_B=(-1/2,0).
    // Relative normal/tangent velocities are point B minus point A.
    let n = [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let t = [0.0, -1.0, -0.5, 0.0, 1.0, -0.5];
    let mut c = contact(1, &n, &t, 0.5);
    c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.25 };
    let result = event(&mut f, &[c]).unwrap();
    vector_close(&result.contact.impulses[0], &[2.5, 1.25, 0.0], 2e-8);
    vector_close(
        f.velocities(),
        &[-0.25, 1.375, -0.625, 0.25, -0.375, -0.625],
        2e-8,
    );
    let v = f.velocities();
    close(2.0 * (v[0] + v[3]), 0.0, 3e-8);
    close(2.0 * (v[1] + v[4]), 2.0, 3e-8);
    // L_z = sum(x*m*v_y - y*m*v_x + I*omega) for this reference inventory.
    close(2.0 * v[4] + v[2] + v[5], -2.0, 3e-8);
    assert_eq!(result.impulse.support_impulse, vec![0.0; 6]);
    close(result.dissipated_energy_j, 285.0 / 64.0, 4e-8);
}

#[test]
fn g1_actual_coulomb_endpoint_sets_finite_force_acceleration_and_newmark_continuation() {
    let mut f = rotary_frame();
    f.static_load(&[0.0, 0.0, 0.0, 8.0, 0.0, 0.0], 1).unwrap();
    let u = f.displacements().to_vec();
    let sections = f.section_deformations(0).to_vec();
    let forces = f.element_forces(0);
    let n = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0];
    let t = [0.0, 0.0, 0.0, 1.0, 0.0, -1.0];
    let mut c = contact(1, &n, &t, 0.25);
    c.obstacle_velocity = [1.0, 3.0];
    let damping = Rayleigh { a0: 0.2, a1: 0.001 };
    let h = 0.001;
    let load = [0.0, 0.0, 0.0, 10.0, 1.0, 2.0];
    let result = f
        .apply_frictional_impacts(
            &[c],
            &load,
            damping,
            h,
            FrictionalContactConfig::default(),
            || false,
        )
        .unwrap();
    vector_close(&result.contact.impulses[0], &[0.8, 0.2, 0.0], 2e-9);
    let v0 = [0.1, 0.4, 0.6];
    let a0 = [0.78, 0.408, 1.88];
    vector_close(&f.velocities()[3..], &v0, 2e-9);
    vector_close(&f.accelerations()[3..], &a0, 2e-8);
    assert_eq!(f.displacements(), u);
    assert_eq!(f.section_deformations(0), sections);
    assert_eq!(f.element_forces(0), forces);
    // Independent fixed-base Euler-Bernoulli M, K and Rayleigh C.
    let m = [2.0, 2.0, 1.0];
    let k = [
        [4000.0, 0.0, 0.0],
        [0.0, 240.0, -120.0],
        [0.0, -120.0, 80.0],
    ];
    let cm = [[4.4, 0.0, 0.0], [0.0, 0.64, -0.12], [0.0, -0.12, 0.28]];
    let mut matrix = [[0.0; 3]; 3];
    let mut rhs = [10.0, 1.0, 2.0];
    for i in 0..3 {
        for j in 0..3 {
            matrix[i][j] =
                (if i == j { m[i] } else { 0.0 }) + 0.5 * h * cm[i][j] + 0.25 * h * h * k[i][j];
            rhs[i] -= k[i][j] * (u[j + 3] + h * v0[j] + 0.25 * h * h * a0[j])
                + cm[i][j] * (v0[j] + 0.5 * h * a0[j]);
        }
    }
    let determinant = matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1];
    let a1 = [
        rhs[0] / matrix[0][0],
        (rhs[1] * matrix[2][2] - matrix[1][2] * rhs[2]) / determinant,
        (matrix[1][1] * rhs[2] - rhs[1] * matrix[2][1]) / determinant,
    ];
    let history = f
        .newmark_with_load(&[0.0, 0.0], h, (0.0, 0.0), damping, &load)
        .unwrap();
    close(
        history.kinetic_energy[0],
        result.impulse.kinetic_energy_after_j,
        2e-9,
    );
    for j in 0..3 {
        close(
            f.displacements()[j + 3],
            u[j + 3] + h * v0[j] + 0.25 * h * h * (a0[j] + a1[j]),
            2e-11,
        );
        close(
            f.velocities()[j + 3],
            v0[j] + 0.5 * h * (a0[j] + a1[j]),
            2e-9,
        );
        close(f.accelerations()[j + 3], a1[j], 2e-7);
    }
}

#[test]
fn g3_zero_friction_agrees_with_the_existing_native_normal_event() {
    let original = point_mass(2.0, [-1.0, 3.0]);
    let n = [1.0, 0.0, 0.0];
    let t = [0.0, 1.0, 0.0];
    let mut f = original.clone();
    let mut c = contact(1, &n, &t, 0.0);
    c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 };
    let coupled = event(&mut f, &[c]).unwrap();
    let mut normal = original;
    let existing = normal
        .apply_normal_impacts(
            &[FrameNormalContact {
                normal_covector: &n,
                gap: 0.0,
                activation_tolerance: 0.0,
                obstacle_normal_speed: 0.0,
                restitution: 0.5,
            }],
            &[0.0; 3],
            UNDAMPED,
            0.01,
            MultiImpactConfig::default(),
            || false,
        )
        .unwrap();
    vector_close(f.velocities(), normal.velocities(), 2e-9);
    vector_close(f.accelerations(), normal.accelerations(), 2e-9);
    close(
        coupled.dissipated_energy_j,
        existing.dissipated_energy_j,
        2e-9,
    );
    close(coupled.contact.impulses[0][1], 0.0, 0.0);
}

#[test]
fn g3_redundant_frictional_contacts_retain_qualified_motion_and_total_work() {
    let mut f = point_mass(2.0, [-1.0, 3.0]);
    let n = [1.0, 0.0, 0.0];
    let t = [0.0, 1.0, 0.0];
    let result = event(&mut f, &[contact(1, &n, &t, 0.5), contact(2, &n, &t, 0.5)]).unwrap();
    close(
        result.contact.impulses.iter().map(|p| p[0]).sum(),
        2.0,
        2e-8,
    );
    close(
        result.contact.impulses.iter().map(|p| p[1]).sum(),
        -1.0,
        2e-8,
    );
    vector_close(f.velocities(), &[0.0, 2.5, 0.0], 2e-8);
    close(result.dissipated_energy_j, 3.75, 3e-8);
    planar_receipt(&result);
}

#[test]
fn g0_global_restitution_energy_failure_is_atomic_and_not_repaired_by_friction() {
    let mut f = coupled_frame();
    let before = format!("{f:?}");
    let n0 = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let t0 = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0];
    let n1 = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let t1 = [0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let mut cs = [contact(1, &n0, &t0, 0.0), contact(2, &n1, &t1, 0.0)];
    for c in &mut cs {
        c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 1.0 };
    }
    match event(&mut f, &cs).unwrap_err() {
        FrameFrictionalImpactError::Contact(FrictionalContactError::EnergyIncrease {
            work,
            ..
        }) => close(work, 3.0 / 16.0, 2e-8),
        other => panic!("wrong refusal: {other:?}"),
    }
    assert_eq!(format!("{f:?}"), before);
}

#[test]
fn g0_late_cancellation_and_evaluation_exhaustion_preserve_the_full_frame_before_retry() {
    let original = point_mass(2.0, [0.0, 0.0]); // No Rayleigh reference captured yet.
    let n = [1.0, 0.0, 0.0];
    let t = [0.0, 1.0, 0.0];
    let mut c = contact(1, &n, &t, 0.5);
    c.obstacle_velocity = [1.0, 3.0];
    let mut expected = original.clone();
    let mut polls = 0;
    let result = expected
        .apply_frictional_impacts(
            &[c],
            &[0.0; 3],
            UNDAMPED,
            0.01,
            FrictionalContactConfig::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 10);
    for cancel_at in [1, 5, polls] {
        let mut f = original.clone();
        let before = format!("{f:?}");
        let mut count = 0;
        let failed = f.apply_frictional_impacts(
            &[c],
            &[0.0; 3],
            UNDAMPED,
            0.01,
            FrictionalContactConfig::default(),
            || {
                count += 1;
                count == cancel_at
            },
        );
        assert!(matches!(
            failed,
            Err(FrameFrictionalImpactError::Contact(
                FrictionalContactError::Cancelled
            ))
        ));
        assert_eq!(format!("{f:?}"), before);
        let retried = event(&mut f, &[c]).unwrap();
        assert_eq!(retried, result);
        assert_eq!(format!("{f:?}"), format!("{expected:?}"));
    }
    let mut f = original;
    let before = format!("{f:?}");
    let config = FrictionalContactConfig {
        maximum_evaluations: 1,
        ..Default::default()
    };
    assert!(matches!(
        f.apply_frictional_impacts(&[c], &[0.0; 3], UNDAMPED, 0.01, config, || false),
        Err(FrameFrictionalImpactError::Contact(
            FrictionalContactError::EvaluationBudget
        ))
    ));
    assert_eq!(format!("{f:?}"), before);
}

#[test]
fn g0_invalid_geometry_limits_and_late_continuation_context_leave_the_frame_unchanged() {
    let original = point_mass(2.0, [0.0, 0.0]);
    let n = [1.0, 0.0, 0.0];
    let t = [0.0, 1.0, 0.0];
    for case in 0..9 {
        let mut f = original.clone();
        let before = format!("{f:?}");
        let mut c = contact(1, &n, &t, 0.5);
        c.obstacle_velocity = [1.0, 3.0];
        let mut config = FrictionalContactConfig::default();
        let mut dt = 0.01;
        match case {
            0 => c.gap = -1e-6,
            1 => c.normal_covector = &[1.0],
            2 => c.tangent_covector = &[0.0, f64::NAN, 0.0],
            3 => c.obstacle_velocity[1] = f64::INFINITY,
            4 => c.friction_coefficient = -0.1,
            5 => c.normal_law = NormalVelocityLaw::Target { speed: -1.0 },
            6 => config.maximum_contacts = 0,
            7 => config.energy_tolerance = f64::NAN,
            _ => dt = 0.0,
        }
        assert!(
            f.apply_frictional_impacts(&[c], &[0.0; 3], UNDAMPED, dt, config, || false)
                .is_err()
        );
        assert_eq!(format!("{f:?}"), before);
    }
    let mut f = original;
    let before = format!("{f:?}");
    let c = contact(1, &n, &t, 0.5);
    assert!(event(&mut f, &[c, c]).is_err());
    assert_eq!(format!("{f:?}"), before);
    let mut massless = Frame2d::new(Geometry::Linear);
    massless.add_node(0.0, 0.0);
    massless.add_mass(0, 2.0); // Free zero-inertia theta.
    let before = format!("{massless:?}");
    assert!(event(&mut massless, &[c]).is_err());
    assert_eq!(format!("{massless:?}"), before);
}

#[test]
fn g0_native_velocity_quantization_cannot_erase_or_change_a_coulomb_kick() {
    let mut f = point_mass(1.0, [1e16, 3.0]);
    let before = format!("{f:?}");
    let n = [1.0, 0.0, 0.0];
    let t = [0.0, 1.0, 0.0];
    let mut c = contact(1, &n, &t, 0.5);
    c.obstacle_velocity = [1e16 + 2.0, 0.0];
    c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 };
    // Normal delta-v=3 is not representable at this stored spacing of 2.
    assert!(matches!(
        event(&mut f, &[c]),
        Err(FrameFrictionalImpactError::Frame(_))
    ));
    assert_eq!(format!("{f:?}"), before);
}
