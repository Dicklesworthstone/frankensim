//! Original cross-mass reactions and work for simultaneous supported Coulomb events.

use super::*;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_solid::dynamics::frictional_impact::{
    SolidFrictionalContact, SolidFrictionalImpactControls,
};
use fs_solid::dynamics::supported::{SupportedSolidError, SupportedSolidFrictionalImpactError};

const ROWS: [[f64; 4]; 3] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, -1.0],
    [0.0, 0.0, 1.0, 0.0],
];
const U: [f64; 4] = [0.05, -0.03, 0.07, 0.02];

fn contact(mu: f64, surface: [f64; 3]) -> SolidFrictionalContact<'static> {
    SolidFrictionalContact {
        id: 1,
        velocity_weights: [&ROWS[0], &ROWS[1], &ROWS[2]],
        gap_m: 0.0,
        activation_tolerance_m: 1e-12,
        prescribed_surface_velocity_mps: surface,
        friction_coefficient: mu,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
    }
}

fn impact_controls() -> SolidFrictionalImpactControls {
    SolidFrictionalImpactControls {
        contact: FrictionalContactConfig {
            velocity_tolerance: 1e-10,
            impulse_tolerance: 1e-10,
            complementarity_tolerance: 1e-10,
            energy_tolerance: 1e-10,
            ..Default::default()
        },
        maximum_impulse_error_kg_m_s: 1e-10,
        maximum_contact_speed_error_mps: 1e-10,
        maximum_work_error_j: 1e-10,
        maximum_workspace_bytes: 1 << 20,
        maximum_response_multiply_terms: 1_000_000,
    }
}

fn assembly() -> TetElasticAssembly {
    let (mut mass, mut stiffness) = (Coo::new(4, 4), Coo::new(4, 4));
    let matrix = [
        [2.0, 0.0, 0.0, 0.35],
        [0.0, 3.0, 0.0, -0.2],
        [0.0, 0.0, 4.0, 0.1],
        [0.35, -0.2, 0.1, 1.5],
    ];
    for (i, row) in matrix.iter().enumerate() {
        for (j, &entry) in row.iter().enumerate() {
            mass.push(i, j, entry);
        }
    }
    for (i, k) in [18.0, 27.0, 36.0].into_iter().enumerate() {
        stiffness.push(i, i, k);
        stiffness.push(i, 3, -k);
        stiffness.push(3, i, -k);
        stiffness.push(3, 3, k);
    }
    TetElasticAssembly {
        stiffness: stiffness.assemble(),
        mass: mass.assemble(),
        free_dofs: vec![0, 1, 2, 3],
        total_mass_kg: 11.0,
        minimum_scaled_jacobian: 1.0,
        element_volumes_m3: vec![1.0],
    }
}

fn make(
    v: [f64; 4],
    basis_scale: f64,
    limits: SupportedSolidControls,
    cx: &Cx<'_>,
) -> SupportedSpectralSolid {
    SupportedSpectralSolid::prepare_with_basis(
        &assembly(),
        &[3],
        &[9.0; 3],
        &[
            basis_scale / 2.0_f64.sqrt(),
            0.0,
            0.0,
            0.0,
            1.0 / 3.0_f64.sqrt(),
            0.0,
            0.0,
            0.0,
            0.5,
        ],
        0.01,
        &U,
        &v,
        &[0.0],
        limits,
        cx,
    )
    .unwrap()
}

fn kinetic(v: &[f64]) -> f64 {
    v[0] * v[0]
        + 1.5 * v[1] * v[1]
        + 2.0 * v[2] * v[2]
        + 0.75 * v[3] * v[3]
        + v[3] * (0.35 * v[0] - 0.2 * v[1] + 0.1 * v[2])
}

#[test]
fn g1_sliding_event_keeps_full_support_reaction_and_native_elastic_continuation() {
    with_cx(false, |cx| {
        let before = [-1.0, 2.4, 0.0, 0.4];
        let mut solid = make(before, 1.0, controls(), cx);
        let report = solid
            .apply_frictional_impacts(&[contact(0.25, [0.0; 3])], impact_controls(), cx)
            .unwrap();
        // W=diag(1/2,1/3,1/4). Restitution gives p_n=3 and sliding gives
        // p_t=-mu*p_n=-.75; the support receives the opposite tangential load.
        let expected = [0.5, 2.15, 0.0, 0.4];
        for i in 0..4 {
            near(solid.velocity_mps()[i], expected[i], 2e-10);
        }
        for (actual, expected) in report.contact.law.impulses[0].iter().zip([3.0, -0.75, 0.0]) {
            near(*actual, expected, 2e-10);
        }
        // M_cf Delta-v = .35*1.5 + (-.2)*(-.25) = .575, while J_c=.75.
        near(report.impulse.momentum_change_kg_m_s[3], 0.575, 2e-10);
        near(report.impulse.support_reaction_kg_m_s[0], -0.175, 2e-10);
        near(report.impulse.energy.support_work_j, -0.07, 2e-10);
        near(report.contact.relative_work_j, -2.15625, 3e-10);
        near(report.contact.dissipated_energy_j, 2.15625, 3e-10);
        near(report.impulse.energy.applied_work_j, -2.15625, 3e-10);
        near(
            kinetic(solid.velocity_mps()) - kinetic(&before),
            -2.22625,
            3e-10,
        );
        assert_eq!(solid.displacement_m(), &U);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(solid.accepted_steps(), 0);
        assert_eq!(solid.velocity_mps()[3].to_bits(), before[3].to_bits());

        // Every free coordinate has omega=3 relative to the translating
        // support. This independent analytical trajectory starts after contact.
        solid
            .advance_affine_supports(&[0.0], &[0.0; 4], &[0.0; 4], cx)
            .unwrap();
        let (s, c) = 0.03_f64.sin_cos();
        let uc = U[3] + 0.01 * before[3];
        for i in 0..3 {
            near(
                solid.displacement_m()[i],
                uc + (U[i] - U[3]) * c + (expected[i] - before[3]) * s / 3.0,
                3e-11,
            );
            near(
                solid.velocity_mps()[i],
                before[3] - 3.0 * (U[i] - U[3]) * s + (expected[i] - before[3]) * c,
                3e-10,
            );
        }
        near(solid.cumulative_applied_work_j(), -2.15625, 3e-10);
    });
}

#[test]
fn g1_sticking_uses_actual_moving_support_velocity_instead_of_zero_slip_frame() {
    with_cx(false, |cx| {
        let mut solid = make([-1.0, 2.4, 0.0, 0.4], 1.0, controls(), cx);
        let report = solid
            .apply_frictional_impacts(&[contact(4.0, [0.0; 3])], impact_controls(), cx)
            .unwrap();
        near(report.contact.law.impulses[0][0], 3.0, 2e-10);
        near(report.contact.law.impulses[0][1], -6.0, 2e-10);
        near(solid.velocity_mps()[1], 0.4, 2e-10);
        near(
            report.contact.relative_velocities_after_mps[0][1],
            0.0,
            2e-10,
        );
        near(report.impulse.support_reaction_kg_m_s[0], -5.075, 3e-10);
        near(report.contact.relative_work_j, -6.75, 3e-10);
        near(report.impulse.energy.support_work_j, -2.03, 3e-10);
        near(report.impulse.kinetic_change_j, -8.78, 5e-10);
    });
}

#[test]
fn g3_obstacle_and_support_work_transform_separately_under_a_common_boost() {
    with_cx(false, |cx| {
        let mut results = Vec::new();
        for boost in [0.0, 2.0] {
            let surface = [0.25 + boost, -0.1, 0.2 + boost];
            let vc = 0.4 + boost;
            let before = [surface[0] - 1.0, vc + surface[1] + 2.0, surface[2], vc];
            let mut solid = make(before, 1.0, controls(), cx);
            let report = solid
                .apply_frictional_impacts(&[contact(0.25, surface)], impact_controls(), cx)
                .unwrap();
            near(report.contact.relative_work_j, -2.15625, 5e-10);
            near(
                report.contact.prescribed_surface_work_j,
                0.825 + 3.0 * boost,
                5e-10,
            );
            near(
                report.impulse.energy.applied_work_j,
                -1.33125 + 3.0 * boost,
                5e-10,
            );
            near(report.impulse.energy.support_work_j, -0.175 * vc, 5e-10);
            near(
                kinetic(solid.velocity_mps()) - kinetic(&before),
                report.impulse.energy.applied_work_j + report.impulse.energy.support_work_j,
                8e-10,
            );
            results.push(report);
        }
        for i in 0..3 {
            near(
                results[0].contact.law.impulses[0][i],
                results[1].contact.law.impulses[0][i],
                2e-10,
            );
        }
        near(
            results[1].impulse.energy.support_work_j - results[0].impulse.energy.support_work_j,
            -0.35,
            5e-10,
        );
    });
}

#[test]
fn g4_refused_frictional_candidates_leave_both_supported_owners_retryable() {
    with_cx(false, |cx| {
        let initial = [-1.0, 2.4, 0.0, 0.4];
        let event = [contact(0.25, [0.0; 3])];
        let mut retry = make(initial, 1.0, controls(), cx);
        let mut reference = make(initial, 1.0, controls(), cx);
        let original = bits(&retry);
        with_cx(true, |cancelled| {
            assert!(
                retry
                    .apply_frictional_impacts(&event, impact_controls(), cancelled)
                    .is_err()
            );
        });
        assert_eq!(bits(&retry), original);
        for budget in 0..2 {
            let mut limits = impact_controls();
            if budget == 0 {
                limits.maximum_workspace_bytes = 1;
            } else {
                limits.maximum_response_multiply_terms = 1;
            }
            assert!(retry.apply_frictional_impacts(&event, limits, cx).is_err());
            assert_eq!(bits(&retry), original);
        }
        // An approaching contact on the prescribed coordinate alone has no
        // normal mobility. Tangential free coordinates cannot repair it.
        let locked = SolidFrictionalContact {
            velocity_weights: [&[0.0, 0.0, 0.0, 1.0], &ROWS[0], &ROWS[2]],
            prescribed_surface_velocity_mps: [1.0, 0.0, 0.0],
            ..event[0]
        };
        assert!(
            retry
                .apply_frictional_impacts(&[locked], impact_controls(), cx)
                .is_err()
        );
        assert_eq!(bits(&retry), original);
        retry
            .apply_frictional_impacts(&event, impact_controls(), cx)
            .unwrap();
        reference
            .apply_frictional_impacts(&event, impact_controls(), cx)
            .unwrap();
        retry
            .advance_affine_supports(&[0.0], &[0.0; 4], &[0.0; 4], cx)
            .unwrap();
        reference
            .advance_affine_supports(&[0.0], &[0.0; 4], &[0.0; 4], cx)
            .unwrap();
        assert_eq!(bits(&retry), bits(&reference));

        // The shared law and free impulse may pass with an admitted rounded
        // basis; the final original full-energy gate must still reject it.
        let mut strict = controls();
        strict.spectral.maximum_relative_operator_error = 0.02;
        strict.spectral.maximum_relative_state_error = 0.02;
        strict.spectral.maximum_work_error_j = 1.0;
        strict.maximum_work_error_j = 1e-16;
        let mut inaccurate = make(initial, 1.001, strict, cx);
        let before = bits(&inaccurate);
        let mut loose = impact_controls();
        loose.maximum_impulse_error_kg_m_s = 1.0;
        assert!(matches!(
            inaccurate.apply_frictional_impacts(&event, loose, cx),
            Err(SupportedSolidFrictionalImpactError::Solid(
                SupportedSolidError::Energy { .. }
            ))
        ));
        assert_eq!(bits(&inaccurate), before);
    });
}
