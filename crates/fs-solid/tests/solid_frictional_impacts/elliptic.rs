//! G1/G3/G4 directional impulses through the actual complete solid owners.

#[path = "elliptic/source.rs"]
mod source;

use super::*;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContactError,
};
use fs_solid::dynamics::frictional_impact::elliptic::SolidEllipticContact;

fn directional(
    c: SolidFrictionalContact<'_>,
    mu: [f64; 2],
    axis: [f64; 2],
) -> SolidEllipticContact<'_> {
    SolidEllipticContact {
        id: c.id,
        velocity_weights: c.velocity_weights,
        gap_m: c.gap_m,
        activation_tolerance_m: c.activation_tolerance_m,
        prescribed_surface_velocity_mps: c.prescribed_surface_velocity_mps,
        normal_law: c.normal_law,
        friction_coefficient: EllipticFriction {
            coefficients_mu: mu,
            principal_axis: axis,
        },
    }
}

#[test]
fn isotropic_limit_preserves_both_original_solid_paths_and_physical_accounts() {
    with_cx(false, |cx| {
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let pencil = assembly(&[0.0; 9], &identity, 3);
        let rows = [&identity[0..3], &identity[3..6], &identity[6..9]];
        for mu in [0.0, 0.2, 2.0] {
            let c = contact(19, rows, [0.1, -0.2, 0.3], mu);
            let elliptic = directional(c, [mu; 2], [3e200, 4e200]);
            let make_dense = || {
                ExactLinearSolid::prepare(
                    &pencil,
                    0.02,
                    &[0.01, -0.02, 0.03],
                    &[-1.1, 1.4, -0.7],
                    dense_controls(),
                    cx,
                )
                .unwrap()
            };
            let make_spectral = || {
                CompleteSpectralSolid::prepare_with_basis(
                    &pencil,
                    &[0.0; 3],
                    &identity,
                    0.02,
                    &[0.01, -0.02, 0.03],
                    &[-1.1, 1.4, -0.7],
                    spectral_controls(),
                    cx,
                )
                .unwrap()
            };
            let (mut dense, mut dense_ellipse) = (make_dense(), make_dense());
            let (mut spectral, mut spectral_ellipse) = (make_spectral(), make_spectral());
            let a = dense
                .apply_frictional_impacts(&[c], impact_controls(), cx)
                .unwrap();
            let b = dense_ellipse
                .apply_elliptic_impacts(&[elliptic], impact_controls(), cx)
                .unwrap();
            let c = spectral
                .apply_frictional_impacts(&[c], impact_controls(), cx)
                .unwrap();
            let d = spectral_ellipse
                .apply_elliptic_impacts(&[elliptic], impact_controls(), cx)
                .unwrap();
            assert_eq!(dense.velocity_mps(), dense_ellipse.velocity_mps());
            assert_eq!(dense.displacement_m(), dense_ellipse.displacement_m());
            assert_eq!(bits(&spectral), bits(&spectral_ellipse));
            for (circular, directional) in [(&a, &b), (&c, &d)] {
                assert_eq!(
                    circular.contact.law.impulses,
                    directional.contact.law.impulses
                );
                assert_eq!(
                    circular.contact.relative_work_j.to_bits(),
                    directional.contact.relative_work_j.to_bits(),
                );
                assert_eq!(
                    circular.contact.prescribed_surface_work_j.to_bits(),
                    directional.contact.prescribed_surface_work_j.to_bits(),
                );
                assert_eq!(
                    circular.impulse.energy.supplied_work_j.to_bits(),
                    directional.impulse.energy.supplied_work_j.to_bits(),
                );
                assert_eq!(directional.impulse.energy.dissipated_j, 0.0);
                assert_eq!(directional.impulse.time_s, 0.0);
            }
            dense.advance(&[0.0; 3], cx).unwrap();
            dense_ellipse.advance(&[0.0; 3], cx).unwrap();
            spectral.advance(&[0.0; 3], cx).unwrap();
            spectral_ellipse.advance(&[0.0; 3], cx).unwrap();
            assert_eq!(dense.velocity_mps(), dense_ellipse.velocity_mps());
            assert_eq!(dense.displacement_m(), dense_ellipse.displacement_m());
            assert_eq!(bits(&spectral), bits(&spectral_ellipse));
        }
    });
}

// Independent support-gradient construction from a chosen physical endpoint;
// it does not reuse the production congruence or scalar-disk solver.
fn sliding_impulse(pn: f64, v: [f64; 2], mu: [f64; 2], axis: [f64; 2]) -> [f64; 3] {
    let [c, s] = axis;
    let material_v = [c * v[0] + s * v[1], -s * v[0] + c * v[1]];
    let support = (mu[0] * material_v[0]).hypot(mu[1] * material_v[1]);
    let p = core::array::from_fn::<_, 2, _>(|a| -pn * mu[a] * mu[a] * material_v[a] / support);
    [pn, c * p[0] - s * p[1], s * p[0] + c * p[1]]
}

#[test]
fn manufactured_stick_slide_keeps_full_normal_shear_cross_contact_mobility_and_modes() {
    with_cx(false, |cx| {
        let phi = [
            1.0, 0.2, 0.0, 0.1, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 1.0, -0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 1.0,
        ];
        let inverse = [
            1.0, -0.2, 0.0, -0.1, -0.015, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.15, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
        ];
        let lambda = [0.0, 4.0, 9.0, 16.0, 25.0, 0.0];
        let (mut mass, mut stiffness, mut mobility) = ([0.0; 36], [0.0; 36], [0.0; 36]);
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    mass[6 * i + j] += inverse[6 * k + i] * inverse[6 * k + j];
                    stiffness[6 * i + j] += inverse[6 * k + i] * lambda[k] * inverse[6 * k + j];
                    mobility[6 * i + j] += phi[6 * i + k] * phi[6 * j + k];
                }
            }
        }
        assert!(mobility[1] != 0.0 && mobility[3] != 0.0 && mobility[22] != 0.0);
        let pencil = assembly(&stiffness, &mass, 6);
        let surface = [0.1, -0.2, 0.3, -0.15, 0.05, -0.1];
        let relative_after = [0.0, 0.4, 0.3, 0.0, 0.0, 0.0];
        let p0 = sliding_impulse(1.2, [0.4, 0.3], [0.3, 0.6], [0.6, 0.8]);
        let expected_p = [p0[0], p0[1], p0[2], 0.8, -0.07, 0.04];
        let delta = action(&mobility, &expected_p);
        let after: [f64; 6] = core::array::from_fn(|i| surface[i] + relative_after[i]);
        let before: [f64; 6] = core::array::from_fn(|i| after[i] - delta[i]);
        let u = [0.03, -0.04, 0.05, 0.06, -0.02, 0.01];
        let rows: [[f64; 6]; 6] =
            core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
        let contacts = [
            directional(
                contact(
                    40,
                    [&rows[0], &rows[1], &rows[2]],
                    [surface[0], surface[1], surface[2]],
                    0.0,
                ),
                [0.3, 0.6],
                [0.6, 0.8],
            ),
            directional(
                contact(
                    7,
                    [&rows[3], &rows[4], &rows[5]],
                    [surface[3], surface[4], surface[5]],
                    0.0,
                ),
                [0.45, 0.2],
                [-4.0, 3.0],
            ),
        ];
        let mut dense =
            ExactLinearSolid::prepare(&pencil, 0.03, &u, &before, dense_controls(), cx).unwrap();
        let make_spectral = || {
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &lambda,
                &phi,
                0.03,
                &u,
                &before,
                spectral_controls(),
                cx,
            )
            .unwrap()
        };
        let (mut spectral, mut permuted) = (make_spectral(), make_spectral());
        let a = dense
            .apply_elliptic_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        let b = spectral
            .apply_elliptic_impacts(&contacts, impact_controls(), cx)
            .unwrap();
        permuted
            .apply_elliptic_impacts(&[contacts[1], contacts[0]], impact_controls(), cx)
            .unwrap();
        assert_eq!(bits(&spectral), bits(&permuted));
        let expected_relative_work: f64 = (0..6)
            .map(|i| expected_p[i] * (relative_after[i] - 0.5 * delta[i]))
            .sum();
        let expected_surface_work: f64 = (0..6).map(|i| expected_p[i] * surface[i]).sum();
        for (velocity, report) in [(dense.velocity_mps(), &a), (spectral.velocity_mps(), &b)] {
            for i in 0..6 {
                near(velocity[i], after[i], 4e-9);
                near(
                    report.contact.law.impulses[i / 3][i % 3],
                    expected_p[i],
                    4e-9,
                );
            }
            near(report.contact.relative_work_j, expected_relative_work, 4e-9);
            near(
                report.contact.prescribed_surface_work_j,
                expected_surface_work,
                4e-9,
            );
            near(
                report.impulse.kinetic_change_j,
                expected_relative_work + expected_surface_work,
                4e-9,
            );
            assert!(report.contact.measured_residuals.sliding_impulse_residual < 1e-9);
            assert!(report.contact.measured_residuals.tangential_duality_gap < 1e-9);
            assert!(
                report
                    .contact
                    .measured_residuals
                    .coordinate_work_residual_j
                    .abs()
                    < 1e-9
            );
        }
        dense.advance(&[0.0; 6], cx).unwrap();
        spectral.advance(&[0.0; 6], cx).unwrap();
        let (q0, dq0) = (action(&inverse, &u), action(&inverse, &after));
        let (mut q, mut dq) = ([0.0; 6], [0.0; 6]);
        for i in 0..6 {
            let omega: f64 = lambda[i].sqrt();
            if omega == 0.0 {
                q[i] = q0[i] + 0.03 * dq0[i];
                dq[i] = dq0[i];
            } else {
                let (s, c) = (0.03 * omega).sin_cos();
                q[i] = q0[i] * c + dq0[i] * s / omega;
                dq[i] = -omega * q0[i] * s + dq0[i] * c;
            }
        }
        let (expected_u, expected_v) = (action(&phi, &q), action(&phi, &dq));
        for (u, v) in [
            (dense.displacement_m(), dense.velocity_mps()),
            (spectral.displacement_m(), spectral.velocity_mps()),
        ] {
            for i in 0..6 {
                near(u[i], expected_u[i], 4e-9);
                near(v[i], expected_v[i], 4e-9);
            }
        }
    });
}

#[test]
fn rotating_material_texture_changes_resistance_and_free_slip_is_frame_covariant() {
    with_cx(false, |cx| {
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let pencil = assembly(&[0.0; 9], &identity, 3);
        for (axis, traction) in [([1.0, 0.0], 0.4), ([0.0, 1.0], 0.1)] {
            let mut solid = ExactLinearSolid::prepare(
                &pencil,
                0.02,
                &[0.0; 3],
                &[-1.0, 2.0, 0.0],
                dense_controls(),
                cx,
            )
            .unwrap();
            let c = directional(
                contact(
                    1,
                    [&identity[0..3], &identity[3..6], &identity[6..9]],
                    [0.0; 3],
                    0.0,
                ),
                [0.4, 0.1],
                axis,
            );
            let report = solid
                .apply_elliptic_impacts(&[c], impact_controls(), cx)
                .unwrap();
            near(report.contact.law.impulses[0][1], -traction, 2e-9);
            near(solid.velocity_mps()[1], 2.0 - traction, 2e-9);
        }
        let surface_world = [0.2, -0.1, 0.4];
        let initial = [-0.8, -0.02, 3.4];
        for angle in [0.0_f64, 0.71, -1.2] {
            let (s, c) = angle.sin_cos();
            let rows = [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]];
            let surface =
                core::array::from_fn(|i| (0..3).map(|a| rows[i][a] * surface_world[a]).sum());
            let contact = directional(
                contact(2, [&rows[0], &rows[1], &rows[2]], surface, 0.0),
                [0.3, 0.0],
                [3.0 * c, -3.0 * s],
            );
            let mut dense =
                ExactLinearSolid::prepare(&pencil, 0.02, &[0.0; 3], &initial, dense_controls(), cx)
                    .unwrap();
            let mut spectral = CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &[0.0; 3],
                &identity,
                0.02,
                &[0.0; 3],
                &initial,
                spectral_controls(),
                cx,
            )
            .unwrap();
            let a = dense
                .apply_elliptic_impacts(&[contact], impact_controls(), cx)
                .unwrap();
            let b = spectral
                .apply_elliptic_impacts(&[contact], impact_controls(), cx)
                .unwrap();
            for (velocity, report) in [(dense.velocity_mps(), &a), (spectral.velocity_mps(), &b)] {
                near(velocity[0], surface_world[0], 3e-9);
                near(velocity[1], surface_world[1], 3e-9);
                near(velocity[2], initial[2], 3e-9);
                let impulse = report.contact.law.impulses[0];
                near(c * impulse[1] - s * impulse[2], -0.08, 3e-9);
                near(s * impulse[1] + c * impulse[2], 0.0, 3e-9);
                near(report.contact.relative_work_j, -0.5032, 3e-9);
                near(report.contact.prescribed_surface_work_j, 0.208, 3e-9);
                near(report.impulse.kinetic_change_j, -0.2952, 3e-9);
                assert!(report.contact.measured_residuals.ellipse_impulse_violation < 1e-9);
            }
        }
    });
}

#[test]
fn late_native_directional_refusal_preserves_state_and_retries_on_the_original_owner() {
    with_cx(false, |cx| {
        let identity = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let pencil = assembly(&[0.0; 9], &identity, 3);
        let make = || {
            let mut controls = spectral_controls();
            controls.maximum_relative_operator_error = 1e-3;
            controls.maximum_relative_state_error = 1e-3;
            controls.maximum_work_error_j = 1e-3;
            CompleteSpectralSolid::prepare_with_basis(
                &pencil,
                &[0.0; 3],
                &[1.0001, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
                0.05,
                &[0.0; 3],
                &[-0.2, 0.3, 0.4],
                controls,
                cx,
            )
            .unwrap()
        };
        let (mut solid, mut reference) = (make(), make());
        let c = directional(
            contact(
                1,
                [&identity[0..3], &identity[3..6], &identity[6..9]],
                [0.0; 3],
                0.0,
            ),
            [0.2, 0.3],
            [1.0, 0.0],
        );
        let before = bits(&solid);
        let mut controls = impact_controls();
        controls.maximum_impulse_error_kg_m_s = 1e-5;
        controls.maximum_work_error_j = 1e-3;
        assert!(matches!(
            solid.apply_elliptic_impacts(&[c], controls, cx),
            Err(SolidFrictionalImpactError::Speed { .. })
        ));
        assert_eq!(bits(&solid), before);
        controls.maximum_contact_speed_error_mps = 1e-3;
        assert!(matches!(
            solid.apply_elliptic_impacts(&[c], controls, cx),
            Err(SolidFrictionalImpactError::EllipticContact(
                EllipticFrictionalContactError::ResidualNotSatisfied(_)
            ))
        ));
        let mut budget = controls;
        budget.maximum_workspace_bytes = 0;
        assert!(matches!(
            solid.apply_elliptic_impacts(&[c], budget, cx),
            Err(SolidFrictionalImpactError::Budget { .. })
        ));
        let mut invalid = c;
        invalid.friction_coefficient.principal_axis = [0.0; 2];
        assert!(
            solid
                .apply_elliptic_impacts(&[invalid], controls, cx)
                .is_err()
        );
        invalid = c;
        invalid.friction_coefficient.coefficients_mu[0] = -0.1;
        assert!(
            solid
                .apply_elliptic_impacts(&[invalid], controls, cx)
                .is_err()
        );
        assert!(solid.apply_elliptic_impacts(&[c, c], controls, cx).is_err());
        with_cx(true, |cancelled| {
            assert!(matches!(
                solid.apply_elliptic_impacts(&[c], controls, cancelled),
                Err(SolidFrictionalImpactError::Cancelled)
            ));
        });
        assert_eq!(bits(&solid), before);
        controls.contact.velocity_tolerance = 1e-3;
        controls.contact.complementarity_tolerance = 1e-4;
        controls.contact.energy_tolerance = 1e-4;
        solid.apply_elliptic_impacts(&[c], controls, cx).unwrap();
        reference
            .apply_elliptic_impacts(&[c], controls, cx)
            .unwrap();
        assert_eq!(bits(&solid), bits(&reference));
        solid.advance(&[0.0; 3], cx).unwrap();
        reference.advance(&[0.0; 3], cx).unwrap();
        assert_eq!(bits(&solid), bits(&reference));
    });
}
