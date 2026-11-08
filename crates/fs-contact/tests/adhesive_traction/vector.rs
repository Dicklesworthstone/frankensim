//! Full-stick shear distribution and work on the original peeled JKR disk.
use super::*;
use core::f64::consts::TAU;
use fs_contact::normal_patch::adhesive::torsion::{JkrTorsionControls, SourceJkrTorsionContact};
use fs_contact::normal_patch::traction::sectors::{
    PolarTransferBudget, RadialMomentBudget, SectorTractionError,
    vector::{SectorAdhesiveVectorTraction, VectorPolarRingScratch},
};
use fs_exec::CancelGate;

fn exact() -> RadialMomentBudget {
    RadialMomentBudget {
        absolute_tolerance_n_m: 0.0,
        relative_tolerance: 0.0,
        maximum_evaluations: 0,
    }
}

fn profile() -> (SourceJkrMixedMode, AdhesiveTractionProfile) {
    let source = mixed_jkr();
    let limit = source.critical_tangential_force_n(0.0).unwrap();
    let profile =
        AdhesiveTractionProfile::jkr_mixed_mode_force(&source, 0.0, [0.3 * limit, -0.4 * limit])
            .unwrap();
    (source, profile)
}

#[test]
fn full_stick_vector_field_and_sector_couples_follow_original_mindlin_solution() {
    let (source, profile) = profile();
    let force = profile.tangential_force_n();
    let limit = source.critical_tangential_force_n(0.0).unwrap();
    let response = source
        .evaluate_force_controlled(0.0, [0.3 * limit, -0.4 * limit])
        .unwrap();
    assert_eq!(response.tangential_force_n, force);
    assert_eq!(response.model_identity(), profile.model_identity());
    let a = profile.contact_radius_m();
    for fraction in [0.0_f64, 0.25, 0.8, 0.999] {
        let actual = profile.traction_pa(fraction * a).unwrap();
        let denominator = 2.0 * PI * a * a * (1.0 - fraction * fraction).sqrt();
        for i in 0..2 {
            let expected = force[i] / denominator;
            close(actual[i], expected, expected.abs(), 2e-13);
        }
        assert_eq!(actual[2], profile.pressure_pa(fraction * a).unwrap());
    }
    assert_eq!(
        profile.traction_pa(a),
        Err(AdhesiveTractionError::SingularEdge)
    );
    assert_eq!(profile.traction_pa(a.next_up()).unwrap(), [0.0; 3]);

    let gate = CancelGate::new();
    let cell = profile
        .vector_sector(0.0, a, 0.0, PI / 2.0, exact(), &gate)
        .unwrap();
    assert_eq!(
        cell.normal,
        profile
            .sector(0.0, a, 0.0, PI / 2.0, exact(), &gate)
            .unwrap()
    );
    for i in 0..2 {
        // Independent quarter-disk integrals of T_i/(2*pi*a*sqrt(a²-r²)).
        close(
            cell.tangential_force_n[i],
            force[i] / 4.0,
            force[i].abs(),
            2e-14,
        );
        for j in 0..2 {
            close(
                cell.tangential_first_moment_n_m[i][j],
                force[i] * a / 8.0,
                force[i].abs() * a,
                2e-14,
            );
        }
    }
    let (resultant, torque) = cell.wrench();
    close(
        torque[2],
        (force[1] - force[0]) * a / 8.0,
        force[0].hypot(force[1]) * a,
        2e-14,
    );
    let translation = [0.21, -0.17, 0.03];
    let omega = [0.7, -0.2, 1.5];
    let gradient = [[0.0, -omega[2]], [omega[2], 0.0], [-omega[1], omega[0]]];
    let (power, estimate) = cell.affine_power(translation, gradient).unwrap();
    let expected = (0..3)
        .map(|i| resultant[i] * translation[i] + torque[i] * omega[i])
        .sum::<f64>();
    close(power, expected, expected.abs(), 2e-14);
    assert_eq!(estimate, 0.0);
    close(
        0.5 * (force[0] * response.tangential_displacement_m[0]
            + force[1] * response.tangential_displacement_m[1]),
        response.tangential_elastic_energy_j,
        response.tangential_elastic_energy_j,
        3e-14,
    );
}

#[test]
fn vector_grid_recovers_all_efforts_affine_work_and_thin_edge_shear() {
    let (source, profile) = profile();
    let a = profile.contact_radius_m();
    let q = profile.tangential_force_n();
    let scale = q[0].hypot(q[1]);
    let gate = CancelGate::new();
    let radii = [0.0, 0.2 * a, 0.63 * a, a.next_down(), a, 1.1 * a];
    let angles = [0.0, 0.2, PI / 2.0, PI, 1.7 * PI, TAU];
    let mut cells = vec![SectorAdhesiveVectorTraction::default(); 25];
    let mut scratch = vec![VectorPolarRingScratch::default(); 5];
    let report = profile
        .transfer_vector_polar_grid(
            &radii,
            &angles,
            &mut cells,
            &mut scratch,
            PolarTransferBudget {
                radial_moment: exact(),
                force_absolute_tolerance_n: 3e-12 * scale,
                moment_absolute_tolerance_n_m: 3e-12 * scale * a,
            },
            &gate,
        )
        .unwrap();
    assert_eq!(report.normal.model_identity, source.identity());
    assert_eq!(
        report.normal.normal_state_identity,
        profile.normal_state_identity()
    );
    assert_eq!(report.normal.authority, InputAuthority::SyntheticFixture);
    assert_eq!(report.normal.radial_evaluations, 0);
    assert_eq!(report.normal.maximum_sector_moment_error_estimate_n_m, 0.0);

    let velocity = [0.21, -0.17, 0.03];
    let gradient = [[2.0, -3.0], [4.0, 1.0], [-1.0, 5.0]];
    let mut force = [0.0; 3];
    let mut moment = [[0.0; 2]; 3];
    let mut power = 0.0;
    for cell in &cells {
        let f = cell.force_n();
        let m = cell.first_moments_n_m();
        for i in 0..3 {
            force[i] += f[i];
            for j in 0..2 {
                moment[i][j] += m[i][j];
            }
        }
        let (work, estimate) = cell.affine_power(velocity, gradient).unwrap();
        assert_eq!(estimate, 0.0);
        power += work;
    }
    for i in 0..2 {
        close(force[i], q[i], scale, 3e-13);
    }
    close(force[2], 0.0, scale, 3e-13);
    for value in moment.into_iter().flatten() {
        close(value, 0.0, scale * a, 3e-13);
    }
    close(power, q[0] * velocity[0] + q[1] * velocity[1], scale, 3e-13);

    let edge = profile
        .vector_sector(a.next_down(), a, 0.0, TAU, exact(), &gate)
        .unwrap();
    let edge_fraction = (((a - a.next_down()) / a) * (1.0 + a.next_down() / a)).sqrt();
    for i in 0..2 {
        assert_ne!(edge.tangential_force_n[i], 0.0);
        close(
            edge.tangential_force_n[i],
            q[i] * edge_fraction,
            (q[i] * edge_fraction).abs(),
            3e-14,
        );
    }
    assert_eq!(edge.tangential_first_moment_n_m, [[0.0; 2]; 2]);

    let whole = profile
        .vector_sector(0.0, 1.1 * a, 0.0, TAU, exact(), &gate)
        .unwrap();
    let (coarse, _) = whole.affine_power(velocity, gradient).unwrap();
    close(coarse, power, scale, 3e-13);
}

#[test]
fn vector_transfer_is_covariant_and_zero_shear_does_not_create_a_gap_force() {
    let (source, first) = profile();
    let limit = source.critical_tangential_force_n(0.0).unwrap();
    let rotated =
        AdhesiveTractionProfile::jkr_mixed_mode_force(&source, 0.0, [0.4 * limit, 0.3 * limit])
            .unwrap();
    let a = first.contact_radius_m();
    assert_eq!(rotated.contact_radius_m(), a);
    let gate = CancelGate::new();
    let one = first
        .vector_sector(0.2 * a, a, 0.3, 2.1, exact(), &gate)
        .unwrap();
    let two = rotated
        .vector_sector(0.2 * a, a, 0.3, 2.1, exact(), &gate)
        .unwrap();
    assert_eq!(one.normal, two.normal);
    assert_eq!(
        two.tangential_force_n,
        [-one.tangential_force_n[1], one.tangential_force_n[0]]
    );
    for j in 0..2 {
        assert_eq!(
            two.tangential_first_moment_n_m[0][j],
            -one.tangential_first_moment_n_m[1][j]
        );
        assert_eq!(
            two.tangential_first_moment_n_m[1][j],
            one.tangential_first_moment_n_m[0][j]
        );
    }

    let pure = AdhesiveTractionProfile::jkr_attached(source.normal_source(), 0.0).unwrap();
    let a = pure.contact_radius_m();
    let cell = pure
        .vector_sector(0.0, a, 0.0, PI / 2.0, exact(), &gate)
        .unwrap();
    assert_eq!(pure.tangential_force_n(), [0.0; 2]);
    assert_eq!(cell.tangential_force_n, [0.0; 2]);
    assert_eq!(cell.tangential_first_moment_n_m, [[0.0; 2]; 2]);
    assert_eq!(
        cell.normal,
        pure.sector(0.0, a, 0.0, PI / 2.0, exact(), &gate).unwrap()
    );

    let source = maugis(0.2);
    let limits = source.noncontact_limits().unwrap().unwrap();
    let delta = limits
        .minimum_indentation_m
        .midpoint(limits.first_contact_indentation_m);
    let state = source.evaluate_noncontact(delta).unwrap();
    let gap = AdhesiveTractionProfile::maugis(&source, &state).unwrap();
    let cell = gap
        .vector_sector(
            0.0,
            gap.interaction_radius_m(),
            0.0,
            PI / 2.0,
            exact(),
            &gate,
        )
        .unwrap();
    assert!(cell.normal.normal_force_n < 0.0);
    assert_eq!(cell.normal.solid_contact_area_m2, 0.0);
    assert_eq!(cell.tangential_force_n, [0.0; 2]);
    assert_eq!(cell.tangential_first_moment_n_m, [[0.0; 2]; 2]);
    assert_eq!(gap.traction_pa(0.0).unwrap()[..2], [0.0; 2]);
}

#[test]
fn vector_partition_cancellation_and_original_source_refusals_are_atomic() {
    let (source, profile) = profile();
    let a = profile.contact_radius_m();
    let gate = CancelGate::new();
    let mut output = [SectorAdhesiveVectorTraction::default(); 2];
    output[0].tangential_force_n = [7.0, -9.0];
    output[1].tangential_first_moment_n_m = [[1.0, 2.0], [3.0, 4.0]];
    let before = output;
    let mut scratch = [VectorPolarRingScratch::default(); 1];
    let budget = PolarTransferBudget {
        radial_moment: exact(),
        force_absolute_tolerance_n: 1e-12,
        moment_absolute_tolerance_n_m: 1e-12,
    };
    for (radii, angles) in [
        ([0.0, 0.9 * a], [0.0, PI, TAU]),
        ([0.0, a], [0.0, PI, 1.9 * PI]),
    ] {
        assert!(
            profile
                .transfer_vector_polar_grid(
                    &radii,
                    &angles,
                    &mut output,
                    &mut scratch,
                    budget,
                    &gate
                )
                .is_err()
        );
        assert_eq!(output, before);
    }
    gate.request();
    assert_eq!(
        profile.transfer_vector_polar_grid(
            &[0.0, a],
            &[0.0, PI, TAU],
            &mut output,
            &mut scratch,
            budget,
            &gate
        ),
        Err(SectorTractionError::Cancelled)
    );
    assert_eq!(output, before);
    let limit = source.critical_tangential_force_n(0.0).unwrap();
    assert!(
        AdhesiveTractionProfile::jkr_mixed_mode_force(&source, 0.0, [limit.next_up(), 0.0],)
            .is_err()
    );
    assert!(profile.traction_pa(f64::NAN).is_err());
}

fn elastic_torsion() -> (SourceJkrTorsionContact, f64, AdhesiveTractionProfile) {
    let normal = jkr();
    let indentation = normal.evaluate_force_controlled(0.0).unwrap().indentation_m;
    let source = SourceJkrTorsionContact::new(
        normal,
        indentation,
        JkrTorsionControls {
            critical_twist_angle_rad: 1e-5,
            maximum_critical_twist_angle_rad: 0.01,
            maximum_mode_iii_energy_release_to_adhesion: 0.01,
        },
    )
    .unwrap();
    let angle = 4e-6;
    let profile = AdhesiveTractionProfile::jkr_torsion_elastic(&source, angle).unwrap();
    (source, angle, profile)
}

#[test]
fn fresh_elastic_torsion_field_retains_original_area_sign_and_refuses_plastic_profiles() {
    let (source, angle, profile) = elastic_torsion();
    let a = source.normal_response().contact_radius_m;
    let g = source.law().reduced_shear_modulus_pa();
    let moment = -(16.0 / 3.0) * g * a.powi(3) * angle;
    close(profile.torsional_moment_n_m(), moment, moment.abs(), 2e-14);
    assert_eq!(profile.contact_radius_m(), a);
    assert_eq!(profile.model_identity(), source.identity());
    assert_eq!(
        profile.normal_state_identity(),
        source.normal_source().normal_state().identity()
    );
    assert_eq!(profile.tangential_force_n(), [0.0; 2]);
    assert_eq!(profile.authority(), InputAuthority::SyntheticFixture);
    let normal = AdhesiveTractionProfile::jkr_attached(
        source.normal_source(),
        source.normal_response().indentation_m,
    )
    .unwrap();
    for s in [0.0_f64, 0.25, 0.8, 0.999] {
        let azimuthal = -4.0 * g * angle * s / (PI * (1.0 - s * s).sqrt());
        for theta in [0.0_f64, 0.7, PI / 2.0, 4.3] {
            let t = profile.traction_polar_pa(s * a, theta).unwrap();
            close(t[0], -azimuthal * theta.sin(), azimuthal.abs(), 3e-13);
            close(t[1], azimuthal * theta.cos(), azimuthal.abs(), 3e-13);
            assert_eq!(t[2], normal.pressure_pa(s * a).unwrap());
        }
    }
    assert!(profile.traction_pa(0.3 * a).is_err());
    assert_eq!(
        profile.traction_polar_pa(a, 0.0),
        Err(AdhesiveTractionError::SingularEdge)
    );
    assert_eq!(
        profile.traction_polar_pa(a.next_up(), 0.0).unwrap(),
        [0.0; 3]
    );
    assert!(profile.traction_polar_pa(0.0, f64::NAN).is_err());
    let critical = source.law().parameters().critical_twist_angle_rad;
    for phi in [critical, -critical, 0.0] {
        let p = AdhesiveTractionProfile::jkr_torsion_elastic(&source, phi).unwrap();
        close(
            p.torsional_moment_n_m(),
            -(16.0 / 3.0) * g * a.powi(3) * phi,
            moment.abs(),
            3e-14,
        );
    }
    for phi in [critical.next_up(), -critical.next_up(), f64::NAN] {
        assert!(AdhesiveTractionProfile::jkr_torsion_elastic(&source, phi).is_err());
    }
}

#[test]
fn torsion_sectors_integrate_edge_singularities_and_tiny_angular_moments() {
    let (source, angle, profile) = elastic_torsion();
    let a = profile.contact_radius_m();
    let m = -(16.0 / 3.0) * source.law().reduced_shear_modulus_pa() * a.powi(3) * angle;
    let gate = CancelGate::new();
    let quarter = profile
        .vector_sector(0.0, a, 0.0, PI / 2.0, exact(), &gate)
        .unwrap();
    // Independent integrals of r²/sqrt(a²-r²) and r³/sqrt(a²-r²).
    close(
        quarter.tangential_force_n[0],
        -3.0 * m / (16.0 * a),
        m.abs() / a,
        3e-14,
    );
    close(
        quarter.tangential_force_n[1],
        3.0 * m / (16.0 * a),
        m.abs() / a,
        3e-14,
    );
    let expected = [[-m / (4.0 * PI), -m / 8.0], [m / 8.0, m / (4.0 * PI)]];
    for (actual, expected) in quarter
        .tangential_first_moment_n_m
        .into_iter()
        .flatten()
        .zip(expected.into_iter().flatten())
    {
        close(actual, expected, m.abs(), 3e-14);
    }
    let edge = profile
        .vector_sector(a.next_down(), a, 0.0, PI / 2.0, exact(), &gate)
        .unwrap();
    let u = (((a - a.next_down()) / a) * (1.0 + a.next_down() / a)).sqrt();
    let edge_force = 3.0 * m / (8.0 * PI * a) * (u.asin() + u * a.next_down() / a);
    close(
        edge.tangential_force_n[0],
        -edge_force,
        edge_force.abs(),
        3e-14,
    );
    close(
        edge.tangential_force_n[1],
        edge_force,
        edge_force.abs(),
        3e-14,
    );
    let edge_torque = 0.375 * m * (u - u.powi(3) / 3.0);
    close(edge.wrench().1[2], edge_torque, edge_torque.abs(), 3e-14);
    assert_ne!(edge_torque, 0.0);

    let width = 1e-8;
    let thin = profile
        .vector_sector(0.0, a, 0.0, width, exact(), &gate)
        .unwrap();
    // ∫sin²(theta)dtheta ~ width³/3: subtracting width-sin(width)
    // directly would erase this nonzero P1 work coefficient.
    let small = -m * width.powi(3) / (6.0 * PI);
    close(
        thin.tangential_first_moment_n_m[0][1],
        small,
        small.abs(),
        3e-14,
    );
    let central = profile
        .vector_sector(0.0, 1e-6 * a, 0.0, TAU, exact(), &gate)
        .unwrap();
    let central_torque = m * 0.375 * 1e-24 * (1.0 + 1e-12 / 3.0);
    close(
        central.wrench().1[2],
        central_torque,
        central_torque.abs(),
        3e-14,
    );
}

#[test]
fn torsion_vector_grid_recovers_source_couple_and_affine_rotation_work() {
    let (_, _, profile) = elastic_torsion();
    let a = profile.contact_radius_m();
    let m = profile.torsional_moment_n_m();
    let whole = profile.annulus(0.0, a).unwrap();
    let scale = whole.hertz_force_n - whole.adhesive_force_n;
    let gate = CancelGate::new();
    let radii = [0.0, 0.2 * a, 0.63 * a, a.next_down(), a, 1.1 * a];
    let angles = [0.0, 0.2, PI / 2.0, PI, 1.7 * PI, TAU];
    let mut cells = vec![SectorAdhesiveVectorTraction::default(); 25];
    let mut scratch = vec![VectorPolarRingScratch::default(); 5];
    let report = profile
        .transfer_vector_polar_grid(
            &radii,
            &angles,
            &mut cells,
            &mut scratch,
            PolarTransferBudget {
                radial_moment: exact(),
                force_absolute_tolerance_n: 3e-12 * scale,
                moment_absolute_tolerance_n_m: 3e-12 * scale * a,
            },
            &gate,
        )
        .unwrap();
    assert_eq!(report.normal.model_identity, profile.model_identity());
    assert_eq!(report.normal.radial_evaluations, 0);
    let omega = 1.7;
    let (mut force, mut torque, mut power) = ([0.0; 2], 0.0, 0.0);
    for cell in cells {
        for i in 0..2 {
            force[i] += cell.tangential_force_n[i];
        }
        torque += cell.wrench().1[2];
        power += cell
            .affine_power([0.0; 3], [[0.0, -omega], [omega, 0.0], [0.0; 2]])
            .unwrap()
            .0;
    }
    for f in force {
        close(f, 0.0, m.abs() / a, 3e-13);
    }
    close(torque, m, m.abs(), 3e-13);
    close(power, m * omega, (m * omega).abs(), 3e-13);
    let full = profile
        .vector_sector(0.0, a, 0.0, TAU, exact(), &gate)
        .unwrap();
    assert_eq!(full.tangential_force_n, [0.0; 2]);
    close(full.wrench().1[2], m, m.abs(), 3e-14);
    let (coarse_power, error) = full
        .affine_power([0.0; 3], [[0.0, -omega], [omega, 0.0], [0.0; 2]])
        .unwrap();
    close(coarse_power, power, power.abs(), 3e-13);
    assert_eq!(error, 0.0);
}
