//! Full-stick shear distribution and work on the original peeled JKR disk.
use super::*;
use core::f64::consts::TAU;
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
