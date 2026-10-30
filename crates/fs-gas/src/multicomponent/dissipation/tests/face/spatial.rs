use super::*;
use crate::multicomponent::{SpatialCell, SpatialPrimitive};

mod reservoir;

fn prepared3(formation: f64) -> SpatialTransport<ConstantTransport> {
    let planar = prepared(formation);
    SpatialTransport::new(planar.gas().clone(), planar.law().clone(), || false).unwrap()
}

fn input3(densities: &[f64], velocity: [f64; 3], pressure: f64) -> SpatialPrimitive<'_> {
    SpatialPrimitive {
        partial_densities: densities,
        normal_velocity: velocity[0],
        transverse_velocity: [velocity[1], velocity[2]],
        pressure,
    }
}

fn cell3(q: &[f64; 6]) -> SpatialCell<'_> {
    SpatialCell {
        partial_densities: &q[..2],
        normal_momentum: q[2],
        transverse_momentum: [q[3], q[4]],
        energy: q[5],
    }
}

#[test]
fn spatial_species_heat_and_full_tensor_work_match_independent_conservative_flux() {
    let mut op = prepared3(1e6);
    let left = input3(&[0.8, 0.2], [0.3, -0.4, 1.2], 1.5);
    let right = input3(&[0.2, 0.8], [1.3, 1.6, -0.8], 2.5);
    let gradients = [[6.0, -3.0, 7.0], [-5.0, 11.0, 4.0]];
    let mut ql = [0.0; 6];
    let mut qr = [0.0; 6];
    op.gas().spatial_conserved_into(left, &mut ql).unwrap();
    op.gas().spatial_conserved_into(right, &mut qr).unwrap();
    let before = (ql, qr, op.workspace_bytes());
    let mut flux = [0.0; 6];
    let report = op
        .face_from_primitives(
            left,
            right,
            0.5,
            gradients,
            MaxwellStefanSolveV1::default(),
            &mut flux,
            || false,
        )
        .unwrap();
    // Equal molecular masses reduce Maxwell--Stefan to Fick's law here.
    // The species still have unequal Cp and large independent formation gauges.
    let j_a = -0.03 * (0.2 - 0.8) / 0.5;
    let tau_nn = (4.0 / 3.0 * 0.07 + 0.02) * 2.0 + (0.02 - 2.0 / 3.0 * 0.07) * 1.0;
    let tau_t1n = 0.07 * (4.0 + 6.0);
    let tau_t2n = 0.07 * (-4.0 - 5.0);
    let h_a = 2.5 * (2.0 - 0.001) + 2.0;
    let h_b = 4.0 * (2.0 - 0.001) + 2.0;
    let expected = [
        j_a,
        -j_a,
        -tau_nn,
        -tau_t1n,
        -tau_t2n,
        j_a * (h_a - h_b) - 0.09 * 2.0 - 0.8 * tau_nn - 0.6 * tau_t1n - 0.2 * tau_t2n,
    ];
    for (a, b) in flux.into_iter().zip(expected) {
        close(a, b);
    }
    assert!(report.maximum_gradient_residual_per_m < 1e-13);
    assert!(report.mass_frame_residual_kg_m2_s.abs() < 1e-14);
    assert!(report.diffusivity_m2_s >= 0.07 * (4.0 / 3.0 + 4.0));
    let mut conservative = [0.0; 6];
    op.face(
        cell3(&ql),
        cell3(&qr),
        0.5,
        gradients,
        MaxwellStefanSolveV1::default(),
        &mut conservative,
        || false,
    )
    .unwrap();
    for (a, b) in conservative.into_iter().zip(flux) {
        close(a, b);
    }
    assert_eq!((ql, qr, op.workspace_bytes()), before);
}

#[test]
fn spatial_flux_rotates_both_tangents_and_preserves_reference_and_galilean_covariance() {
    let mut op = prepared3(0.0);
    let mut shifted = prepared3(1e6);
    let left_velocity = [0.3, -0.4, 1.2];
    let right_velocity = [1.3, 1.6, -0.8];
    let gradients = [[6.0, -3.0, 7.0], [-5.0, 11.0, 4.0]];
    let mut flux = [0.0; 6];
    op.face_from_primitives(
        input3(&[0.8, 0.2], left_velocity, 1.5),
        input3(&[0.2, 0.8], right_velocity, 2.5),
        0.5,
        gradients,
        MaxwellStefanSolveV1::default(),
        &mut flux,
        || false,
    )
    .unwrap();
    let mut gauge = [0.0; 6];
    shifted
        .face_from_primitives(
            input3(&[0.8, 0.2], left_velocity, 1.5),
            input3(&[0.2, 0.8], right_velocity, 2.5),
            0.5,
            gradients,
            MaxwellStefanSolveV1::default(),
            &mut gauge,
            || false,
        )
        .unwrap();
    for (a, b) in gauge.into_iter().zip(flux) {
        close(a, b);
    }
    let boost = [2.1, -1.7, 3.2];
    let mut boosted = [0.0; 6];
    op.face_from_primitives(
        input3(
            &[0.8, 0.2],
            core::array::from_fn(|i| left_velocity[i] + boost[i]),
            1.5,
        ),
        input3(
            &[0.2, 0.8],
            core::array::from_fn(|i| right_velocity[i] + boost[i]),
            2.5,
        ),
        0.5,
        gradients,
        MaxwellStefanSolveV1::default(),
        &mut boosted,
        || false,
    )
    .unwrap();
    for i in 0..5 {
        close(boosted[i], flux[i]);
    }
    close(
        boosted[5],
        flux[5] + (0..3).map(|i| boost[i] * flux[i + 2]).sum::<f64>(),
    );

    // Rotate the complete derivative tensor, including cross derivatives that
    // do not enter this face before rotation. The normal direction stays fixed.
    let basis = [[1.0, 0.0, 0.0], [0.0, 0.6, 0.8], [0.0, -0.8, 0.6]];
    let rotate = |u: [f64; 3]| -> [f64; 3] {
        core::array::from_fn(|a| (0..3).map(|i| basis[a][i] * u[i]).sum())
    };
    let tensor = [[2.0, 4.0, -4.0], gradients[0], gradients[1]];
    let mut rotated_gradient = [[0.0; 3]; 2];
    for a in 0..2 {
        for b in 0..3 {
            for i in 0..3 {
                for j in 0..3 {
                    rotated_gradient[a][b] += basis[a + 1][i] * tensor[i][j] * basis[b][j];
                }
            }
        }
    }
    let mut rotated = [0.0; 6];
    op.face_from_primitives(
        input3(&[0.8, 0.2], rotate(left_velocity), 1.5),
        input3(&[0.2, 0.8], rotate(right_velocity), 2.5),
        0.5,
        rotated_gradient,
        MaxwellStefanSolveV1::default(),
        &mut rotated,
        || false,
    )
    .unwrap();
    for i in [0, 1, 2, 5] {
        close(rotated[i], flux[i]);
    }
    let rotated_momentum = rotate([flux[2], flux[3], flux[4]]);
    for i in 0..3 {
        close(rotated[i + 2], rotated_momentum[i]);
    }
}

#[test]
fn spatial_zero_second_transverse_limit_reuses_planar_flux_bit_for_bit() {
    let mut spatial = prepared3(0.0);
    let mut planar = prepared(0.0);
    let mut expected = [0.0; 5];
    planar
        .face_from_primitives(
            input(&[0.8, 0.2], [0.3, -0.4], 1.5),
            input(&[0.2, 0.8], [1.3, 1.6], 2.5),
            0.5,
            [6.0, -3.0],
            MaxwellStefanSolveV1::default(),
            &mut expected,
            || false,
        )
        .unwrap();
    let mut actual = [0.0; 6];
    spatial
        .face_from_primitives(
            input3(&[0.8, 0.2], [0.3, -0.4, 0.0], 1.5),
            input3(&[0.2, 0.8], [1.3, 1.6, 0.0], 2.5),
            0.5,
            [[6.0, -3.0, 0.0], [0.0; 3]],
            MaxwellStefanSolveV1::default(),
            &mut actual,
            || false,
        )
        .unwrap();
    assert_eq!(
        [actual[0], actual[1], actual[2], actual[3], actual[5]].map(f64::to_bits),
        expected.map(f64::to_bits)
    );
    assert_eq!(actual[4], 0.0);
}

#[test]
fn spatial_slip_wall_retains_both_transverse_dilatations_without_mass_shear_or_work() {
    let mut op = prepared3(0.0);
    let trace = input3(&[1.0, 0.0], [0.3, 1.6, -0.8], 2.0);
    let mut q = [0.0; 6];
    op.gas().spatial_conserved_into(trace, &mut q).unwrap();
    let mut out = [f64::NAN; 6];
    let report = op
        .wall_from_primitive(trace, 2.0, [-3.0, 4.0], &mut out, || false)
        .unwrap();
    close(
        out[2],
        -((4.0 / 3.0 * 0.07 + 0.02) * 2.0 + (0.02 - 2.0 / 3.0 * 0.07) * 1.0),
    );
    for i in [0, 1, 3, 4, 5] {
        assert_eq!(out[i], 0.0);
    }
    let mut conservative = [0.0; 6];
    let conservative_report = op
        .wall(cell3(&q), 2.0, [-3.0, 4.0], &mut conservative, || false)
        .unwrap();
    for (a, b) in conservative.into_iter().zip(out) {
        close(a, b);
    }
    close(
        report.diffusivity_m2_s,
        conservative_report.diffusivity_m2_s,
    );
    assert_eq!(report.mass_frame_residual_kg_m2_s, 0.0);
}

#[test]
fn spatial_original_diffusion_preserves_resolved_trace_and_exact_absent_species() {
    let g = MixtureGas::new(
        bank(vec![
            source("A", R, 3.5, 0.0, 0.0),
            source("B", R, 3.5, 0.0, 0.0),
            source("Trace", R, 3.5, 0.0, 0.0),
            source("Absent", R, 3.5, 0.0, 0.0),
        ]),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let flow = field(g.clone(), 2, Boundaries::Periodic, |_| {
        (vec![0.5, 0.5, 0.0, 0.0], 0.0, 2.0)
    });
    let mut op = SpatialTransport::new(g, law(&flow, 0.0, 0.0, 0.0, 0.03), || false).unwrap();
    let mut flux = [0.0; 8];
    op.face_from_primitives(
        input3(&[0.8, 0.2, 1e-180, 0.0], [0.0; 3], 2.0),
        input3(&[0.2, 0.8, 2e-180, 0.0], [0.0; 3], 2.0),
        0.5,
        [[0.0; 3]; 2],
        MaxwellStefanSolveV1::default(),
        &mut flux,
        || false,
    )
    .unwrap();
    let expected_trace = -0.03 * 1e-180 / 0.5;
    assert!(flux[2] < 0.0);
    assert!((flux[2] / expected_trace - 1.0).abs() < 1e-12);
    assert_eq!(flux[3], 0.0);
    close(flux[0], 0.03 * 0.6 / 0.5);
    close(flux[1], -flux[0]);
}

#[test]
fn spatial_faces_refuse_wrong_sources_unresolved_energy_and_invalid_axes_then_replay() {
    let mut op = prepared3(0.0);
    let other = prepared3(1e6);
    assert!(SpatialTransport::new(op.gas().clone(), other.law().clone(), || false).is_err());
    let left = input3(&[0.8, 0.2], [0.3, -0.4, 1.2], 1.5);
    let right = input3(&[0.2, 0.8], [1.3, 1.6, -0.8], 2.5);
    let gradients = [[6.0, -3.0, 7.0], [-5.0, 11.0, 4.0]];
    let mut out = [0.0; 6];
    for (distance, gradient) in [
        (0.0, gradients),
        (0.5, [[6.0, -3.0, 7.0], [-5.0, f64::NAN, 4.0]]),
    ] {
        assert!(
            op.face_from_primitives(
                left,
                right,
                distance,
                gradient,
                MaxwellStefanSolveV1::default(),
                &mut out,
                || false,
            )
            .is_err()
        );
    }
    assert!(
        op.face_from_primitives(
            left,
            right,
            0.5,
            gradients,
            MaxwellStefanSolveV1::default(),
            &mut out[..5],
            || false,
        )
        .is_err()
    );
    assert!(
        op.face_from_primitives(
            SpatialPrimitive {
                transverse_velocity: [1.0, 1e12],
                ..left
            },
            right,
            0.5,
            gradients,
            MaxwellStefanSolveV1::default(),
            &mut out,
            || false,
        )
        .is_err()
    );
    assert!(
        op.wall_from_primitive(left, 2.0, [-3.0, f64::NAN], &mut out, || false)
            .is_err()
    );
    let mut polls = 0;
    let report = op
        .face_from_primitives(
            left,
            right,
            0.5,
            gradients,
            MaxwellStefanSolveV1::default(),
            &mut out,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let expected = out;
    let capacity = op.workspace_bytes();
    let mut observed = 0;
    assert_eq!(
        op.face_from_primitives(
            left,
            right,
            0.5,
            gradients,
            MaxwellStefanSolveV1::default(),
            &mut out,
            || {
                observed += 1;
                observed == polls
            },
        ),
        Err(DissipationError::Cancelled)
    );
    let replay = op
        .face_from_primitives(
            left,
            right,
            0.5,
            gradients,
            MaxwellStefanSolveV1::default(),
            &mut out,
            || false,
        )
        .unwrap();
    assert_eq!(
        (out, replay, op.workspace_bytes()),
        (expected, report, capacity)
    );
}


#[test]
fn spatial_no_slip_wall_carries_two_shears_heat_and_wall_work() {
    let mut op = prepared3(0.0);
    let trace = input3(&[1.0, 0.0], [0.3, 1.6, -0.8], 2.0);
    let mut q = [0.0; 6];
    let primitive = op.gas().spatial_conserved_into(trace, &mut q).unwrap();
    let mut out = [f64::NAN; 6];
    op.no_slip_wall_from_primitive(
        trace,
        0.25,
        [0.4, -0.2],
        WallThermalCondition::Isothermal(Temperature::new(primitive.temperature - 0.5)),
        2.0,
        [-3.0, 4.0],
        &mut out,
        || false,
    )
    .unwrap();
    let cross = 0.02 - 2.0 / 3.0 * 0.07;
    let normal = (4.0 / 3.0 * 0.07 + 0.02) * 2.0 + cross;
    let shear1 = 0.07 * (0.4 - 1.6) / 0.25;
    let shear2 = 0.07 * (-0.2 + 0.8) / 0.25;
    close(out[2], -normal);
    close(out[3], -shear1);
    close(out[4], -shear2);
    close(
        out[5],
        0.4 * (-shear1) - 0.2 * (-shear2) - 0.09 * (-0.5 / 0.25),
    );
    assert_eq!(out[0], 0.0);
    assert_eq!(out[1], 0.0);

    let mut stationary = [0.0; 6];
    op.no_slip_wall_from_primitive(
        trace,
        0.25,
        [0.0; 2],
        WallThermalCondition::Adiabatic,
        2.0,
        [-3.0, 4.0],
        &mut stationary,
        || false,
    )
    .unwrap();
    assert_ne!(stationary[3], 0.0);
    assert_ne!(stationary[4], 0.0);
    assert_eq!(stationary[5], 0.0);
}
