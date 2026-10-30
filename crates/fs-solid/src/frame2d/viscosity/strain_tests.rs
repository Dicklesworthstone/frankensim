//! Independent integrated axial/Hermite-curvature work on original geometry.
use super::*;

fn frame(length: f64, angle: f64) -> Frame2d {
    let mut frame = Frame2d::new(Geometry::Linear);
    frame.add_node(0.0, 0.0);
    frame.add_node(length * angle.cos(), length * angle.sin());
    frame.add_element(0, 1, 3, 1.0, &|| SectionModel::Elastic { ea: 7.0, ei: 3.0 }).unwrap();
    frame
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-12 * (1.0 + a.abs() + b.abs()), "{a:e} != {b:e}");
}

#[test]
fn strain_factors_recover_independently_integrated_full_stiffness() {
    let points = [0.5 - (3.0_f64 / 5.0).sqrt() / 2.0, 0.5,
        0.5 + (3.0_f64 / 5.0).sqrt() / 2.0];
    let weights = [5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0];
    for length in [0.3, 1.0, 2.7] {
        for angle in [0.0_f64, 0.37, 1.2] {
            let frame = frame(length, angle);
            let (dofs, b, l) = frame.reference_strain_traces(0).unwrap();
            assert_eq!(dofs, [0, 1, 2, 3, 4, 5]);
            let (c, s) = (angle.cos(), angle.sin());
            let axial = [-c / l, -s / l, 0.0, c / l, s / l, 0.0];
            let mut integrated = [[0.0; 6]; 6];
            for (xi, weight) in points.into_iter().zip(weights) {
                // Second derivatives of the original cubic transverse trace,
                // not the three exported factor columns.
                let hi = (12.0 * xi - 6.0) / l / l;
                let hj = -hi;
                let curvature = [-s * hi, c * hi, (6.0 * xi - 4.0) / l,
                    -s * hj, c * hj, (6.0 * xi - 2.0) / l];
                for i in 0..6 {
                    for j in 0..6 {
                        integrated[i][j] += l * weight
                            * (7.0 * axial[i] * axial[j] + 3.0 * curvature[i] * curvature[j]);
                    }
                }
            }
            for i in 0..6 {
                for j in 0..6 {
                    let factored = 7.0 / l * b[0][i] * b[0][j]
                        + 3.0 / l / l / l * (3.0 * b[1][i] * b[1][j] + b[2][i] * b[2][j]);
                    near(factored, integrated[i][j]);
                }
            }
        }
    }
}

#[test]
fn original_strain_maps_annihilate_rigid_motion_and_keep_signed_moment_work() {
    let angle = 0.47_f64;
    let length = 2.3;
    let frame = frame(length, angle);
    let (_, b, l) = frame.reference_strain_traces(0).unwrap();
    near(l, length);
    let theta = 0.03;
    for rigid in [[0.2, -0.4, 0.0, 0.2, -0.4, 0.0],
        [0.0, 0.0, theta, -angle.sin() * length * theta, angle.cos() * length * theta, theta]] {
        for column in b { near(column.iter().zip(rigid).map(|(a, v)| a * v).sum(), 0.0); }
    }
    assert_eq!(b[2][2], l);
    assert_eq!(b[2][5], -l);
    let velocity = [0.2, -0.1, 0.03, -0.2, 0.3, -0.04];
    let forces = [-0.7, 0.2, -0.4];
    let nodal: [f64; 6] = core::array::from_fn(|i|
        (0..3).map(|k| forces[k] * b[k][i]).sum());
    near(nodal.iter().zip(velocity).map(|(f, v)| f * v).sum(),
        (0..3).map(|k| forces[k] * b[k].iter().zip(velocity).map(|(a, v)| a * v).sum::<f64>()).sum());
}

#[test]
fn strain_queries_refuse_absent_or_nonlinear_geometry_without_mutation() {
    let mut frame = frame(1.0, 0.0);
    let original = frame.displacements().to_vec();
    assert!(frame.reference_strain_traces(1).is_err());
    assert_eq!(frame.displacements(), original.as_slice());
    frame.geometry = Geometry::Corotational;
    assert!(frame.reference_strain_traces(0).is_err());
    assert_eq!(frame.displacements(), original.as_slice());
}
