//! Physical frame shapes, participation and signed prestressed modes.

use fs_solid::frame2d::{FrameMassFormulation, FrameModeStatus};
use fs_solid::{Frame2d, Geometry, SectionModel};

const E: f64 = 200e9;
const A: f64 = 0.01;
const I: f64 = 8e-5;
const L: f64 = 3.0;

fn section() -> SectionModel {
    SectionModel::Elastic {
        ea: E * A,
        ei: E * I,
    }
}

fn cantilever(
    elements: usize,
    geometry: Geometry,
    line_mass: f64,
    tip_mass: f64,
) -> (Frame2d, usize) {
    let mut frame =
        Frame2d::with_mass_formulation(geometry, FrameMassFormulation::ConsistentReference);
    let mut previous = frame.add_node(0.0, 0.0);
    frame.fix(previous, [true; 3]);
    for i in 1..=elements {
        let next = frame.add_node(L * i as f64 / elements as f64, 0.0);
        frame
            .add_element(previous, next, 5, line_mass, &section)
            .unwrap();
        previous = next;
    }
    frame.add_mass(previous, tip_mass);
    (frame, previous)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

#[test]
fn tip_mass_modes_restore_the_exact_massless_beam_shape() {
    let mass = 2000.0;
    let (frame, _) = cantilever(4, Geometry::Linear, 0.0, mass);
    let analysis = frame.modal_analysis().unwrap();
    assert_eq!(analysis.modes.len(), 2);
    let bending = &analysis.modes[0];
    let exact_lambda = 3.0 * E * I / (mass * L.powi(3));
    assert!((bending.eigenvalue_per_s2 / exact_lambda - 1.0).abs() < 1e-10);
    assert_eq!(bending.status(), FrameModeStatus::Oscillatory);
    assert_eq!(&bending.displacement[..3], &[0.0; 3]);
    for node in 1..=4 {
        let x = L * node as f64 / 4.0;
        let displacement = x * x * (3.0 * L - x) / (2.0 * L.powi(3) * mass.sqrt());
        let rotation = 3.0 * x * (2.0 * L - x) / (2.0 * L.powi(3) * mass.sqrt());
        assert!((bending.displacement[3 * node + 1] - displacement).abs() < 1e-11);
        assert!((bending.displacement[3 * node + 2] - rotation).abs() < 1e-11);
        assert!(bending.displacement[3 * node].abs() < 1e-12);
    }
    assert!(bending.base_participation[0].abs() < 1e-10);
    assert!((bending.base_participation[1] / mass.sqrt() - 1.0).abs() < 1e-11);
    assert!((bending.effective_mass_kg[1] / mass - 1.0).abs() < 1e-11);
    assert!((analysis.modes[1].effective_mass_kg[0] / mass - 1.0).abs() < 1e-11);
    assert!(bending.relative_residual < 1e-11);
}

#[test]
fn consistent_mass_modes_are_orthonormal_complete_and_work_conjugate() {
    let (frame, _) = cantilever(4, Geometry::Linear, 78.5, 17.0);
    let analysis = frame.modal_analysis().unwrap();
    assert_eq!(analysis.modes.len(), frame.ndof() - 3);
    let mass_shapes: Vec<Vec<f64>> = analysis
        .modes
        .iter()
        .map(|m| frame.inertial_forces(&m.displacement).unwrap())
        .collect();
    let influence: [Vec<f64>; 2] = [0, 1].map(|axis| {
        let r: Vec<f64> = (0..frame.ndof())
            .map(|d| if d % 3 == axis { 1.0 } else { 0.0 })
            .collect();
        frame.inertial_forces(&r).unwrap()
    });
    for (i, mode) in analysis.modes.iter().enumerate() {
        assert_eq!(mode.status(), FrameModeStatus::Oscillatory);
        assert!(mode.relative_residual < 1e-10);
        for (j, mass_shape) in mass_shapes.iter().enumerate() {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert!((dot(&mode.displacement, mass_shape) - expected).abs() < 1e-10);
        }
        for axis in 0..2 {
            assert!(
                (dot(&mode.displacement, &influence[axis]) - mode.base_participation[axis]).abs()
                    < 1e-10
            );
        }
    }
    let state: Vec<f64> = (0..frame.ndof())
        .map(|d| {
            if d < 3 {
                0.0
            } else {
                1e-3 * (d as f64 * 0.37).sin()
            }
        })
        .collect();
    let momentum = frame.inertial_forces(&state).unwrap();
    let mut reconstructed = vec![0.0; frame.ndof()];
    for mode in &analysis.modes {
        let coordinate = dot(&mode.displacement, &momentum);
        for (q, phi) in reconstructed.iter_mut().zip(&mode.displacement) {
            *q += coordinate * phi;
        }
    }
    for (q, original) in reconstructed.iter().zip(&state) {
        assert!((q - original).abs() < 1e-12);
    }
    for axis in 0..2 {
        let effective: f64 = analysis
            .modes
            .iter()
            .map(|m| m.effective_mass_kg[axis])
            .sum();
        assert!(effective > 0.5 * analysis.base_influence_mass_kg[axis]);
        assert!(effective <= (1.0 + 1e-10) * analysis.base_influence_mass_kg[axis]);
    }
}

#[test]
fn free_body_modes_keep_rigid_motion_unresolved_instead_of_clipping_signs() {
    let mut frame =
        Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    let a = frame.add_node(0.0, 0.0);
    let b = frame.add_node(L, 0.0);
    frame.add_element(a, b, 5, 78.5, &section).unwrap();
    let analysis = frame.modal_analysis().unwrap();
    assert_eq!(
        analysis
            .modes
            .iter()
            .filter(|m| m.status() == FrameModeStatus::Unresolved)
            .count(),
        3
    );
    assert_eq!(
        analysis
            .modes
            .iter()
            .filter(|m| m.status() == FrameModeStatus::Oscillatory)
            .count(),
        3
    );
    assert!(
        analysis
            .modes
            .iter()
            .all(|m| m.status() != FrameModeStatus::Unstable)
    );
    assert!(frame.natural_frequencies().is_err());
    for mode in &analysis.modes[..3] {
        assert!(mode.eigenvalue_per_s2.abs() <= mode.eigenvalue_resolution_per_s2);
        assert!(mode.angular_frequency().is_none());
        assert!(mode.growth_rate().is_none());
    }
}

#[test]
fn compressed_column_retains_the_signed_buckling_mode_and_growth_rate() {
    let mass = 2000.0;
    let p_euler = std::f64::consts::PI.powi(2) * E * I / (4.0 * L * L);
    let lambda_zero = 3.0 * E * I / (mass * L.powi(3));
    for ratio in [0.5, 1.2] {
        let (mut frame, tip) = cantilever(8, Geometry::Corotational, 0.0, mass);
        let mut load = vec![0.0; frame.ndof()];
        load[3 * tip] = -ratio * p_euler;
        frame.static_load(&load, 3).unwrap();
        let analysis = frame.modal_analysis().unwrap();
        let first = &analysis.modes[0];
        assert!(first.relative_residual < 1e-10);
        if ratio < 1.0 {
            assert_eq!(first.status(), FrameModeStatus::Oscillatory);
            assert!((first.eigenvalue_per_s2 / lambda_zero - 0.5).abs() < 0.03);
            assert!(frame.natural_frequencies().is_ok());
        } else {
            assert!(first.eigenvalue_per_s2 < -0.1 * lambda_zero);
            assert_eq!(first.status(), FrameModeStatus::Unstable);
            assert!(first.angular_frequency().is_none());
            assert_eq!(first.growth_rate(), Some((-first.eigenvalue_per_s2).sqrt()));
            assert!(frame.natural_frequencies().is_err());
            assert!((first.effective_mass_kg[1] / mass - 1.0).abs() < 1e-10);
        }
    }
}
