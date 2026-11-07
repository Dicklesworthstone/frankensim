//! G1/G3: spatial wall/pressure correlation, source work and original curved maps.

use super::super::time::exact::ExactTimeLimits;
use super::super::time::{HeterogeneousTimeLimits, ResistiveBoundary};
use super::super::{BoundaryAdmittance, HarmonicLimits};
use super::*;
use crate::acoustics::helmholtz::{HelmholtzFe, Medium, QuadraticGeometryControls};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

// G1: exact integral of (1+x) lambda_i lambda_j over the unit right triangle.
// This oracle does not call the FE surface or material quadrature.
fn resistive_power(pressure: &[f64]) -> f64 {
    let gram = [
        [1.0 / 10.0, 7.0 / 120.0, 1.0 / 20.0],
        [7.0 / 120.0, 2.0 / 15.0, 7.0 / 120.0],
        [1.0 / 20.0, 7.0 / 120.0, 1.0 / 10.0],
    ];
    (0..3)
        .flat_map(|i| (0..3).map(move |j| pressure[i] * gram[i][j] * pressure[j]))
        .sum()
}

#[test]
fn spatial_resistance_retains_pressure_correlation_and_midpoint_physical_heat() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let calls = Cell::new(0);
    let boundary = model
        .spatial_resistive_boundary(
            &[FACE],
            4,
            &mut |_, point| {
                calls.set(calls.get() + 1);
                Ok::<_, AcousticFeError>(1.0 + point.point.position_m[0])
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(boundary.faces(), &[FACE]);
    assert_eq!(boundary.quadrature_order(), 4);
    assert_eq!(calls.get(), 16);
    let harmonic_boundary = boundary.at_frequency(FREQUENCY).unwrap();
    let harmonic = model
        .prepare_with_spatial_admittances(
            FREQUENCY,
            &[],
            &[&harmonic_boundary],
            &[],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let harmonic_field = harmonic
        .solve(&[(FACE, [C64::ONE; 3])], &mut || false)
        .unwrap();
    let harmonic_heat = 0.5
        * (resistive_power(
            &harmonic_field
                .pressure_pa
                .iter()
                .map(|p| p.re)
                .collect::<Vec<_>>(),
        ) + resistive_power(
            &harmonic_field
                .pressure_pa
                .iter()
                .map(|p| p.im)
                .collect::<Vec<_>>(),
        ));
    close(harmonic_field.boundary_dissipation_w, harmonic_heat, 1e-12);
    close(harmonic_field.power_balance_residual_w(), 0.0, 1e-12);
    assert!(boundary.at_frequency(0.0).is_err());
    assert!(boundary.at_frequency(f64::NAN).is_err());
    let h = 0.004;
    let owner = model
        .prepare_time_with_spatial_boundaries(
            h,
            &[],
            &[&boundary],
            &[],
            HeterogeneousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let averaged = model
        .prepare_time(
            h,
            &[ResistiveBoundary {
                faces: vec![FACE],
                admittance_m_pa_s: 4.0 / 3.0,
            }],
            HeterogeneousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = owner.state(&[0.0; 4], &[1.0, -0.4, 2.0, 0.3]).unwrap();
    let average_step = averaged.step(&initial, &[], &mut || false).unwrap();
    let mut state = initial;
    let mut dissipated = 0.0;
    let initial_energy = owner.energy(&state).unwrap().total_j();
    for index in 0..30 {
        let step = owner.step(&state, &[], &mut || false).unwrap();
        let midpoint = state
            .pressure_pa()
            .iter()
            .zip(step.state.pressure_pa())
            .map(|(a, b)| 0.5 * (a + b))
            .collect::<Vec<_>>();
        let heat = h * resistive_power(&midpoint);
        close(
            step.work.boundary_dissipation_j,
            heat,
            2e-13 * initial_energy,
        );
        assert!(heat > 0.0);
        close(step.work.balance_residual_j(), 0.0, 2e-12 * initial_energy);
        if index == 0 {
            assert!(
                step.state
                    .pressure_pa()
                    .iter()
                    .zip(average_step.state.pressure_pa())
                    .any(|(a, b)| (a - b).abs() > 1e-3),
                "face averaging must change the solved wave"
            );
        }
        dissipated += heat;
        state = step.state;
    }
    close(
        owner.energy(&state).unwrap().total_j() + dissipated,
        initial_energy,
        2e-11 * initial_energy,
    );
    assert_eq!(calls.get(), 16, "warm steps must not resample the wall");
}

#[test]
fn spatial_resistance_exact_evolution_matches_midpoint_convergence_and_heat() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let boundary = model
        .spatial_resistive_boundary(
            &[FACE],
            4,
            &mut |_, p| Ok::<_, AcousticFeError>(1.0 + p.point.position_m[0]),
            &mut || false,
        )
        .unwrap();
    let duration = 0.04;
    let exact = model
        .prepare_exact_time_with_spatial_boundaries(
            duration,
            &[],
            &[&boundary],
            &[],
            &[],
            &[],
            &[],
            ExactTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = exact.state(&[0.0; 4], &[1.0, -0.4, 2.0, 0.3]).unwrap();
    let expected = exact.step(&initial, &[], &mut || false).unwrap();
    let energy = exact.energy(&initial).unwrap().total_j();
    assert!(expected.work.boundary_dissipation_j > 0.1 * energy);
    close(expected.work.balance_residual_j(), 0.0, 2e-10 * energy);
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let h = duration / steps as f64;
        let midpoint = model
            .prepare_time_with_spatial_boundaries(
                h,
                &[],
                &[&boundary],
                &[],
                HeterogeneousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut state = initial.clone();
        let mut heat = 0.0;
        for _ in 0..steps {
            let result = midpoint.step(&state, &[], &mut || false).unwrap();
            let p = state
                .pressure_pa()
                .iter()
                .zip(result.state.pressure_pa())
                .map(|(a, b)| 0.5 * (a + b))
                .collect::<Vec<_>>();
            heat += h * resistive_power(&p);
            state = result.state;
        }
        let field_error = state
            .pressure_pa()
            .iter()
            .zip(expected.state.pressure_pa())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push(field_error);
        close(
            heat,
            expected.work.boundary_dissipation_j,
            2.0 * energy / (steps * steps) as f64,
        );
    }
    assert!(errors[0] > 1e-5, "time refinement must be live");
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.28 * pair[0],
            "second-order field ladder: {errors:?}"
        );
    }
}

#[test]
fn spatial_resistance_constant_limit_and_original_constraints() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    for order in 1..=3 {
        let geometry = HelmholtzFe::new(&mesh, &POSITIONS, order, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let boundary = model
            .spatial_resistive_boundary(
                &[FACE],
                order + 2,
                &mut |_, _| Ok::<_, AcousticFeError>(0.7),
                &mut || false,
            )
            .unwrap();
        let uniform = [ResistiveBoundary {
            faces: vec![FACE],
            admittance_m_pa_s: 0.7,
        }];
        let spatial = model
            .prepare_time_with_spatial_boundaries(
                0.01,
                &[],
                &[&boundary],
                &[],
                HeterogeneousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let constant = model
            .prepare_time(
                0.01,
                &uniform,
                HeterogeneousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let p = (0..model.ndof())
            .map(|i| 1.0 / (i + 1) as f64)
            .collect::<Vec<_>>();
        let state = spatial.state(&vec![0.0; model.ndof()], &p).unwrap();
        let a = spatial.step(&state, &[], &mut || false).unwrap();
        let b = constant.step(&state, &[], &mut || false).unwrap();
        for (p, q) in a.state.pressure_pa().iter().zip(b.state.pressure_pa()) {
            close(*p, *q, 1e-9);
        }
        close(
            a.work.boundary_dissipation_j,
            b.work.boundary_dissipation_j,
            1e-12,
        );
        assert!(
            model
                .prepare_time_with_spatial_boundaries(
                    0.01,
                    &uniform,
                    &[&boundary],
                    &[],
                    HeterogeneousTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
        assert!(
            model
                .prepare_time_with_spatial_boundaries(
                    0.01,
                    &[],
                    &[&boundary, &boundary],
                    &[],
                    HeterogeneousTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
        assert!(
            model
                .prepare_time_with_spatial_boundaries(
                    0.01,
                    &[],
                    &[&boundary],
                    &[FACE],
                    HeterogeneousTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
        let other_model = super::tests::uniform(&geometry);
        assert!(
            other_model
                .prepare_time_with_spatial_boundaries(
                    0.01,
                    &[],
                    &[&boundary],
                    &[],
                    HeterogeneousTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
        assert!(
            model
                .prepare_exact_time_with_spatial_boundaries(
                    0.01,
                    &[],
                    &[&boundary],
                    &[],
                    &[],
                    &[],
                    &[FACE],
                    ExactTimeLimits::default(),
                    &mut || false,
                )
                .is_err()
        );
    }
}

#[test]
fn spatial_resistance_rejects_active_samples_and_cancels_without_publication() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .spatial_resistive_boundary(
                    &[FACE],
                    4,
                    &mut |_, _| Ok::<_, AcousticFeError>(value),
                    &mut || false,
                )
                .is_err()
        );
    }
    let calls = Cell::new(0);
    let result = model.spatial_resistive_boundary(
        &[FACE],
        4,
        &mut |_, _| {
            calls.set(calls.get() + 1);
            Ok::<_, AcousticFeError>(0.4)
        },
        &mut || calls.get() >= 3,
    );
    assert!(matches!(result, Err(AcousticFeError::Cancelled)));
    assert_eq!(calls.get(), 3);
    let boundary = model
        .spatial_resistive_boundary(
            &[FACE],
            4,
            &mut |_, _| Ok::<_, AcousticFeError>(0.4),
            &mut || false,
        )
        .unwrap();
    assert!(matches!(
        model.prepare_time_with_spatial_boundaries(
            0.01,
            &[],
            &[&boundary],
            &[],
            HeterogeneousTimeLimits::default(),
            &mut || true,
        ),
        Err(HeterogeneousError::Cancelled)
    ));
    assert!(matches!(
        model.prepare_exact_time_with_spatial_boundaries(
            0.01,
            &[],
            &[&boundary],
            &[],
            &[],
            &[],
            &[],
            ExactTimeLimits::default(),
            &mut || true,
        ),
        Err(HeterogeneousError::Cancelled)
    ));
}

#[test]
fn spatial_resistance_and_causal_walls_close_distinct_heat_and_exchange_channels() {
    use super::super::time::mechanical::exact::ExactMechanicalTimeLimits;
    use super::super::time::mechanical::{MechanicalBoundary, MechanicalTimeLimits};
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let boundary = model
        .spatial_resistive_boundary(
            &[FACE],
            4,
            &mut |_, p| Ok::<_, AcousticFeError>(1.0 + p.point.position_m[0]),
            &mut || false,
        )
        .unwrap();
    let walls = [MechanicalBoundary {
        faces: vec![[0, 1, 3]],
        wall: fs_phs::WallPin {
            surface_density: 0.4,
            stiffness_per_area: 1.7,
            resistance: 0.3,
        },
        relaxation: vec![fs_phs::acoustic_outlet::AcousticWallRelaxation {
            stiffness_per_area: 0.8,
            relaxation_time_s: 0.05,
        }],
    }];
    let duration = 0.04;
    let exact = model
        .prepare_exact_time_with_walls_and_spatial_boundaries(
            duration,
            &[],
            &[&boundary],
            &walls,
            &[],
            &[],
            &[],
            &[],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let pressure = [1.0, -0.4, 2.0, 0.3];
    let initial = exact.state(&[0.0; 4], &pressure).unwrap();
    let result = exact.step(&initial, &[], &mut || false).unwrap();
    let energy = exact.energy(&initial).unwrap().total_j();
    close(result.work.balance_residual_j(), 0.0, 1e-9 * energy);
    close(
        result.work.acoustic_balance_residual_j(),
        0.0,
        1e-9 * energy,
    );
    close(result.work.wall_balance_residual_j(), 0.0, 1e-9 * energy);
    assert!(result.work.boundary_dissipation_j > 0.1 * energy);
    assert!(result.work.wall_resistance_dissipation_j > 0.0);
    assert!(result.work.wall_relaxation_dissipation_j > 0.0);
    assert!(result.work.wall_work_j > 0.0);
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let h = duration / steps as f64;
        let midpoint = model
            .prepare_time_with_walls_and_spatial_boundaries(
                h,
                &[],
                &[&boundary],
                &walls,
                &[],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut state = midpoint.state(&[0.0; 4], &pressure).unwrap();
        let mut heat = 0.0;
        for _ in 0..steps {
            let next = midpoint.step(&state, &[], &mut || false).unwrap();
            let p = state
                .pressure_pa()
                .iter()
                .zip(next.state.pressure_pa())
                .map(|(a, b)| 0.5 * (a + b))
                .collect::<Vec<_>>();
            let independent_heat = h * resistive_power(&p);
            close(
                next.work.boundary_dissipation_j,
                independent_heat,
                1e-12 * energy,
            );
            close(next.work.acoustic_balance_residual_j(), 0.0, 1e-10 * energy);
            close(next.work.wall_balance_residual_j(), 0.0, 1e-10 * energy);
            heat += independent_heat;
            state = next.state;
        }
        errors.push(
            state
                .pressure_pa()
                .iter()
                .zip(result.state.pressure_pa())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
        );
        close(
            heat,
            result.work.boundary_dissipation_j,
            2.0 * energy / (steps * steps) as f64,
        );
    }
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.28 * pair[0],
            "coupled time refinement: {errors:?}"
        );
    }
    let conflicting = [MechanicalBoundary {
        faces: vec![FACE],
        ..walls[0].clone()
    }];
    assert!(
        model
            .prepare_time_with_walls_and_spatial_boundaries(
                duration,
                &[],
                &[&boundary],
                &conflicting,
                &[],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
    assert!(
        model
            .prepare_exact_time_with_walls_and_spatial_boundaries(
                duration,
                &[],
                &[&boundary],
                &conflicting,
                &[],
                &[],
                &[],
                &[],
                ExactMechanicalTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
}

const FREQUENCY: f64 = 0.4;
const MEDIUM: Medium = Medium {
    density_kg_m3: 2.0,
    sound_speed_m_s: 3.0,
};
const POSITIONS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const FACE: [u32; 3] = [0, 1, 2];

fn uniform<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn spatial<'m, 'f, 'c>(
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    faces: &[[u32; 3]],
    order: usize,
    law: impl Fn(GeometryFacePoint) -> C64,
) -> SpatialBoundaryAdmittance<'m, 'f, 'c> {
    model
        .spatial_boundary_admittance(
            FREQUENCY,
            faces,
            order,
            &mut |_, _, point| Ok::<_, AcousticFeError>(law(point)),
            &mut || false,
        )
        .unwrap()
}

fn entry(matrix: &Csr, row: usize, col: usize) -> f64 {
    let (cols, values) = matrix.row(row);
    cols.iter()
        .position(|&c| c == col)
        .map_or(0.0, |i| values[i])
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

// Exact barycentric triangle moment, independent of FE or quadrature code.
fn moment(exponents: [usize; 3]) -> f64 {
    let factorial = |n: usize| (1..=n).fold(1.0, |a, b| a * b as f64);
    exponents.into_iter().map(factorial).product::<f64>()
        / factorial(exponents.iter().sum::<usize>() + 2)
}

#[test]
fn affine_linear_admittance_matches_exact_triangle_moments_and_power() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let law = |point: GeometryFacePoint| {
        let [x, y, _] = point.point.position_m;
        C64::new(1.0 + x, 0.25 - y)
    };
    let field = spatial(&model, &[FACE], 4, law);
    let reversed = spatial(&model, &[[2, 0, 1]], 4, law);
    let expected = [
        [1.0 / 10.0, 7.0 / 120.0, 1.0 / 20.0],
        [7.0 / 120.0, 2.0 / 15.0, 7.0 / 120.0],
        [1.0 / 20.0, 7.0 / 120.0, 1.0 / 10.0],
    ];
    for i in 0..4 {
        for j in 0..4 {
            assert_eq!(
                entry(&field.conductance, i, j),
                entry(&reversed.conductance, i, j)
            );
            assert_eq!(
                entry(&field.susceptance, i, j),
                entry(&reversed.susceptance, i, j)
            );
            if i < 3 && j < 3 {
                close(entry(&field.conductance, i, j), expected[i][j], 3e-15);
                let mut exponent = [0; 3];
                exponent[i] += 1;
                exponent[j] += 1;
                let constant = moment(exponent);
                exponent[2] += 1;
                close(
                    entry(&field.susceptance, i, j),
                    0.25 * constant - moment(exponent),
                    3e-15,
                );
            } else {
                assert_eq!(entry(&field.conductance, i, j), 0.0);
            }
        }
    }
    let p = [C64::new(1.0, 0.3), C64::new(-0.4, 0.8), C64::new(2.0, -0.2)];
    let mut power = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            power +=
                0.5 * entry(&field.conductance, i, j) * (p[i].re * p[j].re + p[i].im * p[j].im);
        }
    }
    close(power, 0.3415, 8e-15);
    assert_eq!(field.faces(), &[FACE]);
    assert_eq!(field.frequency_hz(), FREQUENCY);
    assert_eq!(field.quadrature_order(), 4);
}

#[test]
fn constant_fields_recover_uniform_walls_at_every_pressure_order_without_resampling() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let y = C64::new(0.3, -0.1);
    for order in 1..=4 {
        let geometry = HelmholtzFe::new(&mesh, &POSITIONS, order, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let faces = geometry.boundary_faces();
        let calls = Cell::new(0);
        let field = model
            .spatial_boundary_admittance(
                FREQUENCY,
                &faces,
                order + 2,
                &mut |frequency, face, _| {
                    assert_eq!(frequency, FREQUENCY);
                    assert!(face[0] < face[1] && face[1] < face[2]);
                    calls.set(calls.get() + 1);
                    Ok::<_, AcousticFeError>(y)
                },
                &mut || false,
            )
            .unwrap();
        let expected_calls = faces.len() * (order + 2).pow(2);
        assert_eq!(calls.get(), expected_calls);
        let old = model
            .prepare(
                FREQUENCY,
                &[BoundaryAdmittance {
                    faces,
                    admittance_m_pa_s: y,
                }],
                HarmonicLimits::default(),
                &mut || false,
            )
            .unwrap();
        let new = model
            .prepare_with_spatial_admittances(
                FREQUENCY,
                &[],
                &[&field],
                &[],
                HarmonicLimits::default(),
                &mut || false,
            )
            .unwrap();
        let drive = [(FACE, [C64::new(1.0, 0.3); 3])];
        let expected = old.solve(&drive, &mut || false).unwrap();
        let actual = new.solve(&drive, &mut || false).unwrap();
        for (&p, &q) in actual.pressure_pa.iter().zip(&expected.pressure_pa) {
            assert!((p - q).abs() < 2e-10 * (1.0 + q.abs()));
        }
        close(
            actual.boundary_dissipation_w,
            expected.boundary_dissipation_w,
            2e-10 * expected.boundary_dissipation_w,
        );
        close(actual.power_balance_residual_w(), 0.0, 2e-11);
        assert_eq!(calls.get(), expected_calls);
    }
}

#[test]
fn graded_wall_solves_a_manufactured_field_and_closes_independent_surface_work() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let omega = core::f64::consts::TAU * FREQUENCY;
    let beta = 1.0 / (MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2));
    let y = |face: [u32; 3], x: [f64; 3]| {
        if face == FACE {
            C64::new(1.0 + x[0], 0.25 - x[1])
        } else {
            C64::new(0.7, 0.1)
        }
    };
    let graded = spatial(&model, &[FACE], 4, |point| y(FACE, point.point.position_m));
    let other_faces: Vec<_> = geometry
        .boundary_faces()
        .into_iter()
        .filter(|face| *face != FACE)
        .collect();
    let uniform_walls = [BoundaryAdmittance {
        faces: other_faces,
        admittance_m_pa_s: C64::new(0.7, 0.1),
    }];
    // Exact p=1+x Pa. Momentum gives v=i grad(p)/(omega rho); hence
    // inward imposed velocity = Y p - v.n. Its real/imaginary traces are
    // quadratic and use the original boundary-source owner.
    let mut real_traces = Vec::new();
    let mut imaginary_traces = Vec::new();
    for face in geometry.boundary_faces() {
        let (vertices, normal, _) = geometry.boundary_geometry(face).unwrap();
        let mut nodes = [[0.0; 3]; 6];
        nodes[..3].copy_from_slice(&vertices);
        for (i, (a, b)) in [(0, 1), (1, 2), (2, 0)].into_iter().enumerate() {
            nodes[i + 3] =
                core::array::from_fn(|axis| f64::midpoint(vertices[a][axis], vertices[b][axis]));
        }
        let values = nodes.map(|x| {
            y(face, x).scale(1.0 + x[0]) - C64::new(0.0, normal[0] / (omega * MEDIUM.density_kg_m3))
        });
        real_traces.push((face, values.map(|z| z.re)));
        imaginary_traces.push((face, values.map(|z| z.im)));
    }
    let real_source = model
        .quadratic_velocity_source(&real_traces, &mut || false)
        .unwrap();
    let imaginary_source = model
        .quadratic_velocity_source(&imaginary_traces, &mut || false)
        .unwrap();
    let volume_source = model
        .volume_velocity_source(
            4,
            &mut |_, x| Ok::<_, HeterogeneousError>(1.0 + x[0]),
            &mut || false,
        )
        .unwrap();
    let prepared = model
        .prepare_with_spatial_admittances(
            FREQUENCY,
            &uniform_walls,
            &[&graded],
            &[],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
    let boundary_sources = [
        (&real_source, C64::from_re(1.0)),
        (&imaginary_source, C64::new(0.0, 1.0)),
    ];
    let volume_sources = [(&volume_source, C64::new(0.0, omega * beta))];
    let actual = prepared
        .solve_with_boundary_sources(&[], &boundary_sources, &volume_sources, &mut || false)
        .unwrap();
    for (i, x) in POSITIONS.iter().enumerate() {
        assert!((actual.pressure_pa[i] - C64::from_re(1.0 + x[0])).abs() < 2e-11);
    }
    for p in &actual.pressure_pa[4..] {
        assert!(p.abs() < 1e-10);
    }
    // Integral over the graded xy face: int (1+x)^3 / 2 = 13/20.
    // The other faces have int (1+x)^2 = 11/12, 1/2, 11 sqrt(3)/12.
    let expected_work = 0.65 + 0.35 * (17.0 + 11.0 * 3.0_f64.sqrt()) / 12.0;
    close(actual.boundary_dissipation_w, expected_work, 2e-11);
    close(actual.input_power_w, expected_work, 2e-11);
    close(actual.volume_source_power_w, 0.0, 2e-11);
    close(actual.power_balance_residual_w(), 0.0, 2e-11);
    // Replacing the actual law by its face mean is a different operator.
    let mut averaged_walls = uniform_walls.to_vec();
    averaged_walls.push(BoundaryAdmittance {
        faces: vec![FACE],
        admittance_m_pa_s: C64::new(4.0 / 3.0, -1.0 / 12.0),
    });
    let averaged = model
        .prepare(
            FREQUENCY,
            &averaged_walls,
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap()
        .solve_with_boundary_sources(&[], &boundary_sources, &volume_sources, &mut || false)
        .unwrap();
    assert!(
        averaged
            .pressure_pa
            .iter()
            .zip(&actual.pressure_pa)
            .any(|(&a, &b)| (a - b).abs() > 1e-4)
    );
}

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig {
        chunk_bytes: 4096,
        max_chunk_bytes: 4096,
        limit_bytes: Some(1 << 20),
        ..fs_alloc::ArenaConfig::default()
    })
    .scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 93,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

#[test]
fn curved_wall_uses_original_positions_normals_and_measure() {
    with_cx(|cx| {
        let a = 0.8;
        let map = |q: [f64; 3]| [q[0], q[1], q[2] + a * q[0] * q[1]];
        for tet in [[0, 1, 2, 3], [3, 1, 0, 2]] {
            let mesh = TetComplex::from_tets(4, vec![tet]);
            let positions: Vec<_> = POSITIONS.iter().copied().map(&map).collect();
            let midsides: Vec<_> = mesh
                .edges
                .iter()
                .map(|edge| {
                    map(core::array::from_fn(|axis| {
                        f64::midpoint(
                            POSITIONS[edge[0] as usize][axis],
                            POSITIONS[edge[1] as usize][axis],
                        )
                    }))
                })
                .collect();
            let geometry = HelmholtzFe::new_quadratic(
                &mesh,
                &positions,
                &midsides,
                1,
                6,
                MEDIUM,
                QuadraticGeometryControls::default(),
                cx,
            )
            .unwrap();
            let model = uniform(&geometry);
            let field = spatial(&model, &[FACE], 6, |sample| {
                let [x, y, z] = sample.point.position_m;
                let n = sample.outward_normal;
                close(z, a * x * y, 4e-15);
                assert!(n[2] < 0.0);
                close(n[0] / -n[2], a * y, 4e-15);
                close(n[1] / -n[2], a * x, 4e-15);
                C64::from_re((1.0 + z) / sample.surface_jacobian_m2)
            });
            for i in 0..3 {
                for j in 0..3 {
                    let mut exponents = [0; 3];
                    exponents[i] += 1;
                    exponents[j] += 1;
                    let flat = moment(exponents);
                    exponents[1] += 1;
                    exponents[2] += 1;
                    close(
                        entry(&field.conductance, i, j),
                        flat + a * moment(exponents),
                        3e-14,
                    );
                }
            }
        }
    });
}

#[test]
fn spatial_wall_errors_cancellation_and_owner_frequency_conflicts_are_preserved() {
    #[derive(Debug, PartialEq)]
    enum SourceError {
        MissingWall,
        Interior(AcousticFeError),
    }
    impl From<AcousticFeError> for SourceError {
        fn from(error: AcousticFeError) -> Self {
            Self::Interior(error)
        }
    }
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let mut passive =
        |_: f64, _: [u32; 3], _: GeometryFacePoint| Ok::<_, AcousticFeError>(C64::from_re(0.2));
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(
            model
                .spatial_boundary_admittance(frequency, &[FACE], 3, &mut passive, &mut || false)
                .is_err()
        );
    }
    for order in [0, 1, MAX_QUADRATURE_ORDER + 1] {
        assert!(
            model
                .spatial_boundary_admittance(FREQUENCY, &[FACE], order, &mut passive, &mut || false)
                .is_err()
        );
    }
    for faces in [vec![], vec![FACE, [2, 1, 0]], vec![[0, 1, 99]]] {
        assert!(
            model
                .spatial_boundary_admittance(FREQUENCY, &faces, 3, &mut passive, &mut || false)
                .is_err()
        );
    }
    for y in [
        C64::from_re(-0.1),
        C64::new(0.1, f64::NAN),
        C64::from_re(f64::INFINITY),
    ] {
        assert!(
            model
                .spatial_boundary_admittance(
                    FREQUENCY,
                    &[FACE],
                    3,
                    &mut |_, _, _| Ok::<_, AcousticFeError>(y),
                    &mut || false,
                )
                .is_err()
        );
    }
    let original_error = model.spatial_boundary_admittance(
        FREQUENCY,
        &[FACE],
        3,
        &mut |_, _, _| Err::<C64, _>(SourceError::MissingWall),
        &mut || false,
    );
    assert_eq!(original_error.err(), Some(SourceError::MissingWall));
    let calls = Cell::new(0);
    let cancelled = model.spatial_boundary_admittance(
        FREQUENCY,
        &[FACE],
        3,
        &mut |_, _, _| {
            calls.set(calls.get() + 1);
            Ok::<_, AcousticFeError>(C64::from_re(0.2))
        },
        &mut || calls.get() != 0,
    );
    assert!(matches!(cancelled, Err(AcousticFeError::Cancelled)));
    assert_eq!(calls.get(), 1);
    let field = spatial(&model, &[FACE], 3, |_| C64::from_re(0.2));
    let other = uniform(&geometry);
    for (owner, frequency, release) in [
        (&other, FREQUENCY, vec![]),
        (&model, FREQUENCY * 2.0, vec![]),
        (&model, FREQUENCY, vec![FACE]),
    ] {
        assert!(
            owner
                .prepare_with_spatial_admittances(
                    frequency,
                    &[],
                    &[&field],
                    &release,
                    HarmonicLimits::default(),
                    &mut || false
                )
                .is_err()
        );
    }
    assert!(
        model
            .prepare_with_spatial_admittances(
                FREQUENCY,
                &[],
                &[&field, &field],
                &[],
                HarmonicLimits::default(),
                &mut || false
            )
            .is_err()
    );
    assert!(
        model
            .prepare_with_spatial_admittances(
                FREQUENCY,
                &[BoundaryAdmittance {
                    faces: vec![FACE],
                    admittance_m_pa_s: C64::from_re(0.1)
                }],
                &[&field],
                &[],
                HarmonicLimits::default(),
                &mut || false
            )
            .is_err()
    );
    // Adjacent faces may share constrained DOFs; only physical-face overlap refuses.
    model
        .prepare_with_spatial_admittances(
            FREQUENCY,
            &[],
            &[&field],
            &[[0, 1, 3]],
            HarmonicLimits::default(),
            &mut || false,
        )
        .unwrap();
}
