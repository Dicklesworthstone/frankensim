use super::super::tests::{boundary_options, options as harmonic_options, source};
use super::*;
use fs_exec::CancelGate;
use fs_phs::RelaxationBranch;
use fs_plate::shell::reduction::ReductionBudget;
use fs_plate::shell::{ShellMesh, ShellSupport, assemble_shell_sections};
use fs_plate::{ModePair, PlateSection};

fn settings(reference: f64) -> OpenScreenResonanceOptions {
    OpenScreenResonanceOptions {
        boundary: ScreenComplexLimits {
            maximum_matrix_entries: 2 * 256 * 256,
            maximum_workspace_bytes: 4_000_000,
            maximum_work: 1_000_000_000,
            maximum_growth_exponent: 1.0,
        },
        maximum_dimension: 64,
        maximum_dense_bytes: 8_000_000,
        maximum_projection_work: 1_000_000,
        maximum_operator_evaluations: 256,
        maximum_result_bytes: 1_000_000,
        reference_angular_frequency_rad_s: reference,
        physical_residual_tolerance: 1e-10,
    }
}
fn pole_options(scale: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: scale,
        step_tolerance: 1e-12,
        residual_tolerance: 1e-12,
        max_iterations: 32,
        max_dimension: 64,
        max_dense_bytes: 8_000_000,
        max_work: 100_000_000,
    }
}
fn contour_options(scale: f64, nodes: usize, dimension: usize) -> BeynOptions {
    BeynOptions {
        nodes,
        probe_columns: dimension,
        rank_relative_tolerance: 1e-8,
        count_tolerance: 0.02,
        quadrature_tolerance: 0.02,
        contour_min_reciprocal_condition: 1e-13,
        linear_residual_tolerance: 1e-11,
        svd_tolerance: 1e-13,
        svd_max_sweeps: 128,
        eigen_max_sweeps: 256,
        max_dimension: 64,
        max_dense_bytes: 16_000_000,
        max_work: 100_000_000,
        refinement: pole_options(scale),
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()) + 1e-18,
        "{a:?} != {b:?}"
    );
}

// One actual planar source eigenmode, selected by its original physical
// translation. This creates no fabricated shell stiffness or BEM trace.
fn planar_source(transverse: bool, density: f64, sound_speed: f64) -> OpenScreenLoadedModes {
    let mesh = ShellMesh::new(
        vec![[0., 0., 0.], [0.1, 0., 0.], [0., 0.1, 0.], [0.08, 0.09, 0.]],
        vec![[0, 1, 3], [0, 3, 2]],
    )
    .unwrap();
    let sections = vec![PlateSection::isotropic(3.0e6, 0.3, 0.001, 1100.).unwrap(); 2];
    let model =
        assemble_shell_sections(&mesh, &sections, &[0, 1, 2], ShellSupport::Clamped).unwrap();
    let n = model.free;
    let mut k = vec![0.; n * n];
    let mut mass = k.clone();
    for i in 0..n {
        for j in 0..n {
            k[i * n + j] = model.k.get(i, j);
            mass[i * n + j] = model.m.get(i, j);
        }
    }
    let all = fs_modal::eigh_gen_dense(&k, &mass, n).unwrap();
    let translation =
        |mode: &ModePair, axis: usize| model.dof_map[18 + axis].map_or(0.0, |i| mode.phi[i]);
    let selected = all
        .into_iter()
        .filter(|mode| transverse || translation(mode, 2) == 0.0)
        .max_by(|a, b| {
            let score = |mode: &ModePair| {
                if transverse {
                    translation(mode, 2).abs()
                } else {
                    translation(mode, 0).abs() + translation(mode, 1).abs()
                }
            };
            score(a).total_cmp(&score(b))
        })
        .unwrap();
    let omega = selected.lambda.sqrt();
    let shell = ShellReduction::new(
        &mesh,
        &sections,
        &model,
        &[selected],
        ReductionBudget {
            max_modes: 4,
            max_facet_modes: 32,
            relative_tolerance: 1e-6,
        },
    )
    .unwrap();
    let damping = ModalPlateViscosity::from_psd_modal_matrix(
        vec![0.04 * omega],
        1,
        1.,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    OpenScreenLoadedModes::from_shell(
        &shell,
        Medium {
            density,
            sound_speed,
        },
        Some(&damping),
        OpenScreenModesOptions {
            maximum_modes: 4,
            source_surface: OpenScreenRadiationBudget {
                maximum_vertices: 64,
                maximum_facets: 32,
                maximum_corner_weights: 96,
            },
            acoustic_refinements: 1,
            maximum_source_bytes: 1_000_000,
            maximum_setup_work: 1_000_000,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn g3_augmented_real_axis_rows_reproduce_the_actual_full_mutual_harmonic_response() {
    let source = source();
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let pencil = prepared.resonance_pencil(settings(450.), || false).unwrap();
    let force = [C64::new(0.7, -0.2), C64::new(-0.4, 0.3)];
    let response = prepared
        .harmonic_response_batch(450., &[&force], harmonic_options(), || false)
        .unwrap();
    let raw = pencil
        .assemble(C64::from_re(450.), &mut 0, &mut || false)
        .unwrap();
    let jump = response.pressure_jump(0, || false).unwrap();
    let mut state = response.responses()[0].modal_displacement.clone();
    let mut interior = vec![C64::ZERO; pencil.pressure_unknowns()];
    for (&p, &dof) in jump.iter().zip(pencil.pressure_dof_map()) {
        if let Some(i) = dof {
            interior[i] = p;
        }
    }
    state.extend(interior);
    let residual = raw
        .physical_residuals(
            2,
            pencil.pressure_unknowns(),
            &state,
            Some(&force),
            1e-10,
            &mut || false,
        )
        .unwrap();
    assert!(residual[0] < 1e-10 && residual[1] < 1e-10);
    let d = pencil.dimension();
    for i in 0..2 {
        for j in 0..pencil.pressure_unknowns() {
            let g = pencil.velocity_pressure_pairing()[i][j];
            assert_eq!(raw.value.matrix[i * d + 2 + j], C64::from_re(g));
            close(
                raw.value.matrix[(2 + j) * d + i],
                C64::from_re(-source.medium().density * 450_f64.powi(2) * g),
                1e-15,
            );
        }
    }
    // Original complete equations, not an H^-1 quotient hidden in the contour.
    assert_eq!(
        pencil.dimension(),
        source.mode_count() + prepared.boundary().unknowns()
    );
}

#[test]
fn g1_open_pencil_analytic_derivative_and_zero_limit_keep_original_work_coordinates() {
    let source = source();
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let pencil = prepared.resonance_pencil(settings(450.), || false).unwrap();
    let w = C64::new(450., -2.);
    let exact = pencil.original_pencil(w, || false).unwrap();
    for step in [C64::from_re(0.01), C64::new(0., 0.01)] {
        let plus = pencil.original_pencil(w + step, || false).unwrap();
        let minus = pencil.original_pencil(w - step, || false).unwrap();
        for i in 0..exact.matrix.len() {
            let observed = (plus.matrix[i] - minus.matrix[i]) / step.scale(2.);
            let roundoff =
                128. * f64::EPSILON * (plus.matrix[i].abs() + minus.matrix[i].abs()) / step.abs();
            assert!(
                (observed - exact.derivative[i]).abs()
                    <= 5e-7 * exact.derivative[i].abs() + roundoff + 1e-18
            );
        }
    }
    let zero = pencil.original_pencil(C64::ZERO, || false).unwrap();
    let d = pencil.dimension();
    for i in 0..2 {
        for j in 0..2 {
            assert_eq!(
                zero.matrix[i * d + j],
                C64::from_re(source.stiffness_per_s2()[i * 2 + j])
            );
            assert_eq!(
                zero.derivative[i * d + j],
                C64::new(0., -source.viscosity_per_s()[i * 2 + j])
            );
        }
    }
    assert!(
        zero.matrix
            .iter()
            .chain(&zero.derivative)
            .all(|v| finite(*v))
    );
}

#[test]
fn g1_in_plane_shell_mode_keeps_exact_decay_and_survives_three_native_contour_rungs() {
    let source = planar_source(false, 1000., 1500.);
    assert!(
        source.normal_velocity_weights()[0]
            .iter()
            .flatten()
            .all(|v| *v == 0.)
    );
    let k = source.stiffness_per_s2()[0];
    let c = source.viscosity_per_s()[0];
    let expected = C64::new((k - 0.25 * c * c).sqrt(), -0.5 * c);
    let prepared = source
        .prepare(2.2 * expected.abs(), boundary_options(), || false)
        .unwrap();
    let pencil = prepared
        .resonance_pencil(settings(expected.re), || false)
        .unwrap();
    let local = pencil
        .resonance_pole(
            expected + C64::new(0.05 * expected.re, -0.01 * expected.re),
            &[C64::ONE],
            &pole_options(expected.re),
            || false,
        )
        .unwrap();
    close(local.angular_frequency_rad_s, expected, 2e-10);
    assert!(local.pressure_jump_pa.iter().all(|v| v.abs() < 1e-18));
    assert!((local.decay_rate_per_s() - 0.5 * c).abs() < 1e-9 * c);
    assert!((local.quality_factor().unwrap() - expected.re / c).abs() < 1e-9 * expected.re / c);
    let region = Ellipse {
        center: expected + C64::from_re(0.2 * expected.re),
        real_radius: 0.7 * expected.re,
        imaginary_radius: 0.6 * expected.re,
    };
    for nodes in [8, 16, 32] {
        let result = pencil
            .resonance_spectrum(
                region,
                &contour_options(expected.re, nodes, pencil.dimension()),
                || false,
            )
            .unwrap();
        assert_eq!(result.argument_count, 1);
        assert_eq!(result.rank, 1);
        assert_eq!(result.poles.len(), 1);
        assert_eq!(result.nodes, nodes);
        close(result.poles[0].angular_frequency_rad_s, expected, 2e-9);
        assert!(
            result.poles[0].mechanical_residual < 1e-10
                && result.poles[0].boundary_residual < 1e-10
        );
        assert!(result.operator_evaluations <= settings(expected.re).maximum_operator_evaluations);
    }
}

#[test]
fn g1_radiating_shell_pole_approaches_the_independent_incompressible_added_mass_equation() {
    // The same finite spatial pressure operator has an exact incompressible
    // limit. This checks reaction sign/units and true frequency shift; it is
    // not a continuum-disc or measured shell-frequency validation.
    let source = planar_source(true, 1000., 1e8);
    let k = source.stiffness_per_s2()[0];
    let c = source.viscosity_per_s()[0];
    let dry = k.sqrt();
    let prepared = source
        .prepare(1.5 * dry, boundary_options(), || false)
        .unwrap();
    let pencil = prepared.resonance_pencil(settings(dry), || false).unwrap();
    assert_eq!(pencil.pressure_unknowns(), 1);
    let static_op = pencil.original_pencil(C64::ZERO, || false).unwrap();
    let h0 = static_op.matrix[3];
    assert_eq!(h0.im, 0.);
    assert!(h0.re < 0.);
    let g = pencil.velocity_pressure_pairing()[0][0];
    let added_mass = -source.medium().density * g * g / h0.re;
    assert!(added_mass > 0.01);
    let total_mass = 1. + added_mass;
    let decay = c / (2. * total_mass);
    let expected = C64::new((k / total_mass - decay * decay).sqrt(), -decay);
    let pole = pencil
        .resonance_pole(
            expected + C64::new(0.01 * expected.re, 0.),
            &[C64::ONE],
            &pole_options(dry),
            || false,
        )
        .unwrap();
    close(pole.angular_frequency_rad_s, expected, 2e-7);
    assert!(pole.angular_frequency_rad_s.re < 0.995 * dry);
    assert!(pole.pressure_jump_pa.iter().any(|v| v.abs() > 0.));
    for (&jump, &dof) in pole.pressure_jump_pa.iter().zip(pencil.pressure_dof_map()) {
        if dof.is_none() {
            assert_eq!(jump, C64::ZERO);
        }
    }
    assert!(pole.mechanical_residual < 1e-10 && pole.boundary_residual < 1e-10);
}

#[test]
fn g0_open_resonance_caps_whole_complex_domains_and_cancelled_queries_publish_nothing() {
    let source = source();
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let mut too_small = settings(450.);
    too_small.maximum_dimension = 2;
    assert!(prepared.resonance_pencil(too_small, || false).is_err());
    let mut too_small = settings(450.);
    too_small.maximum_dense_bytes = 0;
    assert!(prepared.resonance_pencil(too_small, || false).is_err());
    let pencil = prepared.resonance_pencil(settings(450.), || false).unwrap();
    let raw = pencil
        .original_pencil(C64::new(450., -1.), || false)
        .unwrap();
    assert!(matches!(
        pencil.original_pencil(C64::new(450., -1.), || true),
        Err(BemModesError::Cancelled)
    ));
    let replay = pencil
        .original_pencil(C64::new(450., -1.), || false)
        .unwrap();
    assert_eq!(raw.matrix, replay.matrix);
    assert_eq!(raw.derivative, replay.derivative);
    let region = Ellipse {
        center: C64::new(790., -1.),
        real_radius: 20.,
        imaginary_radius: 1.,
    };
    assert!(matches!(
        pencil.resonance_spectrum(
            region,
            &contour_options(450., 8, pencil.dimension()),
            || false
        ),
        Err(BemModesError::ScreenBoundary(_))
    ));
    let mut tiny_growth = settings(450.);
    tiny_growth.boundary.maximum_growth_exponent = 1e-6;
    let chart = prepared.resonance_pencil(tiny_growth, || false).unwrap();
    assert!(matches!(
        chart.original_pencil(C64::new(450., -10.), || false),
        Err(BemModesError::ScreenBoundary(_))
    ));
    let mut one = settings(450.);
    one.maximum_operator_evaluations = 1;
    let one = prepared.resonance_pencil(one, || false).unwrap();
    assert!(matches!(
        one.resonance_pole(
            C64::new(450., -1.),
            &[C64::ONE, C64::ZERO],
            &pole_options(450.),
            || false
        ),
        Err(BemModesError::EvaluationBudget { maximum: 1 })
    ));
    let mut output = settings(450.);
    output.maximum_result_bytes =
        16 * (source.mode_count() + source.surface().vertices().len()) + 256 + 16 + 512 - 1;
    let output = prepared.resonance_pencil(output, || false).unwrap();
    assert!(
        output
            .resonance_pole(
                C64::new(450., -1.),
                &[C64::ONE, C64::ZERO],
                &pole_options(450.),
                || false
            )
            .is_err()
    );
    assert!(
        pencil
            .resonance_pole(
                C64::new(450., -1.),
                &[C64::ZERO, C64::ZERO],
                &pole_options(450.),
                || false
            )
            .is_err()
    );
}

fn with_memory(
    source: OpenScreenLoadedModes,
    branches: &[RelaxationBranch],
) -> OpenScreenLoadedModes {
    source
        .with_relaxation_branches(
            branches,
            OpenScreenMemoryLimits {
                maximum_branches: 8,
                maximum_payload_bytes: 256_000,
                maximum_work: 1_000_000,
            },
            || false,
        )
        .unwrap()
}

fn signed_memory_source() -> OpenScreenLoadedModes {
    with_memory(
        source(),
        &[
            RelaxationBranch {
                projection: vec![1., -0.3],
                stiffness: 40_000.,
                relaxation_time_s: 0.001,
            },
            RelaxationBranch {
                projection: vec![-0.2, 1.],
                stiffness: 20_000.,
                relaxation_time_s: 0.004,
            },
        ],
    )
}

#[test]
fn g1_retained_screen_memory_matches_driven_rows_and_analytic_frequency_derivatives() {
    let source = signed_memory_source();
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let pencil = prepared.resonance_pencil(settings(450.), || false).unwrap();
    assert_eq!(pencil.memory_count(), 2);
    assert_eq!(
        pencil.dimension(),
        source.mode_count() + pencil.pressure_unknowns() + 2
    );
    let force = [C64::new(0.7, -0.2), C64::new(-0.4, 0.3)];
    let driven = prepared
        .harmonic_response_batch(450., &[&force], harmonic_options(), || false)
        .unwrap();
    let response = &driven.responses()[0];
    let jump = driven.pressure_jump(0, || false).unwrap();
    let mut interior = vec![C64::ZERO; pencil.pressure_unknowns()];
    for (&value, &dof) in jump.iter().zip(pencil.pressure_dof_map()) {
        if let Some(i) = dof {
            interior[i] = value;
        }
    }
    let mut state = response.modal_displacement.clone();
    state.extend(interior);
    state.extend_from_slice(&response.elastic_memory_sqrt_j);
    let raw = pencil
        .assemble(C64::from_re(450.), &mut 0, &mut || false)
        .unwrap();
    let residuals = raw
        .physical_residuals(
            2,
            pencil.pressure_unknowns(),
            &state,
            Some(&force),
            1e-10,
            &mut || false,
        )
        .unwrap();
    assert!(residuals.iter().all(|r| *r < 1e-10));
    assert!(response.power.relaxation_w > 0.);
    // Check the original first-order material equation independently of the
    // coupled matrix assembly, including the signs of both source projections.
    let s = C64::new(0., -450.);
    for (branch, (&tau, row)) in source
        .relaxation_times_s()
        .iter()
        .zip(source.relaxation_rows_per_s())
        .enumerate()
    {
        let strain = row
            .iter()
            .zip(&response.modal_displacement)
            .fold(C64::ZERO, |sum, (&c, &q)| sum + q.scale(c));
        close(
            (C64::ONE + s.scale(tau)) * response.elastic_memory_sqrt_j[branch],
            s.scale(tau) * strain,
            2e-11,
        );
    }
    // Include a bare relaxation rate: the retained pencil and its derivative
    // stay finite there even though a rationally eliminated material would not.
    for w in [C64::new(450., -2.), C64::new(0., -250.)] {
        let exact = pencil.original_pencil(w, || false).unwrap();
        for step in [C64::from_re(0.001), C64::new(0., 0.001)] {
            let plus = pencil.original_pencil(w + step, || false).unwrap();
            let minus = pencil.original_pencil(w - step, || false).unwrap();
            for i in 0..exact.matrix.len() {
                let observed = (plus.matrix[i] - minus.matrix[i]) / step.scale(2.);
                let roundoff = 128. * f64::EPSILON
                    * (plus.matrix[i].abs() + minus.matrix[i].abs())
                    / step.abs();
                assert!(
                    (observed - exact.derivative[i]).abs()
                        <= 5e-7 * exact.derivative[i].abs() + roundoff + 1e-18
                );
            }
        }
    }
}

#[test]
fn g1_original_shell_with_memory_recovers_all_material_roots_and_crosses_the_bare_rate() {
    let base = planar_source(false, 1000., 1e8);
    assert!(
        base.normal_velocity_weights()[0]
            .iter()
            .flatten()
            .all(|v| *v == 0.)
    );
    let k = base.stiffness_per_s2()[0];
    let c = base.viscosity_per_s()[0];
    let a = k.sqrt();
    let rate = 1.5 * a;
    let b = rate + c - a;
    let d = k * rate / a;
    let excess = (rate - a) * (k / a + a - c);
    // For this actual in-plane shell mode, choose a positive Maxwell arm so
    // the independent characteristic polynomial factors as
    // (s+a)(s^2+b*s+d). The original source stiffness/viscosity are untouched.
    assert!(excess > 0. && d > 0.25 * b * b);
    let source = with_memory(
        base,
        &[RelaxationBranch {
            projection: vec![1.],
            stiffness: excess,
            relaxation_time_s: 1. / rate,
        }],
    );
    let prepared = source.prepare(2.5 * a, boundary_options(), || false).unwrap();
    let pencil = prepared.resonance_pencil(settings(a), || false).unwrap();
    for expected in [
        C64::new((d - 0.25 * b * b).sqrt(), -0.5 * b),
        C64::new(-(d - 0.25 * b * b).sqrt(), -0.5 * b),
        C64::new(0., -a),
    ] {
        let pole = pencil
            .resonance_pole(
                expected + C64::new(0.01 * a, 0.01 * a),
                &[C64::ONE],
                &pole_options(a),
                || false,
            )
            .unwrap();
        close(pole.angular_frequency_rad_s, expected, 3e-9);
        assert_eq!(pole.elastic_memory_sqrt_j.len(), 1);
        assert_eq!(pole.pressure_jump_pa.len(), source.surface().vertices().len());
        assert!(pole.pressure_jump_pa.iter().all(|v| v.abs() < 1e-14));
        assert!(
            pole.mechanical_residual < 1e-10
                && pole.boundary_residual < 1e-10
                && pole.memory_residual < 1e-10
        );
    }
    let region = Ellipse {
        center: C64::new(0., -a),
        real_radius: 0.2 * a,
        imaginary_radius: 0.7 * a,
    };
    let spectrum = pencil
        .resonance_spectrum(
            region,
            &contour_options(a, 64, pencil.dimension()),
            || false,
        )
        .unwrap();
    // The ellipse contains the real coupled -i*a pole AND the bare -i*rate
    // material rate. The latter is regular and must not alter the root count.
    assert_eq!(spectrum.argument_count, 1);
    assert_eq!(spectrum.rank, 1);
    assert_eq!(spectrum.poles.len(), 1);
    close(
        spectrum.poles[0].angular_frequency_rad_s,
        C64::new(0., -a),
        3e-9,
    );
    assert_eq!(spectrum.poles[0].elastic_memory_sqrt_j.len(), 1);
    let bare = pencil.original_pencil(C64::new(0., -rate), || false).unwrap();
    assert!(bare.matrix.iter().chain(&bare.derivative).all(|z| finite(*z)));
}

#[test]
fn g0_equal_time_screen_memories_keep_mechanically_hidden_relaxation_states() {
    let source = with_memory(
        source(),
        &[
            RelaxationBranch {
                projection: vec![1., -2.],
                stiffness: 9.,
                relaxation_time_s: 0.01,
            },
            RelaxationBranch {
                projection: vec![2., -4.],
                stiffness: 9.,
                relaxation_time_s: 0.01,
            },
        ],
    );
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let pencil = prepared.resonance_pencil(settings(450.), || false).unwrap();
    let bare = pencil.original_pencil(C64::new(0., -100.), || false).unwrap();
    let first = pencil.mode_count() + pencil.pressure_unknowns();
    let n = pencil.dimension();
    assert_eq!(n, first + 2);
    // q=jump=0, e=(2,-1) has zero net mechanical force but nonzero stored
    // energy. It is an actual retained relaxation mode at the common rate.
    for row in 0..n {
        let residual = bare.matrix[row * n + first].scale(2.)
            - bare.matrix[row * n + first + 1];
        assert!(residual.abs() < 1e-12, "row={row}: {residual:?}");
    }
    // At an ordinary spectral point the independent Schur image keeps all
    // signed cross terms; the production contour never forms this quotient.
    let w = C64::new(450., -2.);
    let raw = pencil.original_pencil(w, || false).unwrap();
    let s_tau = C64::new(w.im, -w.re).scale(0.01);
    let factor = s_tau / (C64::ONE + s_tau);
    let c = [3., -6.];
    for i in 0..2 {
        for j in 0..2 {
            let observed = -(raw.matrix[i * n + first] * raw.matrix[first * n + j]
                / raw.matrix[first * n + first]
                + raw.matrix[i * n + first + 1] * raw.matrix[(first + 1) * n + j]
                    / raw.matrix[(first + 1) * n + first + 1]);
            close(observed, factor.scale(5. * c[i] * c[j]), 2e-12);
        }
    }
}

#[test]
fn g0_memory_pole_admission_counts_each_state_and_cancellation_preserves_replay() {
    let source = signed_memory_source();
    let prepared = source.prepare(800., boundary_options(), || false).unwrap();
    let m = source.mode_count();
    let p = prepared.boundary().unknowns();
    let dry_dimension = m + p;
    let mut dimension = settings(450.);
    dimension.maximum_dimension = dry_dimension;
    assert!(prepared.resonance_pencil(dimension, || false).is_err());
    let mut payload = settings(450.);
    payload.maximum_dense_bytes = 224 * dry_dimension * dry_dimension
        + 512 * dry_dimension
        + 8 * m * p
        + 16 * source.surface().vertices().len()
        + 24 * m;
    assert!(prepared.resonance_pencil(payload, || false).is_err());
    let mut work = settings(450.);
    work.maximum_projection_work = (16 * dry_dimension * dry_dimension + 8 * m * p) as u64;
    assert!(prepared.resonance_pencil(work, || false).is_err());
    let mut output = settings(450.);
    output.maximum_result_bytes =
        16 * (m + source.surface().vertices().len()) + 256 + 16 + 512;
    let output = prepared.resonance_pencil(output, || false).unwrap();
    assert!(
        output
            .resonance_pole(
                C64::new(450., -1.),
                &[C64::ONE, C64::ZERO],
                &pole_options(450.),
                || false,
            )
            .is_err()
    );
    let pencil = prepared.resonance_pencil(settings(450.), || false).unwrap();
    let w = C64::new(450., -2.);
    let original = pencil.original_pencil(w, || false).unwrap();
    assert!(matches!(
        pencil.original_pencil(w, || true),
        Err(BemModesError::Cancelled)
    ));
    let replay = pencil.original_pencil(w, || false).unwrap();
    assert_eq!(original.matrix, replay.matrix);
    assert_eq!(original.derivative, replay.derivative);
    assert!(matches!(
        pencil.resonance_pole(
            C64::new(0., -250.),
            &[C64::ONE, C64::ZERO],
            &pole_options(450.),
            || false,
        ),
        Err(BemModesError::SeedMemorySingular { branch: 1 })
    ));
}
