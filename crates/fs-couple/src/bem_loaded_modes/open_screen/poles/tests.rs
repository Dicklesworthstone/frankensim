use super::super::tests::{boundary_options, options as harmonic_options, source};
use super::*;
use fs_exec::CancelGate;
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
