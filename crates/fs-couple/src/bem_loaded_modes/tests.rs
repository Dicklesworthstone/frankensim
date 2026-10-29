use super::*;

pub(super) fn oscillator(cross_damping: f64) -> BemLoadedModes {
    // Unit-mass radial breathing plus an internal zero-opening coordinate.
    let vertices = [
        [1., 0., 0.],
        [-1., 0., 0.],
        [0., 1., 0.],
        [0., -1., 0.],
        [0., 0., 1.],
        [0., 0., -1.],
    ];
    let faces = [
        [0, 2, 4],
        [2, 1, 4],
        [1, 3, 4],
        [3, 0, 4],
        [2, 0, 5],
        [1, 2, 5],
        [3, 1, 5],
        [0, 3, 5],
    ];
    let surface = SpherePanels::from_triangles(
        faces
            .iter()
            .map(|f| f.map(|i| vertices[i].map(|v| 0.1 * v)))
            .collect(),
    )
    .unwrap();
    BemLoadedModes::build(
        surface,
        vec![vec![1.0; 8], vec![0.0; 8]],
        vec![300.0 * 300.0, 0., 0., 510.0 * 510.0],
        vec![0.8, cross_damping, cross_damping, 0.5],
        Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        BemModesOptions {
            formulation: Formulation::BurtonMiller,
            boundary: ComplexOptions {
                max_panels: 8,
                max_bytes: 64 * 64,
                max_work: 512 * 64,
            },
            panel_refinements: 0,
            reference_angular_frequency_rad_s: 400.,
            physical_residual_tolerance: 1e-10,
            max_dimension: 10,
            max_dense_bytes: 100_000,
            max_operator_evaluations: 200,
            max_projection_work: 100_000,
        },
        &mut || false,
    )
    .unwrap()
}
pub(super) fn options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 400.,
        step_tolerance: 1e-12,
        residual_tolerance: 1e-12,
        max_iterations: 20,
        max_dimension: 10,
        max_dense_bytes: 100_000,
        max_work: 1_000_000,
    }
}

#[test]
fn g1_complete_boundary_pencil_retains_hidden_mode_and_full_material_coupling() {
    let model = oscillator(0.2);
    let opts = BeynOptions {
        nodes: 128,
        probe_columns: 3,
        rank_relative_tolerance: 1e-8,
        count_tolerance: 1e-6,
        quadrature_tolerance: 1e-5,
        contour_min_reciprocal_condition: 1e-12,
        linear_residual_tolerance: 1e-11,
        svd_tolerance: 1e-13,
        svd_max_sweeps: 64,
        eigen_max_sweeps: 128,
        max_dimension: 10,
        max_dense_bytes: 1_000_000,
        max_work: 100_000_000,
        refinement: options(),
    };
    let spectrum = model
        .resonance_spectrum(
            Ellipse {
                center: C64::new(405., -2.),
                real_radius: 140.,
                imaginary_radius: 60.,
            },
            &opts,
            || false,
        )
        .unwrap();
    assert_eq!(spectrum.argument_count, 2);
    assert_eq!(spectrum.poles.len(), 2);
    for pole in &spectrum.poles {
        assert!(pole.angular_frequency_rad_s.im < 0.);
        assert!(pole.mechanical_residual < 1e-10 && pole.boundary_residual < 1e-10);
        assert_eq!(pole.modal_displacement.len(), 2);
        assert_eq!(pole.surface_pressure_pa.len(), 8);
    }
    let seed = [C64::ZERO, C64::ONE];
    let hidden = model
        .resonance_pole(C64::new(510., -0.25), &seed, &options(), || false)
        .unwrap();
    assert!(hidden.modal_displacement[0].abs() > 1e-5);
    assert!(hidden.surface_pressure_pa.iter().any(|p| p.abs() > 1e-5));
    let decoupled = oscillator(0.0)
        .resonance_pole(C64::new(510., -0.25), &seed, &options(), || false)
        .unwrap();
    assert!((hidden.angular_frequency_rad_s - decoupled.angular_frequency_rad_s).abs() > 1e-6);
    let reversed = oscillator(-0.2)
        .resonance_pole(C64::new(510., -0.25), &seed, &options(), || false)
        .unwrap();
    assert!((reversed.angular_frequency_rad_s - hidden.angular_frequency_rad_s).abs() < 1e-9);
}

#[test]
fn g1_full_coupled_analytic_derivative_matches_both_complex_directions() {
    let model = oscillator(0.2);
    let omega = C64::new(380., -4.);
    let value = model.original_pencil(omega, || false).unwrap();
    for direction in [C64::ONE, C64::new(0., 1.)] {
        let h = direction.scale(0.002);
        let plus = model.original_pencil(omega + h, || false).unwrap();
        let minus = model.original_pencil(omega - h, || false).unwrap();
        for (j, &exact) in value.derivative.iter().enumerate() {
            let difference = (plus.matrix[j] - minus.matrix[j]) / (h.scale(2.));
            assert!(
                (difference - exact).abs() < 2e-7 * (1. + exact.abs()),
                "entry {j}: {difference:?} versus {exact:?}"
            );
        }
    }
}
