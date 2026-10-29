//! Physical checks of source-owned, two-sided BEM-loaded shell resonances.
//!
//! The retained structural basis and the panel-constant boundary velocity field
//! are declared approximations. Acoustic panel refinement holds those source
//! fields fixed; the modal ladder holds the structural mesh and acoustic rule
//! fixed. No fitted resistance, acoustic pole or prescribed quality factor is
//! used. These low-frequency contours avoid the closed boundary's interior
//! acoustic resonances; they do not qualify arbitrary complex BIE contours.
#![cfg(not(target_arch = "wasm32"))]

use fs_bem::helmholtz::laplace::complex::ComplexOptions;
use fs_bem::helmholtz::{Formulation, Medium, solve_radiation};
use fs_couple::bem_loaded_modes::{BemLoadedModes, BemModePole, BemModesError, BemModesOptions};
use fs_la::nonlinear_eig::beyn::{BeynOptions, Ellipse};
use fs_la::nonlinear_eig::{NonlinearEigError, PoleOptions, newton_pole_with_cancel};
use fs_math::c64::C64;
use fs_plate::shell::reduction::radiation::RadiationSurfaceBudget;

#[path = "support/material_packs.rs"]
mod material_packs;
#[path = "bem_loaded_plate_modes/source.rs"]
mod source;

fn settings(omega: f64, refinements: usize) -> BemModesOptions {
    BemModesOptions {
        formulation: Formulation::PlainCbie,
        boundary: ComplexOptions {
            max_panels: 512,
            max_bytes: 32 * 1024 * 1024,
            max_work: 512 * 512 * 512,
        },
        panel_refinements: refinements,
        reference_angular_frequency_rad_s: omega,
        physical_residual_tolerance: 1e-10,
        max_dimension: 512,
        max_dense_bytes: 64 * 1024 * 1024,
        max_operator_evaluations: 128,
        max_projection_work: 100_000_000,
    }
}

fn pole_options(omega: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: omega,
        step_tolerance: 1e-12,
        residual_tolerance: 1e-12,
        max_iterations: 16,
        max_dimension: 512,
        max_dense_bytes: 64 * 1024 * 1024,
        max_work: 4_000_000_000,
    }
}

fn spectrum_options(omega: f64) -> BeynOptions {
    BeynOptions {
        nodes: 64,
        probe_columns: 2,
        rank_relative_tolerance: 1e-8,
        count_tolerance: 1e-6,
        quadrature_tolerance: 1e-5,
        contour_min_reciprocal_condition: 1e-12,
        linear_residual_tolerance: 1e-11,
        svd_tolerance: 1e-13,
        svd_max_sweeps: 64,
        eigen_max_sweeps: 128,
        max_dimension: 512,
        max_dense_bytes: 64 * 1024 * 1024,
        max_work: 20_000_000_000,
        refinement: pole_options(omega),
    }
}

fn surface_budget() -> RadiationSurfaceBudget {
    RadiationSurfaceBudget {
        max_panels: 512,
        max_panel_modes: 4096,
    }
}

fn loaded(plate: &source::Plate, medium: Medium, refinements: usize) -> BemLoadedModes {
    BemLoadedModes::from_shell(
        &plate.reduction,
        &plate.thickness,
        surface_budget(),
        medium,
        None,
        settings(plate.reduction.omegas()[0], refinements),
        || false,
    )
    .expect("same source shell, mass and complete exterior boundary")
}

fn first_pole(model: &BemLoadedModes, omega: f64) -> BemModePole {
    let mut seed = vec![C64::ZERO; model.mode_count()];
    seed[0] = C64::ONE;
    model
        .resonance_pole(
            C64::new(omega, -1e-4 * omega),
            &seed,
            &pole_options(omega),
            || false,
        )
        .expect("isolated first structural resonance")
}

fn assert_physical(pole: &BemModePole) {
    assert!(pole.angular_frequency_rad_s.re > 0.0);
    assert!(pole.angular_frequency_rad_s.im < 0.0, "{pole:?}");
    assert!(pole.mechanical_residual <= 1e-10, "{pole:?}");
    assert!(pole.boundary_residual <= 1e-10, "{pole:?}");
    assert!(pole.quality_factor().is_some_and(|q| q > 0.0));
    assert!(pole.surface_pressure_pa.iter().any(|p| p.abs() > 0.0));
}

#[test]
fn g2_source_owned_plate_pole_and_contour_agree_in_original_equations() {
    let plate = source::plate(3, 3, 3);
    let omega = plate.reduction.omegas()[0];
    let model = loaded(&plate, source::air(1.0), 0);
    let seeded = first_pole(&model, omega);
    assert_physical(&seeded);
    let contour = Ellipse {
        center: C64::new(omega, -0.01 * omega),
        real_radius: 0.08 * omega,
        imaginary_radius: 0.06 * omega,
    };
    let spectrum = model
        .resonance_spectrum(contour, &spectrum_options(omega), || false)
        .expect("original boundary/structure contour");
    assert_eq!(spectrum.argument_count, 1);
    assert_eq!(spectrum.rank, 1);
    assert_eq!(spectrum.poles.len(), 1);
    let discovered = &spectrum.poles[0];
    assert_physical(discovered);
    assert!(
        (discovered.angular_frequency_rad_s - seeded.angular_frequency_rad_s).abs() / omega < 1e-10
    );
    assert_eq!(seeded.modal_displacement.len(), 3);
    assert_eq!(seeded.surface_pressure_pa.len(), model.panel_count());
}

#[test]
fn g1_light_fluid_shift_has_second_order_remainder() {
    let plate = source::plate(3, 3, 3);
    let omega = plate.reduction.omegas()[0];
    let mut errors = Vec::new();
    for fraction in [0.25, 0.5, 1.0] {
        let model = loaded(&plate, source::air(fraction), 0);
        let pole = first_pole(&model, omega);
        assert_physical(&pole);
        // Independent original REAL-frequency BEM: one physical unit modal
        // velocity, with the same area-weighted pressure-work transpose.
        let velocity: Vec<_> = model.normal_velocity_weights()[0]
            .iter()
            .map(|&w| C64::from_re(w))
            .collect();
        let real = solve_radiation(
            model.surface(),
            omega / model.medium().sound_speed,
            model.medium(),
            &velocity,
            Formulation::PlainCbie,
        )
        .expect("independent real-frequency radiation response");
        let z = real
            .pressure
            .iter()
            .zip(&model.normal_velocity_weights()[0])
            .zip(model.surface().areas())
            .fold(C64::ZERO, |sum, ((&p, &w), &area)| sum + p.scale(w * area));
        assert!(
            z.re > 0.0,
            "passive real-frequency radiation must remove energy"
        );
        // K-omega^2 I-i omega Z(omega) gives delta omega = -i Z_rr/2
        // for an isolated unit-mass vacuum mode, to first order in density.
        let first_order = C64::new(z.im * 0.5, -z.re * 0.5);
        let error = (pole.angular_frequency_rad_s - C64::from_re(omega) - first_order).abs();
        errors.push(error);
        println!(
            "density_fraction={fraction} omega={:?} first_order={first_order:?} remainder={error:.12e}",
            pole.angular_frequency_rad_s
        );
    }
    for pair in errors.windows(2) {
        let order = (pair[1] / pair[0]).log2();
        assert!(
            order > 1.8 && order < 2.2,
            "density remainder order={order}, errors={errors:?}"
        );
    }
}

#[test]
fn g1_reversed_pressure_work_is_detected_as_radiation_growth() {
    let plate = source::plate(3, 3, 3);
    let omega = plate.reduction.omegas()[0];
    let model = loaded(&plate, source::air(1.0), 0);
    let good = first_pole(&model, omega);
    assert_physical(&good);
    let m = model.mode_count();
    let n = model.dimension();
    let mut seed = good.modal_displacement.clone();
    seed.extend_from_slice(&good.surface_pressure_pa);
    let mutant = newton_pole_with_cancel(
        n,
        good.angular_frequency_rad_s.conj(),
        &seed,
        &pole_options(omega),
        |w, cancel| {
            let mut value = model
                .original_pencil(w, cancel)
                .map_err(|_| NonlinearEigError::FunctionRefused { at: w })?;
            // Mutation ONLY: pressure pushes outward instead of inward. The
            // acoustic field equation and material source remain identical.
            for i in 0..m {
                for j in m..n {
                    value.matrix[i * n + j] = -value.matrix[i * n + j];
                    value.derivative[i * n + j] = -value.derivative[i * n + j];
                }
            }
            Ok(value)
        },
        || false,
    )
    .expect("sign mutation must resolve its own growing root");
    assert!(
        mutant.value.im > 0.0,
        "reversed pressure work should add energy: {mutant:?}"
    );
    assert!(mutant.residual < 1e-10);
}

#[test]
fn g1_panel_and_source_modal_axes_are_independent_and_live() {
    let plate = source::plate(3, 3, 3);
    let omega = plate.reduction.omegas()[0];
    let mut panel_values = Vec::new();
    for refinement in 0..3 {
        let model = loaded(&plate, source::air(1.0), refinement);
        let pole = first_pole(&model, omega);
        assert_physical(&pole);
        assert_eq!(model.panel_count(), 60 * (1 << refinement));
        println!(
            "axis=acoustic_panels panels={} retained_modes=3 omega={:?}",
            model.panel_count(),
            pole.angular_frequency_rad_s
        );
        panel_values.push(pole.angular_frequency_rad_s);
    }
    let panel_steps: Vec<_> = panel_values
        .windows(2)
        .map(|p| (p[1] - p[0]).abs())
        .collect();
    assert!(panel_steps[0] > 1e-4 && panel_steps[1] > 1e-4);
    assert!(panel_steps[1] < 0.8 * panel_steps[0], "{panel_values:?}");

    let mut modal_values = Vec::new();
    for retained in 1..=3 {
        let source = source::plate(3, 3, retained);
        let model = loaded(&source, source::air(1.0), 0);
        let pole = first_pole(&model, source.reduction.omegas()[0]);
        assert_physical(&pole);
        assert_eq!(model.panel_count(), 60);
        println!(
            "axis=source_modes panels=60 retained_modes={retained} omega={:?}",
            pole.angular_frequency_rad_s
        );
        modal_values.push(pole.angular_frequency_rad_s);
    }
    let modal_steps: Vec<_> = modal_values
        .windows(2)
        .map(|p| (p[1] - p[0]).abs())
        .collect();
    assert!(modal_steps[0] > 1e-8 && modal_steps[1] > 1e-8);
    assert!(modal_steps[1] < 0.5 * modal_steps[0], "{modal_values:?}");
    // These observed ladders expose both axes. Three rungs do not establish
    // continuum convergence of the structural mesh or the source velocity map.
}

#[test]
fn g0_original_dimension_cap_and_cancellation_leave_model_reusable() {
    let plate = source::plate(3, 3, 3);
    let omega = plate.reduction.omegas()[0];
    let mut limits = settings(omega, 0);
    limits.max_dimension = 63;
    let make = |limits| {
        BemLoadedModes::from_shell(
            &plate.reduction,
            &plate.thickness,
            surface_budget(),
            source::air(1.0),
            None,
            limits,
            || false,
        )
    };
    let model = make(limits).expect("complete original dimension at its exact cap");
    limits.max_dimension -= 1;
    assert!(
        make(limits).is_err(),
        "the last pressure unknown may not be discarded"
    );
    let first = first_pole(&model, omega);
    let seed = [C64::ONE, C64::ZERO, C64::ZERO];
    assert!(matches!(
        model.resonance_pole(C64::from_re(omega), &seed, &pole_options(omega), || true),
        Err(BemModesError::Cancelled)
    ));
    let repeated = first_pole(&model, omega);
    let bits = |z: C64| (z.re.to_bits(), z.im.to_bits());
    assert_eq!(
        bits(first.angular_frequency_rad_s),
        bits(repeated.angular_frequency_rad_s)
    );
    for (&a, &b) in first
        .modal_displacement
        .iter()
        .chain(&first.surface_pressure_pa)
        .zip(
            repeated
                .modal_displacement
                .iter()
                .chain(&repeated.surface_pressure_pa),
        )
    {
        assert_eq!(bits(a), bits(b));
    }
    let mut exhausted = settings(omega, 0);
    exhausted.max_operator_evaluations = 1;
    let limited = make(exhausted).expect("one construction assembly is admitted");
    assert!(matches!(
        limited.resonance_pole(C64::from_re(omega), &seed, &pole_options(omega), || false),
        Err(BemModesError::EvaluationBudget { maximum: 1 })
    ));
}
