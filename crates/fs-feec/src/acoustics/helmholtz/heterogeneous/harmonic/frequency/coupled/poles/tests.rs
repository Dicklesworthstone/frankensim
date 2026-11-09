//! G1/G3/G4: independent coupled cavity equations, hereditary and dark modes.
use super::*;
use super::super::{MechanicalHarmonicLimits, Wall};
use fs_sparse::{Coo, Csr};

fn scalar(value: f64) -> Csr {
    let mut matrix = Coo::new(1, 1);
    matrix.push(0, 0, value);
    matrix.assemble()
}
fn prepare(omega: f64, stiffness: f64, damping: f64, arms: Vec<(f64, f64)>) -> CoupledPencil {
    let fluid = [scalar(0.0), scalar(1.0), scalar(0.0), scalar(0.0)];
    CoupledPencil::prepare(omega, fluid.each_ref(), &[0], vec![Wall {
        dofs: vec![0], mass: scalar(1.0), area_m2: 1.0, surface_density: 1.0,
        stiffness, resistance: damping, arms,
    }], MechanicalHarmonicLimits::default(), &mut || false).unwrap()
}
fn options() -> MechanicalPoleOptions {
    MechanicalPoleOptions {
        spectral_scale: 3.0, step_tolerance: 2e-13, residual_tolerance: 2e-13,
        max_iterations: 40, max_dimension: 16, max_dense_bytes: 2 << 20,
        max_work: 1_000_000_000,
    }
}
fn limits() -> MechanicalPoleLimits {
    MechanicalPoleLimits { max_evaluations: 512, max_payload_bytes: 2 << 20,
        physical_residual_tolerance: 1e-11 }
}
fn close(actual: C64, expected: C64) {
    assert!((actual-expected).abs() < 2e-10*expected.abs().max(1.0), "{actual:?} != {expected:?}");
}

#[test]
fn g1_coupled_decay_uses_fluid_compliance_and_full_material_memory() {
    let s = C64::new(-0.5, 2.0);
    let rate = 3.0;
    let arm_stiffness: f64 = 2.0;
    // Independent cavity compression p=-w, wall force p=(s²+r*s+k)*w
    // plus the ORIGINAL Maxwell force k_a*s/(s+a)*w.
    let dynamic = s*s + (s/(s+C64::from_re(rate))).scale(arm_stiffness);
    let damping = -dynamic.im/s.im;
    let stiffness = -(dynamic+s.scale(damping)).re-1.0;
    assert!(damping > 0.0 && stiffness > 0.0);
    let coupling = arm_stiffness.sqrt();
    let extension = (s/(s+C64::from_re(rate))).scale(coupling);
    let seed = [C64::from_re(-1.0), C64::ONE, s, extension];
    for reference in [0.7, 1.7, 10.0] {
        let source = prepare(reference, stiffness, damping, vec![(coupling, rate)]);
        let pole = source.resonance_pole(C64::new(-0.45, 1.9), &seed, &options(), limits(), &mut || false).unwrap();
        close(pole.rate_per_s, s);
        let wall = &pole.walls[0];
        let w = wall.displacement_m[0];
        close(pole.pressure_pa[0], -w);
        close(wall.normal_velocity_m_s[0], s*w);
        close(wall.relaxation_extension[0][0], extension*w);
        close(wall.relaxation_memory[0][0], w.scale(coupling)-extension*w);
        assert!(pole.physical_residual < 1e-12);
    }
}

#[test]
fn g1_an_exact_material_relaxation_pole_retains_dark_original_states() {
    // Equal-time arms retain two original memories. Their antisymmetric
    // free decay is invisible in pressure and wall motion, but is a real mode.
    let source = prepare(0.7, 3.0, 0.2, vec![(2.0, 3.0), (2.0, 3.0)]);
    let pole = source.resonance_pole(C64::new(-2.9, 0.0),
        &[C64::ZERO, C64::ZERO, C64::ZERO, C64::ONE, C64::from_re(-1.0)],
        &options(), limits(), &mut || false).unwrap();
    close(pole.rate_per_s, C64::from_re(-3.0));
    assert!(pole.pressure_pa[0].abs() < 1e-13);
    close(pole.walls[0].relaxation_extension[0][0], -pole.walls[0].relaxation_extension[1][0]);
    close(pole.walls[0].relaxation_memory[0][0], -pole.walls[0].relaxation_extension[0][0]);
    assert!(pole.walls[0].relaxation_memory[0][0].abs() > 0.9);
}

#[test]
fn g1_contour_finds_the_loaded_mode_without_a_force_or_receiver() {
    let source = prepare(0.7, 3.0, 0.2, vec![]);
    let expected = C64::new(-0.1, (4.0_f64-0.01).sqrt());
    let region = MechanicalResonanceRegion { center: expected, real_radius: 0.3, imaginary_radius: 0.5 };
    let options = MechanicalSpectrumOptions {
        nodes: 128, probe_columns: 3, rank_relative_tolerance: 1e-9,
        count_tolerance: 1e-6, quadrature_tolerance: 1e-6,
        contour_min_reciprocal_condition: 1e-14, linear_residual_tolerance: 1e-12,
        svd_tolerance: 1e-13, svd_max_sweeps: 128, eigen_max_sweeps: 128,
        max_dimension: 16, max_dense_bytes: 2 << 20, max_work: 2_000_000_000,
        refinement: options(),
    };
    let result = source.resonance_spectrum(region, &options, limits(), &mut || false).unwrap();
    assert_eq!((result.argument_count, result.rank, result.poles.len()), (1,1,1));
    close(result.poles[0].rate_per_s, expected);
    assert!(result.poles[0].physical_residual < 1e-11);
    assert!(matches!(source.resonance_spectrum(MechanicalResonanceRegion {
        center: C64::ZERO, ..region }, &options, limits(), &mut || false), Err(MechanicalPoleError::Invalid(_))));
}

#[test]
fn g4_callback_budgets_and_cancellation_return_no_partial_mode() {
    let source = prepare(0.7, 3.0, 0.2, vec![]);
    let s = C64::new(-0.1, 2.0);
    let seed = [C64::from_re(-1.0), C64::ONE, s];
    assert!(matches!(source.resonance_pole(s, &seed, &options(),
        MechanicalPoleLimits { max_evaluations: 1, ..limits() }, &mut || false), Err(MechanicalPoleError::Budget { .. })));
    assert!(matches!(source.resonance_pole(s, &seed, &options(),
        MechanicalPoleLimits { max_payload_bytes: 1, ..limits() }, &mut || false), Err(MechanicalPoleError::Budget { .. })));
    assert!(matches!(source.resonance_pole(s, &seed, &options(), limits(), &mut || true), Err(MechanicalPoleError::Cancelled)));
    let mut calls=0;
    assert!(matches!(source.resonance_pole(s, &seed, &options(), limits(), &mut || {calls+=1; calls>6}), Err(MechanicalPoleError::Cancelled)));
    assert!(source.resonance_pole(C64::ZERO, &seed, &options(), limits(), &mut || false).is_err());
    // A cancelled read-only query leaves the original harmonic owner usable.
    assert!(source.solve(&[C64::ONE], &mut || false).unwrap().backward_error < 1e-12);
}
