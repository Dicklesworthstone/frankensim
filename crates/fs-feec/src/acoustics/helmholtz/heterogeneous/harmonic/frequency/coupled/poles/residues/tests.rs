//! Independent scalar and hereditary cavity equations for absolute residues.
use super::*;
use super::super::super::{MechanicalHarmonicLimits, Wall};
use super::super::MechanicalPoleOptions;
use fs_sparse::{Coo, Csr};

fn scalar(value: f64) -> Csr {
    let mut a = Coo::new(1, 1);
    a.push(0, 0, value);
    a.assemble()
}
fn pole_options() -> MechanicalPoleOptions {
    MechanicalPoleOptions {
        spectral_scale: 3.0, step_tolerance: 2e-13, residual_tolerance: 2e-13,
        max_iterations: 40, max_dimension: 16, max_dense_bytes: 2 << 20,
        max_work: 1_000_000_000,
    }
}
fn response_options() -> MechanicalPoleResponseOptions {
    MechanicalPoleResponseOptions {
        spectral_scale: 3.0, pole_residual_tolerance: 2e-12,
        pole_step_tolerance: 2e-12, linear_residual_tolerance: 2e-12,
        max_dimension: 16, max_vectors: 1, max_dense_bytes: 2 << 20,
        max_work: 1_000_000_000,
    }
}
fn limits() -> MechanicalPoleLimits {
    MechanicalPoleLimits {
        max_evaluations: 128, max_payload_bytes: 2 << 20,
        physical_residual_tolerance: 1e-11,
    }
}
fn close(actual: C64, expected: C64) {
    assert!((actual-expected).abs() < 3e-10*expected.abs().max(1.0),
        "{actual:?} != {expected:?}");
}
fn cavity(omega: f64, k: f64, c: f64, arms: Vec<(f64, f64)>) -> CoupledPencil {
    let fluid = [scalar(0.0), scalar(1.0), scalar(0.0), scalar(0.0)];
    CoupledPencil::prepare(omega, fluid.each_ref(), &[0], vec![Wall {
        dofs: vec![0], mass: scalar(1.0), surface_density: 1.0,
        stiffness: k, resistance: c, arms,
    }], MechanicalHarmonicLimits::default(), &mut || false).unwrap()
}

#[test]
fn scalar_velocity_source_residue_retains_absolute_units_and_complex_phase() {
    let s = C64::new(-0.2, (4.0_f64-0.04).sqrt());
    let b = C64::new(0.7, -0.3);
    for reference in [0.7, 7.0] {
        let fluid = [scalar(4.0), scalar(1.0), scalar(0.4), scalar(0.0)];
        let owner = CoupledPencil::prepare(reference, fluid.each_ref(), &[0], vec![],
            MechanicalHarmonicLimits::default(), &mut || false).unwrap();
        let pole = owner.resonance_pole(s, &[C64::ONE], &pole_options(),
            limits(), &mut || false).unwrap();
        let prepared = owner.prepare_resonance_response(&pole, response_options(),
            limits(), &mut || false).unwrap();
        let actual = prepared.forcing_residue(&[b], &[], &mut || false).unwrap();
        let expected = s*b/(C64::from_re(0.4)+s.scale(2.0));
        close(actual.pressure_pa_per_s[0], expected);
        close(actual.coordinates[0], expected);
        assert!(actual.walls.is_empty());
        assert!(actual.bordered_relative_residual < 2e-12);
        // Independent Laurent limit of the COMPLETE scalar transfer s*b/T(s).
        let delta = C64::new(1e-5, -2e-5);
        let at = s+delta;
        let limit = delta*at*b/(C64::from_re(4.0)+at.scale(0.4)+at*at);
        assert!((limit-expected).abs() < 2e-5);
        assert!(prepared.forcing_residue(&[b], &[], &mut || true).is_err());
        close(prepared.forcing_residue(&[b], &[], &mut || false).unwrap()
            .pressure_pa_per_s[0], expected);
        let mut wrong = pole.clone();
        wrong.rate_per_s = C64::new(-0.7, 2.7);
        assert!(owner.prepare_resonance_response(&wrong, response_options(),
            limits(), &mut || false).is_err());
        let mut capped = limits();
        capped.max_evaluations = 1;
        assert!(matches!(owner.prepare_resonance_response(&pole, response_options(),
            capped, &mut || false), Err(MechanicalPoleError::Budget { .. })));
    }
}

#[test]
fn coupled_maxwell_residue_keeps_fluid_traction_interference_and_every_memory() {
    let s = C64::new(-0.5, 2.0);
    let rate = 3.0;
    let ka = 2.0_f64;
    let dynamic = s*s+(s/(s+C64::from_re(rate))).scale(ka);
    let c = -dynamic.im/s.im;
    let k = -(dynamic+s.scale(c)).re-1.0;
    let coupling = ka.sqrt();
    let extension = (s/(s+C64::from_re(rate))).scale(coupling);
    let b = C64::new(0.4, -0.7);
    let t = C64::new(-0.2, 0.3);
    for reference in [0.7, 10.0] {
        let owner = cavity(reference, k, c, vec![(coupling, rate)]);
        let seed = [C64::from_re(-1.0), C64::ONE, s, extension];
        let pole = owner.resonance_pole(s, &seed, &pole_options(),
            limits(), &mut || false).unwrap();
        let prepared = owner.prepare_resonance_response(&pole, response_options(),
            limits(), &mut || false).unwrap();
        let traction = owner.prepare_traction(0, &[1.0], &mut || false).unwrap();
        let value = prepared.forcing_residue(&[b], &[(&traction, t)],
            &mut || false).unwrap();
        // Eliminate only in this independent scalar oracle:
        // L(s)w=b/s+t; L=1+k+c*s+s^2+ka*s/(s+rate).
        let d = s+C64::from_re(rate);
        let derivative = s.scale(2.0)+C64::from_re(c)+C64::from_re(ka*rate)/d/d;
        let w = (b/s+t)/derivative;
        close(value.pressure_pa_per_s[0], -w);
        let wall = &value.walls[0];
        close(wall.displacement_m_per_s[0], w);
        close(wall.normal_velocity_m_s2[0], s*w);
        close(wall.relaxation_extension_per_s[0][0], extension*w);
        close(wall.relaxation_memory_per_s[0][0], w.scale(coupling)-extension*w);
        // Coherent forcing can cancel participation while each source is nonzero.
        let cancelled = prepared.forcing_residue(&[b], &[(&traction, -b/pole.rate_per_s)],
            &mut || false).unwrap();
        assert!(cancelled.coordinates.iter().all(|z| z.abs() < 2e-12));
    }
}

#[test]
fn an_exact_dark_relaxation_pole_has_zero_acoustic_and_traction_participation() {
    let owner = cavity(0.7, 3.0, 0.2, vec![(2.0, 3.0), (2.0, 3.0)]);
    let s = C64::from_re(-3.0);
    let seed = [C64::ZERO, C64::ZERO, C64::ZERO, C64::ONE, C64::from_re(-1.0)];
    let pole = owner.resonance_pole(s, &seed, &pole_options(),
        limits(), &mut || false).unwrap();
    let prepared = owner.prepare_resonance_response(&pole, response_options(),
        limits(), &mut || false).unwrap();
    let traction = owner.prepare_traction(0, &[1.0], &mut || false).unwrap();
    let residue = prepared.forcing_residue(&[C64::new(0.3, 0.7)],
        &[(&traction, C64::new(-0.2, 0.4))], &mut || false).unwrap();
    assert!(prepared.mode().walls[0].relaxation_extension[0][0].abs() > 0.9);
    assert!(residue.coordinates.iter().all(|z| z.abs() < 2e-12));
}
