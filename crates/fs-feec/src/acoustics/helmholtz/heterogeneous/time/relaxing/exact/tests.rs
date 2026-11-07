//! Original-state evolution, independent moments/ODEs and physical trajectory work.

use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use crate::acoustics::helmholtz::heterogeneous::time::relaxing::RelaxingTimeLimits;
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

const FACE: [u32; 3] = [0, 1, 2];
const RELEASE: [u32; 3] = [0, 1, 3];
const POSITIONS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
];
const MEDIUM: Medium = Medium { density_kg_m3: 2.0, sound_speed_m_s: 3.0 };

fn material<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(geometry, 4, &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM), &mut || false).unwrap()
}

fn field<'m, 'f, 'c>(model: &'m HeterogeneousHelmholtzFe<'f, 'c>, law: impl Fn([f64; 3]) -> f64) -> SpatialResistiveBoundary<'m, 'f, 'c> {
    model.spatial_resistive_boundary(&[FACE], 4,
        &mut |_, sample| Ok::<_, AcousticFeError>(law(sample.point.position_m)),
        &mut || false).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual-expected).abs()<=tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}");
}

#[test]
fn exact_trajectory_matches_independent_original_moments_and_physical_integrals() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let c0 = field(&model, |_| 0.1);
    let residue = field(&model, |x| 1.0+x[0]);
    let law = model.spatial_relaxing_boundary(Some(&c0), &[(&residue, 0.04)], &mut || false).unwrap();
    let ports = [ExactVelocityPort { profile: vec![(FACE, [1.0; 3])] }];
    let duration = 0.025;
    let owner = model.prepare_exact_time_with_relaxing_boundaries(
        duration, &[], &[], &[&law], &ports, &[], &[], &[RELEASE],
        ExactRelaxingTimeLimits::default(), &mut || false).unwrap();
    // One free pressure vertex. Exact unit-tetra/triangle moments give
    // M22=1/1080, K22=1/12, C0,22=1/120, Cj,22=1/10, b2=0.06/6.
    // Integrate this independent physical ODE and its four separate accounts.
    let derivative = |x: [f64; 7]| {
        let [psi, p, q, _, _, _, _] = x;
        [p, 1080.0*(0.01-psi/12.0-p/120.0-q/10.0), (p-q)/0.04,
            p*p/120.0, q*q/10.0, p*q/10.0, 0.01*p]
    };
    let mut reference = [0.0, 1.0, 0.2, 0.0, 0.0, 0.0, 0.0];
    let dt = duration/4096.0;
    for _ in 0..4096 {
        let a = derivative(reference);
        let b = derivative(core::array::from_fn(|i| reference[i]+0.5*dt*a[i]));
        let c = derivative(core::array::from_fn(|i| reference[i]+0.5*dt*b[i]));
        let d = derivative(core::array::from_fn(|i| reference[i]+dt*c[i]));
        reference = core::array::from_fn(|i| reference[i]+dt*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0);
    }
    let initial = owner.state_with_memories(&[0.0; 4], &[0.0, 0.0, 1.0, 0.0], &[&[0.0, 0.0, 0.2]]).unwrap();
    let result = owner.step(&initial, &[0.06], &mut || false).unwrap();
    close(result.state.pressure_potential_pa_s()[2], reference[0], 2e-11);
    close(result.state.pressure_pa()[2], reference[1], 2e-10);
    close(result.state.memory_pressure_pa(0).unwrap()[2], reference[2], 2e-10);
    for (actual, expected) in [
        (result.work.instantaneous_dissipation_j, reference[3]),
        (result.work.relaxation_dissipation_j, reference[4]),
        (result.work.memory_work_j, reference[5]),
        (result.work.input_work_j, reference[6]),
    ] { close(actual, expected, 2e-12); }
    close(result.integrated_port_pressure_pa_m2_s[0]*0.06, reference[6], 2e-12);
    close(result.work.balance_residual_j(), 0.0, 2e-12);
    close(result.work.acoustic_balance_residual_j(), 0.0, 2e-12);
    close(result.work.memory_balance_residual_j(), 0.0, 2e-12);
    // A second cache with half the clock must follow the same held trajectory.
    let half = model.prepare_exact_time_with_relaxing_boundaries(
        duration/2.0, &[], &[], &[&law], &ports, &[], &[], &[RELEASE],
        ExactRelaxingTimeLimits::default(), &mut || false).unwrap();
    let s0 = half.state_with_memories(initial.pressure_potential_pa_s(), initial.pressure_pa(), &[initial.memory_pressure_pa(0).unwrap()]).unwrap();
    let s1 = half.step(&s0, &[0.06], &mut || false).unwrap();
    let s2 = half.step(&s1.state, &[0.06], &mut || false).unwrap();
    for (partitioned, whole) in [
        (s2.state.pressure_potential_pa_s(), result.state.pressure_potential_pa_s()),
        (s2.state.pressure_pa(), result.state.pressure_pa()),
        (s2.state.memory_pressure_pa(0).unwrap(), result.state.memory_pressure_pa(0).unwrap()),
    ] {
        for (&actual, &expected) in partitioned.iter().zip(whole) {
            close(actual, expected, 3e-11);
        }
    }
    for (partitioned, whole) in [
        (s1.work.instantaneous_dissipation_j+s2.work.instantaneous_dissipation_j, result.work.instantaneous_dissipation_j),
        (s1.work.relaxation_dissipation_j+s2.work.relaxation_dissipation_j, result.work.relaxation_dissipation_j),
        (s1.work.memory_work_j+s2.work.memory_work_j, result.work.memory_work_j),
        (s1.work.input_work_j+s2.work.input_work_j, result.work.input_work_j),
        (s1.integrated_port_pressure_pa_m2_s[0]+s2.integrated_port_pressure_pa_m2_s[0], result.integrated_port_pressure_pa_m2_s[0]),
    ] {
        close(partitioned, whole, 2e-12);
    }
}

#[test]
fn complete_spatial_midpoint_fields_refine_to_the_exact_causal_trajectory() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let c0 = field(&model, |_| 0.08);
    let cx = field(&model, |x| 0.6*(1.0+x[0]));
    let cy = field(&model, |x| 0.4*(1.0+x[1]));
    let law = model.spatial_relaxing_boundary(Some(&c0), &[(&cx, 0.017), (&cy, 0.09)], &mut || false).unwrap();
    let duration = 0.04;
    let exact = model.prepare_exact_time_with_relaxing_boundaries(
        duration, &[], &[], &[&law], &[], &[], &[], &[],
        ExactRelaxingTimeLimits::default(), &mut || false).unwrap();
    let pressure = [1.0, -0.4, 2.0, 0.3];
    let memories: [&[f64]; 2] = [&[0.2, -0.1, 0.4], &[-0.3, 0.2, 0.1]];
    let initial = exact.state_with_memories(&[0.0; 4], &pressure, &memories).unwrap();
    let expected = exact.step(&initial, &[], &mut || false).unwrap();
    assert!(expected.work.instantaneous_dissipation_j>0.0);
    assert!(expected.work.relaxation_dissipation_j>0.0);
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let midpoint = model.prepare_time_with_relaxing_boundaries(
            duration/steps as f64, &[], &[], &[&law], &[], RelaxingTimeLimits::default(), &mut || false).unwrap();
        let mut state = midpoint.state_with_memories(&[0.0; 4], &pressure, &memories).unwrap();
        let mut heat = 0.0;
        for _ in 0..steps {
            let result = midpoint.step(&state, &[], &mut || false).unwrap();
            heat += result.work.relaxation_dissipation_j;
            state = result.state;
        }
        let mut error = state.pressure_pa().iter().zip(expected.state.pressure_pa())
            .map(|(a,b)| (a-b).abs()).fold(0.0, f64::max);
        for arm in 0..2 {
            error = state.memory_pressure_pa(arm).unwrap().iter()
                .zip(expected.state.memory_pressure_pa(arm).unwrap())
                .map(|(a,b)| (a-b).abs()).fold(error, f64::max);
        }
        errors.push(error);
        close(heat, expected.work.relaxation_dissipation_j,
            3.0*exact.energy(&initial).unwrap().total_j()/(steps*steps) as f64);
    }
    assert!(errors[0]>1e-5, "nonvacuous full-state refinement: {errors:?}");
    for pair in errors.windows(2) { assert!(pair[1]<0.28*pair[0], "second-order full-state refinement: {errors:?}"); }
}

#[test]
fn zero_residues_keep_all_release_shared_memories_and_decay_on_the_exact_clock() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let zero = field(&model, |_| 0.0);
    let tau = 0.0002;
    let duration = 0.01;
    let law = model.spatial_relaxing_boundary(None, &[(&zero, tau)], &mut || false).unwrap();
    let owner = model.prepare_exact_time_with_relaxing_boundaries(
        duration, &[], &[], &[&law], &[], &[], &[], &[RELEASE],
        ExactRelaxingTimeLimits::default(), &mut || false).unwrap();
    assert_eq!(owner.memory_dofs(0), Some([0, 1, 2].as_slice()));
    assert_eq!(owner.memory_patch_and_term(0), Some((0, 0)));
    assert_eq!(owner.memory_relaxation_time_s(0), Some(tau));
    let initial = owner.state_with_memories(&[0.0; 4], &[0.0; 4], &[&[1.0, -2.0, 3.0]]).unwrap();
    let result = owner.step(&initial, &[], &mut || false).unwrap();
    let decay = det::exp(-duration/tau);
    for (&actual, authored) in result.state.memory_pressure_pa(0).unwrap().iter().zip([1.0, -2.0, 3.0]) {
        close(actual, authored*decay, 2e-10*decay);
    }
    for &value in result.state.pressure_pa() { close(value, 0.0, 1e-28); }
    assert_eq!(owner.energy(&result.state).unwrap().total_j(), 0.0);
    assert_eq!(result.work.relaxation_dissipation_j, 0.0);
    assert_eq!(result.work.memory_work_j, 0.0);
}

#[test]
fn original_held_sources_share_one_trajectory_and_retain_signed_work_without_resampling() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let calls = Cell::new(0);
    let residue = field(&model, |x| { calls.set(calls.get()+1); 1.0+x[0] });
    let law = model.spatial_relaxing_boundary(None, &[(&residue, 0.05)], &mut || false).unwrap();
    let ports = [ExactVelocityPort { profile: vec![(FACE, [1.0; 3])] }];
    let boundary = model.quadratic_velocity_source(&[(FACE, [1.0; 6])], &mut || false).unwrap();
    let volume = model.volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || false).unwrap();
    let owner = model.prepare_exact_time_with_relaxing_boundaries(
        0.003, &[], &[], &[&law], &ports, &[&boundary], &[&volume], &[],
        ExactRelaxingTimeLimits::default(), &mut || false).unwrap();
    let initial = owner.state(&[0.0; 4], &[0.0; 4]).unwrap();
    let result = owner.step_with_volume_sources(&initial, &[0.2, -0.05], &[0.03], &mut || false).unwrap();
    let face = &result.integrated_port_pressure_pa_m2_s;
    let volume = &result.integrated_volume_source_pressure_pa_s;
    assert_eq!(face.len(), 2); assert_eq!(volume.len(), 1);
    close(face[0], face[1], 2e-12);
    assert!(result.work.boundary_source_work_j<0.0);
    assert!(result.work.volume_source_work_j>0.0);
    close(result.work.boundary_source_work_j, -0.05*face[1], 2e-13);
    close(result.work.volume_source_work_j, 0.03*volume[0], 2e-13);
    close(result.work.input_work_j, 0.2*face[0]-0.05*face[1]+0.03*volume[0], 2e-13);
    close(result.work.balance_residual_j(), 0.0, 2e-12);
    owner.step_with_volume_sources(&result.state, &[-0.1, 0.02], &[-0.01], &mut || false).unwrap();
    assert_eq!(calls.get(), 16);
}

#[test]
fn exact_memory_admission_limits_ownership_and_cancellation_preserve_original_inputs() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let residue = field(&model, |x| 1.0+x[0]);
    let law = model.spatial_relaxing_boundary(None, &[(&residue, 0.1)], &mut || false).unwrap();
    for duration in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        assert!(model.prepare_exact_time_with_relaxing_boundaries(
            duration, &[], &[], &[&law], &[], &[], &[], &[],
            ExactRelaxingTimeLimits::default(), &mut || false).is_err());
    }
    assert!(matches!(model.prepare_exact_time_with_relaxing_boundaries(
        0.01, &[], &[], &[&law], &[], &[], &[], &[],
        ExactRelaxingTimeLimits { max_states: 10, ..ExactRelaxingTimeLimits::default() },
        &mut || false), Err(HeterogeneousError::Limit { .. })));
    assert!(matches!(model.prepare_exact_time_with_relaxing_boundaries(
        0.01, &[], &[], &[&law], &[], &[], &[], &[],
        ExactRelaxingTimeLimits::default(), &mut || true), Err(HeterogeneousError::Cancelled)));
    let other_model = material(&geometry);
    assert!(other_model.prepare_exact_time_with_relaxing_boundaries(
        0.01, &[], &[], &[&law], &[], &[], &[], &[],
        ExactRelaxingTimeLimits::default(), &mut || false).is_err());
    assert!(model.prepare_exact_time_with_relaxing_boundaries(
        0.01, &[], &[], &[&law], &[], &[], &[], &[FACE],
        ExactRelaxingTimeLimits::default(), &mut || false).is_err());
    let owner = model.prepare_exact_time_with_relaxing_boundaries(
        0.01, &[], &[], &[&law], &[], &[], &[], &[],
        ExactRelaxingTimeLimits::default(), &mut || false).unwrap();
    let other = model.prepare_exact_time_with_relaxing_boundaries(
        0.01, &[], &[], &[&law], &[], &[], &[], &[],
        ExactRelaxingTimeLimits::default(), &mut || false).unwrap();
    let initial = owner.state_with_memories(&[0.0; 4], &[1.0, 2.0, 3.0, 4.0], &[&[0.2, -0.1, 0.3]]).unwrap();
    assert!(other.step(&initial, &[], &mut || false).is_err());
    assert!(owner.step(&initial, &[1.0], &mut || false).is_err());
    let mut polls = 0;
    assert!(matches!(owner.step(&initial, &[], &mut || { polls+=1; polls==3 }), Err(HeterogeneousError::Cancelled)));
    assert_eq!(initial.pressure_pa(), &[1.0, 2.0, 3.0, 4.0]);
    assert_eq!(initial.memory_pressure_pa(0).unwrap(), &[0.2, -0.1, 0.3]);
}
