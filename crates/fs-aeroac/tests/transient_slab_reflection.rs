//! G1/G3/G4: complete causal slab returns, source images and original clocks.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::{CompactForceObserver, IdealAcousticPlane, IdealPlaneBoundary};
use fs_aeroac::curle3d::quadrupole::transient::PreparedTransientQuadrupoleObserver;
use fs_aeroac::curle3d::slab::{
    CausalSlabWindow, IdealAcousticSlab, PreparedSlabSourceObserver, TransientAcousticSource,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::dd::Dd;
use std::f64::consts::PI;

type Matrix = [[f64; 3]; 3];
const SOURCE: [f64; 3] = [0.0, 0.0, 0.6];
const RECEIVER: [f64; 3] = [0.4, 0.0, 0.9];
const WIDTH: f64 = 2.0;
const RIGID: IdealPlaneBoundary = IdealPlaneBoundary::Rigid;
const RELEASE: IdealPlaneBoundary = IdealPlaneBoundary::PressureRelease;

fn medium() -> FluidAcousticState { FluidAcousticState::try_new(1.3, 10.0).unwrap() }

fn slab(boundaries: [IdealPlaneBoundary; 2]) -> IdealAcousticSlab {
    let lower = IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], boundaries[0]).unwrap();
    IdealAcousticSlab::new(lower, WIDTH, boundaries[1]).unwrap()
}

fn prepare(
    slab: IdealAcousticSlab,
    receiver: [f64; 3],
    sources: &[[f64; 3]],
    epoch: f64,
    duration: f64,
    maximum_paths: usize,
) -> PreparedSlabSourceObserver {
    PreparedSlabSourceObserver::new(slab, medium(), receiver, sources,
        CausalSlabWindow { source_epoch_s: epoch, duration_s: duration }, maximum_paths, || false).unwrap()
}

fn source(time: f64, complete: bool) -> TransientAcousticSource {
    let mut result = TransientAcousticSource::default();
    if time <= 0.0 { return result; }
    let d = 0.1;
    let t = time / d;
    result.compact.volume_flow_m3_s = 0.01 * t.powi(3);
    result.compact.volume_flow_rate_m3_s2 = 0.03 * t.powi(2) / d;
    if complete {
        let f = [0.3, -0.2, 0.7];
        result.force_impulse_n_s = f.map(|v| v * d * t.powi(4) / 4.0);
        result.compact.force_on_fluid_n = f.map(|v| v * t.powi(3));
        result.compact.force_rate_on_fluid_n_s = f.map(|v| v * 3.0 * t.powi(2) / d);
        let q = [0.7, -0.3, 0.9, 0.23, -0.11, 0.4];
        result.stress.stress_integral_n_m_s = q.map(|v| v * d * t.powi(4) / 4.0);
        result.stress.stress_moment_n_m = q.map(|v| v * t.powi(3));
        result.stress.stress_rate_n_m_s = q.map(|v| v * 3.0 * t.powi(2) / d);
        result.stress.stress_acceleration_n_m_s2 = q.map(|v| v * 6.0 * t / (d*d));
    }
    result
}

fn vector(v: [f64; 3], h: Matrix, sign: f64) -> [f64; 3] {
    core::array::from_fn(|row| sign * (0..3).map(|col| h[row][col] * v[col]).sum::<f64>())
}

fn tensor(q: [f64; 6], h: Matrix, sign: f64) -> [f64; 6] {
    let [xx, yy, zz, xy, xz, yz] = q;
    let q = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
    [(0,0), (1,1), (2,2), (0,1), (0,2), (1,2)].map(|(i,j)| {
        let mut value = 0.0;
        for a in 0..3 { for b in 0..3 { value += h[i][a] * q[a][b] * h[j][b]; } }
        sign * value
    })
}

fn transformed(mut value: TransientAcousticSource, h: Matrix, sign: f64) -> TransientAcousticSource {
    value.compact.volume_flow_m3_s *= sign;
    value.compact.volume_flow_rate_m3_s2 *= sign;
    value.compact.force_on_fluid_n = vector(value.compact.force_on_fluid_n, h, sign);
    value.compact.force_rate_on_fluid_n_s = vector(value.compact.force_rate_on_fluid_n_s, h, sign);
    value.force_impulse_n_s = vector(value.force_impulse_n_s, h, sign);
    value.stress.stress_integral_n_m_s = tensor(value.stress.stress_integral_n_m_s, h, sign);
    value.stress.stress_moment_n_m = tensor(value.stress.stress_moment_n_m, h, sign);
    value.stress.stress_rate_n_m_s = tensor(value.stress.stress_rate_n_m_s, h, sign);
    value.stress.stress_acceleration_n_m_s2 = tensor(value.stress.stress_acceleration_n_m_s2, h, sign);
    value
}

fn free_field(receiver: [f64; 3], position: [f64; 3], value: TransientAcousticSource) -> [f64; 4] {
    let d: [f64; 3] = core::array::from_fn(|i| receiver[i] - position[i]);
    let r = d.iter().map(|v| v*v).sum::<f64>().sqrt();
    let rho = medium().density_kg_m3();
    let c = medium().sound_speed_m_s();
    let green = 1.0 / (4.0 * PI * r);
    let mut field = [0.0; 4];
    field[0] = rho * green * value.compact.volume_flow_rate_m3_s2;
    for i in 0..3 {
        field[i+1] = d[i] / r * green * (value.compact.volume_flow_m3_s / r
            + value.compact.volume_flow_rate_m3_s2 / c);
    }
    let force = CompactForceObserver::new(medium(), receiver, position).unwrap()
        .field_at_retarded(value.force_impulse_n_s, value.compact.force_on_fluid_n,
            value.compact.force_rate_on_fluid_n_s).unwrap();
    let stress = PreparedTransientQuadrupoleObserver::new(medium(), receiver, &[position], 1).unwrap()
        .evaluate_at_retarded(&[value.stress]).unwrap();
    field[0] += force.pressure_pa + stress.pressure_pa;
    for i in 0..3 { field[i+1] += force.particle_velocity_m_s[i] + stress.particle_velocity_m_s[i]; }
    field
}

// Independent PHYSICAL source images. Production instead transforms virtual
// receivers and reuses original source tensors. The much larger lattice also
// checks that the finite horizon omits no active late return.
fn reference(receiver: [f64; 3], time: f64, boundaries: [IdealPlaneBoundary; 2], complete: bool) -> [f64; 4] {
    let mut sum = [0.0; 4];
    for k in -30_i32..=30 {
        for mirrored in [false, true] {
            let mixed = boundaries[0] != boundaries[1];
            let mut sign = if mixed && k % 2 != 0 { -1.0 } else { 1.0 };
            if mirrored && boundaries[0] == RELEASE { sign = -sign; }
            let flip = if mirrored { -1.0 } else { 1.0 };
            let h = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,flip]];
            let position = [SOURCE[0], SOURCE[1], 2.0 * f64::from(k) * WIDTH + flip * SOURCE[2]];
            let r = (0..3).map(|i| (receiver[i] - position[i]).powi(2)).sum::<f64>().sqrt();
            let emission = time - r / medium().sound_speed_m_s();
            if emission < 0.0 { continue; }
            let field = free_field(receiver, position, transformed(source(emission, complete), h, sign));
            for i in 0..4 { sum[i] += field[i]; }
        }
    }
    sum
}

fn close(actual: [f64; 4], expected: [f64; 4]) {
    for i in 0..4 {
        assert!((actual[i] - expected[i]).abs() < 3e-11 * expected[i].abs().max(1.0),
            "channel={i}, actual={:e}, expected={:e}", actual[i], expected[i]);
    }
}

#[test]
fn g1_first_and_repeated_monopole_returns_are_causal_and_complete() {
    let owner = prepare(slab([RIGID, RIGID]), RECEIVER, &[SOURCE], 0.0, 1.0, 100);
    assert_eq!(owner.source_count(), 1);
    assert_eq!(owner.source_position_m(0), Some(SOURCE));
    assert_eq!(owner.source_position_m(1), None);
    assert!(owner.path_count() > 4);
    for time in [0.0, 0.03, 0.1, 0.2, 0.3, 0.5, 1.0] {
        let mut calls = Vec::new();
        let result = owner.evaluate_initially_resting(time, |index, t| {
            assert_eq!(index, 0);
            assert!(!t.lt(Dd::ZERO));
            calls.push(t);
            Ok(source(t.to_f64(), false))
        }, || false).unwrap();
        let actual = [result.pressure_pa, result.particle_velocity_m_s[0],
            result.particle_velocity_m_s[1], result.particle_velocity_m_s[2]];
        close(actual, reference(RECEIVER, time, [RIGID, RIGID], false));
        let active = (0..owner.path_count()).filter(|&i| owner.path_delay_s(i).unwrap() <= time).count();
        assert_eq!(calls.len(), active + 1); // one initial-data check
        for i in 0..3 {
            assert_eq!(result.instantaneous_intensity_w_m2[i], result.pressure_pa * result.particle_velocity_m_s[i]);
        }
        if time < 0.05 { assert_eq!(actual, [0.0; 4]); }
    }
    let direct = free_field(RECEIVER, SOURCE, source(0.1 - 0.05, false));
    close(reference(RECEIVER, 0.1, [RIGID, RIGID], false), direct);
    assert!(reference(RECEIVER, 0.5, [RIGID, RIGID], false)[0]
        > free_field(RECEIVER, SOURCE, source(0.5 - 0.05, false))[0]);
}

#[test]
fn g1_all_wall_signs_and_full_force_tensor_images_enforce_both_boundaries() {
    for boundaries in [[RIGID,RIGID], [RELEASE,RELEASE], [RIGID,RELEASE], [RELEASE,RIGID]] {
        for receiver in [RECEIVER, [0.4,0.0,0.0], [0.4,0.0,WIDTH]] {
            let owner = prepare(slab(boundaries), receiver, &[SOURCE], 0.0, 1.0, 100);
            for time in [0.2, 0.5, 0.9, 1.0] {
                let field = owner.evaluate_initially_resting(time, |_, t| Ok(source(t.to_f64(), true)), || false).unwrap();
                close([field.pressure_pa, field.particle_velocity_m_s[0], field.particle_velocity_m_s[1],
                    field.particle_velocity_m_s[2]], reference(receiver, time, boundaries, true));
                let face = if receiver[2] == 0.0 { Some(0) } else if receiver[2] == WIDTH { Some(1) } else { None };
                if let Some(face) = face {
                    if boundaries[face] == RIGID { assert!(field.particle_velocity_m_s[2].abs() < 3e-11); }
                    else { assert!(field.pressure_pa.abs() < 3e-10); }
                }
            }
        }
    }
}

#[test]
fn g3_original_source_clock_and_rotated_slab_preserve_the_complete_field() {
    let h = [[0.8,0.0,0.6], [0.0,1.0,0.0], [-0.6,0.0,0.8]];
    let origin = [12.0,-5.0,7.0];
    let position = |p| {
        let p = vector(p, h, 1.0);
        core::array::from_fn(|i| p[i] + origin[i])
    };
    let lower = IdealAcousticPlane::new(origin, [0.6,0.0,0.8], RIGID).unwrap();
    let rotated_slab = IdealAcousticSlab::new(lower, WIDTH, RELEASE).unwrap();
    let epoch = 2.0_f64.powi(60);
    let owner = prepare(rotated_slab, position(RECEIVER), &[position(SOURCE)], epoch, 1.0, 100);
    for time in [0.1, 0.5, 1.0] {
        let field = owner.evaluate_initially_resting(time, |_, t| {
            Ok(transformed(source((t - Dd::from_f64(epoch)).to_f64(), true), h, 1.0))
        }, || false).unwrap();
        let expected = reference(RECEIVER, time, [RIGID,RELEASE], true);
        let velocity = vector([expected[1],expected[2],expected[3]], h, 1.0);
        close([field.pressure_pa, field.particle_velocity_m_s[0],field.particle_velocity_m_s[1],
            field.particle_velocity_m_s[2]], [expected[0],velocity[0],velocity[1],velocity[2]]);
    }
}

#[test]
fn g3_window_extension_and_coherent_physical_sources_preserve_early_sound() {
    let a = prepare(slab([RIGID,RELEASE]), RECEIVER, &[SOURCE], 0.0, 0.5, 100);
    let b = prepare(slab([RIGID,RELEASE]), RECEIVER, &[SOURCE], 0.0, 1.0, 100);
    assert!(b.path_count() > a.path_count());
    let f = |owner: &PreparedSlabSourceObserver| owner.evaluate_initially_resting(0.5,
        |_, t| Ok(source(t.to_f64(), true)), || false).unwrap();
    assert_eq!(f(&a), f(&b));
    let pair = prepare(slab([RIGID,RELEASE]), RECEIVER, &[SOURCE, SOURCE], 0.0, 1.0, 200);
    let identity = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]];
    let silent = pair.evaluate_initially_resting(0.9, |i,t| {
        Ok(transformed(source(t.to_f64(), true), identity, if i==0 { 1.0 } else { -1.0 }))
    }, || false).unwrap();
    close([silent.pressure_pa, silent.particle_velocity_m_s[0], silent.particle_velocity_m_s[1],
        silent.particle_velocity_m_s[2]], [0.0;4]);
    assert!(silent.instantaneous_intensity_w_m2.iter().all(|v| v.abs() < 1e-22));
    let empty = prepare(slab([RIGID,RIGID]), RECEIVER, &[], 0.0, 1.0, 0);
    assert_eq!(empty.path_count(), 0);
    assert_eq!(empty.evaluate_initially_resting(1.0, |_,_| panic!("no sources"), || false).unwrap().pressure_pa, 0.0);
}

#[test]
fn g0_g4_complete_window_geometry_onset_and_work_boundaries_refuse_without_a_field() {
    let domain = slab([RIGID,RELEASE]);
    let owner = prepare(domain, RECEIVER, &[SOURCE], 0.0, 0.5, 100);
    let window = owner.window();
    assert!(PreparedSlabSourceObserver::new(domain, medium(), RECEIVER, &[SOURCE], window,
        owner.path_count()-1, || false).is_err());
    assert!(PreparedSlabSourceObserver::new(domain, medium(), RECEIVER, &[SOURCE], window,
        owner.path_count(), || false).is_ok());
    for position in [[0.0,0.0,0.0], [0.0,0.0,WIDTH], [0.0,0.0,WIDTH+0.1]] {
        assert!(PreparedSlabSourceObserver::new(domain, medium(), RECEIVER, &[position], window,
            100, || false).is_err());
    }
    assert!(PreparedSlabSourceObserver::new(domain, medium(), [0.0,0.0,-0.1], &[SOURCE], window,
        100, || false).is_err());
    for time in [-0.01, 0.5_f64.next_up(), f64::NAN] {
        assert!(owner.evaluate_initially_resting(time, |_,_| panic!("invalid window"), || false).is_err());
    }
    assert!(owner.evaluate_initially_resting(0.5, |_,_| Ok(source(0.01, true)), || false).is_err());
    assert!(matches!(PreparedSlabSourceObserver::new(domain, medium(), RECEIVER, &[SOURCE], window,
        100, || true), Err(AeroacError::Cancelled { .. })));
    assert!(matches!(owner.evaluate_initially_resting(0.5, |_,_| panic!("cancelled"), || true),
        Err(AeroacError::Cancelled { .. })));
    assert!(owner.evaluate_initially_resting(0.5, |_,t| Ok(source(t.to_f64(), true)), || false).is_ok());
}
