//! G1/G3/G4: directional interference of physical volume, force and stress sources.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::power::{
    CompactMultipolePhasor, PreparedMultipoleRadiatedMomentum, PreparedMultipoleRadiatedPower,
};
use fs_aeroac::curle3d::quadrupole::QuadrupolePhasor;
use fs_aeroac::curle3d::{
    CompactRadiatedMomentum, CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary,
    PreparedCompactRadiatedMomentum, PreparedCompactRadiatedPower, PreparedCompactSourceObserver,
    RadiatedPowerBudget, RadiatedPowerModeBudget,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 331.0).unwrap()
}

fn budget(n: usize) -> RadiatedPowerBudget {
    RadiatedPowerBudget { max_sources: n, max_source_pairs: n * n }
}

fn modes_budget(n: usize, rank: usize) -> RadiatedPowerModeBudget {
    RadiatedPowerModeBudget {
        max_modes: rank,
        max_source_pair_mode_products: n * n * rank,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "actual={a:e}, expected={b:e}, tolerance={tolerance:e}");
}

fn physical(result: CompactRadiatedMomentum, acoustic: FluidAcousticState) {
    let m = result.power_first_moment_w;
    assert!(result.power_w >= 0.0);
    assert!(m[0].hypot(m[1]).hypot(m[2]) <= result.power_w * (1.0 + 3e-12), "{result:?}");
    for axis in 0..3 {
        close(
            result.wave_momentum_flux_n[axis] * acoustic.sound_speed_m_s(),
            m[axis],
            result.power_w * 2e-14,
        );
    }
}

fn same(actual: CompactRadiatedMomentum, expected: CompactRadiatedMomentum, relative: f64) {
    let tolerance = expected.power_w * relative;
    close(actual.power_w, expected.power_w, tolerance);
    for axis in 0..3 {
        close(actual.power_first_moment_w[axis], expected.power_first_moment_w[axis], tolerance);
    }
}

fn scale(mut source: CompactMultipolePhasor, factor: C64) -> CompactMultipolePhasor {
    source.compact.volume_flow_m3_s = source.compact.volume_flow_m3_s * factor;
    source.compact.force_on_fluid_n = source.compact.force_on_fluid_n.map(|v| v * factor);
    source.quadrupole.stress_moment_n_m = source.quadrupole.stress_moment_n_m.map(|v| v * factor);
    source
}

fn add(a: CompactMultipolePhasor, b: CompactMultipolePhasor) -> CompactMultipolePhasor {
    CompactMultipolePhasor {
        compact: CompactSourcePhasor {
            volume_flow_m3_s: a.compact.volume_flow_m3_s + b.compact.volume_flow_m3_s,
            force_on_fluid_n: core::array::from_fn(|i| {
                a.compact.force_on_fluid_n[i] + b.compact.force_on_fluid_n[i]
            }),
        },
        quadrupole: QuadrupolePhasor {
            stress_moment_n_m: core::array::from_fn(|i| {
                a.quadrupole.stress_moment_n_m[i] + b.quadrupole.stress_moment_n_m[i]
            }),
        },
    }
}

fn fixture() -> ([[f64; 3]; 3], [CompactMultipolePhasor; 3]) {
    let positions = [[-0.25, 0.125, 0.0], [0.375, -0.25, 0.125], [0.0, 0.25, -0.375]];
    let sources = core::array::from_fn(|i| {
        let t = (i + 1) as f64;
        CompactMultipolePhasor {
            compact: CompactSourcePhasor {
                volume_flow_m3_s: C64::new((t - 1.5) * 3e-4, (2.5 - t) * 2e-4),
                force_on_fluid_n: [
                    C64::new(0.3 * t, 0.1),
                    C64::new(-0.2, 0.4 * t),
                    C64::new(0.1 / t, -0.3),
                ],
            },
            quadrupole: QuadrupolePhasor {
                stress_moment_n_m: [
                    C64::new(0.11 * t, 0.05), C64::new(-0.07, 0.13 * t),
                    C64::new(0.04 / t, -0.12), C64::new(0.15 * t, -0.06),
                    C64::new(-0.04, -0.09 * t), C64::new(0.07 * t, 0.03 / t),
                ],
            },
        }
    });
    (positions, sources)
}

fn legendre(n: usize, x: f64) -> (f64, f64) {
    let (mut p, mut q) = (1.0, 0.0);
    for j in 1..=n {
        let next = ((2 * j - 1) as f64 * x * p - (j - 1) as f64 * q) / j as f64;
        q = p;
        p = next;
    }
    (p, q)
}

// Independent angular integration of the physical far-pressure amplitude.
// Uses neither production moments nor the traceless-tensor contractions.
fn quadrature(
    acoustic: FluidAcousticState,
    frequency: f64,
    positions: &[[f64; 3]],
    sources: &[CompactMultipolePhasor],
) -> (f64, [f64; 3]) {
    let order = 64;
    let azimuths = 128;
    let k = TAU * frequency / acoustic.sound_speed_m_s();
    let z = acoustic.characteristic_impedance_pa_s_m();
    let mut out = [0.0; 4];
    for index in 0..order {
        let mut mu = (PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..20 {
            let (p, q) = legendre(order, mu);
            let derivative = order as f64 * (mu * p - q) / (mu * mu - 1.0);
            let step = p / derivative;
            mu -= step;
            if step.abs() < 1e-15 { break; }
        }
        let (p, q) = legendre(order, mu);
        let derivative = order as f64 * (mu * p - q) / (mu * mu - 1.0);
        let weight = 2.0 / ((1.0 - mu * mu) * derivative * derivative);
        for j in 0..azimuths {
            let phi = TAU * (j as f64 + 0.5) / azimuths as f64;
            let r = (1.0 - mu * mu).sqrt();
            let n = [r * phi.cos(), r * phi.sin(), mu];
            let mut amplitude = C64::ZERO;
            for (s, source) in positions.iter().zip(sources) {
                let angle = -k * n.iter().zip(s).map(|(a, b)| a * b).sum::<f64>();
                let t = source.quadrupole.stress_moment_n_m;
                let tensor = [[t[0], t[3], t[4]], [t[3], t[1], t[5]], [t[4], t[5], t[2]]];
                let mut value = source.compact.volume_flow_m3_s.scale(z);
                for row in 0..3 {
                    value = value + source.compact.force_on_fluid_n[row].scale(n[row]);
                    for col in 0..3 {
                        value = value + C64::new(0.0, -k) * tensor[row][col].scale(n[row] * n[col]);
                    }
                }
                amplitude = amplitude + value * C64::new(angle.cos(), angle.sin());
            }
            let power = k * k / (8.0 * PI * z) * amplitude.norm_sq()
                * weight / (2.0 * azimuths as f64);
            out[0] += power;
            for axis in 0..3 { out[axis + 1] += n[axis] * power; }
        }
    }
    (out[0], [out[1], out[2], out[3]])
}

#[test]
fn colocated_force_stress_phase_reverses_direction_at_the_same_power() {
    let frequency = 137.0;
    let k = TAU * frequency / medium().sound_speed_m_s();
    let gain2 = k * k / (8.0 * PI * medium().characteristic_impedance_pa_s_m());
    let prepared = PreparedMultipoleRadiatedMomentum::new(medium(), frequency, &[[0.0; 3]], budget(1)).unwrap();
    for angle in [0.0, PI / 3.0, PI / 2.0, PI, -PI / 3.0] {
        let mut source = CompactMultipolePhasor::default();
        source.compact.force_on_fluid_n[0] = C64::ONE;
        // -i k Txx = exp(i angle): the radiative order difference is physical.
        source.quadrupole.stress_moment_n_m[0] = C64::new(-angle.sin(), angle.cos()).scale(1.0 / k);
        let result = prepared.evaluate(&[source]).unwrap();
        close(result.power_w, gain2 * (1.0 / 3.0 + 1.0 / 5.0), gain2 * 3e-14);
        close(result.power_first_moment_w[0], gain2 * (2.0 / 5.0) * angle.cos(), gain2 * 3e-14);
        assert_eq!(result.power_first_moment_w[1..], [0.0; 2]);
        physical(result, medium());
        source.compact = CompactSourcePhasor::default();
        assert_eq!(prepared.evaluate(&[source]).unwrap().power_first_moment_w, [0.0; 3]);
    }
}

#[test]
fn separated_all_mixed_terms_match_independent_physical_angular_quadrature() {
    let (positions, all) = fixture();
    for frequency in [13.0, 217.0, 1800.0] {
        let prepared = PreparedMultipoleRadiatedMomentum::new(medium(), frequency, &positions, budget(3)).unwrap();
        for subset in 0..4 {
            let sources = all.map(|mut s| {
                if subset == 1 || subset == 3 { s.compact.volume_flow_m3_s = C64::ZERO; }
                if subset == 2 || subset == 3 { s.compact.force_on_fluid_n = [C64::ZERO; 3]; }
                s
            });
            let result = prepared.evaluate(&sources).unwrap();
            let (power, moment) = quadrature(medium(), frequency, &positions, &sources);
            close(result.power_w, power, power * 2e-10);
            for axis in 0..3 { close(result.power_first_moment_w[axis], moment[axis], power * 2e-10); }
            close(result.power_w, prepared.power_operator().evaluate(&sources).unwrap().power_w, power * 2e-14);
            physical(result, medium());
        }
    }
}

#[test]
fn weak_tensor_interference_and_tiny_separation_quotients_are_retained() {
    let frequency = 137.0;
    let k = TAU * frequency / medium().sound_speed_m_s();
    let gain2 = k.powi(4) / (8.0 * PI * medium().characteristic_impedance_pa_s_m());
    let mut source = CompactMultipolePhasor::default();
    source.quadrupole.stress_moment_n_m[0] = C64::ONE;
    let x = 1e-10;
    let delta = x * (5.0_f64 / 7.0).sqrt();
    let prepared = PreparedMultipoleRadiatedMomentum::new(medium(), frequency, &[[0.0; 3], [x / k, 0.0, 0.0]], budget(2)).unwrap();
    let result = prepared.evaluate(&[source, scale(source, C64::new(-1.0, -delta))]).unwrap();
    let expected_power = gain2 * (delta * delta / 5.0 + x * x / 7.0);
    close(result.power_w, expected_power, expected_power * 3e-11);
    close(result.power_first_moment_w[0], -gain2 * 2.0 * x * delta / 7.0, expected_power * 3e-11);
    physical(result, medium());

    // j3(x) underflows here, although j3(x)/x^2 and the directional
    // moment remain representable. Quotient series must retain the effect.
    let x = 1e-120;
    let prepared = PreparedMultipoleRadiatedMomentum::new(medium(), frequency, &[[0.0; 3], [x / k, 0.0, 0.0]], budget(2)).unwrap();
    for sign in [-1.0, 1.0] {
        let result = prepared.evaluate(&[source, scale(source, C64::new(0.0, sign))]).unwrap();
        close(result.power_w, gain2 * 2.0 / 5.0, gain2 * 2e-13);
        close(result.power_first_moment_w[0], sign * gain2 * 2.0 * x / 7.0, gain2 * x * 2e-12);
        assert_eq!(result.power_first_moment_w[1..], [0.0; 2]);
    }
}

#[test]
fn trace_cancellation_and_original_compact_results_survive_the_extension() {
    let acoustic = FluidAcousticState::try_new(1.0, 1.0).unwrap();
    let frequency = 1.0 / TAU;
    assert_eq!(TAU * frequency, 1.0);
    let prepared = PreparedMultipoleRadiatedMomentum::new(acoustic, frequency, &[[0.0; 3]; 2], budget(2)).unwrap();
    let mut cancelled = CompactMultipolePhasor::default();
    cancelled.compact.volume_flow_m3_s = C64::new(0.0, 1.0);
    cancelled.quadrupole.stress_moment_n_m[..3].fill(C64::ONE);
    let eta = 1e-10;
    let mut remainder = CompactMultipolePhasor::default();
    remainder.compact.volume_flow_m3_s = C64::from_re(eta);
    remainder.compact.force_on_fluid_n[0] = C64::from_re(eta);
    let result = prepared.evaluate(&[cancelled, remainder]).unwrap();
    let gain2 = 1.0 / (8.0 * PI);
    close(result.power_w, gain2 * eta * eta * 4.0 / 3.0, result.power_w * 2e-13);
    close(result.power_first_moment_w[0], gain2 * eta * eta * 2.0 / 3.0, result.power_w * 2e-13);
    let dark = prepared.evaluate(&[cancelled, CompactMultipolePhasor::default()]).unwrap();
    assert!(dark.power_w < 1e-60);
    assert_eq!(dark.wave_momentum_flux_n, [0.0; 3]);

    let (positions, all) = fixture();
    let compact = all.map(|s| s.compact);
    let multipoles = compact.map(CompactMultipolePhasor::from);
    let original = PreparedCompactRadiatedMomentum::new(medium(), 137.0, &positions, budget(3)).unwrap();
    let extended = PreparedMultipoleRadiatedMomentum::new(medium(), 137.0, &positions, budget(3)).unwrap();
    same(extended.evaluate(&multipoles).unwrap(), original.evaluate(&compact).unwrap(), 3e-13);
}

#[test]
fn translation_rotation_phase_and_original_material_transform_the_vector_physically() {
    let (positions, sources) = fixture();
    let frequency = 217.0;
    let owner = PreparedMultipoleRadiatedMomentum::new(medium(), frequency, &positions, budget(3)).unwrap();
    assert_eq!(owner.power_operator().compact_power().acoustic_state(), medium());
    let result = owner.evaluate(&sources).unwrap();
    let shifted = positions.map(|p| [p[2] + 8.0, p[0] - 4.0, p[1] + 2.0]);
    let changed = sources.map(|mut s| {
        let f = s.compact.force_on_fluid_n;
        s.compact.force_on_fluid_n = [f[2], f[0], f[1]];
        let t = s.quadrupole.stress_moment_n_m;
        s.quadrupole.stress_moment_n_m = [t[2], t[0], t[1], t[4], t[5], t[3]];
        scale(s, C64::new(0.6, 0.8))
    });
    let transformed = PreparedMultipoleRadiatedMomentum::new(medium(), frequency, &shifted, budget(3)).unwrap().evaluate(&changed).unwrap();
    close(transformed.power_w, result.power_w, result.power_w * 3e-13);
    let expected = [result.power_first_moment_w[2], result.power_first_moment_w[0], result.power_first_moment_w[1]];
    for axis in 0..3 { close(transformed.power_first_moment_w[axis], expected[axis], result.power_w * 3e-13); }
    physical(transformed, medium());

    let without_volume = sources.map(|mut s| { s.compact.volume_flow_m3_s = C64::ZERO; s });
    let reference = owner.evaluate(&without_volume).unwrap();
    let acoustic = FluidAcousticState::try_new(2.0 * medium().density_kg_m3(), 3.0 * medium().sound_speed_m_s()).unwrap();
    let changed = PreparedMultipoleRadiatedMomentum::new(acoustic, 3.0 * frequency, &positions, budget(3)).unwrap().evaluate(&without_volume).unwrap();
    close(changed.power_w * 6.0, reference.power_w, reference.power_w * 3e-13);
    for axis in 0..3 {
        close(changed.wave_momentum_flux_n[axis] * 18.0, reference.wave_momentum_flux_n[axis], reference.power_w * 2e-14);
    }
    physical(changed, acoustic);
}

#[test]
fn independent_and_coherence_modes_preserve_intrinsic_force_stress_direction() {
    let (_, sources) = fixture();
    let owner = PreparedMultipoleRadiatedMomentum::new(medium(), 137.0, &[[0.0; 3]; 3], budget(3)).unwrap();
    let sparse: [[CompactMultipolePhasor; 3]; 3] = core::array::from_fn(|i| {
        let mut out = [CompactMultipolePhasor::default(); 3]; out[i] = sources[i]; out
    });
    let independent = owner.evaluate_independent_sources(&sources).unwrap();
    let factors = [&sparse[0][..], &sparse[1][..], &sparse[2][..]];
    same(owner.evaluate_coherence_modes(&factors, modes_budget(3, 3)).unwrap(), independent, 3e-13);
    assert!(independent.power_first_moment_w.iter().any(|v| v.abs() > 0.0));
    let a = sources;
    let b = [sources[2], scale(sources[0], C64::new(0.0, 1.0)), sources[1]];
    let norm = C64::from_re(1.0 / 2.0_f64.sqrt());
    let plus: [_; 3] = core::array::from_fn(|i| scale(add(a[i], b[i]), norm));
    let minus: [_; 3] = core::array::from_fn(|i| scale(add(a[i], scale(b[i], -C64::ONE)), norm));
    let original = owner.evaluate_coherence_modes(&[&a, &b], modes_budget(3, 2)).unwrap();
    let mixed = owner.evaluate_coherence_modes(&[&plus, &minus], modes_budget(3, 2)).unwrap();
    same(mixed, original, 3e-13);
    physical(mixed, medium());
}

#[test]
fn budgets_invalid_stresses_cancellation_and_wall_domain_refuse_without_partial_results() {
    let (positions, sources) = fixture();
    assert!(PreparedMultipoleRadiatedMomentum::new(medium(), 100.0, &positions, budget(2)).is_err());
    let base = PreparedMultipoleRadiatedPower::new(medium(), 100.0, &positions, budget(3)).unwrap();
    let mut prepare_polls = 0;
    let prepared = base.clone().with_momentum_cancellable(|| { prepare_polls += 1; false }).unwrap();
    for stop in [1, prepare_polls / 2, prepare_polls] {
        let mut polls = 0;
        assert!(matches!(base.clone().with_momentum_cancellable(|| { polls += 1; polls == stop }), Err(AeroacError::Cancelled { .. })));
    }
    let mut application_polls = 0;
    let expected = prepared.evaluate_cancellable(&sources, || { application_polls += 1; false }).unwrap();
    for stop in [1, application_polls / 2, application_polls] {
        let mut polls = 0;
        assert!(matches!(prepared.evaluate_cancellable(&sources, || { polls += 1; polls == stop }), Err(AeroacError::Cancelled { .. })));
        assert_eq!(prepared.evaluate(&sources).unwrap(), expected);
    }
    assert!(prepared.evaluate(&sources[..2]).is_err());
    let mut bad = sources;
    bad[1].quadrupole.stress_moment_n_m[4].im = f64::INFINITY;
    assert!(matches!(prepared.evaluate(&bad), Err(AeroacError::NonFinite { .. })));
    assert!(prepared.evaluate_coherence_modes(&[&sources], modes_budget(3, 0)).is_err());
    assert!(prepared.evaluate_coherence_modes(&[&sources[..2]], modes_budget(3, 1)).is_err());
    assert!(matches!(prepared.evaluate_independent_sources_cancellable(&sources, || true), Err(AeroacError::Cancelled { .. })));
    assert!(matches!(prepared.evaluate_coherence_modes_cancellable(&[&sources], modes_budget(3, 1), || true), Err(AeroacError::Cancelled { .. })));
    let empty = PreparedMultipoleRadiatedMomentum::new(medium(), 100.0, &[], budget(0)).unwrap().evaluate(&[]).unwrap();
    assert_eq!(empty, prepared.evaluate_coherence_modes(&[], modes_budget(3, 0)).unwrap());
    assert_eq!(empty.power_w, 0.0);
    assert_eq!(empty.power_first_moment_w, [0.0; 3]);
    let plane = IdealAcousticPlane::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], IdealPlaneBoundary::Rigid).unwrap();
    let observer = PreparedCompactSourceObserver::new(medium(), 100.0, [0.0, 0.0, 2.0], &positions, 3).unwrap().with_plane(plane, 6).unwrap();
    let images = PreparedCompactRadiatedPower::from_half_space(&observer, budget(6)).unwrap().with_quadrupoles().unwrap();
    assert!(matches!(images.with_momentum(), Err(AeroacError::InvalidParameter { .. })));
}
