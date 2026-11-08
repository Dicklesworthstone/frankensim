//! G1/G3/G4: finite-time outgoing radiation of the original physical multipoles.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::power::{
    CompactMultipolePhasor, MultipoleRadiationLine, PreparedMultipoleRadiatedEnergy,
    PreparedMultipoleRadiatedPower, RadiatedEnergyBudget, WindowedRadiatedEnergy,
};
use fs_aeroac::curle3d::quadrupole::QuadrupolePhasor;
use fs_aeroac::curle3d::{
    CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary,
    PreparedCompactRadiatedPower, PreparedCompactSourceObserver, RadiatedPowerBudget,
};
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn medium() -> FluidAcousticState { FluidAcousticState::try_new(1.3, 331.0).unwrap() }
fn power_budget(n: usize) -> RadiatedPowerBudget {
    RadiatedPowerBudget { max_sources: n, max_source_pairs: n * n }
}
fn budget(lines: usize, sources: usize) -> RadiatedEnergyBudget {
    RadiatedEnergyBudget {
        max_lines: lines,
        max_source_terms: sources,
        max_pair_work: lines * lines + sources * sources,
        max_workspace_bytes: 1 << 20,
    }
}
fn owner(f: f64, positions: &[[f64; 3]]) -> PreparedMultipoleRadiatedPower {
    PreparedMultipoleRadiatedPower::new(medium(), f, positions, power_budget(positions.len())).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "actual={a:e}, expected={b:e}, tolerance={tolerance:e}");
}
fn physical(value: WindowedRadiatedEnergy) {
    assert!(value.radiated_energy_j >= 0.0);
    let m = value.energy_first_moment_j;
    assert!(m[0].hypot(m[1]).hypot(m[2]) <= value.radiated_energy_j * (1.0 + 3e-11), "{value:?}");
    close(value.mean_power_w * value.window.duration_s(), value.radiated_energy_j, value.radiated_energy_j * 3e-14);
    for i in 0..3 {
        close(value.wave_momentum_n_s[i] * medium().sound_speed_m_s(), m[i], value.radiated_energy_j * 3e-14);
    }
}
fn same(a: WindowedRadiatedEnergy, b: WindowedRadiatedEnergy, relative: f64) {
    close(a.radiated_energy_j, b.radiated_energy_j, b.radiated_energy_j * relative);
    for i in 0..3 {
        close(a.energy_first_moment_j[i], b.energy_first_moment_j[i], b.radiated_energy_j * relative);
    }
}
fn scale(mut s: CompactMultipolePhasor, factor: C64) -> CompactMultipolePhasor {
    s.compact.volume_flow_m3_s = s.compact.volume_flow_m3_s * factor;
    s.compact.force_on_fluid_n = s.compact.force_on_fluid_n.map(|v| v * factor);
    s.quadrupole.stress_moment_n_m = s.quadrupole.stress_moment_n_m.map(|v| v * factor);
    s
}
fn fixture(offset: f64) -> ([[f64; 3]; 3], [CompactMultipolePhasor; 3]) {
    let positions = [[-0.25, 0.125, 0.0], [0.375, -0.25, 0.125], [0.0, 0.25, -0.375]];
    let sources = core::array::from_fn(|i| {
        let t = i as f64 + 1.0 + offset;
        CompactMultipolePhasor {
            compact: CompactSourcePhasor {
                volume_flow_m3_s: C64::new((t - 1.5) * 3e-4, (2.5 - t) * 2e-4),
                force_on_fluid_n: [C64::new(0.3 * t, 0.1), C64::new(-0.2, 0.4 * t), C64::new(0.1 / t, -0.3)],
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
fn monopole(value: f64) -> CompactMultipolePhasor {
    CompactSourcePhasor::monopole(C64::from_re(value)).into()
}
fn lines<'p>(
    owners: &'p [PreparedMultipoleRadiatedPower],
    banks: &[&'p [CompactMultipolePhasor]],
) -> Vec<MultipoleRadiationLine<'p>> {
    owners.iter().zip(banks).map(|(operator, &sources)| MultipoleRadiationLine { operator, sources }).collect()
}
fn gauss(n: usize) -> Vec<(f64, f64)> {
    let polynomial = |x| {
        let (mut a, mut b) = (1.0, x);
        for j in 2..=n {
            let next = ((2 * j - 1) as f64 * x * b - (j - 1) as f64 * a) / j as f64;
            a = b; b = next;
        }
        (b, n as f64 * (x * b - a) / (x * x - 1.0))
    };
    (0..n).map(|i| {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..30 {
            let (p, derivative) = polynomial(x);
            let change = p / derivative;
            x -= change;
            if change.abs() < 2e-16 { break; }
        }
        let (_, derivative) = polynomial(x);
        (x, 2.0 / ((1.0 - x * x) * derivative * derivative))
    }).collect()
}

// Independent angular AND temporal quadrature of the actual real far field.
// No production harmonic-window products, spherical moments, source trace
// decomposition or analytic interference contractions enter this oracle.
fn direct_energy(
    frequencies: &[f64],
    positions: &[[f64; 3]],
    banks: &[&[CompactMultipolePhasor]],
    window: HarmonicWindow,
) -> [f64; 4] {
    let z = medium().characteristic_impedance_pa_s_m();
    let time = gauss(64);
    let mut result = [0.0; 4];
    for (mu, weight) in gauss(32) {
        let radius = (1.0 - mu * mu).sqrt();
        for azimuth in 0..64 {
            let phi = TAU * (azimuth as f64 + 0.5) / 64.0;
            let n = [radius * phi.cos(), radius * phi.sin(), mu];
            let amplitudes: Vec<_> = frequencies.iter().zip(banks).map(|(&f, sources)| {
                let k = TAU * f / medium().sound_speed_m_s();
                let mut pressure = C64::ZERO;
                for (position, source) in positions.iter().zip(*sources) {
                    let phase = -k * position.iter().zip(n).map(|(x, n)| x * n).sum::<f64>();
                    let t = source.quadrupole.stress_moment_n_m;
                    let tensor = [[t[0], t[3], t[4]], [t[3], t[1], t[5]], [t[4], t[5], t[2]]];
                    let mut amplitude = source.compact.volume_flow_m3_s.scale(z);
                    for i in 0..3 {
                        amplitude = amplitude + source.compact.force_on_fluid_n[i].scale(n[i]);
                        for j in 0..3 {
                            amplitude = amplitude + C64::new(0.0, -k) * tensor[i][j].scale(n[i] * n[j]);
                        }
                    }
                    pressure = pressure + amplitude * C64::new(phase.cos(), phase.sin())
                        * C64::new(0.0, -k / (8.0 * PI * z).sqrt());
                }
                pressure
            }).collect();
            let mut energy = 0.0;
            for &(point, temporal_weight) in &time {
                let tau = window.start_s() + window.duration_s() * (1.0 + point) / 2.0;
                let real = amplitudes.iter().zip(frequencies).map(|(a, f)| {
                    let phase = TAU * f * tau;
                    a.re * phase.cos() + a.im * phase.sin()
                }).sum::<f64>();
                energy += real * real * temporal_weight * window.duration_s();
            }
            energy *= weight / 128.0;
            result[0] += energy;
            for i in 0..3 { result[i + 1] += n[i] * energy; }
        }
    }
    result
}

#[test]
fn complete_mixed_spectrum_matches_independent_space_time_flux_integration() {
    let (positions, a) = fixture(0.0);
    let (_, b) = fixture(0.5);
    let frequencies = [137.0, 251.0];
    let owners = frequencies.map(|f| owner(f, &positions));
    for subset in 0..4 {
        let select = |mut s: CompactMultipolePhasor| {
            if subset == 1 || subset == 3 { s.compact.volume_flow_m3_s = C64::ZERO; }
            if subset == 2 || subset == 3 { s.compact.force_on_fluid_n = [C64::ZERO; 3]; }
            s
        };
        let a = a.map(select);
        let b = b.map(select);
        let banks = [&a[..], &b[..]];
        let prepared = PreparedMultipoleRadiatedEnergy::new(&lines(&owners, &banks), [0.0; 3], budget(2, 6)).unwrap();
        for (start, duration) in [(0.0, 0.003), (0.001, 0.0002), (-0.002, 0.011)] {
            let window = HarmonicWindow::new(start, duration).unwrap();
            let expected = direct_energy(&frequencies, &positions, &banks, window);
            let actual = prepared.observe(window).unwrap();
            close(actual.radiated_energy_j, expected[0], expected[0] * 2e-10);
            for axis in 0..3 {
                close(actual.energy_first_moment_j[axis], expected[axis + 1], expected[0] * 2e-10);
            }
            physical(actual);
        }
    }
}

fn one_minus_sinc(x: f64) -> f64 {
    let mut term = -x * x / 6.0;
    let mut sum = term;
    for n in 2..=32 { term *= -x * x / ((2 * n) * (2 * n + 1)) as f64; sum += term; }
    -sum
}

#[test]
fn short_centered_windows_retain_the_real_monopole_sum_frequency_term() {
    let frequency = 37.0;
    let operator = owner(frequency, &[[0.0; 3]]);
    let source = [monopole(1.0)];
    let prepared = PreparedMultipoleRadiatedEnergy::new(
        &[MultipoleRadiationLine { operator: &operator, sources: &source }], [0.0; 3], budget(1, 1),
    ).unwrap();
    let cycle_power = operator.evaluate(&source).unwrap().power_w;
    for periods in [1e-8, 0.2, 1.0] {
        let duration = periods / frequency;
        let result = prepared.observe(HarmonicWindow::new(-duration / 2.0, duration).unwrap()).unwrap();
        let expected = cycle_power * duration * one_minus_sinc(TAU * frequency * duration);
        close(result.radiated_energy_j, expected, expected * 3e-13);
        assert_eq!(result.energy_first_moment_j, [0.0; 3]);
        physical(result);
    }
}

#[test]
fn the_explicit_sphere_center_preserves_physical_direction_dependent_arrival_times() {
    let frequency = 37.0;
    let k = TAU * frequency / medium().sound_speed_m_s();
    let displacement = 0.7 / k;
    let operator = owner(frequency, &[[displacement, 0.0, 0.0]]);
    let source = [monopole(1.0)];
    let bank = [MultipoleRadiationLine { operator: &operator, sources: &source }];
    let prepared = PreparedMultipoleRadiatedEnergy::new(&bank, [0.0; 3], budget(1, 1)).unwrap();
    let window = HarmonicWindow::new(0.13 / frequency, 0.17 / frequency).unwrap();
    let cycle_power = operator.evaluate(&source).unwrap().power_w;
    let x = 2.0 * k * displacement;
    let j0 = x.sin() / x;
    let j1 = x.sin() / (x * x) - x.cos() / x;
    let half_span = TAU * frequency * window.duration_s();
    let sensitivity = half_span.sin() / half_span;
    let phase = 2.0 * TAU * frequency * (window.start_s() + window.duration_s() / 2.0);
    let actual = prepared.observe(window).unwrap();
    close(actual.mean_power_w, cycle_power * (1.0 - j0 * sensitivity * phase.cos()), cycle_power * 3e-13);
    close(actual.mean_power_first_moment_w[0], cycle_power * j1 * sensitivity * phase.sin(), cycle_power * 3e-13);
    assert_eq!(actual.mean_power_first_moment_w[1..], [0.0; 2]);
    let centered = PreparedMultipoleRadiatedEnergy::new(&bank, [displacement, 0.0, 0.0], budget(1, 1)).unwrap().observe(window).unwrap();
    assert_eq!(centered.energy_first_moment_j, [0.0; 3]);
    assert!((centered.radiated_energy_j - actual.radiated_energy_j).abs() > 0.01 * actual.radiated_energy_j);
    physical(actual);
}

#[test]
fn coherent_beats_survive_without_rounding_a_short_window_into_the_source_epoch() {
    let frequencies = [1000.0, 1001.0];
    let owners = frequencies.map(|f| owner(f, &[[0.0; 3]]));
    let a = [monopole(1.0)];
    let b = [monopole(-frequencies[0] / frequencies[1])];
    let banks = [&a[..], &b[..]];
    let prepared = PreparedMultipoleRadiatedEnergy::new(&lines(&owners, &banks), [0.0; 3], budget(2, 2)).unwrap();
    let short = HarmonicWindow::new(0.0, 0.0002).unwrap();
    let expected = direct_energy(&frequencies, &[[0.0; 3]], &banks, short);
    let actual = prepared.observe(short).unwrap();
    close(actual.radiated_energy_j, expected[0], expected[0] * 3e-9);
    let incoherent = (owners[0].evaluate(&a).unwrap().power_w + owners[1].evaluate(&b).unwrap().power_w) * short.duration_s();
    assert!(actual.radiated_energy_j < 1e-5 * incoherent);
    // Integer-Hz sources have exactly the same phase after an integer number
    // of seconds. The duration is kept independently of the huge timestamp.
    let shifted = prepared.observe(HarmonicWindow::new(1e12, short.duration_s()).unwrap()).unwrap();
    same(shifted, actual, 3e-12);
    physical(actual);
}

#[test]
fn interval_additivity_and_long_period_limits_match_original_radiation_owners() {
    let (positions, a) = fixture(0.0);
    let (_, b) = fixture(0.5);
    let owners = [owner(137.0, &positions), owner(251.0, &positions)];
    let prepared = PreparedMultipoleRadiatedEnergy::new(&lines(&owners, &[&a, &b]), [0.0; 3], budget(2, 6)).unwrap();
    let (start, first, second) = (0.0013, 0.0011, 0.0023);
    let whole = prepared.observe(HarmonicWindow::new(start, first + second).unwrap()).unwrap();
    let left = prepared.observe(HarmonicWindow::new(start, first).unwrap()).unwrap();
    let right = prepared.observe(HarmonicWindow::new(start + first, second).unwrap()).unwrap();
    close(whole.radiated_energy_j, left.radiated_energy_j + right.radiated_energy_j, whole.radiated_energy_j * 3e-12);
    for axis in 0..3 {
        close(whole.energy_first_moment_j[axis], left.energy_first_moment_j[axis] + right.energy_first_moment_j[axis], whole.radiated_energy_j * 3e-12);
    }
    let single = PreparedMultipoleRadiatedEnergy::new(
        &[MultipoleRadiationLine { operator: &owners[0], sources: &a }], [0.0; 3], budget(1, 3),
    ).unwrap();
    let cycle = owners[0].clone().with_momentum().unwrap().evaluate(&a).unwrap();
    let averaged = single.observe(HarmonicWindow::new(start, 3.0 / 137.0).unwrap()).unwrap();
    close(averaged.mean_power_w, cycle.power_w, cycle.power_w * 3e-12);
    for axis in 0..3 {
        close(averaged.mean_power_first_moment_w[axis], cycle.power_first_moment_w[axis], cycle.power_w * 3e-12);
    }
    physical(whole);
}

#[test]
fn physical_source_coalescing_and_joint_frame_translation_preserve_weak_interference() {
    let (_, sample) = fixture(0.0);
    let operator = owner(137.0, &[[0.0; 3]]);
    let weak = [sample[0]];
    let large = [scale(weak[0], C64::from_re(1e30))];
    let opposite = [scale(large[0], -C64::ONE)];
    let window = HarmonicWindow::new(0.001, 0.002).unwrap();
    let baseline = PreparedMultipoleRadiatedEnergy::new(
        &[MultipoleRadiationLine { operator: &operator, sources: &weak }], [0.0; 3], budget(1, 1),
    ).unwrap().observe(window).unwrap();
    let input = [&large[..], &weak[..], &opposite[..]].map(|sources| MultipoleRadiationLine { operator: &operator, sources });
    let combined = PreparedMultipoleRadiatedEnergy::new(&input, [0.0; 3], budget(3, 3)).unwrap();
    assert_eq!(combined.frequency_count(), 1);
    same(combined.observe(window).unwrap(), baseline, 3e-13);
    assert!(core::ptr::eq(combined.source_lines()[0].operator, &operator));
    let dark = PreparedMultipoleRadiatedEnergy::new(&[input[0], input[2]], [0.0; 3], budget(2, 2)).unwrap().observe(window).unwrap();
    assert_eq!(dark.radiated_energy_j, 0.0);
    assert_eq!(dark.wave_momentum_n_s, [0.0; 3]);

    // A subwavelength cancelling pair retains the O((k d)^2) real-field
    // energy; subtracting rounded cycle powers would erase it.
    let k = TAU * 137.0 / medium().sound_speed_m_s();
    let x = 1e-10;
    let pair_owner = owner(137.0, &[[0.0; 3], [x / k, 0.0, 0.0]]);
    let pair_sources = [monopole(1.0), monopole(-1.0)];
    let duration = 1e-6 / 137.0;
    let pair_window = HarmonicWindow::new(-duration / 2.0, duration).unwrap();
    let radiated = PreparedMultipoleRadiatedEnergy::new(
        &[MultipoleRadiationLine { operator: &pair_owner, sources: &pair_sources }],
        [0.0; 3], budget(1, 2),
    ).unwrap().observe(pair_window).unwrap();
    let single_power = operator.evaluate(&[monopole(1.0)]).unwrap().power_w;
    let expected = single_power * x * x * duration / 3.0
        * (2.0 - one_minus_sinc(TAU * 137.0 * duration));
    close(radiated.radiated_energy_j, expected, expected * 3e-12);
    physical(radiated);

    let (positions, a) = fixture(0.0);
    let (_, b) = fixture(0.5);
    let frequencies = [137.0, 251.0];
    let owners = frequencies.map(|f| owner(f, &positions));
    let original = PreparedMultipoleRadiatedEnergy::new(&lines(&owners, &[&a, &b]), [0.0; 3], budget(2, 6)).unwrap().observe(window).unwrap();
    let translation = [8.0, -4.0, 2.0];
    let shifted_positions = positions.map(|p| core::array::from_fn(|i| p[i] + translation[i]));
    let moved = frequencies.map(|f| owner(f, &shifted_positions));
    let shifted = PreparedMultipoleRadiatedEnergy::new(&lines(&moved, &[&a, &b]), translation, budget(2, 6)).unwrap().observe(window).unwrap();
    same(shifted, original, 3e-13);
}

#[test]
fn exact_budgets_common_medium_uniform_windows_and_cancellation_are_enforced() {
    let (positions, a) = fixture(0.0);
    let (_, b) = fixture(0.5);
    let owners = [owner(137.0, &positions), owner(251.0, &positions)];
    let banks = lines(&owners, &[&a, &b]);
    let mut prepare_polls = 0;
    let prepared = PreparedMultipoleRadiatedEnergy::new_cancellable(&banks, [0.0; 3], budget(2, 6), || { prepare_polls += 1; false }).unwrap();
    assert_eq!(prepared.preparation_pair_work(), 40);
    let exact = RadiatedEnergyBudget { max_workspace_bytes: prepared.workspace_bytes(), ..budget(2, 6) };
    assert!(PreparedMultipoleRadiatedEnergy::new(&banks, [0.0; 3], exact).is_ok());
    for which in 0..4 {
        let mut short = exact;
        match which { 0 => short.max_lines -= 1, 1 => short.max_source_terms -= 1, 2 => short.max_pair_work -= 1, _ => short.max_workspace_bytes -= 1 }
        assert!(PreparedMultipoleRadiatedEnergy::new(&banks, [0.0; 3], short).is_err());
    }
    for stop in [1, prepare_polls / 2, prepare_polls] {
        let mut polls = 0;
        assert!(matches!(PreparedMultipoleRadiatedEnergy::new_cancellable(&banks, [0.0; 3], exact, || { polls += 1; polls == stop }), Err(AeroacError::Cancelled { .. })));
    }
    let window = HarmonicWindow::new(0.001, 0.002).unwrap();
    let mut application_polls = 0;
    let expected = prepared.observe_cancellable(window, || { application_polls += 1; false }).unwrap();
    for stop in [1, application_polls / 2, application_polls] {
        let mut polls = 0;
        assert!(matches!(prepared.observe_cancellable(window, || { polls += 1; polls == stop }), Err(AeroacError::Cancelled { .. })));
        assert_eq!(prepared.observe(window).unwrap(), expected);
    }
    let triangular = HarmonicWindow::with_profile(0.0, 0.01, HarmonicWindowProfile::Triangular).unwrap();
    assert!(prepared.observe(triangular).is_err());
    assert!(PreparedMultipoleRadiatedEnergy::new(&[], [0.0; 3], budget(0, 0)).is_err());
    assert!(PreparedMultipoleRadiatedEnergy::new(&banks, [f64::NAN; 3], exact).is_err());
    for acoustic in [
        FluidAcousticState::try_new(2.6, 331.0).unwrap(),
        FluidAcousticState::try_new(1.3, 332.0).unwrap(),
    ] {
        let wrong_medium = PreparedMultipoleRadiatedPower::new(acoustic, 251.0, &positions, power_budget(3)).unwrap();
        assert!(PreparedMultipoleRadiatedEnergy::new(
            &[banks[0], MultipoleRadiationLine { operator: &wrong_medium, sources: &b }], [0.0; 3], exact,
        ).is_err());
    }
    assert!(PreparedMultipoleRadiatedEnergy::new(
        &[MultipoleRadiationLine { operator: &owners[0], sources: &a[..2] }], [0.0; 3], exact,
    ).is_err());
    let mut invalid = a;
    invalid[0].quadrupole.stress_moment_n_m[4].im = f64::INFINITY;
    assert!(PreparedMultipoleRadiatedEnergy::new(
        &[MultipoleRadiationLine { operator: &owners[0], sources: &invalid }], [0.0; 3], exact,
    ).is_err());
    let plane = IdealAcousticPlane::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], IdealPlaneBoundary::Rigid).unwrap();
    let observer = PreparedCompactSourceObserver::new(medium(), 137.0, [0.0, 0.0, 2.0], &positions, 3).unwrap().with_plane(plane, 6).unwrap();
    let images = PreparedCompactRadiatedPower::from_half_space(&observer, power_budget(6)).unwrap().with_quadrupoles().unwrap();
    assert!(PreparedMultipoleRadiatedEnergy::new(
        &[MultipoleRadiationLine { operator: &images, sources: &a }], [0.0; 3], exact,
    ).is_err());
    let empty_owner = owner(137.0, &[]);
    let dark = PreparedMultipoleRadiatedEnergy::new(
        &[MultipoleRadiationLine { operator: &empty_owner, sources: &[] }], [0.0; 3], budget(1, 0),
    ).unwrap().observe(window).unwrap();
    assert_eq!(dark.radiated_energy_j, 0.0);
    assert_eq!(dark.energy_first_moment_j, [0.0; 3]);
}
