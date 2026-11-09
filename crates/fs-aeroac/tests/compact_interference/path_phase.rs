//! G1/G3: retain short path differences on many-wavelength acoustic paths.

use super::*;
use fs_aeroac::observation::{HarmonicObservationBudget, PreparedHarmonicField};
use fs_math::dd::Dd;

fn axial_stress(amplitude: C64) -> QuadrupolePhasor {
    QuadrupolePhasor {
        stress_moment_n_m: [
            amplitude,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
            C64::ZERO,
        ],
    }
}

// Factor exp(ikr)/(4*pi*r) out before comparing the sources. The remaining
// pressure derivatives are local rational polynomials, so this reference
// never subtracts two rounded many-wavelength carrier phases or ranges.
fn axial_envelopes(r: f64, k: f64) -> [C64; 3] {
    [
        C64::ONE,
        C64::new(1.0 / r, -k),
        C64::new(2.0 / (r * r) - k * k, -2.0 * k / r),
    ]
}

fn check_path_difference(rho: f64, c: f64, frequency: f64, range: f64, separation: f64) {
    let acoustic = medium(rho, c);
    let positions = [[0.0; 3], [separation, 0.0, 0.0]];
    let receiver = [range, 0.0, 0.0];
    let compact =
        PreparedCompactSourceObserver::new(acoustic, frequency, receiver, &positions, 2).unwrap();
    let stress =
        PreparedQuadrupoleObserver::new(acoustic, frequency, receiver, &positions, 2).unwrap();
    let k = compact.wavenumber_per_m().re;
    let first = axial_envelopes(range, k);
    let second = axial_envelopes(range - separation, k);
    let relative_phase = phase(-k * separation);
    let inactive = CompactSourcePhasor::default();
    let source = |component, amplitude| {
        if component == 0 {
            CompactSourcePhasor::monopole(amplitude)
        } else {
            CompactSourcePhasor::force_dipole([amplitude, C64::ZERO, C64::ZERO])
        }
    };
    for component in 0..3 {
        let at = |index| {
            if component == 2 {
                let mut sources = [axial_stress(C64::ZERO); 2];
                sources[index] = axial_stress(C64::ONE);
                stress.evaluate(&sources).unwrap()
            } else {
                let mut sources = [inactive; 2];
                sources[index] = source(component, C64::ONE);
                compact.evaluate(&sources).unwrap()
            }
        };
        let left = at(0);
        let right = at(1);
        let expected = relative_phase
            * (second[component] / first[component]).scale(range / (range - separation));
        let actual = right.pressure_pa / left.pressure_pa;
        assert!(
            (actual - expected).abs() < 3e-15,
            "component={component}, range={range:e}, actual={actual:?}, expected={expected:?}"
        );

        // Independent source control: compensate geometric spreading and the
        // derivative envelope, then apply the small source-path phase once.
        let cancelling_gain = -(C64::ONE / expected);
        let null = if component == 2 {
            stress
                .evaluate(&[axial_stress(C64::ONE), axial_stress(cancelling_gain)])
                .unwrap()
        } else {
            compact
                .evaluate(&[
                    source(component, C64::ONE),
                    source(component, cancelling_gain),
                ])
                .unwrap()
        };
        assert!(null.pressure_pa.abs() < 4e-15 * left.pressure_pa.abs());
    }
}

#[test]
fn g1_ultrasonic_path_phase_preserves_a_designed_spatial_null() {
    // A 1 MHz wave in a declared water-like lossless medium, 100 m from a
    // submillimetre source pair. This checks the propagation model, not a
    // physical no-absorption claim for real water or source calibration.
    check_path_difference(1000.0, 1500.0, 1.0e6, 100.0, 0.0003);
}

#[test]
fn g3_long_path_retains_source_separation_below_the_radius_ulp() {
    // An arithmetic limiting case: binary64 x-y loses the half-metre source
    // separation entirely. All three source mechanisms must retain its phase.
    check_path_difference(1.2, 343.0, 343.0, 2.0_f64.powi(54), 0.5);
}

#[test]
fn g1_diagonal_paths_retain_the_pythagorean_path_difference() {
    let range = 1e5;
    let separation = 1e-4;
    let compact = PreparedCompactSourceObserver::new(
        medium(1000.0, 1500.0),
        1e6,
        [3.0 * range, 4.0 * range, 0.0],
        &[[0.0; 3], [3.0 * separation, 4.0 * separation, 0.0]],
        2,
    )
    .unwrap();
    let source = CompactSourcePhasor::monopole(C64::ONE);
    let silent = CompactSourcePhasor::default();
    let left = compact.evaluate(&[source, silent]).unwrap();
    let right = compact.evaluate(&[silent, source]).unwrap();
    let ratio = phase(-5.0 * compact.wavenumber_per_m().re * separation)
        .scale(range / (range - separation));
    assert!((right.pressure_pa / left.pressure_pa - ratio).abs() < 3e-15);
}

#[test]
fn g1_mixed_volume_and_stress_sources_keep_the_same_short_path_clock() {
    let rho = 1000.0;
    let frequency = 1.0e6;
    let acoustic = medium(rho, 1500.0);
    let range = 100.0;
    let separation = 0.0003;
    let receiver = [range, 0.0, 0.0];
    let compact =
        PreparedCompactSourceObserver::new(acoustic, frequency, receiver, &[[0.0; 3]], 1).unwrap();
    let stress = PreparedQuadrupoleObserver::new(
        acoustic,
        frequency,
        receiver,
        &[[separation, 0.0, 0.0]],
        1,
    )
    .unwrap();
    let omega = TAU * frequency;
    let k = compact.wavenumber_per_m().re;
    let stress_gain = C64::new(0.0, rho * omega) / axial_envelopes(range - separation, k)[2]
        * phase(k * separation).scale((range - separation) / range);
    let volume_line = compact
        .evaluate_harmonic_line(&[CompactSourcePhasor::monopole(C64::ONE)])
        .unwrap();
    let stress_line = stress
        .evaluate_harmonic_line(&[axial_stress(stress_gain)])
        .unwrap();
    let magnitude = volume_line.pressure_pa.abs();
    let field = PreparedHarmonicField::new(
        &[volume_line, stress_line],
        HarmonicObservationBudget {
            max_lines: 2,
            max_line_pairs: 4,
        },
    )
    .unwrap();
    for time in [0.0, 0.13 / frequency, 0.37 / frequency] {
        let sample = field.field_at_clock_dd(Dd::from_f64(time)).unwrap();
        assert!(sample.pressure_pa.abs() < 5e-15 * magnitude);
    }
    // A pressure cancellation does not erase the different near-field
    // gradients. The original velocity readout survives the same composition.
    assert!(
        field
            .field_at_clock_dd(Dd::ZERO)
            .unwrap()
            .particle_velocity_m_s[0]
            .abs()
            > 0.0
    );
}
