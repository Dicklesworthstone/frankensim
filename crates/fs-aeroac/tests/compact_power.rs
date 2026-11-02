//! Analytic compact-array power against independent angular and existing-field flux oracles.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::{
    CompactSourcePhasor, PreparedCompactRadiatedPower, PreparedCompactSourceObserver,
    RadiatedPowerBudget, RadiatedPowerModeBudget,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{RelaxationAbsorption, RelaxingMedium, SpeciesRelaxation};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;
use std::f64::consts::{PI, TAU};

fn budget(n: usize) -> RadiatedPowerBudget {
    RadiatedPowerBudget { max_sources: n, max_source_pairs: n * n }
}

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.3, 331.0).unwrap()
}

fn phase(angle: f64) -> C64 {
    let (sin, cos) = angle.sin_cos();
    C64::new(cos, sin)
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 5e-11 * expected.abs().max(1e-12),
        "actual={actual:e}, expected={expected:e}");
}

fn scale(source: CompactSourcePhasor, factor: C64) -> CompactSourcePhasor {
    CompactSourcePhasor {
        volume_flow_m3_s: source.volume_flow_m3_s * factor,
        force_on_fluid_n: source.force_on_fluid_n.map(|value| value * factor),
    }
}

fn isolated(acoustic: FluidAcousticState, frequency: f64, source: CompactSourcePhasor) -> f64 {
    let rho = acoustic.density_kg_m3();
    let c = acoustic.sound_speed_m_s();
    let omega2 = (TAU * frequency).powi(2);
    rho * omega2 * source.volume_flow_m3_s.norm_sq() / (8.0 * PI * c)
        + omega2 * source.force_on_fluid_n.iter().map(|value| value.norm_sq()).sum::<f64>()
            / (24.0 * PI * rho * c.powi(3))
}

fn fixture() -> ([[f64; 3]; 3], [CompactSourcePhasor; 3]) {
    let positions = [[-0.25, 0.125, 0.0], [0.375, -0.25, 0.125], [0.0, 0.25, -0.375]];
    let sources = [
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(4e-4, 2e-4),
            force_on_fluid_n: [C64::new(0.3, 0.1), C64::new(-0.2, 0.4), C64::new(0.1, -0.3)],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(-1e-4, 3e-4),
            force_on_fluid_n: [C64::new(-0.1, 0.5), C64::new(0.2, -0.3), C64::new(0.4, 0.1)],
        },
        CompactSourcePhasor {
            volume_flow_m3_s: C64::new(2e-4, -1e-4),
            force_on_fluid_n: [C64::new(0.2, -0.1), C64::new(-0.3, 0.2), C64::new(0.1, 0.3)],
        },
    ];
    (positions, sources)
}

fn legendre(n: usize, x: f64) -> (f64, f64) {
    let (mut previous, mut value) = (0.0, 1.0);
    for degree in 1..=n {
        let next = ((2 * degree - 1) as f64 * x * value - (degree - 1) as f64 * previous) / degree as f64;
        previous = value;
        value = next;
    }
    (value, previous)
}

// Gauss-Legendre in polar cosine, periodic trapezoid in azimuth. This
// directly integrates a caller's physical field, not the production moments.
fn sphere_integral(mut value: impl FnMut([f64; 3]) -> f64) -> f64 {
    let order = 32;
    let azimuths = 64;
    let mut integral = 0.0;
    for index in 0..order {
        let mut mu = (PI * (index as f64 + 0.75) / (order as f64 + 0.5)).cos();
        for _ in 0..20 {
            let (p, previous) = legendre(order, mu);
            let derivative = order as f64 * (mu * p - previous) / (mu * mu - 1.0);
            let change = p / derivative;
            mu -= change;
            if change.abs() < 1e-15 { break; }
        }
        let (p, previous) = legendre(order, mu);
        let derivative = order as f64 * (mu * p - previous) / (mu * mu - 1.0);
        let weight = 2.0 / ((1.0 - mu * mu) * derivative * derivative);
        let radial = (1.0 - mu * mu).sqrt();
        for azimuth in 0..azimuths {
            let angle = TAU * azimuth as f64 / azimuths as f64;
            integral += weight * (TAU / azimuths as f64)
                * value([radial * angle.cos(), radial * angle.sin(), mu]);
        }
    }
    integral
}

fn far_field_quadrature(
    acoustic: FluidAcousticState,
    frequency: f64,
    positions: &[[f64; 3]],
    sources: &[CompactSourcePhasor],
) -> f64 {
    let k = TAU * frequency / acoustic.sound_speed_m_s();
    let impedance = acoustic.characteristic_impedance_pa_s_m();
    let integral = sphere_integral(|direction| {
        let mut field = C64::ZERO;
        for (position, source) in positions.iter().zip(sources) {
            let dot = (0..3).map(|axis| direction[axis] * position[axis]).sum::<f64>();
            let dipole = (0..3).fold(C64::ZERO, |sum, axis| sum + source.force_on_fluid_n[axis].scale(direction[axis]));
            field = field + (source.volume_flow_m3_s.scale(impedance) + dipole) * phase(-k * dot);
        }
        field.norm_sq()
    });
    k * k * integral / (32.0 * PI * PI * impedance)
}

#[test]
fn single_source_si_power_and_original_observer_binding() {
    let acoustic = medium();
    let frequency = 413.0;
    let (_, sources) = fixture();
    let position = [[0.125, -0.25, 0.375]];
    let observer = PreparedCompactSourceObserver::new(acoustic, frequency, [0.0, 0.0, 2.0], &position, 1).unwrap();
    let prepared = PreparedCompactRadiatedPower::from_observer(&observer, budget(1)).unwrap();
    assert_eq!(prepared.source_count(), 1);
    assert_eq!(prepared.position_count(), 1);
    assert_eq!(prepared.frequency_hz(), frequency);
    assert_eq!(prepared.acoustic_state().density_kg_m3(), acoustic.density_kg_m3());
    assert_eq!(prepared.acoustic_state().sound_speed_m_s(), acoustic.sound_speed_m_s());
    for source in [
        sources[0],
        CompactSourcePhasor::monopole(sources[0].volume_flow_m3_s),
        CompactSourcePhasor::force_dipole(sources[0].force_on_fluid_n),
    ] {
        close(prepared.evaluate(&[source]).unwrap().power_w, isolated(acoustic, frequency, source));
    }
}

#[test]
fn separated_monopoles_and_quadrature_monopole_dipole_keep_cross_terms() {
    let acoustic = medium();
    let frequency = 180.0;
    let k = TAU * frequency / acoustic.sound_speed_m_s();
    let impedance = acoustic.characteristic_impedance_pa_s_m();
    let z = 1.3;
    let prepared = PreparedCompactRadiatedPower::new(acoustic, frequency, &[[0.0; 3], [0.0, 0.0, z / k]], budget(2)).unwrap();
    let monopole = CompactSourcePhasor::monopole(C64::from_re(1.0 / impedance));
    let unit_power = k * k / (8.0 * PI * impedance);
    let j0 = z.sin() / z;
    let j1 = z.sin() / (z * z) - z.cos() / z;
    for sign in [-1.0, 1.0] {
        close(prepared.evaluate(&[monopole, scale(monopole, C64::from_re(sign))]).unwrap().power_w,
            unit_power * (2.0 + 2.0 * sign * j0));
        let dipole = CompactSourcePhasor::force_dipole([C64::ZERO, C64::ZERO, C64::new(0.0, sign)]);
        // The force is at +z: +i gives constructive integrated Q/F interference.
        close(prepared.evaluate(&[monopole, dipole]).unwrap().power_w,
            unit_power * (4.0 / 3.0 + 2.0 * sign * j1));
    }
}

#[test]
fn total_power_matches_existing_signed_near_field_flux_on_enclosing_spheres() {
    let acoustic = medium();
    let frequency = 180.0;
    let (positions, sources) = fixture();
    let prepared = PreparedCompactRadiatedPower::new(acoustic, frequency, &positions, budget(3)).unwrap();
    let power = prepared.evaluate(&sources).unwrap().power_w;
    close(power, far_field_quadrature(acoustic, frequency, &positions, &sources));
    for radius in [1.4, 3.7] {
        let flux = sphere_integral(|normal| {
            let observer = PreparedCompactSourceObserver::new(acoustic, frequency,
                normal.map(|value| value * radius), &positions, 3).unwrap();
            let field = observer.evaluate(&sources).unwrap();
            radius * radius * (0..3).map(|axis| normal[axis] * field.active_intensity_w_m2[axis]).sum::<f64>()
        });
        close(flux, power);
    }
}

#[test]
fn exact_position_coalescing_preserves_a_weak_residual_and_declared_independence() {
    let acoustic = medium();
    let frequency = 180.0;
    let (_, sources) = fixture();
    let small = sources[0];
    let large = scale(small, C64::from_re(1e16));
    let prepared = PreparedCompactRadiatedPower::new(acoustic, frequency,
        &[[0.25, -0.125, 0.375]; 3], budget(3)).unwrap();
    assert_eq!(prepared.position_count(), 1);
    close(prepared.evaluate(&[large, small, scale(large, -C64::ONE)]).unwrap().power_w,
        isolated(acoustic, frequency, small));
    let pair = PreparedCompactRadiatedPower::new(acoustic, frequency, &[[0.0; 3]; 2], budget(2)).unwrap();
    let opposite = [small, scale(small, -C64::ONE)];
    assert_eq!(pair.evaluate(&opposite).unwrap().power_w, 0.0);
    close(pair.evaluate_independent_sources(&opposite).unwrap().power_w,
        2.0 * isolated(acoustic, frequency, small));
}

#[test]
fn weak_subwavelength_monopole_and_parallel_perpendicular_force_pairs_survive() {
    let acoustic = medium();
    let frequency = 180.0;
    let k = TAU * frequency / acoustic.sound_speed_m_s();
    let impedance = acoustic.characteristic_impedance_pa_s_m();
    let unit_power = k * k / (8.0 * PI * impedance);
    let cases = [
        (CompactSourcePhasor::monopole(C64::from_re(1.0 / impedance)), 3.0),
        (CompactSourcePhasor::force_dipole([C64::ZERO, C64::ZERO, C64::ONE]), 5.0),
        (CompactSourcePhasor::force_dipole([C64::ONE, C64::ZERO, C64::ZERO]), 15.0),
    ];
    for z in [2.0_f64.powi(-40), 2.0_f64.powi(-80)] {
        let prepared = PreparedCompactRadiatedPower::new(acoustic, frequency,
            &[[0.0; 3], [0.0, 0.0, z / k]], budget(2)).unwrap();
        assert_eq!(prepared.position_count(), 2);
        for (source, denominator) in cases {
            let actual = prepared.evaluate(&[source, scale(source, -C64::ONE)]).unwrap().power_w;
            let expected = unit_power * z * z / denominator;
            assert!((actual / expected - 1.0).abs() < 2e-11,
                "z={z:e}, actual={actual:e}, expected={expected:e}");
        }
    }
}

#[test]
fn partial_coherence_matches_angular_quadrature_and_unitary_factor_changes() {
    let acoustic = medium();
    let frequency = 180.0;
    let (positions, sources) = fixture();
    let prepared = PreparedCompactRadiatedPower::new(acoustic, frequency, &positions, budget(3)).unwrap();
    let second: [CompactSourcePhasor; 3] = core::array::from_fn(|index|
        scale(sources[index], phase(0.6 + 0.8 * index as f64).scale(0.7)));
    let mode_budget = RadiatedPowerModeBudget { max_modes: 2, max_source_pair_mode_products: 18 };
    let expected = far_field_quadrature(acoustic, frequency, &positions, &sources)
        + far_field_quadrature(acoustic, frequency, &positions, &second);
    let power = prepared.evaluate_coherence_modes(&[&sources, &second], mode_budget).unwrap().power_w;
    close(power, expected);
    let mix = |a: CompactSourcePhasor, b: CompactSourcePhasor| {
        let imaginary = C64::new(0.0, 1.0);
        scale(CompactSourcePhasor {
            volume_flow_m3_s: a.volume_flow_m3_s + b.volume_flow_m3_s * imaginary,
            force_on_fluid_n: core::array::from_fn(|axis| a.force_on_fluid_n[axis] + b.force_on_fluid_n[axis] * imaginary),
        }, C64::from_re(0.5_f64.sqrt()))
    };
    let mixed_a: [CompactSourcePhasor; 3] = core::array::from_fn(|index| mix(sources[index], second[index]));
    let mixed_b: [CompactSourcePhasor; 3] = core::array::from_fn(|index| mix(second[index], sources[index]));
    close(prepared.evaluate_coherence_modes(&[&mixed_a, &mixed_b], mode_budget).unwrap().power_w, power);

    let sparse: [[CompactSourcePhasor; 3]; 3] = core::array::from_fn(|row|
        core::array::from_fn(|column| if row == column { sources[row] } else { CompactSourcePhasor::default() }));
    let independent = prepared.evaluate_independent_sources(&sources).unwrap().power_w;
    close(independent, sources.iter().map(|source| isolated(acoustic, frequency, *source)).sum());
    close(prepared.evaluate_coherence_modes(&[&sparse[0], &sparse[1], &sparse[2]],
        RadiatedPowerModeBudget { max_modes: 3, max_source_pair_mode_products: 27 }).unwrap().power_w,
        independent);
}

#[test]
fn translation_rotation_and_common_emission_phase_preserve_total_power() {
    let acoustic = medium();
    let frequency = 180.0;
    let (positions, sources) = fixture();
    let expected = PreparedCompactRadiatedPower::new(acoustic, frequency, &positions, budget(3)).unwrap()
        .evaluate(&sources).unwrap().power_w;
    let shifted = positions.map(|x| [-x[1] + 0.5, x[0] + 1.0, x[2] - 0.25]);
    let rotated = sources.map(|source| scale(CompactSourcePhasor {
        force_on_fluid_n: [-source.force_on_fluid_n[1], source.force_on_fluid_n[0], source.force_on_fluid_n[2]],
        ..source
    }, phase(0.71)));
    let actual = PreparedCompactRadiatedPower::new(acoustic, frequency, &shifted, budget(3)).unwrap()
        .evaluate(&rotated).unwrap().power_w;
    close(actual, expected);
}

#[test]
fn preparation_and_application_refusals_publish_no_partial_power() {
    let acoustic = medium();
    let frequency = 180.0;
    let (positions, sources) = fixture();
    for insufficient in [
        RadiatedPowerBudget { max_sources: 2, ..budget(3) },
        RadiatedPowerBudget { max_source_pairs: 8, ..budget(3) },
    ] {
        assert!(matches!(PreparedCompactRadiatedPower::new(acoustic, frequency, &positions, insufficient),
            Err(AeroacError::InvalidParameter { .. })));
    }
    for invalid_frequency in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(PreparedCompactRadiatedPower::new(acoustic, invalid_frequency, &positions, budget(3)).is_err());
    }
    assert!(PreparedCompactRadiatedPower::new(acoustic, frequency, &[[f64::NAN, 0.0, 0.0]], budget(1)).is_err());
    let mut preparation_polls = 0;
    let prepared = PreparedCompactRadiatedPower::new_cancellable(acoustic, frequency, &positions, budget(3), || {
        preparation_polls += 1;
        false
    }).unwrap();
    for stop in [1, preparation_polls / 2, preparation_polls] {
        let mut count = 0;
        assert!(matches!(PreparedCompactRadiatedPower::new_cancellable(acoustic, frequency, &positions, budget(3), || {
            count += 1;
            count == stop
        }), Err(AeroacError::Cancelled { .. })));
        assert_eq!(count, stop);
    }
    assert!(prepared.evaluate(&sources[..2]).is_err());
    let mut bad = sources;
    bad[1].volume_flow_m3_s.re = f64::NAN;
    assert!(matches!(prepared.evaluate(&bad), Err(AeroacError::NonFinite { .. })));
    let mut application_polls = 0;
    let expected = prepared.evaluate_cancellable(&sources, || { application_polls += 1; false }).unwrap();
    for stop in [1, application_polls / 2, application_polls] {
        let mut count = 0;
        assert!(matches!(prepared.evaluate_cancellable(&sources, || { count += 1; count == stop }),
            Err(AeroacError::Cancelled { .. })));
        assert_eq!(prepared.evaluate(&sources).unwrap(), expected);
    }
    for insufficient in [
        RadiatedPowerModeBudget { max_modes: 1, max_source_pair_mode_products: 18 },
        RadiatedPowerModeBudget { max_modes: 2, max_source_pair_mode_products: 17 },
    ] {
        assert!(prepared.evaluate_coherence_modes(&[&sources, &sources], insufficient).is_err());
    }
    let zero = PreparedCompactRadiatedPower::new(acoustic, frequency, &[], budget(0)).unwrap();
    assert_eq!(zero.evaluate(&[]).unwrap().power_w, 0.0);
    assert_eq!(prepared.evaluate_coherence_modes(&[],
        RadiatedPowerModeBudget { max_modes: 0, max_source_pair_mode_products: 0 }).unwrap().power_w, 0.0);
}

#[test]
fn observer_adapter_retains_model_boundary_and_supports_preparation_cancellation() {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let species = SpeciesId::new("synthetic").unwrap();
    // Even the empty-reservoir molecular variant keeps its original owner;
    // the adapter must not silently replace that propagation model.
    let state = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(), processes: vec![],
    }]).unwrap().evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap()).unwrap();
    let molecular = PreparedCompactSourceObserver::new_molecular(state, 180.0, [0.0, 0.0, 2.0], &[[0.0; 3]], 1).unwrap();
    assert!(matches!(PreparedCompactRadiatedPower::from_observer(&molecular, budget(1)),
        Err(AeroacError::InvalidParameter { .. })));
    let observer = PreparedCompactSourceObserver::new(medium(), 180.0, [0.0, 0.0, 2.0], &[[0.0; 3]], 1).unwrap();
    assert!(matches!(PreparedCompactRadiatedPower::from_observer_cancellable(&observer, budget(1), || true),
        Err(AeroacError::Cancelled { .. })));
}
