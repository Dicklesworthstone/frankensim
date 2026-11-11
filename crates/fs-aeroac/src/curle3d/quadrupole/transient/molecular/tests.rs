use super::super::{
    TransientQuadrupoleSource,
    history::{PreparedQuadrupoleHistory, QuadrupoleHistoryBudget, QuadrupoleHistoryEnds},
};
use super::*;
use crate::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use crate::curle3d::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

const RELAXATION_TIME: f64 = 0.0008;
const DURATION: f64 = 0.04;
const RATE: f64 = 800.0;
const TENSOR: [f64; 6] = [0.7, -0.3, 0.9, 0.23, -0.11, 0.4];
const RECEIVER: [f64; 3] = [0.7, -0.4, 0.3];

fn medium(active: bool) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            gamma: 1.25,
            ..GasSpec::dry_air_ussa1976()
        },
        293.15,
        101325.0,
    )
    .unwrap();
    let species = SpeciesId::new("transient-caloric-oracle").unwrap();
    let processes = if active {
        vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: gas.pressure * RELAXATION_TIME,
                t_ref_k: gas.temperature,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 500.0,
            },
        }]
    } else {
        vec![]
    };
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes,
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
    .unwrap()
}

fn config(samples: usize) -> MolecularTransientConfig {
    MolecularTransientConfig {
        timestep_s: DURATION / (samples - 1) as f64,
        samples,
        transform_len: 4 * samples.next_power_of_two(),
        contour_exponent: 24.0,
    }
}

fn frozen(material: &RelaxationState) -> FluidAcousticState {
    FluidAcousticState::try_new(material.gas().density, material.frozen_sound_speed_m_s()).unwrap()
}

fn observer(active: bool, samples: usize) -> PreparedMolecularQuadrupoleObserver {
    PreparedMolecularQuadrupoleObserver::new(
        medium(active),
        RECEIVER,
        &[[0.0; 3]],
        config(samples),
        MolecularTransientBudget::default(),
        || false,
    )
    .unwrap()
}

fn clock(config: MolecularTransientConfig, epoch: f64) -> QuadrupoleSampleClock {
    QuadrupoleSampleClock::new(epoch, config.timestep_s, 0.0).unwrap()
}

// Smooth, causal manufactured stress, with independently analytic derivatives
// and integral. Its transform is RATE^6/(s+RATE)^7, used by the separate oracle.
fn pulse(time: f64) -> TransientQuadrupoleSource {
    if time <= 0.0 {
        return TransientQuadrupoleSource::default();
    }
    let x = RATE * time;
    let exponential = (-x).exp();
    let value = x.powi(6) * exponential / 720.0;
    let first = RATE * (6.0 * x.powi(5) - x.powi(6)) * exponential / 720.0;
    let second =
        RATE * RATE * (30.0 * x.powi(4) - 12.0 * x.powi(5) + x.powi(6)) * exponential / 720.0;
    let mut sum = 1.0;
    let mut term = 1.0;
    for order in 1..=6 {
        term *= x / order as f64;
        sum += term;
    }
    let integral = (1.0 - exponential * sum) / RATE;
    TransientQuadrupoleSource {
        stress_moment_n_m: TENSOR.map(|q| q * value),
        stress_rate_n_m_s: TENSOR.map(|q| q * first),
        stress_acceleration_n_m_s2: TENSOR.map(|q| q * second),
        stress_integral_n_m_s: TENSOR.map(|q| q * integral),
    }
}

fn field_values(field: &RetardedQuadrupoleField) -> [f64; 4] {
    [
        field.pressure_pa,
        field.particle_velocity_m_s[0],
        field.particle_velocity_m_s[1],
        field.particle_velocity_m_s[2],
    ]
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.abs().max(1e-25),
        "{actual:e} != {expected:e}, scale {scale:e}, tolerance {tolerance:e}"
    );
}

// Independent single-reservoir elimination and direct Cartesian Green
// derivatives. No production Laplace/CQ symbol or transient helper is used.
fn original_laplace(material: &RelaxationState, s: C64, x: [f64; 3], q: [f64; 6]) -> [C64; 4] {
    let gas = material.gas();
    let cv0 = gas.specific_heat_cp - gas.specific_gas_constant;
    let cvi = cv0 - gas.specific_gas_constant;
    let cp0 = cv0 + gas.specific_gas_constant;
    let cpi = cvi + gas.specific_gas_constant;
    let bulk = ((C64::from_re(cp0) + s.scale(RELAXATION_TIME * cpi))
        / (C64::from_re(cv0) + s.scale(RELAXATION_TIME * cvi)))
    .scale(gas.pressure);
    let kappa = s * bulk.recip().scale(gas.density).sqrt();
    let range = x.iter().map(|v| v * v).sum::<f64>().sqrt();
    let n = x.map(|v| v / range);
    let phase = -kappa.im * range;
    let green = C64::new(phase.cos(), phase.sin())
        .scale((-kappa.re * range).exp() / (4.0 * core::f64::consts::PI * range));
    let [xx, yy, zz, xy, xz, yz] = q;
    let vector = [
        xx * n[0] + xy * n[1] + xz * n[2],
        xy * n[0] + yy * n[1] + yz * n[2],
        xz * n[0] + yz * n[1] + zz * n[2],
    ];
    let radial = (0..3).map(|i| n[i] * vector[i]).sum::<f64>();
    let trace = xx + yy + zz;
    let a = kappa * kappa + kappa.scale(3.0 / range) + C64::from_re(3.0 / range.powi(2));
    let b = kappa.scale(-1.0 / range) - C64::from_re(1.0 / range.powi(2));
    let c = (kappa * kappa * kappa).scale(-1.0)
        - (kappa * kappa).scale(6.0 / range)
        - kappa.scale(15.0 / range.powi(2))
        - C64::from_re(15.0 / range.powi(3));
    let mut field = [C64::ZERO; 4];
    field[0] = green * (a.scale(radial) + b.scale(trace));
    for axis in 0..3 {
        let gradient = green
            * (c.scale(radial * n[axis]) + a.scale((trace * n[axis] + 2.0 * vector[axis]) / range));
        field[axis + 1] = (gradient / s).scale(-1.0 / gas.density);
    }
    field
}

#[test]
fn all_four_radial_filters_match_original_green_derivatives_and_harmonic_convention() {
    let material = medium(true);
    let paths =
        PreparedTransientQuadrupoleObserver::new(frozen(&material), RECEIVER, &[[0.0; 3]], 1)
            .unwrap();
    let path = &paths.transfers[0];
    for s in [
        C64::new(100.0, 700.0),
        C64::new(0.0, -1700.0),
        C64::new(7000.0, -3000.0),
    ] {
        let response = material.laplace_response(s).unwrap();
        let range = path.delay_s * material.frozen_sound_speed_m_s();
        let e = response.excess_propagation_rate_per_m.scale(-range);
        let mut filter = C64::new(e.im.cos(), e.im.sin()).scale(e.re.exp());
        let slowness = response
            .slowness_s_m
            .scale(material.frozen_sound_speed_m_s());
        let mut actual = [C64::ZERO; 4];
        let factors = [s.recip(), C64::ONE, s, s * s];
        for order in 0..4 {
            for (quantity, factor) in factors.into_iter().enumerate() {
                let mut source = TransientQuadrupoleSource::default();
                match quantity {
                    0 => source.stress_integral_n_m_s = TENSOR,
                    1 => source.stress_moment_n_m = TENSOR,
                    2 => source.stress_rate_n_m_s = TENSOR,
                    _ => source.stress_acceleration_n_m_s2 = TENSOR,
                }
                let projected = observe::inputs(path, source).unwrap();
                for field in 0..4 {
                    actual[field] = actual[field] + filter * factor.scale(projected[order][field]);
                }
            }
            filter = filter * slowness;
        }
        let delay = s.scale(-path.delay_s);
        let delay = C64::new(delay.im.cos(), delay.im.sin()).scale(delay.re.exp());
        let expected = original_laplace(&material, s, RECEIVER, TENSOR);
        for field in 0..4 {
            assert!(
                (actual[field] * delay - expected[field]).abs()
                    < 3e-13 * expected[field].abs().max(1e-20)
            );
        }
        if s.re == 0.0 {
            let harmonic = PreparedQuadrupoleObserver::new_molecular(
                material.clone(),
                -s.im / core::f64::consts::TAU,
                RECEIVER,
                &[[0.0; 3]],
                1,
            )
            .unwrap()
            .evaluate(&[QuadrupolePhasor {
                stress_moment_n_m: TENSOR.map(C64::from_re),
            }])
            .unwrap();
            let harmonic = [
                harmonic.pressure_pa,
                harmonic.particle_velocity_m_s[0],
                harmonic.particle_velocity_m_s[1],
                harmonic.particle_velocity_m_s[2],
            ];
            for field in 0..4 {
                assert!((harmonic[field] - expected[field]).abs() < 3e-13 * expected[field].abs());
            }
        }
    }
}

#[test]
fn empty_reservoir_full_field_matches_retarded_observer_on_large_original_clock() {
    let owner = observer(false, 129);
    let epoch = 1e9;
    let clock = clock(owner.config(), epoch);
    let record = owner
        .observe_initially_resting(
            clock,
            |_, t| Ok(pulse((t - Dd::from_f64(epoch)).to_f64())),
            || false,
        )
        .unwrap();
    for (index, actual) in record.fields().iter().enumerate() {
        let time = record.time_at_clock_dd(index).unwrap();
        let expected = owner
            .paths()
            .evaluate_at_clock_dd(time, |_, t| Ok(pulse((t - Dd::from_f64(epoch)).to_f64())))
            .unwrap();
        let scale = field_values(&expected);
        for ((actual, expected), scale) in field_values(actual).into_iter().zip(scale).zip(scale) {
            close(actual, expected, scale.abs().max(1e-12), 2e-12);
        }
        if (time - Dd::from_f64(epoch)).to_f64() < owner.paths().path_delay_s(0).unwrap() {
            assert_eq!(field_values(actual), [0.0; 4]);
        }
        for axis in 0..3 {
            assert_eq!(
                actual.instantaneous_intensity_w_m2[axis],
                actual.pressure_pa * actual.particle_velocity_m_s[axis]
            );
        }
    }
}

// Direct Bromwich inversion, independent of CQ and source interpolation. The
// smooth source has an analytic transform; the separate original Green law
// includes the physical delay and all spatial derivatives. Frequency extent
// is independently refined; dw=20 gives a 2*pi/dw period far beyond the
// compared window, with exp(-sigma*period) below 3e-14 for its periodic aliases.
fn bromwich(material: &RelaxationState, times: &[f64], extent: usize) -> Vec<[f64; 4]> {
    let sigma = 100.0;
    let dw = 20.0;
    let mut result = vec![[0.0; 4]; times.len()];
    for bin in 0..=extent {
        let omega = bin as f64 * dw;
        let s = C64::new(sigma, omega);
        let denominator = s + C64::from_re(RATE);
        let mut power = C64::ONE;
        for _ in 0..7 {
            power = power * denominator;
        }
        let source = power.recip().scale(RATE.powi(6));
        let field = original_laplace(material, s, RECEIVER, TENSOR).map(|h| h * source);
        let weight = if bin == 0 || bin == extent { 0.5 } else { 1.0 };
        for (&time, output) in times.iter().zip(&mut result) {
            let phase = C64::new((omega * time).cos(), (omega * time).sin());
            for axis in 0..4 {
                output[axis] += weight * (field[axis] * phase).re;
            }
        }
    }
    for (&time, output) in times.iter().zip(&mut result) {
        for value in output {
            *value *= dw * (sigma * time).exp() / core::f64::consts::PI;
        }
    }
    result
}

#[test]
fn bdf2_refines_full_pressure_and_velocity_to_independent_molecular_time_solution() {
    let material = medium(true);
    let times = [0.01, 0.015, 0.02, 0.025, 0.03];
    let reference = bromwich(&material, &times, 20000);
    let lower_extent = bromwich(&material, &times, 10000);
    let impedance = material.gas().density * material.frozen_sound_speed_m_s();
    let scale = reference.iter().map(|f| f[0].abs()).fold(0.0, f64::max);
    for (a, b) in reference.iter().zip(lower_extent) {
        for axis in 0..4 {
            close(
                a[axis],
                b[axis],
                scale / if axis == 0 { 1.0 } else { impedance },
                2e-7,
            );
        }
    }
    let mut errors = Vec::new();
    for samples in [257, 513, 1025] {
        let owner = observer(true, samples);
        let result = owner
            .observe_initially_resting(
                clock(owner.config(), 0.0),
                |_, t| Ok(pulse(t.to_f64())),
                || false,
            )
            .unwrap();
        let mut error: f64 = 0.0;
        for (&time, expected) in times.iter().zip(&reference) {
            let index = (time / owner.config().timestep_s).round() as usize;
            let actual = field_values(&result.fields()[index]);
            for axis in 0..4 {
                error = error.max(
                    (actual[axis] - expected[axis]).abs() * if axis == 0 { 1.0 } else { impedance },
                );
            }
        }
        errors.push(error / scale);
    }
    assert!(
        errors[0] > 1e-5,
        "nonvacuous molecular refinement: {errors:?}"
    );
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.32 * pair[0],
            "BDF2 full-field refinement: {errors:?}"
        );
    }
    assert!(errors[2] < 0.002, "resolved molecular field: {errors:?}");
}

#[test]
fn ideal_images_retain_coherence_and_full_rigid_or_release_boundary_condition() {
    let material = medium(true);
    let source = [0.1, 0.2, 0.3];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], boundary).unwrap();
        let paths = PreparedTransientQuadrupoleObserver::new(
            frozen(&material),
            [0.6, -0.3, 0.0],
            &[source],
            1,
        )
        .unwrap()
        .with_orthogonal_planes(&[plane], 2)
        .unwrap();
        let owner = PreparedMolecularQuadrupoleObserver::from_paths(
            material.clone(),
            paths,
            config(129),
            MolecularTransientBudget::default(),
            || false,
        )
        .unwrap();
        assert_eq!(owner.paths().path_count(), 2);
        let result = owner
            .observe_initially_resting(
                clock(owner.config(), 0.0),
                |_, t| Ok(pulse(t.to_f64())),
                || false,
            )
            .unwrap();
        let peak = result
            .fields()
            .iter()
            .map(|f| {
                f.pressure_pa.abs()
                    + material.gas().density
                        * material.frozen_sound_speed_m_s()
                        * f.particle_velocity_m_s[2].abs()
            })
            .fold(0.0, f64::max);
        assert!(peak > 1e-5);
        for field in result.fields() {
            let value = match boundary {
                IdealPlaneBoundary::Rigid => {
                    field.particle_velocity_m_s[2]
                        * material.gas().density
                        * material.frozen_sound_speed_m_s()
                }
                IdealPlaneBoundary::PressureRelease => field.pressure_pa,
            };
            close(value, 0.0, peak, 2e-12);
            close(
                field.instantaneous_intensity_w_m2[2],
                0.0,
                peak * peak / (material.gas().density * material.frozen_sound_speed_m_s()),
                2e-12,
            );
        }
    }
    let owner = PreparedMolecularQuadrupoleObserver::new(
        material,
        RECEIVER,
        &[[0.0; 3]; 2],
        config(129),
        MolecularTransientBudget::default(),
        || false,
    )
    .unwrap();
    let record = owner
        .observe_initially_resting(
            clock(owner.config(), 0.0),
            |source, t| {
                let mut state = pulse(t.to_f64());
                if source == 1 {
                    state.stress_moment_n_m = state.stress_moment_n_m.map(|v| -v);
                    state.stress_rate_n_m_s = state.stress_rate_n_m_s.map(|v| -v);
                    state.stress_acceleration_n_m_s2 = state.stress_acceleration_n_m_s2.map(|v| -v);
                    state.stress_integral_n_m_s = state.stress_integral_n_m_s.map(|v| -v);
                }
                Ok(state)
            },
            || false,
        )
        .unwrap();
    for field in record.fields() {
        assert_eq!(field_values(field), [0.0; 4]);
        assert_eq!(field.instantaneous_intensity_w_m2, [0.0; 3]);
    }
}

#[test]
fn original_record_binding_contour_refinement_and_no_circular_future_wrap() {
    let owner = observer(true, 257);
    let clock = clock(owner.config(), 0.0);
    let times = (0..owner.config().samples)
        .map(|i| clock.time_at(i).unwrap())
        .collect::<Vec<_>>();
    let moments = times
        .iter()
        .map(|&t| pulse(t).stress_moment_n_m)
        .collect::<Vec<_>>();
    let history = PreparedQuadrupoleHistory::new(
        &[[0.0; 3]],
        clock,
        &times,
        &moments,
        &[QuadrupoleHistoryEnds {
            initial_rate_n_m_s: [0.0; 6],
            final_rate_n_m_s: pulse(DURATION).stress_rate_n_m_s,
            initial_integral_n_m_s: [0.0; 6],
        }],
        QuadrupoleHistoryBudget {
            maximum_sources: 1,
            maximum_samples_per_source: 1025,
            maximum_workspace_bytes: 1 << 20,
        },
        || false,
    )
    .unwrap();
    let recorded = owner.observe_history(&history, || false).unwrap();
    let analytic = owner
        .observe_initially_resting(clock, |_, t| Ok(pulse(t.to_f64())), || false)
        .unwrap();
    let peak = analytic
        .fields()
        .iter()
        .map(|f| f.pressure_pa.abs())
        .fold(0.0, f64::max);
    let error = recorded
        .fields()
        .iter()
        .zip(analytic.fields())
        .map(|(a, b)| (a.pressure_pa - b.pressure_pa).abs())
        .fold(0.0, f64::max);
    assert!(
        error < 0.002 * peak,
        "recorded interpolation {error:e}, peak {peak:e}"
    );
    let larger = PreparedMolecularQuadrupoleObserver::new(
        medium(true),
        RECEIVER,
        &[[0.0; 3]],
        MolecularTransientConfig {
            transform_len: 2 * owner.config().transform_len,
            ..owner.config()
        },
        MolecularTransientBudget::default(),
        || false,
    )
    .unwrap();
    let refined = larger
        .observe_initially_resting(clock, |_, t| Ok(pulse(t.to_f64())), || false)
        .unwrap();
    let changed = owner
        .observe_initially_resting(
            clock,
            |_, t| {
                let mut a = pulse(t.to_f64());
                let b = pulse(t.to_f64() - 0.02);
                for i in 0..6 {
                    a.stress_moment_n_m[i] += 10.0 * b.stress_moment_n_m[i];
                    a.stress_rate_n_m_s[i] += 10.0 * b.stress_rate_n_m_s[i];
                    a.stress_acceleration_n_m_s2[i] += 10.0 * b.stress_acceleration_n_m_s2[i];
                    a.stress_integral_n_m_s[i] += 10.0 * b.stress_integral_n_m_s[i];
                }
                Ok(a)
            },
            || false,
        )
        .unwrap();
    for (index, ((base, refined), changed)) in analytic
        .fields()
        .iter()
        .zip(refined.fields())
        .zip(changed.fields())
        .enumerate()
    {
        close(base.pressure_pa, refined.pressure_pa, peak, 3e-8);
        if index as f64 * clock.timestep_s() < 0.02 + owner.paths().path_delay_s(0).unwrap() {
            close(base.pressure_pa, changed.pressure_pa, peak, 3e-11);
        }
    }
    assert!(
        changed
            .fields()
            .iter()
            .zip(analytic.fields())
            .any(|(a, b)| (a.pressure_pa - b.pressure_pa).abs() > peak)
    );
}

#[test]
fn onset_front_ownership_budget_and_cancellation_are_explicit() {
    let owner = observer(true, 65);
    let clock = clock(owner.config(), 0.0);
    // A finite acceleration jump is physical; no delta onset terms are missing.
    let jump = owner
        .observe_initially_resting(
            clock,
            |_, t| {
                let t = t.to_f64();
                Ok(TransientQuadrupoleSource {
                    stress_moment_n_m: TENSOR.map(|q| q * t * t / 2.0),
                    stress_rate_n_m_s: TENSOR.map(|q| q * t),
                    stress_acceleration_n_m_s2: TENSOR,
                    stress_integral_n_m_s: TENSOR.map(|q| q * t * t * t / 6.0),
                })
            },
            || false,
        )
        .unwrap();
    assert!(jump.fields().iter().any(|f| f.pressure_pa != 0.0));
    for (i, f) in jump.fields().iter().enumerate() {
        if i as f64 * clock.timestep_s() < owner.paths().path_delay_s(0).unwrap() {
            assert_eq!(field_values(f), [0.0; 4]);
        }
    }
    for kind in 0..4 {
        assert!(
            owner
                .observe_initially_resting(
                    clock,
                    |_, t| {
                        let mut q = pulse(t.to_f64());
                        match kind {
                            0 => q.stress_moment_n_m[0] = 1e-30,
                            1 => q.stress_rate_n_m_s[0] = 1e-30,
                            2 => q.stress_integral_n_m_s[0] = 1e-30,
                            _ => q.stress_acceleration_n_m_s2[0] = f64::NAN,
                        };
                        Ok(q)
                    },
                    || false
                )
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        owner.observe_initially_resting(
            clock,
            |_, t| Ok(pulse(t.to_f64())),
            || {
                polls += 1;
                polls == 3
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(
        owner
            .observe_initially_resting(
                QuadrupoleSampleClock::new(0.0, 2.0 * clock.timestep_s(), 0.0).unwrap(),
                |_, t| Ok(pulse(t.to_f64())),
                || false
            )
            .is_err()
    );
    for budget in [
        MolecularTransientBudget {
            maximum_workspace_bytes: 1,
            ..MolecularTransientBudget::default()
        },
        MolecularTransientBudget {
            maximum_transform_work: 1,
            ..MolecularTransientBudget::default()
        },
        MolecularTransientBudget {
            maximum_paths: 0,
            ..MolecularTransientBudget::default()
        },
    ] {
        assert!(
            PreparedMolecularQuadrupoleObserver::new(
                medium(true),
                RECEIVER,
                &[[0.0; 3]],
                config(65),
                budget,
                || false
            )
            .is_err()
        );
    }
    let wrong = PreparedTransientQuadrupoleObserver::new(
        FluidAcousticState::try_new(1.2, 343.0).unwrap(),
        RECEIVER,
        &[[0.0; 3]],
        1,
    )
    .unwrap();
    assert!(
        PreparedMolecularQuadrupoleObserver::from_paths(
            medium(true),
            wrong,
            config(65),
            MolecularTransientBudget::default(),
            || false
        )
        .is_err()
    );
}
