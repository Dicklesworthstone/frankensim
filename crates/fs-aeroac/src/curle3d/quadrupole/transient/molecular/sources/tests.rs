//! G1/G3 complete source Green fields, causal refinement and coherent images;
//! G4 physical onset, numerical extent and cancellation admission.

use super::*;
use crate::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use crate::curle3d::{
    CompactForceObserver, CompactSourcePhasor, IdealAcousticPlane, IdealPlaneBoundary,
    PreparedCompactSourceObserver,
};
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;

const TAU: f64 = 0.0008;
const RATE: f64 = 800.0;
const DURATION: f64 = 0.04;
const POSITION: [f64; 3] = [0.7, -0.4, 0.3];
const VOLUME: f64 = 0.0005;
const FORCE: [f64; 3] = [0.3, -0.7, 0.2];
const STRESS: [f64; 6] = [0.7, -0.3, 0.9, 0.23, -0.11, 0.4];

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
    let species = SpeciesId::new("complete-source-caloric-oracle").unwrap();
    let processes = if active {
        vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: gas.pressure * TAU,
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

fn clock(config: MolecularTransientConfig, epoch: f64) -> MolecularSampleClock {
    MolecularSampleClock::new(epoch, config.timestep_s, 0.0).unwrap()
}

fn observer(active: bool, samples: usize) -> PreparedMolecularSourceObserver {
    PreparedMolecularSourceObserver::new(
        medium(active),
        POSITION,
        &[[0.0; 3]],
        config(samples),
        MolecularTransientBudget::default(),
        || false,
    )
    .unwrap()
}

fn values(field: &RetardedSourceField) -> [f64; 4] {
    let [x, y, z] = field.particle_velocity_m_s;
    [field.pressure_pa, x, y, z]
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.abs().max(1e-24),
        "{actual:e} != {expected:e}, scale={scale:e}, tolerance={tolerance:e}"
    );
}

fn select(mut source: TransientAcousticSource, mechanism: usize) -> TransientAcousticSource {
    if mechanism != 0 && mechanism != 3 {
        source.compact.volume_flow_m3_s = 0.0;
        source.compact.volume_flow_rate_m3_s2 = 0.0;
    }
    if mechanism != 1 && mechanism != 3 {
        source.compact.force_on_fluid_n = [0.0; 3];
        source.compact.force_rate_on_fluid_n_s = [0.0; 3];
        source.force_impulse_n_s = [0.0; 3];
    }
    if mechanism != 2 && mechanism != 3 {
        source.stress = TransientQuadrupoleSource::default();
    }
    source
}

fn scaled(mut source: TransientAcousticSource, factor: f64) -> TransientAcousticSource {
    source.compact.volume_flow_m3_s *= factor;
    source.compact.volume_flow_rate_m3_s2 *= factor;
    source.compact.force_on_fluid_n = source.compact.force_on_fluid_n.map(|v| v * factor);
    source.compact.force_rate_on_fluid_n_s =
        source.compact.force_rate_on_fluid_n_s.map(|v| v * factor);
    source.force_impulse_n_s = source.force_impulse_n_s.map(|v| v * factor);
    source.stress.stress_integral_n_m_s = source.stress.stress_integral_n_m_s.map(|v| v * factor);
    source.stress.stress_moment_n_m = source.stress.stress_moment_n_m.map(|v| v * factor);
    source.stress.stress_rate_n_m_s = source.stress.stress_rate_n_m_s.map(|v| v * factor);
    source.stress.stress_acceleration_n_m_s2 =
        source.stress.stress_acceleration_n_m_s2.map(|v| v * factor);
    source
}

// One causal analytic pulse, with matching derivatives/integrals for every
// mechanism. Its transform is RATE^6/(s+RATE)^7.
fn pulse(time: f64) -> TransientAcousticSource {
    if time <= 0.0 {
        return TransientAcousticSource::default();
    }
    let x = RATE * time;
    let exponential = (-x).exp();
    let p = x.powi(6) * exponential / 720.0;
    let dp = RATE * (6.0 * x.powi(5) - x.powi(6)) * exponential / 720.0;
    let ddp =
        RATE.powi(2) * (30.0 * x.powi(4) - 12.0 * x.powi(5) + x.powi(6)) * exponential / 720.0;
    let mut sum = 1.0;
    let mut term = 1.0;
    for order in 1..=6 {
        term *= x / order as f64;
        sum += term;
    }
    let integral = (1.0 - exponential * sum) / RATE;
    TransientAcousticSource {
        compact: CompactSourceStrength {
            volume_flow_m3_s: VOLUME * p,
            volume_flow_rate_m3_s2: VOLUME * dp,
            force_on_fluid_n: FORCE.map(|v| v * p),
            force_rate_on_fluid_n_s: FORCE.map(|v| v * dp),
        },
        force_impulse_n_s: FORCE.map(|v| v * integral),
        stress: TransientQuadrupoleSource {
            stress_integral_n_m_s: STRESS.map(|v| v * integral),
            stress_moment_n_m: STRESS.map(|v| v * p),
            stress_rate_n_m_s: STRESS.map(|v| v * dp),
            stress_acceleration_n_m_s2: STRESS.map(|v| v * ddp),
        },
    }
}

// Independently eliminated single reservoir and Cartesian Green derivatives.
// This oracle does not use production molecular responses, radial projections,
// or convolution-quadrature coefficients.
fn original_laplace(material: &RelaxationState, s: C64, mechanism: usize) -> [C64; 4] {
    let gas = material.gas();
    let cv0 = gas.specific_heat_cp - gas.specific_gas_constant;
    let cvi = cv0 - gas.specific_gas_constant;
    let cp0 = cv0 + gas.specific_gas_constant;
    let cpi = cvi + gas.specific_gas_constant;
    let bulk = ((C64::from_re(cp0) + s.scale(TAU * cpi))
        / (C64::from_re(cv0) + s.scale(TAU * cvi)))
    .scale(gas.pressure);
    let k = s * bulk.recip().scale(gas.density).sqrt();
    let range = POSITION.iter().map(|x| x * x).sum::<f64>().sqrt();
    let n = POSITION.map(|x| x / range);
    let e = k.scale(-range);
    let green =
        C64::new(e.im.cos(), e.im.sin()).scale(e.re.exp() / (4.0 * core::f64::consts::PI * range));
    let q = if mechanism == 0 || mechanism == 3 {
        VOLUME
    } else {
        0.0
    };
    let f = if mechanism == 1 || mechanism == 3 {
        FORCE
    } else {
        [0.0; 3]
    };
    let tensor = if mechanism == 2 || mechanism == 3 {
        STRESS
    } else {
        [0.0; 6]
    };
    let [xx, yy, zz, xy, xz, yz] = tensor;
    let vector = [
        xx * n[0] + xy * n[1] + xz * n[2],
        xy * n[0] + yy * n[1] + yz * n[2],
        xz * n[0] + yz * n[1] + zz * n[2],
    ];
    let radial = (0..3).map(|i| n[i] * vector[i]).sum::<f64>();
    let rf = (0..3).map(|i| n[i] * f[i]).sum::<f64>();
    let trace = xx + yy + zz;
    let first = k + C64::from_re(1.0 / range);
    let a = k * k + k.scale(3.0 / range) + C64::from_re(3.0 / range.powi(2));
    let b = k.scale(-1.0 / range) - C64::from_re(1.0 / range.powi(2));
    let third = (k * k * k).scale(-1.0)
        - (k * k).scale(6.0 / range)
        - k.scale(15.0 / range.powi(2))
        - C64::from_re(15.0 / range.powi(3));
    let mut field = [C64::ZERO; 4];
    field[0] =
        green * (s.scale(gas.density * q) + first.scale(rf) + a.scale(radial) + b.scale(trace));
    for axis in 0..3 {
        let force = a.scale(n[axis] * rf) + b.scale(f[axis]);
        let stress_gradient =
            third.scale(radial * n[axis]) + a.scale((trace * n[axis] + 2.0 * vector[axis]) / range);
        field[axis + 1] = green
            * (first.scale(q * n[axis]) + ((force - stress_gradient) / s).scale(1.0 / gas.density));
    }
    field
}

#[test]
fn original_volume_force_stress_green_fields_and_harmonic_convention() {
    let owner = observer(true, 65);
    let path = &owner.propagation.geometry.transfers[0];
    let material = owner.molecular_state();
    let c = material.frozen_sound_speed_m_s();
    for s in [
        C64::new(120.0, 700.0),
        C64::new(0.0, -1700.0),
        C64::new(7000.0, -3000.0),
    ] {
        let response = material.laplace_response(s).unwrap();
        let e = response
            .excess_propagation_rate_per_m
            .scale(-path.delay_s * c);
        let first_filter = C64::new(e.im.cos(), e.im.sin()).scale(e.re.exp());
        for mechanism in 0..4 {
            let mut actual = [C64::ZERO; 4];
            let mut filter = first_filter;
            for order in 0..4 {
                for (quantity, factor) in [s.recip(), C64::ONE, s, s * s].into_iter().enumerate() {
                    let mut source = TransientAcousticSource::default();
                    match quantity {
                        0 => {
                            source.force_impulse_n_s = FORCE;
                            source.stress.stress_integral_n_m_s = STRESS;
                        }
                        1 => {
                            source.compact.volume_flow_m3_s = VOLUME;
                            source.compact.force_on_fluid_n = FORCE;
                            source.stress.stress_moment_n_m = STRESS;
                        }
                        2 => {
                            source.compact.volume_flow_rate_m3_s2 = VOLUME;
                            source.compact.force_rate_on_fluid_n_s = FORCE;
                            source.stress.stress_rate_n_m_s = STRESS;
                        }
                        _ => source.stress.stress_acceleration_n_m_s2 = STRESS,
                    }
                    let input = owner.compact[0]
                        .inputs(path, select(source, mechanism))
                        .unwrap();
                    for field in 0..4 {
                        actual[field] = actual[field] + filter * factor.scale(input[order][field]);
                    }
                }
                filter = filter * response.slowness_s_m.scale(c);
            }
            let delay = s.scale(-path.delay_s);
            let delay = C64::new(delay.im.cos(), delay.im.sin()).scale(delay.re.exp());
            let expected = original_laplace(material, s, mechanism);
            for field in 0..4 {
                assert!(
                    (actual[field] * delay - expected[field]).abs()
                        < 4e-13 * expected[field].abs().max(1e-20)
                );
            }
        }
        if s.re == 0.0 {
            let frequency = -s.im / core::f64::consts::TAU;
            let compact = PreparedCompactSourceObserver::new_molecular(
                material.clone(),
                frequency,
                POSITION,
                &[[0.0; 3]],
                1,
            )
            .unwrap()
            .evaluate(&[CompactSourcePhasor {
                volume_flow_m3_s: C64::from_re(VOLUME),
                force_on_fluid_n: FORCE.map(C64::from_re),
            }])
            .unwrap();
            let stress = PreparedQuadrupoleObserver::new_molecular(
                material.clone(),
                frequency,
                POSITION,
                &[[0.0; 3]],
                1,
            )
            .unwrap()
            .evaluate(&[QuadrupolePhasor {
                stress_moment_n_m: STRESS.map(C64::from_re),
            }])
            .unwrap();
            let actual = [
                compact.pressure_pa + stress.pressure_pa,
                compact.particle_velocity_m_s[0] + stress.particle_velocity_m_s[0],
                compact.particle_velocity_m_s[1] + stress.particle_velocity_m_s[1],
                compact.particle_velocity_m_s[2] + stress.particle_velocity_m_s[2],
            ];
            let expected = original_laplace(material, s, 3);
            for field in 0..4 {
                assert!((actual[field] - expected[field]).abs() < 4e-13 * expected[field].abs());
            }
        }
    }
}

#[test]
fn lossless_complete_fields_recover_native_sources_and_original_stress_cache() {
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
    let force = CompactForceObserver::new(owner.paths().acoustic(), POSITION, [0.0; 3]).unwrap();
    let range = POSITION.iter().map(|x| x * x).sum::<f64>().sqrt();
    let n = POSITION.map(|x| x / range);
    let gas = owner.molecular_state().gas();
    let c = owner.molecular_state().frozen_sound_speed_m_s();
    for (i, actual) in record.fields().iter().enumerate() {
        let t = record.time_at_clock_dd(i).unwrap();
        let elapsed = t - Dd::from_f64(epoch);
        let source = pulse((elapsed - Dd::from_f64(force.delay_s())).to_f64());
        let dipole = force
            .field_at_retarded(
                source.force_impulse_n_s,
                source.compact.force_on_fluid_n,
                source.compact.force_rate_on_fluid_n_s,
            )
            .unwrap();
        let stress = owner
            .paths()
            .evaluate_at_clock_dd(t, |_, t| {
                Ok(pulse((t - Dd::from_f64(epoch)).to_f64()).stress)
            })
            .unwrap();
        let g = 1.0 / (4.0 * core::f64::consts::PI * range);
        let expected = [
            gas.density * source.compact.volume_flow_rate_m3_s2 * g
                + dipole.pressure_pa
                + stress.pressure_pa,
            0.0,
            0.0,
            0.0,
        ];
        close(
            actual.pressure_pa,
            expected[0],
            expected[0].abs().max(1e-12),
            3e-11,
        );
        for axis in 0..3 {
            let expected = n[axis]
                * g
                * (source.compact.volume_flow_m3_s / range
                    + source.compact.volume_flow_rate_m3_s2 / c)
                + dipole.particle_velocity_m_s[axis]
                + stress.particle_velocity_m_s[axis];
            close(
                actual.particle_velocity_m_s[axis],
                expected,
                expected.abs().max(1e-14),
                3e-11,
            );
            assert_eq!(
                actual.instantaneous_intensity_w_m2[axis],
                actual.pressure_pa * actual.particle_velocity_m_s[axis]
            );
        }
        if elapsed.to_f64() < force.delay_s() {
            assert_eq!(values(actual), [0.0; 4]);
        }
    }
    let original = owner
        .stress_observer()
        .observe_initially_resting(
            clock,
            |_, t| Ok(pulse((t - Dd::from_f64(epoch)).to_f64()).stress),
            || false,
        )
        .unwrap();
    let only_stress = owner
        .observe_initially_resting(
            clock,
            |_, t| Ok(select(pulse((t - Dd::from_f64(epoch)).to_f64()), 2)),
            || false,
        )
        .unwrap();
    for (original, composed) in original.fields().iter().zip(only_stress.fields()) {
        assert_eq!(values(original), values(composed));
    }
}

// Direct Bromwich inversion of the independent original material/Green law.
// Frequency extent is separately doubled. At dw=20, the alias period greatly
// exceeds the compared interval and exp(-sigma*period) is below 3e-14.
fn bromwich(material: &RelaxationState, times: &[f64], extent: usize) -> Vec<[f64; 4]> {
    let sigma = 100.0;
    let dw = 20.0;
    let mut output = vec![[0.0; 4]; times.len()];
    for bin in 0..=extent {
        let omega = bin as f64 * dw;
        let s = C64::new(sigma, omega);
        let mut denominator = C64::ONE;
        for _ in 0..7 {
            denominator = denominator * (s + C64::from_re(RATE));
        }
        let source = denominator.recip().scale(RATE.powi(6));
        let field = original_laplace(material, s, 3).map(|v| v * source);
        let weight = if bin == 0 || bin == extent { 0.5 } else { 1.0 };
        for (&time, row) in times.iter().zip(&mut output) {
            let phase = C64::new((omega * time).cos(), (omega * time).sin());
            for axis in 0..4 {
                row[axis] += weight * (field[axis] * phase).re;
            }
        }
    }
    for (&time, row) in times.iter().zip(&mut output) {
        for value in row {
            *value *= dw * (sigma * time).exp() / core::f64::consts::PI;
        }
    }
    output
}

#[test]
fn complete_molecular_time_fields_refine_against_independent_bromwich_solution() {
    let material = medium(true);
    let times = [0.01, 0.015, 0.02, 0.025, 0.03];
    let reference = bromwich(&material, &times, 20000);
    let lower_extent = bromwich(&material, &times, 10000);
    let impedance = material.gas().density * material.frozen_sound_speed_m_s();
    let scale = reference.iter().map(|f| f[0].abs()).fold(0.0, f64::max);
    for (a, b) in reference.iter().zip(lower_extent) {
        for field in 0..4 {
            close(
                a[field],
                b[field],
                scale / if field == 0 { 1.0 } else { impedance },
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
            let i = (time / owner.config().timestep_s).round() as usize;
            for (field, actual) in values(&result.fields()[i]).into_iter().enumerate() {
                error = error.max(
                    (actual - expected[field]).abs() * if field == 0 { 1.0 } else { impedance },
                );
            }
        }
        errors.push(error / scale);
    }
    assert!(errors[0] > 1e-7, "nonvacuous molecular field: {errors:?}");
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.4 * pair[0],
            "complete-field BDF2 refinement: {errors:?}"
        );
    }
    assert!(errors[2] < 2e-4, "resolved complete field: {errors:?}");
}

#[test]
fn complete_sources_obey_ideal_boundaries_and_coherent_interference() {
    let material = medium(true);
    let frozen =
        FluidAcousticState::try_new(material.gas().density, material.frozen_sound_speed_m_s())
            .unwrap();
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], boundary).unwrap();
        let paths = PreparedTransientQuadrupoleObserver::new(
            frozen,
            [0.6, -0.3, 0.0],
            &[[0.1, 0.2, 0.3]],
            1,
        )
        .unwrap()
        .with_orthogonal_planes(&[plane], 2)
        .unwrap();
        let owner = PreparedMolecularSourceObserver::from_paths(
            material.clone(),
            paths,
            config(129),
            MolecularTransientBudget::default(),
            || false,
        )
        .unwrap();
        let result = owner
            .observe_initially_resting(
                clock(owner.config(), 0.0),
                |_, t| Ok(pulse(t.to_f64())),
                || false,
            )
            .unwrap();
        let z = material.gas().density * material.frozen_sound_speed_m_s();
        let peak = result
            .fields()
            .iter()
            .map(|f| f.pressure_pa.abs() + z * f.particle_velocity_m_s[2].abs())
            .fold(0.0, f64::max);
        assert!(peak > 1e-5);
        for field in result.fields() {
            let residual = match boundary {
                IdealPlaneBoundary::Rigid => z * field.particle_velocity_m_s[2],
                IdealPlaneBoundary::PressureRelease => field.pressure_pa,
            };
            close(residual, 0.0, peak, 3e-12);
            close(
                field.instantaneous_intensity_w_m2[2],
                0.0,
                peak * peak / z,
                3e-12,
            );
        }
    }
    let owner = observer(true, 129);
    let records = (0..4)
        .map(|mechanism| {
            owner
                .observe_initially_resting(
                    clock(owner.config(), 0.0),
                    |_, t| Ok(select(pulse(t.to_f64()), mechanism)),
                    || false,
                )
                .unwrap()
        })
        .collect::<Vec<_>>();
    let peak = records[3]
        .fields()
        .iter()
        .map(|f| f.pressure_pa.abs())
        .fold(0.0, f64::max);
    let z = material.gas().density * material.frozen_sound_speed_m_s();
    let mut cross_term: f64 = 0.0;
    for i in 0..owner.config().samples {
        let total = &records[3].fields()[i];
        for field in 0..4 {
            let expected = (0..3)
                .map(|j| values(&records[j].fields()[i])[field])
                .sum::<f64>();
            close(
                values(total)[field],
                expected,
                peak / if field == 0 { 1.0 } else { z },
                3e-11,
            );
        }
        for axis in 0..3 {
            let separate = (0..3)
                .map(|j| records[j].fields()[i].instantaneous_intensity_w_m2[axis])
                .sum::<f64>();
            cross_term =
                cross_term.max((total.instantaneous_intensity_w_m2[axis] - separate).abs());
            assert_eq!(
                total.instantaneous_intensity_w_m2[axis],
                total.pressure_pa * total.particle_velocity_m_s[axis]
            );
        }
    }
    assert!(cross_term > 1e-5 * peak * peak / z);
    let opposite = PreparedMolecularSourceObserver::new(
        material,
        POSITION,
        &[[0.0; 3]; 2],
        config(129),
        MolecularTransientBudget::default(),
        || false,
    )
    .unwrap();
    let result = opposite
        .observe_initially_resting(
            clock(opposite.config(), 0.0),
            |index, t| {
                Ok(scaled(
                    pulse(t.to_f64()),
                    if index == 0 { 1.0 } else { -1.0 },
                ))
            },
            || false,
        )
        .unwrap();
    for field in result.fields() {
        assert_eq!(values(field), [0.0; 4]);
    }
}

#[test]
fn source_clock_contour_refinement_and_future_pulse_causality() {
    let owner = observer(true, 257);
    let clock = clock(owner.config(), 1e9);
    let epoch = Dd::from_f64(clock.epoch_s());
    let base = owner
        .observe_initially_resting(
            clock,
            |index, t| {
                assert_eq!(index, 0);
                let elapsed = (t - epoch).to_f64();
                assert!((0.0..=DURATION).contains(&elapsed));
                Ok(pulse(elapsed))
            },
            || false,
        )
        .unwrap();
    let larger = PreparedMolecularSourceObserver::new(
        medium(true),
        POSITION,
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
        .observe_initially_resting(clock, |_, t| Ok(pulse((t - epoch).to_f64())), || false)
        .unwrap();
    let changed = owner
        .observe_initially_resting(
            clock,
            |_, t| {
                let time = (t - epoch).to_f64();
                // A later physical source pulse with the same derivatives/integrals.
                let a = pulse(time);
                let b = scaled(pulse(time - 0.02), 10.0);
                Ok(TransientAcousticSource {
                    compact: CompactSourceStrength {
                        volume_flow_m3_s: a.compact.volume_flow_m3_s + b.compact.volume_flow_m3_s,
                        volume_flow_rate_m3_s2: a.compact.volume_flow_rate_m3_s2
                            + b.compact.volume_flow_rate_m3_s2,
                        force_on_fluid_n: core::array::from_fn(|i| {
                            a.compact.force_on_fluid_n[i] + b.compact.force_on_fluid_n[i]
                        }),
                        force_rate_on_fluid_n_s: core::array::from_fn(|i| {
                            a.compact.force_rate_on_fluid_n_s[i]
                                + b.compact.force_rate_on_fluid_n_s[i]
                        }),
                    },
                    force_impulse_n_s: core::array::from_fn(|i| {
                        a.force_impulse_n_s[i] + b.force_impulse_n_s[i]
                    }),
                    stress: TransientQuadrupoleSource {
                        stress_integral_n_m_s: core::array::from_fn(|i| {
                            a.stress.stress_integral_n_m_s[i] + b.stress.stress_integral_n_m_s[i]
                        }),
                        stress_moment_n_m: core::array::from_fn(|i| {
                            a.stress.stress_moment_n_m[i] + b.stress.stress_moment_n_m[i]
                        }),
                        stress_rate_n_m_s: core::array::from_fn(|i| {
                            a.stress.stress_rate_n_m_s[i] + b.stress.stress_rate_n_m_s[i]
                        }),
                        stress_acceleration_n_m_s2: core::array::from_fn(|i| {
                            a.stress.stress_acceleration_n_m_s2[i]
                                + b.stress.stress_acceleration_n_m_s2[i]
                        }),
                    },
                })
            },
            || false,
        )
        .unwrap();
    let peak = base
        .fields()
        .iter()
        .map(|f| f.pressure_pa.abs())
        .fold(0.0, f64::max);
    let z =
        owner.molecular_state().gas().density * owner.molecular_state().frozen_sound_speed_m_s();
    for i in 0..owner.config().samples {
        let time = i as f64 * clock.timestep_s();
        if time < owner.paths().path_delay_s(0).unwrap() {
            assert_eq!(values(&base.fields()[i]), [0.0; 4]);
        }
        for field in 0..4 {
            let scale = peak / if field == 0 { 1.0 } else { z };
            close(
                values(&base.fields()[i])[field],
                values(&refined.fields()[i])[field],
                scale,
                3e-8,
            );
            if time < 0.02 + owner.paths().path_delay_s(0).unwrap() {
                close(
                    values(&base.fields()[i])[field],
                    values(&changed.fields()[i])[field],
                    scale,
                    4e-11,
                );
            }
        }
    }
    assert!(
        base.fields()
            .iter()
            .zip(changed.fields())
            .any(|(a, b)| (a.pressure_pa - b.pressure_pa).abs() > peak)
    );
}

#[test]
fn mechanism_specific_initial_rest_storage_limits_and_cancellation() {
    let owner = observer(true, 65);
    let clock = clock(owner.config(), 0.0);
    for kind in 0..7 {
        assert!(
            owner
                .observe_initially_resting(
                    clock,
                    |_, t| {
                        let mut source = pulse(t.to_f64());
                        match kind {
                            0 => source.compact.volume_flow_m3_s = 1e-30,
                            1 => source.compact.force_on_fluid_n[1] = 1e-30,
                            2 => source.force_impulse_n_s[2] = 1e-30,
                            3 => source.stress.stress_moment_n_m[0] = 1e-30,
                            4 => source.stress.stress_rate_n_m_s[0] = 1e-30,
                            5 => source.stress.stress_integral_n_m_s[0] = 1e-30,
                            _ => source.compact.volume_flow_rate_m3_s2 = f64::NAN,
                        }
                        Ok(source)
                    },
                    || false
                )
                .is_err()
        );
    }
    let jump = owner
        .observe_initially_resting(
            clock,
            |_, t| {
                let t = t.to_f64();
                Ok(TransientAcousticSource {
                    compact: CompactSourceStrength {
                        volume_flow_m3_s: VOLUME * t,
                        volume_flow_rate_m3_s2: VOLUME,
                        force_on_fluid_n: FORCE.map(|v| v * t),
                        force_rate_on_fluid_n_s: FORCE,
                    },
                    force_impulse_n_s: FORCE.map(|v| v * t * t / 2.0),
                    ..TransientAcousticSource::default()
                })
            },
            || false,
        )
        .unwrap();
    assert!(jump.fields().iter().any(|f| f.pressure_pa != 0.0));
    assert_eq!(
        owner.workspace_bytes() - owner.stress_observer().workspace_bytes(),
        owner.paths().path_count() * core::mem::size_of::<CompactTransfer>()
    );
    for bytes in [
        owner.workspace_bytes() - 1,
        owner.stress_observer().workspace_bytes(),
    ] {
        assert!(
            PreparedMolecularSourceObserver::new(
                medium(true),
                POSITION,
                &[[0.0; 3]],
                owner.config(),
                MolecularTransientBudget {
                    maximum_workspace_bytes: bytes,
                    ..MolecularTransientBudget::default()
                },
                || false
            )
            .is_err()
        );
    }
    assert!(
        PreparedMolecularSourceObserver::new(
            medium(true),
            POSITION,
            &[[0.0; 3]],
            owner.config(),
            MolecularTransientBudget {
                maximum_workspace_bytes: owner.workspace_bytes(),
                ..MolecularTransientBudget::default()
            },
            || false
        )
        .is_ok()
    );
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
                MolecularSampleClock::new(0.0, 2.0 * clock.timestep_s(), 0.0).unwrap(),
                |_, t| Ok(pulse(t.to_f64())),
                || false
            )
            .is_err()
    );
    assert!(
        owner
            .observe_initially_resting(clock, |_, t| Ok(pulse(t.to_f64())), || false)
            .is_ok()
    );
    let empty = PreparedMolecularSourceObserver::new(
        medium(true),
        POSITION,
        &[],
        config(65),
        MolecularTransientBudget::default(),
        || false,
    )
    .unwrap();
    let silence = empty
        .observe_initially_resting(clock, |_, _| panic!("empty source queried"), || false)
        .unwrap();
    for field in silence.fields() {
        assert_eq!(values(field), [0.0; 4]);
    }
}
