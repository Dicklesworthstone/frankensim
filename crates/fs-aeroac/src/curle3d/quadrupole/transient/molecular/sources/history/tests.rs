use super::*;
use crate::curle3d::quadrupole::transient::history::PreparedQuadrupoleHistory;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

const DT: f64 = 1.0 / 4096.0;
const SAMPLES: usize = 64;
const SCALE: f64 = 0.01;
const COMPACT: [[f64; 3]; 1] = [[0.0, 0.0, 0.0]];
const STRESS: [[f64; 3]; 1] = [[0.1, -0.03, 0.02]];
const Q: f64 = 0.001;
const FORCE: [f64; 3] = [0.2, -0.4, 0.3];
const TENSOR: [f64; 6] = [0.1, -0.04, 0.2, 0.03, -0.02, 0.01];

fn medium(active: bool) -> RelaxationState {
    // The same admitted analytic caloric model as the original source tests.
    let gas = GasState::try_new(
        &GasSpec {
            gamma: 1.25,
            ..GasSpec::dry_air_ussa1976()
        },
        293.15,
        101325.0,
    )
    .unwrap();
    let species = SpeciesId::new("joint-history-caloric-oracle").unwrap();
    let processes = if active {
        vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: gas.pressure * 0.0008,
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
fn acoustic() -> FluidAcousticState {
    medium(false).gas().acoustic_state().unwrap()
}
fn budget(samples: usize) -> AcousticSourceHistoryBudget {
    AcousticSourceHistoryBudget {
        maximum_compact_sources: 1,
        maximum_quadrupole_sources: 1,
        maximum_samples_per_source: samples,
        maximum_component_samples: 10 * samples,
        maximum_workspace_bytes: 1 << 20,
    }
}
fn clock(epoch: f64) -> MolecularSampleClock {
    MolecularSampleClock::new(epoch, DT, 0.0).unwrap()
}
fn config() -> MolecularTransientConfig {
    MolecularTransientConfig {
        timestep_s: DT,
        samples: SAMPLES,
        transform_len: 4 * SAMPLES,
        contour_exponent: 24.0,
    }
}
fn polynomial(t: f64) -> (f64, f64, f64, f64) {
    let x = t / SCALE;
    (
        x * x * x,
        3.0 * x * x / SCALE,
        6.0 * x / (SCALE * SCALE),
        SCALE * x * x * x * x / 4.0,
    )
}
fn ends(impulse: [f64; 3], integral: [f64; 6]) -> (CompactHistoryEnds, QuadrupoleHistoryEnds) {
    let rate = polynomial((SAMPLES - 1) as f64 * DT).1;
    (
        CompactHistoryEnds {
            initial_volume_flow_rate_m3_s2: 0.0,
            final_volume_flow_rate_m3_s2: Q * rate,
            initial_force_rate_n_s: [0.0; 3],
            final_force_rate_n_s: FORCE.map(|x| x * rate),
            initial_force_impulse_n_s: impulse,
        },
        QuadrupoleHistoryEnds {
            initial_rate_n_m_s: [0.0; 6],
            final_rate_n_m_s: TENSOR.map(|x| x * rate),
            initial_integral_n_m_s: integral,
        },
    )
}
fn recorder(epoch: f64) -> AcousticSourceHistoryRecorder {
    let mut record = AcousticSourceHistoryRecorder::new(
        acoustic(),
        &COMPACT,
        &STRESS,
        clock(epoch),
        SAMPLES,
        budget(SAMPLES),
        || false,
    )
    .unwrap();
    for i in 0..SAMPLES {
        let v = polynomial(i as f64 * DT).0;
        record
            .push_sample(
                &COMPACT,
                &STRESS,
                clock(epoch).time_at(i).unwrap(),
                &[[Q * v, FORCE[0] * v, FORCE[1] * v, FORCE[2] * v]],
                &[TENSOR.map(|x| x * v)],
                || false,
            )
            .unwrap();
    }
    record
}
fn close(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * scale.abs().max(1e-18),
        "{actual:e} != {expected:e}, scale {scale:e}"
    );
}

#[test]
fn shared_cubic_preserves_physical_derivatives_integrals_and_original_stress_api() {
    let epoch = 2_f64.powi(40);
    let record = recorder(epoch);
    let impulse = [0.02, -0.03, 0.04];
    let integral = [0.1, -0.2, 0.3, 0.04, -0.05, 0.06];
    let (ce, se) = ends(impulse, integral);
    let history = record
        .prepare_clamped_cubic(&[ce], &[se], 100.0, || false)
        .unwrap();
    let original_stress = PreparedQuadrupoleHistory::new(
        &STRESS,
        clock(epoch),
        record.sample_times_s(),
        record.stress_samples_n_m(),
        &[se],
        QuadrupoleHistoryBudget {
            maximum_sources: 1,
            maximum_samples_per_source: SAMPLES,
            maximum_workspace_bytes: 1 << 20,
        },
        || false,
    )
    .unwrap();
    assert_eq!(
        history.initial_source(0).unwrap().force_impulse_n_s,
        impulse
    );
    assert_eq!(
        history
            .initial_source(1)
            .unwrap()
            .stress
            .stress_integral_n_m_s,
        integral
    );
    for i in 1..(SAMPLES - 1) * 3 {
        let elapsed = DT * i as f64 / 3.0;
        let time = Dd::from_f64(epoch) + Dd::from_f64(elapsed);
        let (value, rate, acceleration, primitive) = polynomial(elapsed);
        let compact = history.source_at_clock_dd(0, time).unwrap();
        close(compact.compact.volume_flow_m3_s, Q * value, Q * value);
        close(compact.compact.volume_flow_rate_m3_s2, Q * rate, Q * rate);
        for a in 0..3 {
            close(
                compact.compact.force_on_fluid_n[a],
                FORCE[a] * value,
                FORCE[a] * value,
            );
            close(
                compact.compact.force_rate_on_fluid_n_s[a],
                FORCE[a] * rate,
                FORCE[a] * rate,
            );
            close(
                compact.force_impulse_n_s[a],
                impulse[a] + FORCE[a] * primitive,
                impulse[a].abs() + (FORCE[a] * primitive).abs(),
            );
        }
        let stress = history.source_at_clock_dd(1, time).unwrap().stress;
        assert_eq!(stress, original_stress.source_at_clock_dd(0, time).unwrap());
        for a in 0..6 {
            close(
                stress.stress_moment_n_m[a],
                TENSOR[a] * value,
                TENSOR[a] * value,
            );
            close(
                stress.stress_rate_n_m_s[a],
                TENSOR[a] * rate,
                TENSOR[a] * rate,
            );
            close(
                stress.stress_acceleration_n_m_s2[a],
                TENSOR[a] * acceleration,
                TENSOR[a] * acceleration,
            );
            close(
                stress.stress_integral_n_m_s[a],
                integral[a] + TENSOR[a] * primitive,
                integral[a].abs() + (TENSOR[a] * primitive).abs(),
            );
        }
    }
    assert!(
        history
            .source_at_clock_dd(0, Dd::from_f64(epoch) - Dd::from_f64(DT))
            .is_err()
    );
    assert!(history.source_at_clock_dd(2, Dd::from_f64(epoch)).is_err());
    let observer = history
        .prepare_molecular_observer(
            medium(true),
            [0.7, -0.4, 0.3],
            config(),
            MolecularTransientBudget::default(),
            || false,
        )
        .unwrap();
    // A supplied prehistory impulse is retained; it cannot masquerade as rest.
    assert!(
        history
            .observe_initially_resting(&observer, || false)
            .is_err()
    );
}

#[test]
fn mixed_history_reaches_the_original_causal_molecular_complete_field() {
    let epoch = 2_f64.powi(40);
    let record = recorder(epoch);
    let (ce, se) = ends([0.0; 3], [0.0; 6]);
    let history = record
        .prepare_clamped_cubic(&[ce], &[se], 100.0, || false)
        .unwrap();
    let mut lossless = None;
    for active in [false, true] {
        let observer = history
            .prepare_molecular_observer(
                medium(active),
                [0.7, -0.4, 0.3],
                config(),
                MolecularTransientBudget::default(),
                || false,
            )
            .unwrap();
        let actual = history
            .observe_initially_resting(&observer, || false)
            .unwrap();
        // Independent analytic cubic values/rates/integrals; no spline queries.
        let expected = observer
            .observe_initially_resting(
                clock(epoch),
                |source, time| {
                    let t = (time - Dd::from_f64(epoch)).to_f64();
                    let (v, d, dd, integral) = polynomial(t);
                    let mut result = TransientAcousticSource::default();
                    if source == 0 {
                        result.compact.volume_flow_m3_s = Q * v;
                        result.compact.volume_flow_rate_m3_s2 = Q * d;
                        result.compact.force_on_fluid_n = FORCE.map(|x| x * v);
                        result.compact.force_rate_on_fluid_n_s = FORCE.map(|x| x * d);
                        result.force_impulse_n_s = FORCE.map(|x| x * integral);
                    } else {
                        result.stress.stress_moment_n_m = TENSOR.map(|x| x * v);
                        result.stress.stress_rate_n_m_s = TENSOR.map(|x| x * d);
                        result.stress.stress_acceleration_n_m_s2 = TENSOR.map(|x| x * dd);
                        result.stress.stress_integral_n_m_s = TENSOR.map(|x| x * integral);
                    }
                    Ok(result)
                },
                || false,
            )
            .unwrap();
        let pscale = expected
            .fields()
            .iter()
            .map(|f| f.pressure_pa.abs())
            .fold(0.0, f64::max);
        let uscale = expected
            .fields()
            .iter()
            .flat_map(|f| f.particle_velocity_m_s)
            .map(f64::abs)
            .fold(0.0, f64::max);
        assert!(pscale > 1e-4 && uscale > 1e-7);
        let first_front = (0..observer.paths().path_count())
            .map(|i| observer.paths().path_delay_s(i).unwrap())
            .fold(f64::INFINITY, f64::min);
        for (i, (a, e)) in actual.fields().iter().zip(expected.fields()).enumerate() {
            if i as f64 * DT < first_front {
                assert_eq!(a.pressure_pa, 0.0);
                assert_eq!(a.particle_velocity_m_s, [0.0; 3]);
            }
            close(a.pressure_pa, e.pressure_pa, pscale);
            for axis in 0..3 {
                close(
                    a.particle_velocity_m_s[axis],
                    e.particle_velocity_m_s[axis],
                    uscale,
                );
                close(
                    a.instantaneous_intensity_w_m2[axis],
                    e.instantaneous_intensity_w_m2[axis],
                    pscale * uscale,
                );
            }
        }
        if let Some(previous) = &lossless {
            let difference = actual
                .fields()
                .iter()
                .zip(previous)
                .map(|(field, pressure): (_, &f64)| (field.pressure_pa - pressure).abs())
                .fold(0.0, f64::max);
            assert!(
                difference > 1e-4 * pscale,
                "molecular material must change the physical history"
            );
        } else {
            lossless = Some(
                actual
                    .fields()
                    .iter()
                    .map(|f| f.pressure_pa)
                    .collect::<Vec<_>>(),
            );
        }
    }
}

#[test]
fn lossless_limit_has_the_analytic_monopole_pressure_and_radial_velocity() {
    let epoch = 0.25;
    let mut record = AcousticSourceHistoryRecorder::new(
        acoustic(),
        &COMPACT,
        &STRESS,
        clock(epoch),
        SAMPLES,
        budget(SAMPLES),
        || false,
    )
    .unwrap();
    for i in 0..SAMPLES {
        record
            .push_sample(
                &COMPACT,
                &STRESS,
                clock(epoch).time_at(i).unwrap(),
                &[[Q * polynomial(i as f64 * DT).0, 0.0, 0.0, 0.0]],
                &[[0.0; 6]],
                || false,
            )
            .unwrap();
    }
    let (mut ce, mut se) = ends([0.0; 3], [0.0; 6]);
    ce.final_force_rate_n_s = [0.0; 3];
    se.final_rate_n_m_s = [0.0; 6];
    let history = record
        .prepare_clamped_cubic(&[ce], &[se], 100.0, || false)
        .unwrap();
    let r = 0.7;
    let observer = history
        .prepare_molecular_observer(
            medium(false),
            [r, 0.0, 0.0],
            config(),
            MolecularTransientBudget::default(),
            || false,
        )
        .unwrap();
    let actual = history
        .observe_initially_resting(&observer, || false)
        .unwrap();
    let rho = acoustic().density_kg_m3();
    let c = acoustic().sound_speed_m_s();
    for (i, field) in actual.fields().iter().enumerate() {
        let t = i as f64 * DT - r / c;
        let (value, rate) = if t < 0.0 {
            (0.0, 0.0)
        } else {
            let (v, d, _, _) = polynomial(t);
            (Q * v, Q * d)
        };
        let p = rho * rate / (4.0 * core::f64::consts::PI * r);
        let ux = (value / (r * r) + rate / (c * r)) / (4.0 * core::f64::consts::PI);
        assert!((field.pressure_pa - p).abs() < 2e-8 * p.abs().max(1e-15));
        assert!((field.particle_velocity_m_s[0] - ux).abs() < 2e-8 * ux.abs().max(1e-15));
        assert_eq!(field.particle_velocity_m_s[1], 0.0);
        assert_eq!(field.particle_velocity_m_s[2], 0.0);
        assert!(
            (field.instantaneous_intensity_w_m2[0] - p * ux).abs()
                < 4e-8 * (p * ux).abs().max(1e-15)
        );
    }
}

#[test]
fn original_prefix_clock_geometry_budgets_and_preparation_are_retryable() {
    let make = |limits| {
        AcousticSourceHistoryRecorder::new(
            acoustic(),
            &COMPACT,
            &STRESS,
            clock(0.25),
            SAMPLES,
            limits,
            || false,
        )
    };
    let mut limits = budget(SAMPLES);
    limits.maximum_component_samples -= 1;
    assert!(make(limits).is_err());
    let required = make(budget(SAMPLES)).unwrap().workspace_bytes();
    limits = budget(SAMPLES);
    limits.maximum_workspace_bytes = required - 1;
    assert!(make(limits).is_err());
    limits.maximum_workspace_bytes = required;
    let mut raw = make(limits).unwrap();
    let row = [[0.001, 0.2, -0.4, 0.3]];
    assert!(
        raw.push_sample(&COMPACT, &STRESS, 0.25 + DT, &row, &[TENSOR], || false)
            .is_err()
    );
    assert!(
        raw.push_sample(&COMPACT, &COMPACT, 0.25, &row, &[TENSOR], || false)
            .is_err()
    );
    let mut reference = make(limits).unwrap();
    let mut polls = 0;
    reference
        .push_sample(&COMPACT, &STRESS, 0.25, &row, &[TENSOR], || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(
        raw.push_sample(&COMPACT, &STRESS, 0.25, &row, &[TENSOR], || {
            calls += 1;
            calls == polls
        })
        .is_err()
    );
    assert_eq!(raw.accepted_sample_count(), 0);
    assert!(raw.volume_force_samples().is_empty() && raw.stress_samples_n_m().is_empty());
    raw.push_sample(&COMPACT, &STRESS, 0.25, &row, &[TENSOR], || false)
        .unwrap();
    assert_eq!(raw.volume_force_samples(), reference.volume_force_samples());
    assert_eq!(raw.stress_samples_n_m(), reference.stress_samples_n_m());
    let record = recorder(0.25);
    let (ce, se) = ends([0.0; 3], [0.0; 6]);
    assert!(
        record
            .prepare_clamped_cubic(&[ce], &[se], 0.5 / DT, || false)
            .is_err()
    );
    assert!(
        record
            .prepare_clamped_cubic(&[ce], &[se], 100.0, || true)
            .is_err()
    );
    let history = record
        .prepare_clamped_cubic(&[ce], &[se], 100.0, || false)
        .unwrap();
    assert_eq!(record.accepted_sample_count(), SAMPLES);
    let observer = PreparedMolecularSourceObserver::new(
        medium(true),
        [0.7, -0.4, 0.3],
        &[STRESS[0], COMPACT[0]],
        config(),
        MolecularTransientBudget::default(),
        || false,
    )
    .unwrap();
    assert!(
        history
            .observe_initially_resting(&observer, || false)
            .is_err()
    );
    let mut too_long = config();
    too_long.samples += 1;
    assert!(
        history
            .prepare_molecular_observer(
                medium(true),
                [0.7, -0.4, 0.3],
                too_long,
                MolecularTransientBudget::default(),
                || false
            )
            .is_err()
    );
    assert_eq!(
        history.initial_source(0).unwrap().compact.volume_flow_m3_s,
        0.0
    );
}
