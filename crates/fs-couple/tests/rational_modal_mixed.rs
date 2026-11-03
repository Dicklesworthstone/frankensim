//! Mixed-rank modal radiation keeps every original state and physical work.
use core::f64::consts::{PI, TAU};
use fs_couple::mm_line::continuous::rational::{
    ModalRationalLoad, RationalModalLine, RationalModalLineBudget,
    analysis::RationalModalAnalysisBudget,
};
use fs_duct::RadiusProfile;
use fs_exec::CancelGate;
use fs_la::eigen_complex::lu_complex;
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget};
use fs_vfit::mimo::{
    MimoModel,
    load::{MimoLoadLimits, PreparedMimoLoad},
};
use fs_vfit::model::PoleTerm;

const RADIUS: f64 = 0.025;
const LENGTH: f64 = 0.035;
fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}
fn reference() -> f64 {
    gas().density * gas().sound_speed / (PI * RADIUS * RADIUS)
}
fn budget() -> RationalModalLineBudget {
    RationalModalLineBudget {
        linear: LinearSplitBudget {
            max_modes: 124,
            max_omega_dt: 1.0,
        },
        exact: DenseExponentialBudget {
            max_dimension: 256,
            max_dense_bytes: 64 << 20,
            max_multiply_terms: 20_000_000_000,
        },
        relative_tolerance: 1e-8,
    }
}
fn model(rank: usize, pair: bool) -> MimoModel {
    let z = reference();
    let c = [0.8, 0.28, 0.28, 0.6];
    let mut poles = vec![PoleTerm::Real {
        pole: -3500.0,
        residue: 0.0,
    }];
    let mut residues = vec![c.iter().map(|v| C64::from_re(v * z * 3500.0)).collect()];
    if pair {
        poles.push(PoleTerm::Pair {
            pole: C64::new(-400.0, 12000.0),
            residue: C64::ZERO,
        });
        residues.push(
            c.iter()
                .map(|v| C64::new(v * z * 400.0, v * z * 400.0 * 400.0 / 12000.0))
                .collect(),
        );
    }
    MimoModel {
        ports: 2,
        poles,
        residues,
        d: match rank {
            0 => vec![0.0; 4],
            1 => vec![0.2 * z; 4],
            2 => c.iter().map(|v| 0.2 * z * v).collect(),
            _ => unreachable!(),
        },
        e: vec![0.0; 4],
        reciprocal: true,
    }
}
fn prepared(source: &MimoModel, dt: f64) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        source,
        dt,
        MimoLoadLimits {
            band_rad_s: [1.0, 4e6],
            sampled_passivity_points: 33,
            max_states: 32,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}
fn attach(
    load: PreparedMimoLoad,
) -> Result<RationalModalLine, fs_couple::mm_line::continuous::rational::RationalModalError> {
    RationalModalLine::new(
        &RadiusProfile::exponential(RADIUS, 0.0, LENGTH).unwrap(),
        &gas(),
        2,
        ModalRationalLoad::new(load, RADIUS, gas()).unwrap(),
        budget(),
        &CancelGate::new_clock_free(),
    )
}
fn line(source: &MimoModel, dt: f64) -> RationalModalLine {
    attach(prepared(source, dt)).unwrap()
}
fn close(a: f64, b: f64, scale: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * scale.max(1e-28),
        "{a:e} != {b:e}, scale={scale:e}"
    );
}
fn snapshot(line: &RationalModalLine) -> Vec<u64> {
    let mut x = line.state().to_vec();
    x.extend_from_slice(line.incoming());
    for p in line.port_accounts() {
        x.extend([
            p.boundary_work,
            p.incident_energy,
            p.outgoing_energy,
            p.outgoing_wave,
        ]);
    }
    for p in line.terminal_accounts() {
        x.extend([
            p.mean_duct_pressure_pa,
            p.mean_load_pressure_pa,
            p.mean_outward_flow_m3_s,
            p.duct_work_j,
            p.load_work_j,
        ]);
    }
    for z in line.load().pole_states() {
        x.extend([z.re, z.im]);
    }
    x.into_iter().map(f64::to_bits).collect()
}

#[test]
fn mixed_rank_exact_subdivision_preserves_original_states_means_and_signed_work() {
    let source = model(1, true);
    let dt = 1.0 / 24_000.0;
    let mut coarse = line(&source, dt);
    let mut half = line(&source, dt / 2.0);
    let mut quarter = line(&source, dt / 4.0);
    assert_eq!(coarse.state().len(), 12 + source.order() * 2);
    assert_eq!(coarse.load().source_model(), &source);
    let gate = CancelGate::new_clock_free();
    let mut mutual_live = false;
    let mut amplitudes = vec![0.0_f64; coarse.state().len()];
    let mut state_errors = amplitudes.clone();
    for step in 0..72 {
        let wave = if step < 18 {
            [0.0, 60.0 * det::sin(0.2 * (step as f64 + 0.5))]
        } else {
            [0.0; 2]
        };
        let accepted = coarse.push(&wave, &gate).unwrap();
        mutual_live |= coarse.incoming()[0].abs() > 1e-7;
        for (refined, divisions) in [(&mut half, 2), (&mut quarter, 4)] {
            let (mut work, mut reflected, mut incident) = (0.0, 0.0, 0.0);
            let mut mean = [0.0; 2];
            for _ in 0..divisions {
                let record = refined.push(&wave, &gate).unwrap();
                work += record.terminal_work_j;
                reflected += record.reflected_energy_j;
                incident += record.incident_energy_j;
                for k in 0..2 {
                    mean[k] += refined.incoming()[k] / divisions as f64;
                }
            }
            // Compare each original coordinate in the independent trajectory
            // infinity norm. Instantaneous ringdown zero crossings are not
            // a meaningful denominator for accumulated floating-point error.
            for i in 0..amplitudes.len() {
                amplitudes[i] = amplitudes[i]
                    .max(coarse.state()[i].abs())
                    .max(refined.state()[i].abs());
                state_errors[i] =
                    state_errors[i].max((coarse.state()[i] - refined.state()[i]).abs());
            }
            let scale = coarse.stored_duct_energy_j()
                + accepted.incident_energy_j
                + accepted.reflected_energy_j;
            close(work, accepted.terminal_work_j, scale, 4e-8);
            close(reflected, accepted.reflected_energy_j, scale, 4e-8);
            close(incident, accepted.incident_energy_j, scale, 4e-8);
            for k in 0..2 {
                close(mean[k], coarse.incoming()[k], 60.0, 5e-9);
            }
        }
        for p in coarse.terminal_accounts() {
            close(
                p.mean_duct_pressure_pa,
                p.mean_load_pressure_pa,
                100.0,
                2e-9,
            );
            close(
                p.duct_work_j,
                p.load_work_j,
                coarse.stored_duct_energy_j()
                    + accepted.incident_energy_j
                    + accepted.reflected_energy_j,
                2e-9,
            );
        }
    }
    for i in 0..amplitudes.len() {
        assert!(
            state_errors[i] <= 2e-9 * amplitudes[i],
            "coordinate {i}: error={}, trajectory amplitude={}",
            state_errors[i],
            amplitudes[i]
        );
    }
    eprintln!(
        "maximum original-coordinate trajectory relative error={:e}",
        state_errors
            .iter()
            .zip(&amplitudes)
            .map(|(error, amplitude)| if *amplitude == 0.0 {
                0.0
            } else {
                error / amplitude
            })
            .fold(0.0, f64::max)
    );
    assert!(mutual_live);
    assert_eq!(coarse.load().source_model(), &source);
}

// Independent full weak harmonic equation. At this test frequency the original
// complete Z(s) is invertible even though D is singular. This is an oracle only;
// the runtime itself never uses a frequency-local impedance inverse.
fn weak_reflection(line: &RationalModalLine, omega: f64, incident: &[C64; 2]) -> [C64; 2] {
    let source = line.source();
    let n = source.coordinates();
    let s = C64::new(0.0, -omega);
    let z = line.load().source_model().eval(s);
    let lu = lu_complex(&z, 2).unwrap();
    let mut inverse = [C64::ZERO; 4];
    for col in 0..2 {
        let mut rhs = [C64::ZERO; 2];
        rhs[col] = C64::ONE;
        lu.solve(&mut rhs);
        for row in 0..2 {
            inverse[row * 2 + col] = rhs[row];
        }
    }
    let mut pencil: Vec<_> = source
        .stiffness()
        .iter()
        .zip(source.mass())
        .map(|(k, m)| C64::from_re(k - omega * omega * m))
        .collect();
    let mut rhs = vec![C64::ZERO; n];
    for i in 0..2 {
        pencil[i * n + i] = pencil[i * n + i] + s.scale(1.0 / source.inlet_reference_impedance());
        rhs[i] = incident[i].scale(2.0 / source.inlet_reference_impedance());
        for j in 0..2 {
            let at = (n - 2 + i) * n + n - 2 + j;
            pencil[at] = pencil[at] + s * inverse[i * 2 + j];
        }
    }
    lu_complex(&pencil, n).unwrap().solve(&mut rhs);
    [s * rhs[0] - incident[0], s * rhs[1] - incident[1]]
}

#[test]
fn mixed_rank_playback_converges_to_independent_continuous_mutual_response() {
    let source = model(1, false);
    let omega = 6000.0;
    let drive = [C64::new(1.0, 0.2), C64::new(-0.4, 0.3)];
    let gate = CancelGate::new_clock_free();
    let mut errors = Vec::new();
    for samples in [48, 96, 192] {
        let dt = TAU / (omega * samples as f64);
        let mut region = line(&source, dt);
        let expected = weak_reflection(&region, omega, &drive);
        let independent = region
            .prepare_analysis(
                RationalModalAnalysisBudget {
                    maximum_dimension: 64,
                    maximum_payload_bytes: 4 << 20,
                    maximum_work: 10_000_000,
                    reference_rate_per_s: omega,
                    relative_residual_tolerance: 1e-10,
                },
                &gate,
            )
            .unwrap()
            .factor(omega, &gate)
            .unwrap()
            .response(&drive, &gate)
            .unwrap();
        for i in 0..2 {
            assert!((independent.reflected_pressure_pa[i] - expected[i]).abs() < 1e-10);
        }
        let mut measured = [C64::ZERO; 2];
        for step in 0..24 * samples {
            let phase = omega * dt * (step as f64 + 0.5);
            let cos = det::cos(phase);
            let sin = det::sin(phase);
            let wave = [
                drive[0].re * cos + drive[0].im * sin,
                drive[1].re * cos + drive[1].im * sin,
            ];
            region.push(&wave, &gate).unwrap();
            if step >= 16 * samples {
                for i in 0..2 {
                    measured[i] = measured[i]
                        + C64::new(cos, sin)
                            .scale(2.0 * region.incoming()[i] / (8 * samples) as f64);
                }
            }
        }
        errors.push(
            (0..2)
                .map(|i| (measured[i] - expected[i]).abs())
                .fold(0.0, f64::max),
        );
    }
    assert!(
        errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 2e-4, "{errors:?}");
    eprintln!("48/96/192-samples-per-cycle reflection errors: {errors:?}");
}

#[test]
fn zero_and_full_rank_routes_keep_original_energy_and_mutual_memory() {
    let gate = CancelGate::new_clock_free();
    for rank in [0, 2] {
        let source = model(rank, true);
        let mut region = line(&source, 1.0 / 96_000.0);
        for i in 0..32 {
            let record = region
                .push(&[12.0, if i < 8 { -7.0 } else { 0.0 }], &gate)
                .unwrap();
            let scale = region.stored_duct_energy_j()
                + record.incident_energy_j
                + record.reflected_energy_j;
            close(record.balance_residual_j, 0.0, scale, 1e-8);
            for p in region.terminal_accounts() {
                close(p.duct_work_j, p.load_work_j, scale, 1e-8);
            }
        }
        assert_eq!(region.load().source_model(), &source);
        assert_eq!(region.load_state().len(), source.order() * 2);
    }
}

#[test]
fn mixed_rank_preserves_compatible_history_and_rejects_inconsistent_or_failed_updates() {
    let source = model(1, false);
    let base = prepared(&source, 1.0 / 96_000.0);
    // C*[0.32,0.52] has equal components, hence lies in D's positive range.
    // Nonzero original load history is compatible with a quiescent duct:
    // the algebraic range flow supplies the required opposing direct pressure.
    let state = [0.32e-8, 0.52e-8];
    let load = base
        .continuous_operator()
        .candidate_from_coordinates(&state, &[0.0; 2])
        .unwrap();
    let mut region = attach(load).unwrap();
    assert_eq!(region.load_state(), &state);
    let bad = base
        .continuous_operator()
        .candidate_from_coordinates(&[1e-8, 0.0], &[0.0; 2])
        .unwrap();
    assert!(attach(bad).is_err());
    let gate = CancelGate::new_clock_free();
    region.push(&[10.0, -2.0], &gate).unwrap();
    let accepted = snapshot(&region);
    assert!(region.push(&[1.0], &gate).is_err());
    assert!(region.push(&[f64::NAN, 0.0], &gate).is_err());
    assert_eq!(snapshot(&region), accepted);
    gate.request();
    assert!(region.push(&[0.0; 2], &gate).is_err());
    assert_eq!(snapshot(&region), accepted);
    let profile = RadiusProfile::exponential(RADIUS, 0.0, LENGTH).unwrap();
    let mut b = budget();
    b.exact.max_dense_bytes = 1;
    assert!(
        RationalModalLine::new(
            &profile,
            &gas(),
            2,
            ModalRationalLoad::new(prepared(&source, 1.0 / 96_000.0), RADIUS, gas()).unwrap(),
            b,
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
}
