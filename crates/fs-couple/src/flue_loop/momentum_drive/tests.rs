//! The momentum drive's law at hand-computed states, its refusals, the same-sample guard with the
//! term on, and the ledger. The card-row tests are in `tests/momentum_drive_card.rs`.

use super::*;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::flue_loop::forced_jet::COARSE_SPEC;
use crate::flue_loop::network::FlueNetworkSpec;
use crate::flue_loop::{FlueControl, RetardedVelocity};
use crate::wind_line::WindBoundary;
use fs_exec::CancelGate;
use fs_material::gas::{GasSpec, GasState};

const FLUE_W: f64 = 0.01;
const FLUE_H: f64 = 0.001;
const CUT_UP: f64 = 0.004;
const BORE_R: f64 = 0.006;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

// The voice of `mouth_radiation/tests.rs`: a 0.3 m pipe of radius 6 mm behind a 10 mm x 1 mm flue.
fn bare() -> FlueVoice {
    let g = gas();
    let card = fs_aeroac::jetcard::mint_tonal_interim_card().unwrap();
    let geometry = FlueGeometry {
        flue_width_m: FLUE_W,
        flue_height_m: FLUE_H,
        cut_up_m: CUT_UP,
        labium_offset_m: 0.0,
        jet_angle_rad: 0.0,
        bore_area_m2: core::f64::consts::PI * BORE_R * BORE_R,
        mouth_loss_rel: 0.0,
    };
    let source = FlueNetworkSpec {
        graph: TubeNetworkSpec {
            nodes: vec![
                NetworkNode::Inlet,
                NetworkNode::Termination { reflection: -0.9 },
            ],
            sections: vec![TubeSection {
                nodes: [0, 1],
                length_m: 0.3,
                radius_m: BORE_R,
                max_length_error_m: 0.005,
            }],
            sound_speed_m_s: g.sound_speed,
            max_wave_memory_bytes: 1_000_000,
        },
        gases: vec![g],
        viscothermal: vec![],
        minimum_jet_speed_m_s: 1.0,
        maximum_jet_mach: 0.3,
        max_jet_history_samples: 4096,
    };
    FlueVoice::from_network(
        Some(&card),
        90.0,
        1e-6,
        geometry,
        source,
        96_000,
        3,
        0,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

fn forced() -> FlueVoice {
    bare().with_forced_jet(COARSE_SPEC, 100.0).unwrap()
}

fn voice(kappa: Option<f64>) -> FlueVoice {
    let v = forced();
    match kappa {
        Some(k) => v.with_momentum_drive(k).unwrap(),
        None => v,
    }
}

fn coefficients(eta: Option<f64>) -> JetStepCoefficients {
    JetStepCoefficients {
        q_max: 2e-4,
        half_width_m: 0.0005,
        offset_m: 0.00025,
        band_alpha: 0.1,
        hp_alpha: 0.7,
        leak_rate: 10.0,
        eta_override: eta,
        at_rest: false,
    }
}

fn invalid(result: Result<FlueVoice, FlueError>) -> &'static str {
    match result {
        Err(FlueError::Invalid { what }) => what,
        Err(other) => panic!("expected an Invalid refusal, got {other:?}"),
        Ok(_) => panic!("expected a refusal"),
    }
}

#[test]
fn only_the_two_printed_weights_are_admitted() {
    assert!(is_printed_weight(1.0) && is_printed_weight(2.6));
    assert!(
        is_printed_weight("1.0".parse().unwrap()) && is_printed_weight("2.60".parse().unwrap())
    );
    for k in [
        0.0,
        -1.0,
        0.5,
        2.0,
        3.0,
        4.0,
        1.0 + f64::EPSILON,
        2.6 + 2.0 * f64::EPSILON,
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
    ] {
        assert!(!is_printed_weight(k), "{k} must not be a printed weight");
        assert_eq!(MomentumDrive::new(k, 1.2, 3e-4), Err(KAPPA_REFUSED), "{k}");
    }
    // The refusal names both printed values and where each is printed.
    assert!(KAPPA_REFUSED.contains("1 (Thwaites and Fletcher 1983 Eq. (1), p. 401)"));
    assert!(KAPPA_REFUSED.contains("2.6 (their Sec. II fit, p. 402)"));
    assert!(MomentumDrive::new(1.0, 0.0, 3e-4).is_err());
    assert!(MomentumDrive::new(1.0, 1.2, 0.0).is_err());
    assert!(MomentumDrive::new(1.0, f64::NAN, 3e-4).is_err());
    assert!(MomentumDrive::new(1.0, 1.2, f64::INFINITY).is_err());
}

#[test]
fn pressure_is_kappa_rho_v0_dq_over_sp_at_a_hand_computed_state() {
    // rho 1.2 kg/m^3, v0 10 m/s, S_p 3e-4 m^2, Q_j 5e-4 m^3/s, Q_rest 3e-4 m^3/s:
    // rho v0 = 12, times dQ = 2e-4 is 2.4e-3, over S_p is 8 Pa; kappa 2.6 gives 20.8 Pa.
    let one = MomentumDrive::new(1.0, 1.2, 3e-4).unwrap();
    assert!((one.pressure_per_flow(10.0) / 40_000.0 - 1.0).abs() < 1e-14);
    assert!((one.pressure_pa(10.0, 5e-4, 3e-4) / 8.0 - 1.0).abs() < 1e-14);
    // A planted wrong sign fails here: a flow below rest delivers a NEGATIVE pressure, above rest a positive one.
    assert!((one.pressure_pa(10.0, 1e-4, 3e-4) / -8.0 - 1.0).abs() < 1e-14);
    assert!(one.pressure_pa(10.0, 5e-4, 3e-4) > 0.0);
    assert_eq!(
        one.pressure_pa(10.0, 3e-4, 3e-4).to_bits(),
        0.0f64.to_bits(),
        "zero at rest, exactly"
    );
    let sec2 = MomentumDrive::new(2.6, 1.2, 3e-4).unwrap();
    assert!((sec2.pressure_pa(10.0, 5e-4, 3e-4) / 20.8 - 1.0).abs() < 1e-14);
    assert!(one.source().contains("Eq. (1), p. 401") && sec2.source().contains("Sec. II, p. 402"));
}

#[test]
fn the_voice_adds_exactly_the_term_above_and_below_rest_and_nothing_at_rest() {
    let g = gas();
    let bore = core::f64::consts::PI * BORE_R * BORE_R;
    let c = coefficients(None);
    let v0 = c.q_max / (FLUE_H * FLUE_W);
    let rest = 0.5 * c.q_max * (1.0 + (-c.offset_m / c.half_width_m).tanh());
    assert!((rest_cut_flow_m3_s(c.q_max, c.offset_m, c.half_width_m) - rest).abs() <= 1e-18);
    let off = voice(None);
    for kappa in [1.0, 2.6] {
        let on = voice(Some(kappa));
        for eta in [0.0004, -0.0004] {
            let c = coefficients(Some(eta));
            let a = off.jet_trial(0.0, 0.0, c).unwrap();
            let b = on.jet_trial(0.0, 0.0, c).unwrap();
            assert_eq!(
                a.flow_m3_s.to_bits(),
                b.flow_m3_s.to_bits(),
                "the cut itself is untouched"
            );
            let expected = kappa * g.density * v0 * (a.flow_m3_s - rest) / bore;
            let added = b.source_pa - a.source_pa;
            assert!(expected.abs() > 1e-3);
            assert!(
                (added / expected - 1.0).abs() < 1e-12,
                "kappa {kappa} eta {eta}: {added} vs {expected}"
            );
            assert_eq!(
                added > 0.0,
                a.flow_m3_s > rest,
                "the sign follows the flow's deviation from rest"
            );
        }
        // At rest (labium displacement zero) the cut flow IS the quiescent flow: no source, bit for bit.
        let c = coefficients(Some(0.0));
        let a = off.jet_trial(0.0, 0.0, c).unwrap();
        let b = on.jet_trial(0.0, 0.0, c).unwrap();
        assert_eq!(a.flow_m3_s.to_bits(), rest.to_bits());
        assert_eq!(a.source_pa.to_bits(), b.source_pa.to_bits());
    }
}

#[test]
fn without_the_row_the_source_is_the_original_expression_bit_for_bit() {
    let v = voice(None);
    assert!(v.momentum_drive().is_none());
    let c = coefficients(None);
    assert_eq!(v.momentum_source_slope(c).to_bits(), 0.0f64.to_bits());
    for flow in [0.0, 1e-5, 5.3e-5, 1.99e-4, -3e-5] {
        let original = v.source_coeff * (flow - v.flow_prev) / v.dt;
        assert_eq!(
            v.jet_source_pa(flow, c).unwrap().to_bits(),
            original.to_bits(),
            "flow {flow}"
        );
    }
}

#[test]
fn every_other_weight_and_every_misuse_is_refused() {
    for k in [
        0.0,
        0.5,
        2.0,
        3.0,
        4.0,
        -1.0,
        1.0 + f64::EPSILON,
        f64::NAN,
        f64::INFINITY,
    ] {
        let what = invalid(forced().with_momentum_drive(k));
        assert_eq!(what, KAPPA_REFUSED, "{k}");
    }
    // It acts on the forced jet's cut flow: no forced jet, no momentum drive.
    assert!(invalid(bare().with_momentum_drive(1.0)).contains("after the forced jet"));
    // Once.
    assert!(invalid(voice(Some(1.0)).with_momentum_drive(1.0)).contains("one cold attachment"));
    // Cold only.
    let mut warm = forced();
    warm.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
    warm.step_block(&mut [0.0; 64]).unwrap();
    assert!(invalid(warm.with_momentum_drive(1.0)).contains("one cold attachment"));
    // The provenance names the weight, where it is printed, and the caveats.
    let on = voice(Some(2.6));
    let provenance = &on.island().provenance;
    assert!(provenance.contains("momentum drive: kappa 2.6"));
    assert!(provenance.contains("Sec. II, p. 402"));
    assert!(provenance.contains("laminar Bickley kernel"));
}

// Smallest weight on the retarded velocity at which `solve_jet_mouth` refuses by the feedback bound.
fn first_refused_weight(mut v: FlueVoice) -> f64 {
    let bulk = (2.0 * 100.0 / v.rho).sqrt();
    assert!(v.forced.as_mut().unwrap().prepare(bulk).unwrap().is_some());
    assert!(
        v.forced.as_ref().unwrap().current().1 > 0.0,
        "Kutta g_0 must feed back positively"
    );
    let boundary = WindBoundary {
        direct_reflection: 0.25,
        history_pressure_pa: 0.0,
    };
    let refused = |w: f64| -> bool {
        let delayed = RetardedVelocity {
            past_m_s: 0.0,
            current_weight: w,
        };
        match v.solve_jet_mouth(boundary, delayed, 0.0, coefficients(None), None, None) {
            Ok(_) => false,
            Err(FlueError::Invalid { what })
                if what.contains("current jet feedback is underresolved") =>
            {
                true
            }
            Err(other) => panic!("weight {w}: unexpected refusal {other:?}"),
        }
    };
    let mut low = 1e-3;
    assert!(!refused(low), "a tiny current weight must be admitted");
    let mut high = low;
    while !refused(high) {
        high *= 2.0;
        assert!(high < 1e30, "the bound never engaged");
        if !refused(high) {
            low = high;
        }
    }
    for _ in 0..60 {
        let mid = 0.5 * (low + high);
        if refused(mid) {
            high = mid;
        } else {
            low = mid;
        }
    }
    high
}

#[test]
fn the_same_sample_guard_sees_the_term_and_is_never_loosened() {
    let off = voice(None);
    let sc_over_dt = off.source_coeff / off.dt;
    let w_off = first_refused_weight(voice(None));
    for kappa in [1.0, 2.6] {
        let on = voice(Some(kappa));
        let bulk = coefficients(None).q_max / (FLUE_H * FLUE_W);
        let m = on.momentum_drive().unwrap().pressure_per_flow(bulk);
        assert!(m > 0.0);
        let w_on = first_refused_weight(voice(Some(kappa)));
        assert!(
            w_on < w_off,
            "kappa {kappa}: the term must tighten the bound ({w_on} vs {w_off})"
        );
        // feedback is linear in the weight and in (source_coeff/dt + m): the thresholds scale inversely.
        let predicted = (sc_over_dt + m) / sc_over_dt;
        assert!(
            ((w_off / w_on) / predicted - 1.0).abs() < 1e-9,
            "kappa {kappa}: {} vs {predicted}",
            w_off / w_on
        );
        println!(
            "kappa {kappa}: refused from weight {w_on:.6e} (without the term {w_off:.6e}); slope ratio {predicted:.6}"
        );
    }
}

#[test]
fn with_the_term_on_the_ledger_closes_and_the_sound_changes() {
    let run = |kappa: Option<f64>| {
        let mut v = voice(kappa);
        v.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        let mut out = Vec::new();
        for _ in 0..6 {
            let mut block = [0.0; 4096];
            v.step_block(&mut block).unwrap();
            out.extend_from_slice(&block);
        }
        for d in v.diagnostics() {
            let scale = d.source_work_j.abs()
                + d.pipe_work_j.abs()
                + d.stored_delta_j.abs()
                + d.radiation_loss_j;
            assert!(
                d.ledger_defect_j.abs() <= 1e-9 * scale.max(f64::MIN_POSITIVE),
                "kappa {kappa:?}: ledger defect {} against scale {scale}",
                d.ledger_defect_j
            );
        }
        assert!(out.iter().all(|x| x.is_finite()));
        out
    };
    let base = run(None);
    let a = run(Some(1.0));
    let b = run(Some(2.6));
    let differs = |x: &[f64], y: &[f64]| x.iter().zip(y).any(|(p, q)| p.to_bits() != q.to_bits());
    assert!(
        differs(&base, &a) && differs(&a, &b),
        "the term must change the loop"
    );
    // Replay is exact.
    let again = run(Some(1.0));
    assert!(
        a.iter()
            .zip(&again)
            .all(|(p, q)| p.to_bits() == q.to_bits())
    );
}

#[test]
fn g1_profile_source_uses_the_original_station_for_volume_and_momentum() {
    let on = forced().with_bickley_momentum_drive().unwrap();
    let off = forced().with_jet_cut_width_station().unwrap();
    assert!(on.jet_cut_width_station());
    assert!(on.momentum_drive().is_none());
    assert!(on.bickley_momentum_drive().is_some());
    assert!(on.island().provenance.contains("momentum drive: bickley-profile"));
    let b = on.cut_half_width_m();
    let u = 20.0;
    for eta in [-b, 0.0, 0.25 * b, b] {
        let c = JetStepCoefficients {
            q_max: 2.0 * b * FLUE_W * u,
            half_width_m: b,
            ..coefficients(Some(eta))
        };
        let baseline = off.jet_trial(0.0, 0.0, c).unwrap();
        let profile = on.jet_trial(0.0, 0.0, c).unwrap();
        assert_eq!(profile.flow_m3_s, baseline.flow_m3_s);
        let s = ((eta-c.offset_m)/b).tanh();
        let rest = (-c.offset_m/b).tanh();
        let integral = (s-rest) - (s.powi(3)-rest.powi(3))/3.0;
        let pressure = on.rho*u*u*b*FLUE_W/on.geometry.bore_area_m2*integral;
        if eta == 0.0 {
            assert_eq!(profile.source_pa.to_bits(), baseline.source_pa.to_bits());
        } else {
            assert!((profile.source_pa-baseline.source_pa-pressure).abs() <= 1e-11*pressure.abs());
        }
    }
    assert!(invalid(bare().with_bickley_momentum_drive()).contains("after the forced jet"));
    assert!(invalid(voice(Some(1.0)).with_bickley_momentum_drive()).contains("excludes"));
    assert!(invalid(forced().with_bickley_momentum_drive().unwrap()
        .with_momentum_drive(1.0)).contains("one cold attachment"));
}

#[test]
fn g0_profile_same_sample_guard_uses_the_exact_derivative_upper_bound() {
    let off = forced().with_jet_cut_width_station().unwrap();
    let source = off.source_coeff/off.dt;
    let on = forced().with_bickley_momentum_drive().unwrap();
    let c = coefficients(None);
    let u = c.q_max/(2.0*c.half_width_m*FLUE_W);
    let expected = on.rho*u/on.geometry.bore_area_m2;
    assert_eq!(on.momentum_source_slope(c), expected);
    let w_off = first_refused_weight(off);
    let w_on = first_refused_weight(on);
    assert!(w_on < w_off);
    assert!(((w_off/w_on)/((source+expected)/source)-1.0).abs() < 1e-9);
}

#[test]
fn g0_g5_profile_source_advances_the_actual_loop_with_closed_work_and_replay() {
    let run = |profile: bool| {
        let mut voice = if profile {
            forced().with_bickley_momentum_drive().unwrap()
        } else {
            forced().with_jet_cut_width_station().unwrap()
        };
        voice.apply(FlueControl::SetBlowingPressure(100.0)).unwrap();
        let mut samples = Vec::new();
        for _ in 0..2 {
            let mut block = [0.0; 4096];
            voice.step_block(&mut block).unwrap();
            samples.extend_from_slice(&block);
        }
        for d in voice.diagnostics() {
            let scale = d.source_work_j.abs() + d.pipe_work_j.abs()
                + d.stored_delta_j.abs() + d.radiation_loss_j;
            assert!(d.ledger_defect_j.abs() <= 1e-9*scale.max(f64::MIN_POSITIVE));
        }
        assert!(samples.iter().all(|v| v.is_finite()));
        samples
    };
    let off = run(false);
    let on = run(true);
    assert!(off.iter().zip(&on).any(|(a,b)| a.to_bits() != b.to_bits()));
    assert_eq!(on, run(true));
}
