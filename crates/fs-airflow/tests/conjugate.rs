//! Conjugate coupling battery: the air path against closed form, the exact
//! segment-refinement falsifier, the coupled fixed point against its
//! manufactured solution, the flux-balance audit's sensitivity, and the
//! cancellation drill.
//!
//! Every physical check here has an oracle that is independent of the code it
//! checks:
//!
//! * the single-segment march is compared against `T_w − (T_w − T_in)e^(−NTU)`
//!   evaluated through `f64::exp`, not through the crate's `det::expm1` path;
//! * segment refinement is an EXACT invariance with no tolerance to tune, and
//!   the rejected arithmetic-mean model is built here and shown to fail it, so
//!   the invariance test is demonstrably sensitive rather than merely
//!   satisfied;
//! * the coupled fixed point over a lumped solid has a closed-form solution
//!   for BOTH the solid temperature and the outlet air temperature, derived
//!   below from a global energy statement the driver never evaluates.

use fs_airflow::AirflowError;
use fs_airflow::conjugate::{
    AirPath, AirSegment, ConjugateConfig, Relaxation, SegmentRefinementTransfer, SolidRegionState,
    seam_effort_dimensions, solve_conjugate, solve_conjugate_from,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ladder::{LadderRegistry, Transfer};
use fs_qty::Temperature;

const INLET_K: f64 = 300.0;
const MASS_FLOW_KG_S: f64 = 0.004;
const CP_J_KG_K: f64 = 1_005.0;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    with_gated_cx(&gate, f)
}

fn with_gated_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 0x0000_C0FF_EE00_0001,
                kernel_id: 71,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn capacity_rate() -> f64 {
    MASS_FLOW_KG_S * CP_J_KG_K
}

fn assert_invalid_field(error: AirflowError, expected: &'static str) {
    match error {
        AirflowError::InvalidConjugateInput { field, .. } => {
            assert_eq!(field, expected);
        }
        other => panic!("expected InvalidConjugateInput({expected}), got {other:?}"),
    }
}

fn assert_non_finite_stage(error: AirflowError, expected: &'static str) {
    match error {
        AirflowError::NonFiniteCoupling { stage, .. } => {
            assert_eq!(stage, expected);
        }
        other => panic!("expected NonFiniteCoupling({expected}), got {other:?}"),
    }
}

/// `N` equal segments spanning the same total area and coefficient.
fn uniform_path(segments: usize, total_area_m2: f64, htc: f64) -> AirPath {
    let per = total_area_m2 / segments as f64;
    let declared = (0..segments)
        .map(|index| AirSegment::new(&format!("channel-{index}"), per, htc).expect("valid segment"))
        .collect();
    AirPath::new(INLET_K, MASS_FLOW_KG_S, CP_J_KG_K, declared).expect("valid path")
}

// ---------------------------------------------------------------------------
// Air path against the closed-form heated channel
// ---------------------------------------------------------------------------

#[test]
fn single_segment_march_matches_the_closed_form_heated_channel() {
    // Constant wall temperature: the exact 1-D solution of
    // ṁ c_p dT/dx = h P (T_w − T) is T(x) = T_w − (T_w − T_in) e^(−h A(x)/ṁc_p).
    let wall = 350.0;
    let area = 0.02;
    let htc = 60.0;
    let path = uniform_path(1, area, htc);
    let march = path.march(&[wall]).expect("march");

    let ntu = htc * area / capacity_rate();
    let expected_outlet = wall - (wall - INLET_K) * (-ntu).exp();
    let expected_heat = capacity_rate() * (expected_outlet - INLET_K);

    assert!(
        (march.outlet_temperature_k - expected_outlet).abs() < 1.0e-12,
        "outlet {} vs closed form {expected_outlet}",
        march.outlet_temperature_k
    );
    assert!(
        (march.total_heat_rate_w - expected_heat).abs() < 1.0e-9,
        "heat {} vs closed form {expected_heat}",
        march.total_heat_rate_w
    );
}

#[test]
fn reference_temperature_reproduces_the_segment_heat_rate_exactly() {
    // T_ref,eff is DEFINED so that h A (T_w − T_ref) ≡ Q. If that identity
    // fails, the Robin row handed to the solid describes a different heat
    // rate than the air actually absorbed and the two sides are solving
    // different problems.
    let wall = 372.5;
    for &htc in &[1.0, 25.0, 300.0, 5_000.0] {
        for &area in &[1.0e-4, 0.01, 0.5] {
            let path = uniform_path(1, area, htc);
            let march = path.march(&[wall]).expect("march");
            let segment = &march.segments[0];
            let robin_rate = htc * area * (wall - segment.reference_temperature_k);
            let relative =
                (robin_rate - segment.heat_rate_w).abs() / segment.heat_rate_w.abs().max(1.0e-30);
            assert!(
                relative < 1.0e-12,
                "h={htc} A={area}: Robin rate {robin_rate} vs marched {}",
                segment.heat_rate_w
            );
        }
    }
}

#[test]
fn air_stays_between_the_inlet_and_the_wall_at_every_ntu() {
    let wall = 400.0;
    // NTU up to 25 covers the deeply saturated regime where an arithmetic
    // mean reference leaves the physical interval entirely.
    for &htc in &[0.1, 1.0, 10.0, 100.0, 1_000.0, 5_000.0] {
        let path = uniform_path(1, 0.05, htc);
        let march = path.march(&[wall]).expect("march");
        let segment = &march.segments[0];
        assert!(
            (INLET_K..=wall).contains(&segment.outlet_temperature_k),
            "outlet {} left [{INLET_K}, {wall}] at NTU {}",
            segment.outlet_temperature_k,
            segment.ntu
        );
        assert!(
            (INLET_K..=wall).contains(&segment.reference_temperature_k),
            "reference {} left [{INLET_K}, {wall}] at NTU {}",
            segment.reference_temperature_k,
            segment.ntu
        );
    }
}

#[test]
fn the_arithmetic_mean_reference_leaves_the_physical_interval() {
    // The model this crate deliberately does NOT use: T_ref = (T_in + T_out)/2
    // with a linear T_out = T_in + Q/(ṁc_p) and Q = h A (T_w − T_ref). Solving
    // that pair gives T_ref = T_in + (T_w − T_in)·NTU/(2 + NTU) ... which is
    // fine; the failure is in the OUTLET it implies:
    // T_out = T_in + (T_w − T_in)·2NTU/(2 + NTU), which exceeds T_w for
    // NTU > 2. Past that point the air leaves hotter than the wall heating it.
    let wall = 400.0;
    let ntu = 4.0;
    let arithmetic_outlet = INLET_K + (wall - INLET_K) * 2.0 * ntu / (2.0 + ntu);
    assert!(
        arithmetic_outlet > wall,
        "the arithmetic-mean model is supposed to be non-physical here: {arithmetic_outlet}"
    );

    // The exponential law at the same NTU stays inside the interval.
    let htc = 1_000.0;
    let area = ntu * capacity_rate() / htc;
    let path = uniform_path(1, area, htc);
    let march = path.march(&[wall]).expect("march");
    assert!((march.segments[0].ntu - ntu).abs() < 1.0e-12);
    assert!(march.outlet_temperature_k < wall);
}

// ---------------------------------------------------------------------------
// The exact falsifier: segment-refinement invariance
// ---------------------------------------------------------------------------

#[test]
fn segment_refinement_leaves_a_uniform_wall_outlet_invariant() {
    // e^(−NTU) = (e^(−NTU/N))^N, so splitting one uniform-wall channel into N
    // equal segments must give the SAME outlet for every N. Exact identity,
    // no tolerance to tune beyond floating-point accumulation.
    let wall = 355.0;
    let total_area = 0.03;
    let htc = 80.0;

    let reference = uniform_path(1, total_area, htc)
        .march(&[wall])
        .expect("march")
        .outlet_temperature_k;

    for &n in &[2_usize, 3, 4, 8, 16, 64, 256] {
        let path = uniform_path(n, total_area, htc);
        let march = path.march(&vec![wall; n]).expect("march");
        assert!(
            (march.outlet_temperature_k - reference).abs() < 1.0e-10,
            "N={n} outlet {} drifted from the N=1 outlet {reference}",
            march.outlet_temperature_k
        );
        // The heat rate is the same statement in energy form.
        let expected_heat = capacity_rate() * (reference - INLET_K);
        assert!(
            (march.total_heat_rate_w - expected_heat).abs() < 1.0e-8,
            "N={n} heat {} vs {expected_heat}",
            march.total_heat_rate_w
        );
    }
}

#[test]
fn an_arithmetic_mean_model_fails_segment_refinement_with_clean_n_dependence() {
    // This is the sensitivity proof for the test above. If the invariance
    // check passed for any reasonable model it would falsify nothing. Build
    // the rejected model here and show refinement moves its answer.
    //
    // Arithmetic-mean march, per segment:
    //   Q     = h A (T_w − (T_in + T_out)/2)
    //   T_out = T_in + Q/(ṁ c_p)
    //   ⇒ T_out = T_in + (T_w − T_in) · 2·NTU/(2 + NTU)
    fn arithmetic_outlet(segments: usize, total_area: f64, htc: f64, wall: f64) -> f64 {
        let per_ntu = htc * (total_area / segments as f64) / capacity_rate();
        let gain = 2.0 * per_ntu / (2.0 + per_ntu);
        let mut t = INLET_K;
        for _ in 0..segments {
            t += (wall - t) * gain;
        }
        t
    }

    let wall = 355.0;
    let htc = 80.0;
    // Sweep NTU by stretching the channel. The arithmetic mean's error is
    // second order in the per-segment NTU, so the defect is small in a mild
    // channel and severe in a saturated one. Both facts matter: the first is
    // why the invariance check needs a tolerance far below the defect rather
    // than a loose one, the second is why the model is unusable in enclosures.
    let mut defects = Vec::new();
    for &total_area in &[0.03_f64, 0.10, 0.30, 1.00] {
        let one = arithmetic_outlet(1, total_area, htc, wall);
        let many = arithmetic_outlet(4_096, total_area, htc, wall);
        let ntu = htc * total_area / capacity_rate();
        let defect = (many - one).abs();

        // The exponential law is the N → ∞ limit at every N, so the refined
        // arithmetic march must approach it — which is what makes the coarse
        // arithmetic answer wrong rather than merely a different convention.
        let exact = uniform_path(1, total_area, htc)
            .march(&[wall])
            .expect("march")
            .outlet_temperature_k;
        assert!(
            (many - exact).abs() < 0.05,
            "NTU {ntu}: refined arithmetic mean {many} should approach {exact}"
        );

        // The defect must exceed the invariance test's 1e-10 K tolerance by
        // orders of magnitude, or that test would not be falsifying anything.
        assert!(
            defect > 1.0e-3,
            "NTU {ntu}: rejected model must be visibly N-dependent, saw {defect} K"
        );
        defects.push((ntu, defect));
    }

    // The defect grows with NTU — the failure mode is not a fixed offset.
    for window in defects.windows(2) {
        assert!(
            window[1].1 > window[0].1,
            "defect did not grow with NTU: {:?} then {:?}",
            window[0],
            window[1]
        );
    }
    let (worst_ntu, worst) = *defects.last().expect("nonempty");
    assert!(
        worst > 10.0,
        "at NTU {worst_ntu} the rejected model should be catastrophically \
         N-dependent, saw {worst} K"
    );
}

#[test]
fn downstream_segments_see_hotter_air_than_upstream_ones() {
    // The property that did not exist before this module: air carries state.
    let path = uniform_path(4, 0.04, 90.0);
    let march = path.march(&[360.0; 4]).expect("march");
    for window in march.segments.windows(2) {
        assert!(
            window[1].inlet_temperature_k > window[0].inlet_temperature_k,
            "air did not heat between segments: {:?}",
            march.segments
        );
        // Each segment sees a smaller driving temperature difference, so it
        // absorbs strictly less heat than the one before it.
        assert!(
            window[1].heat_rate_w < window[0].heat_rate_w,
            "downstream segment absorbed at least as much heat as upstream"
        );
    }
}

// ---------------------------------------------------------------------------
// The coupled fixed point against its manufactured solution
// ---------------------------------------------------------------------------

/// A lumped solid at one uniform temperature absorbing `power_w` and losing it
/// only through the declared segments:
/// `P = Σ_j h_j A_j (T_s − T_ref,j)`.
///
/// Its coupled fixed point is known in closed form. Every watt injected must
/// leave through the air, so `T_out = T_in + P/(ṁ c_p)`; and with one uniform
/// wall temperature the whole path is a single exponential of the summed NTU,
/// giving `T_s = T_in + P / (ṁ c_p (1 − e^(−ΣNTU)))`. Neither expression is
/// evaluated anywhere in the driver.
fn lumped_solid(
    regions: &[(String, f64, f64)],
    power_w: f64,
    references: &[f64],
) -> Vec<SolidRegionState> {
    let conductance: f64 = regions.iter().map(|(_, area, htc)| area * htc).sum();
    let weighted: f64 = regions
        .iter()
        .zip(references)
        .map(|((_, area, htc), t_ref)| area * htc * t_ref)
        .sum();
    let solid_temperature = (power_w + weighted) / conductance;
    regions
        .iter()
        .zip(references)
        .map(|((region, area, htc), t_ref)| SolidRegionState {
            region: region.clone(),
            area_m2: *area,
            mean_wall_temperature_k: solid_temperature,
            heat_rate_w: area * htc * (solid_temperature - t_ref),
            // Claim the reference actually used, so every lumped test also
            // exercises the reference-wiring cross-check live.
            mean_reference_temperature_k: Some(*t_ref),
        })
        .collect()
}

fn lumped_fixture(segments: usize) -> (AirPath, Vec<(String, f64, f64)>) {
    let total_area = 0.03;
    let htc = 85.0;
    let path = uniform_path(segments, total_area, htc);
    let regions = path
        .segments()
        .iter()
        .map(|s| (s.region().to_string(), s.area_m2(), s.htc_w_per_m2_k()))
        .collect();
    (path, regions)
}

#[test]
fn the_coupled_fixed_point_reproduces_its_closed_form_solution() {
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };

    let solution = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("the lumped conjugate exchange converges")
    });

    // Oracle 1 — global energy. Every injected watt leaves through the air.
    let expected_outlet = INLET_K + power_w / capacity_rate();
    assert!(
        (solution.march.outlet_temperature_k - expected_outlet).abs() < 1.0e-9,
        "outlet {} vs energy-conservation oracle {expected_outlet}",
        solution.march.outlet_temperature_k
    );

    // Oracle 2 — the closed-form solid temperature.
    let total_ntu: f64 = path
        .segments()
        .iter()
        .map(|s| s.ntu(path.capacity_rate_w_per_k()))
        .sum();
    let expected_solid = INLET_K + power_w / (capacity_rate() * (1.0 - (-total_ntu).exp()));
    let solved_solid = solution.solid[0].mean_wall_temperature_k;
    assert!(
        (solved_solid - expected_solid).abs() < 1.0e-9,
        "solid {solved_solid} vs closed form {expected_solid}"
    );

    // The air actually heated: the exchange is conjugate, not a declared
    // ambient.
    assert!(solution.march.outlet_temperature_k > INLET_K + 1.0);
    assert!(
        solution.reference_temperatures_k[3] > solution.reference_temperatures_k[0],
        "downstream reference must exceed upstream: {:?}",
        solution.reference_temperatures_k
    );
}

#[test]
fn the_history_is_a_complete_monotone_replay_record() {
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let solution = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("converges")
    });

    assert_eq!(solution.history.len(), solution.iterations);
    for (index, record) in solution.history.iter().enumerate() {
        assert_eq!(record.iteration, index);
        assert_eq!(record.reference_temperatures_k.len(), path.segments().len());
    }
    // The first iteration starts from the inlet everywhere.
    assert!(
        solution.history[0]
            .reference_temperatures_k
            .iter()
            .all(|t| *t == INLET_K)
    );
    // Contraction: this seam's composite gain is (1 − ε/NTU) times the solid's
    // gain, both strictly below one, so plain staggering converges without
    // relaxation and the residual falls monotonically. That is a property of
    // the seam, not a tuned schedule.
    for window in solution.history.windows(2) {
        assert!(
            window[1].max_reference_change_k < window[0].max_reference_change_k,
            "residual did not contract: {} then {}",
            window[0].max_reference_change_k,
            window[1].max_reference_change_k
        );
    }
    assert!(
        solution
            .history
            .last()
            .expect("nonempty")
            .max_reference_change_k
            <= 1.0e-12
    );
}

#[test]
fn aitken_relaxation_reaches_the_same_fixed_point() {
    // fs-couple's relaxer is scalar; this checks it is wired correctly and
    // lands on the same answer. It is NOT a claim that relaxation rescues a
    // divergent case — this seam is a contraction, so there is no divergent
    // case here to rescue.
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let plain = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let aitken = ConjugateConfig {
        relaxation: Relaxation::Aitken {
            omega_init: 1.0,
            omega_max: 2.0,
        },
        ..plain
    };

    let run = |config: ConjugateConfig| {
        with_cx(|cx| {
            solve_conjugate(cx, &path, &config, |_, references| {
                Ok(lumped_solid(&regions, power_w, references))
            })
            .expect("converges")
        })
    };
    let a = run(plain);
    let b = run(aitken);
    assert!(
        (a.march.outlet_temperature_k - b.march.outlet_temperature_k).abs() < 1.0e-9,
        "relaxation changed the fixed point: {} vs {}",
        a.march.outlet_temperature_k,
        b.march.outlet_temperature_k
    );
}

#[test]
fn exhausting_the_iteration_budget_refuses_with_a_diagnosis() {
    // No silent max-iter answer: the temperatures reached are not returned.
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-14,
        max_iterations: 2,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect_err("two iterations cannot reach 1e-14 K")
    });
    match error {
        AirflowError::ConjugateNotConverged {
            iterations,
            max_change_bits,
            tolerance_bits,
        } => {
            assert_eq!(iterations, 2);
            assert!(f64::from_bits(max_change_bits) > f64::from_bits(tolerance_bits));
        }
        other => panic!("expected a typed non-convergence refusal, got {other:?}"),
    }
}

#[test]
fn cancellation_stops_at_an_outer_iteration_boundary() {
    // G4 drill. The gate is tripped from inside the solid callback during the
    // first exchange; the driver must refuse at the NEXT checkpoint rather
    // than finishing the loop or returning a partial answer.
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let gate = CancelGate::new();
    let mut calls = 0_usize;
    let error = with_gated_cx(&gate, |cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            calls += 1;
            gate.request();
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect_err("a cancelled exchange must refuse")
    });
    let AirflowError::Cancelled {
        iteration,
        references_k,
    } = &error
    else {
        panic!("expected a cancellation refusal, got {error:?}")
    };
    assert_eq!(*iteration, 1);
    assert_eq!(calls, 1, "the solid side ran again after cancellation");
    // The refusal must carry enough state to resume from, not just an index.
    assert_eq!(references_k.len(), 4);
    assert!(references_k.iter().all(|t| t.is_finite()));
}

// ---------------------------------------------------------------------------
// The flux-balance audit and its sensitivity
// ---------------------------------------------------------------------------

#[test]
fn the_per_region_balance_closes_at_the_fixed_point() {
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let solution = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("converges")
    });

    assert_eq!(solution.balance.regions.len(), 4);
    assert!(
        solution.balance.max_region_imbalance_w < 1.0e-8,
        "per-region imbalance {} W",
        solution.balance.max_region_imbalance_w
    );
    assert!(solution.balance.interface_imbalance_w.abs() < 1.0e-8);
    // The solid delivered exactly the injected power.
    assert!((solution.balance.solid_total_w - power_w).abs() < 1.0e-8);
    assert!(!solution.worst_recorded_imbalance_w.is_nan());
}

#[test]
fn the_per_region_balance_is_sensitive_to_a_wiring_fault() {
    // Sensitivity proof for the audit. A solid side that reports one region's
    // heat against HALF its true area — the "dropped face" fault the audit
    // exists to catch — must show up as a per-region imbalance. Without this,
    // "the balance closed" would be an untested claim.
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            let mut states = lumped_solid(&regions, power_w, references);
            // Region 2 loses half of its faces on the way out of the solver.
            states[2].heat_rate_w *= 0.5;
            Ok(states)
        })
        .expect_err("a dropped-face fault must REFUSE, not return an answer")
    });

    // The critical property: the fault is INVISIBLE to the temperature
    // criterion. The reference temperatures still reach their fixed point, so
    // a driver gated only on kelvin publishes this run as converged. A small
    // residual in K cannot bound a heat rate in W without a conductance, so
    // the watt-denominated gate is a separate obligation, not a restatement.
    match error {
        AirflowError::ConjugateBalanceUnclosed {
            max_region_imbalance_bits,
            tolerance_bits,
            ..
        } => {
            let imbalance = f64::from_bits(max_region_imbalance_bits);
            assert!(
                imbalance > 1.0,
                "the dropped-face imbalance should be large, saw {imbalance} W"
            );
            assert!(imbalance > f64::from_bits(tolerance_bits));
        }
        other => panic!("expected a balance refusal, got {other:?}"),
    }
}

#[test]
fn the_watt_gate_scales_with_the_interface_not_an_absolute_constant() {
    // The SAME dropped-face fault as above, at microwatt scale. Under an
    // absolute-only 1e-6 W gate this mis-wired run RETURNED Ok: the fault's
    // whole signal sits below a constant chosen for watt-scale runs, and a
    // mW-vs-W unit mistake shrinks signal and gate sensitivity together.
    // Keying the gate to the interface's own heat-rate magnitude makes the
    // O(scale) fault refuse at every wattage.
    let power_w = 2.0e-6;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            let mut states = lumped_solid(&regions, power_w, references);
            states[2].heat_rate_w *= 0.5;
            Ok(states)
        })
        .expect_err("the microwatt dropped-face fault must refuse too")
    });
    match error {
        AirflowError::ConjugateBalanceUnclosed {
            max_region_imbalance_bits,
            tolerance_bits,
            ..
        } => {
            let imbalance = f64::from_bits(max_region_imbalance_bits);
            assert!(
                imbalance < 1.0e-6,
                "the whole point: this signal hides under a watt-scale absolute \
                 gate, saw {imbalance} W"
            );
            assert!(imbalance > 0.0, "a zero imbalance proves nothing");
            assert!(imbalance > f64::from_bits(tolerance_bits));
        }
        other => panic!("expected a balance refusal, got {other:?}"),
    }
}

#[test]
fn the_relative_term_admits_legitimate_floating_point_residue_at_scale() {
    // Force admission to come from the relative term ALONE: an absolute floor
    // of 1e-300 W admits no real residue by itself, while a legitimate 45 W
    // fixed point closes its interface only to floating-point residue — above
    // any such floor, far below one part per million of scale. If the
    // relative term were dead, this run would refuse.
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-300,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let solution = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("the relative term admits scale-proportional residue")
    });
    assert!(
        solution.balance.max_region_imbalance_w > 1.0e-300,
        "residue {} W should exceed the absolute floor, or this test proves \
         nothing about the relative term",
        solution.balance.max_region_imbalance_w
    );
}

#[test]
fn inert_and_oversized_relaxation_factors_refuse_at_admission() {
    // ω = 0 never moves the reference: admitting it burns max_iterations of
    // solid solves and then blames convergence for a configuration fault.
    let (path, regions) = lumped_fixture(2);
    for relaxation in [
        Relaxation::Fixed { omega: 0.0 },
        Relaxation::Fixed { omega: -0.5 },
        Relaxation::Fixed { omega: 2.5 },
        Relaxation::Aitken {
            omega_init: 0.0,
            omega_max: 2.0,
        },
        Relaxation::Aitken {
            omega_init: 3.0,
            omega_max: 2.0,
        },
    ] {
        let config = ConjugateConfig {
            relaxation,
            ..ConjugateConfig::default()
        };
        let error = with_cx(|cx| {
            solve_conjugate(cx, &path, &config, |_, references| {
                Ok(lumped_solid(&regions, 10.0, references))
            })
            .expect_err("a stalling relaxation must refuse at admission")
        });
        assert!(
            matches!(
                error,
                AirflowError::InvalidConjugateInput {
                    field: "conjugate relaxation omega",
                    ..
                }
            ),
            "{relaxation:?} must refuse as a config fault, got {error:?}"
        );
    }
}

#[test]
fn a_relaxation_overflow_is_attributed_to_the_relaxation_not_the_solid() {
    // An admissible-but-reckless Aitken magnitude cap combined with a far-off
    // resumed reference overflows the relaxed vector. The refusal must name
    // the relaxation stage: the solid answered honestly, and blaming it sends
    // the caller debugging the wrong component.
    let (path, regions) = lumped_fixture(2);
    let config = ConjugateConfig {
        relaxation: Relaxation::Aitken {
            omega_init: 1.0e300,
            omega_max: 1.0e300,
        },
        ..ConjugateConfig::default()
    };
    let initial = vec![1.0e30; 2];
    let error = with_cx(|cx| {
        solve_conjugate_from(cx, &path, &config, &initial, |_, references| {
            Ok(lumped_solid(&regions, 10.0, references))
        })
        .expect_err("an overflowed relaxation must refuse")
    });
    assert!(
        matches!(
            error,
            AirflowError::NonFiniteCoupling {
                stage: "relaxed reference temperature",
                ..
            }
        ),
        "got {error:?}"
    );
}

#[test]
fn a_mis_reported_reference_temperature_refuses() {
    // The watt gate cannot see a reference-wiring error smaller than hA times
    // its tolerance; the kelvin-denominated cross-check catches it directly.
    let (path, regions) = lumped_fixture(3);
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &ConjugateConfig::default(), |_, references| {
            let mut states = lumped_solid(&regions, 30.0, references);
            states[1].mean_reference_temperature_k = Some(references[1] + 5.0);
            Ok(states)
        })
        .expect_err("a solid that solved a different reference must refuse")
    });
    match error {
        AirflowError::ReferenceTemperatureMismatch {
            index,
            region,
            sent_bits,
            reported_bits,
        } => {
            assert_eq!(index, 1);
            assert_eq!(region, path.segments()[1].region());
            let drift = f64::from_bits(reported_bits) - f64::from_bits(sent_bits);
            assert!((drift - 5.0).abs() < 1.0e-9, "drift {drift}");
        }
        other => panic!("expected a reference mismatch, got {other:?}"),
    }

    // A response that makes no claim is admitted: `None` is a documented
    // no-claim boundary, not a wiring proof.
    let solution = with_cx(|cx| {
        solve_conjugate(cx, &path, &ConjugateConfig::default(), |_, references| {
            let mut states = lumped_solid(&regions, 30.0, references);
            for state in &mut states {
                state.mean_reference_temperature_k = None;
            }
            Ok(states)
        })
        .expect("a no-claim response is admitted")
    });
    assert_eq!(solution.balance.regions.len(), 3);
}

#[test]
fn a_non_finite_solid_heat_rate_cannot_reach_a_closed_balance() {
    // `f64::max` SILENTLY DROPS NaN, so a plain fold over per-region
    // imbalances would report 0 W -- a perfectly closed interface -- for a
    // coupling that blew up. The response gate stops it first, and the fold
    // poisons as a second line.
    let (path, regions) = lumped_fixture(3);
    let config = ConjugateConfig::default();
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            let mut states = lumped_solid(&regions, 30.0, references);
            states[1].heat_rate_w = f64::NAN;
            Ok(states)
        })
        .expect_err("a NaN heat rate must refuse")
    });
    assert!(
        matches!(
            error,
            AirflowError::InvalidConjugateInput {
                field: "solid region heat rate",
                ..
            }
        ),
        "got {error:?}"
    );
}

#[test]
fn a_solid_area_that_disagrees_with_the_declared_path_refuses() {
    // Right region, right count, wrong surface: the two sides would exchange
    // heat across different areas and every watt downstream is mis-scaled.
    let (path, regions) = lumped_fixture(3);
    let config = ConjugateConfig::default();
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            let mut states = lumped_solid(&regions, 30.0, references);
            states[0].area_m2 *= 1.5;
            Ok(states)
        })
        .expect_err("an area mismatch must refuse")
    });
    assert!(
        matches!(error, AirflowError::SegmentAreaMismatch { index: 0, .. }),
        "got {error:?}"
    );
}

#[test]
fn resuming_from_an_actual_cancellation_completes_the_exchange() {
    // The drill that matters: cancel for real, then resume from the refusal's
    // own payload -- not from a history the caller might not have kept.
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let uninterrupted = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("converges")
    });

    let gate = CancelGate::new();
    let mut seen = 0_usize;
    let error = with_gated_cx(&gate, |cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            seen += 1;
            if seen == 3 {
                gate.request();
            }
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect_err("cancelled")
    });
    let AirflowError::Cancelled { references_k, .. } = error else {
        panic!("expected cancellation")
    };

    let resumed = with_cx(|cx| {
        solve_conjugate_from(cx, &path, &config, &references_k, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("resumes from the cancellation payload")
    });
    for (a, b) in resumed
        .reference_temperatures_k
        .iter()
        .zip(&uninterrupted.reference_temperatures_k)
    {
        assert_eq!(a.to_bits(), b.to_bits(), "resumed answer drifted");
    }
}

#[test]
fn the_decomposition_cross_check_catches_an_unowned_robin_face() {
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig::default();
    let solution = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("converges")
    });

    // A whole-domain Robin total that matches the declared decomposition.
    let clean = solution
        .clone()
        .with_decomposition_cross_check(solution.balance.solid_total_w);
    assert!(
        clean
            .balance
            .decomposition_residual_w
            .expect("present")
            .abs()
            < 1.0e-9
    );

    // A whole-domain total carrying 3 W the declared regions never claimed:
    // some Robin face is outside every declared region.
    let dirty = solution
        .clone()
        .with_decomposition_cross_check(solution.balance.solid_total_w + 3.0);
    assert!((dirty.balance.decomposition_residual_w.expect("present") + 3.0).abs() < 1.0e-9);

    // This method predates the structured refusal surface and is infallible.
    // Invalid arithmetic must therefore poison its optional diagnostic rather
    // than leave a citable infinity that looks like an ordinary residual.
    let mut finite_overflow = solution;
    finite_overflow.balance.solid_total_w = f64::MAX;
    let poisoned = finite_overflow.with_decomposition_cross_check(-f64::MAX);
    assert!(
        poisoned
            .balance
            .decomposition_residual_w
            .expect("present")
            .is_nan()
    );
}

// ---------------------------------------------------------------------------
// Port wiring and structural refusals
// ---------------------------------------------------------------------------

#[test]
fn the_seam_port_effort_is_a_temperature() {
    // The coupling variable is a temperature, and the fs-couple port kind this
    // seam declares must agree — otherwise the "port" is decoration.
    assert_eq!(seam_effort_dimensions(), Temperature::DIMS);
}

#[test]
fn a_path_refuses_structurally_invalid_declarations() {
    assert_eq!(
        AirPath::new(INLET_K, MASS_FLOW_KG_S, CP_J_KG_K, Vec::new()).expect_err("empty"),
        AirflowError::EmptyAirPath
    );

    let twice = vec![
        AirSegment::new("duct", 0.01, 50.0).expect("valid"),
        AirSegment::new("duct", 0.01, 50.0).expect("valid"),
    ];
    assert_eq!(
        AirPath::new(INLET_K, MASS_FLOW_KG_S, CP_J_KG_K, twice).expect_err("duplicate"),
        AirflowError::DuplicateAirSegment {
            region: "duct".to_string()
        }
    );

    assert!(matches!(
        AirSegment::new("duct", -1.0, 50.0).expect_err("negative area"),
        AirflowError::InvalidConjugateInput { .. }
    ));
    assert!(matches!(
        AirSegment::new("duct", 0.01, f64::NAN).expect_err("NaN htc"),
        AirflowError::InvalidConjugateInput { .. }
    ));
    assert_eq!(
        AirSegment::new("  ", 0.01, 50.0).expect_err("blank region"),
        AirflowError::EmptyElementName
    );

    let path = uniform_path(2, 0.02, 40.0);
    assert_eq!(
        path.march(&[300.0]).expect_err("arity"),
        AirflowError::SolidResponseArity {
            expected: 2,
            found: 1
        }
    );
}

#[test]
fn derived_path_quantities_refuse_at_the_exact_admission_stage() {
    let ordinary = AirSegment::new("ordinary", 1.0, 1.0).expect("segment");
    let capacity_underflow = AirPath::new(
        INLET_K,
        f64::MIN_POSITIVE,
        f64::MIN_POSITIVE,
        vec![ordinary.clone()],
    )
    .expect_err("positive factors whose product is zero must refuse");
    assert_invalid_field(capacity_underflow, "air capacity rate");

    let capacity_overflow = AirPath::new(INLET_K, f64::MAX, 2.0, vec![ordinary])
        .expect_err("positive finite factors whose product is infinite must refuse");
    assert_invalid_field(capacity_overflow, "air capacity rate");

    let conductance_underflow = AirSegment::new("tiny-conductance", 1.0, f64::from_bits(1))
        .expect("the smallest positive conductance is representable");
    let ntu_underflow = AirPath::new(INLET_K, 1.0, f64::MAX, vec![conductance_underflow])
        .expect_err("a positive conductance/capacity ratio rounded to zero must refuse");
    assert_invalid_field(ntu_underflow, "air segment NTU");

    let conductance_overflow = AirSegment::new("conductance-overflow", f64::MAX, 2.0)
        .expect_err("an unrepresentable hA product must refuse");
    assert_invalid_field(conductance_overflow, "air segment conductance");

    let huge_conductance =
        AirSegment::new("huge-conductance", 1.0, f64::MAX).expect("finite conductance");
    let ntu_overflow = AirPath::new(INLET_K, 1.0, f64::from_bits(1), vec![huge_conductance])
        .expect_err("an unrepresentable conductance/capacity ratio must refuse");
    assert_invalid_field(ntu_overflow, "air segment NTU");

    let wide_a = AirSegment::new("wide-a", f64::MAX, 1.0e-308).expect("finite conductance");
    let wide_b = AirSegment::new("wide-b", f64::MAX, 1.0e-308).expect("finite conductance");
    let total_area_overflow = AirPath::new(INLET_K, 1.0, 1.0, vec![wide_a, wide_b])
        .expect_err("an unrepresentable total wetted area must refuse");
    assert_invalid_field(total_area_overflow, "total air segment area");
}

#[test]
fn march_refuses_each_non_finite_heat_aggregate_at_its_producing_stage() {
    let ordinary = AirPath::new(
        f64::MAX,
        1.0,
        1.0,
        vec![AirSegment::new("ordinary", 1.0, 1.0).expect("segment")],
    )
    .expect("path");
    let driving_overflow = ordinary
        .march(&[-f64::MAX])
        .expect_err("the wall/inlet difference must be checked before products");
    assert_non_finite_stage(driving_overflow, "segment temperature difference");

    let large_segment = AirSegment::new("large", 1.0, 1.0e308).expect("finite conductance");
    let large_path = AirPath::new(INLET_K, 1.0e154, 1.0e154, vec![large_segment]).expect("path");
    let segment_heat_overflow = large_path
        .march(&[INLET_K + 2.0])
        .expect_err("an unrepresentable segment heat rate must refuse");
    assert_non_finite_stage(segment_heat_overflow, "segment heat rate");

    let two_large = vec![
        AirSegment::new("large-a", 1.0, 1.0e308).expect("segment"),
        AirSegment::new("large-b", 1.0, 1.0e308).expect("segment"),
    ];
    let aggregate_path = AirPath::new(INLET_K, 1.0e154, 1.0e154, two_large).expect("path");
    let target_driving = 1.45;
    let independent_effectiveness = 1.0 - (-1.0_f64).exp();
    let walls = [
        INLET_K + target_driving,
        INLET_K + target_driving * independent_effectiveness + target_driving,
    ];
    let aggregate_heat_overflow = aggregate_path
        .march(&walls)
        .expect_err("finite row heat rates with an unrepresentable sum must refuse");
    assert_non_finite_stage(aggregate_heat_overflow, "total air heat rate");
}

#[test]
fn coupled_aggregate_overflows_keep_their_exact_attribution() {
    let (path, regions) = lumped_fixture(2);
    let solid_total_overflow = with_cx(|cx| {
        solve_conjugate(cx, &path, &ConjugateConfig::default(), |_, references| {
            let mut states = lumped_solid(&regions, 10.0, references);
            for state in &mut states {
                state.heat_rate_w = f64::MAX;
            }
            Ok(states)
        })
        .expect_err("finite solid rows with an unrepresentable sum must refuse")
    });
    assert_non_finite_stage(solid_total_overflow, "total solid heat rate");

    let interface_path = AirPath::new(
        INLET_K,
        1.0e154,
        1.0e154,
        vec![AirSegment::new("interface", 1.0, 1.0e308).expect("segment")],
    )
    .expect("path");
    let interface_overflow = with_cx(|cx| {
        solve_conjugate(
            cx,
            &interface_path,
            &ConjugateConfig::default(),
            |_, references| {
                Ok(vec![SolidRegionState {
                    region: "interface".to_string(),
                    area_m2: 1.0,
                    mean_wall_temperature_k: INLET_K - 1.5,
                    heat_rate_w: f64::MAX,
                    mean_reference_temperature_k: Some(references[0]),
                }])
            },
        )
        .expect_err("an unrepresentable air/solid total difference must refuse")
    });
    assert_non_finite_stage(
        interface_overflow,
        "air-solid interface heat-rate imbalance",
    );
}

#[test]
fn area_weighted_residual_overflows_keep_term_and_accumulator_attribution() {
    let contribution_path = AirPath::new(
        INLET_K,
        1.0,
        1.0,
        vec![AirSegment::new("wide", 1.0e308, 1.0e-308).expect("segment")],
    )
    .expect("path");
    let contribution_overflow = with_cx(|cx| {
        solve_conjugate(
            cx,
            &contribution_path,
            &ConjugateConfig::default(),
            |_, references| {
                Ok(vec![SolidRegionState {
                    region: "wide".to_string(),
                    area_m2: 1.0e308,
                    mean_wall_temperature_k: INLET_K + 5.0,
                    heat_rate_w: 0.0,
                    mean_reference_temperature_k: Some(references[0]),
                }])
            },
        )
        .expect_err("an unrepresentable area times residual must refuse")
    });
    assert_non_finite_stage(
        contribution_overflow,
        "area-weighted reference residual contribution",
    );

    let accumulation_path = AirPath::new(
        INLET_K,
        1.0,
        1.0,
        vec![
            AirSegment::new("wide-a", 4.0e307, 2.5e-308).expect("segment"),
            AirSegment::new("wide-b", 4.0e307, 2.5e-308).expect("segment"),
        ],
    )
    .expect("path");
    let accumulation_overflow = with_cx(|cx| {
        solve_conjugate(
            cx,
            &accumulation_path,
            &ConjugateConfig::default(),
            |_, references| {
                Ok(accumulation_path
                    .segments()
                    .iter()
                    .zip(references)
                    .map(|(segment, reference)| SolidRegionState {
                        region: segment.region().to_string(),
                        area_m2: segment.area_m2(),
                        mean_wall_temperature_k: INLET_K + 5.0,
                        heat_rate_w: 0.0,
                        mean_reference_temperature_k: Some(*reference),
                    })
                    .collect())
            },
        )
        .expect_err("finite weighted terms with an unrepresentable sum must refuse")
    });
    assert_non_finite_stage(
        accumulation_overflow,
        "area-weighted reference residual accumulation",
    );
}

#[test]
fn converged_region_differences_refuse_before_a_balance_is_published() {
    let capacity = 1.0e308;
    let segments = vec![
        AirSegment::new("row-a", 1.0, capacity).expect("segment"),
        AirSegment::new("row-b", 1.0, capacity).expect("segment"),
    ];
    let path = AirPath::new(INLET_K, 1.0e154, 1.0e154, segments).expect("path");
    let independent_effectiveness = 1.0 - (-1.0_f64).exp();
    let first_wall = INLET_K - 1.5;
    let second_wall = INLET_K - 1.5 * independent_effectiveness + 1.5;
    let config = ConjugateConfig {
        temperature_tolerance_k: f64::MAX,
        ..ConjugateConfig::default()
    };
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(vec![
                SolidRegionState {
                    region: "row-a".to_string(),
                    area_m2: 1.0,
                    mean_wall_temperature_k: first_wall,
                    heat_rate_w: f64::MAX,
                    mean_reference_temperature_k: Some(references[0]),
                },
                SolidRegionState {
                    region: "row-b".to_string(),
                    area_m2: 1.0,
                    mean_wall_temperature_k: second_wall,
                    heat_rate_w: -f64::MAX,
                    mean_reference_temperature_k: Some(references[1]),
                },
            ])
        })
        .expect_err("an unrepresentable per-region difference must refuse")
    });
    assert_non_finite_stage(error, "converged region interface heat-rate imbalance");
}

#[test]
fn a_solid_response_out_of_declared_order_refuses() {
    // The wiring fault the type system cannot catch: right count, wrong
    // association. Silently accepting it would attribute one surface's heat to
    // another surface's air.
    let (path, regions) = lumped_fixture(3);
    let config = ConjugateConfig::default();
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            let mut states = lumped_solid(&regions, 30.0, references);
            states.swap(0, 2);
            Ok(states)
        })
        .expect_err("swapped regions must refuse")
    });
    assert_eq!(
        error,
        AirflowError::SegmentRegionMismatch {
            index: 0,
            expected: "channel-0".to_string(),
            found: "channel-2".to_string(),
        }
    );
}

#[test]
fn a_zero_iteration_budget_refuses() {
    let (path, regions) = lumped_fixture(2);
    let config = ConjugateConfig {
        max_iterations: 0,
        ..ConjugateConfig::default()
    };
    let error = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, 10.0, references))
        })
        .expect_err("zero iterations is not a stop rule")
    });
    assert!(matches!(
        error,
        AirflowError::InvalidConjugateInput {
            field: "conjugate max iterations",
            ..
        }
    ));
}

// ---------------------------------------------------------------------------
// The CHT ladder's real correlation-rung transfer
// ---------------------------------------------------------------------------

#[test]
fn the_refinement_transfer_round_trips_on_the_coarse_space() {
    // The declared approximation property, same shape as Refine1d's — but see
    // the next test for the one Refine1d cannot state.
    let path = uniform_path(3, 0.03, 70.0);
    let transfer = SegmentRefinementTransfer::new(path, 4).expect("valid factor");
    let coarse = vec![340.0, 350.0, 360.0];
    let fine = transfer.prolongate(&coarse);
    assert_eq!(fine.len(), 12);
    for (restricted, original) in transfer.restrict(&fine).iter().zip(&coarse) {
        assert!((restricted - original).abs() < 1.0e-9);
    }
}

#[test]
fn the_refinement_transfer_preserves_the_air_side_qoi_exactly() {
    // This is the property that makes the transfer REAL rather than a generic
    // vector operation. Prolongate a wall state onto the refined path, march
    // the air there, and the outlet temperature and total heat must match the
    // coarse march bit for bit up to summation order — because
    // e^(−NTU) = (e^(−NTU/N))^N. Refine1d moves anonymous scalars and can make
    // no such statement.
    let coarse_path = uniform_path(3, 0.03, 70.0);
    let walls = vec![340.0, 350.0, 360.0];
    let coarse_march = coarse_path.march(&walls).expect("coarse march");

    for &factor in &[2_usize, 5, 16] {
        let transfer =
            SegmentRefinementTransfer::new(coarse_path.clone(), factor).expect("valid factor");
        let fine_path = transfer.refined_path().expect("refined path");
        let fine_walls = transfer.prolongate(&walls);
        let fine_march = fine_path.march(&fine_walls).expect("fine march");

        assert_eq!(fine_path.segments().len(), 3 * factor);
        assert!(
            (fine_march.outlet_temperature_k - coarse_march.outlet_temperature_k).abs() < 1.0e-10,
            "factor {factor}: outlet {} vs coarse {}",
            fine_march.outlet_temperature_k,
            coarse_march.outlet_temperature_k
        );
        assert!(
            (fine_march.total_heat_rate_w - coarse_march.total_heat_rate_w).abs() < 1.0e-8,
            "factor {factor}: heat {} vs coarse {}",
            fine_march.total_heat_rate_w,
            coarse_march.total_heat_rate_w
        );

        // Restricting the refined wall state returns the coarse one, so the
        // round trip is closed on the physical state too. This is an identity
        // in REAL arithmetic; the restriction evaluates a ratio of geometric
        // sums, so in f64 it lands within a few ULP rather than bit-exact.
        for (restricted, original) in transfer.restrict(&fine_walls).iter().zip(&walls) {
            assert!(
                (restricted - original).abs() < 1.0e-9,
                "factor {factor}: restricted {restricted} vs {original}"
            );
        }
    }
}

#[test]
fn the_cht_ladder_carries_the_injected_transfer() {
    let path = uniform_path(3, 0.03, 70.0);
    let registry = LadderRegistry::cht_with(
        Box::new(SegmentRefinementTransfer::new(path.clone(), 4).expect("valid")),
        Box::new(SegmentRefinementTransfer::new(path, 2).expect("valid")),
    );
    let ladder = registry.ladder("cht").expect("registered");
    assert_eq!(ladder.len(), 3);
    assert_eq!(ladder.bottom().name, "correlation-Nu");

    // Rung 0 -> 1 now runs the conjugate refinement, not Refine1d. Refine1d
    // would map 3 coarse values to 5 fine ones (2n-1 linear interpolation);
    // the conjugate transfer maps them to 3*4 = 12 sub-segments.
    //
    // WHAT THIS DOES NOT SHOW: the edges are NAMED correlation->RANS and
    // RANS->LES because those are the rungs the CHT ladder declares, and
    // installing a correlation self-refinement on them does not make either
    // fine side a RANS or LES state. No RANS rung exists (bead f85xj.5.8).
    // This asserts the INJECTION SEAM carries a caller-supplied operator
    // instead of the demonstrator; the real cross-fidelity transfer is
    // tracked as the remainder of clause (c).
    let fine = ladder
        .prolongate(0, &[340.0, 350.0, 360.0])
        .expect("rung 0");
    assert_eq!(fine.len(), 12);
    for (restricted, original) in ladder
        .restrict(1, &fine)
        .expect("rung 1")
        .iter()
        .zip(&[340.0, 350.0, 360.0])
    {
        assert!((restricted - original).abs() < 1.0e-9);
    }
}

#[test]
fn the_refinement_transfer_restricts_a_nonuniform_state_outlet_preserving() {
    // The round-trip test above only ever restricts states `prolongate`
    // emitted, whose blocks are UNIFORM -- and a plain mean is correct exactly
    // there. That makes the round trip structurally blind to how a genuinely
    // refined fine state comes back. Air entering a sub-segment has already
    // been warmed upstream, so a hot sub-segment near the block outlet moves
    // the exit temperature more than an equally hot one near the inlet.
    let factor = 4;
    let coarse_path = uniform_path(1, 0.03, 80.0);
    let transfer =
        SegmentRefinementTransfer::new(coarse_path.clone(), factor).expect("valid factor");
    let fine_path = transfer.refined_path().expect("refined path");

    // Deliberately nonuniform: the hot end is downstream.
    let fine_walls = vec![310.0, 330.0, 360.0, 400.0];
    let fine_outlet = fine_path
        .march(&fine_walls)
        .expect("fine march")
        .outlet_temperature_k;

    let coarse_walls = transfer.restrict(&fine_walls);
    assert_eq!(coarse_walls.len(), 1);
    let coarse_outlet = coarse_path
        .march(&coarse_walls)
        .expect("coarse march")
        .outlet_temperature_k;
    assert!(
        (coarse_outlet - fine_outlet).abs() < 1.0e-9,
        "restriction did not preserve the outlet: coarse {coarse_outlet} vs fine {fine_outlet}"
    );

    // Sensitivity: the plain mean this replaced is wrong by tens of kelvin
    // here, so the check above falsifies a real alternative.
    let plain_mean = fine_walls.iter().sum::<f64>() / factor as f64;
    let plain_outlet = coarse_path
        .march(&[plain_mean])
        .expect("plain march")
        .outlet_temperature_k;
    assert!(
        (plain_outlet - fine_outlet).abs() > 1.0,
        "the plain mean should be visibly wrong here, saw {plain_outlet} vs {fine_outlet}"
    );

    // And it still reduces to the plain mean on a uniform block, so it does
    // not disturb the states `prolongate` produces.
    let uniform = vec![355.0; factor];
    assert!((transfer.restrict(&uniform)[0] - 355.0).abs() < 1.0e-9);
}

#[test]
fn a_mis_sized_transfer_state_poisons_instead_of_truncating() {
    // The infallible `Transfer` trait cannot refuse, so a state whose arity
    // disagrees with the bound path must come back all-NaN at the correct
    // output arity. Silent alternatives are worse: zip/chunks would drop the
    // surplus downstream entries (the ones restrict weights most heavily) or
    // restrict a partial block with full-block weights.
    let factor = 4;
    let transfer =
        SegmentRefinementTransfer::new(uniform_path(2, 0.03, 80.0), factor).expect("valid factor");

    // One entry short, one surplus, and an entirely wrong shape.
    for bad_fine in [vec![300.0; 7], vec![300.0; 9], vec![300.0; 1]] {
        let coarse = transfer.restrict(&bad_fine);
        assert_eq!(coarse.len(), 2);
        assert!(
            coarse.iter().all(|wall| wall.is_nan()),
            "mis-sized fine state ({} entries) must poison, saw {coarse:?}",
            bad_fine.len()
        );
    }
    for bad_coarse in [vec![300.0; 1], vec![300.0; 3]] {
        let fine = transfer.prolongate(&bad_coarse);
        assert_eq!(fine.len(), 2 * factor);
        assert!(
            fine.iter().all(|wall| wall.is_nan()),
            "mis-sized coarse state ({} entries) must poison, saw {fine:?}",
            bad_coarse.len()
        );
    }

    // Correct arities are untouched by the guard.
    let ok_fine = vec![340.0; 2 * factor];
    assert_eq!(transfer.prolongate(&[340.0, 350.0]).len(), 2 * factor);
    assert_eq!(transfer.restrict(&ok_fine).len(), 2);
}

#[test]
fn a_zero_refinement_factor_refuses() {
    let path = uniform_path(2, 0.02, 50.0);
    assert!(matches!(
        SegmentRefinementTransfer::new(path, 0).expect_err("zero factor"),
        AirflowError::InvalidConjugateInput {
            field: "segment refinement factor",
            ..
        }
    ));
}

// ---------------------------------------------------------------------------
// Pause / resume
// ---------------------------------------------------------------------------

#[test]
fn resuming_from_a_checkpoint_reproduces_the_uninterrupted_tail_bitwise() {
    // G4: every outer iteration is a checkpoint boundary and its recorded
    // reference vector is the complete carried state. Resuming from record k
    // must land on exactly the same answer, bit for bit, in exactly the
    // remaining number of iterations.
    let power_w = 45.0;
    let (path, regions) = lumped_fixture(4);
    let config = ConjugateConfig {
        temperature_tolerance_k: 1.0e-12,
        max_iterations: 200,
        balance_tolerance_w: 1.0e-6,
        relaxation: Relaxation::Fixed { omega: 1.0 },
        ..ConjugateConfig::default()
    };
    let full = with_cx(|cx| {
        solve_conjugate(cx, &path, &config, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("converges")
    });
    assert!(full.iterations > 6, "need a tail worth resuming");

    let checkpoint = 5_usize;
    let seed = full.history[checkpoint].reference_temperatures_k.clone();
    let resumed = with_cx(|cx| {
        solve_conjugate_from(cx, &path, &config, &seed, |_, references| {
            Ok(lumped_solid(&regions, power_w, references))
        })
        .expect("resumes")
    });

    assert_eq!(resumed.iterations, full.iterations - checkpoint);
    for (a, b) in resumed
        .reference_temperatures_k
        .iter()
        .zip(&full.reference_temperatures_k)
    {
        assert_eq!(a.to_bits(), b.to_bits(), "resumed answer drifted");
    }
    for (index, record) in resumed.history.iter().enumerate() {
        let original = &full.history[checkpoint + index];
        assert_eq!(
            record.max_reference_change_k.to_bits(),
            original.max_reference_change_k.to_bits(),
            "resumed residual history drifted at {index}"
        );
    }
}

#[test]
fn a_resume_vector_of_the_wrong_length_refuses() {
    let (path, regions) = lumped_fixture(3);
    let config = ConjugateConfig::default();
    let error = with_cx(|cx| {
        solve_conjugate_from(cx, &path, &config, &[300.0, 301.0], |_, references| {
            Ok(lumped_solid(&regions, 20.0, references))
        })
        .expect_err("wrong arity")
    });
    assert_eq!(
        error,
        AirflowError::SolidResponseArity {
            expected: 3,
            found: 2
        }
    );
}
