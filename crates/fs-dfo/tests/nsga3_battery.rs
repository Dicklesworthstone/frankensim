//! NSGA-III battery (vcia many-objective lane): Das–Dennis direction
//! counts + simplex membership; DTLZ2(3) convergence to the known
//! unit-sphere-octant front with reference-direction COVERAGE
//! (diversity measured by association); the many-objective claim at
//! m = 5 — NSGA-III beats NSGA-II on MC-estimated hypervolume at
//! matched budget; bitwise replay; golden.

use fs_dfo::{NsgaParams, das_dennis, hypervolume, mc_hypervolume, nsga2, nsga3};
use fs_obs::ident::ReplayIdentity;

const SUITE: &str = "fs-dfo-nsga3";
const FIXED_INPUT_SEED: u64 = 0;
const DTLZ2_M3_INPUT_SEED: u64 = 17;
const M5_OPT_INPUT_SEED: u64 = 23;
const M5_MC_INPUT_SEED: u64 = 99;
const GOLDEN_INPUT_SEED: u64 = 3;
const MOEAD_ZDT_INPUT_SEED: u64 = 29;
const MOEAD_DTLZ_INPUT_SEED: u64 = 31;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("NSGA-III verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("NSGA-III verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn measurement(case: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new(SUITE, format!("{case}/measurement"));
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: case.to_string(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("NSGA-III measurement must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("NSGA-III measurement must use the fs-obs wire schema");
    println!("{line}");
}

/// DTLZ2 with m objectives, n = m − 1 + k variables in [0,1].
fn dtlz2(x: &[f64], m: usize) -> Vec<f64> {
    let k = x.len() - (m - 1);
    let g: f64 = x[m - 1..]
        .iter()
        .map(|v| (v - 0.5) * (v - 0.5))
        .sum::<f64>();
    let _ = k;
    let half_pi = core::f64::consts::FRAC_PI_2;
    (0..m)
        .map(|i| {
            let mut f = 1.0 + g;
            for &xj in &x[..m - 1 - i] {
                f *= fs_math::det::cos(xj * half_pi);
            }
            if i > 0 {
                f *= fs_math::det::sin(x[m - 1 - i] * half_pi);
            }
            f
        })
        .collect()
}

#[test]
fn das_dennis_counts_and_simplex() {
    // C(p+m−1, m−1): m=3, p=12 → C(14,2) = 91; m=5, p=4 → C(8,4) = 70.
    let d3 = das_dennis(3, 12);
    assert_eq!(d3.len(), 91);
    let d5 = das_dennis(5, 4);
    assert_eq!(d5.len(), 70);
    for dir in d3.iter().chain(&d5) {
        let s: f64 = dir.iter().sum();
        assert!(
            (s - 1.0).abs() < 1e-12,
            "direction off the simplex: {dir:?}"
        );
        assert!(dir.iter().all(|&v| v >= 0.0));
    }
    verdict(
        "das-dennis",
        true,
        "91 @ (3,12), 70 @ (5,4), on-simplex",
        FIXED_INPUT_SEED,
    );
}

fn guard_nsga_params() -> NsgaParams {
    NsgaParams {
        pop: 4,
        generations: 0,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 0.25,
        seed: 101,
    }
}

#[test]
#[should_panic(expected = "Das-Dennis directions need at least one objective")]
fn das_dennis_rejects_zero_objectives() {
    let _ = das_dennis(0, 1);
}

#[test]
#[should_panic(expected = "Das-Dennis directions need at least one division")]
fn das_dennis_rejects_zero_divisions() {
    let _ = das_dennis(3, 0);
}

#[test]
#[should_panic(expected = "NSGA-III reference-direction dimension must match objective dimension")]
fn nsga3_rejects_reference_direction_dimension_mismatch() {
    let nsga_params = guard_nsga_params();
    let bad_dirs = vec![vec![1.0]];
    let mut f_bad_dir = |_: &[f64]| vec![0.0, 1.0];
    let _ = nsga3(&mut f_bad_dir, 2, (0.0, 1.0), &bad_dirs, &nsga_params);
}

#[test]
#[should_panic(expected = "NSGA-III reference directions vectors must be finite")]
fn nsga3_rejects_all_zero_reference_direction() {
    let zero_dir = vec![vec![0.0, 0.0]];
    let mut f_zero_dir = |_: &[f64]| vec![0.0, 1.0];
    let _ = nsga3(
        &mut f_zero_dir,
        2,
        (0.0, 1.0),
        &zero_dir,
        &guard_nsga_params(),
    );
}

fn guard_moead_params(neighbors: usize) -> fs_dfo::MoeadParams {
    fs_dfo::MoeadParams {
        neighbors,
        max_replace: 1,
        generations: 0,
        eta_c: 20.0,
        eta_m: 20.0,
        p_mut: 0.25,
        seed: 102,
    }
}

#[test]
#[should_panic(expected = "MOEA/D neighborhood size must be positive")]
fn moead_rejects_zero_neighborhood() {
    let weights = das_dennis(2, 1);
    let mut f_zero_neighbors = |_: &[f64]| vec![0.0, 1.0];
    let _ = fs_dfo::moead(
        &mut f_zero_neighbors,
        2,
        (0.0, 1.0),
        &weights,
        &guard_moead_params(0),
    );
}

#[test]
#[should_panic(expected = "MOEA/D weight vectors must not be empty")]
fn moead_rejects_empty_weights() {
    let empty_weights: Vec<Vec<f64>> = Vec::new();
    let mut f_empty_weights = |_: &[f64]| vec![0.0, 1.0];
    let _ = fs_dfo::moead(
        &mut f_empty_weights,
        2,
        (0.0, 1.0),
        &empty_weights,
        &guard_moead_params(1),
    );
}

#[test]
fn dtlz2_m3_convergence_and_coverage() {
    let m = 3usize;
    let dirs = das_dennis(m, 12);
    // Standard DTLZ2 budgets (~250 generations); 150 left the worst
    // straggler at 0.0515 against the 0.05 gate.
    let params = NsgaParams {
        pop: 92,
        generations: 260,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 1.0 / 7.0,
        seed: DTLZ2_M3_INPUT_SEED,
    };
    let mut f = |x: &[f64]| dtlz2(x, m);
    let front = nsga3(&mut f, 7, (0.0, 1.0), &dirs, &params);
    // Convergence: the true front is ‖f‖₂ = 1.
    let mut worst_norm = 0.0f64;
    for ind in &front {
        let n2: f64 = ind.f.iter().map(|v| v * v).sum();
        worst_norm = worst_norm.max((fs_math::det::sqrt(n2) - 1.0).abs());
    }
    assert!(
        worst_norm < 0.05,
        "DTLZ2 front not converged: worst | ||f||-1 | = {worst_norm:.4}"
    );
    // Coverage: fraction of reference directions holding an associate.
    let covered = {
        let mut hit = vec![false; dirs.len()];
        for ind in &front {
            let mut best = (0usize, f64::INFINITY);
            for (k, dir) in dirs.iter().enumerate() {
                let dd: f64 = dir.iter().map(|d| d * d).sum();
                let t: f64 = ind.f.iter().zip(dir).map(|(a, b)| a * b).sum::<f64>() / dd;
                let d2: f64 = ind
                    .f
                    .iter()
                    .zip(dir)
                    .map(|(a, b)| {
                        let r = t.mul_add(-b, *a);
                        r * r
                    })
                    .sum();
                if d2 < best.1 {
                    best = (k, d2);
                }
            }
            hit[best.0] = true;
        }
        hit.iter().filter(|&&h| h).count() as f64 / dirs.len() as f64
    };
    assert!(
        covered > 0.6,
        "reference-direction coverage too low: {covered:.2}"
    );
    verdict(
        "dtlz2-m3",
        true,
        &format!(
            "worst norm dev {worst_norm:.4}, coverage {covered:.2}, front {}; input seed \
             {DTLZ2_M3_INPUT_SEED}",
            front.len()
        ),
        DTLZ2_M3_INPUT_SEED,
    );
}

#[test]
fn many_objective_m5_beats_nsga2_on_hv() {
    let m = 5usize;
    let dirs = das_dennis(m, 4);
    let params = NsgaParams {
        pop: 70,
        generations: 120,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 1.0 / 9.0,
        seed: M5_OPT_INPUT_SEED,
    };
    let mut f3 = |x: &[f64]| dtlz2(x, m);
    let front3 = nsga3(&mut f3, 9, (0.0, 1.0), &dirs, &params);
    let mut f2 = |x: &[f64]| dtlz2(x, m);
    let front2 = nsga2(&mut f2, 9, (0.0, 1.0), &params);
    let reference = vec![1.5f64; m];
    let pts3: Vec<Vec<f64>> = front3.iter().map(|i| i.f.clone()).collect();
    let pts2: Vec<Vec<f64>> = front2.iter().map(|i| i.f.clone()).collect();
    let (hv3, _) = mc_hypervolume(&pts3, &reference, 200_000, M5_MC_INPUT_SEED);
    let (hv2, _) = mc_hypervolume(&pts2, &reference, 200_000, M5_MC_INPUT_SEED);
    assert!(
        hv3 > hv2,
        "NSGA-III should beat NSGA-II at m=5: {hv3:.4} vs {hv2:.4}"
    );
    // Bitwise replay of NSGA-III.
    let mut fr = |x: &[f64]| dtlz2(x, m);
    let ra = nsga3(&mut fr, 9, (0.0, 1.0), &dirs, &params);
    let mut fr2 = |x: &[f64]| dtlz2(x, m);
    let rb = nsga3(&mut fr2, 9, (0.0, 1.0), &dirs, &params);
    assert_eq!(ra.len(), rb.len());
    for (p, q) in ra.iter().zip(&rb) {
        assert!(
            p.f.iter()
                .zip(&q.f)
                .all(|(u, v)| u.to_bits() == v.to_bits())
        );
    }
    verdict(
        "m5-vs-nsga2",
        true,
        &format!(
            "HV nsga3 {hv3:.4} vs nsga2 {hv2:.4} at matched budget, replay bitwise; \
             optimizer input seed {M5_OPT_INPUT_SEED}, MC-HV input seed \
             {M5_MC_INPUT_SEED}; composite aggregate seed zero"
        ),
        FIXED_INPUT_SEED,
    );
}

// The v1 hash named the pre-extension maxima-normalized lane. V2 ratified the
// typed normalization policy. V3 also binds the index-neutral equal-rank
// mating policy and must be frozen only after the complete central selector
// union passes against one source snapshot.
const GOLDEN_HASH_V3: u64 = 0x0282_b7f4_f7eb_9210;

fn golden_feed_bytes(accumulator: &mut u64, bytes: &[u8]) {
    for &byte in bytes {
        *accumulator ^= u64::from(byte);
        *accumulator = accumulator.wrapping_mul(0x0000_0100_0000_01b3);
    }
}

fn golden_feed_u64(accumulator: &mut u64, value: u64) {
    golden_feed_bytes(accumulator, &value.to_le_bytes());
}

fn golden_feed_str(accumulator: &mut u64, value: &str) {
    golden_feed_u64(
        accumulator,
        u64::try_from(value.len()).expect("policy string length fits u64"),
    );
    golden_feed_bytes(accumulator, value.as_bytes());
}

fn golden_feed_policy_identity(accumulator: &mut u64, field: &str, identity: &ReplayIdentity) {
    golden_feed_str(accumulator, field);
    golden_feed_u64(accumulator, u64::from(identity.version()));
    golden_feed_u64(accumulator, identity.root());
}

#[test]
fn nsga3_golden_preimage_consumes_shared_normalization_policy_root() {
    let policy = fs_dfo::moo::NSGA3_NORMALIZATION_POLICY;
    let current = policy.replay_identity();
    let mut mutant_policy = policy;
    mutant_policy.span_floor *= 2.0;
    let mutant = mutant_policy.replay_identity();
    assert_ne!(current.root(), mutant.root());

    let mut current_accumulator = 0xcbf2_9ce4_8422_2325;
    golden_feed_policy_identity(&mut current_accumulator, "normalization-policy", &current);
    let mut mutant_accumulator = 0xcbf2_9ce4_8422_2325;
    golden_feed_policy_identity(&mut mutant_accumulator, "normalization-policy", &mutant);
    assert_ne!(
        current_accumulator, mutant_accumulator,
        "the retained golden preimage must consume the shared typed policy root"
    );
}

#[test]
fn nsga3_golden_preimage_consumes_shared_mating_policy_root() {
    let policy = fs_dfo::moo::NSGA3_MATING_POLICY;
    let current = policy.replay_identity();
    let mut mutant_policy = policy;
    mutant_policy.equal_rank_policy = "lower-live-index-mutant";
    let mutant = mutant_policy.replay_identity();
    assert_ne!(current.root(), mutant.root());

    let mut current_accumulator = 0xcbf2_9ce4_8422_2325;
    golden_feed_policy_identity(&mut current_accumulator, "mating-policy", &current);
    let mut mutant_accumulator = 0xcbf2_9ce4_8422_2325;
    golden_feed_policy_identity(&mut mutant_accumulator, "mating-policy", &mutant);
    assert_ne!(
        current_accumulator, mutant_accumulator,
        "the retained golden preimage must consume the mating-policy root"
    );
}

#[test]
fn nsga3_golden_hash() {
    let normalization = fs_dfo::moo::NSGA3_NORMALIZATION_POLICY;
    let normalization_identity = normalization.replay_identity();
    let mating = fs_dfo::moo::NSGA3_MATING_POLICY;
    let mating_identity = mating.replay_identity();
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    golden_feed_str(&mut acc, "fs-dfo-nsga3-golden-v3");
    golden_feed_policy_identity(&mut acc, "normalization-policy", &normalization_identity);
    golden_feed_policy_identity(&mut acc, "mating-policy", &mating_identity);
    let dirs = das_dennis(3, 6);
    for d in dirs.iter().step_by(5) {
        for v in d {
            golden_feed_u64(&mut acc, v.to_bits());
        }
    }
    let params = NsgaParams {
        pop: 28,
        generations: 30,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 0.2,
        seed: GOLDEN_INPUT_SEED,
    };
    let mut f = |x: &[f64]| dtlz2(x, 3);
    let front = nsga3(&mut f, 5, (0.0, 1.0), &dirs, &params);
    for ind in front.iter().take(10) {
        for v in &ind.f {
            golden_feed_u64(&mut acc, v.to_bits());
        }
    }
    measurement(
        "nsga3-golden",
        format!(
            "{{\"identity_schema\":3,\"actual\":\"{acc:#018x}\",\
             \"expected\":\"{GOLDEN_HASH_V3:#018x}\",\
             \"input_seed\":{GOLDEN_INPUT_SEED},\"normalization_variant\":\"{}\",\
             \"normalization_policy_schema\":{},\"normalization_identity_version\":{},\
             \"normalization_identity_root\":\"0x{:016x}\",\"mating_variant\":\"{}\",\
             \"mating_policy_schema\":{},\"mating_identity_version\":{},\
             \"mating_identity_root\":\"0x{:016x}\"}}",
            normalization.variant,
            normalization.schema_version,
            normalization_identity.version(),
            normalization_identity.root(),
            mating.variant,
            mating.schema_version,
            mating_identity.version(),
            mating_identity.root(),
        ),
    );
    assert_eq!(
        acc, GOLDEN_HASH_V3,
        "nsga3 bits changed: {acc:#018x} vs {GOLDEN_HASH_V3:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
fn moead_zdt1_and_dtlz2_competitive() {
    use fs_dfo::{MoeadParams, moead};
    // ZDT1 convergence + spread (the NSGA-II gates, decomposition path).
    fn zdt1(x: &[f64]) -> Vec<f64> {
        let f1 = x[0];
        let g = 1.0 + 9.0 * x[1..].iter().sum::<f64>() / (x.len() - 1) as f64;
        vec![f1, g * (1.0 - fs_math::det::sqrt(f1 / g))]
    }
    let weights2 = das_dennis(2, 79); // 80 subproblems
    let params = MoeadParams {
        neighbors: 12,
        max_replace: 2,
        generations: 220,
        eta_c: 20.0,
        eta_m: 20.0,
        p_mut: 1.0 / 8.0,
        seed: MOEAD_ZDT_INPUT_SEED,
    };
    let mut f = |x: &[f64]| zdt1(x);
    let front = moead(&mut f, 8, (0.0, 1.0), &weights2, &params);
    let mean_gap: f64 = front
        .iter()
        .map(|ind| (ind.f[1] - (1.0 - fs_math::det::sqrt(ind.f[0]))).abs())
        .sum::<f64>()
        / front.len() as f64;
    assert!(mean_gap < 0.05, "MOEA/D ZDT1 not converged: {mean_gap:.4}");
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for ind in &front {
        lo = lo.min(ind.f[0]);
        hi = hi.max(ind.f[0]);
    }
    assert!(
        hi - lo > 0.7,
        "MOEA/D diversity collapsed: [{lo:.3},{hi:.3}]"
    );
    // DTLZ2(m=3): competitive with NSGA-III on hypervolume at matched
    // budget (within 10% — both are legitimate; numbers ledgered).
    let m = 3usize;
    let dirs = das_dennis(m, 12);
    let params3 = MoeadParams {
        neighbors: 15,
        max_replace: 2,
        generations: 260,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 1.0 / 7.0,
        seed: MOEAD_DTLZ_INPUT_SEED,
    };
    let mut fd = |x: &[f64]| dtlz2(x, m);
    let front_md = moead(&mut fd, 7, (0.0, 1.0), &dirs, &params3);
    let nsga_params = NsgaParams {
        pop: 92,
        generations: 260,
        eta_c: 30.0,
        eta_m: 20.0,
        p_mut: 1.0 / 7.0,
        seed: MOEAD_DTLZ_INPUT_SEED,
    };
    let mut fn3 = |x: &[f64]| dtlz2(x, m);
    let front_n3 = nsga3(&mut fn3, 7, (0.0, 1.0), &dirs, &nsga_params);
    let reference = vec![1.5f64; m];
    let pts_md: Vec<Vec<f64>> = front_md.iter().map(|i| i.f.clone()).collect();
    let pts_n3: Vec<Vec<f64>> = front_n3.iter().map(|i| i.f.clone()).collect();
    let hv_md = hypervolume(&pts_md, &reference);
    let hv_n3 = hypervolume(&pts_n3, &reference);
    assert!(
        hv_md > 0.9 * hv_n3,
        "MOEA/D should be competitive with NSGA-III: {hv_md:.4} vs {hv_n3:.4}"
    );
    // Bitwise replay.
    let mut fr = |x: &[f64]| zdt1(x);
    let ra = moead(&mut fr, 8, (0.0, 1.0), &weights2, &params);
    let mut fr2 = |x: &[f64]| zdt1(x);
    let rb = moead(&mut fr2, 8, (0.0, 1.0), &weights2, &params);
    assert_eq!(ra.len(), rb.len());
    for (p, q) in ra.iter().zip(&rb) {
        assert!(
            p.f.iter()
                .zip(&q.f)
                .all(|(u, v)| u.to_bits() == v.to_bits())
        );
    }
    verdict(
        "moead",
        true,
        &format!(
            "ZDT1 gap {mean_gap:.4} spread {:.2}; DTLZ2 HV {hv_md:.4} vs NSGA-III \
             {hv_n3:.4}; input seeds {MOEAD_ZDT_INPUT_SEED} and \
             {MOEAD_DTLZ_INPUT_SEED}; composite aggregate seed zero",
            hi - lo
        ),
        FIXED_INPUT_SEED,
    );
}
