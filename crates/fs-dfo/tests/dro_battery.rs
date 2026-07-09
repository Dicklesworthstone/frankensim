//! Wasserstein-DRO inner-sup battery: exact kink recovery, large-radius
//! saturation, and fail-fast public contracts.

use fs_dfo::wasserstein_worst_case;
use std::panic::{AssertUnwindSafe, catch_unwind};

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-dfo-dro\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

#[test]
fn kink_recovers_fractional_distribution() {
    // One empirical point, two candidate support points. The budget can
    // move exactly half the mass to the high-loss point:
    // sup q*1 subject to 2q <= 1 is 0.5.
    let report = wasserstein_worst_case(&[0.0, 1.0], &[0.0, 2.0], 1, 1.0);
    assert!(
        (report.worst_case - 0.5).abs() < 1e-9,
        "worst case {}",
        report.worst_case
    );
    assert!(
        (report.q[0] - 0.5).abs() < 1e-8 && (report.q[1] - 0.5).abs() < 1e-8,
        "kink must split mass across both active supports: {:?}",
        report.q
    );
    let q_expectation = report.q[1];
    assert!(
        (q_expectation - report.worst_case).abs() < 1e-8,
        "reported q must realize the reported worst-case value"
    );
    log(
        "kink-fractional-q",
        "pass",
        &format!("lambda {:.6}, q {:?}", report.lambda, report.q),
    );
}

#[test]
fn tiny_scale_kink_uses_scale_relative_recovery() {
    let losses = [0.0, 1.0e-20];
    let report = wasserstein_worst_case(&losses, &[0.0, 1.0], 1, 0.25);
    assert!(
        (report.worst_case - 0.25e-20).abs() < 1e-30,
        "tiny-scale worst case {}",
        report.worst_case
    );
    assert!(
        (report.q[0] - 0.75).abs() < 1e-8 && (report.q[1] - 0.25).abs() < 1e-8,
        "tiny-scale kink must not saturate from an absolute lambda cutoff: {:?}",
        report.q
    );
    let q_expectation = losses
        .iter()
        .zip(&report.q)
        .map(|(loss, mass)| loss * mass)
        .sum::<f64>();
    assert!(
        (q_expectation - report.worst_case).abs() < 1e-30,
        "reported q must realize the reported tiny-scale worst-case value"
    );
    log(
        "tiny-scale-kink",
        "pass",
        &format!("lambda {:.6e}, q {:?}", report.lambda, report.q),
    );
}

#[test]
fn large_radius_saturates_at_max_loss() {
    let losses = [0.0, 3.0, 1.0];
    let costs = [
        0.0, 1.0, 2.0, //
        2.0, 0.0, 1.0,
    ];
    let report = wasserstein_worst_case(&losses, &costs, 2, 10.0);
    assert!(
        (report.worst_case - 3.0).abs() < 1e-12,
        "large radius should saturate at max loss"
    );
    assert!(
        (report.q[1] - 1.0).abs() < 1e-12,
        "all mass can move to the max-loss support: {:?}",
        report.q
    );
    log("large-radius", "pass", "worst case saturates at max loss");
}

#[test]
fn public_contract_guards_fail_fast() {
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = wasserstein_worst_case(&[], &[], 1, 0.0);
        }))
        .is_err(),
        "empty losses must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = wasserstein_worst_case(&[0.0], &[0.0], 0, 0.0);
        }))
        .is_err(),
        "zero empirical sample count must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = wasserstein_worst_case(&[0.0], &[0.0], 1, f64::NAN);
        }))
        .is_err(),
        "non-finite radius must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = wasserstein_worst_case(&[f64::NAN], &[0.0], 1, 0.0);
        }))
        .is_err(),
        "non-finite losses must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = wasserstein_worst_case(&[0.0], &[-1.0], 1, 0.0);
        }))
        .is_err(),
        "negative costs must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = wasserstein_worst_case(&[0.0, 1.0], &[0.1, 0.2], 1, 0.0);
        }))
        .is_err(),
        "each row needs a zero-cost stay-put support"
    );
    log(
        "contract-guards",
        "pass",
        "invalid losses, sample counts, radii, costs, and stay-put rows fail fast",
    );
}
