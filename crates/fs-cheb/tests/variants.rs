//! fs-cheb variants battery (bead kw89): colleague-matrix roots
//! (including the even-multiplicity class the subdivision scanner
//! provably misses), interval-Newton certification, 2D low-rank cross
//! approximation, and the Fourier-periodic representation.

use fs_cheb::Cheb1;
use fs_cheb::budget::{ChebBudget, ChebError};
use fs_cheb::cheb2::Cheb2;
use fs_cheb::colleague::{
    ColleaguePolicy, RootPolynomialClass, certified_roots, certified_roots_bounded,
    certified_roots_with_checkpoint, colleague_roots,
};
use fs_cheb::fourier::FourierSeries;
use fs_ivl::{RootBox, RootSearchConfig, RootSearchTermination};
use fs_math::c64::C64;
use std::panic::{AssertUnwindSafe, catch_unwind};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

/// cheb-101: colleague agrees with the subdivision scanner on
/// simple-root fixtures (and with the analytic roots).
#[test]
fn cheb_101_colleague_simple_roots() {
    let p = Cheb1::build(&|x: f64| (x - 0.3) * (x + 0.62) * (x - 0.91), -1.0, 1.0, 16);
    let mut sub = p.roots();
    sub.sort_by(f64::total_cmp);
    let col = colleague_roots(&p, ColleaguePolicy::default());
    let want = [-0.62, 0.3, 0.91];
    let mut worst = 0.0f64;
    for (got, w) in col.iter().zip(&want) {
        worst = worst.max((got - w).abs());
    }
    let agree =
        col.len() == 3 && sub.len() == 3 && col.iter().zip(&sub).all(|(a, b)| (a - b).abs() < 1e-8);
    verdict(
        "cheb-101-simple-roots",
        agree && worst < 1e-10,
        &format!("colleague {col:?} vs analytic (worst dev {worst:.1e}); subdivision agrees"),
    );
}

/// cheb-102: EVEN-MULTIPLICITY recovery — (x−r)²(x−s) has no sign
/// change at r, so the scanner misses it (the documented v1 no-claim);
/// the colleague matrix finds it.
#[test]
fn cheb_102_even_multiplicity() {
    let (r, s) = (0.25f64, -0.7f64);
    let p = Cheb1::build(&|x: f64| (x - r) * (x - r) * (x - s), -1.0, 1.0, 16);
    let sub = p.roots();
    let scanner_misses = !sub.iter().any(|x| (x - r).abs() < 1e-6);
    let col = colleague_roots(&p, ColleaguePolicy::default());
    let colleague_finds =
        col.iter().any(|x| (x - r).abs() < 1e-6) && col.iter().any(|x| (x - s).abs() < 1e-8);
    verdict(
        "cheb-102-even-multiplicity",
        scanner_misses && colleague_finds && col.len() == 2,
        &format!(
            "scanner roots {sub:?} (double root at {r} MISSED as documented); colleague {col:?} finds both"
        ),
    );
}

/// cheb-103: clustered simple roots stay separated by the colleague
/// matrix down to 1e-3 spacing on a degree-24 interpolant.
#[test]
fn cheb_103_clustered_roots() {
    let roots = [0.100, 0.101, -0.4];
    let p = Cheb1::build(
        &|x: f64| roots.iter().map(|r| x - r).product::<f64>(),
        -1.0,
        1.0,
        24,
    );
    let col = colleague_roots(&p, ColleaguePolicy::default());
    let mut worst = f64::INFINITY;
    let found_all = roots.iter().all(|&r| {
        let d = col
            .iter()
            .map(|x| (x - r).abs())
            .fold(f64::INFINITY, f64::min);
        worst = worst.min(d);
        d < 1e-7
    });
    verdict(
        "cheb-103-clustered",
        found_all && col.len() == 3,
        &format!("colleague separates the 1e-3 cluster: {col:?}"),
    );
}

/// cheb-104: interval-Newton certification — every simple root comes
/// back CERTIFIED with a tight enclosure containing the analytic
/// root; the double-root fixture honestly reports Possible boxes (its
/// derivative encloses zero, as it must).
#[test]
fn cheb_104_certification() {
    let p = Cheb1::build(&|x: f64| (x - 0.3) * (x + 0.62) * (x - 0.91), -1.0, 1.0, 16);
    let boxes = certified_roots(&p, 1e-12);
    let want = [-0.62f64, 0.3, 0.91];
    let certified_boxes: Vec<_> = boxes
        .iter()
        .filter_map(|root| match root {
            RootBox::Certified(interval) => Some(interval),
            RootBox::Possible(_) => None,
        })
        .collect();
    let mut widths = 0.0f64;
    for interval in &certified_boxes {
        widths = widths.max(interval.hi() - interval.lo());
    }
    let simple_roots_distinct = want.iter().all(|&root| {
        certified_boxes
            .iter()
            .filter(|interval| interval.contains(root))
            .count()
            == 1
    });
    let simple_boxes_valid = boxes.iter().all(|root| {
        let interval = root.interval();
        interval.lo().is_finite()
            && interval.hi().is_finite()
            && interval.lo() >= -1.0
            && interval.hi() <= 1.0
    });
    // The double root: certification must NOT claim it.
    let pd = Cheb1::build(&|x: f64| (x - 0.25) * (x - 0.25) * (x + 0.7), -1.0, 1.0, 16);
    let dboxes = certified_roots(&pd, 1e-9);
    let double_certified = dboxes
        .iter()
        .any(|b| matches!(b, RootBox::Certified(iv) if iv.lo() <= 0.25 && 0.25 <= iv.hi()));
    let double_possible = dboxes
        .iter()
        .any(|b| matches!(b, RootBox::Possible(iv) if iv.contains(0.25)));
    let double_boxes_valid = dboxes.iter().all(|root| {
        let interval = root.interval();
        interval.lo().is_finite()
            && interval.hi().is_finite()
            && interval.lo() >= -1.0
            && interval.hi() <= 1.0
    });
    let shifted = Cheb1::build(&|x: f64| (x - 2.4) * (x - 3.25) * (x - 4.7), 2.0, 5.0, 32);
    let shifted_boxes = certified_roots(&shifted, 1e-12);
    let shifted_want = [2.4, 3.25, 4.7];
    let shifted_contains = shifted_want.iter().all(|&r| {
        shifted_boxes
            .iter()
            .any(|b| matches!(b, RootBox::Certified(iv) if iv.lo() <= r && r <= iv.hi()))
    });
    let shifted_in_domain = shifted_boxes
        .iter()
        .all(|b| b.interval().lo() >= 2.0 && b.interval().hi() <= 5.0);
    verdict(
        "cheb-104-certification",
        boxes.len() == 3
            && certified_boxes.len() == 3
            && simple_roots_distinct
            && simple_boxes_valid
            && widths < 1e-10
            && !double_certified
            && double_possible
            && double_boxes_valid
            && shifted_contains
            && shifted_in_domain,
        &format!(
            "3/3 simple roots CERTIFIED (max width {widths:.1e}); double root honestly uncertified ({} boxes); shifted-domain boxes {shifted_boxes:?}",
            dboxes.len(),
        ),
    );
}

#[test]
fn certified_search_preserves_exhaustion_and_physical_coverage() {
    let p = Cheb1::from_coeffs(2.0, 6.0, vec![0.0, 0.0, 1.0]);
    let config = RootSearchConfig {
        min_width: 1e-10,
        max_boxes: 1,
    };
    let partial = certified_roots_bounded(&p, config, &ChebBudget::default()).unwrap();
    assert!(!partial.domain_exhausted);
    assert_eq!(partial.boxes_examined, 1);
    assert_eq!(partial.reference_width_trajectory.len(), 1);
    assert!(partial.reference_width_trajectory[0] >= 2.0);
    assert!(partial.isolated_unique_roots().is_none());
    for root in [4.0 - 2.0_f64.sqrt(), 4.0 + 2.0_f64.sqrt()] {
        assert!(partial.roots.iter().any(|r| r.interval().contains(root)));
    }
    assert!(
        partial
            .roots
            .iter()
            .all(|r| r.interval().lo() >= 2.0 && r.interval().hi() <= 6.0)
    );
    let finished = certified_roots_bounded(
        &p,
        RootSearchConfig {
            max_boxes: 1024,
            ..config
        },
        &ChebBudget::default(),
    )
    .unwrap();
    assert!(finished.domain_exhausted);
    assert_eq!(
        finished.polynomial_class,
        RootPolynomialClass::Nonconstant { degree: 2 }
    );
    assert_eq!(finished.isolated_unique_roots().unwrap().len(), 2);
    assert_eq!(
        finished.boxes_examined,
        finished.reference_width_trajectory.len()
    );
    for root in [4.0 - 2.0_f64.sqrt(), 4.0 + 2.0_f64.sqrt()] {
        assert_eq!(
            finished
                .roots
                .iter()
                .filter(|r| r.interval().contains(root))
                .count(),
            1
        );
    }
}

#[test]
fn certified_search_classifies_zero_constant_and_repeated_roots() {
    let config = RootSearchConfig {
        min_width: 1e-6,
        max_boxes: 4096,
    };
    let budget = ChebBudget::default();
    let zero = Cheb1::from_coeffs(-1.0, 1.0, vec![0.0, -0.0, 0.0]);
    let report = certified_roots_bounded(&zero, config, &budget).unwrap();
    assert_eq!(
        report.polynomial_class,
        RootPolynomialClass::IdenticallyZero
    );
    assert_eq!(report.boxes_examined, 0);
    assert!(report.domain_exhausted);
    assert_eq!(report.roots, vec![RootBox::Possible(report.domain)]);
    assert!(report.isolated_unique_roots().is_none());
    // The stored c0/2 is positive in exact arithmetic even though a scalar
    // binary64 halving of this coefficient would underflow to zero.
    let constant = Cheb1::from_coeffs(-1.0, 1.0, vec![f64::from_bits(1), 0.0]);
    let report = certified_roots_bounded(&constant, config, &budget).unwrap();
    assert_eq!(
        report.polynomial_class,
        RootPolynomialClass::NonzeroConstant
    );
    assert_eq!(report.isolated_unique_roots(), Some([].as_slice()));
    let repeated = Cheb1::from_coeffs(-1.0, 1.0, vec![1.0, 0.0, 0.5]);
    let report = certified_roots_bounded(&repeated, config, &budget).unwrap();
    assert!(report.domain_exhausted);
    assert!(report.isolated_unique_roots().is_none());
    assert!(
        report
            .roots
            .iter()
            .any(|r| matches!(r, RootBox::Possible(iv) if iv.contains(0.0)))
    );
}

#[test]
fn certified_search_merges_unresolved_affine_clusters() {
    let a = 1.0;
    let b = f64::from_bits(1.0_f64.to_bits() + 1);
    let p = Cheb1::from_coeffs(a, b, vec![0.0, 0.0, 1.0]);
    let report = certified_roots_bounded(
        &p,
        RootSearchConfig {
            min_width: 1e-10,
            max_boxes: 1024,
        },
        &ChebBudget::default(),
    )
    .unwrap();
    assert!(report.domain_exhausted);
    assert_eq!(report.roots, vec![RootBox::Possible(report.domain)]);
    assert!(report.isolated_unique_roots().is_none());
}

#[test]
fn certified_search_cancellation_never_becomes_a_complete_root_list() {
    let p = Cheb1::from_coeffs(2.0, 6.0, vec![0.0, 0.0, 1.0]);
    let config = RootSearchConfig {
        min_width: 1e-10,
        max_boxes: 1024,
    };
    let budget = ChebBudget::default();
    let mut polls = 0;
    let complete = certified_roots_with_checkpoint(&p, config, &budget, &mut || {
        polls += 1;
        Ok::<(), usize>(())
    })
    .unwrap();
    assert_eq!(complete.termination, RootSearchTermination::Complete);
    for stop in 0..polls {
        let mut count = 0;
        let outcome = certified_roots_with_checkpoint(&p, config, &budget, &mut || {
            let current = count;
            count += 1;
            if current == stop { Err(stop) } else { Ok(()) }
        })
        .unwrap();
        assert_eq!(outcome.termination, RootSearchTermination::Cancelled(stop));
        assert!(!outcome.report.domain_exhausted);
        assert!(outcome.report.isolated_unique_roots().is_none());
        assert_eq!(
            outcome.report.boxes_examined,
            outcome.report.reference_width_trajectory.len()
        );
        for root in [4.0 - 2.0_f64.sqrt(), 4.0 + 2.0_f64.sqrt()] {
            assert!(
                outcome
                    .report
                    .roots
                    .iter()
                    .any(|r| r.interval().contains(root)),
                "root {root} lost at checkpoint {stop}: {:?}",
                outcome.report
            );
        }
    }
}

#[test]
fn certified_search_admission_refuses_before_polynomial_work() {
    let p = Cheb1::from_coeffs(-1.0, 1.0, vec![0.0, 1.0]);
    let config = RootSearchConfig {
        min_width: 1e-10,
        max_boxes: usize::MAX,
    };
    let result = certified_roots_with_checkpoint(
        &p,
        config,
        &ChebBudget::default(),
        &mut || -> Result<(), ()> { panic!("inadmissible work must not run") },
    );
    assert!(matches!(
        result,
        Err(ChebError::CapExceeded { .. } | ChebError::Overflow { .. })
    ));
    let wide = Cheb1::from_coeffs(-1.0, 1.0, vec![f64::from_bits(1), f64::MAX]);
    let result = certified_roots_bounded(
        &wide,
        RootSearchConfig {
            max_boxes: 32,
            ..config
        },
        &ChebBudget::default(),
    );
    assert!(matches!(result, Err(ChebError::Numerical { .. })));
}

/// cheb-105: 2D low-rank — separable functions captured at exact
/// rank; a non-separable smooth function converges to 1e-9 at modest
/// rank; the separable integral matches the analytic product.
#[test]
fn cheb_105_cheb2() {
    let dom = (-1.0, 1.0, -1.0, 1.0);
    let sep = Cheb2::build(
        &|x: f64, y: f64| x.sin() * (2.0 * y).cos(),
        dom,
        1e-12,
        8,
        64,
    );
    let sum2 = Cheb2::build(
        &|x: f64, y: f64| x.sin() * y.cos() + (x * x) * y,
        dom,
        1e-12,
        8,
        64,
    );
    let smooth = Cheb2::build(
        &|x: f64, y: f64| 1.0 / (1.0 + x * x + 2.0 * y * y),
        dom,
        1e-11,
        24,
        64,
    );
    let mut worst = 0.0f64;
    for i in 0..7 {
        for j in 0..7 {
            let x = -0.9 + 0.3 * f64::from(i);
            let y = -0.9 + 0.3 * f64::from(j);
            let want = 1.0 / (2.0f64.mul_add(y * y, 1.0 + x * x));
            worst = worst.max((smooth.eval(x, y) - want).abs());
        }
    }
    // ∫∫ sin(x)cos(2y) over [−1,1]² = 0 (odd in x).
    let int_sep = sep.integral();
    verdict(
        "cheb-105-cheb2",
        sep.rank() == 1
            && sum2.rank() == 2
            && smooth.rank() <= 24
            && worst < 1e-9
            && int_sep.abs() < 1e-12,
        &format!(
            "ranks: separable {} (exact 1), sum {} (exact 2), smooth {} (residual {:.1e}); eval dev {worst:.1e}; odd integral {int_sep:.1e}",
            sep.rank(),
            sum2.rank(),
            smooth.rank(),
            smooth.residual
        ),
    );
}

/// cheb-106: Fourier-periodic — exact trig recovery, the spectral
/// derivative identity d/dθ sin = cos pointwise, spectral tail decay
/// on a smooth periodic function, and the integral off c₀.
#[test]
fn cheb_106_fourier() {
    let f = FourierSeries::build(&|t: f64| 3.0f64.mul_add((2.0 * t).cos(), t.sin()), 32);
    let mut worst = 0.0f64;
    for k in 0..64 {
        let t = std::f64::consts::TAU * f64::from(k) / 64.0;
        let want = 3.0f64.mul_add((2.0 * t).cos(), t.sin());
        worst = worst.max((f.eval(t) - want).abs());
    }
    let d = FourierSeries::build(&|t: f64| t.sin(), 32).differentiate();
    let mut worst_d = 0.0f64;
    for k in 0..64 {
        let t = std::f64::consts::TAU * f64::from(k) / 64.0;
        worst_d = worst_d.max((d.eval(t) - t.cos()).abs());
    }
    let smooth = FourierSeries::build(&|t: f64| (t.sin()).exp(), 64);
    let tail = smooth.tail_magnitude(24);
    // ∫ exp(sin θ) dθ = 2π·I₀(1) ≈ 7.95492652101284.
    let want_int = 7.954_926_521_012_845;
    let int_dev = (smooth.integral() - want_int).abs();
    verdict(
        "cheb-106-fourier",
        worst < 1e-12 && worst_d < 1e-12 && tail < 1e-12 && int_dev < 1e-10,
        &format!(
            "trig recovery {worst:.1e}; derivative identity {worst_d:.1e}; tail(24+) {tail:.1e}; Bessel integral dev {int_dev:.1e}"
        ),
    );
}

/// cheb-107: determinism — colleague, cross, and Fourier replays are
/// bitwise identical.
#[test]
fn cheb_107_determinism() {
    let p = Cheb1::build(&|x: f64| (x - 0.3) * (x + 0.62) * (x - 0.91), -1.0, 1.0, 16);
    let a = colleague_roots(&p, ColleaguePolicy::default());
    let b = colleague_roots(&p, ColleaguePolicy::default());
    let dom = (-1.0, 1.0, -1.0, 1.0);
    let g = |x: f64, y: f64| 1.0 / (1.0 + x * x + y * y);
    let c1 = Cheb2::build(&g, dom, 1e-10, 16, 64);
    let c2 = Cheb2::build(&g, dom, 1e-10, 16, 64);
    let f1 = FourierSeries::build(&|t: f64| (t.sin()).exp(), 64);
    let f2 = FourierSeries::build(&|t: f64| (t.sin()).exp(), 64);
    let bitwise = a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits())
        && c1.eval(0.37, -0.21).to_bits() == c2.eval(0.37, -0.21).to_bits()
        && f1.eval(1.234).to_bits() == f2.eval(1.234).to_bits();
    verdict(
        "cheb-107-determinism",
        bitwise && a.len() == b.len() && c1.rank() == c2.rank(),
        "colleague, ACA, and Fourier replays bitwise identical",
    );
}

/// cheb-108: public contracts fail fast at the boundary instead of
/// letting malformed spectral objects emit nonsensical claims.
#[test]
fn cheb_108_contract_guards() {
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = Cheb1::from_coeffs(0.0, 1.0, vec![f64::NAN]);
        }))
        .is_err(),
        "non-finite Cheb1 coefficients must fail fast"
    );
    let p = Cheb1::build(&|x: f64| x - 0.25, -1.0, 1.0, 16);
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = colleague_roots(
                &p,
                ColleaguePolicy {
                    cluster_tol: 0.0,
                    ..ColleaguePolicy::default()
                },
            );
        }))
        .is_err(),
        "invalid colleague policy must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = certified_roots(&p, 0.0);
        }))
        .is_err(),
        "non-positive certification widths must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = Cheb2::build(&|_, _| 1.0, (0.0, 0.0, -1.0, 1.0), 1e-8, 1, 16);
        }))
        .is_err(),
        "degenerate Cheb2 domains must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = Cheb2::build(&|_, _| f64::NAN, (-1.0, 1.0, -1.0, 1.0), 1e-8, 1, 16);
        }))
        .is_err(),
        "non-finite Cheb2 samples must fail fast"
    );
    let bad_cheb2 = Cheb2 {
        cols: Vec::new(),
        rows: Vec::new(),
        inv_pivots: vec![1.0],
        residual: 0.0,
    };
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = bad_cheb2.eval(0.0, 0.0);
        }))
        .is_err(),
        "malformed hand-built Cheb2 objects must fail fast"
    );
    let zero_pivot_cheb2 = Cheb2 {
        cols: vec![Cheb1::from_coeffs(0.0, 1.0, vec![2.0])],
        rows: vec![Cheb1::from_coeffs(-1.0, 1.0, vec![2.0])],
        inv_pivots: vec![0.0],
        residual: 0.0,
    };
    assert!(
        catch_unwind(AssertUnwindSafe(|| zero_pivot_cheb2.rank())).is_err(),
        "zero inverse pivots must fail fast"
    );
    let mismatched_column_domains = Cheb2 {
        cols: vec![
            Cheb1::from_coeffs(0.0, 1.0, vec![2.0]),
            Cheb1::from_coeffs(0.0, 2.0, vec![2.0]),
        ],
        rows: vec![
            Cheb1::from_coeffs(-1.0, 1.0, vec![2.0]),
            Cheb1::from_coeffs(-1.0, 1.0, vec![2.0]),
        ],
        inv_pivots: vec![1.0, 1.0],
        residual: 0.0,
    };
    assert!(
        catch_unwind(AssertUnwindSafe(|| mismatched_column_domains.rank())).is_err(),
        "Cheb2 columns with different x-domains must fail fast"
    );
    let mismatched_row_domains = Cheb2 {
        cols: vec![
            Cheb1::from_coeffs(0.0, 1.0, vec![2.0]),
            Cheb1::from_coeffs(0.0, 1.0, vec![2.0]),
        ],
        rows: vec![
            Cheb1::from_coeffs(-1.0, 1.0, vec![2.0]),
            Cheb1::from_coeffs(-2.0, 1.0, vec![2.0]),
        ],
        inv_pivots: vec![1.0, 1.0],
        residual: 0.0,
    };
    assert!(
        catch_unwind(AssertUnwindSafe(|| mismatched_row_domains.rank())).is_err(),
        "Cheb2 rows with different y-domains must fail fast"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = FourierSeries::build(&|_| 0.0, 6);
        }))
        .is_err(),
        "Fourier sample counts must match the radix-2 backend contract"
    );
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = FourierSeries::build(&|_| f64::NAN, 8);
        }))
        .is_err(),
        "non-finite Fourier samples must fail fast"
    );
    let bad_fourier = FourierSeries {
        coeffs: vec![C64::ZERO],
        n: 8,
    };
    assert!(
        catch_unwind(AssertUnwindSafe(|| {
            let _ = bad_fourier.eval(0.0);
        }))
        .is_err(),
        "malformed hand-built FourierSeries objects must fail fast"
    );
    verdict(
        "cheb-108-contract-guards",
        true,
        "invalid policies, domains, samples, widths, and public spectral structs fail fast",
    );
}
