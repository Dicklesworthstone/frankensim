//! G0/G1/G2/G3/G5 complex-kernel checks, within the existing battery target.

use fs_fmm::{
    ComplexKernel, FmmC64, FmmError, HelmholtzAdjointDouble, HelmholtzDouble, HelmholtzGradient,
    HelmholtzHypersingular, HelmholtzSingle, MAX_BB_HELMHOLTZ_KD, MAX_COEFFICIENT_SLOTS,
    MAX_COMPLEX_COEFFICIENT_SLOTS, helmholtz_order_for,
};
use fs_math::{c64::C64, det};
use std::time::Instant;

const SEED: u64 = 0x1070_2026_1016_0001;
const TOL: f64 = 1.0e-6;
const ORDERS: [usize; 5] = [4, 6, 8, 10, 12];

// A small nonuniform cloud keeps the unfused p^6 M2L battery bounded. Three
// occupied leaves (x indices 0,2,3) activate both first-separated M2L pairs
// and distinct adjacent-leaf P2P pairs; all coordinates are dyadic. These
// are operator tests, not a large-cloud performance qualification.
fn fixture() -> (Vec<[f64; 3]>, Vec<C64>, Vec<[f64; 3]>) {
    let mut rng = super::Lcg(SEED);
    let mut points = Vec::new();
    for base in [0.0, 0.5, 0.9375] {
        for j in 0..8 {
            points.push([
                base + (j as f64 / 128.0),
                ((j * 3) % 8) as f64 / 128.0,
                ((j * 5) % 8) as f64 / 128.0,
            ]);
        }
    }
    points[23][0] = 1.0;
    let charges = (0..points.len())
        .map(|_| C64::new(rng.unit() - 0.5, rng.unit() - 0.5))
        .collect();
    let normals = (0..points.len())
        .map(|j| match j % 4 {
            0 => [0.6, 0.0, 0.8],
            1 => [0.0, -0.8, 0.6],
            2 => [1.0, 0.0, 0.0],
            _ => [0.0, 1.0, 0.0],
        })
        .collect();
    (points, charges, normals)
}

fn rel_l2(got: &[C64], oracle: &[C64]) -> f64 {
    let error = got
        .iter()
        .zip(oracle)
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum::<f64>();
    let scale = oracle.iter().map(|&z| z.norm_sq()).sum::<f64>();
    det::sqrt(error / scale)
}

fn bits(values: &[C64]) -> Vec<(u64, u64)> {
    values
        .iter()
        .map(|z| (z.re.to_bits(), z.im.to_bits()))
        .collect()
}

fn row(case: &str, axis: &str, rung: usize, kd: f64, err: f64, tol: f64, elapsed: f64, pass: bool) {
    let (p, oracle) = match axis {
        "scalar_then_double" => (8, "FMM on dyadically translated cloud"),
        "Cartesian_component" => (12, "direct complex point sum"),
        "bitwise_replay" => (
            rung,
            "independent reconstruction and repeated FMM bit patterns",
        ),
        _ => (rung, "direct complex point sum"),
    };
    super::measurement(
        case,
        "helmholtz_check",
        format!(
            "{{\"suite\":\"fs-fmm/battery\",\"case\":\"{case}\",\"inputs\":{{\"fixture\":\"dyadic-three-leaf-cloud-24\",\"provenance\":\"analytic outgoing Green kernel exp(ikr)/(4pi r)\",\"seed\":{SEED}}},\"axis\":\"{axis}\",\"rung\":{rung},\"n\":24,\"kd\":{kd},\"p\":{p},\"observed\":{},\"rel_l2\":{},\"oracle\":\"{oracle}\",\"err\":{},\"tol\":{tol},\"verdict\":\"{}\",\"elapsed_s\":{elapsed}}}",
            super::finite_json(err),
            super::finite_json(err),
            super::finite_json(err),
            if pass { "pass" } else { "fail" }
        ),
    );
}

#[test]
fn fmm_helmholtz_accuracy_vs_direct() {
    let (points, charges, _) = fixture();
    for kd in [1.0, 5.0, 10.0] {
        let kernel = HelmholtzSingle {
            k: kd / (1.0 + 1.0e-9),
        };
        let oracle = FmmC64::new(&kernel, points.clone(), 4, 2)
            .unwrap()
            .direct(&charges)
            .unwrap();
        let mut previous = f64::INFINITY;
        let mut pass = true;
        let mut final_error = previous;
        for p in ORDERS {
            let start = Instant::now();
            let fmm = FmmC64::new(&kernel, points.clone(), p, 2).unwrap();
            let got = fmm.potentials(&charges).unwrap();
            let err = rel_l2(&got, &oracle);
            let monotone = err < previous;
            pass &= monotone;
            row(
                "fmm-hz-single",
                "interpolation_order",
                p,
                kd,
                err,
                TOL,
                start.elapsed().as_secs_f64(),
                monotone && (p < 12 || err <= TOL),
            );
            previous = err;
            final_error = err;
        }
        super::verdict(
            "fmm-hz-single",
            pass && final_error <= TOL,
            "G1/G2: complex single-layer order sweep, kD in {1,5,10}",
            &format!("kD={kd}, final relative L2={final_error}; monotone={pass}"),
            SEED,
        );
        assert_eq!(helmholtz_order_for(kd, TOL).unwrap(), 12);
    }
}

#[test]
fn fmm_helmholtz_double_layer_vs_direct() {
    let (points, charges, normals) = fixture();
    for kd in [1.0, 5.0, 10.0] {
        let k = kd / (1.0 + 1.0e-9);
        let oracle =
            FmmC64::new_double(HelmholtzDouble { k }, points.clone(), normals.clone(), 4, 2)
                .unwrap()
                .direct(&charges)
                .unwrap();
        let mut previous = f64::INFINITY;
        let mut pass = true;
        for p in ORDERS {
            let start = Instant::now();
            let fmm =
                FmmC64::new_double(HelmholtzDouble { k }, points.clone(), normals.clone(), p, 2)
                    .unwrap();
            let err = rel_l2(&fmm.potentials(&charges).unwrap(), &oracle);
            pass &= err < previous;
            row(
                "fmm-hz-double",
                "interpolation_order",
                p,
                kd,
                err,
                TOL,
                start.elapsed().as_secs_f64(),
                err < previous && (p < 12 || err <= TOL),
            );
            previous = err;
        }
        super::verdict(
            "fmm-hz-double",
            pass && previous <= TOL,
            "G2: varying source normals factored through P2M, complex direct oracle",
            &format!("kD={kd}, final relative L2={previous}; monotone={pass}"),
            SEED,
        );
    }
}

#[test]
fn fmm_helmholtz_adjoint_gradient_and_hypersingular_vs_direct() {
    let (points, charges, normals) = fixture();
    let k = 10.0 / (1.0 + 1.0e-9);
    let adjoint = FmmC64::new_adjoint_double(
        HelmholtzAdjointDouble { k },
        points.clone(),
        normals.clone(),
        12,
        2,
    )
    .unwrap();
    let mixed =
        FmmC64::new_hypersingular(HelmholtzHypersingular { k }, points.clone(), normals, 12, 2)
            .unwrap();
    for (case, fmm) in [("fmm-hz-adjoint", adjoint), ("fmm-hz-mixed", mixed)] {
        let start = Instant::now();
        let err = rel_l2(
            &fmm.potentials(&charges).unwrap(),
            &fmm.direct(&charges).unwrap(),
        );
        row(
            case,
            "interpolation_order",
            12,
            10.0,
            err,
            TOL,
            start.elapsed().as_secs_f64(),
            err <= TOL,
        );
        assert!(err <= TOL, "{case}: {err}");
    }
    for axis in 0..3 {
        let kernel = HelmholtzGradient::new(k, axis).unwrap();
        let fmm = FmmC64::new(&kernel, points.clone(), 12, 2).unwrap();
        let start = Instant::now();
        let err = rel_l2(
            &fmm.potentials(&charges).unwrap(),
            &fmm.direct(&charges).unwrap(),
        );
        row(
            "fmm-hz-gradient",
            "Cartesian_component",
            axis,
            10.0,
            err,
            TOL,
            start.elapsed().as_secs_f64(),
            err <= TOL,
        );
        assert!(err <= TOL, "gradient {axis}: {err}");
    }
}

#[test]
fn fmm_helmholtz_refuses_beyond_bb_regime() {
    let start = Instant::now();
    let kd = 1000.0;
    let refusal = helmholtz_order_for(kd, TOL).unwrap_err();
    assert!(matches!(
        refusal,
        FmmError::HelmholtzRegime {
            route: "M2.5-W5b",
            ..
        }
    ));
    assert!(refusal.to_string().contains("M2.5-W5b"));
    let (points, charges, _) = fixture();
    let physical = HelmholtzSingle {
        k: kd / (1.0 + 1.0e-9),
    };
    assert!(matches!(
        FmmC64::new(&physical, points.clone(), 12, 2),
        Err(FmmError::HelmholtzRegime {
            route: "M2.5-W5b",
            ..
        })
    ));
    assert!(helmholtz_order_for(MAX_BB_HELMHOLTZ_KD, TOL).is_ok());
    assert!(helmholtz_order_for(MAX_BB_HELMHOLTZ_KD.next_up(), TOL).is_err());
    assert!(helmholtz_order_for(1.0, TOL.next_down()).is_err());
    // Test-only unguarded kernel: measure why p=12 cannot be assumed accurate
    // at this kD. Production HelmholtzSingle reports its k and refuses above.
    struct Diagnostic(HelmholtzSingle);
    impl ComplexKernel for Diagnostic {
        fn eval(&self, x: [f64; 3], y: [f64; 3]) -> C64 {
            self.0.eval(x, y)
        }
    }
    let raw = Diagnostic(physical);
    let fmm = FmmC64::new(&raw, points, 12, 2).unwrap();
    let err = rel_l2(
        &fmm.potentials(&charges).unwrap(),
        &fmm.direct(&charges).unwrap(),
    );
    super::measurement(
        "fmm-hz-refused",
        "helmholtz_regime_refusal",
        format!(
            "{{\"case\":\"fmm-hz-refused\",\"inputs\":{{\"fixture\":\"dyadic-three-leaf-cloud-24\",\"provenance\":\"analytic outgoing Green kernel\",\"seed\":{SEED}}},\"axis\":\"kD\",\"rung\":1000,\"observed\":{kd},\"err\":{},\"oracle\":\"test-only unguarded p12 vs direct\",\"tol\":{TOL},\"verdict\":\"refused\",\"route\":\"M2.5-W5b\",\"elapsed_s\":{}}}",
            super::finite_json(err),
            start.elapsed().as_secs_f64()
        ),
    );
    assert!(
        err > TOL,
        "the refused high-kD fixture must expose interpolation failure: {err}"
    );
}

#[test]
fn fmm_helmholtz_translation_invariance() {
    let (points, charges, normals) = fixture();
    let moved = points
        .iter()
        .map(|p| [p[0] + 17.25, p[1] - 4.5, p[2] + 9.75])
        .collect::<Vec<_>>();
    let k = 5.0;
    let kernel = HelmholtzSingle { k };
    let engines = [
        (
            FmmC64::new(&kernel, points.clone(), 8, 2).unwrap(),
            FmmC64::new(&kernel, moved.clone(), 8, 2).unwrap(),
        ),
        (
            FmmC64::new_double(HelmholtzDouble { k }, points, normals.clone(), 8, 2).unwrap(),
            FmmC64::new_double(HelmholtzDouble { k }, moved, normals, 8, 2).unwrap(),
        ),
    ];
    for (rung, (a, b)) in engines.into_iter().enumerate() {
        let start = Instant::now();
        let err = rel_l2(
            &a.potentials(&charges).unwrap(),
            &b.potentials(&charges).unwrap(),
        );
        row(
            "fmm-hz-translation",
            "scalar_then_double",
            rung,
            k * a.root_box_size(),
            err,
            1.0e-12,
            start.elapsed().as_secs_f64(),
            err < 1.0e-12,
        );
        assert!(err < 1.0e-12, "translation relative L2: {err}");
    }
}

#[test]
fn fmm_c64_caps() {
    let start = Instant::now();
    assert_eq!(2 * MAX_COMPLEX_COEFFICIENT_SLOTS, MAX_COEFFICIENT_SLOTS);
    FmmC64::validate_coefficient_slots(MAX_COMPLEX_COEFFICIENT_SLOTS).unwrap();
    let error = FmmC64::validate_coefficient_slots(MAX_COMPLEX_COEFFICIENT_SLOTS + 1).unwrap_err();
    assert!(
        matches!(error, FmmError::WorkEnvelopeExceeded { requested, max, .. }
        if requested == MAX_COEFFICIENT_SLOTS + 2 && max == MAX_COEFFICIENT_SLOTS)
    );
    assert!(FmmC64::validate_coefficient_slots(usize::MAX).is_err());
    super::measurement(
        "fmm-hz-caps",
        "complex_coefficient_cap",
        format!(
            "{{\"case\":\"fmm-hz-caps\",\"inputs\":{{\"provenance\":\"analytic sizeof(C64)=2*sizeof(f64)\"}},\"axis\":\"complex_coefficient_slots\",\"rung\":{},\"observed\":{},\"oracle\":{},\"err\":0,\"tol\":0,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
            MAX_COMPLEX_COEFFICIENT_SLOTS,
            MAX_COMPLEX_COEFFICIENT_SLOTS,
            MAX_COEFFICIENT_SLOTS / 2,
            start.elapsed().as_secs_f64()
        ),
    );
}

#[test]
fn fmm_c64_bitwise_repeat() {
    let start = Instant::now();
    let (points, charges, normals) = fixture();
    let kernel = HelmholtzSingle { k: 5.0 };
    let a = FmmC64::new(&kernel, points.clone(), 6, 2).unwrap();
    let b = FmmC64::new(&kernel, points.clone(), 6, 2).unwrap();
    assert_eq!(
        bits(&a.potentials(&charges).unwrap()),
        bits(&a.potentials(&charges).unwrap())
    );
    assert_eq!(
        bits(&a.potentials(&charges).unwrap()),
        bits(&b.potentials(&charges).unwrap())
    );
    let a = FmmC64::new_double(
        HelmholtzDouble { k: 5.0 },
        points.clone(),
        normals.clone(),
        6,
        2,
    )
    .unwrap();
    let b = FmmC64::new_double(HelmholtzDouble { k: 5.0 }, points, normals, 6, 2).unwrap();
    assert_eq!(
        bits(&a.potentials(&charges).unwrap()),
        bits(&b.potentials(&charges).unwrap())
    );
    row(
        "fmm-hz-repeat",
        "bitwise_replay",
        6,
        5.0,
        0.0,
        0.0,
        start.elapsed().as_secs_f64(),
        true,
    );
}

#[test]
fn laplace_bitwise_unchanged() {
    // Freeze the original Laplace kernel, interpolation and real engine's
    // operation tree at ee8e613806eaa9ac3979a8b139af5e1528b699d5. This is a
    // source guard, not an invented historical floating-point measurement.
    // The unmodified fmm-001..003 run alongside this test on the native host.
    let start = Instant::now();
    let original = include_str!("../../src/lib.rs")
        .split_once("impl Kernel for Laplace3d")
        .unwrap()
        .1;
    let fingerprint = original.bytes().fold(0xcbf2_9ce4_8422_2325u64, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3)
    });
    assert_eq!(
        fingerprint, 0x38b0_4e0d_5a5d_1fd2,
        "real operation tree changed; requires a deliberate golden event"
    );
    super::measurement(
        "fmm-laplace-preserved",
        "laplace_source_operation_tree",
        format!(
            "{{\"case\":\"fmm-laplace-preserved\",\"inputs\":{{\"provenance\":\"committed source ee8e613806eaa9ac3979a8b139af5e1528b699d5\"}},\"axis\":\"real_operation_tree\",\"rung\":0,\"observed\":\"{fingerprint:016x}\",\"oracle\":\"38b04e0d5a5d1fd2\",\"err\":0,\"tol\":0,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
            start.elapsed().as_secs_f64()
        ),
    );
}
