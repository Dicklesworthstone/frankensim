//! Independent analytic resistance fixtures for the all-frequency gate.

use fs_math::c64::C64;
use fs_vfit::passivity::{
    CertificateClass, GlobalPassivityBudget, GlobalPassivityUnresolved, GlobalPassivityVerdict,
    PassivityError, check_global_passivity, check_passivity, repair_passivity_global,
};
use fs_vfit::{PoleTerm, RationalModel};

fn real(pole: f64, residue: f64) -> PoleTerm {
    PoleTerm::Real { pole, residue }
}

fn certify(model: &RationalModel) -> fs_vfit::passivity::GlobalPassivityReport {
    check_global_passivity(model, GlobalPassivityBudget::default()).unwrap()
}

#[test]
fn zero_direct_term_and_lossless_improper_term_need_no_artificial_resistance() {
    let models = [
        RationalModel {
            terms: vec![real(-1.0, 1.0), real(-3.0, 0.25)],
            d: 0.0,
            e: 0.0,
        },
        RationalModel {
            terms: Vec::new(),
            d: 0.0,
            e: 2.0,
        },
        RationalModel {
            terms: Vec::new(),
            d: 0.0,
            e: 0.0,
        },
        RationalModel {
            terms: vec![PoleTerm::Pair {
                pole: C64::new(-0.25, 2.0),
                residue: C64::new(0.5, 0.0),
            }],
            d: 0.0,
            e: 0.125,
        },
    ];
    for model in models {
        let proof = certify(&model);
        assert_eq!(proof.verdict, GlobalPassivityVerdict::Certified);
        let (unchanged, repair) =
            repair_passivity_global(&model, (0.0, 100.0), GlobalPassivityBudget::default())
                .unwrap();
        assert_eq!(unchanged, model);
        assert_eq!(repair.rounds, 0);
        assert_eq!(repair.relative_perturbation, 0.0);
        assert_eq!(
            repair.certificate.verdict,
            GlobalPassivityVerdict::Certified
        );
    }
}

#[test]
fn exact_dc_zero_is_preserved_without_roundoff_clipping() {
    // Z(s)=s/(s+1): Re Z(iw)=w²/(1+w²), including exact zero at DC.
    let model = RationalModel {
        terms: vec![real(-1.0, -1.0)],
        d: 1.0,
        e: 0.0,
    };
    let proof = certify(&model);
    assert_eq!(proof.verdict, GlobalPassivityVerdict::Certified);
    assert_eq!(proof.subdivisions, 0);
    // One ulp BELOW this resistance changes the actual stored model at DC.
    // It must not be rounded back onto the passive boundary.
    let active = RationalModel {
        d: fs_math::next_down(1.0),
        ..model
    };
    assert_ne!(certify(&active).verdict, GlobalPassivityVerdict::Certified);
}

#[test]
fn dc_violation_outside_the_legacy_grid_is_repaired_with_poles_and_asymptotes_fixed() {
    // Re Z(0)=-2+1=-1. Above 0.1 rad/s the very slow negative branch
    // is invisible to the old requested-band grid, but can feed a slow loop.
    for gain in [2.0_f64.powi(-40), 1.0, 2.0_f64.powi(40)] {
        let model = RationalModel {
            terms: vec![real(-1.0e-8, -2.0e-8 * gain), real(-1.0, gain)],
            d: 0.0,
            e: 0.125 * gain,
        };
        let legacy = check_passivity(&model, (1.0, 100.0)).unwrap();
        assert!(legacy.passive);
        assert_eq!(legacy.class, CertificateClass::GridOnly);
        let proof = certify(&model);
        assert_eq!(proof.verdict, GlobalPassivityVerdict::Violated);
        assert_eq!(proof.violation_frequency_rad_s, Some(0.0));
        let (repaired, report) =
            repair_passivity_global(&model, (1.0, 100.0), GlobalPassivityBudget::default())
                .unwrap();
        assert_eq!(repaired.poles_expanded(), model.poles_expanded());
        assert_eq!((repaired.d, repaired.e), (model.d, model.e));
        assert!(report.rounds > 0);
        assert!(report.relative_perturbation > 0.0 && report.relative_perturbation < 1e-6);
        assert_eq!(
            report.certificate.verdict,
            GlobalPassivityVerdict::Certified
        );
        assert!(repaired.eval_iw(0.0).re > 0.0);
    }
}

#[test]
fn negative_infinite_tail_is_found_even_when_the_entire_requested_band_is_positive() {
    // Re Z(iw) = (4-2w²)/((1+w²)(4+w²)): positive below sqrt(2),
    // negative above. d=0 must not hide the negative approach to zero.
    let model = RationalModel {
        terms: vec![real(-1.0, 2.0), real(-2.0, -2.0)],
        d: 0.0,
        e: 0.0,
    };
    assert!(check_passivity(&model, (1e-4, 1e-3)).unwrap().passive);
    let proof = certify(&model);
    assert_eq!(proof.verdict, GlobalPassivityVerdict::Violated);
    let w = proof
        .violation_frequency_rad_s
        .expect("finite negative-tail witness");
    assert!(w > 2.0_f64.sqrt());
    assert!((4.0 - 2.0 * w * w) / ((1.0 + w * w) * (4.0 + w * w)) < 0.0);
}

#[test]
fn subdivision_proves_compensated_residues_and_preserves_frequency_and_impedance_scaling() {
    // Re Z(iw)=(w^4-4w²+7)/((w²+1)(w²+4)).
    // The numerator is (w²-2)²+3 > 0, although its middle coefficient
    // and one residue are negative. Termwise signs alone cannot prove it.
    for frequency in [2.0_f64.powi(-20), 1.0, 2.0_f64.powi(20)] {
        for gain in [2.0_f64.powi(-200), 1.0, 2.0_f64.powi(200)] {
            let model = RationalModel {
                terms: vec![
                    real(-frequency, 4.0 * gain * frequency),
                    real(-2.0 * frequency, -6.5 * gain * frequency),
                ],
                d: gain,
                e: 0.125 * gain / frequency,
            };
            let proof = certify(&model);
            assert_eq!(
                proof.verdict,
                GlobalPassivityVerdict::Certified,
                "{frequency} {gain}: {proof:?}"
            );
            assert!(proof.subdivisions > 0);
        }
    }
}

#[test]
fn a_narrow_active_resonance_cannot_hide_between_grid_nodes() {
    let model = RationalModel {
        terms: vec![
            real(-100.0, 100.0),
            PoleTerm::Pair {
                pole: C64::new(-1e-5, 2.345678),
                residue: C64::new(-1e-4, 0.0),
            },
        ],
        d: 0.0,
        e: 0.0,
    };
    assert!(check_passivity(&model, (1.0, 10.0)).unwrap().passive);
    assert!(model.eval_iw(2.345678).re < -8.0);
    let proof = certify(&model);
    assert_eq!(proof.verdict, GlobalPassivityVerdict::Violated, "{proof:?}");
    let w = proof
        .violation_frequency_rad_s
        .expect("finite resonance witness");
    assert!(model.eval_iw(w).re < 0.0);
}

#[test]
fn touching_zero_and_exhausted_budgets_remain_unresolved_and_cannot_trigger_a_repair() {
    // Exact Re Z=(w²-2)²/((w²+1)(w²+4)); the isolated zero is not a
    // negative band and cannot be replaced by an invented loss tolerance.
    let model = RationalModel {
        terms: vec![real(-1.0, 3.0), real(-2.0, -6.0)],
        d: 1.0,
        e: 0.0,
    };
    let budget = GlobalPassivityBudget {
        max_depth: 24,
        max_subdivisions: 256,
        ..Default::default()
    };
    let proof = check_global_passivity(&model, budget).unwrap();
    assert_eq!(proof.verdict, GlobalPassivityVerdict::Unresolved);
    assert!(proof.unresolved.is_some());
    assert!(matches!(
        repair_passivity_global(&model, (0.1, 10.0), budget),
        Err(PassivityError::GlobalUnresolved { .. })
    ));
    let limited = GlobalPassivityBudget {
        max_order: 1,
        ..Default::default()
    };
    let proof = check_global_passivity(&model, limited).unwrap();
    assert_eq!(
        proof.unresolved,
        Some(GlobalPassivityUnresolved::OrderLimit)
    );
    let no_subdivisions = GlobalPassivityBudget {
        max_subdivisions: 0,
        ..Default::default()
    };
    let proof = check_global_passivity(&model, no_subdivisions).unwrap();
    assert_eq!(
        proof.unresolved,
        Some(GlobalPassivityUnresolved::SubdivisionLimit)
    );
}

#[test]
fn malformed_and_nonrepresentable_models_never_receive_a_certificate() {
    let inertance = RationalModel {
        terms: Vec::new(),
        d: 0.0,
        e: -1.0,
    };
    assert_eq!(
        check_global_passivity(&inertance, Default::default()),
        Err(PassivityError::NegativeAsymptote)
    );
    let unstable = RationalModel {
        terms: vec![real(0.0, 1.0)],
        d: 0.0,
        e: 0.0,
    };
    assert_eq!(
        check_global_passivity(&unstable, Default::default()),
        Err(PassivityError::Unstable)
    );
    for model in [
        RationalModel {
            terms: vec![real(-1.0, f64::INFINITY)],
            d: 1.0,
            e: 0.0,
        },
        RationalModel {
            terms: Vec::new(),
            d: f64::NAN,
            e: 0.0,
        },
        RationalModel {
            terms: vec![PoleTerm::Pair {
                pole: C64::new(-1.0, 0.0),
                residue: C64::new(1.0, 0.0),
            }],
            d: 0.0,
            e: 0.0,
        },
    ] {
        assert!(matches!(
            check_global_passivity(&model, Default::default()),
            Err(PassivityError::InvalidModel { .. })
        ));
    }
    let overflow = RationalModel {
        terms: vec![real(-1e-300, 1e300)],
        d: 0.0,
        e: 0.0,
    };
    assert_eq!(
        certify(&overflow).unresolved,
        Some(GlobalPassivityUnresolved::NumericRange)
    );
}
