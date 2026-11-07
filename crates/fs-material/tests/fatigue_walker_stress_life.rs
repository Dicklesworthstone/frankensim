//! Analytical Walker stress-life checks; synthetic inputs are not measured life.
#![deny(unsafe_code)]

use std::sync::Arc;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fatigue::source::{SnLifeMeasure, SnSourceRequest, SnStressMeasure};
use fs_material::fatigue::walker::{
    WalkerParameters, resolve_walker_sn, walker_gamma_property_key,
};
use fs_material::fatigue::{
    FatigueError, FatigueScope, MeanStressCorrection, RainflowLimits, RainflowMiner, SnCurve,
    SnPoint,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn limits() -> RainflowLimits {
    RainflowLimits {
        maximum_samples: 100,
        maximum_residue: 32,
    }
}
fn parameters(gamma: f64) -> WalkerParameters {
    WalkerParameters {
        gamma,
        minimum_load_ratio: -1.0,
        maximum_load_ratio: 0.9,
    }
}
fn scope() -> FatigueScope {
    FatigueScope {
        source_reference: "synthetic Walker oracle".into(),
        material: "synthetic alloy; solid".into(),
        process: "fixed synthetic process".into(),
        surface_and_geometry: "smooth axial coupon".into(),
        environment: "fixed synthetic air".into(),
        load_basis: "signed normal stress amplitude and full cycles".into(),
        temperature_k: 300.0,
    }
}
fn curve(reference_ratio: f64, scale: f64, gamma: f64) -> Arc<SnCurve> {
    // One physical Walker family: N=1e6/(fully-reversed amplitude/scale)^2.
    let factor = ((1.0 - reference_ratio) / 2.0).powf(1.0 - gamma);
    Arc::new(
        SnCurve::new_at_stress_ratio(
            scope(),
            vec![
                SnPoint {
                    amplitude_pa: scale * factor,
                    cycles_to_failure: 1e6,
                },
                SnPoint {
                    amplitude_pa: 100.0 * scale * factor,
                    cycles_to_failure: 100.0,
                },
            ],
            None,
            reference_ratio,
        )
        .unwrap(),
    )
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 3e-12,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn analytic_maximum_swt_and_mean_independent_limits_preserve_reference_r() {
    for (gamma, at_reversed, at_positive) in [
        (0.0, 8.0, 2.0),
        (0.5, (24.0_f64).sqrt(), (6.0_f64).sqrt()),
        (1.0, 3.0, 3.0),
    ] {
        let model = parameters(gamma);
        close(
            model.equivalent_amplitude(3.0, 5.0, -1.0).unwrap(),
            at_reversed,
        );
        close(
            model.equivalent_amplitude(3.0, 5.0, 0.5).unwrap(),
            at_positive,
        );
        close(
            MeanStressCorrection::Walker(model)
                .corrected_amplitude(3.0, 5.0)
                .unwrap(),
            at_reversed,
        );
        for reference in [-1.0, 0.0, 0.1, 0.75] {
            let mean = 3.0 * (1.0 + reference) / (1.0 - reference);
            assert_eq!(
                model.equivalent_amplitude(3.0, mean, reference).unwrap(),
                3.0
            );
        }
    }
    let model = parameters(0.37);
    let independent = 3.0 * (8.0_f64 / 3.0).powf(0.63);
    close(
        model.equivalent_amplitude(3.0, 5.0, -1.0).unwrap(),
        independent,
    );
    // No additional life benefit is inferred for a compressive mean, R=-3.
    assert!(matches!(
        parameters(0.5).equivalent_amplitude(4.0, -2.0, -1.0),
        Err(FatigueError::StressRatioOutsideWalkerDomain { .. })
    ));
    // Near R=1, a small absolute ratio error can mean a large amplitude error.
    let near_unity = WalkerParameters {
        maximum_load_ratio: 1.0 - 0.5 * f64::EPSILON,
        ..parameters(0.5)
    };
    let reference = 1.0 - f64::EPSILON;
    let expected = (0.5 * (1.0 - reference) * (1e15_f64 + 1.0)).sqrt();
    close(
        near_unity
            .equivalent_amplitude(1.0, 1e15, reference)
            .unwrap(),
        expected,
    );
    let narrower = WalkerParameters {
        maximum_load_ratio: 1.0 - 4.0 * f64::EPSILON,
        ..parameters(0.5)
    };
    assert!(matches!(
        narrower.equivalent_amplitude(2.0, 1e16, -1.0),
        Err(FatigueError::StressRatioOutsideWalkerDomain { .. })
    ));
}

#[test]
fn reference_curve_conversion_scale_and_every_chunk_boundary_give_same_damage() {
    for scale in [1e-200, 1.0, 1e200] {
        for reference in [-1.0, 0.0, 0.25] {
            let calibration = curve(reference, scale, 0.5);
            let make = || {
                RainflowMiner::new(
                    calibration.clone(),
                    MeanStressCorrection::Walker(parameters(0.5)),
                    limits(),
                )
                .unwrap()
            };
            // R=0, a=4*scale, sigma_max=8*scale, equivalent R=-1 amplitude=sqrt(32)*scale.
            let history = [0.0, 8.0 * scale, 0.0, 8.0 * scale, 0.0];
            let mut whole = make();
            whole.push_block(&history, &mut || false).unwrap();
            let preview = whole.preview_finish().unwrap();
            close(preview.miner_damage, 64.0 / 1e6);
            close(
                preview.maximum_corrected_amplitude_pa / scale,
                (16.0 * (1.0 - reference)).sqrt(),
            );
            let summary = whole.finish().unwrap();
            assert_eq!(summary, preview);
            for split in 0..=history.len() {
                let mut streamed = make();
                streamed
                    .push_block(&history[..split], &mut || false)
                    .unwrap();
                let mut restarted = streamed.clone();
                restarted
                    .push_block(&history[split..], &mut || false)
                    .unwrap();
                assert_eq!(restarted.finish().unwrap(), summary);
                assert_eq!(restarted.finish().unwrap(), summary);
            }
        }
    }
}

#[test]
fn late_ratio_failure_and_partial_endpoint_refusal_leave_live_history_unchanged() {
    let model = WalkerParameters {
        minimum_load_ratio: -1.0,
        maximum_load_ratio: 0.5,
        ..parameters(0.5)
    };
    let make = || {
        RainflowMiner::new(
            curve(-1.0, 1.0, 0.5),
            MeanStressCorrection::Walker(model),
            limits(),
        )
        .unwrap()
    };
    let mut counter = make();
    counter
        .push_block(&[-4.0, 4.0, -4.0], &mut || false)
        .unwrap();
    let before = counter.clone();
    // Earlier R=-1 and R=-1/2 cycles are valid; the nested [6,8] cycle has R=3/4.
    assert!(matches!(
        counter.push_block(&[4.0, -4.0, 6.0, 8.0, 6.0, 8.0, 6.0], &mut || false),
        Err(FatigueError::StressRatioOutsideWalkerDomain { .. })
    ));
    assert_eq!(counter, before);
    let mut partial = make();
    partial.push_block(&[0.0, 8.0, 6.0], &mut || false).unwrap();
    let before = partial.clone();
    assert!(matches!(
        partial.preview_finish(),
        Err(FatigueError::StressRatioOutsideWalkerDomain { .. })
    ));
    assert!(matches!(
        partial.finish(),
        Err(FatigueError::StressRatioOutsideWalkerDomain { .. })
    ));
    assert_eq!(partial, before);
    // Finishing the original reversal, rather than closing its callback chunk,
    // makes the endpoint an R=0 cycle within the same calibration.
    partial.push(0.0).unwrap();
    close(partial.finish().unwrap().miner_damage, 32.0 / 1e6);
}

#[test]
fn unsupported_models_nontensile_cycles_and_corrected_extrapolation_refuse() {
    for model in [
        WalkerParameters {
            gamma: -0.01,
            ..parameters(0.5)
        },
        WalkerParameters {
            gamma: 1.01,
            ..parameters(0.5)
        },
        WalkerParameters {
            gamma: f64::NAN,
            ..parameters(0.5)
        },
        WalkerParameters {
            minimum_load_ratio: f64::NEG_INFINITY,
            ..parameters(0.5)
        },
        WalkerParameters {
            minimum_load_ratio: -1.1,
            ..parameters(0.5)
        },
        WalkerParameters {
            maximum_load_ratio: 1.0,
            ..parameters(0.5)
        },
        WalkerParameters {
            minimum_load_ratio: 0.9,
            ..parameters(0.5)
        },
    ] {
        assert!(
            RainflowMiner::new(
                curve(-1.0, 1.0, 0.5),
                MeanStressCorrection::Walker(model),
                limits()
            )
            .is_err()
        );
    }
    let model = WalkerParameters {
        minimum_load_ratio: 0.0,
        ..parameters(0.5)
    };
    assert!(
        RainflowMiner::new(
            curve(-1.0, 1.0, 0.5),
            MeanStressCorrection::Walker(model),
            limits()
        )
        .is_err()
    );
    assert!(matches!(
        parameters(0.5).equivalent_amplitude(1.0, -2.0, -1.0),
        Err(FatigueError::NoTensileStress { .. })
    ));
    assert!(matches!(
        parameters(0.5).equivalent_amplitude(1e308, 1e308, -1.0),
        Err(FatigueError::Arithmetic)
    ));
    assert_eq!(
        parameters(0.5)
            .equivalent_amplitude(0.0, -100.0, -1.0)
            .unwrap(),
        0.0
    );
    let mut counter = RainflowMiner::new(
        curve(-1.0, 1.0, 0.0),
        MeanStressCorrection::Walker(WalkerParameters {
            maximum_load_ratio: 0.95,
            ..parameters(0.0)
        }),
        limits(),
    )
    .unwrap();
    // Raw amplitude 20 is in the table, but its maximum-stress equivalent 400 is not.
    counter.push_block(&[360.0, 400.0], &mut || false).unwrap();
    let before = counter.clone();
    assert!(matches!(
        counter.finish(),
        Err(FatigueError::AmplitudeOutsideCurve { .. })
    ));
    assert_eq!(counter, before);
}

fn request(stress_measure: SnStressMeasure, life_measure: SnLifeMeasure) -> SnSourceRequest {
    SnSourceRequest {
        stress_measure,
        life_measure,
        temperature_axis: "T".into(),
        stress_ratio_axis: "R".into(),
        surface_and_geometry: "synthetic smooth coupon".into(),
        environment: "fixed synthetic humidity".into(),
        load_basis: "signed axial amplitude at reference R=0.1".into(),
        maximum_points: 3,
    }
}
fn source_rows(request: &SnSourceRequest) -> (PropertyClaim, PropertyClaim, QueryPoint) {
    let count = if request.life_measure == SnLifeMeasure::Cycles {
        1.0
    } else {
        2.0
    };
    let stress = if request.stress_measure == SnStressMeasure::Amplitude {
        1.0
    } else {
        2.0
    };
    let domain = ValidityDomain::unconstrained()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            293.15,
            313.15,
        )
        .with_quantity("R", QuantitySpec::dimensional(Dims::NONE), 0.1, 0.1)
        .with_quantity(
            request.life_measure.axis_name(),
            QuantitySpec::dimensional(Dims::NONE),
            100.0 * count,
            1e6 * count,
        )
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0);
    let sn = PropertyClaim {
        key: request.stress_measure.property_key(),
        value: PropertyValue::Curve {
            abscissa: request.life_measure.axis_name().into(),
            abscissa_dims: Dims::NONE,
            knots: [(100.0, 100.0), (1e4, 10.0), (1e6, 1.0)]
                .map(|(n, a)| (n * count, a * stress))
                .into(),
            dims: Pressure::DIMS,
        },
        validity: domain.clone(),
        interpolation: InterpolationPolicy::TabulatedOnly,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic common Walker population".into(),
            license: "test-only".into(),
            artifact: None,
        },
    };
    let gamma = PropertyClaim {
        key: walker_gamma_property_key(),
        value: PropertyValue::Scalar {
            value: 0.5,
            dims: Dims::NONE,
        },
        validity: domain.with_quantity("R", QuantitySpec::dimensional(Dims::NONE), -1.0, 0.8),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: sn.provenance.clone(),
    };
    let point = QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 303.15)
        .unwrap()
        .with_quantity("R", QuantitySpec::dimensional(Dims::NONE), 0.1)
        .unwrap()
        .with_quantity(
            request.life_measure.axis_name(),
            QuantitySpec::dimensional(Dims::NONE),
            1e4 * count,
        )
        .unwrap()
        .with_quantity("humidity", QuantitySpec::dimensional(Dims::NONE), 0.2)
        .unwrap();
    (sn, gamma, point)
}
fn card(rows: Vec<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for row in rows {
        claims.insert_claim(row).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic alloy".into(),
            phase: "solid".into(),
            process: "fixed synthetic machining".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn source_fit_preserves_knots_units_receipts_and_binds_the_exact_policy() {
    let mut expected = None;
    for (stress, life) in [
        (SnStressMeasure::Amplitude, SnLifeMeasure::Cycles),
        (SnStressMeasure::Range, SnLifeMeasure::Reversals),
    ] {
        let request = request(stress, life);
        let (sn, gamma, point) = source_rows(&request);
        let material = card(vec![sn.clone(), gamma.clone()]);
        let (curve, correction) =
            resolve_walker_sn(&material, &point, &request, SINGLE, SINGLE).unwrap();
        assert_eq!(curve.reference_stress_ratio(), 0.1);
        assert_eq!(curve.source().unwrap().claim(), &sn);
        assert_eq!(curve.source().unwrap().query(), &point);
        assert_eq!(curve.walker_source().unwrap().claim(), &gamma);
        assert_eq!(curve.walker_source().unwrap().selection(), SINGLE);
        material
            .claims()
            .verify_receipt(&curve.walker_source().unwrap().answer().receipt)
            .unwrap();
        for (_, answer) in curve.source().unwrap().observations() {
            material.claims().verify_receipt(&answer.receipt).unwrap();
        }
        let original = curve.walker_source().unwrap().parameters();
        for changed in [
            MeanStressCorrection::Walker(WalkerParameters {
                gamma: 0.6,
                ..original
            }),
            MeanStressCorrection::Walker(WalkerParameters {
                maximum_load_ratio: 0.9,
                ..original
            }),
            MeanStressCorrection::MatchCalibration,
        ] {
            assert!(RainflowMiner::new(Arc::new(curve.clone()), changed, limits()).is_err());
        }
        let mut counter = RainflowMiner::new(Arc::new(curve), correction, limits()).unwrap();
        let maximum = 20.0 / (0.9_f64).sqrt();
        counter
            .push_block(&[0.0, maximum, 0.0, maximum, 0.0], &mut || false)
            .unwrap();
        let result = counter.finish().unwrap();
        close(result.maximum_corrected_amplitude_pa, 10.0);
        close(result.miner_damage, 2.0 / 1e4);
        if let Some(prior) = expected {
            assert_eq!(result, prior);
        }
        expected = Some(result);
    }
}

#[test]
fn source_temperature_interpolation_is_retained_and_foreign_or_incomplete_fits_refuse() {
    let request = request(SnStressMeasure::Amplitude, SnLifeMeasure::Cycles);
    let (sn, gamma, point) = source_rows(&request);
    let mut thermal = gamma.clone();
    thermal.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(293.15, 0.4), (313.15, 0.6)],
        dims: Dims::NONE,
    };
    thermal.interpolation = InterpolationPolicy::LinearInside;
    let material = card(vec![sn.clone(), thermal.clone()]);
    let (resolved, _) = resolve_walker_sn(&material, &point, &request, SINGLE, SINGLE).unwrap();
    assert_eq!(resolved.walker_source().unwrap().claim(), &thermal);
    close(resolved.walker_source().unwrap().parameters().gamma, 0.5);
    assert!(
        resolve_walker_sn(
            &material,
            &point.clone().with("T", 320.0).unwrap(),
            &request,
            SINGLE,
            SINGLE
        )
        .is_err()
    );
    for change in 0..6 {
        let mut wrong = gamma.clone();
        match change {
            0 => wrong.provenance.source.push_str(" other population"),
            1 => {
                wrong.validity = wrong.validity.with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    310.0,
                )
            }
            2 => {
                wrong.validity = wrong.validity.with_quantity(
                    request.life_measure.axis_name(),
                    QuantitySpec::dimensional(Dims::NONE),
                    1000.0,
                    1e6,
                )
            }
            3 => wrong.key = PropertyKey::new("walker-gamma", Dims::NONE),
            4 => {
                wrong.value = PropertyValue::Scalar {
                    value: 1.2,
                    dims: Dims::NONE,
                }
            }
            5 => {
                wrong.value = PropertyValue::Curve {
                    abscissa: request.life_measure.axis_name().into(),
                    abscissa_dims: Dims::NONE,
                    knots: vec![(100.0, 0.4), (1e4, 0.5), (1e6, 0.6)],
                    dims: Dims::NONE,
                };
                wrong.interpolation = InterpolationPolicy::TabulatedOnly;
            }
            _ => unreachable!(),
        }
        let material = card(vec![sn.clone(), wrong]);
        assert!(
            resolve_walker_sn(&material, &point, &request, SINGLE, SINGLE).is_err(),
            "case {change}"
        );
    }
}
