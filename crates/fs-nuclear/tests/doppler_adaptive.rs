//! G1/G3 manufactured laws through the public ENDF, Doppler and collapse APIs.
//! No evaluated-data accuracy or interval-certificate claim is made.

use std::sync::Arc;

use fs_nuclear::FissionError;
use fs_nuclear::cross_sections::{CrossSectionLimits, NeutronReaction};
use fs_nuclear::endf::doppler::adaptive::{AdaptiveDopplerTabulation, DopplerGridAccuracy};
use fs_nuclear::endf::doppler::{
    DopplerAccuracy, DopplerLimits, DopplerTail, DopplerTails, PreparedDopplerBroadening,
};
use fs_nuclear::endf::spectrum::{self, FluxSpectrum, FluxSpectrumData, FluxSpectrumLimits};
use fs_nuclear::endf::{
    EndfCollapseLimits, EndfError, EndfInterpolation, EndfInterpolationRegion, EndfReadLimits,
    EndfSelection, PointwiseCrossSection,
};
use fs_nuclear::irradiation::{NeutronGroups, NuclideId};

const EV: f64 = 1.602_176_634e-19;
const BARN: f64 = 1e-28;
const ZERO: DopplerTails = DopplerTails {
    below: DopplerTail::Zero,
    above: DopplerTail::Zero,
};
const INVERSE: DopplerTails = DopplerTails {
    below: DopplerTail::InverseSpeed,
    above: DopplerTail::InverseSpeed,
};

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let text = format!("{value:.precision$e}");
        let (mantissa, exponent) = text.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let field = format!("{mantissa}{exponent:+}");
        if field.len() <= 11 {
            return format!("{field:>11}");
        }
    }
    panic!("manufactured ENDF field exceeds eleven columns");
}

fn record(fields: [String; 6], mf: i32, mt: i32) -> String {
    assert!(fields.iter().all(|s| s.len() == 11));
    format!("{}{:4}{mf:2}{mt:3}{:5}\n", fields.concat(), 9237, 1)
}

fn control(c1: f64, c2: f64, fields: [i32; 4], mf: i32, mt: i32) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", fields[0]),
            format!("{:11}", fields[1]),
            format!("{:11}", fields[2]),
            format!("{:11}", fields[3]),
        ],
        mf,
        mt,
    )
}

/// Minimal complete selected MF1/MF3 sections; energies in eV, values in barns.
fn source(points: &[(f64, f64)], law: EndfInterpolation) -> PointwiseCrossSection {
    let mut text = [
        control(92235.0, 233.0, [2, 1, 0, 7], 1, 451),
        control(0.0, 0.0, [0, 0, 0, 6], 1, 451),
        control(1.0, 2e7, [3, 0, 10, 8], 1, 451),
        control(300.0, 0.0, [1, 0, 0, 0], 1, 451),
        control(0.0, 0.0, [0; 4], 1, 0),
        control(92235.0, 233.0, [0; 4], 3, 18),
        control(2e8, -1e6, [0, 0, 1, points.len() as i32], 3, 18),
    ]
    .concat();
    let law = match law {
        EndfInterpolation::Histogram => 1,
        EndfInterpolation::Linear => 2,
    };
    text.push_str(&record(
        [
            format!("{:11}", points.len()),
            format!("{law:11}"),
            format!("{:11}", 0),
            format!("{:11}", 0),
            format!("{:11}", 0),
            format!("{:11}", 0),
        ],
        3,
        18,
    ));
    for chunk in points.chunks(3) {
        let mut fields = std::array::from_fn(|_| real(0.0));
        for (i, &(energy, value)) in chunk.iter().enumerate() {
            fields[2 * i] = real(energy);
            fields[2 * i + 1] = real(value);
        }
        text.push_str(&record(fields, 3, 18));
    }
    text.push_str(&control(0.0, 0.0, [0; 4], 3, 0));
    PointwiseCrossSection::parse(
        "manufactured-adaptive-source".into(),
        EndfSelection {
            material_number: 9237,
            nuclide: NuclideId::new(92, 235, 0).unwrap(),
            reaction: NeutronReaction::Fission,
            temperature_k: 300.0,
        },
        &text,
        EndfReadLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn prepare(
    curve: &PointwiseCrossSection,
    temperature: f64,
    tails: DopplerTails,
) -> PreparedDopplerBroadening {
    PreparedDopplerBroadening::new(curve, temperature, tails, DopplerLimits::default(), || {
        Ok(())
    })
    .unwrap()
}

fn point_accuracy() -> DopplerAccuracy {
    DopplerAccuracy {
        relative_tolerance: 1e-10,
        absolute_tolerance_m2: 1e-13 * BARN,
        ..DopplerAccuracy::default()
    }
}

fn grid_accuracy(relative: f64) -> DopplerGridAccuracy {
    DopplerGridAccuracy {
        relative_tolerance: relative,
        absolute_tolerance_m2: 1e-8 * BARN,
        ..DopplerGridAccuracy::default()
    }
}

fn adapt(
    operator: &PreparedDopplerBroadening,
    interval: [f64; 2],
    relative: f64,
) -> AdaptiveDopplerTabulation {
    operator
        .tabulate_adaptive(
            "manufactured-adaptive-result".into(),
            interval,
            point_accuracy(),
            grid_accuracy(relative),
            DopplerLimits::default(),
            || Ok(()),
        )
        .unwrap()
}

fn inverse_source() -> PointwiseCrossSection {
    // These endpoint rates agree within 4e-13. The very narrow linear middle
    // interval adds still smaller error; both exterior laws are exact 1/v.
    source(
        &[(1.0, 2.0), (1.000001, 1.999999)],
        EndfInterpolation::Linear,
    )
}

#[test]
fn identity_adaptation_retains_the_original_domain_jumps_and_laws_without_work() {
    let curve = source(
        &[(0.0, 5.0), (1.0, 4.0), (1.0, 7.0), (2.0, 7.0)],
        EndfInterpolation::Histogram,
    );
    let operator = prepare(&curve, 300.0, INVERSE);
    let no_work = DopplerLimits {
        max_segment_visits: 0,
        max_kernel_evaluations: 0,
        ..DopplerLimits::default()
    };
    let result = operator
        .tabulate_adaptive(
            "identity".into(),
            [0.0, 2.0 * EV],
            point_accuracy(),
            grid_accuracy(1e-3),
            no_work,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(result.tabulation.curve.energies_j(), curve.energies_j());
    assert_eq!(
        result.tabulation.curve.cross_sections_m2(),
        curve.cross_sections_m2()
    );
    assert_eq!(
        result.tabulation.curve.interpolation_regions(),
        curve.interpolation_regions()
    );
    assert_eq!(result.tabulation.curve.metadata(), curve.metadata());
    assert_eq!(result.point_evaluations, 0);
    assert_eq!(result.maximum_estimated_interpolation_error_m2, 0.0);
    assert_eq!(result.tabulation.source_data_id.as_ref(), curve.data_id());
    assert_eq!(result.tabulation.curve.data_id(), "identity");
    assert!(
        operator
            .tabulate_adaptive(
                "crop".into(),
                [EV, 2.0 * EV],
                point_accuracy(),
                grid_accuracy(1e-3),
                no_work,
                || Ok(())
            )
            .is_err()
    );
}

#[test]
fn source_feature_guard_finds_an_off_center_resonance_invisible_to_initial_probes() {
    let center = 1.17321;
    let curve = source(
        &[
            (0.01, 0.0),
            (center - 0.0002, 0.0),
            (center, 100.0),
            (center + 0.0002, 0.0),
            (4.0, 0.0),
        ],
        EndfInterpolation::Linear,
    );
    let operator = prepare(&curve, 306.0, ZERO);
    for energy in [0.8, 0.95, 1.1, 1.25, 1.4] {
        let sample = operator
            .at_energy(
                energy * EV,
                point_accuracy(),
                DopplerLimits::default(),
                || Ok(()),
            )
            .unwrap();
        assert!(sample.cross_section_m2 < 1e-12 * BARN);
    }
    let result = adapt(&operator, [0.8 * EV, 1.4 * EV], 2e-3);
    let hot = &result.tabulation.curve;
    assert!(hot.cross_section_at_energy(center * EV).unwrap() > BARN);
    assert!(hot.energies_j().len() > curve.energies_j().len());
    assert_eq!(hot.energies_j()[0], 0.8 * EV);
    assert_eq!(*hot.energies_j().last().unwrap(), 1.4 * EV);
    assert_eq!(hot.metadata().temperature_k, 306.0);
    assert_eq!(hot.metadata().nuclide, curve.metadata().nuclide);
    assert_eq!(hot.reaction_q_j(), curve.reaction_q_j());
    assert_eq!(result.tabulation.source_data_id.as_ref(), curve.data_id());
    // Nondyadic probes also inspect the wings between retained output knots.
    for i in 0..137 {
        let energy = (center - 0.016 + (f64::from(i) + 0.37) * 0.032 / 137.0) * EV;
        let direct = operator
            .at_energy(
                energy,
                point_accuracy(),
                DopplerLimits::default(),
                || Ok(()),
            )
            .unwrap()
            .cross_section_m2;
        let sampled = hot.cross_section_at_energy(energy).unwrap();
        assert!(
            (sampled - direct).abs() <= 4e-3 * direct + 3e-8 * BARN,
            "{energy} {sampled} {direct}"
        );
    }
}

#[test]
fn geometric_refinement_resolves_the_inverse_speed_tail_across_energy_decades() {
    let curve = inverse_source();
    let operator = prepare(&curve, 900.0, INVERSE);
    let result = adapt(&operator, [1e-6 * EV, 100.0 * EV], 2e-3);
    for i in 0..117 {
        let energy_ev = 10.0_f64.powf(-6.0 + (f64::from(i) + 0.37) * 8.0 / 117.0);
        let exact = 2.0 * BARN / energy_ev.sqrt();
        let actual = result
            .tabulation
            .curve
            .cross_section_at_energy(energy_ev * EV)
            .unwrap();
        assert!((actual - exact).abs() <= 2.2e-3 * exact);
    }
    assert!(result.maximum_estimated_interpolation_error_m2.is_finite());
    assert!(result.tabulation.curve.energies_j().len() < 10_000);
    assert_eq!(operator.source(), &curve);
    assert_eq!(result, adapt(&operator, [1e-6 * EV, 100.0 * EV], 2e-3));
}

#[test]
fn refined_native_curves_converge_in_flat_and_spectrum_weighted_reaction_rates() {
    let curve = inverse_source();
    let operator = prepare(&curve, 900.0, INVERSE);
    let interval = [0.01 * EV, 16.0 * EV];
    let groups = Arc::new(
        NeutronGroups::new(
            "adaptive-reaction-groups".into(),
            [0.01, 0.25, 1.0, 4.0, 16.0].map(|e| e * EV).to_vec(),
        )
        .unwrap(),
    );
    let flux_scale = 1e22 / EV;
    let flux = FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "manufactured-affine-spectrum".into(),
            energies_j: interval.to_vec(),
            flux_density_m2_s_j: interval.map(|e| flux_scale * (1.0 + e / EV)).to_vec(),
            interpolation_regions: vec![EndfInterpolationRegion {
                end_point: 2,
                interpolation: EndfInterpolation::Linear,
            }],
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap();
    let [low, high] = interval;
    let rate = 2.0 * BARN * EV.sqrt();
    let exact = flux_scale
        * rate
        * (2.0 * (high.sqrt() - low.sqrt())
            + 2.0 / (3.0 * EV) * (high * high.sqrt() - low * low.sqrt()));
    let mut errors = Vec::new();
    for tolerance in [2e-2, 1e-4] {
        let hot = adapt(&operator, interval, tolerance).tabulation.curve;
        let flat = hot
            .collapse_flat_in_energy(&groups, EndfCollapseLimits::default(), || Ok(()))
            .unwrap();
        for (i, &value) in flat.iter().enumerate() {
            let edges = groups.edges_j();
            let expected = 2.0 * rate / (edges[i + 1].sqrt() + edges[i].sqrt());
            assert!((value - expected).abs() <= 2.0 * tolerance * expected);
        }
        let (table, group_flux) = spectrum::collapse_temperature_table(
            "adaptive-native-table".into(),
            &[hot],
            &flux,
            Arc::clone(&groups),
            EndfCollapseLimits::default(),
            CrossSectionLimits::default(),
            || Ok(()),
        )
        .unwrap();
        let actual = table
            .rate_constant(&group_flux, 900.0, groups.len(), || Ok(()))
            .unwrap();
        errors.push((actual - exact).abs() / exact);
    }
    assert!(
        errors[1] < errors[0],
        "coarse {} fine {}",
        errors[0],
        errors[1]
    );
    assert!(errors[1] < 2e-4);
}

#[test]
fn adaptive_work_and_output_limits_are_shared_and_cancellation_is_atomic() {
    let curve = source(&[(0.0, 2.0), (100.0, 2.0)], EndfInterpolation::Linear);
    let operator = prepare(&curve, 900.0, ZERO);
    let theta = 1.380_649e-23 * 600.0 / 233.0;
    let interval = [theta, 25.0 * theta];
    let bounded = DopplerLimits {
        max_kernel_evaluations: 500,
        ..DopplerLimits::default()
    };
    assert!(
        operator
            .at_energy(theta, point_accuracy(), bounded, || Ok(()))
            .is_ok()
    );
    assert!(matches!(
        operator.tabulate_adaptive(
            "shared-work".into(),
            interval,
            point_accuracy(),
            grid_accuracy(1e-3),
            bounded,
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
    for limits in [
        DopplerLimits {
            max_output_points: 2,
            ..DopplerLimits::default()
        },
        DopplerLimits {
            max_source_points: 1,
            ..DopplerLimits::default()
        },
        DopplerLimits {
            max_segment_visits: 0,
            ..DopplerLimits::default()
        },
    ] {
        assert!(matches!(
            operator.tabulate_adaptive(
                "bounded-output".into(),
                interval,
                point_accuracy(),
                grid_accuracy(1e-3),
                limits,
                || Ok(())
            ),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    let mut calls = 0;
    let baseline = operator
        .tabulate_adaptive(
            "cancelled-output".into(),
            interval,
            point_accuracy(),
            grid_accuracy(1e-3),
            DopplerLimits::default(),
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut call = 0;
        let result = operator.tabulate_adaptive(
            "cancelled-output".into(),
            interval,
            point_accuracy(),
            grid_accuracy(1e-3),
            DopplerLimits::default(),
            || {
                call += 1;
                if call == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(result, Err(EndfError::Physics(FissionError::Cancelled)));
        assert_eq!(operator.source(), &curve);
    }
    assert_eq!(
        baseline,
        operator
            .tabulate_adaptive(
                "cancelled-output".into(),
                interval,
                point_accuracy(),
                grid_accuracy(1e-3),
                DopplerLimits::default(),
                || Ok(())
            )
            .unwrap()
    );
}

#[test]
fn invalid_accuracy_and_domains_are_refused_and_zero_curves_stay_compact() {
    let curve = source(
        &[(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (10.0, 0.0)],
        EndfInterpolation::Linear,
    );
    let operator = prepare(&curve, 900.0, ZERO);
    for domain in [
        [0.0, EV],
        [EV, EV],
        [2.0 * EV, EV],
        [f64::NAN, EV],
        [EV, f64::INFINITY],
    ] {
        assert!(
            operator
                .tabulate_adaptive(
                    "invalid".into(),
                    domain,
                    point_accuracy(),
                    grid_accuracy(1e-3),
                    DopplerLimits::default(),
                    || Ok(())
                )
                .is_err()
        );
    }
    for accuracy in [
        DopplerGridAccuracy {
            relative_tolerance: -1.0,
            ..grid_accuracy(1e-3)
        },
        DopplerGridAccuracy {
            relative_tolerance: 0.0,
            absolute_tolerance_m2: 0.0,
            ..grid_accuracy(1e-3)
        },
        DopplerGridAccuracy {
            absolute_tolerance_m2: f64::NAN,
            ..grid_accuracy(1e-3)
        },
        DopplerGridAccuracy {
            maximum_thermal_span: 0.0,
            ..grid_accuracy(1e-3)
        },
        DopplerGridAccuracy {
            maximum_thermal_span: 2.0,
            ..grid_accuracy(1e-3)
        },
    ] {
        assert!(
            operator
                .tabulate_adaptive(
                    "invalid".into(),
                    [EV, 2.0 * EV],
                    point_accuracy(),
                    accuracy,
                    DopplerLimits::default(),
                    || Ok(())
                )
                .is_err()
        );
    }
    let zero = adapt(&operator, [1e-8 * EV, 100.0 * EV], 1e-3);
    assert_eq!(zero.tabulation.curve.energies_j().len(), 2);
    assert_eq!(zero.tabulation.curve.cross_sections_m2(), &[0.0, 0.0]);
    assert_eq!(zero.maximum_estimated_interpolation_error_m2, 0.0);
}
