//! Manufactured laws and independent thermal-velocity averages. No evaluated
//! nuclear-data accuracy, shielding or bound-atom scattering claim is made.

use super::*;
use crate::cross_sections::{CrossSectionLimits, NeutronReaction};
use crate::endf::{EV_J, EndfCollapseLimits, EndfMetadata};
use crate::irradiation::{NeutronGroups, NuclideId};

const BARN: f64 = 1.0e-28;
const ZERO: DopplerTails = DopplerTails {
    below: DopplerTail::Zero,
    above: DopplerTail::Zero,
};
const INVERSE: DopplerTails = DopplerTails {
    below: DopplerTail::InverseSpeed,
    above: DopplerTail::InverseSpeed,
};

fn source(
    energies_j: Vec<f64>,
    cross_sections_m2: Vec<f64>,
    law: EndfInterpolation,
) -> PointwiseCrossSection {
    let count = energies_j.len();
    PointwiseCrossSection {
        data_id: "manufactured-doppler-source".into(),
        metadata: EndfMetadata {
            material_number: 9237,
            nuclide: NuclideId::new(92, 235, 0).unwrap(),
            temperature_k: 300.0,
            atomic_weight_ratio: 233.0,
            library_id: 0,
            library_version: 1,
            library_release: 0,
            material_revision: 2,
            resonance_flag: 2,
            derivation_flag: 1,
        },
        reaction: NeutronReaction::Fission,
        energies_j: energies_j.into(),
        cross_sections_m2: cross_sections_m2.into(),
        interpolation_regions: vec![EndfInterpolationRegion {
            end_point: count,
            interpolation: law,
        }]
        .into(),
        mass_difference_q_j: 1.0e-12,
        reaction_q_j: 2.0e-12,
    }
}

fn prepare(
    source: &PointwiseCrossSection,
    temperature: f64,
    tails: DopplerTails,
) -> PreparedDopplerBroadening {
    PreparedDopplerBroadening::new(source, temperature, tails, DopplerLimits::default(), || {
        Ok(())
    })
    .unwrap()
}

fn tight() -> DopplerAccuracy {
    DopplerAccuracy {
        relative_tolerance: 1.0e-10,
        ..DopplerAccuracy::default()
    }
}

fn sample(operator: &PreparedDopplerBroadening, energy: f64) -> DopplerSample {
    operator
        .at_energy(energy, tight(), DopplerLimits::default(), || Ok(()))
        .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual {actual:.17e}, expected {expected:.17e}, relative {}",
        (actual - expected) / expected
    );
}

#[test]
fn zero_temperature_increment_preserves_original_laws_jumps_and_domain_exactly() {
    let curve = source(
        vec![0.0, EV_J, EV_J, 2.0 * EV_J],
        vec![5.0 * BARN, 4.0 * BARN, 7.0 * BARN, 7.0 * BARN],
        EndfInterpolation::Histogram,
    );
    let before = curve.clone();
    let no_work = DopplerLimits {
        max_segment_visits: 0,
        max_kernel_evaluations: 0,
        ..DopplerLimits::default()
    };
    for tails in [ZERO, INVERSE] {
        let operator = prepare(&curve, 300.0, tails);
        for energy in [0.0, 0.5 * EV_J, EV_J, 2.0 * EV_J] {
            let result = operator
                .at_energy(energy, tight(), no_work, || Ok(()))
                .unwrap();
            assert_eq!(
                result.cross_section_m2,
                curve.cross_section_at_energy(energy).unwrap()
            );
            assert_eq!(result.estimated_quadrature_error_m2, 0.0);
            assert_eq!(result.omitted_kernel_tail_bound_m2, 0.0);
        }
        assert!(
            operator
                .at_energy(3.0 * EV_J, tight(), no_work, || Ok(()))
                .is_err()
        );
        let tabulation = operator
            .tabulate(
                "identity-on-original-grid".into(),
                curve.energies_j().to_vec(),
                tight(),
                no_work,
                || Ok(()),
            )
            .unwrap();
        assert_eq!(tabulation.curve.metadata(), curve.metadata());
        assert_eq!(tabulation.curve.energies_j(), curve.energies_j());
        assert_eq!(
            tabulation.curve.cross_sections_m2(),
            curve.cross_sections_m2()
        );
        assert_eq!(
            tabulation.curve.interpolation_regions(),
            curve.interpolation_regions()
        );
        assert_eq!(tabulation.source_data_id.as_ref(), curve.data_id());
        assert!(
            operator
                .tabulate(
                    "changed-grid".into(),
                    vec![EV_J, 2.0 * EV_J],
                    tight(),
                    no_work,
                    || Ok(())
                )
                .is_err()
        );
    }
    assert_eq!(curve, before);
}

#[test]
fn constant_sigma_matches_mean_relative_speed_and_develops_the_low_energy_tail() {
    let theta = BOLTZMANN_J_K * 600.0 / 233.0;
    let curve = source(
        vec![0.0, 400.0 * theta],
        vec![2.0 * BARN; 2],
        EndfInterpolation::Linear,
    );
    let operator = prepare(&curve, 900.0, ZERO);
    for y in [1e-7, 0.01, 0.2, 1.0, 4.0, 10.0] {
        // E[|v-V|]/v for an isotropic Maxwell gas. The original constant law
        // extends beyond every significant thermal velocity for these probes.
        let factor = (1.0 + 0.5 / (y * y)) * det::erf(y) + det::exp(-y * y) * INV_SQRT_PI / y;
        let result = sample(&operator, theta * y * y);
        close(result.cross_section_m2, 2.0 * BARN * factor, 2.0e-11);
        assert!(
            result.estimated_quadrature_error_m2 + result.omitted_kernel_tail_bound_m2
                <= tight().relative_tolerance * result.cross_section_m2
        );
    }
    let low = sample(&operator, theta * 1e-10).cross_section_m2;
    let four_times_energy = sample(&operator, theta * 4e-10).cross_section_m2;
    close(low / four_times_energy, 2.0, 2e-10);
}

#[test]
fn inverse_speed_tails_preserve_rate_from_tiny_to_large_incident_speed() {
    let theta = BOLTZMANN_J_K * 600.0 / 233.0;
    let rate = 2.0 * BARN * theta.sqrt();
    // A broad linear-in-energy mesh is NOT an exact 1/v law. Confine its
    // interpolation to this tiny interval; the two exact 1/v exterior laws
    // cover the rest. The middle interpolation error is below test tolerance.
    let energies = vec![theta * (1.0 - 1e-8), theta * (1.0 + 1e-8)];
    let values = energies.iter().map(|e| rate / e.sqrt()).collect();
    let operator = prepare(
        &source(energies, values, EndfInterpolation::Linear),
        900.0,
        INVERSE,
    );
    for y in [1e-10, 1e-4, 0.3, 1.0, 8.0, 1e6] {
        let energy = theta * y * y;
        close(
            sample(&operator, energy).cross_section_m2,
            rate / energy.sqrt(),
            2e-11,
        );
    }
}

// Independent radial target-speed integration. The angular integral of
// r*(1+b*r²) is analytic in the physical relative speed r. This does not
// use production speed-difference kernels, source cursors or quadrature nodes.
fn target_velocity_average(y: f64, b: f64) -> f64 {
    let integrand = |u: f64| {
        if u == 0.0 {
            return 0.0;
        }
        let high = y + u;
        let low = (y - u).abs();
        let angular = ((high.powi(3) - low.powi(3)) / 3.0 + b * (high.powi(5) - low.powi(5)) / 5.0)
            / (2.0 * y * u);
        4.0 * INV_SQRT_PI * u * u * (-u * u).exp() * angular / y
    };
    let panels = 32_768;
    let step = 12.0 / f64::from(panels);
    let mut sum = integrand(0.0) + integrand(12.0);
    for index in 1..panels {
        sum += if index % 2 == 0 { 2.0 } else { 4.0 } * integrand(f64::from(index) * step);
    }
    sum * step / 3.0
}

#[test]
fn linear_energy_law_matches_independent_velocity_average_and_uses_delta_temperature_and_awr() {
    let theta = BOLTZMANN_J_K * 600.0 / 233.0;
    let slope = 0.4;
    let curve = source(
        vec![0.0, 400.0 * theta],
        vec![3.0 * BARN, 3.0 * BARN * (1.0 + slope * 400.0)],
        EndfInterpolation::Linear,
    );
    let operator = prepare(&curve, 900.0, ZERO);
    for y in [0.2, 1.0, 3.5] {
        close(
            sample(&operator, theta * y * y).cross_section_m2,
            3.0 * BARN * target_velocity_average(y, slope),
            2e-10,
        );
    }
    let mut warmer_source = curve.clone();
    warmer_source.metadata.temperature_k = 400.0;
    let same_increment = prepare(&warmer_source, 1000.0, ZERO);
    let mut heavier_source = curve.clone();
    heavier_source.metadata.atomic_weight_ratio *= 4.0;
    let same_thermal_width = prepare(&heavier_source, 2700.0, ZERO);
    assert_eq!(sample(&operator, theta), sample(&same_increment, theta));
    assert_eq!(sample(&operator, theta), sample(&same_thermal_width, theta));
    assert!(
        sample(&prepare(&curve, 1000.0, ZERO), theta).cross_section_m2
            > sample(&operator, theta).cross_section_m2
    );
}

#[test]
fn heated_resonance_wings_feed_original_spectrum_collapse_burnup_and_heat() {
    use crate::EnergyChannels;
    use crate::burnup::{
        BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
        BurnupProduct, BurnupRate, BurnupReaction,
    };
    use crate::endf::spectrum::{self, FluxSpectrum, FluxSpectrumData, FluxSpectrumLimits};
    use fs_qty::Frequency;
    use fs_thermochem::decay::DecayLawV1;

    let curve = source(
        [0.01, 0.96, 1.0, 1.04, 4.0].map(|e| e * EV_J).to_vec(),
        vec![0.0, 0.0, 100.0 * BARN, 0.0, 0.0],
        EndfInterpolation::Linear,
    );
    let operator = prepare(&curve, 1500.0, ZERO);
    assert!(sample(&operator, EV_J).cross_section_m2 < 0.6 * 100.0 * BARN);
    assert!(sample(&operator, 1.06 * EV_J).cross_section_m2 > BARN);
    assert_eq!(curve.cross_section_at_energy(1.06 * EV_J).unwrap(), 0.0);
    let energies: Vec<f64> = (0..=400)
        .map(|i| (0.8 + f64::from(i) * 0.001) * EV_J)
        .collect();
    let hot = operator
        .tabulate(
            "manufactured-heated-resonance-1meV-grid".into(),
            energies,
            tight(),
            DopplerLimits::default(),
            || Ok(()),
        )
        .unwrap();
    assert_eq!(hot.curve.metadata().temperature_k, 1500.0);
    assert_eq!(
        hot.curve.metadata().atomic_weight_ratio,
        curve.metadata().atomic_weight_ratio
    );
    assert_eq!(hot.curve.metadata().nuclide, curve.metadata().nuclide);
    assert_eq!(
        hot.curve.metadata().material_revision,
        curve.metadata().material_revision
    );
    assert_eq!(hot.curve.reaction_q_j(), curve.reaction_q_j());
    assert_eq!(hot.source_data_id.as_ref(), curve.data_id());
    assert_eq!(hot.source_temperature_k, 300.0);
    let grid = Arc::new(
        NeutronGroups::new("resonance-window".into(), vec![0.94 * EV_J, 1.06 * EV_J]).unwrap(),
    );
    let flux = FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "manufactured-sloping-flux".into(),
            energies_j: vec![0.8 * EV_J, 1.2 * EV_J],
            flux_density_m2_s_j: vec![0.8e27 / EV_J, 1.2e27 / EV_J],
            interpolation_regions: vec![EndfInterpolationRegion {
                end_point: 2,
                interpolation: EndfInterpolation::Linear,
            }],
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap();
    let (table, group_flux) = spectrum::collapse_temperature_table(
        "heated-evaluation-under-authored-spectrum".into(),
        &[curve.clone(), hot.curve.clone()],
        &flux,
        Arc::clone(&grid),
        EndfCollapseLimits::default(),
        CrossSectionLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let hot_rate = table
        .rate_constant(&group_flux, 1500.0, 1, || Ok(()))
        .unwrap();
    let cold_rate = table
        .rate_constant(&group_flux, 300.0, 1, || Ok(()))
        .unwrap();
    assert!(hot_rate > 0.0 && hot_rate < cold_rate);
    let direct = hot
        .curve
        .collapse_with_spectrum(
            &flux,
            Arc::clone(&grid),
            EndfCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    close(
        hot_rate,
        direct.cross_sections_m2[0] * direct.flux.values_m2_s()[0],
        3e-14,
    );
    let flat = hot
        .curve
        .collapse_flat_in_energy(&grid, EndfCollapseLimits::default(), || Ok(()))
        .unwrap();
    assert!(flat[0] > 0.0);

    let stable = DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap();
    let limits = BurnupLimits::default();
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "doppler-to-original-burnup".into(),
                groups: grid,
                temperature_min_k: 300.0,
                temperature_max_k: 1500.0,
                nuclides: vec![
                    BurnupNuclide {
                        id: curve.metadata().nuclide,
                        law: stable,
                    },
                    BurnupNuclide {
                        id: NuclideId::new(40, 95, 0).unwrap(),
                        law: stable,
                    },
                ],
                reactions: vec![BurnupReaction {
                    data_id: "explicit-fragments-and-event-heat".into(),
                    parent: 0,
                    rate: BurnupRate::TabulatedFission {
                        cross_sections: Arc::new(table),
                    },
                    products: vec![BurnupProduct {
                        nuclide: 1,
                        yield_per_event: 2.0,
                    }],
                    energy_per_event: EnergyChannels {
                        deposited_j: 1e-12,
                        radiation_j: 0.0,
                        escaped_j: 0.0,
                    },
                }],
            },
            limits,
        )
        .unwrap(),
    );
    let mut inventory = BurnupInventory::new(model, vec![100.0, 0.0]).unwrap();
    let report = inventory
        .advance(
            2.0,
            &group_flux,
            1500.0,
            BurnupBudget {
                limits,
                nuclei_tolerance: 1e-9,
                energy_tolerance_j: 1e-18,
            },
            || Ok(()),
        )
        .unwrap();
    let events = 100.0 * -(-2.0 * hot_rate).exp_m1();
    close(report.fissions, events, 2e-12);
    close(inventory.populations()[1], 2.0 * events, 2e-12);
    close(report.released.deposited_j, events * 1e-12, 2e-12);
}

#[test]
fn histogram_jump_is_integrated_with_both_sides_and_a_finite_kernel_tail_bound() {
    let theta = BOLTZMANN_J_K * 600.0 / 233.0;
    let discontinuous = source(
        vec![0.0, theta, theta, 400.0 * theta],
        vec![BARN, BARN, 3.0 * BARN, 3.0 * BARN],
        EndfInterpolation::Histogram,
    );
    let implicit_jump = source(
        vec![0.0, theta, 400.0 * theta],
        vec![BARN, 3.0 * BARN, 3.0 * BARN],
        EndfInterpolation::Histogram,
    );
    let operator = prepare(&discontinuous, 900.0, ZERO);
    let implicit = prepare(&implicit_jump, 900.0, ZERO);
    for y in [0.1, 1.0, 3.0] {
        assert_eq!(
            sample(&operator, theta * y * y),
            sample(&implicit, theta * y * y)
        );
    }
    let truncated = operator
        .at_energy(
            theta,
            DopplerAccuracy {
                relative_tolerance: 1e-4,
                kernel_half_widths: 4.0,
                ..DopplerAccuracy::default()
            },
            DopplerLimits::default(),
            || Ok(()),
        )
        .unwrap();
    let full = sample(&operator, theta);
    assert!(
        (full.cross_section_m2 - truncated.cross_section_m2).abs()
            <= truncated.omitted_kernel_tail_bound_m2 + truncated.estimated_quadrature_error_m2
    );
    assert!(matches!(
        operator.at_energy(
            theta,
            DopplerAccuracy {
                relative_tolerance: 1e-12,
                kernel_half_widths: 4.0,
                ..DopplerAccuracy::default()
            },
            DopplerLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
}

#[test]
fn heating_refuses_cooling_invalid_queries_work_exhaustion_and_unresolved_grid_requests() {
    let curve = source(vec![0.0, EV_J], vec![BARN; 2], EndfInterpolation::Linear);
    for temperature in [0.0, 299.0, f64::NAN, f64::INFINITY] {
        assert!(
            PreparedDopplerBroadening::new(
                &curve,
                temperature,
                ZERO,
                DopplerLimits::default(),
                || Ok(())
            )
            .is_err()
        );
    }
    assert!(
        PreparedDopplerBroadening::new(
            &curve,
            900.0,
            ZERO,
            DopplerLimits {
                max_source_points: 1,
                ..DopplerLimits::default()
            },
            || Ok(())
        )
        .is_err()
    );
    let operator = prepare(&curve, 900.0, ZERO);
    for energy in [0.0, -EV_J, f64::NAN, f64::INFINITY] {
        assert!(
            operator
                .at_energy(energy, tight(), DopplerLimits::default(), || Ok(()))
                .is_err()
        );
    }
    for limits in [
        DopplerLimits {
            max_segment_visits: 0,
            ..DopplerLimits::default()
        },
        DopplerLimits {
            max_kernel_evaluations: 0,
            ..DopplerLimits::default()
        },
    ] {
        assert!(matches!(
            operator.at_energy(0.1 * EV_J, tight(), limits, || Ok(())),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    for accuracy in [
        DopplerAccuracy {
            relative_tolerance: -1.0,
            ..tight()
        },
        DopplerAccuracy {
            absolute_tolerance_m2: f64::NAN,
            ..tight()
        },
        DopplerAccuracy {
            relative_tolerance: 0.0,
            absolute_tolerance_m2: 0.0,
            ..tight()
        },
        DopplerAccuracy {
            kernel_half_widths: 3.0,
            ..tight()
        },
        DopplerAccuracy {
            kernel_half_widths: 27.0,
            ..tight()
        },
    ] {
        assert!(
            operator
                .at_energy(0.1 * EV_J, accuracy, DopplerLimits::default(), || Ok(()))
                .is_err()
        );
    }
    for grid in [vec![0.0, EV_J], vec![EV_J, EV_J], vec![EV_J, 0.5 * EV_J]] {
        assert!(
            operator
                .tabulate(
                    "invalid-grid".into(),
                    grid,
                    tight(),
                    DopplerLimits::default(),
                    || Ok(())
                )
                .is_err()
        );
    }
    assert!(
        operator
            .tabulate(
                "oversized-grid".into(),
                vec![0.1 * EV_J, EV_J],
                tight(),
                DopplerLimits {
                    max_output_points: 1,
                    ..DopplerLimits::default()
                },
                || Ok(())
            )
            .is_err()
    );
}

#[test]
fn source_and_output_remain_immutable_on_cancelled_or_exhausted_tabulation() {
    let curve = source(vec![0.0, EV_J], vec![BARN; 2], EndfInterpolation::Linear);
    let operator = prepare(&curve, 900.0, ZERO);
    let grid = vec![0.01 * EV_J, 0.1 * EV_J, 0.5 * EV_J];
    let mut calls = 0;
    let baseline = operator
        .tabulate(
            "test-output".into(),
            grid.clone(),
            tight(),
            DopplerLimits::default(),
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut call = 0;
        let result = operator.tabulate(
            "test-output".into(),
            grid.clone(),
            tight(),
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
    // Determine the actual per-point allowance through the public API. A
    // guessed panel count can reject a valid single query before this test
    // reaches its intended invariant: all output knots share one work counter.
    let mut per_point_limit = 0;
    for &energy in &grid {
        let mut lower = 0;
        let mut upper = DopplerLimits::default().max_kernel_evaluations;
        while lower < upper {
            let middle = lower + (upper - lower) / 2;
            match operator.at_energy(
                energy,
                tight(),
                DopplerLimits {
                    max_kernel_evaluations: middle,
                    ..DopplerLimits::default()
                },
                || Ok(()),
            ) {
                Ok(_) => upper = middle,
                Err(EndfError::Physics(FissionError::Budget("Doppler kernel evaluations"))) => {
                    lower = middle + 1;
                }
                other => panic!("unexpected point-budget probe: {other:?}"),
            }
        }
        per_point_limit = per_point_limit.max(upper);
    }
    assert!(per_point_limit > 0);
    for &energy in &grid {
        let one = operator
            .at_energy(
                energy,
                tight(),
                DopplerLimits {
                    max_kernel_evaluations: per_point_limit,
                    ..DopplerLimits::default()
                },
                || Ok(()),
            )
            .unwrap();
        assert!(one.cross_section_m2 > 0.0);
    }
    assert!(matches!(
        operator.tabulate(
            "bounded-output".into(),
            grid.clone(),
            tight(),
            DopplerLimits {
                max_kernel_evaluations: per_point_limit,
                ..DopplerLimits::default()
            },
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
    assert_eq!(
        operator
            .tabulate(
                "test-output".into(),
                grid,
                tight(),
                DopplerLimits::default(),
                || Ok(())
            )
            .unwrap(),
        baseline
    );
}
