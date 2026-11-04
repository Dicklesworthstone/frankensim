//! G1/G3/G4: independently soluble field modes, actual material heat and bias.
//! Manufactured coefficients exercise the owner; they are not device validation.

#[path = "nonlinear_winding_circuit/common.rs"]
mod common;

use fs_em::{
    planar::{
        eddy::{
            EddyConductor,
            harmonic::{HarmonicWindingControls, HarmonicWindingError, PlanarHarmonicWindings},
        },
        nonlinear::PlanarNonlinearMagnetostatics,
    },
    winding::{
        WindingConfig,
        incremental::{IncrementalWindingConfig, PlanarIncrementalWindingBasis},
        ohmic::WindingConductor,
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_math::c64::C64;

fn controls() -> HarmonicWindingControls {
    HarmonicWindingControls {
        relative_tolerance: 1e-11,
        absolute_residual_a_m: 1e-12,
        restart: 8,
        max_cycles: 50,
        max_krylov_scalars: 10_000,
        max_field_scalars: 10_000,
        impedance_tolerance_ohm: 1e-8,
        power_tolerance_w: 1e-8,
        current_balance_tolerance_a: 1e-10,
        max_cell_skin_ratio: 1000.0,
    }
}

// Two disjoint stranded windings; each has its own outgoing and return cell.
// The second winding is reversed and has half the first winding's turn density.
// The remaining half of each square is a separate floating conducting region.
fn profiles() -> Vec<Vec<f64>> {
    vec![
        vec![24.0, 0.0, 0.0, 0.0, -24.0, 0.0, 0.0, 0.0],
        vec![0.0, -12.0, 0.0, 0.0, 0.0, 12.0, 0.0, 0.0],
    ]
}
fn incremental_config() -> IncrementalWindingConfig {
    IncrementalWindingConfig {
        winding: WindingConfig {
            max_windings: 4,
            max_basis_values: 10_000,
            pairing_tolerance_h: 1e-10,
            current_balance_tolerance_turns: 1e-10,
            inverse_residual_tolerance: 1e-10,
        },
        minimum_smooth_field_margin_t: 1e-8,
    }
}
fn basis<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    bias: f64,
) -> PlanarIncrementalWindingBasis<'m> {
    PlanarIncrementalWindingBasis::new(
        problem,
        profiles(),
        &common::boundary(),
        vec![bias, 0.0],
        vec![0.02, 0.03],
        &[0.0; 8],
        incremental_config(),
        common::material::controls().nonlinear,
        &mut || false,
    )
    .unwrap()
}
fn eddies(card: &MaterialCard) -> Vec<EddyConductor> {
    let conductor = common::material::conductor(card, 300.0, ElectricalPropertyBasis::Resistivity);
    [2, 3, 6, 7]
        .into_iter()
        .map(|triangle| EddyConductor {
            triangle,
            region: if triangle < 4 { 7 } else { 9 },
            conductor: conductor.clone(),
        })
        .collect()
}
fn copper(card: &MaterialCard) -> Vec<WindingConductor> {
    common::conductors(card, 300.0)
        .into_iter()
        .filter(|s| [0, 1, 4, 5].contains(&s.triangle))
        .collect()
}
fn prepared<'b, 'm>(
    basis: &'b PlanarIncrementalWindingBasis<'m>,
    card: &MaterialCard,
) -> PlanarHarmonicWindings<'b, 'm> {
    PlanarHarmonicWindings::new(
        basis,
        eddies(card),
        &[0, 1, 4, 5],
        copper(card),
        controls(),
        &mut || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn complex_near(actual: C64, expected: C64) {
    near(actual.re, expected.re);
    near(actual.im, expected.im);
}

#[test]
fn bias_frequency_conductivity_depth_and_signed_mutual_impedance_follow_field_modes() {
    for depth in [0.3, 0.8] {
        for sigma in [18.0, 36.0] {
            let source = common::source(true, sigma);
            let problem = common::problem(&source, depth);
            for bias in [0.0, 0.4, 2.0, -2.0] {
                let basis = basis(&problem, bias);
                let model = prepared(&basis, &source);
                let slope = if bias.abs() > 1.0 { 4.0 } else { 1.0 };
                near(model.resistance_ohm()[0], 288.0 * depth / sigma);
                near(model.resistance_ohm()[1], 72.0 * depth / sigma);
                let mut previous_l = f64::INFINITY;
                for omega in [1e-5, 0.5, 8.0, 100.0, 1000.0] {
                    let frequency = model.at_frequency(omega, &mut || false).unwrap();
                    // Independent one-mode elimination in each square:
                    // K=4 depth dH/dB; C= sigma depth/36. The latter follows
                    // from exact P1 mass MINUS the floating-region mean.
                    // Omitting that projection triples C and fails this oracle.
                    let denominator = C64::new(slope, -omega * sigma / 144.0);
                    let linkage = C64::from_re(2.0 * depth) / denominator;
                    let active_z = C64::new(omega * linkage.im, -omega * linkage.re);
                    let turns = [1.0, -0.5];
                    for i in 0..2 {
                        for j in 0..2 {
                            let mut expected = active_z.scale(turns[i] * turns[j]);
                            if i == j {
                                expected = expected
                                    + C64::new(
                                        model.resistance_ohm()[i],
                                        -omega * basis.end_inductance_h()[i],
                                    );
                            }
                            complex_near(frequency.impedance_ohm()[2 * i + j], expected);
                            near(
                                frequency.effective_inductance_h()[2 * i + j],
                                -expected.im / omega,
                            );
                        }
                    }
                    assert!(frequency.impedance_ohm()[0].re > model.resistance_ohm()[0]);
                    assert!(frequency.impedance_ohm()[0].im < 0.0);
                    assert!(frequency.impedance_ohm()[1].re < 0.0); // allowed signed mutual loss
                    assert!(frequency.effective_inductance_h()[0] < previous_l);
                    previous_l = frequency.effective_inductance_h()[0];
                    assert!(frequency.maximum_reciprocity_defect_ohm() < 1e-8);
                    assert!(frequency.maximum_power_pairing_defect_ohm() < 1e-8);
                    assert!(
                        frequency
                            .unit_responses()
                            .iter()
                            .all(|r| r.report().converged_euclidean())
                    );
                    if omega < 1e-4 {
                        near(
                            frequency.effective_inductance_h()[0],
                            basis.inductance_h()[0],
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn arbitrary_phases_close_spatial_heat_storage_and_terminal_power_without_double_counting() {
    let depth = 0.7;
    let sigma = 18.0;
    let source = common::source(true, sigma);
    let problem = common::problem(&source, depth);
    let basis = basis(&problem, 2.0);
    let model = prepared(&basis, &source);
    let omega = 40.0;
    let frequency = model.at_frequency(omega, &mut || false).unwrap();
    let currents = [C64::new(0.2, -0.1), C64::new(-0.4, 0.3)];
    let response = frequency.response(&currents, &mut || false).unwrap();
    let effective_current = currents[0] - currents[1].scale(0.5);
    let a = effective_current.scale(0.5) / C64::new(4.0, -omega * sigma / 144.0);
    complex_near(response.potential_wb_m[6], a);
    complex_near(response.potential_wb_m[7], -a);
    assert!(response.potential_wb_m[..6].iter().all(|&a| a == C64::ZERO));
    let expected_eddy = omega * omega * sigma * depth / 36.0 * a.norm_sq();
    let expected_copper = 0.5
        * (model.resistance_ohm()[0] * currents[0].norm_sq()
            + model.resistance_ohm()[1] * currents[1].norm_sq());
    let expected_storage = 2.0 * depth * 4.0 * a.norm_sq()
        + 0.25 * (0.02 * currents[0].norm_sq() + 0.03 * currents[1].norm_sq());
    near(response.mean_eddy_power_w, expected_eddy);
    near(response.mean_winding_power_w, expected_copper);
    near(
        response.mean_supplied_power_w,
        expected_eddy + expected_copper,
    );
    near(response.mean_magnetic_energy_j, expected_storage);
    near(response.reactive_power_var, -2.0 * omega * expected_storage);
    near(response.power_residual_w, 0.0);
    near(response.reactive_power_residual_var, 0.0);
    assert_eq!(response.eddy_cells.len(), 4);
    assert_eq!(response.winding_cells.len(), 4);
    for cell in &response.eddy_cells {
        let field = if cell.source.triangle < 4 { a } else { -a };
        for (k, fraction) in [-1.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0].iter().enumerate() {
            let e = C64::new(-omega * field.im, omega * field.re).scale(*fraction);
            complex_near(cell.electric_field_v_m[k], e);
            complex_near(cell.current_density_a_m2[k], e.scale(sigma));
        }
        near(
            cell.mean_power_w,
            sigma * (depth / 4.0) * omega * omega * a.norm_sq() / 36.0,
        );
        assert!(
            !response
                .winding_cells
                .iter()
                .any(|w| w.source.triangle == cell.source.triangle)
        );
    }
    for &(_, current) in &response.region_currents_a {
        complex_near(current, C64::ZERO);
    }
    for cell in &response.winding_cells {
        assert_eq!(
            cell.source.conductor,
            model
                .winding_conductors()
                .iter()
                .find(|s| s.triangle == cell.source.triangle)
                .unwrap()
                .conductor
        );
    }
    let zero = frequency.response(&[C64::ZERO; 2], &mut || false).unwrap();
    assert_eq!(zero.mean_supplied_power_w, 0.0);
    assert_eq!(zero.mean_magnetic_energy_j, 0.0);
}

#[test]
fn no_eddy_regions_reuse_dc_derivatives_and_preserve_explicit_end_storage() {
    let source = common::source(true, 18.0);
    let problem = common::problem(&source, 0.6);
    let basis = basis(&problem, 2.0);
    let model = PlanarHarmonicWindings::new(
        &basis,
        vec![],
        &(0..8).collect::<Vec<_>>(),
        copper(&source),
        controls(),
        &mut || false,
    )
    .unwrap();
    let frequency = model.at_frequency(57.0, &mut || false).unwrap();
    for i in 0..2 {
        for j in 0..2 {
            let z = frequency.impedance_ohm()[2 * i + j];
            near(
                z.re,
                if i == j {
                    model.resistance_ohm()[i]
                } else {
                    0.0
                },
            );
            near(z.im, -57.0 * basis.inductance_h()[2 * i + j]);
        }
        assert_eq!(
            frequency.unit_responses()[i].report().iters,
            basis.unit_responses()[i].report().iters
        );
    }
    let response = frequency
        .response(&[C64::ONE, C64::new(0.0, 0.2)], &mut || false)
        .unwrap();
    assert!(response.eddy_cells.is_empty());
    assert_eq!(response.mean_eddy_power_w, 0.0);
    assert_eq!(frequency.maximum_cell_skin_ratio(), 0.0);
    near(
        response.mean_supplied_power_w,
        response.mean_winding_power_w,
    );
}

#[test]
fn affine_dc_bias_keeps_transverse_reluctivity_gauge_invariance_and_skin_bound() {
    let source = common::source(true, 18.0);
    let depth = 0.7;
    let problem = common::problem(&source, depth);
    let omega = 40.0;
    let mut previous = None;
    for gauge in [0.0, 7.0] {
        let initial: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|point| gauge + 1.5 * point[0])
            .collect();
        let boundary: Vec<_> = (0..6).map(|v| (v, initial[v])).collect();
        let basis = PlanarIncrementalWindingBasis::new(
            &problem,
            profiles(),
            &boundary,
            vec![0.0; 2],
            vec![0.02, 0.03],
            &initial,
            incremental_config(),
            common::material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
        let model = prepared(&basis, &source);
        let frequency = model.at_frequency(omega, &mut || false).unwrap();
        // Uniform |B|=1.5 gives H/B=2 and dH/dB=4. The free hat samples
        // each equally, so K=2 depth (2+4), rather than 4 depth * 4.
        let linkage = C64::from_re(2.0 * depth) / C64::new(3.0, -omega * 18.0 / 144.0);
        complex_near(
            frequency.impedance_ohm()[0],
            C64::new(
                model.resistance_ohm()[0] + omega * linkage.im,
                -omega * (linkage.re + 0.02),
            ),
        );
        // The smaller transverse reluctivity controls the smallest skin depth.
        near(
            frequency.maximum_cell_skin_ratio(),
            (omega * 18.0 / 4.0).sqrt(),
        );
        assert!(frequency.unit_responses().iter().all(|r| {
            r.potential_per_current_h_m()[..6]
                .iter()
                .all(|&a| a == C64::ZERO)
        }));
        if let Some(old) = previous {
            for (a, b) in frequency.impedance_ohm().iter().zip(old) {
                complex_near(*a, b);
            }
        }
        previous = Some(frequency.impedance_ohm().to_vec());
    }
}

#[test]
fn refuses_overlap_bad_regions_skin_underresolution_failed_solve_and_cancellation() {
    let source = common::source(true, 18.0);
    let problem = common::problem(&source, 0.6);
    let basis = basis(&problem, 2.0);
    let original = basis.operating_field().potential_wb_m().to_vec();
    let mut overlapping = eddies(&source);
    overlapping[0].triangle = 0;
    overlapping[1].triangle = 1;
    assert!(
        PlanarHarmonicWindings::new(
            &basis,
            overlapping,
            &[2, 3, 4, 5],
            copper(&source),
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut split = eddies(&source);
    split[1].region = 22;
    assert!(
        PlanarHarmonicWindings::new(
            &basis,
            split,
            &[0, 1, 4, 5],
            copper(&source),
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut cfg = controls();
    cfg.max_field_scalars = 1;
    assert!(
        PlanarHarmonicWindings::new(
            &basis,
            eddies(&source),
            &[0, 1, 4, 5],
            copper(&source),
            cfg,
            &mut || false
        )
        .is_err()
    );
    cfg = controls();
    cfg.max_krylov_scalars = 1;
    assert!(
        PlanarHarmonicWindings::new(
            &basis,
            eddies(&source),
            &[0, 1, 4, 5],
            copper(&source),
            cfg,
            &mut || false
        )
        .is_err()
    );
    cfg = controls();
    cfg.max_cell_skin_ratio = 0.1;
    let underresolved = PlanarHarmonicWindings::new(
        &basis,
        eddies(&source),
        &[0, 1, 4, 5],
        copper(&source),
        cfg,
        &mut || false,
    )
    .unwrap();
    assert!(underresolved.at_frequency(100.0, &mut || false).is_err());
    cfg = controls();
    cfg.restart = 1;
    cfg.max_cycles = 1;
    let stalled = PlanarHarmonicWindings::new(
        &basis,
        eddies(&source),
        &[0, 1, 4, 5],
        copper(&source),
        cfg,
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        stalled.at_frequency(100.0, &mut || false),
        Err(HarmonicWindingError::Nonconvergence(_))
    ));
    let model = prepared(&basis, &source);
    for frequency in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(model.at_frequency(frequency, &mut || false).is_err());
    }
    let mut polls = 0;
    assert!(matches!(
        model.at_frequency(100.0, &mut || {
            polls += 1;
            polls > 30
        }),
        Err(HarmonicWindingError::Cancelled)
    ));
    let frequency = model.at_frequency(100.0, &mut || false).unwrap();
    assert!(frequency.response(&[C64::ONE], &mut || false).is_err());
    assert!(matches!(
        frequency.response(&[C64::ONE; 2], &mut || true),
        Err(HarmonicWindingError::Cancelled)
    ));
    assert_eq!(basis.operating_field().potential_wb_m(), original);
}

#[cfg(feature = "harmonic-exposure")]
mod exposure {
    use super::*;
    use fs_em::planar::eddy::harmonic::{
        SolidWindingRegion,
        exposure::{HarmonicWindingExposureLimits, HarmonicWindingHeatSource, HarmonicWindingTone},
    };
    use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};

    fn limits() -> HarmonicWindingExposureLimits {
        HarmonicWindingExposureLimits {
            maximum_lines: 8,
            maximum_pair_products: 100_000,
            maximum_workspace_bytes: 100_000,
            power_tolerance_w: 1e-8,
        }
    }
    fn signal(amplitude: C64, omega: f64, time: f64) -> f64 {
        amplitude.re * (omega * time).cos() + amplitude.im * (omega * time).sin()
    }
    fn velocity(amplitude: C64, omega: f64, time: f64) -> f64 {
        omega * (-amplitude.re * (omega * time).sin() + amplitude.im * (omega * time).cos())
    }
    fn integrate(window: HarmonicWindow, quantity: impl Fn(f64) -> f64) -> f64 {
        let n = 4000;
        let dt = window.duration_s() / f64::from(n);
        let mut sum = 0.0;
        for step in 0..=n {
            let weight = if step == 0 || step == n {
                1.0
            } else if step % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight * quantity(window.start_s() + f64::from(step) * dt);
        }
        sum / (3.0 * f64::from(n))
    }

    #[test]
    fn beat_and_carrier_heat_close_independent_real_field_energy_and_terminal_work() {
        let depth = 0.6;
        let sigma = 18.0;
        let source = common::source(true, sigma);
        let problem = common::problem(&source, depth);
        let basis = basis(&problem, 0.4);
        let model = prepared(&basis, &source);
        let omegas = [13.0, 17.0];
        let frequencies = omegas.map(|w| model.at_frequency(w, &mut || false).unwrap());
        let currents = [
            [C64::new(0.2, -0.07), C64::new(-0.03, 0.09)],
            [C64::new(-0.13, 0.08), C64::new(0.04, 0.06)],
        ];
        let tones = [0, 1].map(|i| HarmonicWindingTone {
            frequency: &frequencies[i],
            currents_a: &currents[i],
        });
        let bank = model
            .prepare_exposure(&tones, limits(), &mut || false)
            .unwrap();
        assert_eq!(bank.distinct_frequencies(), 2);
        // Independent weak one-mode elimination in each square, K=4*depth,
        // C=sigma*depth/36 and F*I=2*depth*(I0-I1/2). The second hat is -a.
        let a = [0, 1].map(|j| {
            (currents[j][0] - currents[j][1].scale(0.5)) / C64::new(2.0, -omegas[j] * sigma / 72.0)
        });
        let resistance = [288.0 * depth / sigma, 72.0 * depth / sigma];
        let end = [0.02, 0.03];
        let amplitude = |time| (0..2).map(|j| signal(a[j], omegas[j], time)).sum::<f64>();
        let derivative = |time| (0..2).map(|j| velocity(a[j], omegas[j], time)).sum::<f64>();
        let current = |winding, time| {
            (0..2)
                .map(|j| signal(currents[j][winding], omegas[j], time))
                .sum::<f64>()
        };
        let stored = |time| {
            4.0 * depth * amplitude(time).powi(2)
                + (0..2)
                    .map(|k| 0.5 * end[k] * current(k, time).powi(2))
                    .sum::<f64>()
        };
        let voltage = |k: usize, time| {
            let turns = [1.0, -0.5][k];
            resistance[k] * current(k, time)
                + 4.0 * depth * turns * derivative(time)
                + end[k]
                    * (0..2)
                        .map(|j| velocity(currents[j][k], omegas[j], time))
                        .sum::<f64>()
        };
        for window in [
            HarmonicWindow::new(0.017, 0.079).unwrap(),
            HarmonicWindow::new(-0.13, 0.51).unwrap(),
        ] {
            let got = bank.observe(window, &mut || false).unwrap();
            let copper = integrate(window, |t| {
                (0..2).map(|k| resistance[k] * current(k, t).powi(2)).sum()
            });
            let eddy = sigma * depth / 18.0 * integrate(window, |t| derivative(t).powi(2));
            near(got.mean_dissipated_power_w, copper + eddy);
            near(got.mean_magnetic_energy_j, integrate(window, stored));
            near(
                got.magnetic_energy_change_j,
                stored(window.start_s() + window.duration_s()) - stored(window.start_s()),
            );
            near(
                got.mean_supplied_power_w,
                integrate(window, |t| {
                    (0..2).map(|k| voltage(k, t) * current(k, t)).sum()
                }),
            );
            near(
                got.supplied_work_j,
                got.dissipated_heat_j + got.magnetic_energy_change_j,
            );
            assert!(got.power_residual_w.abs() < 1e-10);
            assert_eq!(got.cells.len(), 8);
            for cell in got.cells {
                let expected = match cell.source {
                    HarmonicWindingHeatSource::Floating(source) => {
                        assert!([2, 3, 6, 7].contains(&source.triangle));
                        eddy / 4.0
                    }
                    HarmonicWindingHeatSource::Stranded { source, winding } => {
                        assert!([0, 1, 4, 5].contains(&source.triangle));
                        resistance[winding] / 2.0
                            * integrate(window, |t| current(winding, t).powi(2))
                    }
                    HarmonicWindingHeatSource::Solid { .. } => {
                        panic!("fixture has no solid winding")
                    }
                };
                near(cell.mean_power_w, expected);
                near(cell.heat_j, expected * window.duration_s());
                near(cell.volume_m3, 0.25 * depth);
            }
        }
    }

    #[test]
    fn complete_common_period_recovers_the_original_cycle_means() {
        let source = common::source(true, 18.0);
        let problem = common::problem(&source, 0.6);
        let basis = basis(&problem, 2.0);
        let model = prepared(&basis, &source);
        let frequencies = [12.0, 36.0].map(|w| model.at_frequency(w, &mut || false).unwrap());
        let currents = [
            [C64::new(0.2, -0.3), C64::new(0.1, 0.2)],
            [C64::new(-0.1, 0.05), C64::new(0.08, -0.03)],
        ];
        let tones = [0, 1].map(|i| HarmonicWindingTone {
            frequency: &frequencies[i],
            currents_a: &currents[i],
        });
        let original = [0, 1].map(|i| {
            frequencies[i]
                .response(&currents[i], &mut || false)
                .unwrap()
        });
        let bank = model
            .prepare_exposure(&tones, limits(), &mut || false)
            .unwrap();
        let got = bank
            .observe(
                HarmonicWindow::new(-0.029, core::f64::consts::TAU / 12.0).unwrap(),
                &mut || false,
            )
            .unwrap();
        near(
            got.mean_dissipated_power_w,
            original.iter().map(|r| r.mean_dissipated_power_w()).sum(),
        );
        near(
            got.mean_magnetic_energy_j,
            original.iter().map(|r| r.mean_magnetic_energy_j).sum(),
        );
        near(
            got.mean_supplied_power_w,
            original.iter().map(|r| r.mean_supplied_power_w).sum(),
        );
        assert!(got.magnetic_energy_change_j.abs() < 1e-12);
    }

    #[test]
    fn equal_frequency_currents_combine_before_large_cancelled_fields_are_recovered() {
        let source = common::source(true, 18.0);
        let problem = common::problem(&source, 0.6);
        let basis = basis(&problem, 0.4);
        let model = prepared(&basis, &source);
        let frequency = model.at_frequency(30.0, &mut || false).unwrap();
        let large = 2.0_f64.powi(54);
        let currents = [
            [C64::from_re(large), C64::ZERO],
            [C64::new(0.125, 0.05), C64::new(0.03, -0.02)],
            [C64::from_re(-large), C64::ZERO],
        ];
        let tones = [0, 1, 2].map(|i| HarmonicWindingTone {
            frequency: &frequency,
            currents_a: &currents[i],
        });
        let bank = model
            .prepare_exposure(&tones, limits(), &mut || false)
            .unwrap();
        let single = model
            .prepare_exposure(&tones[1..2], limits(), &mut || false)
            .unwrap();
        assert_eq!(bank.distinct_frequencies(), 1);
        let window = HarmonicWindow::new(0.0, 1e-8).unwrap();
        let got = bank.observe(window, &mut || false).unwrap();
        let expected = single.observe(window, &mut || false).unwrap();
        assert_eq!(
            got.mean_dissipated_power_w,
            expected.mean_dissipated_power_w
        );
        assert_eq!(got.mean_magnetic_energy_j, expected.mean_magnetic_energy_j);
        assert_eq!(got.supplied_work_j, expected.supplied_work_j);
        // A subcycle current quadrature is not automatically its cycle average.
        let cycle = frequency.response(&currents[1], &mut || false).unwrap();
        assert!((got.mean_dissipated_power_w - cycle.mean_dissipated_power_w()).abs() > 1e-3);
    }

    #[test]
    fn solid_winding_total_current_heat_retains_cross_frequency_terms() {
        let sigma = 18.0;
        let source = common::source(false, sigma);
        let problem = common::problem(&source, 0.6);
        let basis = PlanarIncrementalWindingBasis::new(
            &problem,
            vec![vec![1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0]],
            &common::boundary(),
            vec![0.0],
            vec![0.03],
            &[0.0; 8],
            incremental_config(),
            common::material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
        let conductor =
            common::material::conductor(&source, 300.0, ElectricalPropertyBasis::Resistivity);
        let eddies = (0..8)
            .map(|triangle| EddyConductor {
                triangle,
                region: if triangle < 4 { 7 } else { 9 },
                conductor: conductor.clone(),
            })
            .collect();
        let model = PlanarHarmonicWindings::new_with_solid_regions(
            &basis,
            eddies,
            &[],
            vec![],
            vec![
                SolidWindingRegion {
                    region: 7,
                    winding: 0,
                    signed_turns: 1.0,
                },
                SolidWindingRegion {
                    region: 9,
                    winding: 0,
                    signed_turns: -1.0,
                },
            ],
            controls(),
            &mut || false,
        )
        .unwrap();
        let omegas = [11.0, 19.0];
        let frequencies = omegas.map(|w| model.at_frequency(w, &mut || false).unwrap());
        let currents = [[C64::new(0.2, -0.1)], [C64::new(-0.1, 0.07)]];
        let tones = [0, 1].map(|i| HarmonicWindingTone {
            frequency: &frequencies[i],
            currents_a: &currents[i],
        });
        let bank = model
            .prepare_exposure(&tones, limits(), &mut || false)
            .unwrap();
        let window = HarmonicWindow::new(0.014, 0.085).unwrap();
        let got = bank.observe(window, &mut || false).unwrap();
        // Independent uniform-current-constrained mode: A=I/(12-i*omega*sigma/6).
        // E at the two exterior vertices and centre is I/sigma + i*omega*(A-A/3).
        let fields = [0, 1].map(|j| {
            let a = currents[j][0] / C64::new(12.0, -omegas[j] * sigma / 6.0);
            let circulation = C64::new(-omegas[j] * a.im, omegas[j] * a.re);
            let drive = currents[j][0].scale(1.0 / sigma);
            [
                drive - circulation.scale(1.0 / 3.0),
                drive - circulation.scale(1.0 / 3.0),
                drive + circulation.scale(2.0 / 3.0),
            ]
        });
        for cell in &got.cells {
            assert!(matches!(
                cell.source,
                HarmonicWindingHeatSource::Solid { .. }
            ));
            // Independent three-point degree-two triangle quadrature at each
            // time, instead of the native mass-matrix product formula.
            let expected = sigma
                * cell.volume_m3
                * integrate(window, |time| {
                    let electric: [f64; 3] = std::array::from_fn(|k| {
                        (0..2).map(|j| signal(fields[j][k], omegas[j], time)).sum()
                    });
                    (0..3)
                        .map(|q| {
                            (0..3)
                                .map(|k| electric[k] * if k == q { 2.0 / 3.0 } else { 1.0 / 6.0 })
                                .sum::<f64>()
                                .powi(2)
                                / 3.0
                        })
                        .sum()
                });
            near(cell.mean_power_w, expected);
        }
        near(
            got.mean_dissipated_power_w,
            got.cells.iter().map(|c| c.mean_power_w).sum(),
        );
        near(
            got.supplied_work_j,
            got.dissipated_heat_j + got.magnetic_energy_change_j,
        );
    }

    #[test]
    fn source_budget_window_and_cancellation_refusals_preserve_original_fields() {
        let source = common::source(true, 18.0);
        let problem = common::problem(&source, 0.6);
        let basis = basis(&problem, 2.0);
        let model = prepared(&basis, &source);
        let other = prepared(&basis, &source);
        let frequency = model.at_frequency(30.0, &mut || false).unwrap();
        let foreign = other.at_frequency(30.0, &mut || false).unwrap();
        let current = [C64::new(0.1, 0.2); 2];
        let tone = HarmonicWindingTone {
            frequency: &frequency,
            currents_a: &current,
        };
        assert!(
            model
                .prepare_exposure(
                    &[HarmonicWindingTone {
                        frequency: &foreign,
                        currents_a: &current
                    }],
                    limits(),
                    &mut || false
                )
                .is_err()
        );
        let mut budget = limits();
        budget.maximum_pair_products = 0;
        assert!(
            model
                .prepare_exposure(&[tone], budget, &mut || false)
                .is_err()
        );
        budget = limits();
        budget.maximum_workspace_bytes = 1;
        assert!(
            model
                .prepare_exposure(&[tone], budget, &mut || false)
                .is_err()
        );
        let bank = model
            .prepare_exposure(&[tone], limits(), &mut || false)
            .unwrap();
        assert!(
            bank.observe(
                HarmonicWindow::with_profile(0.0, 0.1, HarmonicWindowProfile::Triangular).unwrap(),
                &mut || false
            )
            .is_err()
        );
        let window = HarmonicWindow::new(0.0, 0.1).unwrap();
        let original = bank.observe(window, &mut || false).unwrap();
        let mut polls = 0;
        assert!(matches!(
            bank.observe(window, &mut || {
                polls += 1;
                polls > 20
            }),
            Err(HarmonicWindingError::Cancelled)
        ));
        let replay = bank.observe(window, &mut || false).unwrap();
        assert_eq!(replay.supplied_work_j, original.supplied_work_j);
        let dark = model
            .prepare_exposure(&[], limits(), &mut || false)
            .unwrap();
        let got = dark.observe(window, &mut || false).unwrap();
        assert_eq!(got.dissipated_heat_j, 0.0); // nonzero DC bias is not extra AC heat
        assert_eq!(got.mean_magnetic_energy_j, 0.0);
        assert_eq!(got.supplied_work_j, 0.0);
        assert_eq!(got.cells.len(), 8);
    }
}
