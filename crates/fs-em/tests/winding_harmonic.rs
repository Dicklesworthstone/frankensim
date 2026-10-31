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
