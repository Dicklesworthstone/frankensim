//! G0/G1/G3/G4: nonlinear power, native spatial transport and admission.

use super::*;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGridBudget, CellPaletteBudget, MassAbsorptionSpectrum};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::sync::Arc;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 19,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn medium<'a>(cx: &Cx<'_>, density: &'a [f64], kappa: f64) -> AbsorbingGrid<'a> {
    AbsorbingGrid::try_new(
        cx,
        VolumeGrid::new(
            [1, 1, density.len()],
            density,
            [0.0; 3],
            [1.0, 1.0, 1.0 / density.len() as f64],
        ),
        MassAbsorptionSpectrum::constant(kappa).unwrap(),
        None,
        AbsorbingGridBudget {
            max_grid_cells: 64,
            max_steps_per_segment: 128,
        },
    )
    .unwrap()
}

fn law(beta: f64) -> TwoPhotonAbsorption {
    TwoPhotonAbsorption::try_new(beta, (380.0, 780.0), f64::MAX, ContentHash([63; 32])).unwrap()
}

fn beam(
    grid: &TwoPhotonAbsorbingGrid<'_, '_>,
    cx: &Cx<'_>,
    power: f64,
    area: f64,
) -> Result<BeamPowerDeposition<1>, TwoPhotonError> {
    grid.deposit_beam_power(
        cx,
        Point3::new(0.5, 0.5, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        1.0,
        550.0,
        power,
        area,
        64,
    )
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-12 * b.abs().max(1e-300),
        "{a:e} != {b:e}"
    );
}

fn balance(report: &BeamPowerDeposition<1>) {
    close(
        report.absorbed_power_w[0],
        report.cells.iter().map(|c| c.absorbed_power_w[0]).sum(),
    );
    close(
        report.transmitted_power_w[0] + report.absorbed_power_w[0],
        report.incident_power_w[0],
    );
    assert_eq!(report.scattered_power_w, [0.0]);
    assert!(
        report
            .cells
            .iter()
            .all(|c| c.absorbed_power_w[0] > 0.0 && c.scattered_power_w == [0.0])
    );
}

#[test]
fn g1_two_photon_absorption_matches_an_independent_runge_kutta_solution() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.5; 4], 0.4);
        let nonlinear = TwoPhotonAbsorbingGrid::try_new(cx, &medium, law(0.35)).unwrap();
        let result = beam(&nonlinear, cx, 5.0, 2.0).unwrap();
        // Integrate dP/ds directly; no effective-length/transmission expression.
        let rhs = |p: f64| -0.6 * p - 0.175 * p * p;
        let mut reference = 5.0;
        let step = 1.0 / 4096.0;
        for _ in 0..4096 {
            let k1 = rhs(reference);
            let k2 = rhs(reference + step * k1 / 2.0);
            let k3 = rhs(reference + step * k2 / 2.0);
            let k4 = rhs(reference + step * k3);
            reference += step * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        }
        close(result.transmitted_power_w[0], reference);
        balance(&result);
    });
}

#[test]
fn g1_two_photon_linear_and_quadratic_limits_preserve_weak_absorption() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.0, 2.0, 0.0, 3.0], 0.7);
        let nonlinear = TwoPhotonAbsorbingGrid::try_new(cx, &medium, law(0.0)).unwrap();
        let actual = beam(&nonlinear, cx, 7.0, 0.2).unwrap();
        let linear = medium
            .deposit_beam_power(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                1.0,
                &[550.0],
                &[7.0],
                64,
            )
            .unwrap();
        assert_eq!(actual.cells, linear.cells);
        assert_eq!(actual.transmitted_power_w, linear.transmitted_power_w);
        assert_eq!(actual.absorbed_power_w, linear.absorbed_power_w);
        balance(&actual);
    });
    let (absorbed, transmitted) = absorb_two_photon(3.0, 2.0, 1.0, 2.0, 0.0, 2.0);
    close(transmitted, 3.0 / 7.0);
    close(absorbed, 18.0 / 7.0);
    let (absorbed, transmitted) = absorb_two_photon(2.0, 1.0, 1.0, 1.0, 1e-20, 1e-20);
    assert_eq!(transmitted, 2.0);
    close(absorbed, 6e-20);
}

#[test]
fn g1_two_photon_preserves_representable_watts_across_extreme_optical_fractions() {
    // q = 1e-350 underflows, but P*q = 1e-50 W is representable.
    let (absorbed, transmitted) = absorb_two_photon(1e300, 1e200, 1.0, 1e-250, 0.0, 1e-200);
    close(absorbed, 1e-50);
    assert_eq!(transmitted, 1e300);
    // q = 1e400 overflows, but P/(1+q) = 1e-100 W remains finite.
    let (absorbed, transmitted) = absorb_two_photon(1e300, 1e200, 1.0, 1e100, 0.0, 1e200);
    assert_eq!(absorbed, 1e300);
    close(transmitted, 1e-100);
    // A linearly opaque cell deposits the whole beam even when alpha*L overflows.
    assert_eq!(
        absorb_two_photon(2.0, 1.0, 1e300, 1e200, 1e200, 1.0),
        (2.0, 0.0)
    );
    // exp(-tau) underflows before applying a large incident power.
    let (_, transmitted) = absorb_two_photon(1e300, 1e300, 1.0, 1.0, 800.0, 0.2);
    close(transmitted, (1e300_f64.ln() - 800.0).exp() / 1.00025);
}

#[test]
fn g3_two_photon_subdivision_area_scaling_and_non_superposition() {
    with_cx(|cx, _| {
        let coarse = medium(cx, &[1.0], 0.4);
        let refined = medium(cx, &[1.0; 16], 0.4);
        let a = TwoPhotonAbsorbingGrid::try_new(cx, &coarse, law(0.3)).unwrap();
        let b = TwoPhotonAbsorbingGrid::try_new(cx, &refined, law(0.3)).unwrap();
        let full = beam(&a, cx, 8.0, 2.0).unwrap();
        let fine = beam(&b, cx, 8.0, 2.0).unwrap();
        close(full.transmitted_power_w[0], fine.transmitted_power_w[0]);
        close(full.absorbed_power_w[0], fine.absorbed_power_w[0]);
        let same_irradiance = beam(&a, cx, 4.0, 1.0).unwrap();
        close(
            full.absorbed_power_w[0],
            2.0 * same_irradiance.absorbed_power_w[0],
        );
        let split = beam(&a, cx, 4.0, 2.0).unwrap();
        assert!(full.absorbed_power_w[0] > 2.0 * split.absorbed_power_w[0]);
        balance(&full);
        balance(&fine);
        assert_eq!(beam(&b, cx, 8.0, 2.0).unwrap(), fine);
    });
}

#[test]
fn g1_two_photon_indexed_sources_preserve_cell_order_vacuum_and_common_support() {
    with_cx(|cx, _| {
        let medium = medium(cx, &[1.0, 0.0, 2.0], 0.3);
        let laws = CellPalette::try_new(
            cx,
            vec![
                TwoPhotonAbsorption::try_new(0.2, (500.0, 700.0), 100.0, ContentHash([1; 32]))
                    .unwrap(),
                TwoPhotonAbsorption::try_new(0.7, (400.0, 600.0), 100.0, ContentHash([2; 32]))
                    .unwrap(),
            ],
            Arc::from([0, 1, 1]),
            CellPaletteBudget {
                max_models: 2,
                max_cells: 3,
            },
        )
        .unwrap();
        let nonlinear = TwoPhotonAbsorbingGrid::try_new_indexed(cx, &medium, laws).unwrap();
        assert_eq!(nonlinear.wavelength_bounds_nm(), (500.0, 600.0));
        let report = beam(&nonlinear, cx, 4.0, 1.0).unwrap();
        assert_eq!(report.medium_identity, nonlinear.content_identity());
        assert_ne!(report.medium_identity, medium.content_identity());
        assert_eq!(
            report
                .cells
                .iter()
                .map(|c| c.cell_index)
                .collect::<Vec<_>>(),
            [0, 2]
        );
        let first = 4.0 * (-0.1_f64).exp() / (1.0 + 0.2 * 4.0 * -(-0.1_f64).exp_m1() / 0.3);
        let last = first * (-0.2_f64).exp() / (1.0 + 0.7 * first * -(-0.2_f64).exp_m1() / 0.6);
        close(report.cells[0].absorbed_power_w[0], 4.0 - first);
        close(report.cells[1].absorbed_power_w[0], first - last);
        balance(&report);
    });
}

#[test]
fn g0_g4_two_photon_refuses_invalid_sources_beams_domains_budgets_and_cancellation() {
    for beta in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            TwoPhotonAbsorption::try_new(beta, (380.0, 780.0), 1.0, ContentHash([1; 32])),
            Err(TwoPhotonError::Source)
        );
    }
    assert!(TwoPhotonAbsorption::try_new(1.0, (380.0, 780.0), 0.0, ContentHash([1; 32])).is_err());
    assert!(TwoPhotonAbsorption::try_new(1.0, (780.0, 380.0), 1.0, ContentHash([1; 32])).is_err());
    assert!(TwoPhotonAbsorption::try_new(1.0, (380.0, 780.0), 1.0, ContentHash([0; 32])).is_err());
    with_cx(|cx, gate| {
        let medium = medium(cx, &[1.0; 2], 0.4);
        let law =
            TwoPhotonAbsorption::try_new(0.2, (500.0, 600.0), 3.0, ContentHash([63; 32])).unwrap();
        let nonlinear = TwoPhotonAbsorbingGrid::try_new(cx, &medium, law).unwrap();
        for area in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert_eq!(
                beam(&nonlinear, cx, 2.0, area),
                Err(TwoPhotonError::BeamArea)
            );
        }
        assert_eq!(
            beam(&nonlinear, cx, 4.0, 1.0),
            Err(TwoPhotonError::Irradiance { cell_index: 0 })
        );
        assert_eq!(
            beam(&nonlinear, cx, -1.0, 1.0),
            Err(VolumeTransportError::BeamPower { lane: 0 }.into())
        );
        let query = |wavelength, limit| {
            nonlinear.deposit_beam_power(
                cx,
                Point3::new(0.5, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                1.0,
                wavelength,
                2.0,
                1.0,
                limit,
            )
        };
        assert_eq!(
            query(650.0, 2),
            Err(VolumeTransportError::Wavelength.into())
        );
        assert_eq!(query(550.0, 1), Err(VolumeTransportError::Budget.into()));
        gate.request();
        assert_eq!(
            beam(&nonlinear, cx, 2.0, 1.0),
            Err(VolumeTransportError::Cancelled.into())
        );
    });
}
