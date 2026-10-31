//! Analytic layered optical composition on the original borrowed cells.

use super::*;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use std::sync::Arc;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 41,
                    kernel_id: 7,
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

fn budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 4,
        max_steps_per_segment: 16,
    }
}

fn palette_budget() -> CellPaletteBudget {
    CellPaletteBudget {
        max_models: 4,
        max_cells: 4,
    }
}

fn palette(
    cx: &Cx<'_>,
    models: Vec<MassAbsorptionSpectrum>,
    cells: &[u32],
) -> CellPalette<MassAbsorptionSpectrum> {
    CellPalette::try_new(cx, models, Arc::from(cells), palette_budget()).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-100),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_indexed_layers_keep_spectral_optical_depth_and_ordered_lte_source() {
    with_cx(|cx, _| {
        let densities = [2.0, 3.0];
        let temperatures = [1200.0, 1800.0];
        let first = MassAbsorptionSpectrum::table(&[(380.0, 0.2), (780.0, 0.6)]).unwrap();
        let second = MassAbsorptionSpectrum::table(&[(380.0, 1.7), (780.0, 0.9)]).unwrap();
        let optics = palette(
            cx,
            vec![
                MassAbsorptionSpectrum::constant(100.0).unwrap(),
                first,
                second,
            ],
            &[1, 2],
        );
        let medium = AbsorbingGrid::try_new_indexed(
            cx,
            VolumeGrid::new([1, 1, 2], &densities, [0.0; 3], [1.0, 1.0, 0.4]),
            optics,
            Some(&temperatures),
            budget(),
        )
        .unwrap();
        assert!(core::ptr::eq(medium.grid.data.as_ptr(), densities.as_ptr()));
        assert!(core::ptr::eq(
            medium.temperature_k.unwrap().as_ptr(),
            temperatures.as_ptr()
        ));
        assert_eq!(medium.absorption_at_cell(0).unwrap(), &first);
        assert_eq!(medium.absorption_at_cell(1).unwrap(), &second);
        assert_eq!(
            medium.absorption_at_cell(2),
            Err(VolumeTransportError::Palette)
        );
        let wavelengths = [380.0, 550.0, 780.0];
        for (origin, direction, order) in [(-1.0, 1.0, [0, 1]), (2.0, -1.0, [1, 0])] {
            let result = medium
                .transport(
                    cx,
                    Point3::new(0.5, 0.5, origin),
                    Vec3::new(0.0, 0.0, direction),
                    f64::INFINITY,
                    &wavelengths,
                )
                .unwrap();
            for (lane, &wavelength) in wavelengths.iter().enumerate() {
                // Independent linear material laws and exact two-layer solve.
                let kappa = [
                    0.2 + 0.001 * (wavelength - 380.0),
                    1.7 - 0.002 * (wavelength - 380.0),
                ];
                let tau = [0.8 * kappa[0], 1.2 * kappa[1]];
                let (near, far) = (order[0], order[1]);
                let source = [
                    planck_radiance_w_m2_sr_nm(1200.0, wavelength).unwrap(),
                    planck_radiance_w_m2_sr_nm(1800.0, wavelength).unwrap(),
                ];
                close(result.transmittance[lane], (-(tau[0] + tau[1])).exp());
                close(
                    result.emission_w_m2_sr_nm[lane],
                    (1.0 - (-tau[near]).exp()) * source[near]
                        + (-tau[near]).exp() * (1.0 - (-tau[far]).exp()) * source[far],
                );
            }
        }
        let partial = medium
            .transmittance(
                cx,
                Point3::new(0.5, 0.5, 0.2),
                Vec3::new(0.0, 0.0, 1.0),
                0.4,
                &[550.0],
            )
            .unwrap();
        close(partial[0], (-0.2_f64 * (2.0 * 0.37 + 3.0 * 1.36)).exp());
    });
}

#[test]
fn g3_uniform_palette_preserves_native_transport_and_map_changes_identity() {
    with_cx(|cx, _| {
        let densities = [0.4, 1.5, 0.8];
        let temperatures = [1400.0, 1400.0, 2000.0];
        let grid = || VolumeGrid::new([1, 1, 3], &densities, [0.0; 3], [1.0, 1.0, 0.5]);
        let law = MassAbsorptionSpectrum::table(&[(380.0, 0.1), (780.0, 0.7)]).unwrap();
        let global =
            AbsorbingGrid::try_new(cx, grid(), law, Some(&temperatures), budget()).unwrap();
        let indexed = AbsorbingGrid::try_new_indexed(
            cx,
            grid(),
            palette(cx, vec![law], &[0, 0, 0]),
            Some(&temperatures),
            budget(),
        )
        .unwrap();
        for (z, d) in [(-1.0, 1.0), (2.0, -1.0)] {
            let evaluate = |medium: &AbsorbingGrid<'_>| {
                medium
                    .transport(
                        cx,
                        Point3::new(0.3, 0.4, z),
                        Vec3::new(0.0, 0.0, d),
                        f64::INFINITY,
                        &[380.0, 510.0, 780.0],
                    )
                    .unwrap()
            };
            assert_eq!(evaluate(&global), evaluate(&indexed));
        }
        let models = vec![law, MassAbsorptionSpectrum::constant(0.2).unwrap()];
        let a = AbsorbingGrid::try_new_indexed(
            cx,
            grid(),
            palette(cx, models.clone(), &[0, 1, 0]),
            None,
            budget(),
        )
        .unwrap();
        let b = AbsorbingGrid::try_new_indexed(
            cx,
            grid(),
            palette(cx, models, &[1, 0, 0]),
            None,
            budget(),
        )
        .unwrap();
        assert_ne!(a.content_identity(), b.content_identity());
        assert_ne!(
            a.coefficients_at_cell(0, &[550.0]).unwrap(),
            b.coefficients_at_cell(0, &[550.0]).unwrap()
        );
    });
}

#[test]
fn g0_vacuum_and_emission_support_follow_mapped_positive_density_cells() {
    with_cx(|cx, _| {
        let zero = MassAbsorptionSpectrum::constant(0.0).unwrap();
        let active = MassAbsorptionSpectrum::constant(1.0).unwrap();
        let temperatures = [1800.0, 1800.0];
        for (density, map, vacuum) in [
            ([1.0, 1.0], [0, 0], true),
            ([1.0, 0.0], [0, 1], true),
            ([1.0, 1.0], [0, 1], false),
        ] {
            let medium = AbsorbingGrid::try_new_indexed(
                cx,
                VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]),
                palette(cx, vec![zero, active], &map),
                Some(&temperatures),
                budget(),
            )
            .unwrap();
            assert_eq!(medium.is_vacuum(), vacuum);
            assert_eq!(medium.has_thermal_emission(), !vacuum);
            assert_eq!(medium.has_emission_proposal(), !vacuum);
            if !vacuum {
                // The conservative proposal may select nonabsorbing cells,
                // but the actual emitting layer must retain positive support.
                let origin = Point3::new(0.5, 0.5, -1.0);
                let direction = Vec3::new(0.0, 0.0, 1.0);
                assert!(
                    medium
                        .emission_direction_pdf(cx, origin, direction)
                        .unwrap()
                        > 0.0
                );
                let transport = medium
                    .transport(cx, origin, direction, f64::INFINITY, &[550.0])
                    .unwrap();
                let source = planck_radiance_w_m2_sr_nm(1800.0, 550.0).unwrap();
                close(
                    transport.emission_w_m2_sr_nm[0],
                    (1.0 - (-1.0_f64).exp()) * source,
                );
            }
        }
    });
}

#[test]
fn g0_palette_sharing_domains_lengths_budgets_and_cancellation_are_explicit() {
    with_cx(|cx, gate| {
        let original = CellPalette::try_new(
            cx,
            vec![10_u8, 20, 30],
            Arc::from([0_u32, 1]),
            palette_budget(),
        )
        .unwrap();
        let shared = original
            .with_values(cx, vec![0.2_f64, 0.4], palette_budget())
            .unwrap();
        assert!(Arc::ptr_eq(original.cell_indices(), shared.cell_indices()));
        assert_eq!(shared.at_cell(1), Some(&0.4));
        assert_eq!(shared.at_cell(2), None);
        assert_eq!(
            original.with_values(cx, vec![0_u8], palette_budget()),
            Err(VolumeTransportError::Palette)
        );
        for (values, indices) in [(vec![], vec![0_u32]), (vec![1], vec![]), (vec![1], vec![1])] {
            assert_eq!(
                CellPalette::try_new(cx, values, Arc::from(indices), palette_budget()),
                Err(VolumeTransportError::Palette)
            );
        }
        assert_eq!(
            CellPalette::try_new(cx, vec![0_u8; 5], Arc::from([0_u32]), palette_budget()),
            Err(VolumeTransportError::Budget)
        );
        let density = [1.0, 1.0];
        let grid = || VolumeGrid::new([1, 1, 2], &density, [0.0; 3], [1.0; 3]);
        let visible = MassAbsorptionSpectrum::constant(1.0).unwrap();
        let other = MassAbsorptionSpectrum::source_table(
            &[(400.0, 1.0), (800.0, 1.0)],
            ContentHash([1; 32]),
        )
        .unwrap();
        // Even an unreferenced model cannot expand or change the common band.
        assert!(matches!(
            AbsorbingGrid::try_new_indexed(
                cx,
                grid(),
                palette(cx, vec![visible, other], &[0, 0]),
                None,
                budget()
            ),
            Err(VolumeTransportError::AbsorptionSpectrum)
        ));
        assert!(matches!(
            AbsorbingGrid::try_new_indexed(
                cx,
                grid(),
                palette(cx, vec![visible], &[0]),
                None,
                budget()
            ),
            Err(VolumeTransportError::Palette)
        ));
        let broad = MassAbsorptionSpectrum::source_table(
            &[(300.0, 0.2), (1200.0, 0.6)],
            ContentHash([2; 32]),
        )
        .unwrap();
        let broad_medium = AbsorbingGrid::try_new_indexed(
            cx,
            grid(),
            palette(cx, vec![broad], &[0, 0]),
            None,
            budget(),
        )
        .unwrap();
        assert_eq!(broad_medium.wavelength_bounds_nm(), (300.0, 1200.0));
        assert_eq!(
            broad_medium.coefficients_at_cell(0, &[299.0]),
            Err(VolumeTransportError::Wavelength)
        );
        close(
            broad_medium
                .transmittance(
                    cx,
                    Point3::new(0.5, 0.5, -1.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    f64::INFINITY,
                    &[1200.0],
                )
                .unwrap()[0],
            (-1.2_f64).exp(),
        );
        let candidate = palette(cx, vec![visible], &[0, 0]);
        gate.request();
        assert_eq!(
            original.with_values(cx, vec![1, 2], palette_budget()),
            Err(VolumeTransportError::Cancelled)
        );
        assert!(matches!(
            CellPalette::try_new(cx, vec![visible], Arc::from([0_u32]), palette_budget()),
            Err(VolumeTransportError::Cancelled)
        ));
        assert!(matches!(
            AbsorbingGrid::try_new_indexed(cx, grid(), candidate, None, budget()),
            Err(VolumeTransportError::Cancelled)
        ));
    });
}
