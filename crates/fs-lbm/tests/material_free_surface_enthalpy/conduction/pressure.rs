//! Pressure, conduction and finite-body storage share the original transaction.

use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceGasPressure, FreeSurfaceThermalBoundaries,
};
use std::cell::Cell;

fn with_pressure_liquid<R>(
    card: &MaterialCard,
    reference: Option<f64>,
    run: impl FnOnce(MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    let mut cells = Vec::new();
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..4 {
                cells.push(MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: if x == 3 { 0.5 } else { 1.0 },
                    velocity_m_s: [0.0; 3],
                });
            }
        }
    }
    with_cx(&gate, |cx| {
        let mut scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        if let Some(reference) = reference {
            scaling = scaling.with_pressure_reference(reference).unwrap();
        }
        let liquid = scaling
            .free_surface3_from_cells(domain(), &cells, cx)
            .unwrap()
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |p| match p {
                    [1, 3, 3] => 100000.0,
                    [2, 3, 3] => 300000.0,
                    _ => 200000.0,
                },
                &gate,
            )
            .unwrap();
        run(liquid, &gate)
    })
}

fn conducting_face() -> FreeSurfaceConductionFace {
    FreeSurfaceConductionFace {
        cells: [[1, 3, 3], [2, 3, 3]],
        ..face()
    }
}

fn boundary<'a>(pressure_pa: &'a dyn Fn([f64; 3]) -> f64) -> FreeSurfaceThermalBoundaries<'a> {
    FreeSurfaceThermalBoundaries {
        gas_pressure: Some(FreeSurfaceGasPressure {
            maximum_density_excursion: 0.05,
            pressure_pa,
        }),
        ..FreeSurfaceThermalBoundaries::default()
    }
}

#[test]
fn g1_pressure_shares_original_conduction_and_finite_body_heat() {
    let card = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_pressure_liquid(&card, Some(100000.0), |mut loaded, gate| {
        let mut ordinary = loaded.clone();
        let mut walls = [body::wall(&solid, &chart, &context, 450.0)];
        let mut original_walls = walls.clone();
        let before = loaded.total_enthalpy_j().unwrap() + walls[0].enthalpy_j();
        let expected = ordinary
            .step_with_body_conduction(
                &[conducting_face()],
                &[],
                controls(),
                &mut original_walls,
                &[body::contact()],
                gate,
            )
            .unwrap();
        let report = loaded
            .step_with_thermal_boundaries(
                &[conducting_face()],
                &[],
                controls(),
                &mut walls,
                &[body::contact()],
                boundary(&|_| 100000.5),
                false,
                gate,
            )
            .unwrap()
            .conduction;

        // The original thermal stage is identical; pressure then changes the
        // SAME native mass routes without adding an unphysical heat source.
        assert_eq!(report.thermal.substeps, expected.thermal.substeps);
        assert_eq!(report.last_applied, expected.last_applied);
        assert_eq!(report.thermal_endpoint, expected.thermal_endpoint);
        assert_eq!(
            walls[0].specific_enthalpy_j_kg(),
            original_walls[0].specific_enthalpy_j_kg()
        );
        assert_eq!(walls[0].material(), original_walls[0].material());
        assert!(report.bodies.as_ref().unwrap().heat_to_liquid_j > 0.0);
        assert_eq!(report.flow.external_heat_j, 0.0);
        assert_eq!(report.flow.viscous_heat_j, 0.0);
        assert!(report.flow.balance_defect_j.abs() < 1e-8);
        near(
            loaded.total_enthalpy_j().unwrap() + walls[0].enthalpy_j(),
            before,
            2e-10,
        );
        near(report.flow.mass_before_kg, report.flow.mass_after_kg, 2e-12);
        let velocity = loaded
            .flow()
            .sample_at_si(3, 3, 3)
            .unwrap()
            .unwrap()
            .velocity_m_s[0];
        assert!(
            velocity < -1e-6,
            "external pressure must push liquid away from gas"
        );
        assert_ne!(
            loaded.flow().lattice().grid().state_bits(),
            ordinary.flow().lattice().grid().state_bits()
        );
        assert_eq!(loaded.physical_time_s(), 0.001);
    });
}

#[test]
fn g0_reference_pressure_and_empty_boundary_preserve_thermal_trajectory() {
    let card = conducting_card(1.0);
    with_pressure_liquid(&card, Some(100000.0), |mut reference, gate| {
        let mut ordinary = reference.clone();
        let mut empty = reference.clone();
        let expected = ordinary
            .step_with_conduction(&[conducting_face()], &[], controls(), gate)
            .unwrap();
        let first = reference
            .step_with_thermal_boundaries(
                &[conducting_face()],
                &[],
                controls(),
                &mut [],
                &[],
                boundary(&|_| 100000.0),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        let second = empty
            .step_with_thermal_boundaries(
                &[conducting_face()],
                &[],
                controls(),
                &mut [],
                &[],
                FreeSurfaceThermalBoundaries::default(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        for (state, report) in [(&reference, &first), (&empty, &second)] {
            assert_eq!(report.flow, expected.flow);
            assert_eq!(report.last_applied, expected.last_applied);
            assert_eq!(
                state.flow().lattice().grid().state_bits(),
                ordinary.flow().lattice().grid().state_bits()
            );
            assert_eq!(
                state.total_enthalpy_j().unwrap(),
                ordinary.total_enthalpy_j().unwrap()
            );
            assert_eq!(state.physical_time_s(), ordinary.physical_time_s());
        }
    });
}

#[test]
fn g4_late_pressure_refusal_discards_completed_conduction_and_wall_heat() {
    let card = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    with_pressure_liquid(&card, Some(100000.0), |mut liquid, gate| {
        let mut walls = [body::wall(&solid, &chart, &context, 450.0)];
        let before_bits = liquid.flow().lattice().grid().state_bits();
        let before_h = liquid.total_enthalpy_j().unwrap();
        let before_wall = walls[0].specific_enthalpy_j_kg();
        let before_source = walls[0].material().clone();
        for invalid in [0.0, f64::NAN, 110000.0] {
            let calls = Cell::new(0);
            let pressure = |_| {
                calls.set(calls.get() + 1);
                if calls.get() < 5 { 100000.5 } else { invalid }
            };
            assert!(
                liquid
                    .step_with_thermal_boundaries(
                        &[conducting_face()],
                        &[],
                        controls(),
                        &mut walls,
                        &[body::contact()],
                        boundary(&pressure),
                        false,
                        gate,
                    )
                    .is_err()
            );
            assert!(
                calls.get() >= 5,
                "refusal must occur after native pressure reconstruction starts"
            );
            assert_eq!(liquid.flow().lattice().grid().state_bits(), before_bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), before_h);
            assert_eq!(walls[0].specific_enthalpy_j_kg(), before_wall);
            assert_eq!(walls[0].material(), &before_source);
            assert_eq!(liquid.physical_time_s(), 0.0);
        }
        let cancelled = CancelGate::new_clock_free();
        let pressure = |_| {
            cancelled.request();
            100000.5
        };
        assert!(
            liquid
                .step_with_thermal_boundaries(
                    &[conducting_face()],
                    &[],
                    controls(),
                    &mut walls,
                    &[body::contact()],
                    boundary(&pressure),
                    false,
                    &cancelled,
                )
                .is_err()
        );
        assert_eq!(liquid.flow().lattice().grid().state_bits(), before_bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), before_h);
        assert_eq!(walls[0].specific_enthalpy_j_kg(), before_wall);
        assert_eq!(walls[0].material(), &before_source);
        liquid
            .step_with_thermal_boundaries(
                &[conducting_face()],
                &[],
                controls(),
                &mut walls,
                &[body::contact()],
                boundary(&|_| 100000.5),
                false,
                gate,
            )
            .unwrap();
    });
    with_pressure_liquid(&card, None, |mut liquid, gate| {
        let calls = Cell::new(0);
        let pressure = |_| {
            calls.set(calls.get() + 1);
            100000.0
        };
        assert!(
            liquid
                .step_with_thermal_boundaries(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    boundary(&pressure),
                    false,
                    gate,
                )
                .is_err()
        );
        assert_eq!(calls.get(), 0);
        assert_eq!(liquid.physical_time_s(), 0.0);
    });
}
