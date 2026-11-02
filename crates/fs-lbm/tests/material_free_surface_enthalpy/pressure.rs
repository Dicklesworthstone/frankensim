//! Prescribed mechanical pressure on the original thermal VOF transaction.

use super::*;
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use std::cell::Cell;

fn with_liquid<R>(
    hot_viscosity: f64,
    reference: Option<f64>,
    run: impl FnOnce(MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = card(hot_viscosity, 1000.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        if let Some(reference) = reference {
            scaling = scaling.with_pressure_reference(reference).unwrap();
        }
        // A liquid slab has gas on both x faces and stationary lateral walls.
        let cells: Vec<_> = (0..8)
            .flat_map(|z| {
                (0..8).flat_map(move |y| {
                    (2..6).map(move |x| MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 1.0,
                        velocity_m_s: [0.0; 3],
                    })
                })
            })
            .collect();
        let flow = scaling
            .free_surface3_from_cells(domain(), &cells, cx)
            .unwrap();
        let liquid = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        run(liquid, &gate)
    })
}

fn heat() -> Vec<FreeSurfaceCellHeat> {
    (0..8)
        .flat_map(|z| {
            (0..8).flat_map(move |y| {
                (2..6).map(move |x| FreeSurfaceCellHeat {
                    coordinates: [x, y, z],
                    energy_j: 0.1,
                })
            })
        })
        .collect()
}

fn mean_velocity(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> f64 {
    let mut momentum = 0.0;
    let mut mass = 0.0;
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                if let Some(sample) = liquid.flow().sample_at_si(x, y, z).unwrap() {
                    let weight = sample.fill_fraction * sample.density_kg_m3;
                    momentum += weight * sample.velocity_m_s[0];
                    mass += weight;
                }
            }
        }
    }
    momentum / mass
}

#[test]
fn g1_pressure_drives_heated_liquid_without_inventing_thermal_energy() {
    with_liquid(0.02, Some(100000.0), |mut hot, gate| {
        let mut cold = hot.clone();
        let mut unloaded = hot.clone();
        let before_h = hot.total_enthalpy_j().unwrap();
        let before_mass = hot.flow().ledger_mass_kg().unwrap();
        let supplied = heat();
        let supplied_j: f64 = supplied.iter().map(|q| q.energy_j).sum();
        let load = |p: [f64; 3]| 100000.0 + if p[0] < 0.004 { 0.5 } else { -0.5 };
        for step in 0..6 {
            let source = if step == 0 { supplied.as_slice() } else { &[] };
            let report = hot
                .step_with_heat_and_gas_pressure(source, gate, 0.05, load)
                .unwrap();
            cold.step_with_heat_and_gas_pressure(&[], gate, 0.05, load)
                .unwrap();
            unloaded.step_with_heat(source, gate).unwrap();
            near(report.enthalpy_after_j, before_h + supplied_j, 2e-10);
            near(report.mass_after_kg, before_mass, 2e-12);
            assert!(report.balance_defect_j.abs() < 1e-8);
            assert_eq!(report.viscous_heat_j, 0.0);
            assert_eq!(report.body_heat_j, 0.0);
            if step == 0 {
                near(report.external_heat_j, supplied_j, 2e-14);
            } else {
                assert_eq!(report.external_heat_j, 0.0);
            }
        }
        assert!(
            mean_velocity(&hot) > 1e-5,
            "liquid must move toward lower pressure"
        );
        assert!(mean_velocity(&unloaded).abs() < 1e-13);
        assert_ne!(
            hot.flow().lattice().grid().state_bits(),
            cold.flow().lattice().grid().state_bits(),
            "the heated source viscosity must reach the pressure-driven trajectory"
        );
        near(cold.total_enthalpy_j().unwrap(), before_h, 2e-10);
        near(
            hot.total_enthalpy_j().unwrap(),
            unloaded.total_enthalpy_j().unwrap(),
            2e-10,
        );
        assert!(
            hot.material_at([3, 3, 3]).unwrap().dynamic_viscosity_pa_s()
                < cold
                    .material_at([3, 3, 3])
                    .unwrap()
                    .dynamic_viscosity_pa_s()
        );
        near(hot.physical_time_s(), 0.006, 2e-14);
    });
}

#[test]
fn g0_reference_pressure_reproduces_the_complete_heat_step() {
    with_liquid(0.02, Some(100000.0), |mut loaded, gate| {
        let mut original = loaded.clone();
        let supplied = heat();
        for step in 0..4 {
            let source = if step == 0 { supplied.as_slice() } else { &[] };
            let expected = original.step_with_heat(source, gate).unwrap();
            let actual = loaded
                .step_with_heat_and_gas_pressure(source, gate, 0.05, |_| 100000.0)
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(
                loaded.flow().lattice().grid().state_bits(),
                original.flow().lattice().grid().state_bits()
            );
            assert_eq!(loaded.physical_time_s(), original.physical_time_s());
            for x in 0..8 {
                assert_eq!(
                    loaded.specific_enthalpy_at([x, 3, 3]),
                    original.specific_enthalpy_at([x, 3, 3])
                );
                assert_eq!(
                    loaded.material_at([x, 3, 3]),
                    original.material_at([x, 3, 3])
                );
            }
        }
    });
}

#[test]
fn g0_pressure_gauge_shift_preserves_thermal_flow_bits() {
    let trajectory = |reference| {
        with_liquid(0.02, Some(reference), |mut liquid, gate| {
            let supplied = heat();
            for step in 0..4 {
                liquid
                    .step_with_heat_and_gas_pressure(
                        if step == 0 { &supplied } else { &[] },
                        gate,
                        0.05,
                        |p| reference + if p[0] < 0.004 { 0.5 } else { -0.5 },
                    )
                    .unwrap();
            }
            (
                liquid.flow().lattice().grid().state_bits(),
                liquid.total_enthalpy_j().unwrap().to_bits(),
            )
        })
    };
    assert_eq!(trajectory(100000.0), trajectory(100500.0));
}

#[test]
fn g4_late_pressure_and_cancellation_refusals_preserve_heat_flow_and_clock() {
    with_liquid(0.02, Some(100000.0), |mut liquid, gate| {
        let bits = liquid.flow().lattice().grid().state_bits();
        let before_h = liquid.total_enthalpy_j().unwrap();
        let before_source = liquid.material_at([3, 3, 3]).unwrap().clone();
        let source_heat = heat();
        for invalid in [0.0, -1.0, f64::NAN, 110000.0] {
            let calls = Cell::new(0);
            assert!(
                liquid
                    .step_with_heat_and_gas_pressure(&source_heat, gate, 0.05, |_| {
                        calls.set(calls.get() + 1);
                        if calls.get() < 5 { 100000.5 } else { invalid }
                    })
                    .is_err()
            );
            assert!(calls.get() >= 5);
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), before_h);
            assert_eq!(liquid.material_at([3, 3, 3]).unwrap(), &before_source);
            assert_eq!(liquid.physical_time_s(), 0.0);
        }
        for bound in [0.0, 1.0, f64::NAN] {
            let calls = Cell::new(0);
            assert!(
                liquid
                    .step_with_heat_and_gas_pressure(&source_heat, gate, bound, |_| {
                        calls.set(calls.get() + 1);
                        100000.0
                    })
                    .is_err()
            );
            assert_eq!(calls.get(), 0);
        }
        let cancelled = CancelGate::new_clock_free();
        assert!(
            liquid
                .step_with_heat_and_gas_pressure(&source_heat, &cancelled, 0.05, |_| {
                    cancelled.request();
                    100000.5
                })
                .is_err()
        );
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), before_h);
        assert_eq!(liquid.material_at([3, 3, 3]).unwrap(), &before_source);
        assert_eq!(liquid.physical_time_s(), 0.0);
        // A refused candidate leaves the ordinary accepted owner reusable.
        liquid
            .step_with_heat_and_gas_pressure(&source_heat, gate, 0.05, |_| 100000.0)
            .unwrap();
    });
    with_liquid(0.02, None, |mut liquid, gate| {
        let calls = Cell::new(0);
        assert!(
            liquid
                .step_with_heat_and_gas_pressure(&[], gate, 0.05, |_| {
                    calls.set(calls.get() + 1);
                    100000.0
                })
                .is_err()
        );
        assert_eq!(calls.get(), 0);
    });
}
