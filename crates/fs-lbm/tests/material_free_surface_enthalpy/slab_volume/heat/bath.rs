//! Analytic fixed-bath heating drives actual native mass/interface updates.
use super::*;
use fs_lbm::material::free_surface::enthalpy::volume::heat::bath::{
    FreeSurfaceSlabBath, FreeSurfaceSlabBathHeatError,
};

const CP: f64 = 2000.0;
const DT: f64 = 0.001;
const HEIGHT: f64 = 0.002;

#[test]
fn g1_finite_bath_heating_and_cooling_pay_native_work_with_both_density_slopes() {
    for face in [Face3::XMin, Face3::ZMax] {
        for hot_density in [980.0, 1020.0] {
            for cooling in [false, true] {
                with_case(
                    Case {
                        face,
                        hot_density,
                        start_hot: cooling,
                        increment_j_kg: if cooling { -200_000.0 } else { 200_000.0 },
                        ..Case::default()
                    },
                    |old, _, _, gate, cx| {
                        let bits = old.flow().lattice().grid().state_bits();
                        let mass = old.flow().ledger_mass_kg().unwrap();
                        let rho0 = if cooling { hot_density } else { 1000.0 };
                        let rho1 = 0.5 * (1000.0 + hot_density);
                        let speed = ((rho0 - rho1) / rho1) * HEIGHT / DT;
                        let b = mass * (1.0 / 3.0 + 1.0 / 24.0);
                        let work = b * speed * speed;
                        // The applied power pays ONLY the independently computed
                        // mechanical work. Remaining thermal exchange has the
                        // exact constant-Cp half-relaxation to the bath.
                        let source = FreeSurfaceSlabBath {
                            reservoir_temperature_k: if cooling { 300.0 } else { 500.0 },
                            conductance_w_k: mass * CP * 2.0_f64.ln() / DT,
                            applied_power_w: work / DT,
                        };
                        let bracket = if cooling {
                            [-400_000.0, 0.0]
                        } else {
                            [0.0, 400_000.0]
                        };
                        let plan = old
                            .prepare_uniform_slab_bath_heat(
                                face,
                                source,
                                heat_controls(bracket),
                                gate,
                                cx,
                            )
                            .unwrap();
                        let initial_t = if cooling { 500.0 } else { 300.0 };
                        let mean = source.reservoir_temperature_k
                            + (initial_t - source.reservoir_temperature_k) * 0.5 / 2.0_f64.ln();
                        near(
                            plan.native().endpoint().state().temperature_k(),
                            400.0,
                            2e-12,
                        );
                        near(plan.thermal_exchange().mean_temperature_k, mean, 2e-12);
                        assert_eq!(
                            plan.thermal_exchange().initial.phase_curve_identity(),
                            plan.native().endpoint().state().phase_curve_identity()
                        );
                        let step = plan.advance(gate, cx).unwrap();
                        near(step.applied_heat_j, work, 2e-12);
                        near(
                            step.reservoir_heat_j,
                            mass * CP * (400.0 - initial_t),
                            2e-12,
                        );
                        near(step.heat.report.pressure_correction_work_j, work, 2e-8);
                        near(
                            step.heat.volume.report.final_height_m,
                            HEIGHT * rho0 / rho1,
                            2e-12,
                        );
                        near(step.heat.volume.report.mass_after_kg, mass, 2e-13);
                        assert!(step.heat_balance_defect_j.abs() <= 5e-11);
                        assert_eq!(step.heat.volume.liquid.accepted_steps(), 1);
                        assert_eq!(old.flow().lattice().grid().state_bits(), bits);
                        assert_eq!(old.accepted_steps(), 0);
                    },
                );
            }
        }
    }
}

#[test]
fn g3_repeated_bath_intervals_keep_incoming_motion_and_do_not_heat_projection_loss() {
    with_case(
        Case {
            increment_j_kg: 200_000.0,
            ..Case::default()
        },
        |old, _, _, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            let m0 = old.flow().scaling().mass_unit_kg();
            let b0 = inertia(HEIGHT, 0.001, m0, 64.0);
            let u0 = ((1000.0 - 990.0) / 990.0) * HEIGHT / DT;
            let source = FreeSurfaceSlabBath {
                reservoir_temperature_k: 500.0,
                conductance_w_k: mass * CP * 2.0_f64.ln() / DT,
                applied_power_w: b0 * u0 * u0 / DT,
            };
            let hot_plan = old
                .prepare_uniform_slab_bath_heat(
                    Face3::XMin,
                    source,
                    heat_controls([0.0, 400_000.0]),
                    gate,
                    cx,
                )
                .unwrap();
            let hot = hot_plan.advance(gate, cx).unwrap();
            let height1 = HEIGHT * 1000.0 / 990.0;
            let b1 = inertia(height1, 0.001, m0 * 990.0 / 1000.0, 64.0);
            let u1 = ((990.0 - 995.0) / 995.0) * height1 / DT;
            let work1 = b1 * u1 * (u1 - u0);
            let cooling = FreeSurfaceSlabBath {
                reservoir_temperature_k: 300.0,
                applied_power_w: work1 / DT,
                ..source
            };
            let cold_plan = hot
                .heat
                .volume
                .liquid
                .prepare_uniform_slab_bath_heat(
                    Face3::XMin,
                    cooling,
                    heat_controls([-150_000.0, 0.0]),
                    gate,
                    cx,
                )
                .unwrap();
            let cold = cold_plan.advance(gate, cx).unwrap();
            near(
                cold_plan.native().endpoint().state().temperature_k(),
                350.0,
                2e-12,
            );
            near(cold.heat.report.pressure_correction_work_j, work1, 2e-7);
            near(
                cold.heat.volume.report.final_height_m,
                HEIGHT * 1000.0 / 995.0,
                2e-12,
            );
            near(cold.heat.volume.report.mass_after_kg, mass, 2e-13);
            assert!(cold.heat_balance_defect_j.abs() <= 5e-11);
            assert_eq!(cold.heat.volume.liquid.accepted_steps(), 2);
            let stop_plan = cold
                .heat
                .volume
                .liquid
                .prepare_uniform_slab_bath_heat(
                    Face3::XMin,
                    FreeSurfaceSlabBath {
                        reservoir_temperature_k: 350.0,
                        conductance_w_k: 0.0,
                        applied_power_w: 0.0,
                    },
                    heat_controls([-10_000.0, 10_000.0]),
                    gate,
                    cx,
                )
                .unwrap();
            let stop = stop_plan.advance(gate, cx).unwrap();
            assert_eq!(stop.heat.report.net_specific_enthalpy_increment_j_kg, 0.0);
            assert_eq!(stop.heat.report.pressure_correction_work_j, 0.0);
            assert!(
                stop.heat
                    .volume
                    .projection
                    .projection
                    .report
                    .projection_loss_j
                    > 0.0
            );
            assert_eq!(stop.applied_heat_j, 0.0);
            assert_eq!(stop.reservoir_heat_j, 0.0);
            assert_eq!(stop.heat.volume.liquid.accepted_steps(), 3);
            near(stop.heat.volume.liquid.physical_time_s(), 0.253, 1e-14);
        },
    );
}

// Exact constant-density limit with the SAME original source card/native owner.
// Vary the native dt, not a separate heat-only integrator.
fn with_timestep(
    dt: f64,
    run: impl FnOnce(&MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate, &Cx<'_>),
) {
    let card = card(0.1, 1000.0);
    let chart = ResolvedHeatCapacityCurve::try_new_exact(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::with_quantity(
                "specific_heat_capacity",
                QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0])),
            ),
            density: &PropertyKey::with_quantity(
                "density",
                QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
            ),
            temperature_axis: "T",
            lower: &point(300.0),
            upper: &point(500.0),
            selection: fs_matdb::SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        0.0,
    )
    .unwrap();
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(300.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let cells: Vec<_> = (0..8)
        .flat_map(|z| {
            (0..8).flat_map(move |y| {
                (0..2).map(move |x| MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    velocity_m_s: [0.0; 3],
                })
            })
        })
        .collect();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, dt, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &cells, cx)
            .unwrap();
        let liquid = flow
            .with_enthalpy(
                &law,
                chart.curve(),
                FreeSurfaceEnthalpyControls::default(),
                |_| 0.0,
                &gate,
            )
            .unwrap();
        run(&liquid, &gate, cx);
    });
}

// Plans own endpoint source bundles. Recursion retains those real borrowed
// native states for this bounded 1/2/4-step refinement, without replacing them.
fn advance_partition(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    source: FreeSurfaceSlabBath,
    mut controls: FreeSurfaceSlabHeatControls,
    count: usize,
    gate: &CancelGate,
    cx: &Cx<'_>,
) -> (f64, f64, f64, f64) {
    controls.net_specific_enthalpy_bracket_j_kg[1] =
        400_000.0 - liquid.specific_enthalpy_at([0, 0, 0]).unwrap();
    let plan = liquid
        .prepare_uniform_slab_bath_heat(Face3::XMin, source, controls, gate, cx)
        .unwrap();
    let step = plan.advance(gate, cx).unwrap();
    assert_eq!(step.heat.report.pressure_correction_work_j, 0.0);
    near(step.heat.volume.report.final_height_m, HEIGHT, 2e-13);
    assert!(step.heat_balance_defect_j.abs() <= controls.maximum_heat_defect_j);
    if count == 1 {
        (
            plan.native().endpoint().state().temperature_k(),
            step.applied_heat_j,
            step.reservoir_heat_j,
            step.heat.volume.liquid.physical_time_s(),
        )
    } else {
        let (t, applied, bath, time) = advance_partition(
            &step.heat.volume.liquid,
            source,
            controls,
            count - 1,
            gate,
            cx,
        );
        (
            t,
            applied + step.applied_heat_j,
            bath + step.reservoir_heat_j,
            time,
        )
    }
}

#[test]
fn g3_timestep_refinement_retains_exact_fixed_bath_solution_in_zero_expansion_limit() {
    for count in [1, 2, 4] {
        let dt = DT / count as f64;
        with_timestep(dt, |old, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            let capacity = mass * CP;
            let source = FreeSurfaceSlabBath {
                reservoir_temperature_k: 450.0,
                conductance_w_k: capacity / DT,
                applied_power_w: capacity * 10.0 / DT,
            };
            let mut limits = heat_controls([0.0, 400_000.0]);
            limits.volume.volume.projection.duration_s = dt;
            let (t, applied, bath, time) = advance_partition(old, source, limits, count, gate, cx);
            let expected = 460.0 + (300.0 - 460.0) * (-1.0_f64).exp();
            near(t, expected, 3e-12);
            near(applied, source.applied_power_w * DT, 2e-13);
            near(bath, capacity * (expected - 300.0) - applied, 3e-12);
            near(time, DT, 1e-14);
        });
    }
}

#[test]
fn g1_zero_conductance_recovers_supplied_heat_and_zero_source_preserves_thermal_state() {
    with_case(
        Case {
            hot_density: 1000.0,
            increment_j_kg: 200_000.0,
            ..Case::default()
        },
        |old, _, _, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            for power in [0.0, mass * 200_000.0 / DT] {
                let plan = old
                    .prepare_uniform_slab_bath_heat(
                        Face3::XMin,
                        FreeSurfaceSlabBath {
                            reservoir_temperature_k: 400.0,
                            conductance_w_k: 0.0,
                            applied_power_w: power,
                        },
                        heat_controls([0.0, 400_000.0]),
                        gate,
                        cx,
                    )
                    .unwrap();
                let step = plan.advance(gate, cx).unwrap();
                assert_eq!(step.reservoir_heat_j, 0.0);
                assert_eq!(step.heat.report.pressure_correction_work_j, 0.0);
                near(
                    step.heat.report.represented_enthalpy_increment_j,
                    power * DT,
                    2e-13,
                );
                near(step.heat.volume.report.final_height_m, HEIGHT, 2e-13);
            }
        },
    );
}

#[test]
fn g4_bath_source_bracket_native_accuracy_and_cancellation_refuse_atomically() {
    with_case(
        Case {
            increment_j_kg: 200_000.0,
            ..Case::default()
        },
        |old, _, _, gate, cx| {
            let bits = old.flow().lattice().grid().state_bits();
            let h = old.total_enthalpy_j().unwrap().to_bits();
            let source = FreeSurfaceSlabBath {
                reservoir_temperature_k: 500.0,
                conductance_w_k: old.flow().ledger_mass_kg().unwrap() * CP / DT,
                applied_power_w: 0.0,
            };
            assert!(matches!(
                old.prepare_uniform_slab_bath_heat(
                    Face3::XMin,
                    source,
                    heat_controls([0.0, 1.0]),
                    gate,
                    cx
                ),
                Err(FreeSurfaceSlabBathHeatError::Native(
                    FreeSurfaceSlabHeatError::Balance { .. }
                ))
            ));
            assert!(
                old.prepare_uniform_slab_bath_heat(
                    Face3::XMin,
                    FreeSurfaceSlabBath {
                        conductance_w_k: -1.0,
                        ..source
                    },
                    heat_controls([0.0, 400_000.0]),
                    gate,
                    cx
                )
                .is_err()
            );
            let mut exhausted = heat_controls([0.0, 400_000.0]);
            exhausted.maximum_source_evaluations = 2;
            assert!(matches!(
                old.prepare_uniform_slab_bath_heat(Face3::XMin, source, exhausted, gate, cx),
                Err(FreeSurfaceSlabBathHeatError::Native(
                    FreeSurfaceSlabHeatError::Balance { evaluations: 2, .. }
                ))
            ));
            let mut strict = heat_controls([0.0, 400_000.0]);
            strict.volume.maximum_kinetic_mapping_change_j = 0.0;
            let rejected = old
                .prepare_uniform_slab_bath_heat(Face3::XMin, source, strict, gate, cx)
                .unwrap();
            assert!(rejected.advance(gate, cx).is_err());
            let plan = old
                .prepare_uniform_slab_bath_heat(
                    Face3::XMin,
                    source,
                    heat_controls([0.0, 400_000.0]),
                    gate,
                    cx,
                )
                .unwrap();
            gate.request();
            assert!(plan.advance(gate, cx).is_err());
            assert!(
                old.prepare_uniform_slab_bath_heat(
                    Face3::XMin,
                    source,
                    heat_controls([0.0, 400_000.0]),
                    gate,
                    cx
                )
                .is_err()
            );
            assert_eq!(old.flow().lattice().grid().state_bits(), bits);
            assert_eq!(old.total_enthalpy_j().unwrap().to_bits(), h);
            assert_eq!(old.accepted_steps(), 0);
        },
    );
}
