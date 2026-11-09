//! Coupled endpoint contact on the original moving slab and finite caloric body.

use super::*;
use fs_lbm::material::free_surface::enthalpy::volume::heat::reservoir::coupled::{
    FreeSurfaceSlabCoupledReservoirHeatControls, FreeSurfaceSlabCoupledReservoirHeatError,
};

fn coupled_limits(bracket: [f64; 2]) -> FreeSurfaceSlabCoupledReservoirHeatControls {
    FreeSurfaceSlabCoupledReservoirHeatControls {
        heat: heat_controls(bracket),
        maximum_reservoir_energy_defect_j: 5e-11,
        maximum_total_energy_defect_j: 1e-10,
        maximum_contact_heat_defect_j: 5e-10,
    }
}

#[test]
fn g1_joint_contact_pays_actual_work_with_heating_cooling_and_either_density_slope() {
    for hot_density in [980.0, 1020.0] {
        for cooling in [false, true] {
            with_case(
                Case {
                    hot_density,
                    start_hot: cooling,
                    increment_j_kg: if cooling { -200_000.0 } else { 200_000.0 },
                    datum: 1e8,
                    ..Case::default()
                },
                |old, _, _, gate, cx| {
                    let mass = old.flow().ledger_mass_kg().unwrap();
                    let t0 = if cooling { 500.0 } else { 300.0 };
                    let tr0 = if cooling { 300.0 } else { 500.0 };
                    let rho0 = if cooling { hot_density } else { 1000.0 };
                    let rho1 = 0.5 * (1000.0 + hot_density);
                    let u = ((rho0 - rho1) / rho1) * HEIGHT / DT;
                    let work = mass * (1.0 / 3.0 + 1.0 / 24.0) * u * u;
                    let heat = mass * CP * (400.0 - t0) + work;
                    let reservoir_capacity = 4.0 * mass * 1000.0;
                    let tr1 = tr0 - heat / reservoir_capacity;
                    // Independent manufactured solution of the JOINT endpoint
                    // equation, including nonzero measured pressure work.
                    let g = heat / (DT * (tr1 - 400.0));
                    let chart = reservoir_curve(1000.0, 2e8);
                    let original = body(
                        &chart,
                        2e8 + 1000.0 * (tr0 - 300.0),
                        4.0 * mass,
                    );
                    let bracket = if cooling {
                        [-400_000.0, 0.0]
                    } else {
                        [0.0, 400_000.0]
                    };
                    let bits = old.flow().lattice().grid().state_bits();
                    let plan = old
                        .prepare_uniform_slab_coupled_reservoir_heat(
                            Face3::XMin, original, g, coupled_limits(bracket), gate, cx,
                        )
                        .unwrap();
                    let step = plan.advance(gate, cx).unwrap();
                    let liquid_t = plan.native().endpoint().state().temperature_k();
                    near(liquid_t, 400.0, 3e-12);
                    near(step.reservoir.initial.temperature_k(), tr1, 3e-12);
                    near(step.report.heat_to_liquid_j, heat, 3e-12);
                    near(step.heat.report.pressure_correction_work_j, work, 3e-8);
                    near(step.heat.volume.report.final_height_m, HEIGHT * rho0 / rho1, 3e-12);
                    near(step.heat.volume.report.mass_after_kg, mass, 3e-13);
                    let contact = g * DT * (step.reservoir.initial.temperature_k() - liquid_t);
                    assert!((step.report.heat_to_liquid_j - contact).abs() <= 5e-10);
                    assert!(step.report.contact_heat_defect_j.abs() <= 5e-10);
                    // Native scalar heat and coupled scalar contact have
                    // different residuals, even though they share one root.
                    let scalar_heat = step.heat.report.scalar_heat_balance_defect_j;
                    let scalar_contact = step.report.scalar_contact_heat_defect_j;
                    let rounding = 8.0 * f64::EPSILON * (heat.abs() + work.abs());
                    assert!(scalar_heat.abs() <= rounding);
                    let planned_contact =
                        g * DT * (plan.reservoir_endpoint().temperature_k() - liquid_t);
                    assert!(
                        (scalar_contact
                            - (step.heat.report.supplied_heat_j + scalar_heat - planned_contact))
                            .abs()
                            <= rounding
                    );
                    assert!(step.report.reservoir_energy_defect_j.abs() <= 5e-11);
                    assert!(step.report.total_energy_defect_j.abs() <= 1e-10);
                    assert!(std::ptr::eq(step.reservoir.curve, original.curve));
                    assert_eq!(step.reservoir.mass_kg.to_bits(), original.mass_kg.to_bits());
                    assert_eq!(step.reservoir.initial, plan.reservoir_endpoint());
                    assert_eq!(old.flow().lattice().grid().state_bits(), bits);
                    assert_eq!(old.accepted_steps(), 0);
                },
            );
        }
    }
}

#[test]
fn g1_original_latent_reservoir_and_zero_contact_keep_the_joint_endpoints() {
    with_case(
        Case { increment_j_kg: 50_000.0, ..Case::default() },
        |old, _, _, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            let chart = EquilibriumEnthalpyPhaseCurve::try_new(
                ContentHash([0x82; 32]),
                [
                    (0.0, 300.0, 0.0),
                    (50_000.0, 350.0, 0.0),
                    (150_000.0, 350.0, 1.0),
                    (300_000.0, 500.0, 1.0),
                ]
                .map(|(h, t, f)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: h,
                    temperature_k: t,
                    liquid_mass_fraction: f,
                    bulk_density_kg_m3: 1200.0 - 100.0 * f,
                })
                .to_vec(),
            )
            .unwrap();
            let original = body(&chart, 100_000.0, 4.0 * mass);
            let rho1 = 997.5;
            let u = (1000.0 - rho1) / rho1 * HEIGHT / DT;
            let work = 0.375 * mass * u * u;
            let heat = mass * CP * 25.0 + work;
            let plan = old
                .prepare_uniform_slab_coupled_reservoir_heat(
                    Face3::XMin,
                    original,
                    heat / (DT * 25.0),
                    coupled_limits([0.0, 100_000.0]),
                    gate,
                    cx,
                )
                .unwrap();
            let step = plan.advance(gate, cx).unwrap();
            near(plan.native().endpoint().state().temperature_k(), 325.0, 3e-12);
            assert_eq!(step.reservoir.initial.temperature_k(), 350.0);
            assert_eq!(step.reservoir.initial.phase(), SolidLiquidPhase::SolidLiquid);
            let fraction = (100_000.0 - heat / (4.0 * mass) - 50_000.0) / 100_000.0;
            near(step.reservoir.initial.liquid_mass_fraction(), fraction, 3e-12);
            near(step.heat.report.pressure_correction_work_j, work, 3e-8);
            assert!(step.report.total_energy_defect_j.abs() <= 1e-10);
            for temperature in [300.0, 500.0] {
                let sensible = reservoir_curve(1000.0, 0.0);
                let reservoir = body(&sensible, 1000.0 * (temperature - 300.0), mass);
                let no_heat = old
                    .prepare_uniform_slab_coupled_reservoir_heat(
                        Face3::XMin,
                        reservoir,
                        if temperature == 300.0 { 1.0 } else { 0.0 },
                        coupled_limits([0.0, 400_000.0]),
                        gate,
                        cx,
                    )
                    .unwrap();
                let stopped = no_heat.advance(gate, cx).unwrap();
                assert_eq!(stopped.report.heat_to_liquid_j, 0.0);
                assert_eq!(stopped.report.endpoint_contact_heat_j, 0.0);
                assert_eq!(stopped.reservoir.initial, reservoir.initial);
                assert_eq!(stopped.heat.report.net_specific_enthalpy_increment_j_kg, 0.0);
                assert_eq!(stopped.heat.report.pressure_correction_work_j, 0.0);
            }
        },
    );
}

// Recursion retains each immutable source-owning plan for the actual native
// liquid borrowed by the next interval. Neither source mass is reconstructed.
#[allow(clippy::too_many_arguments)]
fn coupled_partitions(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    reservoir: FinitePhaseBody<'_>,
    g: f64,
    cp: f64,
    count: usize,
    mut controls: FreeSurfaceSlabCoupledReservoirHeatControls,
    gate: &CancelGate,
    cx: &Cx<'_>,
) -> (f64, f64, f64, f64) {
    let dt = controls.heat.volume.volume.projection.duration_s;
    let h = liquid.specific_enthalpy_at([0, 0, 0]).unwrap();
    // Authored test bracket limits a trial's interface speed and reservoir
    // debit. A trial outside the ORIGINAL support would correctly refuse.
    controls.heat.net_specific_enthalpy_bracket_j_kg =
        [0.0, (cp * 200.0 - h).min(cp * 200_000.0 * dt)];
    let plan = liquid
        .prepare_uniform_slab_coupled_reservoir_heat(
            Face3::XMin, reservoir, g, controls, gate, cx,
        )
        .unwrap();
    let step = plan.advance(gate, cx).unwrap();
    assert!(std::ptr::eq(step.reservoir.curve, reservoir.curve));
    assert_eq!(step.reservoir.mass_kg, reservoir.mass_kg);
    assert!(step.report.total_energy_defect_j.abs() <= controls.maximum_total_energy_defect_j);
    assert!(step.report.contact_heat_defect_j.abs() <= controls.maximum_contact_heat_defect_j);
    if count == 1 {
        (
            plan.native().endpoint().state().temperature_k(),
            step.reservoir.initial.temperature_k(),
            step.heat.volume.liquid.physical_time_s(),
            step.heat.report.pressure_correction_work_j,
        )
    } else {
        let (t, r, time, w) = coupled_partitions(
            &step.heat.volume.liquid,
            step.reservoir,
            g,
            cp,
            count - 1,
            controls,
            gate,
            cx,
        );
        (t, r, time, w + step.heat.report.pressure_correction_work_j)
    }
}

#[test]
fn g3_zero_expansion_matches_the_backward_euler_finite_pair_and_refines() {
    let equilibrium = (300.0 + 2.0 * 500.0) / 3.0;
    let exact = equilibrium + (300.0 - equilibrium) * (-1.5_f64).exp();
    let mut errors = Vec::new();
    for count in [1, 2, 4, 8] {
        let dt = DT / count as f64;
        with_parameters(dt, 1000.0, CP, |old, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            let chart = reservoir_curve(CP, 0.0);
            let reservoir = body(&chart, CP * 200.0, 2.0 * mass);
            let mut controls = coupled_limits([0.0, CP * 200.0]);
            controls.heat.volume.volume.projection.duration_s = dt;
            let (t, r, time, w) = coupled_partitions(
                old, reservoir, mass * CP / DT, CP, count, controls, gate, cx,
            );
            let attenuation = (1.0 + 1.5 / count as f64).powi(-(count as i32));
            near(t, equilibrium + (300.0 - equilibrium) * attenuation, 3e-12);
            near(r, equilibrium + (500.0 - equilibrium) * attenuation, 3e-12);
            near(time, DT, 1e-14);
            assert_eq!(w, 0.0);
            errors.push((t - exact).abs());
        });
    }
    assert!(errors.windows(2).all(|pair| pair[1] < pair[0]));
    assert!(errors[3] < 0.2 * errors[0]);
}

#[test]
fn g3_nonzero_work_and_finite_reservoir_feedback_refine_together() {
    let cp = 2e-6;
    // This independent continuous native-face-work reference is shared with
    // the split-entry test, not the production endpoint/contact/root kernel.
    let reference = continuous_reference(cp, 32768);
    let coarse_reference = continuous_reference(cp, 16384);
    assert!((reference[0] - coarse_reference[0]).abs() < 1e-5);
    assert!((reference[1] - coarse_reference[1]).abs() < 1e-5);
    let mut errors = Vec::new();
    for count in [2, 4, 8] {
        let dt = DT / count as f64;
        with_parameters(dt, 980.0, cp, |old, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            let chart = reservoir_curve(cp, 0.0);
            let reservoir = body(&chart, cp * 200.0, 2.0 * mass);
            let mut controls = coupled_limits([0.0, cp * 200.0]);
            controls.heat.volume.volume.projection.duration_s = dt;
            controls.heat.maximum_heat_defect_j = 1e-16;
            controls.maximum_reservoir_energy_defect_j = 1e-16;
            controls.maximum_total_energy_defect_j = 2e-16;
            controls.maximum_contact_heat_defect_j = 2e-16;
            let (t, r, time, w) = coupled_partitions(
                old, reservoir, mass * cp / DT, cp, count, controls, gate, cx,
            );
            near(time, DT, 1e-14);
            assert!(w > 5e-9 && t > 360.0 && r < 450.0);
            let defect = mass * cp * (t - 300.0) + 2.0 * mass * cp * (r - 500.0) + w;
            assert!(defect.abs() < 2e-15);
            errors.push((t - reference[0]).abs() + (r - reference[1]).abs());
        });
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1]);
    assert!(errors[2] < 0.55 * errors[0]);
}

#[test]
fn g4_joint_source_native_and_cancellation_refusals_preserve_both_owners() {
    with_case(
        Case { increment_j_kg: 200_000.0, ..Case::default() },
        |old, _, _, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            let chart = reservoir_curve(CP, 0.0);
            let reservoir = body(&chart, CP * 200.0, 2.0 * mass);
            let original_state = reservoir.initial;
            let bits = old.flow().lattice().grid().state_bits();
            let g = mass * CP / DT;
            let limits = coupled_limits([0.0, 400_000.0]);
            assert!(matches!(
                old.prepare_uniform_slab_coupled_reservoir_heat(
                    Face3::XMin, reservoir, g, coupled_limits([0.0, 1.0]), gate, cx,
                ),
                Err(FreeSurfaceSlabCoupledReservoirHeatError::Native(
                    FreeSurfaceSlabHeatError::Balance { .. }
                ))
            ));
            let mut short = limits;
            short.heat.maximum_source_evaluations = 2;
            assert!(old.prepare_uniform_slab_coupled_reservoir_heat(
                Face3::XMin, reservoir, g, short, gate, cx,
            ).is_err());
            let unsupported = body(&chart, -50.0 * CP, 2.0 * mass);
            assert!(matches!(
                old.prepare_uniform_slab_coupled_reservoir_heat(
                    Face3::XMin, unsupported, g, limits, gate, cx,
                ),
                Err(FreeSurfaceSlabCoupledReservoirHeatError::Caloric(_))
            ));
            let other = reservoir_curve(CP, 1.0);
            let mismatch = FinitePhaseBody { curve: &other, ..reservoir };
            assert!(matches!(
                old.prepare_uniform_slab_coupled_reservoir_heat(
                    Face3::XMin, mismatch, 0.0, limits, gate, cx,
                ),
                Err(FreeSurfaceSlabCoupledReservoirHeatError::Caloric(_))
            ));
            assert!(old.prepare_uniform_slab_coupled_reservoir_heat(
                Face3::XMin, FinitePhaseBody { mass_kg: 0.0, ..reservoir },
                g, limits, gate, cx,
            ).is_err());
            assert!(old.prepare_uniform_slab_coupled_reservoir_heat(
                Face3::XMin, reservoir, -g, limits, gate, cx,
            ).is_err());
            let mut strict = limits;
            strict.heat.volume.maximum_kinetic_mapping_change_j = 0.0;
            let rejected = old.prepare_uniform_slab_coupled_reservoir_heat(
                Face3::XMin, reservoir, g, strict, gate, cx,
            ).unwrap();
            assert!(rejected.advance(gate, cx).is_err());
            let plan = old.prepare_uniform_slab_coupled_reservoir_heat(
                Face3::XMin, reservoir, g, limits, gate, cx,
            ).unwrap();
            gate.request();
            assert!(plan.advance(gate, cx).is_err());
            assert!(old.prepare_uniform_slab_coupled_reservoir_heat(
                Face3::XMin, reservoir, g, limits, gate, cx,
            ).is_err());
            assert_eq!(reservoir.initial, original_state);
            assert_eq!(old.flow().lattice().grid().state_bits(), bits);
            assert_eq!(old.accepted_steps(), 0);
        },
    );
}
