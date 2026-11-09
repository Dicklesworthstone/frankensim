//! G1/G3/G4: original coupled Strang-midpoint physics and publication boundaries.

use super::super::super::tests::{
    Fixture, alternating, budget, limits, material, near, reference, snapshot, temperatures,
};
use super::super::adaptive::{ThermalJetAdaptiveControls, ThermalJetAdaptiveTolerance};
use super::super::*;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn controls() -> ThermalJetImplicitControls {
    ThermalJetImplicitControls {
        velocity_tolerance_m_s: 1e-11,
        cell_energy_tolerance_j: 1e-11,
        solver: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 1e-12,
            ..NewtonKrylovConfig::default()
        },
        ..ThermalJetImplicitControls::default()
    }
}

#[test]
fn g1_original_coupled_midpoint_refines_at_second_order_without_be_heat_loss() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let conduction = fixture.conduction(&caloric);
    for heat_transfer in [None, Some(&conduction)] {
        let expected = reference(if heat_transfer.is_some() { 1e-5 } else { 0.0 });
        let mut errors = Vec::new();
        for dt in [0.004, 0.002, 0.001] {
            let mut jet = ThermalCarreauJet::new(
                alternating(&gate),
                &caloric,
                heat_transfer,
                &temperatures(),
                100.0,
                &gate,
            )
            .unwrap();
            let report = jet
                .advance_implicit_with_method(
                    0.02,
                    budget(dt),
                    controls(),
                    ThermalJetImplicitMethod::StrangMidpoint,
                    &gate,
                )
                .unwrap();
            assert_eq!(report.method, ThermalJetImplicitMethod::StrangMidpoint);
            assert_eq!(report.outcome, SlenderJetOutcome::Completed);
            assert_eq!(report.end_time_s, 0.02);
            let actual = [
                1.0 - jet.flow().radius_squared_m2()[0],
                jet.flow().velocities_m_s().next().unwrap(),
                jet.cell_internal_energy_from_reference_j()[0],
                jet.cell_internal_energy_from_reference_j()[1],
                report.viscous_heat_j,
            ];
            errors.push(std::array::from_fn::<_, 5, _>(|i| {
                (actual[i] - expected[i]).abs()
            }));
            assert_eq!(report.implicit_kinetic_damping_j, 0.0);
            near(report.upwind_dissipation_j, 0.0, 1e-12);
            near(report.energy_balance_residual_j, 0.0, 1e-10);
            near(report.combined_energy_residual_j, 0.0, 3e-10);
            near(report.mass_balance_residual_kg, 0.0, 2e-13);
            assert!(report.maximum_velocity_residual_m_s <= controls().velocity_tolerance_m_s);
            assert!(report.maximum_cell_energy_residual_j <= controls().cell_energy_tolerance_j);
            let heat = jet.flow().cell_viscous_heat().unwrap();
            near(heat.total_heat_j, report.viscous_heat_j, 0.0);
            near(heat.cell_heat_j.iter().sum(), report.viscous_heat_j, 1e-10);
            assert!(heat.cell_heat_j.iter().all(|value| *value > 0.0));
            let sample = jet.flow().carreau_sample_at(0).unwrap().unwrap();
            assert_eq!(
                sample.temperature().temperature_k(),
                jet.temperatures_k().next().unwrap()
            );
            assert_eq!(sample.law(), &material().rheology);
            if dt == 0.001 {
                // Independent three-variable midpoint root between conservative
                // SSPRK2 half transports; no native jet/material code in oracle.
                let values = if heat_transfer.is_some() {
                    [
                        0.016757447056814544,
                        0.0858915171071192,
                        301.5957785200866,
                        330.5790043364528,
                        0.10364299253055925,
                    ]
                } else {
                    [
                        0.016757317860028802,
                        0.08588927717786052,
                        301.5645093370957,
                        330.6067070906854,
                        0.10364782777877382,
                    ]
                };
                near(actual[0], values[0], 2e-10);
                near(actual[1], values[1], 2e-10);
                let t: Vec<_> = jet.temperatures_k().collect();
                near(t[0], values[2], 2e-7);
                near(t[1], values[3], 2e-7);
                near(report.viscous_heat_j, values[4], 2e-9);
            }
        }
        for pair in errors.windows(2) {
            for i in 0..5 {
                assert!(pair[1][i] < 0.30 * pair[0][i], "errors={errors:?}");
            }
        }
    }
}

#[test]
fn g4_supported_midpoint_does_not_authorize_out_of_source_reflection() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1.0; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let conduction = fixture.conduction(&caloric);
    let radii: Vec<_> = (0..16)
        .map(|i| {
            if i % 2 == 0 {
                0.9_f64.sqrt()
            } else {
                1.1_f64.sqrt()
            }
        })
        .collect();
    let initial: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 359.9 } else { 280.1 })
        .collect();
    let make_flow = || {
        PeriodicSlenderJet::new_carreau(material(), 4.0, &radii, &[0.0; 16], limits(), &gate)
            .unwrap()
    };
    // Zero velocity/zero surface tension removes transport. The BE half solve
    // is exactly the supported caloric midpoint used by the full method.
    let mut midpoint = ThermalCarreauJet::new(
        make_flow(),
        &caloric,
        Some(&conduction),
        &initial,
        100.0,
        &gate,
    )
    .unwrap();
    midpoint
        .advance_implicit(0.01, budget(0.01), controls(), &gate)
        .unwrap();
    assert!(
        midpoint
            .temperatures_k()
            .all(|t| (280.0..=360.0).contains(&t))
    );
    let mut jet = ThermalCarreauJet::new(
        make_flow(),
        &caloric,
        Some(&conduction),
        &initial,
        100.0,
        &gate,
    )
    .unwrap();
    let before = snapshot(&jet);
    // Independent nonlinear-caloric reflection gives T_even about 271.85 K:
    // the original 280 K source boundary must refuse, without clipping.
    assert!(matches!(
        jet.advance_implicit_with_method(
            0.02,
            budget(0.02),
            controls(),
            ThermalJetImplicitMethod::StrangMidpoint,
            &gate,
        ),
        Err(ThermalJetImplicitError::Thermal(
            ThermalJetError::Thermal { .. }
        ))
    ));
    assert_eq!(snapshot(&jet), before);
    assert_eq!(
        jet.flow().cell_viscous_heat().unwrap().cell_heat_j,
        &[0.0; 16]
    );
}

#[test]
fn g4_midpoint_late_work_and_cancellation_preserve_heat_and_allow_clean_retry() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let conduction = fixture.conduction(&caloric);
    let mut jet = ThermalCarreauJet::new(
        alternating(&gate),
        &caloric,
        Some(&conduction),
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    let mut clean = ThermalCarreauJet::new(
        alternating(&gate),
        &caloric,
        Some(&conduction),
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    for flow in [&mut jet, &mut clean] {
        flow.advance_implicit_with_method(
            0.001,
            budget(0.0005),
            controls(),
            ThermalJetImplicitMethod::StrangMidpoint,
            &gate,
        )
        .unwrap();
    }
    let before = snapshot(&jet);
    let heat = jet.flow().cell_viscous_heat().unwrap().cell_heat_j.to_vec();
    let mut limited = budget(0.0005);
    limited.flow.maximum_substeps = 2;
    assert!(matches!(
        jet.advance_implicit_with_method(
            0.02,
            limited,
            controls(),
            ThermalJetImplicitMethod::StrangMidpoint,
            &gate,
        ),
        Err(ThermalJetImplicitError::Thermal(ThermalJetError::Flow(
            NonlinearJetError::SubstepBudget
        )))
    ));
    assert_eq!(snapshot(&jet), before);
    assert_eq!(jet.flow().cell_viscous_heat().unwrap().cell_heat_j, heat);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        jet.advance_implicit_with_method(
            0.001,
            budget(0.0005),
            controls(),
            ThermalJetImplicitMethod::StrangMidpoint,
            &cancelled,
        ),
        Err(ThermalJetImplicitError::Thermal(ThermalJetError::Flow(
            NonlinearJetError::Cancelled
        )))
    ));
    assert_eq!(snapshot(&jet), before);
    assert_eq!(jet.flow().cell_viscous_heat().unwrap().cell_heat_j, heat);
    for flow in [&mut jet, &mut clean] {
        flow.advance_implicit_with_method(
            0.002,
            budget(0.0005),
            controls(),
            ThermalJetImplicitMethod::StrangMidpoint,
            &gate,
        )
        .unwrap();
    }
    assert_eq!(snapshot(&jet), snapshot(&clean));
    assert_eq!(
        jet.flow().cell_viscous_heat().unwrap().cell_heat_j,
        clean.flow().cell_viscous_heat().unwrap().cell_heat_j,
    );
}

#[test]
fn g4_geometry_stop_restores_the_start_of_the_unfinished_split() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let velocity: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    let limited = SlenderJetLimits {
        minimum_radius_m: 0.998,
        ..limits()
    };
    let flow =
        PeriodicSlenderJet::new_carreau(material(), 4.0, &[1.0; 16], &velocity, limited, &gate)
            .unwrap();
    let mut jet =
        ThermalCarreauJet::new(flow, &caloric, None, &temperatures(), 100.0, &gate).unwrap();
    let report = jet
        .advance_implicit_with_method(
            0.02,
            budget(0.0005),
            controls(),
            ThermalJetImplicitMethod::StrangMidpoint,
            &gate,
        )
        .unwrap();
    assert!(matches!(report.outcome, SlenderJetOutcome::DomainLimit(_)));
    assert!(report.accepted_substeps > 0 && report.end_time_s < 0.02);
    assert_eq!(jet.flow().time_s(), report.end_time_s);
    near(report.energy_balance_residual_j, 0.0, 1e-10);
    near(report.combined_energy_residual_j, 0.0, 3e-10);
    near(
        jet.flow()
            .cell_viscous_heat()
            .unwrap()
            .cell_heat_j
            .iter()
            .sum(),
        report.viscous_heat_j,
        1e-10,
    );
    for ((energy, mass), temperature) in jet
        .cell_internal_energy_from_reference_j()
        .iter()
        .zip(jet.cell_masses_kg())
        .zip(jet.temperatures_k())
    {
        assert_eq!(
            temperature,
            caloric
                .temperature_at_relative_internal_energy(
                    jet.reference_temperature_k(),
                    energy / mass,
                )
                .unwrap(),
        );
    }
}

fn stiff_material() -> CarreauSlenderJetParameters {
    let mut p = material();
    p.rheology = ThermalCarreauYasuda::new(
        CarreauYasuda::new(80.0, 2.0, 10.0, 0.5, 2.0, [0.0, 1000.0]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::WilliamsLandelFerry {
                c1: 1.0,
                c2_k: 100.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    );
    p
}

fn adaptive_controls(scale: f64) -> ThermalJetAdaptiveControls {
    let absolute = |value| ThermalJetAdaptiveTolerance {
        absolute: value,
        relative: 0.0,
    };
    ThermalJetAdaptiveControls {
        area_m2: absolute(1.0),
        radius_m: absolute(1.0),
        velocity_m_s: absolute(scale * 1e-4),
        temperature_k: absolute(100.0),
        specific_energy_j_kg: absolute(1.0),
        cell_heat_j: absolute(scale * 1e-6),
        total_heat_j: absolute(scale * 1e-5),
        absolute_numerical_damping_j: 0.0,
        relative_numerical_damping: 0.0,
        ..ThermalJetAdaptiveControls::default()
    }
}

#[test]
fn g1_adaptive_midpoint_resolves_the_stiff_transient_that_energy_balance_alone_misses() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1.0; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let conduction = fixture.conduction(&caloric);
    let velocity: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    let make_flow = || {
        PeriodicSlenderJet::new_carreau(
            stiff_material(),
            4.0,
            &[1.0; 16],
            &velocity,
            limits(),
            &gate,
        )
        .unwrap()
    };
    let mut step_budget = budget(0.02);
    step_budget.energy_tolerance_j = 1e-8;
    let mut coarse = ThermalCarreauJet::new(
        make_flow(),
        &caloric,
        Some(&conduction),
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    let unresolved = coarse
        .advance_implicit_with_method(
            0.02,
            step_budget,
            controls(),
            ThermalJetImplicitMethod::StrangMidpoint,
            &gate,
        )
        .unwrap();
    assert_eq!(unresolved.accepted_substeps, 1);
    assert_eq!(unresolved.implicit_kinetic_damping_j, 0.0);
    near(unresolved.combined_energy_residual_j, 0.0, 1e-8);
    // A-stable midpoint is NOT L-stable. Its conserved coarse step can reverse
    // this stiff decaying mode while producing only a small fraction of its heat.
    assert!(coarse.flow().velocities_m_s().next().unwrap() < -0.12);
    assert!(unresolved.viscous_heat_j < 0.014);

    let mut resolved = ThermalCarreauJet::new(
        make_flow(),
        &caloric,
        Some(&conduction),
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    let report = resolved
        .advance_implicit_adaptive_with_method(
            0.02,
            step_budget,
            controls(),
            adaptive_controls(1.0),
            ThermalJetImplicitMethod::StrangMidpoint,
            &gate,
        )
        .unwrap();
    assert_eq!(
        report.evolution.method,
        ThermalJetImplicitMethod::StrangMidpoint
    );
    assert_eq!(report.evolution.outcome, SlenderJetOutcome::Completed);
    assert_eq!(report.evolution.end_time_s, 0.02);
    assert!(report.rejected_intervals > 0);
    assert!(report.accepted_intervals > 1 && report.accepted_intervals < 1000);
    assert_eq!(
        report.evolution.accepted_substeps,
        2 * report.accepted_intervals
    );
    assert!(report.evolution.attempted_substeps >= 3 * report.accepted_intervals);
    assert_eq!(report.evolution.implicit_kinetic_damping_j, 0.0);
    assert_eq!(report.numerical_damping_allowance_j, 0.0);
    assert_eq!(report.maximum_accepted_damping_ratio, 0.0);
    assert!(report.maximum_accepted_local.scaled_maximum <= 1.0);
    // Independently converged continuous ODE and equilibrium caloric balance.
    near(report.evolution.viscous_heat_j, 0.196349540849355, 2e-8);
    for temperature in resolved.temperatures_k() {
        near(temperature, 317.824637182205, 2e-5);
    }
    near(
        1.0 - resolved.flow().radius_squared_m2()[0],
        0.00041857904457331,
        1e-7,
    );
    near(report.evolution.energy_balance_residual_j, 0.0, 1e-8);
    near(report.evolution.combined_energy_residual_j, 0.0, 2e-8);
    near(
        resolved
            .flow()
            .cell_viscous_heat()
            .unwrap()
            .cell_heat_j
            .iter()
            .sum(),
        report.evolution.viscous_heat_j,
        1e-8,
    );

    // Inspect the transient before equilibrium: conservation by itself cannot
    // force temperature, velocity, or the spatial heat distribution to be right.
    let expected = [
        0.00016815639350325925,
        0.08670192754843184,
        314.6885328069432,
        318.5785118867534,
        0.1018853150505657,
    ];
    let mut errors = Vec::new();
    for scale in [1.0, 0.1] {
        let mut jet = ThermalCarreauJet::new(
            make_flow(),
            &caloric,
            Some(&conduction),
            &temperatures(),
            100.0,
            &gate,
        )
        .unwrap();
        let report = jet
            .advance_implicit_adaptive_with_method(
                0.0002,
                step_budget,
                controls(),
                adaptive_controls(scale),
                ThermalJetImplicitMethod::StrangMidpoint,
                &gate,
            )
            .unwrap();
        let t: Vec<_> = jet.temperatures_k().collect();
        let actual = [
            1.0 - jet.flow().radius_squared_m2()[0],
            jet.flow().velocities_m_s().next().unwrap(),
            t[0],
            t[1],
            report.evolution.viscous_heat_j,
        ];
        errors.push(std::array::from_fn::<_, 5, _>(|i| {
            (actual[i] - expected[i]).abs()
        }));
        assert!(report.rejected_intervals > 0);
        assert!(report.maximum_accepted_local.scaled_maximum <= 1.0);
        assert_eq!(report.evolution.implicit_kinetic_damping_j, 0.0);
        near(report.evolution.combined_energy_residual_j, 0.0, 2e-8);
        let heat = jet.flow().cell_viscous_heat().unwrap();
        near(heat.total_heat_j, report.evolution.viscous_heat_j, 0.0);
        near(
            heat.cell_heat_j.iter().sum(),
            report.evolution.viscous_heat_j,
            1e-8,
        );
        let sample = jet.flow().carreau_sample_at(0).unwrap().unwrap();
        assert_eq!(sample.temperature().temperature_k(), t[0]);
        assert_eq!(sample.law(), &stiff_material().rheology);
    }
    for i in 0..5 {
        assert!(errors[1][i] < 0.4 * errors[0][i], "errors={errors:?}");
    }
}

