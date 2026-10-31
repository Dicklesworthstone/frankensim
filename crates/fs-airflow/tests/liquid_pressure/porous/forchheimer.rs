//! G1/G3: one porous calibration feeds pressure, liquid properties and heat.
use super::*;
use fs_material::porous::permeability::ForchheimerCalibration;

fn inertial_skeleton(viscous: f64, beta: f64, domain: PorousGasSupport) -> DarcyGasBarrier {
    let original = skeleton(viscous, support());
    let calibration = ForchheimerCalibration::try_new(
        beta,
        [1.0, 0.0, 0.0],
        DarcyFlowConstraint::NormalPressureGradient,
        domain,
        "explicit synthetic liquid superficial-velocity calibration",
    )
    .unwrap();
    original
        .with_forchheimer_calibrations(&vec![calibration; original.layers().len()])
        .unwrap()
}

#[test]
fn calibrated_pore_inertia_changes_the_original_liquid_operating_point() {
    for viscosity_scale in [1.0, 2.0] {
        with_path(
            false,
            viscosity_scale,
            LiquidPropertyTemperature::BulkMean,
            |cx, path| {
                for beta in [0.0, 2000.0, 4000.0] {
                    let barrier = inertial_skeleton(5e9, beta, support());
                    let result = path
                        .solve_pressure_driven(
                            cx,
                            &[300.0, 300.0],
                            &mixed(&barrier),
                            Pressure::new(50.0),
                            local(),
                            pressure_budget(),
                        )
                        .unwrap();
                    // Same fluid in both sections; only the first has porous inertia.
                    // dp = (mu/rho)(a_p+a_channel)m + (b_p/rho)m².
                    let a: f64 = 0.002 * viscosity_scale * 1e10 / 1000.0;
                    let b = beta * 0.01 / (area().value() * area().value() * 1000.0);
                    let expected = 100.0 / (a + (a * a + 200.0 * b).sqrt());
                    close(result.mass_flow_kg_per_s, expected);
                    close(result.hydraulics.total_pressure_drop_pa, 50.0);
                    close(
                        result.hydraulics.total_dissipated_power_w,
                        50.0 * expected / 1000.0,
                    );
                    let pore = &result.hydraulics.regions[0];
                    let receipt = pore.porous.as_ref().unwrap();
                    assert_eq!(receipt.barrier, barrier);
                    close(
                        receipt.inlet_pressure_pa - receipt.outlet_pressure_pa,
                        0.002 * viscosity_scale * 5e9 * expected / 1000.0 + b * expected * expected,
                    );
                    close(
                        pore.loss.dissipated_power_w,
                        pore.loss.pressure_drop_pa * pore.loss.volume_flow_m3_s,
                    );
                    if beta > 0.0 {
                        assert!(expected < 50.0 / a);
                    }
                }
            },
        );
    }
}

#[test]
fn viscous_and_inertial_pore_heat_share_one_enthalpy_update() {
    let beta = 1e9;
    let barrier = inertial_skeleton(5e15, beta, support());
    let law = porous(&barrier, 100000.0);
    let sections = [
        LiquidResistanceSection {
            region: "first",
            resistance: law,
        },
        LiquidResistanceSection {
            region: "last",
            resistance: law,
        },
    ];
    for variable in [false, true] {
        with_path(
            variable,
            1.0,
            LiquidPropertyTemperature::BulkMean,
            |cx, path| {
                let report = path
                    .march_with_dissipation(cx, &[300.0, 300.0], &sections, local())
                    .unwrap();
                let mass: f64 = 0.001;
                let q = mass / 1000.0;
                let ua: f64 = 3.66 * 0.5 / 0.01 * 0.003;
                let e = (-ua / (mass * 2000.0)).exp();
                let f = 1.0 - e;
                let viscous_power_factor = 5e15 * q * q;
                let inertial_geometry = beta * 0.01 / area().value() / area().value();
                let inertial_power = 1000.0 * inertial_geometry * q * q * q;
                let (a, b) = if variable {
                    (0.005, 1e-5)
                } else {
                    (0.002, 0.0)
                };
                let (mut inlet, mut pressure, mut power) = (300.0, 0.0, 0.0);
                for _ in 0..2 {
                    // Independent exact regional exponential balance. The original
                    // affine mu(T_bulk) contributes once, and inertial power once.
                    let outlet = (e * inlet
                        + f * 300.0
                        + f / ua * (viscous_power_factor * (a - 0.5 * b * inlet) + inertial_power))
                        / (1.0 + 0.5 * f * viscous_power_factor * b / ua);
                    let viscosity = a - 0.5 * b * (inlet + outlet);
                    pressure += 5e15 * viscosity * q + 1000.0 * inertial_geometry * q * q;
                    power += viscous_power_factor * viscosity + inertial_power;
                    inlet = outlet;
                }
                close(report.thermal.regions[1].outlet.temperature_k(), inlet);
                close(report.hydraulics.total_pressure_drop_pa, pressure);
                close(report.hydraulics.total_dissipated_power_w, power);
                close(report.enthalpy_gain_w, mass * 2000.0 * (inlet - 300.0));
                close(
                    report.enthalpy_gain_w,
                    report.thermal.total_heat_rate_w + power,
                );
                close(report.balance_residual_w, 0.0);
                for row in &report.hydraulics.regions {
                    assert_eq!(row.porous.as_ref().unwrap().barrier, barrier);
                }
            },
        );
    }
}

#[test]
fn inertial_calibration_limits_remain_active_inside_liquid_pressure_trials() {
    let barrier = inertial_skeleton(
        5e9,
        2000.0,
        PorousGasSupport {
            maximum_superficial_speed_m_s: 0.05,
            ..support()
        },
    );
    with_path(
        false,
        1.0,
        LiquidPropertyTemperature::BulkMean,
        |cx, path| {
            // The bracket's upper trial exceeds the beta calibration, even though
            // the original k source and the eventual root admit that speed.
            assert!(matches!(
                path.solve_pressure_driven(
                    cx,
                    &[300.0, 300.0],
                    &mixed(&barrier),
                    Pressure::new(50.0),
                    local(),
                    pressure_budget(),
                ),
                Err(PressureDrivenError::Resistance { .. })
            ));
        },
    );
}
