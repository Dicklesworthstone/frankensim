//! Exact original equipment/material ingress for the native radial head law.

use super::*;
use fs_airflow::graph::rheology::transient::head::FlowHeadLaw;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{
    power::{PumpDriveCurve, PumpLossRouting, PumpPowerPoint}, radial::RadialPumpHead,
};
use fs_material::liquid::ResolvedLiquidState;
use fs_qty::Power;
use std::sync::Arc;

fn inlet() -> Arc<ResolvedLiquidState> {
    let table = source::LiquidTable::parse(FIXED).unwrap();
    let law = LiquidTemperatureLaw::try_new(
        &table.card,
        &table.point(300.0).unwrap(),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    Arc::new(law.at_temperature(Temperature::new(300.0)).unwrap())
}
fn knot_curve() -> PumpHeadCurve {
    PumpHeadCurve::try_new(
        [(0.0, 0.02), (5e-6, 0.018), (1e-5, 0.0)]
            .map(|(q, h)| PumpHeadPoint {
                flow: VolumetricFlowRate::new(q),
                head: Length::new(h),
            })
            .to_vec(),
        10.0,
        (1e-6, 4e-6),
        SourceProvenance::new("authored numerical pump fixture", "radial-pump-knots"),
    )
    .unwrap()
}

#[test]
fn original_pump_and_material_are_retained_with_exact_knot_tangents() {
    let original = Arc::new(knot_curve());
    let material = inlet();
    let pump = Arc::new(RadialPumpHead::try_new(original.clone(), material.clone()).unwrap());
    assert!(std::ptr::eq(pump.curve(), original.as_ref()));
    assert!(std::ptr::eq(pump.inlet_material(), material.as_ref()));
    assert_eq!(
        pump.reference_density_kg_m3(),
        Some(material.density_kg_m3())
    );
    assert_eq!(pump.flow_interval_m3_s(), [0.0, 1e-5]);
    for (q, pressure, slope) in [
        (0.0, 200.0, -4e6),
        (2.5e-6, 190.0, -4e6),
        (5e-6, 180.0, -36e6),
        (7.5e-6, 90.0, -36e6),
        (1e-5, 0.0, -36e6),
    ] {
        let response = FlowHeadLaw::evaluate(pump.as_ref(), VolumetricFlowRate::new(q)).unwrap();
        close(response.pressure.value(), pressure, 1e-11);
        close(response.pressure_slope_pa_s_m3, slope, 1e-7);
        assert_eq!(
            pump.supply(VolumetricFlowRate::new(q)).unwrap(),
            original
                .evaluate(&material, VolumetricFlowRate::new(q))
                .unwrap()
        );
    }
    assert_eq!(pump.source(), pump.source());
    assert!(pump.evaluate(VolumetricFlowRate::new(-1e-12)).is_err());
    assert!(
        pump.evaluate(VolumetricFlowRate::new(1e-5_f64.next_up()))
            .is_err()
    );
    let unsupported = Arc::new(
        PumpHeadCurve::try_new(
            original.points().to_vec(),
            10.0,
            (1e-9, 2e-9),
            original.source().clone(),
        )
        .unwrap(),
    );
    assert!(RadialPumpHead::try_new(unsupported, material).is_err());
}

#[test]
fn optional_original_drive_support_and_power_refusals_are_not_bypassed() {
    let drive = PumpDriveCurve::try_new(
        [2e-6, 8e-6]
            .map(|q| PumpPowerPoint {
                flow: VolumetricFlowRate::new(q),
                shaft_input: Power::new(0.1),
                electrical_input: Power::new(0.2),
            })
            .to_vec(),
        (900.0, 1100.0),
        (1e-6, 4e-6),
        PumpLossRouting {
            pump_to_reservoir: 0.25,
            drive_to_reservoir: 0.5,
        },
        SourceProvenance::new(
            "separate authored input power fixture",
            "radial-original-drive",
        ),
    )
    .unwrap();
    let original = Arc::new(knot_curve().with_drive_curve(drive));
    let material = inlet();
    let pump = RadialPumpHead::try_new(original.clone(), material.clone()).unwrap();
    assert_eq!(pump.flow_interval_m3_s(), [2e-6, 8e-6]);
    assert!(pump.evaluate(VolumetricFlowRate::new(0.0)).is_err());
    for q in [2e-6, 5e-6, 8e-6] {
        let actual = pump.supply(VolumetricFlowRate::new(q)).unwrap();
        assert!(actual.drive_power.is_some());
        assert_eq!(
            actual,
            original
                .evaluate(&material, VolumetricFlowRate::new(q))
                .unwrap()
        );
        close(
            pump.evaluate(VolumetricFlowRate::new(q))
                .unwrap()
                .pressure
                .value(),
            actual.available_pressure_pa,
            0.0,
        );
    }
    let insufficient = PumpDriveCurve::try_new(
        [2e-6, 8e-6]
            .map(|q| PumpPowerPoint {
                flow: VolumetricFlowRate::new(q),
                shaft_input: Power::new(1e-4),
                electrical_input: Power::new(1e-4),
            })
            .to_vec(),
        (900.0, 1100.0),
        (1e-6, 4e-6),
        PumpLossRouting {
            pump_to_reservoir: 0.0,
            drive_to_reservoir: 0.0,
        },
        original.drive_curve().unwrap().source().clone(),
    )
    .unwrap();
    let limited = RadialPumpHead::try_new(
        Arc::new(knot_curve().with_drive_curve(insufficient)),
        material,
    )
    .unwrap();
    // At this original point hydraulic output is 180 Pa * 5e-6 m3/s,
    // exceeding the independently declared shaft input. No head-only fallback.
    assert!(limited.supply(VolumetricFlowRate::new(5e-6)).is_err());
    assert!(limited.evaluate(VolumetricFlowRate::new(5e-6)).is_err());
}
