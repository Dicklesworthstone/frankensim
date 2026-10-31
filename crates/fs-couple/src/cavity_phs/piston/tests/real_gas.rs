//! Independent real-fluid motion and original-law trajectory energy.
//! Synthetic source parameters; no experimental gas or actuator calibration.

use super::*;
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::fluid_transport::FluidTransport;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    energy::parcel::IsochoricRealGasParcel,
};
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

pub(in crate::cavity_phs::piston) fn real_parcel(t: f64, p: f64) -> IsochoricRealGasParcel {
    let model = ThermochemicalPengRobinsonModel::new(
        gas(0.75, 0.001).model().clone(),
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), p)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    IsochoricRealGasParcel::new(model, point, 0.001, PengRobinsonPhaseV1::Vapor).unwrap()
}

pub(in crate::cavity_phs::piston) fn real_controls() -> PistonSolve {
    PistonSolve {
        volume_bracket_m3: [0.00095, 0.00105],
        gas: GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 450.0,
            energy_tolerance_j_kg: 1e-7,
            maximum_iterations: 128,
        },
        momentum_tolerance_n_s: 1e-10,
        energy_tolerance_j: 1e-7,
        maximum_iterations: 128,
    }
}

fn input(heat_j: f64) -> PistonInput {
    PistonInput {
        duration_s: 0.00005,
        heat_j,
        exterior_pressure_pa: 4.5e6,
        external_force_n: -100.0,
    }
}

fn piston() -> GasPiston<IsochoricRealGasParcel> {
    GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, 0.4).unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_dense_gas_motion_matches_independent_joint_energy_and_impulse_roots() {
    // Independent 75-digit PR/NASA solve of the two coupled equations:
    // um-u0+pm*DeltaV/(2m)=Q/(2m), M*(v1-v0)=dt*(A*(pm-pe)+F).
    // The endpoint separately satisfies u1-u0=2*(um-u0).
    for (q, volume, velocity, tm, pm, t1, p1, work) in [
        (
            100.0,
            0.001_000_231_967_737_106_3,
            0.527_870_948_425_014_1,
            301.072_785_646_489_55,
            5_021_483.793_700_056,
            302.145_285_152_974_2,
            5_042_952.004_407_114,
            1.164_822_232_540_327_3,
        ),
        (
            0.0,
            0.001_000_230_571_603_423_2,
            0.522_286_413_692_779_8,
            299.984_540_921_716_86,
            4_999_145.654_771_119,
            299.969_082_338_507_8,
            4_998_291.517_258_860_5,
            1.152_661_029_366_674_7,
        ),
    ] {
        let mut wall = piston();
        let initial = wall.clone();
        let report = wall.step(input(q), real_controls(), || false).unwrap();
        near(wall.gas().volume_m3(), volume, 2e-15);
        near(wall.velocity_m_s(), velocity, 2e-9);
        near(report.midpoint.temperature_k(), tm, 2e-8);
        near(report.midpoint.pressure_pa(), pm, 5e-4);
        near(wall.gas().state().temperature_k(), t1, 2e-8);
        near(wall.gas().state().pressure_pa(), p1, 5e-4);
        near(report.gas.work_by_gas_j, work, 2e-7);
        near(report.gas.stored_energy_change_j, q - work, 2e-7);
        assert_eq!(wall.gas().mass_kg(), initial.gas().mass_kg());
        assert_eq!(
            wall.gas().density_kg_m3(),
            wall.gas().mass_kg() / wall.gas().volume_m3()
        );
        near(
            report.midpoint.density_kg_m3(),
            wall.gas().mass_kg() / f64::midpoint(initial.gas().volume_m3(), volume),
            2e-10,
        );
        assert_eq!(wall.gas().point().axes().get("specimen-lot"), Some(&42.0));
        assert_eq!(wall.gas().model(), initial.gas().model());
        assert_eq!(
            report.midpoint.transport_state().unwrap().density_kg_m3(),
            report.midpoint.acoustic_state().unwrap().density_kg_m3()
        );
        assert!(report.momentum_residual_n_s.abs() <= real_controls().momentum_tolerance_n_s);
        assert!(report.total_energy_residual_j.abs() <= real_controls().energy_tolerance_j);
        assert!(report.mechanical_energy_residual_j.abs() <= real_controls().energy_tolerance_j);
    }
}

#[test]
fn g3_real_gas_midpoint_reverses_and_closes_a_continuing_energy_account() {
    let initial = piston();
    let mut forward = initial.clone();
    forward
        .step(input(100.0), real_controls(), || false)
        .unwrap();
    let mut reverse = GasPiston::new(
        forward.gas().clone(),
        forward.mass_kg(),
        forward.area_m2(),
        -forward.velocity_m_s(),
    )
    .unwrap();
    reverse
        .step(input(-100.0), real_controls(), || false)
        .unwrap();
    near(reverse.gas().volume_m3(), initial.gas().volume_m3(), 5e-15);
    near(reverse.velocity_m_s(), -initial.velocity_m_s(), 1e-8);
    near(reverse.gas().state().temperature_k(), 300.0, 5e-8);
    near(reverse.gas().state().pressure_pa(), 5e6, 1e-3);

    let mut wall = initial.clone();
    let mut heat = 0.0;
    for i in 0..16 {
        let q = if i == 0 { 100.0 } else { 0.0 };
        wall.step(input(q), real_controls(), || false).unwrap();
        heat += q;
    }
    let stored = wall.gas().mass_kg()
        * wall
            .gas()
            .model()
            .specific_internal_energy_change_between_densities(
                initial.gas().density_kg_m3(),
                initial.gas().state().temperature_k(),
                wall.gas().density_kg_m3(),
                wall.gas().state().temperature_k(),
            )
            .unwrap();
    let kinetic = 0.5
        * wall.mass_kg()
        * (wall.velocity_m_s() - initial.velocity_m_s())
        * (wall.velocity_m_s() + initial.velocity_m_s());
    let dv = wall.gas().volume_m3() - initial.gas().volume_m3();
    let outside =
        input(0.0).exterior_pressure_pa * dv - input(0.0).external_force_n * (dv / wall.area_m2());
    near(
        stored + kinetic + outside,
        heat,
        16.0 * real_controls().energy_tolerance_j,
    );
    assert_eq!(wall.gas().mass_kg(), initial.gas().mass_kg());
}
