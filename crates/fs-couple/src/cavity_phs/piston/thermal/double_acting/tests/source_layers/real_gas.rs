//! Independent real-gas sources, one finite thermal store and one momentum.

use super::{body, exchange, layer, material, point};
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput,
    real_gas::{real_controls, real_pair},
};
use crate::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
use crate::cavity_phs::piston::thermal::{HeatExchangeSolve, ThermalPistonError};
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::UniformEnthalpyStepInput;
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;
use fs_ports::{ModuleWindowAudit, PortedModule};

fn input() -> DoublePistonInput {
    DoublePistonInput {
        duration_s: 2e-5,
        heat_j: [0.03, 0.02],
        external_force_n: 30.0,
    }
}

fn heat_controls() -> DoubleHeatExchangeSolve {
    DoubleHeatExchangeSolve {
        faces: exchange().faces.map(|face| HeatExchangeSolve {
            heat_bracket_j: [-0.001, 0.001],
            energy_tolerance_j: 2e-6,
            ..face
        }),
        maximum_joint_trials: 128,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.15e} != {expected:.15e}, tolerance {tolerance:.3e}"
    );
}

fn gas_storage(
    initial: &DoubleActingGasPiston<IsochoricRealGasParcel>,
    next: &DoubleActingGasPiston<IsochoricRealGasParcel>,
) -> f64 {
    [
        (initial.primary().gas(), next.primary().gas()),
        (initial.opposing(), next.opposing()),
    ]
    .into_iter()
    .map(|(a, b)| {
        a.mass_kg()
            * a.model()
                .specific_internal_energy_change_between_densities(
                    a.density_kg_m3(),
                    a.state().temperature_k(),
                    b.density_kg_m3(),
                    b.state().temperature_k(),
                )
                .unwrap()
    })
    .sum()
}

fn kinetic_change(
    initial: &DoubleActingGasPiston<IsochoricRealGasParcel>,
    next: &DoubleActingGasPiston<IsochoricRealGasParcel>,
) -> f64 {
    let (v0, v1) = (
        initial.primary().velocity_m_s(),
        next.primary().velocity_m_s(),
    );
    0.5 * initial.primary().mass_kg() * (v1 - v0) * (v1 + v0)
}

#[test]
fn asymmetric_real_gases_take_opposite_heat_from_one_sourced_body_and_close_energy() {
    let initial = real_pair();
    let saved = initial.clone();
    let (card, key, chart) = material(0.2);
    let context = point(300.0);
    let thermal = UniformEnthalpyStepInput {
        duration_s: input().duration_s,
        internal_heat_j: 0.015,
        ..body(chart.curve(), 500.0 * (325.0 - 200.0))
    };
    let candidate = initial
        .thermal_layers_candidate(
            input(),
            thermal,
            real_controls(),
            heat_controls(),
            [
                layer(&card, &key, &context, 0.012),
                layer(&card, &key, &context, 0.008),
            ],
            || false,
        )
        .unwrap();
    let report = candidate.step();
    let next = candidate.piston();
    assert!(report.exchanged_heat_j[0] > 0.0);
    assert!(report.exchanged_heat_j[1] < 0.0);
    assert_ne!(initial.primary().gas().model(), initial.opposing().model());
    assert_ne!(
        report.piston.midpoint[0].density_kg_m3(),
        report.piston.midpoint[1].density_kg_m3()
    );
    for (i, (old, gas)) in [
        (initial.primary().gas(), next.primary().gas()),
        (initial.opposing(), next.opposing()),
    ]
    .into_iter()
    .enumerate()
    {
        assert_eq!(gas.mass_kg(), old.mass_kg());
        assert_eq!(gas.model(), old.model());
        assert_eq!(
            gas.point().axes().get("specimen-lot"),
            old.point().axes().get("specimen-lot")
        );
        let midpoint = &report.piston.midpoint[i];
        near(
            midpoint.density_kg_m3(),
            old.mass_kg() / f64::midpoint(old.volume_m3(), gas.volume_m3()),
            2e-10,
        );
        assert_eq!(
            midpoint.transport_state().unwrap().density_kg_m3(),
            midpoint.acoustic_state().unwrap().density_kg_m3()
        );
        let source = layer(&card, &key, &context, [0.012, 0.008][i])
            .resolve(midpoint, report.thermal_midpoint)
            .unwrap();
        assert_eq!(candidate.thermal().report[i], source.source);
        assert_eq!(report.conductance_w_k[i], source.conductance_w_k);
        near(
            report.exchanged_heat_j[i],
            input().duration_s
                * source.conductance_w_k
                * (report.thermal_midpoint.temperature_k() - midpoint.temperature_k()),
            heat_controls().faces[i].source_tolerance_j,
        );
    }
    let stored_body = thermal.mass_kg
        * (candidate.thermal().state.specific_enthalpy_j_kg()
            - thermal.initial.specific_enthalpy_j_kg());
    let displacement = (next.primary().gas().volume_m3() - initial.primary().gas().volume_m3())
        / initial.primary().area_m2();
    near(
        gas_storage(&initial, next) + kinetic_change(&initial, next) + stored_body
            - input().heat_j.iter().sum::<f64>()
            - thermal.internal_heat_j
            - input().external_force_n * displacement,
        0.0,
        heat_controls().faces[0].energy_tolerance_j,
    );
    assert_eq!(
        candidate.thermal().external_heat_j,
        -report.exchanged_heat_j.iter().sum::<f64>()
    );
    assert!(report.thermal_energy_residual_j.abs() <= heat_controls().faces[0].energy_tolerance_j);
    assert!(report.total_energy_residual_j.abs() <= heat_controls().faces[0].energy_tolerance_j);
    assert!(report.joint_trials > 1 && report.joint_trials <= heat_controls().maximum_joint_trials);
    assert_eq!(initial, saved);
}

#[test]
fn continuing_real_pair_preserves_both_sources_and_one_complete_energy_account() {
    let initial = real_pair();
    let mut piston = initial.clone();
    let (card, key, chart) = material(0.2);
    let context = point(300.0);
    let drive = DoublePistonInput {
        heat_j: [0.4, 0.2],
        ..input()
    };
    let mut thermal = UniformEnthalpyStepInput {
        duration_s: drive.duration_s,
        internal_heat_j: 0.1,
        ..body(chart.curve(), 500.0 * (325.0 - 200.0))
    };
    let initial_enthalpy = thermal.initial.specific_enthalpy_j_kg();
    let mut supplied = 0.0;
    for _ in 0..3 {
        let candidate = piston
            .thermal_layers_candidate(
                drive,
                thermal,
                real_controls(),
                heat_controls(),
                [
                    layer(&card, &key, &context, 0.012),
                    layer(&card, &key, &context, 0.008),
                ],
                || false,
            )
            .unwrap();
        let (next, body, report) = candidate.into_parts();
        supplied += drive.heat_j.iter().sum::<f64>()
            + thermal.internal_heat_j
            + report.piston.external_force_work_j;
        for (a, b) in [
            (initial.primary().gas(), next.primary().gas()),
            (initial.opposing(), next.opposing()),
        ] {
            assert_eq!(a.mass_kg(), b.mass_kg());
            assert_eq!(a.model(), b.model());
        }
        assert!(report.exchanged_heat_j[0] > 0.0 && report.exchanged_heat_j[1] < 0.0);
        assert!(
            report.total_energy_residual_j.abs() <= heat_controls().faces[0].energy_tolerance_j
        );
        assert_eq!(next.last_window_duration_s(), Some(drive.duration_s));
        assert!(matches!(
            next.audit_last_window().unwrap(),
            ModuleWindowAudit::Available(_)
        ));
        piston = next;
        thermal.initial = body.state;
    }
    let stored_body =
        thermal.mass_kg * (thermal.initial.specific_enthalpy_j_kg() - initial_enthalpy);
    near(
        gas_storage(&initial, &piston) + kinetic_change(&initial, &piston) + stored_body,
        supplied,
        3.0 * heat_controls().faces[0].energy_tolerance_j,
    );
    assert!(piston.primary().gas().volume_m3() > initial.primary().gas().volume_m3());
    assert!(piston.opposing().volume_m3() < initial.opposing().volume_m3());
}

#[test]
fn real_pair_insulation_reproduces_the_original_motion_owner() {
    let piston = real_pair();
    let (card, key, chart) = material(0.2);
    let context = point(300.0);
    let thermal = UniformEnthalpyStepInput {
        duration_s: input().duration_s,
        ..body(chart.curve(), 500.0 * (325.0 - 200.0))
    };
    let mut expected = piston.clone();
    let native = expected.step(input(), real_controls(), || false).unwrap();
    let candidate = piston
        .thermal_layers_candidate(
            input(),
            thermal,
            real_controls(),
            heat_controls(),
            [
                layer(&card, &key, &context, 0.0),
                layer(&card, &key, &context, 0.0),
            ],
            || false,
        )
        .unwrap();
    assert_eq!(candidate.piston(), &expected);
    assert_eq!(
        candidate.piston().last_window_transfers(),
        expected.last_window_transfers()
    );
    assert_eq!(candidate.step().piston, native);
    assert_eq!(candidate.step().exchanged_heat_j, [0.0; 2]);
    assert_eq!(candidate.thermal().state, thermal.initial);
}

#[test]
fn opposing_source_refusal_and_final_cancel_preserve_both_gases_body_and_ports() {
    let mut piston = real_pair();
    piston.step(input(), real_controls(), || false).unwrap();
    let saved = piston.clone();
    let (card, key, chart) = material(0.2);
    let (wrong, _, _) = material(2.0);
    let context = point(300.0);
    let thermal = UniformEnthalpyStepInput {
        duration_s: input().duration_s,
        ..body(chart.curve(), 500.0 * (325.0 - 200.0))
    };
    let state = thermal.initial;
    assert!(matches!(
        piston.thermal_layers_candidate(
            input(),
            thermal,
            real_controls(),
            heat_controls(),
            [
                layer(&card, &key, &context, 0.012),
                layer(&wrong, &key, &context, 0.0)
            ],
            || false,
        ),
        Err(ThermalPistonError::Interface(_))
    ));
    assert_eq!(piston, saved);
    assert_eq!(thermal.initial, state);

    let mut polls = 0;
    piston
        .thermal_layers_candidate(
            input(),
            thermal,
            real_controls(),
            heat_controls(),
            [
                layer(&card, &key, &context, 0.012),
                layer(&card, &key, &context, 0.008),
            ],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 2);
    let mut count = 0;
    let cancelled = piston.thermal_layers_candidate(
        input(),
        thermal,
        real_controls(),
        heat_controls(),
        [
            layer(&card, &key, &context, 0.012),
            layer(&card, &key, &context, 0.008),
        ],
        || {
            count += 1;
            count == polls
        },
    );
    assert!(matches!(cancelled, Err(ThermalPistonError::Cancelled)));
    assert_eq!(piston, saved);
    assert_eq!(thermal.initial, state);
    assert_eq!(
        piston.last_window_duration_s(),
        saved.last_window_duration_s()
    );
    assert_eq!(
        piston.last_window_transfers(),
        saved.last_window_transfers()
    );
    assert_eq!(
        piston.audit_last_window().unwrap(),
        saved.audit_last_window().unwrap()
    );
}
