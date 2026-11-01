//! Different-mixture PR metering, original motion and finite-body heat together.

use super::*;
use crate::cavity_phs::piston::GasPiston;
use crate::cavity_phs::piston::double_acting::real_gas::{real_controls, real_pair};
use crate::cavity_phs::piston::material::{MaterialDoubleActingGasPiston, tests::body as specimen};
use crate::cavity_phs::piston::thermal::double_acting::DoubleHeatExchangeSolve;
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::{EquilibriumPhaseState, UniformEnthalpyStepInput};
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::nonideal::{
    ResolvedThermochemicalRealGasState, energy::parcel::IsochoricRealGasParcel,
    energy::parcel::transfer::RealGasMassTransferSolve,
    energy::parcel::transfer::nozzle::RealGasNozzleStepSolve,
};
use fs_ports::{ModuleWindowAudit, PortedModule};

fn flow_solve() -> RealGasNozzleStepSolve {
    RealGasNozzleStepSolve {
        entropy: GasEntropySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 450.0,
            entropy_tolerance_j_kg_k: 1e-9,
            maximum_iterations: 128,
        },
        transfer: RealGasMassTransferSolve {
            donor: real_controls().piston.gas,
            receiver: real_controls().opposing_gas,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-7,
        },
        // The expanded branch stays inside the original 200 K source floor.
        minimum_density_fraction: 0.5,
        pressure_tolerance_pa: 1e-4,
        sonic_tolerance_m2_s2: 1e-4,
        maximum_iterations: 128,
        maximum_donor_fraction: 0.05,
    }
}
fn geometry() -> RectangularMeteringValve {
    RectangularMeteringValve {
        reference_primary_volume_m3: 0.001,
        edge_height_at_reference_m: 0.001,
        width_m: 0.0002,
        window_height_m: 0.002,
        opens_on_primary_expansion: true,
        discharge_coefficient: 0.8,
    }
}
fn input() -> DoublePistonInput {
    DoublePistonInput {
        duration_s: 20e-6,
        heat_j: [0.02, -0.01],
        external_force_n: 30.0,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e} (allowance {tolerance:.3e})"
    );
}
fn unchanged(
    actual: &DoubleActingGasPiston<IsochoricRealGasParcel>,
    expected: &DoubleActingGasPiston<IsochoricRealGasParcel>,
) {
    assert_eq!(actual, expected);
    assert_eq!(actual.ports, expected.ports);
    assert_eq!(actual.primary().ports, expected.primary().ports);
}
fn species_masses(parcel: &IsochoricRealGasParcel) -> [f64; 2] {
    let x = parcel
        .model()
        .source()
        .thermodynamics()
        .mole_composition()
        .fractions();
    let molar_mass = x[0] * 0.028 + x[1] * 0.004;
    [
        parcel.mass_kg() * x[0] * 0.028 / molar_mass,
        parcel.mass_kg() * x[1] * 0.004 / molar_mass,
    ]
}
fn energy_change(
    old: &DoubleActingGasPiston<IsochoricRealGasParcel>,
    new: &DoubleActingGasPiston<IsochoricRealGasParcel>,
) -> f64 {
    // These fixtures have zero NASA formation constants. The native endpoint
    // U values retain PR departure and changed composition; no transfer or
    // piston report enters this independent extensive-storage reconstruction.
    let gas: f64 = [
        (old.primary().gas(), new.primary().gas()),
        (old.opposing(), new.opposing()),
    ]
    .into_iter()
    .map(|(a, b)| {
        b.mass_kg() * b.state().specific_internal_energy_j_kg()
            - a.mass_kg() * a.state().specific_internal_energy_j_kg()
    })
    .sum();
    let (v0, v1) = (old.primary().velocity_m_s(), new.primary().velocity_m_s());
    gas + 0.5 * old.primary().mass_kg() * (v1 - v0) * (v1 + v0)
}

#[test]
fn both_pressure_directions_replay_native_pr_nozzle_then_closed_motion() {
    for reverse in [false, true] {
        let original = real_pair();
        let original = if reverse {
            DoubleActingGasPiston::new(
                GasPiston::new(original.opposing().clone(), 2.0, 0.01, 0.4).unwrap(),
                original.primary().gas().clone(),
                0.008,
            )
            .unwrap()
        } else {
            original
        };
        let old = [original.primary().gas(), original.opposing()];
        let donor = usize::from(reverse);
        let opening = geometry().opening(&original).unwrap();
        let (transfer, flow) = old[donor]
            .through_nozzle_to(
                old[1 - donor],
                GasNozzle {
                    area_m2: opening.area_m2,
                    discharge_coefficient: geometry().discharge_coefficient,
                },
                input().duration_s,
                flow_solve(),
                || false,
            )
            .unwrap();
        let native_transfer = transfer.step().clone();
        let (donor_next, receiver_next, _) = transfer.into_parts();
        let mut replay = original.clone();
        if reverse {
            replay.primary.gas = receiver_next;
            replay.opposing = donor_next;
        } else {
            replay.primary.gas = donor_next;
            replay.opposing = receiver_next;
        }
        let native_motion = replay.step(input(), real_controls(), || false).unwrap();
        let candidate = original
            .valve_candidate(
                geometry(),
                input(),
                flow_solve(),
                real_controls(),
                2e-6,
                || false,
            )
            .unwrap();
        let next = candidate.piston();
        let report = candidate.step();
        assert_eq!(report.flow, flow);
        assert_eq!(report.transfer, native_transfer);
        assert_eq!(report.piston, native_motion);
        assert_eq!(next, &replay);
        assert_eq!(
            report.direction,
            Some(if reverse {
                ValveDirection::OpposingToPrimary
            } else {
                ValveDirection::PrimaryToOpposing
            })
        );
        assert!(report.transfer.transferred_mass_kg > 0.0);
        assert!(!report.flow.choked);
        near(
            report.flow.throat().pressure_pa(),
            old[1 - donor].state().pressure_pa(),
            flow_solve().pressure_tolerance_pa,
        );
        near(
            report.flow.mass_flow_kg_s,
            geometry().discharge_coefficient
                * opening.area_m2
                * report.flow.throat().density_kg_m3()
                * report.flow.velocity_m_s,
            1e-15,
        );
        near(
            report.transfer.transferred_mass_kg,
            input().duration_s * report.flow.mass_flow_kg_s,
            1e-18,
        );
        near(
            report.transfer.transferred_enthalpy_j,
            report.transfer.transferred_mass_kg * old[donor].state().specific_enthalpy_j_kg(),
            1e-12,
        );
        assert_ne!(
            report.flow.throat().density_kg_m3(),
            report.flow.throat().source().state().density
        );
        let new = [next.primary().gas(), next.opposing()];
        for species in 0..2 {
            let transported = report.transfer.transferred_mass_kg
                * species_masses(old[donor])[species]
                / old[donor].mass_kg();
            near(report.transfer.species[species].mass_kg, transported, 1e-14);
            near(
                species_masses(old[donor])[species] - species_masses(new[donor])[species],
                transported,
                1e-14,
            );
            near(
                species_masses(new[1 - donor])[species] - species_masses(old[1 - donor])[species],
                transported,
                1e-14,
            );
            near(
                species_masses(new[0])[species] + species_masses(new[1])[species],
                species_masses(old[0])[species] + species_masses(old[1])[species],
                1e-14,
            );
        }
        assert_eq!(new[donor].model(), old[donor].model());
        assert_ne!(
            new[1 - donor]
                .model()
                .source()
                .thermodynamics()
                .mole_composition(),
            old[1 - donor]
                .model()
                .source()
                .thermodynamics()
                .mole_composition()
        );
        for gas in new {
            assert_eq!(gas.density_kg_m3(), gas.mass_kg() / gas.volume_m3());
            assert_eq!(gas.point().axes().get("specimen-lot"), Some(&42.0));
        }
        near(
            energy_change(&original, next),
            input().heat_j.iter().sum::<f64>()
                + input().external_force_n * report.piston.displacement_m,
            2e-6,
        );
        near(report.total_mass_residual_kg, 0.0, 1e-14);
        assert!(report.total_energy_residual_j.abs() <= 2e-6);
        assert_eq!(next.ports.next_tick(), original.ports.next_tick() + 1);
        assert!(matches!(
            next.audit_last_window().unwrap(),
            ModuleWindowAudit::Unavailable { .. }
        ));
    }
}

#[test]
fn wall_position_opens_transport_and_closed_intervals_keep_native_clock_history() {
    let opening = RectangularMeteringValve {
        edge_height_at_reference_m: 0.0,
        ..geometry()
    };
    let mut live = real_pair();
    let initial = live.clone();
    let mut native = live.clone();
    let native_step = native.step(input(), real_controls(), || false).unwrap();
    let first = live
        .step_with_valve(
            opening,
            input(),
            flow_solve(),
            real_controls(),
            2e-6,
            || false,
        )
        .unwrap();
    unchanged(&live, &native);
    assert_eq!(first.piston, native_step);
    assert_eq!(first.direction, None);
    assert_eq!(first.transfer.transferred_mass_kg, 0.0);
    assert!(first.endpoint_opening.area_m2 > 0.0);
    let mut previous = first.endpoint_opening;
    for index in 1..3u64 {
        let report = live
            .step_with_valve(
                opening,
                input(),
                flow_solve(),
                real_controls(),
                2e-6,
                || false,
            )
            .unwrap();
        assert_eq!(report.opening, previous);
        assert_eq!(report.endpoint_opening, opening.opening(&live).unwrap());
        assert_eq!(report.direction, Some(ValveDirection::PrimaryToOpposing));
        assert!(report.transfer.transferred_mass_kg > 0.0);
        assert!(report.endpoint_opening.area_m2 > report.opening.area_m2);
        assert_eq!(live.ports.next_tick(), index + 2);
        assert_eq!(live.primary().ports.next_tick(), index + 2);
        assert!(matches!(
            live.audit_last_window().unwrap(),
            ModuleWindowAudit::Unavailable { .. }
        ));
        previous = report.endpoint_opening;
    }
    for species in 0..2 {
        near(
            species_masses(live.primary().gas())[species]
                + species_masses(live.opposing())[species],
            species_masses(initial.primary().gas())[species]
                + species_masses(initial.opposing())[species],
            1e-14,
        );
    }
    let closed = RectangularMeteringValve {
        discharge_coefficient: 0.0,
        ..opening
    };
    let mut native = live.clone();
    let native_step = native.step(input(), real_controls(), || false).unwrap();
    let report = live
        .step_with_valve(closed, input(), flow_solve(), real_controls(), 2e-6, || {
            false
        })
        .unwrap();
    unchanged(&live, &native);
    assert_eq!(report.piston, native_step);
    assert_eq!(report.direction, None);
    assert_eq!(report.transfer.transferred_mass_kg, 0.0);
    assert_eq!(live.ports.next_tick(), 5);
}

type HeatSource = (ResolvedThermochemicalRealGasState, EquilibriumPhaseState);
fn links(
    gases: &[ResolvedThermochemicalRealGasState; 2],
    body: EquilibriumPhaseState,
) -> Result<[ConductanceSample<HeatSource>; 2], fs_material::MaterialError> {
    let transport = [gases[0].transport_state()?, gases[1].transport_state()?];
    Ok(std::array::from_fn(|i| ConductanceSample {
        // Native gas-layer conduction: authored A/L in meters, original k.
        conductance_w_k: [10.0, 8.0][i] * transport[i].thermal_conductivity_w_m_k(),
        source: (gases[i].clone(), body),
    }))
}
fn heat_solve() -> DoubleHeatExchangeSolve {
    let mut controls = thermal::tests::exchange();
    for face in &mut controls.faces {
        face.heat_bracket_j = [-1e-4, 1e-4];
    }
    controls.maximum_joint_trials = 128;
    controls
}

#[test]
fn finite_body_heat_and_original_solid_density_share_the_metered_real_gas_trajectory() {
    let gas = real_pair();
    let chart = thermal::tests::curve(false);
    let forcing = DoublePistonInput {
        duration_s: 2e-6,
        ..input()
    };
    // This separate finite heat store has its own chart, mass and geometry.
    // The translating source-bound specimen retains its original rigid inertia.
    let body = UniformEnthalpyStepInput {
        initial: chart.state_at_specific_enthalpy(62500.0).unwrap(),
        duration_s: forcing.duration_s,
        internal_heat_j: 1e-5,
        ..thermal::tests::body(&chart, false)
    };
    let run = |density| {
        let solid = specimen(density);
        let original = MaterialDoubleActingGasPiston::new(
            gas.primary().gas().clone(),
            gas.opposing().clone(),
            solid.clone(),
            [0.01, 0.008],
            0.4,
        )
        .unwrap();
        let native = original
            .valve_candidate(
                geometry(),
                forcing,
                flow_solve(),
                real_controls(),
                2e-6,
                || false,
            )
            .unwrap();
        let mut polls = 0;
        let candidate = original
            .thermal_valve_candidate(
                geometry(),
                forcing,
                body,
                real_controls(),
                flow_solve(),
                heat_solve(),
                2e-6,
                links,
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        let report = candidate.step();
        assert_eq!(report.valve.flow, native.step().flow);
        assert_eq!(report.valve.transfer, native.step().transfer);
        assert_eq!(candidate.piston().body(), &solid);
        assert_eq!(
            candidate.piston().piston().primary().mass_kg(),
            solid.inertia().mass()
        );
        assert_eq!(candidate.piston().piston().ports.next_tick(), 2);
        assert!(report.exchanged_heat_j[0] > 0.0 && report.exchanged_heat_j[1] < 0.0);
        let samples = links(&report.valve.piston.midpoint, report.thermal_midpoint).unwrap();
        for i in 0..2 {
            assert_eq!(candidate.thermal().report[i], samples[i].source);
            near(report.conductance_w_k[i], samples[i].conductance_w_k, 1e-14);
            near(
                report.exchanged_heat_j[i],
                forcing.duration_s
                    * samples[i].conductance_w_k
                    * (report.thermal_midpoint.temperature_k()
                        - report.valve.piston.midpoint[i].temperature_k()),
                heat_solve().faces[i].source_tolerance_j,
            );
        }
        let stored = body.mass_kg
            * (candidate.thermal().state.specific_enthalpy_j_kg()
                - body.initial.specific_enthalpy_j_kg());
        near(
            stored + report.exchanged_heat_j.iter().sum::<f64>(),
            body.internal_heat_j,
            1e-9,
        );
        near(
            energy_change(original.piston(), candidate.piston().piston()) + stored,
            forcing.heat_j.iter().sum::<f64>()
                + body.internal_heat_j
                + forcing.external_force_n * report.valve.piston.displacement_m,
            2e-6,
        );
        assert!(report.total_energy_residual_j.abs() <= 2e-6);
        if density == 1000.0 {
            let saved = original.clone();
            let mut count = 0;
            assert!(matches!(
                original.thermal_valve_candidate(
                    geometry(),
                    forcing,
                    body,
                    real_controls(),
                    flow_solve(),
                    heat_solve(),
                    2e-6,
                    links,
                    || {
                        count += 1;
                        count == polls
                    },
                ),
                Err(ThermalPistonError::Cancelled)
            ));
            assert_eq!(count, polls);
            assert_eq!(original.body(), saved.body());
            unchanged(original.piston(), saved.piston());
        }
        let (next, _, _) = candidate.into_parts();
        next
    };
    let light = run(1000.0);
    let dense = run(8000.0);
    assert!(light.piston().primary().velocity_m_s() > dense.piston().primary().velocity_m_s());
    assert!(light.gases()[0].density_kg_m3() < dense.gases()[0].density_kg_m3());
    assert!(light.gases()[1].density_kg_m3() > dense.gases()[1].density_kg_m3());
    assert_eq!(light.gases()[0].mass_kg(), dense.gases()[0].mass_kg());
    assert_eq!(light.gases()[1].mass_kg(), dense.gases()[1].mass_kg());
}

#[test]
fn post_transfer_motion_refusal_and_final_step_cancel_publish_no_owner() {
    let mut live = real_pair();
    live.step_with_valve(
        geometry(),
        input(),
        flow_solve(),
        real_controls(),
        2e-6,
        || false,
    )
    .unwrap();
    let saved = live.clone();
    let mut invalid = real_controls();
    invalid.piston.volume_bracket_m3 = [0.0011, 0.0012];
    assert!(
        live.step_with_valve(geometry(), input(), flow_solve(), invalid, 2e-6, || false,)
            .is_err()
    );
    unchanged(&live, &saved);
    let mut probe = saved.clone();
    let mut polls = 0;
    let accepted = probe
        .step_with_valve(
            geometry(),
            input(),
            flow_solve(),
            real_controls(),
            2e-6,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(accepted.transfer.transferred_mass_kg > 0.0);
    let mut count = 0;
    assert_eq!(
        live.step_with_valve(
            geometry(),
            input(),
            flow_solve(),
            real_controls(),
            2e-6,
            || {
                count += 1;
                count == polls
            },
        ),
        Err(PistonError::Cancelled)
    );
    assert_eq!(count, polls);
    unchanged(&live, &saved);
    assert_eq!(
        geometry().opening(&live).unwrap(),
        geometry().opening(&saved).unwrap()
    );
}
