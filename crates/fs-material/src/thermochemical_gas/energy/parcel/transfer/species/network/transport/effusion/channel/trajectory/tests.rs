//! Adaptive held-wall channels against an independent closed-form heating law.
//! Synthetic source curves test kinetics and accounts, not measured accuracy.
use super::super::super::{
    EffusionNetworkLink,
    tests::{controls, energy, gas, interface, near},
    trajectory::{EffusionAdvanceControls, EffusionAdvanceError},
};
use super::super::{DiffuseChannelNetworkLink, step_candidate as held_step};
use super::{advance_candidate, advance_candidate_with_row_limit};
use crate::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::ClosedGasParcel,
    parcel::transfer::species::{
        effusion::channel::{CylindricalMolecularGeometry, DiffuseMolecularChannel},
        network::{
            SpeciesNetworkError,
            transport::{GasTransportError, GasTransportPaths, GasTransportSolve},
        },
    },
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use std::collections::BTreeMap;

const M: f64 = 0.028;
const VOLUME: f64 = 0.001;
const AREA: f64 = 1e-6;
const T0: f64 = 400.0;
const WALL_T: f64 = 900.0;
const PRESSURE: f64 = 100_000.0;
const DURATION: f64 = 2.0;

fn gases(datum: f64) -> BTreeMap<u64, ClosedGasParcel> {
    BTreeMap::from([(0, gas(T0, PRESSURE, datum)), (1, gas(T0, PRESSURE, datum))])
}
fn channel(datum: f64) -> DiffuseMolecularChannel {
    DiffuseMolecularChannel::try_new(
        interface(datum),
        CylindricalMolecularGeometry::try_new(
            1e-8,
            1e-8,
            1e-4,
            "synthetic diffuse circular channels; no measured collision authority",
        )
        .unwrap(),
    )
    .unwrap()
}
fn link(channel: &DiffuseMolecularChannel) -> DiffuseChannelNetworkLink<'_> {
    DiffuseChannelNetworkLink {
        id: 7,
        parcels: [0, 1],
        channel,
        area_m2: AREA,
        wall_temperature_k: WALL_T,
    }
}
fn endpoints() -> BTreeMap<u64, GasEnergySolve> {
    let mut solve = energy();
    solve.energy_tolerance_j_kg = 1e-8;
    BTreeMap::from([(0, solve), (1, solve)])
}
fn physical() -> GasTransportSolve {
    let mut solve = controls();
    solve.exchange.energy_tolerance_j = 1e-7;
    solve
}
fn adaptive(relative_tolerance: f64) -> EffusionAdvanceControls {
    EffusionAdvanceControls {
        initial_step_s: DURATION,
        minimum_step_s: 1e-7,
        maximum_step_s: DURATION,
        maximum_attempts: 12_000,
        maximum_accepted_intervals: 2_000,
        relative_tolerance,
        temperature_tolerance_k: 1e-7,
        pressure_tolerance_pa: 1e-7,
        species_mass_tolerance_kg: relative_tolerance * 1e-5,
        stream_energy_tolerance_j: relative_tolerance * 10.0,
    }
}
fn temperature(gases: &BTreeMap<u64, ClosedGasParcel>) -> f64 {
    gases[&0].state().state().temperature
}
fn analytic_temperature(duration_s: f64) -> f64 {
    // Equal reservoirs have constant n. Each wall-hit monatomic molecule
    // exchanges 2R(Tw-T), while stored heat uses Cv=3R/2.
    // T'=(4/3)(A/V)(1-Fdirect)*sqrt(R/(2*pi*M))*sqrt(T)*(Tw-T).
    // The substitution y=sqrt(T) has the separable tanh solution below.
    // Equal-disk view factor for L/D=1 is independent of the channel adapter.
    let direct = (2.0_f64.sqrt() - 1.0).powi(2);
    let k = (2.0 / 3.0)
        * (AREA / VOLUME)
        * (1.0 - direct)
        * (R / (2.0 * core::f64::consts::PI * M)).sqrt();
    let r = (T0 / WALL_T).sqrt();
    let z = (k * WALL_T.sqrt() * duration_s).tanh();
    WALL_T * ((r + z) / (1.0 + r * z)).powi(2)
}

#[test]
fn symmetric_finite_gases_refine_to_the_exact_wall_heating_curve() {
    let gases = gases(0.0);
    let original = gases.clone();
    let law = channel(0.0);
    let expected = analytic_temperature(DURATION);
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for tolerance in [1e-2, 1e-5] {
        let result = advance_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law)],
            DURATION,
            &endpoints(),
            physical(),
            adaptive(tolerance),
            || false,
        )
        .unwrap();
        let t = temperature(result.parcels());
        let report = result.step();
        let balance = &report.balance;
        errors.push((t - expected).abs());
        work.push(report.attempted_steps);
        assert!(t > T0 && t < WALL_T);
        near(t, result.parcels()[&1].state().state().temperature, 1e-9);
        for id in [0, 1] {
            near(result.parcels()[&id].mass_kg(), gases[&id].mass_kg(), 1e-15);
            near(
                result.parcels()[&id].state().state().pressure,
                PRESSURE * t / T0,
                1e-5,
            );
        }
        let total_cv = 2.0 * gases[&0].mass_kg() / M * 1.5 * R;
        near(balance.heat_into_walls_j, -total_cv * (t - T0), 1e-7);
        near(
            balance.gas_entropy_change_j_k,
            total_cv * ((t - T0) / T0).ln_1p(),
            1e-9,
        );
        near(
            balance.wall_entropy_change_j_k,
            balance.heat_into_walls_j / WALL_T,
            1e-10,
        );
        assert!(balance.total_entropy_change_j_k > 0.0);
        assert!(balance.energy_residual_j.abs() <= physical().exchange.energy_tolerance_j);
        assert!(
            balance.accumulated_energy_residual_j.abs() <= physical().exchange.energy_tolerance_j
        );
        for residual in [
            balance.total_mass_residual_kg,
            balance.maximum_parcel_species_mass_residual_kg,
            balance.maximum_species_mass_residual_kg,
        ] {
            assert!(residual.abs() <= physical().exchange.mass_tolerance_kg);
        }
        assert_eq!(report.duration_s, DURATION);
        assert_eq!(report.accepted_steps.len(), 2 * report.accepted_intervals);
        assert_eq!(report.walls.len(), 1);
        assert!(
            report.rejected_intervals > 0,
            "the original 2 s withdrawal is too large"
        );
        assert!(report.maximum_accepted_error <= 1.0);
        near(
            report.accepted_steps.iter().map(|s| s.duration_s).sum(),
            DURATION,
            1e-12,
        );
        near(
            report
                .accepted_steps
                .iter()
                .map(|s| s.heat_into_walls_j)
                .sum(),
            balance.heat_into_walls_j,
            1e-9,
        );
    }
    assert!(errors[1] < 0.05, "fine heating error is {} K", errors[1]);
    assert!(
        errors[1] < 0.2 * errors[0],
        "tightening time accuracy must improve the physical curve: {errors:?}"
    );
    assert!(
        work[1] > work[0],
        "tighter accuracy must resolve more of the trajectory: {work:?}"
    );
    assert_eq!(gases, original);
}

#[test]
fn parallel_native_aperture_preserves_the_curve_and_shares_incident_inventory() {
    let gases = gases(0.0);
    let original = gases.clone();
    let law = channel(0.0);
    let channels = [link(&law)];
    let apertures = [EffusionNetworkLink {
        id: 21,
        parcels: [0, 1],
        interface: law.mouth(),
        area_m2: AREA,
    }];
    let paths = GasTransportPaths {
        effusion: &apertures,
        ..Default::default()
    };
    // Each path alone spends about 6.87% over 0.5 s. Together they exceed
    // the same original 10% allowance, even though their net mass is zero.
    held_step(
        &gases,
        GasTransportPaths::default(),
        &channels,
        0.5,
        &endpoints(),
        physical(),
        || false,
    )
    .unwrap();
    held_step(&gases, paths, &[], 0.5, &endpoints(), physical(), || false).unwrap();
    assert!(matches!(
        held_step(
            &gases,
            paths,
            &channels,
            0.5,
            &endpoints(),
            physical(),
            || false
        ),
        Err(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        ))
    ));
    let result = advance_candidate(
        &gases,
        paths,
        &channels,
        DURATION,
        &endpoints(),
        physical(),
        adaptive(1e-5),
        || false,
    )
    .unwrap();
    near(
        temperature(result.parcels()),
        analytic_temperature(DURATION),
        0.05,
    );
    assert!(
        result
            .step()
            .accepted_steps
            .iter()
            .all(|step| step.transport.effusion.len() == 1 && step.channels.len() == 1)
    );
    for id in [0, 1] {
        near(result.parcels()[&id].mass_kg(), gases[&id].mass_kg(), 1e-15);
    }
    assert_eq!(gases, original);
}

#[test]
fn huge_formation_datums_do_not_scale_physical_time_error_or_wall_accounts() {
    let mut baseline: Option<[f64; 7]> = None;
    for datum in [0.0, 1e24, -1e24] {
        let gases = gases(datum);
        let original = gases.clone();
        let law = channel(datum);
        let result = advance_candidate(
            &gases,
            GasTransportPaths::default(),
            &[link(&law)],
            0.5,
            &endpoints(),
            physical(),
            adaptive(1e-5),
            || false,
        )
        .unwrap();
        let balance = &result.step().balance;
        let actual = [
            temperature(result.parcels()),
            result.parcels()[&0].state().state().pressure,
            result.parcels()[&0].mass_kg(),
            balance.heat_into_walls_j,
            balance.gas_entropy_change_j_k,
            balance.wall_entropy_change_j_k,
            balance.total_entropy_change_j_k,
        ];
        if let Some(expected) = baseline {
            for ((x, y), tolerance) in actual
                .into_iter()
                .zip(expected)
                .zip([2e-6, 5e-4, 1e-15, 2e-6, 1e-8, 1e-8, 1e-8])
            {
                near(x, y, tolerance);
            }
        } else {
            baseline = Some(actual);
        }
        near(actual[0], analytic_temperature(0.5), 0.05);
        assert!(balance.energy_residual_j.abs() <= physical().exchange.energy_tolerance_j);
        assert!(
            balance.accumulated_energy_residual_j.abs() <= physical().exchange.energy_tolerance_j
        );
        assert_eq!(gases, original);
    }
}

#[test]
fn zero_area_never_queries_a_dormant_wall_temperature() {
    let gases = gases(0.0);
    let law = channel(0.0);
    let dormant = DiffuseChannelNetworkLink {
        area_m2: 0.0,
        wall_temperature_k: f64::NAN,
        ..link(&law)
    };
    let result = advance_candidate(
        &gases,
        GasTransportPaths::default(),
        &[dormant],
        DURATION,
        &endpoints(),
        physical(),
        adaptive(1e-5),
        || false,
    )
    .unwrap();
    assert_eq!(result.parcels(), &gases);
    assert_eq!(result.step().accepted_intervals, 1);
    assert_eq!(result.step().accepted_steps.len(), 2);
    assert_eq!(result.step().balance.heat_into_walls_j, 0.0);
    assert_eq!(result.step().balance.energy_residual_j, 0.0);
    assert_eq!(result.step().balance.gas_entropy_change_j_k, 0.0);
    assert_eq!(result.step().balance.wall_entropy_change_j_k, 0.0);
    assert_eq!(result.step().balance.total_entropy_change_j_k, 0.0);
}

#[test]
fn work_history_and_late_cancellation_refuse_the_whole_window_atomically() {
    let gases = gases(0.0);
    let original = gases.clone();
    let law = channel(0.0);
    let original_law = law.clone();
    let channels = [link(&law)];
    let solves = endpoints();
    for budget in [
        EffusionAdvanceControls {
            maximum_attempts: 3,
            ..adaptive(1e-4)
        },
        EffusionAdvanceControls {
            maximum_accepted_intervals: 1,
            ..adaptive(1e-4)
        },
    ] {
        assert!(matches!(
            advance_candidate(
                &gases,
                GasTransportPaths::default(),
                &channels,
                DURATION,
                &solves,
                physical(),
                budget,
                || false
            ),
            Err(EffusionAdvanceError::WorkBudget { .. })
        ));
    }
    assert!(matches!(
        advance_candidate_with_row_limit(
            &gases,
            GasTransportPaths::default(),
            &channels,
            DURATION,
            &solves,
            physical(),
            adaptive(1e-4),
            1,
            || false,
        ),
        Err(EffusionAdvanceError::WorkBudget { .. })
    ));
    let mut polls = 0;
    let completed = advance_candidate(
        &gases,
        GasTransportPaths::default(),
        &channels,
        DURATION,
        &solves,
        physical(),
        adaptive(1e-4),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(completed.step().accepted_intervals > 1);
    assert!(polls > 2);
    // Cancel at the last deterministic checkpoint, after native candidate work.
    let mut seen = 0;
    assert!(matches!(
        advance_candidate(
            &gases,
            GasTransportPaths::default(),
            &channels,
            DURATION,
            &solves,
            physical(),
            adaptive(1e-4),
            || {
                seen += 1;
                seen >= polls
            }
        ),
        Err(EffusionAdvanceError::Cancelled)
    ));
    assert_eq!(gases, original);
    assert_eq!(law, original_law);
}

#[test]
fn unresolved_heating_at_the_minimum_interval_is_an_explicit_accuracy_refusal() {
    let gases = gases(0.0);
    let original = gases.clone();
    let law = channel(0.0);
    let strict = EffusionAdvanceControls {
        initial_step_s: 0.5,
        minimum_step_s: 0.5,
        maximum_step_s: 0.5,
        ..adaptive(1e-12)
    };
    match advance_candidate(
        &gases,
        GasTransportPaths::default(),
        &[link(&law)],
        0.5,
        &endpoints(),
        physical(),
        strict,
        || false,
    ) {
        Err(EffusionAdvanceError::MinimumStep {
            elapsed_s,
            duration_s,
            normalized_error,
        }) => {
            assert_eq!(elapsed_s, 0.0);
            assert_eq!(duration_s, 0.5);
            assert!(normalized_error > 1.0);
        }
        other => panic!("the unresolved full/half heating discrepancy must refuse: {other:?}"),
    }
    assert_eq!(gases, original);
}

#[test]
fn opposite_half_step_heat_directions_cannot_cancel_the_accuracy_estimate() {
    use crate::thermochemical_gas::energy::parcel::transfer::species::network::{
        SpeciesNetworkExternalHeat, SpeciesNetworkHeatTransfer,
    };

    let gases = gases(0.0);
    let law = channel(0.0);
    let closed = DiffuseChannelNetworkLink {
        area_m2: 0.0,
        ..link(&law)
    };
    let coarse = held_step(
        &gases,
        GasTransportPaths::default(),
        &[closed],
        1.0,
        &endpoints(),
        physical(),
        || false,
    )
    .unwrap()
    .step()
    .clone();
    assert!(coarse.exchange.heat.is_empty());
    assert!(coarse.exchange.external_heat.is_empty());
    assert_eq!(coarse.channels[0].wall.heat_into_wall_j, 0.0);
    let accuracy = EffusionAdvanceControls {
        relative_tolerance: 0.0,
        stream_energy_tolerance_j: 0.25,
        ..adaptive(1e-5)
    };

    // Isolate each thermal namespace. The signed physical net remains zero,
    // while the accuracy estimate must retain both resolved heat directions.
    for mechanism in ["native", "wall", "accommodation"] {
        let mut first = coarse.clone();
        let mut second = coarse.clone();
        first.duration_s = 0.5;
        second.duration_s = 0.5;
        for (receipt, heat) in [(&mut first, 1.0), (&mut second, -1.0)] {
            match mechanism {
                "native" => receipt.exchange.heat.push(SpeciesNetworkHeatTransfer {
                    parcels: [0, 1],
                    heat_left_to_right_j: heat,
                }),
                "wall" => receipt.channels[0].wall.heat_into_wall_j = heat,
                "accommodation" => {
                    receipt
                        .exchange
                        .external_heat
                        .push(SpeciesNetworkExternalHeat {
                            id: 7,
                            parcel: 0,
                            heat_into_gas_j: heat,
                        });
                }
                _ => unreachable!(),
            }
        }
        let signed_net = match mechanism {
            "native" => {
                first.exchange.heat[0].heat_left_to_right_j
                    + second.exchange.heat[0].heat_left_to_right_j
            }
            "wall" => {
                first.channels[0].wall.heat_into_wall_j + second.channels[0].wall.heat_into_wall_j
            }
            "accommodation" => {
                first.exchange.external_heat[0].heat_into_gas_j
                    + second.exchange.external_heat[0].heat_into_gas_j
            }
            _ => unreachable!(),
        };
        assert_eq!(signed_net, 0.0);
        let error =
            super::accounts::stream_error(&coarse, &first, &second, accuracy, &mut || false)
                .unwrap();
        assert!(
            error > 1.0,
            "{mechanism} heat reversal must exceed the directional accuracy budget: {error}"
        );
    }
}
