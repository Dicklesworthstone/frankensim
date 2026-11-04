//! G0/G1/G3/G4: adaptive physical trajectories through the original gas owners.
//! Synthetic caloric/transport inputs verify equations, not measured accuracy.
use super::*;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::transfer::species::{
        effusion::{
            EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
        },
        network::{
            SpeciesNetworkError, SpeciesNetworkSolve,
            transport::effusion::{EffusionNetworkLink, reservoir},
        },
    },
    tests::{point, species_with_energy_reference, transport_in_domain},
};
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};

const M: f64 = 0.028;
const VOLUME: f64 = 0.001;
const AREA: f64 = 1e-7;

fn model(rows: &[(&str, f64, f64)], fractions: Vec<f64>, datum: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            rows.iter()
                .map(|&(id, m, cp)| species_with_energy_reference(id, m, cp, 0.0, datum))
                .collect(),
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        rows.iter()
            .map(|&(id, _, _)| transport_in_domain(id, 1.7e-5, 0.025, [200.0, 2000.0], [1.0, 1e7]))
            .collect(),
        None,
    )
    .unwrap()
}

fn pure(cp_r: f64, datum: f64, t: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        model(&[("A", M, cp_r)], vec![1.0], datum),
        point(t, p),
        VOLUME,
    )
    .unwrap()
}

fn interface(gas: &ThermochemicalGasModel, minimum_pressure: f64) -> EffusiveInterface {
    let rows = gas
        .thermodynamics()
        .components()
        .iter()
        .map(|s| {
            EffusiveSpeciesTransmission::try_new(
                s.species().clone(),
                1.0,
                "synthetic reciprocal energy-independent aperture",
            )
            .unwrap()
        })
        .collect();
    EffusiveInterface::try_new(
        gas.clone(),
        rows,
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (minimum_pressure, 1e7),
        },
        EffusiveRegime::try_new(
            1e-8,
            1e-11,
            1e-4,
            "synthetic thin-array collisionless support; no measured regime claim",
        )
        .unwrap(),
    )
    .unwrap()
}

fn endpoints(gases: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    gases
        .keys()
        .map(|&id| {
            (
                id,
                GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
            )
        })
        .collect()
}

fn physical() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}

fn adaptive(relative_tolerance: f64) -> EffusionAdvanceControls {
    EffusionAdvanceControls {
        initial_step_s: 5.0,
        minimum_step_s: 1e-7,
        maximum_step_s: 5.0,
        maximum_attempts: 12_000,
        maximum_accepted_intervals: 2_000,
        relative_tolerance,
        temperature_tolerance_k: 1e-7,
        pressure_tolerance_pa: 1e-7,
        species_mass_tolerance_kg: relative_tolerance * 1e-6,
        // This fixture isolates state-error refinement; caloric closure still
        // uses the independent 1e-8 J physical budget on every native step.
        stream_energy_tolerance_j: relative_tolerance,
    }
}

fn outlet(interface: &EffusiveInterface, area_m2: f64) -> VacuumEffusionLink<'_> {
    VacuumEffusionLink {
        id: 7,
        parcel: 0,
        interface,
        area_m2,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g0_window_energy_audit_rejects_changed_storage_with_valid_native_residuals() {
    let gases = BTreeMap::from([(0, pure(2.5, 1e12, 400.0, 200.0))]);
    let law = interface(gases[&0].model(), 1.0);
    let native = reservoir::step_candidate(
        &gases,
        GasTransportPaths::default(),
        &[outlet(&law, AREA)],
        &[],
        0.5,
        &endpoints(&gases),
        physical(),
        || false,
    )
    .unwrap();
    let record = || {
        let mut accepted = accounts::Accepted::default();
        accepted.record(&gases, &native, &mut || false).unwrap();
        accepted
    };
    let (_, actual) = record()
        .finish(&gases, native.parcels(), physical(), &mut || false)
        .unwrap();
    assert!(actual.accumulated_energy_residual_j.abs() <= physical().exchange.energy_tolerance_j);
    assert!(actual.energy_residual_j.abs() <= physical().exchange.energy_tolerance_j);

    // Change only the accepted endpoint's thermal energy. Scaling pressure with
    // temperature preserves density and all species masses under the same EOS.
    // Native receipts and their small original closure residuals are untouched.
    let gas = &native.parcels()[&0];
    let state = gas.state().state();
    let warmer = state.temperature + 0.25;
    let changed = ClosedGasParcel::new(
        gas.model().clone(),
        point(warmer, state.pressure * (warmer / state.temperature)),
        gas.volume_m3(),
    )
    .unwrap();
    assert!((changed.mass_kg() - gas.mass_kg()).abs() < physical().exchange.mass_tolerance_kg);
    let added_energy = gas.mass_kg() * 1.5 * R / M * 0.25;
    assert!(added_energy > 100.0 * physical().exchange.energy_tolerance_j);
    let altered = BTreeMap::from([(0, changed)]);
    assert!(matches!(
        record().finish(&gases, &altered, physical(), &mut || false),
        Err(Error::Invalid(
            "complete effusion window exceeded its fixed species/mass/energy budget"
        ))
    ));
}

#[test]
fn g1_adiabatic_blowdown_refines_to_monatomic_and_diatomic_solutions() {
    for cp_r in [2.5, 3.5] {
        let gases = BTreeMap::from([(0, pure(cp_r, 0.0, 400.0, 200.0))]);
        let original = gases.clone();
        let law = interface(gases[&0].model(), 1.0);
        let vents = [outlet(&law, AREA)];
        let duration = 20.0;
        let k0 = AREA / VOLUME * (R * 400.0 / (2.0 * core::f64::consts::PI * M)).sqrt();
        // d(n Cv T)/dn = (Cv+R/2)T, so T/T0=(n/n0)^(R/(2Cv)).
        // Integrating dn/dt=-k0*n*sqrt(T/T0) gives this independent oracle.
        let exponent = 4.0 * (cp_r - 1.0);
        let denominator = 1.0 + k0 * duration / exponent;
        let expected_mass = gases[&0].mass_kg() / denominator.powf(exponent);
        let expected_t = 400.0 / denominator.powi(2);
        let expected_p = 200.0 / denominator.powf(exponent + 2.0);
        let mut errors = Vec::new();
        let mut work = Vec::new();
        for tolerance in [2e-3, 2e-4, 2e-5] {
            let result = advance_candidate(
                &gases,
                GasTransportPaths::default(),
                &vents,
                &[],
                &[],
                duration,
                &endpoints(&gases),
                physical(),
                adaptive(tolerance),
                || false,
            )
            .unwrap();
            let gas = &result.parcels()[&0];
            errors.push(
                (gas.mass_kg() / expected_mass - 1.0)
                    .abs()
                    .max((gas.state().state().temperature / expected_t - 1.0).abs())
                    .max((gas.state().state().pressure / expected_p - 1.0).abs()),
            );
            let report = result.step();
            work.push(report.attempted_steps);
            assert_eq!(report.duration_s, duration);
            assert_eq!(report.accepted_steps.len(), 2 * report.accepted_intervals);
            assert!(report.maximum_accepted_error <= 1.0);
            near(
                report.accepted_steps.iter().map(|s| s.duration_s).sum(),
                duration,
                1e-12,
            );
            let escaped: f64 = report
                .accepted_steps
                .iter()
                .flat_map(|s| &s.exchange.boundaries)
                .map(|b| b.transfer.outgoing_mass_kg)
                .sum();
            near(gases[&0].mass_kg() - gas.mass_kg(), escaped, 1e-18);
            assert_eq!(gases, original);
        }
        assert!(errors[2] < 3e-3, "Cp/R={cp_r}: {errors:?}");
        for pair in errors.windows(2) {
            assert!(pair[1] < 0.65 * pair[0], "Cp/R={cp_r}: {errors:?}");
        }
        assert!(
            work[2] > work[0],
            "the temporal tolerance must change actual work: {work:?}"
        );
    }
}

#[test]
fn g3_gross_withdrawal_retries_are_deterministic_and_caloric_datum_invariant() {
    let mut endpoints_by_datum = Vec::new();
    for datum in [0.0, 1e12] {
        let gases = BTreeMap::from([(0, pure(2.5, datum, 400.0, 200.0))]);
        let original = gases.clone();
        let law = interface(gases[&0].model(), 1.0);
        let vents = [outlet(&law, 10.0 * AREA)];
        let mut controls = physical();
        controls.exchange.maximum_species_fraction = 0.01;
        let solves = endpoints(&gases);
        assert!(matches!(
            reservoir::step_candidate_with_drift(
                &gases,
                GasTransportPaths::default(),
                &vents,
                &[],
                &[],
                2.0,
                &solves,
                controls,
                || false,
            ),
            Err(GasTransportError::Exchange(
                SpeciesNetworkError::Withdrawal { .. }
            ))
        ));
        let run = || {
            advance_candidate(
                &gases,
                GasTransportPaths::default(),
                &vents,
                &[],
                &[],
                2.0,
                &solves,
                controls,
                adaptive(1e-4),
                || false,
            )
            .unwrap()
        };
        let result = run();
        let repeat = run();
        assert_eq!(result.parcels(), repeat.parcels());
        assert_eq!(result.step(), repeat.step());
        assert!(result.step().rejected_intervals > 0);
        assert!(result.step().accepted_intervals > 1);
        assert!(result.parcels()[&0].mass_kg() < gases[&0].mass_kg());
        assert!(result.parcels()[&0].state().state().temperature < 400.0);
        let gas = &result.parcels()[&0];
        endpoints_by_datum.push((
            gas.mass_kg(),
            gas.state().state().temperature,
            gas.state().state().pressure,
        ));
        assert_eq!(gases, original);
    }
    near(endpoints_by_datum[0].0, endpoints_by_datum[1].0, 2e-16);
    near(endpoints_by_datum[0].1, endpoints_by_datum[1].1, 2e-6);
    near(endpoints_by_datum[0].2, endpoints_by_datum[1].2, 2e-6);
}

#[test]
fn g3_new_species_crosses_finite_chambers_and_fractionates_through_vacuum() {
    let gases = BTreeMap::from([
        (0, pure(2.5, 0.0, 400.0, 200.0)),
        (1, pure(2.5, 0.0, 300.0, 100.0)),
    ]);
    let original = gases.clone();
    let source = ClosedGasParcel::new(
        model(&[("A", M, 2.5), ("B", 0.004, 2.5)], vec![0.5, 0.5], 0.0),
        point(500.0, 400.0),
        VOLUME,
    )
    .unwrap();
    let original_source = source.clone();
    let law = interface(source.model(), 1.0);
    let internal = [EffusionNetworkLink {
        id: 3,
        parcels: [0, 1],
        interface: &law,
        area_m2: AREA,
    }];
    let paths = GasTransportPaths {
        effusion: &internal,
        ..Default::default()
    };
    let vents = [VacuumEffusionLink {
        parcel: 1,
        ..outlet(&law, AREA)
    }];
    let supply = [EffusionReservoirLink {
        id: 9,
        parcel: 0,
        reservoir: &source,
        interface: &law,
        area_m2: AREA,
    }];
    let result = advance_candidate(
        &gases,
        paths,
        &vents,
        &supply,
        &[],
        10.0,
        &endpoints(&gases),
        physical(),
        adaptive(1e-4),
        || false,
    )
    .unwrap();
    let b = SpeciesId::new("B").unwrap();
    assert_eq!(gases[&0].species_mass_kg(&b).unwrap(), 0.0);
    assert!(result.parcels()[&0].species_mass_kg(&b).unwrap() > 0.0);
    assert!(result.parcels()[&1].species_mass_kg(&b).unwrap() > 0.0);
    for name in ["A", "B"] {
        let id = SpeciesId::new(name).unwrap();
        let total = |p: &BTreeMap<u64, ClosedGasParcel>| -> f64 {
            p.values().map(|g| g.species_mass_kg(&id).unwrap()).sum()
        };
        let supplied: f64 = result
            .step()
            .accepted_steps
            .iter()
            .flat_map(|s| &s.exchange.boundaries)
            .filter(|b| b.transfer.species == id)
            .map(|b| b.transfer.net_mass_into_gas_kg)
            .sum();
        near(total(result.parcels()) - total(&gases), supplied, 2e-17);
    }
    let first = &result.step().accepted_steps[0].reservoirs[0].flux.species;
    near(
        first[1].right_to_left.molar_rate_mol_s / first[0].right_to_left.molar_rate_mol_s,
        (M / 0.004).sqrt(),
        1e-12,
    );
    let last = result.step().accepted_steps.last().unwrap();
    let beam = &last.vacuum_outlets[0].flux.species;
    assert_eq!(beam.len(), 2);
    let flux_ratio = beam[1].outward.molar_rate_mol_s / beam[0].outward.molar_rate_mol_s;
    // Partial pressures are sampled from the same original finite chamber:
    // ratio of effused fluxes / ratio of parent mole fractions = sqrt(M_A/M_B).
    let original_moles = |name: &str, molar_mass: f64| {
        let id = SpeciesId::new(name).unwrap();
        let internal_delta: f64 = last
            .exchange
            .fluxes
            .iter()
            .filter(|f| f.transfer.species == id)
            .map(|f| {
                if f.transfer.receiver == 1 {
                    f.transfer.mass_kg
                } else if f.transfer.donor == 1 {
                    -f.transfer.mass_kg
                } else {
                    0.0
                }
            })
            .sum();
        let boundary_delta: f64 = last
            .exchange
            .boundaries
            .iter()
            .filter(|b| b.parcel == 1 && b.transfer.species == id)
            .map(|b| b.transfer.net_mass_into_gas_kg)
            .sum();
        (result.parcels()[&1].species_mass_kg(&id).unwrap() - internal_delta - boundary_delta)
            / molar_mass
    };
    near(
        flux_ratio,
        original_moles("B", 0.004) / original_moles("A", M) * (M / 0.004).sqrt(),
        1e-10,
    );
    let energy = |p: &BTreeMap<u64, ClosedGasParcel>| -> f64 {
        p.values()
            .map(|g| g.mass_kg() * g.state().specific_internal_energy_j_kg())
            .sum()
    };
    let energy_in: f64 = result
        .step()
        .accepted_steps
        .iter()
        .flat_map(|s| &s.exchange.boundaries)
        .map(|b| b.energy_into_gas_j)
        .sum();
    near(energy(result.parcels()) - energy(&gases), energy_in, 1e-8);
    assert_eq!(gases, original);
    assert_eq!(source, original_source);
}

#[test]
fn g4_source_refusal_work_limits_and_late_cancellation_publish_no_endpoint() {
    let gases = BTreeMap::from([(0, pure(2.5, 0.0, 400.0, 200.0))]);
    let original = gases.clone();
    let solves = endpoints(&gases);
    let narrow = interface(gases[&0].model(), 199.0);
    let narrow_vents = [outlet(&narrow, AREA)];
    // A sufficiently short physical step is supported, but a source-domain
    // failure must propagate instead of silently changing the requested model.
    reservoir::step_candidate_with_drift(
        &gases,
        GasTransportPaths::default(),
        &narrow_vents,
        &[],
        &[],
        0.01,
        &solves,
        physical(),
        || false,
    )
    .unwrap();
    let error = advance_candidate(
        &gases,
        GasTransportPaths::default(),
        &narrow_vents,
        &[],
        &[],
        1.0,
        &solves,
        physical(),
        EffusionAdvanceControls {
            maximum_attempts: 3,
            ..adaptive(1e-4)
        },
        || false,
    )
    .unwrap_err();
    match error {
        EffusionAdvanceError::Step {
            elapsed_s,
            duration_s,
            error: GasTransportError::Effusion(_),
        } => {
            assert_eq!(elapsed_s, 0.0);
            assert_eq!(duration_s, 1.0);
        }
        other => panic!("source refusal must not trigger temporal retries: {other}"),
    }
    let law = interface(gases[&0].model(), 1.0);
    let vents = [outlet(&law, AREA)];
    let mut polls = 0;
    let result = advance_candidate(
        &gases,
        GasTransportPaths::default(),
        &vents,
        &[],
        &[],
        20.0,
        &solves,
        physical(),
        adaptive(1e-4),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(result.step().accepted_intervals > 1);
    let cutoff = polls / 2;
    let mut seen = 0;
    let error = advance_candidate(
        &gases,
        GasTransportPaths::default(),
        &vents,
        &[],
        &[],
        20.0,
        &solves,
        physical(),
        adaptive(1e-4),
        || {
            seen += 1;
            seen >= cutoff
        },
    )
    .unwrap_err();
    assert_eq!(error, EffusionAdvanceError::Cancelled);
    assert!(matches!(
        advance_candidate(
            &gases,
            GasTransportPaths::default(),
            &vents,
            &[],
            &[],
            20.0,
            &solves,
            physical(),
            EffusionAdvanceControls {
                maximum_accepted_intervals: 1,
                ..adaptive(1e-4)
            },
            || false,
        ),
        Err(EffusionAdvanceError::WorkBudget {
            accepted_intervals: 1,
            ..
        })
    ));
    assert_eq!(gases, original);
}

#[test]
fn g1_drifting_reservoir_refines_with_shifted_datum_window_energy() {
    use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::{
        drift::EffusiveDriftInterface, momentum::EffusiveMomentumInterface,
    };
    let gases = BTreeMap::from([(0, pure(2.5, 1e12, 400.0, 200.0))]);
    let original = gases.clone();
    let source = pure(2.5, 1e12, 500.0, 300.0);
    let saved_source = source.clone();
    let law = EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            interface(source.model(), 1.0),
            "synthetic uncollimated nondeflecting angular source",
        )
        .unwrap(),
        1000.0,
        "synthetic collisionless moving incidence over original gas support",
    )
    .unwrap();
    let drifting = [EffusionDriftingReservoirLink {
        id: 11,
        parcel: 0,
        reservoir: &source,
        interface: &law,
        area_m2: AREA,
        relative_bulk_velocity_m_s: [-120.0, 180.0, 0.0],
        outward_unit_normal: [1.0, 0.0, 0.0],
    }];
    let solves = endpoints(&gases);
    let duration = 8.0;
    let native = |steps: u32| {
        let mut current = gases.clone();
        for _ in 0..steps {
            current = reservoir::step_candidate_with_drift(
                &current,
                GasTransportPaths::default(),
                &[],
                &[],
                &drifting,
                duration / f64::from(steps),
                &solves,
                physical(),
                || false,
            )
            .unwrap()
            .into_parts()
            .0;
        }
        current
    };
    let coordinates = |parcels: &BTreeMap<u64, ClosedGasParcel>| {
        let gas = &parcels[&0];
        [
            gas.mass_kg(),
            gas.state().state().temperature,
            gas.state().state().pressure,
        ]
    };
    let relative_error = |a: [f64; 3], b: [f64; 3]| {
        a.into_iter()
            .zip(b)
            .map(|(x, y)| (x / y - 1.0).abs())
            .fold(0.0, f64::max)
    };
    // This reference uses the SAME original drift/energy transaction; a
    // second refinement checks its numerical resolution for this comparison.
    let coarse_reference = coordinates(&native(512));
    let reference = coordinates(&native(1024));
    let reference_difference = relative_error(coarse_reference, reference);
    assert!(reference_difference < 2e-5, "{reference_difference}");
    let mut errors = Vec::new();
    for tolerance in [5e-4, 5e-5, 5e-6] {
        let result = advance_candidate(
            &gases,
            GasTransportPaths::default(),
            &[],
            &[],
            &drifting,
            duration,
            &solves,
            physical(),
            adaptive(tolerance),
            || false,
        )
        .unwrap();
        let report = result.step();
        assert!(report.accepted_intervals > 1);
        assert_eq!(report.boundaries.len(), 1);
        let boundary = &report.boundaries[0];
        assert!(boundary.incoming_mass_kg > 0.0 && boundary.outgoing_mass_kg > 0.0);
        near(
            result.parcels()[&0].mass_kg() - gases[&0].mass_kg(),
            boundary.net_mass_into_gas_kg,
            physical().exchange.mass_tolerance_kg,
        );
        assert!(report.balance.energy_residual_j.abs() <= physical().exchange.energy_tolerance_j);
        assert!(
            report.balance.accumulated_energy_residual_j.abs()
                <= physical().exchange.energy_tolerance_j
        );
        assert!(
            report
                .accepted_steps
                .iter()
                .all(|s| s.drifting_reservoirs.len() == 1)
        );
        errors.push(relative_error(coordinates(result.parcels()), reference));
    }
    assert!(errors[2] < 3e-4, "{errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.6 * pair[0], "{errors:?}");
    }
    assert_eq!(gases, original);
    assert_eq!(source, saved_source);
}
