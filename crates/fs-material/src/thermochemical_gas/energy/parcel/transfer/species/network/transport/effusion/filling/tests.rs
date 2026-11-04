//! First-fill equations use synthetic original sources, not measured apertures.
use super::*;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::transfer::species::{
        effusion::{EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport},
        network::{SpeciesNetworkError, SpeciesNetworkSolve, transport},
    },
    tests::{point, species_with_energy_reference, transport_in_domain},
};
use fs_matdb::QueryPoint;
use fs_thermochem::{
    Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    mixture::FrozenIdealGasMixtureModelV1,
};

const M: f64 = 0.028;
const VOLUME: f64 = 0.001;
pub(super) const AREA: f64 = 1e-7;
const DURATION: f64 = 0.25;

fn gas(rows: &[(&str, f64, f64, f64)], fractions: Vec<f64>, datum: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            rows.iter()
                .map(|&(id, mass, cp, slope)| {
                    species_with_energy_reference(id, mass, cp, slope, datum)
                })
                .collect(),
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        rows.iter()
            .map(|&(id, _, _, _)| {
                transport_in_domain(id, 1.7e-5, 0.025, [200.0, 2000.0], [1e-8, 1e7])
            })
            .collect(),
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, point(400.0, 200.0), VOLUME).unwrap()
}

pub(super) fn pure(cp: f64, slope: f64, datum: f64) -> ClosedGasParcel {
    gas(&[("A", M, cp, slope)], vec![1.0], datum)
}

pub(super) fn interface(
    gas: &ClosedGasParcel,
    transmission: f64,
    minimum_pressure: f64,
) -> EffusiveInterface {
    EffusiveInterface::try_new(
        gas.model().clone(),
        gas.model()
            .thermodynamics()
            .components()
            .iter()
            .map(|source| {
                EffusiveSpeciesTransmission::try_new(
                    source.species().clone(),
                    transmission,
                    "synthetic independent reciprocal energy-independent apertures",
                )
                .unwrap()
            })
            .collect(),
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (minimum_pressure, 1e7),
        },
        EffusiveRegime::try_new(
            1e-8,
            1e-11,
            1e-4,
            "synthetic collisionless thin-aperture support; no measured regime claim",
        )
        .unwrap(),
    )
    .unwrap()
}

pub(super) fn link(interface: &EffusiveInterface, area_m2: f64) -> EffusionFillLink<'_> {
    EffusionFillLink {
        id: 17,
        donor: 5,
        receiver: 9,
        interface,
        area_m2,
    }
}

pub(super) fn solve() -> GasEnergySolve {
    GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7,
        maximum_iterations: 128,
    }
}

pub(super) fn controls() -> GasTransportSolve {
    GasTransportSolve {
        exchange: SpeciesNetworkSolve {
            maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        entropy_tolerance_j_k: 1e-10,
    }
}

pub(super) fn empty(volume_m3: f64) -> EvacuatedGasVolume {
    EvacuatedGasVolume::new(volume_m3, QueryPoint::new().with("lot", 42.0).unwrap()).unwrap()
}

fn near(actual: f64, expected: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_captured_molecular_energy_sets_monatomic_diatomic_and_variable_cp_temperature() {
    for (cp, slope) in [(2.5, 0.0), (3.5, 0.0), (2.9, 0.001)] {
        let original = pure(cp, slope, 0.0);
        let receiver = empty(0.002);
        let law = interface(&original, 1.0, 1e-8);
        let candidate = step_candidate(
            &original,
            &receiver,
            link(&law, AREA),
            DURATION,
            solve(),
            solve(),
            controls(),
            || false,
        )
        .unwrap();
        let report = candidate.step();
        let filled = candidate.receiver().gas().unwrap();
        let captured =
            DURATION * AREA * 200.0 * (M / (2.0 * core::f64::consts::PI * R * 400.0)).sqrt();
        near(report.captured_mass_kg, captured, 1e-20);
        near(filled.mass_kg(), captured, 1e-20);
        // u/R = (Cp0/R-1)T + slope*T^2/2 + datum. A captured molecular
        // beam carries u(Td) + R*Td/2 per mole. Invert this independently.
        let cv = cp - 1.0;
        let target = (cv + 0.5) * 400.0 + 0.5 * slope * 400.0_f64.powi(2);
        let expected_t = if slope == 0.0 {
            target / cv
        } else {
            2.0 * target / (cv + (cv * cv + 2.0 * slope * target).sqrt())
        };
        near(filled.state().state().temperature, expected_t, 2e-8);
        near(
            filled.state().state().pressure,
            captured * R / M * expected_t / 0.002,
            1e-9,
        );
        assert!(filled.state().state().temperature > original.state().state().temperature);
        assert!(
            candidate.donor().state().state().temperature < original.state().state().temperature
        );
        assert_eq!(filled.point().axes().get("lot"), Some(&42.0));
        assert_eq!(filled.volume_m3(), receiver.volume_m3());
        assert_eq!(report.captured_species_kg.len(), 1);
        assert_eq!(
            report.fill.as_ref().unwrap().received_mass_kg,
            filled.mass_kg()
        );
        assert!(
            report.balance.total_mass_residual_kg.abs() <= controls().exchange.mass_tolerance_kg
        );
        assert!(
            report.balance.maximum_species_mass_residual_kg
                <= controls().exchange.mass_tolerance_kg
        );
        assert!(report.balance.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
        assert!(report.balance.entropy_change_j_k > 0.0);
        assert!(receiver.context().axes().get("temperature").is_none());
        assert!(receiver.context().axes().get("pressure").is_none());
    }
}

#[test]
fn g3_graham_fractionation_hands_both_real_parcels_to_reciprocal_transport() {
    let original = gas(
        &[("A", M, 3.5, 0.0), ("B", 0.004, 2.5, 0.0)],
        vec![0.75, 0.25],
        0.0,
    );
    let law = interface(&original, 1.0, 1e-8);
    let result = step_candidate(
        &original,
        &empty(0.002),
        link(&law, AREA),
        DURATION,
        solve(),
        solve(),
        controls(),
        || false,
    )
    .unwrap();
    let captured = result.receiver().gas().unwrap();
    let incoming_b = (0.25 / 0.004_f64.sqrt()) / (0.75 / M.sqrt() + 0.25 / 0.004_f64.sqrt());
    near(
        captured
            .model()
            .thermodynamics()
            .mole_composition()
            .fractions()[1],
        incoming_b,
        1e-14,
    );
    assert!(incoming_b > 0.25);
    assert!(
        result
            .donor()
            .model()
            .thermodynamics()
            .mole_composition()
            .fractions()[1]
            < 0.25
    );
    let cv = (1.0 - incoming_b) * 2.5 + incoming_b * 1.5;
    near(
        captured.state().state().temperature,
        400.0 * (cv + 0.5) / cv,
        2e-8,
    );

    let (donor, receiver, _) = result.into_parts();
    let EffusionReceiver::Filled(receiver) = receiver else {
        panic!("positive aperture must fill")
    };
    let gases = BTreeMap::from([(5, donor), (9, receiver)]);
    let paths = [super::super::EffusionNetworkLink {
        id: 17,
        parcels: [5, 9],
        interface: &law,
        area_m2: AREA,
    }];
    let next = transport::step_candidate_with_paths(
        &gases,
        GasTransportPaths {
            effusion: &paths,
            ..Default::default()
        },
        0.01,
        &BTreeMap::from([(5, solve()), (9, solve())]),
        controls(),
        || false,
    )
    .unwrap();
    assert!(next.step().effusion[0].flux.species.iter().all(|row| {
        row.left_to_right.mass_rate_kg_s > 0.0 && row.right_to_left.mass_rate_kg_s > 0.0
    }));
    let balance = audit_pair(
        &original,
        (5, &next.parcels()[&5]),
        Some((9, &next.parcels()[&9])),
        controls(),
        &mut || false,
    )
    .unwrap();
    assert!(balance.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    assert!(balance.entropy_change_j_k > 0.0);
}

#[test]
fn g3_first_fill_is_deterministic_and_invariant_to_large_caloric_reference() {
    let run = |datum| {
        let donor = pure(2.5, 0.0, datum);
        let law = interface(&donor, 1.0, 1e-8);
        step_candidate(
            &donor,
            &empty(0.002),
            link(&law, AREA),
            DURATION,
            solve(),
            solve(),
            controls(),
            || false,
        )
        .unwrap()
    };
    let original = run(0.0);
    assert_eq!(original, run(0.0));
    let shifted = run(1e12);
    for (a, b) in [
        (original.donor(), shifted.donor()),
        (
            original.receiver().gas().unwrap(),
            shifted.receiver().gas().unwrap(),
        ),
    ] {
        assert_eq!(a.mass_kg(), b.mass_kg());
        near(
            a.state().state().temperature,
            b.state().state().temperature,
            1e-9,
        );
        near(a.state().state().pressure, b.state().state().pressure, 1e-9);
    }
    assert_eq!(
        original.step().captured_species_kg,
        shifted.step().captured_species_kg
    );
    near(
        original.step().reference_relative_energy_j,
        shifted.step().reference_relative_energy_j,
        1e-14,
    );
    assert!(
        shifted.step().balance.energy_residual_j.abs() <= controls().exchange.energy_tolerance_j
    );
}

#[test]
fn g0_inactive_geometry_and_zero_transmission_leave_exact_vacuum() {
    let donor = pure(2.5, 0.0, 0.0);
    let receiver = empty(0.002);
    let foreign = gas(&[("other", 0.044, 3.5, 0.0)], vec![1.0], 0.0);
    let unrelated = interface(&foreign, 1.0, 1e-8);
    let blocked = interface(&donor, 0.0, 1e-8);
    for (law, area) in [(&unrelated, 0.0), (&blocked, AREA)] {
        let result = step_candidate(
            &donor,
            &receiver,
            link(law, area),
            DURATION,
            solve(),
            solve(),
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(result.donor(), &donor);
        assert_eq!(
            result.receiver(),
            &EffusionReceiver::Evacuated(receiver.clone())
        );
        assert!(result.step().fill.is_none());
        assert!(result.step().captured_species_kg.is_empty());
        assert_eq!(result.step().captured_mass_kg, 0.0);
        assert_eq!(result.step().reference_relative_energy_j, 0.0);
        assert_eq!(result.step().balance.energy_residual_j, 0.0);
        assert_eq!(result.step().balance.entropy_change_j_k, 0.0);
    }
}

#[test]
fn g3_small_receiver_rejects_entropy_losing_overfill_and_accepts_refinement() {
    let donor = pure(2.5, 0.0, 0.0);
    let receiver = empty(1e-9);
    let law = interface(&donor, 1.0, 1e-8);
    let error = step_candidate(
        &donor,
        &receiver,
        link(&law, AREA),
        0.1,
        solve(),
        solve(),
        controls(),
        || false,
    )
    .unwrap_err();
    assert!(matches!(error, Error::EntropyDecrease { change_j_k, .. } if change_j_k < 0.0));
    let refined = step_candidate(
        &donor,
        &receiver,
        link(&law, AREA),
        1e-5,
        solve(),
        solve(),
        controls(),
        || false,
    )
    .unwrap();
    assert!(refined.receiver().gas().is_some());
    assert!(refined.step().balance.entropy_change_j_k > 0.0);
}

#[test]
fn g4_first_fill_respects_real_source_domains_gross_inventory_and_atomic_cancellation() {
    let donor = pure(2.5, 0.0, 0.0);
    let saved_donor = donor.clone();
    let receiver = empty(0.002);
    let saved_receiver = receiver.clone();
    let law = interface(&donor, 1.0, 1e-8);
    assert!(matches!(
        step_candidate(
            &donor,
            &receiver,
            link(&law, AREA),
            50.0,
            solve(),
            solve(),
            controls(),
            || false
        ),
        Err(Error::Transport(GasTransportError::Exchange(
            SpeciesNetworkError::Withdrawal { .. }
        )))
    ));
    let narrow = interface(&donor, 1.0, 100.0);
    assert!(matches!(
        step_candidate(
            &donor,
            &receiver,
            link(&narrow, AREA),
            DURATION,
            solve(),
            solve(),
            controls(),
            || false
        ),
        Err(Error::Transport(GasTransportError::Effusion(_)))
    ));
    assert!(matches!(
        step_candidate(
            &donor,
            &receiver,
            link(&law, AREA),
            DURATION,
            solve(),
            GasEnergySolve {
                maximum_temperature_k: 500.0,
                ..solve()
            },
            controls(),
            || false
        ),
        Err(Error::Receiver(GasEnergyError::OutsideBracket { .. }))
    ));
    let mut polls = 0;
    step_candidate(
        &donor,
        &receiver,
        link(&law, AREA),
        DURATION,
        solve(),
        solve(),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    let result = step_candidate(
        &donor,
        &receiver,
        link(&law, AREA),
        DURATION,
        solve(),
        solve(),
        controls(),
        || {
            seen += 1;
            seen >= polls
        },
    );
    assert_eq!(result.unwrap_err(), Error::Cancelled);
    assert_eq!(donor, saved_donor);
    assert_eq!(receiver, saved_receiver);
}
