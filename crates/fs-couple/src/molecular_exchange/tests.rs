//! G0/G1/G3/G4: original finite gases, molecular moments and full U + K closure.
//! Synthetic NASA/transport sources are analytic fixtures, not measured gases.
mod trajectory;
use super::{
    MolecularExchangeCandidate, MolecularExchangeControls, MolecularExchangeError,
    MolecularExchangeSurface, step_candidate,
};
use crate::liquid_gas::evaporation::momentum::{MovingGasParcel, boundary::MovingGasBoundarySolve};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergySolve,
        parcel::{
            ClosedGasParcel,
            transfer::species::{
                boundary::GasSpeciesBoundarySolve,
                effusion::{
                    EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission,
                    EffusiveSupport, drift::EffusiveDriftInterface,
                    momentum::EffusiveMomentumInterface,
                },
            },
        },
    },
};
use fs_thermochem::{
    Composition, CompositionBasis, Nasa9RegionV1, Nasa9StandardStateModelV1, SpeciesId,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, mixture::FrozenIdealGasMixtureModelV1,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;

const PURE: &str = "A 1 0.028 2.5 0 1.7e-5 0.025";
const MIXTURE: &str = "A 0.5 0.028 2.5 0 1.7e-5 0.025\nB 0.5 0.004 2.5 0 2e-5 0.15";
const M: f64 = 0.028;
const T: f64 = 400.0;
const P: f64 = 100_000.0;
const VOLUME: f64 = 0.001;
const AREA: f64 = 1e-7;
const DT: f64 = 0.1;

fn catalog(text: &str, datum_j_kg: f64) -> ThermochemicalGasModel {
    let original = source::load(text, "finite molecular exchange analytic fixture").unwrap();
    if datum_j_kg == 0.0 {
        return original;
    }
    let species = original
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            let regions = species
                .regions()
                .iter()
                .map(|region| {
                    let mut card = region.card().clone();
                    card.parameters.get_mut("a7").unwrap().value +=
                        datum_j_kg * species.molar_mass().value() / R;
                    Nasa9RegionV1::from_card(card).unwrap()
                })
                .collect();
            Nasa9StandardStateModelV1::new(
                species.species().clone(),
                species.molar_mass(),
                species.convention().clone(),
                regions,
            )
            .unwrap()
        })
        .collect();
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            species,
            original.thermodynamics().mole_composition().clone(),
        )
        .unwrap(),
        original.transport_materials().to_vec(),
        original.water_species().cloned(),
    )
    .unwrap()
}

fn gas(model: ThermochemicalGasModel, t: f64, p: f64, velocity: [f64; 3]) -> MovingGasParcel {
    MovingGasParcel::new(
        ClosedGasParcel::new(model, source::point(t, p).unwrap(), VOLUME).unwrap(),
        velocity,
    )
    .unwrap()
}

fn interface(model: &ThermochemicalGasModel, maximum_speed: f64) -> EffusiveDriftInterface {
    let transmission = model
        .thermodynamics()
        .components()
        .iter()
        .map(|s| {
            EffusiveSpeciesTransmission::try_new(
                s.species().clone(),
                1.0,
                "synthetic reciprocal energy-independent molecular transmission",
            )
            .unwrap()
        })
        .collect();
    EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            EffusiveInterface::try_new(
                model.clone(),
                transmission,
                EffusiveSupport {
                    temperature_k: (200.0, 1800.0),
                    pressure_pa: (1000.0, 1e7),
                },
                EffusiveRegime::try_new(
                    1e-9,
                    0.0,
                    1e-7,
                    "synthetic collisionless thin-aperture qualification",
                )
                .unwrap(),
            )
            .unwrap(),
            "synthetic uncollimated nondeflecting half-space transmission",
        )
        .unwrap(),
        maximum_speed,
        "synthetic local Maxwellians over the authored relative-speed ball",
    )
    .unwrap()
}

fn surface() -> MolecularExchangeSurface {
    MolecularExchangeSurface {
        area_m2: AREA,
        normal_left_to_right: [1.0, 0.0, 0.0],
        velocity_world_m_s: [0.0; 3],
    }
}

fn controls() -> MolecularExchangeControls {
    let boundary = MovingGasBoundarySolve {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 1800.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-8,
        },
        momentum_absolute_n_s: 1e-12,
        momentum_relative_tolerance: 1e-12,
    };
    MolecularExchangeControls {
        gas: [boundary; 2],
        maximum_duration_s: 1.0,
        maximum_species_fraction: 0.1,
        mass_absolute_kg: 1e-13,
        momentum_absolute_n_s: 1e-12,
        energy_absolute_j: 5e-8,
        relative_tolerance: 1e-10,
    }
}

fn advance(
    gases: &[MovingGasParcel; 2],
    law: &EffusiveDriftInterface,
    aperture: MolecularExchangeSurface,
) -> MolecularExchangeCandidate {
    step_candidate(
        [&gases[0], &gases[1]],
        law,
        aperture,
        DT,
        controls(),
        || false,
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e} (allowance {tolerance:.3e})"
    );
}
fn vector_near(actual: [f64; 3], expected: [f64; 3], tolerance: f64) {
    for (a, b) in actual.into_iter().zip(expected) {
        near(a, b, tolerance);
    }
}
fn internal(gas: &MovingGasParcel) -> f64 {
    gas.parcel().mass_kg() * gas.parcel().state().specific_internal_energy_j_kg()
}
fn temperature(gas: &MovingGasParcel) -> f64 {
    gas.parcel().state().state().temperature
}
fn crossing_mass(t: f64, p: f64, molar_mass: f64) -> f64 {
    DT * AREA * p * (molar_mass / (2.0 * core::f64::consts::PI * R * t)).sqrt()
}
fn accounts(before: &[MovingGasParcel; 2], result: &MolecularExchangeCandidate) {
    let after = result.gases();
    near(
        after.iter().map(|s| s.parcel().mass_kg()).sum::<f64>(),
        before.iter().map(|s| s.parcel().mass_kg()).sum::<f64>(),
        2e-16,
    );
    for d in 0..3 {
        near(
            after.iter().map(|s| s.momentum_kg_m_s()[d]).sum::<f64>(),
            before.iter().map(|s| s.momentum_kg_m_s()[d]).sum::<f64>(),
            2e-12,
        );
    }
    let change: f64 = (0..2)
        .map(|i| {
            internal(&after[i]) - internal(&before[i]) + after[i].kinetic_energy_j()
                - before[i].kinetic_energy_j()
        })
        .sum();
    near(change, 0.0, 5e-8);
    near(
        result.step().mass_residual_kg,
        0.0,
        controls().mass_absolute_kg,
    );
    vector_near(
        result.step().momentum_residual_n_s,
        [0.0; 3],
        controls().momentum_absolute_n_s,
    );
    near(
        result.step().energy_residual_j,
        0.0,
        controls().energy_absolute_j,
    );
}

#[test]
fn unequal_finite_sources_transfer_original_mass_momentum_and_complete_energy() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = [
        gas(model.clone(), 600.0, 180_000.0, [0.0; 3]),
        gas(model, T, 90_000.0, [0.0; 3]),
    ];
    let saved = gases.clone();
    let result = advance(&gases, &law, surface());
    let mlr = crossing_mass(600.0, 180_000.0, M);
    let mrl = crossing_mass(T, 90_000.0, M);
    let impulse = DT * AREA * (180_000.0 + 90_000.0) / 2.0;
    let energy = mlr * 2.0 * R * 600.0 / M - mrl * 2.0 * R * T / M;
    let row = &result.step().exchange.species[0];
    near(row.mass_kg[0], mlr, 2e-18);
    near(row.mass_kg[1], mrl, 2e-18);
    near(row.net_mass_left_to_right_kg, mlr - mrl, 2e-18);
    for (side, sign) in [(0, -1.0), (1, 1.0)] {
        let after = &result.gases()[side];
        let mass = gases[side].parcel().mass_kg() + sign * (mlr - mrl);
        near(after.parcel().mass_kg(), mass, 2e-16);
        vector_near(after.momentum_kg_m_s(), [sign * impulse, 0.0, 0.0], 2e-15);
        let expected = (internal(&gases[side]) + sign * energy - impulse * impulse / (2.0 * mass))
            / (mass * 1.5 * R / M);
        near(temperature(after), expected, 2e-7);
        assert_eq!(after.parcel().volume_m3(), VOLUME);
    }
    accounts(&gases, &result);
    assert_eq!(gases, saved);
}

#[test]
fn zero_net_mass_counterflow_still_transfers_full_normal_stress() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = [
        gas(model.clone(), T, P, [0.0; 3]),
        gas(model, T, P, [0.0; 3]),
    ];
    let result = advance(&gases, &law, surface());
    let row = &result.step().exchange.species[0];
    assert!(row.mass_kg[0] > 0.0);
    assert_eq!(row.mass_kg[0], row.mass_kg[1]);
    assert_eq!(row.net_mass_left_to_right_kg, 0.0);
    // Both oppositely directed particle streams give the right gas +p A dt.
    // The aperture exchange is separate from pressure on the remaining walls.
    let impulse = P * AREA * DT;
    for (side, sign) in [(0, -1.0), (1, 1.0)] {
        let mass = gases[side].parcel().mass_kg();
        let after = &result.gases()[side];
        near(after.parcel().mass_kg(), mass, 2e-16);
        vector_near(after.momentum_kg_m_s(), [sign * impulse, 0.0, 0.0], 2e-15);
        let recoil_energy = impulse * impulse / (2.0 * mass);
        near(
            temperature(after),
            T - recoil_energy / (mass * 1.5 * R / M),
            2e-7,
        );
        assert!(temperature(after) < T);
    }
    accounts(&gases, &result);
}

#[test]
fn tangential_counterflow_thermalizes_its_actual_bulk_kinetic_loss_once() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let speed = 200.0;
    let gases = [
        gas(model.clone(), T, P, [0.0, speed, 0.0]),
        gas(model, T, P, [0.0, -speed, 0.0]),
    ];
    let result = advance(&gases, &law, surface());
    let dm = crossing_mass(T, P, M);
    let normal_impulse = P * AREA * DT;
    for (side, sign) in [(0, 1.0), (1, -1.0)] {
        let mass = gases[side].parcel().mass_kg();
        let py = sign * speed * (mass - 2.0 * dm);
        vector_near(
            result.gases()[side].momentum_kg_m_s(),
            [-sign * normal_impulse, py, 0.0],
            2e-15,
        );
        // Equal gross crossing energies cancel. Thermal storage is therefore
        // the actual bulk loss, including the small normal-recoil energy cost.
        let heat =
            0.5 * mass * speed * speed - (normal_impulse * normal_impulse + py * py) / (2.0 * mass);
        assert!(heat > 0.0);
        near(
            temperature(&result.gases()[side]),
            T + heat / (mass * 1.5 * R / M),
            2e-7,
        );
        near(result.step().exchange.species[0].mass_kg[side], dm, 2e-18);
    }
    accounts(&gases, &result);
}

#[test]
fn common_gas_and_aperture_boost_and_large_formation_datum_preserve_thermal_physics() {
    let base_velocity = [[0.0, 120.0, -30.0], [0.0, -40.0, 60.0]];
    let aperture_velocity = [0.0, 20.0, 10.0];
    let make_pair = |datum, boost: [f64; 3]| {
        let model = catalog(PURE, datum);
        [
            gas(
                model.clone(),
                500.0,
                160_000.0,
                core::array::from_fn(|d| base_velocity[0][d] + boost[d]),
            ),
            gas(
                model,
                700.0,
                90_000.0,
                core::array::from_fn(|d| base_velocity[1][d] + boost[d]),
            ),
        ]
    };
    let original = make_pair(0.0, [0.0; 3]);
    let law = interface(original[0].parcel().model(), 2000.0);
    let reference = advance(
        &original,
        &law,
        MolecularExchangeSurface {
            velocity_world_m_s: aperture_velocity,
            ..surface()
        },
    );
    accounts(&original, &reference);
    for datum in [0.0, 1e12] {
        for boost in [[0.0; 3], [8.0, -16.0, 32.0]] {
            let gases = make_pair(datum, boost);
            let law = interface(gases[0].parcel().model(), 2000.0);
            let result = advance(
                &gases,
                &law,
                MolecularExchangeSurface {
                    velocity_world_m_s: core::array::from_fn(|d| aperture_velocity[d] + boost[d]),
                    ..surface()
                },
            );
            for side in 0..2 {
                near(
                    temperature(&result.gases()[side]),
                    temperature(&reference.gases()[side]),
                    3e-7,
                );
                near(
                    result.gases()[side].parcel().mass_kg(),
                    reference.gases()[side].parcel().mass_kg(),
                    2e-16,
                );
                vector_near(
                    result.gases()[side].velocity_m_s(),
                    core::array::from_fn(|d| reference.gases()[side].velocity_m_s()[d] + boost[d]),
                    3e-11,
                );
                near(
                    result.step().exchange.species[0].mass_kg[side],
                    reference.step().exchange.species[0].mass_kg[side],
                    2e-18,
                );
            }
            // The datum is deliberately too large for a subtraction of
            // absolute U totals to be a useful test of the small thermal work.
            near(
                result.step().energy_residual_j,
                0.0,
                controls().energy_absolute_j,
            );
            vector_near(
                result.step().momentum_residual_n_s,
                [0.0; 3],
                controls().momentum_absolute_n_s,
            );
            if datum == 0.0 {
                accounts(&gases, &result);
            }
        }
    }
}

#[test]
fn gross_withdrawal_cannot_be_funded_by_counterflow_and_new_species_reach_both_gases() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = [
        gas(model.clone(), T, P, [0.0; 3]),
        gas(model, T, P, [0.0; 3]),
    ];
    let saved = gases.clone();
    let fraction = crossing_mass(T, P, M) / gases[0].parcel().mass_kg();
    assert!(matches!(
        step_candidate(
            [&gases[0], &gases[1]],
            &law,
            surface(),
            DT,
            MolecularExchangeControls {
                maximum_species_fraction: 0.5 * fraction,
                ..controls()
            },
            || false
        ),
        Err(MolecularExchangeError::Withdrawal { .. })
    ));
    assert!(
        step_candidate(
            [&gases[0], &gases[1]],
            &law,
            surface(),
            DT,
            MolecularExchangeControls {
                maximum_species_fraction: 2.0 * fraction,
                ..controls()
            },
            || false
        )
        .is_ok()
    );
    assert_eq!(gases, saved);

    let model = catalog(MIXTURE, 0.0);
    let law = interface(&model, 2000.0);
    let pure = |fractions| {
        model
            .with_composition(Composition::new(CompositionBasis::MoleFraction, fractions).unwrap())
            .unwrap()
    };
    let gases = [
        gas(pure(vec![1.0, 0.0]), T, P, [0.0; 3]),
        gas(pure(vec![0.0, 1.0]), T, P, [0.0; 3]),
    ];
    let a = SpeciesId::new("A").unwrap();
    let b = SpeciesId::new("B").unwrap();
    assert_eq!(gases[0].parcel().species_mass_kg(&b).unwrap(), 0.0);
    assert_eq!(gases[1].parcel().species_mass_kg(&a).unwrap(), 0.0);
    let result = advance(&gases, &law, surface());
    let ma = crossing_mass(T, P, M);
    let mb = crossing_mass(T, P, 0.004);
    for (side, species, expected) in [
        (0, &a, gases[0].parcel().mass_kg() - ma),
        (0, &b, mb),
        (1, &a, ma),
        (1, &b, gases[1].parcel().mass_kg() - mb),
    ] {
        near(
            result.gases()[side]
                .parcel()
                .species_mass_kg(species)
                .unwrap(),
            expected,
            2e-16,
        );
    }
    assert_eq!(result.step().exchange.species.len(), 2);
    assert!(
        result
            .gases()
            .iter()
            .all(|s| s.parcel().model().thermodynamics().components().len() == 2)
    );
    accounts(&gases, &result);
}

#[test]
fn source_endpoint_drift_and_terminal_cancellation_preserve_both_original_owners() {
    let model = catalog(PURE, 0.0);
    let law = interface(&model, 2000.0);
    let gases = [
        gas(model.clone(), T, P, [0.0; 3]),
        gas(model.clone(), T, P, [0.0; 3]),
    ];
    let saved = gases.clone();
    let bounded = interface(&model, 0.0);
    // Both initial drifts are zero; real transferred normal momentum puts
    // both endpoints outside a qualification that admits only zero drift.
    assert!(
        step_candidate(
            [&gases[0], &gases[1]],
            &bounded,
            surface(),
            DT,
            controls(),
            || false
        )
        .is_err()
    );
    let foreign = gas(catalog("A 1 0.028 3.5 0 1.7e-5 0.025", 0.0), T, P, [0.0; 3]);
    let foreign_saved = foreign.clone();
    assert!(
        step_candidate(
            [&gases[0], &foreign],
            &law,
            surface(),
            DT,
            controls(),
            || false
        )
        .is_err()
    );
    let closed = step_candidate(
        [&gases[0], &foreign],
        &bounded,
        MolecularExchangeSurface {
            area_m2: 0.0,
            velocity_world_m_s: [f64::NAN, 0.0, 0.0],
            ..surface()
        },
        DT,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(closed.gases(), &[gases[0].clone(), foreign.clone()]);

    let mut polls = 0;
    let accepted = step_candidate(
        [&gases[0], &gases[1]],
        &law,
        surface(),
        DT,
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        step_candidate(
            [&gases[0], &gases[1]],
            &law,
            surface(),
            DT,
            controls(),
            || {
                count += 1;
                count == polls
            }
        ),
        Err(MolecularExchangeError::Cancelled)
    ));
    assert_eq!(advance(&gases, &law, surface()).gases(), accepted.gases());
    assert_eq!(gases, saved);
    assert_eq!(foreign, foreign_saved);
}
