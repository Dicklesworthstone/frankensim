//! Synthetic original NASA/material sources exercise the aperture crosswalk.
use super::*;
use crate::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::parcel::transfer::species::effusion::{
        EffusiveInterface, EffusiveRegime, EffusiveSpeciesTransmission, EffusiveSupport,
    },
    tests::{point, species_with_energy_reference, transport},
};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const AREA: f64 = 1e-8;
const NORMAL: [f64; 3] = [1.0, 0.0, 0.0];

fn gas(datum: f64, temperature: f64) -> ClosedGasParcel {
    let mixture = FrozenIdealGasMixtureModelV1::new(
        vec![
            species_with_energy_reference("A", 0.004, 2.5, 0.0, datum),
            species_with_energy_reference("B", 0.028, 3.5, 0.001, datum),
        ],
        Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
    )
    .unwrap();
    ClosedGasParcel::new(
        ThermochemicalGasModel::new(
            mixture,
            vec![transport("A", 1.7e-5, 0.025), transport("B", 1.7e-5, 0.025)],
            None,
        )
        .unwrap(),
        point(temperature, 100_000.0),
        0.001,
    )
    .unwrap()
}
fn qualified(gas: &ClosedGasParcel, tau: [f64; 2]) -> EffusiveDriftInterface {
    let rows = gas
        .model()
        .thermodynamics()
        .components()
        .iter()
        .zip(tau)
        .map(|(source, probability)| {
            EffusiveSpeciesTransmission::try_new(
                source.species().clone(),
                probability,
                "synthetic energy/internal-state-independent transmission",
            )
            .unwrap()
        })
        .collect();
    let interface = EffusiveInterface::try_new(
        gas.model().clone(),
        rows,
        EffusiveSupport {
            temperature_k: (200.0, 1800.0),
            pressure_pa: (1000.0, 1e7),
        },
        EffusiveRegime::try_new(
            1e-9,
            1e-12,
            1e-7,
            "synthetic collisionless thin aperture across the original gas support",
        )
        .unwrap(),
    )
    .unwrap();
    EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            interface,
            "synthetic angle-independent nondeflecting transmission",
        )
        .unwrap(),
        100_000.0,
        "synthetic exposed local Maxwellian and collisionless aperture over the full speed/T/p support",
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-12 * actual.abs().max(expected.abs()).max(1e-250),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn vector_near(actual: [f64; 3], expected: [f64; 3]) {
    for (a, b) in actual.into_iter().zip(expected) {
        near(a, b);
    }
}
fn sample(
    gas: &ClosedGasParcel,
    interface: &EffusiveDriftInterface,
    velocity: [f64; 3],
    normal: [f64; 3],
) -> DriftingEffusiveFlux {
    gas.drifting_effusive_flux(interface, AREA, velocity, normal, || false)
        .unwrap()
}

#[test]
fn zero_drift_preserves_original_beams_and_half_pressure_momentum() {
    let gas = gas(0.0, 400.0);
    let interface = qualified(&gas, [0.8, 0.3]);
    let original = gas
        .effusive_flux_to_vacuum(interface.momentum_interface().interface(), AREA, || false)
        .unwrap();
    let normal = [0.6, 0.8, 0.0];
    let flux = sample(&gas, &interface, [0.0; 3], normal);
    assert_eq!(flux.interface, interface);
    for ((row, baseline), (mass, fraction, tau)) in flux
        .species
        .iter()
        .zip(original.species)
        .zip([(0.004, 0.25, 0.8), (0.028, 0.75, 0.3)])
    {
        assert_eq!(row.species, baseline.species);
        assert_eq!(row.baseline_incoming, baseline.outward);
        assert_eq!(row.incoming, baseline.outward);
        assert_eq!(row.incoming_mass_rate_change_from_zero_drift_kg_s, 0.0);
        let speed = (core::f64::consts::PI * R * 400.0 / (2.0 * mass)).sqrt();
        vector_near(
            row.mean_relative_velocity_m_s.unwrap(),
            normal.map(|n| -n * speed),
        );
        vector_near(
            row.impulse_for_mass_kg(row.incoming.mass_rate_kg_s)
                .unwrap(),
            normal.map(|n| -n * AREA * tau * fraction * 100_000.0 / 2.0),
        );
    }
}

#[test]
fn tangential_incidence_retains_mass_and_original_variable_heat_capacity_energy() {
    let gas = gas(0.0, 400.0);
    let interface = qualified(&gas, [0.8, 0.3]);
    let velocity = [0.0, 180.0, -75.0];
    let flux = sample(&gas, &interface, velocity, NORMAL);
    for (row, (mass, cp_r, slope)) in flux
        .species
        .iter()
        .zip([(0.004, 2.5, 0.0), (0.028, 3.5, 0.001)])
    {
        assert_eq!(
            row.incoming.mass_rate_kg_s,
            row.baseline_incoming.mass_rate_kg_s
        );
        assert_eq!(
            row.incoming.molar_rate_mol_s,
            row.baseline_incoming.molar_rate_mol_s
        );
        assert_eq!(row.incoming_mass_rate_change_from_zero_drift_kg_s, 0.0);
        let sigma_squared = R * 400.0 / mass;
        let correction = 0.5 * (180.0 * 180.0 + 75.0 * 75.0 - sigma_squared);
        let original_h = R * (cp_r * 400.0 + 0.5 * slope * 400.0 * 400.0) / mass;
        near(row.incoming.source_enthalpy_j_kg.unwrap(), original_h);
        near(
            row.incoming.additional_specific_energy_j_kg.unwrap(),
            correction,
        );
        near(
            row.incoming.additional_specific_energy_j_kg.unwrap(),
            row.baseline_incoming
                .additional_specific_energy_j_kg
                .unwrap()
                + 0.5 * (180.0 * 180.0 + 75.0 * 75.0),
        );
        near(
            row.incoming.crossing_energy_j_kg.unwrap(),
            original_h + correction,
        );
        near(
            row.incoming.energy_rate_w,
            row.incoming.mass_rate_kg_s * (original_h + correction),
        );
        vector_near(
            row.mean_relative_velocity_m_s.unwrap(),
            [
                -(core::f64::consts::PI * sigma_squared / 2.0).sqrt(),
                180.0,
                -75.0,
            ],
        );
        // The vector is weighted by accepted incoming mass, never by a speed
        // inferred from the complete internal/formation-energy coordinate.
        let mass = 2.37e-10;
        vector_near(
            row.impulse_for_mass_kg(mass).unwrap(),
            row.mean_relative_velocity_m_s.unwrap().map(|v| mass * v),
        );
    }
}

#[test]
fn opposite_half_spaces_recover_full_species_momentum_and_enthalpy_flux() {
    let gas = gas(1e5, 400.0);
    let interface = qualified(&gas, [0.8, 0.3]);
    for velocity in [[-600.0, 70.0, -20.0], [140.0, -60.0, 90.0]] {
        let incoming = sample(&gas, &interface, velocity, NORMAL);
        let outgoing = sample(&gas, &interface, velocity, [-1.0, 0.0, 0.0]);
        for ((i, o), (mass, fraction, tau)) in incoming
            .species
            .iter()
            .zip(&outgoing.species)
            .zip([(0.004, 0.25, 0.8), (0.028, 0.75, 0.3)])
        {
            let rho = fraction * 100_000.0 * mass / (R * 400.0);
            let net = AREA * tau * rho * velocity[0];
            near(o.incoming.mass_rate_kg_s - i.incoming.mass_rate_kg_s, net);
            let pi = i.impulse_for_mass_kg(i.incoming.mass_rate_kg_s).unwrap();
            let po = o.impulse_for_mass_kg(o.incoming.mass_rate_kg_s).unwrap();
            for d in 0..3 {
                let pressure = if d == 0 {
                    AREA * tau * fraction * 100_000.0
                } else {
                    0.0
                };
                near(po[d] - pi[d], net * velocity[d] + pressure);
            }
            let squared = velocity.into_iter().map(|v| v * v).sum::<f64>();
            near(
                o.incoming.energy_rate_w - i.incoming.energy_rate_w,
                net * (i.incoming.source_enthalpy_j_kg.unwrap() + 0.5 * squared),
            );
        }
    }
}

#[test]
fn direct_tail_and_tiny_signed_drift_survive_baseline_cancellation() {
    let gas = gas(0.0, 400.0);
    let interface = qualified(&gas, [0.0, 0.6]);
    let sigma = (R * 400.0 / 0.028).sqrt();
    let rho = 0.75 * 100_000.0 * 0.028 / (R * 400.0);
    let tail = sample(&gas, &interface, [10.0 * sigma, 0.0, 0.0], NORMAL);
    let row = &tail.species[1];
    // Independent positive integrals at inward drift a=-10. Alternating
    // Laplace expansions bound the last retained relative terms below 5e-17.
    let h1 = 7.474_560_254_589_328e-25;
    let h2_over_h1 = 0.194_383_033_412_553_3;
    let h3_over_h1 = 0.056_169_665_874_466_94;
    near(row.incoming.mass_rate_kg_s, AREA * 0.6 * rho * sigma * h1);
    near(
        row.mean_relative_velocity_m_s.unwrap()[0],
        -sigma * h2_over_h1,
    );
    near(
        row.incoming.additional_specific_energy_j_kg.unwrap(),
        0.5 * sigma * sigma * (h3_over_h1 - 3.0),
    );
    assert!(row.incoming.mass_rate_kg_s > 0.0);
    let reconstructed =
        row.baseline_incoming.mass_rate_kg_s + row.incoming_mass_rate_change_from_zero_drift_kg_s;
    assert!(
        (reconstructed - row.incoming.mass_rate_kg_s).abs() > 0.01 * row.incoming.mass_rate_kg_s
    );

    let headwind = sample(&gas, &interface, [-20.0 * sigma, 0.0, 0.0], NORMAL);
    let row = &headwind.species[1];
    near(row.incoming.mass_rate_kg_s, AREA * 0.6 * rho * 20.0 * sigma);
    near(
        row.mean_relative_velocity_m_s.unwrap()[0],
        -(20.0 + 1.0 / 20.0) * sigma,
    );
    near(
        row.incoming.additional_specific_energy_j_kg.unwrap(),
        0.5 * (20.0 * sigma).powi(2),
    );

    for sign in [-1.0, 1.0] {
        let speed = sign * sigma * 1e-20;
        let tiny = sample(&gas, &interface, [speed, 0.0, 0.0], NORMAL);
        let row = &tiny.species[1];
        near(
            row.incoming_mass_rate_change_from_zero_drift_kg_s,
            -0.5 * AREA * 0.6 * rho * speed,
        );
        assert!(
            row.incoming_mass_rate_change_from_zero_drift_kg_s.abs()
                < f64::EPSILON * row.baseline_incoming.mass_rate_kg_s
        );
    }
    assert!(
        gas.drifting_effusive_flux(&interface, AREA, [50.0 * sigma, 0.0, 0.0], NORMAL, || false)
            .is_err()
    );
}

#[test]
fn source_and_drift_qualification_distinguish_closed_area_from_zero_transmission() {
    let gas = gas(0.0, 400.0);
    let before = gas.clone();
    let interface = qualified(&gas, [0.0, 0.0]);
    let closed = sample(&gas, &interface, [100.0, 50.0, 0.0], NORMAL);
    for row in closed.species {
        assert_eq!(row.incoming, EffusiveBeam::default());
        assert_eq!(row.mean_relative_velocity_m_s, None);
        assert_eq!(row.impulse_for_mass_kg(0.0).unwrap(), [0.0; 3]);
        assert!(row.impulse_for_mass_kg(1e-12).is_err());
    }
    let foreign = self::gas(1.0, 400.0);
    assert!(
        foreign
            .drifting_effusive_flux(&interface, AREA, [0.0; 3], NORMAL, || false)
            .is_err()
    );
    assert!(
        self::gas(0.0, 1900.0)
            .drifting_effusive_flux(&interface, AREA, [0.0; 3], NORMAL, || false)
            .is_err()
    );
    assert!(
        gas.drifting_effusive_flux(&interface, AREA, [200_000.0, 0.0, 0.0], NORMAL, || false)
            .is_err()
    );
    assert!(
        gas.drifting_effusive_flux(&interface, AREA, [0.0; 3], [2.0, 0.0, 0.0], || false)
            .is_err()
    );
    let inactive = foreign
        .drifting_effusive_flux(&interface, 0.0, [200_000.0, 0.0, 0.0], [0.0; 3], || false)
        .unwrap();
    assert!(inactive.species.is_empty());
    for speed in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            EffusiveDriftInterface::try_new(
                interface.momentum_interface().clone(),
                speed,
                "synthetic drift qualification"
            )
            .is_err()
        );
    }
    assert!(
        EffusiveDriftInterface::try_new(interface.momentum_interface().clone(), 1.0, "").is_err()
    );
    for limit in [0.0, f64::from_bits(1), 1e300] {
        let bound = EffusiveDriftInterface::try_new(
            interface.momentum_interface().clone(),
            limit,
            "synthetic explicit speed support",
        )
        .unwrap();
        bound.admit_relative_velocity([limit, 0.0, 0.0]).unwrap();
        if limit == 0.0 {
            assert!(
                bound
                    .admit_relative_velocity([f64::from_bits(1), 0.0, 0.0])
                    .is_err()
            );
        } else {
            assert!(bound.admit_relative_velocity([limit, limit, 0.0]).is_err());
        }
    }
    let active = qualified(&gas, [0.8, 0.3]);
    let mut count = 0;
    gas.drifting_effusive_flux(&active, AREA, [100.0, 50.0, 0.0], NORMAL, || {
        count += 1;
        false
    })
    .unwrap();
    let mut calls = 0;
    assert!(matches!(
        gas.drifting_effusive_flux(&active, AREA, [100.0, 50.0, 0.0], NORMAL, || {
            calls += 1;
            calls >= count
        }),
        Err(GasEnergyError::Cancelled)
    ));
    assert_eq!(gas, before);
}

#[test]
fn rotations_area_and_caloric_datums_preserve_original_conditional_moments() {
    let gas = gas(0.0, 400.0);
    let interface = qualified(&gas, [0.8, 0.3]);
    let velocity = [125.0, -90.0, 80.0];
    let flux = sample(&gas, &interface, velocity, NORMAL);
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let rotated = sample(&gas, &interface, rotate(velocity), rotate(NORMAL));
    let double = gas
        .drifting_effusive_flux(&interface, 2.0 * AREA, velocity, NORMAL, || false)
        .unwrap();
    let shifted_gas = self::gas(1e24, 400.0);
    let shifted = sample(
        &shifted_gas,
        &qualified(&shifted_gas, [0.8, 0.3]),
        velocity,
        NORMAL,
    );
    for (((row, rotated), double), shifted) in flux
        .species
        .iter()
        .zip(&rotated.species)
        .zip(&double.species)
        .zip(&shifted.species)
    {
        vector_near(
            rotated.mean_relative_velocity_m_s.unwrap(),
            rotate(row.mean_relative_velocity_m_s.unwrap()),
        );
        near(rotated.incoming.mass_rate_kg_s, row.incoming.mass_rate_kg_s);
        near(
            rotated.incoming.additional_specific_energy_j_kg.unwrap(),
            row.incoming.additional_specific_energy_j_kg.unwrap(),
        );
        near(
            double.incoming.mass_rate_kg_s,
            2.0 * row.incoming.mass_rate_kg_s,
        );
        near(
            double.incoming.energy_rate_w,
            2.0 * row.incoming.energy_rate_w,
        );
        assert_eq!(shifted.incoming.mass_rate_kg_s, row.incoming.mass_rate_kg_s);
        assert_eq!(
            shifted.mean_relative_velocity_m_s,
            row.mean_relative_velocity_m_s
        );
        assert_eq!(
            shifted.incoming.additional_specific_energy_j_kg,
            row.incoming.additional_specific_energy_j_kg
        );
        near(
            shifted.incoming.source_enthalpy_j_kg.unwrap(),
            row.incoming.source_enthalpy_j_kg.unwrap() + 1e24,
        );
        near(
            shifted.incoming.energy_rate_w,
            row.incoming.energy_rate_w + row.incoming.mass_rate_kg_s * 1e24,
        );
        assert!(row.impulse_for_mass_kg(-1.0).is_err());
        assert!(row.impulse_for_mass_kg(f64::NAN).is_err());
    }
}
