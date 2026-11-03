//! Source-bound variable-mass and inertial-frame energy checks.
#[path = "../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

use fs_couple::liquid_gas::evaporation::momentum::{
    MovingGasParcel, boundary::MovingGasBoundarySolve,
};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergySolve,
        parcel::{
            ClosedGasParcel,
            transfer::species::boundary::{
                GasSpeciesBoundarySolve, SpeciesBoundaryRelativeTransfer,
            },
        },
    },
};
use fs_thermochem::SpeciesId;

fn model() -> ThermochemicalGasModel {
    source::load(
        "A 0.75 0.028 3.5 0.001 2e-5 0.03\nB 0.25 0.004 2.5 0.0003 1.4e-5 0.021\n",
        "manufactured finite molecular energy",
    )
    .unwrap()
}
fn gas(model: &ThermochemicalGasModel, velocity: [f64; 3]) -> MovingGasParcel {
    MovingGasParcel::new(
        ClosedGasParcel::new(model.clone(), source::point(300.0, 100000.0).unwrap(), 0.01).unwrap(),
        velocity,
    )
    .unwrap()
}
fn controls() -> MovingGasBoundarySolve {
    MovingGasBoundarySolve {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 1000.0,
                energy_tolerance_j_kg: 1e-8,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        momentum_absolute_n_s: 1e-11,
        momentum_relative_tolerance: 1e-11,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}");
}

#[test]
fn species_and_impulse_share_one_original_internal_energy_update() {
    let model = model();
    let before = gas(&model, [20.0, -4.0, 0.0]);
    let b = SpeciesId::new("B").unwrap();
    for fraction in [-0.1, 0.1] {
        let dm = before.parcel().species_mass_kg(&b).unwrap() * fraction;
        let donor = [-10.0, 8.0, 0.0];
        let reference = before.velocity_m_s();
        let impulse = donor.map(|v| dm * v);
        let relative_speed_squared = (0..3)
            .map(|d| (donor[d] - reference[d]).powi(2))
            .sum::<f64>();
        let source_change = model
            .species_specific_enthalpy_change_j_kg(&b, 300.0, 350.0)
            .unwrap();
        let row = SpeciesBoundaryRelativeTransfer {
            species: b.clone(),
            mass_into_gas_kg: dm,
            relative_energy_into_gas_j: dm * (source_change + 0.5 * relative_speed_squared),
        };
        let result = before
            .exchange_species_momentum_relative(
                &model,
                &[row.clone()],
                impulse,
                reference,
                controls(),
                || false,
            )
            .unwrap();
        let after = result.gas();
        near(
            after.parcel().mass_kg() - before.parcel().mass_kg(),
            dm,
            1e-14,
        );
        for d in 0..3 {
            near(
                after.momentum_kg_m_s()[d] - before.momentum_kg_m_s()[d],
                impulse[d],
                1e-14,
            );
            near(result.step().momentum_residual_n_s[d], 0.0, 1e-14);
        }
        let energy = |s: &MovingGasParcel| {
            0.5 * s.parcel().mass_kg()
                * (0..3)
                    .map(|d| (s.velocity_m_s()[d] - reference[d]).powi(2))
                    .sum::<f64>()
        };
        near(
            result.step().kinetic_energy_change_j,
            energy(after) - energy(&before),
            1e-13,
        );
        near(
            result.step().caloric.additional_internal_energy_j,
            -result.step().kinetic_energy_change_j,
            0.0,
        );
        near(result.step().total_energy_residual_j, 0.0, 1e-8);
        assert_eq!(after.parcel().volume_m3(), before.parcel().volume_m3());
        assert_eq!(
            after.parcel().model().transport_materials(),
            before.parcel().model().transport_materials()
        );
    }
}

#[test]
fn translating_energy_frame_preserves_calorics_under_a_common_boost() {
    let model = model();
    let original = gas(&model, [20.0, -4.0, 0.0]);
    let id = SpeciesId::new("B").unwrap();
    let dm = 0.1 * original.parcel().species_mass_kg(&id).unwrap();
    let row = SpeciesBoundaryRelativeTransfer {
        species: id,
        mass_into_gas_kg: dm,
        relative_energy_into_gas_j: 0.2,
    };
    let base = original
        .exchange_species_momentum_relative(
            &model,
            &[row.clone()],
            [-10.0 * dm, 8.0 * dm, 0.0],
            [20.0, -4.0, 0.0],
            controls(),
            || false,
        )
        .unwrap();
    let boost = [10000.0, -20000.0, 0.0];
    let moved = gas(&model, [10020.0, -20004.0, 0.0]);
    let shifted = moved
        .exchange_species_momentum_relative(
            &model,
            &[row],
            [(-10.0 + boost[0]) * dm, (8.0 + boost[1]) * dm, 0.0],
            [10020.0, -20004.0, 0.0],
            controls(),
            || false,
        )
        .unwrap();
    near(
        base.gas().parcel().state().state().temperature,
        shifted.gas().parcel().state().state().temperature,
        1e-8,
    );
    near(
        base.step().kinetic_energy_change_j,
        shifted.step().kinetic_energy_change_j,
        1e-9,
    );
    assert_ne!(
        base.step().laboratory_kinetic_energy_change_j,
        shifted.step().laboratory_kinetic_energy_change_j
    );
}

#[test]
fn zero_exchange_is_exact_and_lost_impulse_caloric_failure_and_cancel_are_atomic() {
    let model = model();
    let before = gas(&model, [20.0, -4.0, 0.0]);
    let same = before
        .exchange_species_momentum_relative(&model, &[], [0.0; 3], [0.0; 3], controls(), || false)
        .unwrap();
    assert_eq!(same.gas(), &before);
    let huge = gas(&model, [1e16, 0.0, 0.0]);
    assert!(
        huge.exchange_species_momentum_relative(
            &model,
            &[],
            [1e-12, 0.0, 0.0],
            [1e16, 0.0, 0.0],
            controls(),
            || false
        )
        .is_err()
    );
    let saved = before.clone();
    let impossible = SpeciesBoundaryRelativeTransfer {
        species: SpeciesId::new("A").unwrap(),
        mass_into_gas_kg: 0.0,
        relative_energy_into_gas_j: 1e100,
    };
    assert!(
        before
            .exchange_species_momentum_relative(
                &model,
                &[impossible],
                [0.0; 3],
                [0.0; 3],
                controls(),
                || false
            )
            .is_err()
    );
    assert!(
        before
            .exchange_species_momentum_relative(&model, &[], [0.0; 3], [0.0; 3], controls(), || {
                true
            })
            .is_err()
    );
    assert_eq!(before, saved);
}
