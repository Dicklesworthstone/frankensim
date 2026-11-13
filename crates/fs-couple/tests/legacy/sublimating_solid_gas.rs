//! Source-backed G0/G1/G3 checks for finite molecular phase/recoil exchange.
#[path = "sublimating_solid_gas/fixture.rs"]
mod fixture;
#[path = "sublimating_solid_gas/adaptive.rs"]
mod adaptive;
use fixture::*;
use fs_couple::solid_gas::kinetic::{
    KineticSublimationError, MovingSolidSurface, advance_solid_with_moving_gas,
};
use fs_material::sublimation::gas::{chamber::GasChamberSurface, momentum::MaxwellReflection};

#[test]
fn sublimation_recoils_both_finite_owners_and_closes_actual_total_energy() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = solid(400.0, caloric, [0.0; 3]);
    let reservoir = gas(&model, 400.0, 0.0, [0.0; 3]);
    let before = total_energy(&body, &reservoir, 0.0);
    let result = advance(
        &body,
        &reservoir,
        &model,
        &binding,
        caloric,
        None,
        DT,
        0.0,
        controls(),
    )
    .unwrap();
    assert!(result.step().outgoing_mass_kg > 0.0);
    assert_eq!(result.step().incoming_mass_kg, 0.0);
    assert!(result.solid().solid().mass_kg < body.solid().mass_kg);
    assert!(result.solid().velocity_m_s()[0] < 0.0);
    assert!(result.gas().velocity_m_s()[0] > 0.0);
    assert!(result.step().surface_temperature_k[1] < 400.0);
    assert!(result.gas().parcel().species_mass_kg(&species()).unwrap() > 0.0);
    near(
        result.solid().solid().mass_kg + result.gas().parcel().mass_kg(),
        body.solid().mass_kg + reservoir.parcel().mass_kg(),
        1e-13,
    );
    for d in 0..3 {
        near(
            result.solid().momentum_kg_m_s()[d] + result.gas().momentum_kg_m_s()[d],
            0.0,
            1e-12,
        );
    }
    near(
        total_energy(result.solid(), result.gas(), 0.0) - before,
        0.0,
        1e-7,
    );
    near(result.step().energy_residual_j, 0.0, 5e-8);
    near(result.step().laboratory_energy_residual_j, 0.0, 5e-8);
}

#[test]
fn deposition_grows_heats_and_preserves_carrier_inventory() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = solid(400.0, caloric, [0.0; 3]);
    let reservoir = gas(&model, 400.0, 5000.0, [0.0; 3]);
    let before = total_energy(&body, &reservoir, 0.0);
    let result = advance(
        &body,
        &reservoir,
        &model,
        &binding,
        caloric,
        None,
        DT,
        0.0,
        controls(),
    )
    .unwrap();
    assert!(result.step().incoming_mass_kg > result.step().outgoing_mass_kg);
    assert!(result.step().transferred_mass_kg < 0.0);
    assert!(result.solid().solid().mass_kg > MASS);
    assert!(result.step().surface_temperature_k[1] > 400.0);
    assert!(result.step().recession_m < 0.0);
    let carrier = fs_thermochem::SpeciesId::new("A").unwrap();
    near(
        result.gas().parcel().species_mass_kg(&carrier).unwrap(),
        reservoir.parcel().species_mass_kg(&carrier).unwrap(),
        1e-13,
    );
    near(
        total_energy(result.solid(), result.gas(), 0.0) - before,
        0.0,
        1e-7,
    );
}

#[test]
fn reflected_molecular_drag_transfers_momentum_and_heats_original_stores_without_phase_mass() {
    let model = catalog(0.0);
    let mut binding = binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let caloric = chart(&model, &binding);
    let body = solid(400.0, caloric, [0.0, 10.0, 0.0]);
    let reservoir = gas(&model, 400.0, 1000.0, [0.0, -5.0, 0.0]);
    let coefficients: Vec<_> = model
        .thermodynamics()
        .components()
        .iter()
        .map(|s| (s.species().clone(), 1.0))
        .collect();
    let reflection =
        MaxwellReflection::from_diffuse_fractions(&model, &coefficients, [200.0, 900.0]).unwrap();
    let before = total_energy(&body, &reservoir, 0.0);
    let result = advance(
        &body,
        &reservoir,
        &model,
        &binding,
        caloric,
        Some(&reflection),
        DT,
        0.0,
        controls(),
    )
    .unwrap();
    assert_eq!(result.step().transferred_mass_kg, 0.0);
    assert_eq!(result.solid().solid().mass_kg, body.solid().mass_kg);
    assert!(result.solid().velocity_m_s()[1] < 10.0);
    assert!(result.gas().velocity_m_s()[1] > -5.0);
    assert!(result.step().surface_temperature_k[1] > 400.0);
    near(
        total_energy(result.solid(), result.gas(), 0.0) - before,
        0.0,
        1e-7,
    );
}

#[test]
fn common_boost_and_consistent_formation_reference_keep_thermal_and_relative_dynamics() {
    let run = |reference: f64, boost: [f64; 3]| {
        let model = catalog(reference);
        let binding = binding(reference);
        let caloric = chart(&model, &binding);
        let body = solid(420.0, caloric, [boost[0] + 2.0, boost[1] + 7.0, boost[2]]);
        let reservoir = gas(
            &model,
            390.0,
            500.0,
            [boost[0] - 3.0, boost[1] - 4.0, boost[2]],
        );
        advance(
            &body,
            &reservoir,
            &model,
            &binding,
            caloric,
            None,
            DT,
            0.0,
            controls(),
        )
        .unwrap()
    };
    let original = run(0.0, [0.0; 3]);
    let boost = [1000.0, -2000.0, 3000.0];
    let changed = run(1e8, boost);
    near(
        changed.step().transferred_mass_kg,
        original.step().transferred_mass_kg,
        1e-17,
    );
    near(
        changed.step().surface_temperature_k[1],
        original.step().surface_temperature_k[1],
        2e-8,
    );
    near(
        changed.step().gas_temperature_k[1],
        original.step().gas_temperature_k[1],
        2e-8,
    );
    for d in 0..3 {
        near(
            changed.solid().velocity_m_s()[d] - boost[d],
            original.solid().velocity_m_s()[d],
            1e-10,
        );
        near(
            changed.gas().velocity_m_s()[d] - boost[d],
            original.gas().velocity_m_s()[d],
            1e-10,
        );
    }
    near(changed.step().laboratory_energy_residual_j, 0.0, 5e-8);
}

#[test]
fn zero_exposure_is_identity_and_external_heat_enters_only_once() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let original = solid(400.0, caloric, [2.0, -3.0, 4.0]);
    let mut thermal = original.solid();
    thermal.area_m2 = 0.0;
    let body = MovingSolidSurface::from_momentum(thermal, original.momentum_kg_m_s()).unwrap();
    let reservoir = gas(&model, 400.0, 0.0, [0.0; 3]);
    let idle = advance(
        &body,
        &reservoir,
        &model,
        &binding,
        caloric,
        None,
        DT,
        0.0,
        controls(),
    )
    .unwrap();
    assert_eq!(idle.solid(), &body);
    assert_eq!(idle.gas(), &reservoir);
    let heat = 0.1;
    let heated = advance(
        &body,
        &reservoir,
        &model,
        &binding,
        caloric,
        None,
        DT,
        heat,
        controls(),
    )
    .unwrap();
    assert_eq!(heated.gas(), &reservoir);
    near(
        heated.step().surface_temperature_k[1],
        400.0 + heat / MASS / caloric.1,
        1e-11,
    );
    near(
        total_energy(heated.solid(), heated.gas(), 0.0) - total_energy(&body, &reservoir, 0.0),
        heat,
        1e-8,
    );
}

#[test]
fn inventory_source_phase_and_late_cancellation_refuse_without_partial_endpoints() {
    let model = catalog(0.0);
    let binding = binding(0.0);
    let caloric = chart(&model, &binding);
    let body = solid(400.0, caloric, [0.0; 3]);
    let reservoir = gas(&model, 400.0, 500.0, [0.0; 3]);
    let original = (body, reservoir.clone());
    let mut small = controls();
    small.maximum_gross_inventory_fraction = 1e-12;
    assert!(
        matches!(advance(&body, &reservoir, &model, &binding, caloric,
        None, DT, 0.0, small), Err(KineticSublimationError::InventoryLimit { maximum_duration_s })
            if maximum_duration_s > 0.0 && maximum_duration_s < DT)
    );
    assert!(matches!(
        advance(
            &body,
            &reservoir,
            &model,
            &binding,
            caloric,
            None,
            1.0,
            0.0,
            controls()
        ),
        Err(KineticSublimationError::Duration { .. })
    ));
    assert!(
        advance(
            &body,
            &reservoir,
            &catalog(1e6),
            &binding,
            caloric,
            None,
            DT,
            0.0,
            controls()
        )
        .is_err()
    );
    let law = phase();
    let callback = |h| temperature(h, caloric);
    let source = GasChamberSurface {
        id: 7,
        law: &law,
        binding: &binding,
        solid: body.solid(),
        temperature_from_enthalpy: &callback,
    };
    assert!(
        advance_solid_with_moving_gas(
            &source,
            &body,
            &reservoir,
            &model,
            frame(),
            [caloric.0 + 1.0, caloric.0 + 2.0],
            None,
            DT,
            0.0,
            controls(),
            || false
        )
        .is_err()
    );
    let mut calls = 0;
    advance_solid_with_moving_gas(
        &source,
        &body,
        &reservoir,
        &model,
        frame(),
        enthalpy_range(caloric),
        None,
        DT,
        0.0,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut count = 0;
    assert!(matches!(
        advance_solid_with_moving_gas(
            &source,
            &body,
            &reservoir,
            &model,
            frame(),
            enthalpy_range(caloric),
            None,
            DT,
            0.0,
            controls(),
            || {
                count += 1;
                count == calls
            }
        ),
        Err(KineticSublimationError::Cancelled)
    ));
    assert_eq!((body, reservoir), original);
}

#[test]
fn specular_pressure_without_mass_transfer_closes_recoil_and_has_quadratic_source_truncation() {
    let model = catalog(0.0);
    let mut binding = binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let caloric = chart(&model, &binding);
    let body = solid(400.0, caloric, [0.0; 3]);
    let reservoir = gas(&model, 400.0, 1000.0, [0.0; 3]);
    let step = |dt| {
        advance(
            &body,
            &reservoir,
            &model,
            &binding,
            caloric,
            None,
            dt,
            0.0,
            controls(),
        )
        .unwrap()
    };
    let full = step(DT);
    let half = step(DT * 0.5);
    assert_eq!(full.step().transferred_mass_kg, 0.0);
    assert_eq!(full.step().incoming_mass_kg, 0.0);
    near(
        full.step().gas.requested_impulse_n_s[0],
        100_000.0 * AREA * DT,
        1e-14,
    );
    assert!(full.step().solid_kinetic_energy_change_j > 0.0);
    assert!(full.step().gas.kinetic_energy_change_j > 0.0);
    assert!(full.step().surface_temperature_k[1] < 400.0);
    assert!(full.step().gas_temperature_k[1] < 400.0);
    // Initially stationary specular walls perform zero instantaneous work.
    // The explicit source kick's O(dt^2) recoil energy is accounted for in
    // original thermal storage, never mislabeled as accommodation heating.
    near(
        full.step().solid_kinetic_energy_change_j,
        4.0 * half.step().solid_kinetic_energy_change_j,
        1e-15,
    );
    near(
        full.step().gas.kinetic_energy_change_j,
        4.0 * half.step().gas.kinetic_energy_change_j,
        1e-15,
    );
    near(
        full.step().solid_thermal_increment_j + full.step().solid_kinetic_energy_change_j,
        0.0,
        1e-9,
    );
    near(
        total_energy(full.solid(), full.gas(), 0.0) - total_energy(&body, &reservoir, 0.0),
        0.0,
        1e-7,
    );
}
