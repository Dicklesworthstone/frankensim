//! G0/G1/G3: actual osmotic particles, finite-rate lag and transactional updates.

use super::*;
use fs_qty::{Amount, AmountConcentration};
use fs_thermochem::solution::solubility::IdealSolubilityV1;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn law(name: &str, concentration: f64, enthalpy: f64) -> IdealSolubilityV1 {
    IdealSolubilityV1::new(
        SpeciesId::new("water").unwrap(),
        SpeciesId::new(name).unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(concentration),
        enthalpy,
        [Temperature::new(280.0), Temperature::new(350.0)],
    )
    .unwrap()
}
fn component(
    name: &str,
    concentration: f64,
    total: f64,
    dissolved: f64,
    conductance: f64,
) -> (SolubilityKineticsV1, SolubilityKineticStateV1) {
    (
        SolubilityKineticsV1::new(
            law(name, concentration, 0.0),
            VolumetricFlowRate::new(conductance),
            PrecipitationPolicyV1::SeededSurface,
        )
        .unwrap(),
        SolubilityKineticStateV1::new(Amount::new(total), Amount::new(dissolved)).unwrap(),
    )
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-13 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}
fn reservoir() -> KineticSolubilitySolution {
    KineticSolubilitySolution::new(
        Temperature::new(300.0),
        18e-6,
        Volume::new(1e-3),
        vec![
            component("A", 100.0, 0.2, 0.0, 1e-3),
            component("B", 50.0, 0.1, 0.1, 5e-4),
        ],
    )
    .unwrap()
}

#[test]
fn actual_dissolved_amounts_drive_osmotic_pressure_during_opposed_phase_transfers() {
    let mut owner = reservoir();
    let phase_buffer = owner.last_steps().as_ptr();
    let dissolved_buffer = owner.solution().solute_moles().as_ptr();
    owner.advance(Time::new(1.0), || false).unwrap();
    let expected = [
        0.1 * (1.0 - (-1.0_f64).exp()),
        0.05 + 0.05 * (-0.5_f64).exp(),
    ];
    for (index, value) in expected.into_iter().enumerate() {
        near(owner.solution().solute_moles()[index], value);
        let state = owner.last_steps()[index].state();
        near(
            state.dissolved().value() + state.solid().value(),
            state.total_solute().value(),
        );
    }
    near(
        owner.solution().osmotic_pressure_pa().unwrap(),
        R * 300.0 * expected.iter().sum::<f64>() / 1e-3,
    );
    near(
        owner.solution().water_activity().unwrap(),
        (-18e-6 * expected.iter().sum::<f64>() / 1e-3).exp(),
    );
    assert!(owner.last_steps()[0].transferred_to_solution().value() > 0.0);
    assert!(owner.last_steps()[1].transferred_to_solution().value() < 0.0);
    assert_eq!(owner.last_steps().as_ptr(), phase_buffer);
    assert_eq!(owner.solution().solute_moles().as_ptr(), dissolved_buffer);
}

#[test]
fn dilution_and_prescribed_heating_preserve_stocks_then_resolve_the_kinetic_lag() {
    let equilibrium = FiniteSolubilitySolution::new(
        Temperature::new(300.0),
        18e-6,
        Volume::new(1e-3),
        vec![(law("A", 100.0, 6000.0), Amount::new(0.2))],
    )
    .unwrap();
    let owner = KineticSolubilitySolution::from_equilibrium(
        &equilibrium,
        &[(
            VolumetricFlowRate::new(1e-3),
            PrecipitationPolicyV1::ExistingSolidOnly,
        )],
    )
    .unwrap();
    assert_eq!(owner.solution(), equilibrium.solution());
    let mut diluted = owner.with_solvent_volume_m3(2e-3).unwrap();
    assert_eq!(
        diluted.last_steps()[0].state(),
        owner.last_steps()[0].state()
    );
    assert!(std::ptr::eq(diluted.models(), owner.models()));
    assert_eq!(diluted.last_steps()[0].equilibrium().dissolved.value(), 0.2);
    assert_eq!(diluted.solution().solute_moles(), &[0.1]);
    diluted.advance(Time::new(1.0), || false).unwrap();
    near(
        diluted.solution().solute_moles()[0],
        0.1 + 0.1 * (1.0 - (-0.5_f64).exp()),
    );
    let mut heated = owner.with_temperature(Temperature::new(330.0)).unwrap();
    assert_eq!(
        heated.last_steps()[0].state(),
        owner.last_steps()[0].state()
    );
    assert!(heated.last_steps()[0].equilibrium().capacity.value() > 0.1);
    heated.advance(Time::new(0.1), || false).unwrap();
    assert!(heated.solution().solute_moles()[0] > 0.1);
    assert!(
        heated.solution().solute_moles()[0]
            < heated.last_steps()[0].equilibrium().dissolved.value()
    );
    assert!(owner.with_temperature(Temperature::new(400.0)).is_err());
    assert_eq!(owner.solution(), equilibrium.solution());
}

#[test]
fn cancellation_and_late_aqueous_refusal_leave_all_component_stocks_unchanged() {
    let mut owner = reservoir();
    let original = owner.clone();
    let mut calls = 0;
    assert_eq!(
        owner.advance(Time::new(1.0), || {
            calls += 1;
            calls == 4
        }),
        Err(KineticSolubilityError::Cancelled)
    );
    assert_eq!(owner, original);
    assert!(owner.advance(Time::new(-1.0), || false).is_err());
    assert_eq!(owner, original);

    let mut excessive = KineticSolubilitySolution::new(
        Temperature::new(300.0),
        18e-6,
        Volume::new(1e-3),
        vec![component("A", 1e9, 1e6, 0.0, 1e-3)],
    )
    .unwrap();
    let original = excessive.clone();
    assert!(excessive.advance(Time::new(1.0), || false).is_err());
    assert_eq!(excessive, original);
}

#[test]
fn layout_and_source_identity_cannot_change_silently() {
    assert!(
        KineticSolubilitySolution::new(Temperature::new(300.0), 18e-6, Volume::new(1e-3), vec![])
            .is_err()
    );
    let duplicate = component("A", 100.0, 0.2, 0.0, 1e-3);
    assert!(
        KineticSolubilitySolution::new(
            Temperature::new(300.0),
            18e-6,
            Volume::new(1e-3),
            vec![duplicate.clone(), duplicate]
        )
        .is_err()
    );
    let equilibrium = FiniteSolubilitySolution::new(
        Temperature::new(300.0),
        18e-6,
        Volume::new(1e-3),
        vec![(law("A", 100.0, 0.0), Amount::new(0.2))],
    )
    .unwrap();
    assert!(KineticSolubilitySolution::from_equilibrium(&equilibrium, &[]).is_err());
}

#[test]
fn dissolved_transport_retains_crystals_and_admits_every_component_before_publication() {
    let mut owner = reservoir();
    owner.advance(Time::new(0.5), || false).unwrap();
    let solid = owner
        .last_steps()
        .iter()
        .map(|step| step.state().solid())
        .collect::<Vec<_>>();
    let phase_buffer = owner.last_steps().as_ptr();
    let aqueous_buffer = owner.solution().solute_moles().as_ptr();
    owner
        .replace_dissolved_inventories(0.00075, &[0.07, 0.02])
        .unwrap();
    for (index, &amount) in [0.07, 0.02].iter().enumerate() {
        let state = owner.last_steps()[index].state();
        assert_eq!(state.solid(), solid[index]);
        assert_eq!(state.dissolved().value(), amount);
        near(state.total_solute().value(), amount + solid[index].value());
        assert_eq!(
            owner.last_steps()[index].transferred_to_solution().value(),
            0.0
        );
    }
    near(
        owner.solution().osmotic_pressure_pa().unwrap(),
        R * 300.0 * 0.09 / 0.00075,
    );
    assert_eq!(owner.last_steps().as_ptr(), phase_buffer);
    assert_eq!(owner.solution().solute_moles().as_ptr(), aqueous_buffer);
    let accepted = owner.clone();
    for amounts in [&[0.03, -0.01][..], &[0.03, f64::NAN], &[0.03]] {
        assert!(owner.replace_dissolved_inventories(0.001, amounts).is_err());
        assert_eq!(owner, accepted);
    }
    // Every phase can be formed, but the resulting activity is inadmissible.
    assert!(
        owner
            .replace_dissolved_inventories(1e-30, &[0.03, 0.04])
            .is_err()
    );
    assert_eq!(owner, accepted);
}

#[test]
fn prepared_transport_copies_preserve_trace_solids_and_reject_different_sources() {
    let state =
        SolubilityKineticStateV1::from_phases(Amount::new(1.0), Amount::new(1e-20)).unwrap();
    let model = SolubilityKineticsV1::new(
        law("A", 100.0, 0.0),
        VolumetricFlowRate::new(1e-3),
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
    .unwrap();
    let mut source = KineticSolubilitySolution::new(
        Temperature::new(300.0),
        18e-6,
        Volume::new(1e-3),
        vec![(model, state)],
    )
    .unwrap();
    let mut prepared = source.clone();
    let buffer = prepared.last_steps().as_ptr();
    source.replace_dissolved_inventories(2e-3, &[0.5]).unwrap();
    assert_eq!(source.last_steps()[0].state().solid().value(), 1e-20);
    assert_eq!(source.last_steps()[0].state().total_solute().value(), 0.5);
    prepared.copy_state_from(&source).unwrap();
    assert_eq!(prepared, source);
    assert_eq!(prepared.last_steps().as_ptr(), buffer);
    assert_ne!(prepared.last_steps().as_ptr(), source.last_steps().as_ptr());
    let before = prepared.clone();
    let incompatible = source.with_temperature(Temperature::new(310.0)).unwrap();
    assert!(prepared.copy_state_from(&incompatible).is_err());
    assert_eq!(prepared, before);
    assert!(prepared.heap_bytes() >= core::mem::size_of_val(prepared.last_steps()));
    assert_eq!(prepared.solvent(), &SpeciesId::new("water").unwrap());
}

fn particle_reservoir(reference_solid: f64) -> KineticSolubilitySolution {
    let model = SolubilityKineticsV1::new(
        law("A", 125.0, 0.0),
        VolumetricFlowRate::new(5e-4),
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
    .unwrap()
    .with_particle_surface(Amount::new(reference_solid))
    .unwrap();
    KineticSolubilitySolution::new(
        Temperature::new(300.0),
        18e-6,
        Volume::new(1e-3),
        vec![(
            model,
            SolubilityKineticStateV1::from_phases(Amount::new(0.125), Amount::new(0.125)).unwrap(),
        )],
    )
    .unwrap()
}

#[test]
fn drying_grows_particle_area_and_dilution_retains_it_through_analytical_redissolution() {
    let original = particle_reservoir(0.125);
    let mut concentrated = original.with_solvent_volume_m3(5e-4).unwrap();
    assert_eq!(
        concentrated.last_steps()[0].state(),
        original.last_steps()[0].state()
    );
    concentrated.advance(Time::new(0.4), || false).unwrap();
    let grown = concentrated.last_steps()[0].state().solid();
    assert!(grown.value() > 0.125);
    assert!(
        concentrated.models()[0]
            .conductance_at(grown)
            .unwrap()
            .value()
            > 5e-4
    );

    let mut diluted = concentrated.with_solvent_volume_m3(2e-3).unwrap();
    assert!(std::ptr::eq(diluted.models(), original.models()));
    assert_eq!(diluted.last_steps()[0].state().solid(), grown);
    let phase_buffer = diluted.last_steps().as_ptr();
    let aqueous_buffer = diluted.solution().solute_moles().as_ptr();
    diluted.advance(Time::new(3.0), || false).unwrap();

    // V*c_sat equals the complete component stock. Integrating the spherical
    // area law then gives s(t)=s0/[1+(2/3)*(Gref/V)*(s0/sref)^(2/3)*t]^(3/2).
    // This independent endpoint also detects accidental reference-area reset
    // during concentration, dilution or aqueous-buffer reconstruction.
    let expected = grown.value()
        / (1.0 + (2.0 / 3.0) * (5e-4 / 2e-3) * (grown.value() / 0.125).powf(2.0 / 3.0) * 3.0)
            .powf(1.5);
    let phase = diluted.last_steps()[0].state();
    assert!((phase.solid().value() / expected - 1.0).abs() < 2e-10);
    near(phase.total_solute().value(), 0.25);
    near(phase.dissolved().value() + phase.solid().value(), 0.25);
    assert!(phase.solid().value() < grown.value());
    assert!(
        diluted.models()[0]
            .conductance_at(phase.solid())
            .unwrap()
            .value()
            < 5e-4
    );
    near(
        diluted.solution().osmotic_pressure_pa().unwrap(),
        R * 300.0 * phase.dissolved().value() / 2e-3,
    );
    assert_eq!(diluted.last_steps().as_ptr(), phase_buffer);
    assert_eq!(diluted.solution().solute_moles().as_ptr(), aqueous_buffer);
}

#[test]
fn dissolved_transport_retains_particle_geometry_and_incompatible_checkpoint_refuses() {
    let mut original = particle_reservoir(0.125);
    original
        .replace_dissolved_inventories(1e-3, &[0.0])
        .unwrap();
    let initial = *original.last_steps()[0].state();
    assert_eq!(initial.solid().value(), 0.125);
    assert_eq!(
        original.models()[0]
            .conductance_at(initial.solid())
            .unwrap()
            .value(),
        5e-4
    );
    original.advance(Time::new(0.5), || false).unwrap();
    let expected_solid = 0.125 / (1.0_f64 + (2.0 / 3.0) * 0.5 * 0.5).powf(1.5);
    assert!(
        (original.last_steps()[0].state().solid().value() / expected_solid - 1.0).abs() < 2e-10
    );
    assert_eq!(
        original.models()[0].conductance(),
        VolumetricFlowRate::new(5e-4)
    );
    assert!(
        original.models()[0]
            .saturation_deficit_half_life(Volume::new(1e-3))
            .is_err()
    );

    let mut other_geometry = particle_reservoir(0.25);
    let before = other_geometry.clone();
    assert!(other_geometry.copy_state_from(&original).is_err());
    assert_eq!(other_geometry, before);

    let mut stopped = original.clone();
    let before = stopped.clone();
    let mut polls = 0;
    assert_eq!(
        stopped.advance(Time::new(0.5), || {
            polls += 1;
            polls == 3
        }),
        Err(KineticSolubilityError::Cancelled)
    );
    assert_eq!(stopped, before);
}
