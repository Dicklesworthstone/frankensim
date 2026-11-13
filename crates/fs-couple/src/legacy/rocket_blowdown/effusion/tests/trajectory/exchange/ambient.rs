//! Fixed-world ambient resampling, physical refinement and original species.
use super::*;
use crate::rocket_blowdown::effusion::exchange::EffusionAmbientPort;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::effusion::drift::EffusiveDriftInterface;

fn ambient_port<'a>(
    source: &'a ClosedGasParcel,
    law: &'a EffusiveDriftInterface,
    velocity: Vec3,
) -> EffusionAmbientPort<'a> {
    EffusionAmbientPort {
        id: 1,
        source,
        velocity_world_m_s: velocity,
        inertial_frame: FRAME,
        interface: law,
        area_m2: AREA,
        outward_normal_world: Vec3::new(1.0, 0.0, 0.0),
    }
}

fn evolve_ambient(
    before: &GasRocketState,
    ports: &[EffusionAmbientPort<'_>],
    duration: f64,
    accuracy: EffusionAdaptiveSolve,
) -> EffusionExchangeTrajectoryCandidate {
    before
        .advance_effusion_exchange_with_ambient_adaptive(
            quiet(),
            &[],
            &[],
            ports,
            duration,
            endpoint(),
            controls(),
            accuracy,
            || false,
        )
        .unwrap()
}

#[test]
fn adaptive_ambient_refreshes_relative_velocity_and_converges_to_refined_physics() {
    let before = initial(PURE, Vec3::new(150.0, -40.0, 20.0));
    let saved = before.clone();
    let source = ClosedGasParcel::new(
        before.reservoir().model().clone(),
        source::point(T0, P0).unwrap(),
        VOLUME,
    )
    .unwrap();
    let saved_source = source.clone();
    let law = EffusiveDriftInterface::try_new(
        interface(before.reservoir().model()),
        1000.0,
        "synthetic local Maxwellian ambient over the authored collisionless drift ball",
    )
    .unwrap();
    let ports = [EffusionAmbientPort {
        area_m2: 8.0 * AREA,
        ..ambient_port(&source, &law, Vec3::ZERO)
    }];
    let duration = 0.5;
    let loose = evolve_ambient(&before, &ports, duration, adaptive(1.0));
    let tight = evolve_ambient(&before, &ports, duration, adaptive(0.01));
    let mut reference = before.clone();
    for _ in 0..512 {
        reference = reference
            .advance_effusion_exchange_with_ambient(
                quiet(),
                &[],
                &[],
                &ports,
                duration / 512.0,
                endpoint(),
                controls(),
                || false,
            )
            .unwrap()
            .into_parts()
            .0;
    }
    // Dimensionless endpoint distance in thermal, inventory and velocity
    // coordinates. The independent fine reference uses only native steps,
    // without the adaptive controller or its acceptance decisions.
    let distance = |state: &GasRocketState| {
        (state.reservoir().state().state().temperature
            - reference.reservoir().state().state().temperature)
            .abs()
            / T0
            + (state.reservoir().mass_kg() - reference.reservoir().mass_kg()).abs()
                / before.reservoir().mass_kg()
            + state
                .vehicle()
                .velocity_m_s
                .sub(reference.vehicle().velocity_m_s)
                .norm_squared()
                .sqrt()
                / (R * T0 / 0.028).sqrt()
    };
    assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
    assert!(tight.step().rejected_intervals > 0);
    assert!(distance(tight.endpoint()) < 0.5 * distance(loose.endpoint()));
    assert_eq!(tight.step().duration_s, duration);
    assert!(tight.step().maximum_accepted_error <= 1.0);

    let mut velocity = before.vehicle().velocity_m_s;
    for accepted in &tight.step().accepted_steps {
        let sample = &accepted.gas.drifting_reservoirs[0].incoming;
        let expected = ports[0].velocity_world_m_s.sub(velocity);
        assert_eq!(
            sample.relative_bulk_velocity_m_s,
            [expected.x, expected.y, expected.z]
        );
        velocity = accepted.mechanics.state_after.velocity_m_s;
    }
    let first = &tight
        .step()
        .accepted_steps
        .first()
        .unwrap()
        .gas
        .drifting_reservoirs[0]
        .incoming;
    let last = &tight
        .step()
        .accepted_steps
        .last()
        .unwrap()
        .gas
        .drifting_reservoirs[0]
        .incoming;
    assert!((last.relative_bulk_velocity_m_s[0] - first.relative_bulk_velocity_m_s[0]).abs() > 1.0);
    assert_eq!(ports[0].velocity_world_m_s, Vec3::ZERO);
    let after = tight.endpoint();
    let receipt = tight.step();
    near(
        after.reservoir().mass_kg() - before.reservoir().mass_kg(),
        receipt.total_net_mass_into_gas_kg,
        1e-13,
    );
    vector_near(
        after
            .vehicle()
            .velocity_m_s
            .scale(mass(after))
            .sub(before.vehicle().velocity_m_s.scale(mass(&before)))
            .sub(receipt.total_net_momentum_into_vessel_world_kg_m_s)
            .sub(receipt.total_external_impulse_n_s),
        Vec3::ZERO,
        1e-10,
    );
    near(
        internal(after) - internal(&before) + kinetic(after)
            - kinetic(&before)
            - receipt.total_net_energy_into_vessel_world_j
            - receipt.total_external_work_j,
        0.0,
        3e-7,
    );
    assert_eq!(
        after.vehicle().propellant_mass_kg,
        after.reservoir().mass_kg()
    );
    assert_eq!(before, saved);
    assert_eq!(source, saved_source);
}

#[test]
fn ambient_species_join_enters_error_axis_and_late_refusal_keeps_original_owners() {
    let catalog = source::load(MIXTURE, "effusive recoil analytic fixture").unwrap();
    let composition = |fractions| {
        catalog
            .with_composition(Composition::new(CompositionBasis::MoleFraction, fractions).unwrap())
            .unwrap()
    };
    let gas = ClosedGasParcel::new(
        composition(vec![1.0, 0.0]),
        source::point(T0, P0).unwrap(),
        VOLUME,
    )
    .unwrap();
    let before = GasRocketState::new(gas, 0.01, Vec3::ZERO, Vec3::ZERO, FRAME).unwrap();
    let saved = before.clone();
    let source = ClosedGasParcel::new(
        composition(vec![0.0, 1.0]),
        source::point(T0, P0).unwrap(),
        VOLUME,
    )
    .unwrap();
    let saved_source = source.clone();
    let a = SpeciesId::new("A").unwrap();
    let b = SpeciesId::new("B").unwrap();
    assert_eq!(
        before
            .reservoir()
            .model()
            .thermodynamics()
            .components()
            .len(),
        1
    );
    assert_eq!(before.reservoir().species_mass_kg(&b).unwrap(), 0.0);
    let original = interface(&catalog);
    let transmissions = catalog
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                if species.species() == &b { 1.0 } else { 0.0 },
                "synthetic B-only molecular aperture",
            )
            .unwrap()
        })
        .collect();
    let law = EffusiveDriftInterface::try_new(
        EffusiveMomentumInterface::try_new(
            EffusiveInterface::try_new(
                catalog,
                transmissions,
                original.interface().support(),
                original.interface().regime().clone(),
            )
            .unwrap(),
            "synthetic uncollimated nondeflecting B transmission",
        )
        .unwrap(),
        1000.0,
        "synthetic B ambient with qualified relative drift",
    )
    .unwrap();
    let ports = [ambient_port(&source, &law, Vec3::new(-40.0, 60.0, 0.0))];
    let accuracy = EffusionAdaptiveSolve {
        initial_step_s: 0.25,
        maximum_step_s: 0.25,
        ..adaptive(1.0)
    };
    let mut polls = 0;
    let result = before
        .advance_effusion_exchange_with_ambient_adaptive(
            quiet(),
            &[],
            &[],
            &ports,
            0.25,
            endpoint(),
            controls(),
            accuracy,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(
        result
            .endpoint()
            .reservoir()
            .model()
            .thermodynamics()
            .components()
            .len(),
        2
    );
    assert!(result.endpoint().reservoir().species_mass_kg(&b).unwrap() > 0.0);
    near(
        result.endpoint().reservoir().species_mass_kg(&a).unwrap(),
        before.reservoir().mass_kg(),
        2e-14,
    );
    assert!(result.step().total_incoming.mass_kg > 0.0);
    assert!(result.step().total_outgoing.mass_kg > 0.0);
    let mut count = 0;
    assert!(matches!(
        before.advance_effusion_exchange_with_ambient_adaptive(
            quiet(),
            &[],
            &[],
            &ports,
            0.25,
            endpoint(),
            controls(),
            accuracy,
            || {
                count += 1;
                count == polls
            },
        ),
        Err(EffusionTrajectoryError::Cancelled)
    ));

    let foreign = initial("A 1 0.028 3.5 0 1.7e-5 0.025", Vec3::ZERO);
    let coast = evolve_ambient(
        &before,
        &[EffusionAmbientPort {
            source: foreign.reservoir(),
            area_m2: 0.0,
            velocity_world_m_s: Vec3::new(f64::NAN, 0.0, 0.0),
            ..ports[0]
        }],
        0.25,
        accuracy,
    );
    assert_eq!(coast.endpoint().reservoir(), before.reservoir());
    assert_eq!(coast.endpoint().vehicle(), before.vehicle());
    assert_eq!(coast.step().total_incoming.mass_kg, 0.0);
    assert_eq!(coast.step().total_outgoing.mass_kg, 0.0);
    assert_eq!(before, saved);
    assert_eq!(source, saved_source);
}
