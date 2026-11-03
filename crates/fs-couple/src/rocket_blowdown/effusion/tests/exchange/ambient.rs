//! Fixed inertial sources: original-source parity, vector moments and rollback.
use super::*;
use crate::rocket_blowdown::effusion::exchange::EffusionAmbientPort;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::effusion::drift::EffusiveDriftInterface;

fn drift_law(model: &ThermochemicalGasModel, maximum: f64) -> EffusiveDriftInterface {
    EffusiveDriftInterface::try_new(
        interface(model),
        maximum,
        "synthetic collisionless nondeflecting aperture over the declared drift ball",
    )
    .unwrap()
}

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

fn advance_ambient(
    before: &GasRocketState,
    ports: &[EffusionAmbientPort<'_>],
) -> EffusionExchangeCandidate {
    before
        .advance_effusion_exchange_with_ambient(
            environment(),
            &[],
            &[],
            ports,
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap()
}

#[test]
fn initially_comoving_ambient_matches_original_ports_and_tangent_slip_heats_once() {
    let before = initial(PURE, Vec3::new(3.0, -2.0, 1.0));
    let source = bath(&before, 600.0, 2.0 * P0);
    let saved_source = source.clone();
    let law = drift_law(before.reservoir().model(), 1000.0);
    let supply = ambient_port(&source, &law, before.vehicle().velocity_m_s);
    let common = EffusionReservoirPort {
        area_m2: 0.5 * AREA,
        ..port(&source, law.momentum_interface())
    };
    let old_supply = EffusionReservoirPort {
        id: supply.id,
        ..port(&source, law.momentum_interface())
    };
    let vacuum = [outlet(law.momentum_interface(), Vec3::new(0.0, 0.0, -1.0))];
    let old = advance(&before, environment(), &vacuum, &[old_supply, common]);
    let new = before
        .advance_effusion_exchange_with_ambient(
            environment(),
            &vacuum,
            &[common],
            &[supply],
            DT,
            endpoint(),
            controls(),
            || false,
        )
        .unwrap();
    near(
        new.endpoint().reservoir().mass_kg(),
        old.endpoint().reservoir().mass_kg(),
        2e-15,
    );
    near(
        new.endpoint().reservoir().state().state().temperature,
        old.endpoint().reservoir().state().state().temperature,
        2e-7,
    );
    vector_near(
        new.endpoint().vehicle().velocity_m_s,
        old.endpoint().vehicle().velocity_m_s,
        2e-12,
    );
    vector_near(
        new.endpoint().vehicle().position_m,
        old.endpoint().vehicle().position_m,
        2e-13,
    );
    for (a, b) in [
        (new.step().incoming, old.step().incoming),
        (new.step().outgoing, old.step().outgoing),
    ] {
        near(a.mass_kg, b.mass_kg, 2e-18);
        vector_near(
            a.relative_momentum_kg_m_s,
            b.relative_momentum_kg_m_s,
            2e-15,
        );
        near(a.relative_energy_j, b.relative_energy_j, 2e-12);
    }
    assert_eq!(
        new.step().ports.iter().map(|p| p.id).collect::<Vec<_>>(),
        vec![1, 4, 9]
    );
    assert_eq!(new.step().gas.drifting_reservoirs.len(), 1);
    assert_eq!(new.step().gas.reservoirs.len(), 1);
    accounts(&before, &new);
    assert_eq!(source, saved_source);

    // Tangential drift leaves the normal arrival count unchanged. Each actual
    // incoming molecule adds its tangential momentum and |u_t|²/2 exactly once.
    let before = initial(PURE, Vec3::ZERO);
    let source = bath(&before, T0, P0);
    let law = drift_law(before.reservoir().model(), 1000.0);
    let tangent = Vec3::new(0.0, 80.0, -50.0);
    let result = advance_ambient(&before, &[ambient_port(&source, &law, tangent)]);
    let dm = escaped_moles(0.028, 1.0) * 0.028;
    let impulse = 0.5 * P0 * AREA * DT;
    let pin = Vec3::new(-impulse, dm * tangent.y, dm * tangent.z);
    let pout = Vec3::new(impulse, 0.0, 0.0);
    near(result.step().incoming.mass_kg, dm, 2e-18);
    near(result.step().outgoing.mass_kg, dm, 2e-18);
    near(result.step().net_mass_into_gas_kg, 0.0, 2e-18);
    vector_near(result.step().incoming.relative_momentum_kg_m_s, pin, 2e-15);
    vector_near(result.step().outgoing.relative_momentum_kg_m_s, pout, 2e-15);
    let heat = 0.5 * dm * tangent.norm_squared();
    near(result.step().net_relative_energy_into_gas_j, heat, 2e-12);
    near(
        result.endpoint().reservoir().state().state().temperature,
        T0 + heat / (before.reservoir().mass_kg() * 1.5 * R / 0.028),
        2e-7,
    );
    vector_near(
        result.endpoint().vehicle().velocity_m_s,
        pin.sub(pout).scale(1.0 / mass(&before)),
        2e-12,
    );
    accounts(&before, &result);
}

#[test]
fn headwind_carries_half_range_energy_and_both_inertial_velocities_boost_together() {
    let velocity = Vec3::new(2.0, -1.0, 0.5);
    let before = initial(PURE, velocity);
    let source = bath(&before, T0, P0);
    let law = drift_law(before.reservoir().model(), 1000.0);
    let sigma2 = R * T0 / 0.028;
    let sigma = sigma2.sqrt();
    let relative = Vec3::new(-sigma, 0.6 * sigma, -0.25 * sigma);
    let ambient = ambient_port(&source, &law, velocity.add(relative));
    let result = advance_ambient(&before, &[ambient]);

    // Independent analytic moments for exactly one thermal-speed headwind:
    // Hk = integral_0^infinity z^k phi(z-1) dz. Phi(1) is a fixed mathematical
    // reference, and these expressions do not call the production drift owner.
    let phi = (-0.5_f64).exp() / (2.0 * core::f64::consts::PI).sqrt();
    let cdf = 0.841_344_746_068_542_9;
    let h1 = phi + cdf;
    let h2 = phi + 2.0 * cdf;
    let h3 = 3.0 * phi + 4.0 * cdf;
    let mi = AREA * DT * (P0 / sigma2) * sigma * h1;
    let mo = escaped_moles(0.028, 1.0) * 0.028;
    let pin = Vec3::new(-AREA * DT * P0 * h2, mi * relative.y, mi * relative.z);
    let pout = Vec3::new(0.5 * AREA * DT * P0, 0.0, 0.0);
    let ei = mi
        * (sigma2 * (1.0 + 0.5 * h3 / h1)
            + 0.5 * (relative.y * relative.y + relative.z * relative.z));
    let eo = mo * 2.0 * sigma2;
    near(result.step().incoming.mass_kg, mi, 2e-18);
    assert!(mi > 2.0 * mo);
    vector_near(result.step().incoming.relative_momentum_kg_m_s, pin, 2e-15);
    vector_near(result.step().outgoing.relative_momentum_kg_m_s, pout, 2e-15);
    near(result.step().incoming.relative_energy_j, ei, 3e-12);
    near(result.step().outgoing.relative_energy_j, eo, 3e-12);
    let dm = mi - mo;
    let expected_temperature =
        (internal(&before) + ei - eo) / ((before.reservoir().mass_kg() + dm) * 1.5 * R / 0.028);
    near(
        result.endpoint().reservoir().state().state().temperature,
        expected_temperature,
        2e-7,
    );
    let increment = pin.sub(pout).scale((dm / mass(&before)).ln_1p() / dm);
    vector_near(
        result.endpoint().vehicle().velocity_m_s,
        velocity.add(increment),
        2e-12,
    );
    accounts(&before, &result);

    let boost = Vec3::new(17.0, -9.0, 23.0);
    let boosted = initial(PURE, velocity.add(boost));
    let shifted = advance_ambient(
        &boosted,
        &[EffusionAmbientPort {
            velocity_world_m_s: ambient.velocity_world_m_s.add(boost),
            ..ambient
        }],
    );
    vector_near(
        shifted.endpoint().vehicle().velocity_m_s,
        result.endpoint().vehicle().velocity_m_s.add(boost),
        3e-12,
    );
    vector_near(
        shifted.endpoint().vehicle().position_m,
        result.endpoint().vehicle().position_m.add(boost.scale(DT)),
        3e-12,
    );
    near(
        shifted.endpoint().reservoir().state().state().temperature,
        result.endpoint().reservoir().state().state().temperature,
        2e-7,
    );
    for (a, b) in [
        (result.step().incoming, shifted.step().incoming),
        (result.step().outgoing, shifted.step().outgoing),
    ] {
        near(a.mass_kg, b.mass_kg, 2e-18);
        vector_near(
            b.momentum_world_kg_m_s,
            a.momentum_world_kg_m_s.add(boost.scale(a.mass_kg)),
            2e-14,
        );
        near(
            b.energy_world_j,
            a.energy_world_j
                + boost.dot(a.momentum_world_kg_m_s)
                + 0.5 * a.mass_kg * boost.norm_squared(),
            3e-10,
        );
    }
    accounts(&boosted, &shifted);
}

#[test]
fn final_drift_frame_source_and_late_cancel_refusals_publish_no_endpoint() {
    let before = initial(PURE, Vec3::new(4.0, -3.0, 2.0));
    let saved = before.clone();
    let source = bath(&before, T0, P0);
    let law = drift_law(before.reservoir().model(), 1000.0);
    let path = ambient_port(&source, &law, before.vehicle().velocity_m_s);
    let bounded = drift_law(before.reservoir().model(), 0.0);
    // The initial drift is exactly zero. The real aperture impulse moves the
    // body, so its endpoint leaves this authored zero-drift qualification.
    assert!(matches!(
        before.advance_effusion_exchange_with_ambient(
            environment(),
            &[],
            &[],
            &[EffusionAmbientPort {
                interface: &bounded,
                ..path
            }],
            DT,
            endpoint(),
            controls(),
            || false,
        ),
        Err(EffusionRocketError::Source { outlet: 1, .. })
    ));
    assert!(matches!(
        before.advance_effusion_exchange_with_ambient(
            environment(),
            &[],
            &[],
            &[EffusionAmbientPort {
                inertial_frame: FRAME + 1,
                ..path
            }],
            DT,
            endpoint(),
            controls(),
            || false,
        ),
        Err(EffusionRocketError::Invalid(_))
    ));
    let foreign = initial("A 1 0.028 3.5 0 1.7e-5 0.025", Vec3::ZERO);
    let foreign_path = EffusionAmbientPort {
        source: foreign.reservoir(),
        ..path
    };
    assert!(
        before
            .advance_effusion_exchange_with_ambient(
                environment(),
                &[],
                &[],
                &[foreign_path],
                DT,
                endpoint(),
                controls(),
                || false,
            )
            .is_err()
    );

    let mut polls = 0;
    let accepted = before
        .advance_effusion_exchange_with_ambient(
            environment(),
            &[],
            &[],
            &[path],
            DT,
            endpoint(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut count = 0;
    assert!(matches!(
        before.advance_effusion_exchange_with_ambient(
            environment(),
            &[],
            &[],
            &[path],
            DT,
            endpoint(),
            controls(),
            || {
                count += 1;
                count == polls
            },
        ),
        Err(EffusionRocketError::Cancelled)
    ));
    assert_eq!(advance_ambient(&before, &[path]), accepted);
    let coast = advance_ambient(
        &before,
        &[EffusionAmbientPort {
            area_m2: 0.0,
            velocity_world_m_s: Vec3::new(f64::NAN, 0.0, 0.0),
            ..foreign_path
        }],
    );
    assert_eq!(coast.endpoint().reservoir(), before.reservoir());
    assert_eq!(
        coast.endpoint().vehicle().velocity_m_s,
        before.vehicle().velocity_m_s
    );
    vector_near(
        coast.endpoint().vehicle().position_m,
        before.vehicle().velocity_m_s.scale(DT),
        2e-14,
    );
    assert_eq!(coast.step().incoming.mass_kg, 0.0);
    assert_eq!(coast.step().outgoing.mass_kg, 0.0);
    assert_eq!(before, saved);
}
