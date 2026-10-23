//! Real source-material evolution on retained native face momentum.
use super::*;
use fs_lbm::conjugate::{FluidProperties, ViscousHeatControls, simple_flow_with_viscosity};

fn cavity(lid_speed: f64) -> SimpleConfig {
    let mut faces = [FvBoundary::Wall { velocity: [0.0; 3] }; 6];
    faces[3] = FvBoundary::Wall {
        velocity: [lid_speed, 0.0, 0.0],
    };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut flow = SimpleConfig::new(faces);
    flow.tolerance = 1e-11;
    flow.pressure_tolerance = 1e-12;
    flow
}

fn velocity_difference(a: &fs_lbm::conjugate::FvFlow, b: &fs_lbm::conjugate::FvFlow) -> f64 {
    a.velocity_m_s
        .iter()
        .zip(&b.velocity_m_s)
        .flat_map(|(a, b)| a.iter().zip(b).map(|(a, b)| (a - b).abs()))
        .fold(0.0, f64::max)
}

#[test]
fn heat_changes_subsequent_momentum_while_old_velocity_alone_advects_h() {
    let card = feedback_card(1.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let flow_config = cavity(0.1);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let heated = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow_config, Heater::default(), &gate, 100)
        .unwrap()
        .with_unsteady_momentum()
        .unwrap();
    let cold = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow_config, (), &gate, 100)
        .unwrap()
        .with_unsteady_momentum()
        .unwrap();
    let (hot1, history1, audit1) = old
        .advance_observed(&heated, 0.125, config(), || Ok(()))
        .unwrap();
    let (cold1, cold_history1, _) = old
        .advance_observed(&cold, 0.125, config(), || Ok(()))
        .unwrap();
    // The first momentum step uses the same h^0 material epoch, even though
    // the thermal endpoints differ. Initial old velocity is exactly at rest.
    assert_eq!(history1.unsteady_flow(), cold_history1.unsteady_flow());
    assert!(
        history1
            .applied()
            .unwrap()
            .flow
            .velocity_m_s
            .iter()
            .all(|v| *v == [0.0; 3])
    );
    assert!(history1.unsteady_step().unwrap().kinetic_energy_j > 0.0);
    near(hot1.cells()[&CellId(1)].specific_h_j_kg, 12.5);
    near(audit1.external_heat_j, 12.5);
    let (hot2, history2, audit2) = hot1
        .advance_observed(&history1, 0.225, config(), || Ok(()))
        .unwrap();
    let (_, cold_history2, _) = cold1
        .advance_observed(&cold_history1, 0.225, config(), || Ok(()))
        .unwrap();
    assert!(
        velocity_difference(
            history2.unsteady_flow().unwrap(),
            cold_history2.unsteady_flow().unwrap()
        ) > 1e-7
    );
    assert_eq!(
        &history2.applied().unwrap().flow,
        history1.unsteady_flow().unwrap()
    );
    assert_eq!(history2.applied().unwrap().time_s, 0.125);
    assert_eq!(history2.unsteady_time_s(), Some(hot2.time_s()));
    assert_eq!(history2.unsteady_step().unwrap().time_s, hot2.time_s());
    let applied = history2.applied().unwrap();
    near(
        applied.samples[1].transport.dynamic_viscosity_pa_s(),
        viscosity(1.0, 12.5),
    );
    for sample in &applied.samples {
        assert_eq!(
            sample.material.material().card_identity(),
            card.content_hash()
        );
    }
    // A newly constructed rest state at the SAME thermal time gives a
    // different result: accepted momentum, not just material data, survived.
    let restarted = input(&ids, &liquid, &caloric)
        .feedback(&hot1, &domain, &flow_config, Heater::default(), &gate, 100)
        .unwrap()
        .with_unsteady_momentum()
        .unwrap();
    assert_eq!(restarted.unsteady_time_s(), Some(0.125));
    let (_, rest2, _) = hot1
        .advance_observed(&restarted, 0.225, config(), || Ok(()))
        .unwrap();
    assert!(
        velocity_difference(
            history2.unsteady_flow().unwrap(),
            rest2.unsteady_flow().unwrap()
        ) > 1e-5
    );
    near(hot2.energy_j() - old.energy_j(), 22.5);
    near(audit2.external_heat_j, 10.0);
    near(history2.audit().external_enthalpy_j, 0.0);
    near(hot2.mass_kg(), old.mass_kg());
    assert_eq!(history2.momentum_solves(), 1);
    assert_eq!(heated.unsteady_time_s(), Some(0.0));
    assert!(heated.unsteady_step().is_none());
}

#[test]
fn native_transient_reaches_existing_full_stress_steady_limit() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let flow_config = cavity(0.1);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow_config, (), &gate, 100)
        .unwrap()
        .with_unsteady_momentum()
        .unwrap();
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 1.0,
    };
    let steady =
        simple_flow_with_viscosity(&domain, &fluid, &[1.0; 4], &flow_config, &gate).unwrap();
    let (first, first_history, _) = old
        .advance_observed(&observer, 0.125, config(), || Ok(()))
        .unwrap();
    assert!(velocity_difference(first_history.unsteady_flow().unwrap(), &steady) > 1e-4);
    let (next, accepted, audit) = first
        .advance_observed(&first_history, 8.0, config(), || Ok(()))
        .unwrap();
    assert!(velocity_difference(accepted.unsteady_flow().unwrap(), &steady) < 1e-8);
    assert!(accepted.unsteady_flow().unwrap().report.max_divergence_m3_s < 1e-12);
    assert_eq!(accepted.unsteady_time_s(), Some(8.0));
    near(next.energy_j(), old.energy_j());
    near(audit.external_heat_j, 0.0);
    assert_eq!(accepted.momentum_solves(), audit.substeps);
    assert!(accepted.momentum_solves() >= 63);
}

#[test]
fn actual_combined_cfl_interval_drives_exactly_one_momentum_step() {
    let card = feedback_card(20.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let flow_config = cavity(10.0);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow_config, Heater::default(), &gate, 100)
        .unwrap()
        .with_unsteady_momentum()
        .unwrap();
    let (next, accepted, audit) = old
        .advance_observed(
            &observer,
            6.0,
            StepConfig {
                max_dt_s: 2.0,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
    let durations = &accepted.inner().durations;
    assert!(durations.len() > 3);
    assert!(durations[1..durations.len() - 1].iter().any(|dt| *dt < 2.0));
    assert_eq!(accepted.momentum_solves(), durations.len());
    assert_eq!(accepted.audit().substeps, durations.len());
    assert_eq!(accepted.unsteady_time_s(), Some(next.time_s()));
    assert_eq!(accepted.unsteady_step().unwrap().time_s, 6.0);
    near(durations.iter().sum(), 6.0);
    near(next.energy_j() - old.energy_j(), 600.0);
    near(audit.external_heat_j, 600.0);
}

#[test]
fn thermal_and_momentum_refusal_cancellation_and_stale_timeline_publish_nothing() {
    let card = feedback_card(1.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let flow_config = cavity(0.1);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(
            &old,
            &domain,
            &flow_config,
            Heater {
                reject_step: Some(2),
                ..Heater::default()
            },
            &gate,
            100,
        )
        .unwrap()
        .with_unsteady_momentum()
        .unwrap();
    let (first, history, _) = old
        .advance_observed(&observer, 0.125, config(), || Ok(()))
        .unwrap();
    let saved_flow = history.unsteady_flow().unwrap().clone();
    let saved_step = *history.unsteady_step().unwrap();
    let saved_inner = history.inner().clone();
    let saved_endpoint = history.endpoint().to_vec();
    assert!(
        first
            .advance_observed(&history, 0.25, config(), || Ok(()))
            .is_err()
    );
    assert_eq!(history.unsteady_flow(), Some(&saved_flow));
    assert_eq!(history.unsteady_step(), Some(&saved_step));
    assert_eq!(history.inner(), &saved_inner);
    assert_eq!(history.endpoint(), saved_endpoint);
    assert_eq!(first.time_s(), 0.125);
    near(first.energy_j(), 12.5);
    // Mismatched clocks refuse before any new material/source history.
    assert!(
        first
            .advance_observed(&observer, 0.25, config(), || Ok(()))
            .is_err()
    );
    assert!(
        old.advance_observed(&history, 0.125, config(), || Ok(()))
            .is_err()
    );
    let (other, _) = old.advance_to(0.125, config(), || Ok(())).unwrap();
    assert!(
        other
            .advance_observed(&history, 0.25, config(), || Ok(()))
            .is_err()
    );
    // Force native momentum refusal after thermal/source work has happened.
    let mut insufficient = flow_config.clone();
    insufficient.max_iterations = 1;
    let refused = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &insufficient, Heater::default(), &gate, 100)
        .unwrap()
        .with_unsteady_momentum()
        .unwrap();
    assert!(matches!(
        old.advance_observed(&refused, 0.125, config(), || Ok(())),
        Err(MaterialFlowError::Momentum(_))
    ));
    assert_eq!(refused.unsteady_time_s(), Some(0.0));
    assert!(refused.unsteady_step().is_none());
    assert_eq!(refused.inner(), &Heater::default());
    // Cancel at the final checkpoint of a proposal which has completed native
    // momentum. The original field, pressure, material and inner history remain.
    let mut checks = 0;
    old.advance_observed(&observer, 0.125, config(), || {
        checks += 1;
        Ok(())
    })
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        old.advance_observed(&observer, 0.125, config(), || {
            seen += 1;
            if seen == checks {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(MaterialFlowError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert!(observer.unsteady_step().is_none());
    assert_eq!(observer.inner().durations.len(), 0);
    assert!(observer.applied().is_none());
    near(old.energy_j(), 0.0);
}

#[test]
fn incompatible_boundaries_and_steady_heat_refuse_in_both_attachment_orders() {
    let card = feedback_card(1.0);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let flow_config = cavity(0.1);
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, 1.0, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &flow_config, (), &gate, 100)
        .unwrap();
    let heat = ViscousHeatControls {
        maximum_reynolds: 0.1,
        power_tolerance_w: 1e-8,
    };
    assert!(
        observer
            .clone()
            .with_unsteady_momentum()
            .unwrap()
            .with_viscous_heating(heat)
            .is_err()
    );
    assert!(
        observer
            .clone()
            .with_viscous_heating(heat)
            .unwrap()
            .with_unsteady_momentum()
            .is_err()
    );
    let mut open = flow_config.clone();
    open.faces[0] = FvBoundary::Outlet;
    assert!(
        input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &open, (), &gate, 100)
            .unwrap()
            .with_unsteady_momentum()
            .is_err()
    );
}
