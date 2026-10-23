//! Latent enthalpy, not temperature-only refresh, drives retained momentum.
use super::*;

fn cavity() -> SimpleConfig {
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall {
        velocity: [0.1, 0.0, 0.0],
    };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut flow = SimpleConfig::new(faces);
    flow.tolerance = 1e-11;
    flow.pressure_tolerance = 1e-12;
    flow
}

#[test]
fn latent_plateau_updates_drag_and_accelerating_flow_without_temperature_change() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let curve = phase_curve(&card, [1.0; 2]);
    let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.1).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let flow_config = cavity();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, &curve, 0.0);
    let resolve = |_: usize, _: &DepositedCell, t: f64| {
        liquid.at_temperature(Temperature::new(t)).map_err(|error| {
            fs_material::MaterialError::State {
                what: error.to_string(),
            }
        })
    };
    let hot = MaterialFlowFeedback::try_new(
        &old,
        &domain,
        input(&ids, &curve),
        &flow_config,
        Heater::default(),
        &gate,
        100,
        resolve,
    )
    .unwrap()
    .with_mushy_mobility(&mobility)
    .unwrap()
    .with_unsteady_momentum()
    .unwrap();
    let cold = MaterialFlowFeedback::try_new(
        &old,
        &domain,
        input(&ids, &curve),
        &flow_config,
        (),
        &gate,
        100,
        resolve,
    )
    .unwrap()
    .with_mushy_mobility(&mobility)
    .unwrap()
    .with_unsteady_momentum()
    .unwrap();
    let (hot1, hot_history1, _) = old
        .advance_observed(&hot, 0.125, config(), || Ok(()))
        .unwrap();
    let (cold1, cold_history1, _) = old
        .advance_observed(&cold, 0.125, config(), || Ok(()))
        .unwrap();
    assert_eq!(hot_history1.unsteady_flow(), cold_history1.unsteady_flow());
    let (hot2, hot_history2, audit) = hot1
        .advance_observed(&hot_history1, 0.225, config(), || Ok(()))
        .unwrap();
    let (_, cold_history2, _) = cold1
        .advance_observed(&cold_history1, 0.225, config(), || Ok(()))
        .unwrap();
    let applied = hot_history2.applied().unwrap();
    assert_eq!(applied.time_s, 0.125);
    assert!(
        hot_history2
            .endpoint()
            .iter()
            .all(|sample| sample.temperature_k == 500.0)
    );
    let phi = 62.5_f64 / 100.0;
    near(
        applied.mushy_drag_pa_s_m2.as_ref().unwrap()[1],
        (1.0 - phi).powi(2) / phi.powi(3),
    );
    assert!(applied.mushy_drag_pa_s_m2.as_ref().unwrap()[1] < 2.0);
    assert!(
        hot_history2.unsteady_step().unwrap().kinetic_energy_j
            > cold_history2.unsteady_step().unwrap().kinetic_energy_j
    );
    assert!(
        hot_history2.unsteady_step().unwrap().kinetic_energy_j
            > hot_history1.unsteady_step().unwrap().kinetic_energy_j
    );
    assert_eq!(&applied.flow, hot_history1.unsteady_flow().unwrap());
    near(hot2.energy_j() - old.energy_j(), 22.5);
    near(audit.external_heat_j, 10.0);
    near(hot2.mass_kg(), old.mass_kg());
    near(hot_history2.audit().external_enthalpy_j, 0.0);
    assert_eq!(hot_history2.unsteady_time_s(), Some(hot2.time_s()));
}

#[test]
fn endpoint_loss_of_mushy_support_rolls_back_heat_and_native_momentum() {
    let card = card(claims(1.0, 1.0, false));
    let liquid = liquid(&card);
    let curve = phase_curve(&card, [1.0; 2]);
    let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.2).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let flow_config = cavity();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, &curve, -400.0);
    let observer = MaterialFlowFeedback::try_new(
        &old,
        &domain,
        input(&ids, &curve),
        &flow_config,
        (),
        &gate,
        100,
        |_, _: &DepositedCell, t| {
            liquid.at_temperature(Temperature::new(t)).map_err(|error| {
                fs_material::MaterialError::State {
                    what: error.to_string(),
                }
            })
        },
    )
    .unwrap()
    .with_mushy_mobility(&mobility)
    .unwrap()
    .with_unsteady_momentum()
    .unwrap();
    assert!(
        old.advance_observed(&observer, 0.1, config(), || Ok(()))
            .is_err()
    );
    assert_eq!(old.cells()[&CellId(1)].specific_h_j_kg, 50.0);
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert!(observer.unsteady_step().is_none());
    assert!(
        observer
            .unsteady_flow()
            .unwrap()
            .velocity_m_s
            .iter()
            .all(|v| *v == [0.0; 3])
    );
    assert!(observer.applied().is_none());
}

#[test]
fn exact_solid_phase_occupancy_refuses_retained_momentum_without_breaking_quasisteady_flow() {
    use super::continuity::{chart, resolve, source};
    use fs_lbm::conjugate::enthalpy::advection::feedback::MaterialFlowError;

    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let flow_config = cavity();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let old = initial(&domain, &curve, 0.0);
    let phase = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        input(&ids, &curve),
        &flow_config,
        (),
        &gate,
        100,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap();
    // Every cell currently flows, but this owner can change native occupancy
    // at a later exact-solid endpoint. It cannot silently discard stored motion.
    assert!(
        phase
            .endpoint()
            .iter()
            .all(|s| s.transport.mobile_transport().is_some())
    );
    assert!(matches!(phase.clone().with_unsteady_momentum(),
        Err(MaterialFlowError::Thermal(DepositionError::InvalidInput(reason)))
            if reason.contains("phase occupancy") && reason.contains("kinetic-energy handoff")));
    assert!(phase.unsteady_flow().is_none());
    assert_eq!(old.time_s(), 0.0);
    assert!(phase.applied().is_none());
    // The same phase/source history still runs through the existing native
    // quasisteady owner, without manufacturing a retained liquid reference.
    let (next, accepted, audit) = old
        .advance_observed(&phase, 0.1, config(), || Ok(()))
        .unwrap();
    assert!(
        accepted
            .applied()
            .unwrap()
            .flow
            .velocity_m_s
            .iter()
            .flatten()
            .any(|u| u.abs() > 0.0)
    );
    assert!(accepted.momentum_solves() > 0);
    assert!(accepted.unsteady_time_s().is_none());
    near(next.energy_j(), old.energy_j());
    near(audit.external_heat_j, 0.0);
    near(next.mass_kg(), old.mass_kg());
}
