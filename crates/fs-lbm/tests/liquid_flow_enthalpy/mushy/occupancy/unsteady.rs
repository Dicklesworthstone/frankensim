//! G1/G3/G4: exact phase occupancy retains native inertia without a capture model.
use super::*;
use fs_lbm::conjugate::ChtError;

#[derive(Clone)]
struct UniformCycle;
impl DepositionObserver for UniformCycle {
    type Error = DepositionError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        let watts = if thermal.time_s() < 0.2 {
            1000.0
        } else if thermal.time_s() < 0.3 {
            -2000.0
        } else {
            2000.0
        };
        Ok(thermal.cells().keys().map(|&cell| (cell, watts)).collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}

#[test]
fn cold_solid_thaws_freezes_at_rest_and_remelts_on_the_original_momentum_clock() {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let old = geometry(&domain, &curve, [-20.0; 9]);
    for heating in [false, true] {
        let mut observer = MaterialFlowFeedback::try_new_phase(
            &old,
            &domain,
            input(&ids, &curve),
            &config,
            UniformCycle,
            &gate,
            1,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_unsteady_momentum()
        .unwrap();
        if heating {
            observer = observer.with_unsteady_viscous_heating(controls()).unwrap();
        }
        assert!(
            observer
                .endpoint()
                .iter()
                .all(|s| s.transport.mobile_transport().is_none())
        );
        assert_eq!(observer.unsteady_time_s(), Some(0.0));
        let (cold, history, audit) = old
            .advance_observed(&observer, 0.01, step(0.01), || Ok(()))
            .unwrap();
        assert_eq!(history.momentum_solves(), 0);
        assert_eq!(history.unsteady_step().unwrap().inner_iterations, 0);
        assert_eq!(history.unsteady_time_s(), Some(cold.time_s()));
        near(audit.external_heat_j, 90.0);
        let (mushy, history, _) = cold
            .advance_observed(&history, 0.1, step(0.1), || Ok(()))
            .unwrap();
        assert_eq!(history.momentum_solves(), 0);
        for s in history.endpoint() {
            near(s.specific_enthalpy_j_kg, 80.0);
            assert!(s.transport.mobile_transport().is_some());
        }
        let (liquid, history, _) = mushy
            .advance_observed(&history, 0.2, step(0.1), || Ok(()))
            .unwrap();
        assert_eq!(history.momentum_solves(), 1);
        let (frozen, history, audit) = liquid
            .advance_observed(&history, 0.3, step(0.1), || Ok(()))
            .unwrap();
        assert_eq!(history.momentum_solves(), 1);
        assert!(
            history
                .endpoint()
                .iter()
                .all(|s| s.transport.mobile_transport().is_none())
        );
        for cell in frozen.cells().values() {
            near(cell.specific_h_j_kg, -20.0);
        }
        near(audit.external_heat_j, -1800.0);
        near(frozen.energy_j(), old.energy_j());
        let (remelted, history, audit) = frozen
            .advance_observed(&history, 0.4, step(0.11), || Ok(()))
            .unwrap();
        assert_eq!(history.momentum_solves(), 0);
        assert!(
            history
                .endpoint()
                .iter()
                .all(|s| s.transport.mobile_transport().is_some())
        );
        near(audit.external_heat_j, 1800.0);
        near(remelted.mass_kg(), old.mass_kg());
        assert_eq!(remelted.cells().len(), 9);
        assert_eq!(history.unsteady_time_s(), Some(remelted.time_s()));
        assert!(
            history
                .unsteady_flow()
                .unwrap()
                .velocity_m_s
                .iter()
                .flatten()
                .all(|v| *v == 0.0)
        );
        if heating {
            let mechanical = history.unsteady_heating_audit().unwrap();
            assert_eq!(mechanical.substeps, 1);
            assert_eq!(mechanical.generated_heat_j, 0.0);
            assert_eq!(mechanical.final_kinetic_energy_j, 0.0);
            assert_eq!(mechanical.backward_euler_loss_j, 0.0);
        }
        // Neither the material query nor the new dormant state advances inputs.
        assert_eq!(observer.unsteady_time_s(), Some(0.0));
        assert!(observer.unsteady_step().is_none());
        assert_eq!(old.cells()[&CellId(0)].specific_h_j_kg, -20.0);
    }
}

#[test]
fn melting_opens_zero_speed_faces_without_restarting_the_moving_mushy_ring() {
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let mut amplitudes = Vec::new();
    for mu in [100.0, 200.0] {
        let card = source(mu);
        let curve = chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let mut h = [50.0; 9];
        h[4] = -20.0;
        let old = geometry(&domain, &curve, h);
        let observer = MaterialFlowFeedback::try_new_phase(
            &old,
            &domain,
            input(&ids, &curve),
            &config,
            CellPower {
                cell: CellId(4),
                watts: 2000.0,
            },
            &gate,
            4,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_unsteady_body_force_heating(corner_force(), controls())
        .unwrap();
        let dt = 0.1;
        let (melted, first, thermal) = old
            .advance_observed(&observer, dt, step(dt), || Ok(()))
            .unwrap();
        // Independent one-coordinate MAC equation: eight unit face masses,
        // 48*mu molecular resistance and 8*beta stationary Darcy resistance.
        let beta = 2.0 * mu;
        let q = 800.0 / (8.0 / dt + 48.0 * mu + 8.0 * beta);
        let mechanical = first.unsteady_mechanical_step().unwrap();
        let flow = first.unsteady_flow().unwrap();
        for x in 1..3 {
            near(flow.field.flux_x(x, 0, 0), -q);
            near(flow.field.flux_x(x, 2, 0), q);
            assert_eq!(flow.field.flux_x(x, 1, 0), 0.0);
        }
        for y in 1..3 {
            near(flow.field.flux_y(0, y, 0), q);
            near(flow.field.flux_y(2, y, 0), -q);
            assert_eq!(flow.field.flux_y(1, y, 0), 0.0);
        }
        near(mechanical.initial_kinetic_energy_j, 0.0);
        near(mechanical.final_kinetic_energy_j, 4.0 * q * q);
        near(mechanical.backward_euler_loss_j, 4.0 * q * q);
        near(mechanical.power.molecular_heat_w, 48.0 * mu * q * q);
        near(mechanical.power.darcy_heat_w, 8.0 * beta * q * q);
        near(mechanical.power.supplied_body_power_w, 800.0 * q);
        near(mechanical.residual_j, 0.0);
        let physical_heat = dt * (48.0 * mu + 8.0 * beta) * q * q;
        near(thermal.external_heat_j, 200.0 + physical_heat);
        near(melted.energy_j() - old.energy_j(), thermal.external_heat_j);
        assert_eq!(mechanical.power.cell_heat_w[4], 0.0);
        near(melted.cells()[&CellId(4)].specific_h_j_kg, 180.0);
        assert!(first.endpoint()[4].transport.mobile_transport().is_some());
        let retained = flow.clone();
        let (next, second, _) = melted
            .advance_observed(&first, 0.11, step(0.01), || Ok(()))
            .unwrap();
        assert_eq!(&second.applied().unwrap().flow, &retained);
        assert_eq!(second.momentum_solves(), 1);
        // Adding zero-speed faces does not alter the original kinetic sum.
        assert_eq!(
            second
                .unsteady_mechanical_step()
                .unwrap()
                .initial_kinetic_energy_j
                .to_bits(),
            mechanical.final_kinetic_energy_j.to_bits(),
        );
        assert!(
            second
                .applied()
                .unwrap()
                .mushy_drag_pa_s_m2
                .as_ref()
                .unwrap()[0]
                < beta
        );
        let reopened = &second.unsteady_flow().unwrap().field;
        assert!(
            [
                reopened.flux_x(1, 1, 0),
                reopened.flux_x(2, 1, 0),
                reopened.flux_y(1, 1, 0),
                reopened.flux_y(1, 2, 0),
            ]
            .iter()
            .any(|q| q.abs() > 1e-9)
        );
        near(next.mass_kg(), old.mass_kg());
        amplitudes.push(q);
    }
    assert!(amplitudes[1] < amplitudes[0]);
}

#[derive(Clone)]
struct CoolingAfterStartup;
impl DepositionObserver for CoolingAfterStartup {
    type Error = DepositionError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(BTreeMap::from([(
            CellId(0),
            if thermal.time_s() < 0.01 {
                0.0
            } else {
                -5000.0
            },
        )]))
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}

#[test]
fn freezing_a_moving_face_refuses_without_losing_kinetic_heat_or_native_history() {
    freezing_a_moving_face_refuses(false);
}

#[test]
fn newly_solid_rheology_uses_zero_rate_before_native_freezing_handoff() {
    freezing_a_moving_face_refuses(true);
}

fn freezing_a_moving_face_refuses(rheology: bool) {
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 9, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let old = geometry(&domain, &curve, [150.0; 9]);
    // Liquid rate may be nonzero, but the original phase source must never be
    // asked for a flowing solid. The physical momentum handoff must still
    // refuse moving faces after this source admission succeeds.
    let rate_source = |_: usize, cell: &DepositedCell, _: f64, rate: f64| {
        if law.classify_at_specific_enthalpy(cell.specific_h_j_kg)?.mobile().is_none()
            && rate != 0.0
        {
            return Err(fs_material::MaterialError::State {
                what: "exact-solid source received a liquid strain rate".into(),
            });
        }
        resolve(&card, &law, cell.specific_h_j_kg)
    };
    let observer = MaterialFlowFeedback::try_new_phase(
        &old,
        &domain,
        input(&ids, &curve),
        &config,
        CoolingAfterStartup,
        &gate,
        4,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_body_force_heating(corner_force(), controls())
    .unwrap();
    let observer = if rheology {
        observer.with_native_rheology(&old, &rate_source).unwrap()
    } else {
        observer
    };
    let (moving, history, _) = old
        .advance_observed(&observer, 0.01, step(0.01), || Ok(()))
        .unwrap();
    let coast = history
        .with_unsteady_body_force_heating(vec![[0.0; 3]; 9], controls())
        .unwrap();
    let native = coast.unsteady_flow().unwrap().clone();
    let mechanical = coast.unsteady_mechanical_step().unwrap().clone();
    let ledger = coast.unsteady_heating_audit();
    assert!(mechanical.final_kinetic_energy_j > 0.0);
    if rheology {
        assert!(coast.endpoint_shear_rate_per_s().unwrap()[0] > 0.0);
    }
    let failed = moving.advance_observed(&coast, 0.05, step(0.04), || Ok(()));
    assert!(matches!(failed,
        Err(MaterialFlowError::Momentum(ChtError::InvalidInput { field: "simple.momentum_state", reason }))
        if reason.contains("freezing a moving native face") && reason.contains("kinetic-energy handoff")
    ));
    assert_eq!(coast.unsteady_flow(), Some(&native));
    assert_eq!(coast.unsteady_mechanical_step(), Some(&mechanical));
    assert_eq!(coast.unsteady_heating_audit(), ledger);
    assert_eq!(coast.unsteady_time_s(), Some(moving.time_s()));
    assert!(
        moving
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg >= 150.0)
    );
    // A shorter mobile-only retry retains the same initial kinetic state.
    let (next, accepted, _) = moving
        .advance_observed(&coast, 0.02, step(0.01), || Ok(()))
        .unwrap();
    assert_eq!(&accepted.applied().unwrap().flow, &native);
    assert_eq!(
        accepted
            .unsteady_mechanical_step()
            .unwrap()
            .initial_kinetic_energy_j
            .to_bits(),
        mechanical.final_kinetic_energy_j.to_bits()
    );
    assert_eq!(accepted.unsteady_time_s(), Some(next.time_s()));
    near(next.mass_kg(), old.mass_kg());
    assert_eq!(coast.unsteady_time_s(), Some(0.01));
}
