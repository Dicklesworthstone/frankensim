//! G1/G3/G4: physical transient loss heats original h; numerical loss does not.
use super::*;
use fs_lbm::conjugate::UnsteadyHeatingAudit;

fn forces(f: f64) -> Vec<[f64; 3]> {
    vec![[-f, f, 0.0], [-f, -f, 0.0], [f, f, 0.0], [f, -f, 0.0]]
}

fn amplitude(old: f64, dt: f64, mu: f64, beta: f64, wall: f64, force: f64) -> f64 {
    // Independently reduced 2x2 MAC system: four unit-volume face masses.
    (4.0 * old / dt + 2.0 * mu * wall + 4.0 * force) / (4.0 / dt + 24.0 * mu + 4.0 * beta)
}

#[test]
fn native_startup_separates_physical_heat_kinetic_storage_and_be_loss() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let old = initial(&domain, 1.0, &caloric);
    for (wall, force) in [(1.0, 0.0), (0.0, 100.0), (1.0, -25.0)] {
        let mut config = driven_box();
        config.faces[3] = FvBoundary::Wall {
            velocity: [wall, 0.0, 0.0],
        };
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_unsteady_body_force_heating(forces(force), controls())
            .unwrap();
        let dt = 0.001;
        let (next, accepted, thermal) = old
            .advance_observed(&observer, dt, thermal_step(), || Ok(()))
            .unwrap();
        let a = amplitude(0.0, dt, 100.0, 0.0, wall, force);
        let step = accepted.unsteady_mechanical_step().unwrap();
        let power = &step.power;
        near(accepted.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a);
        near(accepted.unsteady_flow().unwrap().field.flux_x(1, 0, 0), -a);
        near(
            power.molecular_heat_w,
            100.0 * (24.0 * a * a - 4.0 * wall * a + 2.0 * wall * wall),
        );
        near(power.supplied_wall_power_w, 200.0 * wall * (wall - a));
        near(power.supplied_body_power_w, 4.0 * force * a);
        near(step.initial_kinetic_energy_j, 0.0);
        near(step.final_kinetic_energy_j, 2.0 * a * a);
        near(step.backward_euler_loss_j, 2.0 * a * a);
        assert!(step.backward_euler_loss_j > 0.0);
        near(step.residual_j, 0.0);
        near(step.duration_s, dt);
        assert_eq!(step.end_s, next.time_s());
        for c in 0..4 {
            let heat = 100.0
                * (6.0 * a * a
                    + if c >= 2 {
                        wall * wall - 2.0 * wall * a
                    } else {
                        0.0
                    });
            near(power.cell_heat_w[c], heat);
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, dt * heat);
        }
        let audit = accepted.unsteady_heating_audit().unwrap();
        near(audit.generated_heat_j, dt * power.heat_w);
        near(thermal.external_heat_j, audit.generated_heat_j);
        near(next.energy_j() - old.energy_j(), audit.generated_heat_j);
        near(
            audit.supplied_wall_work_j + audit.supplied_body_work_j,
            audit.generated_heat_j + audit.final_kinetic_energy_j + audit.backward_euler_loss_j,
        );
        assert_eq!(audit.substeps, 1);
        assert_eq!(accepted.momentum_solves(), 1);
        assert!(
            accepted
                .applied()
                .unwrap()
                .flow
                .velocity_m_s
                .iter()
                .all(|v| *v == [0.0; 3])
        );
        assert!(accepted.applied().unwrap().viscous_power.is_none());
        near(accepted.audit().external_enthalpy_j, 0.0);
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(
            observer.unsteady_heating_audit(),
            Some(UnsteadyHeatingAudit::default())
        );
        if force < 0.0 {
            assert!(audit.supplied_body_work_j < 0.0);
        }
    }
}

#[test]
fn actuator_shutoff_preserves_motion_and_does_not_thermalize_numerical_loss() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let old = initial(&domain, 1.0, &caloric);
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let (driven, history, _) = old
        .advance_observed(&observer, 0.001, thermal_step(), || Ok(()))
        .unwrap();
    let a0 = amplitude(0.0, 0.001, 100.0, 0.0, 0.0, 100.0);
    let old_flow = history.unsteady_flow().unwrap().clone();
    let coast = history
        .with_unsteady_body_force_heating(forces(0.0), controls())
        .unwrap();
    assert_eq!(coast.unsteady_flow(), Some(&old_flow));
    let (next, accepted, thermal) = driven
        .advance_observed(&coast, 0.002, thermal_step(), || Ok(()))
        .unwrap();
    let a1 = amplitude(a0, 0.001, 100.0, 0.0, 0.0, 0.0);
    let audit = accepted.unsteady_heating_audit().unwrap();
    near(accepted.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a1);
    assert!(a1 > 0.0 && a1 < a0);
    near(audit.initial_kinetic_energy_j, 2.0 * a0 * a0);
    near(audit.final_kinetic_energy_j, 2.0 * a1 * a1);
    near(audit.backward_euler_loss_j, 2.0 * (a1 - a0).powi(2));
    near(audit.generated_heat_j, 0.001 * 2400.0 * a1 * a1);
    near(audit.supplied_wall_work_j, 0.0);
    near(audit.supplied_body_work_j, 0.0);
    near(
        audit.initial_kinetic_energy_j - audit.final_kinetic_energy_j,
        audit.generated_heat_j + audit.backward_euler_loss_j,
    );
    assert!(audit.generated_heat_j < audit.initial_kinetic_energy_j - audit.final_kinetic_energy_j);
    near(next.energy_j() - driven.energy_j(), audit.generated_heat_j);
    near(thermal.external_heat_j, audit.generated_heat_j);
    assert_eq!(&accepted.applied().unwrap().flow, &old_flow);
}

#[test]
fn dissipated_heat_changes_source_viscosity_and_subsequent_native_motion() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let mut outcomes = Vec::new();
    for (mu, thinning) in [(100.0, false), (200.0, false), (100.0, true)] {
        let card = source(mu, thinning);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_unsteady_body_force_heating(forces(100.0), controls())
            .unwrap();
        let (warm, first, _) = old
            .advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
            .unwrap();
        let a0 = amplitude(0.0, 0.5, mu, 0.0, 0.0, 100.0);
        let heat0 = 0.5 * 24.0 * mu * a0 * a0;
        near(warm.energy_j() - old.energy_j(), heat0);
        let temperature = 300.0 + heat0 / (4.0 * 100.0);
        let mu1 = if thinning {
            mu * (1.0 - 0.09 * (temperature - 300.0))
        } else {
            mu
        };
        let (next, second, thermal) = warm
            .advance_observed(&first, 1.0, thermal_step(), || Ok(()))
            .unwrap();
        let a1 = amplitude(a0, 0.5, mu1, 0.0, 0.0, 100.0);
        near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a1);
        for sample in &second.applied().unwrap().samples {
            near(sample.temperature_k, temperature);
            near(sample.transport.dynamic_viscosity_pa_s(), mu1);
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
        }
        near(thermal.external_heat_j, 0.5 * 24.0 * mu1 * a1 * a1);
        near(next.energy_j() - warm.energy_j(), thermal.external_heat_j);
        near(next.mass_kg(), old.mass_kg());
        outcomes.push((a1, next.energy_j()));
    }
    assert!(outcomes[1].0 < outcomes[0].0 && outcomes[1].1 < outcomes[0].1);
    assert!(outcomes[2].0 > outcomes[0].0 && outcomes[2].1 > outcomes[0].1);
}

#[derive(Clone, Default)]
struct Intervals(Vec<f64>);
impl DepositionObserver for Intervals {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        dt: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.0.push(dt);
        Ok(())
    }
}

#[test]
fn cfl_shortened_intervals_feed_one_native_solve_and_one_heat_deposit_each() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let card = source(1000.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let observer = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, Intervals::default(), &gate, 20)
        .unwrap()
        .with_unsteady_body_force_heating(forces(1000.0), controls())
        .unwrap();
    let (next, accepted, thermal) = old
        .advance_observed(&observer, 20.0, thermal_step(), || Ok(()))
        .unwrap();
    let intervals = &accepted.inner().0;
    assert!(intervals.len() > 2);
    assert_eq!(intervals[0], 10.0);
    assert!(
        intervals[1..intervals.len() - 1]
            .iter()
            .any(|dt| *dt < 10.0)
    );
    let (mut a, mut heat, mut work, mut numerical) = (0.0, 0.0, 0.0, 0.0);
    for &dt in intervals {
        let next_a = amplitude(a, dt, 1000.0, 0.0, 0.0, 1000.0);
        heat += dt * 24000.0 * next_a * next_a;
        work += dt * 4000.0 * next_a;
        numerical += 2.0 * (next_a - a).powi(2);
        a = next_a;
    }
    let audit = accepted.unsteady_heating_audit().unwrap();
    near(audit.duration_s, 20.0);
    near(audit.generated_heat_j, heat);
    near(audit.supplied_body_work_j, work);
    near(audit.backward_euler_loss_j, numerical);
    near(audit.final_kinetic_energy_j, 2.0 * a * a);
    near(next.energy_j() - old.energy_j(), heat);
    near(thermal.external_heat_j, heat);
    assert_eq!(accepted.momentum_solves(), intervals.len());
    assert_eq!(audit.substeps, intervals.len());
    assert_eq!(accepted.unsteady_time_s(), Some(20.0));
}

#[test]
fn late_refusal_rolls_back_native_heat_work_history_and_instantaneous_query_refuses() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let observer = input(&ids, &liquid, &caloric)
        .feedback(
            &old,
            &domain,
            &config,
            History {
                steps: 0,
                reject: Some(2),
            },
            &gate,
            10,
        )
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
    let flow = observer.unsteady_flow().unwrap().clone();
    assert!(observer.requires_step_interval());
    assert!(old.power_balance_observed(&observer, || Ok(())).is_err());
    assert!(
        old.advance_observed(
            &observer,
            0.002,
            StepConfig {
                max_dt_s: 0.001,
                ..thermal_step()
            },
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(observer.unsteady_flow(), Some(&flow));
    assert_eq!(observer.unsteady_time_s(), Some(0.0));
    assert!(observer.unsteady_mechanical_step().is_none());
    assert_eq!(
        observer.unsteady_heating_audit(),
        Some(UnsteadyHeatingAudit::default())
    );
    assert_eq!(observer.inner().steps, 0);
    assert!(observer.applied().is_none());
    near(old.energy_j(), 0.0);
    let (next, accepted, thermal) = old
        .advance_observed(&observer, 0.001, thermal_step(), || Ok(()))
        .unwrap();
    let a = amplitude(0.0, 0.001, 100.0, 0.0, 0.0, 100.0);
    near(next.energy_j(), 0.001 * 2400.0 * a * a);
    near(thermal.external_heat_j, next.energy_j());
    assert_eq!(accepted.inner().steps, 1);
    assert_eq!(accepted.unsteady_heating_audit().unwrap().substeps, 1);
}
