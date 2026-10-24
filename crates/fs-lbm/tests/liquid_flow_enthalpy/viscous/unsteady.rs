//! G1/G3/G4: physical transient loss heats original h; numerical loss does not.
use super::*;
use fs_lbm::conjugate::UnsteadyHeatingAudit;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;

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

fn expanding_source(expansion: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut density = scalar(DENSITY_PROPERTY, Density::DIMS, 1.0);
    density.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(300.0, 1.0), (810.0, 1.0 - 510.0 * expansion)],
        dims: Density::DIMS,
    };
    density.interpolation = InterpolationPolicy::LinearInside;
    claims.insert_claim(density).unwrap();
    for (name, dims, value) in [
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            100.0,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            0.01,
        ),
        (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 100.0),
    ] {
        claims.insert_claim(scalar(name, dims, value)).unwrap();
    }
    card(claims)
}

fn hot_column<'a>(
    domain: &VoxelDomain,
    caloric: &'a ResolvedHeatCapacityCurve,
) -> DepositionSimulation<'a> {
    let seed = initial(domain, 1.0, caloric);
    let events = seed
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            // Left column at 310 K, right at 300 K; no independent heating source.
            cell.specific_h_j_kg = if cell.id.0 % 2 == 0 { 1000.0 } else { 0.0 };
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MAT, caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 0,
            max_events: 4,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_step(), || Ok(()))
    .unwrap()
    .0
}

fn material_gravity() -> BoussinesqConfig<3> {
    // Manufactured acceleration keeps work well above roundoff while density
    // departures remain at most 2%. This is no gravitational-potential model.
    BoussinesqConfig {
        gravity_m_s2: [0.0, -10000.0, 0.0],
        max_relative_density_departure: 0.02,
        max_relative_specific_heat_departure: 0.0,
    }
}

#[test]
fn source_density_changes_retained_motion_and_only_physical_loss_heats_enthalpy() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let dt = 0.1;
    let mut outcomes = Vec::new();
    for expansion in [0.0, 0.0005, 0.001] {
        let card = expanding_source(expansion);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = hot_column(&domain, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback_boussinesq(&old, &domain, &config, (), &gate, 10, material_gravity())
            .unwrap()
            .with_unsteady_viscous_heating(controls())
            .unwrap();
        let (warm, first, thermal) = old
            .advance_observed(&observer, dt, thermal_step(), || Ok(()))
            .unwrap();
        // Two hot-column forces B act on cell velocities a/2, giving B*a
        // total power. Four native face masses give inertia 4*da/dt.
        let b = 100000.0 * expansion;
        let a0 = amplitude(0.0, dt, 100.0, 0.0, 0.0, b / 4.0);
        let heat0 = dt * 2400.0 * a0 * a0;
        let step = first.unsteady_mechanical_step().unwrap();
        near(first.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a0);
        near(step.supplied_buoyancy_power_w, b * a0);
        near(step.power.supplied_body_power_w, b * a0);
        near(step.power.molecular_heat_w, 2400.0 * a0 * a0);
        near(step.power.darcy_heat_w, 0.0);
        near(step.final_kinetic_energy_j, 2.0 * a0 * a0);
        near(step.backward_euler_loss_j, 2.0 * a0 * a0);
        near(step.power.pressure_constraint_power_w, 0.0);
        near(dt * b * a0, heat0 + 4.0 * a0 * a0);
        near(warm.energy_j() - old.energy_j(), heat0);
        near(thermal.external_heat_j, heat0);
        for c in 0..4 {
            near(
                step.cell_buoyancy_power_w.as_ref().unwrap()[c],
                if c % 2 == 0 { b * a0 / 2.0 } else { 0.0 },
            );
            near(
                warm.cells()[&CellId(c as u64)].specific_h_j_kg
                    - old.cells()[&CellId(c as u64)].specific_h_j_kg,
                heat0 / 4.0,
            );
        }
        let (next, second, thermal) = warm
            .advance_observed(&first, 2.0 * dt, thermal_step(), || Ok(()))
            .unwrap();
        // The first physical heat is uniform: it changes source density but
        // its uniform force offset is hydrostatic, leaving the same B drive.
        let a1 = amplitude(a0, dt, 100.0, 0.0, 0.0, b / 4.0);
        let heat1 = dt * 2400.0 * a1 * a1;
        let audit = second.unsteady_heating_audit().unwrap();
        near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a1);
        near(audit.initial_kinetic_energy_j, 2.0 * a0 * a0);
        near(audit.final_kinetic_energy_j, 2.0 * a1 * a1);
        near(audit.backward_euler_loss_j, 2.0 * (a1 - a0).powi(2));
        near(audit.generated_heat_j, heat1);
        near(audit.supplied_buoyancy_work_j, dt * b * a1);
        near(audit.supplied_body_work_j, audit.supplied_buoyancy_work_j);
        near(
            dt * b * a1,
            heat1 + 2.0 * (a1 * a1 - a0 * a0) + 2.0 * (a1 - a0).powi(2),
        );
        near(next.energy_j() - warm.energy_j(), heat1);
        near(thermal.external_heat_j, heat1);
        near(second.audit().external_enthalpy_j, 0.0);
        near(next.mass_kg(), old.mass_kg());
        assert_eq!(
            &second.applied().unwrap().flow,
            first.unsteady_flow().unwrap()
        );
        near(second.applied().unwrap().time_s, dt);
        for sample in &second.applied().unwrap().samples {
            near(
                sample.temperature_k,
                warm.temperature_k(sample.cell).unwrap(),
            );
            near(
                sample.material.resolved_transport().density_kg_m3(),
                1.0 - expansion * (sample.temperature_k - 300.0),
            );
            near(sample.transport.density_kg_m3(), 1.0);
            assert_eq!(
                sample.material.material().material().card_identity(),
                card.content_hash()
            );
        }
        if expansion > 0.0 {
            assert!(a1 > a0 && a0 > 0.0);
            assert!(audit.backward_euler_loss_j > 0.0);
            assert!(audit.generated_heat_j < audit.supplied_body_work_j);
        }
        assert_eq!(
            observer.unsteady_heating_audit(),
            Some(UnsteadyHeatingAudit::default())
        );
        outcomes.push((a1, heat1));
    }
    near(outcomes[0].0, 0.0);
    near(outcomes[0].1, 0.0);
    near(outcomes[2].0, 2.0 * outcomes[1].0);
    near(outcomes[2].1, 4.0 * outcomes[1].1);
}

#[test]
fn simultaneous_buoyancy_and_absorbing_actuators_match_the_same_prescribed_total_force() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    let wall = 0.01;
    config.faces[3] = FvBoundary::Wall {
        velocity: [wall, 0.0, 0.0],
    };
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let card = expanding_source(0.001);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = hot_column(&domain, &caloric);
    let force = -10.0;
    let observer = input(&ids, &liquid, &caloric)
        .feedback_boussinesq(&old, &domain, &config, (), &gate, 10, material_gravity())
        .unwrap()
        .with_unsteady_body_force_heating(forces(force), controls())
        .unwrap();
    let dt = 0.1;
    let (next, accepted, thermal) = old
        .advance_observed(&observer, dt, thermal_step(), || Ok(()))
        .unwrap();
    let b = 100.0;
    let a = amplitude(0.0, dt, 100.0, 0.0, wall, force + b / 4.0);
    let heat = dt * 100.0 * (24.0 * a * a - 4.0 * wall * a + 2.0 * wall * wall);
    let wall_work = dt * 200.0 * wall * (wall - a);
    let external_work = dt * 4.0 * force * a;
    let buoyancy_work = dt * b * a;
    let step = accepted.unsteady_mechanical_step().unwrap();
    let audit = accepted.unsteady_heating_audit().unwrap();
    near(accepted.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a);
    near(audit.supplied_buoyancy_work_j, buoyancy_work);
    near(
        audit.supplied_body_work_j - audit.supplied_buoyancy_work_j,
        external_work,
    );
    near(audit.supplied_wall_work_j, wall_work);
    assert!(wall_work < 0.0 && external_work < 0.0 && buoyancy_work > 0.0);
    near(audit.generated_heat_j, heat);
    near(audit.final_kinetic_energy_j, 2.0 * a * a);
    near(audit.backward_euler_loss_j, 2.0 * a * a);
    near(
        wall_work + external_work + buoyancy_work,
        heat + 4.0 * a * a,
    );
    near(next.energy_j() - old.energy_j(), heat);
    near(thermal.external_heat_j, heat);
    for c in 0..4 {
        let buoyant = if c % 2 == 0 { b * a / 2.0 } else { 0.0 };
        near(step.cell_buoyancy_power_w.as_ref().unwrap()[c], buoyant);
        near(
            step.power.cell_body_power_w.as_ref().unwrap()[c] - buoyant,
            force * a,
        );
    }
    assert_eq!(
        accepted.applied().unwrap().external_force_n_m3.as_deref(),
        Some(forces(force).as_slice())
    );

    // A constant-density source with the same mu/cp/k and total authored force
    // must follow the identical momentum and heat equations, without buoyancy.
    let plain_card = source(100.0, false);
    let plain_liquid = super::liquid(&plain_card);
    let plain_caloric = liquid_caloric(&plain_card);
    let plain_old = hot_column(&domain, &plain_caloric);
    let total: Vec<_> = forces(force)
        .into_iter()
        .enumerate()
        .map(|(c, mut f)| {
            if c % 2 == 0 {
                f[1] += b;
            }
            f
        })
        .collect();
    let prescribed = input(&ids, &plain_liquid, &plain_caloric)
        .feedback(&plain_old, &domain, &config, (), &gate, 10)
        .unwrap()
        .with_unsteady_body_force_heating(total, controls())
        .unwrap();
    let (plain_next, plain, _) = plain_old
        .advance_observed(&prescribed, dt, thermal_step(), || Ok(()))
        .unwrap();
    near(plain.unsteady_flow().unwrap().field.flux_x(1, 1, 0), a);
    near(
        plain.unsteady_heating_audit().unwrap().generated_heat_j,
        heat,
    );
    near(
        plain
            .unsteady_heating_audit()
            .unwrap()
            .supplied_buoyancy_work_j,
        0.0,
    );
    assert!(
        plain
            .unsteady_mechanical_step()
            .unwrap()
            .cell_buoyancy_power_w
            .is_none()
    );
    for c in 0..4 {
        near(
            plain_next.cells()[&CellId(c)].specific_h_j_kg,
            next.cells()[&CellId(c)].specific_h_j_kg,
        );
    }
}

#[test]
fn late_buoyant_heat_refusal_and_cancellation_preserve_accepted_motion_and_all_work() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut config = driven_box();
    config.faces[3] = FvBoundary::wall();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let gate = CancelGate::new();
    let card = expanding_source(0.001);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = hot_column(&domain, &caloric);
    let observer = input(&ids, &liquid, &caloric)
        .feedback_boussinesq(
            &old,
            &domain,
            &config,
            History {
                steps: 0,
                reject: Some(2),
            },
            &gate,
            10,
            material_gravity(),
        )
        .unwrap()
        .with_unsteady_body_force_heating(forces(-10.0), controls())
        .unwrap();
    let (warm, accepted, _) = old
        .advance_observed(&observer, 0.1, thermal_step(), || Ok(()))
        .unwrap();
    let flow = accepted.unsteady_flow().unwrap().clone();
    let mechanical = accepted.unsteady_mechanical_step().unwrap().clone();
    let audit = accepted.unsteady_heating_audit().unwrap();
    let applied = accepted.applied().unwrap().clone();
    let cells = warm.cells().clone();
    assert!(audit.generated_heat_j > 0.0 && audit.supplied_buoyancy_work_j > 0.0);
    assert!(matches!(
        warm.advance_observed(&accepted, 0.2, thermal_step(), || Ok(())),
        Err(MaterialFlowError::Inner(DepositionError::InvalidInput(_)))
    ));
    assert_eq!(accepted.unsteady_flow(), Some(&flow));
    assert_eq!(accepted.unsteady_mechanical_step(), Some(&mechanical));
    assert_eq!(accepted.unsteady_heating_audit(), Some(audit));
    assert_eq!(accepted.applied(), Some(&applied));
    assert_eq!(accepted.unsteady_time_s(), Some(0.1));
    assert_eq!(accepted.inner().steps, 1);
    assert_eq!(warm.cells(), &cells);
    near(warm.time_s(), 0.1);
    gate.request();
    assert!(
        warm.advance_observed(&accepted, 0.2, thermal_step(), || Ok(()))
            .is_err()
    );
    assert_eq!(accepted.unsteady_flow(), Some(&flow));
    assert_eq!(accepted.unsteady_mechanical_step(), Some(&mechanical));
    assert_eq!(accepted.unsteady_heating_audit(), Some(audit));
    assert_eq!(accepted.applied(), Some(&applied));
    assert_eq!(accepted.inner().steps, 1);
    assert_eq!(warm.cells(), &cells);
}
