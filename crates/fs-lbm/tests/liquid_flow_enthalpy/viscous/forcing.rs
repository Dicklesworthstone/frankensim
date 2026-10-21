//! Prescribed body-force work enters the original momentum and caloric owners.
use super::*;
use fs_lbm::conjugate::creeping_flow_with_body_force_heat;

fn forces(f: f64) -> Vec<[f64; 3]> {
    vec![[-f, f, 0.0], [-f, -f, 0.0], [f, f, 0.0], [f, -f, 0.0]]
}
fn stationary() -> SimpleConfig {
    let mut c = driven_box();
    c.faces[3] = FvBoundary::wall();
    c
}

#[test]
fn native_body_work_closes_with_walls_darcy_and_an_absorbing_actuator() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let gate = CancelGate::new();
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 100.0,
        conductivity_w_m_k: 0.01,
        kinematic_viscosity_m2_s: 100.0,
    };
    for (wall, f, beta) in [
        (0.0, 100.0, 0.0),
        (0.0, 100.0, 200.0),
        (1.0, 100.0, 200.0),
        (1.0, -25.0, 0.0),
    ] {
        let mut config = driven_box();
        config.faces[3] = FvBoundary::Wall {
            velocity: [wall, 0.0, 0.0],
        };
        let (flow, p) = creeping_flow_with_body_force_heat(
            &domain,
            &fluid,
            &[100.0; 4],
            &forces(f),
            Some(&[beta; 4]),
            &config,
            controls(),
            &gate,
        )
        .unwrap();
        // Independent one-coordinate MAC solution: four fluxes are +/-a,
        // D=mu*(24a²-4Ua+2U²)+4 beta a²; actuator work is 2mu*U(U-a)+4F*a.
        let a = (200.0 * wall + 4.0 * f) / (2400.0 + 4.0 * beta);
        near(flow.field.flux_x(1, 1, 0), a);
        near(flow.field.flux_x(1, 0, 0), -a);
        near(p.supplied_body_power_w, 4.0 * f * a);
        near(p.supplied_wall_power_w, 200.0 * wall * (wall - a));
        near(
            p.molecular_heat_w,
            100.0 * (24.0 * a * a - 4.0 * wall * a + 2.0 * wall * wall),
        );
        near(p.darcy_heat_w, 4.0 * beta * a * a);
        near(p.heat_w, p.supplied_wall_power_w + p.supplied_body_power_w);
        for power in p.cell_body_power_w.unwrap() {
            near(power, f * a);
        }
        near(p.pressure_constraint_power_w, 0.0);
        near(p.residual_w, 0.0);
        if f < 0.0 {
            assert!(p.supplied_body_power_w < 0.0);
            assert!(p.heat_w >= 0.0);
        }
    }
    // A conservative constant force is balanced by pressure in a closed box.
    // It must not produce an invented f² heat source.
    let (flow, p) = creeping_flow_with_body_force_heat(
        &domain,
        &fluid,
        &[100.0; 4],
        &[[100.0, 0.0, 0.0]; 4],
        None,
        &stationary(),
        controls(),
        &gate,
    )
    .unwrap();
    for v in flow.velocity_m_s.iter().flatten() {
        near(*v, 0.0);
    }
    near(p.heat_w, 0.0);
    near(p.supplied_body_power_w, 0.0);
}

#[test]
fn changing_material_changes_forced_motion_heating_and_the_next_material_epoch() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = stationary();
    let gate = CancelGate::new();
    let ids = [
        Some(CellId(0)),
        Some(CellId(1)),
        Some(CellId(2)),
        Some(CellId(3)),
    ];
    let mut results = Vec::new();
    for (mu, thinning) in [(100.0, false), (200.0, false), (100.0, true)] {
        let card = source(mu, thinning);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let observer = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 10)
            .unwrap()
            .with_body_force_heating(forces(100.0), controls())
            .unwrap();
        let (warm, first, first_thermal) = old
            .advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
            .unwrap();
        let a = 100.0 / (6.0 * mu);
        let h = 0.5 * 100.0 * a;
        let temperature = 300.0 + h / 100.0;
        near(first.applied().unwrap().flow.field.flux_x(1, 1, 0), a);
        assert_eq!(
            first.applied().unwrap().external_force_n_m3.as_deref(),
            Some(forces(100.0).as_slice())
        );
        for c in 0..4 {
            near(warm.temperature_k(CellId(c)).unwrap(), temperature);
        }
        let work = first.viscous_heating_audit().unwrap();
        near(work.supplied_wall_work_j, 0.0);
        near(work.supplied_body_work_j, 4.0 * h);
        near(warm.energy_j() - old.energy_j(), work.supplied_body_work_j);
        near(first_thermal.external_heat_j, work.generated_heat_j);
        let (next, second, audit) = warm
            .advance_observed(&first, 1.0, thermal_step(), || Ok(()))
            .unwrap();
        let mu1 = if thinning {
            mu * (1.0 - 0.09 * (temperature - 300.0))
        } else {
            mu
        };
        let a1 = 100.0 / (6.0 * mu1);
        let applied = second.applied().unwrap();
        near(applied.flow.field.flux_x(1, 1, 0), a1);
        for sample in &applied.samples {
            near(sample.transport.dynamic_viscosity_pa_s(), mu1);
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
        }
        let work = second.viscous_heating_audit().unwrap();
        near(work.supplied_body_work_j, 0.5 * 400.0 * a1);
        near(next.energy_j() - warm.energy_j(), work.generated_heat_j);
        near(audit.external_heat_j, work.generated_heat_j);
        near(work.residual_j, 0.0);
        near(next.mass_kg(), old.mass_kg());
        results.push((a1, next.energy_j()));
    }
    near(results[1].0, 0.5 * results[0].0);
    near(results[1].1, 0.5 * results[0].1);
    assert!(results[2].0 > results[0].0);
    assert!(results[2].1 > results[0].1);
}

#[test]
fn rejected_forced_candidates_leave_the_heat_work_and_source_epoch_untouched() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = stationary();
    let gate = CancelGate::new();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let ids = [
        Some(CellId(0)),
        Some(CellId(1)),
        Some(CellId(2)),
        Some(CellId(3)),
    ];
    let base = input(&ids, &liquid, &caloric)
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
        .unwrap();
    assert!(
        base.clone()
            .with_body_force_heating(vec![], controls())
            .is_err()
    );
    assert!(
        base.clone()
            .with_body_force_heating(vec![[f64::NAN; 3]; 4], controls())
            .is_err()
    );
    let observer = base
        .with_body_force_heating(forces(100.0), controls())
        .unwrap();
    assert!(
        old.advance_observed(
            &observer,
            1.0,
            StepConfig {
                max_dt_s: 0.5,
                ..thermal_step()
            },
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(observer.inner().steps, 0);
    assert!(observer.applied().is_none());
    assert_eq!(
        observer.viscous_heating_audit(),
        Some(ViscousHeatingAudit::default())
    );
    near(old.energy_j(), 0.0);
    near(old.time_s(), 0.0);
    let (retry, accepted, _) = old
        .advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
        .unwrap();
    near(retry.energy_j(), 100.0 / 3.0);
    near(
        accepted
            .viscous_heating_audit()
            .unwrap()
            .supplied_body_work_j,
        retry.energy_j(),
    );
    gate.request();
    assert!(
        old.advance_observed(&observer, 0.5, thermal_step(), || Ok(()))
            .is_err()
    );
    assert!(observer.applied().is_none());
}
