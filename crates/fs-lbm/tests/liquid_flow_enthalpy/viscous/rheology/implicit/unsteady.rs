//! G1/G3/G4: original inertia and material response share the thermal endpoint.
use super::*;
use fs_lbm::conjugate::UnsteadyHeatingAudit;
use fs_lbm::conjugate::enthalpy::advection::feedback::implicit::{
    ImplicitUnsteadyFlowControls, ImplicitUnsteadyFlowReport,
};

fn coupled_controls(rheology: bool) -> ImplicitUnsteadyFlowControls {
    ImplicitUnsteadyFlowControls {
        thermal: implicit_controls(),
        rheology: rheology.then(rheology_controls),
    }
}

fn circulation_root(h0: f64, q0: f64, dt: f64, drive: f64, st: f64, sg: f64) -> (f64, f64, f64) {
    // Independently reduce the four free MAC faces: M=4, damping=24mu.
    // Dot momentum by q: h-h0=dt*drive*q+q0*q-q² per unit-mass cell.
    // Source mu=100+st*h/100+sg*sqrt(6)*q closes a scalar cubic.
    let (mut lo, mut hi) = (0.0, q0.max(drive / 600.0));
    for _ in 0..80 {
        let q = 0.5 * (lo + hi);
        let h = h0 + dt * drive * q + q0 * q - q * q;
        let mu = 100.0 + st * h / 100.0 + sg * 6.0_f64.sqrt() * q;
        if (4.0 / dt + 24.0 * mu) * q > 4.0 * drive + 4.0 * q0 / dt {
            hi = q;
        } else {
            lo = q;
        }
    }
    let q = 0.5 * (lo + hi);
    let h = h0 + dt * drive * q + q0 * q - q * q;
    (q, h, 100.0 + st * h / 100.0 + sg * 6.0_f64.sqrt() * q)
}

fn check_energy(
    report: &ImplicitUnsteadyFlowReport,
    old_q: f64,
    q: f64,
    heat: f64,
    drive: f64,
    dt: f64,
) {
    let m = report.mechanical;
    near(m.initial_kinetic_energy_j, 2.0 * old_q * old_q);
    near(m.final_kinetic_energy_j, 2.0 * q * q);
    near(m.backward_euler_loss_j, 2.0 * (q - old_q).powi(2));
    near(m.supplied_body_work_j, dt * 4.0 * drive * q);
    near(m.generated_heat_j, heat);
    near(report.deposited_heat_j, heat);
    near(report.thermal.stored_energy_change_j, heat);
    near(report.net_stream_enthalpy_j, 0.0);
    near(report.other_thermal_input_j, 0.0);
    assert!(report.combined_residual_j.abs() < 1e-6);
    assert!(report.nonlinear.converged);
    assert_eq!(m.substeps, 1);
}

#[test]
fn thermal_rate_feedback_changes_startup_and_retains_native_inertia_after_drive_off() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let mut first_speeds = Vec::new();
    for (st, sg) in [(0.0, 0.0), (40.0, 0.0), (40.0, 100.0)] {
        let card = rheological_card(st, sg, 2.0);
        let liquid = law(&card);
        let caloric = caloric(&card);
        let old = initial(&domain, 1.0, &caloric);
        let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
        let flow = source
            .feedback(&old, &domain, &config, (), &gate, 8192)
            .unwrap()
            .with_unsteady_body_force_heating(force(), controls())
            .unwrap();
        let dt = 2.0;
        let (q1, h1, mu1) = circulation_root(0.0, 0.0, dt, 100.0, st, sg);
        let (warm, first, report) = flow
            .step_unsteady_implicit(&old, dt, coupled_controls(true), &[1.0; 4], || Ok(()))
            .unwrap();
        let epoch = first.applied().unwrap();
        near(epoch.flow.field.flux_x(1, 1, 0), q1);
        assert_eq!(Some(&epoch.flow), first.unsteady_flow());
        assert_eq!(epoch.time_s, warm.time_s());
        for (i, sample) in epoch.samples.iter().enumerate() {
            near(sample.specific_enthalpy_j_kg, h1);
            near(sample.temperature_k, 300.0 + h1 / 100.0);
            near(sample.material.dynamic_viscosity_pa_s(), mu1);
            near(
                epoch.shear_rate_per_s.as_ref().unwrap()[i],
                6.0_f64.sqrt() * q1,
            );
            near(
                first.endpoint_shear_rate_per_s().unwrap()[i],
                6.0_f64.sqrt() * q1,
            );
            assert_eq!(
                sample.material.material().card_identity(),
                card.content_hash()
            );
        }
        check_energy(&report, 0.0, q1, 4.0 * h1, 100.0, dt);
        assert!(report.momentum_solves > report.rheology.unwrap().iterations);
        assert_eq!(first.momentum_solves(), report.momentum_solves);
        first_speeds.push(q1);

        let coast = first
            .with_unsteady_body_force_heating(vec![[0.0; 3]; 4], controls())
            .unwrap();
        // This short coast resolves substantial native inertia, not equilibrium.
        let dt2 = 0.001;
        let (q2, h2, mu2) = circulation_root(h1, q1, dt2, 0.0, st, sg);
        let (next, second, report2) = coast
            .step_unsteady_implicit(&warm, dt + dt2, coupled_controls(true), &[0.001; 4], || {
                Ok(())
            })
            .unwrap();
        near(second.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q2);
        assert!(q2 > 0.0 && q2 < q1);
        for sample in &second.applied().unwrap().samples {
            near(sample.material.dynamic_viscosity_pa_s(), mu2);
            near(next.cells()[&sample.cell].specific_h_j_kg, h2);
        }
        check_energy(&report2, q1, q2, 4.0 * (h2 - h1), 0.0, dt2);
        assert!(report2.mechanical.backward_euler_loss_j > 0.0);
        assert!(report2.deposited_heat_j < 2.0 * (q1 * q1 - q2 * q2));
        assert!(
            coast
                .step_unsteady_implicit(
                    &old,
                    dt + dt2,
                    coupled_controls(true),
                    &[0.0; 4],
                    || Ok(())
                )
                .is_err()
        );
        assert_eq!(flow.unsteady_time_s(), Some(0.0));
        assert_eq!(
            flow.unsteady_heating_audit(),
            Some(UnsteadyHeatingAudit::default())
        );
        assert_eq!(old.energy_j(), 0.0);
    }
    assert!(first_speeds[0] > first_speeds[1] && first_speeds[1] > first_speeds[2]);
}

#[test]
fn endpoint_velocity_drives_directed_enthalpy_and_the_same_material_face_conduction() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let cycle = [0usize, 2, 3, 1];
    for (mu, k) in [(100.0, 1.0), (200.0, 10.0)] {
        let original = source(mu, false);
        let mut claims = ClaimSet::new();
        for (_, value) in original.claims().claims_ordered() {
            let mut value = value.clone();
            if value.key.name() == THERMAL_CONDUCTIVITY_PROPERTY {
                value.value = PropertyValue::Scalar {
                    value: k,
                    dims: THERMAL_CONDUCTIVITY_DIMS,
                };
            }
            claims.insert_claim(value).unwrap();
        }
        let card = card(claims);
        let liquid = liquid(&card);
        let caloric = liquid_caloric(&card);
        let links: Vec<_> = (0..4)
            .map(|i| ThermalLink {
                id: 100 + i as u64,
                a: CellId(cycle[i] as u64),
                b: CellId(cycle[(i + 1) % 4] as u64),
                area_m2: 1.0,
                conductance_w_k: 999.0,
            })
            .collect();
        let old = thermal_state(
            &domain,
            caloric.curve(),
            &[200.0, 100.0, 100.0, 0.0],
            &links,
        );
        let flow = input(&ids, &liquid, &caloric)
            .feedback(&old, &domain, &config, (), &gate, 512)
            .unwrap()
            .with_unsteady_body_force_heating(force(), controls())
            .unwrap()
            .with_fluid_conduction(&old, &[100, 101, 102, 103])
            .unwrap();
        let dt = 8.0;
        let q = 400.0 / (4.0 / dt + 24.0 * mu);
        let lambda = dt * q;
        let eta = dt * k / 100.0;
        // Original C4 donor cycle plus Fourier links. Its first discrete mode
        // has L(cos)=d*cos-lambda*sin and L(sin)=d*sin+lambda*cos.
        let d = 1.0 + lambda + 2.0 * eta;
        let a = 100.0 * d / (d * d + lambda * lambda);
        let b = 100.0 * lambda / (d * d + lambda * lambda);
        let heat_per_cell = dt * 6.0 * mu * q * q;
        let expected = [
            100.0 + heat_per_cell + a,
            100.0 + heat_per_cell + b,
            100.0 + heat_per_cell - a,
            100.0 + heat_per_cell - b,
        ];
        let (next, accepted, report) = flow
            .step_unsteady_implicit(&old, dt, coupled_controls(false), &[1.0; 4], || Ok(()))
            .unwrap();
        near(accepted.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q);
        for (i, c) in cycle.into_iter().enumerate() {
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, expected[i]);
            assert_eq!(next.cells()[&CellId(c as u64)].boundary.source_w, 0.0);
        }
        for link in next.links().values() {
            assert_eq!(link.conductance_w_k, 999.0);
        }
        check_energy(&report, 0.0, q, 4.0 * heat_per_cell, 100.0, dt);
        near(next.mass_kg(), old.mass_kg());
    }
}

#[test]
fn latent_endpoint_mobility_changes_retained_momentum_without_changing_temperature() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [(0.0, 0.0), (100.0, 1.0)]
            .map(|(h, f)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: 500.0,
                liquid_mass_fraction: f,
                bulk_density_kg_m3: 1.0,
            })
            .to_vec(),
    )
    .unwrap();
    let old = thermal_state(&domain, &curve, &[50.0; 4], &[]);
    let mut speeds = Vec::new();
    for permeability in [0.5, 2.0] {
        let mobility = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, permeability, 0.2).unwrap();
        let flow = MaterialFlowFeedback::try_new(
            &old,
            &domain,
            FlowEnthalpyInput {
                cell_ids: &ids,
                material_id: MAT,
                caloric: &curve,
                inlet_specific_h_j_kg: [Some(50.0); 6],
                max_cells: 4,
                max_connections: 64,
            },
            &config,
            (),
            &gate,
            512,
            |_, _: &DepositedCell, t| {
                liquid.at_temperature(Temperature::new(t)).map_err(|e| {
                    fs_material::MaterialError::State {
                        what: e.to_string(),
                    }
                })
            },
        )
        .unwrap()
        .with_mushy_mobility(&mobility)
        .unwrap()
        .with_unsteady_body_force_heating(force(), controls())
        .unwrap();
        let dt = 0.1;
        let drag = |h: f64| {
            let f = h / 100.0;
            100.0 * (1.0 - f).powi(2) / (permeability * f.powi(3))
        };
        // Retained inertia adds 4/dt. The caloric endpoint supplies beta(h),
        // and h-50=dt*(600+beta(h))*q(h)^2 follows actual physical heating.
        let speed = |h| 400.0 / (4.0 / dt + 2400.0 + 4.0 * drag(h));
        let (mut lo, mut hi) = (50.0, 52.0);
        for _ in 0..80 {
            let h = 0.5 * (lo + hi);
            if h - 50.0 > dt * (600.0 + drag(h)) * speed(h).powi(2) {
                hi = h;
            } else {
                lo = h;
            }
        }
        let h = 0.5 * (lo + hi);
        let q = speed(h);
        let (next, accepted, report) = flow
            .step_unsteady_implicit(&old, dt, coupled_controls(false), &[1.0; 4], || Ok(()))
            .unwrap();
        near(accepted.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q);
        for sample in &accepted.applied().unwrap().samples {
            near(sample.specific_enthalpy_j_kg, h);
            near(sample.temperature_k, 500.0);
            near(next.cells()[&sample.cell].specific_h_j_kg, h);
        }
        near(report.mechanical.darcy_heat_j, dt * 4.0 * drag(h) * q * q);
        near(report.mechanical.molecular_heat_j, dt * 2400.0 * q * q);
        check_energy(&report, 0.0, q, 4.0 * (h - 50.0), 100.0, dt);
        speeds.push(q);
    }
    assert!(speeds[1] > speeds[0]);
}

#[test]
fn constant_material_preserves_native_startup_and_can_resume_the_original_observer() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = source(100.0, false);
    let liquid = liquid(&card);
    let caloric = liquid_caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let flow = input(&ids, &liquid, &caloric)
        .feedback(&old, &domain, &config, (), &gate, 512)
        .unwrap()
        .with_unsteady_body_force_heating(force(), controls())
        .unwrap();
    let dt = 0.001;
    let (explicit, original, _) = old
        .advance_observed(&flow, dt, thermal_step(), || Ok(()))
        .unwrap();
    let (implicit, accepted, report) = flow
        .step_unsteady_implicit(&old, dt, coupled_controls(false), &[0.0; 4], || Ok(()))
        .unwrap();
    assert_eq!(accepted.unsteady_flow(), original.unsteady_flow());
    assert_eq!(
        accepted.unsteady_mechanical_step(),
        original.unsteady_mechanical_step()
    );
    assert_eq!(Some(report.mechanical), original.unsteady_heating_audit());
    for (id, cell) in implicit.cells() {
        near(cell.specific_h_j_kg, explicit.cells()[id].specific_h_j_kg);
    }
    let (_, continued, _) = implicit
        .advance_observed(&accepted, 2.0 * dt, thermal_step(), || Ok(()))
        .unwrap();
    let q0 = 400.0 / (4.0 / dt + 2400.0);
    let q1 = (400.0 + 4.0 * q0 / dt) / (4.0 / dt + 2400.0);
    near(continued.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q1);
    near(
        continued
            .unsteady_heating_audit()
            .unwrap()
            .initial_kinetic_energy_j,
        2.0 * q0 * q0,
    );
}

#[test]
fn nonlinear_work_source_and_late_cancel_refusals_preserve_both_accepted_states() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let config = closed();
    let gate = CancelGate::new();
    let ids: Vec<_> = (0..4).map(|c| Some(CellId(c))).collect();
    let card = rheological_card(40.0, 100.0, 2.0);
    let liquid = law(&card);
    let caloric = caloric(&card);
    let old = initial(&domain, 1.0, &caloric);
    let source = LiquidFlowRheology::try_new(input(&ids, &liquid, &caloric), "rate").unwrap();
    let make = |budget| {
        source
            .feedback(&old, &domain, &config, (), &gate, budget)
            .unwrap()
            .with_unsteady_body_force_heating(force(), controls())
            .unwrap()
    };
    assert!(matches!(
        make(1).step_unsteady_implicit(&old, 0.5, coupled_controls(true), &[1.0; 4], || Ok(())),
        Err(ImplicitFlowError::Flow(MaterialFlowError::Thermal(
            DepositionError::Budget(_)
        )))
    ));
    let flow = make(8192);
    let before = flow.unsteady_flow().unwrap().clone();
    assert!(matches!(
        flow.step_unsteady_implicit(&old, 0.5, coupled_controls(false), &[1.0; 4], || Ok(())),
        Err(ImplicitFlowError::Invalid(_))
    ));
    // An equilibrium rate seed is not an admitted retained velocity history.
    let stationary_source = source
        .quasisteady_feedback(
            &old,
            &domain,
            &config,
            (),
            &gate,
            8192,
            vec![0.1; 4],
            rheology_controls(),
        )
        .unwrap()
        .with_unsteady_body_force_heating(force(), controls())
        .unwrap();
    assert!(matches!(
        stationary_source.step_unsteady_implicit(
            &old,
            0.5,
            coupled_controls(true),
            &[1.0; 4],
            || Ok(())
        ),
        Err(ImplicitFlowError::Invalid(_))
    ));
    let mut one = coupled_controls(true);
    one.rheology.as_mut().unwrap().max_iterations = 1;
    assert!(matches!(
        flow.step_unsteady_implicit(&old, 0.5, one, &[1.0; 4], || Ok(())),
        Err(ImplicitFlowError::Flow(MaterialFlowError::Thermal(
            DepositionError::Budget(_)
        )))
    ));
    let mut polls = 0;
    let (complete, first, report) = flow
        .step_unsteady_implicit(&old, 0.5, coupled_controls(true), &[1.0; 4], || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    let mut repeat = 0;
    assert!(matches!(
        flow.step_unsteady_implicit(&old, 0.5, coupled_controls(true), &[1.0; 4], || {
            repeat += 1;
            if repeat == polls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ImplicitFlowError::Flow(MaterialFlowError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    let (retry, again, report2) = flow
        .step_unsteady_implicit(&old, 0.5, coupled_controls(true), &[1.0; 4], || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), complete.cells());
    assert_eq!(again.unsteady_flow(), first.unsteady_flow());
    assert_eq!(report2.mechanical, report.mechanical);
    assert_eq!(report2.momentum_solves, report.momentum_solves);
    assert_eq!(flow.unsteady_flow(), Some(&before));
    assert_eq!(flow.unsteady_time_s(), Some(0.0));
    assert_eq!(
        flow.unsteady_heating_audit(),
        Some(UnsteadyHeatingAudit::default())
    );
    assert_eq!(flow.momentum_solves(), 0);
    assert!(flow.applied().is_none());
    assert_eq!(old.energy_j(), 0.0);

    let bounded_card = rheological_card(40.0, 100.0, 0.01);
    let bounded_liquid = law(&bounded_card);
    let bounded_caloric = super::super::caloric(&bounded_card);
    let bounded_old = initial(&domain, 1.0, &bounded_caloric);
    let bounded_source =
        LiquidFlowRheology::try_new(input(&ids, &bounded_liquid, &bounded_caloric), "rate")
            .unwrap();
    let bounded = bounded_source
        .feedback(&bounded_old, &domain, &config, (), &gate, 8192)
        .unwrap()
        .with_unsteady_body_force_heating(force(), controls())
        .unwrap();
    assert!(matches!(
        bounded.step_unsteady_implicit(
            &bounded_old,
            0.5,
            coupled_controls(true),
            &[1.0; 4],
            || Ok(())
        ),
        Err(ImplicitFlowError::Flow(MaterialFlowError::Material { .. }))
    ));
    assert_eq!(bounded.unsteady_time_s(), Some(0.0));
    assert_eq!(bounded_old.energy_j(), 0.0);
}
