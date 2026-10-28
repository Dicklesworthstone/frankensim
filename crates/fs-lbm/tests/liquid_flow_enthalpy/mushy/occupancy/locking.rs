//! Independent four-face cavity: supported freezing conserves energy and torque.
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::locking::{
    StationaryPhaseSupport, SupportedFreezingControls,
};

fn capture() -> SupportedFreezingControls {
    SupportedFreezingControls {
        cell_energy_tolerance_j: 1e-10,
        energy_tolerance_j: 1e-9,
    }
}

fn forces(f: f64) -> Vec<[f64; 3]> {
    vec![[-f, f, 0.0], [-f, -f, 0.0], [f, f, 0.0], [f, -f, 0.0]]
}

fn initial<'a>(
    domain: &VoxelDomain,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    h: [f64; 4],
    power: [f64; 4],
    changes: &[(f64, [f64; 4])],
) -> DepositionSimulation<'a> {
    let boundary = |source_w| Boundary {
        ambient_k: 500.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    };
    let mut events = Vec::new();
    for c in 0..4 {
        let [x, y, z] = domain.coords(c);
        let mut start = domain.center(x, y, z);
        let mut end = start;
        start[0] -= 0.5;
        end[0] += 0.5;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(c as u64),
                segment_id: c as u64,
                body_id: 1,
                material_id: MAT,
                mass_kg: 1.0,
                specific_h_j_kg: h[c],
                geometry: CellGeometry {
                    start_m: start,
                    end_m: end,
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: boundary(power[c]),
            }),
        });
    }
    for &(time_s, powers) in changes {
        for (c, power) in powers.into_iter().enumerate() {
            events.push(TimedEvent {
                time_s,
                event: Event::Boundary(CellId(c as u64), boundary(power)),
            });
        }
    }
    let max_events = events.len();
    DepositionSimulation::new(
        events,
        &[(MAT, curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 16,
            max_events,
            max_substeps: 32,
        },
    )
    .unwrap()
    .advance_to(0.0, step(0.1), || Ok(()))
    .unwrap()
    .0
}

fn caloric(card: &MaterialCard, latent: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-100.0, 450.0, 0.0),
            (0.0, 500.0, 0.0),
            (latent, 500.0, 1.0),
            (latent + 200.0, 600.0, 1.0),
        ]
        .into_iter()
        .map(|(h, t, phi)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: phi,
            bulk_density_kg_m3: 1.0,
        })
        .collect(),
    )
    .unwrap()
}

#[test]
fn supported_cavity_freezes_with_native_kinetic_heat_and_frame_angular_impulse_then_remelts() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 4, _>(|c| Some(CellId(c as u64)));
    let supported = [CellId(0), CellId(1), CellId(2), CellId(3)];
    let gate = CancelGate::new();
    let config = closed();
    let mut energies = Vec::new();
    for mu in [100.0, 200.0] {
        let card = source(mu);
        let curve = chart(&card);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let original = initial(
            &domain,
            &curve,
            [150.0; 4],
            [-1000.0; 4],
            &[(0.2, [2000.0; 4])],
        );
        let support = StationaryPhaseSupport::try_new(&original, &supported).unwrap();
        let flow = MaterialFlowFeedback::try_new_phase(
            &original,
            &domain,
            input(&ids, &curve),
            &config,
            (),
            &gate,
            20,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
        let (moving, history, support, startup) = flow
            .advance_with_supported_freezing(
                &original,
                &support,
                0.01,
                step(0.01),
                capture(),
                || Ok(()),
            )
            .unwrap();
        // Four native unit face masses; symmetric stress has resistance 24*mu.
        // Each cell receives 6*mu*q², independently of the source work split.
        let q0 = 100.0 / (100.0 + 6.0 * mu);
        near(history.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q0);
        assert_eq!(startup.capture_heat_j, 0.0);
        let coast = history
            .with_unsteady_body_force_heating(forces(0.0), controls())
            .unwrap();
        assert!(
            moving
                .advance_observed(&coast, 0.2, step(0.2), || Ok(()))
                .is_err()
        );
        let (frozen, history, support, report) = coast
            .advance_with_supported_freezing(
                &moving,
                &support,
                0.2,
                step(0.2),
                capture(),
                || Ok(()),
            )
            .unwrap();
        let dt = 0.2 - 0.01;
        let q1 = q0 / (1.0 + 6.0 * mu * dt);
        let kinetic = 2.0 * q1 * q1;
        let physical = dt * 24.0 * mu * q1 * q1;
        let numerical = 2.0 * (q1 - q0).powi(2);
        near(report.capture_heat_j, kinetic);
        near(report.initial_capture_kinetic_energy_j, kinetic);
        assert_eq!(report.final_capture_kinetic_energy_j, 0.0);
        assert_eq!(report.captured_faces, 4);
        assert_eq!(history.unsteady_time_s(), Some(frozen.time_s()));
        assert_eq!(support.time_s(), frozen.time_s());
        near(support.capture_heat_j(), kinetic);
        for p in report.support_impulse_n_s {
            near(p, 0.0);
        }
        near(report.support_angular_impulse_n_m_s[2], -2.0 * q1);
        near(support.angular_impulse_n_m_s()[2], -2.0 * q1);
        for c in 0..4 {
            near(report.cell_capture_heat_j[&CellId(c)], 0.5 * q1 * q1);
            assert!(
                history.endpoint()[c as usize]
                    .transport
                    .mobile_transport()
                    .is_none()
            );
            assert_eq!(
                history.unsteady_flow().unwrap().velocity_m_s[c as usize],
                [0.0; 3]
            );
            near(
                frozen.cells()[&CellId(c)].specific_h_j_kg,
                moving.cells()[&CellId(c)].specific_h_j_kg - 1000.0 * dt
                    + physical / 4.0
                    + kinetic / 4.0,
            );
        }
        let mechanical = history.unsteady_heating_audit().unwrap();
        near(mechanical.molecular_heat_j, physical);
        near(mechanical.support_capture_heat_j, kinetic);
        near(mechanical.generated_heat_j, physical + kinetic);
        near(mechanical.backward_euler_loss_j, numerical);
        assert_eq!(mechanical.final_kinetic_energy_j, 0.0);
        near(
            report.thermal.external_heat_j,
            -4000.0 * dt + physical + kinetic,
        );
        near(
            frozen.energy_j() - moving.energy_j() - 2.0 * q0 * q0,
            -4000.0 * dt - numerical,
        );
        near(frozen.mass_kg(), original.mass_kg());
        let (remelted, history, support2, report) = history
            .advance_with_supported_freezing(&frozen, &support, 0.3, step(0.11), capture(), || {
                Ok(())
            })
            .unwrap();
        assert_eq!(report.capture_heat_j, 0.0);
        assert_eq!(support2.impulse_n_s(), support.impulse_n_s());
        assert_eq!(
            support2.angular_impulse_n_m_s(),
            support.angular_impulse_n_m_s()
        );
        assert!(
            history
                .endpoint()
                .iter()
                .all(|s| s.transport.mobile_transport().is_some())
        );
        let driven = history
            .with_unsteady_body_force_heating(forces(100.0), controls())
            .unwrap();
        let (_, again, _, _) = driven
            .advance_with_supported_freezing(
                &remelted,
                &support2,
                0.31,
                step(0.02),
                capture(),
                || Ok(()),
            )
            .unwrap();
        near(again.unsteady_flow().unwrap().field.flux_x(1, 1, 0), q0);
        energies.push(kinetic);
    }
    // Changing only the source viscosity changes the retained motion and heat.
    assert!(energies[1] < energies[0]);
}

#[test]
fn changing_only_latent_storage_changes_whether_the_same_cooled_material_locks() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 4, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let card = source(100.0);
    let mut captures = Vec::new();
    for latent in [100.0, 200.0] {
        let curve = caloric(&card, latent);
        let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let original = initial(&domain, &curve, [latent + 50.0; 4], [-1000.0; 4], &[]);
        near(original.temperature_k(CellId(0)).unwrap(), 525.0);
        let support = StationaryPhaseSupport::try_new(
            &original,
            &[CellId(0), CellId(1), CellId(2), CellId(3)],
        )
        .unwrap();
        let flow = MaterialFlowFeedback::try_new_phase(
            &original,
            &domain,
            input(&ids, &curve),
            &config,
            (),
            &gate,
            20,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
        let (moving, history, support, _) = flow
            .advance_with_supported_freezing(
                &original,
                &support,
                0.01,
                step(0.01),
                capture(),
                || Ok(()),
            )
            .unwrap();
        let coast = history
            .with_unsteady_body_force_heating(forces(0.0), controls())
            .unwrap();
        let (_, history, _, report) = coast
            .advance_with_supported_freezing(
                &moving,
                &support,
                0.2,
                step(0.2),
                capture(),
                || Ok(()),
            )
            .unwrap();
        captures.push(report.capture_heat_j);
        assert_eq!(
            history.endpoint()[0].transport.mobile_transport().is_none(),
            latent == 100.0
        );
    }
    assert!(captures[0] > 0.0);
    assert_eq!(captures[1], 0.0);
}

#[test]
fn capture_heating_that_remelts_the_proposed_solid_refuses_without_clipping() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 4, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let q0 = 100.0 / 700.0;
    let q1 = q0 / 7.0;
    let warm_h = 150.0 + 0.01 * 600.0 * q0 * q0;
    // The native thermal step ends at -q1²/4, then capture would add q1²/2.
    let cooling = -(warm_h + 0.01 * 600.0 * q1 * q1 + 0.25 * q1 * q1) / 0.01;
    let original = initial(
        &domain,
        &curve,
        [150.0; 4],
        [0.0; 4],
        &[(0.01, [cooling; 4])],
    );
    let support =
        StationaryPhaseSupport::try_new(&original, &[CellId(0), CellId(1), CellId(2), CellId(3)])
            .unwrap();
    let flow = MaterialFlowFeedback::try_new_phase(
        &original,
        &domain,
        input(&ids, &curve),
        &config,
        (),
        &gate,
        20,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_body_force_heating(forces(100.0), controls())
    .unwrap();
    let (moving, history, support, _) = flow
        .advance_with_supported_freezing(
            &original,
            &support,
            0.01,
            step(0.01),
            capture(),
            || Ok(()),
        )
        .unwrap();
    let coast = history
        .with_unsteady_body_force_heating(forces(0.0), controls())
        .unwrap();
    let old_flow = coast.unsteady_flow().unwrap().clone();
    let error = coast
        .advance_with_supported_freezing(&moving, &support, 0.02, step(0.01), capture(), || Ok(()))
        .err()
        .unwrap()
        .to_string();
    assert!(
        error.contains("capture heat changes the final phase topology"),
        "{error}"
    );
    assert_eq!(coast.unsteady_flow(), Some(&old_flow));
    near(moving.cells()[&CellId(0)].specific_h_j_kg, warm_h);
    assert_eq!(support.capture_heat_j(), 0.0);
    assert_eq!(support.time_s(), 0.01);
}

#[test]
fn unsupported_interface_early_lock_and_late_cancellation_preserve_every_owner() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 4, _>(|c| Some(CellId(c as u64)));
    let gate = CancelGate::new();
    let config = closed();
    let card = source(100.0);
    let curve = chart(&card);
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    for partial in [false, true] {
        let power = if partial {
            [-1000.0, 0.0, 0.0, 0.0]
        } else {
            [-1000.0; 4]
        };
        let original = initial(&domain, &curve, [150.0; 4], power, &[]);
        let support = StationaryPhaseSupport::try_new(
            &original,
            &[CellId(0), CellId(1), CellId(2), CellId(3)],
        )
        .unwrap();
        let flow = MaterialFlowFeedback::try_new_phase(
            &original,
            &domain,
            input(&ids, &curve),
            &config,
            (),
            &gate,
            20,
            &law,
            |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
        )
        .unwrap()
        .with_unsteady_body_force_heating(forces(100.0), controls())
        .unwrap();
        let (moving, history, support, _) = flow
            .advance_with_supported_freezing(
                &original,
                &support,
                0.01,
                step(0.01),
                capture(),
                || Ok(()),
            )
            .unwrap();
        let coast = history
            .with_unsteady_body_force_heating(forces(0.0), controls())
            .unwrap();
        let old_cells = moving.cells().clone();
        let old_flow = coast.unsteady_flow().unwrap().clone();
        let old_endpoint = coast.endpoint().to_vec();
        let old_heat = coast.unsteady_heating_audit();
        if partial {
            let error = coast
                .advance_with_supported_freezing(
                    &moving,
                    &support,
                    0.2,
                    step(0.2),
                    capture(),
                    || Ok(()),
                )
                .err()
                .unwrap()
                .to_string();
            assert!(error.contains("pressure-impulse redistribution"), "{error}");
        } else {
            let missing = StationaryPhaseSupport::try_new(&moving, &[CellId(0)]).unwrap();
            let error = coast
                .advance_with_supported_freezing(
                    &moving,
                    &missing,
                    0.2,
                    step(0.2),
                    capture(),
                    || Ok(()),
                )
                .err()
                .unwrap()
                .to_string();
            assert!(error.contains("explicitly stationary support"), "{error}");
            let early = coast
                .advance_with_supported_freezing(
                    &moving,
                    &support,
                    0.2,
                    step(0.05),
                    capture(),
                    || Ok(()),
                )
                .err()
                .unwrap()
                .to_string();
            assert!(early.contains("handoff"), "{early}");
            let mut calls = 0usize;
            let accepted = coast
                .advance_with_supported_freezing(
                    &moving,
                    &support,
                    0.2,
                    step(0.2),
                    capture(),
                    || {
                        calls += 1;
                        Ok(())
                    },
                )
                .unwrap();
            let mut checked = 0usize;
            assert!(
                coast
                    .advance_with_supported_freezing(
                        &moving,
                        &support,
                        0.2,
                        step(0.2),
                        capture(),
                        || {
                            checked += 1;
                            if checked == calls - 1 {
                                Err(DepositionError::Cancelled)
                            } else {
                                Ok(())
                            }
                        }
                    )
                    .is_err()
            );
            let retry = coast
                .advance_with_supported_freezing(
                    &moving,
                    &support,
                    0.2,
                    step(0.2),
                    capture(),
                    || Ok(()),
                )
                .unwrap();
            assert_eq!(retry.0.cells(), accepted.0.cells());
            assert_eq!(retry.1.unsteady_flow(), accepted.1.unsteady_flow());
            assert_eq!(retry.2, accepted.2);
            assert_eq!(retry.3, accepted.3);
        }
        assert_eq!(moving.cells(), &old_cells);
        assert_eq!(coast.unsteady_flow(), Some(&old_flow));
        assert_eq!(coast.endpoint(), old_endpoint);
        assert_eq!(coast.unsteady_heating_audit(), old_heat);
        assert_eq!(support.capture_heat_j(), 0.0);
    }
}
