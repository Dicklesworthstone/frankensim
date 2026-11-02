//! G1/G3/G4: original finite-body adhesive trajectories and finite heat storage.
use super::*;
use fs_couple::friction_heat::adhesive::transient::*;
use fs_mbd::tangential_friction::adhesive::transient::{
    AdhesiveTransientControls, PreparedRigidAdhesiveTransient,
};

fn phase_controls() -> RigidAdhesiveTransientHeatControls {
    let old = controls();
    let mut mechanical = old.mechanical.mechanical;
    mechanical.energy_tolerance_j = 1e-8;
    mechanical.contact.energy_tolerance_j = 2e-9;
    RigidAdhesiveTransientHeatControls {
        mechanical,
        transient: AdhesiveTransientControls {
            maximum_step_s: 0.05,
            minimum_step_s: 1e-16,
            maximum_attempts: 100_000,
            maximum_evaluations: 400_000,
            scaled_rate_tolerance_m_s: 2e-8,
            scaled_impulse_tolerance_ns: 2e-8,
            scaled_travel_tolerance_m: 2e-8,
            heat_tolerance_j: 1e-9,
            stopping_time_tolerance_s: 2e-7,
        },
        maximum_cells: 2,
        maximum_unknowns: 3,
        maximum_trials: old.maximum_trials,
        maximum_newton_iterations: old.maximum_newton_iterations,
        difference_step_j: old.difference_step_j,
        temperature_scale_j_per_k: old.temperature_scale_j_per_k,
        heat: old.heat,
        combined_energy_tolerance_j: 2e-8,
        newton: old.newton,
        audit_entropy: true,
    }
}

fn pair_thermal(initial: [f64; 2], area: f64) -> DepositionSimulation<'static> {
    let mut events = Vec::new();
    for i in 0..2 {
        let centre = [0.0, -2.0 * HALF_LENGTH][i];
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: CAPACITY[i],
                specific_h_j_kg: initial[i] - 300.0,
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, centre - HALF_LENGTH],
                    end_m: [0.0, 0.0, centre + HALF_LENGTH],
                    volume_m3: area * 2.0 * HALF_LENGTH,
                    surface_area_m2: 2.0 * area + 8.0 * area.sqrt() * HALF_LENGTH,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK[0],
            a: CellId(0),
            b: CellId(1),
            area_m2: area,
            conductance_w_k: series_conductance_w_k(area, [HALF_LENGTH; 2], [K[0], K[1]]).unwrap(),
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &CALORIC as &dyn ThermalLaw),
            (MaterialId(1), &CALORIC as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn contact<'a>(
    sources: &'a Sources,
    route: &'a FrictionHeatRoute,
) -> RigidAdhesiveTransientHeatContact<'a> {
    RigidAdhesiveTransientHeatContact {
        contact: RigidAdhesiveContact {
            binding: RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: Vec3::new(0.0, 0.0, -HALF_LENGTH),
                arm_b_body: Vec3::new(0.0, 0.0, HALF_LENGTH),
            },
            frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
            source: sources.patches()[0],
        },
        route,
        temperature_axis: "T",
        temperature_guess_k: 300.0,
    }
}

fn pair_transport(
    state: &DepositionSimulation<'_>,
    interface: f64,
) -> Result<RigidAdhesiveTransientHeatSource<[f64; 3]>, &'static str> {
    Ok(RigidAdhesiveTransientHeatSource {
        conductivity_w_m_k: [K[0], K[1]],
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: [
                state
                    .temperature_k(CellId(0))
                    .map_err(|_| "temperature A")?,
                state
                    .temperature_k(CellId(1))
                    .map_err(|_| "temperature B")?,
                interface,
            ],
        },
    })
}

// Independent elimination: two linear finite-cell balances leave one scalar
// interface balance. Analytic disk torque supplies a bounded impulse through
// exact spin arrest; no production mechanics/Newton/thermal solver is used.
fn spin_reference(
    initial: [f64; 2],
    area: f64,
    radius: f64,
    dt: f64,
    spin: f64,
) -> (f64, [f64; 2], f64, f64, f64) {
    let reduced = 1.0 / (1.0 / INERTIA[0] + 1.0 / INERTIA[1]);
    let d = [
        dt * K[0] * area / HALF_LENGTH,
        dt * K[1] * area / HALF_LENGTH,
    ];
    let evaluate = |face: f64| {
        let torque = 2.0 * (STRENGTH[0] - SLOPE[0] * (face - 300.0)) * area * radius / 3.0;
        let impulse = (dt * torque).min(reduced * spin);
        let heat = impulse * spin - 0.5 * impulse * impulse / reduced;
        let cells = std::array::from_fn::<_, 2, _>(|i| {
            (CAPACITY[i] * initial[i] + d[i] * face) / (CAPACITY[i] + d[i])
        });
        let residual = d[0] * (face - cells[0]) + d[1] * (face - cells[1]) - heat;
        (residual, cells, heat, impulse, reduced * spin / torque)
    };
    let mut lo = 280.0;
    let mut hi = 320.0;
    assert!(evaluate(lo).0 < 0.0 && evaluate(hi).0 > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if evaluate(mid).0 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    let face = 0.5 * (lo + hi);
    let (_, cells, heat, impulse, stop) = evaluate(face);
    (face, cells, heat, impulse, stop)
}

#[test]
fn pure_spin_heats_finite_cells_through_exact_arrest_and_signed_cooling() {
    let sources = Sources::new();
    let area = sources.areas()[0];
    let radius = sources.radii()[0];
    let initial = [318.0, 292.0];
    let mut heat = pair_thermal(initial, area);
    let route = FrictionHeatRoute::from_link(&heat, LINK[0], [K[0], K[1]]).unwrap();
    let mut island = Island::new(true);
    let original = island.bodies;
    let dt = 4.0;
    let expected = spin_reference(initial, area, radius, dt, SPIN[0] - SPIN[1]);
    assert!(expected.4 < dt);
    let result = apply_rigid_adhesive_transient_heat(
        &mut island.bodies,
        &island.properties,
        contact(&sources, &route),
        &mut heat,
        0.0,
        dt,
        phase_controls(),
        &[0.0; 2],
        pair_transport,
        &CancelGate::new(),
    )
    .unwrap();
    let path = result.mechanical.transient.as_ref().unwrap();
    assert!(path.exact_ray);
    let stopped = path.stopping_time_bracket_s.unwrap();
    close(stopped[0], expected.4, 2e-8);
    close(stopped[1], expected.4, 2e-8);
    close(result.heat.interface_temperature_k, expected.0, 2e-7);
    close(result.heat.dissipated_j, expected.2, 1e-9);
    close(
        result.mechanical.spin_angular_impulse_on_a_world_nms.z,
        -expected.3,
        1e-10,
    );
    for i in 0..2 {
        close(
            heat.temperature_k(CellId(i as u64)).unwrap(),
            expected.1[i],
            2e-7,
        );
        close(result.source[i], expected.1[i], 2e-7);
        close(
            result.heat.body_heat_j[i],
            CAPACITY[i] * (expected.1[i] - initial[i]),
            1e-9,
        );
        close(
            island.bodies[i].angular_momentum_body().z / INERTIA[i],
            1.6,
            1e-9,
        );
        assert_eq!(island.bodies[i].pose(), original[i].pose());
    }
    assert_eq!(island.bodies[2], original[2]);
    assert!(result.heat.body_heat_j[0] < 0.0 && result.heat.body_heat_j[1] > 0.0);
    assert!(result.heat.entropy.is_some());
    close(result.combined_energy_residual_j, 0.0, 2e-9);
    close(result.summed_energy_residual_j, 0.0, 2e-9);
    close(heat.time_s(), dt, 0.0);
}

#[test]
fn thermal_strength_feedback_matches_the_independent_partial_spin_root() {
    let sources = Sources::new();
    let initial = [318.0, 292.0];
    let heat = pair_thermal(initial, sources.areas()[0]);
    let route = FrictionHeatRoute::from_link(&heat, LINK[0], [K[0], K[1]]).unwrap();
    let island = Island::new(true);
    let dt = 0.25;
    let expected = spin_reference(initial, sources.areas()[0], sources.radii()[0], dt, 4.0);
    assert!(expected.4 > dt);
    let candidate = propose_rigid_adhesive_transient_heat(
        &island.bodies,
        &island.properties,
        contact(&sources, &route),
        &heat,
        0.0,
        dt,
        phase_controls(),
        &[0.0; 2],
        pair_transport,
        &CancelGate::new(),
    )
    .unwrap();
    let report = &candidate.report;
    assert!(
        report
            .mechanical
            .transient
            .as_ref()
            .unwrap()
            .stopping_time_bracket_s
            .is_none()
    );
    close(report.heat.interface_temperature_k, expected.0, 2e-7);
    close(report.heat.dissipated_j, expected.2, 2e-9);
    close(
        report.mechanical.spin_angular_impulse_on_a_world_nms.z,
        -expected.3,
        2e-10,
    );
    close(
        report
            .friction_source
            .law()
            .kinetic_patch()
            .parameters()
            .shear_strength_pa,
        STRENGTH[0] - SLOPE[0] * (expected.0 - 300.0),
        5e-6,
    );
    assert!(report.heat.interface_temperature_k > 305.0);
    for i in 0..2 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            expected.1[i],
            2e-7,
        );
    }
    // The immutable proposal did not advance either accepted owner.
    close(heat.time_s(), 0.0, 0.0);
    close(
        island.bodies[0].angular_momentum_body().z,
        INERTIA[0] * SPIN[0],
        0.0,
    );
}

#[test]
fn curved_slip_spin_uses_the_same_native_trajectory_and_one_heat_owner() {
    let sources = Sources::new();
    let heat = pair_thermal([302.0, 292.0], sources.areas()[0]);
    let route = FrictionHeatRoute::from_link(&heat, LINK[0], [K[0], K[1]]).unwrap();
    let mut island = Island::new(true);
    island.bodies[0] = RigidBodyState::new(
        island.bodies[0].pose(),
        Vec3::new(0.4 * MASS[0], -0.2 * MASS[0], 0.0),
        island.bodies[0].angular_momentum_body(),
    )
    .unwrap();
    let contact = contact(&sources, &route);
    let config = phase_controls();
    let dt = 0.04;
    let candidate = propose_rigid_adhesive_transient_heat(
        &island.bodies,
        &island.properties,
        contact,
        &heat,
        0.0,
        dt,
        config,
        &[0.0; 2],
        pair_transport,
        &CancelGate::new(),
    )
    .unwrap();
    let native = PreparedRigidAdhesiveTransient::new(
        &island.bodies,
        &island.properties,
        contact.contact,
        dt,
        config.mechanical,
        config.transient,
        &CancelGate::new(),
    )
    .unwrap()
    .advance_thermal(&candidate.report.friction_source, &CancelGate::new())
    .unwrap();
    assert_eq!(candidate.bodies_after, native.bodies_after);
    assert_eq!(candidate.report.mechanical, native.report);
    let path = candidate.report.mechanical.transient.as_ref().unwrap();
    assert!(!path.exact_ray && path.accepted_steps > 1);
    close(
        candidate.report.heat.dissipated_j,
        path.dissipated_energy_j,
        0.0,
    );
    close(
        energy(&candidate.thermal) - energy(&heat),
        path.dissipated_energy_j,
        2e-9,
    );
    close(candidate.report.combined_energy_residual_j, 0.0, 2e-8);
}

#[test]
fn arrest_continues_as_conduction_and_late_refusals_publish_nothing() {
    let sources = Sources::new();
    let mut heat = pair_thermal([318.0, 292.0], sources.areas()[0]);
    let route = FrictionHeatRoute::from_link(&heat, LINK[0], [K[0], K[1]]).unwrap();
    let mut island = Island::new(false);
    let contact = contact(&sources, &route);
    let saved = island.bodies;
    let before = energy(&heat);
    let baseline = propose_rigid_adhesive_transient_heat(
        &island.bodies,
        &island.properties,
        contact,
        &heat,
        0.0,
        0.2,
        phase_controls(),
        &[0.0; 2],
        pair_transport,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(baseline.bodies_after, [saved[0], saved[1]]);
    close(baseline.report.heat.dissipated_j, 0.0, 0.0);
    assert!(baseline.report.heat.body_heat_j[0] < 0.0 && baseline.report.heat.body_heat_j[1] > 0.0);
    close(energy(&baseline.thermal), before, 2e-9);
    let gate = CancelGate::new();
    let calls = Cell::new(0usize);
    let refused = apply_rigid_adhesive_transient_heat(
        &mut island.bodies,
        &island.properties,
        contact,
        &mut heat,
        0.0,
        0.2,
        phase_controls(),
        &[0.0; 2],
        |state, temperature| {
            calls.set(calls.get() + 1);
            if calls.get() == baseline.report.trials {
                gate.request();
            }
            pair_transport(state, temperature)
        },
        &gate,
    );
    assert!(matches!(
        refused,
        Err(RigidAdhesiveTransientHeatError::Thermal(
            DepositionError::Cancelled
        ))
    ));
    assert_eq!(calls.get(), baseline.report.trials);
    assert_eq!(island.bodies, saved);
    close(heat.time_s(), 0.0, 0.0);
    close(energy(&heat), before, 0.0);
    let mut wrong = contact;
    wrong.temperature_guess_k = 330.0;
    assert!(matches!(
        apply_rigid_adhesive_transient_heat(
            &mut island.bodies,
            &island.properties,
            wrong,
            &mut heat,
            0.0,
            0.2,
            phase_controls(),
            &[0.0; 2],
            pair_transport,
            &CancelGate::new(),
        ),
        Err(RigidAdhesiveTransientHeatError::FrictionSource(_))
    ));
    let mut bounded = phase_controls();
    bounded.maximum_trials = 1;
    assert!(matches!(
        apply_rigid_adhesive_transient_heat(
            &mut island.bodies,
            &island.properties,
            contact,
            &mut heat,
            0.0,
            0.2,
            bounded,
            &[0.0; 2],
            pair_transport,
            &CancelGate::new(),
        ),
        Err(RigidAdhesiveTransientHeatError::TrialBudget)
    ));
    assert_eq!(island.bodies, saved);
    close(heat.time_s(), 0.0, 0.0);
    close(energy(&heat), before, 0.0);
    let accepted = apply_rigid_adhesive_transient_heat(
        &mut island.bodies,
        &island.properties,
        contact,
        &mut heat,
        0.0,
        0.2,
        phase_controls(),
        &[0.0; 2],
        pair_transport,
        &CancelGate::new(),
    )
    .unwrap();
    close(
        accepted.heat.interface_temperature_k,
        baseline.report.heat.interface_temperature_k,
        0.0,
    );
    assert_eq!(island.bodies, saved);
    close(heat.time_s(), 0.2, 0.0);
}
