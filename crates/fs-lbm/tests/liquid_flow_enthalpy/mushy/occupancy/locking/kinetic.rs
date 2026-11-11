//! A small freezing circulation retains its energy beside a large live cavity.
use super::*;
use fs_lbm::conjugate::Voxel;

#[test]
fn accepted_freezing_retains_small_kinetic_heat_beside_large_surviving_circulation() {
    // Two separate 2x2 cavities have the same velocity stencil and viscous
    // power. Scaling the second square increases only its face masses by L².
    // Its unchanged kinetic storage is 2^59 J; the first cavity loses 49/32 J.
    let scale = 2.0_f64.powi(30);
    let mu = 2.0_f64.powi(36);
    let dt = 1.0 / mu;
    let mut domain = VoxelDomain::graded([
        vec![1.0, 1.0, 1.0, scale, scale],
        vec![1.0, 1.0, 1.0, scale, scale],
        vec![1.0],
    ])
    .unwrap();
    for y in 0..5 {
        for x in 0..5 {
            if !((x < 2 && y < 2) || (x >= 3 && y >= 3)) {
                domain.set(x, y, 0, Voxel::Solid(0));
            }
        }
    }
    let card = source(mu);
    // A zero liquid enthalpy datum lets the original caloric inventory retain
    // small physical heating even in the much more massive surviving cavity.
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (-200.0, 450.0, 0.0),
            (-100.0, 500.0, 0.0),
            (0.0, 500.0, 1.0),
            (200.0, 600.0, 1.0),
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
    .unwrap();
    let law = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let boundary = |source_w| Boundary {
        ambient_k: 500.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    };
    let mut events = Vec::new();
    // `DepositionSimulation::new` refuses events that are not ordered in absolute seconds (38f3c07f01): every cell's Birth is at
    // t = 0 and the small cells' Boundary at t = dt, so the Boundary events are collected here and appended after all Births.
    let mut boundary_events = Vec::new();
    let mut ids = vec![None; domain.cell_count()];
    let mut supported = Vec::new();
    let mut drive = vec![[0.0; 3]; domain.cell_count()];
    let mut hold = drive.clone();
    for c in 0..domain.cell_count() {
        if !domain.is_fluid(c) {
            continue;
        }
        let [x, y, z] = domain.coords(c);
        let small = x < 2;
        let id = CellId(c as u64);
        ids[c] = Some(id);
        let mut start = domain.center(x, y, z);
        let mut end = start;
        let width = if small { 1.0 } else { scale };
        start[0] -= 0.5 * width;
        end[0] += 0.5 * width;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: c as u64,
                body_id: 1,
                material_id: MAT,
                mass_kg: width * width,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: start,
                    end_m: end,
                    volume_m3: width * width,
                    surface_area_m2: 2.0 * width * width + 4.0 * width,
                },
                boundary: boundary(0.0),
            }),
        });
        let local_x = if small { x } else { x - 3 };
        let local_y = if small { y } else { y - 3 };
        let circulation = |force: f64| {
            [
                if local_y == 0 { -force } else { force },
                if local_x == 0 { force } else { -force },
                0.0,
            ]
        };
        drive[c] = circulation(if small { (49.0 / 8.0) * mu } else { 0.5 * mu });
        if small {
            supported.push(id);
            boundary_events.push(TimedEvent {
                time_s: dt,
                event: Event::Boundary(id, boundary(-128.0 * mu)),
            });
        } else {
            hold[c] = circulation(3.0 * mu / (scale * scale));
        }
    }
    events.extend(boundary_events);
    let original = DepositionSimulation::new(
        events,
        &[(MAT, &curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 8,
            max_links: 32,
            max_events: 12,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(0.0, step(dt), || Ok(()))
    .unwrap()
    .0;
    let gate = CancelGate::new_clock_free();
    let mut config = closed();
    config.velocity_relaxation = 0.5;
    let input = FlowEnthalpyInput {
        cell_ids: &ids,
        material_id: MAT,
        caloric: &curve,
        inlet_specific_h_j_kg: [None; 6],
        max_cells: 25,
        max_connections: 64,
    };
    let flow = MaterialFlowFeedback::try_new_phase(
        &original,
        &domain,
        input,
        &config,
        (),
        &gate,
        4,
        &law,
        |_, cell: &DepositedCell, _| resolve(&card, &law, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_body_force_heating(
        drive,
        ViscousHeatControls {
            maximum_reynolds: 0.1,
            // Priming the massive cavity uses a scale-aware roundoff allowance.
            // The measured decay/freezing interval below uses an absolute 1e-8 J.
            power_tolerance_w: 16.0 * f64::EPSILON * scale * scale / dt,
        },
    )
    .unwrap();
    let support = StationaryPhaseSupport::try_new(&original, &supported).unwrap();
    let (moving, history, support, _) = flow
        .advance_with_supported_freezing(&original, &support, dt, step(dt), capture(), || Ok(()))
        .unwrap();
    near(
        history.unsteady_flow().unwrap().field.flux_x(1, 1, 0),
        7.0 / 8.0,
    );
    near(
        history.unsteady_flow().unwrap().field.flux_x(4, 4, 0) / scale,
        0.5,
    );
    let coast = history
        .with_unsteady_body_force_heating(
            hold,
            ViscousHeatControls {
                maximum_reynolds: 0.1,
                power_tolerance_w: 1e-8 / dt,
            },
        )
        .unwrap();
    // The interval is dt (2 dt less the dt already advanced) and a circulation is present. Measured on ts2 (2026-11-11): with a
    // step limit of dt this call is refused ("freezing a moving native face requires an explicit momentum/kinetic-energy
    // handoff", the harmonic advection CFL of 52457eb93a cuts a step limit equal to the interval); with 10 dt the audit
    // reports one substep of duration dt and one momentum solve. A smaller multiple above 1 was not tried.
    let (frozen, accepted, _, report) = coast
        .advance_with_supported_freezing(
            &moving,
            &support,
            2.0 * dt,
            step(10.0 * dt),
            capture(),
            || Ok(()),
        )
        .unwrap();
    let native = accepted.unsteady_mechanical_step().unwrap();
    near(native.kinetic_energy_change_j, -1.5);
    near(native.backward_euler_loss_j, 9.0 / 8.0);
    assert_eq!(
        native.initial_kinetic_energy_j,
        native.final_kinetic_energy_j
    );
    assert_eq!(
        report.initial_capture_kinetic_energy_j,
        report.final_capture_kinetic_energy_j
    );
    near(report.capture_kinetic_energy_change_j, -1.0 / 32.0);
    near(report.capture_heat_j, 1.0 / 32.0);
    assert_eq!(report.captured_faces, 4);
    let work = accepted.unsteady_heating_audit().unwrap();
    // The hold is a force density [N/m3] (`with_unsteady_body_force_heating`), 3 mu / scale^2 on each of the FOUR large
    // cells (a 2 x 2 cavity of width scale and thickness 1, kinetic energy 2^59 J at u = 0.5). Per large cell its work is
    // (3 mu / scale^2) x 0.5 x scale^2 x dt = 1.5 J, which is also the cell's viscous heat (the hold balances it), so the
    // large cavity supplies 4 x 1.5 = 6 J and dissipates 6 J, its kinetic energy unchanged. The small cavity's four cells
    // dissipate 4 x 3/32 = 3/8 J. Molecular heat = 3/8 + 6 = 51/8 J.
    near(work.kinetic_energy_change_j, -49.0 / 32.0);
    near(work.molecular_heat_j, 51.0 / 8.0);
    near(work.support_capture_heat_j, 1.0 / 32.0);
    near(work.supplied_body_work_j, 6.0);
    near(work.residual_j, 0.0);
    for id in supported {
        near(report.cell_capture_heat_j[&id], 1.0 / 128.0);
        near(frozen.cells()[&id].specific_h_j_kg, -15783.0 / 128.0);
    }
    near(
        accepted.unsteady_flow().unwrap().field.flux_x(4, 4, 0) / scale,
        0.5,
    );
    // Stored energy change = -512 (unchanged) + the heat deposited, 51/8 + 1/32 = 205/32.
    near(report.thermal.stored_energy_change_j, -512.0 + 205.0 / 32.0);
    near(
        report.thermal.stored_energy_change_j
            + work.kinetic_energy_change_j
            + work.backward_euler_loss_j
            - work.supplied_body_work_j,
        -512.0,
    );
}
