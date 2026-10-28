//! Actual moving modal consumer; retain the stationary fixture regressions too.
include!("modal_friction_heat.rs");
use fs_couple::friction_heat::moving as swept;
use fs_couple::render::schedule::force::coupled::contact::friction::heat::moving::{
    MovingModalHeatSource, MovingModalPatchSource,
};

const MOVING_K: [f64; 2] = [100.0, 200.0];
fn moving_fixture() -> (
    DepositionSimulation<'static>,
    swept::MovingFrictionRoute,
    swept::MovingContactPath,
) {
    let mut events = Vec::new();
    let mut faces: [Vec<swept::StripSurfaceCell>; 2] = std::array::from_fn(|_| Vec::new());
    let depth = 0.0005;
    let width = 0.001;
    let dx = 0.001;
    for side in 0..2 {
        for i in 0..2 {
            let lo = (i as f64 - 1.0) * dx;
            let hi = i as f64 * dx;
            let xc = 0.5 * (lo + hi);
            let (a, b) = if side == 0 {
                ([-depth, xc, 0.0], [0.0, xc, 0.0])
            } else {
                ([0.0, xc, 0.0], [depth, xc, 0.0])
            };
            let id = CellId((2 * side + i) as u64);
            faces[side].push(swept::StripSurfaceCell {
                cell: id,
                span_m: [lo, hi],
            });
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: id.0,
                    body_id: side as u64,
                    material_id: MaterialId(0),
                    mass_kg: if i == 0 {
                        2e-6 * (side + 1) as f64
                    } else {
                        8e-6 * (side + 1) as f64
                    },
                    specific_h_j_kg: if i == 0 { 0.0 } else { 3.0 },
                    geometry: CellGeometry {
                        start_m: a,
                        end_m: b,
                        volume_m3: dx * depth * width,
                        surface_area_m2: 2.0 * (dx * depth + dx * width + depth * width),
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
    }
    for side in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: side as u64,
                a: CellId((2 * side) as u64),
                b: CellId((2 * side + 1) as u64),
                area_m2: depth * width,
                conductance_w_k: MOVING_K[side] * depth * width / dx,
            }),
        });
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 2,
            max_events: 6,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let route = swept::MovingFrictionRoute::new(
        &thermal,
        swept::StripFrame {
            origin_m: [0.0; 3],
            tangent: [0.0, 1.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            width_m: width,
        },
        faces,
        [-0.00075, 0.00075],
        swept::MovingFrictionBudget {
            maximum_surface_cells: 4,
            maximum_overlap_tests: 8,
            maximum_patches: 4,
        },
        || Ok(()),
    )
    .unwrap();
    // Prescribed geometry moves continuously across the central material-cell
    // boundary. Motion refinement remains separate from thermal Newton accuracy.
    let path = swept::MovingContactPath {
        start: swept::MovingContactPosition {
            body_offsets_m: [0.0; 2],
            footprint_offset_m: 0.0,
        },
        end: swept::MovingContactPosition {
            body_offsets_m: [DT, 0.0],
            footprint_offset_m: 0.5 * DT,
        },
    };
    (thermal, route, path)
}
fn moving_controls() -> ModalFrictionHeatSolve {
    let mut c = controls();
    c.maximum_cells = 4;
    c.maximum_trials = 512;
    c.maximum_newton_iterations = 24;
    c.difference_step_j = 1e-9;
    // A numerical Jacobian cannot satisfy analytic-Jacobian forcing below its
    // subtraction floor. Keep the strict outer/local energy gates unchanged.
    c.newton.forcing_minimum = 1e-6;
    c.newton.linear_restart = 8;
    c
}
fn moving_source(
    state: &DepositionSimulation<'_>,
    patches: &[swept::SweptContactPatch],
    t: &[f64],
) -> Result<MovingModalHeatSource<Vec<f64>>, DepositionError> {
    let mut k = BTreeMap::new();
    for p in patches {
        for id in p.cells {
            k.insert(id, MOVING_K[state.cells()[&id].body_id as usize]);
        }
    }
    Ok(MovingModalHeatSource {
        // Explicit uniform pressure is a property of this analytical fixture.
        patches: patches
            .iter()
            .zip(t)
            .map(|(p, &temperature)| MovingModalPatchSource {
                normal_force_n: p.mean_area_m2 / 1.5e-6,
                law: law(temperature),
            })
            .collect(),
        conductivity_w_m_k: k,
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: t.to_vec(),
        },
    })
}
#[test]
fn g1_moving_modal_heat_uses_each_local_temperature_and_exact_original_patch_work() {
    let mut mechanical = mechanics();
    let (mut thermal, route, path) = moving_fixture();
    let count = swept::MovingFrictionHeatStep::new(&thermal, DT, &route, path, || Ok(()))
        .unwrap()
        .patches()
        .len();
    assert!(
        count >= 3,
        "cell crossing must create the intermediate overlap"
    );
    let before = mechanical.total_energy_j().unwrap() + thermal.energy_j();
    let result = mechanical
        .step_with_moving_friction_heat(
            &mut thermal,
            &route,
            path,
            &[0.0; 2],
            drive(),
            moving_controls(),
            &[0.0; 4],
            &vec![300.0; count],
            moving_source,
            &CancelGate::new(),
        )
        .unwrap();
    let mut total = 0.0;
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    for (p, &temperature) in result.heat.patches.iter().zip(&result.source) {
        assert_eq!(p.interface_temperature_k, temperature);
        min = min.min(temperature);
        max = max.max(temperature);
        let mu = law(temperature)
            .kinetic_coefficient((result.mechanical.slip_distance_m / DT).abs())
            .unwrap();
        let force = mu * p.mean_area_m2 / 1.5e-6;
        total += force;
        close(
            p.dissipated_j,
            force * result.mechanical.slip_distance_m.abs(),
            1e-14,
        );
        close(p.body_heat_j.iter().sum(), p.dissipated_j, 1e-12);
    }
    assert!(
        max - min > 0.1,
        "distinct material temperatures must not be globally averaged"
    );
    close(result.mechanical.traction_n, -total, 1e-11);
    close(
        result.heat.dissipated_j,
        result.mechanical.friction_dissipation_j,
        1e-14,
    );
    close(
        mechanical.total_energy_j().unwrap() + thermal.energy_j(),
        before,
        1e-12,
    );
    assert_eq!(result.heat.position, path.end);
    assert_eq!(mechanical.samples_rendered(), 1);
    assert_eq!(thermal.time_s(), DT);
}

#[test]
fn g3_moving_modal_refusal_after_staging_preserves_original_history_and_exact_retry() {
    let mut mechanical = mechanics();
    let (mut thermal, route, path) = moving_fixture();
    let count = swept::MovingFrictionHeatStep::new(&thermal, DT, &route, path, || Ok(()))
        .unwrap()
        .patches()
        .len();
    let initial = states(&mechanical);
    let heat = thermal.cells().clone();
    let calls = Cell::new(0);
    let refuse = |state: &DepositionSimulation<'_>, p: &[swept::SweptContactPatch], t: &[f64]| {
        calls.set(calls.get() + 1);
        if calls.get() == 3 {
            Err(DepositionError::MaterialDomain)
        } else {
            moving_source(state, p, t)
        }
    };
    assert!(matches!(
        mechanical.step_with_moving_friction_heat(
            &mut thermal,
            &route,
            path,
            &[0.0; 2],
            drive(),
            moving_controls(),
            &[0.0; 4],
            &vec![300.0; count],
            refuse,
            &CancelGate::new()
        ),
        Err(ModalFrictionHeatError::Source(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(states(&mechanical), initial);
    assert_eq!(thermal.cells(), &heat);
    assert_eq!(mechanical.port().law, law(300.0));
    assert_eq!(mechanical.samples_rendered(), 0);
    assert_eq!(thermal.time_s(), 0.0);
    let result = mechanical
        .step_with_moving_friction_heat(
            &mut thermal,
            &route,
            path,
            &[0.0; 2],
            drive(),
            moving_controls(),
            &[0.0; 4],
            &vec![300.0; count],
            moving_source,
            &CancelGate::new(),
        )
        .unwrap();
    let mut direct = mechanics();
    let (mut direct_thermal, direct_route, direct_path) = moving_fixture();
    let expected = direct
        .step_with_moving_friction_heat(
            &mut direct_thermal,
            &direct_route,
            direct_path,
            &[0.0; 2],
            drive(),
            moving_controls(),
            &[0.0; 4],
            &vec![300.0; count],
            moving_source,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(result.mechanical, expected.mechanical);
    assert_eq!(states(&mechanical), states(&direct));
    assert_eq!(thermal.cells(), direct_thermal.cells());
}
