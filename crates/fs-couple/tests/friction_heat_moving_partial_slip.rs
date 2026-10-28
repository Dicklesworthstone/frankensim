//! G0/G1/G3/G4/G5: original material memory, brush transport and finite heat.
//! Coefficients are analytical fixtures, not measured material calibration.
use fs_couple::friction_heat::moving::partial_slip::*;
use fs_couple::friction_heat::{FrictionHeatTolerances, moving::*};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::partial_slip::*;
use std::cell::Cell;
use std::collections::BTreeMap;
use std::convert::Infallible;

struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() && h > -300.0 {
            Ok(300.0 + h)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static CALORIC: Caloric = Caloric;
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:.16e} != {b:.16e}, tol {t:.3e}");
}
fn position(a: f64, b: f64, f: f64) -> MovingContactPosition {
    MovingContactPosition {
        body_offsets_m: [a, b],
        footprint_offset_m: f,
    }
}
fn frame() -> TangentFrame {
    TangentFrame::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]).unwrap()
}
fn law(mu: f64) -> PartialSlipLaw {
    PartialSlipLaw::new(
        "uniform-brush/analytical",
        "synthetic/two-spring-density",
        PartialSlipParameters {
            static_mu: mu,
            kinetic_mu: mu,
            tangential_stiffness_n_per_m: 100.0,
            torsional_stiffness_nm_per_rad: 10.0,
            torsional_capacity_factor: 0.5,
            partial_slip_onset_fraction: 0.5,
            partial_slip_hardening_fraction: 0.4,
        },
    )
    .unwrap()
}
fn initial(strain: [f64; 2]) -> PartialSlipState {
    if strain == [0.0; 2] {
        PartialSlipState::zero()
    } else {
        PartialSlipState::from_region_checkpoint(
            [[strain[0], 0.0], [strain[1], 0.0]],
            [0.0; 2],
            [false; 2],
            Some(frame()),
            7,
        )
        .unwrap()
    }
}
fn reference(mu: f64, state: PartialSlipState) -> PartialSlipCheckpoint {
    PartialSlipCheckpoint::new(
        NormalPatchView::new(
            "brush",
            "uniform pressure",
            "synthetic/normal",
            NormalPatchAuthority::SyntheticFixture,
            1.0,
            1.0,
            1.0 / core::f64::consts::PI,
            0.1,
        )
        .unwrap(),
        PartialSlipInterface::new(
            "A->B",
            "material-memory",
            "synthetic/interface",
            NormalPatchAuthority::SyntheticFixture,
        )
        .unwrap(),
        law(mu),
        state,
    )
    .unwrap()
}
fn work() -> GeneralizedWorkOwnership {
    GeneralizedWorkOwnership::new("brush", "interval", "x", "y", "spin").unwrap()
}
fn controls() -> MovingPartialSlipControls {
    MovingPartialSlipControls {
        maximum_regions: 256,
        maximum_cells: 32,
        maximum_trials: 500,
        maximum_newton_iterations: 30,
        difference_step_j: 1e-4,
        temperature_scale_j_per_k: 1.0,
        heat: FrictionHeatTolerances {
            cell_j: 2e-10,
            thermal_j: 3e-10,
            interface_j: 2e-10,
            total_j: 3e-10,
        },
        energy_tolerance_j: 5e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 16,
            ..Default::default()
        },
    }
}
fn fixture(
    edges: &[f64],
    mu: f64,
    strain: [f64; 2],
) -> (DepositionSimulation<'static>, MovingPartialSlip) {
    fixture_state(edges, mu, initial(strain))
}
fn fixture_state(
    edges: &[f64],
    mu: f64,
    state: PartialSlipState,
) -> (DepositionSimulation<'static>, MovingPartialSlip) {
    let n = edges.len() - 1;
    let mut events = Vec::new();
    let mut surfaces: [Vec<StripSurfaceCell>; 2] = std::array::from_fn(|_| Vec::new());
    for side in 0..2 {
        for i in 0..n {
            let id = CellId((side * n + i) as u64);
            let span = [edges[i], edges[i + 1]];
            let length = span[1] - span[0];
            let centre = f64::midpoint(span[0], span[1]);
            let y = if side == 0 { [-1.0, 0.0] } else { [0.0, 1.0] };
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id,
                    segment_id: id.0,
                    body_id: side as u64,
                    material_id: MaterialId(0),
                    mass_kg: 100.0 * length,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [centre, y[0], 0.0],
                        end_m: [centre, y[1], 0.0],
                        volume_m3: length,
                        surface_area_m2: 4.0 * length + 2.0,
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
            surfaces[side].push(StripSurfaceCell {
                cell: id,
                span_m: span,
            });
        }
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2 * n,
            max_links: 1,
            max_events: 2 * n,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let route = MovingFrictionRoute::new(
        &thermal,
        StripFrame {
            origin_m: [0.0; 3],
            tangent: [1.0, 0.0, 0.0],
            normal: [0.0, 1.0, 0.0],
            width_m: 1.0,
        },
        surfaces,
        [0.0, 1.0],
        MovingFrictionBudget {
            maximum_surface_cells: 2 * n,
            maximum_overlap_tests: 4 * n * n,
            maximum_patches: n * n,
        },
        || Ok(()),
    )
    .unwrap();
    let contact = MovingPartialSlip::new(
        &thermal,
        route,
        reference(mu, state),
        position(0.0, 0.0, 0.0),
        || Ok(()),
    )
    .unwrap();
    (thermal, contact)
}
fn source(
    thermal: &DepositionSimulation<'_>,
    regions: &[MovingSlipRegion],
    mu: f64,
) -> MovingPartialSlipSource<()> {
    MovingPartialSlipSource {
        laws: regions.iter().map(|_| law(mu)).collect(),
        conductivity_w_m_k: thermal
            .cells()
            .iter()
            .map(|(&id, c)| (id, if c.body_id == 0 { 1.0 } else { 3.0 }))
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    }
}
fn step(
    contact: &MovingPartialSlip,
    thermal: &DepositionSimulation<'static>,
    target: f64,
    end: MovingContactPosition,
    mu: f64,
) -> MovingPartialSlipCandidate<'static, ()> {
    contact
        .step_candidate(
            thermal,
            target,
            end,
            &work(),
            controls(),
            &vec![0.0; thermal.cells().len()],
            300.0,
            |t, r, _, _| Ok::<_, Infallible>(source(t, r, mu)),
            || false,
        )
        .unwrap()
}

#[test]
fn g0_material_preload_crosses_cells_without_reset_or_fictitious_heat() {
    let (thermal, contact) = fixture(&[-3.0, 0.0, 1.0, 4.0], 1000.0, [0.01, 0.02]);
    let before = contact.parcels().to_vec();
    let candidate = step(&contact, &thermal, 0.1, position(0.0, 0.0, 0.25), 1000.0);
    let r = candidate.report();
    near(r.storage_before_j, 0.011, 1e-14);
    near(r.elastic_outflow_j, 0.00275, 1e-14);
    near(r.storage_after_j, 0.00825, 1e-14);
    near(r.heat.dissipated_j, 0.0, 1e-14);
    near(r.heat.thermal.stored_energy_change_j, 0.0, 1e-14);
    near(r.mechanical_work_j, 0.0, 1e-14);
    near(r.numerical_dissipation_j, 0.0, 1e-14);
    assert_eq!(contact.parcels(), before);
    let retained = &candidate.contact().parcels()[0];
    assert_eq!(retained.span_m, [0.25, 1.0]);
    assert_eq!(
        retained.state().region_elastic_displacement_m(),
        [[0.01, 0.0], [0.02, 0.0]]
    );
    assert_eq!(retained.state().accepted_steps(), 8);
    assert_eq!(
        candidate.contact().parcels()[1]
            .state()
            .region_elastic_displacement_m(),
        [[0.0; 2]; 2]
    );
    // Returning the footprint cannot resurrect a departed material checkpoint.
    // The open boundary explicitly supplied unloaded entering material.
    let (next, heat, _) = candidate.into_parts();
    let returned = step(&next, &heat, 0.2, position(0.0, 0.0, 0.0), 1000.0);
    near(returned.report().storage_after_j, 0.00825, 1e-14);
    assert_eq!(returned.contact().parcels()[0].span_m, [0.0, 0.25]);
    assert_eq!(
        returned.contact().parcels()[0]
            .state()
            .region_elastic_displacement_m(),
        [[0.0; 2]; 2]
    );
}

#[test]
fn g1_native_two_region_work_heats_the_complete_swept_original_network_once() {
    let (thermal, contact) = fixture(&[-3.0, 0.0, 1.0, 4.0], 1.0, [0.0; 2]);
    // A and the footprint move together: the same A material retains its law,
    // while that contact visits different original B cells along the path.
    let candidate = step(&contact, &thermal, 0.1, position(0.02, 0.0, 0.02), 1.0);
    let r = candidate.report();
    near(r.force_on_a_n[0], -1.0, 1e-12);
    near(r.mechanical_work_j, 0.02, 1e-13);
    near(r.storage_after_j, 0.006875, 1e-13);
    near(r.numerical_dissipation_j, 0.006875, 1e-13);
    near(r.heat.dissipated_j, 0.00625, 1e-13);
    near(r.heat.thermal.stored_energy_change_j, 0.00625, 3e-10);
    near(r.elastic_outflow_j, 0.0, 1e-14);
    assert!(r.heat.patches.len() > 1);
    assert!(r.heat.patches.iter().all(|p| p.dissipated_j > 0.0));
    near(
        r.heat.patches.iter().map(|p| p.dissipated_j).sum(),
        0.00625,
        1e-13,
    );
    near(r.total_energy_residual_j, 0.0, 5e-10);
    // Different mechanical histories survive a reversal without being replaced
    // by the aggregate displacement or by a freshly initialized brush.
    let (c, t, _) = candidate.into_parts();
    let reversed = step(&c, &t, 0.2, position(0.019, 0.0, 0.019), 1.0);
    assert!(reversed.report().storage_after_j < 0.006875);
    near(reversed.report().heat.dissipated_j, 0.0, 1e-14);
}

#[test]
fn g3_common_translation_preserves_material_memory_and_distributed_shear_moment() {
    let state = PartialSlipState::from_region_checkpoint(
        [[0.0, 0.01]; 2],
        [0.0; 2],
        [false; 2],
        Some(frame()),
        4,
    )
    .unwrap();
    let (thermal, contact) = fixture_state(&[-3.0, 5.0], 1000.0, state.clone());
    let stationary = step(&contact, &thermal, 0.1, position(0.0, 0.0, 0.0), 1000.0);
    let translated = step(&contact, &thermal, 0.1, position(2.0, 2.0, 2.0), 1000.0);
    assert_eq!(
        stationary.contact().parcels(),
        translated.contact().parcels()
    );
    near(translated.report().force_on_a_n[2], 1.0, 1e-13);
    near(
        stationary.report().moment_on_a_about_origin_nm[1],
        -0.5,
        1e-13,
    );
    near(
        translated.report().moment_on_a_about_origin_nm[1],
        -1.5,
        1e-13,
    );
    near(translated.report().heat.dissipated_j, 0.0, 1e-14);
    near(translated.report().elastic_outflow_j, 0.0, 1e-14);
    near(translated.report().storage_after_j, 0.005, 1e-14);
    near(translated.report().total_energy_residual_j, 0.0, 1e-14);
    assert_eq!(contact.parcels()[0].state(), &state);
}

#[test]
fn g1_local_material_temperature_feedback_matches_independent_scalar_reference() {
    let (thermal, contact) = fixture(&[-3.0, 5.0], 1.0, [0.0; 2]);
    let dt = 0.02;
    let capacity = 800.0;
    let beta = 10000.0;
    let qa = 2.0 / (1.0 + dt * 2.0 / capacity);
    let qb = 6.0 / (1.0 + dt * 6.0 / capacity);
    let share = qa / (qa + qb);
    let mut lo = 0.0;
    let mut hi = 1.0;
    for _ in 0..100 {
        let mu = f64::midpoint(lo, hi);
        let d = 0.1 * mu - 0.01375 * mu * mu;
        let residual = mu * (1.0 + beta * share * d / capacity) - 1.0;
        if residual > 0.0 {
            hi = mu;
        } else {
            lo = mu;
        }
    }
    let mu = f64::midpoint(lo, hi);
    let heat = 0.1 * mu - 0.01375 * mu * mu;
    let candidate = contact
        .step_candidate(
            &thermal,
            dt,
            position(0.1, 0.0, 0.0),
            &work(),
            controls(),
            &[0.0; 2],
            300.0,
            |t, regions, _, _| {
                let mut s = source(t, regions, 1.0);
                s.laws = regions
                    .iter()
                    .map(|r| {
                        let temperature = t.temperature_k(r.surface_a_cell).unwrap();
                        law(1.0 / (1.0 + beta * (temperature - 300.0)))
                    })
                    .collect();
                Ok::<_, Infallible>(s)
            },
            || false,
        )
        .unwrap();
    near(candidate.report().force_on_a_n[0], -mu, 2e-8);
    near(candidate.report().heat.dissipated_j, heat, 3e-9);
    near(
        candidate.thermal().temperature_k(CellId(0)).unwrap(),
        300.0 + share * heat / capacity,
        3e-11,
    );
    assert!(mu < 0.9, "feedback must materially change friction");
    assert!(
        candidate
            .report()
            .regions
            .iter()
            .all(|r| r.region.surface_a_cell == CellId(0))
    );
}

#[test]
fn g1_material_brush_transport_converges_to_independent_steady_and_transient_integrals() {
    let mut errors = Vec::new();
    // Uniform material transport U=1, actual slip v=.1, reference length=1.
    // Each branch has z(age)=min(v*age,C/k). Age is uniform on [0,1].
    let branches = [(0.3, 60.0), (0.7, 40.0)];
    let exact_force: f64 = branches
        .iter()
        .map(|&(c, k)| c * (1.0 - 0.5 * c / (k * 0.1)))
        .sum();
    let exact_storage: f64 = branches
        .iter()
        .map(|&(c, k)| 0.5 * c * c / k * (1.0 - 2.0 / 3.0 * c / (k * 0.1)))
        .sum();
    let exact_heat: f64 = branches
        .iter()
        .map(|&(c, k)| {
            let tau = c / (k * 0.1);
            0.1 * c * (1.0 - tau) * (1.0 - tau)
        })
        .sum();
    for n in [16, 32, 64] {
        let (mut thermal, mut contact) = fixture(&[-3.0, 5.0], 1.0, [0.0; 2]);
        let mut generated = 0.0;
        let mut numerical = 0.0;
        let mut force = 0.0;
        for i in 1..=n {
            let t = i as f64 / n as f64;
            let candidate = step(&contact, &thermal, t, position(0.0, -0.1 * t, t), 1.0);
            generated += candidate.report().heat.dissipated_j;
            numerical += candidate.report().numerical_dissipation_j;
            force = -candidate.report().force_on_a_n[0];
            let (c, h, _) = candidate.into_parts();
            contact = c;
            thermal = h;
        }
        let storage = contact.stored_energy_j().unwrap();
        let deposited: f64 = thermal
            .cells()
            .values()
            .map(|c| c.mass_kg * c.specific_h_j_kg)
            .sum();
        near(deposited, generated, 1e-8);
        let error = [
            (force - exact_force).abs(),
            (storage - exact_storage).abs(),
            (generated - exact_heat).abs(),
            numerical,
        ];
        eprintln!("moving brush n={n} error={error:?}");
        errors.push(error);
        assert!(contact.parcels().len() <= n + 1);
    }
    for axis in [0, 1, 3] {
        for pair in errors.windows(2) {
            let ratio = pair[0][axis] / pair[1][axis];
            assert!((1.55..2.5).contains(&ratio), "axis {axis}: ratio {ratio}");
        }
    }
    // Integrated heat has an O(dt^2) envelope on this constant-coefficient
    // fixture; its observed ratio varies with the two yield times' location
    // inside a step. Do not turn that special cancellation into a general
    // second-order temporal claim for transported history or thermal feedback.
    for (i, n) in [16, 32, 64].into_iter().enumerate() {
        assert!(errors[i][2] <= 0.03 / (n * n) as f64);
        if i > 0 {
            assert!(errors[i][2] < errors[i - 1][2]);
        }
    }
    assert!(errors[2][0] < 0.008);
    assert!(errors[2][1] < 0.00012);
    assert!(errors[2][2] < 1e-5);
}

#[test]
fn g4_g5_late_failure_keeps_all_material_and_thermal_history_retryable() {
    let (thermal, contact) = fixture(&[-3.0, 0.0, 1.0, 4.0], 1.0, [0.0; 2]);
    let snapshot = contact.parcels().to_vec();
    let end = position(0.02, 0.0, 0.02);
    let reference = step(&contact, &thermal, 0.1, end, 1.0);
    let calls = Cell::new(0);
    let error = contact.step_candidate(
        &thermal,
        0.1,
        end,
        &work(),
        controls(),
        &vec![0.0; thermal.cells().len()],
        300.0,
        |t, r, _, _| {
            calls.set(calls.get() + 1);
            if calls.get() >= 4 {
                Err("late pinned-source refusal")
            } else {
                Ok(source(t, r, 1.0))
            }
        },
        || false,
    );
    assert!(matches!(
        error,
        Err(MovingPartialSlipError::Source("late pinned-source refusal"))
    ));
    assert_eq!(contact.parcels(), snapshot);
    assert_eq!(contact.time_s(), 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert!(thermal.cells().values().all(|c| c.specific_h_j_kg == 0.0));
    let retried = step(&contact, &thermal, 0.1, end, 1.0);
    assert_eq!(reference.contact().parcels(), retried.contact().parcels());
    assert_eq!(
        reference.report().heat.dissipated_j,
        retried.report().heat.dissipated_j
    );
    for (a, b) in reference
        .thermal()
        .cells()
        .values()
        .zip(retried.thermal().cells().values())
    {
        assert_eq!(a.specific_h_j_kg, b.specific_h_j_kg);
    }
    let mut cap = controls();
    cap.maximum_regions = 1;
    let limited = contact.step_candidate(
        &thermal,
        0.1,
        position(0.0, 0.0, 0.25),
        &work(),
        cap,
        &vec![0.0; thermal.cells().len()],
        300.0,
        |t, r, _, _| Ok::<_, Infallible>(source(t, r, 1.0)),
        || false,
    );
    assert!(matches!(limited, Err(MovingPartialSlipError::RegionBudget)));
    let polls = Cell::new(0);
    let cancelled = contact.step_candidate(
        &thermal,
        0.1,
        end,
        &work(),
        controls(),
        &vec![0.0; thermal.cells().len()],
        300.0,
        |t, r, _, _| Ok::<_, Infallible>(source(t, r, 1.0)),
        || {
            polls.set(polls.get() + 1);
            polls.get() >= 70
        },
    );
    assert!(matches!(
        cancelled,
        Err(MovingPartialSlipError::Thermal(DepositionError::Cancelled))
    ));
    let changed = contact.step_candidate(
        &thermal,
        0.1,
        end,
        &work(),
        controls(),
        &vec![0.0; thermal.cells().len()],
        300.0,
        |t, r, _, _| {
            let mut s = source(t, r, 1.0);
            for law in &mut s.laws {
                let mut p = law.parameters();
                p.tangential_stiffness_n_per_m *= 2.0;
                *law = PartialSlipLaw::new(law.model_id(), law.source_id(), p).unwrap();
            }
            Ok::<_, Infallible>(s)
        },
        || false,
    );
    assert!(matches!(changed, Err(MovingPartialSlipError::Invalid(_))));
    assert_eq!(contact.parcels(), snapshot);
}
