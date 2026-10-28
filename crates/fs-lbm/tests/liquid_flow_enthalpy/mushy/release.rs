//! Original latent heat and MAC momentum continue into actual native VOF.
//! Manufactured data isolate conservation and data-only material response.
use super::continuity::{chart, resolve, source};
use super::*;
use fs_lbm::conjugate::enthalpy::advection::feedback::{
    PhaseFlowTransport,
    release::{LiquidReleaseControls, LiquidReleaseError, StationaryConstraintRelease},
};
use fs_lbm::conjugate::{FvFlow, ViscousHeatControls, Voxel};
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::material::free_surface::enthalpy::conduction::FreeSurfaceConductionControls;
use fs_lbm::material::free_surface::enthalpy::{FreeSurfaceCellHeat, FreeSurfaceEnthalpyControls};
use fs_lbm::material::free_surface::{MaterialFreeSurfaceDomain3, plan_free_surface_from_material};
use fs_material::MaterialError;
use fs_material::liquid::surface::{LiquidCapillarity, resolve_liquid_surface_state};
use fs_material::phase::transport::ResolvedPhaseTransportState;

fn domain() -> VoxelDomain {
    VoxelDomain::from_fn(4, 4, 4, 1.0, |[x, y, z]| {
        if (1.0..3.0).contains(&x) && (1.0..3.0).contains(&y) && z < 2.0 {
            Voxel::Fluid
        } else {
            Voxel::Solid(0)
        }
    })
    .unwrap()
}
fn config() -> SimpleConfig {
    let mut config = SimpleConfig::new([FvBoundary::wall(); 6]);
    config.tolerance = 1e-12;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-4;
    config
}
fn viscous() -> ViscousHeatControls {
    ViscousHeatControls {
        maximum_reynolds: 0.1,
        power_tolerance_w: 1e-8,
    }
}
fn step() -> StepConfig {
    StepConfig {
        max_dt_s: 0.01,
        energy_tolerance_j: 1e-8,
    }
}
fn limits() -> LiquidReleaseControls {
    LiquidReleaseControls {
        mass_tolerance_kg: 1e-12,
        enthalpy_tolerance_j: 1e-9,
        momentum_tolerance_kg_m_s: 1e-10,
        angular_momentum_tolerance_kg_m2_s: 1e-10,
        maximum_kinetic_projection_deficit_j: 1e-3,
        kinetic_roundoff_tolerance_j: 1e-12,
    }
}
fn receiving(gravity: [f64; 3]) -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [4; 3],
        gravity_m_s2: gravity,
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 64,
    }
}
fn ids(domain: &VoxelDomain) -> Vec<Option<CellId>> {
    (0..domain.cell_count())
        .map(|c| domain.is_fluid(c).then_some(CellId(c as u64)))
        .collect()
}
fn initial<'a>(
    domain: &VoxelDomain,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    h: f64,
    watts: f64,
    heat_until_s: f64,
) -> DepositionSimulation<'a> {
    let boundary = |source_w| Boundary {
        ambient_k: 500.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w,
    };
    let mut events = Vec::new();
    for c in 0..domain.cell_count() {
        if !domain.is_fluid(c) {
            continue;
        }
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
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: start,
                    end_m: end,
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: boundary(watts),
            }),
        });
    }
    if watts != 0.0 {
        for id in ids(domain).into_iter().flatten() {
            events.push(TimedEvent {
                time_s: heat_until_s,
                event: Event::Boundary(id, boundary(0.0)),
            });
        }
    }
    DepositionSimulation::new(
        events,
        &[(MAT, curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 64,
            max_links: 128,
            max_events: 32,
            max_substeps: 256,
        },
    )
    .unwrap()
    .advance_to(0.0, step(), || Ok(()))
    .unwrap()
    .0
}
fn bound<'a>(
    thermal: &DepositionSimulation<'_>,
    domain: &'a VoxelDomain,
    ids: &'a [Option<CellId>],
    config: &'a SimpleConfig,
    gate: &'a CancelGate,
    card: &'a MaterialCard,
    mushy: &'a MushyPhaseLaw<'a>,
) -> MaterialFlowFeedback<
    'a,
    ResolvedPhaseTransportState,
    impl Clone
    + FnMut(usize, &DepositedCell, f64) -> Result<ResolvedPhaseTransportState, MaterialError>
    + 'a,
    (),
    PhaseFlowTransport,
> {
    MaterialFlowFeedback::try_new_phase(
        thermal,
        domain,
        FlowEnthalpyInput {
            cell_ids: ids,
            material_id: MAT,
            caloric: mushy.curve(),
            inlet_specific_h_j_kg: [None; 6],
            max_cells: 64,
            max_connections: 1024,
        },
        config,
        (),
        gate,
        256,
        mushy,
        move |_, cell: &DepositedCell, _| resolve(card, mushy, cell.specific_h_j_kg),
    )
    .unwrap()
    .with_unsteady_viscous_heating(viscous())
    .unwrap()
}
fn point_liquid(t: f64) -> QueryPoint {
    point(t, 100000.0).with("liquid", 1.0).unwrap()
}
fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&fs_exec::Cx<'_>) -> T) -> T {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, Cx, ExecMode, StreamKey};
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} versus {b}; tolerance {tolerance}"
    );
}

#[test]
fn actual_latent_heating_releases_into_native_gravity_flow_and_preserves_the_clock() {
    let card = source(1.0);
    let curve = chart(&card);
    let mushy = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = domain();
    let ids = ids(&domain);
    let config = config();
    let gate = CancelGate::new_clock_free();
    let original = initial(&domain, &curve, -100.0, 2000.0, 0.11);
    let history = bound(&original, &domain, &ids, &config, &gate, &card, &mushy);
    let (plateau, history, _) = original
        .advance_observed(&history, 0.075, step(), || Ok(()))
        .unwrap();
    for cell in plateau.cells().values() {
        close(cell.specific_h_j_kg, 50.0, 1e-10);
        close(plateau.temperature_k(cell.id).unwrap(), 500.0, 0.0);
    }
    let point = point_liquid(510.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 1.0, 0.01, 0.2).unwrap();
        assert!(matches!(
            history.release_fully_liquid(
                &plateau,
                scaling,
                receiving([0.2, 0.0, 0.0]),
                &law,
                &curve,
                StationaryConstraintRelease::Withdraw,
                limits(),
                FreeSurfaceEnthalpyControls::default(),
                cx
            ),
            Err(LiquidReleaseError::Invalid(_))
        ));
        let (molten, history, _) = plateau
            .advance_observed(&history, 0.11, step(), || Ok(()))
            .unwrap();
        close(molten.energy_j() - original.energy_j(), 8.0 * 220.0, 1e-9);
        assert_eq!(history.unsteady_time_s(), Some(0.11));
        let preserved = (
            molten.cells().clone(),
            history.unsteady_flow().unwrap().clone(),
        );
        let (mut released, report) = history
            .release_fully_liquid(
                &molten,
                scaling,
                receiving([0.2, 0.0, 0.0]),
                &law,
                &curve,
                StationaryConstraintRelease::Withdraw,
                limits(),
                FreeSurfaceEnthalpyControls::default(),
                cx,
            )
            .unwrap();
        close(released.physical_time_s(), 0.11, 0.0);
        close(report.source.mass_kg, 8.0, 0.0);
        close(report.target_enthalpy_j, molten.energy_j(), 1e-9);
        close(report.kinetic_projection_deficit_j, 0.0, 1e-20);
        let pristine = released.flow().lattice().grid().state_bits();
        let mut coarse_clock = released.clone().with_time_origin(1e300).unwrap();
        assert!(coarse_clock.step(&gate).is_err());
        assert_eq!(coarse_clock.accepted_steps(), 0);
        assert_eq!(coarse_clock.physical_time_s(), 1e300);
        assert_eq!(coarse_clock.flow().lattice().grid().state_bits(), pristine);
        close(
            coarse_clock.total_enthalpy_j().unwrap(),
            report.target_enthalpy_j,
            0.0,
        );
        for _ in 0..8 {
            let report = released.step(&gate).unwrap();
            close(report.mass_after_kg, 8.0, 1e-11);
            close(report.enthalpy_after_j, molten.energy_j(), 1e-7);
            close(report.external_heat_j, 0.0, 0.0);
        }
        assert_eq!(released.accepted_steps(), 8);
        close(released.physical_time_s(), 0.19, 1e-15);
        assert!(released.clone().with_time_origin(0.0).is_err());
        assert!(
            released
                .flow()
                .sample_at_si(1, 1, 0)
                .unwrap()
                .unwrap()
                .velocity_m_s[0]
                > 0.0
        );
        assert_eq!(molten.cells(), &preserved.0);
        assert_eq!(history.unsteady_flow().unwrap(), &preserved.1);
        let before = (
            released.physical_time_s(),
            released.total_enthalpy_j().unwrap(),
            released.flow().lattice().grid().state_bits(),
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(released.step(&cancelled).is_err());
        assert_eq!(released.accepted_steps(), 8);
        close(released.physical_time_s(), before.0, 0.0);
        close(released.total_enthalpy_j().unwrap(), before.1, 0.0);
        assert_eq!(released.flow().lattice().grid().state_bits(), before.2);
        let impossible = [FreeSurfaceCellHeat {
            coordinates: [1, 1, 0],
            energy_j: -1e9,
        }];
        assert!(released.step_with_heat(&impossible, &gate).is_err());
        assert_eq!(released.accepted_steps(), 8);
        close(released.physical_time_s(), before.0, 0.0);
        assert_eq!(released.flow().lattice().grid().state_bits(), before.2);
    });
}

#[test]
fn nonzero_native_circulation_retains_both_momenta_and_separate_projection_loss() {
    let card = source(100.0);
    let curve = chart(&card);
    let mushy = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = domain();
    let ids = ids(&domain);
    let config = config();
    let gate = CancelGate::new_clock_free();
    let original = initial(&domain, &curve, 200.0, 0.0, 0.0);
    let force = (0..domain.cell_count())
        .map(|c| {
            if !domain.is_fluid(c) {
                return [0.0; 3];
            }
            let [x, y, z] = domain.coords(c);
            let p = domain.center(x, y, z);
            [-(p[1] - 2.0), p[0] - 2.0, 0.0]
        })
        .collect();
    let history = bound(&original, &domain, &ids, &config, &gate, &card, &mushy)
        .with_unsteady_body_force_heating(force, viscous())
        .unwrap();
    let (thermal, history, _) = original
        .advance_observed(&history, 0.01, step(), || Ok(()))
        .unwrap();
    let flow: &FvFlow = history.unsteady_flow().unwrap();
    let speed = flow.field.flux_x(2, 1, 0);
    assert!(speed > 0.0);
    // Two identical 2x2 closed circulation layers: 8 moving MAC faces of speed
    // s give K=4s², while 8 cell velocities (+/-s/2,+/-s/2,0) give K=2s².
    // Their angular momentum is 4s about z in BOTH representations; P is zero.
    for z in 0..2 {
        close(flow.field.flux_x(2, 1, z), speed, 1e-10);
        close(flow.field.flux_x(2, 2, z), -speed, 1e-10);
        close(flow.field.flux_y(1, 2, z), -speed, 1e-10);
        close(flow.field.flux_y(2, 2, z), speed, 1e-10);
    }
    let point = point_liquid(550.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, 1.0, 0.0001, 0.2).unwrap();
        let unchanged = (
            thermal.cells().clone(),
            history.unsteady_flow().unwrap().clone(),
        );
        let mut strict = limits();
        strict.maximum_kinetic_projection_deficit_j = 0.0;
        assert!(
            history
                .release_fully_liquid(
                    &thermal,
                    scaling,
                    receiving([0.0; 3]),
                    &law,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    strict,
                    FreeSurfaceEnthalpyControls::default(),
                    cx
                )
                .is_err()
        );
        let (mut released, report) = history
            .release_fully_liquid(
                &thermal,
                scaling,
                receiving([0.0; 3]),
                &law,
                &curve,
                StationaryConstraintRelease::Withdraw,
                limits(),
                FreeSurfaceEnthalpyControls::default(),
                cx,
            )
            .unwrap();
        close(report.source.kinetic_energy_j, 4.0 * speed * speed, 1e-10);
        close(report.target.kinetic_energy_j, 2.0 * speed * speed, 1e-10);
        close(
            report.kinetic_projection_deficit_j,
            2.0 * speed * speed,
            1e-10,
        );
        for state in [report.source, report.target] {
            for value in state.momentum_kg_m_s {
                close(value, 0.0, 1e-10);
            }
            close(state.angular_momentum_kg_m2_s[2], 4.0 * speed, 1e-10);
        }
        close(report.target_enthalpy_j, thermal.energy_j(), 1e-9);
        let next = released.step(&gate).unwrap();
        close(next.external_heat_j, 0.0, 0.0);
        close(next.viscous_heat_j, 0.0, 0.0);
        close(next.enthalpy_after_j, thermal.energy_j(), 1e-7);
        let heated = released.step_with_viscous_heating(&[], &gate).unwrap();
        assert!(heated.flow.viscous_heat_j > 0.0);
        close(
            heated.flow.enthalpy_after_j - next.enthalpy_after_j,
            heated.flow.viscous_heat_j,
            1e-8,
        );
        let conducted = released
            .step_with_conduction_and_viscous_heating(
                &[],
                &[],
                FreeSurfaceConductionControls {
                    maximum_faces: 64,
                    maximum_substeps: 64,
                    maximum_dt_s: 0.0001,
                    energy_tolerance_j: 1e-8,
                },
                &gate,
            )
            .unwrap();
        close(
            conducted.conduction.flow.enthalpy_after_j - heated.flow.enthalpy_after_j,
            conducted.conduction.flow.viscous_heat_j,
            1e-8,
        );
        close(conducted.conduction.flow.external_heat_j, 0.0, 0.0);
        assert_eq!(released.accepted_steps(), 3);
        close(released.physical_time_s(), 0.0103, 1e-15);
        assert_eq!(thermal.cells(), &unchanged.0);
        assert_eq!(history.unsteady_flow().unwrap(), &unchanged.1);
    });
}

#[test]
fn data_only_viscosity_changes_native_motion_after_the_same_melt() {
    let mut velocities = Vec::new();
    for mu in [1.0, 4.0] {
        let card = source(mu);
        let curve = chart(&card);
        let mushy = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
        let domain = domain();
        let ids = ids(&domain);
        let config = config();
        let gate = CancelGate::new_clock_free();
        let original = initial(&domain, &curve, -100.0, 2000.0, 0.11);
        let history = bound(&original, &domain, &ids, &config, &gate, &card, &mushy);
        let (molten, history, _) = original
            .advance_observed(&history, 0.11, step(), || Ok(()))
            .unwrap();
        let point = point_liquid(510.0);
        let law = LiquidTemperatureLaw::try_new(
            &card,
            &point,
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let surface = resolve_liquid_surface_state(
            &card,
            &point,
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Neglected,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        with_cx(&gate, |cx| {
            let scaling = plan_free_surface_from_material(&surface, 1.0, 0.01, 0.2).unwrap();
            let (mut released, _report) = history
                .release_fully_liquid(
                    &molten,
                    scaling,
                    receiving([0.2, 0.0, 0.0]),
                    &law,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    limits(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx,
                )
                .unwrap();
            for _ in 0..12 {
                released.step(&gate).unwrap();
            }
            close(
                released
                    .material_at([1, 1, 0])
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                mu,
                0.0,
            );
            close(released.total_enthalpy_j().unwrap(), 8.0 * 120.0, 1e-7);
            velocities.push(
                released
                    .flow()
                    .sample_at_si(1, 1, 0)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s[0],
            );
        });
    }
    assert!(velocities.iter().all(|v| *v > 0.0));
    assert!(
        (velocities[0] - velocities[1]).abs() > 1e-8,
        "native velocities {velocities:?}"
    );
}

#[test]
fn source_clock_caloric_identity_geometry_and_cancellation_refuse_without_publication() {
    let card = source(1.0);
    let other_card = source(2.0);
    let curve = chart(&card);
    let copied_curve = chart(&card);
    let mushy = MushyPhaseLaw::try_new_with_immobile_solid(&curve, 1.0, 1.0, 1.0).unwrap();
    let domain = domain();
    let ids = ids(&domain);
    let config = config();
    let gate = CancelGate::new_clock_free();
    let original = initial(&domain, &curve, 200.0, 0.0, 0.0);
    let history = bound(&original, &domain, &ids, &config, &gate, &card, &mushy);
    let wrong_clock = original.advance_to(0.01, step(), || Ok(())).unwrap().0;
    let point = point_liquid(550.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let other_law = LiquidTemperatureLaw::try_new(
        &other_card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let before = (
        original.cells().clone(),
        history.unsteady_flow().unwrap().clone(),
    );
    let scaling = plan_free_surface_from_material(&surface, 1.0, 0.01, 0.2).unwrap();
    with_cx(&gate, |cx| {
        assert!(matches!(
            history.release_fully_liquid(
                &wrong_clock,
                scaling,
                receiving([0.0; 3]),
                &law,
                &curve,
                StationaryConstraintRelease::Withdraw,
                limits(),
                FreeSurfaceEnthalpyControls::default(),
                cx
            ),
            Err(LiquidReleaseError::Source(_))
        ));
        assert!(
            history
                .release_fully_liquid(
                    &original,
                    scaling,
                    receiving([0.0; 3]),
                    &law,
                    &copied_curve,
                    StationaryConstraintRelease::Withdraw,
                    limits(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx
                )
                .is_err()
        );
        assert!(
            history
                .release_fully_liquid(
                    &original,
                    scaling,
                    receiving([0.0; 3]),
                    &other_law,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    limits(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx
                )
                .is_err()
        );
        let mut moved = receiving([0.0; 3]);
        moved.cells[0] = 8;
        assert!(
            history
                .release_fully_liquid(
                    &original,
                    scaling,
                    moved,
                    &law,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    limits(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx
                )
                .is_err()
        );
        let mut budget = FreeSurfaceEnthalpyControls::default();
        budget.max_cells = 1;
        assert!(
            history
                .release_fully_liquid(
                    &original,
                    scaling,
                    receiving([0.0; 3]),
                    &law,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    limits(),
                    budget,
                    cx
                )
                .is_err()
        );
    });
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    with_cx(&cancelled, |cx| {
        assert!(
            history
                .release_fully_liquid(
                    &original,
                    scaling,
                    receiving([0.0; 3]),
                    &law,
                    &curve,
                    StationaryConstraintRelease::Withdraw,
                    limits(),
                    FreeSurfaceEnthalpyControls::default(),
                    cx
                )
                .is_err()
        );
    });
    assert_eq!(original.cells(), &before.0);
    assert_eq!(history.unsteady_flow().unwrap(), &before.1);
    assert_eq!(history.unsteady_time_s(), Some(0.0));
}
