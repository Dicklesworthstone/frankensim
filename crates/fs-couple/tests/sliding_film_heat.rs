//! Combined moving-wall dissipation heats the original finite thermal owners.
//! Synthetic constant layer viscosities and caloric charts are physical oracles,
//! not measured material claims. No viscosity/temperature feedback is inferred.

use fs_blake3::ContentHash;
use fs_couple::render::plate::impact::material::heating::network::ReservoirCell;
use fs_couple::render::plate::impact::material::heating::{
    CaloricReservoir, HeatedImpactSystem, SlidingFilmLayerHeatRoute, SqueezeFilmLayerHeatRoute,
};
use fs_couple::render::plate::impact::squeeze::sliding::{ShearedFilm, WallVelocityPorts};
use fs_couple::render::plate::impact::squeeze::{
    FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm,
};
use fs_couple::render::plate::impact::{ImpactBody, ImpactConfig, ImpactSystem};
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use fs_tribo::resistive_film::stratified::{GapViscosityProfile, StratifiedFilm, ViscosityLayer};

const M1: f64 = 0.2;
const M2: f64 = 0.4;
const GAP: f64 = 0.001;
const WIDTH: f64 = 0.02;
const LENGTH: f64 = 0.03;
const DT: f64 = 0.01;
const LAUNCH: f64 = 0.01;
const STORE_MASS: f64 = 0.0006;
const STORE_VOLUME: f64 = 0.0000006;
const STORE_AREA: f64 = 0.0012;

fn bare() -> ImpactSystem {
    let (a, _) = ImpactBody::free_mass(M1, 0.0, LAUNCH).unwrap();
    let (b, _) = ImpactBody::free_mass(M2, 0.0, 0.0).unwrap();
    ImpactSystem::new(
        vec![a, b],
        vec![],
        vec![],
        vec![],
        ImpactConfig {
            dt_s: DT,
            max_steps: 256,
            maximum_energy_j: 1.0,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8,
            maximum_generalized_force: 1.0,
        },
    )
    .unwrap()
}
fn profile() -> GapViscosityProfile {
    GapViscosityProfile::new(vec![
        ViscosityLayer {
            upper_fraction: 0.5,
            dynamic_viscosity_pa_s: 0.01,
        },
        ViscosityLayer {
            upper_fraction: 1.0,
            dynamic_viscosity_pa_s: 0.04,
        },
    ])
    .unwrap()
}
fn graph() -> ResistiveFilm {
    let gap = GapPort {
        reference_m: GAP,
        closure: vec![0.0; 2],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: WIDTH * LENGTH,
            gap: gap.clone(),
        }],
        (0..2)
            .map(|_| FilmChannel {
                from: 0,
                to: None,
                width_m: WIDTH,
                length_m: LENGTH,
                gap: gap.clone(),
            })
            .collect(),
        2,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 0.01,
            maximum_pressure_pa: 1e6,
        },
    )
    .unwrap()
}
fn mechanics() -> ImpactSystem {
    let rows = [1.0, -1.0]
        .into_iter()
        .map(|sign| WallVelocityPorts {
            lower: vec![sign / M1.sqrt(), 0.0],
            upper: vec![0.0, sign / M2.sqrt()],
        })
        .collect();
    let film = ShearedFilm::from_stratified(
        StratifiedFilm::new(&graph(), vec![profile(); 2]).unwrap(),
        rows,
    )
    .unwrap();
    bare().with_sliding_film(film).unwrap()
}
fn curve(upper_h: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([73; 32]),
        SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: upper_h,
                temperature_k: 300.0 + upper_h / 0.2,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap()
}
fn stores(c: &EquilibriumEnthalpyPhaseCurve) -> Vec<CaloricReservoir<'_>> {
    ["lower", "upper"]
        .into_iter()
        .map(|key| {
            CaloricReservoir::new(
                key.into(),
                c,
                c.state_at_specific_enthalpy(0.0).unwrap(),
                STORE_MASS,
                STORE_VOLUME,
                STORE_AREA,
            )
            .unwrap()
        })
        .collect()
}
fn routes() -> Vec<SlidingFilmLayerHeatRoute> {
    (0..2)
        .map(|channel| SlidingFilmLayerHeatRoute {
            film: 0,
            channel,
            layer_destinations: vec!["lower".into(), "upper".into()],
        })
        .collect()
}
fn heated(c: &EquilibriumEnthalpyPhaseCurve) -> HeatedImpactSystem<'_> {
    HeatedImpactSystem::new_with_sliding_losses(
        mechanics().prepare_analytic().unwrap(),
        stores(c),
        vec![],
        vec![],
        vec![],
        vec![],
        routes(),
        vec![],
        1e-12,
    )
    .unwrap()
}
fn momentum(x: &[f64]) -> f64 {
    M1.sqrt() * x[1] + M2.sqrt() * x[3]
}
fn relative(x: &[f64]) -> f64 {
    x[1] / M1.sqrt() - x[3] / M2.sqrt()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn actual_layer_heat_closes_native_momentum_dg_work_and_enthalpy_for_both_newton_images() {
    // G1/G3: the opposed open channels have p=0. Pure Couette layer work is
    // proportional to thickness/mu: (0.5/.01):(0.5/.04) = 0.8:0.2.
    // The same harmonic viscosity sets the independently known body drag.
    let c = curve(20.0);
    let harmonic = 1.0 / (0.5 / 0.01 + 0.5 / 0.04);
    let drag = 2.0 * WIDTH * LENGTH * harmonic / GAP;
    let rate = drag * (1.0 / M1 + 1.0 / M2);
    let ratio = (1.0 - 0.5 * DT * rate) / (1.0 + 0.5 * DT * rate);
    for analytic in [false, true] {
        let mut system = heated(&c);
        system.set_analytic_newton(analytic);
        let mut native = mechanics().prepare_analytic().unwrap();
        native.set_analytic_newton(analytic);
        let initial = system.mechanics().stored_energy_j();
        let initial_momentum = momentum(system.mechanics().state());
        let mut expected_relative = LAUNCH;
        let mut total = 0.0;
        for tick in 1..=64 {
            let before = relative(system.mechanics().state());
            let frame = system
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            let plain = native
                .step(&[0.0; 2], &CancelGate::new_clock_free())
                .unwrap();
            let after = relative(system.mechanics().state());
            expected_relative *= ratio;
            let expected_heat = DT * drag * (0.5 * (before + after)).powi(2);
            assert_eq!(system.mechanics().state(), native.state());
            assert_eq!(
                frame.mechanical.dissipated_energy_j,
                plain.dissipated_energy_j
            );
            near(after, expected_relative, 2e-11);
            near(
                momentum(system.mechanics().state()),
                initial_momentum,
                2e-11,
            );
            near(frame.deposited_heat_j, expected_heat, 2e-12);
            near(frame.remaining_mechanical_loss_j, 0.0, 2e-12);
            near(frame.thermal_energy_change_j, frame.deposited_heat_j, 1e-12);
            assert_eq!(frame.mechanical.sample, tick);
            near(frame.mechanical.time_s, tick as f64 * DT, 1e-14);
            total += expected_heat;
        }
        near(system.reservoirs()[0].received_heat_j(), 0.8 * total, 2e-12);
        near(system.reservoirs()[1].received_heat_j(), 0.2 * total, 2e-12);
        assert!(
            system.reservoirs()[0].state().temperature_k()
                > system.reservoirs()[1].state().temperature_k()
        );
        near(
            system.mechanics().stored_energy_j()
                + system
                    .reservoirs()
                    .iter()
                    .map(CaloricReservoir::storage_change_j)
                    .sum::<f64>(),
            initial,
            2e-11,
        );
        assert!(system.layer_heat_profile(0, 0).is_none());
    }
}

#[test]
fn selected_channels_keep_unrouted_work_and_normal_indices_separate() {
    // One of two identical sliding channels is routed. A zero-work normal
    // pressure island deliberately has the same film/channel index.
    let c = curve(20.0);
    let normal = StratifiedFilm::new(&graph(), vec![profile(); 2]).unwrap();
    let mechanical = mechanics()
        .with_squeeze_film(normal.film().clone())
        .unwrap()
        .prepare_analytic()
        .unwrap();
    let mut selected = routes();
    selected.truncate(1);
    let mut system = HeatedImpactSystem::new_with_sliding_losses(
        mechanical,
        stores(&c),
        vec![],
        vec![],
        vec![],
        vec![SqueezeFilmLayerHeatRoute {
            film: 0,
            channel: 0,
            profile: profile(),
            layer_destinations: vec!["lower".into(), "upper".into()],
        }],
        selected,
        vec![],
        1e-12,
    )
    .unwrap();
    assert_eq!(system.layer_heat_profile(0, 0), Some(&profile()));
    let initial = system.mechanics().stored_energy_j();
    let mut outside = 0.0;
    for _ in 0..32 {
        let frame = system
            .step(&[0.0; 2], &CancelGate::new_clock_free())
            .unwrap();
        near(
            frame.deposited_heat_j,
            0.5 * frame.mechanical.dissipated_energy_j,
            2e-12,
        );
        near(
            frame.remaining_mechanical_loss_j,
            frame.deposited_heat_j,
            2e-12,
        );
        outside += frame.remaining_mechanical_loss_j;
    }
    assert!(outside > 0.0);
    near(
        system.mechanics().stored_energy_j()
            + system
                .reservoirs()
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>()
            + outside,
        initial,
        2e-11,
    );
}

fn network(c: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for id in [1, 2] {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(1),
                mass_kg: STORE_MASS,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64 + 0.01, 0.0, 0.0],
                    volume_m3: STORE_VOLUME,
                    surface_area_m2: STORE_AREA,
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
            id: 7,
            a: CellId(1),
            b: CellId(2),
            area_m2: 1e-6,
            conductance_w_k: 0.001,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), c as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT / 4.0,
            energy_tolerance_j: 1e-12,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn bindings() -> [ReservoirCell; 2] {
    [
        ReservoirCell {
            reservoir: "lower".into(),
            cell: CellId(1),
        },
        ReservoirCell {
            reservoir: "upper".into(),
            cell: CellId(2),
        },
    ]
}

#[test]
fn sliding_heat_and_conduction_use_one_original_network_transaction() {
    let c = curve(20.0);
    let mut system = heated(&c);
    let mut thermal = network(&c);
    let mut observer = ();
    let initial = system.mechanics().stored_energy_j();
    let mut total = 0.0;
    let mut substeps = 0;
    for tick in 1..=16 {
        let frame = system
            .step_with_thermal_network(
                &mut thermal,
                &mut observer,
                &bindings(),
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-12,
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None),
            )
            .unwrap();
        substeps += frame.network.substeps;
        total += frame.impact.thermal.impact.deposited_heat_j;
        near(frame.external_heat_j, 0.0, 2e-12);
        near(thermal.time_s(), tick as f64 * DT, 1e-14);
        for (i, binding) in bindings().iter().enumerate() {
            assert_eq!(
                system.reservoirs()[i].state().specific_enthalpy_j_kg(),
                thermal.cells()[&binding.cell].specific_h_j_kg
            );
        }
        assert!(
            thermal
                .cells()
                .values()
                .all(|cell| cell.boundary.source_w == 0.0)
        );
    }
    assert!(substeps >= 64);
    near(system.reservoirs()[0].received_heat_j(), 0.8 * total, 2e-12);
    near(system.reservoirs()[1].received_heat_j(), 0.2 * total, 2e-12);
    assert!(system.reservoirs()[0].boundary_heat_j() < 0.0);
    assert!(system.reservoirs()[1].boundary_heat_j() > 0.0);
    near(
        system.mechanics().stored_energy_j()
            + system
                .reservoirs()
                .iter()
                .map(CaloricReservoir::storage_change_j)
                .sum::<f64>(),
        initial,
        2e-11,
    );
}

#[test]
fn invalid_coverage_cancellation_and_late_chart_failure_publish_nothing() {
    let c = curve(20.0);
    let first = routes().remove(0);
    for invalid in [
        vec![SlidingFilmLayerHeatRoute {
            layer_destinations: vec!["lower".into()],
            ..first.clone()
        }],
        vec![SlidingFilmLayerHeatRoute {
            layer_destinations: vec!["lower".into(), "missing".into()],
            ..first.clone()
        }],
        vec![SlidingFilmLayerHeatRoute {
            channel: 2,
            ..first.clone()
        }],
        vec![SlidingFilmLayerHeatRoute {
            film: 1,
            ..first.clone()
        }],
        vec![first.clone(), first.clone()],
    ] {
        assert!(
            HeatedImpactSystem::new_with_sliding_losses(
                mechanics().prepare_analytic().unwrap(),
                stores(&c),
                vec![],
                vec![],
                vec![],
                vec![],
                invalid,
                vec![],
                1e-12
            )
            .is_err()
        );
    }
    // Several physical layers may share one explicitly named finite store.
    assert!(
        HeatedImpactSystem::new_with_sliding_losses(
            mechanics().prepare_analytic().unwrap(),
            stores(&c),
            vec![],
            vec![],
            vec![],
            vec![],
            vec![SlidingFilmLayerHeatRoute {
                layer_destinations: vec!["lower".into(), "lower".into()],
                ..first
            }],
            vec![],
            1e-12
        )
        .is_ok()
    );

    let mut system = heated(&c);
    let before = system.mechanics().state().to_vec();
    let original_stores = system.reservoirs().to_vec();
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(system.step(&[0.0; 2], &stopped).is_err());
    assert_eq!(system.mechanics().state(), before);
    assert_eq!(system.mechanics().samples(), 0);
    assert_eq!(system.reservoirs(), original_stores);
    let retried = system
        .step(&[0.0; 2], &CancelGate::new_clock_free())
        .unwrap();
    let mut clean = heated(&c);
    let expected = clean
        .step(&[0.0; 2], &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(system.mechanics().state(), clean.mechanics().state());
    assert_eq!(system.reservoirs(), clean.reservoirs());
    assert_eq!(retried.deposited_heat_j, expected.deposited_heat_j);

    let narrow = curve(1e-12);
    let mut late = heated(&narrow);
    let motion = late.mechanics().state().to_vec();
    let thermal_stores = late.reservoirs().to_vec();
    assert!(late.step(&[0.0; 2], &CancelGate::new_clock_free()).is_err());
    assert_eq!(late.mechanics().state(), motion);
    assert_eq!(late.mechanics().samples(), 0);
    assert_eq!(late.reservoirs(), thermal_stores);

    let mut joint = heated(&narrow);
    let mut thermal = network(&narrow);
    let mut observer = ();
    assert!(
        joint
            .step_with_thermal_network(
                &mut thermal,
                &mut observer,
                &bindings(),
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-12
                },
                &[0.0; 2],
                &CancelGate::new_clock_free(),
                |_, _| Ok(None)
            )
            .is_err()
    );
    assert_eq!(joint.mechanics().state(), motion);
    assert_eq!(joint.mechanics().samples(), 0);
    assert_eq!(joint.reservoirs(), thermal_stores);
    assert_eq!(thermal.time_s(), 0.0);
    assert!(
        thermal
            .cells()
            .values()
            .all(|cell| cell.specific_h_j_kg == 0.0)
    );
    assert!(
        thermal
            .cells()
            .values()
            .all(|cell| cell.boundary.source_w == 0.0)
    );
}
