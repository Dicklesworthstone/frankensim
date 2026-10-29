//! Original Maxwell work reaches finite enthalpy and the existing heat network.
//! Synthetic constitutive data qualify coupling, not a measured material.
use fs_blake3::ContentHash;
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::material::heating::network::ReservoirCell;
use fs_couple::render::plate::impact::material::heating::{
    CaloricReservoir, HeatedImpactSystem, RelaxationHeatRoute,
};
use fs_couple::render::plate::impact::relaxation::InitialMemory;
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, ImpactConfig, ImpactSystem};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_phs::RelaxationBranch;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

const DT: f64 = 0.003;
const MASS: f64 = 0.01;
const TAU: [f64; 2] = [0.2, 0.07];

fn mechanics(damping: f64) -> ImpactSystem {
    let body = ImpactBody {
        potential: BodyPotential::Linear(vec![3.0]),
        initial: vec![ModalAcousticState {
            displacement_m_sqrt_kg: 0.3,
            velocity_m_sqrt_kg_per_s: 0.4,
        }],
        damping_per_s: vec![damping],
    };
    // The original contact owner supplies H_contact=5*[q]_+^4. This exercises
    // a genuine nonlinear Gonzalez correction in every material coordinate.
    let contact = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        20.0,
        3.0,
        "synthetic quartic obstacle".into(),
    )
    .unwrap();
    ImpactSystem::new(
        vec![body],
        vec![contact],
        vec![],
        vec![],
        ImpactConfig {
            dt_s: DT,
            max_steps: 1000,
            maximum_energy_j: 20.0,
            energy_absolute_tolerance_j: 1e-11,
            energy_relative_tolerance: 1e-8,
            maximum_generalized_force: 10.0,
        },
    )
    .unwrap()
    .with_relaxation_branches(
        vec![
            RelaxationBranch {
                projection: vec![1.0, 0.0],
                stiffness: 4.0,
                relaxation_time_s: TAU[0],
            },
            RelaxationBranch {
                projection: vec![-1.0, 0.0],
                stiffness: 1.0,
                relaxation_time_s: TAU[1],
            },
        ],
        InitialMemory::Unrelaxed,
        2,
    )
    .unwrap()
}

fn curve(upper_h: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([91; 32]),
        SolidLiquidPhase::Solid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: upper_h,
                temperature_k: 300.0 + upper_h / 10.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
        ],
    )
    .unwrap()
}
fn stores(c: &EquilibriumEnthalpyPhaseCurve) -> Vec<CaloricReservoir<'_>> {
    ["body", "support"]
        .map(|key| {
            CaloricReservoir::new(
                key.into(),
                c,
                c.state_at_specific_enthalpy(0.0).unwrap(),
                MASS,
                1e-5,
                0.001,
            )
            .unwrap()
        })
        .to_vec()
}
fn route(branch: usize, to: &str) -> RelaxationHeatRoute {
    RelaxationHeatRoute {
        branch,
        to: to.into(),
    }
}
fn heated(
    c: &EquilibriumEnthalpyPhaseCurve,
    damping: f64,
    routes: Vec<RelaxationHeatRoute>,
) -> HeatedImpactSystem<'_> {
    HeatedImpactSystem::new_with_losses(
        mechanics(damping).prepare_analytic().unwrap(),
        stores(c),
        vec![],
        vec![],
        routes,
        1e-10,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} within {tolerance:e}"
    );
}

// Independent polynomial potential and gradient, without an fs-phs storage,
// stepper or dissipation helper. The two arm efforts include the correction
// from the ORIGINAL nonlinear contact, even though their own laws are linear.
fn independent_heat(a: &[f64], b: &[f64]) -> ([f64; 2], [f64; 2]) {
    let midpoint: [f64; 4] = std::array::from_fn(|i| 0.5 * (a[i] + b[i]));
    let delta: [f64; 4] = std::array::from_fn(|i| b[i] - a[i]);
    assert!(a[0] > 0.0 && b[0] > 0.0);
    let remainder = 5.0 * midpoint[0] * delta[0].powi(3);
    let correction = remainder / delta.iter().map(|d| d * d).sum::<f64>();
    let uncorrected = [midpoint[2] - 2.0 * midpoint[0], midpoint[3] + midpoint[0]];
    let actual = std::array::from_fn(|i| {
        let effort = uncorrected[i] + correction * delta[2 + i];
        DT * effort * effort / TAU[i]
    });
    (actual, uncorrected.map(|e| e * e))
}

#[test]
fn exact_material_work_heats_distinct_stores_without_changing_nonlinear_motion() {
    let c = curve(1000.0);
    let gate = CancelGate::new_clock_free();
    for analytic in [false, true] {
        let mut model = heated(&c, 0.4, vec![route(0, "body"), route(1, "support")]);
        model.set_analytic_newton(analytic);
        let mut plain = mechanics(0.4).prepare().unwrap();
        plain.set_analytic_newton(analytic);
        let initial = plain.stored_energy_j();
        let mut external_loss = 0.0;
        let mut correction_witness = 0.0;
        for _ in 0..80 {
            let before = plain.state().to_vec();
            let inventory: Vec<_> = model
                .reservoirs()
                .iter()
                .map(|s| s.received_heat_j())
                .collect();
            let expected = plain.step(&[0.0], &gate).unwrap();
            let frame = model.step(&[0.0], &gate).unwrap();
            assert_eq!(plain.state(), model.mechanics().state());
            assert_eq!(
                expected.dissipated_energy_j,
                frame.mechanical.dissipated_energy_j
            );
            let (heat, midpoint_power) = independent_heat(&before, plain.state());
            for i in 0..2 {
                near(
                    model.reservoirs()[i].received_heat_j() - inventory[i],
                    heat[i],
                    2e-14,
                );
                correction_witness += (heat[i] - DT * midpoint_power[i] / TAU[i]).abs();
            }
            near(frame.deposited_heat_j, heat.iter().sum(), 2e-14);
            near(frame.combined_energy_residual_j, 0.0, 1e-10);
            external_loss += frame.remaining_mechanical_loss_j;
        }
        assert!(
            correction_witness > 1e-10,
            "test must distinguish corrected effort from midpoint power"
        );
        assert!(
            external_loss > 1e-3,
            "unselected viscous loss must remain outside the thermal stores"
        );
        let thermal: f64 = model
            .reservoirs()
            .iter()
            .map(|s| s.storage_change_j())
            .sum();
        near(
            model.mechanics().stored_energy_j() + thermal + external_loss,
            initial,
            1e-9,
        );
        for store in model.reservoirs() {
            assert!(store.state().temperature_k() > 300.0);
            near(
                store.state().temperature_k() - 300.0,
                store.received_heat_j() / (MASS * 10.0),
                1e-12,
            );
        }
    }
}

#[test]
fn selected_branches_are_canonical_and_retained_parts_do_not_replay_old_heat() {
    let c = curve(1000.0);
    let gate = CancelGate::new_clock_free();
    let routes = || vec![route(0, "body"), route(1, "body")];
    let mut a = heated(&c, 0.0, routes());
    let mut b = heated(&c, 0.0, vec![route(1, "body"), route(0, "body")]);
    let mut partial = heated(&c, 0.0, vec![route(1, "body")]);
    for _ in 0..30 {
        let fa = a.step(&[0.0], &gate).unwrap();
        b.step(&[0.0], &gate).unwrap();
        let fp = partial.step(&[0.0], &gate).unwrap();
        near(
            fa.deposited_heat_j,
            fa.mechanical.dissipated_energy_j,
            2e-15,
        );
        assert!(fp.remaining_mechanical_loss_j > 0.0);
        assert!(fp.deposited_heat_j < fa.deposited_heat_j);
        assert_eq!(a.mechanics().state(), partial.mechanics().state());
    }
    assert_eq!(a.reservoirs(), b.reservoirs());
    assert_eq!(a.reservoirs()[1].received_heat_j(), 0.0);
    let (mechanics, stores) = a.into_parts();
    let mut a =
        HeatedImpactSystem::new_with_losses(mechanics, stores, vec![], vec![], routes(), 1e-10)
            .unwrap();
    for _ in 0..20 {
        a.step(&[0.0], &gate).unwrap();
        b.step(&[0.0], &gate).unwrap();
    }
    assert_eq!(a.reservoirs(), b.reservoirs());
    assert_eq!(a.mechanics().state(), b.mechanics().state());
}

#[test]
fn route_caloric_and_solver_refusals_preserve_material_history_and_all_stores() {
    let c = curve(1000.0);
    for routes in [
        vec![route(0, "body"), route(0, "support")],
        vec![route(2, "body")],
        vec![route(0, "missing")],
    ] {
        assert!(
            HeatedImpactSystem::new_with_losses(
                mechanics(0.0).prepare_analytic().unwrap(),
                stores(&c),
                vec![],
                vec![],
                routes,
                1e-10
            )
            .is_err()
        );
    }
    let small = curve(1e-4);
    let mut bounded = heated(&small, 0.0, vec![route(0, "body")]);
    let gate = CancelGate::new_clock_free();
    let state = bounded.mechanics().state().to_vec();
    let thermal = bounded.reservoirs().to_vec();
    assert!(bounded.step(&[0.0], &gate).is_err());
    assert_eq!(bounded.mechanics().state(), state);
    assert_eq!(bounded.reservoirs(), thermal);
    assert_eq!(bounded.mechanics().samples(), 0);
    let mut valid = heated(&c, 0.0, vec![route(0, "body")]);
    let mut clean = heated(&c, 0.0, vec![route(0, "body")]);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(valid.step(&[0.0], &cancelled).is_err());
    valid.set_iteration_limit(0).unwrap();
    assert!(valid.step(&[0.0], &gate).is_err());
    assert_eq!(valid.mechanics().state(), clean.mechanics().state());
    assert_eq!(valid.reservoirs(), clean.reservoirs());
    valid.set_iteration_limit(20).unwrap();
    valid.step(&[0.0], &gate).unwrap();
    clean.step(&[0.0], &gate).unwrap();
    assert_eq!(valid.mechanics().state(), clean.mechanics().state());
    assert_eq!(valid.reservoirs(), clean.reservoirs());
}

fn network(c: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events = [9, 3]
        .map(|id| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: id,
                material_id: MaterialId(1),
                mass_kg: MASS,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [id as f64, 0.0, 0.0],
                    end_m: [id as f64 + 0.01, 0.0, 0.0],
                    volume_m3: 1e-5,
                    surface_area_m2: 0.001,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .to_vec();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(9),
            b: CellId(3),
            area_m2: 1e-4,
            conductance_w_k: 0.2,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), c as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 64,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: DT / 4.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
fn original_enthalpy_network_conducts_relaxation_heat_on_the_mechanical_clock() {
    let c = curve(1000.0);
    let gate = CancelGate::new_clock_free();
    let mut model = heated(&c, 0.0, vec![route(0, "body"), route(1, "body")]);
    let initial = model.mechanics().stored_energy_j();
    let mut thermal = network(&c);
    let bindings = [
        ReservoirCell {
            reservoir: "support".into(),
            cell: CellId(3),
        },
        ReservoirCell {
            reservoir: "body".into(),
            cell: CellId(9),
        },
    ];
    for sample in 1..=50 {
        let frame = model
            .step_with_thermal_network(
                &mut thermal,
                &mut (),
                &bindings,
                StepConfig {
                    max_dt_s: DT / 4.0,
                    energy_tolerance_j: 1e-10,
                },
                &[0.0],
                &gate,
                |_, _| Ok(None),
            )
            .unwrap();
        near(frame.external_heat_j, 0.0, 1e-14);
        near(
            frame.impact.thermal.impact.combined_energy_residual_j,
            0.0,
            1e-10,
        );
        assert_eq!(thermal.time_s(), f64::from(sample) * DT);
        assert_eq!(model.mechanics().samples(), sample as u64);
        for (index, id) in [9, 3].into_iter().enumerate() {
            assert_eq!(
                model.reservoirs()[index].state().specific_enthalpy_j_kg(),
                thermal.cells()[&CellId(id)].specific_h_j_kg
            );
        }
    }
    assert!(thermal.temperature_k(CellId(3)).unwrap() > 300.0);
    assert!(thermal.temperature_k(CellId(9)).unwrap() > thermal.temperature_k(CellId(3)).unwrap());
    assert_eq!(model.reservoirs()[1].received_heat_j(), 0.0);
    assert!(model.reservoirs()[1].boundary_heat_j() > 0.0);
    near(
        model.mechanics().stored_energy_j() + thermal.energy_j(),
        initial,
        1e-9,
    );
}
