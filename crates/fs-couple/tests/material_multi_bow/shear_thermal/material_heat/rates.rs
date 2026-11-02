//! Manufactured rate/heat composition checks. These coefficients are analytical
//! fixtures, not identified horsehair, rosin or acoustic calibration data.
use super::*;
use fs_couple::modal_acoustic_time::{ModalAcousticMode, ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::{
    material::MaterialRelaxationHeatDestination, rates::MaterialRelaxationTemperature,
};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    FrictionModalSystem, ModalFrictionPort,
    multiple::{FrictionPortFrame, MultiFrictionModalSystem},
};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem, ModalAttachment};
use fs_math::c64::C64;
use fs_phs::RelaxationBranch;

fn budget() -> MaterialHeatBudget {
    MaterialHeatBudget {
        max_groups: 1,
        max_dimension: 256,
        max_dense_bytes: 32 * 1024 * 1024,
        max_multiply_terms: 10_000_000_000,
    }
}

fn caloric_rate_source() -> (
    MaterialCard,
    ResolvedHeatCapacityCurve,
    ResolvedMaterialStatePoint,
) {
    caloric_rate_source_with_knots(vec![(300.0, 0.01), (320.0, 0.03)])
}

fn caloric_rate_source_with_knots(knots: Vec<(f64, f64)>) -> (
    MaterialCard, ResolvedHeatCapacityCurve, ResolvedMaterialStatePoint,
) {
    let (original, _) = coating(false, false);
    let mut claims = original.claims().clone();
    // Explicitly use a nonstandard selector: native memory must follow the
    // physical producer's mapping, not infer a property named "tau".
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("retardation_clock", Time::DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots,
                dims: Time::DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, 320.0)
                .with("omega", 0.0, 1e6),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured analytical caloric Maxwell relaxation".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "analytic-rate-heat".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(300.0),
            upper: &point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-8,
    )
    .unwrap();
    let requirements = [
        ("density", Density::DIMS),
        ("retardation_clock", Time::DIMS),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    })
    .collect::<Vec<_>>();
    let source = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (card, caloric, source)
}

fn isolated_body(caloric: &EquilibriumEnthalpyPhaseCurve, volume: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: BODY,
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: DENSITY * volume,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.001, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 1e-6,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(0), caloric as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 4,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1e-5,
            energy_tolerance_j: 1e-14,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn isolated_memory(rate: u32, gate: &CancelGate) -> MultiFrictionModalSystem {
    let oscillator = ModalAcousticTimeModel::try_new(
        rate,
        vec![ModalAcousticMode {
            angular_frequency_rad_s: 20.0,
            damping_ratio: 0.0,
            pressure_per_modal_velocity: C64::ZERO,
        }],
        ModalAcousticTimeBudget::audible_reference(),
    )
    .unwrap();
    let network =
        CoupledModalSystem::new(vec![oscillator], vec![], options().network, gate).unwrap();
    FrictionModalSystem::new(
        network,
        ModalFrictionPort {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0],
            },
            right: None,
            law: law(),
        },
        options().friction,
        gate,
    )
    .unwrap()
    .with_material_relaxation(
        vec![RelaxationBranch {
            projection: vec![0.0; 2],
            stiffness: 100.0,
            relaxation_time_s: 0.01,
        }],
        Some(vec![0.1]),
        gate,
    )
    .unwrap()
    .with_additional_ports(vec![], joint_limits(), gate)
    .unwrap()
}

#[test]
fn material_rate_heat_matches_caloric_relaxation_and_first_order_time_refinement() {
    let gate = CancelGate::new_clock_free();
    let (card, caloric, source) = caloric_rate_source();
    let volume = 1e-9;
    let capacity = DENSITY * volume * CP;
    let initial_energy = 0.5 * 0.1_f64.powi(2);
    let maximum_rise = initial_energy / capacity;
    // Independent continuous ODE: C Tdot = 2(E0-C dT)/tau(T), with
    // tau(T)=0.01(1+0.1 dT). Integrate dT analytically, then invert once.
    let elapsed = |rise: f64| {
        0.005 * (-(1.0 + 0.1 * maximum_rise) * (1.0 - rise / maximum_rise).ln() - 0.1 * rise)
    };
    let (mut low, mut high) = (0.0, maximum_rise);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if elapsed(mid) < 0.032 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let exact = 300.0 + 0.5 * (low + high);
    let mut errors = Vec::new();
    for rate in [1000, 2000, 4000] {
        let mut system = isolated_memory(rate, &gate);
        let mut heat = isolated_body(caloric.curve(), volume);
        let mut routes = [MaterialRelaxationHeatDestination {
            cell: BODY,
            material: &card,
            caloric: caloric.curve(),
            temperature_axis: "T",
            density_property: "density",
            volume_m3: volume,
            source: source.clone(),
            channel: Some(
                system
                    .prepare_material_heat(&[0..1], budget(), &gate)
                    .unwrap(),
            ),
            nonlinear_channels: vec![],
        }];
        let binding = [MaterialRelaxationTemperature {
            cell: BODY,
            times: vec![(0, "retardation_clock".into())],
        }];
        let mut received = 0.0;
        for _ in 0..rate * 32 / 1000 {
            let old_memory = system.material_memory()[0];
            let old_tau = system.material_relaxation_branches()[0].relaxation_time_s;
            let (frame, change) = system
                .step_with_material_rate_thermal_feedback(
                    &mut heat,
                    &[],
                    &[],
                    &[],
                    &mut routes,
                    &binding,
                    controls(1),
                    &[0.0],
                    &[0.0],
                    &[drive(0.0, 0.0)],
                    &gate,
                    |s, _| other(s),
                )
                .unwrap();
            // Only the accepted rate drives THIS interval: no midpoint or
            // endpoint power may replace the exact physical loss integral.
            let expected_memory = old_memory * (-1.0 / (f64::from(rate) * old_tau)).exp();
            assert!((system.material_memory()[0] - expected_memory).abs() < 1e-13);
            assert_eq!(change.changes.len(), 1);
            let endpoint_energy = 0.5 * expected_memory.powi(2);
            let changed = change.changes[0];
            assert!((changed.stored_energy_j - endpoint_energy).abs() < 1e-13);
            assert!((changed.power_before_w - 2.0 * endpoint_energy / old_tau).abs() < 1e-12);
            assert!(
                (changed.power_after_w - 2.0 * endpoint_energy / changed.relaxation_time_s).abs()
                    < 1e-12
            );
            received += frame.routed_material_heat_j;
            assert!((received + system.total_energy_j().unwrap() - initial_energy).abs() < 2e-11);
            assert!((received - total_heat(&heat)).abs() < 2e-11);
            assert_eq!(
                routes[0].source.state_coordinate("T"),
                Some(heat.temperature_k(BODY).unwrap())
            );
            assert_eq!(
                system.material_relaxation_branches()[0].relaxation_time_s,
                routes[0]
                    .source
                    .property("retardation_clock")
                    .unwrap()
                    .value_si()
            );
            assert!(frame.total_energy_residual_j.abs() <= controls(1).total_energy_tolerance_j);
        }
        errors.push((heat.temperature_k(BODY).unwrap() - exact).abs());
    }
    assert!(
        errors[0] / errors[1] > 1.8 && errors[1] / errors[2] > 1.8,
        "{errors:?}"
    );
    assert!(errors[2] < 0.0004, "{errors:?}");
}

fn branch_snapshot(system: &MultiFrictionModalSystem) -> Vec<(f64, f64, Vec<f64>)> {
    system
        .material_relaxation_branches()
        .iter()
        .map(|b| (b.stiffness, b.relaxation_time_s, b.projection.clone()))
        .collect()
}

#[test]
fn moving_material_rate_source_refusal_restores_warm_geometry_heat_and_history() {
    use fs_couple::render::schedule::force::coupled::contact::friction::multiple::surfaces::geometry::{
        FrictionGeometry, FrictionEndpointGeometry, FrictionGeometryConfig,
    };
    let gate = CancelGate::new_clock_free();
    let volume = 1e-9;
    let dt = 0.001;
    let fraction = -(-2.0_f64 * dt / 0.01).exp_m1();
    let first_rise = 0.5 * 0.1_f64.powi(2) * fraction / (DENSITY * volume * CP);
    // The first interval has a valid constant rate. The second's physical heat
    // reaches a positive source time whose reciprocal is unrepresentable. The
    // caloric solve itself remains valid; refusal belongs to endpoint mobility.
    let (card, caloric, source) = caloric_rate_source_with_knots(vec![
        (300.0, 0.01), (300.0 + 1.5 * first_rise, 0.01),
        (300.0 + 1.6 * first_rise, 1e-310), (320.0, 1e-310),
    ]);
    let mut system = isolated_memory(1000, &gate);
    let mut heat = isolated_body(caloric.curve(), volume);
    let mut routes = [MaterialRelaxationHeatDestination {
        cell: BODY, material: &card, caloric: caloric.curve(), temperature_axis: "T",
        density_property: "density", volume_m3: volume, source,
        channel: Some(system.prepare_material_heat(&[0..1], budget(), &gate).unwrap()),
        nonlinear_channels: vec![],
    }];
    let rates = [MaterialRelaxationTemperature { cell: BODY, times: vec![(0, "retardation_clock".into())] }];
    let geometry = |shape| FrictionGeometry {
        left: ModalAttachment { component: 0, shapes: vec![shape] }, right: None,
    };
    let endpoint = |a, b| FrictionEndpointGeometry { start: geometry(a), end: geometry(b),
        start_speed_m_s: 0.02, end_speed_m_s: 0.02 };
    let config = FrictionGeometryConfig { max_update_terms: 100_000 };
    system.step_with_material_rate_thermal_geometry(&mut heat, &[], &[], &[], &mut routes,
        &rates, controls(1), &[0.0], &[0.0], &[drive(0.02, 0.001)], &[geometry(1.5)],
        Some(&[endpoint(1.0, 2.0)]), config, &gate, |s, _| other(s)).unwrap();
    assert!((heat.temperature_k(BODY).unwrap() - 300.0 - first_rise).abs() < 1e-8);
    let states = system.components()[0].states().to_vec();
    let memory = system.material_memory().to_vec();
    let branches = branch_snapshot(&system);
    let source = routes[0].source.clone();
    let temperature = heat.temperature_k(BODY).unwrap();
    let energy = total_heat(&heat);
    let clock = heat.time_s();
    let failure = system.step_with_material_rate_thermal_geometry(&mut heat, &[], &[], &[],
        &mut routes, &rates, controls(1), &[energy], &[0.0], &[drive(0.02, 0.001)],
        &[geometry(2.5)], Some(&[endpoint(2.0, 3.0)]), config, &gate, |s, _| other(s))
        .unwrap_err();
    assert!(matches!(failure, LayerThermalError::Invalid(
        "source relaxation time and rate must be positive representable values")));
    assert_eq!(system.port(0).unwrap().left, geometry(1.5).left);
    assert_eq!(system.components()[0].states(), states);
    assert_eq!(system.material_memory(), memory);
    assert_eq!(branch_snapshot(&system), branches);
    assert_eq!(routes[0].source, source);
    assert_eq!(system.samples_rendered(), 1);
    assert_eq!(heat.time_s(), clock);
    assert_eq!(heat.temperature_k(BODY).unwrap(), temperature);
    assert_eq!(total_heat(&heat), energy);
    system.admit_material_rate_thermal_feedback(&heat, &routes, &rates, &gate).unwrap();
    // An explicit external cooling source can retain the admitted temperature
    // on retry. It removes exactly this independent uncoupled arm's heat; no
    // geometry work or artificial energy reset is used to clear the refusal.
    let cooling = 0.5 * memory[0].powi(2) * fraction / dt;
    system.step_with_material_rate_thermal_geometry(&mut heat, &[], &[], &[], &mut routes,
        &rates, controls(1), &[energy], &[0.0], &[drive(0.02, 0.001)], &[geometry(2.5)],
        Some(&[endpoint(2.0, 3.0)]), config, &gate, |_, _| Ok::<_, Infallible>(StepCoefficients {
            sources_w: BTreeMap::from([(BODY, -cooling)]), conductances_w_k: BTreeMap::new(), payload: (),
        })).unwrap();
    assert_eq!(system.samples_rendered(), 2);
    assert_eq!(system.port(0).unwrap().left, geometry(2.5).left);
    assert!((system.material_memory()[0] - memory[0] * (-dt/0.01).exp()).abs() < 1e-13);
}

fn warm_body(
    _: &DepositionSimulation<'_>,
    _: &[FrictionPortFrame],
) -> Result<StepCoefficients<()>, Infallible> {
    // Explicit outside thermal work on the ACTUAL source mass. It accelerates
    // the manufactured experiment; no tiny surrogate receiving mass is used.
    Ok(StepCoefficients {
        sources_w: BTreeMap::from([(BODY, 5000.0)]),
        conductances_w_k: BTreeMap::new(),
        payload: (),
    })
}

#[test]
fn material_rate_heat_refreshes_fiber_contact_tail_and_rolls_back_source_and_memory() {
    with_cx(|cx, gate| {
        let (film_card, film_caloric) = coating(true, false);
        let (body_card, body_caloric, ribbon) = fiber(gate, true);
        let (fixed_card, fixed_caloric, fixed_ribbon) = fiber(gate, false);
        let make = |ribbon| {
            base(cx, gate, 7800.0)
                .with_residual_flexibility(gate)
                .unwrap()
                .with_tensioned_fiber_ribbon(ribbon, 0.4, ribbon_limits(), gate)
                .unwrap()
                .with_bow_shear_layers(vec![(0, layer(&film_card))], gate)
                .unwrap()
                .with_bow_residual_flexibility(gate)
                .unwrap()
        };
        let (mut bank, mut fixed) = (make(ribbon), make(fixed_ribbon));
        let (route, rate) = bank
            .fiber_material_rate_heat_destination(
                0,
                0,
                BODY,
                &body_card,
                body_caloric.curve(),
                "T",
                budget(),
                gate,
            )
            .unwrap();
        assert_eq!(rate.times, vec![(0, "tau".into()), (1, "tau".into())]);
        let (reference, fixed_rate) = fixed
            .fiber_material_rate_heat_destination(
                0,
                0,
                BODY,
                &fixed_card,
                fixed_caloric.curve(),
                "T",
                budget(),
                gate,
            )
            .unwrap();
        let (mut heat, mut fixed_heat) = (
            finite_bodies(film_caloric.curve(), body_caloric.curve(), route.volume_m3),
            finite_bodies(
                film_caloric.curve(),
                fixed_caloric.curve(),
                reference.volume_m3,
            ),
        );
        let mut routes = [route];
        let mut references = [reference];
        let bound = binding(&film_card, film_caloric.curve());
        let receiver = bank.system().port(0).unwrap().right.clone();
        assert!(receiver.is_some());
        let first = bank.source().fiber_ribbon().unwrap().component_start();
        let mode = bank.system().components()[..first]
            .iter()
            .map(|m| m.modes().len())
            .sum::<usize>();
        let mut force = vec![0.0; bank.system().mode_count()];
        force[mode] = 10.0;
        let drives = [drive(0.001, 0.1)];
        // Missing and aliased source roles cannot be guessed, even when an
        // unselected property would happen to have an admissible numeric value.
        for times in [
            vec![(0, "tau".into())],
            vec![(0, "density".into()), (1, "tau".into())],
        ] {
            let invalid = [MaterialRelaxationTemperature { cell: BODY, times }];
            assert!(
                bank.step_with_material_rate_thermal_feedback(
                    &mut heat,
                    &bound,
                    &[],
                    &[],
                    &mut routes,
                    &invalid,
                    controls(2),
                    &[0.0; 2],
                    &force,
                    &drives,
                    gate,
                    warm_body
                )
                .is_err()
            );
        }
        assert_eq!(bank.system().samples_rendered(), 0);
        assert_eq!(total_heat(&heat), 0.0);
        let (mut pressure_difference, mut work, mut outside_loss, mut external_heat) =
            (0.0, 0.0, 0.0, 0.0);
        for step in 0..12 {
            if step == 4 {
                let qv = states(&bank);
                let memory = bank.system().material_memory().to_vec();
                let branches = branch_snapshot(bank.system());
                let film = bank.system().material_shear_state(0).unwrap().clone();
                let source = routes[0]
                    .source
                    .properties()
                    .iter()
                    .map(|p| p.value_si())
                    .collect::<Vec<_>>();
                let source_temperature = routes[0].source.state_coordinate("T");
                let cells = heat
                    .cells()
                    .iter()
                    .map(|(id, c)| (*id, c.mass_kg, c.specific_h_j_kg))
                    .collect::<Vec<_>>();
                let clock = heat.time_s();
                assert!(
                    bank.step_with_material_rate_thermal_feedback(
                        &mut heat,
                        &bound,
                        &[],
                        &[],
                        &mut routes,
                        std::slice::from_ref(&rate),
                        controls(2),
                        &[0.0; 2],
                        &force,
                        &drives,
                        gate,
                        |_, _| Err::<StepCoefficients<()>, _>("late material-rate thermal refusal")
                    )
                    .is_err()
                );
                assert_eq!(states(&bank), qv);
                assert_eq!(bank.system().material_memory(), memory);
                assert_eq!(branch_snapshot(bank.system()), branches);
                assert_eq!(bank.system().material_shear_state(0).unwrap(), &film);
                assert_eq!(
                    routes[0]
                        .source
                        .properties()
                        .iter()
                        .map(|p| p.value_si())
                        .collect::<Vec<_>>(),
                    source
                );
                assert_eq!(routes[0].source.state_coordinate("T"), source_temperature);
                assert_eq!(
                    heat.cells()
                        .iter()
                        .map(|(id, c)| (*id, c.mass_kg, c.specific_h_j_kg))
                        .collect::<Vec<_>>(),
                    cells
                );
                assert_eq!(heat.time_s(), clock);
                assert_eq!(bank.system().samples_rendered(), step);
            }
            let (frame, rates) = bank
                .step_with_material_rate_thermal_feedback(
                    &mut heat,
                    &bound,
                    &[],
                    &[],
                    &mut routes,
                    std::slice::from_ref(&rate),
                    controls(2),
                    &[0.0; 2],
                    &force,
                    &drives,
                    gate,
                    warm_body,
                )
                .unwrap();
            let (reference, fixed_changes) = fixed
                .step_with_material_rate_thermal_feedback(
                    &mut fixed_heat,
                    &bound,
                    &[],
                    &[],
                    &mut references,
                    std::slice::from_ref(&fixed_rate),
                    controls(2),
                    &[0.0; 2],
                    &force,
                    &drives,
                    gate,
                    warm_body,
                )
                .unwrap();
            assert_eq!(rates.changes.len(), 2);
            assert!(fixed_changes.changes.is_empty());
            assert_eq!(bank.system().port(0).unwrap().right, receiver);
            assert_eq!(
                routes[0].source.state_coordinate("T"),
                Some(heat.temperature_k(BODY).unwrap())
            );
            let tau = routes[0].source.property("tau").unwrap().value_si();
            assert!(
                bank.system()
                    .material_relaxation_branches()
                    .iter()
                    .all(|b| b.relaxation_time_s.to_bits() == tau.to_bits())
            );
            let material_heat = bank
                .system()
                .coupled_frame()
                .unwrap()
                .material_dissipation_j;
            assert!((frame.routed_material_heat_j - material_heat).abs() < 2e-11);
            pressure_difference += (frame.mechanics.observer_pressure_pa
                - reference.mechanics.observer_pressure_pa)
                .abs();
            work += frame.mechanics.external_work_j + frame.mechanics.drive_work_j;
            external_heat += 5000.0 * bank.system().sample_period_s();
            outside_loss += frame.mechanics.network_dissipation_j - frame.routed_material_heat_j
                + frame.mechanics.friction_residual_loss_j;
            // Film modulus is constant; its storage has no reversible T*S term.
            let balance =
                bank.system().total_energy_j().unwrap() + total_heat(&heat) + outside_loss
                    - work
                    - external_heat;
            assert!(balance.abs() < 2e-9, "{balance}");
            assert!(frame.total_energy_residual_j.abs() <= controls(2).total_energy_tolerance_j);
        }
        assert!(pressure_difference > 1e-12, "{pressure_difference}");
        assert!(bank.system().material_relaxation_branches()[0].relaxation_time_s < 0.0005);
    });
}
