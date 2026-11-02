//! Source-backed aiding flow on the native finite wall/fluid heat link.
use super::*;
use fs_airflow::graph::thermal::medium::convection::AidingVerticalPlateFrame;
use fs_airflow::graph::thermal::transport::transient::material::convection::GraphConvectionError;

fn plate() -> GraphConvectiveLink {
    GraphConvectiveLink {
        link: 9,
        fluid_cell: CellId(2),
        characteristic_length: Length::new(0.025),
        property_temperature: ConvectionPropertyTemperature::FilmMean,
        law: GraphConvectionLaw::AidingVerticalPlate {
            branch: 0,
            plate_up: [0.0, 0.0, 1.0],
            maximum_relative_density_departure: 0.1,
        },
    }
}

// rho(300 K)=1000; a positive slope gives negative thermal expansion. Both
// signs are authored original material claims, not gas-law substitutions.
fn density_card(slope: f64) -> MaterialCard {
    let base = card(false, 1.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "density" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (250.0, 1000.0 - 50.0 * slope),
                    (450.0, 1000.0 + 150.0 * slope),
                ],
                dims: Dims([-3, 1, 0, 0, 0, 0]),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

#[test]
fn authored_frame_preserves_rotations_signed_flow_and_pure_limits() {
    for up in [[0.0, 0.0, 1.0], [2.0 / 7.0, 3.0 / 7.0, 6.0 / 7.0]] {
        for gravity in [0.0, 9.81] {
            let frame =
                AidingVerticalPlateFrame::try_new(up.map(|x| 7.0 * x), up.map(|x| -gravity * x))
                    .unwrap();
            for speed in [-0.02, 0.0, 0.02] {
                let settings = frame.settings(up.map(|x| speed * x), 0.1).unwrap();
                near(settings.speed_up_m_s, speed);
                near(settings.gravity_m_s2, gravity);
            }
        }
    }
    let frame = AidingVerticalPlateFrame::try_new([0.0, 0.0, 1.0], [0.0, 0.0, -9.81]).unwrap();
    // A vertical projection cannot admit a transverse stream or tilted plate.
    assert!(frame.settings([0.01, 0.0, 0.02], 0.1).is_err());
    assert!(AidingVerticalPlateFrame::try_new([0.1, 0.0, 1.0], [0.0, 0.0, -9.81],).is_err());
    assert!(AidingVerticalPlateFrame::try_new([0.0, 0.0, -1.0], [0.0, 0.0, -9.81],).is_err());
}

#[test]
fn graph_mixed_heat_retains_signed_source_film_and_debits_each_owner_once() {
    context(|cx| {
        for slope in [-0.5_f64, 0.5] {
            let card = density_card(slope);
            let (law, caloric) = source_pair(&card);
            let law = law
                .with_density_expansion(Temperature::new(250.0), Temperature::new(450.0))
                .unwrap();
            for (wall, gravity, q) in [
                (310.0, 9.81, -slope.signum() * 1e-5),
                (290.0, 9.81, slope.signum() * 1e-5),
                (310.0, 9.81, 0.0),
                (310.0, 0.0, 1e-5),
                (300.0, 9.81, 0.0),
                (300.0, 9.81, 1e-5),
            ] {
                let initial = thermal(&caloric, [300.0; 6], Some(wall), false);
                let graph = graph();
                let inventory = bind_inventory(cx, &graph, &initial);
                let laws = laws(2300.0);
                let owner = LiquidGraphSources::try_new(
                    cx,
                    &law,
                    &caloric,
                    &initial,
                    SourceGraphBinding {
                        graph: &graph,
                        inventory: &inventory,
                        laws: &laws,
                        acceleration_m_s2: [0.0, 0.0, -gravity],
                    },
                    caps(),
                )
                .unwrap();
                let resolved = owner.resolve(cx, &initial, &[q, 0.0, q]).unwrap();
                let report = native_report(
                    cx,
                    &resolved.graph,
                    &inventory,
                    [q, 0.0, q],
                    [1e5 * q, 0.0, 1e5 * q],
                    [q, 0.0, q],
                    0.1,
                );
                let heat = owner
                    .evaluate_links(cx, &initial, &resolved.sources, &report, &[plate()], 16)
                    .unwrap();
                let sample = &heat.convection.samples()[0];
                let film = 0.5 * (wall + 300.0);
                let rho = 1000.0 + slope * (film - 300.0);
                let beta = -slope / rho;
                near(sample.property_temperature_k, film);
                near(sample.material.density_kg_m3(), rho);
                near(sample.material.dynamic_viscosity_pa_s(), 0.002);
                near(
                    sample.material.material().query().axes()["pressure"],
                    100000.0,
                );
                near(sample.thermal.volumetric_expansion_per_k().unwrap(), beta);
                if wall == 300.0 && q == 0.0 {
                    assert!(sample.evaluation.is_none());
                    assert_eq!(sample.conductance_w_k, 0.0);
                } else {
                    let evaluation = sample.evaluation.as_ref().unwrap();
                    assert_eq!(
                        evaluation.nusselt.card().id,
                        CorrelationId::ManningQureshiAidingVerticalPlate
                    );
                    assert_eq!(
                        evaluation.nusselt.card().source.identifier,
                        "doi:10.1115/1.2910940"
                    );
                    // Independent dimensional closure using the original
                    // film density, viscosity, conductivity and Cp, while
                    // Q/A is the native reference-volume speed.
                    near(
                        evaluation.nusselt.groups()["Re"],
                        rho * (q / PIPE_AREA).abs() * 0.025 / 0.002,
                    );
                    near(evaluation.nusselt.groups()["Pr"], 8.0);
                    near(
                        evaluation.nusselt.groups()["Ra"],
                        gravity * (beta * (wall - 300.0)).abs() * 0.025_f64.powi(3)
                            / ((0.002 / rho) * (0.5 / (rho * 2000.0))),
                    );
                    near(
                        sample.conductance_w_k,
                        0.002 * evaluation.nusselt.evidence().value * 0.5 / 0.025,
                    );
                }
                near(
                    sample.heat_into_fluid_w,
                    sample.conductance_w_k * (wall - 300.0),
                );
                let mut coefficients = StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                };
                heat.convection
                    .apply_to(&mut coefficients.conductances_w_k, || false)
                    .unwrap();
                let problem = BackwardEulerStep::new(&initial, 0.1, || Ok(())).unwrap();
                let trial = problem
                    .evaluate(
                        &[0.0; 7],
                        |_| Ok::<_, DepositionError>(coefficients),
                        || Ok(()),
                    )
                    .unwrap();
                near(trial.residual_j()[2], -0.1 * sample.heat_into_fluid_w);
                near(trial.residual_j()[6], 0.1 * sample.heat_into_fluid_w);
                near(trial.residual_j().iter().sum(), 0.0);
            }
        }
    });
}

#[test]
fn mixed_graph_stop_is_continuous_but_reversal_cannot_hide_opposing_flow() {
    context(|cx| {
        let card = buoyant_card();
        let (law, caloric) = source_pair(&card);
        let law = law
            .with_density_expansion(Temperature::new(250.0), Temperature::new(450.0))
            .unwrap();
        let initial = thermal(&caloric, [300.0; 6], Some(310.0), false);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0, 0.0, -9.81],
            },
            caps(),
        )
        .unwrap();
        let mut stopped_conductance = Vec::new();
        for q in [0.0, 1e-15] {
            let resolved = owner.resolve(cx, &initial, &[q, 0.0, q]).unwrap();
            let report = native_report(
                cx,
                &resolved.graph,
                &inventory,
                [q, 0.0, q],
                [1e5 * q, 0.0, 1e5 * q],
                [q, 0.0, q],
                0.1,
            );
            let heat = owner
                .evaluate_links(cx, &initial, &resolved.sources, &report, &[plate()], 16)
                .unwrap();
            stopped_conductance.push(heat.convection.samples()[0].conductance_w_k);
        }
        assert!(stopped_conductance[0] > 0.0);
        near(stopped_conductance[0], stopped_conductance[1]);
        for (old, new) in [(1e-5, 0.0), (0.0, 2e-5), (-1e-5, 3e-5)] {
            let average = 0.5 * (old + new);
            let resolved = owner
                .resolve(cx, &initial, &[average, 0.0, average])
                .unwrap();
            let drive = 1e6 * (new - old) / 0.1 + 1e5 * average;
            let report = native_report(
                cx,
                &resolved.graph,
                &inventory,
                [old, 0.0, old],
                [drive, 0.0, drive],
                [average, 0.0, average],
                0.1,
            );
            owner
                .admit_actual(cx, &initial, &resolved.sources, &report)
                .unwrap();
            let heat =
                owner.evaluate_links(cx, &initial, &resolved.sources, &report, &[plate()], 16);
            if old < 0.0 {
                // Average upward speed is aiding, but the independently
                // admitted old downward endpoint opposes the hot wall.
                assert!(matches!(
                    heat,
                    Err(GraphConvectionError::Envelope { link: 9, .. })
                ));
            } else {
                assert!(heat.unwrap().convection.samples()[0].heat_into_fluid_w > 0.0);
            }
        }
    });
}

#[test]
fn local_cell_crossflow_cannot_borrow_a_vertical_whole_branch_axis() {
    context(|cx| {
        let card = buoyant_card();
        let (law, caloric) = source_pair(&card);
        let law = law
            .with_density_expansion(Temperature::new(250.0), Temperature::new(450.0))
            .unwrap();
        let bent = thermal(&caloric, [300.0; 6], Some(310.0), true);
        let mut schedule: Vec<_> = bent
            .cells()
            .values()
            .cloned()
            .map(|mut cell| {
                if cell.id == CellId(2) || cell.id == CellId(3) {
                    let dx: [f64; 3] = core::array::from_fn(|axis| {
                        cell.geometry.end_m[axis] - cell.geometry.start_m[axis]
                    });
                    let length = dx[0].hypot(dx[1]).hypot(dx[2]);
                    cell.geometry.volume_m3 = PIPE_AREA * length;
                    cell.mass_kg = 1000.0 * cell.geometry.volume_m3;
                }
                TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell),
                }
            })
            .collect();
        schedule.extend(bent.links().values().map(|link| TimedEvent {
            time_s: 0.0,
            event: Event::Contact(*link),
        }));
        let initial = DepositionSimulation::new(
            schedule,
            &[(MaterialId(1), caloric.curve() as &dyn ThermalLaw)],
            DepositionBudget {
                max_cells: 16,
                max_links: 16,
                max_events: 16,
                max_substeps: 16,
            },
        )
        .unwrap()
        .advance_to(
            0.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap()
        .0;
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(2300.0);
        let owner = LiquidGraphSources::try_new(
            cx,
            &law,
            &caloric,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0, 0.0, -9.81],
            },
            caps(),
        )
        .unwrap();
        for q in [0.0, 1e-5] {
            let resolved = owner.resolve(cx, &initial, &[q, 0.0, q]).unwrap();
            let report = native_report(
                cx,
                &resolved.graph,
                &inventory,
                [q, 0.0, q],
                [1e5 * q, 0.0, 1e5 * q],
                [q, 0.0, q],
                0.1,
            );
            let heat =
                owner.evaluate_links(cx, &initial, &resolved.sources, &report, &[plate()], 16);
            if q == 0.0 {
                assert!(heat.unwrap().convection.samples()[0].heat_into_fluid_w > 0.0);
            } else {
                assert!(matches!(
                    heat,
                    Err(GraphConvectionError::Geometry { link: 9, .. })
                ));
            }
        }
    });
}
