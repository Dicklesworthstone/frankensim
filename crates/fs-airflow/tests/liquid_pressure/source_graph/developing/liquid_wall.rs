//! Original liquid wall Prandtl correction on accepted developing graph flow.
use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::developing::DevelopingSectionEvaluation;

fn wall_property_card(ceiling: f64) -> MaterialCard {
    let base = card_with_transport_ceiling(true, 1.0, ceiling);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(250.0, 0.25), (450.0, 1.25)],
                dims: Dims([1, 1, -3, -1, 0, 0]),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, vec![]).unwrap()
}

fn corrected_sections() -> [GraphDevelopingSection; 2] {
    let mut bindings = sections();
    for binding in &mut bindings {
        binding.duct = DevelopingDuct::TurbulentCircularLiquidWall {
            diameter: Length::new(diameter()),
        };
        binding.property_temperature = ConvectionPropertyTemperature::Reservoir;
    }
    bindings
}

fn flow_rate() -> f64 {
    20_000.0 * 0.002 * PIPE_AREA / (1000.0 * diameter())
}

fn wall_prandtl(temperature: f64) -> f64 {
    2000.0 * (0.005 - 1e-5 * temperature) / (0.005 * temperature - 1.0)
}

// Independent corrected-friction and rationalized section equation. All
// transport curves are the original synthetic source, not production queries.
fn expected_h(rate: f64, bulk: f64, wall: f64, start: f64) -> f64 {
    let mu = 0.005 - 1e-5 * bulk;
    let conductivity = 0.005 * bulk - 1.0;
    let reynolds = 1000.0 * rate.abs() * diameter() / (PIPE_AREA * mu);
    let prandtl = 2000.0 * mu / conductivity;
    let friction = (1.8 * reynolds.log10() - 1.5).powi(-2);
    let developed = friction / 8.0 * (reynolds - 1000.0) * prandtl
        / (1.0 + 12.7 * (friction / 8.0).sqrt() * (prandtl.powf(2.0 / 3.0) - 1.0));
    let a = (start / diameter()).cbrt();
    let b = ((start + 0.5) / diameter()).cbrt();
    developed * (1.0 + 1.0 / (a * a + a * b + b * b))
        * (prandtl / wall_prandtl(wall)).powf(0.11) * conductivity / diameter()
}

#[test]
fn g1_wall_prandtl_preserves_bulk_source_signed_entrance_and_actual_section_evidence() {
    context(|cx| {
        let card = wall_property_card(450.0);
        let (law, caloric) = source_pair(&card);
        let initial = section_fixture(&caloric, 400.0);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(1e6);
        let owner = LiquidGraphSources::try_new(
            cx, &law, &caloric, &initial, SourceGraphBinding {
                graph: &graph, inventory: &inventory, laws: &laws,
                acceleration_m_s2: [0.0; 3],
            }, caps(),
        ).unwrap();
        let bindings = corrected_sections();
        for sign in [-1.0, 1.0] {
            let q = sign * flow_rate();
            let resolved = owner.resolve(cx, &initial, &[q, q, -q]).unwrap();
            let hydraulics = resolved.graph.solve(&[
                FixedPressure { node: 0, pressure: Pressure::new(A * 0.002 * q) },
                FixedPressure { node: 1, pressure: Pressure::new(0.0) },
            ], hydraulic_controls(), cx).unwrap();
            let heat = owner.evaluate_quasi_steady_developing_sections(
                cx, &initial, &resolved.sources, &hydraulics, 0.01, &bindings, 2,
            ).unwrap();
            let mut uncorrected = bindings;
            for binding in &mut uncorrected {
                binding.duct = DevelopingDuct::TurbulentCircular {
                    diameter: Length::new(diameter()),
                };
            }
            let base = owner.evaluate_quasi_steady_developing_sections(
                cx, &initial, &resolved.sources, &hydraulics, 0.01, &uncorrected, 2,
            ).unwrap();
            for ((sample, baseline), binding) in heat.convection.samples().iter()
                .zip(base.convection.samples()).zip(bindings)
            {
                let start = if sign > 0.0 { binding.forward_start } else { binding.reverse_start };
                let expected = expected_h(q, 300.0, 400.0, start.value());
                near(sample.coefficient.value.value(), expected);
                near(sample.section.nusselt() * 0.5 / diameter(), expected);
                near(sample.conductance_w_k, area() * expected);
                near(sample.heat_into_fluid_w, 100.0 * area() * expected);
                assert_eq!(sample.thermal.transport(), baseline.thermal.transport());
                assert_eq!(sample.section.downstream().groups(), baseline.section.downstream().groups());
                near(sample.property_temperature_k, 300.0);
                near(sample.section.downstream().groups()["Re"], 20_000.0);
                near(sample.section.downstream().groups()["Pr"], 8.0);
                assert!(baseline.surface.is_none());
                let wall = sample.surface.as_ref().unwrap();
                near(wall.thermal.temperature_k(), 400.0);
                near(wall.thermal.transport().prandtl(), 2.0);
                assert_eq!(wall.material.material().card_identity(), card.content_hash());
                near(wall.material.material().query().axes()["pressure"], 100000.0);
                near(sample.coefficient.sensitivity.d_qoi["Pr_wall"], -0.11 * expected / 2.0);
                assert!(sample.coefficient.model.discrepancy_rel.is_infinite());
                let DevelopingSectionEvaluation::TurbulentCircular(section) = &sample.section else {
                    panic!("original turbulent section was not retained");
                };
                assert!(section.property_correction().is_some());
                // Wall conductivity changes too: a viscosity-only ratio gives
                // 2^0.11, whereas this original wall Prandtl ratio is 4.
                near(sample.conductance_w_k / baseline.conductance_w_k, 4.0_f64.powf(0.11));
            }
            let mut links = BTreeMap::new();
            let mut boundaries = BoundaryConductances::default();
            heat.convection.apply_to(&mut links, &mut boundaries, || false).unwrap();
            near(links[&9], heat.convection.samples()[0].conductance_w_k);
            near(boundaries.reservoir_w_k[&CellId(3)], heat.convection.samples()[1].conductance_w_k);
            let mut film = bindings;
            film[0].property_temperature = ConvectionPropertyTemperature::FilmMean;
            assert!(owner.evaluate_quasi_steady_developing_sections(
                cx, &initial, &resolved.sources, &hydraulics, 0.01, &film, 2,
            ).is_err());
        }

        // K=1 needs only admitted bulk transport. Enabling the correction
        // must query the real wall, never extend a source's temperature support.
        let limited_card = wall_property_card(350.0);
        let (limited, limited_caloric) = source_pair(&limited_card);
        let limited_initial = section_fixture(&limited_caloric, 400.0);
        let limited_inventory = bind_inventory(cx, &graph, &limited_initial);
        let limited_owner = LiquidGraphSources::try_new(
            cx, &limited, &limited_caloric, &limited_initial, SourceGraphBinding {
                graph: &graph, inventory: &limited_inventory, laws: &laws,
                acceleration_m_s2: [0.0; 3],
            }, caps(),
        ).unwrap();
        let q = flow_rate();
        let source = limited_owner.resolve(cx, &limited_initial, &[q, q, -q]).unwrap();
        let hydraulics = source.graph.solve(&[
            FixedPressure { node: 0, pressure: Pressure::new(A * 0.002 * q) },
            FixedPressure { node: 1, pressure: Pressure::new(0.0) },
        ], hydraulic_controls(), cx).unwrap();
        assert!(limited_owner.evaluate_quasi_steady_developing_sections(
            cx, &limited_initial, &source.sources, &hydraulics, 0.01, &bindings, 2,
        ).is_err());
        let mut baseline = bindings;
        for binding in &mut baseline {
            binding.duct = DevelopingDuct::TurbulentCircular { diameter: Length::new(diameter()) };
        }
        assert!(limited_owner.evaluate_quasi_steady_developing_sections(
            cx, &limited_initial, &source.sources, &hydraulics, 0.01, &baseline, 2,
        ).is_ok());
    });
}

#[test]
fn g3_corrected_sections_close_native_momentum_advection_and_one_wall_heat_budget() {
    use fs_airflow::graph::thermal::transport::transient::inertial::{
        InertialGraphHeat, InertialThermalControls, InertialThermalProblem,
        ResolvedInertialGraph, advance_inertial_graph_distributed,
    };
    use fs_solver::NewtonKrylovConfig;
    context(|cx| {
        let card = wall_property_card(450.0);
        let (law, caloric) = source_pair(&card);
        let initial = section_fixture(&caloric, 400.0);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &initial);
        let laws = laws(1e6);
        let owner = LiquidGraphSources::try_new(
            cx, &law, &caloric, &initial, SourceGraphBinding {
                graph: &graph, inventory: &inventory, laws: &laws,
                acceleration_m_s2: [0.0; 3],
            }, caps(),
        ).unwrap();
        let bindings = corrected_sections();
        let q = flow_rate();
        let rates = [q, 0.0, q];
        let first = owner.resolve(cx, &initial, &rates).unwrap();
        let branches: Vec<_> = inventory.inertial_branches().iter().map(|b| b.unwrap()).collect();
        let pressure = vec![FixedPressure { node: 0, pressure: Pressure::new(0.0) }];
        let momentum = InertialGraph::try_new(
            cx, &first.graph, Density::new(1000.0), &branches, &pressure, &rates, 0.0, 32,
        ).unwrap();
        let problem = InertialThermalProblem {
            momentum: &momentum, node_cells: inventory.node_cells(), maximum_connections: 32,
        };
        let controls = InertialThermalControls {
            solver: NewtonKrylovConfig {
                absolute_tolerance: 1e-9, relative_tolerance: 0.0, ..Default::default()
            },
            maximum_cells: 16, maximum_newton_iterations: 32, maximum_endpoint_evaluations: 2048,
            energy_scale_j: 1.0, flow_scale_m3_s: q, finite_difference_step: 1e-4,
            cell_tolerance_j: 1e-6, total_thermal_tolerance_j: 1e-6,
            impulse_tolerance_pa_s: 1e-6, mechanical_tolerance_j: 1e-8,
            average_flow_tolerance_m3_s: 1e-12,
            momentum: InertialGraphControls {
                maximum_drag_step: 2.0, maximum_neglected_power_w: 100.0,
                maximum_neglected_energy_j: 100.0,
            },
        };
        let dt = 0.01;
        let (next_momentum, next, report) = advance_inertial_graph_distributed(
            cx, &problem, &inventory, &initial, dt, controls, &[0.0; 7], &rates,
            |cx, candidate, flows| {
                let resolved = owner.resolve(cx, candidate, flows).map_err(|e| e.to_string())?;
                Ok::<_, String>(ResolvedInertialGraph {
                    graph: resolved.graph, pressure_boundaries: pressure.clone(),
                    driving_pressure_pa: vec![A * 0.002 * q, 0.0, A * 0.002 * q],
                    payload: resolved.sources,
                })
            },
            |cx, candidate, sources, hydraulics| {
                let heat = owner.evaluate_developing_sections(
                    cx, candidate, sources, hydraulics, &bindings, 2,
                ).map_err(|e| e.to_string())?;
                let mut links = BTreeMap::new();
                let mut boundaries = BoundaryConductances::default();
                heat.convection.apply_to(&mut links, &mut boundaries, || cx.checkpoint().is_err())
                    .map_err(|e| e.to_string())?;
                Ok::<_, String>(InertialGraphHeat {
                    mass_boundaries: vec![],
                    heat: BoundaryStepCoefficients {
                        bulk: StepCoefficients {
                            sources_w: BTreeMap::new(), conductances_w_k: links, payload: heat,
                        }, boundaries,
                    },
                })
            },
        ).unwrap();
        let samples = report.heat.convection.samples();
        let rate = report.hydraulics.average_volume_flows_m3_s[0];
        for (sample, start) in samples.iter().zip([0.0, 0.5]) {
            let bulk = next.temperature_k(sample.fluid_cell).unwrap();
            let expected = expected_h(rate, bulk, sample.wall_temperature_k, start);
            near(sample.coefficient.value.value(), expected);
            near(sample.property_temperature_k, bulk);
            near(sample.surface.as_ref().unwrap().thermal.temperature_k(), sample.wall_temperature_k);
        }
        let external = dt * samples[1].heat_into_fluid_w;
        let stored: f64 = next.cells().iter().map(|(id, cell)| {
            cell.mass_kg * (cell.specific_h_j_kg - initial.cells()[id].specific_h_j_kg)
        }).sum();
        assert!((stored - external).abs() < 2e-6);
        assert!((report.thermal.external_heat_j - external).abs() < 2e-6);
        assert!(next.temperature_k(CellId(2)).unwrap() > 300.0);
        assert!(next.temperature_k(CellId(3)).unwrap() > 300.0);
        assert!(next.temperature_k(CellId(6)).unwrap() < 400.0);
        near(next.time_s(), dt);
        near(next_momentum.accepted_time_s(), dt);
        assert!(report.maximum_cell_residual_j <= controls.cell_tolerance_j);
        assert_eq!(initial.time_s(), 0.0);
        assert_eq!(initial.links()[&9].conductance_w_k, 0.0);
    });
}
