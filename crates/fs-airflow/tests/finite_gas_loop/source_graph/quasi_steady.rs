//! Actual solved gas flow, original film properties and endpoint Mach admission.
use super::*;
use fs_airflow::graph::GraphSolveConfig;
use fs_airflow::graph::thermal::medium::convection::ForcedConvectionSettings;
use fs_airflow::graph::thermal::transport::transient::material::convection::GraphConvectionError;
use fs_convection::CorrelationId;
use fs_qty::VolumetricFlowRate;

#[test]
fn solved_gas_endpoint_uses_original_film_and_cannot_hide_mach_behind_source_probe() {
    context(|cx| {
        let law = gas_law();
        let initial = thermal(&law, [300.0; 6], Some(330.0));
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        let maximum_mach = 1e-3;
        let owner = GasGraphSources::try_new(
            cx,
            &law,
            &initial,
            SourceGraphBinding {
                graph: &graph,
                inventory: &inventory,
                laws: &laws,
                acceleration_m_s2: [0.0; 3],
            },
            source_limits(),
            maximum_mach,
        )
        .unwrap();
        // Zero is only the material-law probe. Each native pressure solve must
        // independently qualify its actual rate before any wall heat is used.
        let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let bulk = &resolved.sources.sections()[0][0].momentum.medium;
        let bulk_transport = bulk.transport_state().unwrap();
        let sound = bulk.state().sound_speed;
        let resistance = VISCOUS_GEOMETRY * bulk_transport.dynamic_viscosity_pa_s();
        let links = [GraphConvectiveLink {
            link: 11,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(0.03),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: GraphConvectionLaw::Forced {
                branch: 0,
                correlation: CorrelationId::CircularDuctLaminarCwt,
                settings: ForcedConvectionSettings {
                    length_over_hydraulic_diameter: Some(100.0),
                    ..ForcedConvectionSettings::default()
                },
            },
        }];
        let duration = 0.2;
        for fraction in [0.5, 2.0] {
            let expected_flow = fraction * maximum_mach * sound * AREA;
            // Independent linear pressure balance: Delta p = mu*geometry*Q.
            // The original resolved graph retains its native gravity receipt.
            let hydraulics = resolved
                .graph
                .solve(
                    &[
                        FixedPressure {
                            node: 0,
                            pressure: Pressure::new(resistance * expected_flow),
                        },
                        FixedPressure {
                            node: 1,
                            pressure: Pressure::new(0.0),
                        },
                    ],
                    GraphSolveConfig {
                        max_sweeps: 256,
                        max_node_iterations: 96,
                        absolute_flow_tolerance: VolumetricFlowRate::new(1e-14),
                        relative_flow_tolerance: 1e-12,
                    },
                    cx,
                )
                .unwrap();
            let actual_flow = hydraulics.branches[0].flow.value();
            near(actual_flow, expected_flow);
            near(hydraulics.branches[1].flow.value(), expected_flow);
            near(hydraulics.branches[2].flow.value(), -expected_flow);
            let admission =
                owner.admit_quasi_steady(cx, &initial, &resolved.sources, &hydraulics, duration);
            let heat = owner.evaluate_quasi_steady_links(
                cx,
                &initial,
                &resolved.sources,
                &hydraulics,
                duration,
                &links,
                16,
            );
            if fraction > 1.0 {
                match admission {
                    Err(SourceGraphError::GasMach {
                        branch: 0,
                        section: 0,
                        value,
                        maximum,
                    }) => {
                        near(value, actual_flow / (AREA * sound));
                        near(maximum, maximum_mach);
                        assert!(value > maximum);
                    }
                    other => panic!("actual solved Mach was not independently refused: {other:?}"),
                }
                assert!(matches!(
                    heat,
                    Err(GraphConvectionError::Source(SourceGraphError::GasMach {
                        branch: 0,
                        section: 0,
                        ..
                    }))
                ));
                continue;
            }
            let admission = admission.unwrap();
            near(admission.actual_losses[0][0].volume_flow_m3_s, actual_flow);
            near(admission.speed_envelope_m3_s[0], actual_flow);
            near(
                admission.neglected_hydraulic_power_w,
                3.0 * resistance * actual_flow.powi(2),
            );
            near(
                admission.neglected_energy_j,
                duration * admission.neglected_power_w,
            );
            let heat = heat.unwrap();
            let sample = &heat.convection.samples()[0];
            near(sample.property_temperature_k, 315.0);
            let direct = law
                .model()
                .evaluate(&point(sample.property_temperature_k))
                .unwrap();
            assert_eq!(sample.material, direct);
            assert_eq!(owner.source().mole_composition().fractions(), &[0.75, 0.25]);
            near(sample.material.state().pressure, PRESSURE);
            let film = sample.material.transport_state().unwrap();
            near(
                film.specific_heat_j_kg_k(),
                R / MOLAR_MASS * (3.25 + 0.75 * SLOPE * sample.property_temperature_k),
            );
            let evaluation = sample.evaluation.as_ref().unwrap();
            near(
                evaluation.nusselt.groups()["Re"],
                film.density_kg_m3() * actual_flow / AREA * 0.03 / film.dynamic_viscosity_pa_s(),
            );
            assert!(evaluation.nusselt.groups()["Re"] > 0.0);
            let conductance = 0.002 * 3.66 * film.thermal_conductivity_w_m_k() / 0.03;
            near(sample.conductance_w_k, conductance);
            near(sample.heat_into_fluid_w, 30.0 * conductance);
        }
    });
}
