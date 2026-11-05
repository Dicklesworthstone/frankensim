//! Native plate strips retain the actual mixture, film state and caloric owner.
use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::{
    PlateFlowDirection, reservoir::GraphConvectiveReservoir,
};
use fs_conduction::transient::deposition::{Reservoir, backward_euler::BoundaryConductances};
use fs_convection::CorrelationId;

#[test]
fn g1_g3_gas_plate_strips_use_original_film_eos_and_cumulative_wall_heat() {
    context(|cx| {
        let law = gas_law();
        let original = thermal(&law, [300.0; 6], Some(360.0));
        let area = original.links()[&11].area_m2;
        let mut schedule: Vec<_> = original
            .cells()
            .values()
            .cloned()
            .map(|mut cell| {
                if cell.id == CellId(3) {
                    cell.boundary.reservoir = Some(Reservoir {
                        area_m2: area,
                        conductance_w_k: 0.0,
                        temperature_k: 360.0,
                    });
                }
                TimedEvent {
                    time_s: 0.0,
                    event: Event::Birth(cell),
                }
            })
            .collect();
        schedule.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(original.links()[&11].clone()),
        });
        let initial = DepositionSimulation::new(
            schedule,
            &[(MaterialId(1), &law as &dyn ThermalLaw)],
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
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
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
            0.1,
        )
        .unwrap();
        let q = 1e-4;
        let flows = [q, 0.0, q];
        let resolved = owner.resolve(cx, &initial, &flows).unwrap();
        let drag = VISCOUS_GEOMETRY * owner.reference().state().dynamic_viscosity;
        let hydraulics = native_report(
            cx,
            &resolved.graph,
            &inventory,
            flows,
            [drag * q, 0.0, drag * q],
            flows,
            0.1,
        );
        let plate = |start| GraphConvectionLaw::ForcedPlateSection {
            branch: 0,
            correlation: CorrelationId::FlatPlateLaminarAverage,
            start: Length::new(start),
            direction: PlateFlowDirection::AlongBranch,
        };
        let links = owner
            .evaluate_links(
                cx,
                &initial,
                &resolved.sources,
                &hydraulics,
                &[GraphConvectiveLink {
                    link: 11,
                    fluid_cell: CellId(2),
                    characteristic_length: Length::new(0.5),
                    property_temperature: ConvectionPropertyTemperature::FilmMean,
                    law: plate(0.0),
                }],
                1,
            )
            .unwrap();
        let reservoirs = owner
            .evaluate_reservoirs(
                cx,
                &initial,
                &resolved.sources,
                &hydraulics,
                &[GraphConvectiveReservoir {
                    fluid_cell: CellId(3),
                    characteristic_length: Length::new(1.0),
                    property_temperature: ConvectionPropertyTemperature::FilmMean,
                    law: plate(0.5),
                }],
                1,
            )
            .unwrap();
        let first = &links.convection.samples()[0];
        let second = &reservoirs.convection.samples()[0];
        let (_, cp_intercept, cp_slope_half) = coefficients();
        for (material, temperature) in [
            (&first.material, first.property_temperature_k),
            (&second.material, second.property_temperature_k),
        ] {
            near(temperature, 330.0);
            near(material.state().pressure, PRESSURE);
            near(
                material.state().density,
                PRESSURE * MOLAR_MASS / (R * temperature),
            );
            near(material.specific_enthalpy_j_kg(), h(temperature));
            near(
                material.transport_state().unwrap().specific_heat_j_kg_k(),
                cp_intercept + 2.0 * cp_slope_half * temperature,
            );
        }
        assert_eq!(owner.source().mole_composition().fractions(), &[0.75, 0.25]);
        let transport = first.thermal.transport();
        assert_eq!(transport, second.thermal.transport());
        // Independent laminar primitive C sqrt(x): the second strip must
        // subtract its upstream transfer instead of restarting a leading edge.
        let c = 0.664
            * transport.thermal_conductivity_w_m_k()
            * (transport.density_kg_m3() * (q / AREA) / transport.dynamic_viscosity_pa_s()).sqrt()
            * transport.prandtl().powf(1.0 / 3.0);
        near(first.conductance_w_k, area * c / 0.5_f64.sqrt());
        near(second.conductance_w_k, area * c / (1.0 + 0.5_f64.sqrt()));
        near(
            first.heat_into_fluid_w + second.heat_into_fluid_w,
            60.0 * 2.0 * area * c,
        );
        assert!(
            first
                .evaluation
                .as_ref()
                .unwrap()
                .section
                .as_ref()
                .unwrap()
                .upstream
                .is_none()
        );
        assert!(
            second
                .evaluation
                .as_ref()
                .unwrap()
                .section
                .as_ref()
                .unwrap()
                .upstream
                .is_some()
        );
        let mut finite = BTreeMap::new();
        let mut prescribed = BoundaryConductances::default();
        links.convection.apply_to(&mut finite, || false).unwrap();
        reservoirs
            .convection
            .apply_to(&mut prescribed, || false)
            .unwrap();
        near(finite[&11], first.conductance_w_k);
        near(prescribed.reservoir_w_k[&CellId(3)], second.conductance_w_k);
        assert!(prescribed.ambient_w_k.is_empty());
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
    });
}
