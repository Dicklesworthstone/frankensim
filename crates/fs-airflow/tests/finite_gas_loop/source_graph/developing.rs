//! Developing-section convection retains the original gas and native flow limits.
//! Synthetic source data verify implementation; they are not calibration data.

use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::developing::{
    DevelopingWall, GraphDevelopingConvectionError, GraphDevelopingSection,
};
use fs_conduction::transient::deposition::{
    Reservoir, backward_euler::BoundaryConductances,
};
use std::collections::BTreeMap;

fn with_developing_reservoir(law: &GasEnthalpyLaw, diameter: f64) -> DepositionSimulation<'_> {
    let original = thermal(law, [300.0; 6], None);
    let schedule = original.cells().values().cloned().map(|mut cell| {
        if cell.id == CellId(2) {
            // The native first section is half a metre long. Its whole
            // circular wetted perimeter has one original reservoir owner.
            cell.boundary.reservoir = Some(Reservoir {
                area_m2: core::f64::consts::PI * diameter * 0.5,
                conductance_w_k: 0.0,
                temperature_k: 360.0,
            });
        }
        TimedEvent { time_s: 0.0, event: Event::Birth(cell) }
    }).collect();
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 16, max_links: 16, max_events: 16, max_substeps: 16 },
    ).unwrap().advance_to(
        0.0,
        StepConfig { max_dt_s: 1.0, energy_tolerance_j: 1e-8 },
        || Ok(()),
    ).unwrap().0
}

#[test]
fn g1_g4_developing_gas_retains_film_calorics_and_native_endpoint_mach() {
    context(|cx| {
        let law = gas_law();
        let diameter = (4.0 * AREA / core::f64::consts::PI).sqrt();
        let initial = with_developing_reservoir(&law, diameter);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        let binding = || SourceGraphBinding {
            graph: &graph, inventory: &inventory, laws: &laws,
            acceleration_m_s2: [0.0; 3],
        };
        let owner = GasGraphSources::try_new(
            cx, &law, &initial, binding(), source_limits(), 0.1,
        ).unwrap();
        let q = 1e-4;
        let rates = [q, 0.0, q];
        let resolved = owner.resolve(cx, &initial, &rates).unwrap();
        let drag = VISCOUS_GEOMETRY * owner.reference().state().dynamic_viscosity;
        let report = native_report(
            cx, &resolved.graph, &inventory, rates, [drag * q, 0.0, drag * q], rates, 0.1,
        );
        let sections = [GraphDevelopingSection {
            fluid_cell: CellId(2), wall: DevelopingWall::Reservoir,
            diameter: Length::new(diameter),
            forward_start: Length::new(1.0), reverse_start: Length::new(0.0),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
        }];
        let heat = owner.evaluate_developing_sections(
            cx, &initial, &resolved.sources, &report, &sections, 1,
        ).unwrap();
        let sample = &heat.convection.samples()[0];
        near(report.average_volume_flows_m3_s[0], q);
        near(sample.wall_temperature_k, 360.0);
        near(sample.fluid_temperature_k, 300.0);
        near(sample.property_temperature_k, 330.0);
        near(sample.material.state().pressure, PRESSURE);
        near(sample.material.state().density, PRESSURE * MOLAR_MASS / (R * 330.0));
        near(sample.material.specific_enthalpy_j_kg(), h(330.0));
        assert_eq!(owner.source().mole_composition().fractions(), &[0.75, 0.25]);
        let transport = sample.material.transport_state().unwrap();
        let (_, a, b) = coefficients();
        near(transport.specific_heat_j_kg_k(), a + 2.0 * b * 330.0);
        assert_eq!(sample.thermal.transport(), transport);
        let reynolds = transport.density_kg_m3() * (q / AREA) * diameter
            / transport.dynamic_viscosity_pa_s();
        let prandtl = transport.specific_heat_j_kg_k() * transport.dynamic_viscosity_pa_s()
            / transport.thermal_conductivity_w_m_k();
        assert!((1.0..2300.0).contains(&reynolds));
        // Independent cumulative mean subtraction. The section begins one
        // metre from the heated entrance, not at a new entrance at every cell.
        let mean = |distance: f64| {
            let gz = reynolds * prandtl * diameter / distance;
            3.66 + 0.0668 * gz / (1.0 + 0.04 * gz.powf(2.0 / 3.0))
        };
        let section_nusselt = (1.5 * mean(1.5) - mean(1.0)) / 0.5;
        let area = core::f64::consts::PI * diameter * 0.5;
        let conductance = area * section_nusselt * transport.thermal_conductivity_w_m_k()
            / diameter;
        near(sample.section.nusselt(), section_nusselt);
        assert!(sample.section.upstream().is_some());
        assert!(sample.section.nusselt() < mean(0.5));
        near(sample.area_m2, area);
        near(sample.length_m, 0.5);
        near(sample.conductance_w_k, conductance);
        near(sample.heat_into_fluid_w, 60.0 * conductance);
        let mut links = BTreeMap::new();
        let mut boundaries = BoundaryConductances::default();
        heat.convection.apply_to(&mut links, &mut boundaries, || false).unwrap();
        assert!(links.is_empty());
        near(boundaries.reservoir_w_k[&CellId(2)], conductance);

        // A zero average cannot hide the native trajectory's Mach violation;
        // the original source/flow owner refuses before forced-card evaluation.
        let bounded = GasGraphSources::try_new(
            cx, &law, &initial, binding(), source_limits(), 1e-4,
        ).unwrap();
        let low_mach = bounded.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let inertia = inventory.reference_density().value() / AREA;
        let drive = -2.0 * inertia * q / 0.1;
        let crossing = native_report(
            cx, &low_mach.graph, &inventory, rates, [drive, 0.0, drive], [0.0; 3], 0.1,
        );
        assert!(crossing.average_volume_flows_m3_s[0].abs() < 1e-12);
        assert!(matches!(bounded.evaluate_developing_sections(
            cx, &initial, &low_mach.sources, &crossing, &sections, 1,
        ), Err(GraphDevelopingConvectionError::Source(SourceGraphError::GasMach { .. }))));
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
    });
}
