//! A prescribed wall uses original gas film properties and native flow limits.
use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::reservoir::{
    GraphConvectiveReservoir, GraphReservoirConvectionError,
};
use fs_conduction::transient::deposition::{
    Reservoir, backward_euler::BoundaryConductances,
};

fn with_reservoir(law: &GasEnthalpyLaw) -> DepositionSimulation<'_> {
    let original = thermal(law, [300.0; 6], None);
    let schedule = original.cells().values().cloned().map(|mut cell| {
        if cell.id == CellId(2) {
            cell.boundary.reservoir = Some(Reservoir {
                area_m2: 0.002, conductance_w_k: 0.0, temperature_k: 360.0,
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
fn g1_g4_gas_reservoir_retains_original_film_eos_calorics_and_endpoint_mach_limit() {
    context(|cx| {
        let law = gas_law();
        let initial = with_reservoir(&law);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        let binding = || SourceGraphBinding {
            graph: &graph, inventory: &inventory, laws: &laws,
            acceleration_m_s2: [0.0, 0.0, -9.81],
        };
        let owner = GasGraphSources::try_new(
            cx, &law, &initial, binding(), source_limits(), 0.1,
        ).unwrap();
        let resolved = owner.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let report = native_report(
            cx, &resolved.graph, &inventory, [0.0; 3], [0.0; 3], [0.0; 3], 0.1,
        );
        let reservoirs = [GraphConvectiveReservoir {
            fluid_cell: CellId(2), characteristic_length: Length::new(0.03),
            property_temperature: ConvectionPropertyTemperature::FilmMean,
            law: GraphConvectionLaw::Natural { geometry: NaturalConvectionGeometry::VerticalPlate },
        }];
        let heat = owner.evaluate_reservoirs(
            cx, &initial, &resolved.sources, &report, &reservoirs, 1,
        ).unwrap();
        let sample = &heat.convection.samples()[0];
        near(sample.wall_temperature_k, 360.0);
        near(sample.fluid_temperature_k, 300.0);
        near(sample.property_temperature_k, 330.0);
        near(sample.material.state().pressure, PRESSURE);
        near(sample.material.state().density, PRESSURE * MOLAR_MASS / (R * 330.0));
        near(sample.material.specific_enthalpy_j_kg(), h(330.0));
        let transport = sample.material.transport_state().unwrap();
        let (_, a, b) = coefficients();
        near(transport.specific_heat_j_kg_k(), a + 2.0 * b * 330.0);
        // Independent Churchill-Chu algebra: full 60 K wall/bulk difference,
        // with beta=1/T supplied by this original fixed-pressure ideal mixture.
        let rayleigh = 9.81 / 330.0 * 60.0 * 0.03_f64.powi(3)
            / (transport.kinematic_viscosity_m2_s() * transport.thermal_diffusivity_m2_s());
        let prandtl = transport.prandtl();
        let nusselt = (0.825 + 0.387 * rayleigh.powf(1.0 / 6.0)
            / (1.0 + (0.492 / prandtl).powf(9.0 / 16.0)).powf(8.0 / 27.0)).powi(2);
        let conductance = 0.002 * nusselt * transport.thermal_conductivity_w_m_k() / 0.03;
        near(sample.conductance_w_k, conductance);
        near(sample.heat_into_fluid_w, 60.0 * conductance);
        assert!(heat.admission.actual_losses.iter().flatten()
            .all(|loss| loss.dissipated_power_w.abs() < 1e-20));
        let mut native = BoundaryConductances::default();
        heat.convection.apply_to(&mut native, || false).unwrap();
        near(native.reservoir_w_k[&CellId(2)], conductance);
        assert!(native.ambient_w_k.is_empty());

        // Average flow is zero, but the native trajectory's endpoints exceed
        // the original Mach ceiling. Boundary evaluation must retain that gate.
        let bounded = GasGraphSources::try_new(
            cx, &law, &initial, binding(), source_limits(), 1e-4,
        ).unwrap();
        let low_mach = bounded.resolve(cx, &initial, &[0.0; 3]).unwrap();
        let inertia = inventory.reference_density().value() / AREA;
        let drive = -2.0 * inertia * 1e-4 / 0.1;
        let crossing = native_report(
            cx, &low_mach.graph, &inventory, [1e-4, 0.0, 1e-4],
            [drive, 0.0, drive], [0.0; 3], 0.1,
        );
        assert!(matches!(bounded.evaluate_reservoirs(
            cx, &initial, &low_mach.sources, &crossing, &reservoirs, 1,
        ), Err(GraphReservoirConvectionError::Source(SourceGraphError::GasMach { .. }))));
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 300.0);
    });
}
