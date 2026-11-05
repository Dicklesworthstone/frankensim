//! Source-backed gas spheres: bulk transport, original wall viscosity, no film substitution.
use super::*;
use fs_airflow::graph::thermal::transport::transient::material::convection::reservoir::{
    GraphConvectiveReservoir, GraphReservoirConvectionError,
};
use fs_conduction::transient::deposition::Reservoir;

fn variable_gas_law() -> GasEnthalpyLaw {
    // Every constituent has the same positive mu(T)/mu(300)=T/300 factor.
    // The original mixture law therefore retains that factor exactly, while
    // density and NASA Cp change independently with the actual query state.
    let transport = |name: &str, mu: f64, conductivity: f64| {
        let mut claims = ClaimSet::new();
        for (key, value, interpolation) in [
            (
                PropertyKey::new("dynamic_viscosity", DynViscosity::DIMS),
                PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0, mu * 250.0 / 300.0), (450.0, mu * 450.0 / 300.0)],
                    dims: DynViscosity::DIMS,
                },
                InterpolationPolicy::LinearInside,
            ),
            (
                PropertyKey::with_quantity(
                    "thermal_conductivity",
                    typed(QuantityKind::ThermalConductivity),
                ),
                PropertyValue::Scalar {
                    value: conductivity,
                    dims: fs_material::state_point::THERMAL_CONDUCTIVITY_DIMS,
                },
                InterpolationPolicy::ConstantWithinValidity,
            ),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key,
                    value,
                    interpolation,
                    validity: ValidityDomain::unconstrained()
                        .with_quantity(
                            "temperature",
                            typed(QuantityKind::AbsoluteTemperature),
                            250.0,
                            450.0,
                        )
                        .with_quantity(
                            "pressure",
                            typed(QuantityKind::Pressure),
                            100000.0,
                            200000.0,
                        ),
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: provenance(),
                })
                .unwrap();
        }
        SpeciesTransportMaterial {
            species: SpeciesId::new(name).unwrap(),
            material: MaterialCard::assemble(
                MaterialStateId {
                    chemistry: name.into(),
                    phase: "gas".into(),
                    process: "synthetic sphere temperature law".into(),
                    revision: 0,
                },
                claims,
                vec![],
            )
            .unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }
    };
    let composition = Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap();
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species("A", 0.028, 3.5, SLOPE, -100000.0),
                species("B", 0.004, 2.5, 0.0, 300000.0),
            ],
            composition.clone(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)],
        None,
    )
    .unwrap();
    GasEnthalpyLaw::from_mole(
        &model,
        composition,
        point(350.0),
        GasEnergySolve {
            minimum_temperature_k: 280.0,
            maximum_temperature_k: 360.0,
            energy_tolerance_j_kg: 1e-8,
            maximum_iterations: 64,
        },
    )
    .unwrap()
}

fn fixture(law: &GasEnthalpyLaw, wall_k: f64) -> DepositionSimulation<'_> {
    let original = thermal(law, [350.0; 6], Some(wall_k));
    let mut schedule: Vec<_> = original
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            if cell.id == CellId(3) {
                cell.boundary.reservoir = Some(Reservoir {
                    area_m2: 0.002,
                    conductance_w_k: 0.0,
                    temperature_k: wall_k,
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
    DepositionSimulation::new(
        schedule,
        &[(MaterialId(1), law as &dyn ThermalLaw)],
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
    .0
}

#[test]
fn g1_g4_gas_spheres_retain_wall_viscosity_without_a_second_storage_density() {
    context(|cx| {
        let law = variable_gas_law();
        let initial = fixture(&law, 280.0);
        let graph = graph();
        let inventory = bind_inventory(cx, &graph, &law, &initial);
        let laws = laws();
        let limits = source_limits();
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
            limits,
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
        let diameter = 0.04;
        let link = GraphConvectiveLink {
            link: 11,
            fluid_cell: CellId(2),
            characteristic_length: Length::new(diameter),
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: GraphConvectionLaw::ForcedSphere { branch: 0 },
        };
        let reservoir = GraphConvectiveReservoir {
            fluid_cell: CellId(3),
            characteristic_length: Length::new(diameter),
            property_temperature: ConvectionPropertyTemperature::Reservoir,
            law: GraphConvectionLaw::ForcedSphere { branch: 0 },
        };
        let links = owner
            .evaluate_links(cx, &initial, &resolved.sources, &hydraulics, &[link], 1)
            .unwrap();
        let reservoirs = owner
            .evaluate_reservoirs(
                cx,
                &initial,
                &resolved.sources,
                &hydraulics,
                &[reservoir],
                1,
            )
            .unwrap();
        let finite = &links.convection.samples()[0];
        let prescribed = &reservoirs.convection.samples()[0];
        for (bulk, thermal, surface, evaluation) in [
            (
                &finite.material,
                &finite.thermal,
                finite.surface.as_ref().unwrap(),
                finite.evaluation.as_ref().unwrap(),
            ),
            (
                &prescribed.material,
                &prescribed.thermal,
                prescribed.surface.as_ref().unwrap(),
                prescribed.evaluation.as_ref().unwrap(),
            ),
        ] {
            assert_eq!(*bulk, law.model().evaluate(&point(350.0)).unwrap());
            assert_eq!(
                surface.material,
                law.model().evaluate(&point(280.0)).unwrap()
            );
            near(thermal.temperature_k(), 350.0);
            near(surface.thermal.temperature_k(), 280.0);
            near(bulk.state().pressure, PRESSURE);
            near(surface.material.state().pressure, PRESSURE);
            near(bulk.specific_enthalpy_j_kg(), h(350.0));
            let transport = thermal.transport();
            let ratio = transport.dynamic_viscosity_pa_s()
                / surface.thermal.transport().dynamic_viscosity_pa_s();
            near(ratio, 1.25);
            near(evaluation.nusselt.groups()["mu_bulk_over_surface"], ratio);
            // This large surface-density departure has no mass inventory. Only
            // its original viscosity enters Whitaker; bulk storage stays fixed.
            let surface_departure =
                surface.material.state().density / owner.reference().state().density - 1.0;
            near(surface_departure, 0.25);
            assert!(surface_departure > limits.maximum_relative_density_departure);
            let re = transport.density_kg_m3() * (q / AREA) * diameter
                / transport.dynamic_viscosity_pa_s();
            let pr = transport.prandtl();
            let nu = 2.0
                + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0)) * pr.powf(0.4) * ratio.powf(0.25);
            near(evaluation.nusselt.groups()["Re"], re);
            near(evaluation.nusselt.groups()["Pr"], pr);
            near(
                evaluation.coefficient.value.value(),
                nu * transport.thermal_conductivity_w_m_k() / diameter,
            );
        }
        assert_eq!(owner.source().mole_composition().fractions(), &[0.75, 0.25]);
        near(finite.property_temperature_k, 350.0);
        near(prescribed.property_temperature_k, 350.0);
        near(finite.conductance_w_k, prescribed.conductance_w_k);
        near(finite.heat_into_fluid_w, -70.0 * finite.conductance_w_k);
        near(prescribed.heat_into_fluid_w, finite.heat_into_fluid_w);
        // Heating this gas from a hotter wall gives mu_bulk/mu_wall < 1.
        // The original card refuses; no unit ratio, film rule or clamp appears.
        let hot = fixture(&law, 360.0);
        let hot_sources = owner.resolve(cx, &hot, &flows).unwrap();
        assert!(matches!(
            owner.evaluate_reservoirs(cx, &hot, &hot_sources.sources, &hydraulics, &[reservoir], 1),
            Err(GraphReservoirConvectionError::Convection(_))
        ));
        near(initial.time_s(), 0.0);
        near(initial.temperature_k(CellId(2)).unwrap(), 350.0);
        near(initial.temperature_k(CellId(6)).unwrap(), 280.0);
    });
}
