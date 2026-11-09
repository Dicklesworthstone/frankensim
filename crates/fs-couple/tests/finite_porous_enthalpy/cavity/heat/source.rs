//! Source transport acts on actual gas, wall and retained compaction histories.
use super::*;
use fs_couple::air_path::input::convection::{
    GasConvectionBoundary, NaturalConvectionGeometry, NaturalGasBoundary,
};
use fs_couple::solid_enthalpy::finite_plastic::porous::cavity::heat::CavityWallSourceContact;
use fs_couple::solid_enthalpy::gas::{
    GasConductionGeometry, SolidGasConductivityProfile, SolidGasContact, SolidGasConvection,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_qty::Temperature;

const AREA_PER_CAPACITY: f64 = 4e-4;
const GAS_DISTANCE: f64 = 2e-5;
const SOLID_DISTANCE: f64 = 1e-3;
const SPHERE_DIAMETER: f64 = 0.03;

fn source_fixture(k300: f64, slope: f64, upper: f64, coupled: bool) -> Fixture {
    let original = if coupled {
        Fixture::new(CP, 31)
    } else {
        thermally_constant_fixture()
    };
    let mut claims = original.card.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (280.0, k300 - 20.0 * slope),
                    (upper, k300 + (upper - 300.0) * slope),
                ],
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "temperature",
                temperature_quantity(),
                280.0,
                upper,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: provenance(),
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic source-conducting porous wall".into(),
            phase: "solid".into(),
            process: "numerical regression".into(),
            revision: 0,
        },
        claims,
        original.card.models().to_vec(),
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP_NAME, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "temperature",
            lower: &source_context(280.0),
            upper: &source_context(420.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * CP,
        1e-6,
    )
    .unwrap();
    Fixture { card, caloric }
}

fn contacts(thermal: &DepositionSimulation<'_>, convection: bool) -> Vec<CavityWallSourceContact> {
    thermal
        .cells()
        .values()
        .enumerate()
        .flat_map(|(i, cell)| {
            [0.4, 0.6]
                .into_iter()
                .enumerate()
                .map(move |(patch, fraction)| {
                    let area = fraction * AREA_PER_CAPACITY * cell.mass_kg * CP;
                    CavityWallSourceContact {
                        cavity: 0,
                        contact: SolidGasContact {
                            id: (1000 + 2 * i + patch) as u64,
                            cell: cell.id,
                            geometry: GasConductionGeometry {
                                area_m2: area,
                                gas_distance_m: if convection { 0.0 } else { GAS_DISTANCE },
                                cell_distance_m: if convection { 0.01 } else { SOLID_DISTANCE },
                            },
                            convection: convection.then_some(SolidGasConvection {
                                boundary: GasConvectionBoundary::Natural(NaturalGasBoundary {
                                    geometry: NaturalConvectionGeometry::Sphere,
                                    characteristic_length_m: SPHERE_DIAMETER,
                                    exposed_area_m2: area,
                                    gravity_m_s2: 0.0,
                                }),
                                solid_profile: SolidGasConductivityProfile::IntegratedPurePhase {
                                    phase: SolidLiquidPhase::Solid,
                                    temperature_range_k: [280.0, 420.0],
                                },
                                maximum_surface_iterations: 96,
                                surface_power_tolerance_w: 1e-11,
                            }),
                            conductivity: PropertyKey::new(
                                "thermal_conductivity",
                                THERMAL_CONDUCTIVITY_DIMS,
                            ),
                            selection: MaterialPropertySelection::SingleClaimOnly,
                        },
                    }
                })
        })
        .collect()
}

fn gas_k(temperature: f64) -> f64 {
    // Independent interpolation of the original synthetic transport claim.
    0.025 * (1.0 + (temperature - 150.0) / 1850.0)
}

#[test]
fn sourced_series_conduction_uses_actual_gas_temperature_and_shared_native_decay() {
    context(|cx| {
        let fixture = source_fixture(1.0, 0.0, 420.0, false);
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let cold = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.0, 1e24);
        let (gas, _) = cold.at_temperature(360.0).unwrap();
        let capacity =
            gas.mass_kg() * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.028 * 2.5;
        let original = gas.clone();
        let mut patches = contacts(&thermal, false);
        patches.reverse();
        let system = ThermalPorousGasSystem::try_new(
            solid,
            &thermal,
            vec![cavity_entry(Arc::new(cavity), gas)],
            links(&thermal, 0.03),
            cx,
        )
        .unwrap()
        .with_source_contacts(&thermal, patches, cx)
        .unwrap();
        let next = system
            .advance(
                &thermal,
                &Times::default(),
                1.0,
                &affine(&nodes, 1.0),
                &[],
                [0.0; 3],
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                &[0.02],
                &[0.0],
                solve_heat(),
                cx,
            )
            .unwrap();
        assert!(next.observer.0.len() > 4);
        let wall_capacity: f64 = thermal.cells().values().map(|cell| cell.mass_kg * CP).sum();
        let wall = |temperature: f64, frozen: bool| {
            let coefficient = gas_k(if frozen { 360.0 } else { temperature });
            let rate = 0.03 + AREA_PER_CAPACITY / (GAS_DISTANCE / coefficient + SOLID_DISTANCE);
            let decay: f64 = next.observer.0.iter().map(|dt| 1.0 - rate * dt).product();
            300.0 + (temperature - 300.0) * (1.0 - decay)
        };
        let root = |frozen: bool| {
            let (mut lo, mut hi) = (300.0, 360.0);
            for _ in 0..80 {
                let t = 0.5 * (lo + hi);
                let residual =
                    capacity * (t - 360.0) + wall_capacity * (wall(t, frozen) - 300.0) - 0.02;
                if residual < 0.0 {
                    lo = t;
                } else {
                    hi = t;
                }
            }
            0.5 * (lo + hi)
        };
        let expected = root(false);
        assert!((root(true) - expected).abs() > 1e-3);
        close(
            next.system.cavities()[0].gas.state().state().temperature,
            expected,
            3e-7,
        );
        for cell in thermal.cells().values() {
            close(
                next.thermal.temperature_k(cell.id).unwrap(),
                wall(expected, false),
                3e-7,
            );
        }
        assert_eq!(next.balance.source_contacts.len(), 2 * tets.len());
        assert!(
            next.balance
                .source_contacts
                .windows(2)
                .all(|p| p[0].id < p[1].id)
        );
        let source_heat: f64 = next
            .balance
            .source_contacts
            .iter()
            .map(|r| r.heat_into_solid_j)
            .sum();
        let fixed_heat: f64 = next.balance.links.iter().map(|r| r.heat_into_solid_j).sum();
        close(
            source_heat + fixed_heat,
            next.balance.gas_to_solid_heat_j,
            2e-10,
        );
        close(
            next.balance.total_gas_energy_change_j,
            capacity * (expected - 360.0),
            3e-9,
        );
        close(
            next.balance.gas_to_solid_heat_j,
            wall_capacity * (wall(expected, false) - 300.0),
            3e-9,
        );
        close(next.balance.external_solid_heat_j, 0.0, 1e-9);
        close(
            next.balance.combined_residual_j,
            0.0,
            solve_heat().energy_tolerance_j,
        );
        for row in &next.balance.source_contacts {
            let patch = contacts(&thermal, false)
                .into_iter()
                .find(|p| p.contact.id == row.id)
                .unwrap();
            let expected_g =
                patch.contact.geometry.area_m2 / (GAS_DISTANCE / gas_k(expected) + SOLID_DISTANCE);
            close(row.last_conductance_w_k, expected_g, 2e-12);
            close(row.last_explicit_decay_conductance_w_k, expected_g, 2e-12);
            close(
                row.last_power_into_solid_w,
                expected_g * (expected - row.last_cell_temperature_k),
                2e-10,
            );
            assert_eq!(row.source.card_identity(), fixture.card.content_hash());
            assert!(row.convection.is_none());
        }
        assert_eq!(system.cavities()[0].gas, original);
        assert_eq!(thermal.time_s(), 0.0);
    });
}

#[test]
fn integrated_source_convection_survives_heated_porous_compression_epochs() {
    context(|cx| {
        let fixture = source_fixture(0.06, 0.001, 420.0, true);
        let (nodes, tets) = shell();
        let mut thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let cold = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.001, 1e24);
        let (gas, _) = cold.at_temperature(340.0).unwrap();
        let original = gas.clone();
        let patches = contacts(&thermal, true);
        let mut system = ThermalPorousGasSystem::try_new(
            solid,
            &thermal,
            vec![cavity_entry(Arc::new(cavity), gas)],
            vec![],
            cx,
        )
        .unwrap()
        .with_source_contacts(&thermal, patches.clone(), cx)
        .unwrap();
        let mut previous_profile = None;
        for (time, stretch) in [(1.0, 0.94), (2.0, 0.90)] {
            let before = system.cavities()[0].gas.clone();
            let old_states = system.solid().mechanics().element_states().to_vec();
            let next = system
                .advance(
                    &thermal,
                    &(),
                    time,
                    &affine(&nodes, stretch),
                    &[],
                    [0.0; 3],
                    None,
                    cavity_settings(),
                    controls(),
                    thermal_config(),
                    config(),
                    CouplingMethod::RelaxedPicard,
                    &[0.015],
                    &[0.0],
                    solve_heat(),
                    cx,
                )
                .unwrap();
            let stage = next.balance.gases[0].heat_stage_temperature_k;
            let pressure =
                before.state().state().pressure * stage / before.state().state().temperature;
            for (patch, row) in patches.iter().zip(&next.balance.source_contacts) {
                let film = row.convection.as_ref().unwrap();
                let cell = row.last_cell_temperature_k;
                let surface = film.film.wall_temperature_k;
                let area = patch.contact.geometry.area_m2;
                let gas_g = area * 2.0 * gas_k(0.5 * (stage + surface)) / SPHERE_DIAMETER;
                // Independent integral of the authored affine k(T), not a
                // frozen cell coefficient and not the native layer solver.
                let solid_power = area / patch.contact.geometry.cell_distance_m
                    * (surface - cell)
                    * (0.06 + 0.0005 * (surface + cell - 600.0));
                close(
                    row.last_power_into_solid_w,
                    gas_g * (stage - surface),
                    2e-10,
                );
                close(row.last_power_into_solid_w, solid_power, 3e-10);
                close(
                    row.last_explicit_decay_conductance_w_k,
                    row.last_power_into_solid_w / (stage - cell),
                    2e-12,
                );
                close(
                    film.film.nusselt.as_ref().unwrap().evidence().value,
                    2.0,
                    0.0,
                );
                close(film.film.property_gas.state().pressure, pressure, 2e-5);
                close(film.film.reservoir_temperature_k, stage, 0.0);
                assert!(film.surface_balance_residual_w.abs() <= 1e-11);
                assert!(film.solid_layer_entropy_generation_w_k >= 0.0);
                assert_eq!(row.source.card_identity(), fixture.card.content_hash());
            }
            let profile = next.balance.source_contacts[0]
                .convection
                .as_ref()
                .unwrap()
                .solid_conductivity_profile
                .as_ref()
                .unwrap();
            if let Some(previous) = &previous_profile {
                assert!(Arc::ptr_eq(previous, profile));
            }
            previous_profile = Some(Arc::clone(profile));
            close(
                next.balance
                    .source_contacts
                    .iter()
                    .map(|row| row.heat_into_solid_j)
                    .sum(),
                next.balance.gas_to_solid_heat_j,
                2e-9,
            );
            close(next.balance.external_solid_heat_j, 0.0, 2e-8);
            close(
                next.balance.combined_residual_j,
                0.0,
                solve_heat().energy_tolerance_j,
            );
            let gas = &next.system.cavities()[0].gas;
            assert_eq!(gas.mass_kg().to_bits(), original.mass_kg().to_bits());
            assert_eq!(gas.model(), original.model());
            assert!(gas.volume_m3() < before.volume_m3());
            for (old, new) in old_states
                .iter()
                .zip(next.system.solid().mechanics().element_states())
            {
                assert!(new.alpha >= old.alpha);
            }
            thermal = next.thermal;
            system = next.system;
        }
        assert_eq!(system.time_s(), 2.0);
        assert!(
            system
                .solid()
                .mechanics()
                .element_states()
                .iter()
                .all(|state| state.alpha > 0.0)
        );
    });
}

#[derive(Clone, Default)]
struct Heater {
    calls: usize,
}
impl DepositionObserver for Heater {
    type Error = DepositionError;
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        self.calls += 1;
        Ok(thermal
            .cells()
            .values()
            .map(|cell| (cell.id, cell.mass_kg * CP))
            .collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        Ok(())
    }
}

#[test]
fn late_original_conductivity_refusal_rolls_back_wall_gas_and_observer() {
    context(|cx| {
        let fixture = source_fixture(1.0, 0.0, 300.1, false);
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
        let (_, cavity) = boundaries(&nodes, &tets, cx);
        let gas = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.0, 0.0);
        let original = gas.clone();
        let system = ThermalPorousGasSystem::try_new(
            solid,
            &thermal,
            vec![cavity_entry(Arc::new(cavity), gas)],
            vec![],
            cx,
        )
        .unwrap()
        .with_source_contacts(&thermal, contacts(&thermal, false), cx)
        .unwrap();
        let old_cells = thermal.cells().clone();
        let old_states = system.solid().mechanics().element_states().to_vec();
        let observer = Heater::default();
        let refused = system.advance(
            &thermal,
            &observer,
            1.0,
            &affine(&nodes, 1.0),
            &[],
            [0.0; 3],
            None,
            cavity_settings(),
            controls(),
            thermal_config(),
            config(),
            CouplingMethod::RelaxedPicard,
            &[0.0],
            &[0.0],
            solve_heat(),
            cx,
        );
        assert!(matches!(
            refused,
            Err(CavityHeatError::Feedback(
                FiniteThermoplasticError::Observer(CavityHeatObserverError::Source(_))
            ))
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(system.cavities()[0].gas, original);
        assert_eq!(system.solid().mechanics().element_states(), old_states);
        assert_eq!(observer.calls, 0);
        let replay = system
            .advance(
                &thermal,
                &(),
                1.0,
                &affine(&nodes, 1.0),
                &[],
                [0.0; 3],
                None,
                cavity_settings(),
                controls(),
                thermal_config(),
                config(),
                CouplingMethod::RelaxedPicard,
                &[0.0],
                &[0.0],
                solve_heat(),
                cx,
            )
            .unwrap();
        close(replay.balance.gas_to_solid_heat_j, 0.0, 0.0);
        assert_eq!(replay.system.cavities()[0].gas, original);
    });
}

#[test]
fn source_patch_area_and_cross_mode_identity_are_admitted_before_stepping() {
    context(|cx| {
        let fixture = source_fixture(1.0, 0.0, 420.0, false);
        let (nodes, tets) = shell();
        let thermal = original_thermal(&fixture, &nodes, &tets);
        let make = || {
            let solid = bind_shell(&fixture, &thermal, &nodes, &tets, cx);
            let (_, cavity) = boundaries(&nodes, &tets, cx);
            let gas = gas_fixture::gas(cavity.reference_volume_m3(), 5e6, 0.0, 0.0);
            ThermalPorousGasSystem::try_new(
                solid,
                &thermal,
                vec![cavity_entry(Arc::new(cavity), gas)],
                links(&thermal, 0.01),
                cx,
            )
            .unwrap()
        };
        let mut patch = contacts(&thermal, false).remove(0);
        patch.contact.id = 0;
        assert!(matches!(
            make().with_source_contacts(&thermal, vec![patch.clone()], cx),
            Err(CavityHeatError::Invalid(_))
        ));
        patch.contact.id = 1000;
        patch.contact.geometry.area_m2 = 2.0 * thermal.exposed_area_m2(patch.contact.cell).unwrap();
        assert!(matches!(
            make().with_source_contacts(&thermal, vec![patch], cx),
            Err(CavityHeatError::Invalid(_))
        ));
        assert_eq!(thermal.time_s(), 0.0);
    });
}
