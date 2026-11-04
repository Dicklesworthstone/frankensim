//! Original finite-solid rows, sourced contacts and one reciprocal atmosphere.
use super::*;
use fs_conduction::enthalpy::CellGeometry;
use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity};
use fs_evidence::Evidence;
use fs_lbm::material::free_surface::enthalpy::conduction::surface::{
    FreeSurfaceImplicitSolve,
    convection::{
        FreeSurfaceConvectionQuery, FreeSurfaceGasConvection, FreeSurfaceHeatTransferCoefficient,
    },
};
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceBodyGasFace, FreeSurfaceBodyLink, FreeSurfaceBodyRadiationFace,
    FreeSurfaceThermalBody,
};
use fs_material::fluid_transport::FluidTransport;
use fs_material::thermal::source::{
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, ThermalMaterialPoint,
};
use fs_qty::Length;
use std::cell::Cell;

const WET: [usize; 3] = [0, 3, 3];
const R_GAS: f64 = 8.314_462_618_153_24 / 0.028;

fn implicit() -> FreeSurfaceInterfaceSolve {
    FreeSurfaceInterfaceSolve {
        implicit: Some(FreeSurfaceImplicitSolve::default()),
        ..FreeSurfaceInterfaceSolve::default()
    }
}

fn solid_face() -> FreeSurfaceBodyGasFace {
    FreeSurfaceBodyGasFace {
        wall: [-1, 0, 3],
        gas: [0, 0, 3],
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }
}

fn emissive_wall() -> MaterialCard {
    let base = body::wall_card(0.5, 500.0);
    let mut claims = base.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.2), (500.0, 0.8)],
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured finite-wall emissivity".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn rear_chart(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    // Cp=50 J/(kg K), ten times smaller than the front wall's original chart.
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 15000.0,
                temperature_k: 600.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2000.0,
            },
        ],
    )
    .unwrap()
}

fn rear<'a>(
    card: &'a MaterialCard,
    chart: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint,
) -> FreeSurfaceThermalBody<'a> {
    FreeSurfaceThermalBody::try_new(
        43,
        ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        chart,
        CellGeometry {
            start_m: [-0.002, 0.0035, 0.0035],
            end_m: [-0.001, 0.0035, 0.0035],
            volume_m3: 1e-9,
            surface_area_m2: 6e-6,
        },
        (480.0 - 300.0) * 50.0,
        &PropertyKey::with_quantity(
            "thermal_conductivity",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
        ),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn solid_link() -> FreeSurfaceBodyLink {
    FreeSurfaceBodyLink {
        other_body_id: 43,
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }
}

fn temperature_of(body: &FreeSurfaceThermalBody<'_>) -> f64 {
    body.material().query().axes()["T"]
}

fn gas_storage(mass: f64, start: f64, end: f64) -> f64 {
    // Independent NASA9 cv/R = 2 + .001 T; absolute formation energy cancels.
    mass * R_GAS * (end - start) * (2.0 + 0.0005 * (end + start))
}

fn contact_card(
    liquid: &MaterialCard,
    solid: &MaterialCard,
    scale: f64,
    liquid_upper_k: f64,
) -> fs_matdb::InterfaceSystemCard {
    let base = body::interface_card(liquid, solid, scale, liquid_upper_k, false);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        // This original resistance curve is queried at a typed mean kelvin
        // state. Its support must declare that quantity, as do both bulk axes.
        claim.validity = claim.validity.with_quantity(
            "meanT",
            QuantitySpec::dimensional(Temperature::DIMS),
            300.0,
            500.0,
        );
        claims.insert_claim(claim).unwrap();
    }
    fs_matdb::InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        base.models().to_vec(),
    )
    .unwrap()
}

#[test]
fn stiff_liquid_solid_gas_network_closes_every_original_row_and_source() {
    let liquid = conducting_card(1.0);
    let front = emissive_wall();
    let back = body::wall_card(20.0, 500.0);
    let front_curve = body::wall_chart(&front);
    let back_curve = rear_chart(&back);
    let front_point = point(450.0);
    let back_point = point(480.0);
    let interface = contact_card(&liquid, &front, 1e-3, 500.0);
    let interface_point = body::interface_point();
    let binding = body::interface(&liquid, &front, &interface, &interface_point);
    let run = |implicit_step| {
        with_exposed_liquid(&liquid, false, |thermal, gate| {
            let mut gas = IsochoricGasParcel::new(
                gas_source::model(1e12, 250000.0, 100.0),
                gas_source::point(350.0, 100000.0),
                1e-10,
            )
            .unwrap();
            let patch = solid_face();
            let mut bodies = [
                body::wall(&front, &front_curve, &front_point, 450.0)
                    .with_conduction_links(vec![solid_link()])
                    .with_gas_faces(vec![patch])
                    .with_radiation_faces(
                        vec![FreeSurfaceBodyRadiationFace {
                            wall: patch.wall,
                            gas: patch.gas,
                            area_m2: patch.area_m2,
                            surroundings_k: 1000.0,
                        }],
                        MaterialPropertySelection::SingleClaimOnly,
                    )
                    .unwrap(),
                rear(&back, &back_curve, &back_point),
            ];
            let liquid_before = thermal.total_enthalpy_j().unwrap();
            let solids_before: f64 = bodies.iter().map(|b| b.enthalpy_j()).sum();
            let authored = 1e-5;
            let report = thermal
                .step_with_finite_gas_surface_exchange(
                    &[],
                    &[FreeSurfaceCellHeat {
                        coordinates: EXPOSED,
                        energy_j: authored,
                    }],
                    FreeSurfaceConductionControls {
                        maximum_dt_s: 0.001,
                        ..controls()
                    },
                    &mut bodies,
                    &[body::contact()],
                    std::slice::from_ref(&binding),
                    &mut gas,
                    &[face(900.0)],
                    gas_controls(),
                    if implicit_step {
                        implicit()
                    } else {
                        FreeSurfaceInterfaceSolve::default()
                    },
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            if !implicit_step {
                return report.thermal.substeps;
            }

            let dt = 0.001;
            let exchange = report.surface.as_ref().unwrap();
            let wall = report.bodies.as_ref().unwrap();
            let film = report.body_environment.as_ref().unwrap();
            let radiation = report.body_radiation.as_ref().unwrap();
            let tw = temperature(&wall.thermal_endpoint[0].liquid);
            let tl = exchange.thermal_endpoint[0].liquid_temperature_k;
            let ta = temperature_of(&bodies[0]);
            let tb = temperature_of(&bodies[1]);
            let tg = gas.state().state().temperature;
            let pressure = 100000.0 * tg / 350.0;
            let kg = 2.6 * (0.5 + (pressure - 1000.0) / 249000.0);
            let ka = 0.5 * (ta - 200.0);
            let kb = 20.0 * (tb - 200.0);
            let resistance = 1e-3 * (1.0 + (0.5 * (tw + ta) - 300.0) / 200.0);
            let qaw = 1e-6 * (ta - tw) / (0.0005 / (tw - 200.0) + 0.0005 / ka + resistance);
            let qba = 1e-6 * (tb - ta) / (0.0005 / ka + 0.0005 / kb);
            let qga = 1e-6 * (tg - ta) / (0.0005 / ka + 0.0005 / kg);
            let qra = (0.2 + 0.003 * (ta - 300.0))
                * 5.670374419e-8
                * 1e-6
                * (1000.0_f64.powi(4) - ta.powi(4));
            let (_, qgl, qrl) = reference_surface(tl, tg, 900.0, 0.8, tl - 200.0, |_| kg);
            let gas_change = gas_storage(gas.mass_kg(), 350.0, tg);
            // Independent physical BE equations for the two liquid nodes, two
            // DISTINCT finite capacities, and the one original caloric gas.
            for defect in [
                0.002 * (tw - 400.0) - dt * qaw,
                0.002 * (tl - 400.0) - authored - dt * (qgl + qrl),
                0.001 * (ta - 450.0) - dt * (qba - qaw + qga + qra),
                0.0001 * (tb - 480.0) + dt * qba,
                gas_change + dt * (qgl + qga),
            ] {
                assert!(defect.abs() < 1e-9, "independent coupled row: {defect}");
            }
            let stored = thermal.total_enthalpy_j().unwrap() - liquid_before
                + bodies.iter().map(|b| b.enthalpy_j()).sum::<f64>()
                - solids_before
                + gas_change;
            assert!((stored - authored - dt * (qra + qrl)).abs() < 1e-9);
            assert!((gas_change + exchange.heat_from_gas_j + film.heat_to_solids_j).abs() < 1e-9);
            assert!(exchange.heat_from_gas_j != 0.0 && film.heat_to_solids_j != 0.0);
            assert!(radiation.heat_to_solids_j > 0.0);
            near(radiation.heat_to_solids_j, dt * qra, 1e-9);
            near(gas.state().state().pressure, pressure, 2e-12);
            assert!(tw > 400.0 && ta > 450.0 && tb < 480.0 && tg > 350.0);
            let contact = &wall.thermal_endpoint[0];
            assert_eq!(&contact.body, bodies[0].material());
            assert_eq!(&film.thermal_endpoint[0].solid, bodies[0].material());
            assert_eq!(
                wall.solid_thermal_endpoint[0].materials[1],
                *bodies[1].material()
            );
            let source = contact.interface.as_ref().unwrap();
            assert_eq!(source.card_identity(), interface.content_hash());
            near(source.state_coordinate("liquidT").unwrap(), tw, 2e-12);
            near(source.state_coordinate("solidT").unwrap(), ta, 2e-12);
            near(
                source.state_coordinate("meanT").unwrap(),
                0.5 * (tw + ta),
                2e-12,
            );
            assert_eq!(
                film.thermal_endpoint[0].gas,
                gas.model()
                    .evaluate(&gas_source::point(
                        0.5 * ta + 0.5 * tg,
                        gas.state().state().pressure,
                    ))
                    .unwrap()
            );
            assert_eq!(wall.last_applied, wall.thermal_endpoint);
            assert_eq!(wall.solid_last_applied, wall.solid_thermal_endpoint);
            let work = exchange.implicit.as_ref().unwrap();
            assert_eq!(work.intervals, 1);
            assert!(
                work.maximum_cell_residual_j
                    <= FreeSurfaceImplicitSolve::default().cell_energy_tolerance_j
            );
            assert!(
                work.maximum_gas_residual_j
                    <= FreeSurfaceImplicitSolve::default().cell_energy_tolerance_j
            );
            assert!(
                report
                    .finite_gas
                    .as_ref()
                    .unwrap()
                    .combined_thermal_balance_defect_j
                    .abs()
                    < 1e-9
            );
            report.thermal.substeps
        })
    };
    assert_eq!(run(true), 1);
    assert!(
        run(false) >= 10,
        "original explicit capacities must require substantially more intervals"
    );
}

#[test]
fn prescribed_surface_keeps_solid_contact_inside_the_original_implicit_owner() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    let interface = contact_card(&liquid, &solid, 1e-5, 500.0);
    let interface_point = body::interface_point();
    let binding = body::interface(&liquid, &solid, &interface, &interface_point);
    let gas = gas_source::model(0.0, 250000.0, 100.0);
    let gas_point = gas_source::point(450.0, 100000.0);
    with_exposed_liquid(&liquid, false, |thermal, gate| {
        let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
        let before = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j();
        let mut insulated = face(400.0);
        insulated.face.area_m2 = 0.0;
        let report = thermal
            .step_with_surface_exchange(
                &[],
                &[],
                FreeSurfaceConductionControls {
                    maximum_dt_s: 0.001,
                    ..controls()
                },
                &mut bodies,
                &[body::contact()],
                std::slice::from_ref(&binding),
                &[FreeSurfaceInterfaceContact {
                    surface: insulated,
                    model: &gas,
                    context: &gas_point,
                }],
                implicit(),
                false,
                gate,
            )
            .unwrap()
            .conduction;
        // Energy conservation eliminates the liquid state; solve the remaining
        // nonlinear source-dependent contact equation independently by bisection.
        let (mut lo, mut hi) = (400.0, 450.0);
        for _ in 0..80 {
            let ts = 0.5 * (lo + hi);
            let tl = 400.0 + 0.5 * (450.0 - ts);
            let resistance = 1e-5 * (1.0 + (0.5 * (tl + ts) - 300.0) / 200.0);
            let q = 1e-6 * (ts - tl)
                / (0.0005 / (tl - 200.0) + 0.0005 / (0.5 * (ts - 200.0)) + resistance);
            if 0.001 * (ts - 450.0) + 0.001 * q < 0.0 {
                lo = ts;
            } else {
                hi = ts;
            }
        }
        let ts = 0.5 * (lo + hi);
        near(temperature_of(&bodies[0]), ts, 1e-9);
        near(
            temperature(thermal.material_at(WET).unwrap()),
            400.0 + 0.5 * (450.0 - ts),
            1e-9,
        );
        assert!(report.finite_gas.is_none());
        assert_eq!(report.thermal.external_heat_j, 0.0);
        assert_eq!(
            report
                .surface
                .as_ref()
                .unwrap()
                .implicit
                .as_ref()
                .unwrap()
                .intervals,
            1
        );
        assert!(
            (thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j() - before).abs() < 1e-9
        );
    });
}

#[derive(Default)]
struct ObservedPlate<'a> {
    maximum_gas_k: Cell<f64>,
    cancel: Option<&'a CancelGate>,
}
impl FreeSurfaceGasConvection for ObservedPlate<'_> {
    fn coefficient(
        &self,
        q: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        self.maximum_gas_k
            .set(self.maximum_gas_k.get().max(q.bulk.state().temperature));
        if q.bulk.state().temperature > 351.0 {
            if let Some(gate) = self.cancel {
                gate.request();
            }
        }
        let state = q
            .film
            .transport_state()
            .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))?;
        let re = state
            .reynolds(10.0, 0.01)
            .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))?;
        fs_convection::evaluate(
            CorrelationId::FlatPlateLaminarAverage,
            CorrelationInputs::forced(re, state.prandtl()),
        )
        .and_then(|nu| {
            nu.heat_transfer_coefficient(
                ThermalConductivity::new(state.thermal_conductivity_w_m_k()),
                Length::new(0.01),
            )
        })
        .map(Some)
        .map_err(|e| FreeSurfaceEnthalpyError::Convection(e.to_string()))
    }
}

fn liquid_states(thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> Vec<Option<f64>> {
    let mut states = Vec::new();
    for z in 1..=6 {
        for y in 1..=6 {
            for x in 0..=3 {
                states.push(thermal.specific_enthalpy_at([x, y, z]));
            }
        }
    }
    states
}

#[test]
fn late_contact_and_cancellation_refusals_preserve_liquid_each_body_and_gas() {
    let liquid = conducting_card(1.0);
    let solid = body::wall_card(0.5, 500.0);
    let back = body::wall_card(20.0, 500.0);
    let curve = body::wall_chart(&solid);
    let back_curve = rear_chart(&back);
    let context = point(450.0);
    let back_point = point(480.0);
    let interface = contact_card(&liquid, &solid, 1e-3, 500.0);
    let narrow = contact_card(&liquid, &solid, 1e-3, 400.004);
    let point = body::interface_point();
    let binding = body::interface(&liquid, &solid, &interface, &point);
    let restricted = body::interface(&liquid, &solid, &narrow, &point);
    let run = |refuse_first| {
        with_exposed_liquid(&liquid, true, |thermal, gate| {
            let mut gas = IsochoricGasParcel::new(
                gas_source::model(0.0, 250000.0, 1.0),
                gas_source::point(350.0, 100000.0),
                1e-10,
            )
            .unwrap();
            let mut bodies = [
                body::wall(&solid, &curve, &context, 450.0)
                    .with_conduction_links(vec![solid_link()])
                    .with_gas_faces(vec![solid_face()]),
                rear(&back, &back_curve, &back_point),
            ];
            if refuse_first {
                let grid = thermal.flow().lattice().grid().state_bits();
                let enthalpies = liquid_states(thermal);
                let original_gas = gas.clone();
                let original_bodies: Vec<_> = bodies
                    .iter()
                    .map(|b| (b.specific_enthalpy_j_kg(), b.material().clone()))
                    .collect();
                let observed = ObservedPlate::default();
                let failure = thermal
                    .step_with_finite_gas_surface_convection(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[body::contact()],
                        std::slice::from_ref(&restricted),
                        &mut gas,
                        &[face(400.0)],
                        gas_controls(),
                        implicit(),
                        &observed,
                        None,
                        false,
                        gate,
                    )
                    .unwrap_err();
                assert!(
                    matches!(failure, FreeSurfaceEnthalpyError::Material(_)),
                    "must refuse the original ordered interface: {failure:?}"
                );
                assert!(
                    observed.maximum_gas_k.get() > 351.0,
                    "an actual coupled trial must precede source refusal"
                );
                assert_eq!(thermal.flow().lattice().grid().state_bits(), grid);
                assert_eq!(liquid_states(thermal), enthalpies);
                assert_eq!(gas, original_gas);
                for (b, (h, material)) in bodies.iter().zip(&original_bodies) {
                    assert_eq!(b.specific_enthalpy_j_kg(), *h);
                    assert_eq!(b.material(), material);
                }
                let cancellation = CancelGate::new_clock_free();
                let cancelling = ObservedPlate {
                    cancel: Some(&cancellation),
                    ..ObservedPlate::default()
                };
                let failure = thermal
                    .step_with_finite_gas_surface_convection(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[body::contact()],
                        std::slice::from_ref(&binding),
                        &mut gas,
                        &[face(400.0)],
                        gas_controls(),
                        implicit(),
                        &cancelling,
                        None,
                        false,
                        &cancellation,
                    )
                    .unwrap_err();
                assert!(cancellation.is_requested());
                assert_eq!(
                    failure,
                    FreeSurfaceEnthalpyError::from(
                        fs_conduction::enthalpy::DepositionError::Cancelled
                    )
                );
                assert_eq!(thermal.flow().lattice().grid().state_bits(), grid);
                assert_eq!(liquid_states(thermal), enthalpies);
                assert_eq!(gas, original_gas);
                for (b, (h, material)) in bodies.iter().zip(&original_bodies) {
                    assert_eq!(b.specific_enthalpy_j_kg(), *h);
                    assert_eq!(b.material(), material);
                }
                assert_eq!(thermal.accepted_steps(), 0);
                assert_eq!(thermal.physical_time_s(), 0.0);
            }
            let report = thermal
                .step_with_finite_gas_surface_convection(
                    &[],
                    &[],
                    controls(),
                    &mut bodies,
                    &[body::contact()],
                    std::slice::from_ref(&binding),
                    &mut gas,
                    &[face(400.0)],
                    gas_controls(),
                    implicit(),
                    &ObservedPlate::default(),
                    None,
                    false,
                    gate,
                )
                .unwrap()
                .conduction;
            (
                thermal.flow().lattice().grid().state_bits(),
                liquid_states(thermal),
                gas,
                bodies
                    .iter()
                    .map(|b| (b.specific_enthalpy_j_kg(), b.material().clone()))
                    .collect::<Vec<_>>(),
                report.surface.unwrap().heat_from_gas_j,
                report.body_environment.unwrap().heat_to_solids_j,
            )
        })
    };
    assert_eq!(run(true), run(false));
}
