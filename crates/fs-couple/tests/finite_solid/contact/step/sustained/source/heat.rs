//! Actual native impulse work and two-body thermal feedback, with an independent
//! scalar elimination. Caloric/source values are manufactured, not measurements.
use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_couple::friction_heat::directional::DirectionalHeatMaterial;
use fs_couple::friction_heat::finite_solid::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_matdb::InterfaceSystemCard;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::{ScalarAdmissibility, ScalarPropertyRequirement};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use fs_tribo::InputAuthority;
use std::collections::BTreeMap;

fn heat_source() -> (
    InterfaceSystemCard,
    SourceThermalEllipticCoulomb,
    QueryPoint,
) {
    let (base, _) = interface_source([1.0, 0.0]);
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut c = claim.clone();
        c.validity = c
            .validity
            .with_quantity("Ta", temperature, 299.0, 301.0)
            .with_quantity("Tb", temperature, 299.0, 301.0);
        if c.key.name() == "kinetic-y" {
            c.value = PropertyValue::Surface2D {
                axes: ["Ta".into(), "Tb".into()],
                axis_dims: [Temperature::DIMS; 2],
                samples: [
                    [299.0, 299.0],
                    [301.0, 299.0],
                    [299.0, 301.0],
                    [301.0, 301.0],
                ]
                .into_iter()
                .map(|p| (p, 0.4 + 0.1 * (p[0] - 300.0) - 0.05 * (p[1] - 300.0)))
                .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: Dims::NONE,
            };
        }
        claims.insert_claim(c).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        base.surface_a().clone(),
        base.surface_b().clone(),
        base.context().clone(),
        claims,
        vec![],
    )
    .unwrap();
    let query = conditions(300.0)
        .query
        .with("path", 0.25)
        .unwrap()
        .with_quantity("Ta", temperature, 300.0)
        .unwrap()
        .with_quantity("Tb", temperature, 300.0)
        .unwrap();
    let coefficient = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    let principal = [0, 1].map(|i| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            query.clone(),
            coefficient(["static-x", "static-y"][i]),
            coefficient(["kinetic-x", "kinetic-y"][i]),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "T".into(),
                nominal_pressure: None,
            },
            [
                ("Ta", FrictionQueryAxisKind::Temperature),
                ("Tb", FrictionQueryAxisKind::Temperature),
                ("speed", FrictionQueryAxisKind::SlidingSpeed),
                ("path", FrictionQueryAxisKind::AccumulatedDistance),
            ]
            .into_iter()
            .map(|(name, kind)| FrictionQueryAxis {
                name: name.into(),
                kind,
            })
            .collect(),
        )
        .unwrap()
    });
    let law = SourceThermalEllipticCoulomb::from_principal_sources(principal, [1.0, 0.0]).unwrap();
    (card, law, query)
}

fn thermal(curves: &[EquilibriumEnthalpyPhaseCurve; 2]) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for i in 0..2 {
        let length = [1.0 / 12.0, 1.0 / 6.0][i];
        let (start, end) = if i == 0 {
            (-length, 0.0)
        } else {
            (0.0, length)
        };
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64 + 1),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(i as u64),
                mass_kg: 1.0,
                specific_h_j_kg: 0.0,
                geometry: CellGeometry {
                    start_m: [0.25, 0.25, start],
                    end_m: [0.25, 0.25, end],
                    volume_m3: length,
                    surface_area_m2: 2.0 + 4.0 * length,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 7,
            a: CellId(1),
            b: CellId(2),
            area_m2: 1.0,
            conductance_w_k: 0.5,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), &curves[0] as &dyn ThermalLaw),
            (MaterialId(1), &curves[1] as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

#[test]
#[allow(clippy::too_many_lines)]
fn sustained_heat_uses_native_tangential_work_and_actual_body_temperatures_atomically() {
    with_cx(|_, cx| {
        let cards = [
            source_card(12.0, 1.0, false, false),
            source_card(6.0, 1.0, false, false),
        ];
        let material = bind(&cards[1]);
        let solid = stored_solid(&material, cx);
        let surface = FiniteSolidContactSurface::new(&solid, 64, cx).unwrap();
        let scene = FiniteSolidSphereScene::new(&[sphere(9, 0, 0.0, 0.0)], &[], 0).unwrap();
        let bodies = [body_at(
            [0.25, 0.25, -0.125 - GAP],
            Vec3::ZERO,
            Vec3::new(0.0, 1.0, 0.0),
        )];
        let history = FiniteSolidSourceTrajectoryHistory::new(&solid, &scene);
        let (interface, law, query) = heat_source();
        let curves = [0, 1].map(|i| {
            EquilibriumEnthalpyPhaseCurve::try_single_phase(
                cards[i].content_hash(),
                SolidLiquidPhase::Solid,
                [280.0, 320.0]
                    .into_iter()
                    .map(|t| EnthalpyPhaseKnot {
                        specific_enthalpy_j_kg: t - 300.0,
                        temperature_k: t,
                        liquid_mass_fraction: 0.0,
                        bulk_density_kg_m3: [12.0, 6.0][i],
                    })
                    .collect(),
            )
            .unwrap()
        });
        let thermal = thermal(&curves);
        let route = FrictionHeatRoute::from_link(&thermal, 7, [1.0 / 24.0, 1.0 / 12.0]).unwrap();
        let detected = scene
            .contacts(&surface, &solid, &bodies, source_controls(STEP).surface, cx)
            .unwrap();
        let mut contact = FiniteSolidHeatContact {
            contact_id: detected.contacts[0].id,
            key: detected.keys[0],
            binding: binding(&law),
            route: &route,
            material: DirectionalHeatMaterial {
                point: ThermalInterfacePoint {
                    interface: &interface,
                    materials: [&cards[0], &cards[1]],
                    context: &query,
                    temperature_axes: ["Ta", "Tb"],
                    liquid_fraction_axes: [None, None],
                    mean_temperature_axis: None,
                },
                caloric: [&curves[0], &curves[1]],
                interface_temperature_axis: Some("T"),
            },
            heat_fraction: 0.75,
            temperature_guess_k: 300.0,
            patch_radius_m: 0.02,
        };
        let motion = FiniteSolidHeatMotion {
            solid: &solid,
            surface: &surface,
            scene: &scene,
            bodies: &bodies,
            properties: &[homogeneous_sphere()],
            nodal_forces_n: &[],
            body_force_n_m3: [0.0; 3],
            wrenches: &[Wrench::ZERO],
            gravity: upward_gravity(),
            origin_m: ORIGIN,
        };
        let tolerance = FrictionHeatTolerances {
            cell_j: 1e-9,
            thermal_j: 1e-9,
            interface_j: 1e-9,
            total_j: 1e-9,
        };
        let controls = FiniteSolidHeatControls {
            mechanical: source_controls(STEP),
            history: Default::default(),
            maximum_cells: 2,
            maximum_interfaces: 1,
            maximum_unknowns: 3,
            maximum_trials: 256,
            maximum_newton_iterations: 20,
            difference_step_j: 1e-6,
            temperature_scale_j_per_k: 1.0,
            geometry_tolerance_m: 2e-8,
            normal_tolerance: 1e-10,
            heat: tolerance,
            combined_energy_tolerance_j: 3e-9,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 1e-11,
                linear_restart: 3,
                ..Default::default()
            },
            audit_entropy: true,
        };
        let transport = |ids: &[u64]| FiniteSolidHeatSource {
            interfaces: ids
                .iter()
                .map(|&link_id| FiniteSolidHeatConductivity {
                    link_id,
                    conductivity_w_m_k: [1.0 / 24.0, 1.0 / 12.0],
                })
                .collect(),
            transport: StepCoefficients {
                sources_w: BTreeMap::from([(CellId(1), 0.2 / STEP), (CellId(2), -0.2 / STEP)]),
                conductances_w_k: BTreeMap::new(),
                payload: (),
            },
        };
        for fraction in [0.0, 0.75] {
            contact.heat_fraction = fraction;
            let next = propose_finite_solid_contact_heat(
                &motion,
                &history,
                &thermal,
                STEP,
                &[contact],
                controls,
                &[0.0; 2],
                |_, ids, _| Ok::<_, &'static str>(transport(ids)),
                cx,
            )
            .unwrap();
            // C_A=C_B=1, G_A=G_B=1. Equal/opposite outside heat gives
            // ta=(dt*ti+0.2)/(1+dt), tb=(dt*ti-0.2)/(1+dt).
            // Native Wtt=7/4, Jn=2*a*dt; Q=f*mu*Jn*(1-7*mu*Jn/8).
            // mu=0.4+0.1*ta-0.05*tb=b0+0.025*Q.
            let jn = 2.0 * ACCELERATION * STEP;
            let b0 = 0.4 + 0.03 / (1.0 + STEP);
            let a = 0.025 * 0.875 * fraction * jn * jn;
            let b = 1.0 - 0.025 * fraction * jn;
            let mu = 2.0 * b0 / (b + (b * b + 4.0 * a * b0).sqrt());
            let work = mu * jn * (1.0 - 0.875 * mu * jn);
            let q = fraction * work;
            let ti = q * (1.0 + STEP) / (2.0 * STEP);
            let ta = (STEP * ti + 0.2) / (1.0 + STEP);
            let tb = (STEP * ti - 0.2) / (1.0 + STEP);
            assert_eq!(next.mechanical.steps.len(), 1);
            let event = next.mechanical.steps[0].event_data.as_ref().unwrap();
            near(event.sources[0].sample.coefficients_mu[1], mu, 2e-8);
            near(event.work[0].friction_heat_candidate_j.unwrap(), work, 2e-9);
            near(next.heat.dissipated_j, q, 2e-9);
            near(
                next.exported_tangential_work_j,
                (1.0 - fraction) * work,
                2e-9,
            );
            near(
                next.unconverted_contact_energy_j,
                (ACCELERATION * STEP).powi(2),
                2e-9,
            );
            near(
                next.heat.interfaces[0].interface_temperature_k,
                300.0 + ti,
                2e-7,
            );
            for (id, t) in [(CellId(1), ta), (CellId(2), tb)] {
                near(next.thermal.temperature_k(id).unwrap(), 300.0 + t, 2e-8);
            }
            for source in event.sources[0].sample.state.friction_sources() {
                for (axis, t) in [("T", ti), ("Ta", ta), ("Tb", tb)] {
                    near(source.state_coordinate(axis).unwrap(), 300.0 + t, 2e-7);
                }
            }
            near(next.combined_energy_residual_j, 0.0, 3e-9);
            let travel = (1.0 - 1.75 * mu * jn) * STEP;
            near(
                next.history
                    .entries()
                    .next()
                    .unwrap()
                    .state
                    .friction_sources()[0]
                    .state_coordinate("path")
                    .unwrap(),
                0.25 + travel,
                2e-10,
            );
            assert_eq!(next.history.time_s(), next.thermal.time_s());
        }
        // Native motion leaves a zero-radius patch; neither owner publishes it.
        let bad = FiniteSolidHeatContact {
            patch_radius_m: 0.0,
            ..contact
        };
        assert!(
            propose_finite_solid_contact_heat(
                &motion,
                &history,
                &thermal,
                STEP,
                &[bad],
                controls,
                &[0.0; 2],
                |_, ids, _| Ok::<_, &'static str>(transport(ids)),
                cx
            )
            .is_err()
        );
        let failed = propose_finite_solid_contact_heat(
            &motion,
            &history,
            &thermal,
            STEP,
            &[contact],
            controls,
            &[0.0; 2],
            |state, ids, _| {
                if state.temperature_k(CellId(1)).unwrap() > 300.1 {
                    Err("late conductivity refusal")
                } else {
                    Ok(transport(ids))
                }
            },
            cx,
        );
        assert!(matches!(
            failed,
            Err(FiniteSolidHeatError::Source("late conductivity refusal"))
        ));
        assert_eq!(history.entries().len(), 0);
        assert_eq!(solid.native_state().t, 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.temperature_k(CellId(1)).unwrap(), 300.0);
        assert_eq!(thermal.temperature_k(CellId(2)).unwrap(), 300.0);
    });
}
