//! G1/G3/G4: thermal material state reaches the actual planar contact boundary.

use super::*;
use fs_conduction::enthalpy::CellGeometry;
use fs_lbm::d3q19::freesurface3::PlanarWall3;
use fs_lbm::freesurface::{
    ADVANCING_CONTACT_ANGLE_PROPERTY, ContactLineRegime2, DYNAMIC_WETTING_SPEED_AXIS,
    DYNAMIC_WETTING_SPEED_DIMS, RECEDING_CONTACT_ANGLE_PROPERTY,
};
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceBodyContact, FreeSurfaceThermalBody,
};
use fs_lbm::material::free_surface::wetting::{FreeSurfaceWettingControls, FreeSurfaceWettingError};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::state_point::resolve_interface_state_point;
use fs_material::thermal::source::ThermalMaterialPoint;

#[path = "wetting/gas.rs"]
mod gas;

const WALL_CELL: [usize; 3] = [0, 1, 3];

fn wetting_controls() -> FreeSurfaceWettingControls {
    FreeSurfaceWettingControls { pinning_speed_m_s: 0.0, maximum_wall_cells: 64 }
}

fn interface(
    liquid: &MaterialCard, solid: &MaterialCard, slope: f64, rate_slope: f64,
    upper_temperature: f64, maximum_speed: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, sign) in [(RECEDING_CONTACT_ANGLE_PROPERTY, -1.0), (ADVANCING_CONTACT_ANGLE_PROPERTY, 1.0)] {
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, Dims::NONE),
            value: PropertyValue::Surface2D {
                axes: ["T".into(), DYNAMIC_WETTING_SPEED_AXIS.into()],
                axis_dims: [Temperature::DIMS, DYNAMIC_WETTING_SPEED_DIMS],
                samples: [[300.0, 0.0], [500.0, 0.0], [300.0, maximum_speed], [500.0, maximum_speed]]
                    .into_iter().map(|p| (p, 1.2 + slope * (p[0] - 400.0) + sign * rate_slope * p[1]))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, upper_temperature)
                .with_quantity(DYNAMIC_WETTING_SPEED_AXIS,
                    QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS), 0.0, maximum_speed),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic thermal contact-angle surface; not measured".into(),
                license: "test-only".into(), artifact: None,
            },
        };
        if slope == 0.0 && rate_slope == 0.0 {
            claim.value = PropertyValue::Scalar { value: 1.2, dims: Dims::NONE };
            claim.interpolation = InterpolationPolicy::ConstantWithinValidity;
        }
        claims.insert_claim(claim.clone()).unwrap();
        if upper_temperature < 500.0 {
            let mut fallback = claim;
            fallback.validity = ValidityDomain::unconstrained()
                .with("T", upper_temperature + 0.000001, 500.0)
                .with_quantity(DYNAMIC_WETTING_SPEED_AXIS,
                    QuantitySpec::dimensional(DYNAMIC_WETTING_SPEED_DIMS), 0.0, maximum_speed);
            claims.insert_claim(fallback).unwrap();
        }
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec { material: solid.id().clone(), texture_frame: "stationary wall finish".into() },
        SurfaceSpec { material: liquid.id().clone(), texture_frame: "liquid".into() },
        SystemContext {
            medium: "declared synthetic liquid".into(), third_body: None,
            environment: "declared synthetic gas".into(), history: "fixed clean surfaces".into(),
        }, claims, vec![],
    ).unwrap()
}

fn with_thermal<R>(
    card: &MaterialCard, interface: &InterfaceSystemCard, shear: bool, capillary: bool,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = if capillary { liquid_source(card, shear) } else {
        ResolvedLiquidSurfaceState::try_from_material(
            bundle(card, false, 400.0), LiquidViscosityBasis::Dynamic, LiquidCapillarity::Neglected,
        ).unwrap()
    };
    let law = LiquidTemperatureLaw::try_from_material(
        card, source.material().clone(), "T", LiquidViscosityBasis::Dynamic,
    ).unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        let thermal = if shear {
            scaling.free_surface3_from_cells_with_rheology(
                domain(), &cells(), &law, &chart, FreeSurfaceEnthalpyControls::default(),
                RheologySolve {
                    maximum_iterations: 128, relaxation: 0.5,
                    absolute_relaxation_tolerance: 1e-13, relative_relaxation_tolerance: 1e-13,
                }, |_| 200000.0, cx, &gate,
            )
        } else {
            scaling.free_surface3_from_cells(domain(), &cells(), cx).unwrap().with_enthalpy(
                &law, &chart, FreeSurfaceEnthalpyControls::default(), |_| 200000.0, &gate,
            )
        }.unwrap();
        let mut thermal = thermal.with_dynamic_planar_wetting(
            PlanarWall3::XMin, interface, MaterialPropertySelection::SingleClaimOnly,
            wetting_controls(), &gate,
        ).unwrap();
        run(&mut thermal, &gate)
    })
}

fn check_samples(thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>, slope: f64, rate_slope: f64) {
    let original = thermal.wetting_source().unwrap();
    for (&position, sample) in thermal.wetting_samples().unwrap() {
        let t = sample.material.query_point().iter().find(|(axis, _)| axis == "T").unwrap().1;
        if let Some(material) = thermal.material_at(position) {
            near(t, temperature(material), 2e-12);
        } else {
            assert_eq!(t, 400.0, "dry numerical ghosts retain explicit reference context");
        }
        let receding = 1.2 + slope * (t - 400.0) - rate_slope * sample.signed_speed_m_s.abs();
        let advancing = 1.2 + slope * (t - 400.0) + rate_slope * sample.signed_speed_m_s.abs();
        for (name, value) in [(RECEDING_CONTACT_ANGLE_PROPERTY, receding), (ADVANCING_CONTACT_ANGLE_PROPERTY, advancing)] {
            let property = sample.material.property(name).unwrap();
            near(property.value_si(), value, 2e-12);
            assert_eq!(property.answer().receipt.selected, original.property(name).unwrap().answer().receipt.selected);
        }
        let (regime, angle) = if sample.signed_speed_m_s > 0.0 {
            (ContactLineRegime2::Advancing, advancing)
        } else if sample.signed_speed_m_s < 0.0 {
            (ContactLineRegime2::Receding, receding)
        } else {
            (ContactLineRegime2::Pinned, sample.geometric_angle_radians.clamp(receding, advancing))
        };
        assert_eq!(sample.regime, regime);
        near(sample.applied_angle_radians, angle, 2e-12);
    }
}

fn wall<'a>(card: &'a MaterialCard, chart: &'a EquilibriumEnthalpyPhaseCurve, context: &'a QueryPoint) -> FreeSurfaceThermalBody<'a> {
    FreeSurfaceThermalBody::try_new(
        42, ThermalMaterialPoint { material: card, context, temperature_axis: "T", liquid_fraction_axis: None },
        chart, CellGeometry {
            start_m: [-0.001, 0.0015, 0.0035], end_m: [0.0, 0.0015, 0.0035],
            volume_m3: 1e-9, surface_area_m2: 6e-6,
        }, 75000.0,
        &PropertyKey::with_quantity("thermal_conductivity", QuantitySpec::semantic(
            SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))),
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

fn contact() -> FreeSurfaceBodyContact {
    FreeSurfaceBodyContact { liquid: WALL_CELL, wall: [-1, 1, 3], ..body::contact() }
}

#[test]
fn g1_constant_temperature_thermal_wetting_preserves_the_original_native_trajectory() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    let source = liquid_source(&liquid, false);
    let gate = CancelGate::new_clock_free();
    let reference = with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap()
            .free_surface3_from_cells(domain(), &cells(), cx).unwrap();
        let mut wetting = flow.with_dynamic_planar_wetting(
            PlanarWall3::XMin, &interface, MaterialPropertySelection::SingleClaimOnly,
            wetting_controls(), &gate,
        ).unwrap();
        (0..2).map(|_| {
            wetting.step_with_gate(&gate).unwrap();
            wetting.flow().lattice().grid().state_bits()
        }).collect::<Vec<_>>()
    });
    with_thermal(&liquid, &interface, false, true, |thermal, gate| {
        let energy = thermal.total_enthalpy_j().unwrap();
        for expected in reference {
            let report = thermal.step(gate).unwrap();
            assert_eq!(thermal.flow().lattice().grid().state_bits(), expected);
            near(thermal.total_enthalpy_j().unwrap(), energy, 2e-11);
            assert_eq!(report.external_heat_j, 0.0);
            assert_eq!(report.viscous_heat_j, 0.0);
            check_samples(thermal, 0.0, 0.0);
        }
    });
}

#[test]
fn g3_temperature_dependent_interface_data_changes_the_same_native_step() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let run = |slope, capillary| {
        let interface = interface(&liquid, &solid, slope, 0.002, 500.0, 100.0);
        with_thermal(&liquid, &interface, false, capillary, |thermal, gate| {
            let before = thermal.total_enthalpy_j().unwrap();
            let heat = cells().iter().map(|cell| FreeSurfaceCellHeat {
                coordinates: cell.coordinates, energy_j: 0.04,
            }).collect::<Vec<_>>();
            let supplied = heat.iter().map(|source| source.energy_j).sum::<f64>();
            let report = thermal.step_with_heat(&heat, gate).unwrap();
            assert!((thermal.total_enthalpy_j().unwrap() - before - supplied).abs() < 1e-9);
            near(report.external_heat_j, supplied, 2e-13);
            assert_eq!(report.viscous_heat_j, 0.0);
            near(temperature(thermal.material_at(WALL_CELL).unwrap()), 420.0, 2e-11);
            check_samples(thermal, slope, 0.002);
            thermal.flow().lattice().grid().state_bits()
        })
    };
    let fixed = run(0.0, true);
    let responsive = run(0.004, true);
    assert_eq!(fixed.len(), responsive.len());
    let delta = fixed.iter().zip(&responsive)
        .map(|(&a, &b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
        .fold(0.0_f64, f64::max);
    assert!(delta > 1e-12, "heated interface data must reach native momentum immediately");
    assert_eq!(run(0.0, false), run(0.004, false), "omitted capillary traction keeps angle changes inert");
}

#[test]
fn g1_g3_finite_wall_and_gas_heat_drive_local_wetting_and_rheological_motion() {
    let liquid = source_card(0.004, true, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let chart = body::wall_chart(&solid);
    let context = point(450.0);
    let run = |slope| {
        let interface = interface(&liquid, &solid, slope, 0.002, 500.0, 100.0);
        with_thermal(&liquid, &interface, true, true, |thermal, gate| {
            let mut gas = gas(100.0);
            let mut bodies = [wall(&solid, &chart, &context)];
            let before = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j()
                + gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
            let report = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[contact()], &[], &mut gas,
                &[gas_face()], gas_controls(), false, gate,
            ).unwrap().conduction;
            // Independent two-layer explicit conduction with the authored
            // liquid/wall heat capacities, not a production thermal helper.
            let (mut tl, mut tw, mut qwall) = (400.0, 450.0, 0.0);
            for _ in 0..4 {
                let g = 1e-6 / (0.0005 / (tl - 200.0) + 0.0005 / (0.5 * (tw - 200.0)));
                let q = 0.00025 * g * (tw - tl);
                tl += q / 0.002;
                tw -= q / 0.001;
                qwall += q;
            }
            let wall_report = report.bodies.as_ref().unwrap();
            near(temperature(&wall_report.thermal_endpoint[0].liquid), tl, 2e-11);
            near(wall_report.heat_to_liquid_j, qwall, 2e-11);
            near(bodies[0].specific_enthalpy_j_kg(), (tw - 300.0) * 500.0, 2e-11);
            assert_eq!(report.thermal.substeps, 4);
            assert_eq!(report.flow.viscous_heat_j, 0.0);
            let after = thermal.total_enthalpy_j().unwrap() + bodies[0].enthalpy_j()
                + gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
            assert!((after - before).abs() < 1e-9);
            assert!(report.finite_gas.as_ref().unwrap().combined_thermal_balance_defect_j.abs() < 1e-9);
            check_samples(thermal, slope, 0.002);
            let state = thermal.material_at(WALL_CELL).unwrap();
            let rate = thermal.rheology_at(WALL_CELL).unwrap().shear_rate_per_s;
            near(state.dynamic_viscosity_pa_s(), 0.1 - 0.0002 * (temperature(state) - 400.0) + 0.002 * rate, 2e-12);
            (thermal.flow().lattice().grid().state_bits(), wall_report.heat_to_liquid_j,
                report.environment.as_ref().unwrap().heat_to_liquid_j)
        })
    };
    let fixed = run(0.0);
    let responsive = run(0.004);
    assert_ne!(fixed.0, responsive.0);
    assert_eq!((fixed.1, fixed.2), (responsive.1, responsive.2), "wetting must not manufacture surface heat");
}

#[test]
fn g4_late_temperature_and_endpoint_speed_refusals_preserve_every_participant() {
    for (upper_t, max_speed) in [(401.2, 100.0), (500.0, 1e-6)] {
        let liquid = source_card(0.004, false, 500.0);
        let solid = body::wall_card(0.5, 500.0);
        let chart = body::wall_chart(&solid);
        let context = point(450.0);
        let interface = interface(&liquid, &solid, 0.004, 0.002, upper_t, max_speed);
        with_thermal(&liquid, &interface, false, true, |thermal, gate| {
            let mut gas = gas(100.0);
            let mut bodies = [wall(&solid, &chart, &context)];
            let bits = thermal.flow().lattice().grid().state_bits();
            let energy = thermal.total_enthalpy_j().unwrap();
            let mass = thermal.flow().ledger_mass_kg().unwrap();
            let samples = thermal.wetting_samples().unwrap().clone();
            let material = thermal.material_at(WALL_CELL).unwrap().clone();
            let body_h = bodies[0].specific_enthalpy_j_kg();
            let body_source = bodies[0].material().clone();
            let gas_state = gas.state().clone();
            let gas_point = gas.point().clone();
            if upper_t < 500.0 {
                let source = thermal.wetting_source().unwrap();
                let requirements = source.properties().iter().map(|p| p.requirement().clone()).collect::<Vec<_>>();
                let fresh = resolve_interface_state_point(
                    &interface, &source.query().with("T", 450.0).unwrap(), &requirements,
                    MaterialPropertySelection::SingleClaimOnly,
                ).unwrap();
                assert_ne!(fresh.property(RECEDING_CONTACT_ANGLE_PROPERTY).unwrap().answer().receipt.selected,
                    source.property(RECEDING_CONTACT_ANGLE_PROPERTY).unwrap().answer().receipt.selected);
            }
            // Temperature crosses 401.2 K on the SECOND thermal substep. The
            // speed case instead passes the whole thermal stage and refuses
            // the contact history produced by the candidate native motion.
            let error = thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[contact()], &[], &mut gas,
                &[gas_face()], gas_controls(), false, gate,
            ).unwrap_err();
            assert!(matches!(error, FreeSurfaceEnthalpyError::Wetting(FreeSurfaceWettingError::Material(_))));
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
            assert_eq!(thermal.flow().ledger_mass_kg().unwrap(), mass);
            assert_eq!(thermal.wetting_samples().unwrap(), &samples);
            assert_eq!(thermal.material_at(WALL_CELL).unwrap(), &material);
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), body_h);
            assert_eq!(bodies[0].material(), &body_source);
            assert_eq!(gas.state(), &gas_state);
            assert_eq!(gas.point(), &gas_point);
            if upper_t < 500.0 {
                thermal.step_with_finite_gas_environment(
                    &[], &[], controls(), &mut bodies, &[], &[], &mut gas,
                    &[gas_face()], gas_controls(), false, gate,
                ).unwrap();
            }
            let accepted = thermal.flow().lattice().grid().state_bits();
            let accepted_samples = thermal.wetting_samples().unwrap().clone();
            let accepted_gas = gas.state().clone();
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(thermal.step_with_finite_gas_environment(
                &[], &[], controls(), &mut bodies, &[], &[], &mut gas,
                &[gas_face()], gas_controls(), false, &cancelled,
            ).is_err());
            assert_eq!(thermal.flow().lattice().grid().state_bits(), accepted);
            assert_eq!(thermal.wetting_samples().unwrap(), &accepted_samples);
            assert_eq!(gas.state(), &accepted_gas);
            assert_eq!(bodies[0].specific_enthalpy_j_kg(), body_h);
        });
    }
}

#[test]
fn g1_cox_contact_angle_uses_local_heated_viscosity_without_inventing_heat() {
    check_heated_cox_viscosity(0.0);
}

#[test]
fn g1_two_viscous_wedges_retain_explicit_ratio_along_local_thermal_motion() {
    check_heated_cox_viscosity(1.0);
}

fn check_heated_cox_viscosity(outer_to_liquid_ratio: f64) {
    use fs_lbm::material::free_surface::wetting::{CoxVoinovControls, CoxVoinovLaw};
    let base = source_card(0.004, false, 500.0);
    let mut claims = ClaimSet::new();
    for (_, original) in base.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == "dynamic_viscosity" {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(300.0, 0.12), (500.0, 0.08)],
                dims: claim.key.dims(),
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    let liquid = MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap();
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    let source = liquid_source(&liquid, false);
    let law = LiquidTemperatureLaw::try_from_material(
        &liquid,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(&liquid);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let initial = cells()
            .into_iter()
            .map(|mut cell| {
                cell.velocity_m_s = [0.0; 3];
                cell
            })
            .collect::<Vec<_>>();
        let flow = plan_free_surface_from_material(&source, 0.001, 0.00002, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial, cx)
            .unwrap();
        let mut thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap()
            .with_dynamic_planar_wetting(
                PlanarWall3::XMin,
                &interface,
                MaterialPropertySelection::SingleClaimOnly,
                wetting_controls(),
                &gate,
            )
            .unwrap()
            .with_cox_voinov_wetting(
                CoxVoinovLaw::try_new_with_viscosity_ratio(
                    CoxVoinovControls {
                        microscopic_cutoff_m: 1e-9,
                        observation_length_m: 0.001,
                        maximum_capillary_number: 0.01,
                        minimum_angle_radians: 0.01,
                        maximum_angle_radians: 3.0,
                    },
                    outer_to_liquid_ratio,
                    &gate,
                )
                .unwrap(),
                &gate,
            )
            .unwrap();
        let before = thermal.total_enthalpy_j().unwrap();
        let heat = initial
            .iter()
            .map(|cell| FreeSurfaceCellHeat {
                coordinates: cell.coordinates,
                energy_j: 0.04,
            })
            .collect::<Vec<_>>();
        let supplied = heat.iter().map(|entry| entry.energy_j).sum::<f64>();
        let report = thermal.step_with_heat(&heat, &gate).unwrap();
        near(thermal.total_enthalpy_j().unwrap() - before, supplied, 2e-9);
        near(report.external_heat_j, supplied, 2e-12);
        assert_eq!(report.viscous_heat_j, 0.0);
        let mut heated_moving_samples = 0;
        for (&position, sample) in thermal.wetting_samples().unwrap() {
            let prediction = sample.hydrodynamic.unwrap();
            assert_eq!(
                prediction.outer_to_liquid_viscosity_ratio,
                outer_to_liquid_ratio
            );
            assert_eq!(
                sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
                Some(0.0)
            );
            if let Some(state) = thermal.material_at(position) {
                let mu = state.dynamic_viscosity_pa_s();
                near(prediction.dynamic_viscosity_pa_s, mu, 1e-14);
                near(mu, 0.1 - 0.0002 * (temperature(state) - 400.0), 2e-12);
                near(
                    prediction.signed_capillary_number,
                    mu * sample.signed_speed_m_s / 0.004,
                    1e-14,
                );
                if sample.signed_speed_m_s != 0.0 && (mu - 0.1).abs() > 1e-4 {
                    heated_moving_samples += 1;
                }
            }
            if sample.regime != ContactLineRegime2::Pinned {
                let expected = thermal
                    .cox_voinov_wetting_law()
                    .unwrap()
                    .apparent_angle_radians(1.2, prediction.signed_capillary_number)
                    .unwrap();
                near(sample.applied_angle_radians, expected, 1e-13);
            }
        }
        assert!(heated_moving_samples > 0);
    });
}
