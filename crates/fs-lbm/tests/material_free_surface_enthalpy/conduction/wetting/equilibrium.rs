//! Thermal adhesion uses the existing accepted heat and native flow history.

use super::*;
use fs_lbm::material::free_surface::wetting::EquilibriumAngleOrigin;
use fs_material::liquid::surface::EQUILIBRIUM_CONTACT_ANGLE_PROPERTY;
use fs_material::liquid::surface::adhesion::EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY;

#[path = "equilibrium_height_graph.rs"]
mod height_graph;

fn equilibrium_interface(
    liquid: &MaterialCard,
    solid: &MaterialCard,
    origin: EquilibriumAngleOrigin,
    slope: f64,
    upper_temperature: f64,
) -> InterfaceSystemCard {
    let original = interface(liquid, solid, 0.0, 0.0, 500.0, 100.0);
    let (name, dims, reference) = match origin {
        EquilibriumAngleOrigin::SourcedAngle => (
            EQUILIBRIUM_CONTACT_ANGLE_PROPERTY,
            Dims::NONE,
            core::f64::consts::FRAC_PI_2,
        ),
        EquilibriumAngleOrigin::YoungDupre => (
            EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY,
            SURFACE_TENSION_DIMS,
            0.004,
        ),
    };
    let value = if slope == 0.0 {
        PropertyValue::Scalar {
            value: reference,
            dims,
        }
    } else {
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![
                (300.0, (reference - 100.0 * slope).max(0.25 * reference)),
                (400.0, reference),
                (500.0, reference + 100.0 * slope),
            ],
            dims,
        }
    };
    let claim = PropertyClaim {
        key: PropertyKey::new(name, dims),
        value,
        validity: ValidityDomain::unconstrained().with("T", 300.0, upper_temperature),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: if slope == 0.0 {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        observations: vec![],
        provenance: Provenance {
            source: "synthetic reversible wetting thermodynamics; not measured".into(),
            license: "test-only".into(),
            artifact: None,
        },
    };
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim.clone()).unwrap();
    if upper_temperature < 500.0 {
        // A different supported source exists beyond the first source's range.
        // A thermal trajectory must refuse, rather than silently switching.
        let mut fallback = claim;
        fallback.validity =
            ValidityDomain::unconstrained().with("T", upper_temperature + 0.000001, 500.0);
        claims.insert_claim(fallback).unwrap();
    }
    InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        vec![],
    )
    .unwrap()
}

fn with_equilibrium<R>(
    card: &MaterialCard,
    interface: &InterfaceSystemCard,
    origin: EquilibriumAngleOrigin,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(card, false);
    let law = LiquidTemperatureLaw::try_from_material(
        card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &cells(), cx)
            .unwrap();
        let thermal = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        let mut thermal = thermal
            .with_equilibrium_planar_wetting(
                PlanarWall3::XMin,
                interface,
                origin,
                MaterialPropertySelection::SingleClaimOnly,
                64,
                &gate,
            )
            .unwrap();
        run(&mut thermal, &gate)
    })
}

fn check_equilibrium(thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>, slope: f64) {
    let original = thermal.wetting_source().unwrap();
    let selected = original
        .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
        .unwrap()
        .answer()
        .receipt
        .selected;
    for (&position, sample) in thermal.wetting_samples().unwrap() {
        let t = sample
            .material
            .query_point()
            .iter()
            .find(|(axis, _)| axis == "T")
            .unwrap()
            .1;
        if let Some(material) = thermal.material_at(position) {
            near(t, temperature(material), 2e-12);
        } else {
            assert_eq!(
                t, 400.0,
                "dry wall ghosts retain the declared reference context"
            );
        }
        let expected_work = 0.004 + slope * (t - 400.0);
        let work = sample
            .material
            .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
            .unwrap();
        near(work.value_si(), expected_work, 2e-12);
        assert_eq!(work.answer().receipt.selected, selected);
        near(
            sample.applied_angle_radians,
            (expected_work / 0.004 - 1.0).acos(),
            2e-12,
        );
        assert_eq!(sample.equilibrium, Some(EquilibriumAngleOrigin::YoungDupre));
        assert_eq!(sample.regime, ContactLineRegime2::Equilibrium);
        assert_eq!(sample.signed_speed_m_s, 0.0);
        assert!(sample.hydrodynamic.is_none());
        assert!(
            sample
                .material
                .property(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY)
                .is_none()
        );
        assert!(
            sample
                .material
                .property(ADVANCING_CONTACT_ANGLE_PROPERTY)
                .is_none()
        );
        assert!(
            !sample
                .material
                .query_point()
                .iter()
                .any(|(axis, _)| axis == DYNAMIC_WETTING_SPEED_AXIS)
        );
    }
}

#[test]
fn g1_constant_adhesion_and_measured_angle_keep_the_same_native_thermal_trajectory() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let run = |origin| {
        let interface = equilibrium_interface(&liquid, &solid, origin, 0.0, 500.0);
        with_equilibrium(&liquid, &interface, origin, |thermal, gate| {
            let original = thermal.wetting_source().unwrap().clone();
            let initial = thermal.wetting_samples().unwrap().clone();
            for sample in initial.values() {
                assert_eq!(sample.material, original);
                assert_eq!(sample.applied_angle_radians, core::f64::consts::FRAC_PI_2);
            }
            let energy = thermal.total_enthalpy_j().unwrap();
            let mut result = Vec::new();
            for _ in 0..3 {
                let report = thermal.step(gate).unwrap();
                near(thermal.total_enthalpy_j().unwrap(), energy, 2e-11);
                assert_eq!(report.external_heat_j, 0.0);
                for sample in thermal.wetting_samples().unwrap().values() {
                    assert_eq!(sample.material, original);
                }
                result.push(thermal.flow().lattice().grid().state_bits());
            }
            result
        })
    };
    assert_eq!(
        run(EquilibriumAngleOrigin::YoungDupre),
        run(EquilibriumAngleOrigin::SourcedAngle)
    );
}

#[test]
fn g1_g3_local_thermal_adhesion_drives_native_momentum_and_retains_original_sources() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let run = |slope| {
        let interface = equilibrium_interface(
            &liquid,
            &solid,
            EquilibriumAngleOrigin::YoungDupre,
            slope,
            500.0,
        );
        with_equilibrium(
            &liquid,
            &interface,
            EquilibriumAngleOrigin::YoungDupre,
            |thermal, gate| {
                let initial_energy = thermal.total_enthalpy_j().unwrap();
                let heat = [FreeSurfaceCellHeat {
                    coordinates: WALL_CELL,
                    energy_j: 0.04,
                }];
                let report = thermal.step_with_heat(&heat, gate).unwrap();
                near(
                    thermal.total_enthalpy_j().unwrap() - initial_energy,
                    0.04,
                    2e-11,
                );
                near(report.external_heat_j, 0.04, 2e-13);
                near(report.mass_before_kg, report.mass_after_kg, 2e-13);
                assert_eq!(report.viscous_heat_j, 0.0);
                check_equilibrium(thermal, slope);
                let boundary = &thermal.wetting_samples().unwrap()[&WALL_CELL];
                assert!(
                    boundary
                        .material
                        .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                        .unwrap()
                        .value_si()
                        >= 0.004
                );
                thermal.flow().lattice().grid().state_bits()
            },
        )
    };
    let fixed = run(0.0);
    let responsive = run(0.00002);
    assert_eq!(fixed.len(), responsive.len());
    let change = fixed
        .iter()
        .zip(&responsive)
        .map(|(&a, &b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
        .fold(0.0_f64, f64::max);
    assert!(
        change > 1e-12,
        "adhesion at the heated wall must reach the same native momentum step"
    );
    assert_eq!(
        responsive,
        run(0.00002),
        "accepted thermal adhesion is deterministic"
    );
}

#[test]
fn g3_g4_work_support_complete_wetting_and_tension_support_refuse_atomically() {
    let solid = body::wall_card(0.5, 500.0);
    for (work_limit, slope, sigma_limit) in [
        (410.0, 0.00002, 500.0),
        (500.0, 0.0003, 500.0),
        (500.0, 0.00002, 410.0),
    ] {
        let liquid = source_card(0.004, false, sigma_limit);
        let interface = equilibrium_interface(
            &liquid,
            &solid,
            EquilibriumAngleOrigin::YoungDupre,
            slope,
            work_limit,
        );
        with_equilibrium(
            &liquid,
            &interface,
            EquilibriumAngleOrigin::YoungDupre,
            |thermal, gate| {
                let populations = thermal.flow().lattice().grid().state_bits();
                let samples = thermal.wetting_samples().unwrap().clone();
                let before_h = cells()
                    .iter()
                    .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                    .collect::<Vec<_>>();
                let energy = thermal.total_enthalpy_j().unwrap();
                let heat = [FreeSurfaceCellHeat {
                    coordinates: WALL_CELL,
                    energy_j: 0.04,
                }];
                assert!(thermal.step_with_heat(&heat, gate).is_err());
                assert_eq!(thermal.flow().lattice().grid().state_bits(), populations);
                assert_eq!(thermal.wetting_samples().unwrap(), &samples);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
                assert_eq!(
                    cells()
                        .iter()
                        .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                        .collect::<Vec<_>>(),
                    before_h
                );
                thermal.step(gate).unwrap();
                check_equilibrium(thermal, slope);
            },
        );
    }
}

#[test]
fn g3_equilibrium_adhesion_cannot_silently_acquire_cox_dynamics() {
    use fs_lbm::material::free_surface::wetting::{CoxVoinovControls, CoxVoinovLaw};
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let interface = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.00002,
        500.0,
    );
    with_equilibrium(
        &liquid,
        &interface,
        EquilibriumAngleOrigin::YoungDupre,
        |thermal, gate| {
            let law = CoxVoinovLaw::try_new(
                CoxVoinovControls {
                    microscopic_cutoff_m: 1e-8,
                    observation_length_m: 0.001,
                    maximum_capillary_number: 0.01,
                    minimum_angle_radians: 0.1,
                    maximum_angle_radians: 3.0,
                },
                gate,
            )
            .unwrap();
            assert!(thermal.clone().with_cox_voinov_wetting(law, gate).is_err());
            thermal.step(gate).unwrap();
            check_equilibrium(thermal, 0.00002);
        },
    );
}

fn varying_tension(slope: f64, upper: f64) -> MaterialCard {
    let base = source_card(0.004, false, upper);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == SURFACE_TENSION_PROPERTY && slope != 0.0 {
            claim.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (300.0, 0.004 - 100.0 * slope),
                    (400.0, 0.004),
                    (500.0, 0.004 + 100.0 * slope),
                ],
                dims: SURFACE_TENSION_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        claims.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn uniform_cells() -> Vec<MaterialFreeSurfaceInitialCell3> {
    (0..4)
        .flat_map(|z| {
            (0..8).flat_map(move |y| {
                (0..8).map(move |x| MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: if z == 3 { 0.5 } else { 1.0 },
                    velocity_m_s: [0.0; 3],
                })
            })
        })
        .collect()
}

fn with_uniform<R>(
    card: &MaterialCard,
    interface: &InterfaceSystemCard,
    initial: &[MaterialFreeSurfaceInitialCell3],
    initial_h: f64,
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(card, false);
    let law = LiquidTemperatureLaw::try_from_material(
        card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut thermal = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells_with_uniform_thermal_capillarity(
                domain(),
                initial,
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                None,
                None,
                |_| initial_h,
                cx,
                &gate,
            )
            .unwrap()
            .with_equilibrium_planar_wetting(
                PlanarWall3::XMin,
                interface,
                EquilibriumAngleOrigin::YoungDupre,
                MaterialPropertySelection::SingleClaimOnly,
                64,
                &gate,
            )
            .unwrap();
        assert!(thermal.uniform_thermal_capillarity());
        run(&mut thermal, &gate)
    })
}

fn homogeneous_heat(thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> Vec<FreeSurfaceCellHeat> {
    uniform_cells()
        .iter()
        .map(|cell| {
            let p = cell.coordinates;
            let mass = thermal
                .flow()
                .region_mass_kg(|x, y, z| [x, y, z] == p)
                .unwrap();
            // A power-of-two specific increment preserves identical inputs exactly.
            FreeSurfaceCellHeat {
                coordinates: p,
                energy_j: mass * 65536.0,
            }
        })
        .collect()
}

#[test]
fn g1_uniform_thermal_work_and_tension_share_native_pressure_and_young_dupre() {
    let liquid = varying_tension(0.00002, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let interface = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.00002,
        500.0,
    );
    let run = |heat| {
        with_uniform(
            &liquid,
            &interface,
            &uniform_cells(),
            200000.0,
            |thermal, gate| {
                let original = thermal
                    .material_at([0, 2, 3])
                    .unwrap()
                    .material()
                    .property(SURFACE_TENSION_PROPERTY)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected;
                let supplied = if heat {
                    homogeneous_heat(thermal)
                } else {
                    vec![]
                };
                let report = thermal.step_with_heat(&supplied, gate).unwrap();
                near(
                    report.external_heat_j,
                    supplied.iter().map(|q| q.energy_j).sum(),
                    2e-12,
                );
                near(report.mass_before_kg, report.mass_after_kg, 2e-12);
                let state = thermal.material_at([0, 2, 3]).unwrap();
                let sigma = state.material().property(SURFACE_TENSION_PROPERTY).unwrap();
                assert_eq!(sigma.answer().receipt.selected, original);
                near(
                    sigma.value_si(),
                    0.004 + 0.00002 * (temperature(state) - 400.0),
                    2e-14,
                );
                if heat {
                    assert!(sigma.value_si() > 0.004);
                }
                // rho=1000, dx=dt=0.001: sigma_native=sigma/(rho*(dx/dt)^2*dx)=sigma.
                near(
                    thermal
                        .flow()
                        .lattice()
                        .surface_tension_at(0, 2, 3)
                        .unwrap(),
                    sigma.value_si(),
                    2e-14,
                );
                let sample = &thermal.wetting_samples().unwrap()[&[0, 2, 3]];
                near(
                    sample.applied_angle_radians,
                    core::f64::consts::FRAC_PI_2,
                    2e-14,
                );
                near(
                    sample
                        .material
                        .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                        .unwrap()
                        .value_si(),
                    sigma.value_si(),
                    2e-14,
                );
                assert_eq!(
                    sample.material.state_coordinate("T"),
                    Some(temperature(state))
                );
                thermal.flow().lattice().grid().state_bits()
            },
        )
    };
    assert_eq!(
        run(true),
        run(false),
        "uniform sigma with W=sigma cannot invent tangential stress or a wall angle change"
    );

    // A separately prepared hot state exercises a nonconstant W/sigma ratio.
    let interface = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.00001,
        500.0,
    );
    with_uniform(
        &liquid,
        &interface,
        &uniform_cells(),
        240000.0,
        |thermal, gate| {
            let sample = &thermal.wetting_samples().unwrap()[&[0, 2, 3]];
            near(
                sample.applied_angle_radians,
                (0.0042_f64 / 0.0044 - 1.0).acos(),
                2e-13,
            );
            near(
                thermal
                    .flow()
                    .lattice()
                    .surface_tension_at(0, 2, 3)
                    .unwrap(),
                0.0044,
                2e-14,
            );
            assert_eq!(sample.material.state_coordinate("T"), Some(420.0));
            // The admitted live coefficient and angle drive the actual native
            // pressure/collision owner, independently of later thermal transport.
            let mut native = thermal.flow().clone();
            let before = native.lattice().grid().state_bits();
            native.step_with_gate(gate).unwrap();
            assert!(
                before
                    .iter()
                    .zip(native.lattice().grid().state_bits())
                    .any(|(&a, b)| { (f64::from_bits(a) - f64::from_bits(b)).abs() > 1e-12 })
            );
        },
    );
}

#[test]
fn g3_uniform_mode_preserves_constant_tension_native_trajectory_and_receipts() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let interface = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.0,
        500.0,
    );
    let run = |thermal: &mut MaterialFreeSurfaceEnthalpy3<'_, '_>, gate: &CancelGate| {
        let original = thermal.wetting_source().unwrap().clone();
        for sample in thermal.wetting_samples().unwrap().values() {
            assert_eq!(sample.material, original);
        }
        let mut history = vec![];
        for _ in 0..3 {
            thermal.step(gate).unwrap();
            for sample in thermal.wetting_samples().unwrap().values() {
                assert_eq!(
                    sample
                        .material
                        .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                        .unwrap()
                        .answer()
                        .receipt
                        .selected,
                    original
                        .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                        .unwrap()
                        .answer()
                        .receipt
                        .selected,
                );
                if sample.material.state_coordinate("T") == original.state_coordinate("T") {
                    assert_eq!(sample.material, original);
                }
            }
            history.push(thermal.flow().lattice().grid().state_bits());
        }
        history
    };
    let ordinary = with_equilibrium(&liquid, &interface, EquilibriumAngleOrigin::YoungDupre, run);
    let uniform = with_uniform(&liquid, &interface, &cells(), 200000.0, run);
    assert_eq!(uniform, ordinary);
}

#[test]
fn g3_g4_uniform_adhesion_rejects_spatial_gradients_source_exit_and_partial_wetting_loss() {
    let solid = body::wall_card(0.5, 500.0);
    for (work_slope, upper) in [(0.00002, 500.0), (0.00002, 410.0), (0.0003, 500.0)] {
        let liquid = varying_tension(0.00002, upper);
        let interface = equilibrium_interface(
            &liquid,
            &solid,
            EquilibriumAngleOrigin::YoungDupre,
            work_slope,
            500.0,
        );
        with_uniform(
            &liquid,
            &interface,
            &uniform_cells(),
            200000.0,
            |thermal, gate| {
                let before = (
                    thermal.flow().lattice().grid().state_bits(),
                    thermal.wetting_samples().unwrap().clone(),
                    thermal.total_enthalpy_j().unwrap(),
                    thermal.physical_time_s(),
                    thermal.flow().lattice().surface_tension_at(0, 2, 3),
                );
                let assert_unchanged = |thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>| {
                    assert_eq!(thermal.flow().lattice().grid().state_bits(), before.0);
                    assert_eq!(thermal.wetting_samples().unwrap(), &before.1);
                    assert_eq!(thermal.total_enthalpy_j().unwrap(), before.2);
                    assert_eq!(thermal.physical_time_s(), before.3);
                    assert_eq!(
                        thermal.flow().lattice().surface_tension_at(0, 2, 3),
                        before.4
                    );
                };
                let local = [FreeSurfaceCellHeat {
                    coordinates: [0, 2, 3],
                    energy_j: 0.02,
                }];
                assert!(thermal.step_with_heat(&local, gate).is_err());
                assert_unchanged(thermal);
                // Rejection also occurs inside the existing thermal substep owner,
                // before a final momentum candidate or partial heat can publish.
                assert!(
                    thermal
                        .step_with_conduction(&[], &local, controls(), gate)
                        .is_err()
                );
                assert_unchanged(thermal);
                if upper < 500.0 || work_slope > 0.00002 {
                    let heat = homogeneous_heat(thermal);
                    assert!(thermal.step_with_heat(&heat, gate).is_err());
                    assert_unchanged(thermal);
                }
                let cancelled = CancelGate::new_clock_free();
                cancelled.request();
                assert!(thermal.step(&cancelled).is_err());
                assert_unchanged(thermal);
                thermal.step(gate).unwrap();
            },
        );
    }
}

// The same native uniform-tension owner must also compose with a MOVING
// contact line. Keep the original source, caloric and geometry fixtures.
fn with_dynamic_capillarity<R>(
    card: &MaterialCard,
    interface: &InterfaceSystemCard,
    uniform: bool,
    initial: &[MaterialFreeSurfaceInitialCell3],
    run: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    use fs_lbm::material::free_surface::wetting::{CoxVoinovControls, CoxVoinovLaw};

    let source = liquid_source(card, false);
    let law = LiquidTemperatureLaw::try_from_material(
        card,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&source, 0.001, 0.00002, 0.1).unwrap();
        let thermal = if uniform {
            scaling.free_surface3_from_cells_with_uniform_thermal_capillarity(
                domain(),
                initial,
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                None,
                None,
                |_| 200000.0,
                cx,
                &gate,
            )
        } else {
            scaling
                .free_surface3_from_cells(domain(), initial, cx)
                .unwrap()
                .with_enthalpy(
                    &law,
                    &chart,
                    FreeSurfaceEnthalpyControls::default(),
                    |_| 200000.0,
                    &gate,
                )
        }
        .unwrap();
        let mut thermal = thermal
            .with_dynamic_planar_wetting(
                PlanarWall3::XMin,
                interface,
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
                    1.0,
                    &gate,
                )
                .unwrap(),
                &gate,
            )
            .unwrap();
        assert_eq!(thermal.uniform_thermal_capillarity(), uniform);
        run(&mut thermal, &gate)
    })
}

#[test]
fn g3_uniform_capillarity_preserves_constant_sigma_dynamic_wetting() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    let initial = cells()
        .into_iter()
        .map(|mut cell| {
            cell.velocity_m_s = [0.0; 3];
            cell
        })
        .collect::<Vec<_>>();
    let run = |uniform| {
        with_dynamic_capillarity(&liquid, &interface, uniform, &initial, |thermal, gate| {
            let mut history = Vec::new();
            for _ in 0..3 {
                thermal.step(gate).unwrap();
                history.push((
                    thermal.flow().lattice().grid().state_bits(),
                    thermal
                        .wetting_samples()
                        .unwrap()
                        .values()
                        .map(|sample| {
                            (
                                sample.signed_speed_m_s,
                                sample.applied_angle_radians,
                                sample.hydrodynamic,
                            )
                        })
                        .collect::<Vec<_>>(),
                ));
            }
            history
        })
    };
    assert_eq!(run(false), run(true));
}

#[test]
fn g1_uniform_heating_refreshes_native_pressure_and_cox_capillary_number_together() {
    let slope = -0.00002;
    let liquid = varying_tension(slope, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let interface = interface(&liquid, &solid, 0.0, 0.0, 500.0, 100.0);
    with_dynamic_capillarity(
        &liquid,
        &interface,
        true,
        &uniform_cells(),
        |thermal, gate| {
            let selected = thermal
                .material_at([0, 2, 3])
                .unwrap()
                .material()
                .property(SURFACE_TENSION_PROPERTY)
                .unwrap()
                .answer()
                .receipt
                .selected;
            let before = thermal.total_enthalpy_j().unwrap();
            let heat = homogeneous_heat(thermal);
            let supplied = heat.iter().map(|q| q.energy_j).sum::<f64>();
            let report = thermal.step_with_heat(&heat, gate).unwrap();
            near(report.external_heat_j, supplied, 2e-12);
            near(thermal.total_enthalpy_j().unwrap() - before, supplied, 2e-9);
            near(report.mass_before_kg, report.mass_after_kg, 2e-12);
            assert_eq!(report.viscous_heat_j, 0.0);
            let mut heated_moving = 0;
            for (&position, sample) in thermal.wetting_samples().unwrap() {
                let prediction = sample.hydrodynamic.unwrap();
                assert_eq!(
                    sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
                    Some(0.0)
                );
                if let Some(state) = thermal.material_at(position) {
                    let sigma = state.material().property(SURFACE_TENSION_PROPERTY).unwrap();
                    assert_eq!(sigma.answer().receipt.selected, selected);
                    near(
                        sigma.value_si(),
                        0.004 + slope * (temperature(state) - 400.0),
                        2e-14,
                    );
                    near(prediction.surface_tension_n_m, sigma.value_si(), 2e-14);
                    near(
                        prediction.dynamic_viscosity_pa_s,
                        state.dynamic_viscosity_pa_s(),
                        2e-14,
                    );
                    near(
                        prediction.signed_capillary_number,
                        state.dynamic_viscosity_pa_s() * sample.signed_speed_m_s / sigma.value_si(),
                        2e-14,
                    );
                    let [x, y, z] = position;
                    near(
                        thermal
                            .flow()
                            .lattice()
                            .surface_tension_at(x, y, z)
                            .unwrap(),
                        sigma.value_si() * 0.00002_f64.powi(2) / (1000.0 * 0.001_f64.powi(3)),
                        2e-14,
                    );
                    if sample.signed_speed_m_s != 0.0 && sigma.value_si() < 0.0039 {
                        heated_moving += 1;
                    }
                }
                if sample.regime != ContactLineRegime2::Pinned {
                    let expected = thermal
                        .cox_voinov_wetting_law()
                        .unwrap()
                        .apparent_angle_radians(1.2, prediction.signed_capillary_number)
                        .unwrap();
                    near(sample.applied_angle_radians, expected, 2e-13);
                }
            }
            assert!(
                heated_moving > 0,
                "heated source tension must reach a moving contact line"
            );
        },
    );
}

#[test]
fn g3_g4_uniform_dynamic_wetting_refuses_sigma_gradients_and_source_exits_atomically() {
    let solid = body::wall_card(0.5, 500.0);
    for (sigma_upper, angle_upper) in [(500.0, 500.0), (410.0, 500.0), (500.0, 410.0)] {
        let liquid = varying_tension(-0.00002, sigma_upper);
        let interface = interface(&liquid, &solid, 0.0, 0.0, angle_upper, 100.0);
        with_dynamic_capillarity(
            &liquid,
            &interface,
            true,
            &uniform_cells(),
            |thermal, gate| {
                let snapshot = |thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>| {
                    (
                        thermal.flow().lattice().grid().state_bits(),
                        thermal.wetting_samples().unwrap().clone(),
                        thermal.total_enthalpy_j().unwrap(),
                        thermal.physical_time_s(),
                        uniform_cells()
                            .iter()
                            .map(|cell| {
                                (
                                    thermal.specific_enthalpy_at(cell.coordinates),
                                    thermal.flow().lattice().surface_tension_at(
                                        cell.coordinates[0],
                                        cell.coordinates[1],
                                        cell.coordinates[2],
                                    ),
                                )
                            })
                            .collect::<Vec<_>>(),
                    )
                };
                let before = snapshot(thermal);
                let local = [FreeSurfaceCellHeat {
                    coordinates: [0, 2, 3],
                    energy_j: 0.02,
                }];
                assert!(thermal.step_with_heat(&local, gate).is_err());
                assert_eq!(snapshot(thermal), before);
                assert!(
                    thermal
                        .step_with_conduction(&[], &local, controls(), gate)
                        .is_err()
                );
                assert_eq!(snapshot(thermal), before);
                if sigma_upper < 500.0 || angle_upper < 500.0 {
                    let heat = homogeneous_heat(thermal);
                    assert!(thermal.step_with_heat(&heat, gate).is_err());
                    assert_eq!(snapshot(thermal), before);
                }
                let cancelled = CancelGate::new_clock_free();
                cancelled.request();
                assert!(thermal.step(&cancelled).is_err());
                assert_eq!(snapshot(thermal), before);
                thermal.step(gate).unwrap();
            },
        );
    }
}
