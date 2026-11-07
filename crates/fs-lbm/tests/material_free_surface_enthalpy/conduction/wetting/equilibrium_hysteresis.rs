//! Independent static thresholds preserve Young adhesion and drive native motion.

use super::super::gas::gas_model;
use super::*;
use fs_lbm::material::free_surface::wetting::{
    CoxVoinovLaw, EquilibriumWettingPatch, EquilibriumWettingRateLaw, FreeSurfaceWettingSample,
    MolecularKineticControls, MolecularKineticLaw,
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

fn with_thresholds(
    original: &InterfaceSystemCard,
    young: f64,
    endpoints: [[f64; 2]; 2],
    upper_temperature: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (_, source) in original.claims().claims_ordered() {
        let mut claim = source.clone();
        let shift = if claim.key.name() == EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY {
            0.004 * (1.0 + young.cos()) - 0.004
        } else {
            young - core::f64::consts::FRAC_PI_2
        };
        match &mut claim.value {
            PropertyValue::Scalar { value, .. } => *value += shift,
            PropertyValue::Curve { knots, .. } => {
                for (_, value) in knots {
                    *value += shift;
                }
            }
            _ => panic!("synthetic equilibrium fixture must be scalar or a temperature curve"),
        }
        claims.insert_claim(claim).unwrap();
    }
    for (name, [reference, slope]) in [
        (RECEDING_CONTACT_ANGLE_PROPERTY, endpoints[0]),
        (ADVANCING_CONTACT_ANGLE_PROPERTY, endpoints[1]),
    ] {
        let mut claim = PropertyClaim {
            key: PropertyKey::new(name, Dims::NONE),
            value: if slope == 0.0 {
                PropertyValue::Scalar {
                    value: reference,
                    dims: Dims::NONE,
                }
            } else {
                PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![
                        (399.0, reference - slope),
                        (500.0, reference + 100.0 * slope),
                    ],
                    dims: Dims::NONE,
                }
            },
            validity: ValidityDomain::unconstrained().with("T", 399.0, upper_temperature),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: if slope == 0.0 {
                InterpolationPolicy::ConstantWithinValidity
            } else {
                InterpolationPolicy::LinearInside
            },
            observations: vec![],
            provenance: Provenance {
                source: "synthetic independently sourced static pinning endpoint; not measured"
                    .into(),
                license: "test-only".into(),
                artifact: None,
            },
        };
        claims.insert_claim(claim.clone()).unwrap();
        if upper_temperature < 500.0 {
            claim.validity =
                ValidityDomain::unconstrained().with("T", upper_temperature + 0.000001, 500.0);
            claims.insert_claim(claim).unwrap();
        }
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

fn with_flow<R>(
    liquid: &MaterialCard,
    interface: &InterfaceSystemCard,
    origin: EquilibriumAngleOrigin,
    initial: &[MaterialFreeSurfaceInitialCell3],
    run: impl FnOnce(MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let source = liquid_source(liquid, false);
    let law = LiquidTemperatureLaw::try_from_material(
        liquid,
        source.material().clone(),
        "T",
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    let chart = caloric(liquid);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let thermal = plan_free_surface_from_material(&source, 0.001, 0.00002, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), initial, cx)
            .unwrap()
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap()
            .with_equilibrium_planar_wetting(
                PlanarWall3::XMin,
                interface,
                origin,
                MaterialPropertySelection::SingleClaimOnly,
                64,
                &gate,
            )
            .unwrap();
        run(thermal, &gate)
    })
}

fn moving_cells() -> Vec<MaterialFreeSurfaceInitialCell3> {
    cells()
        .into_iter()
        .map(|mut cell| {
            cell.velocity_m_s = [0.0, 0.00004, 0.0];
            cell
        })
        .collect()
}

fn heat(initial: &[MaterialFreeSurfaceInitialCell3]) -> Vec<FreeSurfaceCellHeat> {
    initial
        .iter()
        .map(|cell| FreeSurfaceCellHeat {
            coordinates: cell.coordinates,
            energy_j: 0.04 * cell.fill_fraction,
        })
        .collect()
}

fn rate<'s, 'm>(
    thermal: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    model: usize,
    gas: &ResolvedThermochemicalGasState,
    gate: &CancelGate,
) -> MaterialFreeSurfaceEnthalpy3<'s, 'm> {
    match model {
        0 => thermal
            .with_dynamic_equilibrium_wetting(
                EquilibriumWettingRateLaw::CoxVoinov(
                    CoxVoinovLaw::try_new(outer_gas::cox_controls(), gate).unwrap(),
                ),
                gate,
            )
            .unwrap(),
        1 => thermal
            .with_dynamic_equilibrium_wetting(
                EquilibriumWettingRateLaw::MolecularKinetic(
                    MolecularKineticLaw::try_new(MolecularKineticControls {
                        jump_length_m: 5e-10,
                        equilibrium_jump_frequency_per_s: 1e9,
                        minimum_temperature_k: 399.0,
                        maximum_temperature_k: 500.0,
                        maximum_speed_m_s: 100.0,
                        minimum_angle_radians: 0.1,
                        maximum_angle_radians: 3.0,
                    })
                    .unwrap(),
                ),
                gate,
            )
            .unwrap(),
        _ => thermal
            .with_dynamic_equilibrium_gas_wetting(outer_gas::cox_controls(), gas, gate)
            .unwrap(),
    }
}

fn branch_angle(sample: &FreeSurfaceWettingSample) -> f64 {
    let h = sample.equilibrium_hysteresis.as_ref().unwrap();
    match sample.regime {
        ContactLineRegime2::Advancing => h.advancing_angle_radians,
        ContactLineRegime2::Receding => h.receding_angle_radians,
        ContactLineRegime2::Pinned => sample.geometric_angle_radians,
        ContactLineRegime2::Equilibrium => h.equilibrium_angle_radians,
    }
}

#[test]
fn g1_sourced_rest_pinning_and_both_depinning_directions_reach_native_wall_mass() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let original = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.0,
        500.0,
    );
    let initial = (0..8)
        .flat_map(|x| {
            (0..4).flat_map(move |y| {
                (0..8).map(move |z| MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: if y == 3 { 0.5 } else { 1.0 },
                    velocity_m_s: [0.0; 3],
                })
            })
        })
        .collect::<Vec<_>>();
    for (young, interval, regime, direction) in [
        (1.3, [1.2, 1.9], ContactLineRegime2::Pinned, 0.0),
        (1.3, [1.2, 1.4], ContactLineRegime2::Advancing, 1.0),
        (1.8, [1.75, 1.9], ContactLineRegime2::Receding, -1.0),
    ] {
        let interface = with_thresholds(
            &original,
            young,
            [[interval[0], 0.0], [interval[1], 0.0]],
            500.0,
        );
        with_flow(
            &liquid,
            &interface,
            EquilibriumAngleOrigin::YoungDupre,
            &initial,
            |thermal, gate| {
                let mut thermal = thermal
                    .with_equilibrium_wetting_hysteresis(
                        MaterialPropertySelection::SingleClaimOnly,
                        1e-10,
                        gate,
                    )
                    .unwrap();
                let mut resolved = 0;
                for sample in thermal.wetting_samples().unwrap().values() {
                    let h = sample.equilibrium_hysteresis.as_ref().unwrap();
                    assert_eq!(sample.signed_speed_m_s, 0.0);
                    assert_eq!(h.model_speed_m_s, 0.0);
                    assert_eq!(h.pinning_dissipation_per_length_w_m, 0.0);
                    if h.contact_line_resolved {
                        near(
                            sample.geometric_angle_radians,
                            core::f64::consts::FRAC_PI_2,
                            2e-13,
                        );
                        assert_eq!(sample.regime, regime);
                        resolved += 1;
                    } else {
                        assert_eq!(sample.regime, ContactLineRegime2::Equilibrium);
                    }
                }
                assert!(resolved > 0);
                let before = thermal
                    .flow()
                    .region_mass_kg(|x, _, z| x == 0 && (1..7).contains(&z))
                    .unwrap();
                let energy = thermal.total_enthalpy_j().unwrap();
                for _ in 0..4 {
                    let report = thermal.step(gate).unwrap();
                    assert_eq!(report.external_heat_j, 0.0);
                    assert_eq!(report.viscous_heat_j, 0.0);
                    near(report.mass_before_kg, report.mass_after_kg, 2e-12);
                }
                let change = thermal
                    .flow()
                    .region_mass_kg(|x, _, z| x == 0 && (1..7).contains(&z))
                    .unwrap()
                    - before;
                if direction == 0.0 {
                    assert!(
                        change.abs() < 1e-15,
                        "a flat pinned interface must remain at rest: {change}"
                    );
                } else {
                    assert!(
                        direction * change > 1e-16,
                        "sourced depinning must move actual wall-adjacent liquid in its declared direction: {change}"
                    );
                }
                near(thermal.total_enthalpy_j().unwrap(), energy, 2e-10);
            },
        );
    }
}

#[test]
fn g1_patch_thresholds_keep_equilibrium_sources_through_thermal_and_rate_feedback() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let base = with_thresholds(
        &equilibrium_interface(
            &liquid,
            &solid,
            EquilibriumAngleOrigin::YoungDupre,
            0.00002,
            500.0,
        ),
        core::f64::consts::FRAC_PI_2,
        [[1.2, 0.001], [1.9, -0.001]],
        500.0,
    );
    let patch_card = with_thresholds(
        &equilibrium_interface(
            &liquid,
            &solid,
            EquilibriumAngleOrigin::SourcedAngle,
            -0.002,
            500.0,
        ),
        core::f64::consts::FRAC_PI_2,
        [[1.3, -0.001], [1.7, 0.001]],
        500.0,
    );
    let patch = EquilibriumWettingPatch {
        cells: (4..8)
            .flat_map(|y| (0..8).map(move |z| [0, y, z]))
            .collect(),
        interface: &patch_card,
        origin: EquilibriumAngleOrigin::SourcedAngle,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    let initial = moving_cells();
    for model in 0..3 {
        let mut gas =
            IsochoricGasParcel::new(gas_model(5000.0), gas_source::point(450.0, 100000.0), 1e-8)
                .unwrap();
        let outer = gas.state().clone();
        with_flow(
            &liquid,
            &base,
            EquilibriumAngleOrigin::YoungDupre,
            &initial,
            |thermal, gate| {
                let patterned = thermal
                    .with_equilibrium_wetting_patches(std::slice::from_ref(&patch), gate)
                    .unwrap();
                let mut reference = rate(patterned.clone(), model, &outer, gate);
                let mut thermal = rate(
                    patterned
                        .with_equilibrium_wetting_hysteresis(
                            MaterialPropertySelection::SingleClaimOnly,
                            0.0,
                            gate,
                        )
                        .unwrap(),
                    model,
                    &outer,
                    gate,
                );
                let originals = thermal.wetting_samples().unwrap().clone();
                if model == 2 {
                    let mut reference_gas = gas.clone();
                    thermal
                        .step_with_finite_gas_environment(
                            &heat(&initial),
                            &[],
                            controls(),
                            &mut [],
                            &[],
                            &[],
                            &mut gas,
                            &[gas_face()],
                            gas_controls(),
                            false,
                            gate,
                        )
                        .unwrap();
                    reference
                        .step_with_finite_gas_environment(
                            &heat(&initial),
                            &[],
                            controls(),
                            &mut [],
                            &[],
                            &[],
                            &mut reference_gas,
                            &[gas_face()],
                            gas_controls(),
                            false,
                            gate,
                        )
                        .unwrap();
                    assert!(gas.state().state().temperature < outer.state().temperature);
                } else {
                    thermal.step_with_heat(&heat(&initial), gate).unwrap();
                    reference.step_with_heat(&heat(&initial), gate).unwrap();
                }
                let mut moving = [0; 2];
                for (&position, sample) in thermal.wetting_samples().unwrap() {
                    let h = sample.equilibrium_hysteresis.as_ref().unwrap();
                    let old_h = originals[&position]
                        .equilibrium_hysteresis
                        .as_ref()
                        .unwrap();
                    assert_eq!(sample.material.properties().len(), 1);
                    assert_eq!(h.thresholds.properties().len(), 2);
                    assert_eq!(
                        h.thresholds.card_identity(),
                        sample.material.card_identity()
                    );
                    assert_eq!(
                        sample.material.properties()[0].answer().receipt.selected,
                        originals[&position].material.properties()[0]
                            .answer()
                            .receipt
                            .selected
                    );
                    for name in [
                        RECEDING_CONTACT_ANGLE_PROPERTY,
                        ADVANCING_CONTACT_ANGLE_PROPERTY,
                    ] {
                        assert_eq!(
                            h.thresholds
                                .property(name)
                                .unwrap()
                                .answer()
                                .receipt
                                .selected,
                            old_h
                                .thresholds
                                .property(name)
                                .unwrap()
                                .answer()
                                .receipt
                                .selected
                        );
                    }
                    assert_eq!(
                        sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
                        None
                    );
                    assert_eq!(
                        h.thresholds.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
                        None
                    );
                    let t = sample.material.state_coordinate("T").unwrap();
                    near(h.thresholds.state_coordinate("T").unwrap(), t, 2e-13);
                    if let Some(state) = thermal.material_at(position) {
                        near(t, temperature(state), 2e-12);
                    }
                    let (receding, advancing) = if position[1] >= 4 {
                        (1.3 - 0.001 * (t - 400.0), 1.7 + 0.001 * (t - 400.0))
                    } else {
                        (1.2 + 0.001 * (t - 400.0), 1.9 - 0.001 * (t - 400.0))
                    };
                    near(h.receding_angle_radians, receding, 2e-13);
                    near(h.advancing_angle_radians, advancing, 2e-13);
                    let branch = branch_angle(sample);
                    let pinning = 0.004 * (h.equilibrium_angle_radians.cos() - branch.cos());
                    assert!((h.pinning_force_per_length_n_m - pinning).abs() < 1e-17);
                    assert!(h.pinning_dissipation_per_length_w_m >= 0.0);
                    if model == 1 {
                        let kinetics = sample.molecular_kinetic.unwrap();
                        assert_eq!(kinetics.signed_speed_m_s, h.model_speed_m_s);
                        let force =
                            2.0 * 1.380_649e-23 * t / 5e-10_f64.powi(2) * h.model_speed_m_s.asinh();
                        near(
                            sample.applied_angle_radians,
                            (branch.cos() - force / 0.004).acos(),
                            2e-12,
                        );
                    } else {
                        let cox = sample.hydrodynamic.unwrap();
                        near(
                            cox.signed_capillary_number,
                            0.1 * sample.signed_speed_m_s / 0.004,
                            2e-12,
                        );
                        near(
                            cox.model_signed_capillary_number,
                            0.1 * h.model_speed_m_s / 0.004,
                            2e-12,
                        );
                        let ratio = if model == 2 {
                            assert_eq!(sample.outer_gas.as_deref(), Some(gas.state()));
                            gas.state().state().dynamic_viscosity / 0.1
                        } else {
                            0.0
                        };
                        near(cox.outer_to_liquid_viscosity_ratio, ratio, 2e-12);
                        let forward =
                            outer_gas::cox_difference(branch, sample.applied_angle_radians, ratio);
                        assert!(
                            (forward - cox.model_signed_capillary_number * 1e6_f64.ln()).abs()
                                < 3e-10
                        );
                    }
                    if h.model_speed_m_s.abs() > 1e-9 {
                        moving[usize::from(position[1] >= 4)] += 1;
                    }
                }
                assert!(moving.into_iter().all(|count| count > 0));
                thermal.step(gate).unwrap();
                reference.step(gate).unwrap();
                assert_ne!(
                    thermal.flow().lattice().grid().state_bits(),
                    reference.flow().lattice().grid().state_bits(),
                    "separate sourced thresholds must affect the actual native flow"
                );
            },
        );
    }
}

#[test]
fn g3_zero_width_thresholds_recover_smooth_kinetics_despite_a_finite_deadband() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let outer = gas_model(5000.0)
        .evaluate(&gas_source::point(450.0, 100000.0))
        .unwrap();
    let initial = moving_cells();
    for origin in [
        EquilibriumAngleOrigin::YoungDupre,
        EquilibriumAngleOrigin::SourcedAngle,
    ] {
        let interface = with_thresholds(
            &equilibrium_interface(&liquid, &solid, origin, 0.0, 500.0),
            core::f64::consts::FRAC_PI_2,
            [[core::f64::consts::FRAC_PI_2, 0.0]; 2],
            500.0,
        );
        for model in 0..3 {
            with_flow(&liquid, &interface, origin, &initial, |thermal, gate| {
                let mut smooth = rate(thermal.clone(), model, &outer, gate);
                // Also prove attachment after rate selection admits the same
                // source/calibration domain without selecting a second model.
                let mut threshold = rate(thermal, model, &outer, gate)
                    .with_equilibrium_wetting_hysteresis(
                        MaterialPropertySelection::SingleClaimOnly,
                        1.0,
                        gate,
                    )
                    .unwrap();
                smooth.step_with_heat(&heat(&initial), gate).unwrap();
                threshold.step_with_heat(&heat(&initial), gate).unwrap();
                assert_eq!(
                    smooth.flow().lattice().grid().state_bits(),
                    threshold.flow().lattice().grid().state_bits()
                );
                for (&position, sample) in threshold.wetting_samples().unwrap() {
                    let h = sample.equilibrium_hysteresis.as_ref().unwrap();
                    assert_eq!(h.pinning_force_per_length_n_m, 0.0);
                    assert_eq!(h.pinning_dissipation_per_length_w_m, 0.0);
                    assert_eq!(
                        sample.applied_angle_radians,
                        smooth.wetting_samples().unwrap()[&position].applied_angle_radians
                    );
                }
            });
        }
    }
}

#[test]
fn g4_patch_threshold_source_exit_or_loss_of_equilibrium_order_rolls_back_the_whole_tick() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let original = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.0,
        500.0,
    );
    let base = with_thresholds(
        &original,
        core::f64::consts::FRAC_PI_2,
        [[1.2, 0.0], [1.9, 0.0]],
        500.0,
    );
    let initial = moving_cells();
    for (upper, advancing_slope) in [(410.0, 0.0), (500.0, -0.02)] {
        let patch_card = with_thresholds(
            &original,
            core::f64::consts::FRAC_PI_2,
            [[1.2, 0.0], [1.7, advancing_slope]],
            upper,
        );
        let patch = EquilibriumWettingPatch {
            cells: vec![WALL_CELL],
            interface: &patch_card,
            origin: EquilibriumAngleOrigin::YoungDupre,
            selection: MaterialPropertySelection::SingleClaimOnly,
        };
        with_flow(
            &liquid,
            &base,
            EquilibriumAngleOrigin::YoungDupre,
            &initial,
            |thermal, gate| {
                let mut thermal = thermal
                    .with_equilibrium_wetting_patches(std::slice::from_ref(&patch), gate)
                    .unwrap()
                    .with_equilibrium_wetting_hysteresis(
                        MaterialPropertySelection::SingleClaimOnly,
                        0.0,
                        gate,
                    )
                    .unwrap();
                assert!(
                    thermal
                        .clone()
                        .with_equilibrium_wetting_hysteresis(
                            MaterialPropertySelection::SingleClaimOnly,
                            0.0,
                            gate,
                        )
                        .is_err()
                );
                assert!(
                    thermal
                        .clone()
                        .with_equilibrium_wetting_patches(&[patch], gate)
                        .is_err()
                );
                let before = (
                    thermal.flow().lattice().grid().state_bits(),
                    thermal.wetting_samples().unwrap().clone(),
                    thermal.total_enthalpy_j().unwrap(),
                    thermal.physical_time_s(),
                    initial
                        .iter()
                        .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                        .collect::<Vec<_>>(),
                );
                let mut retry = thermal.clone();
                assert!(thermal.step_with_heat(&heat(&initial), gate).is_err());
                let cancelled = CancelGate::new_clock_free();
                cancelled.request();
                assert!(thermal.step(&cancelled).is_err());
                assert_eq!(thermal.flow().lattice().grid().state_bits(), before.0);
                assert_eq!(thermal.wetting_samples().unwrap(), &before.1);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), before.2);
                assert_eq!(thermal.physical_time_s(), before.3);
                assert_eq!(
                    initial
                        .iter()
                        .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                        .collect::<Vec<_>>(),
                    before.4
                );
                thermal.step(gate).unwrap();
                retry.step(gate).unwrap();
                assert_eq!(
                    thermal.flow().lattice().grid().state_bits(),
                    retry.flow().lattice().grid().state_bits()
                );
                assert_eq!(thermal.wetting_samples(), retry.wetting_samples());
            },
        );
    }
}
