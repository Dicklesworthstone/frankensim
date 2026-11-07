//! Substrate chemistry changes the existing capillary boundary and native flow.

use super::super::gas::gas_model;
use super::*;
use fs_lbm::material::free_surface::wetting::{
    EquilibriumWettingPatch, EquilibriumWettingRateLaw, FreeSurfaceWettingSample,
    MolecularKineticControls, MolecularKineticLaw,
};

fn patch<'m>(
    interface: &'m InterfaceSystemCard,
    origin: EquilibriumAngleOrigin,
    rows: std::ops::Range<usize>,
) -> EquilibriumWettingPatch<'m> {
    EquilibriumWettingPatch {
        cells: rows.flat_map(|y| (0..8).map(move |z| [0, y, z])).collect(),
        interface,
        origin,
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}

fn heat() -> Vec<FreeSurfaceCellHeat> {
    cells()
        .iter()
        .map(|cell| FreeSurfaceCellHeat {
            coordinates: cell.coordinates,
            energy_j: 0.04,
        })
        .collect()
}

fn microscopic(sample: &FreeSurfaceWettingSample) -> f64 {
    match sample.equilibrium.unwrap() {
        EquilibriumAngleOrigin::SourcedAngle => sample
            .material
            .property(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY)
            .unwrap()
            .value_si(),
        EquilibriumAngleOrigin::YoungDupre => (sample
            .material
            .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
            .unwrap()
            .value_si()
            / 0.004
            - 1.0)
            .acos(),
    }
}

#[test]
fn g1_g3_patterned_adhesion_drives_native_flow_and_is_patch_order_independent() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let base = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.0,
        500.0,
    );
    let work = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.00002,
        500.0,
    );
    let angle = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::SourcedAngle,
        -0.002,
        500.0,
    );
    with_equilibrium(
        &liquid,
        &base,
        EquilibriumAngleOrigin::YoungDupre,
        |initial, gate| {
            let patches = [
                patch(&work, EquilibriumAngleOrigin::YoungDupre, 1..4),
                patch(&angle, EquilibriumAngleOrigin::SourcedAngle, 4..7),
            ];
            let mut patterned = initial
                .clone()
                .with_equilibrium_wetting_patches(&patches, gate)
                .unwrap();
            let mut reordered = initial
                .clone()
                .with_equilibrium_wetting_patches(&[patches[1].clone(), patches[0].clone()], gate)
                .unwrap();
            let mut homogeneous = initial
                .clone()
                .with_equilibrium_wetting_patches(
                    &[patch(&base, EquilibriumAngleOrigin::YoungDupre, 1..7)],
                    gate,
                )
                .unwrap();
            let originals = patterned.wetting_samples().unwrap().clone();
            let before = patterned.total_enthalpy_j().unwrap();
            let supplied = heat().iter().map(|h| h.energy_j).sum::<f64>();
            let report = patterned.step_with_heat(&heat(), gate).unwrap();
            reordered.step_with_heat(&heat(), gate).unwrap();
            homogeneous.step_with_heat(&heat(), gate).unwrap();
            initial.step_with_heat(&heat(), gate).unwrap();
            near(
                patterned.total_enthalpy_j().unwrap() - before,
                supplied,
                2e-9,
            );
            assert_eq!(report.viscous_heat_j, 0.0);
            assert_eq!(patterned.wetting_samples(), reordered.wetting_samples());
            assert_eq!(
                patterned.flow().lattice().grid().state_bits(),
                reordered.flow().lattice().grid().state_bits()
            );
            assert_eq!(homogeneous.wetting_samples(), initial.wetting_samples());
            assert_eq!(
                homogeneous.flow().lattice().grid().state_bits(),
                initial.flow().lattice().grid().state_bits()
            );
            let mut heated_patches = [0; 2];
            for (&position, sample) in patterned.wetting_samples().unwrap() {
                let expected_card = match position[1] {
                    1..=3 => work.content_hash(),
                    4..=6 => angle.content_hash(),
                    _ => base.content_hash(),
                };
                assert_eq!(sample.material.card_identity(), expected_card);
                assert_eq!(
                    sample.material.properties()[0].answer().receipt.selected,
                    originals[&position].material.properties()[0]
                        .answer()
                        .receipt
                        .selected
                );
                assert_eq!(
                    sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
                    None
                );
                let t = sample.material.state_coordinate("T").unwrap();
                if let Some(state) = patterned.material_at(position) {
                    near(t, temperature(state), 2e-12);
                }
                let expected = match position[1] {
                    1..=3 => {
                        assert_eq!(sample.equilibrium, Some(EquilibriumAngleOrigin::YoungDupre));
                        if t > 401.0 {
                            heated_patches[0] += 1;
                        }
                        ((0.004 + 0.00002 * (t - 400.0)) / 0.004 - 1.0).acos()
                    }
                    4..=6 => {
                        assert_eq!(
                            sample.equilibrium,
                            Some(EquilibriumAngleOrigin::SourcedAngle)
                        );
                        if t > 401.0 {
                            heated_patches[1] += 1;
                        }
                        core::f64::consts::FRAC_PI_2 - 0.002 * (t - 400.0)
                    }
                    _ => core::f64::consts::FRAC_PI_2,
                };
                near(sample.applied_angle_radians, expected, 2e-12);
            }
            assert!(heated_patches.into_iter().all(|count| count > 0));
            let change = patterned
                .flow()
                .lattice()
                .grid()
                .state_bits()
                .iter()
                .zip(initial.flow().lattice().grid().state_bits())
                .map(|(&a, b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
                .fold(0.0_f64, f64::max);
            assert!(
                change > 1e-12,
                "spatial adhesion must reach native capillary momentum: {change}"
            );
        },
    );
}

#[test]
fn g1_patterned_molecular_wetting_uses_each_original_equilibrium_and_temperature() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let angle = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::SourcedAngle,
        -0.002,
        500.0,
    );
    dynamic::with_source(
        EquilibriumAngleOrigin::YoungDupre,
        500.0,
        |initial, gate| {
            let mut thermal = initial
                .with_equilibrium_wetting_patches(
                    &[patch(&angle, EquilibriumAngleOrigin::SourcedAngle, 4..8)],
                    gate,
                )
                .unwrap()
                .with_dynamic_equilibrium_wetting(
                    EquilibriumWettingRateLaw::MolecularKinetic(
                        MolecularKineticLaw::try_new(MolecularKineticControls {
                            jump_length_m: 5e-10,
                            equilibrium_jump_frequency_per_s: 1e9,
                            minimum_temperature_k: 300.0,
                            maximum_temperature_k: 500.0,
                            maximum_speed_m_s: 100.0,
                            minimum_angle_radians: 0.1,
                            maximum_angle_radians: 3.0,
                        })
                        .unwrap(),
                    ),
                    gate,
                )
                .unwrap();
            let originals = thermal.wetting_samples().unwrap().clone();
            let report = thermal.step_with_heat(&heat(), gate).unwrap();
            assert_eq!(report.viscous_heat_j, 0.0);
            let mut moving = [0; 2];
            for (&position, sample) in thermal.wetting_samples().unwrap() {
                assert_eq!(
                    sample.material.properties()[0].answer().receipt.selected,
                    originals[&position].material.properties()[0]
                        .answer()
                        .receipt
                        .selected
                );
                let prediction = sample.molecular_kinetic.unwrap();
                assert!(sample.hydrodynamic.is_none());
                near(
                    prediction.temperature_k,
                    sample.material.state_coordinate("T").unwrap(),
                    2e-12,
                );
                let expected_force = 2.0 * 1.380_649e-23 * prediction.temperature_k
                    / 5e-10_f64.powi(2)
                    * sample.signed_speed_m_s.asinh();
                near(
                    prediction.driving_force_per_length_n_m,
                    expected_force,
                    2e-12,
                );
                near(
                    sample.applied_angle_radians,
                    (microscopic(sample).cos() - expected_force / 0.004).acos(),
                    2e-12,
                );
                if sample.signed_speed_m_s.abs() > 1e-9 {
                    moving[usize::from(position[1] >= 4)] += 1;
                }
            }
            assert!(moving.into_iter().all(|count| count > 0));
        },
    );
}

#[test]
fn g1_patterned_cox_wetting_refreshes_the_same_evolving_outer_gas_on_every_patch() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let angle = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::SourcedAngle,
        -0.002,
        500.0,
    );
    let mut atmosphere =
        IsochoricGasParcel::new(gas_model(5000.0), gas_source::point(450.0, 100000.0), 1e-8)
            .unwrap();
    let original_gas = atmosphere.state().clone();
    dynamic::with_source(
        EquilibriumAngleOrigin::YoungDupre,
        500.0,
        |initial, gate| {
            let mut thermal = initial
                .with_equilibrium_wetting_patches(
                    &[patch(&angle, EquilibriumAngleOrigin::SourcedAngle, 4..8)],
                    gate,
                )
                .unwrap()
                .with_dynamic_equilibrium_gas_wetting(
                    outer_gas::cox_controls(),
                    &original_gas,
                    gate,
                )
                .unwrap();
            let originals = thermal.wetting_samples().unwrap().clone();
            thermal
                .step_with_finite_gas_environment(
                    &[],
                    &[],
                    controls(),
                    &mut [],
                    &[],
                    &[],
                    &mut atmosphere,
                    &[gas_face()],
                    gas_controls(),
                    false,
                    gate,
                )
                .unwrap();
            let evolved = atmosphere.state();
            assert!(evolved.state().temperature < original_gas.state().temperature);
            assert!(evolved.state().dynamic_viscosity < original_gas.state().dynamic_viscosity);
            let mut moving = [0; 2];
            for (&position, sample) in thermal.wetting_samples().unwrap() {
                assert_eq!(sample.outer_gas.as_deref(), Some(evolved));
                assert_eq!(
                    sample.material.properties()[0].answer().receipt.selected,
                    originals[&position].material.properties()[0]
                        .answer()
                        .receipt
                        .selected
                );
                assert_eq!(
                    sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
                    None
                );
                let prediction = sample.hydrodynamic.unwrap();
                let ratio = evolved.state().dynamic_viscosity / 0.1;
                near(prediction.outer_to_liquid_viscosity_ratio, ratio, 2e-12);
                let observed = outer_gas::cox_difference(
                    microscopic(sample),
                    sample.applied_angle_radians,
                    ratio,
                );
                let expected = 0.1 * sample.signed_speed_m_s / 0.004 * 1e6_f64.ln();
                assert!((observed - expected).abs() < 3e-10);
                if sample.signed_speed_m_s.abs() > 1e-9 {
                    moving[usize::from(position[1] >= 4)] += 1;
                }
            }
            assert!(moving.into_iter().all(|count| count > 0));
        },
    );
}

#[test]
fn g3_g4_pattern_geometry_and_patch_source_exits_refuse_without_partial_updates() {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let base = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.0,
        500.0,
    );
    // A distinct fallback above 410 K exists but cannot replace the original.
    let limited = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.00002,
        410.0,
    );
    let foreign_liquid = InterfaceSystemCard::assemble(
        limited.surface_a().clone(),
        limited.surface_a().clone(),
        limited.context().clone(),
        limited.claims().clone(),
        vec![],
    )
    .unwrap();
    let mut foreign_context = limited.context().clone();
    foreign_context.environment = "synthetic incompatible surrounding gas".into();
    let foreign_gas = InterfaceSystemCard::assemble(
        limited.surface_a().clone(),
        limited.surface_b().clone(),
        foreign_context,
        limited.claims().clone(),
        vec![],
    )
    .unwrap();
    with_equilibrium(
        &liquid,
        &base,
        EquilibriumAngleOrigin::YoungDupre,
        |initial, gate| {
            let valid = patch(&limited, EquilibriumAngleOrigin::YoungDupre, 1..4);
            assert!(
                initial
                    .clone()
                    .with_equilibrium_wetting_patches(&[valid.clone(), valid.clone()], gate,)
                    .is_err()
            );
            for coordinates in [vec![], vec![[1, 2, 2]], vec![[0, 8, 2]]] {
                let mut invalid = valid.clone();
                invalid.cells = coordinates;
                assert!(
                    initial
                        .clone()
                        .with_equilibrium_wetting_patches(&[invalid], gate)
                        .is_err()
                );
            }
            for foreign in [&foreign_liquid, &foreign_gas] {
                assert!(
                    initial
                        .clone()
                        .with_equilibrium_wetting_patches(
                            &[patch(foreign, EquilibriumAngleOrigin::YoungDupre, 1..4)],
                            gate,
                        )
                        .is_err()
                );
            }
            let mut thermal = initial
                .clone()
                .with_equilibrium_wetting_patches(std::slice::from_ref(&valid), gate)
                .unwrap();
            assert!(
                thermal
                    .clone()
                    .with_equilibrium_wetting_patches(std::slice::from_ref(&valid), gate,)
                    .is_err()
            );
            let before = (
                thermal.flow().lattice().grid().state_bits(),
                thermal.wetting_samples().unwrap().clone(),
                thermal.total_enthalpy_j().unwrap(),
                thermal.physical_time_s(),
                cells()
                    .iter()
                    .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                    .collect::<Vec<_>>(),
            );
            let mut retry = thermal.clone();
            assert!(thermal.step_with_heat(&heat(), gate).is_err());
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(thermal.step(&cancelled).is_err());
            assert_eq!(thermal.flow().lattice().grid().state_bits(), before.0);
            assert_eq!(thermal.wetting_samples().unwrap(), &before.1);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), before.2);
            assert_eq!(thermal.physical_time_s(), before.3);
            assert_eq!(
                cells()
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
            initial.step(gate).unwrap();
            assert!(
                initial
                    .clone()
                    .with_equilibrium_wetting_patches(&[valid], gate)
                    .is_err()
            );
        },
    );
}
