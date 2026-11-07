//! Source-resolved outer gas and original equilibrium adhesion share one trajectory.

use super::super::gas::gas_model;
use super::*;
use fs_lbm::material::free_surface::wetting::{
    CoxVoinovControls, CoxVoinovLaw, EquilibriumWettingRateLaw,
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

pub(super) fn cox_controls() -> CoxVoinovControls {
    CoxVoinovControls {
        microscopic_cutoff_m: 1e-9,
        observation_length_m: 0.001,
        maximum_capillary_number: 0.01,
        minimum_angle_radians: 0.01,
        maximum_angle_radians: 3.0,
    }
}

fn fixed_ratio(
    outer: &ResolvedThermochemicalGasState,
    gate: &CancelGate,
) -> EquilibriumWettingRateLaw {
    EquilibriumWettingRateLaw::CoxVoinov(
        CoxVoinovLaw::try_new_with_viscosity_ratio(
            cox_controls(),
            outer.state().dynamic_viscosity / 0.1,
            gate,
        )
        .unwrap(),
    )
}

// Direct Simpson integration of the original two-viscous-wedge equation,
// independent of the production angular table, local mapping and inverse.
pub(super) fn cox_difference(from: f64, to: f64, ratio: f64) -> f64 {
    let intervals = 256;
    let h = (to - from) / intervals as f64;
    let slope = |theta: f64| {
        let other = core::f64::consts::PI - theta;
        let sin = theta.sin();
        let cos = theta.cos();
        let liquid = theta * theta - sin * sin;
        let outer = other * other - sin * sin;
        (ratio * liquid * (other + sin * cos) + outer * (theta - sin * cos))
            / (2.0
                * sin
                * (ratio * ratio * liquid + 2.0 * ratio * (theta * other + sin * sin) + outer))
    };
    let mut sum = slope(from) + slope(to);
    for i in 1..intervals {
        sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * slope(from + i as f64 * h);
    }
    sum * h / 3.0
}

fn check_sources_and_forward_law(
    thermal: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    outer: &ResolvedThermochemicalGasState,
    origin: EquilibriumAngleOrigin,
) -> usize {
    assert_eq!(thermal.wetting_outer_gas(), Some(outer));
    let original = thermal.wetting_source().unwrap();
    let mu_g = 5000.0 * 1.8e-5 * (0.5 + (outer.state().pressure - 1000.0) / 249000.0);
    near(outer.state().dynamic_viscosity, mu_g, 2e-13);
    let mut moving = 0;
    for (&position, sample) in thermal.wetting_samples().unwrap() {
        assert_eq!(sample.equilibrium, Some(origin));
        assert_eq!(sample.outer_gas.as_deref(), Some(outer));
        assert_eq!(sample.material.properties().len(), 1);
        assert_eq!(
            sample.material.state_coordinate(DYNAMIC_WETTING_SPEED_AXIS),
            None
        );
        assert_eq!(
            sample.material.properties()[0].answer().receipt.selected,
            original.properties()[0].answer().receipt.selected
        );
        let t = sample.material.state_coordinate("T").unwrap();
        if let Some(liquid) = thermal.material_at(position) {
            near(t, temperature(liquid), 2e-12);
        }
        let microscopic = match origin {
            EquilibriumAngleOrigin::YoungDupre => {
                ((0.004 + 0.00002 * (t - 400.0)) / 0.004 - 1.0).acos()
            }
            EquilibriumAngleOrigin::SourcedAngle => {
                core::f64::consts::FRAC_PI_2 - 0.002 * (t - 400.0)
            }
        };
        let prediction = sample.hydrodynamic.unwrap();
        let ratio = mu_g / 0.1;
        let capillary = 0.1 * sample.signed_speed_m_s / 0.004;
        near(prediction.dynamic_viscosity_pa_s, 0.1, 2e-13);
        near(prediction.surface_tension_n_m, 0.004, 2e-13);
        near(prediction.outer_to_liquid_viscosity_ratio, ratio, 2e-13);
        near(prediction.signed_capillary_number, capillary, 2e-13);
        let expected = capillary * 1e6_f64.ln();
        let observed = cox_difference(microscopic, sample.applied_angle_radians, ratio);
        assert!(
            (observed - expected).abs() < 3e-10,
            "two-phase forward equation: {observed:e} != {expected:e}"
        );
        if sample.signed_speed_m_s.abs() > 1e-9 {
            moving += 1;
            assert_eq!(
                sample.regime,
                if sample.signed_speed_m_s > 0.0 {
                    ContactLineRegime2::Advancing
                } else {
                    ContactLineRegime2::Receding
                }
            );
        }
    }
    moving
}

#[test]
fn g1_sourced_outer_gas_retains_equilibrium_claims_and_fixed_snapshot_parity() {
    for origin in [
        EquilibriumAngleOrigin::YoungDupre,
        EquilibriumAngleOrigin::SourcedAngle,
    ] {
        let outer = gas_model(5000.0)
            .evaluate(&gas_source::point(450.0, 100000.0))
            .unwrap();
        dynamic::with_source(origin, 500.0, |initial, gate| {
            let mut reference = initial
                .clone()
                .with_dynamic_equilibrium_wetting(fixed_ratio(&outer, gate), gate)
                .unwrap();
            let mut thermal = initial
                .with_dynamic_equilibrium_gas_wetting(cox_controls(), &outer, gate)
                .unwrap();
            let before = thermal.total_enthalpy_j().unwrap();
            let heat = cells()
                .iter()
                .map(|cell| FreeSurfaceCellHeat {
                    coordinates: cell.coordinates,
                    energy_j: 0.04,
                })
                .collect::<Vec<_>>();
            let supplied = heat.iter().map(|entry| entry.energy_j).sum::<f64>();
            let report = thermal.step_with_heat(&heat, gate).unwrap();
            reference.step_with_heat(&heat, gate).unwrap();
            near(thermal.total_enthalpy_j().unwrap() - before, supplied, 2e-9);
            near(report.external_heat_j, supplied, 2e-12);
            assert_eq!(report.viscous_heat_j, 0.0);
            assert!(check_sources_and_forward_law(&thermal, &outer, origin) > 0);
            assert_eq!(
                thermal.wetting_outer_gas().unwrap().state().temperature,
                450.0
            );
            for sample in thermal.wetting_samples().unwrap().values() {
                assert_ne!(sample.material.state_coordinate("T"), Some(450.0));
            }
            assert_eq!(
                thermal.flow().lattice().grid().state_bits(),
                reference.flow().lattice().grid().state_bits()
            );
            thermal.step(gate).unwrap();
            reference.step(gate).unwrap();
            assert_eq!(
                thermal.flow().lattice().grid().state_bits(),
                reference.flow().lattice().grid().state_bits(),
                "an unchanged sourced ratio must preserve the original fixed-ratio physics"
            );
        });
    }
}

#[test]
fn g1_finite_gas_changes_equilibrium_adhesion_angles_and_following_native_flow() {
    for origin in [
        EquilibriumAngleOrigin::YoungDupre,
        EquilibriumAngleOrigin::SourcedAngle,
    ] {
        let mut atmosphere =
            IsochoricGasParcel::new(gas_model(5000.0), gas_source::point(450.0, 100000.0), 1e-8)
                .unwrap();
        let original_gas = atmosphere.state().clone();
        dynamic::with_source(origin, 500.0, |initial, gate| {
            let mut reference = initial
                .clone()
                .with_dynamic_equilibrium_wetting(fixed_ratio(&original_gas, gate), gate)
                .unwrap();
            let mut reference_gas = atmosphere.clone();
            let mut thermal = initial
                .with_dynamic_equilibrium_gas_wetting(cox_controls(), &original_gas, gate)
                .unwrap();
            let before = thermal.total_enthalpy_j().unwrap()
                + atmosphere.mass_kg() * atmosphere.state().specific_internal_energy_j_kg();
            let report = thermal
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
                .unwrap()
                .conduction;
            reference
                .step_with_finite_gas_environment(
                    &[],
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
            let evolved = atmosphere.state();
            assert!(evolved.state().temperature < original_gas.state().temperature);
            assert!(evolved.state().pressure < original_gas.state().pressure);
            assert!(evolved.state().dynamic_viscosity < original_gas.state().dynamic_viscosity);
            assert_eq!(reference_gas.state(), evolved);
            assert!(check_sources_and_forward_law(&thermal, evolved, origin) > 0);
            let after = thermal.total_enthalpy_j().unwrap()
                + atmosphere.mass_kg() * evolved.specific_internal_energy_j_kg();
            near(after, before, 2e-10);
            assert!(
                report
                    .finite_gas
                    .unwrap()
                    .combined_thermal_balance_defect_j
                    .abs()
                    < 1e-9
            );
            assert_eq!(
                thermal.flow().lattice().grid().state_bits(),
                reference.flow().lattice().grid().state_bits(),
                "initial zero-speed contact angles and the same thermal exchange coincide"
            );
            let angle_difference = thermal
                .wetting_samples()
                .unwrap()
                .iter()
                .map(|(position, sample)| {
                    (sample.applied_angle_radians
                        - reference.wetting_samples().unwrap()[position].applied_angle_radians)
                        .abs()
                })
                .fold(0.0_f64, f64::max);
            assert!(
                angle_difference > 1e-12,
                "the original gas viscosity must alter the moving boundary: {angle_difference}"
            );
            thermal.step(gate).unwrap();
            reference.step(gate).unwrap();
            let difference = thermal
                .flow()
                .lattice()
                .grid()
                .state_bits()
                .iter()
                .zip(reference.flow().lattice().grid().state_bits())
                .map(|(&a, b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
                .fold(0.0_f64, f64::max);
            assert!(
                difference > 1e-13,
                "updated gas-controlled adhesion must reach native momentum: {difference}"
            );
        });
    }
}

#[test]
fn g3_g4_equilibrium_gas_refusals_preserve_the_original_thermal_and_flow_state() {
    let outer = gas_model(5000.0)
        .evaluate(&gas_source::point(450.0, 100000.0))
        .unwrap();
    dynamic::with_source(
        EquilibriumAngleOrigin::YoungDupre,
        500.0,
        |initial, gate| {
            assert!(
                initial
                    .clone()
                    .with_cox_voinov_gas_wetting(cox_controls(), &outer, gate)
                    .is_err()
            );
            let mut thermal = initial
                .with_dynamic_equilibrium_gas_wetting(cox_controls(), &outer, gate)
                .unwrap();
            assert!(
                thermal
                    .clone()
                    .with_dynamic_equilibrium_gas_wetting(cox_controls(), &outer, gate)
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
            let mut foreign =
                IsochoricGasParcel::new(gas_model(1.0), gas_source::point(450.0, 100000.0), 1e-8)
                    .unwrap();
            let gas_before = foreign.state().clone();
            assert!(
                thermal
                    .step_with_finite_gas_environment(
                        &[],
                        &[],
                        controls(),
                        &mut [],
                        &[],
                        &[],
                        &mut foreign,
                        &[gas_face()],
                        gas_controls(),
                        false,
                        gate,
                    )
                    .is_err()
            );
            assert_eq!(foreign.state(), &gas_before);
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(thermal.step(&cancelled).is_err());
            assert_eq!(
                (
                    thermal.flow().lattice().grid().state_bits(),
                    thermal.wetting_samples().unwrap().clone(),
                    thermal.total_enthalpy_j().unwrap(),
                    thermal.physical_time_s(),
                    cells()
                        .iter()
                        .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                        .collect::<Vec<_>>(),
                ),
                before
            );
            assert_eq!(thermal.wetting_outer_gas(), Some(&outer));
            thermal.step(gate).unwrap();
        },
    );
}
