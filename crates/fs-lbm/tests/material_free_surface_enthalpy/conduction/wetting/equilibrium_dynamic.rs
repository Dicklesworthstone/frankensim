//! Original reversible adhesion drives explicit moving-contact-line kinetics.

use super::*;
use fs_lbm::material::free_surface::wetting::{
    CoxVoinovControls, CoxVoinovLaw, EquilibriumWettingRateLaw, MolecularKineticControls,
    MolecularKineticError, MolecularKineticLaw,
};

fn rate_law(
    molecular: bool,
    maximum_temperature_k: f64,
    gate: &CancelGate,
) -> EquilibriumWettingRateLaw {
    if molecular {
        EquilibriumWettingRateLaw::MolecularKinetic(
            MolecularKineticLaw::try_new(MolecularKineticControls {
                jump_length_m: 5e-10,
                equilibrium_jump_frequency_per_s: 1e9,
                minimum_temperature_k: 300.0,
                maximum_temperature_k,
                maximum_speed_m_s: 100.0,
                minimum_angle_radians: 0.1,
                maximum_angle_radians: 3.0,
            })
            .unwrap(),
        )
    } else {
        EquilibriumWettingRateLaw::CoxVoinov(
            CoxVoinovLaw::try_new(
                CoxVoinovControls {
                    microscopic_cutoff_m: 1e-8,
                    observation_length_m: 0.001,
                    maximum_capillary_number: 0.01,
                    minimum_angle_radians: 0.1,
                    maximum_angle_radians: 3.0,
                },
                gate,
            )
            .unwrap(),
        )
    }
}

fn with_source<R>(
    origin: EquilibriumAngleOrigin,
    upper_temperature: f64,
    run: impl FnOnce(MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let liquid = source_card(0.004, false, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    let slope = match origin {
        EquilibriumAngleOrigin::YoungDupre => 0.00002,
        EquilibriumAngleOrigin::SourcedAngle => -0.002,
    };
    let interface = equilibrium_interface(&liquid, &solid, origin, slope, upper_temperature);
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
                // A small physical translation along the wall resolves motion
                // independently of contact-angle heating; Ca remains below 0.01.
                cell.velocity_m_s = [0.0, 0.00004, 0.0];
                cell
            })
            .collect::<Vec<_>>();
        let thermal = plan_free_surface_from_material(&source, 0.001, 0.00002, 0.1)
            .unwrap()
            .free_surface3_from_cells(domain(), &initial, cx)
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
                &interface,
                origin,
                MaterialPropertySelection::SingleClaimOnly,
                64,
                &gate,
            )
            .unwrap();
        run(thermal, &gate)
    })
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

// Independent full-angle inviscid-outer Cox integral. Simpson quadrature here
// checks the physical forward equation without the production inverse/table.
fn cox_difference(a: f64, b: f64) -> f64 {
    let intervals = 256;
    let h = (b - a) / intervals as f64;
    let integrand = |theta: f64| (theta - theta.sin() * theta.cos()) / (2.0 * theta.sin());
    let mut sum = integrand(a) + integrand(b);
    for i in 1..intervals {
        sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(a + i as f64 * h);
    }
    h * sum / 3.0
}

fn residual(actual: f64, expected: f64, absolute_tolerance: f64) {
    assert!(
        (actual - expected).abs() <= absolute_tolerance,
        "forward-equation residual: {actual} != {expected} (tolerance {absolute_tolerance})",
    );
}

#[test]
fn g1_dynamic_equilibrium_uses_original_adhesion_temperature_and_forward_rate_laws() {
    for origin in [
        EquilibriumAngleOrigin::YoungDupre,
        EquilibriumAngleOrigin::SourcedAngle,
    ] {
        for molecular in [false, true] {
            with_source(origin, 500.0, |thermal, gate| {
                let mut reference = thermal.clone();
                let mut thermal = thermal
                    .with_dynamic_equilibrium_wetting(rate_law(molecular, 500.0, gate), gate)
                    .unwrap();
                let original = thermal.wetting_source().unwrap().clone();
                assert_eq!(original.properties().len(), 1);
                for sample in thermal.wetting_samples().unwrap().values() {
                    assert_eq!(sample.regime, ContactLineRegime2::Equilibrium);
                    assert_eq!(sample.signed_speed_m_s, 0.0);
                    near(
                        sample.applied_angle_radians,
                        core::f64::consts::FRAC_PI_2,
                        2e-14,
                    );
                }
                let initial_h = thermal.total_enthalpy_j().unwrap();
                let heat = heat();
                let supplied = heat.iter().map(|entry| entry.energy_j).sum::<f64>();
                let report = thermal.step_with_heat(&heat, gate).unwrap();
                reference.step_with_heat(&heat, gate).unwrap();
                near(
                    thermal.total_enthalpy_j().unwrap() - initial_h,
                    supplied,
                    2e-9,
                );
                near(report.external_heat_j, supplied, 2e-12);
                assert_eq!(report.viscous_heat_j, 0.0);
                let mut moving = 0;
                let mut heated_moving = 0;
                for (&position, sample) in thermal.wetting_samples().unwrap() {
                    assert_eq!(sample.equilibrium, Some(origin));
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
                    if let Some(state) = thermal.material_at(position) {
                        near(t, temperature(state), 2e-12);
                    }
                    let theta = match origin {
                        EquilibriumAngleOrigin::YoungDupre => {
                            let w = 0.004 + 0.00002 * (t - 400.0);
                            (w / 0.004 - 1.0).acos()
                        }
                        EquilibriumAngleOrigin::SourcedAngle => {
                            core::f64::consts::FRAC_PI_2 - 0.002 * (t - 400.0)
                        }
                    };
                    if molecular {
                        let prediction = sample.molecular_kinetic.unwrap();
                        near(prediction.temperature_k, t, 2e-12);
                        near(prediction.signed_speed_m_s, sample.signed_speed_m_s, 2e-14);
                        let argument = 0.004
                            * 5e-10_f64.powi(2)
                            * (theta.cos() - sample.applied_angle_radians.cos())
                            / (2.0 * 1.380_649e-23 * t);
                        // Absolute residuals remain meaningful at the stationary
                        // samples, where subtracting two cosines loses digits.
                        residual(
                            2.0 * 1e9 * 5e-10 * argument.sinh(),
                            sample.signed_speed_m_s,
                            1e-13,
                        );
                        let force = 0.004 * (theta.cos() - sample.applied_angle_radians.cos());
                        residual(prediction.driving_force_per_length_n_m, force, 2e-15);
                        residual(
                            prediction.dissipation_per_length_w_m,
                            force * sample.signed_speed_m_s,
                            2e-15 * sample.signed_speed_m_s.abs(),
                        );
                        assert!(prediction.dissipation_per_length_w_m >= 0.0);
                        assert!(sample.hydrodynamic.is_none());
                    } else {
                        let prediction = sample.hydrodynamic.unwrap();
                        near(
                            prediction.signed_capillary_number,
                            prediction.dynamic_viscosity_pa_s * sample.signed_speed_m_s / 0.004,
                            2e-14,
                        );
                        residual(
                            cox_difference(theta, sample.applied_angle_radians),
                            prediction.signed_capillary_number * 100000.0_f64.ln(),
                            3e-11,
                        );
                        assert!(sample.molecular_kinetic.is_none());
                    }
                    if sample.signed_speed_m_s.abs() > 1e-9 {
                        moving += 1;
                        if t > 401.0 {
                            heated_moving += 1;
                        }
                        assert_eq!(
                            sample.regime,
                            if sample.signed_speed_m_s > 0.0 {
                                ContactLineRegime2::Advancing
                            } else {
                                ContactLineRegime2::Receding
                            }
                        );
                        assert!((sample.applied_angle_radians - theta).abs() > 1e-9);
                        assert_eq!(
                            (sample.applied_angle_radians - theta).is_sign_positive(),
                            sample.signed_speed_m_s.is_sign_positive()
                        );
                    }
                }
                assert!(
                    moving > 0,
                    "the real accepted VOF history must drive a moving line"
                );
                assert!(
                    heated_moving > 0,
                    "the moving line must use the actual heated adhesion source"
                );
                thermal.step(gate).unwrap();
                reference.step(gate).unwrap();
                let populations = thermal.flow().lattice().grid().state_bits();
                let reference_populations = reference.flow().lattice().grid().state_bits();
                assert_eq!(populations.len(), reference_populations.len());
                let difference = populations
                    .iter()
                    .zip(&reference_populations)
                    .map(|(&a, &b)| (f64::from_bits(a) - f64::from_bits(b)).abs())
                    .fold(0.0_f64, f64::max);
                assert!(
                    difference > 1e-12,
                    "rate-dependent angles must change native momentum on the next step: {difference}"
                );
            });
        }
    }
}

#[test]
fn g1_dynamic_equilibrium_uses_same_evolving_tension_as_native_pressure() {
    let liquid = varying_tension(0.00002, 500.0);
    let solid = body::wall_card(0.5, 500.0);
    // W(T) = sigma(T) isolates the tension in the rate response from a
    // changing Young angle. The planar interface translates along the wall.
    let interface = equilibrium_interface(
        &liquid,
        &solid,
        EquilibriumAngleOrigin::YoungDupre,
        0.00002,
        500.0,
    );
    let initial = uniform_cells()
        .into_iter()
        .map(|mut cell| {
            cell.velocity_m_s = [0.0, 0.0, 0.00002];
            cell
        })
        .collect::<Vec<_>>();
    for molecular in [false, true] {
        with_uniform(&liquid, &interface, &initial, 200000.0, |initial, gate| {
            let mut thermal = initial
                .clone()
                .with_dynamic_equilibrium_wetting(rate_law(molecular, 500.0, gate), gate)
                .unwrap();
            let heat = homogeneous_heat(&thermal);
            thermal.step_with_heat(&heat, gate).unwrap();
            let position = [0, 2, 3];
            let state = thermal.material_at(position).unwrap();
            let sigma = state
                .material()
                .property(SURFACE_TENSION_PROPERTY)
                .unwrap()
                .value_si();
            assert!(sigma > 0.004);
            // Here dx = dt = 0.001 and rho = 1000, so native sigma equals SI sigma.
            near(
                thermal
                    .flow()
                    .lattice()
                    .surface_tension_at(0, 2, 3)
                    .unwrap(),
                sigma,
                2e-14,
            );
            let sample = &thermal.wetting_samples().unwrap()[&position];
            assert!(sample.signed_speed_m_s.abs() > 1e-9);
            let t = temperature(state);
            assert_eq!(sample.material.state_coordinate("T"), Some(t));
            let work = sample
                .material
                .property(EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY)
                .unwrap()
                .value_si();
            near(work, sigma, 2e-14);
            if molecular {
                let prediction = sample.molecular_kinetic.unwrap();
                near(prediction.surface_tension_n_m, sigma, 2e-14);
                near(prediction.temperature_k, t, 2e-14);
                let argument = -sigma * 5e-10_f64.powi(2) * sample.applied_angle_radians.cos()
                    / (2.0 * 1.380_649e-23 * t);
                residual(
                    2.0 * 1e9 * 5e-10 * argument.sinh(),
                    sample.signed_speed_m_s,
                    1e-13,
                );
            } else {
                let prediction = sample.hydrodynamic.unwrap();
                near(prediction.surface_tension_n_m, sigma, 2e-14);
                let capillary = 0.1 * sample.signed_speed_m_s / sigma;
                near(prediction.signed_capillary_number, capillary, 2e-14);
                residual(
                    cox_difference(core::f64::consts::FRAC_PI_2, sample.applied_angle_radians),
                    capillary * 100000.0_f64.ln(),
                    3e-11,
                );
            }
        });
    }
}

#[test]
fn g3_g4_dynamic_equilibrium_retains_source_and_atomic_refusals() {
    for (source_limit, calibration_limit) in [(410.0, 500.0), (500.0, 410.0)] {
        with_source(
            EquilibriumAngleOrigin::YoungDupre,
            source_limit,
            |thermal, gate| {
                let mut thermal = thermal
                    .with_dynamic_equilibrium_wetting(rate_law(true, calibration_limit, gate), gate)
                    .unwrap();
                let before = (
                    thermal.flow().lattice().grid().state_bits(),
                    thermal.wetting_samples().unwrap().clone(),
                    thermal.total_enthalpy_j().unwrap(),
                    thermal.physical_time_s(),
                    thermal.material_at(WALL_CELL).unwrap().clone(),
                    cells()
                        .iter()
                        .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                        .collect::<Vec<_>>(),
                );
                assert!(
                    thermal
                        .clone()
                        .with_dynamic_equilibrium_wetting(rate_law(false, 500.0, gate), gate,)
                        .is_err()
                );
                let error = thermal.step_with_heat(&heat(), gate).unwrap_err();
                if source_limit == 410.0 {
                    assert!(matches!(
                        error,
                        FreeSurfaceEnthalpyError::Wetting(FreeSurfaceWettingError::Material(_))
                    ));
                } else {
                    assert!(matches!(
                        error,
                        FreeSurfaceEnthalpyError::Wetting(
                            FreeSurfaceWettingError::MolecularKinetic(
                                MolecularKineticError::OutsideDomain("temperature")
                            )
                        )
                    ));
                }
                let cancelled = CancelGate::new_clock_free();
                cancelled.request();
                assert!(thermal.step(&cancelled).is_err());
                assert_eq!(
                    (
                        thermal.flow().lattice().grid().state_bits(),
                        thermal.wetting_samples().unwrap().clone(),
                        thermal.total_enthalpy_j().unwrap(),
                        thermal.physical_time_s(),
                        thermal.material_at(WALL_CELL).unwrap().clone(),
                        cells()
                            .iter()
                            .map(|cell| thermal.specific_enthalpy_at(cell.coordinates))
                            .collect::<Vec<_>>()
                    ),
                    before
                );
                thermal.step(gate).unwrap();
            },
        );
    }
}
