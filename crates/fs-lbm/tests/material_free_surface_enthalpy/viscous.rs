//! Molecular stress work reaches the same liquid h, viscosity and buoyancy owners.

use super::*;
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;

fn with_liquid<R>(
    card: &MaterialCard,
    shear: f64,
    buoyant: bool,
    initial_h: impl Fn([usize; 3]) -> f64,
    exercise: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        card,
        &point(300.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    let mut cells = initial(true);
    for cell in &mut cells {
        cell.velocity_m_s = [shear * (cell.coordinates[2] as f64 - 3.5), 0.0, 0.0];
    }
    with_cx(&gate, |cx| {
        let gravity = if buoyant { [0.0, 0.0, -0.1] } else { [0.0; 3] };
        let flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(
                MaterialFreeSurfaceDomain3 {
                    gravity_m_s2: gravity,
                    ..domain()
                },
                &cells,
                cx,
            )
            .unwrap();
        let mut liquid = if buoyant {
            flow.with_boussinesq_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                BoussinesqConfig {
                    gravity_m_s2: gravity,
                    max_relative_density_departure: 0.02,
                    max_relative_specific_heat_departure: 0.0,
                },
                initial_h,
                &gate,
            )
        } else {
            flow.with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                initial_h,
                &gate,
            )
        }
        .unwrap();
        exercise(&mut liquid, &gate)
    })
}

#[test]
fn g1_molecular_stress_heats_finite_liquid_and_refreshes_source_viscosity() {
    let source = card(0.02, 1000.0);
    with_liquid(
        &source,
        0.01,
        false,
        |_| 200000.0,
        |liquid, gate| {
            // Local equilibria do not yet carry resolved viscous stress.
            // Native collision/stream creates the non-equilibrium shear field.
            liquid.step(gate).unwrap();
            let mu_before = liquid
                .material_at([3, 3, 3])
                .unwrap()
                .dynamic_viscosity_pa_s();
            let report = liquid.step_with_viscous_heating(&[], gate).unwrap();
            assert!(report.flow.viscous_heat_j > 1e-12);
            assert_eq!(report.flow.external_heat_j, 0.0);
            let independently_summed: f64 = report
                .applied
                .iter()
                .map(|cell| {
                    assert_eq!(
                        cell.material.material().card_identity(),
                        source.content_hash()
                    );
                    assert!(cell.liquid_volume_m3 > 0.0);
                    cell.material.dynamic_viscosity_pa_s()
                        * cell.shear_rate_per_s.powi(2)
                        * cell.liquid_volume_m3
                        * 0.001
                })
                .sum();
            near(report.flow.viscous_heat_j, independently_summed, 1e-12);
            near(
                report
                    .applied
                    .iter()
                    .map(|cell| cell.liquid_volume_m3)
                    .sum::<f64>(),
                report.flow.mass_before_kg / 1000.0,
                1e-12,
            );
            assert!(
                (report.flow.enthalpy_after_j
                    - report.flow.enthalpy_before_j
                    - report.flow.viscous_heat_j)
                    .abs()
                    < 1e-9
            );
            assert!(report.flow.enthalpy_after_j > report.flow.enthalpy_before_j);
            assert!(
                liquid
                    .material_at([3, 3, 3])
                    .unwrap()
                    .dynamic_viscosity_pa_s()
                    < mu_before
            );
            near(report.flow.mass_after_kg, report.flow.mass_before_kg, 1e-12);
        },
    );
}

#[test]
fn g3_rest_adds_no_resolved_viscous_heat_or_change_to_original_flow() {
    let source = card(0.02, 1000.0);
    let ordinary = with_liquid(
        &source,
        0.0,
        false,
        |_| 200000.0,
        |liquid, gate| {
            liquid.step(gate).unwrap();
            (
                liquid.flow().lattice().grid().state_bits(),
                liquid.total_enthalpy_j().unwrap(),
            )
        },
    );
    let heated = with_liquid(
        &source,
        0.0,
        false,
        |_| 200000.0,
        |liquid, gate| {
            let report = liquid.step_with_viscous_heating(&[], gate).unwrap();
            // Floating reduction of equilibrium may leave roundoff stress, which
            // must not become resolved heat or alter the native trajectory.
            assert!(report.flow.viscous_heat_j < 1e-25);
            (
                liquid.flow().lattice().grid().state_bits(),
                liquid.total_enthalpy_j().unwrap(),
            )
        },
    );
    assert_eq!(ordinary, heated);
}

#[test]
fn g4_source_refusal_and_cancellation_preserve_native_motion_and_enthalpy() {
    let source = card(0.02, 1000.0);
    with_liquid(
        &source,
        0.01,
        false,
        |_| 200000.0,
        |liquid, gate| {
            liquid.step(gate).unwrap();
            let bits = liquid.flow().lattice().grid().state_bits();
            let enthalpy = liquid.total_enthalpy_j().unwrap();
            let material = liquid.material_at([3, 3, 3]).unwrap().clone();
            assert!(
                liquid
                    .step_with_viscous_heating(
                        &[FreeSurfaceCellHeat {
                            coordinates: [3, 3, 3],
                            energy_j: 1000.0
                        }],
                        gate,
                    )
                    .is_err()
            );
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), enthalpy);
            assert_eq!(liquid.material_at([3, 3, 3]).unwrap(), &material);
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(liquid.step_with_viscous_heating(&[], &cancelled).is_err());
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), enthalpy);
            assert_eq!(liquid.material_at([3, 3, 3]).unwrap(), &material);
            assert!(
                liquid
                    .step_with_viscous_heating(&[], gate)
                    .unwrap()
                    .flow
                    .viscous_heat_j
                    > 0.0
            );
        },
    );
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn g1_conduction_buoyancy_and_molecular_heat_share_one_liquid_interval() {
    use fs_lbm::material::free_surface::enthalpy::conduction::{
        FreeSurfaceConductionControls, FreeSurfaceConductionFace,
    };
    let source = card(0.02, 990.0);
    with_liquid(
        &source,
        0.01,
        true,
        |p| 200000.0 + 1000.0 * p[0] as f64,
        |liquid, gate| {
            liquid.step(gate).unwrap();
            let report = liquid
                .step_with_conduction_and_viscous_heating(
                    &[FreeSurfaceConductionFace {
                        cells: [[3, 3, 3], [4, 3, 3]],
                        area_m2: 1e-6,
                        distances_m: [0.0005; 2],
                    }],
                    &[FreeSurfaceCellHeat {
                        coordinates: [3, 3, 3],
                        energy_j: 1e-6,
                    }],
                    FreeSurfaceConductionControls {
                        maximum_faces: 1,
                        maximum_substeps: 20,
                        maximum_dt_s: 0.001,
                        energy_tolerance_j: 1e-9,
                    },
                    gate,
                )
                .unwrap();
            let flow = report.conduction.flow;
            assert!(flow.viscous_heat_j > 0.0);
            assert_eq!(flow.external_heat_j, 1e-6);
            assert!(
                (report.conduction.thermal.external_heat_j
                    - flow.external_heat_j
                    - flow.viscous_heat_j)
                    .abs()
                    < 1e-9
            );
            assert!(
                (flow.enthalpy_after_j
                    - flow.enthalpy_before_j
                    - flow.external_heat_j
                    - flow.viscous_heat_j)
                    .abs()
                    < 1e-9
            );
            assert!(report.conduction.last_applied[0].conductance_w_k > 0.0);
            let force = liquid.flow().lattice().force_density_at(3, 3, 3).unwrap();
            assert_ne!(force[2], -0.0001);
            near(flow.mass_after_kg, flow.mass_before_kg, 1e-12);
        },
    );
}
