//! Native heat-driven volume evolution, with independently calculated face inertia.
#[path = "heat/bath.rs"]
mod bath;

#[path = "heat/reservoir.rs"]
mod reservoir;

use super::*;
use fs_lbm::material::free_surface::enthalpy::volume::{
    FreeSurfaceSlabHeatControls, FreeSurfaceSlabHeatError,
};

fn heat_controls(bracket: [f64; 2]) -> FreeSurfaceSlabHeatControls {
    FreeSurfaceSlabHeatControls {
        volume: controls(),
        net_specific_enthalpy_bracket_j_kg: bracket,
        maximum_heat_defect_j: 5e-11,
        maximum_pressure_work_defect_j: 1e-16,
        maximum_source_evaluations: 80,
    }
}

// Each internal normal face owns its two ORIGINAL fluid half cells;
// the opening owns only its final liquid half cell. No gas mass is added.
fn inertia(height: f64, dx: f64, cell_mass: f64, columns: f64) -> f64 {
    let cells = (height / dx).ceil() as usize;
    let terminal = height / dx - (cells - 1) as f64;
    let mut b = 0.0;
    for face in 1..=cells {
        let x = if face == cells {
            height
        } else {
            face as f64 * dx
        };
        let mass = if face == cells {
            0.5 * cell_mass * terminal
        } else if face + 1 == cells {
            0.5 * cell_mass * (1.0 + terminal)
        } else {
            cell_mass
        };
        b += mass * (x / height).powi(2);
    }
    columns * b
}

#[test]
fn g1_heat_pays_original_pressure_work_on_all_slab_orientations_and_both_density_slopes() {
    for face in Face3::ALL {
        for hot_density in [980.0, 1020.0] {
            for datum in [0.0, 1e8] {
                with_case(
                    Case {
                        face,
                        hot_density,
                        datum,
                        increment_j_kg: 200_000.0,
                        ..Case::default()
                    },
                    |old, _, _, gate, cx| {
                        let bits = old.flow().lattice().grid().state_bits();
                        let cell_mass = old.flow().scaling().mass_unit_kg();
                        let mass = 128.0 * cell_mass;
                        let b = mass * (1.0 / 3.0 + 1.0 / 24.0);
                        let rho1 = 0.5 * (1000.0 + hot_density);
                        let speed = ((1000.0 - rho1) / rho1) * 0.002 / 0.001;
                        let work = b * speed * speed;
                        let supplied = mass * 200_000.0 + work;
                        let plan = old
                            .prepare_uniform_slab_heat(
                                face,
                                supplied,
                                heat_controls([0.0, 400_000.0]),
                                gate,
                                cx,
                            )
                            .unwrap();
                        let stepped = plan.advance(gate, cx).unwrap();
                        let r = stepped.report;
                        near(r.net_specific_enthalpy_increment_j_kg, 200_000.0, 1e-11);
                        near(r.surface_speed_inertia_kg, b, 2e-13);
                        near(r.pressure_correction_work_j, work, 2e-8);
                        near(plan.endpoint().state().temperature_k(), 400.0, 1e-12);
                        assert!(r.pressure_correction_work_j > 0.0);
                        assert!(r.heat_balance_defect_j.abs() <= 5e-11);
                        assert!(r.pressure_response_defect_j.abs() <= 1e-16);
                        near(stepped.volume.report.mass_after_kg, mass, 2e-13);
                        near(
                            stepped.volume.report.final_height_m,
                            0.002 * 1000.0 / rho1,
                            2e-12,
                        );
                        let mut location = [3; 3];
                        location[face as usize / 2] = if face as usize % 2 == 0 { 0 } else { 7 };
                        near(
                            stepped
                                .volume
                                .liquid
                                .specific_enthalpy_at(location)
                                .unwrap(),
                            datum + 200_000.0,
                            1e-11,
                        );
                        assert_eq!(stepped.volume.liquid.accepted_steps(), 1);
                        assert_eq!(old.flow().lattice().grid().state_bits(), bits);
                        assert_eq!(old.accepted_steps(), 0);
                    },
                );
            }
        }
    }
}

#[test]
fn g3_partial_and_subcell_heat_cooling_sequences_keep_actual_incoming_motion_and_clock() {
    for (depth_cells, terminal_fill) in [(3, 0.4), (1, 0.3)] {
        for face in [Face3::XMin, Face3::ZMax] {
            with_case(
                Case {
                    face,
                    depth_cells,
                    terminal_fill,
                    increment_j_kg: 200_000.0,
                    ..Case::default()
                },
                |old, _, _, gate, cx| {
                    let dx = 0.001;
                    let dt = 0.001;
                    let height0 = ((depth_cells - 1) as f64 + terminal_fill) * dx;
                    let m0 = old.flow().scaling().mass_unit_kg();
                    let mass = 64.0 * ((depth_cells - 1) as f64 + terminal_fill) * m0;
                    let b0 = inertia(height0, dx, m0, 64.0);
                    let u0 = ((1000.0 - 990.0) / 990.0) * height0 / dt;
                    let heat = mass * 200_000.0 + b0 * u0 * u0;
                    let hot_plan = old
                        .prepare_uniform_slab_heat(
                            face,
                            heat,
                            heat_controls([0.0, 400_000.0]),
                            gate,
                            cx,
                        )
                        .unwrap();
                    let hot = hot_plan.advance(gate, cx).unwrap();
                    let height1 = height0 * 1000.0 / 990.0;
                    let b1 = inertia(height1, dx, m0 * 990.0 / 1000.0, 64.0);
                    let u1 = ((990.0 - 995.0) / 995.0) * height1 / dt;
                    let cooling = mass * (-100_000.0) + b1 * u1 * (u1 - u0);
                    let cold_plan = hot
                        .volume
                        .liquid
                        .prepare_uniform_slab_heat(
                            face,
                            cooling,
                            heat_controls([-150_000.0, 0.0]),
                            gate,
                            cx,
                        )
                        .unwrap();
                    let cold = cold_plan.advance(gate, cx).unwrap();
                    near(
                        cold.report.net_specific_enthalpy_increment_j_kg,
                        -100_000.0,
                        1e-10,
                    );
                    near(
                        cold.report.pressure_correction_work_j,
                        b1 * u1 * (u1 - u0),
                        2e-7,
                    );
                    near(cold.volume.report.mass_after_kg, mass, 2e-13);
                    near(
                        cold.volume.report.final_height_m,
                        height0 * 1000.0 / 995.0,
                        2e-12,
                    );
                    assert!(cold.report.heat_balance_defect_j.abs() <= 5e-11);
                    assert_eq!(cold.volume.liquid.accepted_steps(), 2);
                    near(cold.volume.liquid.physical_time_s(), 0.252, 1e-14);
                    // Removing affine divergence at zero source spends the
                    // original numerical projection loss; it cannot heat h.
                    let zero_plan = cold
                        .volume
                        .liquid
                        .prepare_uniform_slab_heat(
                            face,
                            0.0,
                            heat_controls([-10_000.0, 10_000.0]),
                            gate,
                            cx,
                        )
                        .unwrap();
                    let zero = zero_plan.advance(gate, cx).unwrap();
                    assert_eq!(zero.report.net_specific_enthalpy_increment_j_kg, 0.0);
                    assert_eq!(zero.report.pressure_correction_work_j, 0.0);
                    assert!(zero.volume.projection.projection.report.projection_loss_j > 0.0);
                    assert_eq!(zero.volume.liquid.accepted_steps(), 3);
                    near(
                        zero.volume.report.final_height_m,
                        cold.volume.report.final_height_m,
                        1e-13,
                    );
                },
            );
        }
    }
}

#[test]
fn g0_bad_heat_bracket_exhaustion_and_native_accuracy_refuse_without_mutation() {
    with_case(
        Case {
            increment_j_kg: 200_000.0,
            ..Case::default()
        },
        |old, _, _, gate, cx| {
            let before = old.flow().lattice().grid().state_bits();
            let h = old.total_enthalpy_j().unwrap().to_bits();
            let q = old.flow().ledger_mass_kg().unwrap() * 200_000.0;
            assert!(matches!(
                old.prepare_uniform_slab_heat(Face3::XMin, q, heat_controls([0.0, 1.0]), gate, cx),
                Err(FreeSurfaceSlabHeatError::Balance { .. })
            ));
            let mut exhausted = heat_controls([0.0, 400_000.0]);
            exhausted.maximum_source_evaluations = 2;
            assert!(matches!(
                old.prepare_uniform_slab_heat(Face3::XMin, q, exhausted, gate, cx),
                Err(FreeSurfaceSlabHeatError::Balance { evaluations: 2, .. })
            ));
            let mut strict = heat_controls([0.0, 400_000.0]);
            strict.volume.maximum_kinetic_mapping_change_j = 0.0;
            let plan = old
                .prepare_uniform_slab_heat(Face3::XMin, q, strict, gate, cx)
                .unwrap();
            assert!(plan.advance(gate, cx).is_err());
            gate.request();
            assert!(
                old.prepare_uniform_slab_heat(
                    Face3::XMin,
                    q,
                    heat_controls([0.0, 400_000.0]),
                    gate,
                    cx
                )
                .is_err()
            );
            assert_eq!(old.flow().lattice().grid().state_bits(), before);
            assert_eq!(old.total_enthalpy_j().unwrap().to_bits(), h);
            assert_eq!(old.accepted_steps(), 0);
        },
    );
}

#[test]
fn g1_constant_density_uses_all_supplied_heat_without_expansion_work() {
    with_case(
        Case {
            hot_density: 1000.0,
            increment_j_kg: 200_000.0,
            ..Case::default()
        },
        |old, _, _, gate, cx| {
            let q = old.flow().ledger_mass_kg().unwrap() * 200_000.0;
            let plan = old
                .prepare_uniform_slab_heat(
                    Face3::XMin,
                    q,
                    heat_controls([0.0, 400_000.0]),
                    gate,
                    cx,
                )
                .unwrap();
            let step = plan.advance(gate, cx).unwrap();
            assert_eq!(step.report.pressure_correction_work_j, 0.0);
            near(step.report.represented_enthalpy_increment_j, q, 1e-13);
            near(step.volume.report.final_height_m, 0.002, 1e-13);
            assert_eq!(step.volume.report.native_kinetic_energy_j, 0.0);
        },
    );
}
