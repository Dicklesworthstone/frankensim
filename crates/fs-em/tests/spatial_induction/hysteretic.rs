//! Actual spatial vector memory, charge continuity and separated physical heat.
#[path = "hysteretic/common.rs"]
mod common;

use super::{cross, grid, near};
use common::*;
use fs_em::spatial::{
    SpatialMagneticError,
    vector_potential::induction::hysteretic::{
        SpatialHystereticInduction, SpatialHystereticInductionError,
    },
};
use fs_material::field_storage::hysteresis::vector::VectorHysteresisError;

fn assert_mode(state: &SpatialHystereticInduction<'_>, h: f64, z: f64, b: f64) {
    let diagonal = state.problem().mesh().edges.binary_search(&[0, 7]).unwrap();
    near(state.potential_wb()[diagonal], 0.5 * b, 3e-9);
    let field = state.field(&mut || false).unwrap();
    near(field.stored_energy_j(), stored(h, z), 3e-8);
    for (cell, direction) in DIRECTIONS.iter().enumerate() {
        let memory = &state.cell_states()[cell];
        assert_eq!(memory.polarizations_t().len(), 1);
        for (axis, direction) in direction.iter().enumerate() {
            near(memory.field_a_m()[axis], h * direction, 3e-8);
            near(memory.flux_density_t()[axis], b * direction, 3e-9);
            near(memory.polarizations_t()[0][axis], z * direction, 3e-9);
            near(field.field_strength_a_m()[cell][axis], h * direction, 3e-8);
            near(field.flux_density_t()[cell][axis], b * direction, 3e-9);
        }
        near(memory.stored_energy_j_m3(), stored(h, z), 3e-8);
    }
}

#[test]
fn one_mode_play_closes_charge_and_energy_and_retains_remanence_after_switch_off() {
    let card = card(500.0);
    let (positions, tets) = grid(1);
    let problem = problem(&card, positions, tets);
    let mut state = initial(&problem, Some(conductor(&card, 300.0)));
    let (mut old_h, mut old_z, mut old_b) = (0.0, 0.0, 0.0);
    let sigma = resistivity(300.0).recip();
    let mass = 3.0 * sigma / 40.0;
    let mut minor_return = None;
    // Real current-driven preparation, saturation, reversal and a closed minor
    // loop precede source-free relaxation. No detached initial magnetic state.
    for index in 0..21 {
        let (h, z, b, current) = if index < 5 {
            let h = [12.0, 0.0, 6.0, 1.0, 6.0][index];
            let z = play(old_z, h);
            let b = MU * h + z;
            (h, z, b, 4.0 * h + mass * (b - old_b) / DT)
        } else {
            let (h, z, b) = endpoint(old_b, old_z, sigma, DT, 0.0);
            (h, z, b, 0.0)
        };
        let candidate = state
            .step_candidate(&[[current; 3]; 6], DT, field_controls(), &mut || false)
            .unwrap();
        assert_mode(candidate.state(), h, z, b);
        let step = candidate.step();
        let joule = mass * (0.5 * (b - old_b)).powi(2) / DT;
        let pinning = PIN * (z - old_z).abs();
        let change = stored(h, z) - stored(old_h, old_z);
        let numerical = h * (b - old_b) - change - pinning;
        near(step.joule_heat_j, joule, 3e-8);
        near(step.hysteresis_heat_j, pinning, 3e-8);
        near(step.heat_j, joule + pinning, 3e-8);
        near(step.stored_change_j, change, 3e-8);
        near(step.numerical_dissipation_j, numerical, 3e-8);
        near(step.impressed_work_j, 0.25 * current * (b - old_b), 3e-8);
        near(step.energy_residual_j, 0.0, 3e-8);
        assert!(step.maximum_charge_residual_a <= 2e-8);
        assert!(step.maximum_ampere_residual_a <= 2e-8);
        assert_eq!(step.cells.len(), 6);
        assert_eq!(step.magnetic_cells.len(), 6);
        for cell in &step.cells {
            let vertices = problem.mesh().tets[cell.cell];
            let projection = |v| match v {
                0 => 0.0,
                7 => 0.5,
                _ => 0.25,
            };
            for i in 1..4 {
                near(
                    cell.potential_v[i] - cell.potential_v[0],
                    -(b - old_b) / (2.0 * DT) * (projection(vertices[i]) - projection(vertices[0])),
                    3e-8,
                );
            }
            near(cell.heat_j, joule / 6.0, 3e-8);
        }
        for cell in &step.magnetic_cells {
            near(cell.irreversible_heat_j, pinning / 6.0, 3e-8);
            near(cell.energy_change_j, change / 6.0, 3e-8);
            near(cell.numerical_dissipation_j, numerical / 6.0, 3e-8);
        }
        if index == 0 {
            assert!(joule > 1.0 && pinning > 1.0 && numerical > 1.0);
        }
        if index == 2 {
            minor_return = Some(candidate.state().cell_states().to_vec());
        } else if index == 4 {
            for (returned, original) in candidate
                .state()
                .cell_states()
                .iter()
                .zip(minor_return.as_ref().unwrap())
            {
                for axis in 0..3 {
                    near(
                        returned.polarizations_t()[0][axis],
                        original.polarizations_t()[0][axis],
                        3e-9,
                    );
                    near(returned.field_a_m()[axis], original.field_a_m()[axis], 3e-8);
                }
            }
        }
        near(state.time_s(), index as f64 * DT, 1e-14);
        state = candidate.into_parts().0;
        (old_h, old_z, old_b) = (h, z, b);
    }
    assert!(old_h.abs() < 3e-6);
    near(old_b, CHI * PIN, 1e-6);
    near(old_z, CHI * PIN, 1e-6);
    assert_mode(&state, old_h, old_z, old_b);
}

fn rotate(v: [f64; 3]) -> [f64; 3] {
    // Proper rotation from quaternion (1,2,3,4)/sqrt(30); no axis is preserved.
    [
        (-10.0 * v[0] + 2.0 * v[1] + 11.0 * v[2]) / 15.0,
        (10.0 * v[0] - 5.0 * v[1] + 10.0 * v[2]) / 15.0,
        (5.0 * v[0] + 14.0 * v[1] + 2.0 * v[2]) / 15.0,
    ]
}
fn rotated_near(actual: [f64; 3], original: [f64; 3]) {
    for (a, b) in actual.into_iter().zip(rotate(original)) {
        near(a, b, 3e-7);
    }
}

#[test]
fn noncollinear_spatial_history_and_charge_are_covariant_under_rigid_rotation() {
    let card = card(500.0);
    let (positions, tets) = grid(2);
    let rotated = problem(
        &card,
        positions.iter().copied().map(rotate).collect(),
        tets.clone(),
    );
    let original = problem(&card, positions, tets);
    let mut a = initial(&original, Some(conductor(&card, 300.0)));
    let mut b = initial(&rotated, Some(conductor(&card, 300.0)));
    let count = original.mesh().tets.len();
    let mut turned = false;
    for current in [[50.0, 50.0, 50.0], [-40.0, 55.0, 15.0], [0.0; 3]] {
        let ca = a
            .step_candidate(&vec![current; count], DT, field_controls(), &mut || false)
            .unwrap();
        let cb = b
            .step_candidate(
                &vec![rotate(current); count],
                DT,
                field_controls(),
                &mut || false,
            )
            .unwrap();
        for ((old, next), rotated) in a
            .cell_states()
            .iter()
            .zip(ca.state().cell_states())
            .zip(cb.state().cell_states())
        {
            rotated_near(rotated.field_a_m(), next.field_a_m());
            rotated_near(rotated.flux_density_t(), next.flux_density_t());
            rotated_near(rotated.polarizations_t()[0], next.polarizations_t()[0]);
            let turn = cross(old.polarizations_t()[0], next.polarizations_t()[0]);
            turned |= turn.into_iter().map(|v| v * v).sum::<f64>() > 1e-7;
        }
        for (ea, eb) in ca
            .state()
            .potential_wb()
            .iter()
            .zip(cb.state().potential_wb())
        {
            near(*ea, *eb, 3e-8);
        }
        for (la, lb) in ca.step().cells.iter().zip(&cb.step().cells) {
            assert_eq!(la.cell, lb.cell);
            near(la.heat_j, lb.heat_j, 3e-8);
            for i in 0..4 {
                near(la.potential_v[i], lb.potential_v[i], 3e-8);
                rotated_near(lb.current_density_a_m2[i], la.current_density_a_m2[i]);
            }
        }
        near(ca.step().stored_after_j, cb.step().stored_after_j, 3e-8);
        near(
            ca.step().hysteresis_heat_j,
            cb.step().hysteresis_heat_j,
            3e-8,
        );
        near(
            ca.step().numerical_dissipation_j,
            cb.step().numerical_dissipation_j,
            3e-8,
        );
        near(ca.step().impressed_work_j, cb.step().impressed_work_j, 3e-8);
        for step in [ca.step(), cb.step()] {
            assert!(step.maximum_charge_residual_a <= 2e-8);
            near(step.energy_residual_j, 0.0, 3e-8);
        }
        a = ca.into_parts().0;
        b = cb.into_parts().0;
    }
    assert!(
        turned,
        "fixture must rotate retained polarization, not just rescale it"
    );
}

#[test]
fn late_cancellation_and_out_of_fit_excitation_preserve_every_accepted_history() {
    let card = card(500.0);
    let (positions, tets) = grid(1);
    let problem = problem(&card, positions, tets);
    let state = initial(&problem, Some(conductor(&card, 300.0)))
        .step_candidate(&[[63.12; 3]; 6], DT, field_controls(), &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let before_a = state.potential_wb().to_vec();
    let before_memory = state.cell_states().to_vec();
    let before_current = state.impressed_current_density_a_m2().to_vec();
    let mut checkpoints = 0;
    let reference = state
        .step_candidate(&[[0.0; 3]; 6], DT, field_controls(), &mut || {
            checkpoints += 1;
            false
        })
        .unwrap();
    assert!(checkpoints > 10);
    let mut reached = 0;
    let cancelled = state.step_candidate(&[[0.0; 3]; 6], DT, field_controls(), &mut || {
        reached += 1;
        reached == checkpoints - 1
    });
    assert!(matches!(
        cancelled,
        Err(SpatialHystereticInductionError::Spatial(
            SpatialMagneticError::Cancelled
        )) | Err(SpatialHystereticInductionError::Material(
            VectorHysteresisError::Cancelled
        ))
    ));
    assert!(
        state
            .step_candidate(&[[1e6; 3]; 6], DT, field_controls(), &mut || false)
            .is_err()
    );
    assert_eq!(state.time_s(), DT);
    assert_eq!(state.potential_wb(), before_a);
    assert_eq!(state.cell_states(), before_memory);
    assert_eq!(state.impressed_current_density_a_m2(), before_current);
    let retry = state
        .step_candidate(&[[0.0; 3]; 6], DT, field_controls(), &mut || false)
        .unwrap();
    assert_eq!(
        retry.state().potential_wb(),
        reference.state().potential_wb()
    );
    assert_eq!(retry.state().cell_states(), reference.state().cell_states());
    assert_eq!(retry.step().heat_j, reference.step().heat_j);
}
