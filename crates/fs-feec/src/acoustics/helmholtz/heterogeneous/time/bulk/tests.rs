//! Original FE volume moments, complete histories and work-conjugate sources.
use super::*;
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use crate::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_rep_mesh::TetComplex;

const X: [[f64; 3]; 5] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, -1.0],
];
const MEDIUM: Medium = Medium {
    density_kg_m3: 2.0,
    sound_speed_m_s: 3.0,
};
fn model<'f, 'c>(g: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        g,
        g.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}
fn arm(strength: f64, time: f64) -> BulkMemoryArm {
    BulkMemoryArm {
        delta_compressibility_pa_inv: strength,
        relaxation_time_s: time,
    }
}
fn moment(a: &[f64], b: &[f64]) -> f64 {
    // Integral of two independent affine fields on either V=1/6 tetrahedron.
    (a.iter().sum::<f64>() * b.iter().sum::<f64>()
        + a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>())
        / 120.0
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:e} != {b:e}, tolerance {tol:e}");
}

#[test]
fn original_consistent_volume_forms_keep_independent_interface_memories() {
    let mesh = TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]);
    let g = HelmholtzFe::new(&mesh, &X, 1, MEDIUM).unwrap();
    let m = model(&g);
    let arms = [arm(0.03, 0.2), arm(0.07, 0.5)];
    let regions = [
        BulkMemoryRegion {
            tetrahedra: vec![0],
            arms: vec![arms[0]],
        },
        BulkMemoryRegion {
            tetrahedra: vec![1],
            arms: vec![arms[1]],
        },
    ];
    let h = 0.01;
    let owner = m
        .prepare_bulk_time(
            h,
            &regions,
            &[],
            &[],
            &[],
            BulkTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(owner.arm_dofs(0), Some([0, 1, 2, 3].as_slice()));
    assert_eq!(owner.arm_dofs(1), Some([0, 1, 2, 4].as_slice()));
    let p = [1.0, -0.2, 0.7, 1.8, -0.4];
    let histories = vec![vec![0.3, -0.1, 0.6, 1.2], vec![-0.8, 0.7, -0.2, 0.4]];
    let initial = owner
        .state(&[0.0; 5], &p, &histories, &mut || false)
        .unwrap();
    let energy = owner.energy(&initial, &mut || false).unwrap();
    let beta = 1.0 / 18.0;
    close(
        energy.frozen_compression_j,
        0.5 * beta
            * (moment(&p[..4], &p[..4])
                + moment(&[p[0], p[1], p[2], p[4]], &[p[0], p[1], p[2], p[4]])),
        2e-16,
    );
    close(
        energy.molecular_storage_j,
        0.5 * (arms[0].delta_compressibility_pa_inv * moment(&histories[0], &histories[0])
            + arms[1].delta_compressibility_pa_inv * moment(&histories[1], &histories[1])),
        2e-16,
    );
    let step = owner
        .step_projected(&initial, &[0.0; 5], &mut || false)
        .unwrap();
    let mut heat = 0.0;
    for (i, a) in arms.iter().enumerate() {
        let q1 = step.state.memory_pressure_pa(i).unwrap();
        let mismatch: Vec<_> = owner
            .arm_dofs(i)
            .unwrap()
            .iter()
            .enumerate()
            .map(|(j, &d)| 0.5 * (p[d] + step.state.pressure_pa()[d] - histories[i][j] - q1[j]))
            .collect();
        let expected =
            h * a.delta_compressibility_pa_inv / a.relaxation_time_s * moment(&mismatch, &mismatch);
        close(step.memories[i].dissipation_j, expected, 2e-16);
        heat += expected;
    }
    close(step.work.molecular_dissipation_j, heat, 2e-16);
    close(
        step.work.final_energy.total_j() + heat,
        energy.total_j(),
        2e-15,
    );
    assert!(owner.acoustic_state(&step.state).unwrap().pressure_pa() == step.state.pressure_pa());
    let duplicate = [regions[0].clone(), regions[0].clone()];
    assert!(
        m.prepare_bulk_time(
            h,
            &duplicate,
            &[],
            &[],
            &[],
            BulkTimeLimits::default(),
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn high_order_pressure_release_preserves_all_volume_memory_coordinates() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let g = HelmholtzFe::new(&mesh, &X[..4], 2, MEDIUM).unwrap();
    let m = model(&g);
    let region = BulkMemoryRegion {
        tetrahedra: vec![0],
        arms: vec![arm(0.04, 0.02)],
    };
    let h = 0.001;
    let face = [0, 1, 2];
    let owner = m
        .prepare_bulk_time(
            h,
            &[region],
            &[],
            &[],
            &[face],
            BulkTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(owner.arm_dofs(0).unwrap().len(), 10);
    let fixed = g
        .pressure_release_constraints(&[face], EmptyFreePressure::Refuse, &mut || false)
        .unwrap();
    assert_eq!(fixed.constrained.iter().filter(|&&v| v).count(), 6);
    let q: Vec<_> = (0..10).map(|j| 0.2 * (j as f64 + 1.0)).collect();
    let initial = owner
        .state(&[0.0; 10], &[0.0; 10], &[q.clone()], &mut || false)
        .unwrap();
    let step = owner
        .step_projected(&initial, &[0.0; 10], &mut || false)
        .unwrap();
    let decay = (0.02 - h / 2.0) / (0.02 + h / 2.0);
    for (j, &d) in owner.arm_dofs(0).unwrap().iter().enumerate() {
        if fixed.constrained[d] {
            assert_eq!(step.state.pressure_pa()[d], 0.0);
            assert_eq!(step.state.pressure_potential_pa_s()[d], 0.0);
            close(
                step.state.memory_pressure_pa(0).unwrap()[j],
                decay * q[j],
                2e-15,
            );
        }
    }
    assert!(step.work.molecular_dissipation_j > 0.0);
    assert!(
        owner
            .step_with_sources(&initial, &[(face, [0.0; 3])], &[], &[], &mut || false)
            .is_err()
    );
}

#[test]
fn original_prepared_volume_and_boundary_sources_drive_the_same_physical_work() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let g = HelmholtzFe::new(&mesh, &X[..4], 1, MEDIUM).unwrap();
    let m = model(&g);
    let boundary = m
        .quadratic_velocity_source(&[([0, 1, 2], [2.0; 6])], &mut || false)
        .unwrap();
    let volume = m
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    let region = BulkMemoryRegion {
        tetrahedra: vec![0],
        arms: vec![arm(0.03, 0.1)],
    };
    let h = 0.002;
    let owner = m
        .prepare_bulk_time(
            h,
            &[region],
            &[],
            &[],
            &[],
            BulkTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = owner
        .state(
            &[0.0; 4],
            &[1.0, -0.2, 0.6, 0.4],
            &[vec![0.0; 4]],
            &mut || false,
        )
        .unwrap();
    let result = owner
        .step_with_sources(
            &initial,
            &[([0, 1, 2], [0.1; 3])],
            &[(&boundary, 2.0)],
            &[(&volume, 3.0)],
            &mut || false,
        )
        .unwrap();
    let p: Vec<_> = initial
        .pressure_pa()
        .iter()
        .zip(result.state.pressure_pa())
        .map(|(a, b)| 0.5 * (a + b))
        .collect();
    let boundary_work = h * 2.0 * p[..3].iter().sum::<f64>() / 3.0;
    let volume_work = h * 3.0 * p.iter().sum::<f64>() / 4.0;
    let linear_work = h * p[..3].iter().sum::<f64>() / 60.0;
    close(result.work.boundary_source_work_j, boundary_work, 2e-16);
    close(result.work.volume_source_work_j, volume_work, 2e-16);
    close(
        result.work.input_work_j,
        boundary_work + volume_work + linear_work,
        2e-16,
    );
    let other = model(&g);
    let wrong = other
        .volume_velocity_source(3, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    assert!(
        owner
            .step_with_sources(&initial, &[], &[], &[(&wrong, 1.0)], &mut || false)
            .is_err()
    );
}
