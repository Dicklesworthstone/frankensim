//! Source-to-native-flow regressions. Authored transport tuples below are
//! synthetic test inputs, not measured material cards or turbulence validation.

use super::*;
use crate::d3q19::sparse::sgs::tests::{ACTIVE, seeded, settings};
use crate::material::plan_scaling_from_material;
use crate::sgs3::{SgsError, SgsModel};
use fs_exec::TilePool;
use std::cell::Cell;

struct Source {
    value: FluidTransportState,
    reads: Cell<usize>,
}

impl FluidTransport for Source {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        Ok(self.value)
    }
}

fn fluid(mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, 1.0, mu).unwrap()
}

#[test]
fn sparse_material_identity_and_heterogeneous_zero_sgs_reach_native_collision() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, 0.01, 0.01, 0.05).unwrap();
    let cold = Source {
        value: fluid(0.4),
        reads: Cell::new(0),
    };
    // The momentum adapter retains this transport response without imposing
    // an unrelated thermal tau bound on its small conductivity.
    let hot = Source {
        value: FluidTransportState::try_new(1000.0, 2000.0, 1e-250, 1.6).unwrap(),
        reads: Cell::new(0),
    };
    let gate = CancelGate::new_clock_free();
    let mut actual = seeded(&ACTIVE);
    let count = ACTIVE.len() * 64;
    let bound = scaling
        .bind_sparse_sgs_materials(&actual, count, &gate, |p| {
            Ok(if p[0] < 4 { &cold } else { &hot })
        })
        .unwrap();
    for sample in bound.samples() {
        let source = if sample.cell[0] < 4 { &cold } else { &hot };
        assert!(core::ptr::eq(sample.material, source));
        assert_eq!(sample.transport, source.value);
        let physical_nu = source.value.dynamic_viscosity_pa_s() / 1000.0;
        let expected = 0.5 + 3.0 * physical_nu * 0.01 / (0.01 * 0.01);
        assert!((sample.molecular_relaxation_time - expected).abs() < 2e-15);
    }
    assert_eq!(cold.reads.get() + hot.reads.get(), count);
    let mut expected = actual.clone();
    let mut wrong_uniform = actual.clone();
    let mut workspace = actual.sgs_workspace(&gate).unwrap();
    let allocation = workspace.collisions().as_ptr();
    for _ in 0..3 {
        expected
            .step_serial_inner(Some(&gate), Some(bound.molecular_collisions()))
            .unwrap();
        bound
            .step_serial(
                &mut actual,
                &gate,
                settings(SgsModel::Wale { cw: 0.0 }),
                &mut workspace,
            )
            .unwrap();
        wrong_uniform.step_serial().unwrap();
    }
    assert_eq!(actual.state_bits(), expected.state_bits());
    assert_ne!(actual.state_bits(), wrong_uniform.state_bits());
    assert_eq!(workspace.collisions(), bound.molecular_collisions());
    assert_eq!(workspace.collisions().as_ptr(), allocation);
    assert_eq!(
        cold.reads.get() + hot.reads.get(),
        count,
        "held material epoch was resampled"
    );
}

#[test]
fn sparse_local_molecular_viscosity_controls_eddy_addition_and_ratio_refusal() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, 0.01, 0.01, 0.05).unwrap();
    let cold = fluid(0.4);
    let hot = fluid(1.6);
    let gate = CancelGate::new_clock_free();
    let mut grid = seeded(&ACTIVE);
    grid.force = [0.0; 3];
    let bound = scaling
        .bind_sparse_sgs_materials(&grid, ACTIVE.len() * 64, &gate, |p| {
            Ok(if p[0] < 4 { &cold } else { &hot })
        })
        .unwrap();
    let mut workspace = grid.sgs_workspace(&gate).unwrap();
    let config = settings(SgsModel::Smagorinsky { cs: 0.18 });
    let first = bound
        .prepare_sgs(&grid, &gate, config, &mut workspace)
        .unwrap();
    // Both interior samples have the same known linear gradient, while their
    // source viscosities differ fourfold. The eddy addition stays the same.
    let expected_eddy =
        0.18 * 0.18 * (0.004_f64.powi(2) + 0.003_f64.powi(2) + 0.002_f64.powi(2)).sqrt();
    for cell in [[1, 1, 1], [5, 1, 1]] {
        let sample = bound.samples().iter().position(|s| s.cell == cell).unwrap();
        let nu = bound.molecular_collisions()[sample].kinematic_viscosity();
        let total = workspace.collisions()[sample].kinematic_viscosity();
        assert!((total - nu - expected_eddy).abs() < 2e-16);
    }
    let expected_ratio = workspace
        .collisions()
        .iter()
        .zip(bound.molecular_collisions())
        .map(|(total, molecular)| {
            (total.kinematic_viscosity() - molecular.kinematic_viscosity())
                / molecular.kinematic_viscosity()
        })
        .fold(0.0_f64, f64::max);
    assert!((first.maximum_eddy_to_molecular - expected_ratio).abs() < 2e-14);
    let prepared = workspace.collisions().to_vec();
    assert_eq!(
        bound
            .prepare_sgs(&grid, &gate, config, &mut workspace)
            .unwrap(),
        first
    );
    assert_eq!(
        workspace.collisions(),
        prepared,
        "SGS viscosity compounded between preparations"
    );
    let before = grid.state_bits();
    let mut limited = config;
    limited.maximum_eddy_to_molecular = 0.5 * first.maximum_eddy_to_molecular;
    assert!(matches!(
        bound.step_serial(&mut grid, &gate, limited, &mut workspace),
        Err(SparseMaterialSgsError3::Sgs(SparseSgsError3::Closure(
            SgsError::ViscosityRatioExceeded { .. }
        )))
    ));
    assert_eq!(grid.state_bits(), before);
    assert_eq!(grid.steps(), 0);
}

#[test]
fn sparse_material_pool_and_epoch_refusals_preserve_native_state() {
    let reference = fluid(1.0);
    let scaling = plan_scaling_from_material(&reference, 0.01, 0.01, 0.05).unwrap();
    let cold = fluid(0.4);
    let hot = fluid(1.6);
    let gate = CancelGate::new_clock_free();
    let tiles: Vec<_> = (0..3)
        .flat_map(|x| (0..2).flat_map(move |y| (0..2).map(move |z| (x, y, z))))
        .filter(|p| *p != (2, 1, 1))
        .collect();
    let initial = seeded(&tiles);
    let count = tiles.len() * 64;
    let bound = scaling
        .bind_sparse_sgs_materials(&initial, count, &gate, |p| {
            Ok(if p[0] < 4 { &cold } else { &hot })
        })
        .unwrap();
    let config = settings(SgsModel::Vreman { cv: 0.07 });
    let mut serial = initial.clone();
    let mut workspace = serial.sgs_workspace(&gate).unwrap();
    for _ in 0..3 {
        bound
            .step_serial(&mut serial, &gate, config, &mut workspace)
            .unwrap();
    }
    let pool = TilePool::for_host(3, 847);
    let mut pooled = initial.clone();
    for _ in 0..3 {
        bound
            .step_pooled(&mut pooled, &pool, &gate, config, &mut workspace)
            .unwrap();
    }
    assert_eq!(pooled.state_bits(), serial.state_bits());
    let reads = Cell::new(0);
    assert!(
        scaling
            .bind_sparse_sgs_materials(&initial, count - 1, &gate, |_| {
                reads.set(reads.get() + 1);
                Ok(&cold)
            })
            .is_err()
    );
    assert_eq!(reads.get(), 0, "budget refusal sampled a source");
    let bad_density = FluidTransportState::try_new(999.0, 1000.0, 1.0, 0.4).unwrap();
    assert!(matches!(
        scaling.bind_sparse_sgs_materials(&initial, count, &gate, |_| Ok(&bad_density)),
        Err(SparseMaterialSgsError3::Source { .. })
    ));
    let cancelled = CancelGate::new_clock_free();
    assert!(matches!(
        scaling.bind_sparse_sgs_materials(&initial, count, &cancelled, |_| {
            cancelled.request();
            Ok(&cold)
        }),
        Err(SparseMaterialSgsError3::Sgs(SparseSgsError3::Flow(
            SparseError3::Cancelled
        )))
    ));
    assert_eq!(initial.steps(), 0);
    let mut changed = initial.clone();
    changed.deactivate_tiles(&[(0, 0, 0)]).unwrap();
    changed.activate_tiles(&[(2, 1, 1)]).unwrap();
    let before = changed.state_bits();
    // A fresh SGS workspace cannot bless the stale material assignment.
    let mut fresh = changed.sgs_workspace(&gate).unwrap();
    assert!(matches!(
        bound.step_serial(&mut changed, &gate, config, &mut fresh),
        Err(SparseMaterialSgsError3::Sgs(
            SparseSgsError3::TopologyChanged
        ))
    ));
    assert_eq!(changed.state_bits(), before);
}
