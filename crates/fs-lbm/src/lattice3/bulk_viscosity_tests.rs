use super::{D3q19, GridBudget3, Lattice3, StepData, ViscousLattice3, coordinates, macro_cell};
use crate::d3q19::{CollisionModel3, Q3, equilibrium3};
use core::marker::PhantomData;

// Exercise the SAME collision and pull-streaming passes as the native serial
// grid, including its 64-cell partition and physical macroscopic admission.
// The fixture owns only initial populations, buffers and Fourier measurements.
fn decay(nx: usize, axis: usize, bulk_tau: f64, longitudinal: bool) -> f64 {
    let mut dims = [2; 3];
    dims[axis] = nx;
    let cells = dims.iter().product();
    let solids = vec![false; cells];
    let data = StepData::<D3q19, Q3> {
        dims,
        periodic: [true; 3],
        solid: &solids,
        collision: CollisionModel3::CentralMomentBulk {
            shear_rate: 1.0 / 0.68,
            bulk_rate: 1.0 / bulk_tau,
            higher_order_rate: 1.0,
        },
        collisions: None,
        force: [0.0; 3],
        forces: None,
        budget: GridBudget3 {
            maximum_cells: cells,
            maximum_mach: 0.2,
        },
        marker: PhantomData,
    };
    let k = core::f64::consts::TAU / nx as f64;
    let cs = D3q19::CS2.sqrt();
    let mut populations: Vec<_> = (0..cells)
        .map(|cell| {
            let perturbation = 1e-6 * (k * coordinates(cell, dims)[axis] as f64).cos();
            let mut velocity = [0.0; 3];
            velocity[if longitudinal { axis } else { (axis + 1) % 3 }] = if longitudinal {
                cs * perturbation
            } else {
                perturbation
            };
            equilibrium3(
                if longitudinal {
                    1.0 + perturbation
                } else {
                    1.0
                },
                velocity,
            )
        })
        .collect();
    let amplitude = |populations: &[[f64; Q3]]| {
        let mut mode = [0.0_f64; 2];
        for (cell, f) in populations.iter().enumerate() {
            let (rho, velocity) = macro_cell::<D3q19, Q3>(f, [0.0; 3], cell, 0.2).unwrap();
            let wave = if longitudinal {
                rho - 1.0 + velocity[axis] / cs
            } else {
                velocity[(axis + 1) % 3]
            };
            let phase = k * coordinates(cell, dims)[axis] as f64;
            mode[0] += wave * phase.cos();
            mode[1] += wave * phase.sin();
        }
        mode[0].hypot(mode[1])
    };
    let mut post = vec![[0.0; Q3]; cells];
    let mut initial = 0.0;
    for step in 1..=320 {
        for (tile, output) in post.chunks_mut(64).enumerate() {
            data.process(&populations, output, tile * 64, false)
                .unwrap();
        }
        for (tile, output) in populations.chunks_mut(64).enumerate() {
            data.process(&post, output, tile * 64, true).unwrap();
        }
        if step == 32 {
            initial = amplitude(&populations);
        }
    }
    let final_amplitude = amplitude(&populations);
    assert!(final_amplitude > 0.0 && final_amplitude < initial);
    (initial / final_amplitude).ln() / 288.0
}

#[test]
fn native_sound_decay_recovers_three_dimensional_bulk_viscosity() {
    let nu = (0.68 - 0.5) / 3.0;
    for (nx, axis) in [(64, 0), (64, 1), (64, 2), (128, 0)] {
        let k = core::f64::consts::TAU / nx as f64;
        for bulk_tau in [0.6, 1.8] {
            let bulk_nu = (2.0 / 9.0) * (bulk_tau - 0.5);
            let expected = 0.5 * ((4.0 / 3.0) * nu + bulk_nu) * k * k;
            let measured = decay(nx, axis, bulk_tau, true);
            let error = (measured / expected - 1.0).abs();
            println!("D3Q19 sound nx={nx} axis={axis} tau_bulk={bulk_tau}: relative error={error}");
            assert!(error < 0.02, "attenuation {measured} != {expected}");
        }
    }
}

#[test]
fn native_shear_decay_is_independent_of_bulk_relaxation() {
    let k = core::f64::consts::TAU / 64.0;
    let expected = (0.68 - 0.5) / 3.0 * k * k;
    for axis in 0..3 {
        for bulk_tau in [0.6, 1.8] {
            let measured = decay(64, axis, bulk_tau, false);
            let error = (measured / expected - 1.0).abs();
            println!("D3Q19 shear axis={axis} tau_bulk={bulk_tau}: relative error={error}");
            assert!(error < 0.015, "attenuation {measured} != {expected}");
        }
    }
}

#[test]
fn material_and_sgs_shear_rebinding_preserve_bulk_and_higher_rates() {
    let original = CollisionModel3::CentralMomentBulk {
        shear_rate: 1.2,
        bulk_rate: 0.7,
        higher_order_rate: 1.6,
    };
    assert!((original.kinematic_bulk_viscosity() - (2.0 / 9.0) * (1.0 / 0.7 - 0.5)).abs() < 1e-15);
    for tau in [0.56, 0.9, 1.7] {
        let changed = D3q19::with_shear_relaxation(original, tau).unwrap();
        assert_eq!(
            changed,
            CollisionModel3::CentralMomentBulk {
                shear_rate: 1.0 / tau,
                bulk_rate: 0.7,
                higher_order_rate: 1.6,
            }
        );
        assert_eq!(
            changed.kinematic_bulk_viscosity().to_bits(),
            original.kinematic_bulk_viscosity().to_bits()
        );
        assert!((changed.kinematic_viscosity() - (tau - 0.5) / 3.0).abs() < 1e-15);
    }
}
