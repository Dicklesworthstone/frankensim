use super::*;
use crate::d3q19::CollisionModel3;
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use crate::lattice3::{D3q19, GridBudget3};

fn controls() -> WrenchImpulseControls3 {
    WrenchImpulseControls3 {
        impulse: ImpulseControls {
            max_cells: 128,
            maximum_speed: 0.1,
            maximum_momentum_defect: 2e-12,
            maximum_work_defect: 2e-12,
        },
        maximum_angular_momentum_defect: 3e-12,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn dense<L: ImpulseLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    periodic: [bool; 3],
) -> Grid3<L, Q> {
    Grid3::try_new(
        [2, 1, 1],
        collision,
        [0.0; 3],
        periodic,
        GridBudget3 {
            maximum_cells: 128,
            maximum_mach: 0.3,
        },
        |_| false,
        |p| L::W.map(|w| w * if p[0] == 0 { 1.0 } else { 3.0 }),
    )
    .unwrap()
}

fn pure_couple<L: ImpulseLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let mut grid = dense::<L, Q>(collision, [false; 3]);
    let mass: Vec<f64> = grid.populations().iter().map(|f| f.iter().sum()).collect();
    let request = WrenchImpulse3 {
        origin: [0.0; 3],
        linear: [0.0; 3],
        angular: [0.0, 0.0, 0.003],
    };
    let report = grid
        .apply_wrench_impulse(&[0, 1], request, controls(), &CancelGate::new_clock_free())
        .unwrap();
    for a in 0..3 {
        close(report.impulse.momentum_change[a], 0.0, 2e-15);
        close(report.angular_momentum_change[a], request.angular[a], 2e-15);
    }
    for (i, expected_y) in [-0.003, 0.001].into_iter().enumerate() {
        let (rho, u) = grid.macroscopic(i).unwrap();
        close(rho, mass[i], 2e-15);
        close(u[1], expected_y, 2e-15);
        close(u[0], 0.0, 2e-15);
        close(u[2], 0.0, 2e-15);
    }
    close(report.impulse.work, 6e-6, 2e-15);
    close(report.impulse.kinetic_energy_change, 6e-6, 2e-15);
    assert_eq!(grid.completed_steps(), 0);
}

#[test]
fn native_wrench_couple_uses_actual_inertia_on_d3q19_and_d3q27() {
    pure_couple::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
    pure_couple::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.8,
            equilibrium: Equilibrium27::Quadratic,
        },
        maximum_mach: 0.3,
    });
}

#[test]
fn native_wrench_retains_spatial_force_epoch_and_observes_physical_torque() {
    let grid = dense::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, [false; 3]);
    let forces = vec![[1e-5, 2e-5, -3e-5], [-2e-5, 1e-5, 3e-5]];
    // This is the admitted initial force-frame convention: rho*u = raw+F/2.
    let mut source = crate::lattice3::sources::ForcedGrid3 { grid, forces };
    let original_forces = source.forces.clone();
    let original: Vec<_> = (0..2).map(|i| source.macroscopic(i).unwrap()).collect();
    let request = WrenchImpulse3 {
        origin: [0.0; 3],
        linear: [0.002, 0.0, 0.0],
        angular: [0.0, 0.001, 0.002],
    };
    let report = source
        .apply_wrench_impulse(&[0, 1], request, controls(), &CancelGate::new_clock_free())
        .unwrap();
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (i, &(rho0, u0)) in original.iter().enumerate() {
        let (rho1, u1) = source.macroscopic(i).unwrap();
        let j: [f64; 3] = core::array::from_fn(|a| rho1 * u1[a] - rho0 * u0[a]);
        let r = [i as f64 + 0.5, 0.5, 0.5];
        let l = [
            r[1] * j[2] - r[2] * j[1],
            r[2] * j[0] - r[0] * j[2],
            r[0] * j[1] - r[1] * j[0],
        ];
        for a in 0..3 {
            linear[a] += j[a];
            angular[a] += l[a];
        }
    }
    for a in 0..3 {
        close(linear[a], request.linear[a], 3e-15);
        close(angular[a], request.angular[a], 3e-15);
        close(report.angular_momentum_change[a], angular[a], 3e-15);
    }
    assert_eq!(source.forces, original_forces);
    assert_eq!(source.completed_steps(), 0);
}

#[test]
fn native_wrench_refuses_unresolved_periodic_solid_and_duplicate_stencils_atomically() {
    let mut grid = dense::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, [false; 3]);
    let before = grid.populations().to_vec();
    let request = WrenchImpulse3 {
        origin: [0.0; 3],
        linear: [0.0; 3],
        angular: [0.0, 0.0, 0.003],
    };
    let gate = CancelGate::new_clock_free();
    for indices in [&[0, 0][..], &[0, 2][..], &[][..]] {
        assert!(
            grid.apply_wrench_impulse(indices, request, controls(), &gate)
                .is_err()
        );
        assert_eq!(grid.populations(), before);
    }
    let axial = WrenchImpulse3 {
        angular: [0.001, 0.0, 0.0],
        ..request
    };
    assert!(
        grid.apply_wrench_impulse(&[0, 1], axial, controls(), &gate)
            .is_err()
    );
    assert_eq!(grid.populations(), before);
    grid.solid[1] = true;
    assert!(
        grid.apply_wrench_impulse(&[0, 1], request, controls(), &gate)
            .is_err()
    );
    assert_eq!(grid.populations(), before);
    grid.solid[1] = false;
    grid.periodic[0] = true;
    assert!(
        grid.apply_wrench_impulse(&[0, 1], request, controls(), &gate)
            .is_err()
    );
    assert_eq!(grid.populations(), before);
    grid.periodic[0] = false;
    gate.request();
    assert_eq!(
        grid.apply_wrench_impulse(&[0, 1], request, controls(), &gate),
        Err(ImpulseError::Cancelled)
    );
    assert_eq!(grid.populations(), before);
    assert_eq!(grid.completed_steps(), 0);
}

#[test]
fn native_wrench_zero_late_torque_refusal_and_final_cancellation_keep_populations() {
    let mut grid = dense::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 }, [false; 3]);
    let before = grid.populations().to_vec();
    let gate = CancelGate::new_clock_free();
    let zero = WrenchImpulse3 {
        origin: [0.0; 3],
        linear: [0.0; 3],
        angular: [0.0; 3],
    };
    let report = grid
        .apply_wrench_impulse(&[0, 1], zero, controls(), &gate)
        .unwrap();
    assert_eq!(report.impulse.affected_cells, 0);
    assert_eq!(grid.populations(), before);

    // A single-cell point force is geometrically exact. The zero angular
    // allowance must nevertheless reject finite-precision population readback,
    // BEFORE publishing any of the staged impulse.
    let request = WrenchImpulse3 {
        linear: [0.001, 0.0, 0.0],
        angular: [0.0, 0.0005, -0.0005],
        ..zero
    };
    let strict = WrenchImpulseControls3 {
        maximum_angular_momentum_defect: 0.0,
        ..controls()
    };
    assert_eq!(
        grid.apply_wrench_impulse(&[0], request, strict, &gate),
        Err(ImpulseError::Balance)
    );
    assert_eq!(grid.populations(), before);
    assert_ne!(grid.next, before);

    let result = grid.apply_impulse_with_receipt(
        controls().impulse,
        None,
        &gate,
        |_| [0.001, 0.0, 0.0],
        |report, _, _| {
            gate.request();
            Ok(report)
        },
    );
    assert_eq!(result, Err(ImpulseError::Cancelled));
    assert_eq!(grid.populations(), before);
    assert_eq!(grid.completed_steps(), 0);
}
