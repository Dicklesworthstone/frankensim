//! G0 native transaction/balance and G1 full-tensor diffusion modes.
use super::*;
use crate::thermal::tensor::DiffusivityTensor;
use fs_math::det::{cos, exp};

fn isotropic_native<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        for axis in 0..3 {
            let make = || {
                let mut velocity = [0.004, -0.003, 0.005];
                velocity[axis] = -0.025; // Reversed flow exercises both reservoir roles.
                fixture::<L, Q>(cx, collision, axis, velocity, [0.0; 3])
            };
            let mut scalar = make();
            let mut tensor = make();
            let cells = scalar.populations().len();
            let initial: Vec<_> = (0..cells).map(|i| 0.7 + 0.002 * i as f64).collect();
            let prepare = |flow: &OpenBoundaryGrid3<L, Q>| {
                flow.prepare_thermal_field(cx, &initial, 0.9, thermal_faces(axis, 0.8), limits(3.0))
                    .unwrap()
            };
            let mut a = prepare(&scalar);
            let mut b = prepare(&tensor);
            let law = DiffusivityTensor::from_lattice([
                [0.1, 0.0, 0.0],
                [0.0, 0.1, 0.0],
                [0.0, 0.0, 0.1],
            ])
            .unwrap();
            let field = vec![law; cells];
            let source: Vec<_> = (0..cells).map(|i| 0.0001 * (i % 7) as f64).collect();
            for step in 0..12 {
                let expected = scalar
                    .step_serial_with_thermal(&mut a, cx, None, &source, None, [0.0; 3])
                    .unwrap();
                // A uniform law and its per-cell representation share semantics.
                let diffusion = if step % 2 == 0 {
                    DiffusivityField::Uniform(&law)
                } else {
                    DiffusivityField::PerCell(&field)
                };
                let actual = tensor
                    .step_serial_with_thermal_tensor_diffusion(
                        &mut b, cx, None, &source, diffusion, [0.0; 3],
                    )
                    .unwrap();
                assert_eq!(scalar.populations(), tensor.populations());
                assert_eq!(expected.momentum, actual.momentum);
                near(actual.heat.defect, 0.0, 2e-12);
                near(actual.heat.source, expected.heat.source, 2e-13);
                for face in 0..6 {
                    near(
                        actual.heat.reservoirs[face].heat_in,
                        expected.heat.reservoirs[face].heat_in,
                        2e-13,
                    );
                }
                for (a, b) in a
                    .populations()
                    .iter()
                    .flatten()
                    .zip(b.populations().iter().flatten())
                {
                    near(*a, *b, 3e-14);
                }
            }
        }
    });
}

#[test]
fn native_tensor_isotropic_limit_preserves_open_flux_and_both_momentum_lattices() {
    isotropic_native::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    isotropic_native::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk {
            tau: 0.9,
            equilibrium: Equilibrium27::Factorized,
        },
        maximum_mach: 0.2,
    });
}

fn mode_error(cx: &Cx<'_>, n: usize, axes: [usize; 2], sign: f64) -> f64 {
    let mut dims = [1; 3];
    for axis in axes {
        dims[axis] = n;
    }
    let cells: usize = dims.into_iter().product();
    let grid = Grid3::<D3q19, 19>::try_new_equilibrium(
        dims,
        CollisionModel3::Bgk { tau: 0.9 },
        [0.0; 3],
        [true; 3],
        GridBudget3 {
            maximum_cells: cells,
            maximum_mach: 0.2,
        },
        |_| false,
        |_| (1.0, [0.0; 3]),
    )
    .unwrap();
    let mut flow =
        OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new([FaceBoundary3::Periodic; 6]))
            .unwrap();
    let basis: Vec<_> = (0..cells)
        .map(|mut cell| {
            let coordinates: [usize; 3] = core::array::from_fn(|axis| {
                let position = cell % dims[axis];
                cell /= dims[axis];
                position
            });
            let phase = std::f64::consts::TAU / n as f64
                * (coordinates[axes[0]] as f64 + sign * coordinates[axes[1]] as f64);
            cos(phase)
        })
        .collect();
    let initial: Vec<_> = basis.iter().map(|&wave| 1.0 + 0.1 * wave).collect();
    let mut budget = limits(3.0);
    budget.maximum_cells = cells;
    let mut heat = flow
        .prepare_thermal_field(
            cx,
            &initial,
            0.9,
            ThermalOpenBoundaries3 {
                faces: [ThermalOpenFace3::Periodic; 6],
                internal_wall: ThermalWall::Adiabatic,
            },
            budget,
        )
        .unwrap();
    let alpha = [
        [0.12, 0.035, 0.025],
        [0.035, 0.09, -0.02],
        [0.025, -0.02, 0.07],
    ];
    let law = DiffusivityTensor::from_lattice(alpha).unwrap();
    let source = vec![0.0; cells];
    let steps = n * n / 16;
    for _ in 0..steps {
        let result = flow
            .step_serial_with_thermal_tensor_diffusion(
                &mut heat,
                cx,
                None,
                &source,
                DiffusivityField::Uniform(&law),
                [0.0; 3],
            )
            .unwrap();
        near(result.heat.stored_change, 0.0, 2e-11);
    }
    let amplitude = 2.0 / cells as f64
        * basis
            .iter()
            .enumerate()
            .map(|(i, wave)| (heat.temperature(i).unwrap() - 1.0) * wave)
            .sum::<f64>();
    let k = std::f64::consts::TAU / n as f64;
    // Independent continuum Fourier symbol k^T alpha k, including cross terms.
    let rate = k
        * k
        * (alpha[axes[0]][axes[0]]
            + alpha[axes[1]][axes[1]]
            + 2.0 * sign * alpha[axes[0]][axes[1]]);
    (amplitude / (0.1 * exp(-rate * steps as f64)) - 1.0).abs()
}

#[test]
fn native_tensor_rotated_fourier_decay_refines_in_every_coordinate_plane() {
    with_cx(|cx| {
        for axes in [[0, 1], [0, 2], [1, 2]] {
            for sign in [-1.0, 1.0] {
                let coarse = mode_error(cx, 16, axes, sign);
                let fine = mode_error(cx, 32, axes, sign);
                assert!(fine < 0.006, "axes={axes:?} sign={sign}: {fine}");
                assert!(
                    fine < 0.32 * coarse,
                    "axes={axes:?} sign={sign}: {coarse} -> {fine}"
                );
            }
        }
    });
}

#[test]
fn native_tensor_refusal_preserves_flow_heat_history_and_scalar_retry() {
    with_cx(|cx| {
        let make = || {
            fixture::<D3q19, 19>(
                cx,
                CollisionModel3::Bgk { tau: 0.9 },
                0,
                [0.025, 0.0, 0.0],
                [0.0; 3],
            )
        };
        let mut flow = make();
        let mut control = make();
        let cells = flow.populations().len();
        let prepare = |native: &OpenBoundaryGrid3<D3q19, 19>| {
            native
                .prepare_thermal_field(
                    cx,
                    &vec![0.4; cells],
                    0.9,
                    thermal_faces(0, 0.4),
                    limits(0.7),
                )
                .unwrap()
        };
        let mut heat = prepare(&flow);
        let mut expected = prepare(&control);
        let law = DiffusivityTensor::from_lattice([
            [0.12, 0.02, 0.0],
            [0.02, 0.09, -0.01],
            [0.0, -0.01, 0.07],
        ])
        .unwrap();
        let before_flow = flow.populations().to_vec();
        let before_heat = heat.populations().to_vec();
        let malformed = vec![law; cells - 1];
        assert!(matches!(
            flow.step_serial_with_thermal_tensor_diffusion(
                &mut heat,
                cx,
                None,
                &vec![0.0; cells],
                DiffusivityField::PerCell(&malformed),
                [0.0; 3],
            ),
            Err(OpenThermalError3::Thermal(ThermalOpenError3::Invalid(_)))
        ));
        // This passes the first source half-kick and transport but exceeds
        // the temperature domain at the final kick, after momentum staging.
        assert!(matches!(
            flow.step_serial_with_thermal_tensor_diffusion(
                &mut heat,
                cx,
                None,
                &vec![0.4; cells],
                DiffusivityField::Uniform(&law),
                [0.0; 3],
            ),
            Err(OpenThermalError3::Thermal(
                ThermalOpenError3::Temperature { .. }
            ))
        ));
        let mut checkpoints = 0;
        let cancelled = heat.begin_step_tensor::<ThermalOpenError3>(
            &vec![0.001; cells],
            DiffusivityField::Uniform(&law),
            || {
                checkpoints += 1;
                if checkpoints == 4 {
                    Err(ThermalOpenError3::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(checkpoints, 4);
        assert!(matches!(cancelled, Err(ThermalOpenError3::Cancelled)));
        assert_eq!(flow.populations(), before_flow);
        assert_eq!(heat.populations(), before_heat);
        assert_eq!(flow.completed_steps(), 0);
        assert_eq!(heat.completed_steps(), 0);
        assert!(heat.last_balance().is_none());
        let source = vec![0.001; cells];
        let actual = flow
            .step_serial_with_thermal(&mut heat, cx, None, &source, None, [0.0; 3])
            .unwrap();
        let reference = control
            .step_serial_with_thermal(&mut expected, cx, None, &source, None, [0.0; 3])
            .unwrap();
        assert_eq!(flow.populations(), control.populations());
        assert_eq!(heat.populations(), expected.populations());
        assert_eq!(actual.momentum, reference.momentum);
        assert_eq!(actual.heat, reference.heat);
    });
}
