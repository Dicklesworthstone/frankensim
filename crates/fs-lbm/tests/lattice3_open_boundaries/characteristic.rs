//! Native outgoing acoustic and tangential modes, source frames and rollback.

use super::*;
use fs_lbm::d3q19::BoundaryProfileError3;
use fs_lbm::lattice3::open_boundaries::CharacteristicOutlet3;
use fs_lbm::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;

fn control(speed: f64) -> CharacteristicOutlet3 {
    CharacteristicOutlet3 {
        reference_normal_velocity_lu: speed,
        incoming_characteristic_lu: 0.0,
        maximum_relative_density_departure: 0.05,
    }
}

fn fixture<L: EquilibriumLattice3<Q>, const Q: usize>(
    cx: &Cx<'_>,
    collision: L::Collision,
    face: Face3,
    force: [f64; 3],
    density: f64,
    velocity: [f64; 3],
) -> OpenBoundaryGrid3<L, Q> {
    let axis = face as usize / 2;
    let mut periodic = [true; 3];
    periodic[axis] = false;
    let grid = Grid3::<L, Q>::try_new_equilibrium(
        [5; 3],
        collision,
        force,
        periodic,
        budget(),
        |_| false,
        |_| (density, velocity),
    )
    .unwrap();
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[face as usize] = FaceBoundary3::Pressure { density: 1.3 };
    faces[face as usize ^ 1] = FaceBoundary3::Velocity { velocity };
    OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
}

fn trace_expected<L: Lattice3<Q>, const Q: usize>(
    face: Face3,
    target: CharacteristicOutlet3,
    density: f64,
    velocity: [f64; 3],
) -> (f64, [f64; 3]) {
    let axis = face as usize / 2;
    let sign = if face as usize % 2 == 0 { -1.0 } else { 1.0 };
    let c = L::CS2.sqrt();
    let plus =
        sign * velocity[axis] - target.reference_normal_velocity_lu + c * (density - 1.3) / 1.3;
    let mut next = velocity;
    next[axis] = sign
        * (target.reference_normal_velocity_lu + (plus + target.incoming_characteristic_lu) / 2.0);
    (
        1.3 * (1.0 + (plus - target.incoming_characteristic_lu) / (2.0 * c)),
        next,
    )
}

fn physical_modes<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        for face in Face3::ALL {
            let axis = face as usize / 2;
            let sign = if face as usize % 2 == 0 { -1.0 } else { 1.0 };
            let mut velocity = [0.007, -0.004, 0.003];
            velocity[axis] = sign * 0.045;
            for force in [[0.0; 3], [0.00011, -0.00007, 0.00009]] {
                let mut flow = fixture::<L, Q>(cx, collision, face, force, 1.3026, velocity);
                let target = CharacteristicOutlet3 {
                    incoming_characteristic_lu: 0.0004,
                    ..control(0.04)
                };
                flow.set_characteristic_outlet(face, target, &CancelGate::new_clock_free())
                    .unwrap();
                let before = totals(&flow);
                let receipt = flow.step_serial(cx, None, [0.0; 3]).unwrap();
                let after = totals(&flow);
                near(after.0 - before.0, receipt.mass_change, 3e-12);
                for a in 0..3 {
                    near(
                        after.1[a] - before.1[a],
                        receipt.impulse[a] + 125.0 * force[a],
                        2e-12,
                    );
                }
                let mut point = [2; 3];
                point[axis] = if sign < 0.0 { 0 } else { 4 };
                let cell = address(point, [5; 3]);
                point[axis] = if sign < 0.0 { 1 } else { 3 };
                let donor = address(point, [5; 3]);
                let (rho_d, u_d) = flow.macroscopic(donor).unwrap();
                let expected = trace_expected::<L, Q>(face, target, rho_d, u_d);
                let actual = flow.macroscopic(cell).unwrap();
                near(actual.0, expected.0, 3e-14);
                for a in 0..3 {
                    near(actual.1[a], expected.1[a], 3e-14);
                }
                let donor_stress = nonequilibrium_stress(&flow, donor);
                let trace_stress = nonequilibrium_stress(&flow, cell);
                for a in 0..3 {
                    for b in 0..3 {
                        near(trace_stress[a][b], donor_stress[a][b], 3e-14);
                    }
                }
                let minus = sign * actual.1[axis]
                    - target.reference_normal_velocity_lu
                    - L::CS2.sqrt() * (actual.0 - 1.3) / 1.3;
                near(minus, target.incoming_characteristic_lu, 3e-14);
            }
        }
    });
}

#[test]
fn native_characteristics_preserve_tangential_modes_stress_and_force_on_every_face() {
    physical_modes::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    physical_modes::<D3q19, 19>(CollisionModel3::Trt {
        tau: 0.9,
        magic: 0.1875,
    });
    physical_modes::<D3q19, 19>(CollisionModel3::CentralMomentBulk {
        shear_rate: 1.1,
        bulk_rate: 0.8,
        higher_order_rate: 1.2,
    });
    physical_modes::<D3q27, 27>(bgk());
    physical_modes::<D3q27, 27>(config(CollisionModel27::Bgk {
        tau: 0.9,
        equilibrium: Equilibrium27::Quadratic,
    }));
    physical_modes::<D3q27, 27>(config(CollisionModel27::Cumulant {
        rates: Rates27 { omega: [1.1; 10] },
    }));
}

fn acceleration_trace<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let face = Face3::XMax;
        let mut native = fixture::<L, Q>(cx, collision, face, [0.0; 3], 1.303, [0.04, 0.005, 0.0]);
        native
            .set_characteristic_outlet(face, control(0.04), &CancelGate::new_clock_free())
            .unwrap();
        let mut flow = ForcedOpenBoundaryGrid3::try_new_with_body_acceleration_preserving_velocity(
            native,
            cx,
            vec![[0.0001, 0.00002, -0.00003]; 125],
        )
        .unwrap();
        let end: Vec<_> = (0..125)
            .map(|cell| [0.00002 * (cell % 5) as f64, -0.00004, 0.00001])
            .collect();
        flow.advance_serial_with_body_acceleration(cx, None, &end, [0.0; 3])
            .unwrap();
        let donor = address([3, 2, 2], [5; 3]);
        let cell = address([4, 2, 2], [5; 3]);
        let (rho_d, u_d) = flow.macroscopic(donor).unwrap();
        let expected = trace_expected::<L, Q>(face, control(0.04), rho_d, u_d);
        let actual = flow.macroscopic(cell).unwrap();
        near(actual.0, expected.0, 3e-14);
        for a in 0..3 {
            near(actual.1[a], expected.1[a], 3e-14);
            near(flow.forces()[cell][a], actual.0 * end[cell][a], 1e-17);
        }
        assert_eq!(flow.characteristic_outlet(face), Some(control(0.04)));
        let saved = flow.populations().to_vec();
        let forces = flow.forces().to_vec();
        let end = vec![[-0.2, 0.0, 0.0]; 125];
        assert!(
            flow.advance_serial_with_body_acceleration(cx, None, &end, [0.0; 3])
                .is_err()
        );
        assert_eq!(flow.populations(), saved);
        assert_eq!(flow.forces(), forces);
        assert_eq!(flow.completed_steps(), 1);
    });
}

#[test]
fn characteristic_trace_uses_actual_density_acceleration_and_preserves_source_epoch_on_refusal() {
    acceleration_trace::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    acceleration_trace::<D3q27, 27>(bgk());
}

#[test]
fn characteristic_control_ownership_limits_and_backflow_refuse_without_publication() {
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let make = || fixture::<D3q27, 27>(cx, bgk(), Face3::XMax, [0.0; 3], 1.3, [0.04, 0.0, 0.0]);
        let mut flow = make();
        let mut equilibrium = make();
        equilibrium
            .set_characteristic_outlet(Face3::XMax, control(0.04), &gate)
            .unwrap();
        for _ in 0..8 {
            equilibrium.step_serial(cx, None, [0.0; 3]).unwrap();
            for cell in 0..125 {
                let (rho, velocity) = equilibrium.macroscopic(cell).unwrap();
                near(rho, 1.3, 3e-13);
                for a in 0..3 {
                    near(velocity[a], if a == 0 { 0.04 } else { 0.0 }, 3e-14);
                }
            }
        }
        assert!(
            flow.set_characteristic_outlet(Face3::XMin, control(0.04), &gate)
                .is_err()
        );
        flow.set_pressure_profile(Face3::XMax, &gate, |_| Ok(1.3))
            .unwrap();
        assert!(
            flow.set_characteristic_outlet(Face3::XMax, control(0.04), &gate)
                .is_err()
        );
        flow.clear_pressure_profile(Face3::XMax);
        flow.set_characteristic_outlet(Face3::XMax, control(0.04), &gate)
            .unwrap();
        let mut called = false;
        assert_eq!(
            flow.set_pressure_profile(Face3::XMax, &gate, |_| {
                called = true;
                Ok(1.3)
            }),
            Err(BoundaryProfileError3::ConflictingControl)
        );
        assert!(!called);
        for invalid in [
            control(0.0),
            control(-0.02),
            control(1.0),
            CharacteristicOutlet3 {
                incoming_characteristic_lu: f64::NAN,
                ..control(0.04)
            },
            CharacteristicOutlet3 {
                maximum_relative_density_departure: 1.0,
                ..control(0.04)
            },
        ] {
            assert!(
                flow.set_characteristic_outlet(Face3::XMax, invalid, &gate)
                    .is_err()
            );
            assert_eq!(flow.characteristic_outlet(Face3::XMax), Some(control(0.04)));
        }
        let mut faces = core::array::from_fn(|a| flow.boundaries().face(Face3::ALL[a]));
        faces[1] = FaceBoundary3::Pressure { density: 1.31 };
        assert!(
            flow.rebind_open_boundaries(BoundarySpec3::new(faces), &gate)
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(
            flow.set_characteristic_outlet(Face3::XMax, control(0.05), &cancelled)
                .is_err()
        );
        assert_eq!(flow.characteristic_outlet(Face3::XMax), Some(control(0.04)));
        // Removing only the control restores the exact original pressure path.
        flow.clear_characteristic_outlet(Face3::XMax);
        let mut original = make();
        for _ in 0..3 {
            flow.step_serial(cx, None, [0.0; 3]).unwrap();
            original.step_serial(cx, None, [0.0; 3]).unwrap();
            assert_eq!(flow.populations(), original.populations());
        }
        for (rho, u) in [(1.3, -0.02), (1.5, 0.04)] {
            let mut rejected =
                fixture::<D3q27, 27>(cx, bgk(), Face3::XMax, [0.0; 3], rho, [u, 0.0, 0.0]);
            rejected
                .set_characteristic_outlet(Face3::XMax, control(0.04), &gate)
                .unwrap();
            let before = rejected.populations().to_vec();
            assert!(rejected.step_serial(cx, None, [0.0; 3]).is_err());
            assert_eq!(rejected.populations(), before);
            assert_eq!(rejected.completed_steps(), 0);
        }
        let grid = Grid3::<D3q27, 27>::try_new_equilibrium(
            [5; 3],
            bgk(),
            [0.0; 3],
            [false, false, true],
            budget(),
            |_| false,
            |_| (1.3, [0.04, 0.0, 0.0]),
        )
        .unwrap();
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity {
            velocity: [0.04, 0.0, 0.0],
        };
        faces[1] = FaceBoundary3::Pressure { density: 1.3 };
        faces[2] = FaceBoundary3::Wall { velocity: [0.0; 3] };
        faces[3] = faces[2];
        let mut walls = OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap();
        assert!(
            walls
                .set_characteristic_outlet(Face3::XMax, control(0.04), &gate)
                .is_err()
        );
        assert_eq!(walls.characteristic_outlet(Face3::XMax), None);
    });
}

fn pulse_reflection<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    upper: bool,
) {
    with_cx(|cx| {
        let dims = [128, 2, 2];
        let sign = if upper { 1.0 } else { -1.0 };
        let face = if upper { Face3::XMax } else { Face3::XMin };
        let c = L::CS2.sqrt();
        let mean = 0.04;
        let make = || {
            let grid = Grid3::<L, Q>::try_new_equilibrium(
                dims,
                collision,
                [0.0; 3],
                [false, true, true],
                budget(),
                |_| false,
                |p| {
                    let outward_x = if upper { p[0] } else { 127 - p[0] } as f64;
                    let perturbation = 0.0001 * (-0.5 * ((outward_x - 64.0) / 8.0).powi(2)).exp();
                    (
                        1.3 * (1.0 + perturbation),
                        [sign * (mean + c * perturbation), 0.0, 0.0],
                    )
                },
            )
            .unwrap();
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[face as usize] = FaceBoundary3::Pressure { density: 1.3 };
            faces[face as usize ^ 1] = FaceBoundary3::Velocity {
                velocity: [sign * mean, 0.0, 0.0],
            };
            OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
        };
        let mut pressure = make();
        let mut characteristic = make();
        characteristic
            .set_characteristic_outlet(face, control(mean), &CancelGate::new_clock_free())
            .unwrap();
        for _ in 0..180 {
            pressure.step_serial(cx, None, [0.0; 3]).unwrap();
            characteristic.step_serial(cx, None, [0.0; 3]).unwrap();
        }
        // The outgoing pulse has left. Compare the returning acoustic J- in
        // the actual original populations before any inlet rereflection.
        let returning = |flow: &OpenBoundaryGrid3<L, Q>| {
            (12..116)
                .map(|x| {
                    let (rho, u) = flow.macroscopic(address([x, 0, 0], dims)).unwrap();
                    (sign * u[0] - mean - c * (rho - 1.3) / 1.3).powi(2)
                })
                .sum::<f64>()
                .sqrt()
        };
        let fixed = returning(&pressure);
        let outgoing = returning(&characteristic);
        eprintln!(
            "D3Q{Q} {face:?}: returning J- norm, pressure={fixed:.9e}, characteristic={outgoing:.9e}, ratio={:.9e}",
            outgoing / fixed,
        );
        assert!(
            fixed > 1e-4,
            "pressure reflection must be resolved: {fixed:e}"
        );
        assert!(
            outgoing < 0.1 * fixed,
            "characteristic={outgoing:e}, pressure={fixed:e}"
        );
    });
}

#[test]
fn outgoing_acoustic_pulse_reflects_less_than_fixed_pressure_on_both_stencils_and_normals() {
    for upper in [false, true] {
        pulse_reflection::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.6 }, upper);
        pulse_reflection::<D3q27, 27>(
            config(CollisionModel27::Bgk {
                tau: 0.6,
                equilibrium: Equilibrium27::Factorized,
            }),
            upper,
        );
    }
}
