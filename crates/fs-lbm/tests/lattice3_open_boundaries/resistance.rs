//! Porous constitutive endpoints on the native open-flow transaction.

use super::*;
use fs_lbm::d3q19::{InletFlux3, ResistanceMomentumExchange3};
use fs_lbm::lattice3::open_boundaries::CharacteristicOutlet3;
use fs_lbm::lattice3::open_boundaries::forcing::ForcedOpenBoundaryGrid3;
use fs_lbm::lattice3::open_boundaries::forcing::resistance::{
    DarcyForchheimer3, LinearResistance3, ResistanceOpenGrid3,
};

const AXES: [[f64; 3]; 3] = [
    [2.0 / 3.0, -2.0 / 3.0, 1.0 / 3.0],
    [1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
    [-2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0],
];

mod endpoint;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}

fn drag(law: DarcyForchheimer3, u: [f64; 3]) -> [f64; 3] {
    let slip = core::array::from_fn(|a| u[a] - law.linear.skeleton_velocity[a]);
    let speed = dot(slip, slip).sqrt();
    core::array::from_fn(|a| {
        -law.quadratic * speed * slip[a]
            - (0..3).map(|k| {
                law.linear.principal[k] * law.linear.axes[k][a] * dot(law.linear.axes[k], slip)
            }).sum::<f64>()
    })
}

fn porous_law() -> DarcyForchheimer3 {
    DarcyForchheimer3 {
        linear: LinearResistance3 {
            principal: [0.04, 0.08, 0.025],
            axes: AXES,
            skeleton_velocity: [0.011, -0.004, 0.006],
            body_force: [0.0; 3],
        },
        quadratic: 1.7,
    }
}

fn pressure_faces(rho: f64) -> BoundarySpec3 {
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[0] = FaceBoundary3::Pressure { density: rho };
    faces[1] = faces[0];
    BoundarySpec3::new(faces)
}

fn source_stress<L: EquilibriumLattice3<Q>, const Q: usize>(
    flow: &ForcedOpenBoundaryGrid3<L, Q>,
    cell: usize,
) -> [[f64; 3]; 3] {
    let (rho, u) = flow.macroscopic(cell).unwrap();
    let force = flow.forces()[cell];
    core::array::from_fn(|a| core::array::from_fn(|b| {
        flow.populations()[cell].iter().zip(L::E)
            .map(|(f, e)| f * f64::from(e[a] * e[b])).sum::<f64>()
            + 0.5 * (u[a] * force[b] + force[a] * u[b])
            - rho * (u[a] * u[b] + if a == b { L::CS2 } else { 0.0 })
    }))
}

fn attachment<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let dims = [6, 3, 3];
        let grid = Grid3::<L, Q>::try_new_equilibrium(
            dims, collision, [0.0002, -0.0001, 0.00005], [false, true, true], budget(),
            |_| false, |p| (1.1 + 0.0002 * p[0] as f64, [0.021, -0.007, 0.011]),
        ).unwrap();
        let mut flow = OpenBoundaryGrid3::try_new(grid, cx, boundaries(0, [0.021, -0.007, 0.011], 1.1)).unwrap();
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        flow.step_serial(cx, None, [0.0; 3]).unwrap();
        let macros: Vec<_> = (0..54).map(|cell| flow.macroscopic(cell).unwrap()).collect();
        let stress: Vec<_> = (0..54).map(|cell| nonequilibrium_stress(&flow, cell)).collect();
        let mut law = porous_law();
        law.linear.body_force = [0.0003, -0.0002, 0.0001];
        let porous = ResistanceOpenGrid3::try_new_preserving_velocity(flow, cx, vec![law; 54], [0.5, -1.0, 0.25]).unwrap();
        assert_eq!(porous.flow().completed_steps(), 2);
        for cell in 0..54 {
            let (rho, u) = porous.flow().macroscopic(cell).unwrap();
            near(rho, macros[cell].0, 8e-14);
            for a in 0..3 {
                near(u[a], macros[cell].1[a], 8e-14);
                near(porous.flow().forces()[cell][a], law.linear.body_force[a] + drag(law, u)[a], 8e-14);
                for b in 0..3 {
                    near(source_stress(porous.flow(), cell)[a][b], stress[cell][a][b], 8e-14);
                }
            }
        }
    });
}

#[test]
fn porous_attachment_preserves_evolved_physical_velocity_and_stress_on_both_stencils() {
    attachment::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    attachment::<D3q27, 27>(bgk());
}

fn prescribed_endpoint<L: EquilibriumLattice3<Q>, const Q: usize>(
    collision: L::Collision,
    boost: [f64; 3],
) {
    with_cx(|cx| {
        let rho = 1.2;
        let u0: [f64; 3] = core::array::from_fn(|a| [0.022, -0.013, 0.018][a] + boost[a]);
        let u1: [f64; 3] = core::array::from_fn(|a| [0.018, -0.008, 0.021][a] + boost[a]);
        let mut law = porous_law();
        for a in 0..3 { law.linear.skeleton_velocity[a] += boost[a]; }
        let (d0, d1) = (drag(law, u0), drag(law, u1));
        // Manufacture the independent load from the physical trapezoidal
        // momentum equation, rather than reproducing the production root solve.
        law.linear.body_force = core::array::from_fn(|a| rho * (u1[a] - u0[a]) - 0.5 * (d0[a] + d1[a]));
        let grid = Grid3::<L, Q>::try_new_equilibrium(
            [5, 3, 3], collision, [0.0; 3], [false, true, true], budget(),
            |_| false, |_| (rho, u0),
        ).unwrap();
        let flow = OpenBoundaryGrid3::try_new(grid, cx, pressure_faces(rho)).unwrap();
        let mut porous = ResistanceOpenGrid3::try_new_preserving_velocity(flow, cx, vec![law; 45], [0.0; 3]).unwrap();
        let exchange = porous.step_serial(cx, None).unwrap();
        for cell in 0..45 {
            let (density, u) = porous.flow().macroscopic(cell).unwrap();
            near(density, rho, 2e-13);
            for a in 0..3 {
                near(u[a], u1[a], 2e-13);
                near(porous.flow().forces()[cell][a], law.linear.body_force[a] + d1[a], 2e-13);
            }
        }
        let midpoint = core::array::from_fn(|a| 0.5 * (u0[a] + u1[a]));
        let reaction = core::array::from_fn(|a| -45.0 * 0.5 * (d0[a] + d1[a]));
        for a in 0..3 {
            near(exchange.external_impulse[a], 45.0 * law.linear.body_force[a], 2e-13);
            near(exchange.skeleton_impulse[a], reaction[a], 2e-13);
        }
        near(exchange.fluid_drag_work, -dot(reaction, midpoint), 2e-13);
        near(exchange.skeleton_work, dot(reaction, law.linear.skeleton_velocity), 2e-13);
        near(exchange.fluid_drag_work + exchange.skeleton_work, -exchange.relative_dissipation, 2e-13);
        assert!(exchange.relative_dissipation > 0.0);
        assert!(exchange.open.open_links > 0);
    });
}

#[test]
fn anisotropic_vector_speed_drag_has_the_physical_endpoint_and_common_boost() {
    for boost in [[0.0; 3], [0.017, -0.009, 0.005]] {
        prescribed_endpoint::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 }, boost);
        prescribed_endpoint::<D3q27, 27>(bgk(), boost);
    }
}

fn slab<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let dims = [24, 3, 3];
        let (left, right) = (1.0003, 0.9997);
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Pressure { density: left };
        faces[1] = FaceBoundary3::Pressure { density: right };
        let law = DarcyForchheimer3 {
            linear: LinearResistance3 { principal: [0.05; 3], ..LinearResistance3::default() },
            quadratic: 5.0,
        };
        let grid = Grid3::<L, Q>::try_new_equilibrium(
            dims, collision, [0.0; 3], [false, true, true], budget(), |_| false,
            |p| (left + (right - left) * p[0] as f64 / 23.0, [0.0; 3]),
        ).unwrap();
        let flow = OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap();
        let mut porous = ResistanceOpenGrid3::try_new_preserving_velocity(flow, cx, vec![law; 216], [0.0; 3]).unwrap();
        for _ in 0..800 { porous.step_serial(cx, None).unwrap(); }
        let drive = L::CS2 * (left - right) / 23.0;
        let expected = 2.0 * drive / (0.05 + (0.05_f64.powi(2) + 20.0 * drive).sqrt());
        let mean = (2..22).map(|x| porous.flow().macroscopic(address([x, 1, 1], dims)).unwrap().1[0]).sum::<f64>() / 20.0;
        // Small density departure and speed bound the compressible/inertial
        // corrections to the continuum Darcy–Forchheimer slab solution.
        near(mean, expected, 0.03 * expected);
        near(porous.flow().macroscopic(address([0, 1, 1], dims)).unwrap().0, left, 2e-13);
        near(porous.flow().macroscopic(address([23, 1, 1], dims)).unwrap().0, right, 2e-13);
        assert!(porous.last_exchange().unwrap().skeleton_impulse[0] > 0.0);
    });
}

#[test]
fn pressure_driven_porous_slab_reaches_darcy_forchheimer_throughflow() {
    slab::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    slab::<D3q27, 27>(bgk());
}

fn zero_parity<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let dims = [6, 5, 3];
        let solid = [3, 2, 1];
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity { velocity: [0.018, 0.0, 0.004] };
        faces[1] = FaceBoundary3::Pressure { density: 1.1 };
        faces[2] = FaceBoundary3::stationary_wall();
        faces[3] = faces[2];
        let make = || {
            let grid = Grid3::<L, Q>::try_new_equilibrium(
                dims, collision, [0.0; 3], [false, false, true], budget(),
                |p| p == solid, |p| (1.1 + 0.0001 * p[0] as f64, [0.012, 0.0, 0.004]),
            ).unwrap();
            OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
        };
        let laws: Vec<_> = (0..90).map(|cell| DarcyForchheimer3::from(LinearResistance3 {
            body_force: [0.00001 * (cell % 6) as f64, -0.00002, 0.00003],
            skeleton_velocity: [0.01, -0.02, 0.015],
            ..LinearResistance3::default()
        })).collect();
        let forces: Vec<_> = laws.iter().enumerate().map(|(cell, law)| {
            if cell == address(solid, dims) { [0.0; 3] } else { law.linear.body_force }
        }).collect();
        let mut explicit = ForcedOpenBoundaryGrid3::try_new_preserving_velocity(make(), cx, forces).unwrap();
        let mut porous = ResistanceOpenGrid3::try_new_preserving_velocity(make(), cx, laws, [0.5, -0.2, 0.75]).unwrap();
        for _ in 0..4 {
            let ordinary = explicit.step_serial(cx, None, [0.5, -0.2, 0.75]).unwrap();
            let coupled = porous.step_serial(cx, None).unwrap();
            assert_eq!(
                porous.flow().populations().iter().flatten().map(|v| v.to_bits()).collect::<Vec<_>>(),
                explicit.populations().iter().flatten().map(|v| v.to_bits()).collect::<Vec<_>>(),
            );
            assert_eq!(porous.flow().forces(), explicit.forces());
            assert_eq!(coupled.open, ordinary);
            assert_eq!(coupled.skeleton_impulse, [0.0; 3]);
            assert_eq!(coupled.relative_dissipation, 0.0);
            assert!(coupled.walls.measured_links > 0);
        }
    });
}

#[test]
fn zero_resistance_preserves_native_force_density_seams_and_voxel_bits() {
    zero_parity::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    zero_parity::<D3q27, 27>(bgk());
}

fn inventory<L: EquilibriumLattice3<Q>, const Q: usize>(
    porous: &ResistanceOpenGrid3<L, Q>,
    solid: Option<[usize; 3]>,
) -> (f64, [f64; 3], [f64; 3]) {
    let dims = porous.flow().dimensions();
    let origin = porous.moment_origin();
    let (mut mass, mut momentum, mut angular) = (0.0, [0.0; 3], [0.0; 3]);
    for cell in 0..porous.flow().populations().len() {
        let p = [cell % dims[0], cell / dims[0] % dims[1], cell / (dims[0] * dims[1])];
        if Some(p) == solid { continue; }
        let (rho, u) = porous.flow().macroscopic(cell).unwrap();
        let arm: [f64; 3] = core::array::from_fn(|a| p[a] as f64 + 0.5 - origin[a]);
        mass += rho;
        for a in 0..3 {
            let (b, c) = ((a + 1) % 3, (a + 2) % 3);
            momentum[a] += rho * u[a];
            angular[a] += arm[b] * rho * u[c] - arm[c] * rho * u[b];
        }
    }
    (mass, momentum, angular)
}

fn balance<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let dims = [7, 5, 4];
        let solid = [3, 2, 1];
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity { velocity: [0.015, 0.0, 0.004] };
        faces[1] = FaceBoundary3::Pressure { density: 1.1 };
        faces[2] = FaceBoundary3::Wall { velocity: [0.0, 0.0, 0.006] };
        faces[3] = faces[2];
        let grid = Grid3::<L, Q>::try_new_equilibrium(
            dims, collision, [0.0; 3], [false, false, true], budget(), |p| p == solid,
            |p| (1.1 + 0.0001 * p[0] as f64, [0.012, 0.001 * (p[1] as f64 - 2.0), 0.003]),
        ).unwrap();
        let flow = OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap();
        let laws: Vec<_> = (0..140).map(|cell| {
            let (x, y) = (cell % 7, cell / 7 % 5);
            let mut law = porous_law();
            law.linear.body_force = [0.00005, -0.00003, 0.00002];
            if (x == 0 || x == 6) && (y == 0 || y == 4) {
                law.linear.principal = [0.0; 3];
                law.quadratic = 0.0;
            }
            law
        }).collect();
        let mut porous = ResistanceOpenGrid3::try_new_preserving_velocity(flow, cx, laws, [1.25, -0.75, 0.5]).unwrap();
        for _ in 0..3 {
            let before = inventory(&porous, Some(solid));
            let r = porous.step_serial(cx, None).unwrap();
            let after = inventory(&porous, Some(solid));
            near(after.0 - before.0, r.open.mass_change + r.walls.fluid_mass_change, 2e-11);
            for a in 0..3 {
                near(after.1[a] - before.1[a],
                    r.external_impulse[a] - r.skeleton_impulse[a]
                        + r.open.impulse[a] + r.walls.fluid_population_impulse[a], 2e-11);
                near(after.2[a] - before.2[a],
                    r.external_angular_impulse[a] - r.skeleton_angular_impulse[a]
                        + r.open.angular_impulse[a] + r.walls.fluid_population_angular_impulse[a]
                        + r.periodic.angular_momentum_change[a], 5e-11);
                near(r.walls.wall_impulse[a] + r.walls.fluid_population_impulse[a],
                    r.walls.wall_velocity_mass_impulse[a], 2e-13);
            }
            near(r.fluid_drag_work + r.skeleton_work, -r.relative_dissipation, 2e-13);
            assert!(r.walls.measured_links > porous.flow().last_wall_exchange().unwrap().measured_links);
            assert!(r.periodic.face_crossings > 0);
            assert!(r.open.open_links > 0);
            assert_eq!(porous.last_exchange(), Some(r));
        }
    });
}

#[test]
fn moving_wall_voxel_reservoir_and_periodic_terms_close_the_actual_momentum_account() {
    balance::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    balance::<D3q27, 27>(bgk());
}

fn traces<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let dims = [6, 4, 3];
        let grid = Grid3::<L, Q>::try_new_equilibrium(
            dims, collision, [0.0; 3], [false, true, true], budget(),
            |_| false, |_| (1.1, [0.025, -0.002, 0.003]),
        ).unwrap();
        let flow = OpenBoundaryGrid3::try_new(grid, cx, boundaries(0, [0.025, -0.002, 0.003], 1.1)).unwrap();
        let law = porous_law();
        let mut porous = ResistanceOpenGrid3::try_new_preserving_velocity(flow, cx, vec![law; 72], [0.0; 3]).unwrap();
        let gate = CancelGate::new_clock_free();
        for tick in 0..2 {
            let profile = |p: [usize; 3]| [0.023 + 0.001 * tick as f64 + 0.0002 * p[1] as f64, -0.002, 0.003 + 0.0001 * p[2] as f64];
            porous.set_velocity_profile(Face3::XMin, InletFlux3::Unconstrained, &gate, |p| Ok(profile(p))).unwrap();
            porous.set_pressure_profile(Face3::XMax, &gate, |p| Ok(1.1 + 0.0001 * p[1] as f64)).unwrap();
            porous.step_serial(cx, None).unwrap();
            for z in 0..3 {
                for y in 0..4 {
                    let inlet = address([0, y, z], dims);
                    let outlet = address([5, y, z], dims);
                    let (_, u) = porous.flow().macroscopic(inlet).unwrap();
                    for a in 0..3 { near(u[a], profile([0, y, z])[a], 2e-13); }
                    near(porous.flow().macroscopic(outlet).unwrap().0, 1.1 + 0.0001 * y as f64, 2e-13);
                }
            }
            for cell in 0..72 {
                let u = porous.flow().macroscopic(cell).unwrap().1;
                let d = drag(law, u);
                for a in 0..3 { near(porous.flow().forces()[cell][a], law.linear.body_force[a] + d[a], 2e-13); }
            }
        }
        porous.clear_pressure_profile(Face3::XMax);
        porous.set_characteristic_outlet(Face3::XMax, CharacteristicOutlet3 {
            reference_normal_velocity_lu: 0.025,
            incoming_characteristic_lu: 0.0,
            maximum_relative_density_departure: 0.05,
        }, &gate).unwrap();
        porous.step_serial(cx, None).unwrap();
        assert!(porous.flow().characteristic_outlet(Face3::XMax).is_some());
        for z in 0..3 {
            for y in 0..4 {
                let (rho, u) = porous.flow().macroscopic(address([5, y, z], dims)).unwrap();
                near((u[0] - 0.025) - L::CS2.sqrt() * (rho - 1.1) / 1.1, 0.0, 2e-13);
            }
        }
    });
}

#[test]
fn updated_profiles_and_characteristic_traces_use_their_own_implicit_source() {
    traces::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    traces::<D3q27, 27>(bgk());
}

fn snapshot<L: EquilibriumLattice3<Q>, const Q: usize>(
    grid: &ResistanceOpenGrid3<L, Q>,
) -> (Vec<u64>, Vec<u64>, usize, Option<ResistanceMomentumExchange3>) {
    (
        grid.flow().populations().iter().flatten().map(|v| v.to_bits()).collect(),
        grid.flow().forces().iter().flatten().map(|v| v.to_bits()).collect(),
        grid.flow().completed_steps(),
        grid.last_exchange(),
    )
}

fn seam_refusal<L: EquilibriumLattice3<Q>, const Q: usize>(collision: L::Collision) {
    with_cx(|cx| {
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[0] = FaceBoundary3::Velocity { velocity: [0.01, 0.0, 0.0] };
        faces[1] = FaceBoundary3::Pressure { density: 1.0 };
        faces[2] = FaceBoundary3::stationary_wall();
        faces[3] = faces[2];
        let make = || {
            let grid = Grid3::<L, Q>::try_new_equilibrium(
                [5, 4, 3], collision, [0.0; 3], [false, false, true], budget(),
                |_| false, |_| (1.0, [0.01, 0.0, 0.0]),
            ).unwrap();
            OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap()
        };
        for law in [
            DarcyForchheimer3::from(LinearResistance3 {
                principal: [0.02; 3], ..LinearResistance3::default()
            }),
            DarcyForchheimer3 { quadratic: 0.5, ..DarcyForchheimer3::default() },
        ] {
            assert!(matches!(
                ResistanceOpenGrid3::try_new_preserving_velocity(make(), cx, vec![law; 60], [0.0; 3]),
                Err(GridError3::MaterialField { cell: 0, .. })
            ));
        }
        // Nonzero laws are admitted in the interior and on complete face
        // traces when only the mixed boundary cells use zero drag.
        let mut laws = vec![porous_law(); 60];
        for cell in 0..60 {
            let (x, y) = (cell % 5, cell / 5 % 4);
            if (x == 0 || x == 4) && (y == 0 || y == 3) {
                laws[cell] = DarcyForchheimer3::default();
            }
        }
        let mut admitted = ResistanceOpenGrid3::try_new_preserving_velocity(make(), cx, laws, [0.0; 3]).unwrap();
        admitted.step_serial(cx, None).unwrap();
    });
}

fn rollback_fixture(cx: &Cx<'_>) -> ResistanceOpenGrid3<D3q27, 27> {
    let mut faces = [FaceBoundary3::Periodic; 6];
    faces[4] = FaceBoundary3::Velocity { velocity: [0.0; 3] };
    faces[5] = faces[4];
    let grid = Grid3::<D3q27, 27>::try_new_equilibrium(
        [8, 8, 5], bgk(), [0.0; 3], [true, true, false], budget(),
        |_| false, |_| (1.0, [0.0; 3]),
    ).unwrap();
    let flow = OpenBoundaryGrid3::try_new(grid, cx, BoundarySpec3::new(faces)).unwrap();
    let law = DarcyForchheimer3 {
        linear: LinearResistance3 { principal: [0.04; 3], ..LinearResistance3::default() },
        quadratic: 1.0,
    };
    ResistanceOpenGrid3::try_new_preserving_velocity(flow, cx, vec![law; 320], [0.5, -0.25, 1.0]).unwrap()
}

#[test]
fn active_mixed_seams_refuse_and_late_endpoint_failure_preserves_all_epochs() {
    seam_refusal::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.9 });
    seam_refusal::<D3q27, 27>(bgk());
    with_cx(|cx| {
        let mut porous = rollback_fixture(cx);
        porous.step_serial(cx, None).unwrap();
        let mut faces = [FaceBoundary3::Periodic; 6];
        faces[4] = FaceBoundary3::Velocity { velocity: [0.0, 0.0, 0.03] };
        faces[5] = faces[4];
        let gate = CancelGate::new_clock_free();
        porous.rebind_open_boundaries(BoundarySpec3::new(faces), &gate).unwrap();
        let before = snapshot(&porous);
        let mut collisions = vec![bgk(); 320];
        collisions[319].maximum_mach = 0.01; // Opening rest passes; the complete trace does not.
        assert!(matches!(
            porous.step_serial(cx, Some(&collisions)),
            Err(GridError3::Mach { cell: 319, .. })
        ));
        assert_eq!(snapshot(&porous), before);
        assert_eq!(porous.flow().boundaries(), BoundarySpec3::new(faces));
        porous.step_serial(cx, None).unwrap();
        assert_eq!(porous.flow().completed_steps(), 2);
    });
}

#[test]
fn pooled_resistance_open_steps_match_serial_and_cancel_after_stream_atomically() {
    struct CancelAfterStream {
        pool: TilePool,
        calls: AtomicUsize,
    }
    impl KernelRunner for CancelAfterStream {
        fn workers(&self) -> usize { self.pool.workers() }
        fn run_with_gate<K: TileKernel>(
            &self, kernel: &K, gate: &CancelGate,
        ) -> (Result<K::Out, RunError>, RunReport) {
            let result = self.pool.run_with_gate(kernel, gate);
            if self.calls.fetch_add(1, Ordering::AcqRel) == 1 { gate.request(); }
            result
        }
    }
    with_cx(|cx| {
        let setup = || {
            let mut owner = rollback_fixture(cx);
            let mut faces = [FaceBoundary3::Periodic; 6];
            faces[4] = FaceBoundary3::Velocity { velocity: [0.004, -0.002, 0.03] };
            faces[5] = faces[4];
            owner.rebind_open_boundaries(BoundarySpec3::new(faces), &CancelGate::new_clock_free()).unwrap();
            owner
        };
        let mut serial = setup();
        let expected = serial.step_serial(cx, None).unwrap();
        for workers in [1, 2, 7] {
            let mut parallel = setup();
            let actual = parallel.step_pooled(
                &TilePool::for_host(workers, SEED), &CancelGate::new_clock_free(), None,
            ).unwrap();
            assert_eq!(snapshot(&parallel), snapshot(&serial));
            assert_eq!(actual, expected);
        }
        let mut cancelled = setup();
        let before = snapshot(&cancelled);
        let runner = CancelAfterStream {
            pool: TilePool::for_host(2, SEED),
            calls: AtomicUsize::new(0),
        };
        assert!(matches!(
            cancelled.step_pooled(&runner, &CancelGate::new_clock_free(), None),
            Err(GridError3::Cancelled)
        ));
        assert_eq!(snapshot(&cancelled), before);
        let retry = cancelled.step_serial(cx, None).unwrap();
        assert_eq!(retry, expected);
        assert_eq!(snapshot(&cancelled), snapshot(&serial));
    });
}
