//! SI constitutive stress on the actual native moving-wall/force stencil.
//! Analytic profiles are prescribed here; this is not a flow-convergence test.

use std::cell::Cell;

use fs_exec::{CancelGate, Cx};
use fs_lbm::d3q19::CollisionModel3;
use fs_lbm::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27};
use fs_lbm::kinematics::{FlowKinematics3, KinematicsWorkspace3};
use fs_lbm::lattice3::sources::{ForcedGrid3, ForcedLattice3};
use fs_lbm::lattice3::{D3q19, Grid3, GridBudget3, Lattice3};
use fs_lbm::material::{MaterialLbmError, plan_scaling_from_material};
use fs_material::MaterialError;
use fs_material::fluid_transport::{
    FluidTransport, FluidTransportState,
    rheology::{CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode},
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

const DX: f64 = 0.01;
const DT: f64 = 0.01;
const MU: f64 = 1.0;

struct SingleRead(Cell<usize>);

impl FluidTransport for SingleRead {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        assert_eq!(self.0.replace(self.0.get() + 1), 0, "source was read twice");
        FluidTransportState::try_new(1000.0, 1000.0, 1000.0, MU)
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e} (tol {tolerance})"
    );
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let arenas = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    arenas.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey { seed: 0x5E4A, kernel_id: 3, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        ))
    })
}

fn equilibrium<L: Lattice3<Q>, const Q: usize>(rho: f64, u: [f64; 3]) -> [f64; Q] {
    let u2: f64 = u.iter().map(|v| v * v).sum();
    core::array::from_fn(|q| {
        let eu: f64 = (0..3).map(|a| f64::from(L::E[q][a]) * u[a]).sum();
        L::W[q] * rho * (1.0 + eu / L::CS2 + 0.5 * eu * eu / L::CS2.powi(2)
            - 0.5 * u2 / L::CS2)
    })
}

fn thermal_law() -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(3.0, 0.2, 4.0, 0.4, 2.0, [0.0, 3.0]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius { activation_energy_j_mol: 30_000.0 },
            300.0,
            [290.0, 330.0],
        ).unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}

fn native_stress<L: Lattice3<Q> + ForcedLattice3<Q>, const Q: usize>(collision: L::Collision) {
    let source = SingleRead(Cell::new(0));
    let scaling = plan_scaling_from_material(&source, DX, DT, 0.05).unwrap();
    let gate = CancelGate::new_clock_free();
    for normal in 0..3 {
        for tangent in (0..3).filter(|&a| a != normal) {
            let mut dims = [3; 3];
            dims[normal] = 7;
            let cells = dims.into_iter().product();
            for (linear, quadratic) in [(0.0, 0.0), (0.002, 0.0), (0.0, 0.0003)] {
                for forced in [false, true] {
                    let force = |p: [usize; 3]| {
                        let mut f = [0.0; 3];
                        if forced { f[tangent] = 0.0004 * (p[normal] as f64 + 0.5); }
                        f
                    };
                    let grid = Grid3::<L, Q>::try_new(
                        dims, collision, [0.0; 3],
                        core::array::from_fn(|a| a != normal),
                        GridBudget3 { maximum_cells: cells, maximum_mach: 0.2 },
                        |_| false,
                        |p| {
                            let y = p[normal] as f64 + 0.5;
                            let rho = 1.4 + 0.002 * y;
                            let mut u = force(p).map(|f| -0.5 * f / rho);
                            u[tangent] += 0.01 + linear * y + quadratic * y * (7.0 - y);
                            equilibrium::<L, Q>(rho, u)
                        },
                    ).unwrap();
                    let before = grid.populations().to_vec();
                    let probe = with_cx(|cx| grid.wall_momentum_probe(cx, [0.0; 3], |_| true).unwrap());
                    let wall: Vec<_> = probe.wall_links().map(|link| {
                        let mut u = [0.0; 3];
                        u[tangent] = 0.01 + if link.exterior_faces[2 * normal + 1] {
                            linear * 7.0
                        } else { 0.0 };
                        u
                    }).collect();
                    let mut workspace = KinematicsWorkspace3::new(
                        dims, KinematicsWorkspace3::required_bytes(dims).unwrap(),
                    ).unwrap();
                    // Retain the current force frame on the original owner.
                    let forces: Vec<_> = (0..cells).map(|i| {
                        force([i % dims[0], i / dims[0] % dims[1], i / dims[0] / dims[1]])
                    }).collect();
                    let grid = with_cx(|cx| ForcedGrid3::try_new(grid, cx, forces.clone()).unwrap());
                    let field = grid.kinematics_with_moving_walls(
                        &gate, &probe, &wall, &mut workspace,
                    ).unwrap();
                    for (i, sample) in field.iter().enumerate() {
                        let p = [i % dims[0], i / dims[0] % dims[1], i / dims[0] / dims[1]];
                        let derivative = (linear + quadratic * (7.0 - 2.0 * (p[normal] as f64 + 0.5))) / DT;
                        let kinematics = sample.unwrap();
                        for bulk in [0.0, 7.0] {
                            let response = scaling.viscous_stress_from_lattice_kinematics(
                                kinematics, bulk,
                            ).unwrap();
                            let mut expected = [0.0; 6];
                            expected[match (normal.min(tangent), normal.max(tangent)) {
                                (0, 1) => 3, (1, 2) => 4, _ => 5,
                            }] = MU * derivative;
                            for (&actual, expected) in response.stress_pa().iter().zip(expected) {
                                near(actual, expected, 3e-12);
                            }
                            near(response.dissipation_w_m3(), MU * derivative.powi(2), 3e-12);
                            for sign in [-1.0, 1.0] {
                                let mut n = [0.0; 3]; n[normal] = sign;
                                let traction = response.traction(n).unwrap();
                                for axis in 0..3 {
                                    near(traction.shear_pa[axis], if axis == tangent {
                                        sign * MU * derivative
                                    } else { 0.0 }, 3e-12);
                                }
                                near(traction.normal_pa, 0.0, 3e-12);
                            }
                        }
                        if forced && linear != 0.0 && p[normal] == 3 {
                            let law = thermal_law();
                            let reference = scaling.thermal_viscous_stress_from_lattice_kinematics(
                                kinematics, &law, 300.0, 0.0,
                            ).unwrap();
                            let expected_mu = 0.2 + 2.8 * (1.0 + (4.0 * derivative).powi(2)).powf(-0.3);
                            near(reference.viscous.shear_viscosity_pa_s(), expected_mu, 3e-12);
                            near(reference.rheology.shear_rate_per_s(), derivative, 3e-12);
                            assert_eq!(reference.rheology.law(), &law);
                            let hot = scaling.thermal_viscous_stress_from_lattice_kinematics(
                                kinematics, &law, 320.0, 0.0,
                            ).unwrap();
                            assert!(hot.viscous.dissipation_w_m3() < reference.viscous.dissipation_w_m3());
                            assert!(scaling.thermal_viscous_stress_from_lattice_kinematics(
                                kinematics, &law, 340.0, 0.0,
                            ).is_err());
                        }
                    }
                    assert_eq!(grid.populations(), before);
                    assert_eq!(grid.forces(), forces);
                    assert_eq!(grid.completed_steps(), 0);
                }
            }
        }
    }
    assert_eq!(source.0.get(), 1);
}

#[test]
fn d3q19_material_traction_and_heat_match_all_axis_couette_and_poiseuille() {
    native_stress::<D3q19, 19>(CollisionModel3::Bgk { tau: 0.8 });
}

#[test]
fn d3q27_material_traction_and_heat_match_all_axis_couette_and_poiseuille() {
    native_stress::<D3q27, 27>(CollisionConfig27 {
        model: CollisionModel27::Bgk { tau: 0.8, equilibrium: Equilibrium27::Factorized },
        maximum_mach: 0.2,
    });
}

#[test]
fn rotation_dilation_and_source_refusals_retain_their_physical_meaning() {
    let source = SingleRead(Cell::new(0));
    let scaling = plan_scaling_from_material(&source, DX, DT, 0.05).unwrap();
    let rotation = FlowKinematics3::from_gradient([
        [0.0, -0.3, 0.1], [0.3, 0.0, -0.2], [-0.1, 0.2, 0.0],
    ]).unwrap();
    let no_stress = scaling.viscous_stress_from_lattice_kinematics(rotation, 4.0).unwrap();
    assert_eq!(no_stress.stress_pa(), &[0.0; 6]);
    assert_eq!(no_stress.dissipation_w_m3(), 0.0);
    let dilation = FlowKinematics3::from_gradient([
        [0.02, -0.3, 0.1], [0.3, 0.02, -0.2], [-0.1, 0.2, 0.02],
    ]).unwrap();
    let response = scaling.thermal_viscous_stress_from_lattice_kinematics(
        dilation, &thermal_law(), 300.0, 4.0,
    ).unwrap();
    assert_eq!(response.rheology.shear_rate_per_s(), 0.0);
    assert_eq!(response.viscous.shear_viscosity_pa_s(), 3.0);
    assert_eq!(response.viscous.deviatoric_stress_pa(), &[0.0; 6]);
    near(response.viscous.dissipation_w_m3(), 4.0 * 6.0_f64.powi(2), 1e-12);
    let n = [0.6, 0.8, 0.0];
    assert_eq!(response.viscous.traction(n).unwrap().shear_pa, [0.0; 3]);
    for bulk in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(scaling.viscous_stress_from_lattice_kinematics(rotation, bulk),
            Err(MaterialLbmError::Transport(_))));
    }
    let excessive_rate = FlowKinematics3::from_gradient([
        [0.0, 0.1, 0.0], [0.0; 3], [0.0; 3],
    ]).unwrap();
    assert!(scaling.thermal_viscous_stress_from_lattice_kinematics(
        excessive_rate, &thermal_law(), 300.0, 0.0,
    ).is_err());
    assert_eq!(source.0.get(), 1);
}
