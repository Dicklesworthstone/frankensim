//! Original gas, material, reference mass and boundary pressure in one solve.
use super::gas_fixture::gas;
use super::*;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_solid::hyper3::gravity::{GravitySample3, ReferenceGravityDensity, StationaryGravity3};
use fs_solid::hyper3::pressure::cavity::AdiabaticCavity;

#[test]
fn g1_adiabatic_cavity_gravity_and_ambient_pressure_share_native_equilibrium_and_reactions() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let (mu, lambda, stretch) = (2.0e6, 3.0e6, 0.94_f64);
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let piola = mu * (stretch - 1.0 / stretch) + 3.0 * lambda * stretch.ln() / stretch;
        let ambient = 2.0e5;
        let target_pressure = -piola / stretch.powi(2);
        // This original NASA-9 fixture has constant Cp/R=3.5, hence gamma=1.4.
        // The independent ideal-gas isentrope supplies the reference pressure,
        // while production still resolves its original source laws at every trial.
        let gamma = 1.4;
        let reference_pressure = target_pressure * stretch.powf(3.0 * gamma);
        let gas = gas(cavity.reference_volume_m3(), reference_pressure, 0.0, 1e24);
        let gas_before = gas.clone();
        let controls = GasEntropySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 600.0,
            entropy_tolerance_j_kg_k: 1e-10,
            maximum_iterations: 128,
        };
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let density: Vec<_> = (0..tets.len()).map(|i| 1500.0 + 250.0 * i as f64).collect();
        let constant = [0.1, -0.2, -9.8];
        let gradient = [-0.05, -0.08, -0.02];
        let acceleration = |x: [f64; 3]| std::array::from_fn(|i| constant[i] + gradient[i] * x[i]);
        let field = |x: [f64; 3]| {
            Ok(GravitySample3 {
                potential_j_kg: 1e20
                    - (0..3)
                        .map(|i| constant[i] * x[i] + 0.5 * gradient[i] * x[i] * x[i])
                        .sum::<f64>(),
                acceleration_m_s2: acceleration(x),
                acceleration_gradient_s2: std::array::from_fn(|i| {
                    std::array::from_fn(|j| if i == j { gradient[i] } else { 0.0 })
                }),
            })
        };
        let increment_calls = std::cell::Cell::new(0);
        let increment = |a: [f64; 3], b: [f64; 3]| {
            increment_calls.set(increment_calls.get() + 1);
            Ok(-(0..3)
                .map(|i| (b[i] - a[i]) * (constant[i] + 0.5 * gradient[i] * (a[i] + b[i])))
                .sum::<f64>())
        };

        // Exact simplex moments of rho*N_a*N_b independently assemble the
        // affine gravity load at the known target. Current-position gravity
        // differs during Newton, so this is not a frozen dead-force substitute.
        let mut gravity_at_target = vec![0.0; 3 * nodes.len()];
        for (tet, &rho) in tets.iter().zip(&density) {
            let x = tet.map(|node| nodes[node]);
            let reference_volume =
                dot(sub(x[1], x[0]), cross(sub(x[2], x[0]), sub(x[3], x[0]))).abs() / 6.0;
            let g = x.map(|point| acceleration(point.map(|value| stretch * value)));
            for a in 0..4 {
                for axis in 0..3 {
                    gravity_at_target[3 * tet[a] + axis] += rho * reference_volume / 20.0
                        * ((0..4).map(|b| g[b][axis]).sum::<f64>() + g[a][axis]);
                }
            }
        }
        let (_, outer_gradient) = tetra_volume_gradient([nodes[0], nodes[1], nodes[2], nodes[3]]);
        let outer_traction = piola + ambient * stretch.powi(2);
        let dead: Vec<_> = (0..3 * nodes.len())
            .filter(|dof| !fixed.iter().any(|&(id, _)| id == *dof))
            .map(|dof| {
                let outer = if dof < 12 {
                    outer_traction * outer_gradient[dof]
                } else {
                    0.0
                };
                (dof, outer - gravity_at_target[dof])
            })
            .collect();
        let problem = HyperTetProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            material: HyperelasticLawField::Uniform(&material),
            prescribed_m: &fixed,
            nodal_forces_n: &dead,
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings {
                load_steps: 8,
                ..HyperTetSettings::default()
            },
        };
        // Ambient pressure acts on the entire closed material boundary. Its
        // inner contribution cancels the SAME ambient counterpressure in the
        // cavity load: gas pressure acts once inside, ambient once outside.
        let solution = problem
            .solve_with_adiabatic_cavity_loads(
                AdiabaticCavity {
                    cavity: &cavity,
                    reference: &gas,
                    exterior_pressure_pa: ambient,
                    controls,
                },
                Some((
                    ReferenceGravityDensity::PerElement(&density),
                    StationaryGravity3::with_potential_increment(&field, &increment),
                )),
                Some((cavity.surface(), TetPressure::Uniform(ambient))),
                cx,
            )
            .unwrap();
        assert!(increment_calls.get() > 0);
        assert_eq!(
            gas, gas_before,
            "equilibrium cannot mutate the source parcel"
        );
        let u: Vec<_> = solution.displacement_m.iter().flatten().copied().collect();
        for (&actual, &reference) in u.iter().zip(nodes.iter().flatten()) {
            near(actual, (stretch - 1.0) * reference, 4e-8);
        }
        near(solution.min_det_f, stretch.powi(3), 4e-8);
        let load = cavity
            .adiabatic_gas_load(&u, &gas, ambient, controls, cx)
            .unwrap();
        near(load.volume_m3, stretch.powi(3) / 6.0, 3e-9);
        near(load.net_pressure_pa + ambient, target_pressure, 0.1);
        near(
            load.gas.state().state().temperature,
            300.0 * stretch.powf(-3.0 * (gamma - 1.0)),
            3e-6,
        );
        assert_eq!(load.gas.mass_kg().to_bits(), gas.mass_kg().to_bits());
        assert_eq!(load.gas.model(), gas.model());
        for &(dof, reaction) in &solution.reactions_n {
            near(
                reaction,
                outer_traction * outer_gradient[dof] - gravity_at_target[dof],
                0.15,
            );
        }
        let gravity = problem.gravity_load(&density, &u, &field, cx).unwrap();
        let pressure_force = cavity
            .surface()
            .force_n(&u, TetPressure::Uniform(ambient), cx)
            .unwrap();
        for axis in 0..3 {
            let reaction: f64 = solution
                .reactions_n
                .iter()
                .filter(|(dof, _)| dof % 3 == axis)
                .map(|(_, value)| value)
                .sum();
            let dead_force: f64 = dead
                .iter()
                .filter(|(dof, _)| dof % 3 == axis)
                .map(|(_, value)| value)
                .sum();
            let external: f64 = (0..nodes.len())
                .map(|node| {
                    let dof = 3 * node + axis;
                    gravity.force_n[dof] + pressure_force[dof] + load.force_n[dof]
                })
                .sum();
            near(reaction + dead_force + external, 0.0, 0.15);
        }
    });
}
