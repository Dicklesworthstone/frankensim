//! Release a sheared, spinning quadratic tetrahedron without artificial damping.
//! Run `cargo run -p fs-solid --example conservative_elastic_motion`.
//!
//! Both runs start from the same volume-preserving deformation and velocity.
//! FullMaterial retains the finite bulk stiffness. ExactIncompressible enforces
//! continuous P1 weak volume constraints and never evaluates its unused bulk
//! modulus. The pressure printed for that run is an interval work multiplier,
//! not an endpoint Cauchy pressure. Weak volume moments are not pointwise J=1.
//!
//! All quantities use SI units: a 1 m reference edge, density 6 kg/m^3, and
//! Mooney-Rivlin coefficients c10=4 Pa, c01=2 Pa give a deliberately soft demo.
//! This conservative, second-order integrator is not frequency matched; audio
//! rendering still requires the separate exact/frequency-matched linear split.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::discrete_gradient::DiscreteGradientSettings;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::mixed::conservative::{ConservativeFormulation, ConservativeTetDefinition};
use fs_solid::hyper3::mixed::{MixedLawField, MixedUpSettings};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::tet10::promote_p1_to_p2;
use fs_solver::nonlinear::NewtonKrylovConfig;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cancel = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &cancel,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 24,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let (nodes, elements) = promote_p1_to_p2(
            &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
            &[[0, 1, 2, 3]],
            None,
            &cx,
        )?;
        // A nonaffine triangular map: its determinant is exactly one.
        let initial_u: Vec<_> = nodes
            .iter()
            .map(|x| [0.25 * x[1] + 0.4 * x[1] * x[2], 0.0, 0.0])
            .collect();
        let initial_v: Vec<_> = nodes
            .iter()
            .zip(&initial_u)
            .map(|(x, u)| {
                let [a, b, c] = core::array::from_fn(|i| x[i] + u[i] - 0.25);
                // Translation + angular velocity (0.3,-0.2,0.5) rad/s +
                // a trace-free strain rate, so div_x(v)=0 initially.
                [
                    0.2 - 0.2 * c - 0.5 * b + 0.08 * a,
                    -0.15 + 0.5 * a - 0.3 * c - 0.06 * b,
                    0.1 + 0.3 * b + 0.2 * a - 0.02 * c,
                ]
            })
            .collect();
        let settings = MixedUpSettings {
            load_steps: 1,
            maximum_newton_iterations: 30,
            newton: NewtonKrylovConfig {
                absolute_tolerance: 2e-11,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            length_scale_m: 1.0,
            stress_scale_pa: 12.0,
            maximum_assembly_bytes: 2 * 1024 * 1024,
            maximum_factor_bytes: 2 * 1024 * 1024,
        };
        let formulations = [
            ("full", ConservativeFormulation::FullMaterial, 60.0),
            (
                "exact",
                ConservativeFormulation::ExactIncompressible {
                    maximum_relative_volume_residual: 1e-9,
                    maximum_initial_relative_volume_rate_per_s: 1e-9,
                    maximum_point_volume_work_defect: 1e-12,
                },
                1e300,
            ),
        ];
        println!(
            "mode,time_s,corner_x_m,stored_j,kinetic_j,energy_drift_j,interval_balance_j,min_det_f,max_weak_volume_m3,max_interval_pressure_pa"
        );
        for (name, formulation, kappa) in formulations {
            let material = Hyperelastic::new(
                HyperelasticModel::MooneyRivlin {
                    c10: 4.0,
                    c01: 2.0,
                    kappa,
                },
                4.0,
            )?;
            let density = [6.0];
            let mut solid = ConservativeTetDefinition {
                nodes_m: &nodes,
                elements: &elements,
                material: MixedLawField::Uniform(&material),
                reference_density_kg_m3: &density,
                formulation,
                budget: TetAssemblyBudget::standard(),
                quadrature_points_per_dir: 4,
                maximum_quadrature_points: 64,
                maximum_mass_assembly_bytes: 128 * 1024,
                material_settings: DiscreteGradientSettings {
                    minimum_relative_strain_increment: 1e-6,
                    absolute_work_tolerance_pa: 1e-10,
                    relative_work_tolerance: 1e-12,
                },
            }
            .prepare(&initial_u, &initial_v, &cx)?;
            let initial_energy_j = solid.total_energy_j();
            for interval in 1..=32 {
                // No support, applied force or body force: only elastic motion.
                // Each accepted interval and the accumulated energy balance
                // must close within 2e-9 J; a refusal leaves solid unchanged.
                let step = solid.advance(0.01, &[], &[], [0.; 3], settings, 2e-9, &cx)?;
                if interval == 1 || interval % 8 == 0 {
                    let exact = !step.pressure_residual_m3.is_empty();
                    let weak_volume = exact.then(|| {
                        step.pressure_residual_m3.iter().map(|v| v.abs()).fold(0.0, f64::max)
                    });
                    let pressure = exact.then(|| {
                        step.solid.algorithmic_pressure_pa().iter().map(|v| v.abs()).fold(0.0, f64::max)
                    });
                    println!(
                        "{name},{:.3},{:.9},{:.9},{:.9},{:.3e},{:.3e},{:.9},{},{}",
                        step.solid.time_s(),
                        nodes[1][0] + step.solid.displacement_m()[1][0],
                        step.solid.stored_energy_j(),
                        step.solid.kinetic_energy_j(),
                        step.solid.total_energy_j() - initial_energy_j,
                        step.energy.balance_error_j,
                        step.min_det_f,
                        weak_volume.map(|v| format!("{v:.3e}")).unwrap_or_default(),
                        pressure.map(|v| format!("{v:.6}")).unwrap_or_default(),
                    );
                }
                solid = step.solid;
            }
        }
        Ok(())
    })
}
