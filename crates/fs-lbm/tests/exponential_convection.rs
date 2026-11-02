//! Exact convection-diffusion profiles through the real CHT assembly/solve.

use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    ConvectionScheme, EnergyConfig, FlowFace, FlowField, FluidProperties, ThermalFace,
    ThermalSetup, VoxelDomain, solve_energy,
};

#[test]
fn exponential_convection_matches_the_boundary_value_problem_on_coarse_graded_grids() {
    let gate = CancelGate::new();
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    };
    let domains = [
        VoxelDomain::graded([vec![0.125; 8], vec![0.2], vec![0.3]]).unwrap(),
        VoxelDomain::graded([vec![0.03, 0.17, 0.06, 0.24, 0.5], vec![0.2], vec![0.3]]).unwrap(),
    ];
    for domain in &domains {
        for speed in [-20.0_f64, -2.0, 0.0, 2.0, 20.0] {
            let mut boundaries = [FlowFace::Wall; 6];
            boundaries[0] = FlowFace::Fixed;
            boundaries[1] = FlowFace::Fixed;
            let (flow, _) = FlowField::from_cell_velocities(
                domain,
                &vec![[speed, 0.0, 0.0]; domain.cell_count()],
                boundaries,
                1e-12,
                &gate,
            )
            .unwrap();
            // The Inflow rule is a diffusive Dirichlet boundary on either
            // flow sign. Prescribe both end temperatures to resolve the
            // downstream thermal boundary layer, rather than zero gradient.
            let mut faces = [ThermalFace::Adiabatic; 6];
            faces[0] = ThermalFace::Inflow { temperature: 300.0 };
            faces[1] = ThermalFace::Inflow { temperature: 400.0 };
            let solution = solve_energy(
                domain,
                &fluid,
                &[],
                &flow,
                &ThermalSetup::new(faces),
                &EnergyConfig {
                    scheme: ConvectionScheme::Exponential,
                    ..EnergyConfig::default()
                },
                &gate,
            )
            .unwrap();
            for cell in 0..domain.cell_count() {
                let [i, j, k] = domain.coords(cell);
                let x = domain.center(i, j, k)[0];
                let profile = if speed == 0.0 {
                    x
                } else {
                    (speed * x).exp_m1() / speed.exp_m1()
                };
                let expected = 300.0 + 100.0 * profile;
                let computed = solution.temperature[cell];
                assert!(
                    (computed - expected).abs() < 2e-7,
                    "speed={speed}, x={x}: {computed} != {expected}"
                );
                assert!((300.0 - 2e-7..=400.0 + 2e-7).contains(&computed));
            }
            assert!(solution.report.balance.residual_w.abs() < 2e-7 * (1.0 + speed.abs()));
        }
    }
}
