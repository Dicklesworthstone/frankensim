use super::*;
use crate::conjugate::simple::{FvBoundary, SimpleConfig, TimeScheme};
use crate::conjugate::{FvSgsConfig, VoxelDomain};
use crate::sgs3::{SgsModel, SgsUse};
use fs_exec::CancelGate;

fn sampler<'a>(
    domain: &'a VoxelDomain,
    config: &'a SimpleConfig,
    density: &[f64],
    molecular: &[f64],
) -> Solver<'a> {
    let gate = CancelGate::new_clock_free();
    let mut solver = Solver::with_molecular_coefficients(domain, density[0], molecular[0], config);
    solver
        .set_material_transport(molecular, None, &gate)
        .unwrap();
    let mass: Vec<_> = density
        .iter()
        .enumerate()
        .map(|(c, rho)| rho * domain.volume(c))
        .collect();
    solver.set_sgs_density(&mass, &gate).unwrap();
    solver.molecular_bulk_viscosity = molecular.iter().map(|mu| (2.0 / 3.0) * mu).collect();
    solver.start_time(0.1, TimeScheme::BackwardEuler);
    solver
        .enable_sgs(
            FvSgsConfig {
                model: SgsModel::Smagorinsky { cs: 0.2 },
                filter_width_factor: 1.0,
                maximum_eddy_to_molecular: 1e6,
                usage: SgsUse::Exploration,
            },
            &gate,
        )
        .unwrap();
    solver
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * actual.abs().max(expected.abs()).max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn heterogeneous_edge_power_uses_one_total_stress_and_actual_material_density() {
    let domain = VoxelDomain::graded([vec![0.1, 0.3], vec![0.2, 0.5], vec![0.4]]).unwrap();
    let config = SimpleConfig::new([FvBoundary::Symmetry; 6]);
    let density = [1.0, 7.0, 2.0, 11.0];
    let molecular = [0.02, 0.09, 0.04, 0.03];
    let mut solver = sampler(&domain, &config, &density, &molecular);
    // A manufactured reciprocal shear sample; no solve or projection changes
    // this observed edge. Its actual native closure supplies the coefficient.
    for index in 0..solver.vel[0].len() {
        let at = solver.comps[0].coords(index);
        solver.vel[0][index] = 3.0 * domain.center(0, at[1], 0)[1];
    }
    solver.update_sgs(&CancelGate::new_clock_free()).unwrap();
    let cells = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]];
    let eddy = solver
        .sgs_edge_eddy([0, 1, 0], [1, 1, 0], 1, false)
        .unwrap();
    assert!(eddy > 0.0);
    let total: [f64; 4] = core::array::from_fn(|c| molecular[c] + density[c] * eddy);
    let volume: [f64; 4] = core::array::from_fn(|c| domain.volume(c) / 4.0);
    let edge_volume: f64 = volume.iter().sum();
    let harmonic = 1.0
        / (0..4)
            .map(|c| volume[c] / edge_volume / total[c])
            .sum::<f64>();
    let native = solver.edge_mu([0, 1, 0], [1, 1, 0], 1, false);
    near(native, harmonic);
    let power = native * edge_volume * 9.0;
    let partition = edge_power(&solver, cells, [0, 1], eddy, power).unwrap();
    let stress = harmonic * 3.0;
    for c in 0..4 {
        // Independent local traction continuity: gamma_c = tau / mu_total,c.
        let local_gradient = stress / total[c];
        near(
            partition[c][0],
            molecular[c] * local_gradient * local_gradient * volume[c],
        );
        near(
            partition[c][1],
            density[c] * eddy * local_gradient * local_gradient * volume[c],
        );
        assert!(partition[c][0] > 0.0 && partition[c][1] > 0.0);
    }
    near(partition.iter().flatten().sum(), power);
}

#[test]
fn pure_volume_rate_keeps_original_molecular_heat_without_modeled_bulk_loss() {
    let domain = VoxelDomain::new(3, 3, 3, 1.0).unwrap();
    let config = SimpleConfig::new([FvBoundary::Symmetry; 6]);
    let mut solver = sampler(&domain, &config, &[3.0; 27], &[0.7; 27]);
    let rate = 0.125;
    for a in 0..3 {
        for index in 0..solver.vel[a].len() {
            solver.vel[a][index] = rate * solver.comps[a].coords(index)[a] as f64;
        }
    }
    solver.update_sgs(&CancelGate::new_clock_free()).unwrap();
    let center = 1 + 3 * (1 + 3);
    // A positive closure coefficient must still contribute no spherical
    // stress. The original 2 mu D molecular response keeps its full loss.
    solver.eddy_viscosity[center] = 0.3;
    let power = solver
        .sample_material_sgs_power(&CancelGate::new_clock_free())
        .unwrap();
    near(power.heat[center], 6.0 * 0.7 * rate * rate);
    assert_eq!(power.sgs_transfer[center], 0.0);
    near(
        power.quadrature_w,
        power.molecular_heat_w + power.sgs_transfer_w,
    );
}

#[test]
fn resolved_sgs_wall_work_and_heating_stay_entirely_molecular() {
    let domain = VoxelDomain::new(3, 3, 3, 0.5).unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall {
        velocity: [1.0, 0.0, 0.0],
    };
    let config = SimpleConfig::new(faces);
    let mut solver = sampler(&domain, &config, &[2.0; 27], &[0.4; 27]);
    // Resting fluid beside a moving wall has positive modeled gradients, but
    // the native wall closure itself remains molecular and supplies all work.
    let original = solver
        .sample_molecular_power::<false>(&CancelGate::new_clock_free())
        .unwrap();
    solver.update_sgs(&CancelGate::new_clock_free()).unwrap();
    assert!(solver.eddy_viscosity.iter().any(|eddy| *eddy > 0.0));
    let power = solver
        .sample_material_sgs_power(&CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(power.heat, original.heat);
    assert_eq!(power.wall_power_w, original.wall_power_w);
    assert_eq!(
        power.wall_tangential_force_n,
        original.wall_tangential_force_n
    );
    assert_eq!(power.quadrature_w, original.quadrature_w);
    assert_eq!(power.sgs_transfer_w, 0.0);
    assert!(power.sgs_transfer.iter().all(|value| *value == 0.0));
}
