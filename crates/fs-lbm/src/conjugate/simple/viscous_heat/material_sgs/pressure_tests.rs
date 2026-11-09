use super::*;
use crate::conjugate::simple::{FvBoundary, SimpleConfig, TimeScheme};
use crate::conjugate::{FvSgsConfig, VoxelDomain};
use crate::sgs3::{SgsModel, SgsUse};
use fs_exec::CancelGate;

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * actual.abs().max(expected.abs()).max(1.0),
        "{actual} != {expected}"
    );
}

fn sampler<'a>(
    domain: &'a VoxelDomain,
    config: &'a SimpleConfig,
    axis: usize,
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
    solver.pressure_work_axis = Some(axis);
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

#[test]
fn graded_mouth_power_uses_actual_density_and_total_traction_continuity_on_each_axis() {
    let density = [1.5, 9.0];
    let molecular = [0.03, 0.11];
    let gate = CancelGate::new_clock_free();
    for a in 0..3 {
        let d = (a + 1) % 3;
        let mut widths = [vec![0.3], vec![0.3], vec![0.3]];
        widths[a] = vec![0.4];
        widths[d] = vec![0.2, 0.6];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mut faces = [FvBoundary::Symmetry; 6];
        faces[2 * a] = FvBoundary::Outlet;
        faces[2 * a + 1] = FvBoundary::Outlet;
        let config = SimpleConfig::new(faces);
        let mut solver = sampler(&domain, &config, a, &density, &molecular);
        let minus = [0, 0, 0];
        let mut plus = minus;
        plus[d] = 1;
        let ld = 0.4;
        let la = 0.2;
        let tangent_face = solver.cell_face(d, minus, true);
        solver.vel[d][tangent_face] = 0.3;
        for (high, slope) in [(false, 0.7), (true, 1.1)] {
            let face = solver.cell_face(a, plus, high);
            solver.vel[a][face] = slope * ld;
        }
        solver.update_sgs(&gate).unwrap();
        let mut expected_heat = [0.0; 2];
        let mut expected_sgs = [0.0; 2];
        let mut expected_power = 0.0;
        for (high, slope, sign) in [(false, 0.7, 1.0), (true, 1.1, -1.0)] {
            let eddy = solver
                .pressure_mouth_sgs_eddy(d, minus, plus, high)
                .unwrap();
            assert!(eddy > 0.0);
            let total: [f64; 2] = core::array::from_fn(|c| molecular[c] + density[c] * eddy);
            let volume: [f64; 2] = core::array::from_fn(|c| domain.volume(c) / 4.0);
            let edge_volume = volume[0] + volume[1];
            let harmonic = edge_volume / (volume[0] / total[0] + volume[1] / total[1]);
            near(solver.pressure_mouth_mu(d, minus, plus, high), harmonic);
            let gamma = slope + sign * 0.3 / la;
            let stress = harmonic * gamma;
            let power = harmonic * edge_volume * gamma * gamma;
            let partition = mouth_power(&solver, [minus, plus], d, eddy, power).unwrap();
            for c in 0..2 {
                // Independent material traction continuity: each subvolume
                // has gamma_c = tau / (mu_c + rho_c nu_t).
                let local_gradient = stress / total[c];
                let heat = molecular[c] * local_gradient * local_gradient * volume[c];
                let modeled = density[c] * eddy * local_gradient * local_gradient * volume[c];
                near(partition[c][0], heat);
                near(partition[c][1], modeled);
                expected_heat[c] += heat;
                expected_sgs[c] += modeled;
            }
            expected_power += power;
        }
        let mut heat = [0.0; 2];
        let mut modeled = [0.0; 2];
        let mut quadrature = 0.0;
        solver
            .pressure_mouth_power::<true>(&mut heat, &mut modeled, &mut quadrature, &gate)
            .unwrap();
        for c in 0..2 {
            near(heat[c], expected_heat[c]);
            near(modeled[c], expected_sgs[c]);
            assert!(heat[c] > 0.0 && modeled[c] > 0.0);
        }
        near(quadrature, expected_power);
        near(heat.iter().chain(&modeled).sum(), quadrature);

        // The complete physical sampler adds original cell normal stresses
        // exactly once. There are no four-cell interior edges or wall losses.
        let full = solver.sample_material_sgs_power(&gate).unwrap();
        for c in 0..2 {
            let at = domain.coords(c);
            let normal: [f64; 3] = core::array::from_fn(|axis| {
                (solver.vel[axis][solver.cell_face(axis, at, true)]
                    - solver.vel[axis][solver.cell_face(axis, at, false)])
                    / domain.width(axis, at[axis])
            });
            let divergence: f64 = normal.iter().sum();
            let molecular_normal =
                2.0 * molecular[c] * domain.volume(c) * normal.iter().map(|g| g * g).sum::<f64>();
            let modeled_normal = 2.0
                * density[c]
                * solver.eddy_viscosity[c]
                * domain.volume(c)
                * normal
                    .iter()
                    .map(|g| (g - divergence / 3.0).powi(2))
                    .sum::<f64>();
            near(full.heat[c], expected_heat[c] + molecular_normal);
            near(full.sgs_transfer[c], expected_sgs[c] + modeled_normal);
        }
        near(
            full.quadrature_w,
            full.molecular_heat_w + full.sgs_transfer_w,
        );
        assert_eq!(full.wall_power_w, [0.0; 6]);
    }
}

#[test]
fn mouth_normal_and_tangent_rates_cancel_with_the_correct_outward_side_sign() {
    let domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let config = SimpleConfig::new(faces);
    let mut solver = sampler(&domain, &config, 0, &[2.0; 2], &[0.4; 2]);
    let gate = CancelGate::new_clock_free();
    let (minus, plus) = ([0, 0, 0], [0, 1, 0]);
    let tangent_face = solver.cell_face(1, minus, true);
    solver.vel[1][tangent_face] = 0.25;
    let low_plus = solver.cell_face(0, plus, false);
    solver.vel[0][low_plus] = -0.5;
    solver.update_sgs(&gate).unwrap();
    let mut heat = [0.0; 2];
    let mut modeled = [0.0; 2];
    let mut total = 0.0;
    solver
        .pressure_mouth_power::<true>(&mut heat, &mut modeled, &mut total, &gate)
        .unwrap();
    // At the low mouth du_n/dt + du_t/dn = -0.5 + 0.25/0.5 = 0.
    // At the high mouth it is 0 - 0.25/0.5 = -0.5. One normal cell
    // shares material IDs between both mouths, but their closures are distinct.
    let high_eddy = solver
        .pressure_mouth_sgs_eddy(1, minus, plus, true)
        .unwrap();
    let high_mu = 0.4 + 2.0 * high_eddy;
    let high_power = high_mu * 0.5 * 0.5_f64.powi(2);
    let expected = mouth_power(&solver, [minus, plus], 1, high_eddy, high_power).unwrap();
    near(total, high_power);
    for c in 0..2 {
        near(heat[c], expected[c][0]);
        near(modeled[c], expected[c][1]);
    }
    assert!(modeled.iter().sum::<f64>() > 0.0);
}

#[test]
fn mouth_stress_and_power_use_the_cached_edge_not_cell_averaged_eddy_viscosity() {
    let domain = VoxelDomain::new(1, 2, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let config = SimpleConfig::new(faces);
    let mut solver = sampler(&domain, &config, 0, &[1.0, 8.0], &[0.2, 0.7]);
    let gate = CancelGate::new_clock_free();
    let tangent_face = solver.cell_face(1, [0, 0, 0], true);
    solver.vel[1][tangent_face] = 0.5;
    solver.update_sgs(&gate).unwrap();
    let sample = |solver: &Solver<'_>| {
        let mut heat = [0.0; 2];
        let mut modeled = [0.0; 2];
        let mut total = 0.0;
        solver
            .pressure_mouth_power::<true>(&mut heat, &mut modeled, &mut total, &gate)
            .unwrap();
        (heat, modeled, total)
    };
    let original = sample(&solver);
    assert!(original.1.iter().sum::<f64>() > 0.0);
    // The normal-stress samples may differ from the native mouth shear sample.
    // Perturb only those cached cell coefficients to expose accidental reuse.
    solver.eddy_viscosity[0] = 100.0;
    solver.eddy_viscosity[1] = 0.0;
    assert_eq!(sample(&solver), original);
}
