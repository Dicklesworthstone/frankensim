//! Boundary-gradient and material-cap checks independent of a momentum solve.

use super::*;
use crate::conjugate::{FluidProperties, FvBoundary, SimpleConfig, TimeScheme, Voxel, VoxelDomain};
use crate::sgs3::SgsUse;

fn fluid(rho: f64) -> FluidProperties {
    FluidProperties {
        density_kg_m3: rho,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.002,
        kinematic_viscosity_m2_s: 0.01,
    }
}

fn controls(normal: usize) -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[2 * normal] = FvBoundary::Outlet;
    faces[2 * normal + 1] = FvBoundary::Outlet;
    SimpleConfig::new(faces)
}

fn config(model: SgsModel) -> super::super::FvSgsConfig {
    super::super::FvSgsConfig {
        model,
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 1.0,
        usage: SgsUse::Exploration,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-13 * expected.abs().max(1.0),
        "actual={actual}, expected={expected}"
    );
}

#[test]
fn reservoir_tangent_trace_changes_only_its_own_cell_gradient_on_graded_axes() {
    for normal in 0..3 {
        let tangent = (normal + 1) % 3;
        let mut widths = [vec![0.5], vec![0.5], vec![0.5]];
        widths[normal] = vec![0.25, 0.75];
        let domain = VoxelDomain::graded_from_fn(widths, |_| Voxel::Fluid).unwrap();
        let controls = controls(normal);
        let fluid = fluid(1.0);
        let mut solver = Solver::new(&domain, &fluid, &controls);
        solver.vel[tangent].fill(0.2);
        let velocity: Vec<_> = (0..domain.cell_count())
            .map(|cell| solver.cell_velocity(cell))
            .collect();
        // Legacy outlets remain zero normal gradient for this constant field.
        for cell in 0..domain.cell_count() {
            assert_eq!(solver.sgs_gradient(cell, &velocity), [[0.0; 3]; 3]);
        }
        solver.pressure_work_axis = Some(normal);
        for cell in 0..domain.cell_count() {
            let position = domain.coords(cell)[normal];
            let mut expected = [[0.0; 3]; 3];
            // Quadratic through the actual 0.125/0.375 m mouth distances
            // and the 0.5 m neighbouring centre separation, with zero trace.
            expected[tangent][normal] = if position == 0 {
                (0.5 / 0.125) * 0.2 / (0.125 + 0.5)
            } else {
                -(0.5 / 0.375) * 0.2 / (0.5 + 0.375)
            };
            let actual = solver.sgs_gradient(cell, &velocity);
            for row in 0..3 {
                for column in 0..3 {
                    close(actual[row][column], expected[row][column]);
                }
            }
        }
    }
}

fn fill_single_layer(solver: &mut Solver<'_>, normal: usize, tangent: usize) {
    // At the low mouth du_n/dt=.8, du_t/dn=.4. At the high mouth
    // du_n/dt=0, du_t/dn=-.4. The two-cell interpolated diagonal
    // entries cancel, leaving these independently known full tensors.
    for face in 0..solver.vel[normal].len() {
        let at = solver.comps[normal].coords(face);
        solver.vel[normal][face] = if at[normal] == 0 {
            if at[tangent] == 0 { -0.1 } else { 0.1 }
        } else {
            0.0
        };
    }
    for face in 0..solver.vel[tangent].len() {
        solver.vel[tangent][face] = if solver.comps[tangent].coords(face)[tangent] == 1 {
            0.05
        } else {
            0.0
        };
    }
}

fn mouth_cells(tangent: usize) -> ([usize; 3], [usize; 3]) {
    let minus = [0; 3];
    let mut plus = minus;
    plus[tangent] = 1;
    (minus, plus)
}

#[test]
fn each_mouth_keeps_its_native_gradient_for_all_models_and_axis_permutations() {
    let h: f64 = 0.25;
    let a: f64 = 0.8;
    let b: f64 = 0.4;
    let strain2 = 0.5 * (a + b).powi(2);
    let deviator2 = (2.0 / 3.0) * (a * b).powi(2);
    let wale = (0.5 * h).powi(2) * deviator2.powf(1.5) / (strain2.powf(2.5) + deviator2.powf(1.25));
    let vreman = 0.1 * h * h * (a * b).abs() / (a * a + b * b).sqrt();
    let models = [
        (SgsModel::Smagorinsky { cs: 0.2 }, 0.003, 0.001),
        (SgsModel::Wale { cw: 0.5 }, wale, 0.0),
        (SgsModel::Vreman { cv: 0.1 }, vreman, 0.0),
    ];
    let gate = CancelGate::new_clock_free();
    for normal in 0..3 {
        for tangent in 0..3 {
            if tangent == normal {
                continue;
            }
            let mut n = [1; 3];
            n[tangent] = 2;
            let domain = VoxelDomain::new(n[0], n[1], n[2], h).unwrap();
            let controls = controls(normal);
            let fluid = fluid(64.0);
            let (minus, plus) = mouth_cells(tangent);
            let molecular = [0.01; 2];
            let masses: Vec<_> = (0..2)
                .map(|cell| (cell + 1) as f64 * domain.volume(cell))
                .collect();
            for (model, expected_low, expected_high) in models {
                let mut solver = Solver::new(&domain, &fluid, &controls);
                solver.pressure_work_axis = Some(normal);
                solver.start_time(0.01, TimeScheme::BackwardEuler);
                fill_single_layer(&mut solver, normal, tangent);
                solver
                    .set_material_transport(&molecular, None, &gate)
                    .unwrap();
                solver.set_sgs_density(&masses, &gate).unwrap();
                solver.enable_sgs(config(model), &gate).unwrap();
                solver.update_sgs(&gate).unwrap();
                for (high, expected) in [(false, expected_low), (true, expected_high)] {
                    close(
                        solver
                            .pressure_mouth_sgs_eddy(tangent, minus, plus, high)
                            .unwrap(),
                        expected,
                    );
                    // Original source density 64 must not enter either total.
                    let harmonic = 2.0 / (1.0 / (0.01 + expected) + 1.0 / (0.01 + 2.0 * expected));
                    close(
                        solver.pressure_mouth_mu(tangent, minus, plus, high),
                        harmonic,
                    );
                }
            }
        }
    }
}

#[test]
fn every_adjacent_material_caps_the_boundary_edge_without_a_cell_cap_false_positive() {
    let domain = VoxelDomain::new(1, 2, 1, 0.25).unwrap();
    let controls = controls(0);
    let fluid = fluid(0.125);
    let gate = CancelGate::new_clock_free();
    let mut solver = Solver::new(&domain, &fluid, &controls);
    solver.pressure_work_axis = Some(0);
    solver.start_time(0.01, TimeScheme::BackwardEuler);
    fill_single_layer(&mut solver, 0, 1);
    let molecular = [0.01; 2];
    let base_mass = vec![domain.volume(0); 2];
    solver
        .set_material_transport(&molecular, None, &gate)
        .unwrap();
    solver.set_sgs_density(&base_mass, &gate).unwrap();
    solver
        .enable_sgs(config(SgsModel::Smagorinsky { cs: 0.2 }), &gate)
        .unwrap();
    solver.update_sgs(&gate).unwrap();
    let accepted = solver.sgs.as_ref().unwrap().mouth_eddy.clone();
    let eddy = solver.eddy_viscosity.clone();
    for cell in 0..2 {
        let mut mass = base_mass.clone();
        mass[cell] *= 6.0;
        solver.set_sgs_density(&mass, &gate).unwrap();
        // Cell cap max ~=.99499, low-mouth edge cap=1.8, high=0.6.
        for (other, &nu_t) in eddy.iter().enumerate() {
            assert!(solver.sgs_density(other) * nu_t / molecular[other] < 1.0);
        }
        for result in [
            solver.admit_sgs_material_endpoint(&molecular, &gate),
            solver.update_sgs(&gate),
        ] {
            assert!(matches!(result,
                Err(ChtError::InvalidInput { field: "simple.sgs", reason })
                    if reason.starts_with(&format!("cell {cell}:"))
            ));
        }
        solver.set_sgs_density(&base_mass, &gate).unwrap();
        solver.update_sgs(&gate).unwrap();
        assert_eq!(solver.sgs.as_ref().unwrap().mouth_eddy, accepted);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        solver.update_sgs(&cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(solver.sgs.as_ref().unwrap().mouth_eddy, accepted);
    solver.update_sgs(&gate).unwrap();
    assert_eq!(solver.sgs.as_ref().unwrap().mouth_eddy, accepted);
}
