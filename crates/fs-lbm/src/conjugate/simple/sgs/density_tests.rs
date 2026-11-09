//! Actual material density in SGS stress and molecular applicability limits.

use super::*;
use crate::conjugate::{FluidProperties, SimpleConfig, TimeScheme, Voxel, VoxelDomain};

fn fluid(rho: f64) -> FluidProperties {
    FluidProperties {
        density_kg_m3: rho,
        specific_heat_j_kg_k: 2.0,
        conductivity_w_m_k: 0.002,
        kinematic_viscosity_m2_s: 0.01,
    }
}

fn config() -> FvSgsConfig {
    FvSgsConfig {
        model: SgsModel::Smagorinsky { cs: 0.2 },
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 1.0,
        usage: SgsUse::Exploration,
    }
}

#[test]
fn actual_density_uses_each_volume_and_refused_updates_preserve_the_field() {
    let domain = VoxelDomain::graded_from_fn(
        [vec![0.125, 0.375, 0.5], vec![0.25, 0.75], vec![0.5]],
        |p| {
            if p[0] > 0.5 && p[1] > 0.25 {
                Voxel::Solid(0)
            } else {
                Voxel::Fluid
            }
        },
    )
    .unwrap();
    let controls = SimpleConfig::new([FvBoundary::Outlet; 6]);
    let fluid = fluid(17.0);
    let gate = CancelGate::new_clock_free();
    let mut solver = Solver::new(&domain, &fluid, &controls);
    assert!(solver.sgs_density_kg_m3.is_empty());
    for cell in 0..domain.cell_count() {
        assert_eq!(solver.sgs_density(cell), 17.0);
    }
    // Binary-exact dimensions and densities distinguish M/V from a shared
    // cell volume, raw mass, and the unrelated flow reference density.
    let expected: Vec<_> = (0..domain.cell_count())
        .map(|cell| {
            if domain.is_fluid(cell) {
                0.5 * (cell + 1) as f64
            } else {
                0.0
            }
        })
        .collect();
    let mass: Vec<_> = expected
        .iter()
        .enumerate()
        .map(|(cell, rho)| rho * domain.volume(cell))
        .collect();
    solver.set_sgs_density(&mass, &gate).unwrap();
    assert_eq!(solver.sgs_density_kg_m3, expected);
    assert_eq!(solver.rho, 17.0);
    for (cell, &rho) in expected.iter().enumerate() {
        assert_eq!(solver.sgs_density(cell), rho);
    }

    let late_fluid = (0..domain.cell_count())
        .rev()
        .find(|&cell| domain.is_fluid(cell))
        .unwrap();
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        let mut proposed = mass.clone();
        proposed[late_fluid] = invalid; // MAX overflows only after M/V.
        assert!(solver.set_sgs_density(&proposed, &gate).is_err());
        assert_eq!(solver.sgs_density_kg_m3, expected);
    }
    assert!(
        solver
            .set_sgs_density(&mass[..mass.len() - 1], &gate)
            .is_err()
    );
    assert_eq!(solver.sgs_density_kg_m3, expected);
    let solid = (0..domain.cell_count())
        .find(|&cell| !domain.is_fluid(cell))
        .unwrap();
    let mut solid_mass = mass.clone();
    solid_mass[solid] = 1.0;
    assert!(solver.set_sgs_density(&solid_mass, &gate).is_err());
    assert_eq!(solver.sgs_density_kg_m3, expected);

    let doubled: Vec<_> = mass.iter().map(|value| 2.0 * value).collect();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        solver.set_sgs_density(&doubled, &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(solver.sgs_density_kg_m3, expected);
    solver.set_sgs_density(&doubled, &gate).unwrap();
    for (cell, &rho) in expected.iter().enumerate() {
        assert_eq!(solver.sgs_density(cell), 2.0 * rho);
    }
}

#[test]
fn affine_shear_uses_actual_density_for_dynamic_stress_and_the_cell_cap() {
    let domain = VoxelDomain::new(3, 4, 2, 0.25).unwrap();
    let controls = SimpleConfig::new([
        FvBoundary::Outlet,
        FvBoundary::Outlet,
        FvBoundary::wall(),
        FvBoundary::Wall {
            velocity: [1.0, 0.0, 0.0],
        },
        FvBoundary::Symmetry,
        FvBoundary::Symmetry,
    ]);
    let gate = CancelGate::new_clock_free();
    let molecular = vec![0.01; domain.cell_count()];
    // u_x=y has |S|=1. Its kinematic closure is independent of any density:
    // nu_t=(0.2*0.25)^2=0.0025, so rho=2 gives mu_total=0.015 and ratio=0.5.
    let expected_eddy = (0.2_f64 * 0.25).powi(2);
    for reference in [0.125, 64.0] {
        let fluid = fluid(reference);
        let mut solver = Solver::new(&domain, &fluid, &controls);
        for face in 0..solver.vel[0].len() {
            let y = solver.comps[0].coords(face)[1];
            solver.vel[0][face] = domain.face_coord(1, y) + 0.5 * domain.width(1, y);
        }
        solver.start_time(0.01, TimeScheme::BackwardEuler);
        solver
            .set_material_transport(&molecular, None, &gate)
            .unwrap();
        solver.enable_sgs(config(), &gate).unwrap();
        let mass: Vec<_> = (0..domain.cell_count())
            .map(|cell| 2.0 * domain.volume(cell))
            .collect();
        solver.set_sgs_density(&mass, &gate).unwrap();
        solver.update_sgs(&gate).unwrap();
        for cell in 0..domain.cell_count() {
            assert!((solver.eddy_viscosity[cell] - expected_eddy).abs() < 2e-15);
            assert!((solver.mu_eff(domain.coords(cell)) - 0.015).abs() < 2e-15);
            assert!((solver.sgs_total_viscosity(cell, expected_eddy) - 0.015).abs() < 2e-15);
        }
        assert_eq!(solver.rho, reference);
        assert_eq!(solver.molecular_viscosity, molecular);

        // Actual rho=8 gives nu_t/(mu/rho)=2 and must refuse even when the
        // scalar reference would accept; rho=2 above must accept even when
        // the scalar reference would refuse.
        let dense: Vec<_> = mass.iter().map(|value| 4.0 * value).collect();
        solver.set_sgs_density(&dense, &gate).unwrap();
        assert!(matches!(
            solver.update_sgs(&gate),
            Err(ChtError::InvalidInput {
                field: "simple.sgs",
                ..
            })
        ));
        assert_eq!(solver.molecular_viscosity, molecular);
        solver.set_sgs_density(&mass, &gate).unwrap();
        solver.update_sgs(&gate).unwrap();
    }
}

#[test]
fn all_four_native_edge_densities_limit_and_transmit_the_total_viscosity() {
    // Exactly one four-cell xy shear edge: a high-density corner cannot be
    // detected accidentally in a different edge's lower-left position.
    let domain = VoxelDomain::new(2, 2, 1, 0.25).unwrap();
    let controls = SimpleConfig::new([FvBoundary::Outlet; 6]);
    let fluid = fluid(1.0);
    let gate = CancelGate::new_clock_free();
    let mut solver = Solver::new(&domain, &fluid, &controls);
    for face in 0..solver.vel[0].len() {
        solver.vel[0][face] = if solver.comps[0].coords(face)[1] == 0 {
            0.1
        } else {
            -0.1
        };
    }
    solver.start_time(0.01, TimeScheme::BackwardEuler);
    let molecular = vec![0.01; domain.cell_count()];
    solver
        .set_material_transport(&molecular, None, &gate)
        .unwrap();
    solver.enable_sgs(config(), &gate).unwrap();
    solver.update_sgs(&gate).unwrap();

    // Native |du_x/dy|=0.8. The external zero-normal-gradient conditions
    // halve each cell-centred derivative, so nu_cell=0.001, nu_edge=0.002.
    let edge_eddy = (0.2_f64 * 0.25).powi(2) * (0.2 / 0.25);
    let cell_eddy = 0.5 * edge_eddy;
    for &eddy in &solver.eddy_viscosity {
        assert!((eddy - cell_eddy).abs() < 2e-15);
    }
    assert!(
        (solver.sgs_edge_eddy([0, 0, 0], [1, 0, 0], 1, true).unwrap() - edge_eddy).abs() < 2e-15
    );
    let fallback_cells = solver.eddy_viscosity.clone();
    let fallback_edges = solver.sgs.as_ref().unwrap().shear_eddy.clone();
    let fallback_stress = solver.edge_mu([0, 0, 0], [1, 0, 0], 1, true);
    let uniform: Vec<_> = (0..domain.cell_count())
        .map(|cell| domain.volume(cell))
        .collect();
    solver.set_sgs_density(&uniform, &gate).unwrap();
    solver.update_sgs(&gate).unwrap();
    assert_eq!(solver.eddy_viscosity, fallback_cells);
    assert_eq!(solver.sgs.as_ref().unwrap().shear_eddy, fallback_edges);
    assert_eq!(
        solver.edge_mu([0, 0, 0], [1, 0, 0], 1, true).to_bits(),
        fallback_stress.to_bits()
    );

    let density = [1.0, 2.0, 3.0, 4.0];
    let mass: Vec<_> = density
        .iter()
        .enumerate()
        .map(|(cell, rho)| rho * domain.volume(cell))
        .collect();
    solver.set_sgs_density(&mass, &gate).unwrap();
    solver.update_sgs(&gate).unwrap();
    // Four equal areas transmit harmonic TOTAL mu_i=0.01+rho_i*0.002,
    // not a harmonic molecular value with one reference-density addition.
    let expected = 4.0
        / density
            .iter()
            .map(|rho| 1.0 / (0.01 + rho * edge_eddy))
            .sum::<f64>();
    let actual = solver.edge_mu([0, 0, 0], [1, 0, 0], 1, true);
    assert!((actual - expected).abs() < 2e-15);
    assert!((actual - (0.01 + fluid.density_kg_m3 * edge_eddy)).abs() > 0.002);
    assert!((actual - (0.01 + 2.5 * edge_eddy)).abs() > 1e-4);
    assert_eq!(actual, solver.edge_mu([0, 0, 0], [0, 1, 0], 0, true));
    assert_eq!(solver.wall_mu([0, 0, 0], [1, 0, 0]), 0.01);
    for (cell, &rho) in density.iter().enumerate() {
        assert!(
            (solver.sgs_total_viscosity(cell, edge_eddy) - (0.01 + rho * edge_eddy)).abs() < 2e-15
        );
    }

    for cell in 0..4 {
        let mut proposed = mass.clone();
        proposed[cell] = 6.0 * domain.volume(cell);
        solver.set_sgs_density(&proposed, &gate).unwrap();
        // Every cell cap still accepts (at most 0.6); only the sole native
        // edge ratio reaches 1.2. Both endpoint admission and refresh must
        // inspect this corner, independently of the other three densities.
        for index in 0..4 {
            assert!(solver.sgs_density(index) * cell_eddy / 0.01 < 1.0);
        }
        assert!(matches!(
            solver.admit_sgs_material_endpoint(&molecular, &gate),
            Err(ChtError::InvalidInput { field: "simple.sgs", reason })
                if reason.starts_with(&format!("cell {cell}:"))
        ));
        assert!(matches!(
            solver.update_sgs(&gate),
            Err(ChtError::InvalidInput { field: "simple.sgs", reason })
                if reason.starts_with(&format!("cell {cell}:"))
        ));
        assert_eq!(solver.molecular_viscosity, molecular);
        solver.set_sgs_density(&mass, &gate).unwrap();
        solver.update_sgs(&gate).unwrap();
        assert_eq!(solver.edge_mu([0, 0, 0], [1, 0, 0], 1, true), actual);
    }
}
