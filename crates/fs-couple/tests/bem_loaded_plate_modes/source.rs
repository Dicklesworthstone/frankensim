//! Source-owned geometry, material, mass and modes for the exterior-load tests.
use fs_bem::helmholtz::Medium;
use fs_material::gas::{GasSpec, GasState};
use fs_plate::shell::reduction::{ReductionBudget, ShellReduction};
use fs_plate::shell::{ShellMesh, ShellSupport, assemble_shell, modes_shell};
use fs_plate::{PlateSection, SliceOptions};

use super::material_packs::{ALUMINUM_2024_T3, MaterialPack};

pub const THICKNESS_M: f64 = 0.002;

pub struct Plate {
    pub reduction: ShellReduction,
    pub thickness: Vec<f64>,
}

/// A mildly skewed, finite-thickness aluminium plate, clamped on its perimeter.
/// The skew removes reflection symmetries that would make a modal-basis ladder
/// appear converged merely because all newly admitted modes are orthogonal to
/// the radiation load. It is declared geometry, never a fitted acoustic input.
pub fn plate(nx: usize, ny: usize, retained: usize) -> Plate {
    assert!(nx >= 3 && ny >= 3);
    let mut nodes = Vec::new();
    let mut boundary = Vec::new();
    for j in 0..=ny {
        for i in 0..=nx {
            let u = i as f64 / nx as f64;
            let v = j as f64 / ny as f64;
            nodes.push([0.4 * (u + 0.12 * u * v), 0.3 * (v + 0.07 * u * v), 0.0]);
            if i == 0 || i == nx || j == 0 || j == ny {
                boundary.push(j * (nx + 1) + i);
            }
        }
    }
    let mut triangles = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let a = j * (nx + 1) + i;
            triangles.push([a, a + 1, a + nx + 2]);
            triangles.push([a, a + nx + 2, a + nx + 1]);
        }
    }
    let mesh = ShellMesh::new(nodes, triangles).expect("source shell mesh");
    let pack = MaterialPack::seed(ALUMINUM_2024_T3);
    let property = |name, dimensions| {
        pack.scalar_as_measured(name, dimensions)
            .expect("measured aluminium property")
    };
    let section = PlateSection::isotropic(
        property("young_modulus", fs_qty::Pressure::DIMS),
        property("poisson_ratio", fs_qty::Dims::NONE),
        THICKNESS_M,
        property("density", fs_qty::Density::DIMS),
    )
    .expect("source aluminium section");
    let model = assemble_shell(&mesh, &section, &boundary, ShellSupport::Clamped)
        .expect("original shell stiffness and mass");
    let options = SliceOptions {
        max_lanczos: 4 * model.free,
        max_restarts: model.free + 2,
        mass_diagonal_equilibration: true,
        ..SliceOptions::default()
    };
    let upper = (2.0 * core::f64::consts::PI * 4000.0).powi(2);
    let report =
        modes_shell(&model, (1.0, upper), &options).expect("certified in-vacuo shell modes");
    assert!(report.modes.len() >= retained);
    let reduction = ShellReduction::new(
        &mesh,
        &vec![section; mesh.tris.len()],
        &model,
        &report.modes[..retained],
        ReductionBudget {
            max_modes: retained,
            max_facet_modes: mesh.tris.len() * retained,
            relative_tolerance: 1e-7,
        },
    )
    .expect("mass-orthonormal reduction of the same shell pencil");
    Plate {
        reduction,
        thickness: vec![THICKNESS_M; mesh.nodes.len()],
    }
}

/// The density ladder changes the source pressure at fixed temperature, so
/// density scales while the ideal-gas sound speed remains unchanged. There is
/// no artificial load multiplier in the production coupled operator.
pub fn air(density_fraction: f64) -> Medium {
    let state = GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        293.15,
        101_325.0 * density_fraction,
    )
    .expect("source dry-air thermodynamic state");
    Medium::from(&state)
}
