//! Executable real-FEM / continuum-mean comparison for a heated slab.
//!
//! cargo run -p fs-verify --features thermal-conduction --example conduction_mean -- 4 8.0 4.0 2.0 1.0
//! Arguments: cells per axis, source W/m^3, optional kx ky kz in W/(m K).
//! Geometry and coefficients are illustrative exact nominal declarations.
use std::error::Error;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_conduction::verification::{MeanSolveConfig, solve_with_mean_bound};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess, Nonlinearity, ScalarField,
    ThermalBc, ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let (n, source, k) = match args.as_slice() {
        [] => (4, 2.0, [1.0; 3]),
        [n, source] => (n.parse::<usize>()?, source.parse::<f64>()?, [1.0; 3]),
        [n, source, kx, ky, kz] => (
            n.parse::<usize>()?,
            source.parse::<f64>()?,
            [kx.parse::<f64>()?, ky.parse::<f64>()?, kz.parse::<f64>()?],
        ),
        _ => {
            return Err("usage: conduction_mean [cells-per-axis source-W-per-m3 [kx ky kz]]".into());
        }
    };
    if ![2, 4, 8, 16].contains(&n)
        || !source.is_finite()
        || !(0.0..=1_000.0).contains(&source)
        || k.iter()
            .any(|v| !v.is_finite() || !(0.01..=1_000.0).contains(v))
    {
        return Err("admitted example: n in {2,4,8,16}, source in [0,1000] W/m^3, each k in [0.01,1000] W/(m K)".into());
    }
    let (complex, vertices) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, vertices)?;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "left",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::dirichlet(300.0)?,
        )?
        .region(
            "right",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::dirichlet(300.0)?,
        )?
        .adiabatic_remainder()
        .finish()?;
    let material =
        ConductivityModel::constant_tensor([[k[0], 0.0, 0.0], [0.0, k[1], 0.0], [0.0, 0.0, k[2]]])?;
    let forcing = ScalarField::Uniform(source);
    let mut config = MeanSolveConfig::default();
    config.primal.nonlinearity = Nonlinearity::FixedPoint {
        relaxation: 1.0,
        max_backtracks: 4,
    };
    config.dual = config.primal.clone();
    config.dual.initial = InitialGuess::Uniform(0.0);
    config.flux.max_cells = 6 * n * n * n;
    let gate = CancelGate::new();
    let result = ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 11,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        solve_with_mean_bound(
            &cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: None,
                source: &forcing,
            },
            config,
        )
    })?;
    let exact = 300.0 + source / (12.0 * k[0]);
    let bound = &result.bound;
    println!(
        concat!(
            "{{\"functional\":\"whole-domain-volume-mean-temperature\",\"unit\":\"K\",",
            "\"cells\":{},\"source_w_m3\":{:.17e},\"conductivity_diagonal_w_m_k\":{:?},",
            "\"mean_lower\":{:.17e},\"mean_upper\":{:.17e},",
            "\"candidate_lower\":{:.17e},\"candidate_upper\":{:.17e},",
            "\"primal_energy_upper\":{:.17e},\"dual_energy_upper\":{:.17e},",
            "\"analytic_mean\":{:.17e},\"analytic_contained\":{},",
            "\"scope\":\"nominal linear PDE on the declared conforming polyhedral domain; not a point maximum, CAD or physical-uncertainty certificate\"}}"
        ),
        mesh.element_count(),
        source,
        k,
        bound.enclosure.lo,
        bound.enclosure.hi,
        bound.candidate_mean.lo,
        bound.candidate_mean.hi,
        bound.integral.primal.energy_error_upper,
        bound.integral.dual.energy_error_upper,
        exact,
        bound.enclosure.lo <= exact && exact <= bound.enclosure.hi
    );
    if exact < bound.enclosure.lo || exact > bound.enclosure.hi {
        return Err("analytic comparison is outside the returned mean enclosure".into());
    }
    Ok(())
}
