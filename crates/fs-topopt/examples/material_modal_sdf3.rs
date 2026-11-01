//! Run with:
//! cargo run -p fs-topopt --features cutfem-marquee --example material_modal_sdf3
//!
//! A fixed implicit geometry, an authored heterogeneous constitutive field and
//! one prepared consistent inertia map drive minimum-mass vibration design.
//! These are numerical fixture inputs, not calibrated data for a named material.
use std::ops::ControlFlow;

use fs_ascent::projected_al::ProjectedAlOptions;
use fs_blake3::ContentHash;
use fs_cutfem::elastic3::{CutElasticity3, ElasticMaterialField3, ElasticityError3};
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_topopt::modal::MatrixFreeEigenOptions;
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::modal::Sdf3EigenfrequencyOptions3;
use fs_topopt::sdf3::modal::design::{
    Sdf3EigenfrequencyDesign3, Sdf3EigenfrequencyDesignOptions3, Sdf3EigenfrequencyGoal3,
};
use fs_topopt::sdf3::response::limits::design::MinimumMaterialObjective3;
use fs_topopt::{SimpParams, SolveBudget, SolveControl};

struct Top;
impl CutSdf3 for Top {
    fn value(&self, point: [f64; 3]) -> f64 {
        point[2] - 0.53
    }

    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(0.53, 0.53)
    }

    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let derivative = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(derivative, derivative)
    }
}

fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}

fn hertz(eigenvalue: f64) -> f64 {
    eigenvalue.sqrt() / (2.0 * std::f64::consts::PI)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let domain = HexCell::try_new([0.0; 3], [1.6, 0.9, 0.7])?;
    let field = |cell: HexCell| {
        let right = cell.lo()[0] >= 0.8;
        let diagonal = [20.0, 17.0, 29.0, 8.0, 13.0, 11.0];
        let tensor = std::array::from_fn(|i| {
            std::array::from_fn(|j| if i == j { diagonal[i] } else { 0.0 })
        });
        ElasticAcousticState::try_new_mandel(
            tensor,
            if right { 5.0 } else { 2.0 },
            ContentHash([if right { 72 } else { 71 }; 32]),
            0.01,
        )
        .map_err(|_| ElasticityError3::Invalid("invalid authored material"))
    };
    let mut callback = |_| poll();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 { depth: 1, ..Default::default() },
        &mut callback,
    )?;
    let operator = CutElasticity3::build_with_material_field(
        domain,
        [2, 2, 1],
        &Top,
        ElasticMaterialField3::ByCell(&field),
        &|point| point[0] == 0.0,
        Default::default(),
        &mut quadrature,
    )?;
    // Geometry/material quadrature is performed once. The same immutable
    // inertia map is reused by every candidate and finite-difference probe.
    let mass = operator.prepare_material_mass(100_000, poll)?;
    let mut study = CutDensityStudy3::new(
        operator,
        0.12,
        SimpParams { e_min: 0.02, penal: 3.0, beta: 2.0, eta: 0.5 },
    );
    let rho = vec![0.7; study.cells()];
    let frequency = Sdf3EigenfrequencyOptions3 {
        solver: MatrixFreeEigenOptions {
            count: 2,
            oversampling: 60,
            max_iterations: 8,
            relative_tolerance: 1e-9,
            linear_tolerance: 1e-12,
            seed: 17,
        },
        aggregation_beta: 30.0,
    };
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let baseline = study.evaluate_eigenfrequency(&mass, &rho, frequency, &mut control)?;
    let required_eigenvalue = 0.7 * baseline.objective.aggregate;
    if required_eigenvalue <= 0.0 {
        return Err("this fixture requires a positive smooth-min baseline".into());
    }
    let goal = Sdf3EigenfrequencyGoal3::MinimumMaterial {
        objective: MinimumMaterialObjective3::MassKg {
            scale_kg: baseline.physical_mass_kg,
        },
        volume_scale: 1.0,
    };
    let settings = Sdf3EigenfrequencyDesignOptions3 {
        frequency,
        eigenvalue_reference: required_eigenvalue,
        density_floor: 0.3,
        optimizer: ProjectedAlOptions {
            max_spectral_step: 0.03,
            ..Default::default()
        },
        ..Default::default()
    };
    let mut design = Sdf3EigenfrequencyDesign3::new(
        &mut study, &mass, &rho, goal, settings, &mut control,
    )?;
    let report = design.run(12)?;
    let best = design
        .recheck_best()?
        .ok_or("no freshly rechecked feasible design")?;
    println!("Authored material fixture; discrete zero-prestress vibration estimates.");
    println!(
        "mass: {:.8} -> {:.8} kg",
        baseline.physical_mass_kg, best.evaluation.physical_mass_kg
    );
    println!(
        "lowest mode: {:.8} -> {:.8} Hz; required {:.8} Hz",
        hertz(baseline.objective.modal.values[0]),
        hertz(best.evaluation.objective.modal.values[0]),
        hertz(required_eigenvalue),
    );
    println!(
        "projected volume fraction: {:.8}",
        best.evaluation.volume_fraction
    );
    println!("gradient check passed: {}", design.gradient_check().passed());
    println!(
        "optimizer stop: {:?}; accepted history rows: {}",
        report.stop,
        design.history().len()
    );
    println!("consistent inertia terms retained: {}", mass.scalar_terms());
    Ok(())
}
