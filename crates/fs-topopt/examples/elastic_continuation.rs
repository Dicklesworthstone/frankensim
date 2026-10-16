//! Run an independently loaded, gradient-checked 3D SIMP continuation study.
//!
//! From the repository root (use the project's RCH lane where required):
//! `cargo run -p fs-topopt --release --example elastic_continuation -- OUTPUT_DIR 4 250000`
//!
//! OUTPUT_DIR must not exist. Optional trailing arguments are accepted updates
//! PER STAGE (default 4) and cumulative Krylov iterations (default 250000).
//! Exports stage_history.csv, gradient_checks.csv, summary.json, and, only if
//! an equilibrium was accepted, design.vtk. Export performs no additional solves.
//! On a numerical, gradient, or resource stop the accepted prefix is still
//! exported and the process exits nonzero. A complete schedule is not optimality.
//!
//! The fixture is a dimensionless unit cube clamped at x=0, with separate
//! unit-resultant y/z loads at x=1 and weights 0.3/0.7. The background is one
//! fixed Kuhn tetrahedral mesh. This is NOT the raw-SDF/CutFEM marquee and makes
//! no continuum, manufacturing, cross-ISA, or experimental-validation claim.

use std::fs::{self, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::ops::ControlFlow;
use std::path::Path;

use fs_rep_mesh::TetComplex;
use fs_topopt::pipeline::LoadCase;
use fs_topopt::{
    ContinuationTermination, DensityElasticity, DensityFilter, DesignPipeline,
    GradientCheckOptions, MultiLoadContinuationReport, MultiLoadGradientCheck,
    MultiLoadOcOptions, MultiLoadOcReport, SimpParams, SolveBudget, SolveControl,
    controlled_gradient_checked_multi_load_continuation,
};

const VOLUME_CAP: f64 = 0.4;

fn solve_fixture(
    updates: usize, linear_iterations: usize,
) -> (TetComplex, Vec<[f64; 3]>, MultiLoadContinuationReport) {
    let (complex, positions) = fs_feec::kuhn_cube(2);
    let mut elasticity = DensityElasticity::new(
        &complex, &positions, 1.0, 0.3, &|p| p[0] < 1e-12,
    );
    let volumes: Vec<f64> = fs_feec::element_geometry(&complex, &positions)
        .vol_signed.iter().map(|v| v.abs()).collect();
    let mut pipeline = DesignPipeline {
        filter: DensityFilter::new(&complex, &positions, 0.15),
        params: SimpParams::default(),
    };
    let loaded = positions.iter().filter(|p| p[0] > 1.0 - 1e-12).count();
    let nodal_force = -1.0 / f64::from(u32::try_from(loaded).expect("fixture face fits u32"));
    let mut force_y = vec![0.0; elasticity.n()];
    let mut force_z = vec![0.0; elasticity.n()];
    for (vertex, p) in positions.iter().enumerate() {
        if p[0] > 1.0 - 1e-12 {
            force_y[3 * vertex + 1] = nodal_force;
            force_z[3 * vertex + 2] = nodal_force;
        }
    }
    let loads = [LoadCase { force: &force_y, weight: 0.3 },
        LoadCase { force: &force_z, weight: 0.7 }];
    let rho = vec![VOLUME_CAP; elasticity.cells()];
    let schedule = [(1.0, 1.0), (2.0, 2.0), (3.0, 4.0), (3.0, 8.0)]
        .map(|(penal, beta)| SimpParams { penal, beta, ..SimpParams::default() });
    let mut callback = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget {
        total_iterations: linear_iterations, ..SolveBudget::default()
    }, &mut callback);
    let report = controlled_gradient_checked_multi_load_continuation(
        &mut pipeline, &mut elasticity, &loads, &rho, &volumes, &schedule,
        MultiLoadOcOptions { volume_fraction: VOLUME_CAP, max_iterations: updates,
            ..MultiLoadOcOptions::default() },
        GradientCheckOptions::default(), &mut control,
    );
    (complex, positions, report)
}

fn write_history(out: &mut impl Write, report: &MultiLoadContinuationReport) -> io::Result<()> {
    writeln!(out, "stage,penal,beta,eta,incoming_volume,restoration_scale,iteration,compliance,case_y,case_z,volume_fraction,max_change,stage_termination")?;
    for stage in &report.stages {
        for row in &stage.history {
            if row.case_compliances.len() != 2 {
                return Err(io::Error::new(io::ErrorKind::InvalidData, "expected two independent loads"));
            }
            writeln!(out, "{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:?}",
                stage.stage, stage.params.penal, stage.params.beta, stage.params.eta,
                stage.incoming_volume_fraction, stage.restoration_scale,
                row.iteration, row.compliance, row.case_compliances[0], row.case_compliances[1],
                row.volume_fraction, row.max_change, stage.termination)?;
        }
    }
    Ok(())
}

fn write_audit(
    out: &mut impl Write, stage: usize, retained_stage: bool, audit: &MultiLoadGradientCheck,
) -> io::Result<()> {
    for probe in &audit.probes {
        writeln!(out, "{},{},{:.17e},{:.17e},{:.17e},{:.17e},{:?},{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{}",
            stage, retained_stage, audit.params.penal, audit.params.beta,
            audit.options.step, audit.options.relative_tolerance, probe.direction,
            probe.active_densities, probe.compliance_analytic, probe.compliance_difference,
            probe.compliance_relative_error, probe.volume_analytic, probe.volume_difference,
            probe.volume_relative_error, audit.passed())?;
    }
    Ok(())
}

fn write_audits(out: &mut impl Write, report: &MultiLoadContinuationReport) -> io::Result<()> {
    writeln!(out, "stage,retained_stage,penal,beta,step,relative_tolerance,direction,active_densities,compliance_analytic,compliance_difference,compliance_relative_error,volume_analytic,volume_difference,volume_relative_error,passed")?;
    for stage in &report.stages {
        if let Some(audit) = &stage.gradient_check {
            write_audit(out, stage.stage, true, audit)?;
        }
    }
    if let (Some(stage), Some(audit)) = (report.stopped_stage, &report.rejected_gradient_check) {
        write_audit(out, stage, false, audit)?;
    }
    Ok(())
}

fn write_summary(out: &mut impl Write, report: &MultiLoadContinuationReport) -> io::Result<()> {
    // Every string is a fixed literal or an enum name; there is no unescaped
    // caller-supplied path/text in the JSON. Numeric model fields are admitted.
    write!(out, "{{\"schema\":1,\"kind\":\"fixed_mesh_multi_load_simp_continuation\",\"termination\":\"{:?}\",\"stages_solved\":{},\"stopped_stage\":",
        report.termination, report.stages.len())?;
    match report.stopped_stage {
        Some(stage) => write!(out, "{stage}")?,
        None => write!(out, "null")?,
    }
    write!(out, ",\"linear_solves\":{},\"linear_iterations\":{},\"has_accepted_design\":{},\"accepted_model\":",
        report.work.linear_solves, report.work.linear_iterations, report.last.is_some())?;
    if report.last.is_some() {
        write!(out, "{{\"penal\":{:.17e},\"beta\":{:.17e},\"eta\":{:.17e},\"e_min\":{:.17e}}}",
            report.params.penal, report.params.beta, report.params.eta, report.params.e_min)?;
    } else {
        write!(out, "null")?;
    }
    writeln!(out, ",\"volume_cap\":{VOLUME_CAP:.17e},\"optimality_certified\":false,\"continuum_certified\":false,\"raw_sdf_cutfem\":false}}")
}

fn write_design(
    out: &mut impl Write, complex: &TetComplex, positions: &[[f64; 3]], report: &MultiLoadOcReport,
) -> io::Result<()> {
    let dofs = positions.len().checked_mul(3)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "dof count overflow"))?;
    if report.history.is_empty() || report.rho.len() != complex.tets.len()
        || report.projected_rho.len() != complex.tets.len()
        || report.displacements.len() != 2
        || report.displacements.iter().any(|u| u.len() != dofs)
        || !report.rho.iter().chain(&report.projected_rho)
            .chain(report.displacements.iter().flatten()).all(|x| x.is_finite()) {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "no aligned, finite accepted physical state"));
    }
    writeln!(out, "# vtk DataFile Version 3.0\nAccepted fixed-mesh SIMP state; not a binary geometry or certificate\nASCII\nDATASET UNSTRUCTURED_GRID")?;
    writeln!(out, "POINTS {} double", positions.len())?;
    for p in positions { writeln!(out, "{:.17e} {:.17e} {:.17e}", p[0], p[1], p[2])?; }
    writeln!(out, "CELLS {} {}", complex.tets.len(), 5 * complex.tets.len())?;
    for tet in &complex.tets { writeln!(out, "4 {} {} {} {}", tet[0], tet[1], tet[2], tet[3])?; }
    writeln!(out, "CELL_TYPES {}", complex.tets.len())?;
    for _ in &complex.tets { writeln!(out, "10")?; }
    writeln!(out, "CELL_DATA {}", complex.tets.len())?;
    for (name, field) in [("raw_density", &report.rho), ("projected_density", &report.projected_rho)] {
        writeln!(out, "SCALARS {name} double 1\nLOOKUP_TABLE default")?;
        for value in field { writeln!(out, "{value:.17e}")?; }
    }
    writeln!(out, "POINT_DATA {}", positions.len())?;
    for (load, displacement) in report.displacements.iter().enumerate() {
        writeln!(out, "VECTORS displacement_{load} double")?;
        for u in displacement.chunks_exact(3) {
            writeln!(out, "{:.17e} {:.17e} {:.17e}", u[0], u[1], u[2])?;
        }
    }
    Ok(())
}

fn new_output(path: &Path) -> io::Result<BufWriter<std::fs::File>> {
    Ok(BufWriter::new(OpenOptions::new().write(true).create_new(true).open(path)?))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() == 1 && args[0] == "--help" {
        println!("elastic_continuation OUTPUT_DIR [UPDATES_PER_STAGE [TOTAL_KRYLOV_ITERATIONS]]\nThe output directory must not exist. Defaults: 4 updates per stage, 250000 Krylov iterations.\nOutputs histories, numerical gradient audits, a status summary, and any accepted design.");
        return Ok(());
    }
    if args.is_empty() || args.len() > 3 {
        return Err("usage: elastic_continuation OUTPUT_DIR [UPDATES_PER_STAGE [TOTAL_KRYLOV_ITERATIONS]]".into());
    }
    let updates = args.get(1).map_or(Ok(4), |value| value.parse::<usize>())?;
    let budget = args.get(2).map_or(Ok(250_000), |value| value.parse::<usize>())?;
    let directory = Path::new(&args[0]);
    fs::create_dir(directory)?;
    let (complex, positions, report) = solve_fixture(updates, budget);
    let mut history = new_output(&directory.join("stage_history.csv"))?;
    write_history(&mut history, &report)?;
    history.flush()?;
    let mut audits = new_output(&directory.join("gradient_checks.csv"))?;
    write_audits(&mut audits, &report)?;
    audits.flush()?;
    let mut summary = new_output(&directory.join("summary.json"))?;
    write_summary(&mut summary, &report)?;
    summary.flush()?;
    if let Some(last) = &report.last {
        let mut design = new_output(&directory.join("design.vtk"))?;
        write_design(&mut design, &complex, &positions, last)?;
        design.flush()?;
    }
    eprintln!("{:?}: {} solved stages, {} Krylov iterations; accepted design: {}",
        report.termination, report.stages.len(), report.work.linear_iterations, report.last.is_some());
    if report.termination != ContinuationTermination::ScheduleComplete {
        return Err(io::Error::other(format!("study stopped: {:?}; detail: {:?}",
            report.termination, report.evaluation_stop)).into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_real_accepted_state_exports_without_new_physics_work() {
        let (complex, positions, report) = solve_fixture(0, 250_000);
        assert_eq!(report.termination, ContinuationTermination::ScheduleComplete, "{report:?}");
        let work = report.work;
        let mut vtk = Vec::new();
        let last = report.last.as_ref().unwrap();
        write_design(&mut vtk, &complex, &positions, last).unwrap();
        let text = String::from_utf8(vtk).unwrap();
        for field in ["raw_density", "projected_density", "displacement_0", "displacement_1"] {
            assert!(text.contains(field));
        }
        let mut rows = Vec::new();
        write_history(&mut rows, &report).unwrap();
        assert_eq!(String::from_utf8(rows).unwrap().lines().count(), 5);
        let mut audits = Vec::new();
        write_audits(&mut audits, &report).unwrap();
        assert_eq!(String::from_utf8(audits).unwrap().lines().count(), 9);
        assert_eq!(report.work, work);
        let mut corrupted = last.clone();
        corrupted.displacements[0].pop();
        assert!(write_design(&mut Vec::new(), &complex, &positions, &corrupted).is_err());
    }

    #[test]
    fn g4_unfunded_study_reports_no_accepted_model_or_design() {
        let (_, _, report) = solve_fixture(4, 0);
        assert_eq!(report.termination, ContinuationTermination::EvaluationStopped);
        assert!(report.last.is_none() && report.stages.is_empty());
        let mut bytes = Vec::new();
        write_summary(&mut bytes, &report).unwrap();
        let text = String::from_utf8(bytes).unwrap();
        assert!(text.contains("\"has_accepted_design\":false"));
        assert!(text.contains("\"accepted_model\":null"));
        assert!(text.contains("\"optimality_certified\":false"));
    }
}
