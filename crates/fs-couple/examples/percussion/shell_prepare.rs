//! Shared cold shell reduction for single and paired cymbals.
//! Preserves the existing source pencil, vertical inertia, mode order and budgets.
use super::{Error, specimen, ShellSupport, modes_shell, ModePair, SliceOptions,
    ShellReduction, ReductionBudget};
use fs_plate::shell::survey::MeshShell;

pub fn prepare(specimen:&specimen::Specimen,dt_s:f64)->Result<(MeshShell,ShellReduction),Error> {
    let shell=specimen.build()?;
    let model=shell.assemble(&[],ShellSupport::Free)?;
    let pi=std::f64::consts::PI;
    let [lower,upper]=specimen.band_hz;
    if !dt_s.is_finite() || dt_s<=0.0 || 2.0*pi*upper*dt_s>=0.9*pi {
        return Err("supplied shell frequency window exceeds the mechanical Nyquist guard".into());
    }
    let report=modes_shell(&model,((2.0*pi*lower).powi(2),(2.0*pi*upper).powi(2)),&SliceOptions::default())?;
    // Keep a true vertical free-translation coordinate for felt mounting.
    // Other rigid rotations/translations are omitted in this bounded example;
    // this is NOT a fully rocking 6-DOF cymbal stand model.
    let mut phi=vec![0.0;model.free];
    for node in 0..shell.mesh.nodes.len() {if let Some(i)=model.dof_map[6*node+2] {phi[i]=1.0/shell.mass_kg.sqrt();}}
    let mut defect=vec![0.0;model.free];model.k.spmv(&phi,&mut defect);
    let residual=defect.iter().enumerate().map(|(i,r)|r*r/model.m.get(i,i)).sum::<f64>().sqrt();
    let mut modes=vec![ModePair{lambda:0.0,phi,residual,interval:(-residual,residual)}];
    modes.extend(report.modes);
    if modes.len()>32 {return Err("splash retains too many modes for this declared reference; narrow the explicit window or increase the host budget".into());}
    let reduction=ShellReduction::new(&shell.mesh,&shell.sections,&model,&modes,
        ReductionBudget{max_modes:32,max_facet_modes:20000,relative_tolerance:1e-5})?;
    Ok((shell,reduction))
}
