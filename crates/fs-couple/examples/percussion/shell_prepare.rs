//! Shared cold shell reduction for single and paired cymbals.
//! Preserves the existing source pencil, vertical inertia, mode order and budgets.
use super::{Error, specimen, ShellSupport, modes_shell, ModePair, SliceOptions,
    ShellReduction, ReductionBudget};
use fs_plate::shell::survey::MeshShell;

pub fn prepare(specimen:&specimen::Specimen,dt_s:f64)->Result<(MeshShell,ShellReduction),Error> {
    build(specimen,dt_s,None)
}

/// Add heave and two small-angle rocking directions about the reference origin.
/// The original shell pencil supplies their full physical angular inertia;
/// mounts, contacts and radiation consume the same resulting modal coordinates.
pub fn prepare_rocking(specimen:&specimen::Specimen,dt_s:f64,maximum_angle_rad:f64)
    ->Result<(MeshShell,ShellReduction),Error>
{
    if !maximum_angle_rad.is_finite() || maximum_angle_rad<=0. || maximum_angle_rad>0.1 {
        return Err("shell rocking needs an explicit angle limit in (0,0.1] rad".into());
    }
    build(specimen,dt_s,Some(maximum_angle_rad))
}

fn build(specimen:&specimen::Specimen,dt_s:f64,rocking:Option<f64>)->Result<(MeshShell,ShellReduction),Error> {
    let shell=specimen.build()?;
    let model=shell.assemble(&[],ShellSupport::Free)?;
    let pi=std::f64::consts::PI;
    let [lower,upper]=specimen.band_hz;
    if !dt_s.is_finite() || dt_s<=0.0 || 2.0*pi*upper*dt_s>=0.9*pi {
        return Err("supplied shell frequency window exceeds the mechanical Nyquist guard".into());
    }
    let report=modes_shell(&model,((2.0*pi*lower).powi(2),(2.0*pi*upper).powi(2)),&SliceOptions::default())?;
    let mut modes=if rocking.is_some() {
        fs_plate::shell::reduction::rigid::mass_normalized_modes(&shell.mesh,&model,&[
            [0.,0.,1.,0.,0.,0.], [0.,0.,0.,1.,0.,0.], [0.,0.,0.,0.,1.,0.],
        ])?
    }else{
        // Preserve the original vertical-only image, including its normalization.
        let mut phi=vec![0.0;model.free];
        for node in 0..shell.mesh.nodes.len() {if let Some(i)=model.dof_map[6*node+2] {phi[i]=1.0/shell.mass_kg.sqrt();}}
        let mut defect=vec![0.0;model.free];model.k.spmv(&phi,&mut defect);
        let residual=defect.iter().enumerate().map(|(i,r)|r*r/model.m.get(i,i)).sum::<f64>().sqrt();
        vec![ModePair{lambda:0.0,phi,residual,interval:(-residual,residual)}]
    };
    modes.extend(report.modes);
    if modes.len()>63 {return Err("splash retains too many modes for this declared reference; narrow the explicit window or increase the host budget".into());}
    let budget=ReductionBudget{max_modes:63,max_facet_modes:65_536,relative_tolerance:1e-5};
    let reduction=match rocking {
        Some(maximum_rotation_rad)=>ShellReduction::new_with_linearized_rigid_motion(
            &shell.mesh,&shell.sections,&model,&modes,budget,
            fs_plate::shell::reduction::RigidMotionLimit{modes:3,maximum_rotation_rad})?,
        None=>ShellReduction::new(&shell.mesh,&shell.sections,&model,&modes,budget)?,
    };
    Ok((shell,reduction))
}
