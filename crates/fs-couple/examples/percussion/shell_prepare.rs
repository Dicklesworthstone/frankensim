//! Shared cold shell reduction for single and paired cymbals.
//! Preserves the existing source pencil, vertical inertia, mode order and budgets.
use super::{Error, specimen, ShellSupport, modes_shell, ModePair, SliceOptions,
    ShellReduction, ReductionBudget};
use fs_plate::shell::survey::MeshShell;

// Cold host ceilings, not a physical bandwidth claim. The previous 32-mode
// cap blocked supplied multi-kHz cymbal bases before the nonlinear shell owner
// could use them. Keep enough headroom for sticks, mounts and other mechanical
// participants under MAX_IMPACT_MODES, and still refuse rather than truncate.
const MAX_RETAINED_SHELL_MODES: usize = 192;
const MAX_SHELL_FACET_MODES: usize = 131_072;

fn reduction_budget(mode_count: usize, facet_count: usize) -> Result<ReductionBudget, Error> {
    if mode_count == 0 || mode_count > MAX_RETAINED_SHELL_MODES {
        return Err(format!("splash shell needs 1..={MAX_RETAINED_SHELL_MODES} retained modes; narrow the explicit frequency window or use a larger host").into());
    }
    let work = facet_count.checked_mul(mode_count).ok_or("splash shell facet-mode work overflow")?;
    if work > MAX_SHELL_FACET_MODES {
        return Err(format!("splash shell reduction needs {work} facet-mode records, above the host ceiling {MAX_SHELL_FACET_MODES}; refine the requested band/mesh together").into());
    }
    Ok(ReductionBudget { max_modes: MAX_RETAINED_SHELL_MODES,
        max_facet_modes: MAX_SHELL_FACET_MODES, relative_tolerance: 1e-5 })
}

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
    let budget=reduction_budget(modes.len(),shell.mesh.tris.len())?;
    let reduction=match rocking {
        Some(maximum_rotation_rad)=>ShellReduction::new_with_linearized_rigid_motion(
            &shell.mesh,&shell.sections,&model,&modes,budget,
            fs_plate::shell::reduction::RigidMotionLimit{modes:3,maximum_rotation_rad})?,
        None=>ShellReduction::new(&shell.mesh,&shell.sections,&model,&modes,budget)?,
    };
    Ok((shell,reduction))
}

#[cfg(test)]
mod budget_tests {
    use super::*;

    #[test]
    fn high_band_shell_budget_accepts_more_than_the_old_32_mode_ceiling_without_truncation() {
        let budget = reduction_budget(96, 512).unwrap();
        assert_eq!(budget.max_modes, MAX_RETAINED_SHELL_MODES);
        assert_eq!(budget.max_facet_modes, MAX_SHELL_FACET_MODES);
        assert!(reduction_budget(33, 512).is_ok());
        assert!(reduction_budget(MAX_RETAINED_SHELL_MODES + 1, 1).is_err());
        assert!(reduction_budget(192, 682).is_ok());
        assert!(reduction_budget(192, 683).is_err());
    }
}
