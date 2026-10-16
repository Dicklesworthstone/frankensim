//! First-order rigid pose separated from the shell's elastic strain.
//! The same physical translational/rotary pencil supplies inertia; all point
//! and skin ports continue to use u = t + omega cross x and theta = omega.
use super::{ShellReduction,ShellMesh,ShellModel,ModePair,PlateError,bad,dot};

/// Explicit validity domain of a linearized rigid-pose image.
#[derive(Debug,Clone,Copy)]
pub struct RigidMotionLimit {
    /// Number of leading, mass-orthonormal geometric rigid modes (1..=6).
    pub modes:usize,
    /// Largest norm of the physical rotation vector [rad], in (0,0.1].
    /// This is a model restriction, not a certified finite-rotation error.
    pub maximum_rotation_rad:f64,
}
#[derive(Debug,Clone)]
pub(super) struct RigidMotion {
    pub(super) rotations:Vec<[f64;3]>,
    maximum_rotation_rad:f64,
}
fn cross(a:[f64;3],b:[f64;3])->[f64;3] {
    [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
}

/// Construct a physical rigid subspace using the supplied lumped shell matrix.
/// Each motion is `[tx,ty,tz,rx,ry,rz]` about the coordinate origin; a shifted
/// pivot is represented by its corresponding translation. Two fixed-order
/// M-Gram--Schmidt passes retain the requested span without inventing diagonal
/// moments or discarding translation/rotation cross inertia. A free shell with
/// positive diagonal nodal mass is required, so the attached M-inverse residual
/// norm uses that same inertia exactly. This is cold preparation, not an
/// eigensolver or time integrator.
/// # Errors
/// Invalid/free-set dimensions, dependent or nonfinite motions, nonpositive
/// physical inertia, or overflowing arithmetic.
pub fn mass_normalized_modes(mesh:&ShellMesh,model:&ShellModel,motions:&[[f64;6]])
    ->Result<Vec<ModePair>,PlateError>
{
    mesh.validate()?;
    if motions.is_empty() || motions.len()>6 || motions.iter().flatten().any(|x|!x.is_finite())
        || model.free!=6*mesh.nodes.len() || model.dof_map.len()!=model.free
        || model.dof_map.iter().enumerate().any(|(i,j)|*j!=Some(i))
        || model.m.nrows()!=model.free || model.m.ncols()!=model.free
        || model.k.nrows()!=model.free || model.k.ncols()!=model.free {
        return Err(bad("rigid shell basis requires a complete free pencil and one to six finite motions"));
    }
    for row in 0..model.free {
        let (columns,values)=model.m.row(row);
        if !model.m.get(row,row).is_finite() || model.m.get(row,row)<=0.
            || columns.iter().zip(values).any(|(&column,value)|!value.is_finite() || column!=row && *value!=0.) {
            return Err(bad("rigid shell residual requires the existing positive diagonal nodal mass"));
        }
    }
    let mut modes:Vec<ModePair>=Vec::with_capacity(motions.len());
    let mut work=vec![0.;model.free];
    for motion in motions {
        let rotation=[motion[3],motion[4],motion[5]];
        let mut phi=vec![0.;model.free];
        for (node,p) in mesh.nodes.iter().enumerate() {
            let u=cross(rotation,*p);
            for c in 0..3 {phi[6*node+c]=motion[c]+u[c];phi[6*node+3+c]=rotation[c];}
        }
        model.m.spmv(&phi,&mut work);
        let original=phi.iter().zip(&work).map(|(a,b)|a*b).sum::<f64>();
        for _ in 0..2 {for prior in &modes {
            model.m.spmv(&phi,&mut work);
            let projection=prior.phi.iter().zip(&work).map(|(a,b)|a*b).sum::<f64>();
            for (v,p) in phi.iter_mut().zip(&prior.phi) {*v-=projection*p;}
        }}
        model.m.spmv(&phi,&mut work);
        let inertia=phi.iter().zip(&work).map(|(a,b)|a*b).sum::<f64>();
        if !original.is_finite() || original<=0. || !inertia.is_finite()
            || inertia<=256.*f64::EPSILON*original {
            return Err(bad("rigid shell motions are dependent or their inertia is not positive finite"));
        }
        for p in &mut phi {*p/=inertia.sqrt();}
        model.k.spmv(&phi,&mut work);
        let residual=work.iter().enumerate().map(|(i,r)|r*r/model.m.get(i,i)).sum::<f64>().sqrt();
        if !residual.is_finite() || phi.iter().any(|x|!x.is_finite()) {
            return Err(bad("rigid shell normalization or residual overflows"));
        }
        modes.push(ModePair{lambda:0.,phi,residual,interval:(-residual,residual)});
    }
    Ok(modes)
}

pub(super) fn admit(mesh:&ShellMesh,model:&ShellModel,modes:&[ModePair],limit:RigidMotionLimit)
    ->Result<RigidMotion,PlateError>
{
    if limit.modes==0 || limit.modes>6 || limit.modes>=modes.len()
        || !limit.maximum_rotation_rad.is_finite() || limit.maximum_rotation_rad<=0.
        || limit.maximum_rotation_rad>0.1 || model.dof_map.iter().any(Option::is_none) {
        return Err(bad("linearized rigid shell needs one to six free rigid modes, elastic modes and an angle cap in (0,0.1] rad"));
    }
    let mut rotations=Vec::with_capacity(limit.modes);
    for mode in &modes[..limit.modes] {
        if mode.lambda!=0. {return Err(bad("linearized shell rigid prefix must have zero eigenvalues"));}
        let at=|node:usize,c:usize|mode.phi[model.dof_map[6*node+c].expect("free pencil admitted")];
        let origin=mesh.nodes[0];
        let translation:[f64;3]=core::array::from_fn(|c|at(0,c));
        let rotation:[f64;3]=core::array::from_fn(|c|at(0,3+c));
        for (node,p) in mesh.nodes.iter().enumerate() {
            let arm=core::array::from_fn(|c|p[c]-origin[c]);let angular=cross(rotation,arm);
            for c in 0..3 {
                let expected=translation[c]+angular[c];
                let scale=translation[c].abs()+rotation.iter().map(|x|x.abs()).sum::<f64>()
                    *arm.iter().map(|x|x.abs()).sum::<f64>();
                if (at(node,c)-expected).abs()>128.*f64::EPSILON*scale.max(1e-300)
                    || (at(node,3+c)-rotation[c]).abs()>128.*f64::EPSILON*rotation[c].abs() {
                    return Err(bad("supplied shell rigid mode is not a geometric infinitesimal rigid motion"));
                }
            }
        }
        rotations.push(rotation);
    }
    Ok(RigidMotion{rotations,maximum_rotation_rad:limit.maximum_rotation_rad})
}

impl ShellReduction {
    /// Validate a complete body displacement and observe rigid rotation [rad].
    /// `None` identifies the original elastic/linear-displacement image. For the
    /// opt-in rigid-pose image, the norm must stay within the supplied restriction.
    /// No allocation, force clipping, potential barrier or state mutation occurs.
    /// The time owner calls this before accepting initial/endpoint states.
    /// # Errors
    /// Wrong dimensions, nonfinite displacement, or excessive rigid rotation.
    pub fn rigid_rotation(&self,q:&[f64])->Result<Option<[f64;3]>,PlateError> {
        if q.len()!=self.mode_count() || q.iter().any(|x|!x.is_finite()) {
            return Err(bad("shell rotation query requires the complete finite body displacement"));
        }
        let Some(rigid)=&self.rigid_motion else{return Ok(None);};
        let mut angle=[0.;3];
        for (rotation,q) in rigid.rotations.iter().zip(q) {for c in 0..3 {angle[c]+=rotation[c]*q;}}
        let norm=dot(angle,angle).sqrt();
        if !norm.is_finite() || norm>rigid.maximum_rotation_rad {
            return Err(bad("shell state exceeds its declared linearized rigid-rotation limit"));
        }
        Ok(Some(angle))
    }
}
