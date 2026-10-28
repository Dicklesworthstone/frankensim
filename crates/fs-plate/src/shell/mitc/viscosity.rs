//! Kelvin--Voigt resultants on the original shear-deformable elastic space.
//!
//! Membrane [N s/m], bending [N m s] and transverse shear [N s/m] are
//! independent symmetric tensors in each facet frame. No drilling penalty is
//! dissipated. Bubble rotations are eliminated using the ELASTIC transformation;
//! the complete physical C is then projected by congruence on the supplied modes.
//! This is a fixed condensed trial space, not a frequency-dependent elimination
//! of dynamic bubble coordinates.

use super::{element, MitcBudget, NodalDirectors, PlateError, ShearShellSection, ShellMesh};
use fs_math::det;

/// Full local material resultants; zero must be explicitly supplied.
#[derive(Clone, Copy, Debug)]
pub struct ShearShellViscosity {
    membrane: [f64; 9],
    bending: [f64; 9],
    shear: [f64; 3],
}
fn bad(what: &'static str) -> PlateError { PlateError::BadSection { what } }
fn admit_spd3(v: &[f64; 9]) -> Result<(), PlateError> {
    if v.iter().any(|x| !x.is_finite())
        || (0..3).any(|i| (0..3).any(|j| v[3*i+j] != v[3*j+i])) {
        return Err(bad("shell viscosity must be finite and symmetric"));
    }
    if v.iter().all(|x| *x == 0.0) { return Ok(()); }
    if [v[0], v[4], v[8]].iter().any(|x| *x <= 0.0) {
        return Err(bad("shell viscosity must be positive definite or explicitly zero"));
    }
    let s = [det::sqrt(v[0]), det::sqrt(v[4]), det::sqrt(v[8])];
    let a = v[1]/s[0]/s[1]; let b = v[2]/s[0]/s[2]; let c = v[5]/s[1]/s[2];
    let pivot = 1.0-a*a;
    if !pivot.is_finite() || pivot <= 0.0 {
        return Err(bad("shell viscosity has an unresolved or nonpositive pivot"));
    }
    let t = (c-a*b)/det::sqrt(pivot);
    if !t.is_finite() || 1.0-b*b-t*t <= 0.0 {
        return Err(bad("shell viscosity has an unresolved or nonpositive pivot"));
    }
    Ok(())
}
fn admit3(v: &[f64; 9]) -> Result<(), PlateError> {
    if admit_spd3(v).is_ok() { return Ok(()); }
    fs_material::plane_stress_viscosity::factor_passive_resultant(v)
        .map(|_| ())
        .map_err(|_| bad("shell viscosity must be symmetric positive semidefinite"))
}
impl ShearShellViscosity {
    /// Matrices already expressed in the facet strain coordinates, with
    /// engineering in-plane shear in slot three. Shear is [13,13/23,23];
    /// any physical shear-correction factor is included by the section compiler.
    /// Directional-only loss and undamped subspaces are admitted without a floor.
    pub fn new(membrane: [f64; 9], bending: [f64; 9], shear: [f64; 3])
        -> Result<Self, PlateError>
    {
        admit3(&membrane)?; admit3(&bending)?;
        let [a,b,c] = shear;
        if shear.iter().any(|x| !x.is_finite())
            || (!shear.iter().all(|x| *x == 0.0)
                && (a <= 0.0 || c <= 0.0 || (b/det::sqrt(a)/det::sqrt(c)).abs() >= 1.0)) {
            admit3(&[a,b,0.0,b,c,0.0,0.0,0.0,0.0])?;
        }
        Ok(Self { membrane, bending, shear })
    }
    /// The three independent physical resultants in the facet frame.
    pub const fn matrices(&self) -> ([f64; 9], [f64; 9], [[f64; 2]; 2]) {
        let [a,b,c] = self.shear;
        (self.membrane, self.bending, [[a,b],[b,c]])
    }
    /// Explicit source zero; no weak damping is inserted to resolve a gauge.
    pub fn is_zero(&self) -> bool {
        self.membrane.iter().chain(&self.bending).chain(&self.shear).all(|x| *x == 0.0)
    }
}

/// Checked conservative operation count for condensation and modal projection.
/// This is cold computation, not a wall-clock promise.
pub fn projection_terms(nodes: usize, triangles: usize, modes: usize) -> Option<usize> {
    triangles.checked_mul(1_000_000usize.checked_add(
        324usize.checked_mul(modes)?)?.checked_add(
        18usize.checked_mul(modes)?.checked_mul(modes)?)?)?
        .checked_add(6usize.checked_mul(nodes)?.checked_mul(modes)?)
}

/// Project complete physical material dissipation on already expanded global
/// displacement/rotation vectors. For M-orthonormal modes the result is C [1/s].
/// Work-conjugate support frames must be expanded by MitcShellModel::expand,
/// never reconstructed using dof_map alone. No modal diagonalization occurs.
///
/// Geometry/section/budget admission precedes temporary allocations; cancellation
/// is polled at each element, quadrature rule and mode row. This function makes
/// no eigenpair, material-source, or truncation-accuracy claim.
/// # Errors
/// Missing coverage, invalid vectors, nonpositive viscosity, excessive work,
/// or a refused elastic enrichment.
#[allow(clippy::too_many_arguments)]
pub fn project_viscosity(
    mesh: &ShellMesh, directors: &NodalDirectors, sections: &[ShearShellSection],
    viscosity: &[ShearShellViscosity], modes: &[Vec<f64>],
    budget: MitcBudget, max_projection_terms: usize,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Vec<f64>, PlateError> {
    if cancelled() { return Err(PlateError::Cancelled); }
    if mesh.nodes.len() > budget.max_nodes || mesh.tris.len() > budget.max_triangles {
        return Err(bad("shell viscosity exceeds its declared mesh budget"));
    }
    directors.validate(mesh)?;
    if sections.len() != mesh.tris.len() || viscosity.len() != mesh.tris.len()
        || modes.is_empty() || modes.len() > 512
        || modes.iter().any(|p| p.len() != 6*mesh.nodes.len() || p.iter().any(|x| !x.is_finite())) {
        return Err(bad("shell viscosity requires complete sections and bounded physical mode vectors"));
    }
    let n = modes.len();
    if projection_terms(mesh.nodes.len(), mesh.tris.len(), n)
        .is_none_or(|work| work > max_projection_terms) {
        return Err(bad("shell viscosity projection exceeds its explicit work allowance"));
    }
    for section in sections { section.validate()?; }
    let mut out = vec![0.0; n*n];
    for (e, tri) in mesh.tris.iter().enumerate() {
        if cancelled() { return Err(PlateError::Cancelled); }
        if viscosity[e].is_zero() { continue; }
        let g = mesh.facet(e)?;
        let c = element::viscosity_matrix(mesh, tri, &g, tri.map(|i| directors.vectors[i]),
            &sections[e], &viscosity[e], &mut cancelled)?;
        let local: Vec<[f64; 18]> = modes.iter().map(|p|
            core::array::from_fn(|i| p[6*tri[i/6]+i%6])).collect();
        // One element action per mode, followed by its Gram projection. The
        // previous pairwise contraction repeated the 18x18 action n times.
        let mut action = Vec::with_capacity(n);
        for vector in &local {
            if cancelled() { return Err(PlateError::Cancelled); }
            let row: [f64; 18] = core::array::from_fn(|a|
                (0..18).map(|b| c[a][b]*vector[b]).sum());
            if row.iter().any(|v| !v.is_finite()) {
                return Err(bad("shell viscosity element action overflow"));
            }
            action.push(row);
        }
        for i in 0..n {
            if cancelled() { return Err(PlateError::Cancelled); }
            for j in i..n {
                let mut value = 0.0;
                for a in 0..18 { value += local[i][a]*action[j][a]; }
                out[n*i+j] += value;
                if !out[n*i+j].is_finite() { return Err(bad("shell modal viscosity overflow")); }
                out[n*j+i] = out[n*i+j];
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shell::mitc::{element, ShellMass};
    fn section() -> ShearShellSection {
        ShearShellSection::from_orthotropic_solid(
            2.0e8,1.1e8,0.25,4.0e7,3.0e7,1.0e7,0.012,800.0,0.31,5.0/6.0).unwrap()
    }
    #[test]
    fn g0_same_elastic_enrichment_preserves_proportional_material_work() {
        let mesh=ShellMesh::new(vec![[0.0,0.0,0.0],[0.4,0.03,0.05],
            [0.07,0.3,0.04]],vec![[0,1,2]]).unwrap();
        let g=mesh.facet(0).unwrap();
        let normals=[g.frame[2];3];
        let section=section();
        let (a,d,s)=section.matrices();
        let tau=0.0017;
        let law=ShearShellViscosity::new(a.map(|v|tau*v),d.map(|v|tau*v),
            [s[0][0]*tau,s[0][1]*tau,s[1][1]*tau]).unwrap();
        let (k,_)=element::matrices(&mesh,&mesh.tris[0],&g,normals,&section,ShellMass::Consistent,&mut ||false).unwrap();
        let c=element::viscosity_matrix(&mesh,&mesh.tris[0],&g,normals,&section,&law,&mut ||false).unwrap();
        let scale=k.iter().flatten().map(|x|x.abs()).fold(0.0_f64,f64::max)*tau;
        for i in 0..18 { for j in 0..18 {
            assert!((c[i][j]-tau*k[i][j]).abs()<2e-11*scale);
        } }
        let directors=NodalDirectors::new(&mesh,normals.to_vec(),0.9).unwrap();
        let vectors=vec![(0..18).map(|i|((i*7+3) as f64).sin()*0.03).collect::<Vec<_>>(),
            (0..18).map(|i|((i*11+2) as f64).cos()*0.07).collect()];
        let budget=MitcBudget {max_nodes:3,max_triangles:1};
        let projected=project_viscosity(&mesh,&directors,&[section],&[law],&vectors,
            budget,2_000_000,||false).unwrap();
        for i in 0..2 {for j in 0..2 {
            let mut expected=0.0;
            for a in 0..18 {for b in 0..18 {expected+=vectors[i][a]*c[a][b]*vectors[j][b];}}
            assert!((projected[2*i+j]-expected).abs()<1e-12*expected.abs().max(1.0));
        }}
        assert_ne!(projected[1],0.0);
        assert!(projected[0]*projected[3]>=projected[1]*projected[1]);
        assert!(project_viscosity(&mesh,&directors,&[section],&[law],&vectors,budget,1,||false).is_err());
        assert!(project_viscosity(&mesh,&directors,&[section],&[law],&vectors,budget,2_000_000,||true).is_err());
    }
    #[test]
    fn g0_explicit_zero_and_indefinite_viscosity_have_no_hidden_regularization() {
        let zero=ShearShellViscosity::new([0.0;9],[0.0;9],[0.0;3]).unwrap();
        assert!(zero.is_zero());
        let mut v=[0.0;9]; v[0]=1.0; v[4]=1.0; v[8]=1.0; v[1]=2.0; v[3]=2.0;
        assert!(ShearShellViscosity::new(v,[0.0;9],[0.0;3]).is_err());
        assert!(ShearShellViscosity::new([0.0;9],[0.0;9],[1.0,2.0,1.0]).is_err());
    }

    #[test]
    fn directional_membrane_loss_matches_affine_channel_energy_and_keeps_null_modes() {
        let mesh=ShellMesh::new(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],
            [0.0,1.0,0.0]],vec![[0,1,2]]).unwrap();
        let directors=NodalDirectors::new(&mesh,vec![[0.0,0.0,1.0];3],0.9).unwrap();
        let a=[1.0,-2.0,0.0,0.0,0.5,3.0,0.0,-1.0];
        let modes:Vec<Vec<f64>>=a.iter().enumerate().map(|(i,slope)| {
            mesh.nodes.iter().flat_map(|p| [slope*p[0]+i as f64*p[1]+0.25,
                0.7*p[0]-0.3*p[1],0.0,0.0,0.0,0.0]).collect()
        }).collect();
        let law=ShearShellViscosity::new([4.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
            [0.0;9],[0.0;3]).unwrap();
        let projected=project_viscosity(&mesh,&directors,&[section()],&[law],&modes,
            MitcBudget{max_nodes:3,max_triangles:1},2_000_000,||false).unwrap();
        for i in 0..a.len() {for j in 0..a.len() {
            // Independent affine membrane energy: area * A11 * du_i/dx * du_j/dx.
            let expected=0.5*4.0*a[i]*a[j];
            assert!((projected[a.len()*i+j]-expected).abs()<2e-12);
        }}
        assert!(ShearShellViscosity::new([0.0;9],[0.0;9],[4.0,6.0,9.0]).is_ok());
        assert!(ShearShellViscosity::new([0.0;9],[0.0;9],[0.0,1.0,0.0]).is_err());
    }
}
