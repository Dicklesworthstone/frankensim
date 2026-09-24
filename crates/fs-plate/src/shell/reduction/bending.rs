//! Physical flexural energy in the original shell coordinates. Drilling is
//! numerical stabilization, not a material modulus or a relaxation channel.
use super::super::{ShellReduction, ShellMesh, PlateError, bad, dot};

impl ShellReduction {
    /// Number of original material facets, in source connectivity order.
    #[must_use]
    pub fn facet_count(&self)->usize {self.facets.len()}

    /// Project the existing DKT bending law, excluding membrane and drilling.
    /// Returns row-major K_b with H_b = q^T K_b q / 2. All modal cross terms and
    /// original physical rotations are retained. `facet_weights` are nonnegative
    /// multipliers of each source facet's material energy, not mode gains.
    ///
    /// This cold query reconstructs D from the SAME retained area-integrated
    /// membrane resultant (A=12D/h^2), at roundoff precision. It does not alter
    /// the equilibrium potential, eigenvectors, mass or any dynamic history.
    /// Constant translation is removed before spatial differencing, exactly as
    /// in the original shell projection; it creates no internal bending loss.
    ///
    /// # Errors
    /// Incomplete/nonfinite/negative weights, exceeded facet-mode-pair budget,
    /// failed original DKT admission or nonfinite derived coefficients.
    pub fn bending_stiffness(&self, facet_weights: &[f64], max_facet_mode_pairs: usize)
        -> Result<Vec<f64>, PlateError>
    {
        let n=self.mode_count();
        if facet_weights.len()!=self.facets.len()
            || facet_weights.iter().any(|x| !x.is_finite() || *x<0.0)
            || n.checked_mul(n).and_then(|x|x.checked_mul(self.facets.len()))
                .is_none_or(|x|x>max_facet_mode_pairs) {
            return Err(bad("shell bending projection weights or work budget are invalid"));
        }
        let mesh=ShellMesh {nodes:self.reference_positions.clone(),tris:self.triangles.clone()};
        let mut matrix=vec![0.0;n*n];
        let mut local=vec![[0.0;9];n];
        for (f,tri) in self.triangles.iter().enumerate() {
            if facet_weights[f]==0.0 {continue;}
            let g=mesh.facet(f)?;let h=self.section_thicknesses[f];
            let inverse=h*h/(12.0*g.area_m2);
            let d=self.facets[f].membrane.map(|a|a*inverse);
            if d.iter().any(|x|!x.is_finite()) {return Err(bad("shell bending section overflow"));}
            let (k,_) = crate::dkt_stiffness(&g.x,&g.y,&d,f)?;
            for (mode,row) in local.iter_mut().enumerate() {
                let reference=self.translations[mode*self.nodes+tri[0]];
                for (a,&node) in tri.iter().enumerate() {
                    let u=self.translations[mode*self.nodes+node];
                    let rotation=self.rotations[mode*self.nodes+node];
                    row[3*a]=dot(g.frame[2],core::array::from_fn(|c|u[c]-reference[c]));
                    row[3*a+1]=-dot(g.frame[1],rotation);
                    row[3*a+2]=dot(g.frame[0],rotation);
                }
            }
            for i in 0..n {for j in 0..=i {
                let mut value=0.0;
                for a in 0..9 {for b in 0..9 {value+=local[i][a]*k[9*a+b]*local[j][b];}}
                matrix[i*n+j]+=facet_weights[f]*value;
                matrix[j*n+i]=matrix[i*n+j];
            }}
        }
        if matrix.iter().any(|x|!x.is_finite()) {return Err(bad("shell bending projection overflow"));}
        Ok(matrix)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ModePair,PlateSection};
    use crate::shell::{assemble_shell_sections,ShellSupport};
    use crate::shell::reduction::ReductionBudget;
    fn fixture()->(ShellMesh,Vec<PlateSection>,crate::shell::ShellModel,Vec<ModePair>,ShellReduction) {
        let mesh=ShellMesh::new(vec![[0.,0.,0.],[0.1,0.,0.01],[0.,0.1,0.02],[0.08,0.09,0.035]],
            vec![[0,1,3],[0,3,2]]).unwrap();
        let sections=vec![PlateSection::isotropic(112.6e9,0.342,0.001,8607.).unwrap(),
            PlateSection::isotropic(90e9,0.31,0.0007,8000.).unwrap()];
        let model=assemble_shell_sections(&mesh,&sections,&[0,1,2],ShellSupport::Clamped).unwrap();
        let n=model.free;let mut k=vec![0.;n*n];let mut m=k.clone();
        for i in 0..n {for j in 0..n {k[i*n+j]=model.k.get(i,j);m[i*n+j]=model.m.get(i,j);}}
        let modes=fs_modal::eigh_gen_dense(&k,&m,n).unwrap();
        let r=ShellReduction::new(&mesh,&sections,&model,&modes,ReductionBudget {
            max_modes:12,max_facet_modes:100,relative_tolerance:1e-6}).unwrap();
        (mesh,sections,model,modes,r)
    }
    #[test]
    fn weighted_flexure_matches_independent_nodal_dkt_work_and_keeps_cross_terms() {
        let (mesh,sections,model,modes,r)=fixture();let n=modes.len();
        let weights=[0.3,1.7];let kb=r.bending_stiffness(&weights,1000).unwrap();
        let q:Vec<_>=(0..n).map(|i|1e-6*(i as f64-2.3)).collect();
        let mut nodal=vec![[0.;6];mesh.nodes.len()];
        for (node,v) in nodal.iter_mut().enumerate() {for c in 0..6 {
            v[c]=model.dof_map[6*node+c].map_or(0.,|free|
                modes.iter().zip(&q).map(|(mode,q)|mode.phi[free]*q).sum());
        }}
        let mut want=0.;
        for (f,tri) in mesh.tris.iter().enumerate() {
            let g=mesh.facet(f).unwrap();let (k,_)=crate::dkt_stiffness(&g.x,&g.y,&sections[f].d,f).unwrap();
            let mut x=[0.;9];
            for (a,&node) in tri.iter().enumerate() {
                x[3*a]=dot(g.frame[2],[nodal[node][0],nodal[node][1],nodal[node][2]]);
                x[3*a+1]=-dot(g.frame[1],[nodal[node][3],nodal[node][4],nodal[node][5]]);
                x[3*a+2]=dot(g.frame[0],[nodal[node][3],nodal[node][4],nodal[node][5]]);
            }
            for i in 0..9 {for j in 0..9 {want+=weights[f]*x[i]*k[9*i+j]*x[j];}}
        }
        let actual:f64=(0..n).map(|i|(0..n).map(|j|q[i]*kb[i*n+j]*q[j]).sum::<f64>()).sum();
        assert!((actual-want).abs()<1e-10*want.abs());assert!(actual>0.);
        assert!((0..n).any(|i|(0..i).any(|j|kb[i*n+j].abs()>1e-6)));
        for i in 0..n {for j in 0..n {assert_eq!(kb[i*n+j],kb[j*n+i]);}}
    }
    #[test]
    fn rigid_translation_and_pure_facet_drilling_have_no_flexural_loss() {
        let (_,_,_,_,mut r)=fixture();
        // Kinematic query fixture, not an eigenpair claim. Replace every shape
        // by a common translation and remove rotation: exact strain-free rows.
        r.translations.fill([3.,-2.,1.]);r.rotations.fill([0.;3]);
        assert!(r.bending_stiffness(&[1.,1.],1000).unwrap().iter().all(|v|*v==0.));
        let normal=ShellMesh{nodes:r.reference_positions.clone(),tris:r.triangles.clone()}.facet(0).unwrap().frame[2];
        r.translations.fill([0.;3]);r.rotations.fill(normal);
        // Isolate the first facet: rotation around ITS normal is only drilling.
        let matrix=r.bending_stiffness(&[1.,0.],1000).unwrap();
        assert!(matrix.iter().all(|v|v.abs()<1e-20));
    }
    #[test]
    fn material_weight_superposition_and_refusal_do_not_change_the_source() {
        let (_,_,_,_,r)=fixture();let n=r.mode_count();let q=vec![1e-7;n];let before=r.potential(&q);
        let a=r.bending_stiffness(&[1.,0.],1000).unwrap();let b=r.bending_stiffness(&[0.,1.],1000).unwrap();
        let c=r.bending_stiffness(&[2.,3.],1000).unwrap();
        for i in 0..c.len() {assert!((c[i]-2.*a[i]-3.*b[i]).abs()<1e-12*c[i].abs().max(1.));}
        assert!(r.bending_stiffness(&[1.],1000).is_err());
        assert!(r.bending_stiffness(&[1.,-1.],1000).is_err());
        assert!(r.bending_stiffness(&[1.,f64::NAN],1000).is_err());
        assert!(r.bending_stiffness(&[1.,1.],2*n*n-1).is_err());
        assert_eq!(r.potential(&q),before);
    }
}
