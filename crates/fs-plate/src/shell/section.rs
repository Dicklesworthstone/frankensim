//! Independent effective membrane A and bending D for the original shell owner.
//!
//! A homogeneous solid plate may use A = 12 D / h^2. A characterized textile
//! must not: its measured N/m and N m quantities are separate capabilities.
use crate::{PlateError, PlateSection};
use fs_material::sheet::SheetMatrix;

/// Complete independent elastic shell section in the local facet frame.
/// Membrane and bending are symmetric positive definite; there is no membrane-
/// bending coupling B, prestrain, preferred curvature or hidden material loss.
/// Inertia uses the original plate layout, including uniform-through-thickness
/// rotary inertia. For an effective sheet, density = reference areal mass / h
/// is an inertia representation only, not a constituent density law.
#[derive(Clone, Copy, Debug)]
pub struct ShellSection {
    plate: PlateSection,
    membrane: SheetMatrix,
}
impl ShellSection {
    /// Independent effective membrane A [N/m] and original bending/inertia data.
    /// Both matrices must already be expressed in the SAME local facet axes.
    /// # Errors
    /// Invalid bending rigidity, thickness or inertia representation.
    pub fn new(plate: PlateSection, membrane: SheetMatrix) -> Result<Self, PlateError> {
        plate.validate()?;
        let out = Self { plate, membrane };
        out.validate_section()?;
        Ok(out)
    }
    /// Original bending/inertia section; no membrane-to-bending conversion.
    pub const fn plate(&self) -> &PlateSection { &self.plate }
    /// Independent effective engineering-strain stiffness [N/m].
    pub const fn membrane(&self) -> &SheetMatrix { &self.membrane }
}

// One private lowering interface for linear assembly, geometric nonlinear
// reduction, static preload and drilling. Public legacy entrypoints retain
// their signatures; a caller never has to infer the type of an empty slice.
pub(crate) trait Section {
    fn plate(&self) -> &PlateSection;
    fn membrane(&self) -> Result<[f64;9], PlateError>;
    fn validate_section(&self) -> Result<(), PlateError> {
        self.plate().validate()?;
        let a = self.membrane()?;
        if a.iter().any(|x| !x.is_finite()) || [a[0],a[4],a[8]].iter().any(|x| *x<=0.) {
            return Err(PlateError::BadSection { what: "unrepresentable shell membrane stiffness" });
        }
        Ok(())
    }
}
impl Section for PlateSection {
    fn plate(&self) -> &PlateSection { self }
    fn membrane(&self) -> Result<[f64;9], PlateError> {
        // Preserve the established homogeneous solid-section arithmetic.
        let a = self.d.map(|d| (12.0/(self.thickness*self.thickness))*d);
        if a.iter().any(|x| !x.is_finite()) {
            return Err(PlateError::BadSection { what: "unrepresentable homogeneous shell membrane stiffness" });
        }
        Ok(a)
    }
}
impl Section for ShellSection {
    fn plate(&self) -> &PlateSection { &self.plate }
    fn membrane(&self) -> Result<[f64;9], PlateError> { Ok(*self.membrane.coefficients()) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shell::{ShellMesh,ShellModel,ShellSupport,assemble_shell_materials,assemble_shell_sections};
    fn mesh() -> ShellMesh {
        ShellMesh::new(vec![[0.,0.,0.],[0.2,0.,0.],[0.,0.1,0.]],vec![[0,1,2]]).unwrap()
    }
    fn section(a:f64,d:f64) -> ShellSection {
        ShellSection::new(PlateSection {d:*SheetMatrix::orthotropic(0.04*d,0.02*d,0.005*d,0.01*d).unwrap().coefficients(),
            thickness:0.001,density:300.},SheetMatrix::orthotropic(600.*a,300.*a,50.*a,100.*a).unwrap()).unwrap()
    }
    fn energy(model:&ShellModel,u:&[f64]) -> f64 {
        let mut f=vec![0.;u.len()];model.k.spmv(u,&mut f);
        0.5*u.iter().zip(f).map(|(a,b)|a*b).sum::<f64>()
    }
    fn near(a:f64,b:f64) {assert!((a-b).abs()<1e-11*a.abs().max(b.abs()).max(1e-15),"{a} != {b}");}
    #[test]
    fn measured_membrane_and_bending_change_only_their_own_energy_channels() {
        let mesh=mesh();let build=|a,d|assemble_shell_materials(&mesh,&[section(a,d)],&[],ShellSupport::Free).unwrap();
        let base=build(1.,1.);let tensile=build(3.,1.);let bending=build(1.,5.);
        let mut extension=vec![0.;18];let mut curvature=extension.clone();
        for (i,p) in mesh.nodes.iter().enumerate() {
            extension[6*i]=0.02*p[0];
            curvature[6*i+2]=0.5*p[0]*p[0];curvature[6*i+4]=-p[0];
        }
        assert!(energy(&base,&extension)>0. && energy(&base,&curvature)>0.);
        near(energy(&tensile,&extension),3.*energy(&base,&extension));
        near(energy(&bending,&extension),energy(&base,&extension));
        near(energy(&bending,&curvature),5.*energy(&base,&curvature));
        near(energy(&tensile,&curvature),energy(&base,&curvature));
        for row in 0..base.free {assert_eq!(base.m.row(row),tensile.m.row(row));assert_eq!(base.m.row(row),bending.m.row(row));}
    }
    #[test]
    fn independent_sections_preserve_legacy_equivalence_and_six_rigid_motions() {
        let mesh=mesh();let plate=PlateSection::isotropic(4e6,0.25,0.001,300.).unwrap();
        let a=Section::membrane(&plate).unwrap();
        let section=ShellSection::new(plate,SheetMatrix::orthotropic(a[0],a[4],a[1],a[8]).unwrap()).unwrap();
        let old=assemble_shell_sections(&mesh,&[plate],&[],ShellSupport::Free).unwrap();
        let new=assemble_shell_materials(&mesh,&[section],&[],ShellSupport::Free).unwrap();
        for i in 0..new.free {for j in 0..new.free {near(new.k.get(i,j),old.k.get(i,j));near(new.m.get(i,j),old.m.get(i,j));}}
        let anisotropic=assemble_shell_materials(&mesh,&[self::section(1.,1.)],&[],ShellSupport::Free).unwrap();
        for axis in 0..6 {
            let mut u=vec![0.;18];
            for (i,p) in mesh.nodes.iter().enumerate() {
                if axis<3 {u[6*i+axis]=1.;} else {
                    let mut omega=[0.;3];omega[axis-3]=1.;
                    u[6*i]=omega[1]*p[2]-omega[2]*p[1];
                    u[6*i+1]=omega[2]*p[0]-omega[0]*p[2];
                    u[6*i+2]=omega[0]*p[1]-omega[1]*p[0];
                    u[6*i+3..6*i+6].copy_from_slice(&omega);
                }
            }
            let mut residual=vec![0.;18];anisotropic.k.spmv(&u,&mut residual);
            assert!(residual.iter().all(|v|v.abs()<1e-9),"rigid axis {axis}: {residual:?}");
        }
        assert!(assemble_shell_materials(&mesh,&[],&[],ShellSupport::Free).is_err());
    }
}
