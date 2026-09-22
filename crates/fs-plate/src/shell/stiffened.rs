//! Bonded eccentric 3-D beam reinforcement of the existing CST/DKT shell.
//!
//! A beam's centroidal endpoints are shell_node + offset. The SAME rigid
//! offset map d_c = d + theta cross offset transforms both stiffness and
//! inertia. This retains axial/bending coupling and the sign of eccentricity;
//! adding EAe^2 to EI alone cannot do that for a deformable shell membrane.
//! Straight Euler-Bernoulli segments, Saint-Venant GJ torsion, perfect bonds,
//! small motions about the supplied geometry, lumped centroidal inertia.
//! No shear deformation, glue slip, stress stiffening or preload equilibrium.

use super::{ShellMesh, ShellModel, ShellSupport, assemble_shell_sections};
use crate::{PlateError, PlateSection};
use fs_sparse::{Coo, Csr};

/// Centroidal section in coherent SI. y/z are the beam's transverse axes.
#[derive(Debug, Clone, Copy)]
pub struct BeamSection {
    pub young_pa: f64,
    pub shear_pa: f64,
    pub density_kg_m3: f64,
    pub area_m2: f64,
    /// Integral z^2 dA; bending in the beam's x/z plane.
    pub iy_m4: f64,
    /// Integral y^2 dA; bending in the beam's x/y plane.
    pub iz_m4: f64,
    /// Saint-Venant torsion constant, NOT the polar area moment Iy+Iz.
    pub torsion_m4: f64,
}
impl BeamSection {
    fn validate(self) -> Result<(), PlateError> {
        if [self.young_pa, self.shear_pa, self.density_kg_m3, self.area_m2,
            self.iy_m4, self.iz_m4, self.torsion_m4].iter()
            .any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(bad("beam section requires finite positive E,G,rho,A,Iy,Iz,J"));
        }
        Ok(())
    }
}

/// Two shell nodes joined by one physical centroidal beam segment.
#[derive(Debug, Clone, Copy)]
pub struct ShellBeam {
    pub nodes: [usize; 2],
    /// Cartesian offsets from each midsurface node to its beam centroid [m].
    pub offsets_m: [[f64; 3]; 2],
    /// Direction defining +z of the cross section. Its projection normal to
    /// the centroidal chord is normalized; collinear directions refuse.
    pub section_z: [f64; 3],
    pub section: BeamSection,
}
fn bad(what: &'static str) -> PlateError { PlateError::BadStiffener { what } }
fn dot(a: [f64;3], b: [f64;3]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
fn cross(a: [f64;3], b: [f64;3]) -> [f64;3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}
fn unit(v: [f64;3]) -> Result<[f64;3], PlateError> {
    let length = dot(v,v).sqrt();
    if !length.is_finite() || length <= 1e-12 { return Err(bad("unresolved beam axis")); }
    Ok(v.map(|x| x/length))
}

/// Assemble reinforced geometry through the existing shell and sparse owners.
/// Each beam is a distinct physical segment; callers must not list it twice.
/// No partially reinforced model escapes on a refusal.
/// # Errors
/// Existing shell refusals, invalid beam geometry/section, or numeric overflow.
pub fn assemble_stiffened_shell(mesh: &ShellMesh, sections: &[PlateSection],
    boundary: &[usize], support: ShellSupport, beams: &[ShellBeam])
    -> Result<ShellModel, PlateError> {
    let mut model = assemble_shell_sections(mesh, sections, boundary, support)?;
    if beams.is_empty() { return Ok(model); }
    let mut k = Coo::new(model.free, model.free);
    let mut m = Coo::new(model.free, model.free);
    copy_sparse(&model.k, &mut k);
    copy_sparse(&model.m, &mut m);
    for beam in beams {
        let (ke, me) = element(mesh, beam)?;
        for i in 0..12 {
            let Some(row) = model.dof_map[6*beam.nodes[i/6]+i%6] else { continue; };
            for j in 0..12 {
                let Some(col) = model.dof_map[6*beam.nodes[j/6]+j%6] else { continue; };
                if ke[12*i+j] != 0.0 { k.push(row,col,ke[12*i+j]); }
                if me[12*i+j] != 0.0 { m.push(row,col,me[12*i+j]); }
            }
        }
    }
    model.k = k.assemble(); model.m = m.assemble();
    if (0..model.free).any(|row| model.k.row(row).1.iter()
        .chain(model.m.row(row).1).any(|x| !x.is_finite())) {
        return Err(bad("reinforced shell accumulation overflow"));
    }
    Ok(model)
}
fn copy_sparse(source: &Csr, target: &mut Coo) {
    for row in 0..source.nrows() {
        let (columns, values) = source.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            target.push(row,column,value);
        }
    }
}

// Element coordinates: [u,v,w,theta_x,theta_y,theta_z] at each shell node.
// Congruence, rather than separately authored force/motion maps, guarantees
// reciprocal work. A beam bending slope is theta_z for v and -theta_y for w.
fn element(mesh: &ShellMesh, beam: &ShellBeam)
    -> Result<([f64;144],[f64;144]),PlateError> {
    beam.section.validate()?;
    if beam.nodes[0] == beam.nodes[1]
        || beam.nodes.iter().any(|&i| i >= mesh.nodes.len())
        || beam.offsets_m.iter().flatten().chain(beam.section_z.iter()).any(|x| !x.is_finite()) {
        return Err(bad("invalid beam nodes, offsets or orientation"));
    }
    let p: [[f64;3];2] = std::array::from_fn(|i|
        std::array::from_fn(|j| mesh.nodes[beam.nodes[i]][j]+beam.offsets_m[i][j]));
    let chord = std::array::from_fn(|i| p[1][i]-p[0][i]);
    let length = dot(chord,chord).sqrt();
    let x = unit(chord)?;
    let hint = unit(beam.section_z)?;
    let projection = dot(hint,x);
    let z = unit(std::array::from_fn(|i| hint[i]-projection*x[i]))?;
    let axes = [x,cross(z,x),z];
    let s = beam.section;
    let mut k = [0.0;144];
    let mut m = [0.0;144];
    for (a,b,coefficient) in [(0,6,s.young_pa*s.area_m2/length),
        (3,9,s.shear_pa*s.torsion_m4/length)] {
        k[12*a+a] += coefficient; k[12*b+b] += coefficient;
        k[12*a+b] -= coefficient; k[12*b+a] -= coefficient;
    }
    let l = length;
    let hermite = [[12.0,6.0*l,-12.0,6.0*l],
        [6.0*l,4.0*l*l,-6.0*l,2.0*l*l],
        [-12.0,-6.0*l,12.0,-6.0*l],
        [6.0*l,2.0*l*l,-6.0*l,4.0*l*l]];
    for (slots,signs,inertia) in [([1,5,7,11],[1.0;4],s.iz_m4),
        ([2,4,8,10],[1.0,-1.0,1.0,-1.0],s.iy_m4)] {
        let coefficient = s.young_pa*inertia/(l*l*l);
        for i in 0..4 { for j in 0..4 {
            k[12*slots[i]+slots[j]] += coefficient*hermite[i][j]*signs[i]*signs[j];
        }}
    }
    // Centroidal nodal masses. Longitudinal extent is represented by the two
    // nodal translational masses, cross-section rotary inertia by Iy and Iz.
    let half = s.density_kg_m3*l/2.0;
    for node in 0..2 { for (dof,inertia) in [s.area_m2,s.area_m2,s.area_m2,
        s.iy_m4+s.iz_m4,s.iy_m4,s.iz_m4].into_iter().enumerate() {
        let i=6*node+dof; m[12*i+i] = half*inertia;
    }}
    let mut transform=[0.0;144];
    for node in 0..2 {
        let base=6*node;
        for (a,axis) in axes.iter().enumerate() {
            for (b,&component) in axis.iter().enumerate() {
                transform[12*(base+a)+base+b]=component;
                transform[12*(base+3+a)+base+3+b]=component;
                let mut basis=[0.0;3]; basis[b]=1.0;
                transform[12*(base+a)+base+3+b]=dot(*axis,cross(basis,beam.offsets_m[node]));
            }
        }
    }
    let k = congruence(&k,&transform);
    let m = congruence(&m,&transform);
    if k.iter().chain(m.iter()).any(|x| !x.is_finite()) {
        return Err(bad("beam stiffness/inertia overflow"));
    }
    Ok((k,m))
}
fn congruence(a: &[f64;144], t: &[f64;144]) -> [f64;144] {
    let mut at=[0.0;144]; let mut out=[0.0;144];
    for i in 0..12 { for j in 0..12 { for k in 0..12 {
        at[12*i+j] += a[12*i+k]*t[12*k+j];
    }}}
    for i in 0..12 { for j in i..12 {
        let value=(0..12).map(|k|t[12*k+i]*at[12*k+j]).sum();
        out[12*i+j]=value; out[12*j+i]=value;
    }}
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    fn mesh() -> ShellMesh {
        ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.],[0.,0.5,0.]],vec![[0,1,2]]).unwrap()
    }
    fn beam() -> ShellBeam {
        ShellBeam { nodes:[0,1], offsets_m:[[0.,0.,0.04];2], section_z:[0.,0.,1.],
            section:BeamSection {young_pa:12e9,shear_pa:0.8e9,density_kg_m3:450.,
                area_m2:0.001,iy_m4:1.3e-7,iz_m4:4.5e-8,torsion_m4:8e-8} }
    }
    fn quadratic(a:&[f64;144],u:&[f64;12])->f64 {
        (0..12).map(|i|u[i]*(0..12).map(|j|a[12*i+j]*u[j]).sum::<f64>()).sum()
    }
    #[test]
    fn offset_beam_preserves_all_six_rigid_motions() {
        let mut mesh=mesh(); mesh.nodes[1]=[0.8,0.3,0.2];
        let b=beam(); let (k,_)=element(&mesh,&b).unwrap();
        let scale=k.iter().fold(0.0_f64,|a,x|a.max(x.abs()));
        for axis in 0..3 { for rotational in [false,true] {
            let mut direction=[0.0;3];direction[axis]=1.0;
            let mut u=[0.0;12];
            for node in 0..2 {
                let displacement=if rotational {cross(direction,mesh.nodes[node])} else {direction};
                u[6*node..6*node+3].copy_from_slice(&displacement);
                if rotational {u[6*node+3..6*node+6].copy_from_slice(&direction);}
            }
            for i in 0..12 {
                let residual=(0..12).map(|j|k[12*i+j]*u[j]).sum::<f64>();
                assert!(residual.abs()<1e-12*scale,"rigid residual {residual}");
            }
        }}
    }
    #[test]
    fn axial_bending_torsion_and_eccentricity_have_analytic_energy() {
        let mesh=mesh(); let b=beam(); let (k,m)=element(&mesh,&b).unwrap();
        let s=b.section;
        let mut axial=[0.0;12];axial[6]=1.0;
        assert!((quadratic(&k,&axial)/(s.young_pa*s.area_m2)-1.0).abs()<1e-13);
        // Constant curvature w=x^2/2, theta_y=-x. Offset produces centroidal
        // axial strain -e and hence precisely EAe^2 in ADDITION to EI.
        let mut bending=[0.0;12];bending[8]=0.5;bending[10]=-1.0;
        let want=s.young_pa*(s.iy_m4+s.area_m2*0.04_f64.powi(2));
        assert!((quadratic(&k,&bending)/want-1.0).abs()<1e-12);
        let mut centered=b;centered.offsets_m=[[0.0;3];2];
        let (kc,_)=element(&mesh,&centered).unwrap();
        let mut torsion=[0.0;12];torsion[9]=1.0;
        assert!((quadratic(&kc,&torsion)/(s.shear_pa*s.torsion_m4)-1.0).abs()<1e-13);
        for axis in 0..3 {
            let mut u=[0.0;12];u[axis]=1.;u[6+axis]=1.;
            assert!((quadratic(&m,&u)/(s.density_kg_m3*s.area_m2)-1.).abs()<1e-13);
        }
        // Changing the sign of the offset must change axial/rotation coupling.
        let mut other=b;other.offsets_m=[[0.,0.,-0.04];2];
        let (negative,_)=element(&mesh,&other).unwrap();
        assert!(k[12*6+10].abs()>1.0);
        assert!((k[12*6+10]+negative[12*6+10]).abs()<1e-12);
    }
    #[test]
    fn reinforced_shell_uses_the_same_supports_and_adds_positive_energy_and_mass() {
        let mesh=mesh(); let b=beam();
        let section=PlateSection::isotropic(12e9,0.3,0.008,450.).unwrap();
        let plain=assemble_shell_sections(&mesh,&[section],&[0],ShellSupport::Clamped).unwrap();
        let reinforced=assemble_stiffened_shell(&mesh,&[section],&[0],ShellSupport::Clamped,&[b]).unwrap();
        assert_eq!(plain.dof_map,reinforced.dof_map);
        let u:Vec<_>=(0..plain.free).map(|i|((i+3) as f64).sin()).collect();
        for (a,b) in [(&plain.k,&reinforced.k),(&plain.m,&reinforced.m)] {
            let mut av=vec![0.;plain.free];let mut bv=av.clone();
            a.spmv(&u,&mut av);b.spmv(&u,&mut bv);
            let extra:f64=u.iter().zip(av.iter().zip(&bv)).map(|(u,(a,b))|u*(b-a)).sum();
            assert!(extra>0.0);
        }
        assert_eq!(reinforced.k,assemble_stiffened_shell(&mesh,&[section],&[0],ShellSupport::Clamped,&[b]).unwrap().k);
    }
    #[test]
    fn invalid_beams_refuse_before_returning_a_partial_model() {
        let mesh=mesh();
        let mut b=beam();b.nodes=[0,3];assert!(element(&mesh,&b).is_err());
        let mut b=beam();b.offsets_m[0][0]=f64::NAN;assert!(element(&mesh,&b).is_err());
        let mut b=beam();b.section_z=[1.,0.,0.];assert!(element(&mesh,&b).is_err());
        let mut b=beam();b.section.iy_m4=-1.;assert!(element(&mesh,&b).is_err());
        let mut b=beam();b.section.area_m2=f64::MAX;assert!(element(&mesh,&b).is_err());
    }
}
