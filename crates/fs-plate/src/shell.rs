//! Flat-facet CST/DKT shells. Nodal coordinates are three displacements and
//! three PHYSICAL axial rotations, not plate slopes. In a local tangent frame
//! the DKT slopes are (w_x,w_y)=(-theta_y,theta_x). Drilling stabilization
//! penalizes rotation relative to membrane spin, preserving rigid motions.
//! Lumped translational/rotary mass and the existing fs-modal solver are used.
//! Moderate rotations, facet-local material axes; no plastic forming model.

use crate::{PlateError, PlateSection, dkt_stiffness};
use fs_modal::{SliceOptions, SliceReport, slice_window};
use fs_sparse::{Coo, Csr};

/// Radially sampled shells, thickness fields and explicit surface indentations.
pub mod profile;
/// Nonlinear membrane energy projected from the same shell geometry.
pub mod reduction;
/// Circular prestressed films through the existing DKT/membrane assembly.
pub mod head;

/// Relative membrane-spin stabilization; numerical, not a material property.
pub const DRILLING_ALPHA: f64 = 1e-3;

/// Triangular midsurface, in metres. Each node has (u,v,w,theta_x,theta_y,theta_z).
#[derive(Debug, Clone, PartialEq)]
pub struct ShellMesh {
    /// Cartesian midsurface positions.
    pub nodes: Vec<[f64; 3]>,
    /// Oriented triangle indices.
    pub tris: Vec<[usize; 3]>,
}
impl ShellMesh {
    /// Admit finite coordinates and nondegenerate, in-range triangles.
    /// # Errors
    /// Empty, nonfinite, invalid-index or degenerate mesh.
    pub fn new(nodes: Vec<[f64; 3]>, tris: Vec<[usize; 3]>) -> Result<Self, PlateError> {
        let mesh = Self { nodes, tris };
        mesh.validate()?;
        Ok(mesh)
    }
    /// Number of positions.
    #[must_use]
    pub fn node_count(&self) -> usize { self.nodes.len() }
    /// Number of facets.
    #[must_use]
    pub fn element_count(&self) -> usize { self.tris.len() }
    fn validate(&self) -> Result<(), PlateError> {
        if self.nodes.len() < 3 || self.tris.is_empty()
            || self.nodes.iter().flatten().any(|x| !x.is_finite()) {
            return Err(PlateError::BadSection { what: "shell requires finite nonempty geometry" });
        }
        for e in 0..self.tris.len() { self.facet(e)?; }
        Ok(())
    }
    /// Local orthonormal frame and constant P1 derivatives for one facet.
    /// This same geometry is usable by modal nonlinear membrane reductions.
    /// # Errors
    /// Out-of-range element/node, nonfinite or degenerate geometry.
    pub fn facet(&self, element: usize) -> Result<FacetGeometry, PlateError> {
        let fail = || PlateError::DegenerateElement { element, twice_area: 0.0 };
        let tri = self.tris.get(element).ok_or_else(fail)?;
        for &n in tri {
            if n >= self.nodes.len() { return Err(PlateError::BadBoundary { node: n, node_count: self.nodes.len() }); }
        }
        let a = sub(self.nodes[tri[1]], self.nodes[tri[0]]);
        let b = sub(self.nodes[tri[2]], self.nodes[tri[0]]);
        let length = norm(a);
        let normal = cross(a, b);
        let twice_area = norm(normal);
        if !length.is_finite() || length < 1e-12 || !twice_area.is_finite() || twice_area < 2e-14 { return Err(fail()); }
        let ex = a.map(|v| v / length);
        let ez = normal.map(|v| v / twice_area);
        let ey = cross(ez, ex);
        let x = [0.0, length, dot(b, ex)];
        let y = [0.0, 0.0, dot(b, ey)];
        let gradient = [
            [(y[1]-y[2])/twice_area, (x[2]-x[1])/twice_area],
            [(y[2]-y[0])/twice_area, (x[0]-x[2])/twice_area],
            [(y[0]-y[1])/twice_area, (x[1]-x[0])/twice_area],
        ];
        Ok(FacetGeometry { frame: [ex, ey, ez], x, y, area_m2: 0.5*twice_area, gradient })
    }
}
fn sub(a: [f64;3], b: [f64;3]) -> [f64;3] { [a[0]-b[0],a[1]-b[1],a[2]-b[2]] }
fn dot(a: [f64;3], b: [f64;3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
fn cross(a: [f64;3], b: [f64;3]) -> [f64;3] { [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]] }
fn norm(a: [f64;3]) -> f64 { dot(a,a).sqrt() }

/// Element geometry shared by stiffness and physical modal strain projections.
#[derive(Debug, Clone, Copy)]
pub struct FacetGeometry {
    /// Rows ex, ey, ez mapping Cartesian vectors to local tangent coordinates.
    pub frame: [[f64;3];3],
    /// Local nodal x coordinates [m].
    pub x: [f64;3],
    /// Local nodal y coordinates [m].
    pub y: [f64;3],
    /// Midsurface facet area [m^2].
    pub area_m2: f64,
    /// P1 shape derivatives, [node][x/y], in [1/m].
    pub gradient: [[f64;2];3],
}

/// Strong essential boundary conditions. Compliant supports belong at ports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShellSupport {
    /// No constrained degrees of freedom (six rigid modes remain).
    Free,
    /// All translations and physical rotations constrained.
    Clamped,
    /// Translations constrained, physical rotations free.
    Pinned,
}
/// Reduced stiffness/mass pencil, with explicit full-to-free coordinate map.
#[derive(Debug, Clone)]
pub struct ShellModel {
    /// Membrane, bending and relative-spin stabilization stiffness.
    pub k: Csr,
    /// Positive lumped translational and rotary inertia.
    pub m: Csr,
    /// Full DOF (6*node+component) to free coordinate.
    pub dof_map: Vec<Option<usize>>,
    /// Free coordinate count.
    pub free: usize,
}

/// Assemble a homogeneous shell using the same element path as thickness fields.
/// # Errors
/// Invalid geometry, section, boundary, mass or unrepresentable assembly.
pub fn assemble_shell(mesh: &ShellMesh, section: &PlateSection,
    boundary_nodes: &[usize], support: ShellSupport) -> Result<ShellModel, PlateError> {
    assemble(mesh, core::slice::from_ref(section), true, boundary_nodes, support)
}
/// Assemble one explicit material/thickness section per facet. D is expressed
/// in that facet's ex/ey axes; rotate anisotropic data into those axes beforehand.
/// There is no thickness averaging across facets and no inferred density.
/// # Errors
/// Section count, physical admission, mesh, boundary or finite-set refusal.
pub fn assemble_shell_sections(mesh: &ShellMesh, sections: &[PlateSection],
    boundary_nodes: &[usize], support: ShellSupport) -> Result<ShellModel, PlateError> {
    if sections.len() != mesh.tris.len() {
        return Err(PlateError::SectionCount { expected: mesh.tris.len(), actual: sections.len() });
    }
    assemble(mesh, sections, false, boundary_nodes, support)
}
fn assemble(mesh: &ShellMesh, sections: &[PlateSection], uniform: bool,
    boundary_nodes: &[usize], support: ShellSupport) -> Result<ShellModel, PlateError> {
    mesh.validate()?;
    for section in sections { section.validate()?; }
    let ndof = mesh.nodes.len().checked_mul(6).ok_or(PlateError::BadSection { what: "shell DOF overflow" })?;
    let mut constrained = vec![false; ndof];
    for &node in boundary_nodes {
        if node >= mesh.nodes.len() { return Err(PlateError::BadBoundary { node, node_count: mesh.nodes.len() }); }
        let count = match support { ShellSupport::Free => 0, ShellSupport::Pinned => 3, ShellSupport::Clamped => 6 };
        for c in 0..count { constrained[6*node+c] = true; }
    }
    let mut free = 0;
    let dof_map: Vec<_> = constrained.iter().map(|fixed| if *fixed { None } else { let i=free; free+=1; Some(i) }).collect();
    if free == 0 { return Err(PlateError::BadSection { what: "shell has no free coordinates" }); }
    let mut k = Coo::new(free,free);
    let mut mass = vec![0.0;free];
    for (e,tri) in mesh.tris.iter().enumerate() {
        let section = &sections[if uniform { 0 } else { e }];
        let g = mesh.facet(e)?;
        let local = local_stiffness(&g,section,e)?;
        // Transform translations and physical axial rotations with the same R.
        for i in 0..18 {
            let Some(ri) = dof_map[6*tri[i/6]+i%6] else { continue; };
            for j in 0..18 {
                let Some(cj) = dof_map[6*tri[j/6]+j%6] else { continue; };
                let bi=6*(i/6)+3*((i%6)/3);
                let bj=6*(j/6)+3*((j%6)/3);
                let mut value=0.0;
                for p in 0..3 { for q in 0..3 {
                    value += g.frame[p][i%3]*local[(bi+p)*18+bj+q]*g.frame[q][j%3];
                }}
                if !value.is_finite() { return Err(PlateError::BadSection { what: "shell stiffness overflow" }); }
                k.push(ri,cj,value);
            }
        }
        let m=section.density*section.thickness*g.area_m2/3.0;
        let jr=m*section.thickness*section.thickness/12.0;
        for &node in tri { for c in 0..6 {
            if let Some(i)=dof_map[6*node+c] { mass[i]+=if c<3 {m} else {jr}; }
        }}
    }
    let mut m = Coo::new(free,free);
    for (i,&value) in mass.iter().enumerate() {
        // An isolated coordinate is a mesh error, not permission to invent mass.
        if !value.is_finite() || value<=0.0 { return Err(PlateError::BadSection { what: "shell contains unrepresented mass or isolated free nodes" }); }
        m.push(i,i,value);
    }
    Ok(ShellModel { k:k.assemble(),m:m.assemble(),dof_map,free })
}
fn local_stiffness(g:&FacetGeometry, s:&PlateSection, element:usize) -> Result<[f64;324],PlateError> {
    let mut k=[0.0;324];
    // Membrane resultant modulus A = 12 D / h^2, including D16/D26.
    let a=s.d.map(|d| (12.0/(s.thickness*s.thickness))*d);
    let mut b=[[0.0;18];3];
    for i in 0..3 {
        let [dx,dy]=g.gradient[i];
        b[0][6*i]=dx; b[1][6*i+1]=dy;
        b[2][6*i]=dy; b[2][6*i+1]=dx;
    }
    for i in 0..18 { for j in 0..18 { for p in 0..3 { for q in 0..3 {
        k[i*18+j]+=g.area_m2*b[p][i]*a[3*p+q]*b[q][j];
    }}}}
    let bending = local_bending_stiffness(g, s, element)?;
    for (value, extra) in k.iter_mut().zip(bending) { *value += extra; }
    Ok(k)
}
// Project this positive remainder separately from membrane energy. Subtracting
// two nearly equal dense stiffnesses would destroy positive small eigenvalues.
pub(super) fn local_bending_stiffness(g:&FacetGeometry, s:&PlateSection, element:usize)
    -> Result<[f64;324],PlateError> {
    let mut k=[0.0;324];
    let a=s.d.map(|d| (12.0/(s.thickness*s.thickness))*d);
    let (kb,_) = dkt_stiffness(&g.x,&g.y,&s.d,element)?;
    let slots=[2,4,3,8,10,9,14,16,15];
    let signs=[1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,-1.0,1.0];
    for i in 0..9 { for j in 0..9 { k[slots[i]*18+slots[j]]+=signs[i]*kb[9*i+j]*signs[j]; }}
    // Positive penalty on theta_z - (v_x-u_y)/2. Its lumped three-point
    // integration also controls nonconstant drilling modes without grounding
    // the mean rigid rotation. The old diagonal theta_z spring did not.
    for node in 0..3 {
        let mut spin=[0.0;18];
        spin[6*node+5]=1.0;
        for i in 0..3 { spin[6*i]=0.5*g.gradient[i][1]; spin[6*i+1]=-0.5*g.gradient[i][0]; }
        let scale=DRILLING_ALPHA*a[8]*g.area_m2/3.0;
        for i in 0..18 { for j in 0..18 { k[18*i+j]+=scale*spin[i]*spin[j]; }}
    }
    Ok(k)
}

/// Certified generalized modes in the angular-frequency-squared window.
/// # Errors
/// The existing fs-modal eigensolver's refusals.
pub fn modes_shell(model:&ShellModel,window:(f64,f64),opts:&SliceOptions)->Result<SliceReport,PlateError> {
    Ok(slice_window(&model.k,&model.m,window,opts)?)
}

/// Existing cylindrical geometry helper. For checked production profile input,
/// use [`profile::revolve`] instead.
#[must_use]
pub fn generate_cylinder_shell(r:f64,h:f64,n_theta:usize,n_z:usize)->ShellMesh {
    let mut nodes=Vec::with_capacity((n_theta+1)*(n_z+1));
    for j in 0..=n_z { let z=j as f64/n_z as f64*h; for i in 0..n_theta {
        let theta=i as f64/n_theta as f64*2.0*core::f64::consts::PI;
        nodes.push([r*theta.cos(),r*theta.sin(),z]);
    }}
    let mut tris=Vec::with_capacity(2*n_theta*n_z);
    for j in 0..n_z { for i in 0..n_theta {
        let ni=(i+1)%n_theta;
        let (a,b,c,d)=(j*n_theta+i,j*n_theta+ni,(j+1)*n_theta+i,(j+1)*n_theta+ni);
        tris.push([a,b,d]);tris.push([a,d,c]);
    }}
    ShellMesh {nodes,tris}
}
/// Revolve the existing crown-to-lip bell profile, without admission changes.
#[must_use]
pub fn generate_bell_shell(profile:&[(f64,f64)],n_theta:usize)->ShellMesh {
    let mut nodes=Vec::with_capacity(profile.len()*n_theta);
    for &(r,z) in profile { for i in 0..n_theta {
        let theta=i as f64/n_theta as f64*2.0*core::f64::consts::PI;
        nodes.push([r*theta.cos(),r*theta.sin(),z]);
    }}
    let mut tris=Vec::with_capacity(2*(profile.len()-1)*n_theta);
    for j in 0..profile.len()-1 { for i in 0..n_theta {
        let ni=(i+1)%n_theta;
        let(a,b,c,d)=(j*n_theta+i,j*n_theta+ni,(j+1)*n_theta+i,(j+1)*n_theta+ni);
        tris.push([a,b,d]);tris.push([a,d,c]);
    }}
    ShellMesh {nodes,tris}
}
/// Existing normalized church-bell research profile; not a measured specimen.
#[must_use]
pub fn canonical_church_bell_profile(scale_m:f64)->Vec<(f64,f64)> {
    [(0.10,1.00),(0.18,0.88),(0.24,0.72),(0.30,0.55),(0.38,0.38),
        (0.50,0.20),(0.65,0.08),(0.75,0.00)].iter()
        .map(|&(r,z)|(r*scale_m,z*scale_m)).collect()
}
