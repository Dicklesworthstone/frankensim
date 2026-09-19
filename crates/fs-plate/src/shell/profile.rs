//! Physical revolved midsurfaces with radial taper and explicit local features.
//! A cymbal, gong, bowl or annular diaphragm is data, never a solver branch.
//! Indentations modify the mechanical mesh; thickness relief modifies D and M.
//! These are stress-free reference shapes, NOT a simulation of hammer forming,
//! residual stress, work hardening, lathing damage or a particular manufacturer.
use super::{ShellMesh, assemble_shell_sections, ShellModel, ShellSupport};
use crate::{PlateError, PlateSection};

/// One measured or explicitly estimated meridian station.
#[derive(Debug, Clone, Copy)]
pub struct ProfileStation {
    /// Nonnegative radius [m], strictly increasing between stations.
    pub radius_m: f64,
    /// Midsurface height [m].
    pub height_m: f64,
    /// Physical thickness [m], positive. It is not a normal-map amplitude.
    pub thickness_m: f64,
}
/// Explicit axisymmetric manufacturing relief (lathe groove or thickness band).
#[derive(Debug, Clone, Copy)]
pub struct AnnularRelief {
    /// Ring center radius [m].
    pub radius_m: f64,
    /// Positive half width [m].
    pub half_width_m: f64,
    /// Signed midsurface shift at the ring center [m].
    pub height_delta_m: f64,
    /// Signed thickness change at the ring center [m].
    pub thickness_delta_m: f64,
}
/// Local C1 compact indentation. Center and support size are explicit;
/// no random feature field or proprietary manufacturing pattern is inferred.
#[derive(Debug, Clone, Copy)]
pub struct SurfaceIndentation {
    /// Planar center [m] in the profile's x/y coordinates.
    pub center_m: [f64;2],
    /// Positive support radius [m].
    pub radius_m: f64,
    /// Signed axial height change at the center [m].
    pub height_delta_m: f64,
    /// Signed thickness change at the center [m].
    pub thickness_delta_m: f64,
}
/// Input mesh-work ceilings. The eigensolve has its own separate budget.
#[derive(Debug, Clone, Copy)]
pub struct ProfileBudget {
    /// Maximum generated nodes.
    pub max_nodes: usize,
    /// Maximum generated facets.
    pub max_triangles: usize,
    /// Maximum node-by-feature evaluations.
    pub max_feature_evaluations: usize,
}
/// Generated shell with one section per facet and inspectable sampled thickness.
#[derive(Debug, Clone)]
pub struct ProfileShell {
    /// Actual deformed reference midsurface, not just a display mesh.
    pub mesh: ShellMesh,
    /// Isotropic sections in element order.
    pub sections: Vec<PlateSection>,
    /// Thickness at each generated vertex [m].
    pub nodal_thickness_m: Vec<f64>,
    /// Center vertex or inner annular ring; no support is imposed automatically.
    pub inner_nodes: Vec<usize>,
    /// Outer rim nodes; free edges remain free unless a caller constrains them.
    pub outer_nodes: Vec<usize>,
    /// Sum rho*h*A over the actual triangulated midsurface [kg].
    pub mass_kg: f64,
    /// Longest facet edge [m]. Compare with the smallest declared feature.
    pub max_edge_m: f64,
    /// Count of features whose support is narrower than twice max_edge_m.
    /// Nonzero is an explicit under-resolution warning, not a quality certificate.
    pub underresolved_features: usize,
}
fn bad(what:&'static str)->PlateError {PlateError::BadSection{what}}
fn bump(s:f64)->f64 { if s<1.0 {(1.0-s*s).powi(2)}else{0.0} }

/// Revolve supplied radial stations, including a unique center vertex when r=0.
/// No duplicate center triangles or seam vertices. Thickness is nodally sampled
/// and averaged per facet before the EXISTING section law constructs stiffness.
/// Increasing radial/azimuthal resolution is necessary for narrow hammer/lathe
/// details; `underresolved_features` reports a conservative warning.
/// # Errors
/// Invalid profile/material/features, nonpositive resulting thickness, budget
/// overflow, allocation failure or invalid generated facets.
#[allow(clippy::too_many_arguments)]
pub fn revolve(stations:&[ProfileStation], azimuths:usize,
    young_pa:f64, poisson:f64, density_kg_m3:f64,
    rings:&[AnnularRelief], dents:&[SurfaceIndentation], budget:ProfileBudget)
    ->Result<ProfileShell,PlateError> {
    if stations.len()<2 || azimuths<3 {return Err(bad("profile needs two stations and three azimuths"));}
    // Reuse section-law material admission, without guessing from alloy labels.
    PlateSection::isotropic(young_pa,poisson,stations[0].thickness_m,density_kg_m3)?;
    for (i,s) in stations.iter().enumerate() {
        if ![s.radius_m,s.height_m,s.thickness_m].iter().all(|x|x.is_finite())
            || s.radius_m<0.0 || s.thickness_m<=0.0 || (i>0 && s.radius_m<=stations[i-1].radius_m) {
            return Err(bad("profile radii must increase and geometry must be finite with positive thickness"));
        }
    }
    for r in rings { if ![r.radius_m,r.half_width_m,r.height_delta_m,r.thickness_delta_m].iter().all(|x|x.is_finite())
        || r.radius_m<0.0 || r.half_width_m<=0.0 { return Err(bad("invalid annular relief")); } }
    for d in dents { if ![d.center_m[0],d.center_m[1],d.radius_m,d.height_delta_m,d.thickness_delta_m].iter().all(|x|x.is_finite())
        || d.radius_m<=0.0 {return Err(bad("invalid surface indentation"));} }
    let center=stations[0].radius_m==0.0;
    let n=stations.len().checked_sub(usize::from(center)).and_then(|v|v.checked_mul(azimuths))
        .and_then(|v|v.checked_add(usize::from(center))).ok_or_else(||bad("profile size overflow"))?;
    let triangles=(stations.len()-1).checked_mul(2).and_then(|v|v.checked_sub(usize::from(center)))
        .and_then(|v|v.checked_mul(azimuths)).ok_or_else(||bad("profile size overflow"))?;
    let features=rings.len().checked_add(dents.len()).ok_or_else(||bad("profile feature overflow"))?;
    let work=n.checked_mul(features).ok_or_else(||bad("profile feature work overflow"))?;
    if n>budget.max_nodes || triangles>budget.max_triangles || work>budget.max_feature_evaluations {
        return Err(bad("profile exceeds its declared mesh/feature budget"));
    }
    let mut nodes=Vec::new();let mut thickness=Vec::new();let mut tris=Vec::new();
    nodes.try_reserve_exact(n).map_err(|_|bad("profile node allocation failed"))?;
    thickness.try_reserve_exact(n).map_err(|_|bad("profile thickness allocation failed"))?;
    tris.try_reserve_exact(triangles).map_err(|_|bad("profile facet allocation failed"))?;
    for (j,s) in stations.iter().enumerate() {
        let count=if center && j==0 {1}else{azimuths};
        for i in 0..count {
            let theta=2.0*core::f64::consts::PI*i as f64/azimuths as f64;
            let (x,y)=(s.radius_m*theta.cos(),s.radius_m*theta.sin());
            let(mut z,mut h)=(s.height_m,s.thickness_m);
            for r in rings {let b=bump((s.radius_m-r.radius_m).abs()/r.half_width_m);z+=b*r.height_delta_m;h+=b*r.thickness_delta_m;}
            for d in dents {let b=bump((x-d.center_m[0]).hypot(y-d.center_m[1])/d.radius_m);z+=b*d.height_delta_m;h+=b*d.thickness_delta_m;}
            if !z.is_finite() || !h.is_finite() || h<=0.0 {return Err(bad("profile feature produced invalid thickness or height"));}
            nodes.push([x,y,z]);thickness.push(h);
        }
    }
    let start=|j:usize|if center {if j==0 {0}else{1+(j-1)*azimuths}}else{j*azimuths};
    for j in 0..stations.len()-1 {for i in 0..azimuths {
        let ni=(i+1)%azimuths;
        if center && j==0 {tris.push([0,start(1)+i,start(1)+ni]);}
        else {let(a,b,c,d)=(start(j)+i,start(j)+ni,start(j+1)+i,start(j+1)+ni);
            tris.push([a,c,d]);tris.push([a,d,b]);}
    }}
    let mesh=ShellMesh::new(nodes,tris)?;
    let mut sections=Vec::new();sections.try_reserve_exact(triangles).map_err(|_|bad("profile section allocation failed"))?;
    let(mut mass,mut max_edge)=(0.0,0.0_f64);
    for (e,tri) in mesh.tris.iter().enumerate() {
        let h=(thickness[tri[0]]+thickness[tri[1]]+thickness[tri[2]])/3.0;
        let section=PlateSection::isotropic(young_pa,poisson,h,density_kg_m3)?;
        mass+=density_kg_m3*h*mesh.facet(e)?.area_m2;
        sections.push(section);
        for edge in 0..3 {let a=mesh.nodes[tri[edge]];let b=mesh.nodes[tri[(edge+1)%3]];
            max_edge=max_edge.max((a[0]-b[0]).hypot(a[1]-b[1]).hypot(a[2]-b[2]));}
    }
    if !mass.is_finite() {return Err(bad("profile mass overflow"));}
    let underresolved_features=rings.iter().filter(|r|r.half_width_m<2.0*max_edge).count()
        +dents.iter().filter(|d|d.radius_m<2.0*max_edge).count();
    Ok(ProfileShell {mesh,sections,nodal_thickness_m:thickness,
        inner_nodes:(0..if center{1}else{azimuths}).collect(),
        outer_nodes:(start(stations.len()-1)..n).collect(),mass_kg:mass,max_edge_m:max_edge,underresolved_features})
}
impl ProfileShell {
    /// Feed the actual variable-thickness geometry into the existing shell FEM.
    /// # Errors
    /// Shell admission/assembly failures; no implicit mounting constraints.
    pub fn assemble(&self,nodes:&[usize],support:ShellSupport)->Result<ShellModel,PlateError> {
        assemble_shell_sections(&self.mesh,&self.sections,nodes,support)
    }
}
