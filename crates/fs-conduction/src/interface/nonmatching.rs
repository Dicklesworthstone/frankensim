//! Common-refinement integration for two independently triangulated planar
//! contact traces. No volume remeshing, node welding, nearest-node matching,
//! temperature averaging or penalty approximation to perfect contact.
//!
//! Integrate `(Ta-Tb)(va-vb)/R` on the triangle intersections. Each overlap
//! triangle has three positive degree-two quadrature weights; both P1 bases
//! are evaluated at the SAME physical points. In exact arithmetic this gives
//! a symmetric positive-semidefinite contact operator with a constant null
//! mode and equal/opposite heat. Floating-point geometry/coverage tests are
//! explicitly tolerance-based, not certified intersection predicates.
use std::collections::{BTreeMap, BTreeSet};
use fs_exec::Cx;
use fs_sparse::Coo;
use crate::{ConductionError, ConductionMesh, ThermalBoundary};
use super::{InterfaceFacePair, InterfaceFlux, InterfaceResistance};

/// Explicit geometry/work policy for one planar contact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NonmatchingOptions {
    /// Maximum distance of any trace vertex from side A's reference plane, m.
    /// May be zero; may not exceed one millionth of the surface coordinate span.
    pub plane_tolerance_m: f64,
    /// Relative per-face coverage/overlap tolerance, in `(0, 1e-6]`.
    pub coverage_relative_tolerance: f64,
    /// Pair tests, including both same-side disjointness checks.
    pub max_pair_tests: usize,
    /// Maximum retained common-refinement triangles, not volume elements.
    pub max_overlap_triangles: usize,
}

/// Explicit, uniformly resistive planar interface. Face slots belong to the
/// supplied conduction mesh; A and B remain independent nodal traces.
#[derive(Debug, Clone, PartialEq)]
pub struct NonmatchingSurface {
    pub(super) name: String,
    pub(super) side_a: Vec<usize>,
    pub(super) side_b: Vec<usize>,
    pub(super) resistance: InterfaceResistance,
    options: NonmatchingOptions,
}
impl NonmatchingSurface {
    /// Declare both complete tessellations and a card-backed constant R''.
    /// Geometry, coverage and resource admission occur during binding.
    pub fn new(name: impl Into<String>, mut side_a: Vec<usize>, mut side_b: Vec<usize>,
        resistance: InterfaceResistance, options: NonmatchingOptions) -> Result<Self, ConductionError> {
        let name = name.into();
        if name.trim().is_empty() || side_a.is_empty() || side_b.is_empty() {
            return Err(error(&name, "name and both contact face sets must be nonempty"));
        }
        side_a.sort_unstable(); side_b.sort_unstable();
        if side_a.windows(2).chain(side_b.windows(2)).any(|p| p[0] == p[1])
            || side_a.iter().any(|i| side_b.binary_search(i).is_ok()) {
            return Err(error(&name, "each boundary face must have exactly one side/owner"));
        }
        if !(options.plane_tolerance_m.is_finite() && options.plane_tolerance_m >= 0.0
            && options.coverage_relative_tolerance.is_finite()
            && options.coverage_relative_tolerance > 0.0 && options.coverage_relative_tolerance <= 1e-6)
            || options.max_pair_tests == 0 || options.max_overlap_triangles == 0 {
            return Err(error(&name, "invalid explicit nonmatching geometry or work policy"));
        }
        let density = 1.0 / resistance.value_m2_k_per_w();
        if !density.is_finite() || density <= 0.0 {
            return Err(error(&name, "contact conductance density is not representable"));
        }
        Ok(Self { name, side_a, side_b, resistance, options })
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Face {
    slot: usize,
    vertices: [usize; 3],
    points: [[f64; 3]; 3],
    normal: [f64; 3],
    area: f64,
}
impl Face {
    fn read(mesh: &ConductionMesh, boundary: &ThermalBoundary, slot: usize,
        name: &str) -> Result<Self, ConductionError> {
        let face = mesh.boundary().get(slot).ok_or_else(|| error(name, "contact face slot is outside the mesh"))?;
        if boundary.condition_for(slot).is_some() {
            return Err(error(name, "a contact trace also carries an external boundary condition"));
        }
        Ok(Self { slot, vertices: face.vertices.map(|i| i as usize),
            points: face.vertices.map(|i| mesh.positions()[i as usize]),
            normal: face.outward_normal, area: face.area })
    }
    fn same_bits(&self, other: &Self) -> bool {
        self.slot == other.slot && self.vertices == other.vertices
            && self.area.to_bits() == other.area.to_bits()
            && self.points.iter().flatten().chain(&self.normal)
                .zip(other.points.iter().flatten().chain(&other.normal))
                .all(|(a,b)| a.to_bits() == b.to_bits())
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Stencil {
    vertices: [usize; 6],
    /// Signed A-minus-B basis evaluations at three degree-two quadrature points.
    basis: [[f64; 6]; 3],
    weight_m2: f64,
    /// Exact matching pairs stay on the unchanged matching assembly path.
    delegated: bool,
}
impl Stencil {
    fn jump(&self, values: &[f64], q: &[f64; 6]) -> f64 {
        let v = self.vertices;
        // Constant traces reproduce their exact jump, even when barycentric
        // coefficient summation would lose low bits of the absolute offset.
        (values[v[0]] - values[v[3]])
            + q[1]*(values[v[1]]-values[v[0]]) + q[2]*(values[v[2]]-values[v[0]])
            + q[4]*(values[v[4]]-values[v[3]]) + q[5]*(values[v[5]]-values[v[3]])
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(super) struct Bound {
    pub(super) name: String,
    resistance: InterfaceResistance,
    faces: Vec<Face>,
    stencils: Vec<Stencil>,
}

impl Bound {
    pub(super) fn build(cx: &Cx<'_>, mesh: &ConductionMesh, boundary: &ThermalBoundary,
        surface: &NonmatchingSurface, exact: &BTreeSet<(usize,usize)>) -> Result<Self, ConductionError> {
        let name = &surface.name;
        poll(cx,0)?;
        let count = surface.side_a.len().checked_add(surface.side_b.len())
            .ok_or_else(|| error(name,"contact face count overflow"))?;
        let tests = count.checked_mul(count.saturating_sub(1)).map(|n| n/2)
            .ok_or_else(|| error(name,"contact pair budget overflow"))?;
        if tests > surface.options.max_pair_tests {
            return Err(error(name,"nonmatching pair-test budget exhausted before geometry work"));
        }
        let a = surface.side_a.iter().map(|&s| Face::read(mesh,boundary,s,name)).collect::<Result<Vec<_>,_>>()?;
        let b = surface.side_b.iter().map(|&s| Face::read(mesh,boundary,s,name)).collect::<Result<Vec<_>,_>>()?;
        let av: BTreeSet<_> = a.iter().flat_map(|f| f.vertices).collect();
        if b.iter().flat_map(|f| f.vertices).any(|v| av.contains(&v)) {
            return Err(error(name,"contact sides must own distinct temperature nodes"));
        }
        let origin = a[0].points[0];
        let normal = a[0].normal;
        let axis = (0..3).max_by(|&i,&j| normal[i].abs().total_cmp(&normal[j].abs())).unwrap_or(0);
        let u = (axis+1)%3; let v = (axis+2)%3;
        let mut scale = 0.0_f64;
        for face in a.iter().chain(&b) { for point in face.points { for d in 0..3 {
            scale = scale.max(finite(point[d]-origin[d],name)?.abs());
        } } }
        if !scale.is_finite() || scale <= 0.0 || surface.options.plane_tolerance_m/scale > 1e-6 {
            return Err(error(name,"invalid planar span or plane tolerance larger than 1e-6 of the span"));
        }
        let area_scale = finite(scale*scale/normal[axis].abs(),name)?;
        if area_scale <= 0.0 { return Err(error(name,"projected contact area underflow")); }
        let project = |face: &Face, sign: f64| -> Result<Triangle,ConductionError> {
            let dot = face.normal.iter().zip(normal).map(|(x,y)| x*y).sum::<f64>();
            if !dot.is_finite() || sign*dot < 1.0-256.0*f64::EPSILON {
                return Err(error(name,"contact facets are not coplanar with consistent opposing normals"));
            }
            let mut points = [[0.0;2];3];
            let mut ids = face.vertices;
            for (i,p) in face.points.iter().enumerate() {
                let delta = [p[0]-origin[0],p[1]-origin[1],p[2]-origin[2]];
                let distance = finite(delta.iter().zip(normal).map(|(x,y)|x*y).sum(),name)?;
                if distance.abs() > surface.options.plane_tolerance_m {
                    return Err(error(name,"contact vertices leave the declared plane tolerance"));
                }
                points[i] = [delta[u]/scale,delta[v]/scale];
            }
            let det = cross(points[0],points[1],points[2]);
            if !det.is_finite() || det == 0.0 { return Err(error(name,"degenerate projected contact triangle")); }
            if det < 0.0 { points.swap(1,2); ids.swap(1,2); }
            Ok(Triangle { points, ids, slot:face.slot, area:face.area })
        };
        let ta = a.iter().map(|f|project(f,1.0)).collect::<Result<Vec<_>,_>>()?;
        let tb = b.iter().map(|f|project(f,-1.0)).collect::<Result<Vec<_>,_>>()?;
        let tolerance = surface.options.coverage_relative_tolerance;
        let mut at = 0;
        for side in [&ta,&tb] { for i in 0..side.len() { for j in i+1..side.len() {
            poll(cx,at)?; at+=1;
            let polygon = intersection(&side[i],&side[j]);
            let area = polygon_area(&polygon)*area_scale;
            if !area.is_finite() || area > tolerance*side[i].area.min(side[j].area) {
                return Err(error(name,"overlapping facets within one contact side"));
            }
        } } }
        let mut coverage = BTreeMap::<usize,f64>::new();
        let mut stencils = Vec::new();
        for left in &ta { for right in &tb {
            poll(cx,at)?; at+=1;
            let polygon = intersection(left,right);
            for i in 1..polygon.len().saturating_sub(1) {
                let points = [polygon[0],polygon[i],polygon[i+1]];
                let area = finite(0.5*cross(points[0],points[1],points[2])*area_scale,name)?;
                if area <= 0.0 { continue; }
                if stencils.len() >= surface.options.max_overlap_triangles {
                    return Err(error(name,"nonmatching overlap-triangle budget exhausted"));
                }
                for slot in [left.slot,right.slot] {
                    let sum = coverage.entry(slot).or_insert(0.0);
                    *sum = finite(*sum+area,name)?;
                }
                let mut basis = [[0.0;6];3];
                for (q,row) in basis.iter_mut().enumerate() {
                    let p = [0,1].map(|d| (2.0/3.0)*points[q][d]
                        +(1.0/6.0)*points[(q+1)%3][d]+(1.0/6.0)*points[(q+2)%3][d]);
                    let l = barycentric(left,p,name)?; let r = barycentric(right,p,name)?;
                    for k in 0..3 { row[k]=l[k]; row[k+3]=-r[k]; }
                    if row.iter().any(|x| x.abs()>1.0+32.0*tolerance) {
                        return Err(error(name,"unstable barycentric evaluation on a contact overlap"));
                    }
                }
                let weight_m2 = finite(area/3.0,name)?;
                if weight_m2 <= 0.0 { return Err(error(name,"contact quadrature weight underflow")); }
                let key = (left.slot.min(right.slot),left.slot.max(right.slot));
                stencils.push(Stencil { vertices:[left.ids[0],left.ids[1],left.ids[2],right.ids[0],right.ids[1],right.ids[2]],
                    basis, weight_m2, delegated:exact.contains(&key) });
            }
        } }
        for face in a.iter().chain(&b) {
            let area = coverage.get(&face.slot).copied().unwrap_or(0.0);
            if (area-face.area).abs() > tolerance*face.area {
                return Err(error(name,format!("contact face {} is not completely and singly covered: overlap {area}, face {} m2",face.slot,face.area)));
            }
        }
        poll(cx,at)?;
        Ok(Self { name:name.clone(), resistance:surface.resistance.clone(),
            faces:a.into_iter().chain(b).collect(), stencils })
    }

    pub(super) fn validate_for(&self, mesh: &ConductionMesh, boundary: &ThermalBoundary) -> Result<(),ConductionError> {
        for face in &self.faces {
            if !face.same_bits(&Face::read(mesh,boundary,face.slot,&self.name)?) {
                return Err(error(&self.name,"nonmatching contact mesh changed after binding"));
            }
        }
        Ok(())
    }
    pub(super) fn assemble_into(&self,cx:&Cx<'_>,coo:&mut Coo)->Result<(),ConductionError> {
        let g=1.0/self.resistance.value_m2_k_per_w();
        for (index,stencil) in self.stencils.iter().enumerate() {
            poll(cx,index)?;
            if stencil.delegated {continue;}
            for q in &stencil.basis { for i in 0..6 { for j in i..6 {
                let value=finite((g*stencil.weight_m2)*(q[i]*q[j]),&self.name)?;
                coo.push(stencil.vertices[i],stencil.vertices[j],value);
                if i!=j {coo.push(stencil.vertices[j],stencil.vertices[i],value);}
            } } }
        }
        Ok(())
    }
    fn check_field(&self,values:&[f64])->Result<(),ConductionError> {
        for face in &self.faces { for vertex in face.vertices {
            if vertex>=values.len() || !values[vertex].is_finite() {
                return Err(error(&self.name,"contact field has missing or nonfinite trace entries"));
            }
        } }
        Ok(())
    }
    /// Only the non-delegated portion; matching flux is merged by the owner.
    pub(super) fn flux(&self,values:&[f64])->Result<Option<InterfaceFlux>,ConductionError> {
        self.check_field(values)?;
        let mut area=0.0; let mut heat=0.0; let g=1.0/self.resistance.value_m2_k_per_w();
        for stencil in &self.stencils {
            if stencil.delegated {continue;}
            for q in &stencil.basis {
                area=finite(area+stencil.weight_m2,&self.name)?;
                heat=finite(heat+g*stencil.weight_m2*stencil.jump(values,q),&self.name)?;
            }
        }
        if area==0.0 {return Ok(None);}
        let conductance=finite(area*g,&self.name)?;
        if conductance<=0.0 {return Err(error(&self.name,"contact conductance underflow"));}
        Ok(Some(InterfaceFlux {interface:self.name.clone(),area_m2:area,
            conductance_w_per_k:conductance,mean_jump_k:finite(heat/conductance,&self.name)?,
            heat_rate_a_to_b_w:heat,card_identity:self.resistance.card_identity(),receipt:self.resistance.receipt().clone()}))
    }
    /// Total `lambda^T K_contact T`, including any exactly matched subpatch.
    pub(super) fn log_resistance_pullback(&self,cx:&Cx<'_>,temperature:&[f64],lambda:&[f64])->Result<f64,ConductionError> {
        self.check_field(temperature)?;self.check_field(lambda)?;
        let mut result=0.0;let g=1.0/self.resistance.value_m2_k_per_w();
        for (index,stencil) in self.stencils.iter().enumerate() {
            poll(cx,index)?;
            for q in &stencil.basis {
                result=finite(result+g*stencil.weight_m2*stencil.jump(temperature,q)*stencil.jump(lambda,q),&self.name)?;
            }
        }
        Ok(result)
    }
    pub(super) fn resistances(&self)->Vec<f64> {
        vec![self.resistance.value_m2_k_per_w();self.stencils.len()]
    }
}

struct Triangle {points:[[f64;2];3],ids:[usize;3],slot:usize,area:f64}
fn cross(a:[f64;2],b:[f64;2],c:[f64;2])->f64 {
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn intersection(a:&Triangle,b:&Triangle)->Vec<[f64;2]> {
    let mut polygon=a.points.to_vec();
    for edge in 0..3 {
        if polygon.is_empty(){break;}
        let u=b.points[edge];let v=b.points[(edge+1)%3];
        let mut next=Vec::with_capacity(6);
        let mut old=*polygon.last().expect("nonempty polygon");
        let mut d_old=cross(u,v,old);
        for &point in &polygon {
            let d=cross(u,v,point);
            if (d>=0.0)!=(d_old>=0.0) {
                let t=d_old/(d_old-d);
                next.push([old[0]+t*(point[0]-old[0]),old[1]+t*(point[1]-old[1])]);
            }
            if d>=0.0 {next.push(point);}
            old=point;d_old=d;
        }
        polygon=next;
    }
    polygon
}
fn polygon_area(points:&[[f64;2]])->f64 {
    (1..points.len().saturating_sub(1)).map(|i|0.5*cross(points[0],points[i],points[i+1])).sum()
}
fn barycentric(triangle:&Triangle,p:[f64;2],name:&str)->Result<[f64;3],ConductionError> {
    let [a,b,c]=triangle.points;let d=cross(a,b,c);
    let l1=finite(cross(a,p,c)/d,name)?;let l2=finite(cross(a,b,p)/d,name)?;
    Ok([finite(1.0-l1-l2,name)?,l1,l2])
}
fn finite(value:f64,name:&str)->Result<f64,ConductionError> {
    if value.is_finite(){Ok(value)}else{Err(error(name,"nonfinite contact integration arithmetic"))}
}
fn poll(cx:&Cx<'_>,at:usize)->Result<(),ConductionError> {
    cx.checkpoint().map_err(|_|ConductionError::Cancelled{stage:"nonmatching-contact",at})
}
fn error(name:&str,what:impl Into<String>)->ConductionError {
    ConductionError::Interface{interface:name.into(),what:what.into(),
        fix:"supply two complete, disjoint, oppositely oriented planar tessellations with independent nodes and adequate explicit geometry budgets".into()}
}
