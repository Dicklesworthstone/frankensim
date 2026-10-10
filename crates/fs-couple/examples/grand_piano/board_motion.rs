//! Cold, full-vector shell/plate motion at an explicitly supplied acoustic skin.
//! P1 translation/rotation or the explicit edge-cubic flat-plate field,
//! with physical axial rotation and u_skin = u + theta x arm.
//! Shapes remain in the SAME bare-board modal basis used by bridge mechanics.
//! A projection is along the actual facet normal, never a nearest-node snap.
use fs_plate::ShellMesh;
use super::MAX_BOARD_MODES;

/// Finite acoustic skin derived from the admitted board's section thicknesses.
#[path = "board_skin.rs"]
pub mod skin;

/// Area-normalized degree-three triangle rule. A cubic displacement plus
/// its quadratic rotation crossed with a linear arm integrates exactly on
/// one structural facet. The signed centroid weight is for linear motion
/// projection, not a positive mass/energy quadrature or a cross-facet rule.
pub const EDGE_CUBIC_QUADRATURE: [([f64; 3], f64); 4] = [
    ([1. / 3.; 3], -27. / 48.),
    ([0.6, 0.2, 0.2], 25. / 48.),
    ([0.2, 0.6, 0.2], 25. / 48.),
    ([0.2, 0.2, 0.6], 25. / 48.),
];

/// A supplied physical force/displacement probe, independent of a modal basis.
/// Preparing the map once permits projecting a certified source slice before
/// reduction without allocating an oversized runtime `MotionSurface`.
#[derive(Clone, Copy, Debug)]
pub struct SourceBridgePort {
    pub triangle: usize,
    pub weights: [f64; 3],
    pub arm_m: [f64; 3],
    pub direction: [f64; 3],
}

/// The two admitted transverse directions of one source string course.
/// The primary must match the geometric board's existing bridge coefficient
/// in every source mode before that mode can be removed by reduction.
#[derive(Clone, Copy, Debug)]
pub struct SourceBridgeFrame {
    pub midi: u8,
    pub primary: SourceBridgePort,
    pub secondary: SourceBridgePort,
}

/// One structural site's exact linear map, shared by source force snapshots
/// and the final runtime motion surface. It only needs three nodal states.
pub struct SiteProjection {
    nodes: [usize; 3],
    port: SourceBridgePort,
    cubic: Option<([f64; 9], [[f64; 9]; 2])>,
}

impl SourceBridgePort {
    pub fn prepare(&self, mesh: &ShellMesh, edge_cubic: bool) -> Result<SiteProjection, String> {
        let &nodes = mesh.tris.get(self.triangle)
            .ok_or("motion projection has no such structural facet")?;
        if self.weights.iter().chain(&self.arm_m).chain(&self.direction).any(|v| !v.is_finite())
            || self.weights.iter().any(|b| *b < -1e-10 || *b > 1. + 1e-10)
            || (self.weights.iter().sum::<f64>() - 1.).abs() > 1e-8
            || (dot(self.direction, self.direction) - 1.).abs() > 1e-8 {
            return Err("motion projection requires in-triangle weights, a finite SI arm and a unit direction".into());
        }
        let cubic = if edge_cubic {
            let [Some(a), Some(b), Some(c)] = nodes.map(|node| mesh.nodes.get(node)) else {
                return Err("motion projection has an invalid structural node".into());
            };
            if a.iter().chain(b).chain(c).any(|v| !v.is_finite()) || a[2] != b[2] || a[2] != c[2] {
                return Err("edge-cubic motion requires a finite flat XY facet".into());
            }
            let x = [a[0], b[0], c[0]]; let y = [a[1], b[1], c[1]];
            Some((fs_plate::edge_cubic_transverse_shape(&x, &y, self.weights),
                fs_plate::edge_cubic_transverse_gradient_shape(&x, &y, self.weights)))
        } else { None };
        Ok(SiteProjection { nodes, port: *self, cubic })
    }
}

impl SiteProjection {
    pub fn nodes(&self) -> [usize; 3] { self.nodes }

    /// Project physical `(u_x,u_y,u_z,theta_x,theta_y,theta_z)` at the three
    /// nodes returned by `nodes()`. The second result bounds absolute terms
    /// for primary-port consistency checks without a unit-sized floor.
    pub fn project_nodal(&self, nodal: [[f64; 6]; 3]) -> Result<(f64, f64), String> {
        if nodal.iter().flatten().any(|v| !v.is_finite()) {
            return Err("motion projection contains a nonfinite nodal coordinate".into());
        }
        let SourceBridgePort { weights: barycentric, arm_m, direction, .. } = self.port;
        let (value, scale) = if let Some((shape, gradient)) = &self.cubic {
            if nodal.iter().any(|q| q[0] != 0. || q[1] != 0. || q[5] != 0.) {
                return Err("edge-cubic motion cannot discard in-plane or drilling coordinates".into());
            }
            let dofs: [f64; 9] = std::array::from_fn(|i| {
                let q = nodal[i / 3]; match i % 3 { 0 => q[2], 1 => -q[4], _ => q[3] }
            });
            let sample = |row: &[f64; 9]| row.iter().zip(dofs).map(|(a, b)| a * b).sum::<f64>();
            let w = sample(shape); let theta = [sample(&gradient[1]), -sample(&gradient[0]), 0.];
            let rotation = cross(theta, arm_m);
            let value = dot(direction, [rotation[0], rotation[1], w + rotation[2]]);
            let dx_scale = (direction[2] * arm_m[0]).abs() + (direction[0] * arm_m[2]).abs();
            let dy_scale = (direction[2] * arm_m[1]).abs() + (direction[1] * arm_m[2]).abs();
            let scale = (0..9).map(|i| dofs[i].abs() * (direction[2].abs() * shape[i].abs()
                + dx_scale * gradient[0][i].abs() + dy_scale * gradient[1][i].abs())).sum();
            (value, scale)
        } else {
            let mut value = 0.; let mut scale = 0.;
            // Preserve the original P1 projection's value arithmetic.
            for (i, q) in nodal.into_iter().enumerate() {
                let rotation = cross([q[3], q[4], q[5]], arm_m);
                let u = std::array::from_fn(|c| q[c] + rotation[c]);
                value += barycentric[i] * dot(direction, u);
                let rotation_scale = [(q[4] * arm_m[2]).abs() + (q[5] * arm_m[1]).abs(),
                    (q[5] * arm_m[0]).abs() + (q[3] * arm_m[2]).abs(),
                    (q[3] * arm_m[1]).abs() + (q[4] * arm_m[0]).abs()];
                scale += barycentric[i].abs() * (0..3).map(|c|
                    direction[c].abs() * (q[c].abs() + rotation_scale[c])).sum::<f64>();
            }
            (value, scale)
        };
        if !value.is_finite() || !scale.is_finite() { return Err("skin motion projection overflow".into()); }
        Ok((value, scale))
    }
}

#[derive(Debug)]
pub struct MotionSurface {
    pub mesh: ShellMesh,
    /// Mode-major, node-major (u_x,u_y,u_z,theta_x,theta_y,theta_z).
    pub shapes: Vec<Vec<[f64; 6]>>,
    edge_cubic: bool,
    projection_tree: Vec<ProjectionNode>,
}
fn dot(a:[f64;3],b:[f64;3])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}
fn cross(a:[f64;3],b:[f64;3])->[f64;3] {
    [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]
}
#[derive(Clone, Copy, Debug)]
struct Bounds { min: [f64; 3], max: [f64; 3] }
impl Bounds {
    fn of_triangle(mesh: &ShellMesh, element: usize) -> Self {
        let tri = mesh.tris[element];
        let min = std::array::from_fn(|c| tri.iter().map(|&i| mesh.nodes[i][c]).fold(f64::INFINITY, f64::min));
        let max = std::array::from_fn(|c| tri.iter().map(|&i| mesh.nodes[i][c]).fold(f64::NEG_INFINITY, f64::max));
        Self { min, max }
    }
    fn union(self, other: Self) -> Self {
        Self { min: std::array::from_fn(|c| self.min[c].min(other.min[c])),
            max: std::array::from_fn(|c| self.max[c].max(other.max[c])) }
    }
    fn contains_offset(self, point: [f64; 3], offset: f64) -> bool {
        (0..3).all(|c| {
            // The exact facet gate admits barycentric roundoff of 1e-10.
            // Enlarge this broad phase by its worst coordinate extrapolation;
            // a candidate accepted by the old scan must never be omitted.
            let reach = offset + 3e-10 * (self.max[c] - self.min[c]) + 1e-12;
            point[c] >= self.min[c] - reach && point[c] <= self.max[c] + reach
        })
    }
}
#[derive(Debug)]
struct ProjectionNode {
    bounds: Bounds,
    children: Option<(usize, usize)>,
    elements: Vec<usize>,
}
fn build_projection_tree(mesh: &ShellMesh, elements: &mut [usize], nodes: &mut Vec<ProjectionNode>) -> usize {
    let bounds = elements.iter().map(|&e| Bounds::of_triangle(mesh, e))
        .reduce(Bounds::union).expect("admitted shell has triangles");
    let index = nodes.len();
    nodes.push(ProjectionNode { bounds, children: None, elements: Vec::new() });
    if elements.len() <= 8 {
        nodes[index].elements.extend_from_slice(elements);
    } else {
        let axis = (0..3).max_by(|&a, &b| (bounds.max[a] - bounds.min[a])
            .total_cmp(&(bounds.max[b] - bounds.min[b]))).unwrap();
        elements.sort_unstable_by(|&a, &b| {
            let midpoint = |e| {
                let box_ = Bounds::of_triangle(mesh, e);
                box_.min[axis] / 2.0 + box_.max[axis] / 2.0
            };
            midpoint(a).total_cmp(&midpoint(b)).then(a.cmp(&b))
        });
        let (left, right) = elements.split_at_mut(elements.len() / 2);
        let left = build_projection_tree(mesh, left, nodes);
        let right = build_projection_tree(mesh, right, nodes);
        nodes[index].children = Some((left, right));
    }
    index
}
impl MotionSurface {
    pub fn new(mesh:ShellMesh,shapes:Vec<Vec<[f64;6]>>)->Result<Self,String> {
        if mesh.nodes.len()>20_000 || mesh.tris.len()>40_000
            || !(1..=MAX_BOARD_MODES).contains(&shapes.len())
            || shapes.iter().any(|m|m.len()!=mesh.nodes.len()
                || m.iter().flatten().any(|v|!v.is_finite())) {
            return Err("acoustic motion must retain every finite nodal coordinate and admitted mode".into());
        }
        // ShellMesh fields are public; validate the actual supplied instance.
        let mesh=ShellMesh::new(mesh.nodes,mesh.tris).map_err(|e|e.to_string())?;
        let mut projection_tree = Vec::new();
        let mut elements: Vec<_> = (0..mesh.tris.len()).collect();
        build_projection_tree(&mesh, &mut elements, &mut projection_tree);
        Ok(Self {mesh,shapes,edge_cubic:false,projection_tree})
    }
    /// Retain the flat plate's existing edge-cubic field in the same modal
    /// basis. Nodal storage still uses physical rotations: wx=-theta_y,
    /// wy=theta_x. Interior rotations come from the analytic cubic gradient,
    /// not P1 interpolation of those nodal rotations. This is not a shell
    /// displacement law and admits no in-plane/drilling DOFs or crown.
    pub fn new_edge_cubic(mesh:ShellMesh,shapes:Vec<Vec<[f64;6]>>)->Result<Self,String> {
        let mut surface=Self::new(mesh,shapes)?;
        if surface.mesh.nodes.iter().any(|p|p[2]!=surface.mesh.nodes[0][2])
            || surface.shapes.iter().flatten().any(|q|q[0]!=0. || q[1]!=0. || q[5]!=0.) {
            return Err("edge-cubic motion requires a flat XY plate with transverse displacement and physical slope rotations".into());
        }
        surface.edge_cubic=true;
        Ok(surface)
    }
    pub fn is_edge_cubic(&self)->bool { self.edge_cubic }
    /// Project a known structural site, including its supplied physical arm.
    /// The first vector is the motion/effort coefficient in each retained mode;
    /// the second is an absolute-term scale for relative roundoff checks, with
    /// no unit-sized floor that could hide disagreement in very small modes.
    /// Both mechanics and skin pressure must use this same linear map.
    pub fn project_at(&self,element:usize,barycentric:[f64;3],arm_m:[f64;3],direction:[f64;3])
        ->Result<(Vec<f64>,Vec<f64>),String> {
        let projection=SourceBridgePort {triangle:element,weights:barycentric,arm_m,direction}
            .prepare(&self.mesh,self.edge_cubic)?;
        let mut result=Vec::with_capacity(self.shapes.len());
        let mut scales=Vec::with_capacity(self.shapes.len());
        for mode in &self.shapes {
            let [Some(a),Some(b),Some(c)]=projection.nodes().map(|node|mode.get(node)) else {
                return Err("motion projection is missing a retained nodal coordinate".into());
            };
            let (value,scale)=projection.project_nodal([*a,*b,*c])?;
            result.push(value);scales.push(scale);
        }
        Ok((result,scales))
    }
    /// Incremental normal motion at a skin point, including through-thickness
    /// rotation. The nearest admissible facet-normal projection must lie inside
    /// a real structural triangle and within the declared offset [m]. No
    /// extrapolation across a hole, invented skin thickness or CAD repair.
    pub fn normal_weights(&self,point:[f64;3],normal:[f64;3],offset_limit_m:f64)
        ->Result<Vec<f64>,String> {
        self.normal_weights_with_work(point, normal, offset_limit_m).map(|(weights, _)| weights)
    }
    /// The same exact projection, with its bounded broad-phase visit count.
    /// A caller preparing many acoustic panels can enforce a real work budget
    /// instead of pessimistically multiplying every panel by every triangle.
    pub fn normal_weights_with_work(&self,point:[f64;3],normal:[f64;3],offset_limit_m:f64)
        ->Result<(Vec<f64>,usize),String> {
        if point.iter().chain(&normal).any(|v|!v.is_finite())
            || (dot(normal,normal)-1.).abs()>1e-8 || !offset_limit_m.is_finite()
            || !(0.0..=0.05).contains(&offset_limit_m) {
            return Err("skin motion requires finite SI point, unit normal and offset in [0,0.05] m".into());
        }
        // Any admissible orthogonal projection lies inside its triangle's
        // bounds enlarged by the allowed distance. Sort back to source order
        // to preserve nearest-distance and equal-distance tie behavior.
        let mut stack = vec![0];
        let mut candidates = Vec::new();
        let mut work = 0usize;
        let margin = offset_limit_m;
        while let Some(index) = stack.pop() {
            work += 1;
            let node = &self.projection_tree[index];
            if !node.bounds.contains_offset(point, margin) { continue; }
            if let Some((left, right)) = node.children {
                stack.push(right); stack.push(left);
            } else {
                for &element in &node.elements {
                    work += 1;
                    if Bounds::of_triangle(&self.mesh, element).contains_offset(point, margin) {
                        candidates.push(element);
                    }
                }
            }
        }
        candidates.sort_unstable();
        work += candidates.len();
        let mut selected=None;let mut best=f64::INFINITY;
        for element in candidates {
            let tri = self.mesh.tris[element];
            let g=self.mesh.facet(element).map_err(|e|e.to_string())?;
            let d=std::array::from_fn(|c|point[c]-self.mesh.nodes[tri[0]][c]);
            let distance=dot(d,g.frame[2]);
            if distance.abs()>offset_limit_m+1e-12 || distance.abs()>=best {continue;}
            let x=dot(d,g.frame[0]);let y=dot(d,g.frame[1]);
            let weights:[f64;3]=std::array::from_fn(|i|
                (if i==0 {1.} else {0.})+g.gradient[i][0]*x+g.gradient[i][1]*y);
            if weights.iter().any(|b|!b.is_finite() || *b < -1e-10 || *b > 1.+1e-10) {continue;}
            best=distance.abs();selected=Some((element,weights,g.frame[2].map(|n|distance*n)));
        }
        let (element,weights,arm)=selected.ok_or("moving acoustic skin has no in-panel structural projection within its declared offset")?;
        let (result,_)=self.project_at(element,weights,arm,normal)?;
        Ok((result, work))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn cubic_mode() -> MotionSurface {
        let mesh=ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],vec![[0,1,2]]).unwrap();
        // Exact existing cubic reconstruction of nodal x^3 values/slopes:
        // w=x^3-(1-x-y)*x*y. Include a tiny second mode for relative checks.
        let q=vec![[0.;6],[0.,0.,1.,0.,-3.,0.],[0.;6]];
        let small=q.iter().map(|a|a.map(|v|v*1e-12)).collect();
        MotionSurface::new_edge_cubic(mesh,vec![q,small]).unwrap()
    }
    #[test]
    fn source_site_map_commutes_with_modal_reduction_and_preserves_force_work() {
        for edge_cubic in [false, true] {
            let motion = if edge_cubic {
                let mut motion = cubic_mode();
                motion.shapes[1] = vec![[0.; 6], [0.,0.,-0.3,0.4,1.2,0.],
                    [0.,0.,0.2,-0.5,0.7,0.]];
                motion
            } else {
                let mesh = ShellMesh::new(vec![[0.,0.,0.], [1.,0.,0.03], [0.,1.,0.02]],
                    vec![[0,1,2]]).unwrap();
                MotionSurface::new(mesh, vec![
                    vec![[1.,2.,3.,4.,5.,6.]; 3],
                    vec![[-0.2,0.3,0.7,-0.8,0.4,0.9]; 3],
                ]).unwrap()
            };
            for direction in [[0.,0.,1.], [0.6,0.8,0.], [0.,0.6,0.8]] {
                let port = SourceBridgePort { triangle: 0, weights: [0.17,0.29,0.54],
                    arm_m: [0.014,-0.009,0.02], direction };
                let projection = port.prepare(&motion.mesh, edge_cubic).unwrap();
                let (runtime, scales) = motion.project_at(port.triangle, port.weights,
                    port.arm_m, port.direction).unwrap();
                let mut source = Vec::new();
                for (i, mode) in motion.shapes.iter().enumerate() {
                    let (value, scale) = projection.project_nodal(projection.nodes().map(|node| mode[node])).unwrap();
                    assert_eq!(value.to_bits(), runtime[i].to_bits());
                    assert_eq!(scale.to_bits(), scales[i].to_bits());
                    source.push(value);
                }
                let coordinates = [0.31,-0.27];
                let mixed = projection.nodes().map(|node| std::array::from_fn(|c|
                    motion.shapes.iter().zip(coordinates).map(|(mode,q)| mode[node][c]*q).sum()));
                let (displacement, scale) = projection.project_nodal(mixed).unwrap();
                let modal: f64 = source.iter().zip(coordinates).map(|(b,q)| b*q).sum();
                assert!((displacement-modal).abs() < 1e-14*scale.max(modal.abs()));
                let force = 0.42;
                let modal_work: f64 = source.iter().zip(coordinates).map(|(b,q)| (force*b)*q).sum();
                assert!((force*displacement-modal_work).abs() < 1e-14*scale.max(modal.abs()));
            }
        }
    }
    #[test]
    fn cubic_motion_uses_analytic_rotations_and_the_actual_attachment_arm() {
        let motion=cubic_mode();
        assert!(motion.is_edge_cubic());
        for [x,y] in [[0.,0.],[1.,0.],[0.37,0.],[0.23,0.41],[0.6,0.4]] {
            let w=x*x*x-(1.-x-y)*x*y;
            let wx=3.*x*x+2.*x*y+y*y-y;let wy=x*x+2.*x*y-x;
            for arm in [[0.;3],[0.,0.,0.012],[0.018,-0.013,0.007]] {
                let u=[-arm[2]*wx,-arm[2]*wy,w+arm[0]*wx+arm[1]*wy];
                for direction in [[0.,0.,1.],[1.,0.,0.],[0.6,0.8,0.]] {
                    let (actual,scale)=motion.project_at(0,[1.-x-y,x,y],arm,direction).unwrap();
                    let expected=dot(direction,u);
                    assert!((actual[0]-expected).abs()<1e-14);
                    assert!((actual[1]-1e-12*expected).abs()<1e-26);
                    assert!((scale[1]-1e-12*scale[0]).abs()<1e-26);
                    for i in 0..2 {assert!(actual[i].abs()<=scale[i]*(1.+1e-14));}
                }
            }
        }
        let at=motion.project_at(0,[0.4,0.23,0.37],[0.,0.,0.009],[1.,0.,0.]).unwrap().0;
        let searched=motion.normal_weights([0.23,0.37,0.009],[1.,0.,0.],0.01).unwrap();
        for (a,b) in at.iter().zip(searched) {assert!((a-b).abs()<1e-14*a.abs().max(1e-30));}
        let linear=MotionSurface::new(motion.mesh.clone(),motion.shapes.clone()).unwrap();
        assert!(!linear.is_edge_cubic());
        let legacy=linear.project_at(0,[0.4,0.23,0.37],[0.,0.,0.009],[1.,0.,0.]).unwrap().0;
        assert!((legacy[0]-at[0]).abs()>1e-4,"interior gradient must not use P1 nodal slopes");
    }
    #[test]
    fn cubic_motion_refuses_crown_or_discarded_dofs_and_invalid_known_sites() {
        let motion=cubic_mode();
        let mut crown=motion.mesh.clone();crown.nodes[2][2]=0.001;
        assert!(MotionSurface::new_edge_cubic(crown,motion.shapes.clone()).is_err());
        for coordinate in [0,1,5] {
            let mut shapes=motion.shapes.clone();shapes[0][0][coordinate]=0.1;
            assert!(MotionSurface::new_edge_cubic(motion.mesh.clone(),shapes).is_err());
        }
        for (element,bary,arm,direction) in [
            (1,[0.4,0.3,0.3],[0.;3],[0.,0.,1.]),
            (0,[-0.1,0.6,0.5],[0.;3],[0.,0.,1.]),
            (0,[0.4,0.3,0.4],[0.;3],[0.,0.,1.]),
            (0,[0.4,0.3,0.3],[f64::NAN,0.,0.],[0.,0.,1.]),
            (0,[0.4,0.3,0.3],[0.;3],[0.,0.,0.]),
        ] {assert!(motion.project_at(element,bary,arm,direction).is_err());}
    }
    #[test]
    fn cubic_triangle_rule_integrates_every_monomial_through_degree_three() {
        let factorial=[1.,1.,2.,6.,24.,120.];
        for a in 0..=3 {for b in 0..=3-a {for c in 0..=3-a-b {
            let actual:f64=EDGE_CUBIC_QUADRATURE.iter().map(|(q,weight)|
                weight*q[0].powi(a as i32)*q[1].powi(b as i32)*q[2].powi(c as i32)).sum();
            let expected=2.*factorial[a]*factorial[b]*factorial[c]/factorial[a+b+c+2];
            assert!((actual-expected).abs()<1e-14);
        }}}
    }
    #[test]
    fn tilted_skin_retains_full_vector_rigid_motion_and_offset_rotation() {
        let mesh=ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.02],[0.,1.,0.01]],vec![[0,1,2]]).unwrap();
        let theta=[0.12,-0.08,0.04];let shift=[0.01,0.02,-0.03];
        let mode=mesh.nodes.iter().map(|&p| {
            let u=cross(theta,p);
            [u[0]+shift[0],u[1]+shift[1],u[2]+shift[2],theta[0],theta[1],theta[2]]
        }).collect();
        let normal=mesh.facet(0).unwrap().frame[2];
        let point=std::array::from_fn(|c|0.2*mesh.nodes[1][c]+0.3*mesh.nodes[2][c]+0.008*normal[c]);
        let surface=MotionSurface::new(mesh,vec![mode]).unwrap();
        let skin_normal=[1.,0.,0.];
        let actual=surface.normal_weights(point,skin_normal,0.01).unwrap()[0];
        assert!((actual-shift[0]-cross(theta,point)[0]).abs()<1e-14);
        let reverse=surface.normal_weights(point,[-1.,0.,0.],0.01).unwrap()[0];
        assert_eq!(actual,-reverse);
    }
    #[test]
    fn nodal_shapes_and_holes_never_receive_a_default_or_nearest_node() {
        let mesh=ShellMesh::new(vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]],vec![[0,1,2]]).unwrap();
        assert!(MotionSurface::new(mesh.clone(),vec![vec![[0.;6];2]]).is_err());
        let s=MotionSurface::new(mesh,vec![vec![[0.,0.,1.,0.,0.,0.];3]]).unwrap();
        for p in [[0.8,0.8,0.],[-0.1,0.1,0.],[0.2,0.2,0.03],[f64::NAN,0.,0.]] {
            assert!(s.normal_weights(p,[0.,0.,1.],0.01).is_err());
        }
        assert_eq!(s.normal_weights([0.2,0.3,0.005],[0.,0.,1.],0.01).unwrap(),vec![1.]);
    }
    #[test]
    fn spatial_projection_finds_a_distant_last_facet_without_all_pairs_work() {
        let mut nodes = Vec::new();
        let mut tris = Vec::new();
        let mut shape = Vec::new();
        for element in 0..128 {
            let x = 2.0 * element as f64;
            let start = nodes.len();
            nodes.extend([[x,0.0,0.0],[x+1.0,0.0,0.0],[x,1.0,0.0]]);
            tris.push([start,start+1,start+2]);
            shape.extend([[0.0,0.0,element as f64,0.0,0.0,0.0];3]);
        }
        let surface = MotionSurface::new(ShellMesh::new(nodes,tris).unwrap(),vec![shape]).unwrap();
        let (weights,work) = surface.normal_weights_with_work([254.2,0.3,0.005],[0.0,0.0,1.0],0.01).unwrap();
        assert!((weights[0]-127.0).abs()<1e-10);
        assert!(work < 128, "spatial projection inspected {work} nodes/candidates");
        assert!(surface.normal_weights([253.5,0.3,0.005],[0.0,0.0,1.0],0.01).is_err());
    }
    #[test]
    fn spatial_projection_keeps_the_nearest_overlapping_structural_face() {
        let mesh = ShellMesh::new(
            vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],
                 [0.0,0.0,0.01],[1.0,0.0,0.01],[0.0,1.0,0.01]],
            vec![[0,1,2],[3,4,5]]).unwrap();
        let mut shape = vec![[0.0;6];6];
        for q in &mut shape[..3] {q[2]=1.0;}
        for q in &mut shape[3..] {q[2]=2.0;}
        let surface = MotionSurface::new(mesh,vec![shape]).unwrap();
        let upper = surface.normal_weights([0.2,0.3,0.006],[0.0,0.0,1.0],0.01).unwrap();
        let lower = surface.normal_weights([0.2,0.3,0.004],[0.0,0.0,1.0],0.01).unwrap();
        assert!((upper[0]-2.0).abs()<1e-12);
        assert!((lower[0]-1.0).abs()<1e-12);
    }
}
