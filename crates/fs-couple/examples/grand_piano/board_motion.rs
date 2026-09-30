//! Cold, full-vector shell/plate motion at an explicitly supplied acoustic skin.
//! P1 translation and physical axial rotation, with u_skin = u + theta x arm.
//! Shapes remain in the SAME bare-board modal basis used by bridge mechanics.
//! A projection is along the actual facet normal, never a nearest-node snap.
use fs_plate::ShellMesh;
use super::MAX_BOARD_MODES;

/// Finite acoustic skin derived from the admitted board's section thicknesses.
#[path = "board_skin.rs"]
pub mod skin;

#[derive(Debug)]
pub struct MotionSurface {
    pub mesh: ShellMesh,
    /// Mode-major, node-major (u_x,u_y,u_z,theta_x,theta_y,theta_z).
    pub shapes: Vec<Vec<[f64; 6]>>,
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
        Ok(Self {mesh,shapes,projection_tree})
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
            best=distance.abs();selected=Some((tri,weights,g.frame[2].map(|n|distance*n)));
        }
        let (tri,weights,arm)=selected.ok_or("moving acoustic skin has no in-panel structural projection within its declared offset")?;
        let mut result=Vec::with_capacity(self.shapes.len());
        for mode in &self.shapes {
            let mut value=0.;
            for (i,node) in tri.into_iter().enumerate() {
                let q=mode[node];let rotation=cross([q[3],q[4],q[5]],arm);
                let u=std::array::from_fn(|c|q[c]+rotation[c]);
                value+=weights[i]*dot(normal,u);
            }
            if !value.is_finite() {return Err("skin motion projection overflow".into());}
            result.push(value);
        }
        Ok((result, work))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
