//! Explicit original-mesh bridge attachment, sharing the shell point owner.
use super::*;
use fs_plate::shell::point::{triangle_offset_point_rows,triangle_point_rows};

#[derive(Clone, Debug)]
pub(super) struct Bridge {
    triangle: Option<usize>,
    nodes: [usize; 3],
    barycentric: [f64; 3],
    pub(super) direction: [f64; 3],
    offset_m: [f64; 3],
}

impl Bridge {
    /// Number of fields excludes the optional assembly component name.
    pub(super) fn record(r: &Records<'_>) -> (&'static str, usize) {
        if r.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
            == Some("bridge-triangle")
        {
            ("bridge-triangle", 10)
        } else {
            ("bridge", 4)
        }
    }

    pub(super) fn read_fields(v: &[&str], node_count: usize, triangles: &[[usize; 3]]) -> Result<Self>
    {
        let bridge = match v.len() {
            4 => {
                let node = index(v[0], node_count)?;
                Self { triangle: None, nodes: [node; 3], barycentric: [1.0,0.0,0.0],
                    direction: unit_direction(&v[1..4])?, offset_m: [0.0;3] }
            }
            10 => {
                let triangle = index(v[0], triangles.len())?;
                Self { triangle: Some(triangle), nodes: triangles[triangle],
                    barycentric: [scalar(v[1])?,scalar(v[2])?,scalar(v[3])?],
                    direction: unit_direction(&v[4..7])?,
                    offset_m: [scalar(v[7])?,scalar(v[8])?,scalar(v[9])?] }
            }
            _ => return Err("shell bridge has an invalid physical footprint record".into()),
        };
        // Direction, weights and the arm are pure input admission. Physical
        // membrane-spin/support admission needs the actual SI facet geometry.
        triangle_point_rows(bridge.barycentric,bridge.direction,[bridge.offset_m;3])?;
        Ok(bridge)
    }

    fn rows(&self, mesh: &ShellMesh, directors: &NodalDirectors,
        supports: &[ShellConstraint]) -> Result<[[f64;6];3]>
    {
        if self.nodes.iter().any(|&node|node>=mesh.nodes.len()) {
            return Err("shell bridge lost its original mesh node".into());
        }
        let rows = if let Some(triangle) = self.triangle {
            if mesh.tris.get(triangle) != Some(&self.nodes) {
                return Err("shell bridge lost its original triangle".into());
            }
            let facet = mesh.facet(triangle)?;
            triangle_offset_point_rows(self.barycentric,self.direction,self.offset_m,
                facet.frame,facet.gradient,self.nodes.map(|node|directors.vectors()[node]))?
        } else {
            triangle_point_rows(self.barycentric,self.direction,[[0.0;3];3])?
        };
        if !(0..3).any(|a| {
            let fixed = supports.iter().find(|s|s.node == self.nodes[a]);
            (0..6).any(|dof| rows[a][dof] != 0.0 && fixed.is_none_or(|s|!s.mask[dof]))
        }) {
            return Err("shell bridge has no unconstrained physical force or moment direction".into());
        }
        Ok(rows)
    }

    pub(super) fn force_port(&self, mesh: &ShellMesh, directors: &NodalDirectors,
        supports: &[ShellConstraint]) -> Result<Vec<f64>>
    {
        let extent = mesh.nodes.len().checked_mul(6).ok_or("shell bridge row extent overflows")?;
        let mut port = vec![0.0;extent];
        for (node, row) in self.nodes.into_iter().zip(self.rows(mesh,directors,supports)?) {
            for dof in 0..6 { port[6*node+dof] += row[dof]; }
        }
        Ok(port)
    }

    pub(super) fn position(&self, shell: &Spec, spec: &input::Spec) -> Result<[f64; 3]> {
        let mut point = [0.0;3];
        for a in 0..3 {
            let p = shell.position(self.nodes[a], spec);
            for i in 0..3 { point[i] += self.barycentric[a]*(p[i]+self.offset_m[i]); }
        }
        if point.iter().any(|v|!v.is_finite()) {
            return Err("shell bridge physical position overflows".into());
        }
        Ok(point)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn geometry() -> (ShellMesh,NodalDirectors) {
        let mesh = ShellMesh::new(vec![[1.0,0.0,0.0],[1.0,1.0,0.0],
            [0.0,0.0,0.0],[0.0,1.0,0.0]],vec![[2,0,3],[0,1,3]]).unwrap();
        let directors = NodalDirectors::new(&mesh,vec![[0.0,0.0,1.0];4],0.8).unwrap();
        (mesh,directors)
    }

    #[test]
    fn interior_and_offset_ports_keep_the_original_three_node_six_dof_rows() {
        let (mesh,directors) = geometry();
        let bridge = Bridge::read_fields(&["0","0.2","0.3","0.5","1","0","0","0","0","0.02"],
            4, &mesh.tris).unwrap();
        let row = bridge.force_port(&mesh,&directors,&[]).unwrap();
        for (node, weight) in [(2,0.2),(0,0.3),(3,0.5)] {
            assert_eq!(row[6*node],weight);
            assert_eq!(row[6*node+4],weight*0.02);
        }
        assert_eq!(&row[6..12],&[0.0;6]);
        let old = Bridge::read_fields(&["2","0","0","1"],4,&mesh.tris).unwrap();
        let mut expected = vec![0.0;24]; expected[14] = 1.0;
        assert_eq!(old.force_port(&mesh,&directors,&[]).unwrap(),expected);
    }

    #[test]
    fn admission_uses_the_entire_physical_row_including_rotational_supports() {
        let mesh = ShellMesh::new(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]],
            vec![[0,1,2]]).unwrap();
        let directors = NodalDirectors::new(&mesh,vec![[0.0,0.0,1.0];3],0.8).unwrap();
        let supports = [ShellConstraint { node: 0, mask: [true,true,true,false,false,false],
            frame: ConstraintFrame::Global }];
        let node = Bridge::read_fields(&["0","1","0","0"],3,&mesh.tris).unwrap();
        assert!(node.force_port(&mesh,&directors,&supports).is_err());
        let offset = ["0","1","0","0","1","0","0","0","0","0.02"];
        let offset = Bridge::read_fields(&offset,3,&mesh.tris).unwrap();
        assert!(offset.force_port(&mesh,&directors,&supports).is_ok());
        let fixed = [ShellConstraint { node: 0, mask: [true;6], frame: ConstraintFrame::Global }];
        assert!(offset.force_port(&mesh,&directors,&fixed).is_err());
        // Unconstrained zero-weight vertices cannot rescue a fixed point.
        let centered = ["0","0.5","0.5","0","1","0","0","0","0","0.02"];
        let centered = Bridge::read_fields(&centered,3,&mesh.tris).unwrap();
        assert!(centered.force_port(&mesh,&directors,&fixed).is_ok());
        for bad in [
            ["0","-0.2","0.7","0.5","1","0","0","0","0","0.02"],
            ["0","0.2","0.3","0.4","1","0","0","0","0","0.02"],
            ["1","0.2","0.3","0.5","1","0","0","0","0","0.02"],
        ] {
            assert!(Bridge::read_fields(&bad,3,&mesh.tris).is_err());
        }
        // A normal-axis arm torque uses membrane spin, never the free drill.
        let spin = Bridge::read_fields(&["0","1","0","0","0","1","0","0.02","0","0"],
            3,&mesh.tris).unwrap();
        let translations_fixed: Vec<_> = (0..3).map(|node| ShellConstraint { node,
            mask:[true,true,true,false,false,false],frame:ConstraintFrame::Global }).collect();
        assert!(spin.force_port(&mesh,&directors,&translations_fixed).is_err());
        let row = spin.force_port(&mesh,&directors,&[]).unwrap();
        assert!((0..3).all(|a|row[6*a+5] == 0.0));
        assert!(row[6]!=0.0 || row[7]!=0.0 || row[12]!=0.0 || row[13]!=0.0);
    }
}
