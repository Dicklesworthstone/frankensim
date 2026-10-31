//! Explicit original-mesh bridge attachment, sharing the shell point owner.
use super::*;
use fs_plate::shell::point::triangle_point_rows;

#[derive(Clone, Debug)]
pub(super) struct Bridge {
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

    pub(super) fn read_fields(v: &[&str], node_count: usize, triangles: &[[usize; 3]],
        supports: &[ShellConstraint]) -> Result<Self>
    {
        let bridge = match v.len() {
            4 => {
                let node = index(v[0], node_count)?;
                Self { nodes: [node; 3], barycentric: [1.0,0.0,0.0],
                    direction: unit_direction(&v[1..4])?, offset_m: [0.0;3] }
            }
            10 => {
                let triangle = index(v[0], triangles.len())?;
                Self { nodes: triangles[triangle],
                    barycentric: [scalar(v[1])?,scalar(v[2])?,scalar(v[3])?],
                    direction: unit_direction(&v[4..7])?,
                    offset_m: [scalar(v[7])?,scalar(v[8])?,scalar(v[9])?] }
            }
            _ => return Err("shell bridge has an invalid physical footprint record".into()),
        };
        let rows = bridge.rows()?;
        if !(0..3).any(|a| {
            let fixed = supports.iter().find(|s|s.node == bridge.nodes[a]);
            (0..6).any(|dof| rows[a][dof] != 0.0 && fixed.is_none_or(|s|!s.mask[dof]))
        }) {
            return Err("shell bridge has no unconstrained physical force or moment direction".into());
        }
        Ok(bridge)
    }

    fn rows(&self) -> Result<[[f64; 6]; 3]> {
        Ok(triangle_point_rows(self.barycentric, self.direction, [self.offset_m;3])?)
    }

    pub(super) fn force_port(&self, node_count: usize) -> Result<Vec<f64>> {
        let extent = node_count.checked_mul(6).ok_or("shell bridge row extent overflows")?;
        let mut port = vec![0.0;extent];
        for (node, row) in self.nodes.into_iter().zip(self.rows()?) {
            if node >= node_count { return Err("shell bridge lost its original mesh node".into()); }
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

    #[test]
    fn interior_and_offset_ports_keep_the_original_three_node_six_dof_rows() {
        let bridge = Bridge::read_fields(&["0","0.2","0.3","0.5","1","0","0","0","0","0.02"],
            4, &[[2,0,3]], &[]).unwrap();
        let row = bridge.force_port(4).unwrap();
        for (node, weight) in [(2,0.2),(0,0.3),(3,0.5)] {
            assert_eq!(row[6*node],weight);
            assert_eq!(row[6*node+4],weight*0.02);
        }
        assert_eq!(&row[6..12],&[0.0;6]);
        let old = Bridge::read_fields(&["2","0","0","1"],4,&[[2,0,3]],&[]).unwrap();
        let mut expected = vec![0.0;24]; expected[14] = 1.0;
        assert_eq!(old.force_port(4).unwrap(),expected);
    }

    #[test]
    fn admission_uses_the_entire_physical_row_including_rotational_supports() {
        let supports = [ShellConstraint { node: 0, mask: [true,true,true,false,false,false],
            frame: ConstraintFrame::Global }];
        assert!(Bridge::read_fields(&["0","1","0","0"],3,&[[0,1,2]],&supports).is_err());
        let offset = ["0","1","0","0","1","0","0","0","0","0.02"];
        assert!(Bridge::read_fields(&offset,3,&[[0,1,2]],&supports).is_ok());
        let fixed = [ShellConstraint { node: 0, mask: [true;6], frame: ConstraintFrame::Global }];
        assert!(Bridge::read_fields(&offset,3,&[[0,1,2]],&fixed).is_err());
        // Unconstrained zero-weight vertices cannot rescue a fixed point.
        let centered = ["0","0.5","0.5","0","1","0","0","0","0","0.02"];
        assert!(Bridge::read_fields(&centered,3,&[[0,1,2]],&fixed).is_ok());
        for bad in [
            ["0","-0.2","0.7","0.5","1","0","0","0","0","0.02"],
            ["0","0.2","0.3","0.4","1","0","0","0","0","0.02"],
            ["1","0.2","0.3","0.5","1","0","0","0","0","0.02"],
        ] {
            assert!(Bridge::read_fields(&bad,3,&[[0,1,2]],&[]).is_err());
        }
    }
}
