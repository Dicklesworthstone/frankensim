use super::*;

pub(super) type Vector = [Iv; 3];
pub(super) fn sub(a: Vector, b: Vector) -> Vector {
    std::array::from_fn(|i| a[i].sub(b[i]))
}
pub(super) fn scale(a: Vector, b: Iv) -> Vector {
    a.map(|x| x.mul(b))
}
pub(super) fn dot(a: Vector, b: Vector) -> Iv {
    (0..3).fold(Iv::zero(), |s, i| s.add(a[i].mul(b[i])))
}
fn cross(a: Vector, b: Vector) -> Vector {
    [a[1].mul(b[2]).sub(a[2].mul(b[1])),
     a[2].mul(b[0]).sub(a[0].mul(b[2])),
     a[0].mul(b[1]).sub(a[1].mul(b[0]))]
}

/// Exact simplex P1 mass moments: measure/(m*(m+1)) *
/// (sum nodal squares + square of nodal sum), where m is the vertex count.
pub(super) fn integral_square(values: &[Iv], measure: Iv, denominator: f64) -> Iv {
    let sum = values.iter().copied().fold(Iv::zero(), Iv::add);
    let sq = values.iter().fold(Iv::zero(), |s, v| s.add(v.sq()));
    measure.mul(sq.add(sum.sq())).div_pos(Iv::point(denominator))
}

#[derive(Debug)]
pub(super) struct Cell {
    pub points: [Vector; 4],
    pub volume: Iv,
    pub gradient: Vector,
    pub faces: [usize; 4],
}
#[derive(Debug)]
pub(super) struct Face {
    pub vertices: [usize; 3],
    pub sides: Vec<(usize, usize)>, // element, opposite local vertex
    pub area: Iv,
    pub normal_area: Vector, // outward from sides[0]
    pub condition: Option<BoundaryCondition>,
}
impl Face {
    pub fn sign(&self, e: usize) -> f64 {
        if self.sides[0].0 == e { 1.0 } else { -1.0 }
    }
}

fn geometry(points: [Vector; 4]) -> Result<(Iv, [Vector; 4]), TetError> {
    let det = dot(sub(points[1], points[0]),
        cross(sub(points[2], points[0]), sub(points[3], points[0])));
    let positive = if det.lo > 0.0 { det } else if det.hi < 0.0 {
        signed(det, -1.0)
    } else { return Err(TetError::Invalid("tet determinant contains zero")); };
    let volume = positive.div_pos(Iv::point(6.0));
    if volume.is_unbounded() || volume.lo <= 0.0 { return Err(TetError::Unbounded); }
    let mut normals = [[Iv::zero(); 3]; 4];
    for (opposite, normal) in normals.iter_mut().enumerate() {
        let local: Vec<_> = (0..4).filter(|&i| i != opposite).collect();
        let [a, b, c] = [points[local[0]], points[local[1]], points[local[2]]];
        let raw = cross(sub(b, a), sub(c, a));
        let side = dot(raw, sub(points[opposite], a));
        let sign = if side.lo > 0.0 { -1.0 } else if side.hi < 0.0 { 1.0 }
            else { return Err(TetError::Invalid("ambiguous outward face orientation")); };
        *normal = scale(raw, Iv::point(sign).div_pos(Iv::point(2.0)));
    }
    Ok((volume, normals))
}

pub(super) fn build(
    problem: &TetProblem<'_>, candidate: &[f64], budget: FluxBudget,
    keep_going: &mut impl FnMut() -> bool,
) -> Result<(Vec<Cell>, Vec<Face>), TetError> {
    poll(keep_going)?;
    let n = problem.tets.len();
    if n == 0 || problem.vertices.is_empty() { return Err(TetError::Invalid("empty domain")); }
    if n > budget.max_cells || problem.vertices.len() > budget.max_cells.saturating_mul(4)
        || budget.max_iterations > 1_000_000 { return Err(TetError::Budget); }
    if candidate.len() != problem.vertices.len() || problem.source.len() != n
        || problem.conductivity.len() != n || problem.boundary.len() > n.saturating_mul(4) {
        return Err(TetError::Invalid("field or boundary length"));
    }
    for chunk in problem.vertices.chunks(256) {
        poll(keep_going)?;
        if chunk.iter().flatten().any(|x| !x.is_finite()) {
            return Err(TetError::Invalid("non-finite coordinates"));
        }
    }
    if candidate.iter().chain(problem.source).any(|v| !v.is_finite())
        || problem.conductivity.iter().any(|k| !k.is_finite() || *k <= 0.0) {
        return Err(TetError::Invalid("finite data and positive conductivity required"));
    }
    let mut declarations = BTreeMap::new();
    for boundary in problem.boundary {
        poll(keep_going)?;
        let mut order = [0, 1, 2];
        order.sort_by_key(|&i| boundary.vertices[i]);
        let key = order.map(|i| boundary.vertices[i]);
        if key[0] == key[1] || key[1] == key[2] || key[2] >= candidate.len() {
            return Err(TetError::Invalid("boundary vertex indices"));
        }
        let condition = match boundary.condition {
            BoundaryCondition::Dirichlet(values) => {
                let values = order.map(|i| values[i]);
                if values.iter().any(|v| !v.is_finite())
                    || (0..3).any(|i| candidate[key[i]] != values[i]) {
                    return Err(TetError::Invalid("candidate does not match prescribed Dirichlet trace"));
                }
                BoundaryCondition::Dirichlet(values)
            }
            BoundaryCondition::Neumann(q) => {
                if !q.is_finite() { return Err(TetError::Invalid("Neumann flux")); }
                BoundaryCondition::Neumann(q)
            }
            BoundaryCondition::Robin { h, reference } => {
                if !h.is_finite() || h <= 0.0 || reference.iter().any(|v| !v.is_finite()) {
                    return Err(TetError::Invalid("positive Robin h and finite reference required"));
                }
                BoundaryCondition::Robin { h, reference: order.map(|i| reference[i]) }
            }
        };
        if declarations.insert(key, condition).is_some() {
            return Err(TetError::Invalid("duplicate boundary declaration"));
        }
    }
    let mut cells = Vec::with_capacity(n);
    let mut faces = Vec::<Face>::new();
    let mut slots = BTreeMap::new();
    let mut seen = BTreeSet::new();
    for (e, tet) in problem.tets.iter().enumerate() {
        poll(keep_going)?;
        let mut key = *tet;
        key.sort_unstable();
        if key[3] >= candidate.len() || key.windows(2).any(|w| w[0] == w[1])
            || !seen.insert(key) { return Err(TetError::Invalid("invalid or duplicate tet")); }
        let points = tet.map(|i| problem.vertices[i].map(Iv::point));
        let (volume, normals) = geometry(points)?;
        let mut gradient = [Iv::zero(); 3];
        // Sum grad(lambda_i)=0 exactly: subtract a common temperature first.
        for i in 1..4 {
            let coefficient = Iv::point(candidate[tet[i]]).sub(Iv::point(candidate[tet[0]]))
                .div_pos(volume.scale_pos(3.0));
            for d in 0..3 { gradient[d] = gradient[d].sub(normals[i][d].mul(coefficient)); }
        }
        let mut cell_faces = [0; 4];
        for opposite in 0..4 {
            let mut key = [0; 3];
            let mut j = 0;
            for (i, &v) in tet.iter().enumerate() {
                if i != opposite { key[j] = v; j += 1; }
            }
            key.sort_unstable();
            let f = if let Some(&f) = slots.get(&key) {
                let face: &mut Face = &mut faces[f];
                if face.sides.len() != 1 || dot(face.normal_area, normals[opposite]).hi >= 0.0 {
                    return Err(TetError::Invalid("nonmanifold face or neighbors on the same side"));
                }
                if face.condition.is_some() {
                    return Err(TetError::Invalid("boundary condition on an interior face"));
                }
                face.sides.push((e, opposite));
                f
            } else {
                let area = dot(normals[opposite], normals[opposite]).sqrt();
                if area.is_unbounded() || area.lo <= 0.0 { return Err(TetError::Unbounded); }
                let f = faces.len();
                faces.push(Face { vertices: key, sides: vec![(e, opposite)], area,
                    normal_area: normals[opposite], condition: declarations.remove(&key) });
                slots.insert(key, f);
                f
            };
            cell_faces[opposite] = f;
        }
        cells.push(Cell { points, volume, gradient, faces: cell_faces });
    }
    if !declarations.is_empty() || faces.iter().any(|f| f.sides.len() == 1 && f.condition.is_none()) {
        return Err(TetError::Invalid("boundary declarations must partition the extracted exterior"));
    }
    Ok((cells, faces))
}
