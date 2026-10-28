//! Conservative capillary forces on fixed-connectivity triangular surfaces.
//!
//! Positions are SI metres and prescribed isotropic surface tension is N/m.
//! Open patches retain their boundary line forces; closed patches assemble the
//! same surface-stress divergence without a curvature reconstruction.
//!
//! The endpoint force is a discrete gradient of triangle area through the Gram
//! invariants of its two edges. It is work-conjugate to the finite area change
//! and conserves linear and midpoint angular momentum. This is the discrete-
//! gradient principle of Gonzalez (1996), doi:10.1007/BF02440162; the triangle
//! formulas below follow directly from the Gram determinant.
//!
//! A prescribed coefficient change supplies free energy A_mean * delta_sigma.
//! It is not heat: a temperature-dependent tension needs a separate surface
//! entropy/internal-energy model. No surfactant transport, bending elasticity,
//! topology changes, intermediate-path validity or time integrator is implied.
//! The triangle leaves allocate nothing. Patch assembly is synchronous and
//! linear in the supplied vertices/faces; callers own tiling and cancellation.

type Vector = [f64; 3];

/// Invalid geometry, coefficient field or unrepresentable arithmetic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CapillarySurfaceError {
    /// A position contains NaN or infinity.
    NonFinitePosition,
    /// A triangle has zero representable area.
    DegenerateTriangle,
    /// Tension must be finite and nonnegative.
    InvalidSurfaceTension,
    /// Endpoint vertex counts or per-face coefficient counts disagree.
    LengthMismatch,
    /// A triangle names a vertex outside the supplied coordinate arrays.
    InvalidVertexIndex {
        /// Face in the supplied connectivity.
        face: usize,
        /// Invalid vertex index.
        vertex: usize,
    },
    /// An energy, force, geometry invariant or sum cannot be represented.
    Unrepresentable,
}

/// Signed energy transfers, with forces doing positive work on the liquid.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SurfaceEnergyChange {
    /// Initial geometric area [m²].
    pub area_before_m2: f64,
    /// Accepted geometric area [m²].
    pub area_after_m2: f64,
    /// Initial prescribed-coefficient free energy [J].
    pub energy_before_j: f64,
    /// Accepted prescribed-coefficient free energy [J].
    pub energy_after_j: f64,
    /// Free-energy change, evaluated without subtracting aggregate energies [J].
    pub energy_change_j: f64,
    /// External free-energy supply from prescribed coefficient changes [J].
    pub coefficient_supply_j: f64,
    /// Work of the returned force over the endpoint displacement [J].
    /// Work + energy_change - coefficient_supply vanishes to roundoff.
    pub mechanical_work_j: f64,
}

/// Allocation-free response for one triangle; local ordering matches its nodes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TriangleCapillaryResponse {
    /// Restoring forces [N] at the three vertices.
    pub forces_n: [Vector; 3],
    /// Geometric surface storage and work.
    pub energy: SurfaceEnergyChange,
}

/// A constant tension or one explicitly prescribed value per triangle.
#[derive(Clone, Copy, Debug)]
pub enum SurfaceTensionField<'a> {
    /// One coefficient [N/m] shared by every face.
    Uniform(f64),
    /// One coefficient [N/m] for each face, in connectivity order.
    PerFace(&'a [f64]),
}

impl SurfaceTensionField<'_> {
    fn validate(self, faces: usize) -> Result<(), CapillarySurfaceError> {
        match self {
            Self::Uniform(value) => tension(value),
            Self::PerFace(values) if values.len() == faces => {
                values.iter().try_for_each(|&value| tension(value))
            }
            Self::PerFace(_) => Err(CapillarySurfaceError::LengthMismatch),
        }
    }
    fn at(self, face: usize) -> f64 {
        match self {
            Self::Uniform(value) => value,
            Self::PerFace(values) => values[face],
        }
    }
}

/// Assembled patch response with one force per supplied vertex.
#[derive(Clone, Debug, PartialEq)]
pub struct CapillaryPatchResponse {
    /// Assembled restoring forces [N]; unused vertices retain zero force.
    pub forces_n: Vec<Vector>,
    /// Deterministically accumulated surface storage and work.
    pub energy: SurfaceEnergyChange,
}

#[derive(Clone, Copy)]
struct Geometry {
    a: Vector,
    b: Vector,
    normal: Vector,
    area: f64,
    gram: [f64; 3],
}

/// Instantaneous conservative force -sigma * grad(area).
///
/// Reversing triangle orientation changes no physical force. Nonmanifold or
/// open connectivity is not certified here; each submitted face contributes
/// once and the caller owns which interfaces belong to the physical surface.
pub fn triangle_surface_tension(
    vertices_m: [Vector; 3],
    surface_tension_n_m: f64,
) -> Result<TriangleCapillaryResponse, CapillarySurfaceError> {
    tension(surface_tension_n_m)?;
    let geometry = geometry(vertices_m)?;
    let ga = scaled(cross(geometry.b, geometry.normal), 0.5);
    let gb = scaled(cross(geometry.normal, geometry.a), 0.5);
    let forces = nodal_forces(ga, gb, surface_tension_n_m);
    let energy = surface_tension_n_m * geometry.area;
    response(
        forces,
        SurfaceEnergyChange {
            area_before_m2: geometry.area,
            area_after_m2: geometry.area,
            energy_before_j: energy,
            energy_after_j: energy,
            ..SurfaceEnergyChange::default()
        },
    )
}

/// Exactly work-conjugate endpoint capillary force for a finite deformation.
///
/// Let a,b be triangle edges; s=a·a, t=b·b, u=a·b, so 4A²=st-u².
/// With endpoint arithmetic means, the area discrete gradient is
/// G_a=(t_mean*a_mean-u_mean*b_mean)/(2(A_before+A_after)),
/// G_b=(s_mean*b_mean-u_mean*a_mean)/(2(A_before+A_after)).
/// The symmetric Gram coefficients make the torque about midpoint positions
/// vanish, including finite rotations. Averaging endpoint unit normals alone
/// would not have that angular-momentum property.
///
/// The force is -sigma_mean*G. The exact product identity splits the free
/// energy into geometric work and A_mean*delta_sigma. Caller-supplied changes
/// in sigma remain external free-energy supply, never an inferred heat source.
pub fn triangle_surface_tension_step(
    before_m: [Vector; 3],
    after_m: [Vector; 3],
    tension_before_n_m: f64,
    tension_after_n_m: f64,
) -> Result<TriangleCapillaryResponse, CapillarySurfaceError> {
    tension(tension_before_n_m)?;
    tension(tension_after_n_m)?;
    let before = geometry(before_m)?;
    let after = geometry(after_m)?;
    let a = mean(before.a, after.a);
    let b = mean(before.b, after.b);
    let [s, t, u] = mean(before.gram, after.gram);
    let mean_area = 0.5 * before.area + 0.5 * after.area;
    let factor = 0.25 / mean_area;
    let ga = core::array::from_fn(|i| product_difference(t, a[i], u, b[i]) * factor);
    let gb = core::array::from_fn(|i| product_difference(s, b[i], u, a[i]) * factor);
    let (ga, gb) = correct_area_work(before, after, ga, gb)?;
    let sigma = 0.5 * tension_before_n_m + 0.5 * tension_after_n_m;
    let forces = nodal_forces(ga, gb, sigma);
    let supply = mean_area * (tension_after_n_m - tension_before_n_m);
    // Remove common translation algebraically, using the zero resultant.
    let work = dot(forces[1], sub(after.a, before.a)) + dot(forces[2], sub(after.b, before.b));
    response(
        forces,
        SurfaceEnergyChange {
            area_before_m2: before.area,
            area_after_m2: after.area,
            energy_before_j: tension_before_n_m * before.area,
            energy_after_j: tension_after_n_m * after.area,
            energy_change_j: sigma * (after.area - before.area) + supply,
            coefficient_supply_j: supply,
            mechanical_work_j: work,
        },
    )
}

/// Assemble a prescribed-tension surface patch at two endpoint geometries.
///
/// Connectivity uses the same vertex-index triples as the existing mesh
/// owners. The two coordinate slices share the same vertex identities.
/// Outputs are allocated once and published only on complete success. The
/// traversal is fixed, with compensated scalar energy/area reductions.
/// Coefficients may vary by face; no smoothing or material-name default occurs.
///
/// This is a force/work operator, not an interface advection or remeshing
/// scheme. Apply its force through the receiving momentum/time owner.
pub fn assemble_surface_tension_step(
    before_m: &[Vector],
    after_m: &[Vector],
    faces: &[[usize; 3]],
    tension_before: SurfaceTensionField<'_>,
    tension_after: SurfaceTensionField<'_>,
) -> Result<CapillaryPatchResponse, CapillarySurfaceError> {
    if before_m.len() != after_m.len() {
        return Err(CapillarySurfaceError::LengthMismatch);
    }
    tension_before.validate(faces.len())?;
    tension_after.validate(faces.len())?;
    if before_m
        .iter()
        .chain(after_m)
        .flatten()
        .any(|x| !x.is_finite())
    {
        return Err(CapillarySurfaceError::NonFinitePosition);
    }
    let mut forces = vec![[0.0; 3]; before_m.len()];
    let mut totals = [Sum::default(); 7];
    for (face, &indices) in faces.iter().enumerate() {
        for vertex in indices {
            if vertex >= before_m.len() {
                return Err(CapillarySurfaceError::InvalidVertexIndex { face, vertex });
            }
        }
        let before = indices.map(|i| before_m[i]);
        let after = indices.map(|i| after_m[i]);
        let local = if before == after && tension_before.at(face) == tension_after.at(face) {
            triangle_surface_tension(before, tension_before.at(face))?
        } else {
            triangle_surface_tension_step(
                before,
                after,
                tension_before.at(face),
                tension_after.at(face),
            )?
        };
        for (node, force) in indices.into_iter().zip(local.forces_n) {
            for (total, value) in forces[node].iter_mut().zip(force) {
                *total += value;
                if !total.is_finite() {
                    return Err(CapillarySurfaceError::Unrepresentable);
                }
            }
        }
        for (total, value) in totals.iter_mut().zip(local.energy.values()) {
            total.add(value);
        }
    }
    let values = totals.map(Sum::value);
    if values.iter().any(|value| !value.is_finite()) {
        return Err(CapillarySurfaceError::Unrepresentable);
    }
    Ok(CapillaryPatchResponse {
        forces_n: forces,
        energy: SurfaceEnergyChange::from_values(values),
    })
}

impl SurfaceEnergyChange {
    fn values(self) -> [f64; 7] {
        [
            self.area_before_m2,
            self.area_after_m2,
            self.energy_before_j,
            self.energy_after_j,
            self.energy_change_j,
            self.coefficient_supply_j,
            self.mechanical_work_j,
        ]
    }
    fn from_values(v: [f64; 7]) -> Self {
        Self {
            area_before_m2: v[0],
            area_after_m2: v[1],
            energy_before_j: v[2],
            energy_after_j: v[3],
            energy_change_j: v[4],
            coefficient_supply_j: v[5],
            mechanical_work_j: v[6],
        }
    }
}

// Gram invariants can lose small terms on skinny faces. Correct only a work
// residual larger than ordinary force-work roundoff, along a symmetric
// invariant direction. Its torque about midpoint edges is identically zero.
fn correct_area_work(
    before: Geometry,
    after: Geometry,
    mut ga: Vector,
    mut gb: Vector,
) -> Result<(Vector, Vector), CapillarySurfaceError> {
    let da = sub(after.a, before.a);
    let db = sub(after.b, before.b);
    let work = dot(ga, da) + dot(gb, db);
    let residual = (after.area - before.area) - work;
    let scale = before.area
        + after.area
        + (0..3)
            .map(|i| (ga[i] * da[i]).abs() + (gb[i] * db[i]).abs())
            .sum::<f64>();
    if residual.abs() > 32.0 * f64::EPSILON * scale {
        let a = mean(before.a, after.a);
        let b = mean(before.b, after.b);
        let ds = 2.0 * dot(a, da);
        let dt = 2.0 * dot(b, db);
        let du = dot(a, db) + dot(b, da);
        let magnitude = ds.abs().max(dt.abs()).max(du.abs());
        if !magnitude.is_finite() || magnitude == 0.0 {
            return Err(CapillarySurfaceError::Unrepresentable);
        }
        let [s, t, u] = [ds / magnitude, dt / magnitude, du / magnitude];
        // C_a=s*a+u*b, C_b=t*b+u*a; C dot displacement is the positive
        // sum (ds²/2+dt²/2+du²)/magnitude, with no rotational component.
        let denominator = 0.5 * s * ds + 0.5 * t * dt + u * du;
        let correction = residual / denominator;
        if !correction.is_finite() {
            return Err(CapillarySurfaceError::Unrepresentable);
        }
        for i in 0..3 {
            ga[i] += correction * s.mul_add(a[i], u * b[i]);
            gb[i] += correction * t.mul_add(b[i], u * a[i]);
        }
    }
    Ok((ga, gb))
}

fn geometry(vertices: [Vector; 3]) -> Result<Geometry, CapillarySurfaceError> {
    if vertices.iter().flatten().any(|x| !x.is_finite()) {
        return Err(CapillarySurfaceError::NonFinitePosition);
    }
    let a = sub(vertices[1], vertices[0]);
    let b = sub(vertices[2], vertices[0]);
    let n = cross(a, b);
    let magnitude = n[0].hypot(n[1]).hypot(n[2]);
    let area = 0.5 * magnitude;
    if area == 0.0 {
        return Err(CapillarySurfaceError::DegenerateTriangle);
    }
    let gram = [dot(a, a), dot(b, b), dot(a, b)];
    if !area.is_finite() || gram.iter().any(|x| !x.is_finite()) {
        return Err(CapillarySurfaceError::Unrepresentable);
    }
    Ok(Geometry {
        a,
        b,
        normal: n.map(|x| x / magnitude),
        area,
        gram,
    })
}
fn tension(value: f64) -> Result<(), CapillarySurfaceError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(CapillarySurfaceError::InvalidSurfaceTension)
    }
}
fn response(
    forces_n: [Vector; 3],
    energy: SurfaceEnergyChange,
) -> Result<TriangleCapillaryResponse, CapillarySurfaceError> {
    if forces_n
        .iter()
        .flatten()
        .chain(energy.values().iter())
        .any(|x| !x.is_finite())
    {
        return Err(CapillarySurfaceError::Unrepresentable);
    }
    Ok(TriangleCapillaryResponse { forces_n, energy })
}
fn nodal_forces(ga: Vector, gb: Vector, sigma: f64) -> [Vector; 3] {
    let a = scaled(ga, -sigma);
    let b = scaled(gb, -sigma);
    [core::array::from_fn(|i| -(a[i] + b[i])), a, b]
}
fn sub(a: Vector, b: Vector) -> Vector {
    core::array::from_fn(|i| a[i] - b[i])
}
fn mean(a: Vector, b: Vector) -> Vector {
    core::array::from_fn(|i| 0.5 * a[i] + 0.5 * b[i])
}
fn scaled(a: Vector, s: f64) -> Vector {
    a.map(|x| x * s)
}
fn dot(a: Vector, b: Vector) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn cross(a: Vector, b: Vector) -> Vector {
    [
        product_difference(a[1], b[2], a[2], b[1]),
        product_difference(a[2], b[0], a[0], b[2]),
        product_difference(a[0], b[1], a[1], b[0]),
    ]
}
fn product_difference(a: f64, b: f64, c: f64, d: f64) -> f64 {
    let product = c * d;
    let product_error = c.mul_add(d, -product);
    a.mul_add(b, -product) - product_error
}
#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }
    fn value(self) -> f64 {
        self.sum + self.correction
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-12 * (1.0 + actual.abs() + expected.abs()),
            "actual={actual:e}, expected={expected:e}"
        );
    }
    fn vec_close(actual: Vector, expected: Vector) {
        for i in 0..3 {
            close(actual[i], expected[i]);
        }
    }
    fn sum(vectors: [Vector; 3]) -> Vector {
        core::array::from_fn(|i| vectors.iter().map(|v| v[i]).sum())
    }
    fn fixture() -> ([Vector; 3], [Vector; 3]) {
        (
            [[0.3, -0.2, 1.1], [2.0, 0.5, 1.3], [0.1, 1.6, 1.9]],
            [[0.8, -0.4, 1.0], [2.4, 1.0, 1.2], [-0.3, 1.8, 2.6]],
        )
    }
    fn balance(response: TriangleCapillaryResponse, before: [Vector; 3], after: [Vector; 3]) {
        vec_close(sum(response.forces_n), [0.0; 3]);
        let torque = core::array::from_fn(|node| {
            cross(mean(before[node], after[node]), response.forces_n[node])
        });
        vec_close(sum(torque), [0.0; 3]);
        let direct_work: f64 = (0..3)
            .map(|node| dot(response.forces_n[node], sub(after[node], before[node])))
            .sum();
        let e = response.energy;
        close(e.mechanical_work_j, direct_work);
        close(e.energy_change_j, e.energy_after_j - e.energy_before_j);
        close(
            e.mechanical_work_j + e.energy_change_j,
            e.coefficient_supply_j,
        );
    }

    #[test]
    fn g1_right_triangle_has_analytic_line_forces_and_area_energy() {
        let vertices = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 4.0, 0.0]];
        let state = triangle_surface_tension(vertices, 2.0).unwrap();
        close(state.energy.area_before_m2, 6.0);
        close(state.energy.energy_before_j, 12.0);
        for (actual, expected) in
            state
                .forces_n
                .into_iter()
                .zip([[4.0, 3.0, 0.0], [-4.0, 0.0, 0.0], [0.0, -3.0, 0.0]])
        {
            vec_close(actual, expected);
        }
        balance(state, vertices, vertices);
        let reversed =
            triangle_surface_tension([vertices[0], vertices[2], vertices[1]], 2.0).unwrap();
        vec_close(reversed.forces_n[0], state.forces_n[0]);
        vec_close(reversed.forces_n[1], state.forces_n[2]);
        vec_close(reversed.forces_n[2], state.forces_n[1]);
        let endpoint = triangle_surface_tension_step(vertices, vertices, 2.0, 2.0).unwrap();
        for (a, b) in endpoint.forces_n.into_iter().zip(state.forces_n) {
            vec_close(a, b);
        }
    }

    #[test]
    fn g1_instantaneous_force_is_negative_energy_gradient() {
        let (vertices, _) = fixture();
        let sigma = 0.072;
        let response = triangle_surface_tension(vertices, sigma).unwrap();
        let h = 1e-5;
        for node in 0..3 {
            for axis in 0..3 {
                let mut plus = vertices;
                let mut minus = vertices;
                plus[node][axis] += h;
                minus[node][axis] -= h;
                let energy = |p: [Vector; 3]| {
                    let a = sub(p[1], p[0]);
                    let b = sub(p[2], p[0]);
                    let n = cross(a, b);
                    0.5 * sigma * dot(n, n).sqrt()
                };
                let negative_gradient = -(energy(plus) - energy(minus)) / (2.0 * h);
                assert!((negative_gradient - response.forces_n[node][axis]).abs() < 2e-11);
            }
        }
    }

    #[test]
    fn g0_finite_deformation_preserves_work_and_midpoint_angular_momentum() {
        let (before, after) = fixture();
        for (sigma_before, sigma_after) in [(0.072, 0.072), (0.08, 0.03), (0.0, 0.04)] {
            let response =
                triangle_surface_tension_step(before, after, sigma_before, sigma_after).unwrap();
            balance(response, before, after);
            close(
                response.energy.coefficient_supply_j,
                0.5 * (response.energy.area_before_m2 + response.energy.area_after_m2)
                    * (sigma_after - sigma_before),
            );
            let reverse =
                triangle_surface_tension_step(after, before, sigma_after, sigma_before).unwrap();
            for (a, b) in response.forces_n.into_iter().zip(reverse.forces_n) {
                vec_close(a, b);
            }
            close(
                response.energy.mechanical_work_j,
                -reverse.energy.mechanical_work_j,
            );
            close(
                response.energy.coefficient_supply_j,
                -reverse.energy.coefficient_supply_j,
            );
        }
    }

    #[test]
    fn g3_translation_rotation_and_unit_scaling_are_objective() {
        let (before, after) = fixture();
        let rotation = |v: Vector| [v[2], v[0], v[1]];
        let translated = |v: Vector| {
            let r = rotation(v);
            [r[0] + 13.0, r[1] - 7.0, r[2] + 9.0]
        };
        let base = triangle_surface_tension_step(before, after, 0.072, 0.05).unwrap();
        let transformed = triangle_surface_tension_step(
            before.map(translated),
            after.map(translated),
            0.072,
            0.05,
        )
        .unwrap();
        for (a, b) in base.forces_n.into_iter().zip(transformed.forces_n) {
            vec_close(rotation(a), b);
        }
        close(
            base.energy.energy_change_j,
            transformed.energy.energy_change_j,
        );
        close(
            base.energy.mechanical_work_j,
            transformed.energy.mechanical_work_j,
        );
        let rigid =
            triangle_surface_tension_step(before, before.map(translated), 0.072, 0.072).unwrap();
        balance(rigid, before, before.map(translated));
        close(rigid.energy.energy_change_j, 0.0);
        close(rigid.energy.mechanical_work_j, 0.0);
        let scaled_response = triangle_surface_tension_step(
            before.map(|p| scaled(p, 3.0)),
            after.map(|p| scaled(p, 3.0)),
            0.072,
            0.05,
        )
        .unwrap();
        for (a, b) in base.forces_n.into_iter().zip(scaled_response.forces_n) {
            vec_close(scaled(a, 3.0), b);
        }
        close(
            scaled_response.energy.energy_change_j,
            9.0 * base.energy.energy_change_j,
        );
        close(
            scaled_response.energy.mechanical_work_j,
            9.0 * base.energy.mechanical_work_j,
        );
    }

    #[test]
    fn g0_prescribed_tension_change_is_external_free_energy_not_motion_work() {
        let (vertices, _) = fixture();
        let response = triangle_surface_tension_step(vertices, vertices, 0.03, 0.07).unwrap();
        assert_eq!(response.energy.mechanical_work_j, 0.0);
        close(
            response.energy.coefficient_supply_j,
            0.04 * response.energy.area_before_m2,
        );
        close(
            response.energy.energy_change_j,
            response.energy.coefficient_supply_j,
        );
        let no_tension = triangle_surface_tension_step(vertices, fixture().1, 0.0, 0.0).unwrap();
        assert_eq!(no_tension.forces_n, [[0.0; 3]; 3]);
        assert_eq!(no_tension.energy.energy_change_j, 0.0);
    }

    #[test]
    fn g1_patch_scatter_preserves_variable_tension_line_forces_and_work() {
        let before = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let faces = [[0, 1, 2], [0, 2, 3]];
        let field = SurfaceTensionField::PerFace(&[1.0, 2.0]);
        let response =
            assemble_surface_tension_step(&before, &before, &faces, field, field).unwrap();
        close(response.energy.area_before_m2, 1.0);
        close(response.energy.energy_before_j, 1.5);
        for (actual, expected) in response.forces_n.into_iter().zip([
            [0.5, 1.0, 0.0],
            [-0.5, 0.5, 0.0],
            [-1.0, -0.5, 0.0],
            [1.0, -1.0, 0.0],
        ]) {
            vec_close(actual, expected);
        }
        let after = [
            [0.1, -0.1, 0.2],
            [1.2, 0.2, 0.1],
            [1.1, 0.9, 0.5],
            [-0.2, 1.1, 0.3],
        ];
        let moved = assemble_surface_tension_step(
            &before,
            &after,
            &faces,
            field,
            SurfaceTensionField::PerFace(&[0.5, 3.0]),
        )
        .unwrap();
        let direct_work: f64 = moved
            .forces_n
            .iter()
            .enumerate()
            .map(|(i, &f)| dot(f, sub(after[i], before[i])))
            .sum();
        close(direct_work, moved.energy.mechanical_work_j);
        close(
            direct_work + moved.energy.energy_change_j,
            moved.energy.coefficient_supply_j,
        );
        let net: Vector = core::array::from_fn(|i| moved.forces_n.iter().map(|f| f[i]).sum());
        vec_close(net, [0.0; 3]);
        let torque: Vector = core::array::from_fn(|axis| {
            moved
                .forces_n
                .iter()
                .enumerate()
                .map(|(i, &f)| cross(mean(before[i], after[i]), f)[axis])
                .sum()
        });
        vec_close(torque, [0.0; 3]);
    }

    #[test]
    fn invalid_geometry_coefficients_and_topology_refuse() {
        let (before, after) = fixture();
        for sigma in [-0.01, f64::NAN, f64::INFINITY] {
            assert_eq!(
                triangle_surface_tension(before, sigma).unwrap_err(),
                CapillarySurfaceError::InvalidSurfaceTension
            );
        }
        assert_eq!(
            triangle_surface_tension([[0.0; 3]; 3], 0.072).unwrap_err(),
            CapillarySurfaceError::DegenerateTriangle
        );
        let mut invalid = before;
        invalid[1][2] = f64::NAN;
        assert_eq!(
            triangle_surface_tension(invalid, 0.072).unwrap_err(),
            CapillarySurfaceError::NonFinitePosition
        );
        let field = SurfaceTensionField::Uniform(0.072);
        assert_eq!(
            assemble_surface_tension_step(&before, &after[..2], &[[0, 1, 2]], field, field)
                .unwrap_err(),
            CapillarySurfaceError::LengthMismatch
        );
        assert_eq!(
            assemble_surface_tension_step(&before, &after, &[[0, 1, 3]], field, field).unwrap_err(),
            CapillarySurfaceError::InvalidVertexIndex { face: 0, vertex: 3 }
        );
        assert_eq!(
            assemble_surface_tension_step(
                &before,
                &after,
                &[[0, 1, 2]],
                field,
                SurfaceTensionField::PerFace(&[])
            )
            .unwrap_err(),
            CapillarySurfaceError::LengthMismatch
        );
    }

    #[test]
    fn skinny_face_work_and_nonzero_collinear_geometry_remain_physical() {
        let edge = [0.3, 0.7, 0.9];
        assert_eq!(
            triangle_surface_tension([[0.0; 3], edge, edge], 0.072).unwrap_err(),
            CapillarySurfaceError::DegenerateTriangle
        );
        for width in [1e-3, 1e-6, 1e-8, 1e-10] {
            let before = [[0.0; 3], [1.0, 0.0, 0.0], [1.0, width, 0.0]];
            let after = [[0.0; 3], [1.1, 0.1, 0.0], [1.1, 0.1 + width, 0.0]];
            let result = triangle_surface_tension_step(before, after, 1.0, 1.0).unwrap();
            let residual = result.energy.mechanical_work_j + result.energy.energy_change_j;
            assert!(
                residual.abs() < 2e-15,
                "width={width:e}, residual={residual:e}"
            );
            balance(result, before, after);
        }
    }
}
