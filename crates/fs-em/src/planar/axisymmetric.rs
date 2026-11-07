//! Axisymmetric poloidal magnetostatics with azimuthal impressed current.
//!
//! Coordinates are `[r,z]` in metres. The unknown is the regularized azimuthal
//! potential `u = A_phi/r` [T], continuous P1 on each meridional triangle:
//! `B_r = -r u_z`, `B_z = 2u + r u_r`. Consequently `A_phi = r u` vanishes
//! on the axis without setting the finite axial field to zero. No `1/r` is
//! evaluated, including in cells adjacent to the axis.
//!
//! The weak equation is `integral H.curl(r v e_phi) dV = integral J_phi r v dV`,
//! with `dV = angular_span r dr dz`. Exact cubic triangle moments integrate
//! the cylindrical curl pairing and impressed-current work. This is NOT the
//! scalar radial Laplacian used by electrostatics. The regularization follows
//! the axisymmetric vector-potential formulation in Miah et al.,
//! <https://doi.org/10.1007/s00366-025-02126-y>.
//!
//! Materials are source-resolved, linear, isotropic and lossless at fixed state.
//! A prescribed trace is `u` [T] on complete exterior edges of positive swept
//! area; every component needs such a trace. Unprescribed physical edges have
//! zero tangential H. Axis edges impose regularity only. A finite outer boundary
//! must be chosen explicitly. The force methods integrate complete-revolution
//! loads; the nonlinear owner supplies reversible saturation and recoil. No
//! infinite-domain boundary, eddy current or arbitrary angular-sector excitation
//! is inferred. Results are finite-element estimates, not continuum enclosures.

mod force;
/// Reversible saturation and bounded recoil on this same cylindrical geometry.
pub mod nonlinear;

use fs_feec::weighted::{WeightedAssemblyLimits, WeightedError, assemble_triplets};
use fs_material::{
    field_storage::linear_magnetic_permeability, state_point::ResolvedMaterialStatePoint,
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, CheckedCgError, CsrOp, SolveReport, checked_cg};
use fs_sparse::{Csr, precond::Precond};

use super::{MagneticError, PlanarMagnetostatics, finite, gradients, poll};

struct Cell {
    gradients: [[f64; 2]; 3],
    // integral N_i N_j dV; the physical radial weight makes these cubic.
    moments: [[f64; 3]; 3],
    // integral r N_i dV, multiplying the azimuthal current density.
    source_moments: [f64; 3],
}

impl Cell {
    fn new(points: [[f64; 2]; 3], angular_span: f64) -> Result<Self, MagneticError> {
        let gradients = gradients(points)?;
        let [a, b, c] = points;
        let area =
            finite(0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs())?;
        let mut moments = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in i..3 {
                let mut mean = 0.0;
                for (k, point) in points.iter().enumerate() {
                    // Exact integral N_i N_j N_k / area is 1/10,
                    // 1/30 or 1/60 according to index multiplicities.
                    let denominator = if i == j && j == k {
                        10.0
                    } else if i == j || i == k || j == k {
                        30.0
                    } else {
                        60.0
                    };
                    mean = finite(mean + point[0] / denominator)?;
                }
                let value = finite(area * angular_span * mean)?;
                moments[i][j] = value;
                moments[j][i] = value;
            }
        }
        let mut source_moments = [0.0; 3];
        for i in 0..3 {
            for (j, point) in points.iter().enumerate() {
                source_moments[i] = finite(point[0].mul_add(moments[i][j], source_moments[i]))?;
            }
        }
        Ok(Self {
            gradients,
            moments,
            source_moments,
        })
    }

    fn flux_vertices(
        &self,
        potential: [f64; 3],
        points: [[f64; 2]; 3],
    ) -> Result<[[f64; 2]; 3], MagneticError> {
        let mut gradient = [0.0; 2];
        for i in 1..3 {
            let delta = finite(potential[i] - potential[0])?;
            for (a, component) in gradient.iter_mut().enumerate() {
                *component = finite(delta.mul_add(self.gradients[i][a], *component))?;
            }
        }
        let mut b = [[0.0; 2]; 3];
        for i in 0..3 {
            let r = points[i][0];
            b[i] = [
                finite(-r * gradient[1])?,
                finite(r.mul_add(gradient[0], 2.0 * potential[i]))?,
            ];
        }
        Ok(b)
    }
}

struct Diagonal(Vec<f64>);
impl Precond for Diagonal {
    fn apply(&self, residual: &[f64], output: &mut [f64]) {
        for ((out, &value), &scale) in output.iter_mut().zip(residual).zip(&self.0) {
            *out = value * scale;
        }
    }
}

fn assembly_error(error: WeightedError) -> MagneticError {
    match error {
        WeightedError::Cancelled => MagneticError::Cancelled,
        error => MagneticError::Assembly(error),
    }
}

/// Prepared cylindrical curl operator and complete immutable material sources.
/// Repeated excitations reuse geometry, material admission and sparse assembly.
pub struct AxisymmetricMagnetostatics {
    mesh: TriComplex2,
    materials: Vec<ResolvedMaterialStatePoint>,
    permeability_h_m: Vec<f64>,
    cells: Vec<Cell>,
    stiffness: Csr,
    angular_span: f64,
}

impl AxisymmetricMagnetostatics {
    /// Assemble the correct cylindrical curl-curl pairing for `u = A_phi/r`.
    /// The supplied work limits are checked before staging cell contributions.
    pub fn new(
        mesh: TriComplex2,
        materials: Vec<ResolvedMaterialStatePoint>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, MagneticError> {
        poll(cancelled)?;
        let Metric2::Axisymmetric { angular_span } = mesh.metric() else {
            return Err(MagneticError::Invalid("A_phi requires axisymmetric metric"));
        };
        let n = mesh.vertices().len();
        let count = mesh.faces().len();
        if count == 0 || materials.len() != count {
            return Err(MagneticError::Invalid(
                "a nonempty mesh and one source material per triangle are required",
            ));
        }
        let triplets = count.checked_mul(9).ok_or(MagneticError::Invalid(
            "axisymmetric assembly extent overflow",
        ))?;
        if count > limits.max_cells || n > limits.max_dofs || triplets > limits.max_triplets {
            return Err(MagneticError::Invalid(
                "mesh exceeds declared assembly limits",
            ));
        }
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(triplets)
            .map_err(|_| MagneticError::Assembly(WeightedError::Allocation))?;
        let mut cells = Vec::with_capacity(count);
        let mut permeability_h_m = Vec::with_capacity(count);
        for (cell, face) in mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let mu =
                linear_magnetic_permeability(&materials[cell]).map_err(MagneticError::Material)?;
            let nu = mu.recip();
            if !nu.is_finite() || nu <= 0.0 {
                return Err(MagneticError::Invalid("unrepresentable reluctivity"));
            }
            let points = face.map(|v| mesh.vertices()[v as usize]);
            let geometry = Cell::new(points, angular_span)?;
            let basis: [[[f64; 2]; 3]; 3] = core::array::from_fn(|i| {
                core::array::from_fn(|a| {
                    [
                        -points[a][0] * geometry.gradients[i][1],
                        points[a][0]
                            .mul_add(geometry.gradients[i][0], if a == i { 2.0 } else { 0.0 }),
                    ]
                })
            });
            for i in 0..3 {
                for j in i..3 {
                    let value = finite(pairing(&geometry.moments, &basis[i], &basis[j])? * nu)?;
                    if i == j && value <= 0.0 {
                        return Err(MagneticError::Invalid(
                            "axisymmetric stiffness diagonal must be positive",
                        ));
                    }
                    entries.push((face[i] as usize, face[j] as usize, value));
                    if i != j {
                        entries.push((face[j] as usize, face[i] as usize, value));
                    }
                }
            }
            permeability_h_m.push(mu);
            cells.push(geometry);
        }
        let stiffness = assemble_triplets(n, entries, cancelled).map_err(assembly_error)?;
        Ok(Self {
            mesh,
            materials,
            permeability_h_m,
            cells,
            stiffness,
            angular_span,
        })
    }

    /// Original meridional geometry, topology and feature identities.
    pub fn mesh(&self) -> &TriComplex2 {
        &self.mesh
    }
    /// Complete original source records at the fixed constitutive state.
    pub fn materials(&self) -> &[ResolvedMaterialStatePoint] {
        &self.materials
    }
    /// Absolute isotropic permeability [H/m], in original triangle order.
    pub fn permeability_h_m(&self) -> &[f64] {
        &self.permeability_h_m
    }
    /// Rotationally invariant sweep angle [rad]; full revolution is 2 pi.
    pub fn angular_span_rad(&self) -> f64 {
        self.angular_span
    }

    /// Solve for cellwise `J_phi` [A/m2] and prescribed nodal `A_phi/r` [T].
    /// Azimuthal current follows positive cylindrical phi. No balanced axial
    /// return is required: each axisymmetric azimuthal current is already a
    /// closed ring. Whole current-carrying rings must lie in the modeled region.
    /// Each connected component needs a complete non-axis Dirichlet edge.
    pub fn solve(
        &self,
        current_density_a_m2: &[f64],
        boundary_potential_over_radius_t: &[(usize, f64)],
        config: CheckedCgConfig,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricMagneticField<'_>, MagneticError> {
        poll(cancelled)?;
        if current_density_a_m2.len() != self.cells.len()
            || current_density_a_m2.iter().any(|v| !v.is_finite())
            || !config.tolerance.is_finite()
            || config.tolerance <= 0.0
            || config.tolerance >= 1.0
        {
            return Err(MagneticError::Invalid(
                "invalid azimuthal current or residual tolerance",
            ));
        }
        let prescribed =
            prescribed_boundary(&self.mesh, boundary_potential_over_radius_t, cancelled)?;
        let n = prescribed.len();
        let rhs = current_load(&self.mesh, &self.cells, current_density_a_m2, cancelled)?;
        let free: Vec<_> = (0..n).filter(|&v| prescribed[v].is_none()).collect();
        let mut index = vec![usize::MAX; n];
        for (i, &v) in free.iter().enumerate() {
            index[v] = i;
        }
        let mut row_ptr = vec![0];
        let mut columns = Vec::new();
        let mut values = Vec::new();
        let mut reduced_rhs = Vec::with_capacity(free.len());
        let mut inverse_diagonal = Vec::with_capacity(free.len());
        for &v in &free {
            poll(cancelled)?;
            let mut load = rhs[v];
            let (cols, vals) = self.stiffness.row(v);
            let diagonal = self.stiffness.get(v, v);
            if diagonal <= 0.0 || !diagonal.recip().is_finite() {
                return Err(MagneticError::Invalid("invalid magnetic Jacobi diagonal"));
            }
            inverse_diagonal.push(diagonal.recip());
            for (&w, &value) in cols.iter().zip(vals) {
                if let Some(potential) = prescribed[w] {
                    load = finite((-value).mul_add(potential, load))?;
                } else {
                    columns.push(index[w]);
                    values.push(value);
                }
            }
            row_ptr.push(columns.len());
            reduced_rhs.push(load);
        }
        let matrix = Csr::from_parts(free.len(), free.len(), row_ptr, columns, values);
        let mut potential: Vec<_> = prescribed.iter().map(|v| v.unwrap_or(0.0)).collect();
        let report = if free.is_empty() {
            None
        } else {
            // u has no constant gauge: shifting it changes B_z by twice the
            // shift. Do not reuse scalar-voltage gauge subtraction here.
            let result = checked_cg(
                &CsrOp::symmetric(matrix),
                &Diagonal(inverse_diagonal),
                &reduced_rhs,
                config,
                |_| poll(cancelled),
            )
            .map_err(|error| match error {
                CheckedCgError::InvalidInput(why) => MagneticError::Invalid(why),
                CheckedCgError::Interrupted(error) => error,
            })?;
            if !result.report.converged_euclidean() {
                return Err(MagneticError::Nonconvergence(result.report));
            }
            for (&v, &value) in free.iter().zip(&result.x) {
                potential[v] = value;
            }
            Some(result.report)
        };
        self.recover(
            potential,
            current_density_a_m2.to_vec(),
            &rhs,
            &prescribed,
            report,
            cancelled,
        )
    }

    fn recover(
        &self,
        potential: Vec<f64>,
        impressed_current_density_a_m2: Vec<f64>,
        rhs: &[f64],
        prescribed: &[Option<f64>],
        solve_report: Option<SolveReport>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<AxisymmetricMagneticField<'_>, MagneticError> {
        let mut flux = Vec::with_capacity(self.cells.len());
        let mut cell_energy_j = Vec::with_capacity(self.cells.len());
        let mut energy_j = 0.0;
        for (cell, face) in self.mesh.faces().iter().enumerate() {
            poll(cancelled)?;
            let u = face.map(|v| potential[v as usize]);
            let points = face.map(|v| self.mesh.vertices()[v as usize]);
            let b = self.cells[cell].flux_vertices(u, points)?;
            let energy = finite(
                0.5 * pairing(&self.cells[cell].moments, &b, &b)? / self.permeability_h_m[cell],
            )?;
            if energy < 0.0 {
                return Err(MagneticError::Invalid("negative recovered magnetic energy"));
            }
            energy_j = finite(energy_j + energy)?;
            cell_energy_j.push(energy);
            flux.push(b);
        }
        let mut impressed_current_work_j = 0.0;
        let mut boundary_work_j = 0.0;
        for (v, &u) in potential.iter().enumerate() {
            poll(cancelled)?;
            impressed_current_work_j = finite(u.mul_add(rhs[v], impressed_current_work_j))?;
            if prescribed[v].is_some() {
                let (cols, vals) = self.stiffness.row(v);
                let mut reaction = -rhs[v];
                for (&w, &value) in cols.iter().zip(vals) {
                    reaction = finite(value.mul_add(potential[w], reaction))?;
                }
                boundary_work_j = finite(u.mul_add(reaction, boundary_work_j))?;
            }
        }
        let work_balance_error_j =
            finite(2.0 * energy_j - impressed_current_work_j - boundary_work_j)?;
        poll(cancelled)?;
        Ok(AxisymmetricMagneticField {
            problem: self,
            potential_over_radius_t: potential,
            impressed_current_density_a_m2,
            cell_vertex_flux_density_t: flux,
            cell_energy_j,
            energy_j,
            impressed_current_work_j,
            boundary_work_j,
            work_balance_error_j,
            solve_report,
        })
    }
}

fn prescribed_boundary(
    mesh: &TriComplex2,
    boundary: &[(usize, f64)],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Option<f64>>, MagneticError> {
    let prescribed = PlanarMagnetostatics::admit_mesh_boundary(mesh, boundary, cancelled)?;
    // A_phi=0 on the axis follows from regularity. Prescribing u there
    // would impose an extra axial magnetic field, not a vector-potential trace.
    let trace = mesh
        .boundary_trace()
        .map_err(|_| MagneticError::Invalid("invalid magnetic boundary trace"))?;
    let mut supported = vec![false; prescribed.len()];
    for edge in &trace.edges {
        poll(cancelled)?;
        let [a, b] = edge.oriented_vertices.map(|v| v as usize);
        if prescribed[a].is_some()
            && prescribed[b].is_some()
            && (mesh.vertices()[a][0] > 0.0 || mesh.vertices()[b][0] > 0.0)
        {
            supported[a] = true;
            supported[b] = true;
        }
    }
    if prescribed
        .iter()
        .zip(supported)
        .any(|(v, support)| v.is_some() && !support)
    {
        return Err(MagneticError::Invalid(
            "axis regularity is not a prescribed magnetic boundary",
        ));
    }
    Ok(prescribed)
}

fn current_load(
    mesh: &TriComplex2,
    cells: &[Cell],
    current_density_a_m2: &[f64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, MagneticError> {
    if current_density_a_m2.len() != cells.len()
        || current_density_a_m2.iter().any(|value| !value.is_finite())
    {
        return Err(MagneticError::Invalid(
            "complete finite azimuthal current required",
        ));
    }
    let mut rhs = vec![0.0; mesh.vertices().len()];
    for (cell, face) in mesh.faces().iter().enumerate() {
        poll(cancelled)?;
        for (i, &v) in face.iter().enumerate() {
            rhs[v as usize] = finite(
                current_density_a_m2[cell].mul_add(cells[cell].source_moments[i], rhs[v as usize]),
            )?;
        }
    }
    Ok(rhs)
}

fn pairing(
    moments: &[[f64; 3]; 3],
    a: &[[f64; 2]; 3],
    b: &[[f64; 2]; 3],
) -> Result<f64, MagneticError> {
    let mut value = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let dot = finite(a[i][0].mul_add(b[j][0], a[i][1] * b[j][1]))?;
            value = finite(moments[i][j].mul_add(dot, value))?;
        }
    }
    Ok(value)
}

/// A local field sample; cylindrical components are `[radial, axial]`.
#[derive(Debug, Clone, Copy)]
pub struct AxisymmetricMagneticSample {
    /// Azimuthal vector potential [Wb/m].
    pub azimuthal_potential_wb_m: f64,
    /// Magnetic induction [T] in the meridional plane.
    pub flux_density_t: [f64; 2],
    /// Magnetic field [A/m] in the meridional plane.
    pub magnetic_field_a_m: [f64; 2],
}

/// Computed estimate borrowing the exact mesh and complete constitutive sources.
/// B/H are affine within each cell; A_phi is quadratic. Material interfaces
/// retain their separate one-sided cell fields. No cross-interface averaging.
pub struct AxisymmetricMagneticField<'a> {
    problem: &'a AxisymmetricMagnetostatics,
    potential_over_radius_t: Vec<f64>,
    impressed_current_density_a_m2: Vec<f64>,
    cell_vertex_flux_density_t: Vec<[[f64; 2]; 3]>,
    cell_energy_j: Vec<f64>,
    energy_j: f64,
    impressed_current_work_j: f64,
    boundary_work_j: f64,
    work_balance_error_j: f64,
    solve_report: Option<SolveReport>,
}

impl<'a> AxisymmetricMagneticField<'a> {
    /// Source owner with the same lifetime as this field.
    pub fn problem(&self) -> &'a AxisymmetricMagnetostatics {
        self.problem
    }
    /// Nodal `u=A_phi/r` [T], including finite values on the axis.
    pub fn potential_over_radius_t(&self) -> &[f64] {
        &self.potential_over_radius_t
    }
    /// Original cellwise azimuthal current excitation [A/m2].
    pub fn impressed_current_density_a_m2(&self) -> &[f64] {
        &self.impressed_current_density_a_m2
    }
    /// Affine B values [T] at each cell's three stored vertices, in `[r,z]` order.
    pub fn cell_vertex_flux_density_t(&self) -> &[[[f64; 2]; 3]] {
        &self.cell_vertex_flux_density_t
    }
    /// Exact discrete integral `B.H/2` over each swept triangle [J].
    pub fn cell_energy_j(&self) -> &[f64] {
        &self.cell_energy_j
    }
    /// Total linear magnetic energy [J].
    pub fn energy_j(&self) -> f64 {
        self.energy_j
    }
    /// `integral J_phi A_phi dV` [J]; twice current-supplied stored energy
    /// under homogeneous prescribed potentials, up to the solved residual.
    pub fn impressed_current_work_j(&self) -> f64 {
        self.impressed_current_work_j
    }
    /// Prescribed-potential reaction pairing [J]. The static work identity is
    /// `2 W = current_work + boundary_work`, up to the algebraic residual.
    pub fn boundary_work_j(&self) -> f64 {
        self.boundary_work_j
    }
    /// Actual `2 W - current_work - boundary_work` [J], not an error enclosure.
    pub fn work_balance_error_j(&self) -> f64 {
        self.work_balance_error_j
    }
    /// Native checked solver report; absent when every coordinate is prescribed.
    pub fn solve_report(&self) -> Option<&SolveReport> {
        self.solve_report.as_ref()
    }

    /// Evaluate the native polynomial field within one original triangle.
    /// The local barycentric coordinates must be finite, nonnegative and sum
    /// to one within 32 machine epsilons. Axis samples require no special case.
    pub fn sample(
        &self,
        cell: usize,
        barycentric: [f64; 3],
    ) -> Result<AxisymmetricMagneticSample, MagneticError> {
        let (potential, b) = sample_components(
            &self.problem.mesh,
            &self.potential_over_radius_t,
            &self.cell_vertex_flux_density_t,
            cell,
            barycentric,
        )?;
        Ok(AxisymmetricMagneticSample {
            azimuthal_potential_wb_m: potential,
            flux_density_t: b,
            magnetic_field_a_m: [
                finite(b[0] / self.problem.permeability_h_m[cell])?,
                finite(b[1] / self.problem.permeability_h_m[cell])?,
            ],
        })
    }
}

fn sample_components(
    mesh: &TriComplex2,
    potential: &[f64],
    flux: &[[[f64; 2]; 3]],
    cell: usize,
    barycentric: [f64; 3],
) -> Result<(f64, [f64; 2]), MagneticError> {
    if cell >= mesh.faces().len()
        || barycentric.iter().any(|v| !v.is_finite() || *v < 0.0)
        || (barycentric.iter().sum::<f64>() - 1.0).abs() > 32.0 * f64::EPSILON
    {
        return Err(MagneticError::Invalid("invalid magnetic cell sample"));
    }
    let mut radius = 0.0;
    let mut u = 0.0;
    let mut b = [0.0; 2];
    for (i, &weight) in barycentric.iter().enumerate() {
        let v = mesh.faces()[cell][i] as usize;
        radius = finite(weight.mul_add(mesh.vertices()[v][0], radius))?;
        u = finite(weight.mul_add(potential[v], u))?;
        for (a, component) in b.iter_mut().enumerate() {
            *component = finite(weight.mul_add(flux[cell][i][a], *component))?;
        }
    }
    Ok((finite(radius * u)?, b))
}
