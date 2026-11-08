//! Full meridional P1 pressure acoustics for a complete axisymmetric body.
//!
//! Coordinates are (r,z) in metres; [`AxisymmetricPressureFe::new`] keeps m=0,
//! while [`AxisymmetricPressureFe::new_fourier`] retains one signed azimuthal
//! Fourier order. The sealed
//! TriComplex2 must carry the full-turn metric 2*pi*r dA. Original nodal
//! pressure-potential and pressure fields retain radial as well as axial
//! variation; m=0 also retains constant compression. No reduced basis is used.
//!
//! With p=psi_dot and v=-grad(psi)/rho, the original physical weak equations are
//! M_beta p_dot + C_Y p + K_rho psi = b_in, where
//! M_ij=int(beta phi_i phi_j dV), K_ij=int(grad phi_i.grad phi_j/rho dV),
//! C_ij=int(Y phi_i phi_j dS), and b_i=int(phi_i v_in dS).
//! The existing weighted2 owner integrates M exactly. Constant P1 gradients
//! and the original swept cell measure integrate K exactly; cubic radial edge
//! moments integrate C and linearly varying normal motion exactly.
//!
//! Public peak phasors use exp(+i omega t), as the other pressure-FE owners.
//! Conjugation at the existing fs-modal harmonic factor boundary converts its
//! exp(-i omega t) convention; its generalized displacement/velocity are psi/p.
//! Reported powers use the original pressure/normal-velocity pair, including
//! the factor 1/2. A prepared frequency reuses its full factor for new drives.
//!
//! Axis edges have zero swept surface measure and receive no boundary data or
//! surface boundary condition. For m=0 their pressure coordinates remain live:
//! regularity is the natural finite-energy weak condition, without division by radius.
//! A finite P1 mesh does not impose a pointwise zero radial derivative in each
//! axis-adjacent cell. Nonzero Fourier orders instead constrain axis pressure
//! exactly to zero and retain their angular kinetic energy. Geometry and h
//! refinement remain necessary.
//!
//! This is a frozen, lossless, isotropic bulk fluid with piecewise constant
//! density/compressibility and real passive local wall admittance. Straight
//! meridian triangles approximate curved profiles. There is no PML, exterior
//! radiation, boundary-layer or mean-flow model here.
//! This delivers the interior Fourier slice, not completion of plan M2.4-W1d.

use crate::weighted::{WeightedAssemblyLimits, WeightedError, WeightedMass, assemble_triplets};
use crate::weighted2::{CellWeight2, weighted_mass_matrix_2d};
use fs_exec::Cx;
use fs_math::{c64::C64, dd::Dd};
use fs_modal::harmonic::{
    HarmonicResponseControls, HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_qty::Dims;
use fs_rep_mesh::{Metric2, TriComplex2, TriComplex2Error};
use fs_sparse::Csr;
use std::collections::BTreeSet;

mod fourier;
pub use fourier::{
    FourierAcousticPoint, FourierQuadrature, MAX_AZIMUTHAL_ORDER, MAX_FOURIER_QUADRATURE_POINTS,
};

/// Frozen physical coefficients for one original meridian triangle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricFluid {
    pub density_kg_m3: f64,
    /// Isentropic compressibility [1/Pa]; c² = 1/(rho beta).
    pub compressibility_pa_inv: f64,
}

/// Real passive local law v_out = Y p - v_imposed on original swept edges.
#[derive(Clone, Debug)]
pub struct AxisymmetricAdmittance {
    /// Unique physical boundary edges; vertex winding does not change the law.
    pub edges: Vec<[u32; 2]>,
    /// Specific admittance [m/(Pa s)], finite and nonnegative.
    pub admittance_m_pa_s: f64,
}

/// Original edge and inward normal velocity [m/s] at its two vertices.
/// Values follow the supplied vertex order; reversing that order with its
/// values preserves the physical drive. Velocity is linear along the edge.
pub type AxisymmetricVelocityTrace = ([u32; 2], [C64; 2]);

/// Physical swept edge, in the original complex's canonical edge order.
/// The rotation axis is excluded because it is not a physical fluid surface.
#[derive(Clone, Debug)]
pub struct AxisymmetricBoundaryEdge {
    pub global_edge: usize,
    pub vertices: [u32; 2],
    pub outward_normal_r_z: [f64; 2],
    pub surface_area_m2: f64,
    mass: [[f64; 2]; 2],
}

#[derive(Debug)]
pub enum AxisymmetricError {
    Invalid(&'static str),
    Cancelled,
    Weighted(WeightedError),
    Mesh(TriComplex2Error),
    Harmonic(HarmonicResponseError),
}
impl core::fmt::Display for AxisymmetricError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-FEEC-AXISYMMETRIC: {self:?}")
    }
}
impl std::error::Error for AxisymmetricError {}
impl From<WeightedError> for AxisymmetricError {
    fn from(error: WeightedError) -> Self {
        match error {
            WeightedError::Cancelled => Self::Cancelled,
            error => Self::Weighted(error),
        }
    }
}
impl From<TriComplex2Error> for AxisymmetricError {
    fn from(error: TriComplex2Error) -> Self {
        Self::Mesh(error)
    }
}
impl From<HarmonicResponseError> for AxisymmetricError {
    fn from(error: HarmonicResponseError) -> Self {
        match error {
            HarmonicResponseError::Cancelled => Self::Cancelled,
            error => Self::Harmonic(error),
        }
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), AxisymmetricError> {
    cx.checkpoint().map_err(|_| AxisymmetricError::Cancelled)
}
fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}

struct Triangle {
    gradients: [[f64; 2]; 3],
    orientation: f64,
}

/// Original meridian mesh, physical coefficients and sparse consistent forms.
pub struct AxisymmetricPressureFe<'m> {
    mesh: &'m TriComplex2,
    materials: Vec<AxisymmetricFluid>,
    triangles: Vec<Triangle>,
    mass: WeightedMass,
    stiffness: Csr,
    boundary: Vec<AxisymmetricBoundaryEdge>,
    boundary_index: Vec<usize>,
    volume_m3: f64,
    fourier: Option<fourier::FourierForm>,
}

/// One full original-space harmonic factor, with optional exact zero-pressure
/// edge restriction. All other vertex coordinates are retained.
pub struct PreparedAxisymmetricHarmonic<'f, 'm> {
    model: &'f AxisymmetricPressureFe<'m>,
    frequency_hz: f64,
    free_dofs: Vec<usize>,
    pressure_release: BTreeSet<usize>,
    conductance: Csr,
    factor: PreparedHarmonicResponse,
}

/// Original field and independent time-averaged physical diagnostics.
pub struct AxisymmetricHarmonicField<'f, 'm> {
    model: &'f AxisymmetricPressureFe<'m>,
    potential: Vec<C64>,
    pressure: Vec<C64>,
    pub input_power_w: f64,
    /// Original real wall loss; signed floating-point roundoff is retained.
    pub boundary_dissipation_w: f64,
    /// Total meridional plus azimuthal kinetic energy, integrated over 2*pi.
    pub kinetic_energy_j: f64,
    /// Angular part of kinetic_energy_j; exactly zero for the original m=0.
    pub azimuthal_kinetic_energy_j: f64,
    pub compressional_energy_j: f64,
    /// Native componentwise backward error; not a spatial or forward-error bound.
    pub relative_residual: f64,
}

/// A cell-side sample. Velocity is the original P1 gradient divided by that
/// cell's density, including at material interfaces; it is constant per cell.
#[derive(Clone, Copy, Debug)]
pub struct AxisymmetricAcousticPoint {
    pub position_r_z_m: [f64; 2],
    pub pressure_pa: C64,
    pub pressure_potential_pa_s: C64,
    pub particle_velocity_r_z_m_s: [C64; 2],
    /// Time-averaged active intensity Re(p conj(v))/2 [W/m²].
    pub active_intensity_r_z_w_m2: [f64; 2],
}

impl<'m> AxisymmetricPressureFe<'m> {
    /// Assemble the exact P1 physical forms on the ORIGINAL full-turn meridian.
    ///
    /// One positive finite material pair is required per stored triangle.
    /// Shared vertices retain continuous pressure across material interfaces;
    /// density and compressibility jumps are integrated before global assembly.
    /// Orphan vertices refuse rather than being silently removed.
    ///
    /// Existing weighted assembly limits bound each bulk form before allocation.
    /// Cached geometry and boundary data additionally occupy O(cells+edges+DOFs).
    /// Cancellation is polled between cells, assembly entries and boundary edges.
    pub fn new(
        mesh: &'m TriComplex2,
        materials: &[AxisymmetricFluid],
        limits: WeightedAssemblyLimits,
        cx: &Cx<'_>,
    ) -> Result<Self, AxisymmetricError> {
        poll(cx)?;
        if !matches!(mesh.metric(), Metric2::Axisymmetric { angular_span }
            if angular_span == core::f64::consts::TAU)
        {
            return Err(AxisymmetricError::Invalid(
                "full 2*pi axisymmetric meridian metric required",
            ));
        }
        let cells = mesh.faces().len();
        let n = mesh.vertices().len();
        let triplets = cells.checked_mul(9).ok_or(WeightedError::SizeOverflow)?;
        limits.check(cells, n, triplets)?;
        if materials.len() != cells {
            return Err(AxisymmetricError::Invalid(
                "one physical fluid pair per original triangle required",
            ));
        }
        let mut weights = Vec::with_capacity(cells);
        for material in materials {
            poll(cx)?;
            if !material.density_kg_m3.is_finite()
                || material.density_kg_m3 <= 0.0
                || !material.compressibility_pa_inv.is_finite()
                || material.compressibility_pa_inv <= 0.0
                || !(1.0 / material.density_kg_m3).is_finite()
            {
                return Err(AxisymmetricError::Invalid(
                    "positive finite density and compressibility required",
                ));
            }
            weights.push(CellWeight2::Scalar(material.compressibility_pa_inv));
        }
        let mass = weighted_mass_matrix_2d(
            mesh,
            0,
            &weights,
            Dims([1, -1, 2, 0, 0, 0]),
            limits,
            &mut || cx.checkpoint().is_err(),
        )?;
        let mut entries = Vec::new();
        entries
            .try_reserve_exact(triplets)
            .map_err(|_| WeightedError::Allocation)?;
        let mut triangles = Vec::with_capacity(cells);
        let mut volume = Dd::ZERO;
        for (cell, &face) in mesh.faces().iter().enumerate() {
            poll(cx)?;
            let [a, b, c] = face.map(|v| mesh.vertices()[v as usize]);
            let determinant = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
            let gradients = [
                [(b[1] - c[1]) / determinant, (c[0] - b[0]) / determinant],
                [(c[1] - a[1]) / determinant, (a[0] - c[0]) / determinant],
                [(a[1] - b[1]) / determinant, (b[0] - a[0]) / determinant],
            ];
            if !determinant.is_finite()
                || determinant == 0.0
                || gradients.iter().flatten().any(|v| !v.is_finite())
            {
                return Err(AxisymmetricError::Invalid(
                    "representable original meridian triangle gradients required",
                ));
            }
            let measure = mesh.face_measure(cell).ok_or(AxisymmetricError::Invalid(
                "missing original swept cell measure",
            ))?;
            let scale = measure / materials[cell].density_kg_m3;
            if !scale.is_finite() || scale <= 0.0 {
                return Err(AxisymmetricError::Invalid(
                    "representable original inverse-density volume measure required",
                ));
            }
            volume = volume + Dd::from_f64(measure);
            for i in 0..3 {
                for j in i..3 {
                    let value = scale
                        * (gradients[i][0] * gradients[j][0] + gradients[i][1] * gradients[j][1]);
                    if !value.is_finite() || (i == j && value <= 0.0) {
                        return Err(AxisymmetricError::Invalid(
                            "representable original meridian stiffness required",
                        ));
                    }
                    entries.push((face[i] as usize, face[j] as usize, value));
                    if i != j {
                        entries.push((face[j] as usize, face[i] as usize, value));
                    }
                }
            }
            triangles.push(Triangle {
                gradients,
                orientation: determinant.signum(),
            });
        }
        let stiffness = assemble_triplets(n, entries, &mut || cx.checkpoint().is_err())?;
        for row in 0..n {
            poll(cx)?;
            if mass.matrix().get(row, row) <= 0.0 {
                return Err(AxisymmetricError::Invalid(
                    "original meridian vertices must belong to physical cells",
                ));
            }
        }
        let trace = mesh.boundary_trace()?;
        let mut boundary = Vec::with_capacity(trace.edges.len());
        let mut boundary_index = vec![usize::MAX; mesh.edges().len()];
        for edge in trace.edges {
            poll(cx)?;
            let vertices = mesh.edges()[edge.global_edge];
            let [a, b] = vertices.map(|v| mesh.vertices()[v as usize]);
            if a[0] == 0.0 && b[0] == 0.0 {
                continue;
            }
            let area = mesh
                .edge_measure(edge.global_edge)
                .ok_or(AxisymmetricError::Invalid(
                    "missing original swept boundary measure",
                ))?;
            let length = (b[0] - a[0]).hypot(b[1] - a[1]);
            let sweep = core::f64::consts::TAU * length;
            // Exact cubic moments int_0^1 r lambda_i lambda_j ds.
            let aa = sweep * (a[0] / 4.0 + b[0] / 12.0);
            let ab = sweep * (a[0] / 12.0 + b[0] / 12.0);
            let bb = sweep * (a[0] / 12.0 + b[0] / 4.0);
            if [area, aa, ab, bb]
                .iter()
                .any(|v| !v.is_finite() || *v <= 0.0)
            {
                return Err(AxisymmetricError::Invalid(
                    "positive representable physical edge moments required",
                ));
            }
            let [from, to] = edge.oriented_vertices.map(|v| mesh.vertices()[v as usize]);
            let sign = triangles[edge.source_face].orientation;
            let normal = [
                sign * (to[1] - from[1]) / length,
                -sign * (to[0] - from[0]) / length,
            ];
            if normal.iter().any(|v| !v.is_finite()) {
                return Err(AxisymmetricError::Invalid(
                    "finite physical boundary normal required",
                ));
            }
            boundary_index[edge.global_edge] = boundary.len();
            boundary.push(AxisymmetricBoundaryEdge {
                global_edge: edge.global_edge,
                vertices,
                outward_normal_r_z: normal,
                surface_area_m2: area,
                mass: [[aa, ab], [ab, bb]],
            });
        }
        let volume_m3 = volume.to_f64();
        if !volume_m3.is_finite() || volume_m3 <= 0.0 {
            return Err(AxisymmetricError::Invalid(
                "representable positive swept volume required",
            ));
        }
        poll(cx)?;
        Ok(Self {
            mesh,
            materials: materials.to_vec(),
            triangles,
            mass,
            stiffness,
            boundary,
            boundary_index,
            volume_m3,
            fourier: None,
        })
    }

    #[must_use]
    pub const fn mesh(&self) -> &TriComplex2 {
        self.mesh
    }
    #[must_use]
    pub fn ndof(&self) -> usize {
        self.mesh.vertices().len()
    }
    #[must_use]
    pub const fn volume_m3(&self) -> f64 {
        self.volume_m3
    }
    #[must_use]
    pub fn mass(&self) -> &Csr {
        self.mass.matrix()
    }
    #[must_use]
    pub const fn stiffness(&self) -> &Csr {
        &self.stiffness
    }
    /// Original physical swept boundary edges; the zero-area axis is excluded.
    #[must_use]
    pub fn boundary_edges(&self) -> &[AxisymmetricBoundaryEdge] {
        &self.boundary
    }

    fn edge(&self, vertices: [u32; 2]) -> Result<&AxisymmetricBoundaryEdge, AxisymmetricError> {
        let index =
            self.mesh
                .edge_index(vertices[0], vertices[1])
                .ok_or(AxisymmetricError::Invalid(
                    "requested original edge is absent",
                ))?;
        let boundary = self.boundary_index[index];
        if boundary == usize::MAX {
            return Err(AxisymmetricError::Invalid(
                "physical boundary edge required; internal and axis edges carry no boundary data",
            ));
        }
        Ok(&self.boundary[boundary])
    }

    /// Prepare one positive physical frequency [Hz] in the full original space.
    ///
    /// Optional pressure-release edges impose psi=p=0 by exact coordinate
    /// restriction, before factorization. No large admittance or penalty is used.
    /// Release and admittance edges must be disjoint and unique; driven edges
    /// may overlap admittance but never release edges. A release is an ideal
    /// prescribed pressure, not an approximation to an open radiating mouth.
    ///
    /// The native factor's coordinate/residual limits apply unchanged. Its
    /// payload budget additionally reserves the three dense input matrices
    /// (24*n_free² bytes); original sparse forms, boundary caches and lifted
    /// O(ndof) fields are additional. No numerical factor is implemented here.
    pub fn prepare<'f>(
        &'f self,
        frequency_hz: f64,
        boundaries: &[AxisymmetricAdmittance],
        pressure_release_edges: &[[u32; 2]],
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<PreparedAxisymmetricHarmonic<'f, 'm>, AxisymmetricError> {
        poll(cx)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        if !frequency_hz.is_finite() || frequency_hz <= 0.0 || !omega.is_finite() {
            return Err(AxisymmetricError::Invalid(
                "positive representable acoustic frequency required",
            ));
        }
        let mut fixed = vec![false; self.ndof()];
        for &vertex in self.axis_pressure_constraints() {
            fixed[vertex] = true;
        }
        let mut pressure_release = BTreeSet::new();
        for &vertices in pressure_release_edges {
            poll(cx)?;
            let edge = self.edge(vertices)?;
            if !pressure_release.insert(edge.global_edge) {
                return Err(AxisymmetricError::Invalid(
                    "pressure-release edges must be unique",
                ));
            }
            for &vertex in &edge.vertices {
                fixed[vertex as usize] = true;
            }
        }
        let mut seen = pressure_release.clone();
        let mut entries = Vec::new();
        for boundary in boundaries {
            poll(cx)?;
            let y = boundary.admittance_m_pa_s;
            if !y.is_finite() || y < 0.0 {
                return Err(AxisymmetricError::Invalid(
                    "finite nonnegative real specific wall admittance required",
                ));
            }
            for &vertices in &boundary.edges {
                poll(cx)?;
                let edge = self.edge(vertices)?;
                if !seen.insert(edge.global_edge) {
                    return Err(AxisymmetricError::Invalid(
                        "admittance and release edges must be unique and disjoint",
                    ));
                }
                for i in 0..2 {
                    for j in 0..2 {
                        entries.push((
                            edge.vertices[i] as usize,
                            edge.vertices[j] as usize,
                            y * edge.mass[i][j],
                        ));
                    }
                }
            }
        }
        let conductance =
            assemble_triplets(self.ndof(), entries, &mut || cx.checkpoint().is_err())?;
        let free_dofs = (0..self.ndof()).filter(|&i| !fixed[i]).collect::<Vec<_>>();
        let n = free_dofs.len();
        if n == 0 {
            return Err(AxisymmetricError::Invalid(
                "pressure-release edges leave no free pressure coordinates",
            ));
        }
        if n > controls.maximum_coordinates {
            return Err(HarmonicResponseError::Budget {
                quantity: "axisymmetric free pressure coordinates",
                requested: n,
                cap: controls.maximum_coordinates,
            }
            .into());
        }
        let input_bytes =
            n.checked_mul(n)
                .and_then(|v| v.checked_mul(24))
                .ok_or(AxisymmetricError::Invalid(
                    "axisymmetric dense input extent overflow",
                ))?;
        let remaining = controls
            .maximum_payload_bytes
            .checked_sub(input_bytes)
            .ok_or(AxisymmetricError::Harmonic(HarmonicResponseError::Budget {
                quantity: "axisymmetric dense input pencil bytes",
                requested: input_bytes,
                cap: controls.maximum_payload_bytes,
            }))?;
        let mut original_to_free = vec![usize::MAX; self.ndof()];
        for (i, &global) in free_dofs.iter().enumerate() {
            original_to_free[global] = i;
        }
        let mass = restrict_dense(self.mass(), &free_dofs, &original_to_free, cx)?;
        let stiffness = restrict_dense(&self.stiffness, &free_dofs, &original_to_free, cx)?;
        let damping = restrict_dense(&conductance, &free_dofs, &original_to_free, cx)?;
        let factor = PreparedHarmonicResponse::try_new(
            &mass,
            &damping,
            &stiffness,
            n,
            omega,
            HarmonicResponseControls {
                maximum_payload_bytes: remaining,
                ..controls
            },
            cx,
        )?;
        poll(cx)?;
        Ok(PreparedAxisymmetricHarmonic {
            model: self,
            frequency_hz,
            free_dofs,
            pressure_release,
            conductance,
            factor,
        })
    }

    fn velocity_load(
        &self,
        traces: &[AxisymmetricVelocityTrace],
        pressure_release: &BTreeSet<usize>,
        cx: &Cx<'_>,
    ) -> Result<Vec<C64>, AxisymmetricError> {
        let mut load = vec![[Dd::ZERO; 2]; self.ndof()];
        let mut seen = BTreeSet::new();
        for &(vertices, velocities) in traces {
            poll(cx)?;
            if velocities.iter().any(|v| !finite(*v)) {
                return Err(AxisymmetricError::Invalid(
                    "finite original edge velocity phasors required",
                ));
            }
            let edge = self.edge(vertices)?;
            if !seen.insert(edge.global_edge) || pressure_release.contains(&edge.global_edge) {
                return Err(AxisymmetricError::Invalid(
                    "drive edges must be unique and cannot be pressure-release edges",
                ));
            }
            let values = if vertices == edge.vertices {
                velocities
            } else {
                [velocities[1], velocities[0]]
            };
            for i in 0..2 {
                let row = edge.vertices[i] as usize;
                for (j, value) in values.iter().enumerate() {
                    load[row][0] =
                        load[row][0] + Dd::from_f64(edge.mass[i][j]) * Dd::from_f64(value.re);
                    load[row][1] =
                        load[row][1] + Dd::from_f64(edge.mass[i][j]) * Dd::from_f64(value.im);
                }
            }
        }
        let result = load
            .into_iter()
            .map(|v| C64::new(v[0].to_f64(), v[1].to_f64()))
            .collect::<Vec<_>>();
        if result.iter().any(|v| !finite(*v)) {
            return Err(AxisymmetricError::Invalid(
                "original swept velocity load is unrepresentable",
            ));
        }
        poll(cx)?;
        Ok(result)
    }

    fn gradient(&self, cell: usize, potential: &[C64]) -> [C64; 2] {
        core::array::from_fn(|axis| {
            let mut value = [Dd::ZERO; 2];
            for (i, &vertex) in self.mesh.faces()[cell].iter().enumerate() {
                let gradient = Dd::from_f64(self.triangles[cell].gradients[i][axis]);
                value[0] = value[0] + gradient * Dd::from_f64(potential[vertex as usize].re);
                value[1] = value[1] + gradient * Dd::from_f64(potential[vertex as usize].im);
            }
            C64::new(value[0].to_f64(), value[1].to_f64())
        })
    }
}

impl<'f, 'm> PreparedAxisymmetricHarmonic<'f, 'm> {
    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    /// Every retained ORIGINAL vertex coordinate, in ascending global order.
    #[must_use]
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Solve another original normal-velocity trace with this prepared factor.
    /// All outputs are peak phasors in exp(+i omega t). Original full nodal
    /// vectors are returned with exactly zero entries on pressure-release traces
    /// and, for a nonzero Fourier order, on the physical axis.
    pub fn solve(
        &self,
        inward_velocity: &[AxisymmetricVelocityTrace],
        cx: &Cx<'_>,
    ) -> Result<AxisymmetricHarmonicField<'f, 'm>, AxisymmetricError> {
        poll(cx)?;
        let load = self
            .model
            .velocity_load(inward_velocity, &self.pressure_release, cx)?;
        let rhs = self
            .free_dofs
            .iter()
            .map(|&i| load[i].conj())
            .collect::<Vec<_>>();
        let solved = self.factor.solve(&rhs, cx)?;
        let mut potential = vec![C64::ZERO; self.model.ndof()];
        let mut pressure = vec![C64::ZERO; self.model.ndof()];
        for (i, &global) in self.free_dofs.iter().enumerate() {
            poll(cx)?;
            potential[global] = solved.displacement[i].conj();
            pressure[global] = solved.velocity[i].conj();
        }
        // Evaluate work from the ORIGINAL full-space pressure/velocity pair.
        // Fixed pressure reactions do no work. No endpoint-energy deficit or
        // extra circumference, density, or angular-frequency factor enters.
        let mut input = Dd::ZERO;
        for (&p, &b) in pressure.iter().zip(&load) {
            input = input
                + Dd::from_f64(p.re) * Dd::from_f64(b.re)
                + Dd::from_f64(p.im) * Dd::from_f64(b.im);
        }
        let boundary_dissipation_w = 0.5 * quadratic(&self.conductance, &pressure, cx)?;
        let compressional_energy_j = 0.25 * quadratic(self.model.mass(), &pressure, cx)?;
        // Original cell gradients avoid subtracting large constant-gauge terms
        // in kinetic storage; the same exact swept measure enters stiffness.
        let mut kinetic = Dd::ZERO;
        for (cell, material) in self.model.materials.iter().enumerate() {
            poll(cx)?;
            let scale = self
                .model
                .mesh
                .face_measure(cell)
                .ok_or(AxisymmetricError::Invalid("missing swept cell measure"))?
                / material.density_kg_m3;
            for value in self.model.gradient(cell, &potential) {
                kinetic = kinetic
                    + Dd::from_f64(scale) * Dd::from_f64(value.re) * Dd::from_f64(value.re)
                    + Dd::from_f64(scale) * Dd::from_f64(value.im) * Dd::from_f64(value.im);
            }
        }
        let input_power_w = 0.5 * input.to_f64();
        let mut kinetic_energy_j = 0.25 * kinetic.to_f64();
        let azimuthal_kinetic_energy_j = if let Some(fourier) = &self.model.fourier {
            let angular = 0.25 * quadratic(&fourier.stiffness, &potential, cx)?;
            kinetic_energy_j += angular;
            angular
        } else {
            0.0
        };
        if [
            input_power_w,
            boundary_dissipation_w,
            kinetic_energy_j,
            azimuthal_kinetic_energy_j,
            compressional_energy_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(AxisymmetricError::Invalid(
                "original acoustic energy or work is unrepresentable",
            ));
        }
        poll(cx)?;
        Ok(AxisymmetricHarmonicField {
            model: self.model,
            potential,
            pressure,
            input_power_w,
            boundary_dissipation_w,
            kinetic_energy_j,
            azimuthal_kinetic_energy_j,
            compressional_energy_j,
            relative_residual: solved.relative_residual,
        })
    }
}

impl AxisymmetricHarmonicField<'_, '_> {
    #[must_use]
    pub fn pressure_pa(&self) -> &[C64] {
        &self.pressure
    }
    #[must_use]
    pub fn pressure_potential_pa_s(&self) -> &[C64] {
        &self.potential
    }
    #[must_use]
    pub fn power_balance_residual_w(&self) -> f64 {
        self.input_power_w - self.boundary_dissipation_w
    }

    /// Interpolate the original cell-side Fourier amplitudes at barycentric
    /// coordinates. No division by r appears in these meridional components,
    /// so evaluation on the axis remains finite. For a nonzero order use
    /// [`Self::probe_fourier`] to include azimuthal velocity and the phase.
    /// The P1 radial derivative is a cell approximation, not an enforced
    /// pointwise regularity condition or a spatial error estimate.
    pub fn probe(
        &self,
        cell: usize,
        barycentric: [f64; 3],
    ) -> Result<AxisymmetricAcousticPoint, AxisymmetricError> {
        if cell >= self.model.mesh.faces().len()
            || barycentric
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0 || *v > 1.0)
            || (barycentric.iter().sum::<f64>() - 1.0).abs() > 32.0 * f64::EPSILON
        {
            return Err(AxisymmetricError::Invalid(
                "original cell and barycentric point required",
            ));
        }
        let face = self.model.mesh.faces()[cell];
        let mut position_r_z_m = [0.0; 2];
        let mut pressure_pa = C64::ZERO;
        let mut pressure_potential_pa_s = C64::ZERO;
        for (i, &vertex) in face.iter().enumerate() {
            let point = self.model.mesh.vertices()[vertex as usize];
            for axis in 0..2 {
                position_r_z_m[axis] += barycentric[i] * point[axis];
            }
            pressure_pa = pressure_pa + self.pressure[vertex as usize].scale(barycentric[i]);
            pressure_potential_pa_s =
                pressure_potential_pa_s + self.potential[vertex as usize].scale(barycentric[i]);
        }
        let inverse_density = -1.0 / self.model.materials[cell].density_kg_m3;
        let particle_velocity_r_z_m_s = self
            .model
            .gradient(cell, &self.potential)
            .map(|value| value.scale(inverse_density));
        let active_intensity_r_z_w_m2 = particle_velocity_r_z_m_s.map(|velocity| {
            0.5 * (Dd::from_f64(pressure_pa.re) * Dd::from_f64(velocity.re)
                + Dd::from_f64(pressure_pa.im) * Dd::from_f64(velocity.im))
            .to_f64()
        });
        if !finite(pressure_pa)
            || !finite(pressure_potential_pa_s)
            || particle_velocity_r_z_m_s.iter().any(|v| !finite(*v))
            || position_r_z_m
                .iter()
                .chain(&active_intensity_r_z_w_m2)
                .any(|v| !v.is_finite())
        {
            return Err(AxisymmetricError::Invalid(
                "original acoustic point is unrepresentable",
            ));
        }
        Ok(AxisymmetricAcousticPoint {
            position_r_z_m,
            pressure_pa,
            pressure_potential_pa_s,
            particle_velocity_r_z_m_s,
            active_intensity_r_z_w_m2,
        })
    }
}

fn restrict_dense(
    matrix: &Csr,
    free_dofs: &[usize],
    original_to_free: &[usize],
    cx: &Cx<'_>,
) -> Result<Vec<f64>, AxisymmetricError> {
    let n = free_dofs.len();
    let mut dense = vec![0.0; n * n];
    for (row, &global) in free_dofs.iter().enumerate() {
        poll(cx)?;
        let (columns, values) = matrix.row(global);
        for (&column, &value) in columns.iter().zip(values) {
            let col = original_to_free[column];
            if col != usize::MAX {
                dense[row * n + col] = value;
            }
        }
    }
    Ok(dense)
}

fn quadratic(matrix: &Csr, field: &[C64], cx: &Cx<'_>) -> Result<f64, AxisymmetricError> {
    let mut value = Dd::ZERO;
    let mut absolute = 0.0;
    for (row, &left) in field.iter().enumerate() {
        poll(cx)?;
        let (columns, entries) = matrix.row(row);
        for (&col, &entry) in columns.iter().zip(entries) {
            let right = field[col];
            let term = Dd::from_f64(entry)
                * (Dd::from_f64(left.re) * Dd::from_f64(right.re)
                    + Dd::from_f64(left.im) * Dd::from_f64(right.im));
            value = value + term;
            absolute += term.to_f64().abs();
        }
    }
    let value = value.to_f64();
    if !value.is_finite() || !absolute.is_finite() || value < -128.0 * f64::EPSILON * absolute {
        return Err(AxisymmetricError::Invalid(
            "negative or unrepresentable original acoustic quadratic form",
        ));
    }
    Ok(value)
}

#[cfg(test)]
mod tests;
