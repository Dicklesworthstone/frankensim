//! Resolved linearized Navier--Stokes acoustics on straight tetrahedra.
//!
//! The velocity/pressure pair is Taylor--Hood P2/P1; temperature is P2. With
//! `exp(+i omega t)` and `theta = T'/T0`, the equations are
//! `i omega rho v = div(-p I + 2 mu dev(eps(v)) + zeta div(v) I)`,
//! `div(v) + i omega (p/p0 - theta) = 0`, and
//! `i omega (rho cp T0 theta - p) - kappa T0 laplacian(theta) = 0`.
//! These are full linear equations at rest, including axial diffusion and
//! dilatational stress; no thin-layer or low-reduced-frequency approximation
//! replaces their volume operators. See Kampinga, *Viscothermal acoustics using
//! finite elements* (2010), section 2.4.
//!
//! Walls are stationary no-slip, with an explicit isothermal subset; remaining
//! walls and all pressure ports have zero normal conductive heat flux. A port
//! shape psi prescribes TOTAL traction `-P psi n`, including viscous stress.
//! Its work-conjugate inward flow is `Q = -integral(psi v.n)`. Shapes reuse
//! [`BoundaryMode`]'s real P1 traces, but they are traction shapes here, rather
//! than prescribed velocity shapes. Uniform shapes have ordinary pressure/flow
//! semantics. Nonuniform coordinates need the caller's dual mapping before
//! connecting to a prescribed-velocity modal basis.
//!
//! The assembled operator is `R + G + i omega H`: R is the viscous/thermal
//! dissipation, G is exactly antisymmetric, and H is real symmetric for a
//! calorically perfect gas. Passive complex cp adds a separately accounted
//! caloric loss. Thus boundary input work closes against resolved volume losses.
//! One cancellable original complex LU serves all simultaneous drives or all
//! columns of a multiport response. Hard limits bound this dense local-domain
//! owner; it is not a large-domain iterative solver or a geometry error bound.

use super::helmholtz::multiport::{BoundaryMode, MAX_PORT_MODES};
use super::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use crate::highorder::simplex::{SimplexSpace, duffy_quadrature};
use crate::whitney::{ElementGeometry, try_element_geometry};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};
use fs_rep_mesh::TetComplex;
use std::collections::BTreeSet;

/// Hard dense-solver cap, independent of the caller's smaller budget.
pub const MAX_LNS_DOFS: usize = 1536;
/// Hard cap on element work in this local-domain owner.
pub const MAX_LNS_ELEMENTS: usize = 8192;

/// Uniform ideal-gas equilibrium and independently specified transport in SI.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermoviscousMedium {
    /// Equilibrium density [kg/m^3].
    pub density_kg_m3: f64,
    /// Equilibrium absolute pressure [Pa].
    pub pressure_pa: f64,
    /// Equilibrium absolute temperature [K].
    pub temperature_k: f64,
    /// Equilibrium constant-pressure heat capacity [J/(kg K)].
    pub cp_j_kg_k: f64,
    /// Dynamic shear viscosity [Pa s].
    pub dynamic_viscosity_pa_s: f64,
    /// Thermal conductivity [W/(m K)].
    pub thermal_conductivity_w_m_k: f64,
    /// Independent dilatational viscosity [Pa s], nonnegative. Do not also put
    /// molecular relaxation represented by complex cp in this coefficient.
    pub bulk_viscosity_pa_s: f64,
}

impl ThermoviscousMedium {
    /// Ideal-mixture gas constant [J/(kg K)].
    #[must_use]
    pub fn gas_constant_j_kg_k(self) -> f64 {
        self.pressure_pa / self.density_kg_m3 / self.temperature_k
    }

    fn validate(self) -> Result<(), ThermoviscousError> {
        let positive = [
            self.density_kg_m3,
            self.pressure_pa,
            self.temperature_k,
            self.cp_j_kg_k,
            self.dynamic_viscosity_pa_s,
            self.thermal_conductivity_w_m_k,
        ];
        if positive.iter().any(|x| !x.is_finite() || *x <= 0.0)
            || !self.bulk_viscosity_pa_s.is_finite()
            || self.bulk_viscosity_pa_s < 0.0
            || !self.gas_constant_j_kg_k().is_finite()
            || self.cp_j_kg_k <= self.gas_constant_j_kg_k()
        {
            return Err(ThermoviscousError::InvalidInput(
                "gas equilibrium/transport",
            ));
        }
        Ok(())
    }
}

/// Boundary partition. Ports may share faces with other port modes, but never
/// with walls. Every exterior face must belong to a wall or a port.
#[derive(Debug, Clone)]
pub struct ThermoviscousBoundary {
    /// Stationary no-slip wall faces.
    pub no_slip_faces: Vec<[u32; 3]>,
    /// Isothermal faces, a subset of `no_slip_faces`; other walls are insulated.
    pub isothermal_faces: Vec<[u32; 3]>,
    /// Pressure-traction shapes. Thermal condition at each port is insulated.
    pub ports: Vec<BoundaryMode>,
}

/// Computational and geometric admission; no acoustic accuracy is inferred.
#[derive(Debug, Clone, Copy)]
pub struct ThermoviscousLimits {
    /// Maximum unconstrained scalar unknowns across all five fields.
    pub max_dofs: usize,
    /// Peak dense arrays budget [bytes], `48 * ndof^2` (two real operators,
    /// one complex system, one complex factor). Small geometry/vectors are extra.
    pub max_dense_bytes: usize,
    /// Maximum element count.
    pub max_elements: usize,
    /// Maximum `longest edge / smallest altitude` per tetrahedron.
    pub max_element_aspect_ratio: f64,
}

impl Default for ThermoviscousLimits {
    fn default() -> Self {
        Self {
            max_dofs: 1024,
            max_dense_bytes: 64 * 1024 * 1024,
            max_elements: 4096,
            max_element_aspect_ratio: 100.0,
        }
    }
}

/// Explicit admission, cancellation and numerical failures.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermoviscousError {
    /// Invalid physical input, boundary partition or response coordinates.
    InvalidInput(&'static str),
    /// Caller or hard resource cap exceeded.
    Budget {
        what: &'static str,
        requested: usize,
        limit: usize,
    },
    /// A tetrahedron exceeds the declared aspect bound.
    Aspect {
        element: usize,
        observed: f64,
        limit: f64,
    },
    /// Invalid volume geometry.
    Geometry,
    /// Cancellation before a complete result exists.
    Cancelled,
    /// Original pressure trace owner refused.
    Acoustic(AcousticFeError),
    /// Original complex factorization/substitution refused.
    Factor(LuComplexError),
    /// Nonfinite assembled coefficient or field.
    NonFinite,
    /// Extended-precision residual still exceeds the numerical solve tolerance.
    Residual { observed: f64, limit: f64 },
}

impl std::fmt::Display for ThermoviscousError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FS-FEEC-LNS: {self:?}")
    }
}
impl std::error::Error for ThermoviscousError {}
impl From<AcousticFeError> for ThermoviscousError {
    fn from(e: AcousticFeError) -> Self {
        Self::Acoustic(e)
    }
}
impl From<LuComplexError> for ThermoviscousError {
    fn from(e: LuComplexError) -> Self {
        Self::Factor(e)
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ThermoviscousError> {
    if cancelled() {
        Err(ThermoviscousError::Cancelled)
    } else {
        Ok(())
    }
}
fn canonical(mut face: [u32; 3]) -> [u32; 3] {
    face.sort_unstable();
    face
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}

/// Full mixed P2/P1/P2 spatial discretization; no modal volume reduction.
pub struct ThermoviscousFe<'c> {
    fe: HelmholtzFe<'c>,
    space: SimplexSpace<'c>,
    geometry: ElementGeometry,
    medium: ThermoviscousMedium,
    velocity_map: Vec<Option<usize>>,
    temperature_map: Vec<Option<usize>>,
    velocity_count: usize,
    pressure_start: usize,
    temperature_start: usize,
    n: usize,
    constant: Vec<f64>,
    mass: Vec<f64>,
    loads: Vec<Vec<f64>>,
    areas: Vec<f64>,
    max_aspect: f64,
}

/// One reusable factorization at a single frequency and caloric state.
pub struct PreparedThermoviscous<'a, 'c> {
    model: &'a ThermoviscousFe<'c>,
    frequency_hz: f64,
    cp: C64,
    factor: LuComplex,
    rows: Vec<f64>,
    columns: Vec<f64>,
}

/// Cycle-averaged work [W]; no clipping of loss or residual is performed.
#[derive(Debug, Clone, Copy)]
pub struct ThermoviscousPower {
    /// Work entering through prescribed total traction.
    pub input_w: f64,
    /// Deviatoric shear-viscous dissipation.
    pub viscous_w: f64,
    /// Independent dilatational-viscous dissipation.
    pub bulk_w: f64,
    /// Resolved heat-conduction dissipation.
    pub conductive_w: f64,
    /// Molecular/caloric relaxation from passive complex heat capacity.
    pub caloric_w: f64,
}
impl ThermoviscousPower {
    /// Sum of physical volume dissipation.
    #[must_use]
    pub fn dissipated_w(self) -> f64 {
        self.viscous_w + self.bulk_w + self.conductive_w + self.caloric_w
    }
    /// Input minus physical volume dissipation [W].
    #[must_use]
    pub fn balance_residual_w(self) -> f64 {
        self.input_w - self.dissipated_w()
    }
}

/// Resolved phasor fields in the original hierarchical basis.
#[derive(Debug, Clone)]
pub struct ThermoviscousField {
    /// Pressure-traction coefficients [Pa], in declared mode order.
    pub pressure_drive_pa: Vec<C64>,
    /// Work-conjugate inward flow coordinates [m^3/s].
    pub inward_flow_m3_s: Vec<C64>,
    /// P1 pressure coefficients [Pa], in vertex order.
    pub pressure_pa: Vec<C64>,
    /// P2 velocity coefficients [m/s], vertices then hierarchical edge bubbles.
    pub velocity_m_s: Vec<[C64; 3]>,
    /// P2 temperature coefficients [K], in the same hierarchical ordering.
    pub temperature_k: Vec<C64>,
    /// Independently integrated resolved losses and boundary input work.
    pub power: ThermoviscousPower,
}
impl ThermoviscousField {
    /// Mean traction work entering the domain [W].
    #[must_use]
    pub const fn input_power_w(&self) -> f64 {
        self.power.input_w
    }
}

impl<'c> ThermoviscousFe<'c> {
    /// Assemble the full mixed system once. P2/P1/P2 is fixed; an unstable
    /// equal-order pair is not an option. Curved geometry is not inferred.
    pub fn new(
        complex: &'c TetComplex,
        positions: &[[f64; 3]],
        medium: ThermoviscousMedium,
        boundary: &ThermoviscousBoundary,
        limits: ThermoviscousLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ThermoviscousError> {
        poll(cancelled)?;
        medium.validate()?;
        if limits.max_dofs == 0
            || limits.max_dense_bytes == 0
            || limits.max_elements == 0
            || !limits.max_element_aspect_ratio.is_finite()
            || limits.max_element_aspect_ratio <= 0.0
            || boundary.ports.is_empty()
            || boundary.ports.len() > MAX_PORT_MODES
            || complex.tets.is_empty()
        {
            return Err(ThermoviscousError::InvalidInput(
                "limits, mesh or port count",
            ));
        }
        let element_limit = limits.max_elements.min(MAX_LNS_ELEMENTS);
        if complex.tets.len() > element_limit {
            return Err(ThermoviscousError::Budget {
                what: "elements",
                requested: complex.tets.len(),
                limit: element_limit,
            });
        }
        let geometry =
            try_element_geometry(complex, positions).map_err(|_| ThermoviscousError::Geometry)?;
        let mut max_aspect: f64 = 0.0;
        for (t, tet) in complex.tets.iter().enumerate() {
            poll(cancelled)?;
            let mut longest: f64 = 0.0;
            for a in 0..4 {
                for b in a + 1..4 {
                    let d = std::array::from_fn::<_, 3, _>(|j| {
                        positions[tet[a] as usize][j] - positions[tet[b] as usize][j]
                    });
                    longest = longest.max(det::sqrt(d.iter().map(|v| v * v).sum()));
                }
            }
            let inverse_altitude = geometry.grads[t]
                .iter()
                .map(|g| det::sqrt(g.iter().map(|v| v * v).sum()))
                .fold(0.0, f64::max);
            let aspect = longest * inverse_altitude;
            if !aspect.is_finite() || aspect > limits.max_element_aspect_ratio {
                return Err(ThermoviscousError::Aspect {
                    element: t,
                    observed: aspect,
                    limit: limits.max_element_aspect_ratio,
                });
            }
            max_aspect = max_aspect.max(aspect);
        }
        let space = SimplexSpace::new(complex, 2);
        let gamma = medium.cp_j_kg_k / (medium.cp_j_kg_k - medium.gas_constant_j_kg_k());
        let fe = HelmholtzFe::new(
            complex,
            positions,
            2,
            Medium {
                density_kg_m3: medium.density_kg_m3,
                sound_speed_m_s: det::sqrt(gamma * medium.pressure_pa / medium.density_kg_m3),
            },
        )?;
        poll(cancelled)?;
        let mut walls = BTreeSet::new();
        let mut thermal = BTreeSet::new();
        for &face in &boundary.no_slip_faces {
            fe.boundary_geometry(face)?;
            if !walls.insert(canonical(face)) {
                return Err(ThermoviscousError::InvalidInput("duplicate wall"));
            }
        }
        for &face in &boundary.isothermal_faces {
            let key = canonical(face);
            if !walls.contains(&key) || !thermal.insert(key) {
                return Err(ThermoviscousError::InvalidInput("isothermal wall subset"));
            }
        }
        let mut covered = walls.clone();
        for mode in &boundary.ports {
            if mode.traces.is_empty() {
                return Err(ThermoviscousError::InvalidInput("empty port"));
            }
            let mut seen = BTreeSet::new();
            for &(face, values) in &mode.traces {
                fe.boundary_geometry(face)?;
                let key = canonical(face);
                if walls.contains(&key)
                    || !seen.insert(key)
                    || values.iter().any(|v| !v.is_finite())
                {
                    return Err(ThermoviscousError::InvalidInput(
                        "port/wall partition or trace",
                    ));
                }
                covered.insert(key);
            }
        }
        if fe.boundary_faces().iter().any(|f| !covered.contains(f)) {
            return Err(ThermoviscousError::InvalidInput("unassigned exterior face"));
        }
        let mask = |faces: &BTreeSet<[u32; 3]>| {
            let mut fixed = vec![false; space.ndof];
            for face in faces {
                for &v in face {
                    fixed[v as usize] = true;
                }
                for a in 0..3 {
                    for b in a + 1..3 {
                        let edge = [face[a], face[b]];
                        let id = complex
                            .edges
                            .binary_search(&edge)
                            .expect("canonical face edge");
                        fixed[space.edge_off + id] = true;
                    }
                }
            }
            let mut count = 0;
            let map = fixed
                .into_iter()
                .map(|fixed| {
                    if fixed {
                        None
                    } else {
                        let id = count;
                        count += 1;
                        Some(id)
                    }
                })
                .collect::<Vec<_>>();
            (map, count)
        };
        let (velocity_map, nv) = mask(&walls);
        let (temperature_map, nt) = mask(&thermal);
        let pressure_start = 3 * nv;
        let temperature_start = pressure_start + complex.vertex_count;
        let n = temperature_start + nt;
        let cap = limits.max_dofs.min(MAX_LNS_DOFS);
        if nv == 0 || nt == 0 {
            return Err(ThermoviscousError::InvalidInput(
                "unresolved velocity/thermal space",
            ));
        }
        if n > cap {
            return Err(ThermoviscousError::Budget {
                what: "mixed unknowns",
                requested: n,
                limit: cap,
            });
        }
        let bytes = n
            .checked_mul(n)
            .and_then(|v| v.checked_mul(48))
            .ok_or(ThermoviscousError::NonFinite)?;
        if bytes > limits.max_dense_bytes {
            return Err(ThermoviscousError::Budget {
                what: "dense bytes",
                requested: bytes,
                limit: limits.max_dense_bytes,
            });
        }
        let mut model = Self {
            fe,
            space,
            geometry,
            medium,
            velocity_map,
            temperature_map,
            velocity_count: nv,
            pressure_start,
            temperature_start,
            n,
            constant: vec![0.0; n * n],
            mass: vec![0.0; n * n],
            loads: Vec::new(),
            areas: Vec::new(),
            max_aspect,
        };
        model.assemble(cancelled)?;
        for mode in &boundary.ports {
            poll(cancelled)?;
            let mut traces = [Vec::new(), Vec::new(), Vec::new()];
            let mut area = 0.0;
            for &(face, values) in &mode.traces {
                let (_, normal, a) = model.fe.boundary_geometry(face)?;
                area += a;
                for c in 0..3 {
                    traces[c].push((face, values.map(|v| -v * normal[c])));
                }
            }
            let mut load = vec![0.0; n];
            for c in 0..3 {
                let original = model.fe.boundary_linear_load(&traces[c])?;
                for (i, &b) in original.iter().enumerate() {
                    if let Some(j) = model.velocity_map[i] {
                        load[c * nv + j] = b;
                    }
                }
            }
            if load.iter().all(|v| *v == 0.0) {
                return Err(ThermoviscousError::InvalidInput(
                    "unresolved/zero traction mode",
                ));
            }
            model.loads.push(load);
            model.areas.push(area);
        }
        poll(cancelled)?;
        Ok(model)
    }

    /// Number of unconstrained mixed unknowns.
    #[must_use]
    pub const fn ndof(&self) -> usize {
        self.n
    }
    /// Original equilibrium coefficients.
    #[must_use]
    pub const fn medium(&self) -> ThermoviscousMedium {
        self.medium
    }
    /// Number of declared traction coordinates.
    #[must_use]
    pub fn port_mode_count(&self) -> usize {
        self.loads.len()
    }
    /// Geometric area of each mode's face set [m^2]; no shape renormalization.
    #[must_use]
    pub fn port_areas_m2(&self) -> &[f64] {
        &self.areas
    }
    /// Largest admitted geometric aspect ratio.
    #[must_use]
    pub const fn max_element_aspect_ratio(&self) -> f64 {
        self.max_aspect
    }

    fn assemble(&mut self, cancelled: &mut impl FnMut() -> bool) -> Result<(), ThermoviscousError> {
        let n = self.n;
        let nv = self.velocity_count;
        let gas = self.medium;
        let mu = gas.dynamic_viscosity_pa_s;
        let lambda = gas.bulk_viscosity_pa_s - 2.0 * mu / 3.0;
        let quad = duffy_quadrature(4);
        for (t, tet) in self.space.complex.tets.iter().enumerate() {
            poll(cancelled)?;
            let dofs = self.space.element_dofs(t);
            let jac = 6.0 * self.geometry.vol_signed[t].abs();
            for &(bary, w) in &quad {
                let weight = jac * w;
                let values = dofs
                    .iter()
                    .map(|(_, f)| f.eval(bary, 2))
                    .collect::<Vec<_>>();
                let gradients = dofs
                    .iter()
                    .map(|(_, f)| {
                        let dl = f.d_lambda(bary, 2);
                        std::array::from_fn::<_, 3, _>(|c| {
                            (0..4)
                                .map(|a| dl[a] * self.geometry.grads[t][a][c])
                                .sum::<f64>()
                        })
                    })
                    .collect::<Vec<_>>();
                for (i, &(gi, _)) in dofs.iter().enumerate() {
                    for (j, &(gj, _)) in dofs.iter().enumerate() {
                        let mass = weight * values[i] * values[j];
                        let grad_dot: f64 = (0..3).map(|c| gradients[i][c] * gradients[j][c]).sum();
                        if let (Some(vi), Some(vj)) = (self.velocity_map[gi], self.velocity_map[gj])
                        {
                            for a in 0..3 {
                                for b in 0..3 {
                                    let row = a * nv + vi;
                                    let col = b * nv + vj;
                                    let strain = mu * gradients[i][b] * gradients[j][a]
                                        + lambda * gradients[i][a] * gradients[j][b]
                                        + if a == b { mu * grad_dot } else { 0.0 };
                                    self.constant[row * n + col] += weight * strain;
                                    if a == b {
                                        self.mass[row * n + col] += gas.density_kg_m3 * mass;
                                    }
                                }
                            }
                        }
                        if let (Some(ti), Some(tj)) =
                            (self.temperature_map[gi], self.temperature_map[gj])
                        {
                            let row = self.temperature_start + ti;
                            let col = self.temperature_start + tj;
                            self.constant[row * n + col] += gas.thermal_conductivity_w_m_k
                                * gas.temperature_k
                                * weight
                                * grad_dot;
                            self.mass[row * n + col] +=
                                gas.density_kg_m3 * gas.cp_j_kg_k * gas.temperature_k * mass;
                        }
                    }
                    for a in 0..4 {
                        let p = self.pressure_start + tet[a] as usize;
                        if let Some(vi) = self.velocity_map[gi] {
                            for c in 0..3 {
                                let v = c * nv + vi;
                                let divergence = weight * bary[a] * gradients[i][c];
                                self.constant[p * n + v] += divergence;
                                self.constant[v * n + p] -= divergence;
                            }
                        }
                        if let Some(ti) = self.temperature_map[gi] {
                            let temp = self.temperature_start + ti;
                            let cross = weight * bary[a] * values[i];
                            self.mass[p * n + temp] -= cross;
                            self.mass[temp * n + p] -= cross;
                        }
                    }
                }
                for a in 0..4 {
                    for b in 0..4 {
                        let row = self.pressure_start + tet[a] as usize;
                        let col = self.pressure_start + tet[b] as usize;
                        self.mass[row * n + col] += weight * bary[a] * bary[b] / gas.pressure_pa;
                    }
                }
            }
        }
        if self
            .constant
            .iter()
            .chain(&self.mass)
            .any(|x| !x.is_finite())
        {
            return Err(ThermoviscousError::NonFinite);
        }
        Ok(())
    }

    fn validate_frequency_cp(&self, frequency_hz: f64, cp: C64) -> Result<f64, ThermoviscousError> {
        let omega = std::f64::consts::TAU * frequency_hz;
        if !omega.is_finite()
            || omega <= 0.0
            || !finite(cp)
            || cp.re <= self.medium.gas_constant_j_kg_k()
            || cp.im > 0.0
        {
            return Err(ThermoviscousError::InvalidInput(
                "frequency/passive exp(+iwt) heat capacity",
            ));
        }
        Ok(omega)
    }

    fn factor(
        &self,
        frequency_hz: f64,
        cp: C64,
        reference: Option<&[f64]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(LuComplex, Vec<f64>, Vec<f64>), ThermoviscousError> {
        poll(cancelled)?;
        let omega = self.validate_frequency_cp(frequency_hz, cp)?;
        if let Some(reference) = reference {
            if reference.len() != self.loads.len()
                || reference.iter().any(|r| !r.is_finite() || *r <= 0.0)
            {
                return Err(ThermoviscousError::InvalidInput(
                    "positive real port references",
                ));
            }
        }
        let n = self.n;
        let mut a = vec![C64::ZERO; n * n];
        for row in 0..n {
            poll(cancelled)?;
            for col in 0..n {
                let i = row * n + col;
                a[i] = self.operator_entry(row, col, omega, cp);
            }
        }
        if let Some(reference) = reference {
            for (load, &r) in self.loads.iter().zip(reference) {
                for i in 0..self.pressure_start {
                    if load[i] != 0.0 {
                        poll(cancelled)?;
                        for j in 0..self.pressure_start {
                            a[i * n + j].re += r * load[i] * load[j];
                        }
                    }
                }
            }
        }
        // Two-sided algebraic equilibration changes neither fields nor work.
        // It prevents Pa, m/s and dimensionless temperature coordinates from
        // creating avoidable scale disparity in the original LU owner.
        let mut rows = vec![0.0; n];
        let mut columns = vec![0.0; n];
        for i in 0..n {
            poll(cancelled)?;
            let max = a[i * n..(i + 1) * n]
                .iter()
                .map(|z| z.abs())
                .fold(0.0, f64::max);
            if !max.is_finite() || max == 0.0 {
                return Err(ThermoviscousError::NonFinite);
            }
            rows[i] = 1.0 / max;
            if !rows[i].is_finite() {
                return Err(ThermoviscousError::NonFinite);
            }
            for j in 0..n {
                a[i * n + j] = a[i * n + j].scale(rows[i]);
            }
        }
        for j in 0..n {
            poll(cancelled)?;
            let max = (0..n).map(|i| a[i * n + j].abs()).fold(0.0, f64::max);
            if !max.is_finite() || max == 0.0 {
                return Err(ThermoviscousError::NonFinite);
            }
            columns[j] = 1.0 / max;
            if !columns[j].is_finite() {
                return Err(ThermoviscousError::NonFinite);
            }
            for i in 0..n {
                a[i * n + j] = a[i * n + j].scale(columns[j]);
            }
        }
        if a.iter().any(|z| !finite(*z)) {
            return Err(ThermoviscousError::NonFinite);
        }
        let factor = lu_complex_with_cancel(&a, n, &mut *cancelled)?;
        poll(cancelled)?;
        Ok((factor, rows, columns))
    }

    fn operator_entry(&self, row: usize, col: usize, omega: f64, cp: C64) -> C64 {
        let i = row * self.n + col;
        let mut value = C64::new(self.constant[i], omega * self.mass[i]);
        if row >= self.temperature_start && col >= self.temperature_start {
            let thermal_mass = self.mass[i] / self.medium.cp_j_kg_k;
            value.re -= omega * cp.im * thermal_mass;
            value.im = omega * cp.re * thermal_mass;
        }
        value
    }

    // Near static compression, a small acoustic flow is the difference between
    // large traction/pressure terms. Correct the ORIGINAL, unscaled equations
    // with the retained LU, using the shared double-double residual owner.
    // Refining only the rounded equilibrated matrix would retain its amplified
    // scaling error. At most two O(n^2) corrections; no extra dense storage.
    fn solve_system(
        &self,
        factor: &LuComplex,
        rows: &[f64],
        columns: &[f64],
        omega: f64,
        cp: C64,
        reference: Option<&[f64]>,
        rhs: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, ThermoviscousError> {
        poll(cancelled)?;
        let n = self.n;
        let mut x = rhs
            .iter()
            .zip(rows)
            .map(|(&b, &r)| b.scale(r))
            .collect::<Vec<_>>();
        factor.solve_with_cancel(&mut x, &mut *cancelled)?;
        for (v, &scale) in x.iter_mut().zip(columns) {
            *v = v.scale(scale);
        }
        for correction in 0..=2 {
            poll(cancelled)?;
            if x.iter().any(|z| !finite(*z)) {
                return Err(ThermoviscousError::NonFinite);
            }
            let flows = if reference.is_some() {
                self.loads
                    .iter()
                    .map(|load| {
                        let mut re = Dd::ZERO;
                        let mut im = Dd::ZERO;
                        for (&b, &v) in load.iter().zip(&x) {
                            if b != 0.0 {
                                re = re + Dd::from_f64(b) * Dd::from_f64(v.re);
                                im = im + Dd::from_f64(b) * Dd::from_f64(v.im);
                            }
                        }
                        (re, im)
                    })
                    .collect::<Vec<_>>()
            } else {
                Vec::new()
            };
            let mut residual = vec![C64::ZERO; n];
            let mut backward_error: f64 = 0.0;
            for i in 0..n {
                poll(cancelled)?;
                let mut re = Dd::from_f64(rhs[i].re);
                let mut im = Dd::from_f64(rhs[i].im);
                let mut scale = rhs[i].abs();
                for (j, &v) in x.iter().enumerate() {
                    let a = self.operator_entry(i, j, omega, cp);
                    if a.re == 0.0 && a.im == 0.0 {
                        continue;
                    }
                    re = re
                        - (Dd::from_f64(a.re) * Dd::from_f64(v.re)
                            - Dd::from_f64(a.im) * Dd::from_f64(v.im));
                    im = im
                        - (Dd::from_f64(a.re) * Dd::from_f64(v.im)
                            + Dd::from_f64(a.im) * Dd::from_f64(v.re));
                    scale += a.abs() * v.abs();
                }
                if let Some(reference) = reference {
                    for (port, &r) in reference.iter().enumerate() {
                        let b = self.loads[port][i];
                        if b != 0.0 {
                            let coefficient = Dd::from_f64(r) * Dd::from_f64(b);
                            re = re - coefficient * flows[port].0;
                            im = im - coefficient * flows[port].1;
                            scale += (r * b).abs()
                                * C64::new(flows[port].0.to_f64(), flows[port].1.to_f64()).abs();
                        }
                    }
                }
                let value = C64::new(re.to_f64(), im.to_f64());
                if !finite(value) || !scale.is_finite() {
                    return Err(ThermoviscousError::NonFinite);
                }
                backward_error = backward_error.max(value.abs() / scale.max(f64::MIN_POSITIVE));
                residual[i] = value.scale(rows[i]);
            }
            if correction == 2 || backward_error < 1e-18 {
                const LIMIT: f64 = 5e-13;
                if backward_error > LIMIT {
                    return Err(ThermoviscousError::Residual {
                        observed: backward_error,
                        limit: LIMIT,
                    });
                }
                return Ok(x);
            }
            factor.solve_with_cancel(&mut residual, &mut *cancelled)?;
            for i in 0..n {
                x[i] = x[i] + residual[i].scale(columns[i]);
            }
        }
        unreachable!("bounded refinement returns on its final residual")
    }

    /// Factor with the equilibrium real heat capacity.
    pub fn prepare<'a>(
        &'a self,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscous<'a, 'c>, ThermoviscousError> {
        self.prepare_with_heat_capacity(
            frequency_hz,
            C64::new(self.medium.cp_j_kg_k, 0.0),
            cancelled,
        )
    }

    /// Factor with source-evaluated cp [J/(kg K)] in exp(+iwt). A passive
    /// relaxing reservoir has Im(cp)<=0. Source validity/band remains the caller's.
    pub fn prepare_with_heat_capacity<'a>(
        &'a self,
        frequency_hz: f64,
        cp: C64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscous<'a, 'c>, ThermoviscousError> {
        let (factor, rows, columns) = self.factor(frequency_hz, cp, None, cancelled)?;
        Ok(PreparedThermoviscous {
            model: self,
            frequency_hz,
            cp,
            factor,
            rows,
            columns,
        })
    }

    /// Full row-major admittance Q=Y P with one factorization, in declared
    /// traction-mode coordinates. No reciprocity averaging is applied.
    pub fn admittance_with_heat_capacity(
        &self,
        frequency_hz: f64,
        cp: C64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, ThermoviscousError> {
        let prepared = self.prepare_with_heat_capacity(frequency_hz, cp, cancelled)?;
        let m = self.loads.len();
        let mut y = vec![C64::ZERO; m * m];
        for col in 0..m {
            let mut drive = vec![C64::ZERO; m];
            drive[col] = C64::ONE;
            let x = prepared.solve(&drive, cancelled)?;
            for row in 0..m {
                y[row * m + col] = self.flow(&x, row);
            }
        }
        poll(cancelled)?;
        Ok(y)
    }

    /// Direct loaded solve in positive-real-reference power coordinates:
    /// a=(P+RQ)/(2 sqrt(R)), b=(P-RQ)/(2 sqrt(R)). Avoids inverting an
    /// intermediate admittance or impedance at a closed-domain pole/zero.
    pub fn scattering_with_heat_capacity(
        &self,
        frequency_hz: f64,
        cp: C64,
        reference: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, ThermoviscousError> {
        let (factor, rows, columns) = self.factor(frequency_hz, cp, Some(reference), cancelled)?;
        let m = self.loads.len();
        let mut s = vec![C64::ZERO; m * m];
        for col in 0..m {
            poll(cancelled)?;
            let rhs = self.loads[col]
                .iter()
                .map(|&b| C64::new(2.0 * det::sqrt(reference[col]) * b, 0.0))
                .collect::<Vec<_>>();
            let solution = self.solve_system(
                &factor,
                &rows,
                &columns,
                std::f64::consts::TAU * frequency_hz,
                cp,
                Some(reference),
                &rhs,
                cancelled,
            )?;
            for row in 0..m {
                s[row * m + col] = C64::new(if row == col { 1.0 } else { 0.0 }, 0.0)
                    - self.flow(&solution, row).scale(det::sqrt(reference[row]));
            }
        }
        if s.iter().any(|z| !finite(*z)) {
            return Err(ThermoviscousError::NonFinite);
        }
        poll(cancelled)?;
        Ok(s)
    }

    fn flow(&self, x: &[C64], port: usize) -> C64 {
        x.iter()
            .zip(&self.loads[port])
            .fold(C64::ZERO, |sum, (&v, &b)| sum + v.scale(b))
    }
}

impl PreparedThermoviscous<'_, '_> {
    fn solve(
        &self,
        drive: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, ThermoviscousError> {
        poll(cancelled)?;
        if drive.len() != self.model.loads.len() || drive.iter().any(|z| !finite(*z)) {
            return Err(ThermoviscousError::InvalidInput("traction drive"));
        }
        let mut rhs = vec![C64::ZERO; self.model.n];
        for (load, &p) in self.model.loads.iter().zip(drive) {
            for i in 0..rhs.len() {
                rhs[i] = rhs[i] + p.scale(load[i]);
            }
        }
        self.model.solve_system(
            &self.factor,
            &self.rows,
            &self.columns,
            std::f64::consts::TAU * self.frequency_hz,
            self.cp,
            None,
            &rhs,
            cancelled,
        )
    }

    /// Resolve arbitrary simultaneous complex traction drives and independently
    /// integrate the physical volume losses on the same FE quadrature.
    pub fn response(
        &self,
        drive: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ThermoviscousField, ThermoviscousError> {
        let x = self.solve(drive, cancelled)?;
        let model = self.model;
        let nv = model.velocity_count;
        let pressure_pa = x[model.pressure_start..model.temperature_start].to_vec();
        let velocity_m_s = model
            .velocity_map
            .iter()
            .map(|id| match id {
                Some(i) => std::array::from_fn(|c| x[c * nv + i]),
                None => [C64::ZERO; 3],
            })
            .collect::<Vec<_>>();
        let temperature_k = model
            .temperature_map
            .iter()
            .map(|id| match id {
                Some(i) => x[model.temperature_start + i].scale(model.medium.temperature_k),
                None => C64::ZERO,
            })
            .collect::<Vec<_>>();
        let inward_flow_m3_s = (0..drive.len())
            .map(|p| model.flow(&x, p))
            .collect::<Vec<_>>();
        let mut power = ThermoviscousPower {
            input_w: drive
                .iter()
                .zip(&inward_flow_m3_s)
                .map(|(&p, &q)| 0.5 * (p.conj() * q).re)
                .sum(),
            viscous_w: 0.0,
            bulk_w: 0.0,
            conductive_w: 0.0,
            caloric_w: 0.0,
        };
        let quad = duffy_quadrature(4);
        for t in 0..model.space.complex.tets.len() {
            poll(cancelled)?;
            let dofs = model.space.element_dofs(t);
            for &(bary, w) in &quad {
                let mut grad_v = [[C64::ZERO; 3]; 3];
                let mut grad_t = [C64::ZERO; 3];
                let mut temp = C64::ZERO;
                for &(i, f) in &dofs {
                    let dl = f.d_lambda(bary, 2);
                    let gradient = std::array::from_fn::<_, 3, _>(|c| {
                        (0..4)
                            .map(|a| dl[a] * model.geometry.grads[t][a][c])
                            .sum::<f64>()
                    });
                    temp = temp + temperature_k[i].scale(f.eval(bary, 2));
                    for a in 0..3 {
                        grad_t[a] = grad_t[a] + temperature_k[i].scale(gradient[a]);
                        for b in 0..3 {
                            grad_v[a][b] = grad_v[a][b] + velocity_m_s[i][a].scale(gradient[b]);
                        }
                    }
                }
                let divergence = grad_v[0][0] + grad_v[1][1] + grad_v[2][2];
                let mut dev_norm = 0.0;
                for a in 0..3 {
                    for b in 0..3 {
                        let dev = (grad_v[a][b] + grad_v[b][a]).scale(0.5)
                            - if a == b {
                                divergence.scale(1.0 / 3.0)
                            } else {
                                C64::ZERO
                            };
                        dev_norm += dev.norm_sq();
                    }
                }
                let volume_weight = 6.0 * model.geometry.vol_signed[t].abs() * w;
                let gas = model.medium;
                power.viscous_w += volume_weight * gas.dynamic_viscosity_pa_s * dev_norm;
                power.bulk_w +=
                    0.5 * volume_weight * gas.bulk_viscosity_pa_s * divergence.norm_sq();
                power.conductive_w += 0.5 * volume_weight * gas.thermal_conductivity_w_m_k
                    / gas.temperature_k
                    * grad_t.iter().map(|z| z.norm_sq()).sum::<f64>();
                power.caloric_w += -0.5
                    * volume_weight
                    * std::f64::consts::TAU
                    * self.frequency_hz
                    * gas.density_kg_m3
                    * self.cp.im
                    / gas.temperature_k
                    * temp.norm_sq();
            }
        }
        if [
            power.input_w,
            power.viscous_w,
            power.bulk_w,
            power.conductive_w,
            power.caloric_w,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(ThermoviscousError::NonFinite);
        }
        poll(cancelled)?;
        Ok(ThermoviscousField {
            pressure_drive_pa: drive.to_vec(),
            inward_flow_m3_s,
            pressure_pa,
            velocity_m_s,
            temperature_k,
            power,
        })
    }
}

#[cfg(test)]
mod tests;
