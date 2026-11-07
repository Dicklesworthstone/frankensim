//! Full original-coordinate tensor-Prony shell response with reciprocal radiation.
//!
//! The immutable source shell supplies its geometry, supports, material law,
//! inertia and physical coordinate map. Its oriented midsurface is the closed
//! exterior acoustic boundary. Uniform planar refinement retains that surface;
//! no thickness offset, cavity, modal basis or missing face is inferred.
//!
//! Every original free translation and physical rotation remains in the loaded
//! solve. The established P1 midsurface translation trace supplies area-mean P0
//! normal velocity on each acoustic child panel. Rotations have no direct
//! normal velocity at this zero-offset trace; their original mechanical coupling
//! and inertia remain. Acoustic reaction is exactly the transpose work map.
//!
//! The original reciprocal P1/P0 owner supplies the full complex self/mutual
//! impedance, coherent field and receiver arrays. The original shell owner
//! supplies frequency-dependent material reduction, one loaded factor, strain
//! admission and intrinsic heat. All phasors are peak values under exp(-i w t).
//! Pointwise forced response establishes no causal realization, loaded-pole
//! certificate, open-screen solution or continuum/discretization error bound.

use fs_bem::helmholtz::{HelmholtzError, MAX_RADIATION_FIELDS_PER_BATCH, Medium};
use fs_bem::near_field::Options as ReceiverOptions;
use fs_bem::panel3d::SpherePanels;
use fs_bem::reciprocal_radiation::{
    Budget as BoundaryBudget, CoherentPortFieldLimits, CoherentRadiationPower,
    CoherentVelocityPortField, MAX_RECIPROCAL_PANELS, Prepared, PreparedReceivers,
    ReceiverGeometry, VelocityPortSolution,
};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicPortControls, HarmonicResponse, HarmonicResponseControls,
};
use fs_plate::shell::viscosity::harmonic::{
    PreparedTensorPronyShell, TensorPronyShell, TensorShellHarmonicError,
};

/// Explicit source-to-boundary geometry and physical-map resources.
#[derive(Clone, Copy, Debug)]
pub struct TensorShellRadiationControls {
    /// Uniform planar quartering levels. Each original facet has 4^levels children.
    pub acoustic_refinements: u32,
    /// Includes clamped faces, which still scatter sound.
    pub maximum_panels: usize,
    /// All original free coordinates; also bounded by the native source-batch cap.
    pub maximum_coordinates: usize,
    /// Geometry, refinement scratch, identity ports and complete velocity rows.
    /// Bound: 8*q*P + 8*q*q + 224*P + 24*q + 2048*(levels+1) + 256 bytes.
    /// P is acoustic panels and q original free coordinates; allocator overhead
    /// and the borrowed shell are excluded. Native static BEM cache is separate.
    pub maximum_setup_bytes: usize,
    /// Conservative geometry/map score: 256*P + 16*q*P + q*q.
    pub maximum_setup_work: u64,
}

/// Additional admission for one frequency; simultaneous preparations add their
/// storage. The native static geometry cache and receiver rows are separate.
#[derive(Clone, Copy, Debug)]
pub struct TensorShellRadiationFrequencyControls {
    /// Original complete shell material, matrix and factor resources.
    pub mechanical: HarmonicResponseControls,
    /// Original full physical-port projection and pointwise PSD admission.
    pub ports: HarmonicPortControls,
    /// Original source columns and native frequency assembly/LU/angular scratch.
    /// Conservative numeric bound: 64*N^2 + 128*q*(N+P) + 32*q^2
    /// + 64*A*N + 256*N + 65536, where A=64*4^angular_refinement.
    /// Mechanical/port payload, coherent results and allocator overhead are separate.
    pub maximum_boundary_payload_bytes: usize,
    /// One coherent result and its streaming outgoing-flux scratch.
    pub field: CoherentPortFieldLimits,
    /// Work for each coherent outgoing-flux assessment, not a cumulative sweep.
    pub maximum_radiated_power_work: u64,
    /// All three independent work defects must lie below this finite (0,1) value.
    pub power_relative_tolerance: f64,
}

/// Native refusals retain their numerical or material source.
#[derive(Debug)]
pub enum TensorShellRadiationError {
    Invalid(&'static str),
    Budget { quantity: &'static str, required: u128, maximum: u128 },
    Geometry(fs_bem::BemError),
    Boundary(HelmholtzError),
    Shell(TensorShellHarmonicError),
    /// Total input/heat/flux, projected port/surface work, and surface/flux defects.
    PowerBalance { coupled: f64, port_work: f64, flux: f64, tolerance: f64 },
    Allocation(&'static str),
    Cancelled,
}
impl core::fmt::Display for TensorShellRadiationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-COUPLE-TENSOR-SHELL-RADIATION: {self:?}")
    }
}
impl std::error::Error for TensorShellRadiationError {}
impl From<fs_bem::BemError> for TensorShellRadiationError {
    fn from(value: fs_bem::BemError) -> Self { Self::Geometry(value) }
}
impl From<HelmholtzError> for TensorShellRadiationError {
    fn from(value: HelmholtzError) -> Self { Self::Boundary(value) }
}
impl From<TensorShellHarmonicError> for TensorShellRadiationError {
    fn from(value: TensorShellHarmonicError) -> Self { Self::Shell(value) }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), TensorShellRadiationError> {
    cx.checkpoint().map_err(|_| TensorShellRadiationError::Cancelled)
}
fn admit(quantity: &'static str, required: u128, maximum: u128)
    -> Result<(), TensorShellRadiationError>
{
    if required > maximum {
        Err(TensorShellRadiationError::Budget { quantity, required, maximum })
    } else {
        Ok(())
    }
}
fn reserved<T>(length: usize, what: &'static str) -> Result<Vec<T>, TensorShellRadiationError> {
    let mut values = Vec::new();
    values.try_reserve_exact(length).map_err(|_| TensorShellRadiationError::Allocation(what))?;
    Ok(values)
}

/// Immutable exact shell source, homogeneous fluid and physical midsurface map.
pub struct TensorShellRadiation<'shell> {
    shell: &'shell TensorPronyShell,
    surface: SpherePanels,
    velocity_rows: Vec<Vec<f64>>,
    identity_ports: Vec<f64>,
    medium: Medium,
}
impl<'shell> TensorShellRadiation<'shell> {
    /// Derive every panel and motion row from the source shell itself.
    /// All facets remain, including those with zero free motion. The native
    /// closed-exterior topology/orientation and source band are admitted by
    /// `prepare`; this constructor never caps, flips or welds a refused surface.
    ///
    /// # Errors
    /// Nonfinite fluid, incomplete resources, native panel geometry, allocation
    /// failure or cooperative cancellation. No coordinate is dropped to fit.
    pub fn new(
        shell: &'shell TensorPronyShell,
        medium: Medium,
        controls: TensorShellRadiationControls,
        cx: &Cx<'_>,
    ) -> Result<Self, TensorShellRadiationError> {
        checkpoint(cx)?;
        if !medium.density.is_finite() || medium.density <= 0.0
            || !medium.sound_speed.is_finite() || medium.sound_speed <= 0.0
        {
            return Err(TensorShellRadiationError::Invalid("positive finite homogeneous medium required"));
        }
        let q = shell.model().free;
        if q == 0 || controls.acoustic_refinements > 6 {
            return Err(TensorShellRadiationError::Invalid("nonempty original coordinates and at most six refinement levels required"));
        }
        let panels = shell.mesh().tris.len().checked_mul(4usize.pow(controls.acoustic_refinements))
            .ok_or(TensorShellRadiationError::Invalid("refined panel count overflow"))?;
        admit("original shell coordinates", q as u128,
            controls.maximum_coordinates.min(MAX_RADIATION_FIELDS_PER_BATCH) as u128)?;
        admit("closed radiation panels", panels as u128,
            controls.maximum_panels.min(MAX_RECIPROCAL_PANELS) as u128)?;
        let (n, p) = (q as u128, panels as u128);
        let bytes = 8*n*p + 8*n*n + 224*p + 24*n
            + 2048*(controls.acoustic_refinements as u128 + 1) + 256;
        let work = 256*p + 16*n*p + n*n;
        admit("shell radiation setup bytes", bytes, controls.maximum_setup_bytes as u128)?;
        admit("shell radiation setup work", work, controls.maximum_setup_work as u128)?;
        let mut triangles = reserved(panels, "refined original shell triangles")?;
        let mut traces = reserved(panels, "original shell P1 trace coordinates")?;
        for (element, nodes) in shell.mesh().tris.iter().enumerate() {
            checkpoint(cx)?;
            append_refined(
                nodes.map(|node| shell.mesh().nodes[node]),
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                controls.acoustic_refinements, element, &mut triangles, &mut traces, cx,
            )?;
        }
        let surface = SpherePanels::from_triangles(triangles)?;
        let mut velocity_rows = reserved(q, "complete shell radiation velocity rows")?;
        for _ in 0..q {
            checkpoint(cx)?;
            let mut row = reserved(panels, "shell radiation velocity row")?;
            row.resize(panels, 0.0);
            velocity_rows.push(row);
        }
        for (panel, ((element, weights), normal)) in traces.iter().zip(surface.normals()).enumerate() {
            checkpoint(cx)?;
            for (corner, &node) in shell.mesh().tris[*element].iter().enumerate() {
                for (component, &direction) in normal.iter().enumerate() {
                    if let Some(free) = shell.model().dof_map[6*node + component] {
                        velocity_rows[free][panel] += weights[corner] * direction;
                    }
                }
            }
        }
        let mut identity_ports = reserved(q*q, "identity physical shell ports")?;
        identity_ports.resize(q*q, 0.0);
        for i in 0..q { identity_ports[i*q+i] = 1.0; }
        checkpoint(cx)?;
        Ok(Self { shell, surface, velocity_rows, identity_ports, medium })
    }

    /// Same material/geometry/support source; no mutable access is exposed.
    pub const fn shell(&self) -> &'shell TensorPronyShell { self.shell }
    /// Exact original planar surface with optional conforming planar subdivisions.
    pub const fn surface(&self) -> &SpherePanels { &self.surface }
    pub const fn medium(&self) -> Medium { self.medium }
    /// Row j converts original generalized velocity j to outward panel velocity.
    /// Translation velocities are m/s, rotations rad/s; zero rotational rows
    /// are retained. The identity generalized effort chart is force/torque.
    pub fn normal_velocity_rows(&self) -> &[Vec<f64>] { &self.velocity_rows }

    /// Admit the entire requested exterior band before frequency source work.
    /// Inherits the original reciprocal owner's closed-surface, wavelength,
    /// angular-resolution and sub-resonance restrictions unchanged.
    ///
    /// # Errors
    /// Native closed geometry/band/cache/work refusal or cancellation.
    pub fn prepare(
        &self,
        highest_omega_rad_s: f64,
        budget: BoundaryBudget,
        cx: &Cx<'_>,
    ) -> Result<PreparedTensorShellRadiation<'_, 'shell>, TensorShellRadiationError> {
        checkpoint(cx)?;
        if !highest_omega_rad_s.is_finite() || highest_omega_rad_s <= 0.0 {
            return Err(TensorShellRadiationError::Invalid("positive finite acoustic source band required"));
        }
        let boundary = Prepared::new(&self.surface, highest_omega_rad_s/self.medium.sound_speed,
            budget, || cx.checkpoint().is_err())?;
        admit("complete original-coordinate BEM batch work",
            boundary.velocity_port_work(self.velocity_rows.len())? as u128,
            budget.maximum_work as u128)?;
        checkpoint(cx)?;
        Ok(PreparedTensorShellRadiation {
            source: self, boundary, highest_omega_rad_s, budget,
        })
    }
}

fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    // Identical shared endpoints produce identical children in either winding.
    core::array::from_fn(|i| 0.5*a[i] + 0.5*b[i])
}
fn append_refined(
    points: [[f64; 3]; 3], barycentric: [[f64; 3]; 3], levels: u32, element: usize,
    triangles: &mut Vec<[[f64; 3]; 3]>, traces: &mut Vec<(usize, [f64; 3])>, cx: &Cx<'_>,
) -> Result<(), TensorShellRadiationError> {
    checkpoint(cx)?;
    if levels == 0 {
        triangles.push(points);
        traces.push((element, core::array::from_fn(|i|
            (barycentric[0][i] + barycentric[1][i] + barycentric[2][i])/3.0)));
        return Ok(());
    }
    let p = [points[0], points[1], points[2],
        midpoint(points[0], points[1]), midpoint(points[1], points[2]),
        midpoint(points[2], points[0])];
    let b = [barycentric[0], barycentric[1], barycentric[2],
        midpoint(barycentric[0], barycentric[1]), midpoint(barycentric[1], barycentric[2]),
        midpoint(barycentric[2], barycentric[0])];
    for indices in [[0, 3, 5], [3, 1, 4], [5, 4, 2], [3, 4, 5]] {
        append_refined(indices.map(|i| p[i]), indices.map(|i| b[i]), levels-1,
            element, triangles, traces, cx)?;
    }
    Ok(())
}

/// One source-bound native geometry/cache, reusable across frequencies and forces.
pub struct PreparedTensorShellRadiation<'geometry, 'shell> {
    source: &'geometry TensorShellRadiation<'shell>,
    boundary: Prepared<'geometry>,
    highest_omega_rad_s: f64,
    budget: BoundaryBudget,
}
impl<'geometry, 'shell> PreparedTensorShellRadiation<'geometry, 'shell> {
    pub const fn source(&self) -> &'geometry TensorShellRadiation<'shell> { self.source }
    /// Native immutable boundary, for original diagnostics and independent checks.
    pub const fn boundary(&self) -> &Prepared<'geometry> { &self.boundary }

    /// Admit an entire 1..=64 exterior receiver set before source solves.
    /// Clearances are actual shortest distances to the supplied triangles.
    pub fn receiver_geometry(
        &self, points: &[[f64; 3]], minimum_clearance_m: f64, cx: &Cx<'_>,
    ) -> Result<ReceiverGeometry<'_, 'geometry>, TensorShellRadiationError> {
        checkpoint(cx)?;
        let geometry = self.boundary.receivers(points, minimum_clearance_m)?;
        checkpoint(cx)?;
        Ok(geometry)
    }

    /// Solve every original unit-velocity acoustic column with one native BEM
    /// factor, then load the full tensor shell through its existing factor owner.
    /// The real generalized force map is identity in the ORIGINAL free coordinate
    /// chart. Panel areas already occur in BEM impedance; none are applied twice.
    /// Dark and rotational columns remain; full self/mutual/reactive entries are
    /// neither clipped nor symmetrized. No unloaded mechanical solve is required.
    ///
    /// # Errors
    /// Out-of-band request, per-frequency budgets, native acoustic/material/
    /// passivity/factor refusal, allocation failure or cancellation.
    pub fn prepare_frequency(
        &self, omega_rad_s: f64, controls: TensorShellRadiationFrequencyControls, cx: &Cx<'_>,
    ) -> Result<PreparedTensorShellRadiationFrequency<'_, 'geometry, 'shell>, TensorShellRadiationError> {
        checkpoint(cx)?;
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 || omega_rad_s > self.highest_omega_rad_s {
            return Err(TensorShellRadiationError::Invalid("frequency lies outside the admitted exterior band"));
        }
        if !controls.power_relative_tolerance.is_finite()
            || controls.power_relative_tolerance <= 0.0 || controls.power_relative_tolerance >= 1.0
        {
            return Err(TensorShellRadiationError::Invalid("finite power tolerance in (0,1) required"));
        }
        let q = self.source.velocity_rows.len();
        let (n, p, ports) = (self.boundary.unknowns() as u128,
            self.source.surface.areas().len() as u128, q as u128);
        let angular = 64_u128 * (1_u128 << (2*self.budget.angular_refinement));
        admit("original harmonic coordinates", ports, controls.mechanical.maximum_coordinates as u128)?;
        admit("original harmonic ports", ports, controls.ports.maximum_ports as u128)?;
        let boundary_bytes = 64*n*n + 128*ports*(n+p) + 32*ports*ports
            + 64*angular*n + 256*n + 65536;
        admit("original BEM frequency payload", boundary_bytes,
            controls.maximum_boundary_payload_bytes as u128)?;
        admit("coherent field payload", 32*n + 32*p + 32*angular + 256,
            controls.field.maximum_payload_bytes as u128)?;
        admit("coherent field work", 64*ports*(n+p) + 128*p,
            controls.field.maximum_work as u128)?;
        admit("coherent outgoing-flux work", 128*angular*(n+16*p),
            controls.maximum_radiated_power_work as u128)?;
        let mut rows = reserved(q, "original shell velocity row references")?;
        rows.extend(self.source.velocity_rows.iter().map(Vec::as_slice));
        let wavenumber = omega_rad_s/self.source.medium.sound_speed;
        let radiation = self.boundary.solve_velocity_ports(wavenumber, self.source.medium,
            &rows, || cx.checkpoint().is_err())?;
        let shell: &'shell TensorPronyShell = self.source.shell;
        let mechanical = shell.prepare_harmonic_response_with_port_impedance(
            omega_rad_s, &self.source.identity_ports, q, radiation.impedance(),
            controls.mechanical, controls.ports, cx,
        )?;
        checkpoint(cx)?;
        Ok(PreparedTensorShellRadiationFrequency {
            prepared: self, radiation, mechanical, controls, wavenumber,
        })
    }
}

/// Reusable loaded original-coordinate factor and SAME bound acoustic sources.
pub struct PreparedTensorShellRadiationFrequency<'p, 'geometry, 'shell> {
    prepared: &'p PreparedTensorShellRadiation<'geometry, 'shell>,
    radiation: VelocityPortSolution<'p, 'geometry>,
    mechanical: PreparedTensorPronyShell<'shell>,
    controls: TensorShellRadiationFrequencyControls,
    wavenumber: f64,
}
impl<'p, 'geometry, 'shell> PreparedTensorShellRadiationFrequency<'p, 'geometry, 'shell> {
    /// Allows native material stress/heat observation of the returned displacement.
    pub const fn shell_preparation(&self) -> &PreparedTensorPronyShell<'shell> { &self.mechanical }
    /// Full physical row-major impedance, including dark original coordinates.
    pub fn impedance(&self) -> &[C64] { self.radiation.impedance() }
    pub const fn wavenumber_per_m(&self) -> f64 { self.wavenumber }

    /// Prepare the original shared pressure/particle-velocity rows at this exact
    /// source frequency and fluid. Observer geometry never changes a termination.
    /// The native quadrature budget is aggregate over the entire receiver array.
    pub fn prepare_receivers(
        &self, points: &[[f64; 3]], minimum_clearance_m: f64, options: ReceiverOptions,
        particle_velocity: bool, cx: &Cx<'_>,
    ) -> Result<PreparedReceivers<'p, 'geometry>, TensorShellRadiationError> {
        checkpoint(cx)?;
        let boundary: &'p Prepared<'geometry> = &self.prepared.boundary;
        Ok(boundary.receivers(points, minimum_clearance_m)?.prepare(
            self.wavenumber, self.prepared.source.medium, options, particle_velocity,
            || cx.checkpoint().is_err(),
        )?)
    }

    /// Solve original force/torque phasors and publish one coherent motion/field.
    /// Intrinsic tensor-Prony/Kelvin heat stays separate from boundary work and
    /// integrated outgoing flux. Source-cycle strain and original physical-row
    /// residual admission run in the native shell solver before observation.
    ///
    /// # Errors
    /// Native force/material/strain refusal, coherent field/work overflow or
    /// budget, any failed work identity, or cancellation. No partial result.
    pub fn solve(
        &self, force: &[C64], cx: &Cx<'_>,
    ) -> Result<TensorShellRadiationResponse<'p, 'geometry>, TensorShellRadiationError> {
        checkpoint(cx)?;
        let motion = self.mechanical.solve(force, cx)?;
        let field = self.radiation.coherent(&motion.velocity, self.controls.field,
            || cx.checkpoint().is_err())?;
        let power = field.power(self.controls.maximum_radiated_power_work,
            || cx.checkpoint().is_err())?;
        let port = motion.port_load.as_ref().ok_or(TensorShellRadiationError::Invalid(
            "native loaded shell omitted its physical port work"))?;
        let port_work = relative_balance(port.average_absorbed_power_w, power.surface_power_w.re, 0.0);
        let coupled = relative_balance(motion.average_input_power_w,
            motion.average_viscous_power_w, power.outgoing_power_w);
        let tolerance = self.controls.power_relative_tolerance;
        if !port_work.is_finite() || !coupled.is_finite()
            || port_work > tolerance || coupled > tolerance || power.relative_flux_defect > tolerance
        {
            return Err(TensorShellRadiationError::PowerBalance {
                coupled, port_work, flux: power.relative_flux_defect, tolerance,
            });
        }
        checkpoint(cx)?;
        Ok(TensorShellRadiationResponse {
            motion, field, power, relative_coupled_power_defect: coupled,
            relative_port_work_defect: port_work,
        })
    }
}

// Scaled real power identity a=b+c; no arbitrary one-watt floor or hidden repair.
fn relative_balance(a: f64, b: f64, c: f64) -> f64 {
    let scale = a.abs().max(b.abs()).max(c.abs());
    if scale == 0.0 { return 0.0; }
    let (a, b, c) = (a/scale, b/scale, c/scale);
    (a-b-c).abs()/(a.abs()+b.abs()+c.abs())
}

/// Paired immutable physical response. Independent repeated forces share the same
/// loaded factor; each result owns its coherent densities and work diagnostics.
pub struct TensorShellRadiationResponse<'p, 'geometry> {
    motion: HarmonicResponse,
    field: CoherentVelocityPortField<'p, 'geometry>,
    power: CoherentRadiationPower,
    relative_coupled_power_defect: f64,
    relative_port_work_defect: f64,
}
impl<'p, 'geometry> TensorShellRadiationResponse<'p, 'geometry> {
    /// Original displacement/rotation, velocity, force work and intrinsic heat.
    pub const fn motion(&self) -> &HarmonicResponse { &self.motion }
    /// Evaluate with the native coherent pressure/particle-velocity receiver APIs.
    pub const fn field(&self) -> &CoherentVelocityPortField<'p, 'geometry> { &self.field }
    pub const fn power(&self) -> CoherentRadiationPower { self.power }
    pub const fn relative_coupled_power_defect(&self) -> f64 { self.relative_coupled_power_defect }
    pub const fn relative_port_work_defect(&self) -> f64 { self.relative_port_work_defect }
}
