//! Penetrable homogeneous-fluid scattering on the original closed BEM surface.
//!
//! Both sides retain their original FluidAcoustics owner and one coherent
//! snapshot at a shared omega. These are stationary, inviscid, nondispersive,
//! positive-density fluids: no solid shear, empirical loss or porous-frame law.
//! Pressure and outward normal velocity are continuous across the interface.
//!
//! With exp(-i omega t), n pointing from the bounded interior to the exterior,
//! q_o=i omega rho_o v, q_i=eta q_o, eta=rho_i/rho_o, the Green trace identities are
//! (D_o-I/2)p-S_o q_o=-p_inc, (D_i+I/2)p-eta S_i q_o=0,
//! N_o p-(D'_o+I/2)q_o=-q_inc, N_i p-eta(D'_i-I/2)q_o=0.
//! Subtracting the exterior identities from the interior identities gives
//! [I+D_i-D_o, S_o-eta S_i; N_i-N_o, (1+eta)I/2+D'_o-eta D'_i].
//! This is the direct transmission combination of Bruno & Perez-Arancibia,
//! Proc. R. Soc. A 473 (2017), 20170161, equations 3.3–3.7, after changing
//! their interior normal derivative to q_o (doi:10.1098/rspa.2017.0161).
//!
//! The solve uses (p, q_o/k_o), so both unknowns have pressure units. It solves
//! directly for scattered traces against the contrast operator M-I, retaining
//! weak reflections that total-minus-incident subtraction would erase.
//! The common hypersingular static kernel cancels analytically. Off-diagonal
//! panels and equivalent-disc self terms retain the original constant-panel
//! Helmholtz discretization; this is not a new Galerkin or accuracy guarantee.
use super::bc::{IncidentField, IncidentTraces};
use super::*;
use crate::near_field;
use fs_exec::CancelGate;
use fs_la::eigen_complex::{
    LuComplexError, lu_complex_blocked_with_cancel, lu_complex_blocked_workspace_bytes,
};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};

/// A transmission solve has two pressure-scaled unknowns per panel.
pub const MAX_TRANSMISSION_PANELS: usize = MAX_DENSE_PANELS / 2;

/// Explicit dense admission, applied before quadratic allocation.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub maximum_panels: usize,
    /// Peak dense matrices, LU descriptors and vectors. Excludes thread/runtime
    /// metadata and allocator overhead; not a process-RSS bound.
    pub maximum_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_panels: MAX_TRANSMISSION_PANELS,
            maximum_bytes: 512 * 1024 * 1024,
        }
    }
}

/// Original material refusal is distinct from geometry/solver refusal.
#[derive(Debug)]
pub enum TransmissionError {
    Material(MaterialError),
    Helmholtz(HelmholtzError),
}
impl From<MaterialError> for TransmissionError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<HelmholtzError> for TransmissionError {
    fn from(error: HelmholtzError) -> Self {
        Self::Helmholtz(error)
    }
}
impl core::fmt::Display for TransmissionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => core::fmt::Display::fmt(e, f),
            Self::Helmholtz(e) => core::fmt::Display::fmt(e, f),
        }
    }
}
impl std::error::Error for TransmissionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(e) => Some(e),
            Self::Helmholtz(e) => Some(e),
        }
    }
}
fn bad(what: &'static str) -> HelmholtzError {
    HelmholtzError::BadParameter { what }
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(gate: &CancelGate) -> Result<(), HelmholtzError> {
    if gate.is_requested() {
        Err(HelmholtzError::Cancelled)
    } else {
        Ok(())
    }
}
fn lu_error(error: LuComplexError) -> HelmholtzError {
    match error {
        LuComplexError::Cancelled => HelmholtzError::Cancelled,
        LuComplexError::Singular { .. } => HelmholtzError::Singular,
    }
}
fn medium(state: FluidAcousticState) -> Medium {
    Medium {
        density: state.density_kg_m3(),
        sound_speed: state.sound_speed_m_s(),
    }
}

/// One admitted fluid inclusion, frequency and reusable coupled factorization.
///
/// The retained triangles must form ONE outward, non-self-intersecting closed
/// component. Original closure checks do not certify global nonintersection.
/// Exterior incident sources must be regular throughout that component.
pub struct FluidInterface<'a> {
    surface: &'a SpherePanels,
    exterior_source: &'a dyn FluidAcoustics,
    interior_source: &'a dyn FluidAcoustics,
    exterior: FluidAcousticState,
    interior: FluidAcousticState,
    omega: f64,
    k_o: f64,
    k_i: f64,
    ppw: f64,
    contrast: Vec<C64>,
    factor: Option<LuComplex>,
    condition_lower_bound: f64,
}
impl<'a> FluidInterface<'a> {
    /// Query each original material exactly once, admit BOTH wavelengths, and
    /// factor the coupled interface. Cancellation is checked during assembly
    /// and by the original cancellable blocked LU owner.
    pub fn new(
        surface: &'a SpherePanels,
        exterior_source: &'a dyn FluidAcoustics,
        interior_source: &'a dyn FluidAcoustics,
        omega_rad_s: f64,
        options: Options,
        gate: &CancelGate,
    ) -> Result<Self, TransmissionError> {
        poll(gate)?;
        if !omega_rad_s.is_finite()
            || omega_rad_s <= 0.0
            || options.maximum_panels == 0
            || options.maximum_panels > MAX_TRANSMISSION_PANELS
        {
            return Err(bad("invalid fluid-interface frequency or dense panel budget").into());
        }
        let n = surface.areas().len();
        if n > options.maximum_panels {
            return Err(HelmholtzError::WorkCap { panels: n }.into());
        }
        surface
            .validate()
            .map_err(|_| bad("invalid fluid-interface surface"))?;
        let triangles = surface
            .triangles()
            .ok_or_else(|| bad("fluid interfaces require retained triangles"))?;
        if near_field::components(triangles)?.len() != 1 {
            return Err(bad("fluid interface requires one closed outward component").into());
        }
        poll(gate)?;
        let exterior = exterior_source.acoustic_state()?;
        poll(gate)?;
        let interior = interior_source.acoustic_state()?;
        poll(gate)?;
        let k_o = omega_rad_s / exterior.sound_speed_m_s();
        let k_i = omega_rad_s / interior.sound_speed_m_s();
        let ppw = admit_resolution(surface, k_o)?.min(admit_resolution(surface, k_i)?);
        let eta = interior.density_kg_m3() / exterior.density_kg_m3();
        let delta_rho =
            (interior.density_kg_m3() - exterior.density_kg_m3()) / exterior.density_kg_m3();
        if [
            eta,
            k_o,
            k_i,
            omega_rad_s * exterior.density_kg_m3(),
            omega_rad_s * interior.density_kg_m3(),
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
            || !delta_rho.is_finite()
        {
            return Err(bad("unrepresentable fluid-interface contrast or frequency scale").into());
        }
        let m = 2 * n;
        let workers = DENSE_WORKERS.with(core::cell::Cell::get);
        let complex_bytes = core::mem::size_of::<C64>();
        let permutation_bytes = m * core::mem::size_of::<usize>();
        let dense_bytes = 2 * m * m * complex_bytes + permutation_bytes;
        let lu_workspace = lu_complex_blocked_workspace_bytes(m, workers)
            .ok_or_else(|| bad("fluid-interface LU workspace byte count overflow"))?;
        let factor_peak = dense_bytes + m * complex_bytes + lu_workspace;
        let solve_peak = dense_bytes + 6 * m * complex_bytes;
        let required_bytes = factor_peak.max(solve_peak);
        if required_bytes > options.maximum_bytes {
            return Err(HelmholtzError::CacheBudget {
                required_bytes,
                maximum_bytes: options.maximum_bytes,
            }
            .into());
        }
        let (contrast, factor, condition_lower_bound) = if medium(exterior) == medium(interior) {
            // Transparent material identity is exact, including finite sources.
            (Vec::new(), None, 1.0)
        } else {
            let mut matrix = assemble(surface, k_o, k_i, eta, delta_rho, gate)?;
            let diagonal: Vec<_> = (0..m).map(|i| matrix[i * m + i]).collect();
            for i in 0..m {
                matrix[i * m + i] = matrix[i * m + i] + C64::ONE;
            }
            let factor =
                lu_complex_blocked_with_cancel(&matrix, m, workers, || gate.is_requested())
                    .map_err(lu_error)?;
            poll(gate)?;
            let condition = condition_lower_bound(&matrix, m, &factor);
            poll(gate)?;
            if !condition.is_finite() {
                return Err(bad("nonfinite fluid-interface condition diagnostic").into());
            }
            // Restore SAVED contrasts, not rounded (1+contrast)-1.
            for i in 0..m {
                matrix[i * m + i] = diagonal[i];
            }
            (matrix, Some(factor), condition)
        };
        poll(gate)?;
        Ok(Self {
            surface,
            exterior_source,
            interior_source,
            exterior,
            interior,
            omega: omega_rad_s,
            k_o,
            k_i,
            ppw,
            contrast,
            factor,
            condition_lower_bound,
        })
    }
    pub fn exterior_source(&self) -> &'a dyn FluidAcoustics {
        self.exterior_source
    }
    pub fn interior_source(&self) -> &'a dyn FluidAcoustics {
        self.interior_source
    }
    pub fn exterior_state(&self) -> FluidAcousticState {
        self.exterior
    }
    pub fn interior_state(&self) -> FluidAcousticState {
        self.interior
    }
    pub fn angular_frequency_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn exterior_wavenumber_rad_m(&self) -> f64 {
        self.k_o
    }
    pub fn interior_wavenumber_rad_m(&self) -> f64 {
        self.k_i
    }

    /// Solve immutable incident traces bound to this exterior medium, surface
    /// and frequency. Repeated calls reuse the same factor; materials are not
    /// re-queried. Caller-authored traces retain their original regularity duty.
    pub fn scatter<'p>(
        &'p self,
        incident: &IncidentTraces,
        gate: &CancelGate,
    ) -> Result<TransmissionSolution<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        let exterior = medium(self.exterior);
        incident.validate_binding(self.surface, self.k_o, exterior)?;
        let n = self.surface.areas().len();
        let m = 2 * n;
        let z_o = self.exterior.characteristic_impedance_pa_s_m();
        let mut incoming = incident.pressure_pa().to_vec();
        incoming.extend(
            incident
                .normal_velocity_m_s()
                .iter()
                .map(|&v| v * C64::new(0.0, z_o)),
        );
        if incoming.iter().any(|&z| !finite(z)) {
            return Err(bad("nonfinite pressure-scaled incident interface field"));
        }
        let mut scattered = vec![C64::ZERO; m];
        if let Some(factor) = &self.factor {
            for (i, value) in scattered.iter_mut().enumerate() {
                poll(gate)?;
                for (j, &b) in incoming.iter().enumerate() {
                    *value = *value - self.contrast[i * m + j] * b;
                }
            }
            if scattered.iter().any(|&z| !finite(z)) {
                return Err(bad("nonfinite fluid-interface contrast drive"));
            }
            factor
                .solve_with_cancel(&mut scattered, || gate.is_requested())
                .map_err(lu_error)?;
        }
        if scattered.iter().any(|&z| !finite(z)) {
            return Err(bad("nonfinite fluid-interface scattered solution"));
        }
        let pressure = scattered[..n].to_vec();
        let transparent = self.factor.is_none();
        let velocity: Vec<_> = if transparent {
            vec![C64::ZERO; n]
        } else {
            scattered[n..]
                .iter()
                .map(|&u| C64::new(u.im / z_o, -u.re / z_o))
                .collect()
        };
        let total_pressure: Vec<_> = if transparent {
            incident.pressure_pa().to_vec()
        } else {
            pressure
                .iter()
                .zip(incident.pressure_pa())
                .map(|(&s, &i)| s + i)
                .collect()
        };
        let total_velocity: Vec<_> = if transparent {
            incident.normal_velocity_m_s().to_vec()
        } else {
            velocity
                .iter()
                .zip(incident.normal_velocity_m_s())
                .map(|(&s, &i)| s + i)
                .collect()
        };
        if total_pressure
            .iter()
            .chain(&total_velocity)
            .chain(&velocity)
            .any(|&z| !finite(z))
        {
            return Err(bad("nonfinite total fluid-interface traces"));
        }
        let (_, work_interval, work) = evaluate_radiated_power(
            &total_pressure,
            &total_velocity,
            self.surface.areas(),
            exterior,
        );
        let radiation = self.radiation(self.k_o, exterior, pressure, velocity)?;
        if !finite(work.surface_power) || !work.apparent_power.is_finite() {
            return Err(bad("nonfinite fluid-interface power diagnostics"));
        }
        poll(gate)?;
        Ok(TransmissionSolution {
            interface: self,
            scattered: radiation,
            total_pressure,
            total_velocity,
            interface_work: work,
            interface_work_roundoff_interval: work_interval,
        })
    }

    /// Sample the original coherent plane-wave/monopole owner and retain that
    /// same field for finite exterior total-pressure observations.
    pub fn scatter_field<'p>(
        &'p self,
        field: IncidentField<'p>,
        minimum_source_clearance_m: f64,
        gate: &CancelGate,
    ) -> Result<FluidScatteredField<'p, 'a>, HelmholtzError> {
        poll(gate)?;
        let traces = IncidentTraces::from_field(
            self.surface,
            self.k_o,
            medium(self.exterior),
            field,
            minimum_source_clearance_m,
        )?;
        let solution = self.scatter(&traces, gate)?;
        Ok(FluidScatteredField {
            incident: field,
            solution,
        })
    }
    fn radiation(
        &self,
        k: f64,
        medium: Medium,
        pressure: Vec<C64>,
        velocity: Vec<C64>,
    ) -> Result<RadiationSolution, HelmholtzError> {
        let (power, interval, diagnostics) =
            evaluate_radiated_power(&pressure, &velocity, self.surface.areas(), medium);
        if !power.is_finite()
            || !finite(diagnostics.surface_power)
            || !diagnostics.apparent_power.is_finite()
        {
            return Err(bad("nonfinite transmission representation power"));
        }
        Ok(RadiationSolution {
            pressure,
            velocity,
            k,
            medium,
            surface_fingerprint: surface_fingerprint(self.surface),
            panels_per_wavelength: self.ppw,
            radiated_power: power,
            radiated_power_roundoff_interval: interval,
            power_diagnostics: diagnostics,
            condition_lower_bound: self.condition_lower_bound,
            dense_cap_utilization: 2.0 * self.surface.areas().len() as f64
                / MAX_DENSE_PANELS as f64,
        })
    }
}

/// Outgoing scattered traces plus the continuous physical interface traces.
/// The prepared interface keeps both original material owners alive.
pub struct TransmissionSolution<'p, 's> {
    interface: &'p FluidInterface<'s>,
    scattered: RadiationSolution,
    total_pressure: Vec<C64>,
    total_velocity: Vec<C64>,
    /// Signed total surface work. Real work is zero in the continuum for this
    /// source-free lossless interior; its numerical defect is NOT absorption.
    /// Imaginary work describes the interior's reactive energy exchange.
    pub interface_work: RadiationPowerDiagnostics,
    /// Floating-point panel-dot enclosure only; excludes discretization error.
    pub interface_work_roundoff_interval: (f64, f64),
}
impl TransmissionSolution<'_, '_> {
    /// Original outgoing representation for existing finite/far-field owners.
    pub fn scattered(&self) -> &RadiationSolution {
        &self.scattered
    }
    pub fn total_pressure_pa(&self) -> &[C64] {
        &self.total_pressure
    }
    /// The SAME normal velocity on both sides, along the outward body normal.
    pub fn total_normal_velocity_m_s(&self) -> &[C64] {
        &self.total_velocity
    }

    /// Interior pressure from S_i q_i-D_i p with the original adaptive triangle
    /// receiver quadrature, exact winding and point-to-triangle clearance.
    pub fn interior_pressure_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
    ) -> Result<near_field::Observation, HelmholtzError> {
        let interface = self.interface;
        let geometry =
            near_field::Geometry::new_interior(interface.surface, points_m, minimum_clearance_m)?;
        // Existing receiver evaluates Dp-Sq. Negating BOTH original traces
        // gives the interior representation; rho_i supplies q_i=i omega rho_i v.
        let traces = interface.radiation(
            interface.k_i,
            medium(interface.interior),
            self.total_pressure.iter().map(|&p| -p).collect(),
            self.total_velocity.iter().map(|&v| -v).collect(),
        )?;
        geometry
            .prepare(interface.k_i, medium(interface.interior), options)?
            .evaluate(&traces)
    }
}

/// One coherent incident field and its solved transmitting body.
pub struct FluidScatteredField<'p, 's> {
    incident: IncidentField<'p>,
    solution: TransmissionSolution<'p, 's>,
}
impl FluidScatteredField<'_, '_> {
    pub fn solution(&self) -> &TransmissionSolution<'_, '_> {
        &self.solution
    }
    /// Incident plus scattered pressure at the original finite receivers.
    /// Quadrature estimates cover the scattered representation only.
    pub fn total_pressure_at(
        &self,
        points_m: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: near_field::Options,
    ) -> Result<near_field::Observation, HelmholtzError> {
        let interface = self.solution.interface;
        let geometry = near_field::Geometry::new(interface.surface, points_m, minimum_clearance_m)?;
        let incident: Vec<_> = points_m
            .iter()
            .map(|&x| {
                self.incident
                    .pressure_velocity_at(interface.k_o, medium(interface.exterior), x)
            })
            .collect::<Result<_, _>>()?;
        let mut observed = geometry
            .prepare(interface.k_o, medium(interface.exterior), options)?
            .evaluate(self.solution.scattered())?;
        for (p, incident) in observed.pressure.iter_mut().zip(incident) {
            *p = *p + incident.pressure_pa;
            if !finite(*p) {
                return Err(bad("nonfinite total transmitted-body exterior field"));
            }
        }
        Ok(observed)
    }
}

// Analytic differences of the ORIGINAL outgoing point kernels. Midpoint
// integration of d[e^(iz)(1-iz)]/dz=z e^(iz) and
// d[e^(iz)(z²+2iz-2)]/dz=i z² e^(iz) removes their common static singularity
// before floating-point subtraction, including low frequency/weak contrast.
fn kernel_difference(
    k_o: f64,
    k_i: f64,
    x: [f64; 3],
    nx: [f64; 3],
    y: [f64; 3],
    ny: [f64; 3],
) -> [C64; 4] {
    let delta = sub(x, y);
    let r = norm(delta);
    let h = 0.5 * (k_i - k_o) * r;
    let mid = k_o * r + h;
    let sine = det::sin(h);
    let h2 = h * h;
    let moment = if h.abs() < 0.01 {
        h * h2 * (1.0 / 3.0 - h2 / 30.0 + h2 * h2 / 840.0 - h2 * h2 * h2 / 45360.0)
    } else {
        sine - h * det::cos(h)
    };
    let e = expik(mid);
    let f = (e * C64::new(mid * sine, moment)).scale(2.0);
    let g = (e * C64::new(-2.0 * mid * moment, (mid * mid + h2) * sine - 2.0 * moment)).scale(2.0);
    let xn = dot(nx, delta) / r;
    let yn = dot(ny, delta) / r;
    let scale = 1.0 / (4.0 * core::f64::consts::PI * r);
    [
        (e * C64::new(0.0, 2.0 * sine)).scale(scale),
        f.scale(scale * yn / r),
        f.scale(-scale * xn / r),
        (f.scale(dot(nx, ny) - xn * yn) + g.scale(xn * yn)).scale(scale / (r * r)),
    ]
}

// Same equivalent-disc self integrals as self_terms, evaluated as factored
// series differences. Existing wavelength admission gives k*a <= sqrt(pi)/3,
// so 24 terms leave a remainder far below f64 roundoff, without 1/k cancellation.
fn self_contrast(k_o: f64, k_i: f64, area: f64) -> (C64, C64, C64) {
    let a = (area / core::f64::consts::PI).sqrt();
    let (z_o, z_i, dz) = (k_o * a, k_i * a, (k_i - k_o) * a);
    let (mut factorial, mut p_o, mut p_i, mut difference) = (1.0, 1.0, 1.0, 0.0);
    let (mut phase, mut single, mut ds, mut dn) =
        (C64::ONE, C64::from_re(0.5 * a), C64::ZERO, C64::ZERO);
    for n in 1..=24 {
        difference = p_i * dz + z_o * difference;
        p_o *= z_o;
        p_i *= z_i;
        factorial *= f64::from(n);
        phase = phase * C64::new(0.0, 1.0);
        single = single + phase.scale(0.5 * a * p_o / (factorial * f64::from(n + 1)));
        ds = ds + phase.scale(0.5 * a * difference / (factorial * f64::from(n + 1)));
        if n >= 2 {
            dn = dn - phase.scale(difference / (2.0 * a * factorial));
        }
    }
    (single, ds, dn)
}

fn assemble(
    surface: &SpherePanels,
    k_o: f64,
    k_i: f64,
    eta: f64,
    delta_rho: f64,
    gate: &CancelGate,
) -> Result<Vec<C64>, HelmholtzError> {
    let n = surface.areas().len();
    let m = 2 * n;
    let mut matrix = Vec::new();
    matrix
        .try_reserve_exact(m * m)
        .map_err(|_| HelmholtzError::AllocationFailed {
            what: "fluid-interface matrix",
        })?;
    matrix.resize(m * m, C64::ZERO);
    for i in 0..n {
        poll(gate)?;
        let (x, nx) = (surface.centroids()[i], surface.normals()[i]);
        // Original D_0 closure and its weighted adjoint use the same geometry.
        let mut d0 = C64::ZERO;
        for j in 0..n {
            if j != i {
                d0 = d0
                    + kernels(0.0, x, nx, surface.centroids()[j], surface.normals()[j])
                        .1
                        .scale(surface.areas()[j]);
            }
        }
        for j in 0..n {
            let area = surface.areas()[j];
            let (s, dp, ds, dd, ddp, dn) = if i == j {
                let (s, ds, dn) = self_contrast(k_o, k_i, area);
                (s, -C64::from_re(0.5) - d0, ds, C64::ZERO, C64::ZERO, dn)
            } else {
                let (y, ny) = (surface.centroids()[j], surface.normals()[j]);
                let (s, _, dp, _) = kernels(k_o, x, nx, y, ny);
                let [ds, dd, ddp, dn] =
                    kernel_difference(k_o, k_i, x, nx, y, ny).map(|z| z.scale(area));
                (s.scale(area), dp.scale(area), ds, dd, ddp, dn)
            };
            matrix[i * m + j] = dd;
            matrix[i * m + n + j] = (-s.scale(delta_rho) - ds.scale(eta)).scale(k_o);
            matrix[(n + i) * m + j] = dn.scale(1.0 / k_o);
            matrix[(n + i) * m + n + j] = -dp.scale(delta_rho) - ddp.scale(eta)
                + C64::from_re(if i == j { 0.5 * delta_rho } else { 0.0 });
        }
    }
    if matrix.iter().any(|&z| !finite(z)) {
        return Err(bad("nonfinite fluid-interface operator"));
    }
    Ok(matrix)
}
