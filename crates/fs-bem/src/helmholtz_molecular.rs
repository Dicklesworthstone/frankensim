//! Molecular bulk relaxation in the exterior reaction and observation problem.
//!
//! A uniform stationary gas obeys `div(v)=i*omega*p/K(omega)` and
//! `grad(p)=i*rho*omega*v` for exp(-i omega t). Thus `k²=rho*omega²/K`,
//! but the velocity Neumann trace remains `q=i*rho*omega*v`: replacing
//! omega by k times the equilibrium sound speed changes the boundary physics.
//! The original material caloric law supplies K; the original complex BEM
//! assembly and Green representation supply the pressure and its gradient.
//!
//! This is the inviscid-momentum, relaxing-compression model. Exterior shear
//! viscosity, heat conduction, wall layers and mean flow are not represented.
//! In particular this does not reinterpret a plane-wave attenuation coefficient
//! as a bulk modulus or add an independent absorption filter after radiation.
//! Boundary work includes molecular absorption in the exterior and outgoing
//! energy; it is not flux measured on a distant sphere in an absorbing medium.

use super::aperture::{self, Aperture, ApertureImpedance};
use super::laplace::complex::{self, ComplexError, ComplexOptions};
use super::laplace::{DampedWaveNumber, PointField, fields_from_traces};
use super::*;
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_material::{MaterialError, relaxation::RelaxationState};
use std::sync::Arc;

/// Material, geometric, work-envelope, solve or cancellation refusal.
#[derive(Debug)]
pub enum MolecularError {
    Material(MaterialError),
    Boundary(HelmholtzError),
    Operator(ComplexError),
    LinearSolve(LuComplexError),
}
impl core::fmt::Display for MolecularError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(e) => e.fmt(f),
            Self::Boundary(e) => e.fmt(f),
            Self::Operator(e) => e.fmt(f),
            Self::LinearSolve(e) => e.fmt(f),
        }
    }
}
impl std::error::Error for MolecularError {}
impl From<MaterialError> for MolecularError {
    fn from(e: MaterialError) -> Self {
        Self::Material(e)
    }
}
impl From<HelmholtzError> for MolecularError {
    fn from(e: HelmholtzError) -> Self {
        Self::Boundary(e)
    }
}
impl From<ComplexError> for MolecularError {
    fn from(e: ComplexError) -> Self {
        Self::Operator(e)
    }
}
impl From<LuComplexError> for MolecularError {
    fn from(e: LuComplexError) -> Self {
        Self::LinearSolve(e)
    }
}

fn bad(what: &'static str) -> MolecularError {
    HelmholtzError::BadParameter { what }.into()
}
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MolecularError> {
    if cancelled() {
        Err(HelmholtzError::Cancelled.into())
    } else {
        Ok(())
    }
}

/// One evaluated material and declared harmonic band, retained by all results.
/// Molecular source evaluation occurs once before BEM assembly at each frequency.
#[derive(Debug, Clone)]
pub struct MolecularExterior {
    source: Arc<RelaxationState>,
    band_rad_s: [f64; 2],
}

/// Immutable solution, binding reaction, finite receivers, gas and geometry.
#[derive(Debug, Clone)]
pub struct MolecularRadiation {
    source: Arc<RelaxationState>,
    omega: f64,
    k: C64,
    fingerprint: u64,
    pressure: Vec<C64>,
    velocity: Vec<C64>,
    boundary_work_w: f64,
}

impl MolecularExterior {
    /// Bind an already evaluated material; no gas composition or reservoir
    /// constants are inferred. The band is a caller declaration, not a source
    /// accuracy certificate or a guarantee of BEM spatial convergence.
    pub fn new(source: RelaxationState, band_rad_s: [f64; 2]) -> Result<Self, MolecularError> {
        source.gas().validate_coherence()?;
        if !band_rad_s.iter().all(|w| w.is_finite() && *w > 0.0) || band_rad_s[1] <= band_rad_s[0] {
            return Err(bad(
                "molecular exterior needs a positive increasing angular-frequency band",
            ));
        }
        Ok(Self {
            source: Arc::new(source),
            band_rad_s,
        })
    }
    pub fn source(&self) -> &RelaxationState {
        &self.source
    }
    pub const fn band_rad_s(&self) -> [f64; 2] {
        self.band_rad_s
    }

    /// Solve a batch using one original A/B assembly and one LU. The values-only
    /// specialization holds two N² matrices, without computing unused frequency
    /// derivatives. The existing triangle CBIE or centroid/disc Burton–Miller
    /// rule is preserved; |k| controls the spatial-resolution admission.
    ///
    /// Cancellation is checked at source admission, each assembly/RHS row and
    /// each batch field. The original cancellable LU polls at its elimination
    /// and substitution boundaries, including its declared-worker block policy.
    pub fn solve_batch(
        &self,
        surface: &SpherePanels,
        omega: f64,
        fields: &[&[C64]],
        formulation: Formulation,
        options: &ComplexOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<MolecularRadiation>, MolecularError> {
        poll(&mut cancelled)?;
        if !omega.is_finite() || omega < self.band_rad_s[0] || omega > self.band_rad_s[1] {
            return Err(bad(
                "molecular exterior frequency is outside its declared band",
            ));
        }
        let n = surface.areas().len();
        if fields.is_empty()
            || fields.len() > MAX_RADIATION_FIELDS_PER_BATCH
            || fields
                .iter()
                .any(|v| v.len() != n || v.iter().any(|&x| !finite(x)))
        {
            return Err(bad(
                "molecular exterior needs a finite, bounded boundary-field batch",
            ));
        }
        let k = self.source.relaxing_wavenumber(omega)?;
        if !finite(k) || k.re <= 0.0 || k.im < 0.0 {
            return Err(bad(
                "molecular exterior wavenumber must be outgoing and passive",
            ));
        }
        let op =
            complex::assemble_values_with_cancel(surface, k, formulation, options, &mut cancelled)?;
        poll(&mut cancelled)?;
        let lu = lu_complex_with_cancel(&op.a, n, &mut cancelled)?;
        // The factor owns its coefficients; retain only B across the RHS batch.
        drop(op.a);
        poll(&mut cancelled)?;
        let q_scale = C64::new(0.0, self.source.gas().density * omega);
        if !finite(q_scale) {
            return Err(bad("molecular exterior Neumann scale is unrepresentable"));
        }
        let fingerprint = surface_fingerprint(surface);
        let mut out = Vec::with_capacity(fields.len());
        for &velocity in fields {
            poll(&mut cancelled)?;
            let mut pressure = vec![C64::ZERO; n];
            for i in 0..n {
                poll(&mut cancelled)?;
                for (j, &v) in velocity.iter().enumerate() {
                    pressure[i] = pressure[i] + op.b[i * n + j] * (q_scale * v);
                }
            }
            lu.solve_with_cancel(&mut pressure, &mut cancelled)?;
            let mut boundary_work_w = 0.0;
            for ((&p, &v), &area) in pressure.iter().zip(velocity).zip(surface.areas()) {
                if !finite(p) {
                    return Err(bad("nonfinite molecular boundary pressure"));
                }
                boundary_work_w += 0.5 * (p * v.conj()).re * area;
            }
            if !boundary_work_w.is_finite() {
                return Err(bad("nonfinite molecular boundary work"));
            }
            out.push(MolecularRadiation {
                source: Arc::clone(&self.source),
                omega,
                k,
                fingerprint,
                pressure,
                velocity: velocity.to_vec(),
                boundary_work_w,
            });
        }
        poll(&mut cancelled)?;
        Ok(out)
    }

    /// Compute all self/mutual loads at one frequency, retaining the same solved
    /// traces for receiver evaluation. `impedance().acoustic(0)` is the original
    /// plane-port Pa·s/m³ matrix consumed by `CoupledTermination`, including
    /// `MolecularPortScattering::response`; use outward mouth flows in that order.
    /// Reference planes and real aperture shapes have the original admission.
    /// No end correction, chimney or duplicate exterior load is added here.
    pub fn aperture_response(
        &self,
        surface: &SpherePanels,
        apertures: &[Aperture],
        omega: f64,
        formulation: Formulation,
        options: &ComplexOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularApertureResponse, MolecularError> {
        poll(&mut cancelled)?;
        aperture::admit_apertures(surface, apertures)?;
        let fields: Vec<Vec<_>> = apertures
            .iter()
            .map(|a| a.shape.iter().map(|&v| C64::from_re(v)).collect())
            .collect();
        let refs: Vec<_> = fields.iter().map(Vec::as_slice).collect();
        let solutions =
            self.solve_batch(surface, omega, &refs, formulation, options, &mut cancelled)?;
        let pressures: Vec<_> = solutions.iter().map(|s| s.pressure.as_slice()).collect();
        let impedance = aperture::project_pressures(surface, apertures, &pressures)?;
        let mut declaration = aperture::mix(
            aperture::aperture_declaration(surface, apertures),
            0x4d4f4c4543554c41,
        );
        let gas = self.source.gas();
        for v in [
            omega,
            solutions[0].k.re,
            solutions[0].k.im,
            gas.density,
            gas.temperature,
            gas.pressure,
            gas.specific_gas_constant,
            gas.specific_heat_cp,
            self.source.frozen_cv_j_kg_k(),
            self.band_rad_s[0],
            self.band_rad_s[1],
        ] {
            declaration = aperture::mix(declaration, v.to_bits());
        }
        for p in self.source.processes() {
            for v in [p.heat_capacity_j_kg_k, p.relaxation_time_s] {
                declaration = aperture::mix(declaration, v.to_bits());
            }
        }
        declaration = aperture::mix(
            declaration,
            match formulation {
                Formulation::PlainCbie => 0,
                Formulation::BurtonMiller => 1,
                Formulation::BurtonMillerWrongAlphaSign => 2,
            },
        );
        let volume_gains = apertures
            .iter()
            .map(|a| {
                a.shape
                    .iter()
                    .zip(surface.areas())
                    .map(|(s, area)| s * area)
                    .sum()
            })
            .collect();
        let n = apertures.len();
        let matrix = ApertureImpedance {
            declaration,
            apertures: n,
            volume_gains,
            rows: vec![aperture::row(omega, impedance, n, "molecular-dense")],
        };
        poll(&mut cancelled)?;
        Ok(MolecularApertureResponse { matrix, solutions })
    }
}

impl MolecularRadiation {
    pub fn source(&self) -> &RelaxationState {
        &self.source
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn wavenumber_per_m(&self) -> C64 {
        self.k
    }
    pub fn pressure_pa(&self) -> &[C64] {
        &self.pressure
    }
    pub fn velocity_m_s(&self) -> &[C64] {
        &self.velocity
    }
    /// Raw discrete mean surface work [W], without positivity clipping.
    /// Spatial discretization can produce a negative numerical defect; callers
    /// must refine or refuse it, not interpret this value as a passive fit.
    pub const fn boundary_work_w(&self) -> f64 {
        self.boundary_work_w
    }

    /// Finite pressure and gradient from the SAME medium, k and boundary solve.
    /// The original centroid Green representation and closed-surface clearance
    /// admission apply; no near-field quadrature convergence is implied.
    pub fn observe(
        &self,
        surface: &SpherePanels,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<PointField>, MolecularError> {
        poll(&mut cancelled)?;
        if self.fingerprint != surface_fingerprint(surface) {
            return Err(bad(
                "molecular solution belongs to a different source surface",
            ));
        }
        crate::near_field::Geometry::new(surface, points, minimum_clearance_m)?;
        let mut out = Vec::with_capacity(points.len());
        for point in points {
            poll(&mut cancelled)?;
            out.extend(fields_from_traces(
                surface,
                &self.pressure,
                &self.velocity,
                DampedWaveNumber {
                    k: self.k,
                    phase_origin_m: 0.0,
                },
                C64::new(0.0, self.source.gas().density * self.omega),
                &[*point],
            )?);
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
}

/// One-frequency aperture matrix plus its immutable boundary solutions.
#[derive(Debug, Clone)]
pub struct MolecularApertureResponse {
    matrix: ApertureImpedance,
    solutions: Vec<MolecularRadiation>,
}
impl MolecularApertureResponse {
    pub fn impedance(&self) -> &ApertureImpedance {
        &self.matrix
    }
    pub fn boundary_solutions(&self) -> &[MolecularRadiation] {
        &self.solutions
    }

    /// Observe solved outward aperture volume flows [m³/s]. Superposition occurs
    /// on the retained boundary traces before one receiver traversal, avoiding
    /// one Npanels×Nreceivers traversal per aperture. This shares both the load's
    /// source gas and all mutual scattering; no independent propagation filter.
    pub fn observe_volume_flows(
        &self,
        surface: &SpherePanels,
        volume_flows: &[C64],
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<PointField>, MolecularError> {
        poll(&mut cancelled)?;
        if volume_flows.len() != self.solutions.len()
            || volume_flows.iter().any(|&v| !finite(v))
            || self
                .matrix
                .volume_gains
                .iter()
                .any(|v| !v.is_finite() || *v == 0.0)
        {
            return Err(bad(
                "molecular aperture observation needs finite matched flows and nonzero volume gains",
            ));
        }
        let first = &self.solutions[0];
        let mut combined = MolecularRadiation {
            source: Arc::clone(&first.source),
            omega: first.omega,
            k: first.k,
            fingerprint: first.fingerprint,
            pressure: vec![C64::ZERO; first.pressure.len()],
            velocity: vec![C64::ZERO; first.velocity.len()],
            boundary_work_w: 0.0,
        };
        for ((solution, &flow), &gain) in self
            .solutions
            .iter()
            .zip(volume_flows)
            .zip(&self.matrix.volume_gains)
        {
            poll(&mut cancelled)?;
            let amplitude = flow.scale(1.0 / gain);
            if !finite(amplitude) {
                return Err(bad(
                    "molecular aperture velocity amplitude is unrepresentable",
                ));
            }
            for i in 0..combined.pressure.len() {
                combined.pressure[i] = combined.pressure[i] + solution.pressure[i] * amplitude;
                combined.velocity[i] = combined.velocity[i] + solution.velocity[i] * amplitude;
                if !finite(combined.pressure[i]) || !finite(combined.velocity[i]) {
                    return Err(bad("nonfinite molecular aperture superposition"));
                }
            }
        }
        combined.observe(surface, points, minimum_clearance_m, cancelled)
    }
}

#[cfg(test)]
#[path = "helmholtz_molecular/tests.rs"]
mod tests;
