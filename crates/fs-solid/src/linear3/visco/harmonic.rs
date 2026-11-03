//! Forced resonance on the original spatial tensor-Prony pencil.
//!
//! At the requested positive frequency, solve
//! (K'(omega) - omega^2 M - i K''(omega)) u = f with peak phasors under
//! exp(-i omega t). The existing harmonic owner receives exactly
//! C_eff = K'' / omega; this is a frequency-domain identity, NOT a
//! memoryless replacement of the hereditary material in time integration.
//! Material dispersion, anisotropy, mutual coupling and consistent inertia
//! remain in the full original coordinates. No modal basis is required.

use super::{ViscoElementError, ViscoElements, checkpoint};
use fs_exec::Cx;
use fs_math::{c64::C64, det};
use fs_modal::harmonic::{
    HarmonicPortControls, HarmonicPortResponse, HarmonicResponseControls, HarmonicResponseError,
    PreparedHarmonicResponse,
};
use fs_sparse::Csr;

/// A source-material refusal or a failure of the shared harmonic solve.
#[derive(Clone, Debug, PartialEq)]
pub enum ViscoHarmonicError {
    /// Original geometry, source frequency band or strain applicability failed.
    Element(ViscoElementError),
    /// Size, numerical payload, cancellation, singularity or residual refusal.
    Harmonic(HarmonicResponseError),
}
impl core::fmt::Display for ViscoHarmonicError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-SOLID-VISCO-HARMONIC: {self:?}")
    }
}
impl std::error::Error for ViscoHarmonicError {}
impl From<ViscoElementError> for ViscoHarmonicError {
    fn from(error: ViscoElementError) -> Self {
        Self::Element(error)
    }
}
impl From<HarmonicResponseError> for ViscoHarmonicError {
    fn from(error: HarmonicResponseError) -> Self {
        Self::Harmonic(error)
    }
}

/// Full, material-admitted periodic particular solution in free-DOF order.
#[derive(Clone, Debug)]
pub struct ViscoHarmonicResponse {
    /// Peak displacement phasors [m], under exp(-i omega t).
    pub displacement_m: Vec<C64>,
    /// Peak velocity phasors [m/s], equal to -i omega u.
    pub velocity_m_s: Vec<C64>,
    /// Maximum rowwise numerical backward error of the dynamic solve.
    /// This is not a forward-error bound near resonance.
    pub relative_residual: f64,
    /// Cycle-averaged applied work rate Re(f^H v)/2 [W].
    pub average_input_power_w: f64,
    /// Cycle-averaged original constitutive loss omega u^H K'' u / 2 [W].
    /// Includes every Prony branch and the original Kelvin viscosity once.
    pub average_material_heat_w: f64,
    /// Original boundary-port flow, opposing effort and absorbed work.
    /// Boundary work is kept separate from intrinsic material heat: its source
    /// determines whether it represents radiation, heat or another subsystem.
    pub port_load: Option<HarmonicPortResponse>,
    /// Input work minus material heat and boundary absorbed work [W]; not clipped.
    pub power_balance_residual_w: f64,
}

/// Immutable factor at one frequency, borrowed from the exact source solid.
/// Reuse it for distinct forcing phases, ports, or load cases at that frequency.
/// A new frequency requires a fresh material evaluation and preparation.
/// Dense factorization is O(n^3); retained storage and each solve are O(n^2),
/// followed by the original quadrature-point strain admission.
pub struct ViscoHarmonicOperator<'a> {
    source: &'a ViscoElements,
    omega_rad_s: f64,
    response: PreparedHarmonicResponse,
}

struct HarmonicBoundary<'a> {
    force_map: &'a [f64],
    ports: usize,
    impedance: &'a [C64],
    controls: HarmonicPortControls,
}

impl ViscoElements {
    /// Prepare the complete dispersive forced response, without truncating or
    /// diagonalizing the mass, storage or loss matrices. Source laws must admit
    /// this positive angular frequency. Stationary homogeneous supports use
    /// the original free-DOF map; moving supports and DC have separate owners.
    ///
    /// controls.maximum_payload_bytes bounds the harmonic owner's dense
    /// numerical payload PLUS the three dense source matrices used here. The
    /// existing sparse material assembly retains its original mesh/assembly
    /// limits; this is not a bound on its COO temporaries or the source solid.
    /// The coordinate and dense-payload gates run before spectral assembly.
    pub fn prepare_harmonic_response(
        &self,
        omega_rad_s: f64,
        controls: HarmonicResponseControls,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicOperator<'_>, ViscoHarmonicError> {
        self.prepare_harmonic_response_impl(omega_rad_s, controls, None, cx)
    }

    /// Prepare the original hereditary solid with a complete physical boundary
    /// impedance evaluated at this same frequency. The real row-major force map
    /// `B` has one row per original free DOF and one column per power-conjugate
    /// port. Its flow is `B^T v`, and the opposing effort is `Z B^T v`.
    ///
    /// The complete pencil is `K' - omega^2 M - i K'' - i omega B Z B^T`.
    /// All reactive and mutual entries are retained. Structure and load enter
    /// one factorization, so a physical boundary can regularize an otherwise
    /// singular bare resonance. The original material spectrum, full-cycle
    /// strain admission and intrinsic heat accounting remain unchanged.
    ///
    /// The shared owner admits the Hermitian part of `Z` as pointwise PSD.
    /// The caller owns the load's frequency domain, physical units and causality;
    /// a matrix at one frequency does not establish broadband passivity. The
    /// structural controls retain their existing source-matrix accounting;
    /// `port_controls` separately bound the additional projection and payload.
    #[allow(clippy::too_many_arguments)] // physical source, ports and separate explicit budgets
    pub fn prepare_harmonic_response_with_port_impedance(
        &self,
        omega_rad_s: f64,
        force_map: &[f64],
        ports: usize,
        impedance: &[C64],
        controls: HarmonicResponseControls,
        port_controls: HarmonicPortControls,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicOperator<'_>, ViscoHarmonicError> {
        self.prepare_harmonic_response_impl(
            omega_rad_s,
            controls,
            Some(HarmonicBoundary {
                force_map,
                ports,
                impedance,
                controls: port_controls,
            }),
            cx,
        )
    }

    fn prepare_harmonic_response_impl(
        &self,
        omega_rad_s: f64,
        controls: HarmonicResponseControls,
        boundary: Option<HarmonicBoundary<'_>>,
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicOperator<'_>, ViscoHarmonicError> {
        checkpoint(cx)?;
        let invalid = |what| ViscoHarmonicError::Harmonic(HarmonicResponseError::Invalid(what));
        if !omega_rad_s.is_finite() || omega_rad_s <= 0.0 {
            return Err(invalid(
                "material harmonic response requires a finite positive frequency",
            ));
        }
        let n = self.equilibrium.free_dofs.len();
        if n == 0
            || !controls.maximum_relative_residual.is_finite()
            || controls.maximum_relative_residual <= 0.0
            || controls.maximum_relative_residual >= 1.0
        {
            return Err(invalid(
                "nonempty pencil and residual budget in (0,1) required",
            ));
        }
        if n > controls.maximum_coordinates {
            return Err(HarmonicResponseError::Budget {
                quantity: "coordinates",
                requested: n,
                cap: controls.maximum_coordinates,
            }
            .into());
        }
        let source_bytes = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(24))
            .ok_or_else(|| invalid("harmonic source payload extent overflow"))?;
        // The shared owner admits 128*n*n + 128*n bytes. Include its complete
        // declared peak before allocating any of the additional dense inputs.
        let total = n
            .checked_mul(n)
            .and_then(|x| x.checked_mul(128))
            .and_then(|x| n.checked_mul(128).and_then(|y| x.checked_add(y)))
            .and_then(|x| x.checked_add(source_bytes))
            .ok_or_else(|| invalid("harmonic total payload extent overflow"))?;
        if total > controls.maximum_payload_bytes {
            return Err(HarmonicResponseError::Budget {
                quantity: "dense numerical payload bytes",
                requested: total,
                cap: controls.maximum_payload_bytes,
            }
            .into());
        }
        let stiffness = self.complex_stiffness(omega_rad_s, cx)?;
        let mass = dense(&self.equilibrium.mass, n, cx)?;
        let storage = dense(&stiffness.real, n, cx)?;
        let mut damping = dense(&stiffness.imag, n, cx)?;
        for value in &mut damping {
            let loss = *value;
            *value /= omega_rad_s;
            if !value.is_finite() || loss != 0.0 && *value == 0.0 {
                return Err(invalid(
                    "effective harmonic material loss is not representable",
                ));
            }
        }
        let solve_controls = HarmonicResponseControls {
            maximum_payload_bytes: controls.maximum_payload_bytes - source_bytes,
            ..controls
        };
        let response = if let Some(boundary) = boundary {
            PreparedHarmonicResponse::try_with_port_impedance(
                &mass,
                &damping,
                &storage,
                n,
                boundary.force_map,
                boundary.ports,
                boundary.impedance,
                omega_rad_s,
                solve_controls,
                boundary.controls,
                cx,
            )?
        } else {
            PreparedHarmonicResponse::try_new(
                &mass,
                &damping,
                &storage,
                n,
                omega_rad_s,
                solve_controls,
                cx,
            )?
        };
        checkpoint(cx)?;
        Ok(ViscoHarmonicOperator {
            source: self,
            omega_rad_s,
            response,
        })
    }
}

impl ViscoHarmonicOperator<'_> {
    /// Original full node-major DOF indices corresponding to the response rows.
    pub fn free_dofs(&self) -> &[usize] {
        &self.source.equilibrium.free_dofs
    }

    /// Requested physical angular frequency [rad/s].
    pub fn omega_rad_s(&self) -> f64 {
        self.omega_rad_s
    }

    /// Solve one peak force phasor [N] per original free coordinate. The cached
    /// factor is reused. Before returning, admit the maximum Mandel strain norm
    /// over the ENTIRE sinusoidal cycle at every original material quadrature
    /// point; checking only real/imaginary endpoint phases is insufficient.
    ///
    /// This is a finite-dimensional numerical particular solution. It is not
    /// a startup transient, modal/spatial convergence certificate, unsampled
    /// strain bound, or experimental validation of the source material.
    /// Cancellation or a source-domain refusal returns no partial response.
    pub fn solve(
        &self,
        force_n: &[C64],
        cx: &Cx<'_>,
    ) -> Result<ViscoHarmonicResponse, ViscoHarmonicError> {
        let result = self.response.solve(force_n, cx)?;
        self.source
            .validate_harmonic_displacement(&result.displacement, cx)?;
        checkpoint(cx)?;
        Ok(ViscoHarmonicResponse {
            displacement_m: result.displacement,
            velocity_m_s: result.velocity,
            relative_residual: result.relative_residual,
            average_input_power_w: result.average_input_power_w,
            average_material_heat_w: result.average_viscous_power_w,
            port_load: result.port_load,
            power_balance_residual_w: result.power_balance_residual_w,
        })
    }
}

impl ViscoElements {
    // One original full-cycle material admission for condensed and complete
    // memory-state forced response. Finite amplitudes cannot use pole-only
    // normalization as a substitute for this physical strain check.
    pub(super) fn validate_harmonic_displacement(
        &self,
        displacement: &[C64],
        cx: &Cx<'_>,
    ) -> Result<(), ViscoElementError> {
        for element in &self.elements {
            checkpoint(cx)?;
            for point in &element.points {
                checkpoint(cx)?;
                let mut re = [0.0; 6];
                let mut im = [0.0; 6];
                for (b, &dof) in point.b.iter().zip(&element.dofs) {
                    if let Ok(row) = self.equilibrium.free_dofs.binary_search(&dof) {
                        let u = displacement[row];
                        for i in 0..6 {
                            re[i] = b[i].mul_add(u.re, re[i]);
                            im[i] = b[i].mul_add(u.im, im[i]);
                        }
                    }
                }
                // Maximize ||re*cos(theta)+im*sin(theta)|| via the largest
                // eigenvector of its 2x2 Gram matrix. Scaling prevents overflow
                // in the phase calculation and preserves very small strains.
                let scale = re
                    .iter()
                    .chain(&im)
                    .map(|v| v.abs())
                    .fold(0.0_f64, f64::max);
                if re.iter().chain(&im).any(|v| !v.is_finite()) {
                    return Err(ViscoElementError::NonFiniteResult);
                }
                if scale == 0.0 {
                    continue;
                }
                let (mut a, mut b, mut c) = (0.0, 0.0, 0.0);
                for i in 0..6 {
                    let (r, s) = (re[i] / scale, im[i] / scale);
                    a += r * r;
                    b += s * s;
                    c += r * s;
                }
                let angle = 0.5 * det::atan2(2.0 * c, a - b);
                let (sin, cos) = det::sin_cos(angle);
                let peak = core::array::from_fn(|i| cos.mul_add(re[i], sin * im[i]));
                element
                    .law
                    .validate_mandel_strain(&peak)
                    .map_err(ViscoElementError::Material)?;
            }
        }
        checkpoint(cx)?;
        Ok(())
    }
}

fn dense(matrix: &Csr, n: usize, cx: &Cx<'_>) -> Result<Vec<f64>, ViscoElementError> {
    let mut values = vec![0.0; n * n];
    for i in 0..n {
        checkpoint(cx)?;
        let (columns, entries) = matrix.row(i);
        for (&j, &value) in columns.iter().zip(entries) {
            values[i * n + j] = value;
        }
    }
    Ok(values)
}
