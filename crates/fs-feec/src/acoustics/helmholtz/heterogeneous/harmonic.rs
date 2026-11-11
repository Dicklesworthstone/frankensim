//! Prepared full-space pressure solves and their original physical work.

use super::super::pressure_release::PressureReleaseConstraints;
use super::super::tri_quadrature;
use super::source::PressureSourceOwner;
use super::{
    BoundaryVelocitySource, HeterogeneousError, HeterogeneousHelmholtzFe,
    SpatialBoundaryAdmittance, VolumeVelocitySource, finite, invalid, poll,
};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeSet;

mod dispersive;
pub use dispersive::BoundaryAdmittanceFrequencySample;
mod equivalent;
pub use equivalent::{
    EquivalentFluidCoefficients, EquivalentFluidControls, EquivalentFluidElementPower,
    EquivalentFluidPressureField, PreparedEquivalentFluid, PreparedEquivalentFluidTraction,
};
mod frequency;
mod mechanical;
pub use mechanical::{
    AntiresonanceControls, AntiresonanceRecovery, DrivenAntiresonance, DrivenResonance,
    MechanicalHarmonicField, MechanicalHarmonicFrequencyResponse, MechanicalHarmonicLimits,
    MechanicalPoleError, MechanicalPoleLimits, MechanicalPoleOptions,
    MechanicalPoleResponseOptions, MechanicalPoleWall, MechanicalPressureObservation,
    MechanicalRelaxationParameterDirection, MechanicalResonanceDrive,
    MechanicalResonanceObservation, MechanicalResonanceParameterDerivative,
    MechanicalResonancePole, MechanicalResonanceRegion, MechanicalResonanceResidue,
    MechanicalResonanceResponseError, MechanicalResonanceSpectrum, MechanicalSpectrumOptions,
    MechanicalWallAccounting, MechanicalWallFrequencyDerivative, MechanicalWallParameterDerivative,
    MechanicalWallParameterDirection, MechanicalWallPhasors, MechanicalWallPower,
    MechanicalWallResidue, PreparedMechanicalHarmonic, PreparedMechanicalResonance,
    ResonanceControls, ResonanceSearchError,
};
mod resonance;
pub use resonance::{ResonanceSurvey, ResonanceSurveyControls, ResonanceSurveyError};
mod source_load;

impl From<frequency::PressureError> for HeterogeneousError {
    fn from(error: frequency::PressureError) -> Self {
        match error {
            frequency::PressureError::Invalid(what) => invalid(what).into(),
            frequency::PressureError::Factor(error) => error.into(),
            frequency::PressureError::Residual { observed, limit } => {
                Self::Residual { observed, limit }
            }
            frequency::PressureError::Limit {
                what,
                requested,
                maximum,
            } => Self::Limit {
                what,
                requested,
                maximum,
            },
            frequency::PressureError::Cancelled => Self::Cancelled,
        }
    }
}

/// Hard limit for the original dense full-pressure factorization.
pub const MAX_HARMONIC_DOFS: usize = 1536;

/// Passive, locally reacting boundary: outward velocity equals `Y p`.
///
/// `Y` is a SPECIFIC admittance [m/(Pa s)], not a volume-flow admittance.
/// Faces must be unique across all groups. An imposed inward velocity may
/// occupy the same face, giving `v_out = Y p - v_imposed` there.
#[derive(Debug, Clone)]
pub struct BoundaryAdmittance {
    /// Original boundary triangles; vertex winding does not change the physics.
    pub faces: Vec<[u32; 3]>,
    /// Finite complex `Y` in `exp(+i omega t)`, with `Re(Y) >= 0`.
    pub admittance_m_pa_s: C64,
}

/// Original face and complex inward normal velocity [m/s] at its three vertices.
/// The values follow the supplied vertex order. This is a piecewise-linear
/// trace, independently refinable from the hierarchical volume pressure.
pub type VelocityTrace = ([u32; 3], [C64; 3]);

/// Explicit limits on full-space harmonic preparation.
#[derive(Debug, Clone, Copy)]
pub struct HarmonicLimits {
    /// Unconstrained pressure DOF cap, also bounded by [`MAX_HARMONIC_DOFS`].
    pub max_dofs: usize,
    /// Peak scalar dense matrices: input plus factor, `32 * ndof^2` bytes.
    /// Sparse operators, O(ndof) vectors, allocator/runtime and an explicitly
    /// configured parallel LU worker's workspace are additional.
    pub max_dense_bytes: usize,
    /// Maximum componentwise backward error of the original physical equations.
    pub residual_tolerance: f64,
}

impl Default for HarmonicLimits {
    fn default() -> Self {
        Self {
            max_dofs: 1024,
            max_dense_bytes: 64 * 1024 * 1024,
            residual_tolerance: 5e-13,
        }
    }
}

/// One factorization of the heterogeneous pressure space at one frequency.
/// Exact homogeneous pressure-release constraints are eliminated; all other
/// original pressure coordinates remain, without a modal reduction.
/// No original medium is queried while solving subsequent right-hand sides.
pub struct PreparedHeterogeneous<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    frequency_hz: f64,
    omega: f64,
    factor: LuComplex,
    scales: Vec<f64>,
    conductance: Csr,
    susceptance: Csr,
    boundary_frequency_derivative: Option<(Csr, Csr)>,
    pressure_release: PressureReleaseConstraints,
    residual_tolerance: f64,
}

/// Full pressure coefficients and time-averaged physical work/energy.
#[derive(Debug, Clone)]
pub struct HeterogeneousPressureField {
    /// Coefficients [Pa] in the ORIGINAL hierarchical H1 basis. Only the vertex
    /// coordinates are point values; use the model's `pressure_at` for others.
    pub pressure_pa: Vec<C64>,
    /// Total boundary and volume-source work [W], with peak phasor factor 1/2.
    pub input_power_w: f64,
    /// Signed work of prepared quadratic boundary sources [W], already included
    /// in input power. The unprepared linear traces are not part of this subtotal.
    /// Coherent sources may individually absorb work.
    pub boundary_source_power_w: f64,
    /// Signed `Re(integral(conj(p) s))/2` [W], already included in input power.
    /// Coherent sources may individually absorb work; this is not bulk heat.
    pub volume_source_power_w: f64,
    /// Passive boundary heat `integral(Re(Y) |p|^2)/2` [W]. Bulk heat is zero.
    pub boundary_dissipation_w: f64,
    /// Time-averaged kinetic energy `p^H K_rho p/(4 omega^2)` [J].
    pub kinetic_energy_j: f64,
    /// Time-averaged compression energy `p^H M_beta p/4` [J].
    pub compressional_energy_j: f64,
    /// Physical componentwise backward error. This is an algebraic solve
    /// diagnostic, not a bound on mesh, source, geometry or quadrature error.
    pub backward_error: f64,
}

/// Original pressure response and its analytic frequency tangent.
#[derive(Debug, Clone)]
pub struct HeterogeneousFrequencyResponse {
    /// Original physical response, including its ordinary work and energy.
    pub response: HeterogeneousPressureField,
    /// Derivative of the original pressure coefficients per rad/s [Pa s/rad].
    /// Pressure-release trace coefficients remain exactly zero. This is a
    /// pressure tangent, not another physical pressure field or energy state.
    pub pressure_derivative_pa_per_rad_s: Vec<C64>,
    /// Componentwise backward error of the original differentiated equations.
    /// The denominator retains separate source and operator-derivative terms.
    /// This is not a bound on resonance, geometry or spatial-discretization error.
    pub derivative_backward_error: f64,
}

impl HeterogeneousPressureField {
    /// Signed input-minus-boundary-loss residual [W]; never forced to zero.
    #[must_use]
    pub fn power_balance_residual_w(&self) -> f64 {
        self.input_power_w - self.boundary_dissipation_w
    }
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare `(K_rho - omega^2 M_beta + i omega C_Y)` in the full pressure
    /// space. The original complex LU uses symmetric diagonal equilibration;
    /// physical residual correction later uses the separate ORIGINAL sparse
    /// terms, retaining compression terms lost to rounding in `K - omega^2 M`.
    pub fn prepare<'m>(
        &'m self,
        frequency_hz: f64,
        admittances: &[BoundaryAdmittance],
        limits: HarmonicLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedHeterogeneous<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_with_pressure_release(frequency_hz, admittances, &[], limits, cancelled)
    }

    /// Prepare the original physical operator with exact `p = 0` on the
    /// supplied boundary faces. Vertex, edge and face trace coefficients are
    /// eliminated before factorization, then lifted as exact zeros on output.
    /// No large-admittance penalty, altered material or modal truncation is used.
    ///
    /// Remaining faces are rigid unless they carry a passive admittance or
    /// an imposed inward velocity. A pressure-release face cannot also own
    /// either of those laws. Adjacent faces may share constrained trace DOFs.
    /// Original volume sources and physical energy/work operators are retained.
    /// Residual checks cover free equations: eliminated rows contain boundary
    /// reactions, not failed equations. A zero-pressure reaction does no work.
    ///
    /// This ideal pressure-soft boundary is not an aperture radiation model:
    /// exterior loading, end corrections and radiated power require their
    /// corresponding coupled boundary owners. Empty faces retain `prepare`'s
    /// original arithmetic. Dense limits count only the free pressure space;
    /// original sparse matrices and lifted O(ndof) fields remain additional.
    pub fn prepare_with_pressure_release<'m>(
        &'m self,
        frequency_hz: f64,
        admittances: &[BoundaryAdmittance],
        pressure_release_faces: &[[u32; 3]],
        limits: HarmonicLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedHeterogeneous<'m, 'f, 'c>, HeterogeneousError> {
        self.prepare_with_spatial_admittances(
            frequency_hz,
            admittances,
            &[],
            pressure_release_faces,
            limits,
            cancelled,
        )
    }

    /// Prepare the SAME full-space operator with pointwise sampled wall laws.
    ///
    /// Uniform groups and spatial fields may coexist on disjoint faces. Each
    /// spatial field must belong to this exact material/geometry owner and
    /// have been sampled at this exact frequency. All admittance faces must
    /// be unique and disjoint from pressure-release faces; imposed velocity
    /// sources may share admittance faces as before.
    ///
    /// The original physical forms receive +i omega C_Y. Boundary heat remains
    /// p^H Re(C_Y) p / 2, retaining spatial wall/pressure correlation. No
    /// material callback runs during factorization or subsequent solves.
    /// Empty spatial fields preserve the existing uniform-law arithmetic.
    pub fn prepare_with_spatial_admittances<'m>(
        &'m self,
        frequency_hz: f64,
        admittances: &[BoundaryAdmittance],
        spatial_admittances: &[&SpatialBoundaryAdmittance<'_, '_, '_>],
        pressure_release_faces: &[[u32; 3]],
        limits: HarmonicLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedHeterogeneous<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        if !omega.is_finite()
            || omega <= 0.0
            || !(omega * omega).is_finite()
            || omega * omega == 0.0
        {
            return Err(invalid("positive representable harmonic frequency required").into());
        }
        if limits.max_dofs == 0
            || limits.max_dense_bytes == 0
            || !limits.residual_tolerance.is_finite()
            || limits.residual_tolerance <= 0.0
        {
            return Err(invalid("positive harmonic limits and residual tolerance required").into());
        }
        let cap = limits.max_dofs.min(MAX_HARMONIC_DOFS);
        // Preserve cheap admission for the original unconstrained path before
        // allocating an identity constraint map over a potentially large mesh.
        if pressure_release_faces.is_empty() && self.ndof() > cap {
            return Err(HeterogeneousError::Limit {
                what: "unconstrained pressure degrees of freedom",
                requested: self.ndof(),
                maximum: cap,
            });
        }
        let pressure_release = self.geometry.pressure_release_constraints(
            pressure_release_faces,
            EmptyFreePressure::Refuse,
            cancelled,
        )?;
        let n = pressure_release.free_dofs.len();
        if n > cap {
            return Err(HeterogeneousError::Limit {
                what: "unconstrained pressure degrees of freedom",
                requested: n,
                maximum: cap,
            });
        }
        let entries = n
            .checked_mul(n)
            .ok_or_else(|| invalid("dense pressure extent overflow"))?;
        let bytes = entries
            .checked_mul(2 * core::mem::size_of::<C64>())
            .ok_or_else(|| invalid("dense pressure byte count overflow"))?;
        if bytes > limits.max_dense_bytes {
            return Err(HeterogeneousError::Limit {
                what: "full pressure dense matrix bytes",
                requested: bytes,
                maximum: limits.max_dense_bytes,
            });
        }
        for field in spatial_admittances {
            poll(cancelled)?;
            field.validate(self, frequency_hz)?;
        }
        let mut seen = BTreeSet::new();
        for faces in admittances
            .iter()
            .map(|boundary| boundary.faces.as_slice())
            .chain(
                spatial_admittances
                    .iter()
                    .map(|field| field.faces.as_slice()),
            )
        {
            for &face in faces {
                poll(cancelled)?;
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(invalid("admittance boundary faces must not overlap").into());
                }
                if pressure_release.faces.contains(&key) {
                    return Err(
                        invalid("pressure-release and admittance faces must not overlap").into(),
                    );
                }
            }
        }
        let (mut conductance, mut susceptance) =
            self.admittance_matrices(admittances, cancelled)?;
        if !spatial_admittances.is_empty() {
            let mut real = Coo::new(self.ndof(), self.ndof());
            let mut imaginary = Coo::new(self.ndof(), self.ndof());
            for (re, im) in core::iter::once((&conductance, &susceptance)).chain(
                spatial_admittances
                    .iter()
                    .map(|field| (&field.conductance, &field.susceptance)),
            ) {
                for row in 0..self.ndof() {
                    poll(cancelled)?;
                    for (coo, matrix) in [(&mut real, re), (&mut imaginary, im)] {
                        let (cols, values) = matrix.row(row);
                        for (&col, &value) in cols.iter().zip(values) {
                            coo.push(row, col, value);
                        }
                    }
                }
            }
            conductance = real.assemble();
            susceptance = imaginary.assemble();
            for row in 0..self.ndof() {
                poll(cancelled)?;
                if conductance
                    .row(row)
                    .1
                    .iter()
                    .chain(susceptance.row(row).1)
                    .any(|value| !value.is_finite())
                {
                    return Err(invalid("boundary admittance assembly overflow").into());
                }
            }
        }
        let mut a = vec![C64::ZERO; entries];
        let mut scales = vec![0.0; n];
        for (row, &original_row) in pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let (cols, vals) = self.stiffness.row(original_row);
            for (&col, &value) in cols.iter().zip(vals) {
                let free_col = pressure_release.original_to_free[col];
                if free_col != usize::MAX {
                    a[row * n + free_col].re += value;
                }
            }
            let (cols, vals) = self.mass.row(original_row);
            for (&col, &value) in cols.iter().zip(vals) {
                let free_col = pressure_release.original_to_free[col];
                if free_col != usize::MAX {
                    a[row * n + free_col].re -= omega * omega * value;
                }
            }
            let (cols, vals) = susceptance.row(original_row);
            for (&col, &value) in cols.iter().zip(vals) {
                let free_col = pressure_release.original_to_free[col];
                if free_col != usize::MAX {
                    a[row * n + free_col].re -= omega * value;
                }
            }
            let (cols, vals) = conductance.row(original_row);
            for (&col, &value) in cols.iter().zip(vals) {
                let free_col = pressure_release.original_to_free[col];
                if free_col != usize::MAX {
                    a[row * n + free_col].im += omega * value;
                }
            }
            let values = &a[row * n..(row + 1) * n];
            if values.iter().any(|&z| !finite(z)) {
                return Err(invalid("harmonic pressure operator overflow").into());
            }
            let max = values.iter().map(|z| z.abs()).fold(0.0, f64::max);
            if max == 0.0 {
                return Err(LuComplexError::Singular { index: row }.into());
            }
            scales[row] = 1.0 / det::sqrt(max);
            if !scales[row].is_finite() || scales[row] <= 0.0 {
                return Err(invalid("unrepresentable pressure equilibration").into());
            }
        }
        for row in 0..n {
            poll(cancelled)?;
            for col in 0..n {
                a[row * n + col] = a[row * n + col].scale(scales[row]).scale(scales[col]);
            }
        }
        if a.iter().any(|&z| !finite(z)) {
            return Err(invalid("pressure equilibration overflow").into());
        }
        let factor = lu_complex_with_cancel(&a, n, &mut *cancelled)?;
        poll(cancelled)?;
        Ok(PreparedHeterogeneous {
            model: self,
            frequency_hz,
            omega,
            factor,
            scales,
            conductance,
            susceptance,
            boundary_frequency_derivative: None,
            pressure_release,
            residual_tolerance: limits.residual_tolerance,
        })
    }

    pub(super) fn admittance_matrices(
        &self,
        admittances: &[BoundaryAdmittance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Csr, Csr), HeterogeneousError> {
        self.boundary_coefficient_matrices(
            admittances
                .iter()
                .map(|boundary| (boundary.faces.as_slice(), boundary.admittance_m_pa_s)),
            true,
            cancelled,
        )
    }

    fn boundary_coefficient_matrices<'a>(
        &self,
        boundaries: impl IntoIterator<Item = (&'a [[u32; 3]], C64)>,
        require_passive: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Csr, Csr), HeterogeneousError> {
        Self::boundary_coefficient_matrices_for(
            self.geometry,
            boundaries,
            require_passive,
            cancelled,
        )
    }

    fn boundary_coefficient_matrices_for<'a>(
        geometry: &super::super::HelmholtzFe<'_>,
        boundaries: impl IntoIterator<Item = (&'a [[u32; 3]], C64)>,
        require_passive: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Csr, Csr), HeterogeneousError> {
        let mut real = Coo::new(geometry.ndof(), geometry.ndof());
        let mut imaginary = Coo::new(geometry.ndof(), geometry.ndof());
        let quad = tri_quadrature(geometry.face_quadrature_order(geometry.order + 1));
        let mut seen = BTreeSet::new();
        for (faces, y) in boundaries {
            if !finite(y) || (require_passive && y.re < 0.0) || faces.is_empty() {
                return Err(
                    invalid("finite passive admittance on nonempty boundary required").into(),
                );
            }
            for &face in faces {
                poll(cancelled)?;
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(invalid("admittance boundary faces must not overlap").into());
                }
                let rule = geometry.face_rule(face, &quad)?;
                for (i, &row) in rule.dofs.iter().enumerate() {
                    for (j, &col) in rule.dofs.iter().enumerate().skip(i) {
                        let value: f64 = rule
                            .weights
                            .iter()
                            .zip(rule.phi[i].iter().zip(&rule.phi[j]))
                            .map(|(w, (a, b))| w * a * b)
                            .sum();
                        for (coo, coefficient) in [(&mut real, y.re), (&mut imaginary, y.im)] {
                            if coefficient != 0.0 {
                                let entry = value * coefficient;
                                if !entry.is_finite() {
                                    return Err(
                                        invalid("boundary admittance assembly overflow").into()
                                    );
                                }
                                coo.push(row, col, entry);
                                if i != j {
                                    coo.push(col, row, entry);
                                }
                            }
                        }
                    }
                }
            }
        }
        poll(cancelled)?;
        Ok((real.assemble(), imaginary.assemble()))
    }

    fn velocity_load(
        &self,
        traces: &[VelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HeterogeneousError> {
        Self::velocity_load_for(self.geometry, traces, cancelled)
    }

    fn velocity_load_for(
        geometry: &super::super::HelmholtzFe<'_>,
        traces: &[VelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HeterogeneousError> {
        let quad = tri_quadrature(geometry.face_quadrature_order(geometry.order + 1));
        let mut load = vec![C64::ZERO; geometry.ndof()];
        let mut seen = BTreeSet::new();
        for &(face, values) in traces {
            poll(cancelled)?;
            let mut key = face;
            key.sort_unstable();
            if !seen.insert(key) || values.iter().any(|&v| !finite(v)) {
                return Err(
                    invalid("velocity trace requires unique faces and finite values").into(),
                );
            }
            let rule = geometry.face_rule(face, &quad)?;
            let mut local_values = [C64::ZERO; 3];
            for (i, vertex) in rule.vertices.iter().enumerate() {
                let input = face
                    .iter()
                    .position(|v| v == vertex)
                    .ok_or_else(|| invalid("velocity trace vertex mismatch"))?;
                local_values[i] = values[input];
            }
            for (q, &(barycentric, _)) in quad.iter().enumerate() {
                let value = local_values
                    .iter()
                    .zip(barycentric)
                    .fold(C64::ZERO, |sum, (&v, b)| sum + v.scale(b));
                for (local, &global) in rule.dofs.iter().enumerate() {
                    load[global] = load[global] + value.scale(rule.weights[q] * rule.phi[local][q]);
                }
            }
        }
        if load.iter().any(|&v| !finite(v)) {
            return Err(invalid("velocity load projection overflow").into());
        }
        Ok(load)
    }
}

impl PreparedHeterogeneous<'_, '_, '_> {
    /// Frequency of this factorization [Hz].
    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }

    /// Solve for arbitrary simultaneous spatial velocity traces using the
    /// retained factor. Unlisted faces remain rigid unless they carry an
    /// admittance or were declared pressure-release during preparation.
    /// Sources and admittances can share a face; pressure release cannot.
    ///
    /// At most two residual corrections use the original separate physical
    /// sparse operators and the same LU. No source is resampled or refactored.
    pub fn solve(
        &self,
        inward_velocity: &[VelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousPressureField, HeterogeneousError> {
        poll(cancelled)?;
        self.admit_velocity_faces(inward_velocity, cancelled)?;
        let load = self.model.velocity_load(inward_velocity, cancelled)?;
        self.solve_projected(&load, cancelled)
    }

    /// Solve pressure and its analytic frequency derivative with the SAME LU.
    ///
    /// Derivatives are per angular frequency [rad/s], although preparation
    /// accepts Hz. Optional velocity traces are `dv_in/domega`, including
    /// coherent source phase changes. Missing derivative faces contribute zero;
    /// their order and face set need not match the primary traces.
    ///
    /// Geometry, density and compressibility remain frozen. Boundary slopes
    /// supplied through `prepare_with_dispersive_admittances` are included;
    /// other preparations hold their admittances fixed. A sampled complex or
    /// spatial wall law never implies an unstated dispersion derivative. No
    /// finite differences, modal truncation or fitted pole is used. Exact
    /// pressure-release traces and their tangents stay zero.
    pub fn solve_with_frequency_derivative(
        &self,
        inward_velocity: &[VelocityTrace],
        inward_velocity_derivative: Option<&[VelocityTrace]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousFrequencyResponse, HeterogeneousError> {
        self.solve_with_sources_frequency_derivative(
            inward_velocity,
            inward_velocity_derivative,
            &[],
            &[],
            cancelled,
        )
    }

    /// Analytic frequency response with every original source family.
    ///
    /// Each source triple is `(original_source, amplitude, d_amplitude/domega)`.
    /// Boundary amplitudes are inward velocity [m/s]; volume amplitudes are
    /// volume flow [m3/s]. Their derivatives are per rad/s. Both complete drives
    /// are projected coherently through the same original source shapes before
    /// solving; zero-amplitude sources still retain their supplied derivatives.
    /// Linear trace slopes follow [`Self::solve_with_frequency_derivative`].
    ///
    /// Source counts, model binding, face conflicts, finite values and the
    /// original residual tolerance apply to BOTH drives. The derivative reuses
    /// the retained factor and at most two corrections with separate original
    /// sparse physical terms in double-double arithmetic. No dense matrix is
    /// allocated or factored again. Cancellation publishes neither a partial
    /// primary response nor a partial derivative and leaves this factor reusable.
    pub fn solve_with_sources_frequency_derivative(
        &self,
        inward_velocity: &[VelocityTrace],
        inward_velocity_derivative: Option<&[VelocityTrace]>,
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousFrequencyResponse, HeterogeneousError> {
        poll(cancelled)?;
        super::boundary_source::admit_count(boundary_sources.len())?;
        super::source::admit_count(volume_sources.len())?;
        let (boundary_values, boundary_slopes): (Vec<_>, Vec<_>) = boundary_sources
            .iter()
            .map(|&(source, amplitude, slope)| ((source, amplitude), (source, slope)))
            .unzip();
        let (volume_values, volume_slopes): (Vec<_>, Vec<_>) = volume_sources
            .iter()
            .map(|&(source, amplitude, slope)| ((source, amplitude), (source, slope)))
            .unzip();
        let load =
            self.project_sources(inward_velocity, &boundary_values, &volume_values, cancelled)?;
        let load_derivative = self.project_sources(
            inward_velocity_derivative.unwrap_or(&[]),
            &boundary_slopes,
            &volume_slopes,
            cancelled,
        )?;
        let pencil = self.pressure_pencil();
        let primary = pencil.solve(&load, cancelled)?;
        let derivative =
            pencil.frequency_derivative(&primary.pressure, &load, &load_derivative, cancelled)?;
        let mut response =
            self.field(primary.pressure, &load, primary.backward_error, cancelled)?;
        Self::source_work(&mut response, &boundary_values, &volume_values, cancelled)?;
        Ok(HeterogeneousFrequencyResponse {
            response,
            pressure_derivative_pa_per_rad_s: derivative.pressure,
            derivative_backward_error: derivative.backward_error,
        })
    }

    /// Drive the SAME full-pressure factor with simultaneous boundary motion
    /// and original-volume sources. Each amplitude is a complex peak volume
    /// flow [m^3/s] in `exp(+i omega t)`. Sources interfere before observation
    /// and work evaluation; no source normalization or propagation is repeated.
    pub fn solve_with_volume_sources(
        &self,
        inward_velocity: &[VelocityTrace],
        sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousPressureField, HeterogeneousError> {
        self.solve_with_boundary_sources(inward_velocity, &[], sources, cancelled)
    }

    /// Drive the retained factor with linear boundary traces, prepared
    /// quadratic boundary motion and original-volume injection together.
    ///
    /// A prepared boundary amplitude is inward normal velocity [m/s]; a volume
    /// amplitude is volume flow [m^3/s]. Both are complex peak phasors in
    /// `exp(+i omega t)`. Their original load columns add before solving, so
    /// overlapping sources interfere without any area or shape normalization.
    /// Source preparation and material sampling are not repeated. Every source
    /// must belong to this model, even at zero amplitude; prepared faces cannot
    /// overlap a pressure-release face. Passive admittance may share a source
    /// face. Signed source work subtotals are already included in total input.
    pub fn solve_with_boundary_sources(
        &self,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousPressureField, HeterogeneousError> {
        let load =
            self.project_sources(inward_velocity, boundary_sources, volume_sources, cancelled)?;
        let mut field = self.solve_projected(&load, cancelled)?;
        Self::source_work(&mut field, boundary_sources, volume_sources, cancelled)?;
        Ok(field)
    }

    fn project_sources(
        &self,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HeterogeneousError> {
        Self::project_sources_for(
            self.model,
            &self.pressure_release,
            inward_velocity,
            boundary_sources,
            volume_sources,
            cancelled,
        )
    }

    pub(super) fn project_sources_for(
        model: &HeterogeneousHelmholtzFe<'_, '_>,
        pressure_release: &PressureReleaseConstraints,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HeterogeneousError> {
        Self::project_sources_for_owner(
            PressureSourceOwner::Heterogeneous(model),
            pressure_release,
            inward_velocity,
            boundary_sources,
            volume_sources,
            cancelled,
        )
    }

    pub(super) fn project_sources_for_owner(
        owner: PressureSourceOwner<'_, '_, '_>,
        pressure_release: &PressureReleaseConstraints,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HeterogeneousError> {
        poll(cancelled)?;
        Self::admit_velocity_faces_for(pressure_release, inward_velocity, cancelled)?;
        Self::project_sources_with_boundary_exclusions(
            owner,
            &pressure_release.faces,
            inward_velocity,
            boundary_sources,
            volume_sources,
            cancelled,
        )
    }

    /// Shared original source projection; exclusions carry the receiving
    /// owner's physical boundary laws without introducing another trace map.
    pub(in crate::acoustics::helmholtz) fn project_sources_with_boundary_exclusions(
        owner: PressureSourceOwner<'_, '_, '_>,
        excluded_faces: &BTreeSet<[u32; 3]>,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HeterogeneousError> {
        poll(cancelled)?;
        super::boundary_source::admit_count(boundary_sources.len())?;
        super::source::admit_count(volume_sources.len())?;
        for &(source, amplitude) in boundary_sources {
            poll(cancelled)?;
            source.validate_owner(owner)?;
            if !finite(amplitude) {
                return Err(invalid("finite harmonic boundary-velocity amplitude required").into());
            }
            if !excluded_faces.is_empty() {
                for face in &source.faces {
                    poll(cancelled)?;
                    if excluded_faces.contains(face) {
                        return Err(invalid(
                            "prepared boundary source overlaps a constrained boundary law",
                        )
                        .into());
                    }
                }
            }
        }
        for &(source, amplitude) in volume_sources {
            poll(cancelled)?;
            source.validate_owner(owner)?;
            if !finite(amplitude) {
                return Err(invalid("finite harmonic volume-flow amplitude required").into());
            }
        }
        let mut load = HeterogeneousHelmholtzFe::velocity_load_for(
            owner.geometry(),
            inward_velocity,
            cancelled,
        )?;
        for (row, value) in load.iter_mut().enumerate() {
            poll(cancelled)?;
            let mut real = Dd::from_f64(value.re);
            let mut imaginary = Dd::from_f64(value.im);
            for &(source, amplitude) in boundary_sources {
                real = real + Dd::from_f64(source.load[row]) * Dd::from_f64(amplitude.re);
                imaginary = imaginary + Dd::from_f64(source.load[row]) * Dd::from_f64(amplitude.im);
            }
            for &(source, amplitude) in volume_sources {
                real = real + Dd::from_f64(source.load[row]) * Dd::from_f64(amplitude.re);
                imaginary = imaginary + Dd::from_f64(source.load[row]) * Dd::from_f64(amplitude.im);
            }
            *value = C64::new(real.to_f64(), imaginary.to_f64());
        }
        Ok(load)
    }

    pub(super) fn source_work(
        field: &mut HeterogeneousPressureField,
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HeterogeneousError> {
        let (boundary, volume) = Self::source_power(
            &field.pressure_pa,
            boundary_sources,
            volume_sources,
            cancelled,
        )?;
        field.boundary_source_power_w = boundary;
        field.volume_source_power_w = volume;
        Ok(())
    }

    pub(in crate::acoustics::helmholtz) fn source_power(
        pressure_pa: &[C64],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, f64), HeterogeneousError> {
        let mut work = Dd::from_f64(0.0);
        for &(source, amplitude) in boundary_sources {
            poll(cancelled)?;
            let pressure = source.harmonic_pressure(pressure_pa);
            work = work
                + Dd::from_f64(pressure.re) * Dd::from_f64(amplitude.re)
                + Dd::from_f64(pressure.im) * Dd::from_f64(amplitude.im);
        }
        let boundary_source_power_w = 0.5 * work.to_f64();
        work = Dd::from_f64(0.0);
        for &(source, amplitude) in volume_sources {
            poll(cancelled)?;
            let pressure = source.harmonic_pressure(pressure_pa);
            work = work
                + Dd::from_f64(pressure.re) * Dd::from_f64(amplitude.re)
                + Dd::from_f64(pressure.im) * Dd::from_f64(amplitude.im);
        }
        let volume_source_power_w = 0.5 * work.to_f64();
        if !boundary_source_power_w.is_finite() || !volume_source_power_w.is_finite() {
            return Err(invalid("harmonic prepared-source work overflow").into());
        }
        poll(cancelled)?;
        Ok((boundary_source_power_w, volume_source_power_w))
    }

    fn admit_velocity_faces(
        &self,
        traces: &[VelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HeterogeneousError> {
        Self::admit_velocity_faces_for(&self.pressure_release, traces, cancelled)
    }

    fn admit_velocity_faces_for(
        pressure_release: &PressureReleaseConstraints,
        traces: &[VelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HeterogeneousError> {
        if pressure_release.faces.is_empty() {
            return Ok(());
        }
        for &(face, _) in traces {
            poll(cancelled)?;
            let mut key = face;
            key.sort_unstable();
            if pressure_release.faces.contains(&key) {
                return Err(invalid(
                    "pressure-release and imposed-velocity faces must not overlap",
                )
                .into());
            }
        }
        Ok(())
    }

    fn solve_projected(
        &self,
        load: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousPressureField, HeterogeneousError> {
        let solved = self.pressure_pencil().solve(load, cancelled)?;
        self.field(solved.pressure, load, solved.backward_error, cancelled)
    }

    fn pressure_pencil(&self) -> frequency::PressurePencil<'_> {
        frequency::PressurePencil {
            omega: self.omega,
            stiffness: &self.model.stiffness,
            mass: &self.model.mass,
            conductance: &self.conductance,
            susceptance: &self.susceptance,
            boundary_frequency_derivative: self
                .boundary_frequency_derivative
                .as_ref()
                .map(|(g, b)| (g, b)),
            factor: &self.factor,
            scales: &self.scales,
            free_dofs: &self.pressure_release.free_dofs,
            residual_tolerance: self.residual_tolerance,
        }
    }

    fn field(
        &self,
        pressure_pa: Vec<C64>,
        load: &[C64],
        backward_error: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousPressureField, HeterogeneousError> {
        Self::field_for(
            self.model,
            self.omega,
            &self.conductance,
            pressure_pa,
            load,
            backward_error,
            cancelled,
        )
    }

    pub(super) fn field_for(
        model: &HeterogeneousHelmholtzFe<'_, '_>,
        omega: f64,
        conductance: &Csr,
        pressure_pa: Vec<C64>,
        load: &[C64],
        backward_error: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HeterogeneousPressureField, HeterogeneousError> {
        let mut work = Dd::from_f64(0.0);
        for (p, v) in pressure_pa.iter().zip(load) {
            work = work
                + Dd::from_f64(p.re) * Dd::from_f64(v.re)
                + Dd::from_f64(p.im) * Dd::from_f64(v.im);
        }
        let input_power_w = 0.5 * work.to_f64();
        let boundary_dissipation_w =
            0.5 * quadratic_nonnegative(conductance, &pressure_pa, cancelled)?;
        let kinetic_energy_j = 0.25
            * quadratic_nonnegative(&model.stiffness, &pressure_pa, cancelled)?
            / omega
            / omega;
        let compressional_energy_j =
            0.25 * quadratic_nonnegative(&model.mass, &pressure_pa, cancelled)?;
        if [
            input_power_w,
            boundary_dissipation_w,
            kinetic_energy_j,
            compressional_energy_j,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(invalid("unrepresentable acoustic energy or work").into());
        }
        poll(cancelled)?;
        Ok(HeterogeneousPressureField {
            pressure_pa,
            input_power_w,
            boundary_source_power_w: 0.0,
            volume_source_power_w: 0.0,
            boundary_dissipation_w,
            kinetic_energy_j,
            compressional_energy_j,
            backward_error,
        })
    }
}

// All three forms are positive semidefinite before floating-point assembly.
// Remove only a negative rounding residue bounded by its absolute summands;
// a materially negative form refuses instead of inventing positive energy.
fn quadratic_nonnegative(
    matrix: &Csr,
    field: &[C64],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, HeterogeneousError> {
    frequency::quadratic_nonnegative(matrix, field, cancelled).map_err(Into::into)
}

#[cfg(test)]
#[path = "harmonic/pressure_release_tests.rs"]
mod pressure_release_tests;
