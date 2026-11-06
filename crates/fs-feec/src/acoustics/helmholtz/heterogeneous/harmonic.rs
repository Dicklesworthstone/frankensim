//! Prepared full-space pressure solves and their original physical work.

use super::super::pressure_release::PressureReleaseConstraints;
use super::super::tri_quadrature;
use super::{
    BoundaryVelocitySource, HeterogeneousError, HeterogeneousHelmholtzFe, SpatialBoundaryAdmittance,
    VolumeVelocitySource, finite, invalid, poll,
};
use fs_la::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeSet;

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
        let pressure_release = self
            .geometry
            .pressure_release_constraints(pressure_release_faces, cancelled)?;
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
            .chain(spatial_admittances.iter().map(|field| field.faces.as_slice()))
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
            pressure_release,
            residual_tolerance: limits.residual_tolerance,
        })
    }

    pub(super) fn admittance_matrices(
        &self,
        admittances: &[BoundaryAdmittance],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Csr, Csr), HeterogeneousError> {
        let mut real = Coo::new(self.ndof(), self.ndof());
        let mut imaginary = Coo::new(self.ndof(), self.ndof());
        let quad = tri_quadrature(self.geometry.face_quadrature_order(self.geometry.order + 1));
        let mut seen = BTreeSet::new();
        for boundary in admittances {
            let y = boundary.admittance_m_pa_s;
            if !finite(y) || y.re < 0.0 || boundary.faces.is_empty() {
                return Err(
                    invalid("finite passive admittance on nonempty boundary required").into(),
                );
            }
            for &face in &boundary.faces {
                poll(cancelled)?;
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(invalid("admittance boundary faces must not overlap").into());
                }
                let rule = self.geometry.face_rule(face, &quad)?;
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
        let quad = tri_quadrature(self.geometry.face_quadrature_order(self.geometry.order + 1));
        let mut load = vec![C64::ZERO; self.ndof()];
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
            let rule = self.geometry.face_rule(face, &quad)?;
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
        poll(cancelled)?;
        self.admit_velocity_faces(inward_velocity, cancelled)?;
        super::boundary_source::admit_count(boundary_sources.len())?;
        super::source::admit_count(volume_sources.len())?;
        for &(source, amplitude) in boundary_sources {
            poll(cancelled)?;
            source.validate(self.model)?;
            if !finite(amplitude) {
                return Err(invalid("finite harmonic boundary-velocity amplitude required").into());
            }
            if !self.pressure_release.faces.is_empty() {
                for face in &source.faces {
                    poll(cancelled)?;
                    if self.pressure_release.faces.contains(face) {
                        return Err(invalid(
                            "pressure-release and prepared boundary-source faces must not overlap",
                        )
                        .into());
                    }
                }
            }
        }
        for &(source, amplitude) in volume_sources {
            poll(cancelled)?;
            source.validate(self.model)?;
            if !finite(amplitude) {
                return Err(invalid("finite harmonic volume-flow amplitude required").into());
            }
        }
        let mut load = self.model.velocity_load(inward_velocity, cancelled)?;
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
        let mut field = self.solve_projected(&load, cancelled)?;
        let mut work = Dd::from_f64(0.0);
        for &(source, amplitude) in boundary_sources {
            poll(cancelled)?;
            let pressure = source.harmonic_pressure(&field.pressure_pa);
            work = work
                + Dd::from_f64(pressure.re) * Dd::from_f64(amplitude.re)
                + Dd::from_f64(pressure.im) * Dd::from_f64(amplitude.im);
        }
        field.boundary_source_power_w = 0.5 * work.to_f64();
        work = Dd::from_f64(0.0);
        for &(source, amplitude) in volume_sources {
            poll(cancelled)?;
            let pressure = source.harmonic_pressure(&field.pressure_pa);
            work = work
                + Dd::from_f64(pressure.re) * Dd::from_f64(amplitude.re)
                + Dd::from_f64(pressure.im) * Dd::from_f64(amplitude.im);
        }
        field.volume_source_power_w = 0.5 * work.to_f64();
        if !field.boundary_source_power_w.is_finite() || !field.volume_source_power_w.is_finite() {
            return Err(invalid("harmonic prepared-source work overflow").into());
        }
        poll(cancelled)?;
        Ok(field)
    }

    fn admit_velocity_faces(
        &self,
        traces: &[VelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), HeterogeneousError> {
        if self.pressure_release.faces.is_empty() {
            return Ok(());
        }
        for &(face, _) in traces {
            poll(cancelled)?;
            let mut key = face;
            key.sort_unstable();
            if self.pressure_release.faces.contains(&key) {
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
        let rhs = load
            .iter()
            .map(|v| C64::new(-self.omega * v.im, self.omega * v.re))
            .collect::<Vec<_>>();
        if rhs.iter().any(|&z| !finite(z)) {
            return Err(invalid("harmonic velocity source overflow").into());
        }
        let mut free_pressure = self
            .pressure_release
            .free_dofs
            .iter()
            .zip(&self.scales)
            .map(|(&original, &s)| rhs[original].scale(s))
            .collect::<Vec<_>>();
        self.factor
            .solve_with_cancel(&mut free_pressure, &mut *cancelled)?;
        let mut pressure = vec![C64::ZERO; self.model.ndof()];
        for ((&original, &p), &scale) in self
            .pressure_release
            .free_dofs
            .iter()
            .zip(&free_pressure)
            .zip(&self.scales)
        {
            pressure[original] = p.scale(scale);
        }
        for correction in 0..=2 {
            poll(cancelled)?;
            if pressure.iter().any(|&z| !finite(z)) {
                return Err(invalid("nonfinite harmonic pressure").into());
            }
            let (mut residual, backward_error) = self.residual(&pressure, &rhs, cancelled)?;
            if correction == 2 || backward_error < 1e-18 {
                if backward_error > self.residual_tolerance {
                    return Err(HeterogeneousError::Residual {
                        observed: backward_error,
                        limit: self.residual_tolerance,
                    });
                }
                return self.field(pressure, load, backward_error, cancelled);
            }
            for (r, &scale) in residual.iter_mut().zip(&self.scales) {
                *r = r.scale(scale);
            }
            self.factor
                .solve_with_cancel(&mut residual, &mut *cancelled)?;
            for (i, &original) in self.pressure_release.free_dofs.iter().enumerate() {
                pressure[original] = pressure[original] + residual[i].scale(self.scales[i]);
            }
        }
        unreachable!("the bounded correction loop returns on its last residual")
    }

    fn residual(
        &self,
        pressure: &[C64],
        rhs: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<C64>, f64), HeterogeneousError> {
        let omega = Dd::from_f64(self.omega);
        let omega2 = omega * omega;
        let mut residual = vec![C64::ZERO; self.pressure_release.free_dofs.len()];
        let mut worst: f64 = 0.0;
        for (free_row, &row) in self.pressure_release.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let mut real = Dd::from_f64(rhs[row].re);
            let mut imaginary = Dd::from_f64(rhs[row].im);
            let mut scale = rhs[row].abs();
            for (matrix, multiplier, rotate) in [
                (&self.model.stiffness, Dd::from_f64(1.0), false),
                (&self.model.mass, -omega2, false),
                (&self.susceptance, -omega, false),
                (&self.conductance, omega, true),
            ] {
                let (cols, vals) = matrix.row(row);
                for (&col, &value) in cols.iter().zip(vals) {
                    let coefficient = multiplier * Dd::from_f64(value);
                    let p = pressure[col];
                    if rotate {
                        real = real + coefficient * Dd::from_f64(p.im);
                        imaginary = imaginary - coefficient * Dd::from_f64(p.re);
                    } else {
                        real = real - coefficient * Dd::from_f64(p.re);
                        imaginary = imaginary - coefficient * Dd::from_f64(p.im);
                    }
                    scale += coefficient.to_f64().abs() * p.abs();
                }
            }
            let value = C64::new(real.to_f64(), imaginary.to_f64());
            if !finite(value) || !scale.is_finite() {
                return Err(invalid("unrepresentable physical pressure residual").into());
            }
            worst = worst.max(value.abs() / scale.max(f64::MIN_POSITIVE));
            residual[free_row] = value;
        }
        Ok((residual, worst))
    }

    fn field(
        &self,
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
            0.5 * quadratic_nonnegative(&self.conductance, &pressure_pa, cancelled)?;
        let kinetic_energy_j = 0.25
            * quadratic_nonnegative(&self.model.stiffness, &pressure_pa, cancelled)?
            / self.omega
            / self.omega;
        let compressional_energy_j =
            0.25 * quadratic_nonnegative(&self.model.mass, &pressure_pa, cancelled)?;
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
    let mut sum = Dd::from_f64(0.0);
    let mut absolute = 0.0;
    for (row, p) in field.iter().enumerate() {
        poll(cancelled)?;
        let (cols, values) = matrix.row(row);
        for (&col, &coefficient) in cols.iter().zip(values) {
            let q = field[col];
            sum = sum
                + Dd::from_f64(coefficient)
                    * (Dd::from_f64(p.re) * Dd::from_f64(q.re)
                        + Dd::from_f64(p.im) * Dd::from_f64(q.im));
            absolute += coefficient.abs() * p.abs() * q.abs();
        }
    }
    let value = sum.to_f64();
    if !value.is_finite() || !absolute.is_finite() || value < -64.0 * f64::EPSILON * absolute {
        return Err(invalid("nonpositive or unrepresentable acoustic energy form").into());
    }
    Ok(value.max(0.0))
}

#[cfg(test)]
#[path = "harmonic/pressure_release_tests.rs"]
mod pressure_release_tests;
