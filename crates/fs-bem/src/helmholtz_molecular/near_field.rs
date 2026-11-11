//! Source-bound adaptive observation of the original molecular boundary field.
use super::*;
use crate::near_field::{Geometry, VectorObservation, complex as receiver};

/// Accuracy/work and retained operator-row storage for molecular receivers.
/// The original boundary discretization and solve retain their own budgets.
#[derive(Debug, Clone, Copy)]
pub struct MolecularNearFieldOptions {
    /// Shared adaptive triangle quadrature, including separate gradient gates.
    pub quadrature: crate::near_field::Options,
    /// Retained row payload and row descriptors. Bounded receiver positions/output buffers,
    /// geometry admission scratch, allocator metadata and stacks are excluded.
    /// Receiver cardinality remains bounded by the original geometry owner.
    pub maximum_row_bytes: usize,
}

impl Default for MolecularNearFieldOptions {
    fn default() -> Self {
        Self {
            quadrature: crate::near_field::Options::default(),
            maximum_row_bytes: 16 * 1024 * 1024,
        }
    }
}

impl MolecularNearFieldOptions {
    pub(crate) fn admit(self, panels: usize, receivers: usize) -> Result<(), HelmholtzError> {
        self.quadrature.validate()?;
        if panels
            .checked_mul(receivers)
            .and_then(|pairs| pairs.checked_mul(80))
            .is_none_or(|work| work > self.quadrature.maximum_kernel_evaluations)
        {
            return Err(HelmholtzError::BadParameter {
                what: "molecular receivers exceed the minimum adaptive quadrature work",
            });
        }
        let required_bytes = receiver::required_row_bytes(panels, receivers).ok_or(
            HelmholtzError::BadParameter {
                what: "molecular receiver row storage overflow",
            },
        )?;
        if required_bytes > self.maximum_row_bytes {
            return Err(HelmholtzError::CacheBudget {
                required_bytes,
                maximum_bytes: self.maximum_row_bytes,
            });
        }
        Ok(())
    }
}

/// One reusable pressure/gradient operator bound to the original molecular
/// source, frequency and boundary geometry. No boundary solve is performed.
pub struct PreparedMolecularNearField<'a> {
    source: Arc<RelaxationState>,
    omega: f64,
    k: C64,
    fingerprint: u64,
    points: Vec<[f64; 3]>,
    rows: receiver::Prepared<'a>,
}

impl MolecularRadiation {
    /// Integrate the same constant pressure/actual-velocity boundary traces on
    /// retained triangles. Pressure and every analytic gradient component use
    /// their own adaptive discrepancy test, with common geometric separation,
    /// complex-wave variation, depth and complete-array quadrature-work gates.
    ///
    /// Reuse the returned rows across the original source solution bank. The
    /// source's real rho*omega supplies Euler momentum independently of its
    /// dispersive complex k. Matching rho/c alone does not authorize reuse with
    /// a different source owner. The ordinary geometry owner admits exact
    /// winding/clearance; cancellation reaches adaptive panels and triangle
    /// tiles. A failure never substitutes centroid observations.
    pub fn prepare_near_field<'a>(
        &self,
        surface: &'a SpherePanels,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: MolecularNearFieldOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedMolecularNearField<'a>, MolecularError> {
        poll(&mut cancelled)?;
        if self.fingerprint != surface_fingerprint(surface) {
            return Err(bad(
                "molecular receiver belongs to a different source surface",
            ));
        }
        options.admit(surface.areas().len(), points.len())?;
        let geometry = Geometry::new(surface, points, minimum_clearance_m)?;
        let rows = receiver::Prepared::new(
            &geometry,
            DampedWaveNumber {
                k: self.k,
                phase_origin_m: 0.0,
            },
            self.source.gas().density * self.omega,
            options.quadrature,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedMolecularNearField {
            source: Arc::clone(&self.source),
            omega: self.omega,
            k: self.k,
            fingerprint: self.fingerprint,
            points: points.to_vec(),
            rows,
        })
    }

    /// Paired adaptive pressure/particle velocity at the original receivers.
    /// Estimates describe quadrature of the fixed solved traces only; they do
    /// not bound boundary discretization, source error or relative error near
    /// a coherent cancellation. Use prepare_near_field to amortize quadrature
    /// across multiple solutions from the same original source owner.
    pub fn observe_near_field(
        &self,
        surface: &SpherePanels,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: MolecularNearFieldOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VectorObservation, MolecularError> {
        self.prepare_near_field(
            surface,
            points,
            minimum_clearance_m,
            options,
            &mut cancelled,
        )?
        .evaluate(self, cancelled)
    }
}

impl PreparedMolecularNearField<'_> {
    /// Original physical receiver positions, in the prepared row order.
    pub fn points_m(&self) -> &[[f64; 3]] {
        &self.points
    }

    /// Actual preparation work, shared by all subsequent source columns.
    pub fn kernel_evaluations(&self) -> usize {
        self.rows.kernel_evaluations()
    }

    fn admit_solution(&self, solution: &MolecularRadiation) -> Result<(), MolecularError> {
        if !Arc::ptr_eq(&self.source, &solution.source)
            || self.omega.to_bits() != solution.omega.to_bits()
            || self.k.re.to_bits() != solution.k.re.to_bits()
            || self.k.im.to_bits() != solution.k.im.to_bits()
            || self.fingerprint != solution.fingerprint
        {
            return Err(bad(
                "molecular receiver rows require the original source, frequency and surface",
            ));
        }
        Ok(())
    }

    /// Apply the original adaptive rows without new quadrature or a new solve.
    /// A rejected/cancelled evaluation returns no partial observation and does
    /// not mutate the prepared rows or the source boundary solution.
    pub fn evaluate(
        &self,
        solution: &MolecularRadiation,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VectorObservation, MolecularError> {
        poll(&mut cancelled)?;
        self.admit_solution(solution)?;
        Ok(self
            .rows
            .evaluate(&solution.pressure, &solution.velocity, &mut cancelled)?)
    }

    pub(crate) fn evaluate_fields(
        &self,
        solution: &MolecularRadiation,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<receiver::Fields, MolecularError> {
        poll(cancelled)?;
        self.admit_solution(solution)?;
        Ok(self
            .rows
            .evaluate_fields(&solution.pressure, &solution.velocity, cancelled)?)
    }
}

impl MolecularApertureResponse {
    /// Observe the same coherently combined outward aperture volume flows
    /// with adaptive paired receivers. Source superposition precedes one row
    /// contraction; every aperture's mutual boundary load and original volume
    /// gain remain intact. No additional boundary solve or end correction.
    pub fn observe_volume_flows_near_field(
        &self,
        surface: &SpherePanels,
        volume_flows: &[C64],
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        options: MolecularNearFieldOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VectorObservation, MolecularError> {
        poll(&mut cancelled)?;
        options.admit(surface.areas().len(), points.len())?;
        self.combined_volume_flows(volume_flows, &mut cancelled)?
            .observe_near_field(surface, points, minimum_clearance_m, options, cancelled)
    }
}

#[cfg(test)]
#[path = "near_field/tests.rs"]
mod tests;
