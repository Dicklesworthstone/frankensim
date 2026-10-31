//! One retained physical boundary solve supplies loading and observation.
//! Material-coated passive bodies participate in the same mutual scattering
//! problem as every finite source; pressure is never fitted or symmetrized.
use super::*;
use crate::helmholtz::aperture::{
    project_real_velocity_rows, reciprocity_and_passivity, reciprocity_defect,
};
use crate::helmholtz::bc::{IncidentSource, Monopole};

pub mod directivity;
pub use directivity::{
    IterativeCoherentDirectivity, IterativeDirectivityOptions, IterativeFarFieldSample,
};

/// Immutable exp(-i omega t) fields at one positive physical frequency.
/// Point-source solutions are scattered fields; finite-source solutions are
/// total outgoing fields. The owning geometry/source/material binding stays
/// borrowed so these fields cannot drift away from their physical meaning.
pub struct IterativeBoundaryResponse<'o, 's> {
    observer: &'o IterativeBoundaryObservation<'s>,
    omega: f64,
    solutions: Vec<helmholtz::RadiationSolution>,
    reports: Vec<IterativeSourceReport>,
}

impl<'s> IterativeBoundaryObservation<'s> {
    /// Retain every source pressure AND normal-velocity array from the same
    /// shared FMM/preconditioner. `maximum_retained_entries` counts complex
    /// boundary coefficients, bounded by MAX_OBSERVATION_ENTRIES, before any
    /// constitutive query or solve. It is not an RSS estimate. The old receiver
    /// path remains streaming and does not inherit this extra retention cap.
    /// This returns unwindowed physical fields; DC or unsupported material
    /// frequencies refuse instead of creating artificial boundary solutions.
    pub fn solve_frequency(
        &self,
        omega: f64,
        maximum_retained_entries: usize,
        gate: &CancelGate,
    ) -> Result<IterativeBoundaryResponse<'_, 's>, IterativeObservationError> {
        poll(gate)?;
        self.observation.admit_frequency(omega)?;
        if omega == 0.0 {
            return Err(TransferError::Invalid(
                "retained boundary fields require positive frequency",
            )
            .into());
        }
        let entries = self
            .sources()
            .len()
            .checked_mul(self.field_geometry().0.areas().len())
            .and_then(|n| n.checked_mul(2))
            .ok_or(TransferError::Invalid("retained boundary extent overflow"))?;
        cap(
            entries,
            maximum_retained_entries.min(MAX_OBSERVATION_ENTRIES),
        )?;
        let (solutions, reports) = self.visit_solutions(
            omega,
            gate,
            || {
                Ok((
                    Vec::with_capacity(self.sources().len()),
                    Vec::with_capacity(self.sources().len()),
                ))
            },
            |(solutions, reports), _, radiation, report| {
                solutions.push(radiation);
                reports.push(report);
                Ok(())
            },
        )?;
        Ok(IterativeBoundaryResponse {
            observer: self,
            omega,
            solutions,
            reports,
        })
    }
}

impl<'o, 's> IterativeBoundaryResponse<'o, 's> {
    pub fn surface(&self) -> &'s SpherePanels {
        self.observer.field_geometry().0
    }
    pub fn medium(&self) -> Medium {
        self.observer.field_geometry().1
    }
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    pub fn sources(&self) -> &[SourceId] {
        self.observer.sources()
    }
    pub fn solutions(&self) -> &[helmholtz::RadiationSolution] {
        &self.solutions
    }
    pub fn reports(&self) -> &[IterativeSourceReport] {
        &self.reports
    }
    pub fn operator_hash(&self) -> ContentHash {
        self.observer.operator_hash()
    }

    /// Exact direct unit source for a compact-source column. Volume flow is
    /// one m³/s; force is one N along the source's declared direction. Its
    /// field must be added to the retained outgoing scattered trace.
    pub fn source_incident(
        &self,
        index: usize,
    ) -> Result<Option<IncidentSource>, IterativeObservationError> {
        if index >= self.solutions.len() {
            return Err(TransferError::Invalid("retained source index is out of range").into());
        }
        Ok(match &self.observer.observation.excitation {
            Excitation::Point(owner) => Some(owner.compact_sources[index].unit_incident()),
            Excitation::Surface { .. } => None,
        })
    }

    /// Volume-only adapter for the original point-source API. Force columns
    /// explicitly refuse instead of silently losing their direct field; use
    /// `source_incident` for a mixed compact-source response.
    pub fn source_monopole(
        &self,
        index: usize,
    ) -> Result<Option<Monopole>, IterativeObservationError> {
        match self.source_incident(index)? {
            Some(IncidentSource::Monopole(source)) => Ok(Some(source)),
            None => Ok(None),
            Some(_) => Err(TransferError::Invalid(
                "retained force source requires the typed source_incident accessor",
            )
            .into()),
        }
    }

    /// Observe the original bound receivers without another source solve.
    /// The same adaptive vector rows and exact direct-field addition as the
    /// streaming path preserve phase, pressure units and reactive velocity.
    pub fn observe_receivers(
        &self,
        gate: &CancelGate,
    ) -> Result<IterativeObservationSample, IterativeObservationError> {
        poll(gate)?;
        let (_, medium, geometry, quadrature) = self.observer.field_geometry();
        let rows = geometry.prepare_velocity_cancellable(
            self.omega / medium.sound_speed,
            medium,
            quadrature,
            gate,
        )?;
        let mut output = IterativeObservationSample {
            omega_rad_s: self.omega,
            operator_hash: self.operator_hash(),
            transfer: Vec::with_capacity(self.observer.observation.channel_count()),
            sources: self.reports.clone(),
        };
        for (index, radiation) in self.solutions.iter().enumerate() {
            poll(gate)?;
            self.observer.append_receiver_field(
                self.omega,
                index,
                &rows,
                radiation,
                &mut output.transfer,
            )?;
        }
        poll(gate)?;
        self.observer.admit_sample(&output)?;
        Ok(output)
    }

    /// Bind the retained fields to the common harmonic/coherence and microphone
    /// owner. Only receiver quadrature is prepared; no boundary solve is repeated.
    pub fn prepare_harmonic(
        &self,
        gate: &CancelGate,
    ) -> Result<PreparedHarmonicObservation, IterativeObservationError> {
        let sample = self.observe_receivers(gate)?;
        Ok(PreparedHarmonicObservation::from_samples(
            self.omega,
            self.medium(),
            self.sources(),
            self.observer.receivers(),
            sample.operator_hash,
            sample.transfer,
        )?)
    }

    /// Full self/mutual reaction of finite physical source coordinates, with
    /// all passive material patches and disconnected bodies present. Z_ij is
    /// integral b_i p_j dA under the original velocity/work map; no transpose
    /// averaging, positive loading or fitted scalar radiation is substituted.
    /// Point-source self-loading has no finite surface work port and refuses.
    /// Matrix and panel-product budgets include every original source pair.
    pub fn velocity_port_impedance(
        &self,
        maximum_matrix_entries: usize,
        maximum_panel_products: u64,
        gate: &CancelGate,
    ) -> Result<BoundaryPortImpedance, IterativeObservationError> {
        poll(gate)?;
        let Excitation::Surface { owner, .. } = &self.observer.observation.excitation else {
            return Err(TransferError::Invalid(
                "point-source self-loading has no finite velocity work port",
            )
            .into());
        };
        let n = self.sources().len();
        let entries = n
            .checked_mul(n)
            .ok_or(TransferError::Invalid("port matrix extent overflow"))?;
        cap(entries, maximum_matrix_entries.min(MAX_OBSERVATION_ENTRIES))?;
        let products = (entries as u64)
            .checked_mul(self.surface().areas().len() as u64)
            .ok_or(TransferError::Invalid("port projection work overflow"))?;
        if products > maximum_panel_products {
            return Err(TransferError::Invalid(
                "port projection exceeds declared panel-product budget",
            )
            .into());
        }
        let pressures: Vec<_> = self
            .solutions
            .iter()
            .map(|solution| solution.pressure.as_slice())
            .collect();
        let impedance = project_real_velocity_rows(
            self.surface(),
            n,
            |row, panel| owner.fields[row][panel].re,
            &pressures,
            Some(gate),
        )?;
        // Scaling only the diagnostic's scratch copy avoids overflow or
        // underflow; the physical matrix is returned byte-for-byte unchanged.
        let scale = impedance
            .iter()
            .fold(0.0_f64, |s, z| s.max(z.re.abs()).max(z.im.abs()));
        poll(gate)?;
        let (reciprocity_defect, passivity_margin) = if scale == 0.0 {
            (0.0, Some(0.0))
        } else {
            let scaled: Vec<_> = impedance
                .iter()
                .map(|z| C64::new(z.re / scale, z.im / scale))
                .collect();
            // Reuse fs-la's service allowance for all 60 Jacobi sweeps.
            // Large physical matrices remain available without claiming an
            // eigenvalue diagnostic that would block cooperative cancellation.
            let dimension = 2 * n;
            let eigen_work = dimension
                .checked_pow(3)
                .and_then(|work| work.checked_mul(1024));
            if eigen_work.is_some_and(|work| work <= fs_la::eigen::MAX_EIGEN_UNPOLLED_SCALAR_WORK)
                && fs_la::eigen::admit_jacobi_eigh(dimension).is_ok()
            {
                let (defect, margin) = reciprocity_and_passivity(&scaled, n);
                (defect, Some(margin))
            } else {
                (reciprocity_defect(&scaled, n), None)
            }
        };
        poll(gate)?;
        if !reciprocity_defect.is_finite()
            || passivity_margin.is_some_and(|margin| !margin.is_finite())
        {
            return Err(TransferError::Invalid("nonfinite port work diagnostic").into());
        }
        Ok(BoundaryPortImpedance {
            omega: self.omega,
            operator_hash: self.operator_hash(),
            sources: self.sources().to_vec(),
            impedance,
            reciprocity_defect,
            passivity_margin,
            projection_panel_products: products,
        })
    }
}

/// Unrepaired frequency-domain work-conjugate reaction matrix. This is sampled
/// loading, not a certified reciprocal/passive rational realization. Different
/// source coordinate kinds retain their SI work-dual effort/rate dimensions:
/// structural velocity uses force, aperture volume flow uses pressure.
#[derive(Clone, Debug)]
pub struct BoundaryPortImpedance {
    omega: f64,
    operator_hash: ContentHash,
    sources: Vec<SourceId>,
    impedance: Vec<C64>,
    reciprocity_defect: f64,
    passivity_margin: Option<f64>,
    projection_panel_products: u64,
}
impl BoundaryPortImpedance {
    pub const fn omega_rad_s(&self) -> f64 {
        self.omega
    }
    pub const fn operator_hash(&self) -> ContentHash {
        self.operator_hash
    }
    pub fn sources(&self) -> &[SourceId] {
        &self.sources
    }
    /// Complete row-major matrix: effort row, source-rate column.
    pub fn impedance(&self) -> &[C64] {
        &self.impedance
    }
    pub const fn reciprocity_defect(&self) -> f64 {
        self.reciprocity_defect
    }
    /// Smallest Hermitian-work eigenvalue divided by the matrix Frobenius norm.
    /// A negative value stays negative; it is never repaired into acceptance.
    /// `None` leaves this diagnostic unresolved when the original eigensolver's
    /// conservative unpolled-work allowance would be exceeded. The complete
    /// physical matrix and coherent source work remain available in that case.
    pub const fn passivity_margin(&self) -> Option<f64> {
        self.passivity_margin
    }
    pub const fn projection_panel_products(&self) -> u64 {
        self.projection_panel_products
    }

    /// Only original unit-volume-flow aperture coordinates have a work-equivalent
    /// pressure/flow matrix [Pa s/m^3]. A structural mode is never relabeled.
    pub fn acoustic_impedance(&self) -> Result<&[C64], IterativeObservationError> {
        if self
            .sources
            .iter()
            .any(|source| !matches!(source, SourceId::ApertureVolumeFlow(_)))
        {
            return Err(TransferError::Invalid(
                "acoustic impedance requires only aperture-volume-flow coordinates",
            )
            .into());
        }
        Ok(&self.impedance)
    }

    /// Z times the original coordinate rates. The fluid reaction ON the
    /// driven structure is the negative of these resisting efforts.
    pub fn resisting_efforts(&self, rates: &[C64]) -> Result<Vec<C64>, IterativeObservationError> {
        let n = self.sources.len();
        if rates.len() != n || rates.iter().any(|&value| !finite(value)) {
            return Err(TransferError::Invalid(
                "port rates must match the finite source coordinates",
            )
            .into());
        }
        let efforts: Vec<_> = self
            .impedance
            .chunks_exact(n)
            .map(|row| {
                row.iter()
                    .zip(rates)
                    .fold(C64::ZERO, |sum, (&z, &rate)| sum + z * rate)
            })
            .collect();
        if efforts.iter().any(|&value| !finite(value)) {
            return Err(TransferError::Invalid("nonfinite port resisting effort").into());
        }
        Ok(efforts)
    }

    /// Average harmonic work delivered to radiation plus passive walls [W].
    /// Coherent cross terms and their phases are retained exactly.
    pub fn input_power_w(&self, rates: &[C64]) -> Result<f64, IterativeObservationError> {
        let efforts = self.resisting_efforts(rates)?;
        let power = 0.5
            * rates
                .iter()
                .zip(efforts)
                .map(|(&v, f)| (v.conj() * f).re)
                .sum::<f64>();
        if !power.is_finite() {
            return Err(TransferError::Invalid("nonfinite source input power").into());
        }
        Ok(power)
    }
}

#[cfg(test)]
mod tests;
