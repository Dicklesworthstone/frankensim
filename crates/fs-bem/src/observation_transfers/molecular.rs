//! One material compression law in both boundary reaction and receiver fields.
//!
//! The existing complex boundary solvers and molecular exterior own the physics. This
//! adapter supplies their complete signed source basis to the common FIR
//! owner. No absorption filter is applied after summing the source fields.
//! One explicitly selected boundary formulation spans the declared band. The
//! triangle-CBIE arm retains its geometric resonance-exclusion bound; the
//! centroid/disc Burton--Miller arm retains the original combined operator.
//! Neither switches discretizations inside a Fourier transfer.

use super::*;
use crate::helmholtz::bc::Admittance;
use crate::helmholtz::laplace::complex::ComplexOptions;
use crate::helmholtz::molecular::{
    MolecularBoundarySolution, MolecularError, MolecularExterior, MolecularNearFieldOptions,
    MolecularRadiation,
};
use fs_material::relaxation::RelaxationState;
use std::cell::RefCell;

/// A single solved source bank for both reaction and finite receivers.
pub mod prepared;
pub use prepared::{MolecularResponseLimits, MolecularResponsePower, PreparedMolecularResponse};

#[path = "molecular/near_field.rs"]
mod adaptive;
pub use adaptive::{MolecularReceiverError, MolecularReceiverQuadrature};

/// The original material/boundary refusal or shared observation refusal.
#[derive(Debug)]
pub enum MolecularObservationError {
    Transfer(TransferError),
    Molecular(MolecularError),
}
impl core::fmt::Display for MolecularObservationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Transfer(error) => error.fmt(f),
            Self::Molecular(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for MolecularObservationError {}
impl From<TransferError> for MolecularObservationError {
    fn from(error: TransferError) -> Self {
        Self::Transfer(error)
    }
}
impl From<MolecularError> for MolecularObservationError {
    fn from(error: MolecularError) -> Self {
        Self::Molecular(error)
    }
}
fn invalid(what: &'static str) -> MolecularObservationError {
    TransferError::Invalid(what).into()
}

/// One original spatial discretization, fixed for the entire frequency band.
/// This choice does not certify spatial convergence or mechanical passivity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolecularObservationFormulation {
    /// Original triangle integration below the geometric CBIE exclusion bound.
    TriangleCbie,
    /// Original centroid/disc combined equation, including its low-frequency
    /// discretization error. No triangle quadrature is implied by this arm.
    CentroidBurtonMiller,
}

impl MolecularObservationFormulation {
    /// Conservative kernel-evaluation reservation of the original assembler.
    pub const fn kernel_terms_per_pair(self) -> u64 {
        match self {
            Self::TriangleCbie => 512,
            Self::CentroidBurtonMiller => 8,
        }
    }

    /// Stable physical/numerical description shared by observer reports.
    pub const fn label(self) -> &'static str {
        match self {
            Self::TriangleCbie => "molecular_triangle_cbie",
            Self::CentroidBurtonMiller => "molecular_centroid_burton_miller",
        }
    }

    const fn boundary_formulation(self) -> helmholtz::Formulation {
        match self {
            Self::TriangleCbie => helmholtz::Formulation::PlainCbie,
            Self::CentroidBurtonMiller => helmholtz::Formulation::BurtonMiller,
        }
    }
}

/// Source-bound finite-distance observation and work-conjugate reaction of a
/// closed body in an arbitrary admitted relaxing gas.
///
/// This inherits the original molecular owner's inviscid momentum equation:
/// molecular compression loss is present; shear viscosity, heat conduction,
/// wall layers and mean flow are not. The gas data are evaluated independently
/// of instrument identity. Mesh, receiver quadrature and Fourier convergence
/// remain separate; this adapter grants no continuum-accuracy certificate.
pub struct MolecularBemObservation<'a> {
    boundary: BemObservation<'a>,
    exterior: MolecularExterior,
    formulation: MolecularObservationFormulation,
    limits: ComplexOptions,
    plain_cbie_limit_rad_per_m: f64,
    local_boundary: Option<LocalBoundary>,
    near_field_options: Option<MolecularNearFieldOptions>,
}

impl<'a> MolecularBemObservation<'a> {
    /// Bind original source traces and receivers to one evaluated molecular
    /// medium. The entire declared band must be below the same closed-geometry
    /// Dirichlet exclusion bound used by the ordinary radiation owner.
    /// The source's equilibrium speed provides a conservative upper bound on
    /// |k| for positive Debye reservoirs; its frozen speed gives the earliest
    /// physical wavefront. Neither changes the actual finite-frequency solve.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        surface: &'a SpherePanels,
        source: RelaxationState,
        sources: &[SourceId],
        traces: &[Vec<f64>],
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
        band_rad_s: [f64; 2],
        limits: ComplexOptions,
    ) -> Result<Self, MolecularObservationError> {
        Self::new_with_formulation(
            surface,
            source,
            sources,
            traces,
            receivers,
            minimum_clearance_m,
            band_rad_s,
            MolecularObservationFormulation::TriangleCbie,
            limits,
        )
    }

    /// Select one original boundary formulation for the entire band, including
    /// reaction estimates and receiver queries below the CBIE exclusion bound.
    /// Burton--Miller permits bands crossing that bound while retaining the
    /// same wavelength, storage and work admission at every actual frequency.
    /// Its low-frequency resistance can be spatially biased; the Fourier
    /// convergence check does not establish boundary-quadrature accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_formulation(
        surface: &'a SpherePanels,
        source: RelaxationState,
        sources: &[SourceId],
        traces: &[Vec<f64>],
        receivers: &[ReceiverSpec],
        minimum_clearance_m: f64,
        band_rad_s: [f64; 2],
        formulation: MolecularObservationFormulation,
        limits: ComplexOptions,
    ) -> Result<Self, MolecularObservationError> {
        let exterior = MolecularExterior::new(source, band_rad_s)?;
        let gas = exterior.source().gas();
        let policy = GeometryPolicy::new(surface).map_err(TransferError::from)?;
        let plain_cbie_limit_rad_per_m = policy.plain_cbie_limit();
        if formulation == MolecularObservationFormulation::TriangleCbie
            && band_rad_s[1] / gas.sound_speed > plain_cbie_limit_rad_per_m
        {
            return Err(invalid(
                "molecular observation band exceeds the triangle-CBIE resonance-exclusion bound",
            ));
        }
        let panels = surface.areas().len();
        let entries = panels
            .checked_mul(panels)
            .ok_or_else(|| invalid("molecular observation matrix extent overflow"))?;
        if limits.max_panels == 0
            || limits.max_panels > helmholtz::MAX_DENSE_PANELS
            || panels > limits.max_panels
            || entries
                .checked_mul(32)
                .is_none_or(|bytes| bytes > limits.max_bytes)
            || (entries as u64)
                .checked_mul(formulation.kernel_terms_per_pair())
                .is_none_or(|work| work > limits.max_work)
        {
            return Err(invalid(
                "molecular observation exceeds its complete boundary storage or work declaration",
            ));
        }
        let mut boundary = BemObservation::new_admitted(
            surface,
            Medium::from(gas),
            sources,
            traces,
            receivers,
            minimum_clearance_m,
            true,
        )?;
        // Only the FIR's exact prefix changes. Each complex Green function
        // already contains the complete travel and dispersive phase, once.
        let speed_ratio = gas.sound_speed / exterior.source().frozen_sound_speed_m_s();
        for delay in &mut boundary.clearances_s {
            *delay *= speed_ratio;
        }
        Ok(Self {
            boundary,
            exterior,
            formulation,
            limits,
            plain_cbie_limit_rad_per_m,
            local_boundary: None,
            near_field_options: None,
        })
    }

    /// Attach the existing frequency-dependent material/mechanical wall laws
    /// to this same relaxing exterior, source basis and receiver geometry.
    /// Every panel supporting a prescribed source must remain explicitly
    /// rigid. Otherwise `v = v_drive - Y p` would silently change that source's
    /// structural normal velocity or prescribed aperture volume flow.
    ///
    /// The immutable wall table is evaluated once per frequency. Its original
    /// specific impedance, source support and exp(-i omega t) convention are
    /// retained; only the exterior's complex-k propagation changes. The real
    /// equilibrium rho*c normalizes Y and is never replaced by a complex wave
    /// impedance. Porous laws retain their original equivalent-fluid model;
    /// no new molecular relaxation law inside their pores is inferred.
    /// A positive material lower cutoff remains a harmonic-only binding here:
    /// this owner's FIR artifact has no explicit lower-band window, so its
    /// bake refuses instead of extrapolating such a wall to zero frequency.
    ///
    /// Surface work and gas absorption remain separate channels. This does
    /// not establish spatial convergence, a passive reduced realization or
    /// reciprocity of the original collocation discretization.
    pub fn with_local_boundary(
        mut self,
        boundary: LocalBoundary,
    ) -> Result<Self, MolecularObservationError> {
        if boundary.panel_count() != self.boundary.surface.areas().len()
            || boundary.medium() != self.boundary.medium
        {
            return Err(invalid(
                "molecular wall boundary must use the original panel count and equilibrium gas",
            ));
        }
        let tau = core::f64::consts::TAU;
        let band = self.exterior.band_rad_s();
        boundary.admit_physical_band_hz([band[0] / tau, band[1] / tau])?;
        for field in &self.boundary.fields {
            if field
                .iter()
                .enumerate()
                .any(|(panel, &value)| value != C64::ZERO && !boundary.rigid_at(panel))
            {
                return Err(invalid(
                    "molecular absorbing/compliant wall overlaps a prescribed source trace",
                ));
            }
        }
        self.local_boundary = Some(boundary);
        Ok(self)
    }

    /// Original immutable wall table, including its material source support.
    pub fn local_boundary(&self) -> Option<&LocalBoundary> {
        self.local_boundary.as_ref()
    }

    pub fn source(&self) -> &RelaxationState {
        self.exterior.source()
    }
    pub const fn formulation(&self) -> MolecularObservationFormulation {
        self.formulation
    }
    /// Geometry diagnostic; this is not a band cap on the Burton--Miller arm.
    pub fn plain_cbie_limit_rad_per_m(&self) -> f64 {
        self.plain_cbie_limit_rad_per_m
    }
    pub fn geometric_delays_s(&self) -> Vec<f64> {
        self.boundary.geometric_delays_s()
    }

    /// Cache identity includes the original gas and EVERY evaluated reservoir.
    /// Equal equilibrium rho/c with different relaxation cannot reuse a bake.
    pub fn operator_hash(&self) -> ContentHash {
        let mut h = DomainHasher::new("fs-bem:observation-boundary:molecular:v2");
        h.update(self.formulation.label().as_bytes());
        h.update(&self.boundary.operator_hash().0);
        let gas = self.source().gas();
        for value in [
            gas.temperature,
            gas.pressure,
            gas.density,
            gas.sound_speed,
            gas.dynamic_viscosity,
            gas.thermal_conductivity,
            gas.gamma,
            gas.specific_gas_constant,
            gas.specific_heat_cp,
            gas.prandtl,
            gas.characteristic_impedance,
            gas.water_mole_fraction,
            self.source().frozen_cv_j_kg_k(),
        ]
        .into_iter()
        .chain(self.exterior.band_rad_s())
        {
            h.update(&value.to_bits().to_le_bytes());
        }
        h.update(&(self.source().processes().len() as u64).to_le_bytes());
        for process in self.source().processes() {
            let name = process.species.as_str();
            h.update(&(name.len() as u64).to_le_bytes());
            h.update(name.as_bytes());
            h.update(&(process.index as u64).to_le_bytes());
            h.update(&process.heat_capacity_j_kg_k.to_bits().to_le_bytes());
            h.update(&process.relaxation_time_s.to_bits().to_le_bytes());
        }
        if let Some(boundary) = &self.local_boundary {
            h.update(b"original-local-boundary-v1");
            h.update(boundary.content_hash().as_bytes());
        }
        if let Some(options) = self.near_field_options {
            h.update(b"original-adaptive-molecular-receivers-v1");
            h.update(
                &options
                    .quadrature
                    .relative_tolerance
                    .to_bits()
                    .to_le_bytes(),
            );
            for value in [
                options.quadrature.maximum_depth,
                options.quadrature.maximum_kernel_evaluations,
                options.maximum_row_bytes,
            ] {
                h.update(&(value as u64).to_le_bytes());
            }
        }
        h.finalize()
    }

    fn boundary_solutions(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<MolecularBoundarySolution>, MolecularObservationError> {
        let coefficients = self
            .local_boundary
            .as_ref()
            .map(|boundary| boundary.normalized_admittance(omega, &mut *cancelled))
            .transpose()?;
        let admittance =
            coefficients
                .as_ref()
                .map_or(Admittance::Rigid, |beta| Admittance::Local {
                    beta,
                    declared_active: false,
                });
        let mut solutions = Vec::with_capacity(self.boundary.fields.len());
        for batch in self
            .boundary
            .fields
            .chunks(helmholtz::MAX_RADIATION_FIELDS_PER_BATCH)
        {
            let fields: Vec<_> = batch.iter().map(Vec::as_slice).collect();
            solutions.extend(self.exterior.solve_boundary_batch(
                self.boundary.surface,
                omega,
                &fields,
                admittance,
                self.formulation.boundary_formulation(),
                &self.limits,
                &mut *cancelled,
            )?);
        }
        Ok(solutions)
    }

    fn solutions(
        &self,
        omega: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<MolecularRadiation>, MolecularObservationError> {
        Ok(self
            .boundary_solutions(omega, cancelled)?
            .into_iter()
            .map(|solution| solution.radiation)
            .collect())
    }

    /// Complete exp(-i omega t) p/gradient rows from the original material
    /// boundary solution. The Neumann datum remains i*rho*omega*v; replacing
    /// omega with k*c would assign a different force to the same source.
    pub fn frequency_sample(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, MolecularObservationError> {
        if cancelled() {
            return Err(TransferError::Cancelled.into());
        }
        if !omega.is_finite() || omega < 0.0 {
            return Err(invalid("invalid molecular observation frequency"));
        }
        if omega == 0.0 {
            if self
                .local_boundary
                .as_ref()
                .is_some_and(|wall| wall.supported_band_hz().0 > 0.0)
            {
                return Err(invalid(
                    "molecular observation wall has no zero-frequency source support",
                ));
            }
            return Ok(vec![
                C64::ZERO;
                4 * self.boundary.points.len()
                    * self.boundary.fields.len()
            ]);
        }
        let solved = self.solutions(omega, &mut cancelled)?;
        Ok(self.observation_samples(&solved, &mut cancelled)?.0)
    }

    /// Prepare material-bound pressure and particle velocity for coherent
    /// fields, source-coherence statistics and microphone cross spectra.
    ///
    /// The original complex-wavenumber boundary solve supplies every pressure
    /// and gradient. Inviscid momentum remains `u = grad(p)/(i*rho*omega)`;
    /// neither the equilibrium `rho*c` nor a complex wave impedance replaces
    /// that relation. Subsequent observations reuse this complete transfer
    /// bank without another boundary solve or an absorption correction.
    ///
    /// Ideal first-order microphones retain their existing real equilibrium
    /// `rho*c` pressure-equivalent calibration. In this dispersive medium that
    /// reference is not the actual frequency-dependent plane-wave impedance,
    /// so unity plane-wave sensitivity is not implied. Spatial uncertainty
    /// and material/band/cancellation refusals remain with the original owners.
    pub fn prepare_harmonic(
        &self,
        omega: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicObservation, MolecularObservationError> {
        let medium = self.boundary.medium;
        harmonic::admit_frequency(omega, medium)?;
        let samples = self.frequency_sample(omega, cancelled)?;
        Ok(PreparedHarmonicObservation::from_samples(
            omega,
            medium,
            &self.boundary.sources,
            &self.boundary.receivers,
            self.operator_hash(),
            samples,
        )?
        .without_radiation_momentum())
    }

    /// One row of the work-conjugate boundary impedance at the same frequency.
    /// Projects each source's pressure against the original receiving trace
    /// and panel areas, preserving mixed structural/volume-flow units.
    /// This is the physical reaction the one-way approximation leaves out;
    /// no sign clipping, symmetrization or fitted loss is performed.
    pub fn impedance_row(
        &self,
        omega: f64,
        receiving_source: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, MolecularObservationError> {
        let test = self
            .boundary
            .fields
            .get(receiving_source)
            .ok_or_else(|| invalid("unknown molecular reaction receiving source"))?;
        let solved = self.solutions(omega, &mut cancelled)?;
        Ok(solved
            .iter()
            .map(|solution| {
                solution
                    .pressure_pa()
                    .iter()
                    .zip(test)
                    .zip(self.boundary.surface.areas())
                    .fold(C64::ZERO, |sum, ((&pressure, &trace), &area)| {
                        sum + pressure * trace.conj().scale(area)
                    })
            })
            .collect())
    }

    /// The original bounded FIR owner, including material/boundary cancellation.
    /// The geometric prefix uses frozen sound speed; the sampled response
    /// already contains dispersion, so there is no second propagation delay.
    pub fn bake(
        &self,
        declaration: &TransferDeclaration,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ObservationTransferBake, MolecularObservationError> {
        if self
            .local_boundary
            .as_ref()
            .is_some_and(|wall| wall.supported_band_hz().0 > 0.0)
        {
            return Err(invalid(
                "positive-band molecular wall needs an explicit lower-window observation owner",
            ));
        }
        let band = self.exterior.band_rad_s();
        let finest = declaration.fft_lengths.last().copied().unwrap_or(0);
        if declaration.sources != self.boundary.sources
            || declaration.receivers != self.boundary.receivers
            || declaration.geometric_delays_s != self.geometric_delays_s()
            || declaration.producer_hash != self.operator_hash()
            || core::f64::consts::TAU * declaration.taper_stop_hz > band[1]
            || finest == 0
            || core::f64::consts::TAU * declaration.fs_hz / (finest as f64) < band[0]
            || matches!(
                declaration.producer_uncertainty,
                ProducerUncertainty::ExactAnalytic
            )
        {
            return Err(invalid(
                "cross-wired or out-of-band molecular observation declaration",
            ));
        }
        let cancel = RefCell::new(cancelled);
        bake_observation_transfers_with_errors(
            |omega| self.frequency_sample(omega, || (*cancel.borrow_mut())()),
            declaration,
            0.0,
            || (*cancel.borrow_mut())(),
        )
    }
}

#[cfg(test)]
mod tests;
