//! Frequency-dependent passive boundaries in the real diffraction observer.
//!
//! Material source support and any lower reconstruction window remain explicit.
//! The constitutive owners, boundary solver, adaptive receiver integration and
//! FIR convergence owner each retain their original physical/numerical roles.
use super::*;
use crate::helmholtz::bc::{
    Admittance, IncidentField, IncidentSource, IncidentTraces, PreparedBoundary,
};
use crate::near_field::Options;
use std::cell::RefCell;

mod boundary;
pub use boundary::{LocalBoundary, LocalWall, LocalWallLayer};
mod iterative;
pub use iterative::{
    BoundaryWorkExposureFactor, BoundaryWorkExposureLimits, BoundaryWorkExposureLine,
    PreparedBoundaryWorkExposure, WindowedBoundaryWork,
    BoundaryPortImpedance, IterativeBoundaryObservation, IterativeBoundaryObservationBake,
    IterativeBoundaryPower, IterativeBoundaryResponse, IterativeCoherentDirectivity,
    IterativeDirectivityOptions, IterativeFarFieldSample, IterativeObservationError,
    IterativeObservationOptions, IterativeObservationSample, IterativeSourceReport,
};

/// Declared C4 lower reconstruction window: zero at/below `stop_hz`, one at/
/// above `pass_hz`. It is a band-limited observation, not an extrapolated wall
/// law. Its zero-phase taper needs the FIR owner's explicit reconstruction
/// latency and must pass that owner's original tail/wrap gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LowerBandWindow {
    stop_hz: f64,
    pass_hz: f64,
}
impl LowerBandWindow {
    pub fn try_new(stop_hz: f64, pass_hz: f64) -> Result<Self, TransferError> {
        if !stop_hz.is_finite() || stop_hz <= 0.0 || !pass_hz.is_finite() || pass_hz <= stop_hz {
            return Err(TransferError::Invalid(
                "lower observation window needs 0 < stop < pass",
            ));
        }
        Ok(Self { stop_hz, pass_hz })
    }
    pub const fn stop_hz(self) -> f64 {
        self.stop_hz
    }
    pub const fn pass_hz(self) -> f64 {
        self.pass_hz
    }

    fn weight(self, frequency_hz: f64) -> f64 {
        if frequency_hz <= self.stop_hz {
            return 0.0;
        }
        if frequency_hz >= self.pass_hz {
            return 1.0;
        }
        let x = (frequency_hz - self.stop_hz) / (self.pass_hz - self.stop_hz);
        if x <= 0.5 {
            observation_transition_low(x)
        } else {
            1.0 - observation_transition_low(1.0 - x)
        }
    }

    fn admit_ladder(self, declaration: &TransferDeclaration) -> Result<(), TransferError> {
        if self.pass_hz >= declaration.pass_band_upper_hz {
            return Err(TransferError::Invalid(
                "lower and upper observation transitions overlap",
            ));
        }
        for &length in &declaration.fft_lengths {
            let spacing = declaration.fs_hz / length as f64;
            let first = (self.pass_hz / spacing).ceil();
            let last = (declaration.pass_band_upper_hz / spacing).floor();
            // A necessary spectral-liveness precondition, not an error bound
            // for unsampled material resonances. Actual FIR gates still apply.
            if last - first < 3.0
                || (self.pass_hz - self.stop_hz) / spacing < 4.0
                || (declaration.taper_stop_hz - declaration.pass_band_upper_hz) / spacing < 4.0
            {
                return Err(TransferError::Invalid(
                    "observation band is unresolved on an FFT rung",
                ));
            }
        }
        Ok(())
    }
}

enum Excitation<'a> {
    Point(ScatteringObservation<'a>),
    Surface {
        owner: BemObservation<'a>,
        geometry: Geometry<'a>,
        quadrature: Options,
    },
    Combined {
        surface: BemObservation<'a>,
        compact: ScatteringObservation<'a>,
        sources: Vec<SourceId>,
    },
}

/// A prescribed finite boundary radiator or point source observed in the same
/// closed-body Robin diffraction solve as passive material/compliant patches.
/// The boundary solve includes mutual scattering between every component.
/// This is a stationary observation, not an inferred in-loop source load.
pub struct LocalBoundaryObservation<'a> {
    excitation: Excitation<'a>,
    boundary: LocalBoundary,
    highest_omega: f64,
    lower_window: Option<LowerBandWindow>,
}

/// The lower band declaration travels with the upper-band FIR artifact and
/// cannot be mistaken for a full-band material response through this owner.
#[derive(Clone, Debug)]
pub struct BoundaryObservationBake {
    transfer: ObservationTransferBake,
    lower_window: Option<LowerBandWindow>,
    boundary_hash: ContentHash,
}
impl BoundaryObservationBake {
    pub fn transfer(&self) -> &ObservationTransferBake {
        &self.transfer
    }
    pub const fn lower_band_window(&self) -> Option<LowerBandWindow> {
        self.lower_window
    }
    pub const fn boundary_hash(&self) -> ContentHash {
        self.boundary_hash
    }
    pub fn observer(&self) -> Result<ObservationFir, TransferError> {
        self.transfer.observer()
    }
}

impl<'a> LocalBoundaryObservation<'a> {
    pub fn for_scattering(
        owner: ScatteringObservation<'a>,
        boundary: LocalBoundary,
        lower_window: Option<LowerBandWindow>,
    ) -> Result<Self, TransferError> {
        admit_binding(
            owner.surface,
            owner.medium,
            owner.highest_omega,
            &boundary,
            lower_window,
        )?;
        Ok(Self {
            highest_omega: owner.highest_omega,
            excitation: Excitation::Point(owner),
            boundary,
            lower_window,
        })
    }

    /// Consume the original admitted finite velocity/volume-flow traces.
    /// Every driven panel must be explicitly rigid: adding a shunt there would
    /// change the actual prescribed volume flow or structural velocity.
    /// Geometry-routed or triangle-ruled source owners refuse because this
    /// closed Robin owner uses centroid collocation; no rule is substituted.
    pub fn for_surface(
        owner: BemObservation<'a>,
        boundary: LocalBoundary,
        highest_omega: f64,
        quadrature: Options,
        lower_window: Option<LowerBandWindow>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        if owner.policy.is_some() || !matches!(owner.panel_rule, helmholtz::PanelRule::Centroid) {
            return Err(TransferError::Invalid(
                "local-wall observation requires the original centroid boundary owner",
            ));
        }
        quadrature.validate()?;
        admit_binding(
            owner.surface,
            owner.medium,
            highest_omega,
            &boundary,
            lower_window,
        )?;
        for field in &owner.fields {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            if field
                .iter()
                .enumerate()
                .any(|(panel, &value)| value != C64::ZERO && !boundary.rigid_at(panel))
            {
                return Err(TransferError::Invalid(
                    "absorbing/compliant wall overlaps a prescribed source trace",
                ));
            }
        }
        let geometry = Geometry::new(owner.surface, &owner.points, owner.minimum_clearance_m)?;
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        Ok(Self {
            excitation: Excitation::Surface {
                owner,
                geometry,
                quadrature,
            },
            boundary,
            highest_omega,
            lower_window,
        })
    }

    /// Combine already-admitted finite wall/aperture and compact flow/force
    /// observations of the SAME boundary. Finite columns precede compact
    /// columns; every source retains its original SI coordinate and phase.
    /// This constructor performs no solve. Each frequency shares one boundary
    /// setup and receiver table across the complete basis. It does not infer
    /// point-source self-loading or solve for prescribed structural motion.
    pub fn combine_surface_and_scattering(
        surface: Self,
        scattering: Self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        let Excitation::Surface {
            owner, quadrature, ..
        } = &surface.excitation
        else {
            return Err(TransferError::Invalid(
                "first observation must contain finite surface sources",
            ));
        };
        let Excitation::Point(compact) = &scattering.excitation else {
            return Err(TransferError::Invalid(
                "second observation must contain compact sources",
            ));
        };
        if !std::ptr::eq(owner.surface, compact.surface)
            || owner.medium != compact.medium
            || owner.receivers != compact.receivers
            || owner.minimum_clearance_m != compact.minimum_clearance_m
            || surface.boundary.content_hash() != scattering.boundary.content_hash()
            || surface.highest_omega != scattering.highest_omega
            || surface.lower_window != scattering.lower_window
            || quadrature.relative_tolerance != compact.quadrature.relative_tolerance
            || quadrature.maximum_depth != compact.quadrature.maximum_depth
            || quadrature.maximum_kernel_evaluations
                != compact.quadrature.maximum_kernel_evaluations
        {
            return Err(TransferError::Invalid(
                "combined sources require identical boundary, medium, receivers, band and quadrature",
            ));
        }
        let count = owner
            .sources
            .len()
            .checked_add(compact.sources().len())
            .ok_or(TransferError::Invalid("combined source extent overflow"))?;
        cap(count, MAX_OBSERVATION_INPUTS)?;
        let entries = count
            .checked_mul(owner.receivers.len())
            .and_then(|n| n.checked_mul(4))
            .ok_or(TransferError::Invalid(
                "combined observation extent overflow",
            ))?;
        cap(entries, MAX_OBSERVATION_ENTRIES)?;
        let mut names = std::collections::BTreeSet::new();
        let mut sources = Vec::with_capacity(count);
        for source in owner.sources.iter().chain(compact.sources()) {
            if cancelled() {
                return Err(TransferError::Cancelled);
            }
            if !names.insert(source.name()) {
                return Err(TransferError::Invalid(
                    "combined source names must be distinct",
                ));
            }
            sources.push(source.clone());
        }
        let (Excitation::Surface { owner, .. }, Excitation::Point(compact)) =
            (surface.excitation, scattering.excitation)
        else {
            unreachable!()
        };
        Ok(Self {
            excitation: Excitation::Combined {
                surface: owner,
                compact,
                sources,
            },
            boundary: surface.boundary,
            highest_omega: surface.highest_omega,
            lower_window: surface.lower_window,
        })
    }

    pub fn sources(&self) -> &[SourceId] {
        match &self.excitation {
            Excitation::Point(owner) => owner.sources(),
            Excitation::Surface { owner, .. } => &owner.sources,
            Excitation::Combined { sources, .. } => sources,
        }
    }
    fn has_force_sources(&self) -> bool {
        self.sources()
            .iter()
            .any(|source| matches!(source, SourceId::PointForce(_)))
    }
    pub fn receivers(&self) -> &[ReceiverSpec] {
        match &self.excitation {
            Excitation::Point(owner) => &owner.receivers,
            Excitation::Surface { owner, .. } => &owner.receivers,
            Excitation::Combined { compact, .. } => &compact.receivers,
        }
    }
    pub fn geometric_delays_s(&self) -> Vec<f64> {
        match &self.excitation {
            Excitation::Point(owner) => owner.geometric_delays_s(),
            Excitation::Surface { owner, .. } => owner.geometric_delays_s(),
            Excitation::Combined {
                surface, compact, ..
            } => {
                let mut delays = surface.geometric_delays_s();
                delays.extend(compact.geometric_delays_s());
                delays
            }
        }
    }
    pub const fn lower_band_window(&self) -> Option<LowerBandWindow> {
        self.lower_window
    }

    fn field_geometry(&self) -> (&'a SpherePanels, Medium, &Geometry<'a>, Options) {
        match &self.excitation {
            Excitation::Point(owner) | Excitation::Combined { compact: owner, .. } => (
                owner.surface,
                owner.medium,
                &owner.receiver_geometry,
                owner.quadrature,
            ),
            Excitation::Surface {
                owner,
                geometry,
                quadrature,
            } => (owner.surface, owner.medium, geometry, *quadrature),
        }
    }
    fn source_incident(&self, index: usize) -> Option<IncidentSource> {
        match &self.excitation {
            Excitation::Point(owner) => owner
                .compact_sources
                .get(index)
                .map(CompactSourceSpec::unit_incident),
            Excitation::Surface { .. } => None,
            Excitation::Combined {
                surface, compact, ..
            } => index
                .checked_sub(surface.fields.len())
                .and_then(|index| compact.compact_sources.get(index))
                .map(CompactSourceSpec::unit_incident),
        }
    }
    fn source_velocity(&self, index: usize) -> Option<&[C64]> {
        match &self.excitation {
            Excitation::Point(_) => None,
            Excitation::Surface { owner, .. } | Excitation::Combined { surface: owner, .. } => {
                owner.fields.get(index).map(Vec::as_slice)
            }
        }
    }

    pub fn operator_hash(&self) -> ContentHash {
        let mut hash = DomainHasher::new("fs-bem:local-boundary-observation:v1");
        hash.update(self.boundary.content_hash().as_bytes());
        hash.update(&self.highest_omega.to_bits().to_le_bytes());
        match &self.excitation {
            Excitation::Point(owner) => {
                hash.update(&[0]);
                hash.update(owner.operator_hash().as_bytes());
            }
            Excitation::Surface {
                owner, quadrature, ..
            } => {
                hash.update(&[1]);
                hash.update(owner.operator_hash().as_bytes());
                hash.update(&quadrature.relative_tolerance.to_bits().to_le_bytes());
                hash.update(&(quadrature.maximum_depth as u64).to_le_bytes());
                hash.update(&(quadrature.maximum_kernel_evaluations as u64).to_le_bytes());
            }
            Excitation::Combined {
                surface, compact, ..
            } => {
                hash.update(&[2]);
                hash.update(surface.operator_hash().as_bytes());
                hash.update(compact.operator_hash().as_bytes());
            }
        }
        match self.lower_window {
            None => hash.update(&[0]),
            Some(window) => {
                hash.update(&[1]);
                hash.update(&window.stop_hz.to_bits().to_le_bytes());
                hash.update(&window.pass_hz.to_bits().to_le_bytes());
            }
        }
        hash.finalize()
    }

    /// Physical, unwindowed source-major p/grad transfer. Unsupported material
    /// queries refuse before boundary preparation, even below a selected window.
    pub fn physical_frequency_sample(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        self.admit_frequency(omega)?;
        if omega == 0.0 {
            if self.has_force_sources() {
                return Err(TransferError::Invalid(
                    "force-dipole observation has no admitted DC solution",
                ));
            }
            if self.boundary.supported_band_hz().0 > 0.0 {
                return Err(TransferError::Invalid(
                    "DC is outside local-wall source support",
                ));
            }
            return Ok(vec![C64::ZERO; self.channel_count()]);
        }
        let beta = self.boundary.normalized_admittance(omega, &mut cancelled)?;
        let law = Admittance::Local {
            beta: &beta,
            declared_active: false,
        };
        let result = match &self.excitation {
            Excitation::Point(owner) => owner
                .frequency_sample_with_admittance(omega, law, &mut cancelled)
                .map_err(|error| match error {
                    HelmholtzError::Cancelled => TransferError::Cancelled,
                    other => TransferError::Helmholtz(other),
                })?,
            Excitation::Surface { .. } | Excitation::Combined { .. } => {
                let (surface, medium, geometry, quadrature) = self.field_geometry();
                let k = omega / medium.sound_speed;
                let boundary = PreparedBoundary::new(
                    surface,
                    k,
                    medium,
                    law,
                    helmholtz::Formulation::BurtonMiller,
                )?;
                if cancelled() {
                    return Err(TransferError::Cancelled);
                }
                let rows = geometry.prepare_velocity(k, medium, quadrature)?;
                let momentum = C64::new(0.0, omega * medium.density);
                let mut output = Vec::with_capacity(self.channel_count());
                // Reuse factorization and receiver quadrature for every source;
                // no matrix or integrated Green row is rebuilt inside this loop.
                for index in 0..self.sources().len() {
                    if cancelled() {
                        return Err(TransferError::Cancelled);
                    }
                    let incident = self.source_incident(index).map(|source| [source]);
                    let radiation = if let Some(sources) = &incident {
                        let traces = IncidentTraces::from_admitted_field(
                            surface,
                            k,
                            medium,
                            IncidentField::Sources(sources),
                            None,
                        )?;
                        boundary.scatter(&traces)?.scattered
                    } else {
                        boundary
                            .solve(
                                self.source_velocity(index)
                                    .ok_or(TransferError::Invalid("missing finite source trace"))?,
                            )?
                            .radiation
                    };
                    let received = rows.evaluate_velocity(&radiation)?;
                    for ((&pressure, velocity), receiver) in received
                        .scalar
                        .pressure
                        .iter()
                        .zip(&received.particle_velocity_m_s)
                        .zip(self.receivers())
                    {
                        if let Some(sources) = &incident {
                            let direct = IncidentField::Sources(sources).pressure_velocity_at(
                                k,
                                medium,
                                receiver.position_m,
                            )?;
                            output.push(pressure + direct.pressure_pa);
                            output.extend(
                                (0..3).map(|i| momentum * (velocity[i] + direct.velocity_m_s[i])),
                            );
                        } else {
                            output.push(pressure);
                            output.extend(velocity.iter().map(|&v| momentum * v));
                        }
                    }
                }
                output
            }
        };
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        if result.len() != self.channel_count() || result.iter().any(|&v| !finite(v)) {
            return Err(TransferError::Invalid(
                "nonfinite or incomplete local-wall observation",
            ));
        }
        Ok(result)
    }

    /// Physical single-frequency observations of the same local material or
    /// compliant boundary. This uses the unwindowed constitutive response;
    /// unsupported material frequencies still refuse before the solve.
    pub fn prepare_harmonic(
        &self,
        omega: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHarmonicObservation, TransferError> {
        let medium = self.field_geometry().1;
        harmonic::admit_frequency(omega, medium)?;
        let samples = self.physical_frequency_sample(omega, cancelled)?;
        PreparedHarmonicObservation::from_samples(
            omega,
            medium,
            self.sources(),
            self.receivers(),
            self.operator_hash(),
            samples,
        )
    }

    /// The explicitly lower-windowed transfer used by `bake`. A zero window
    /// skips constitutive queries and solves, never cancellation or band checks.
    pub fn band_limited_frequency_sample(
        &self,
        omega: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, TransferError> {
        if cancelled() {
            return Err(TransferError::Cancelled);
        }
        self.admit_frequency(omega)?;
        let weight = self
            .lower_window
            .map_or(1.0, |window| window.weight(omega / core::f64::consts::TAU));
        if weight == 0.0 {
            return Ok(vec![C64::ZERO; self.channel_count()]);
        }
        let mut field = self.physical_frequency_sample(omega, &mut cancelled)?;
        if weight != 1.0 {
            for value in &mut field {
                *value = value.scale(weight);
            }
        }
        Ok(field)
    }

    pub fn bake(
        &self,
        declaration: &TransferDeclaration,
        cancelled: impl FnMut() -> bool,
    ) -> Result<BoundaryObservationBake, TransferError> {
        self.admit_declaration(declaration, self.operator_hash())?;
        let cancel = RefCell::new(cancelled);
        let transfer = bake_observation_transfers_with_errors(
            |omega| self.band_limited_frequency_sample(omega, || (*cancel.borrow_mut())()),
            declaration,
            self.lower_window.map_or(0.0, |window| window.pass_hz),
            || (*cancel.borrow_mut())(),
        )?;
        Ok(BoundaryObservationBake {
            transfer,
            lower_window: self.lower_window,
            boundary_hash: self.boundary.content_hash(),
        })
    }

    fn admit_declaration(
        &self,
        declaration: &TransferDeclaration,
        producer_hash: ContentHash,
    ) -> Result<(), TransferError> {
        // The common FIR owner sets bin zero without querying this producer.
        // Only a caller-declared lower window can remove force-source DC.
        if self.has_force_sources() && self.lower_window.is_none() {
            return Err(TransferError::Invalid(
                "force-dipole FIR requires an explicit lower band window",
            ));
        }
        declaration.admit_work()?;
        if declaration.sources != self.sources()
            || declaration.receivers != self.receivers()
            || declaration.geometric_delays_s != self.geometric_delays_s()
            || declaration.producer_hash != producer_hash
            || core::f64::consts::TAU * declaration.taper_stop_hz > self.highest_omega
            || matches!(
                declaration.producer_uncertainty,
                ProducerUncertainty::ExactAnalytic
            )
        {
            return Err(TransferError::Invalid(
                "cross-wired local-boundary observation declaration",
            ));
        }
        if let Some(window) = self.lower_window {
            window.admit_ladder(declaration)?;
        }
        Ok(())
    }

    fn channel_count(&self) -> usize {
        4 * self.sources().len() * self.receivers().len()
    }
    fn admit_frequency(&self, omega: f64) -> Result<(), TransferError> {
        if !omega.is_finite() || omega < 0.0 || omega > self.highest_omega {
            return Err(TransferError::Invalid(
                "local-boundary observation frequency outside declared band",
            ));
        }
        Ok(())
    }
}

fn admit_binding(
    surface: &SpherePanels,
    medium: Medium,
    highest_omega: f64,
    boundary: &LocalBoundary,
    lower: Option<LowerBandWindow>,
) -> Result<(), TransferError> {
    admit_binding_for_route(surface, medium, highest_omega, boundary, lower, false)
}

fn admit_binding_for_route(
    surface: &SpherePanels,
    medium: Medium,
    highest_omega: f64,
    boundary: &LocalBoundary,
    lower: Option<LowerBandWindow>,
    iterative: bool,
) -> Result<(), TransferError> {
    if boundary.panel_count() != surface.areas().len()
        || boundary.medium() != medium
        || !highest_omega.is_finite()
        || highest_omega <= 0.0
        || !(highest_omega * medium.density).is_finite()
    {
        return Err(TransferError::Invalid(
            "local wall is bound to another medium or panel shape",
        ));
    }
    boundary.admit_band(highest_omega, lower)?;
    if iterative {
        helmholtz::admit_iterative_resolution(surface, highest_omega / medium.sound_speed)?;
    } else {
        helmholtz::admit_resolution(surface, highest_omega / medium.sound_speed)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests;
