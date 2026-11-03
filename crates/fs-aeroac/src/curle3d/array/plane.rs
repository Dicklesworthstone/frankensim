//! Exact ideal-plane Green functions for prescribed compact acoustic sources.
//!
//! For a plane through a with inward unit normal n, reflect y to
//! y' = y - 2 n (n.(y-a)). Its half-space Green function is G(x,y) + s G(x,y'),
//! with s=+1 for a rigid plane and s=-1 for pressure release. Differentiating
//! with respect to the ORIGINAL source position gives image force s R F,
//! R=I-2nn^T, while the image volume injection is s Q. Reflecting a dipole's
//! position without its orientation violates the wall boundary condition.
//!
//! This is the ideal-boundary specialization of the image construction in
//! Allen & Berkley, JASA 65 (1979), 943–950, doi:10.1121/1.382599. No constant
//! reflection coefficient is substituted for a finite-impedance wall: that
//! problem has an angle-dependent, generally nonlocal reflected Green function.

use super::*;
use fs_math::dd::Dd;

const LOSSLESS_SCOPE: &str = "Prescribed compact 3D Q/F sources in a homogeneous stationary \
lossless half-space with one infinite ideal rigid or pressure-release plane; direct and \
reflected pressure/velocity interfere coherently. No finite source-size validation, wall \
impedance, edge diffraction, motion, source feedback or measured-SPL qualification.";
const MOLECULAR_SCOPE: &str = "Prescribed compact 3D Q/F sources in a homogeneous stationary \
half-space with material-owned molecular bulk relaxation, inviscid momentum and one infinite \
ideal rigid or pressure-release plane; direct and reflected fields interfere coherently. \
No shear/conduction boundary layers, finite source-size validation, finite wall impedance, \
edge diffraction, motion, source feedback or measured-SPL qualification.";

/// Exact homogeneous boundary condition on an infinite stationary plane.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IdealPlaneBoundary {
    /// Zero normal particle velocity; the pressure Green function is even.
    Rigid,
    /// Zero acoustic pressure; the pressure Green function is odd.
    PressureRelease,
}

impl IdealPlaneBoundary {
    const fn sign(self) -> f64 {
        match self {
            Self::Rigid => 1.0,
            Self::PressureRelease => -1.0,
        }
    }
}

/// A plane and its explicitly selected physical half-space, in metres.
///
/// The normal points INTO the fluid. Sources must lie strictly in that
/// half-space; receivers may also lie on the plane. No geometry is snapped.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IdealAcousticPlane {
    point_m: [f64; 3],
    normal: [f64; 3],
    boundary: IdealPlaneBoundary,
}

impl IdealAcousticPlane {
    /// Admit finite geometry and normalize an arbitrary nonzero inward normal.
    pub fn new(
        point_m: [f64; 3],
        inward_normal: [f64; 3],
        boundary: IdealPlaneBoundary,
    ) -> Result<Self, AeroacError> {
        if point_m.iter().chain(&inward_normal).any(|v| !v.is_finite()) {
            return Err(nonfinite("acoustic reflection plane geometry"));
        }
        let scale = inward_normal.iter().fold(0.0_f64, |s, v| s.max(v.abs()));
        if scale == 0.0 {
            return Err(invalid("acoustic reflection plane needs a nonzero normal"));
        }
        let scaled = inward_normal.map(|v| v / scale);
        let length = det::hypot(det::hypot(scaled[0], scaled[1]), scaled[2]);
        Ok(Self {
            point_m,
            normal: scaled.map(|v| v / length),
            boundary,
        })
    }

    pub const fn point_m(&self) -> [f64; 3] {
        self.point_m
    }
    pub const fn inward_normal(&self) -> [f64; 3] {
        self.normal
    }
    pub const fn boundary(&self) -> IdealPlaneBoundary {
        self.boundary
    }

    fn distance(&self, position: [f64; 3]) -> Result<Dd, AeroacError> {
        let mut distance = Dd::ZERO;
        for axis in 0..3 {
            distance = distance
                + (Dd::from_f64(position[axis]) - Dd::from_f64(self.point_m[axis]))
                    * Dd::from_f64(self.normal[axis]);
        }
        if !distance.to_f64().is_finite() {
            return Err(nonfinite("distance to acoustic reflection plane"));
        }
        Ok(distance)
    }

    fn image(&self, source: [f64; 3]) -> Result<[f64; 3], AeroacError> {
        let distance = self.distance(source)?;
        if distance.to_f64() <= 0.0 {
            return Err(invalid(
                "compact sources must lie strictly inside the fluid half-space",
            ));
        }
        let image = core::array::from_fn(|axis| {
            (Dd::from_f64(source[axis]) - distance * Dd::from_f64(2.0 * self.normal[axis])).to_f64()
        });
        if image.iter().any(|v| !v.is_finite()) || self.distance(image)?.to_f64() >= 0.0 {
            return Err(invalid(
                "reflected compact source position is not representable",
            ));
        }
        Ok(image)
    }
}

/// Prepared direct and reflected fields on the original source clock.
///
/// Both paths of each physical source remain coherent, including when distinct
/// sources are declared independent. Preparation folds the dipole reflection
/// into the cached transfer, so evaluation has no allocation, trigonometry or
/// source-vector reflection. The original gas (including molecular relaxation)
/// and all sixteen pressure/velocity coefficients per path are retained.
#[derive(Clone, Debug)]
pub struct PreparedHalfSpaceObserver {
    direct: PreparedCompactSourceObserver,
    reflected: PreparedCompactSourceObserver,
    plane: IdealAcousticPlane,
}

impl PreparedCompactSourceObserver {
    /// Add one exact ideal reflecting plane to this original propagation owner.
    ///
    /// `maximum_source_paths` bounds twice the physical source count, including
    /// the existing direct paths. The plane applies to the whole source array.
    pub fn with_plane(
        self,
        plane: IdealAcousticPlane,
        maximum_source_paths: usize,
    ) -> Result<PreparedHalfSpaceObserver, AeroacError> {
        self.with_plane_cancellable(plane, maximum_source_paths, || false)
    }

    /// Poll at each source during geometry, kernel and reflection preparation.
    pub fn with_plane_cancellable(
        self,
        plane: IdealAcousticPlane,
        maximum_source_paths: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedHalfSpaceObserver, AeroacError> {
        poll(&mut cancelled)?;
        if self
            .source_count()
            .checked_mul(2)
            .is_none_or(|n| n > maximum_source_paths)
        {
            return Err(invalid("half-space source paths exceed preparation budget"));
        }
        if plane.distance(self.observer_m)?.to_f64() < 0.0 {
            return Err(invalid("compact receiver is outside the fluid half-space"));
        }
        let mut positions = Vec::new();
        positions
            .try_reserve_exact(self.source_count())
            .map_err(|_| invalid("reflected source allocation exceeds available memory"))?;
        for source in &self.transfers {
            poll(&mut cancelled)?;
            positions.push(plane.image(source.source_m)?);
        }
        let mut reflected = Self::prepare_cancellable(
            self.acoustic,
            self.frequency_hz,
            self.observer_m,
            &positions,
            self.source_count(),
            self.wavenumber_per_m,
            self.molecular.clone(),
            &mut cancelled,
        )?;
        for transfer in &mut reflected.transfers {
            poll(&mut cancelled)?;
            for row in &mut transfer.field {
                let original = *row;
                row[0] = original[0].scale(plane.boundary.sign());
                for col in 0..3 {
                    let mut sum = ComplexSum::default();
                    for axis in 0..3 {
                        let identity = if axis == col { 1.0 } else { 0.0 };
                        let coefficient = plane.boundary.sign()
                            * (identity - 2.0 * plane.normal[axis] * plane.normal[col]);
                        sum.add(original[axis + 1].scale(coefficient))?;
                    }
                    row[col + 1] = sum.value();
                    if !finite(row[col + 1]) {
                        return Err(nonfinite("reflected compact force transfer"));
                    }
                }
            }
        }
        poll(&mut cancelled)?;
        Ok(PreparedHalfSpaceObserver {
            direct: self,
            reflected,
            plane,
        })
    }
}

impl PreparedHalfSpaceObserver {
    /// Original free-field owner, with unchanged geometry and material identity.
    pub const fn direct(&self) -> &PreparedCompactSourceObserver {
        &self.direct
    }
    pub const fn plane(&self) -> IdealAcousticPlane {
        self.plane
    }
    pub fn source_count(&self) -> usize {
        self.direct.source_count()
    }
    pub const fn scope(&self) -> &'static str {
        if self.direct.molecular.is_some() {
            MOLECULAR_SCOPE
        } else {
            LOSSLESS_SCOPE
        }
    }
    pub fn image_position_m(&self, source: usize) -> Option<[f64; 3]> {
        self.reflected.source_position_m(source)
    }
    /// Reflected path equilibrium travel time; dispersive paths are not pure delays.
    pub fn reflected_delay_s(&self, source: usize) -> Option<f64> {
        self.reflected.source_delay_s(source)
    }
    /// Reflected-path phase delay, already included in its complex transfer.
    pub fn reflected_phase_delay_s(&self, source: usize) -> Option<f64> {
        self.reflected.source_phase_delay_s(source)
    }

    /// Combine all paths before forming signed active and reactive intensity.
    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }

    pub fn evaluate_cancellable(
        &self,
        sources: &[CompactSourcePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(sources, &mut cancelled)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("half-space field intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: self.scope(),
        })
    }

    /// Export to the existing finite-window/ensemble observers, with both path
    /// phases present exactly once. Frequency and source clock remain original.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<crate::observation::HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(sources, &mut || false)?;
        Ok(crate::observation::HarmonicFieldLine {
            frequency_hz: self.direct.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }

    /// Explicit source coherence factors; `max_source_mode_pairs` counts BOTH
    /// physical paths per factor. Reflections never become independent sources.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.evaluate_coherence_modes_cancellable(modes, budget, || false)
    }

    pub fn evaluate_coherence_modes_cancellable(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        if modes.len() > budget.max_modes
            || modes
                .len()
                .checked_mul(self.source_count())
                .and_then(|n| n.checked_mul(2))
                .is_none_or(|n| n > budget.max_source_mode_pairs)
        {
            return Err(invalid("half-space source coherence exceeds work budget"));
        }
        let mut statistics = StatisticsSum::default();
        for mode in modes {
            statistics.add(self.amplitudes(mode, &mut cancelled)?)?;
        }
        poll(&mut cancelled)?;
        statistics.finish(self.scope())
    }

    /// Mutually independent PHYSICAL sources, with each source's direct/image
    /// interference retained before adding its power. Allocation-free O(N).
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.evaluate_independent_sources_cancellable(sources, || false)
    }

    /// Bounded cancellation without changing the prepared geometry or material.
    pub fn evaluate_independent_sources_cancellable(
        &self,
        sources: &[CompactSourcePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        self.direct.check_count(sources)?;
        let mut statistics = StatisticsSum::default();
        for ((direct, image), source) in self
            .direct
            .transfers
            .iter()
            .zip(&self.reflected.transfers)
            .zip(sources)
        {
            poll(&mut cancelled)?;
            let mut sums = [ComplexSum::default(); 4];
            direct.accumulate(*source, &mut sums)?;
            image.accumulate(*source, &mut sums)?;
            statistics.add(sums.map(ComplexSum::value))?;
        }
        poll(&mut cancelled)?;
        statistics.finish(self.scope())
    }

    fn amplitudes(
        &self,
        sources: &[CompactSourcePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[C64; 4], AeroacError> {
        poll(cancelled)?;
        self.direct.check_count(sources)?;
        let mut sums = [ComplexSum::default(); 4];
        for ((direct, image), source) in self
            .direct
            .transfers
            .iter()
            .zip(&self.reflected.transfers)
            .zip(sources)
        {
            poll(cancelled)?;
            direct.accumulate(*source, &mut sums)?;
            image.accumulate(*source, &mut sums)?;
        }
        let values = sums.map(ComplexSum::value);
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("accumulated half-space field"));
        }
        poll(cancelled)?;
        Ok(values)
    }
}

pub(super) fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AeroacError> {
    if cancelled() {
        Err(AeroacError::Cancelled {
            what: "compact acoustic field",
        })
    } else {
        Ok(())
    }
}
