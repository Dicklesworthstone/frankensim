//! Complete driven cavity response with causal Maxwell material/wall memory.
//!
//! The existing mechanical pHS owner retains q=[u,psi], momentum M*q_dot
//! and every viscous-strain state. The shared harmonic owner factors that full
//! system with peak exp(-i*omega*t) phasors. K is RELAXED. Branch stiffness
//! must not already be included in K; there is no sampled damping surrogate.

use super::StructuralInterior;
pub use super::resonance::{
    StructuralInteriorParameterDirection, StructuralInteriorSourceParameterJet,
};
use super::resonance::{StructuralInteriorReceiver, StructuralInteriorSources};
use fs_exec::Cx;
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::HarmonicTransferSample;
pub use fs_modal::relaxation_modes::MaxwellDampingBranch;
pub use fs_phs::harmonic::resonance::{
    DrivenResonance, ResonanceControls, ResonanceSearchError, ResonanceSurvey,
    ResonanceSurveyControls, ResonanceSurveyError,
};
use fs_phs::harmonic::resonance::{search_resonance, survey_resonances};
pub use fs_phs::mechanical::MechanicalRelaxationDirection as StructuralInteriorRelaxationBranchDirection;

mod fields;
mod power;
mod sources;
mod transfer;
use fields::StructuralInteriorRelaxationHarmonicError as Error;
pub use fields::{
    StructuralInteriorRelaxationCurvatureResponse, StructuralInteriorRelaxationFixedReceiver,
    StructuralInteriorRelaxationFrequencyParameterResponse,
    StructuralInteriorRelaxationFrequencyResponse, StructuralInteriorRelaxationHarmonicError,
    StructuralInteriorRelaxationHarmonicLimits, StructuralInteriorRelaxationMemoryCurvature,
    StructuralInteriorRelaxationMemoryFrequency, StructuralInteriorRelaxationMemoryParameterJet,
    StructuralInteriorRelaxationMemoryResponse, StructuralInteriorRelaxationParameterDirection,
    StructuralInteriorRelaxationPowerCurvature, StructuralInteriorRelaxationPowerFrequency,
    StructuralInteriorRelaxationPowerParameterJet, StructuralInteriorRelaxationReceiverJet,
    StructuralInteriorRelaxationReceiverSpectrum, StructuralInteriorRelaxationResponse,
};
pub use power::StructuralInteriorRelaxationPowerChannel;
pub use sources::{
    StructuralInteriorRelaxationSourceOrder, StructuralInteriorRelaxationSourceSpectrum,
};

/// Source-bound complete linear image; all frequencies share native assembly.
pub struct PreparedStructuralInteriorRelaxation<'m, 'a, 'f, 'c> {
    model: &'m StructuralInterior<'a, 'f, 'c>,
    inner: fields::Prepared<'m>,
}

/// One original full-state factor, reusable across coherent force/flow drives.
pub struct PreparedStructuralInteriorRelaxationHarmonic<'p, 'm> {
    inner: fields::Factor<'p, 'm>,
}

/// One physical parameter direction, prepared once for reuse across frequencies.
/// It is bound to the original immutable mechanical/cavity preparation.
pub struct PreparedStructuralInteriorRelaxationDirection<'p, 'm> {
    inner: fields::PreparedDirection<'p, 'm>,
}

impl<'m, 'a, 'f, 'c> PreparedStructuralInteriorRelaxation<'m, 'a, 'f, 'c> {
    /// Prepare the complete original M/C/relaxed-K, reciprocal B and Maxwell
    /// direction using the already admitted mechanical inverse and arm map.
    ///
    /// One supplied B derivative enters both opposite spatial blocks and both
    /// force/flow observations. Each arm supplies its stiffness, relaxation-time
    /// and original strain-projection derivative; equal-time arms remain
    /// separate. Source derivatives belong to the later response query.
    ///
    /// No factorization, material fit or finite model perturbation occurs.
    /// The same prepared direction can be reused by factors of this model at
    /// any admitted positive frequency. A factor from another preparation
    /// refuses it. The existing adapter limits bound complete retained
    /// directions, simultaneous native jets and independent physical recovery.
    pub fn prepare_parameter_direction<'p>(
        &'p self,
        direction: StructuralInteriorRelaxationParameterDirection<'_>,
        cx: &Cx<'_>,
    ) -> Result<PreparedStructuralInteriorRelaxationDirection<'p, 'm>, Error> {
        Ok(PreparedStructuralInteriorRelaxationDirection {
            inner: self
                .inner
                .prepare_parameter_direction(direction, &mut || cx.checkpoint().is_err())?,
        })
    }
}

impl<'a, 'f, 'c> StructuralInterior<'a, 'f, 'c> {
    /// Bind original material channels before solving any driven frequency.
    ///
    /// Branch strains follow [structural displacement, pressure potential].
    /// A potential-only channel realizes k*tau/(1+s*tau)*b*b^T wall admittance;
    /// structural-only and mixed channels retain their declared constitutive law.
    /// The caller supplies that law and its valid material/geometry regime.
    /// Every original coordinate, including constant potential, survives.
    ///
    /// Native mechanical admission bounds 2*n+b states and performs no modal
    /// truncation. Adapter bounds precede physical assembly; the native harmonic
    /// owner separately admits its complete image and factors. Cancellation is
    /// checked around bounded mechanical assembly and throughout
    /// harmonic factorization, response, original-row and physical work recovery.
    pub fn prepare_relaxation_response<'m>(
        &'m self,
        branches: &'m [MaxwellDampingBranch],
        limits: StructuralInteriorRelaxationHarmonicLimits,
        cx: &Cx<'_>,
    ) -> Result<PreparedStructuralInteriorRelaxation<'m, 'a, 'f, 'c>, Error> {
        let inner = fields::Prepared::new(
            fields::Source {
                structural: self.structural_dofs,
                coordinates: self.structural_dofs + self.pressure_dofs(),
                mass: &self.mass,
                damping: &self.damping,
                stiffness: &self.stiffness,
                coupling: &self.coupling,
            },
            branches,
            limits,
            &mut || cx.checkpoint().is_err(),
        )?;
        Ok(PreparedStructuralInteriorRelaxation { model: self, inner })
    }
}

impl PreparedStructuralInteriorRelaxation<'_, '_, '_, '_> {
    /// Factor at one strictly positive angular frequency [rad/s].
    /// A true singular coupled response refuses; no loss floor is introduced.
    pub fn factor(
        &self,
        angular_frequency_rad_s: f64,
        cx: &Cx<'_>,
    ) -> Result<PreparedStructuralInteriorRelaxationHarmonic<'_, '_>, Error> {
        Ok(PreparedStructuralInteriorRelaxationHarmonic {
            inner: self
                .inner
                .factor(angular_frequency_rad_s, &mut || cx.checkpoint().is_err())?,
        })
    }

    /// Resolve a selected receiver's local peak and, optionally, its actual
    /// squared-amplitude half-height crossings. These depend on coherent drive
    /// and observation; their bandwidth Q is not a pole-decay Q. The interval
    /// must isolate the requested peak. Shared search owns all query/refinement
    /// limits and refuses missing crossings, singularities and cancellation.
    pub fn driven_resonance(
        &self,
        interval_rad_s: (f64, f64),
        sources: StructuralInteriorSources<'_>,
        receiver: StructuralInteriorReceiver<'_>,
        controls: ResonanceControls,
        cx: &Cx<'_>,
    ) -> Result<DrivenResonance, ResonanceSearchError<Error>> {
        search_resonance(
            interval_rad_s,
            controls,
            || cx.checkpoint().is_err(),
            |omega, _, _| self.sample(omega, sources, receiver, cx),
        )
    }

    /// Survey the driven peaks resolved by the supplied positive rad/s grid.
    /// This is not a complete free-mode census. Each query reuses assembly and
    /// one native factor for the complete response and its analytic tangent.
    pub fn survey_resonances(
        &self,
        frequencies_rad_s: &[f64],
        sources: StructuralInteriorSources<'_>,
        receiver: StructuralInteriorReceiver<'_>,
        controls: ResonanceSurveyControls,
        cx: &Cx<'_>,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<Error>> {
        survey_resonances(
            frequencies_rad_s,
            controls,
            || cx.checkpoint().is_err(),
            |omega, _, _| self.sample(omega, sources, receiver, cx),
        )
    }

    fn sample(
        &self,
        omega: f64,
        sources: StructuralInteriorSources<'_>,
        receiver: StructuralInteriorReceiver<'_>,
        cx: &Cx<'_>,
    ) -> Result<(C64, C64, f64), Error> {
        let result = self.factor(omega, cx)?.solve_with_frequency_derivative(
            sources.structural_force,
            sources.inward_flow_load,
            None,
            cx,
        )?;
        let sample = match receiver {
            StructuralInteriorReceiver::StructuralDisplacement(weights) => {
                HarmonicTransferSample::project_real(
                    weights,
                    &result.response.structural_displacement,
                    &result.structural_displacement_derivative_per_rad_s,
                )
                .map_err(|_| Error::Invalid("invalid original displacement receiver"))?
            }
            StructuralInteriorReceiver::StructuralVelocity(weights) => {
                HarmonicTransferSample::project_real(
                    weights,
                    &result.response.structural_velocity,
                    &result.structural_velocity_derivative_per_rad_s,
                )
                .map_err(|_| Error::Invalid("invalid original velocity receiver"))?
            }
            StructuralInteriorReceiver::PressurePoint {
                tetrahedron,
                barycentric,
            } => HarmonicTransferSample {
                amplitude: self
                    .model
                    .interior
                    .pressure_at(&result.response.pressure_pa, tetrahedron, barycentric)
                    .map_err(|_| Error::Invalid("invalid original pressure receiver"))?,
                angular_derivative: self
                    .model
                    .interior
                    .pressure_at(
                        &result.pressure_derivative_pa_per_rad_s,
                        tetrahedron,
                        barycentric,
                    )
                    .map_err(|_| Error::Invalid("invalid original pressure receiver"))?,
            },
        };
        Ok((
            sample.amplitude,
            sample.angular_derivative,
            result
                .response
                .relative_residual
                .max(result.relative_derivative_residual),
        ))
    }
}

impl PreparedStructuralInteriorRelaxationHarmonic<'_, '_> {
    /// Apply original peak generalized force and inward pressure-basis flow.
    /// This is the complete steady driven solution, including all backgrounds
    /// and memory channels. Initial transient history is not an extra AC source.
    pub fn solve(
        &self,
        structural_force: &[C64],
        inward_flow_load: &[C64],
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorRelaxationResponse, Error> {
        self.inner
            .solve(structural_force, inward_flow_load, &mut || {
                cx.checkpoint().is_err()
            })
    }

    /// Resolve the complete original frequency/parameter response with this
    /// SAME factor, including material k/tau/strain and geometry coupling terms.
    ///
    /// Five native state solves retain every memory state. The original
    /// parameter and mixed momentum/branch equations are checked separately.
    /// Changing sqrt(k) is included when canonical branch effort becomes
    /// physical stress. Independent material heat includes the direct decay-
    /// rate derivative; instantaneous wall/structural loss includes original
    /// damping changes. Interface work retains both B and B_theta contributions.
    ///
    /// Full coherent source amplitude/phase jets are supplied in original
    /// force/inward-flow order. The result exposes fields and signed power jets
    /// for the existing receiver, extremum and bandwidth query owners.
    /// Coordinates, topology and arm count stay fixed. No finite parameter
    /// change, missing constitutive derivative or transient change work is inferred.
    pub fn solve_with_frequency_parameter_derivative(
        &self,
        sources: StructuralInteriorSourceParameterJet<'_>,
        direction: &PreparedStructuralInteriorRelaxationDirection<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorRelaxationFrequencyParameterResponse, Error> {
        self.inner
            .solve_with_frequency_parameter_derivative(sources, &direction.inner, &mut || {
                cx.checkpoint().is_err()
            })
    }

    /// Recover every original field, material spring and heat channel through
    /// its second frequency derivative from the same complete-state factor.
    ///
    /// Both optional source jets are partial derivatives per rad/s at fixed
    /// original M, C, relaxed K, geometry and branch laws. They retain coherent
    /// source phase and amplitude variation. None supplies a zero derivative.
    /// All memory states remain dynamic; no impedance is frozen at this
    /// frequency. Native and original second-order equations are checked
    /// independently, and physical source, interface and heat work are exposed
    /// with their complete frequency curvature.
    pub fn solve_with_frequency_curvature(
        &self,
        structural_force: &[C64],
        inward_flow_load: &[C64],
        source_derivative: Option<StructuralInteriorSources<'_>>,
        source_second_derivative: Option<StructuralInteriorSources<'_>>,
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorRelaxationCurvatureResponse, Error> {
        self.inner.solve_with_frequency_curvature(
            structural_force,
            inward_flow_load,
            source_derivative.map(|d| (d.structural_force, d.inward_flow_load)),
            source_second_derivative.map(|d| (d.structural_force, d.inward_flow_load)),
            &mut || cx.checkpoint().is_err(),
        )
    }

    /// Reuse this full factor for every original field's derivative per rad/s.
    /// Explicit source derivatives cover both force and flow; None fixes both.
    /// Pressure comes from the canonical velocity and includes the full
    /// derivative of -i*omega*psi. No finite-difference solve or pole fit is used.
    pub fn solve_with_frequency_derivative(
        &self,
        structural_force: &[C64],
        inward_flow_load: &[C64],
        source_derivative: Option<StructuralInteriorSources<'_>>,
        cx: &Cx<'_>,
    ) -> Result<StructuralInteriorRelaxationFrequencyResponse, Error> {
        self.inner.solve_with_frequency_derivative(
            structural_force,
            inward_flow_load,
            source_derivative.map(|d| (d.structural_force, d.inward_flow_load)),
            &mut || cx.checkpoint().is_err(),
        )
    }
}
