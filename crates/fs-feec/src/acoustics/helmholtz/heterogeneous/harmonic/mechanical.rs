//! Complete pressure and causal mechanical-wall harmonic fields.
//!
//! Every original wall displacement, velocity and Maxwell memory is retained.
//! Reciprocal work uses the same original consistent face mass in both fluid
//! and wall equations. No reciprocal wall impedance or modal reduction occurs.

use super::super::time::{
    ResistiveBoundary,
    mechanical::{MechanicalBoundary, traction::WallTractionProfile},
};
use super::{
    BoundaryAdmittance, BoundaryVelocitySource, C64, Coo, Csr, HeterogeneousError,
    HeterogeneousHelmholtzFe, HeterogeneousPressureField, MAX_HARMONIC_DOFS, PreparedHeterogeneous,
    PressureReleaseConstraints, ResonanceSurvey, ResonanceSurveyControls, ResonanceSurveyError,
    VelocityTrace, VolumeVelocitySource, frequency, invalid, poll,
};
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use frequency::coupled::{CoupledPencil, Wall, WallTraction};
pub use frequency::coupled::{
    MechanicalHarmonicLimits, MechanicalWallAccounting, MechanicalWallFrequencyDerivative,
    MechanicalWallPhasors, MechanicalWallPower,
};
use fs_math::dd::Dd;
use fs_phs::RelaxationBranch;
use fs_phs::harmonic::resonance::survey_resonances;

mod driven;
pub use driven::{
    AntiresonanceControls, AntiresonanceRecovery, DrivenAntiresonance, DrivenResonance,
    MechanicalPressureObservation, MechanicalResonanceDrive, MechanicalResonanceObservation,
    ResonanceControls, ResonanceSearchError,
};

mod poles;
pub use poles::{
    MechanicalPoleError, MechanicalPoleLimits, MechanicalPoleOptions,
    MechanicalPoleResponseOptions, MechanicalPoleWall, MechanicalRelaxationParameterDirection,
    MechanicalResonanceParameterDerivative, MechanicalResonancePole, MechanicalResonanceRegion,
    MechanicalResonanceResidue, MechanicalResonanceResponseError, MechanicalResonanceSpectrum,
    MechanicalSpectrumOptions, MechanicalWallParameterDerivative, MechanicalWallParameterDirection,
    MechanicalWallResidue, PreparedMechanicalResonance,
};

/// Original acoustic pressure together with every distributed causal wall state.
#[derive(Clone, Debug)]
pub struct MechanicalHarmonicField {
    /// Original source work, fluid storage and real-Robin heat. Its backward
    /// error covers the COMPLETE coupled system. Its pressure-only power
    /// residual excludes mechanical transfer; use this wrapper's balances.
    pub acoustic: HeterogeneousPressureField,
    pub walls: Vec<MechanicalWallPhasors>,
    pub wall_accounting: MechanicalWallAccounting,
}

impl MechanicalHarmonicField {
    #[must_use]
    pub fn total_energy_j(&self) -> f64 {
        self.acoustic.kinetic_energy_j
            + self.acoustic.compressional_energy_j
            + self.wall_accounting.total_energy_j()
    }
    /// Fluid and applied mechanical work minus real-Robin and wall heat.
    #[must_use]
    pub fn power_balance_residual_w(&self) -> f64 {
        self.acoustic.input_power_w + self.wall_accounting.wall_traction_power_w
            - self.acoustic.boundary_dissipation_w
            - self.wall_accounting.dissipation_w()
    }
    /// Source work minus Robin heat and signed fluid-to-wall transfer.
    #[must_use]
    pub fn fluid_power_balance_residual_w(&self) -> f64 {
        self.acoustic.input_power_w
            - self.acoustic.boundary_dissipation_w
            - self.wall_accounting.fluid_to_wall_power_w
    }
}

/// Complete response and analytic derivative per fundamental angular frequency.
#[derive(Clone, Debug)]
pub struct MechanicalHarmonicFrequencyResponse {
    pub response: MechanicalHarmonicField,
    pub pressure_derivative_pa_per_rad_s: Vec<C64>,
    pub wall_derivatives: Vec<MechanicalWallFrequencyDerivative>,
    pub derivative_backward_error: f64,
}

/// One retained LU of the complete original pressure/wall/memory equations.
pub struct PreparedMechanicalHarmonic<'m, 'f, 'c> {
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    frequency_hz: f64,
    omega: f64,
    pressure_release: PressureReleaseConstraints,
    conductance: Csr,
    boundaries: Vec<MechanicalBoundary>,
    pencil: CoupledPencil,
}

/// One reusable applied traction shape on one complete harmonic wall owner.
/// Peak amplitudes are Pa, positive outward like interior pressure. The cached
/// consistent force uses the ORIGINAL hierarchical trace, with no area or
/// impedance normalization. Repeated solves reuse both this column and the LU.
pub struct PreparedMechanicalHarmonicTraction<'p, 'm, 'f, 'c> {
    owner: &'p PreparedMechanicalHarmonic<'m, 'f, 'c>,
    wall: usize,
    coefficients: Vec<f64>,
    column: WallTraction,
}

impl PreparedMechanicalHarmonicTraction<'_, '_, '_, '_> {
    #[must_use]
    pub const fn wall(&self) -> usize {
        self.wall
    }

    /// Dimensionless coefficients in the owner's original `wall_dofs` order.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.coefficients
    }
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Prepare full harmonic pressure, wall motion and every original memory.
    ///
    /// Original uniform memoryless Robin, mechanical and pressure-release
    /// faces must be disjoint; coherent velocity forcing may share any
    /// non-release face. Every active original wall trace coordinate survives,
    /// including coordinates adjacent to a pressure-release face. Adjacent
    /// mechanical patches retain their separate local wall states.
    ///
    /// A lossless LOCAL wall-impedance zero causes no division or refusal by
    /// itself. A genuinely singular COUPLED operator still fails its native LU
    /// or residual gate. No damping is added and no dark wall mode is removed.
    /// Dense limits count ALL free pressure, displacement, velocity and Maxwell
    /// coordinates, with the existing MAX_HARMONIC_DOFS hard cap. The original
    /// sparse operators and complete lifted fields are additional storage.
    ///
    /// Bulk medium/geometry and original wall coefficients remain fixed. This
    /// is the existing locally reacting wall law, without lateral bending,
    /// exterior radiation, nonmatching mortar or inferred material dispersion.
    #[allow(clippy::too_many_arguments)]
    pub fn prepare_with_coupled_mechanical_walls<'m>(
        &'m self,
        frequency_hz: f64,
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        mut limits: MechanicalHarmonicLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalHarmonic<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        if !omega.is_finite()
            || omega <= 0.0
            || !(omega * omega).is_finite()
            || omega * omega == 0.0
            || limits.max_states == 0
            || limits.max_dense_bytes == 0
            || !limits.residual_tolerance.is_finite()
            || limits.residual_tolerance <= 0.0
        {
            return Err(invalid("positive complete harmonic frequency and limits required").into());
        }
        limits.max_states = limits.max_states.min(MAX_HARMONIC_DOFS);
        if pressure_release_faces.is_empty() && self.ndof() > limits.max_states {
            return Err(HeterogeneousError::Limit {
                what: "complete pressure and mechanical wall states",
                requested: self.ndof(),
                maximum: limits.max_states,
            });
        }
        let pressure_release = self.geometry.pressure_release_constraints(
            pressure_release_faces,
            if boundaries.is_empty() {
                EmptyFreePressure::Refuse
            } else {
                EmptyFreePressure::AdmitWithWallStates
            },
            cancelled,
        )?;
        let states = pressure_release.free_dofs.len();
        if states > limits.max_states {
            return Err(HeterogeneousError::Limit {
                what: "complete pressure and mechanical wall states",
                requested: states,
                maximum: limits.max_states,
            });
        }
        let mut seen = pressure_release.faces.clone();
        for faces in resistive
            .iter()
            .map(|b| &b.faces)
            .chain(boundaries.iter().map(|b| &b.faces))
        {
            for &face in faces {
                poll(cancelled)?;
                let mut key = face;
                key.sort_unstable();
                if !seen.insert(key) {
                    return Err(invalid(
                        "mechanical, real Robin and pressure-release faces must be disjoint",
                    )
                    .into());
                }
            }
        }
        let robin: Vec<_> = resistive
            .iter()
            .map(|b| BoundaryAdmittance {
                faces: b.faces.clone(),
                admittance_m_pa_s: C64::from_re(b.admittance_m_pa_s),
            })
            .collect();
        let (conductance, susceptance) = self.admittance_matrices(&robin, cancelled)?;
        let (walls, _) = mechanical_walls_for(
            self.geometry,
            boundaries,
            states,
            limits.max_states,
            cancelled,
        )?;
        let pencil = CoupledPencil::prepare(
            omega,
            [&self.stiffness, &self.mass, &conductance, &susceptance],
            &pressure_release.free_dofs,
            walls,
            limits,
            cancelled,
        )?;
        Ok(PreparedMechanicalHarmonic {
            model: self,
            frequency_hz,
            omega,
            pressure_release,
            conductance,
            boundaries: boundaries.to_vec(),
            pencil,
        })
    }

    /// Discover grid-resolved driven pressure peaks of the COMPLETE coupled
    /// wall/pressure operator, including through a regular local Z=0 point.
    /// Original Robin/wall laws and coherent source phasors remain fixed;
    /// the observation is an original tetrahedron and barycentric point.
    /// The shared search owns refinement and budgets. Every query retains
    /// all physical coordinates and reuses its LU for the analytic tangent.
    /// Genuine coupled singularities and cancellation refuse transactionally.
    /// No complete free-mode census, half-power width or decay Q is claimed.
    #[allow(clippy::too_many_arguments)]
    pub fn survey_pressure_resonances_with_coupled_mechanical_walls(
        &self,
        frequencies_rad_s: &[f64],
        resistive: &[ResistiveBoundary],
        boundaries: &[MechanicalBoundary],
        pressure_release_faces: &[[u32; 3]],
        inward_velocity: &[VelocityTrace],
        observation: (usize, [f64; 4]),
        limits: MechanicalHarmonicLimits,
        controls: ResonanceSurveyControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResonanceSurvey, ResonanceSurveyError<HeterogeneousError>> {
        poll(cancelled).map_err(ResonanceSearchError::Source)?;
        self.validate_sample(self.ndof(), observation.0, observation.1)
            .map_err(HeterogeneousError::from)
            .map_err(ResonanceSearchError::Source)?;
        survey_resonances(
            frequencies_rad_s,
            controls,
            cancelled,
            |omega, _, cancelled| {
                let result = self
                    .prepare_with_coupled_mechanical_walls(
                        omega / core::f64::consts::TAU,
                        resistive,
                        boundaries,
                        pressure_release_faces,
                        limits,
                        cancelled,
                    )?
                    .solve_with_frequency_derivative(inward_velocity, None, cancelled)?;
                Ok((
                    self.pressure_at(
                        &result.response.acoustic.pressure_pa,
                        observation.0,
                        observation.1,
                    )?,
                    self.pressure_at(
                        &result.pressure_derivative_pa_per_rad_s,
                        observation.0,
                        observation.1,
                    )?,
                    result
                        .response
                        .acoustic
                        .backward_error
                        .max(result.derivative_backward_error),
                ))
            },
        )
    }
}

/// The original distributed wall trace and material-memory preparation, shared
/// by real and complex bulk fluids. Pressure release does not remove wall DOFs.
pub(super) fn mechanical_walls_for(
    geometry: &super::super::super::HelmholtzFe<'_>,
    boundaries: &[MechanicalBoundary],
    mut states: usize,
    max_states: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(Vec<Wall>, usize), HeterogeneousError> {
    let mut walls = Vec::with_capacity(boundaries.len());
    for boundary in boundaries {
        poll(cancelled)?;
        fs_phs::wall_specific_impedance(&boundary.wall, 1.0).map_err(HeterogeneousError::Wall)?;
        let mut arms = Vec::with_capacity(boundary.relaxation.len());
        for arm in &boundary.relaxation {
            poll(cancelled)?;
            let original = RelaxationBranch {
                projection: vec![1.0],
                stiffness: arm.stiffness_per_area,
                relaxation_time_s: arm.relaxation_time_s,
            };
            let coupling = original
                .scalar_relaxed_memory(1.0)
                .map_err(HeterogeneousError::Wall)?;
            arms.push((coupling, 1.0 / arm.relaxation_time_s));
        }
        // The same high-order/curved face Gram matrix as the time-domain
        // owner, without face averaging or an inverse wall impedance.
        let (global, _) = HeterogeneousHelmholtzFe::boundary_coefficient_matrices_for(
            geometry,
            [(boundary.faces.as_slice(), C64::ONE)],
            true,
            cancelled,
        )?;
        let mut dofs = Vec::new();
        for row in 0..geometry.ndof() {
            poll(cancelled)?;
            let (columns, values) = global.row(row);
            if columns.binary_search(&row).is_ok_and(|i| values[i] > 0.0) {
                dofs.push(row);
            }
        }
        let count = arms
            .len()
            .checked_add(2)
            .and_then(|n| dofs.len().checked_mul(n))
            .ok_or_else(|| invalid("complete wall state extent overflow"))?;
        states = states
            .checked_add(count)
            .ok_or_else(|| invalid("complete harmonic state extent overflow"))?;
        if states > max_states {
            return Err(HeterogeneousError::Limit {
                what: "complete pressure and mechanical wall states",
                requested: states,
                maximum: max_states,
            });
        }
        let mut local = Coo::new(dofs.len(), dofs.len());
        let mut area = Dd::ZERO;
        for (row, &global_row) in dofs.iter().enumerate() {
            poll(cancelled)?;
            let (cols, values) = global.row(global_row);
            for (&col, &value) in cols.iter().zip(values) {
                if value != 0.0 {
                    // The original hierarchical H1 constant trace has one on
                    // vertex DOFs and zero on all edge/face bubbles. Thus its
                    // quadratic form is the physical area, including curved
                    // face weights; summing every matrix entry would be wrong.
                    if global_row < geometry.space.edge_off && col < geometry.space.edge_off {
                        area = area + Dd::from_f64(value);
                    }
                    let col = dofs
                        .binary_search(&col)
                        .map_err(|_| invalid("incomplete original mechanical trace support"))?;
                    local.push(row, col, value);
                }
            }
        }
        let area_m2 = area.to_f64();
        let mass_kg = area_m2 * boundary.wall.surface_density;
        if !area_m2.is_finite() || area_m2 <= 0.0 || !mass_kg.is_finite() || mass_kg <= 0.0 {
            return Err(
                invalid("original mechanical wall area and mass must be representable").into(),
            );
        }
        walls.push(Wall {
            dofs,
            mass: local.assemble(),
            area_m2,
            surface_density: boundary.wall.surface_density,
            stiffness: boundary.wall.stiffness_per_area,
            resistance: boundary.wall.resistance,
            arms,
        });
    }
    Ok((walls, states))
}

impl<'m, 'f, 'c> PreparedMechanicalHarmonic<'m, 'f, 'c> {
    #[must_use]
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    #[must_use]
    pub fn wall_count(&self) -> usize {
        self.boundaries.len()
    }
    #[must_use]
    pub fn wall(&self, wall: usize) -> Option<&MechanicalBoundary> {
        self.boundaries.get(wall)
    }
    #[must_use]
    pub fn wall_dofs(&self, wall: usize) -> Option<&[usize]> {
        self.pencil.wall_dofs(wall)
    }

    /// Independent positive resistance/Maxwell absorption on each original wall.
    /// The caller binds these coherent phasors to this owner and frequency.
    /// Exact native relaxation extensions are required; pressure/traction work
    /// and mutable aggregate heat receipts are not inputs to this calculation.
    pub fn wall_dissipation(
        &self,
        walls: &[MechanicalWallPhasors],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<MechanicalWallPower>, HeterogeneousError> {
        Ok(self.pencil.wall_dissipation(walls, cancelled)?)
    }

    /// Project an original wall shape once. All mechanical coordinates survive
    /// pressure release on an adjacent fluid face. Storage is O(wall DOFs).
    pub fn prepare_wall_traction(
        &self,
        profile: &WallTractionProfile,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedMechanicalHarmonicTraction<'_, 'm, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let dofs = self
            .wall_dofs(profile.wall)
            .ok_or_else(|| invalid("traction wall index"))?;
        let coefficients = profile.local_coefficients(dofs, cancelled)?;
        let column = self
            .pencil
            .prepare_traction(profile.wall, &coefficients, cancelled)?;
        Ok(PreparedMechanicalHarmonicTraction {
            owner: self,
            wall: profile.wall,
            coefficients,
            column,
        })
    }

    fn traction_columns<'p>(
        &self,
        tractions: &[(&'p PreparedMechanicalHarmonicTraction<'_, 'm, 'f, 'c>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<(&'p WallTraction, C64)>, HeterogeneousError> {
        poll(cancelled)?;
        super::super::source::admit_count(tractions.len())?;
        tractions
            .iter()
            .map(|&(source, amplitude)| {
                poll(cancelled)?;
                if !core::ptr::eq(source.owner, self)
                    || !amplitude.re.is_finite()
                    || !amplitude.im.is_finite()
                {
                    return Err(invalid(
                        "finite traction amplitude and original prepared harmonic owner required",
                    )
                    .into());
                }
                Ok((&source.column, amplitude))
            })
            .collect()
    }

    pub fn solve(
        &self,
        inward_velocity: &[VelocityTrace],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalHarmonicField, HeterogeneousError> {
        self.solve_with_boundary_sources(inward_velocity, &[], &[], cancelled)
    }

    /// Coherent original linear traces, prepared boundary motion and volume
    /// injection. Source units, binding, release-face refusal and signed source
    /// work are exactly those of the original full-pressure harmonic owner.
    pub fn solve_with_boundary_sources(
        &self,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalHarmonicField, HeterogeneousError> {
        self.solve_with_wall_tractions(
            inward_velocity,
            boundary_sources,
            volume_sources,
            &[],
            cancelled,
        )
    }

    /// Coherent fluid sources and independently applied wall tractions [Pa].
    /// Tractions force momentum through W*shape, not prescribed wall motion.
    /// Fluid input power, applied mechanical power and internal pressure work
    /// remain separate signed accounts; no incoherent source-power sum is used.
    pub fn solve_with_wall_tractions(
        &self,
        inward_velocity: &[VelocityTrace],
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64)],
        tractions: &[(&PreparedMechanicalHarmonicTraction<'_, 'm, 'f, 'c>, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalHarmonicField, HeterogeneousError> {
        let columns = self.traction_columns(tractions, cancelled)?;
        let load = PreparedHeterogeneous::project_sources_for(
            self.model,
            &self.pressure_release,
            inward_velocity,
            boundary_sources,
            volume_sources,
            cancelled,
        )?;
        let (solution, traction_power) = self
            .pencil
            .solve_with_tractions(&load, &columns, cancelled)?;
        let mut field = self.field(solution, &load, cancelled)?;
        field.wall_accounting.wall_traction_power_w = traction_power;
        PreparedHeterogeneous::source_work(
            &mut field.acoustic,
            boundary_sources,
            volume_sources,
            cancelled,
        )?;
        Ok(field)
    }

    /// Complete analytic tangent per rad/s using the SAME full-state LU.
    /// Optional trace slopes retain coherent source phase changes. Neither a
    /// wall impedance inverse nor a finite-difference frequency step is used.
    pub fn solve_with_frequency_derivative(
        &self,
        inward_velocity: &[VelocityTrace],
        inward_velocity_derivative: Option<&[VelocityTrace]>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalHarmonicFrequencyResponse, HeterogeneousError> {
        self.solve_with_sources_frequency_derivative(
            inward_velocity,
            inward_velocity_derivative,
            &[],
            &[],
            cancelled,
        )
    }

    /// Every source triple is (original source, amplitude, d_amplitude/domega).
    /// Derivatives retain every pressure, displacement, velocity and original
    /// Maxwell memory coordinate. Exact pressure-release derivatives stay zero.
    pub fn solve_with_sources_frequency_derivative(
        &self,
        inward_velocity: &[VelocityTrace],
        inward_velocity_derivative: Option<&[VelocityTrace]>,
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64, C64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalHarmonicFrequencyResponse, HeterogeneousError> {
        self.solve_with_wall_tractions_frequency_derivative(
            inward_velocity,
            inward_velocity_derivative,
            boundary_sources,
            volume_sources,
            &[],
            cancelled,
        )
    }

    /// Each source triple is (original source, peak amplitude, d/domega).
    /// Traction derivatives are Pa/(rad/s), supplied independently of the
    /// acoustic m/s or m^3/s spectra. The exact RHS derivative and complete
    /// operator tangent use the same LU, including every Maxwell coordinate.
    pub fn solve_with_wall_tractions_frequency_derivative(
        &self,
        inward_velocity: &[VelocityTrace],
        inward_velocity_derivative: Option<&[VelocityTrace]>,
        boundary_sources: &[(&BoundaryVelocitySource<'_, '_, '_>, C64, C64)],
        volume_sources: &[(&VolumeVelocitySource<'_, '_, '_>, C64, C64)],
        tractions: &[(
            &PreparedMechanicalHarmonicTraction<'_, 'm, 'f, 'c>,
            C64,
            C64,
        )],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalHarmonicFrequencyResponse, HeterogeneousError> {
        poll(cancelled)?;
        super::super::boundary_source::admit_count(boundary_sources.len())?;
        super::super::source::admit_count(volume_sources.len())?;
        super::super::source::admit_count(tractions.len())?;
        let (traction_values, traction_slopes): (Vec<_>, Vec<_>) = tractions
            .iter()
            .map(|&(source, value, slope)| ((source, value), (source, slope)))
            .unzip();
        let columns = self.traction_columns(&traction_values, cancelled)?;
        let column_slopes = self.traction_columns(&traction_slopes, cancelled)?;
        let (boundary_values, boundary_slopes): (Vec<_>, Vec<_>) = boundary_sources
            .iter()
            .map(|&(source, value, slope)| ((source, value), (source, slope)))
            .unzip();
        let (volume_values, volume_slopes): (Vec<_>, Vec<_>) = volume_sources
            .iter()
            .map(|&(source, value, slope)| ((source, value), (source, slope)))
            .unzip();
        let load = PreparedHeterogeneous::project_sources_for(
            self.model,
            &self.pressure_release,
            inward_velocity,
            &boundary_values,
            &volume_values,
            cancelled,
        )?;
        let slope = PreparedHeterogeneous::project_sources_for(
            self.model,
            &self.pressure_release,
            inward_velocity_derivative.unwrap_or(&[]),
            &boundary_slopes,
            &volume_slopes,
            cancelled,
        )?;
        let (primary, traction_power) = self
            .pencil
            .solve_with_tractions(&load, &columns, cancelled)?;
        let mut derivative = self.pencil.frequency_derivative_with_tractions(
            &primary.pressure,
            &load,
            &slope,
            &column_slopes,
            cancelled,
        )?;
        let wall_derivatives =
            self.pencil
                .wall_derivatives(&primary.pressure, &derivative.pressure, cancelled)?;
        derivative.pressure.truncate(self.model.ndof());
        let mut response = self.field(primary, &load, cancelled)?;
        response.wall_accounting.wall_traction_power_w = traction_power;
        PreparedHeterogeneous::source_work(
            &mut response.acoustic,
            &boundary_values,
            &volume_values,
            cancelled,
        )?;
        Ok(MechanicalHarmonicFrequencyResponse {
            response,
            pressure_derivative_pa_per_rad_s: derivative.pressure,
            wall_derivatives,
            derivative_backward_error: derivative.backward_error,
        })
    }

    fn field(
        &self,
        mut solution: frequency::PressureSolution,
        load: &[C64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MechanicalHarmonicField, HeterogeneousError> {
        let (walls, wall_accounting) = self.pencil.wall_fields(&solution.pressure, cancelled)?;
        solution.pressure.truncate(self.model.ndof());
        let acoustic = PreparedHeterogeneous::field_for(
            self.model,
            self.omega,
            &self.conductance,
            solution.pressure,
            load,
            solution.backward_error,
            cancelled,
        )?;
        Ok(MechanicalHarmonicField {
            acoustic,
            walls,
            wall_accounting,
        })
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod traction_tests;
