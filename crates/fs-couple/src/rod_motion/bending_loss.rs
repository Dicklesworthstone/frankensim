//! Source-bound bending viscosity on the native finite-rotation rod trajectory.
//!
//! For each circular section, the existing material owner supplies C_b=sum(eta I).
//! With phi=Log(R_i^-1 R_j), its Rayleigh dissipation is
//! C_b/(2L) * (phi_dot_y^2 + phi_dot_z^2). The shared angular-strain pullback
//! uses the same left/right SO(3) Jacobians as the elastic Cosserat owner.
//! A common rigid spin therefore cannot become material bending loss.
//!
//! Each frozen-pose rank-one flow uses exact exponential relaxation and its
//! positive integrated work. A palindromic sweep composes the two bending
//! directions and adjacent segments. Optional torsional viscosity is a distinct
//! symmetric flow with a separate heat account. Elastic motion still requires
//! timestep refinement; this does not turn coarse splitting into exact dynamics.
//! Sources are frozen, small-strain Kelvin laws. No axial/shear damping,
//! thermal feedback, fitted loss ratio or hereditary spectrum is inferred.

use fs_exec::CancelGate;
use fs_material::section::ConcentricKelvinVoigtSection;
use fs_solid::rod::MAX_MATERIAL_ROD_SEGMENTS;

use super::torsional_loss::{TwistMode, angular_modes, relax_angular, same_elastic_source};
use super::{
    KickAccounts, PreparedRodLoads, PreparedRodMotion, PreparedRodTorsionalLoss, RodMotion,
    RodMotionControls, RodMotionError, RodMotionInertia, RodMotionStep, RodNodalWrench,
    RodTorsionalHeat, StageLoads, checkpoint, finite,
};

/// Explicit original bending-viscosity sections, one per physical rod segment.
#[derive(Clone, Debug)]
pub struct PreparedRodBendingLoss {
    sections: Vec<ConcentricKelvinVoigtSection>,
    resistance_nm_s: Vec<f64>,
    segment_lengths_m: Vec<f64>,
}

/// Physical bending heat, distinct from torsion and the native integration defect.
#[derive(Clone, Debug, PartialEq)]
pub struct RodBendingHeat {
    /// Sum over both bending directions and every original segment [J].
    pub total_heat_j: f64,
    /// Original source-segment order, summing its two bending directions [J].
    pub segment_heat_j: Vec<f64>,
}

/// Independently selected angular material losses on the same physical step.
/// `None` leaves that mechanism unmodeled; it does not assert a lossless material.
#[derive(Clone, Copy, Debug)]
pub struct RodAngularLosses<'a> {
    /// Explicit viscosity of both material bending directions.
    pub bending: Option<&'a PreparedRodBendingLoss>,
    /// Independently sourced circular-shaft shear viscosity.
    pub torsional: Option<&'a PreparedRodTorsionalLoss>,
}

impl PreparedRodBendingLoss {
    /// Bind C_b/L to each original segment without modifying its elastic law.
    ///
    /// Each loss section must retain the original geometry, material identities,
    /// complete query and selected mechanical properties. It may additionally
    /// resolve the explicit solid bending-viscosity property. Every source band
    /// must support [0, f64::MAX] rad/s: a finite-band loss sample cannot authorize
    /// this transient Kelvin law. Explicit zero is supported; missing data is not.
    pub fn new<I>(
        model: &PreparedRodMotion<I>,
        sections: Vec<ConcentricKelvinVoigtSection>,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        checkpoint(gate)?;
        if sections.is_empty()
            || sections.len() > MAX_MATERIAL_ROD_SEGMENTS
            || sections.len() != model.source().sections().len()
        {
            return Err(RodMotionError::Invalid(
                "bending loss needs matching bounded source sections",
            ));
        }
        let segment_lengths_m = model.source().reference_mass().segment_lengths_m().to_vec();
        let mut resistance_nm_s = Vec::with_capacity(sections.len());
        for (section, &segment_length_m) in sections.iter().zip(&segment_lengths_m) {
            checkpoint(gate)?;
            if section.omega_band_rad_s() != (0.0, f64::MAX) {
                return Err(RodMotionError::Invalid(
                    "transient bending viscosity requires explicit broadband source support",
                ));
            }
            let value = section.viscous_stiffness_n_m2_s() / segment_length_m;
            if !value.is_finite()
                || value < 0.0
                || (section.viscous_stiffness_n_m2_s() > 0.0 && value == 0.0)
            {
                return Err(RodMotionError::Invalid(
                    "unrepresentable bending resistance",
                ));
            }
            resistance_nm_s.push(value);
        }
        let loss = Self {
            sections,
            resistance_nm_s,
            segment_lengths_m,
        };
        loss.admit(model, gate)?;
        Ok(loss)
    }

    /// Original circular geometry, source selectors and resolved material states.
    pub fn sections(&self) -> &[ConcentricKelvinVoigtSection] {
        &self.sections
    }

    /// C_b/L multiplying each of the two bending log-strain rates [N m s].
    pub fn resistance_nm_s(&self) -> &[f64] {
        &self.resistance_nm_s
    }

    pub(super) fn admit<I>(
        &self,
        model: &PreparedRodMotion<I>,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        let mechanical = model.source().sections();
        if self.sections.len() != mechanical.len()
            || self.segment_lengths_m != model.source().reference_mass().segment_lengths_m()
        {
            return Err(RodMotionError::Invalid(
                "bending loss must match the reference rod segments",
            ));
        }
        for (source, current) in self.sections.iter().zip(mechanical) {
            checkpoint(gate)?;
            if !same_elastic_source(current.elastic(), source.elastic()) {
                return Err(RodMotionError::Invalid(
                    "bending loss source is stale or differs from the actual rod material state",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn heat(&self) -> RodBendingHeat {
        RodBendingHeat {
            total_heat_j: 0.0,
            segment_heat_j: vec![0.0; self.sections.len()],
        }
    }

    fn flow<I: RodMotionInertia>(
        &self,
        state: &mut RodMotion<I>,
        dt: f64,
        accounts: &mut KickAccounts,
        heat: &mut RodBendingHeat,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        // One bounded stack array, one log/Jacobian evaluation per active pair.
        // No global damping matrix or per-edge allocation in the loss flow.
        let mut modes = [TwistMode::ZERO; 2 * MAX_MATERIAL_ROD_SEGMENTS];
        let count = 2 * self.sections.len();
        for segment in 0..self.sections.len() {
            checkpoint(gate)?;
            if self.resistance_nm_s[segment] != 0.0 {
                let pair = angular_modes(state, segment, [1, 2])?;
                modes[2 * segment..2 * segment + 2].copy_from_slice(&pair);
            }
        }
        let mut relax = |index: usize, interval: f64| -> Result<(), RodMotionError> {
            let segment = index / 2;
            let dissipated = relax_angular(
                state,
                segment,
                modes[index],
                self.resistance_nm_s[segment],
                interval,
                accounts,
                gate,
            )?;
            heat.segment_heat_j[segment] += dissipated;
            heat.total_heat_j += dissipated;
            finite(&[heat.segment_heat_j[segment], heat.total_heat_j])
        };
        for index in 0..count - 1 {
            relax(index, 0.5 * dt)?;
        }
        relax(count - 1, dt)?;
        for index in (0..count - 1).rev() {
            relax(index, 0.5 * dt)?;
        }
        Ok(())
    }
}

impl RodAngularLosses<'_> {
    pub(super) fn admit<I>(
        self,
        model: &PreparedRodMotion<I>,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        if let Some(law) = self.torsional {
            law.admit(model, gate)?;
        }
        if let Some(law) = self.bending {
            law.admit(model, gate)?;
        }
        Ok(())
    }

    pub(super) fn flow<I: RodMotionInertia>(
        self,
        state: &mut RodMotion<I>,
        dt: f64,
        accounts: &mut KickAccounts,
        torsional_heat: &mut Option<RodTorsionalHeat>,
        bending_heat: &mut Option<RodBendingHeat>,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        // T(dt/2) B(dt) T(dt/2) preserves second order when both mechanisms
        // share inertia. An explicitly zero bending law must not subdivide
        // the existing torsional sweep: adjacent pair flows need not commute.
        let active_bending = self
            .bending
            .is_some_and(|law| law.resistance_nm_s.iter().any(|&c| c != 0.0));
        let torsional_dt = if active_bending { 0.5 * dt } else { dt };
        if let (Some(law), Some(heat)) = (self.torsional, torsional_heat.as_mut()) {
            law.flow(state, torsional_dt, accounts, heat, gate)?;
        }
        if active_bending && let (Some(law), Some(heat)) = (self.bending, bending_heat.as_mut()) {
            law.flow(state, dt, accounts, heat, gate)?;
        }
        if active_bending && let (Some(law), Some(heat)) = (self.torsional, torsional_heat.as_mut())
        {
            law.flow(state, torsional_dt, accounts, heat, gate)?;
        }
        Ok(())
    }
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Native motion with separately sourced bending and/or torsional damping.
    /// Frozen material coefficients and both positive heat accounts are retained.
    /// An omitted mechanism is not inferred from another viscosity or from Q.
    pub fn step_damped_candidate(
        &self,
        dt: f64,
        losses: RodAngularLosses<'_>,
        loads: &[RodNodalWrench],
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodMotionStep), RodMotionError> {
        self.step_with_angular_losses(dt, StageLoads::World(loads), losses, controls, gate)
    }

    /// The same material losses under held material-point forces and couples.
    pub fn step_spatial_damped_candidate(
        &self,
        dt: f64,
        losses: RodAngularLosses<'_>,
        loads: &PreparedRodLoads,
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodMotionStep), RodMotionError> {
        self.step_with_angular_losses(dt, StageLoads::Spatial(loads), losses, controls, gate)
    }
}
