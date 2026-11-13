//! Source-bound transverse shear dissipation on the native Cosserat trajectory.
//!
//! `e = L0 [gamma_y, gamma_z]` uses the original geodesic-midpoint frame.
//! Each direction has resistance `C_s/L0`, and its full translation/rotation
//! pullback supplies both velocity and reaction. Frozen-pose rank-one flows
//! use the SAME exact exponential impulse as axial/torsional viscosity.
//! Palindromic direction, segment and mechanism sweeps preserve second order.
//! Rigid translations and rotations do not become material heat. Shear heat is
//! separate from axial, bending, torsional and numerical loss; no unbound
//! thermal reservoir, hereditary spectrum or source update is inferred.

use fs_exec::CancelGate;
use fs_material::section::transverse::ConcentricTransverseViscosity;
use fs_solid::rod::MAX_MATERIAL_ROD_SEGMENTS;

use super::axial_loss::{DisplacementMode, displacement_mode, relax_displacement};
use super::torsional_loss::same_mechanics;
use super::{
    KickAccounts, PreparedRodAxialLoss, PreparedRodLoads, PreparedRodMotion, RodAngularLosses,
    RodMotion, RodMotionControls, RodMotionError, RodMotionInertia, RodMotionStep, RodNodalWrench,
    StageLoads, checkpoint, finite,
};

/// Explicit shear-viscosity sections, one per original reference segment.
#[derive(Clone, Debug)]
pub struct PreparedRodTransverseLoss {
    sections: Vec<ConcentricTransverseViscosity>,
    resistance_n_s_m: Vec<f64>,
    segment_lengths_m: Vec<f64>,
}

/// Physical transverse shear heat, summed over both material directions.
#[derive(Clone, Debug, PartialEq)]
pub struct RodTransverseHeat {
    pub total_heat_j: f64,
    /// Original reference-segment order [J].
    pub segment_heat_j: Vec<f64>,
}

/// Independently selected material loss mechanisms on one native rod history.
/// An omitted mechanism stays unmodeled; it is not inferred from another law.
#[derive(Clone, Copy, Debug)]
pub struct RodMaterialLosses<'a> {
    pub angular: RodAngularLosses<'a>,
    pub axial: Option<&'a PreparedRodAxialLoss>,
    pub transverse: Option<&'a PreparedRodTransverseLoss>,
}

impl PreparedRodTransverseLoss {
    /// Retain the SAME geometry, shear correction, material identities, query
    /// and elastic properties as the rod. Additional explicit viscosity is
    /// allowed. A finite-band coupon cannot authorize this transient Kelvin
    /// law: every source section must support [0, f64::MAX] rad/s.
    pub fn new<I>(
        model: &PreparedRodMotion<I>,
        sections: Vec<ConcentricTransverseViscosity>,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        checkpoint(gate)?;
        if sections.is_empty()
            || sections.len() > MAX_MATERIAL_ROD_SEGMENTS
            || sections.len() != model.source().sections().len()
        {
            return Err(RodMotionError::Invalid(
                "transverse loss needs matching bounded source sections",
            ));
        }
        let segment_lengths_m = model.source().reference_mass().segment_lengths_m().to_vec();
        let mut resistance_n_s_m = Vec::with_capacity(sections.len());
        for (section, &segment_length_m) in sections.iter().zip(&segment_lengths_m) {
            checkpoint(gate)?;
            if section.omega_band_rad_s() != (0.0, f64::MAX) {
                return Err(RodMotionError::Invalid(
                    "transient transverse viscosity requires explicit broadband source support",
                ));
            }
            let value = section.viscous_rigidity_n_s() / segment_length_m;
            if !value.is_finite()
                || value < 0.0
                || (section.viscous_rigidity_n_s() > 0.0 && value == 0.0)
            {
                return Err(RodMotionError::Invalid(
                    "unrepresentable transverse shear resistance",
                ));
            }
            resistance_n_s_m.push(value);
        }
        let loss = Self {
            sections,
            resistance_n_s_m,
            segment_lengths_m,
        };
        loss.admit(model, gate)?;
        Ok(loss)
    }

    pub fn sections(&self) -> &[ConcentricTransverseViscosity] {
        &self.sections
    }
    /// C_s/L0 multiplying each physical transverse displacement rate [N s/m].
    pub fn resistance_n_s_m(&self) -> &[f64] {
        &self.resistance_n_s_m
    }

    pub(super) fn admit<I>(
        &self,
        model: &PreparedRodMotion<I>,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        if self.sections.len() != model.source().sections().len()
            || self.segment_lengths_m != model.source().reference_mass().segment_lengths_m()
        {
            return Err(RodMotionError::Invalid(
                "transverse loss must match the reference rod segments",
            ));
        }
        for (source, current) in self.sections.iter().zip(model.source().sections()) {
            checkpoint(gate)?;
            if !same_mechanics(current, source.shear()) {
                return Err(RodMotionError::Invalid(
                    "transverse loss source is stale or differs from the actual rod material state",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn heat(&self) -> RodTransverseHeat {
        RodTransverseHeat {
            total_heat_j: 0.0,
            segment_heat_j: vec![0.0; self.sections.len()],
        }
    }
    pub(super) fn active(&self) -> bool {
        self.resistance_n_s_m.iter().any(|&c| c != 0.0)
    }

    pub(super) fn flow<I: RodMotionInertia>(
        &self,
        state: &mut RodMotion<I>,
        dt: f64,
        origin: [f64; 3],
        accounts: &mut KickAccounts,
        heat: &mut RodTransverseHeat,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        let n = self.sections.len();
        let mut gradients = [[[[0.0; 6]; 2]; 2]; MAX_MATERIAL_ROD_SEGMENTS];
        state.model.transverse_displacement_gradients_into(
            &state.pose,
            &mut gradients[..n],
            gate,
        )?;
        let mut modes = [DisplacementMode::ZERO; 2 * MAX_MATERIAL_ROD_SEGMENTS];
        for segment in 0..n {
            checkpoint(gate)?;
            if self.resistance_n_s_m[segment] != 0.0 {
                for axis in 0..2 {
                    modes[2 * segment + axis] =
                        displacement_mode(state, segment, gradients[segment][axis])?;
                }
            }
        }
        let mut relax = |index: usize, interval: f64| -> Result<(), RodMotionError> {
            let segment = index / 2;
            let heat_j = relax_displacement(
                state,
                segment,
                modes[index],
                self.resistance_n_s_m[segment],
                interval,
                origin,
                accounts,
                gate,
            )?;
            heat.segment_heat_j[segment] += heat_j;
            heat.total_heat_j += heat_j;
            finite(&[heat.segment_heat_j[segment], heat.total_heat_j])
        };
        for index in 0..2 * n - 1 {
            relax(index, 0.5 * dt)?;
        }
        relax(2 * n - 1, dt)?;
        for index in (0..2 * n - 1).rev() {
            relax(index, 0.5 * dt)?;
        }
        Ok(())
    }
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Native motion with any explicit axial, transverse, bending and torsional
    /// sources. Physical heat, support impulses and the integration defect are
    /// checked before publishing the candidate; refusal preserves this state.
    pub fn step_material_damped_candidate(
        &self,
        dt: f64,
        losses: RodMaterialLosses<'_>,
        loads: &[RodNodalWrench],
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodMotionStep), RodMotionError> {
        self.step_with_all_material_losses(dt, StageLoads::World(loads), losses, controls, gate)
    }

    /// The same material-loss composition with prepared material-point loads.
    pub fn step_spatial_material_damped_candidate(
        &self,
        dt: f64,
        losses: RodMaterialLosses<'_>,
        loads: &PreparedRodLoads,
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodMotionStep), RodMotionError> {
        self.step_with_all_material_losses(dt, StageLoads::Spatial(loads), losses, controls, gate)
    }
}
