//! Intrinsic axial viscosity on the original finite-rotation rod trajectory.
//!
//! The existing section supplies C_a = sum(eta_i A_i). For extension
//! e = L0*gamma_x the resistance is C_a/L0. Its velocity and reaction use the
//! SAME geodesic-midpoint strain derivative as the native elastic rod. Axial
//! viscosity therefore includes rotating-frame terms in a bent configuration;
//! a common rigid translation or spin cannot become material damping.
//!
//! Frozen-pose rank-one flows use the original exact viscous impulse. Symmetric
//! segment and mechanism sweeps retain second-order splitting. Reference mass,
//! full central rotary tensors, eccentric centroids and support reactions are
//! preserved. Reported physical heat is not numerical loss or an implicit
//! thermal deposit. Small-strain Kelvin coefficients remain frozen at their
//! admitted source state; transverse shear and hereditary loss are not inferred.

use fs_exec::CancelGate;
use fs_material::section::axial::ConcentricAxialViscosity;
use fs_solid::rod::MAX_MATERIAL_ROD_SEGMENTS;

use super::torsional_loss::{rank_one_relaxation, same_elastic_source};
use super::{
    KickAccounts, PreparedRodLoads, PreparedRodMotion, RodAngularLosses, RodMotion,
    RodMotionControls, RodMotionError, RodMotionInertia, RodMotionStep, RodNodalWrench, StageLoads,
    add, checkpoint, cross, finite, rotate, sub,
};

/// Original axial-viscosity section for each reference rod segment.
#[derive(Clone, Debug)]
pub struct PreparedRodAxialLoss {
    sections: Vec<ConcentricAxialViscosity>,
    resistance_n_s_m: Vec<f64>,
    segment_length_m: f64,
}

/// Physical axial heat in original segment order, separate from angular loss.
#[derive(Clone, Debug, PartialEq)]
pub struct RodAxialHeat {
    pub total_heat_j: f64,
    pub segment_heat_j: Vec<f64>,
}

impl PreparedRodAxialLoss {
    /// Bind explicit broadband axial viscosity to the SAME elastic sources,
    /// held material states and reference geometry as the dynamic rod.
    /// Missing viscosity is not zero, and a finite-frequency sample cannot
    /// authorize a transient Kelvin law. Additional loss properties are allowed
    /// in each source bundle; replacing an original mechanical property is not.
    pub fn new<I>(
        model: &PreparedRodMotion<I>,
        sections: Vec<ConcentricAxialViscosity>,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        checkpoint(gate)?;
        if sections.is_empty()
            || sections.len() > MAX_MATERIAL_ROD_SEGMENTS
            || sections.len() != model.source().sections().len()
        {
            return Err(RodMotionError::Invalid(
                "axial loss needs matching bounded source sections",
            ));
        }
        let segment_length_m = model.source().reference_mass().segment_length_m();
        let mut resistance_n_s_m = Vec::with_capacity(sections.len());
        for section in &sections {
            checkpoint(gate)?;
            if section.omega_band_rad_s() != (0.0, f64::MAX) {
                return Err(RodMotionError::Invalid(
                    "transient axial viscosity requires explicit broadband source support",
                ));
            }
            let value = section.viscous_rigidity_n_s() / segment_length_m;
            if !value.is_finite()
                || value < 0.0
                || (section.viscous_rigidity_n_s() > 0.0 && value == 0.0)
            {
                return Err(RodMotionError::Invalid("unrepresentable axial resistance"));
            }
            resistance_n_s_m.push(value);
        }
        let loss = Self {
            sections,
            resistance_n_s_m,
            segment_length_m,
        };
        loss.admit(model, gate)?;
        Ok(loss)
    }

    /// Original ordered layers, viscosity selectors and material receipts.
    pub fn sections(&self) -> &[ConcentricAxialViscosity] {
        &self.sections
    }
    /// C_a/L0 multiplying the actual axial extension rate [N s/m].
    pub fn resistance_n_s_m(&self) -> &[f64] {
        &self.resistance_n_s_m
    }

    pub(super) fn admit<I>(
        &self,
        model: &PreparedRodMotion<I>,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        if self.sections.len() != model.source().sections().len()
            || self.segment_length_m != model.source().reference_mass().segment_length_m()
        {
            return Err(RodMotionError::Invalid(
                "axial loss must match the reference rod segments",
            ));
        }
        for (source, current) in self.sections.iter().zip(model.source().sections()) {
            checkpoint(gate)?;
            if !same_elastic_source(current.elastic(), source.elastic()) {
                return Err(RodMotionError::Invalid(
                    "axial loss source is stale or differs from the actual rod material state",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn heat(&self) -> RodAxialHeat {
        RodAxialHeat {
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
        heat: &mut RodAxialHeat,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        let n = self.sections.len();
        let mut gradients = [[[0.0; 6]; 2]; MAX_MATERIAL_ROD_SEGMENTS];
        state
            .model
            .axial_extension_gradients_into(&state.pose, &mut gradients[..n], gate)?;
        // Convert material-node covectors into the SAME COM/spin coordinates
        // used by native motion. This preserves all translation/rotation terms.
        let mut modes = [AxialMode::ZERO; MAX_MATERIAL_ROD_SEGMENTS];
        for (segment, mode) in modes[..n].iter_mut().enumerate() {
            checkpoint(gate)?;
            if self.resistance_n_s_m[segment] == 0.0 {
                continue;
            }
            let mut free = false;
            for side in 0..2 {
                let node = segment + side;
                let k = state.kinetics[node];
                let g = gradients[segment][side];
                let force = [g[0], g[1], g[2]];
                let torque = [g[3], g[4], g[5]];
                mode.force[side] = force;
                mode.spin[side] = if k.eccentric() {
                    sub(
                        torque,
                        cross(k.center, rotate(state.pose.frames[node].inverse(), force)?),
                    )
                } else {
                    torque
                };
                for (c, &g) in force.iter().enumerate() {
                    if !state.model.fixed_positions()[node][c] {
                        free |= g != 0.0;
                        mode.mobility += g * (g / k.mass);
                    }
                }
                if !state.model.fixed_frames()[node] {
                    free |= mode.spin[side].iter().any(|&g| g != 0.0);
                    mode.mobility += dot(mode.spin[side], k.omega(mode.spin[side])?);
                }
            }
            finite(&[mode.mobility])?;
            if mode.mobility < 0.0 || (free && mode.mobility == 0.0) {
                return Err(RodMotionError::Invalid(
                    "unrepresentable positive axial-strain mobility",
                ));
            }
        }
        let mut relax = |segment: usize, h: f64| -> Result<(), RodMotionError> {
            checkpoint(gate)?;
            let mode = modes[segment];
            let mut rate = 0.0;
            for side in 0..2 {
                let node = segment + side;
                let k = state.kinetics[node];
                rate += dot(mode.force[side], state.linear[node].map(|p| p / k.mass))
                    + dot(mode.spin[side], k.omega(state.spin[node])?);
            }
            finite(&[rate])?;
            let (impulse, dissipated) =
                rank_one_relaxation(rate, mode.mobility, self.resistance_n_s_m[segment], h)?;
            if impulse == 0.0 {
                return Ok(());
            }
            for side in 0..2 {
                let node = segment + side;
                let frame = state.pose.frames[node];
                let k = state.kinetics[node];
                let mut support_force = [0.0; 3];
                let mut support_torque = [0.0; 3];
                for c in 0..3 {
                    if state.model.fixed_positions()[node][c] {
                        support_force[c] = impulse * mode.force[side][c];
                    } else {
                        state.linear[node][c] =
                            (-impulse).mul_add(mode.force[side][c], state.linear[node][c]);
                    }
                    if state.model.fixed_frames()[node] {
                        support_torque[c] = impulse * mode.spin[side][c];
                    } else {
                        state.spin[node][c] =
                            (-impulse).mul_add(mode.spin[side][c], state.spin[node][c]);
                    }
                }
                if k.eccentric() && state.model.fixed_frames()[node] {
                    support_torque = add(
                        support_torque,
                        cross(k.center, rotate(frame.inverse(), support_force)?),
                    );
                }
                let support_torque = rotate(frame, support_torque)?;
                accounts.linear = add(accounts.linear, support_force);
                accounts.angular = add(
                    accounts.angular,
                    add(
                        support_torque,
                        cross(sub(state.pose.positions_m[node], origin), support_force),
                    ),
                );
                accounts.supports[node].linear_n_s =
                    add(accounts.supports[node].linear_n_s, support_force);
                accounts.supports[node].angular_nm_s =
                    add(accounts.supports[node].angular_nm_s, support_torque);
                finite(&state.linear[node])?;
                finite(&state.spin[node])?;
                finite(&accounts.supports[node].linear_n_s)?;
                finite(&accounts.supports[node].angular_nm_s)?;
            }
            heat.segment_heat_j[segment] += dissipated;
            heat.total_heat_j += dissipated;
            finite(&[heat.segment_heat_j[segment], heat.total_heat_j])
        };
        for i in 0..n - 1 {
            relax(i, 0.5 * dt)?;
        }
        relax(n - 1, dt)?;
        for i in (0..n - 1).rev() {
            relax(i, 0.5 * dt)?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct AxialMode {
    force: [[f64; 3]; 2],
    spin: [[f64; 3]; 2],
    mobility: f64,
}
impl AxialMode {
    const ZERO: Self = Self {
        force: [[0.0; 3]; 2],
        spin: [[0.0; 3]; 2],
        mobility: 0.0,
    };
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Native motion with explicit axial viscosity and independently selected
    /// angular losses. Omitted angular mechanisms are not inferred from eta_a.
    /// Frozen-pose loss substeps surround the same elastic motion. The public
    /// state is published only after source, cancellation and energy gates pass.
    pub fn step_axially_damped_candidate(
        &self,
        dt: f64,
        axial: &PreparedRodAxialLoss,
        angular: RodAngularLosses<'_>,
        loads: &[RodNodalWrench],
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodMotionStep), RodMotionError> {
        self.step_with_material_losses(
            dt,
            StageLoads::World(loads),
            angular,
            Some(axial),
            controls,
            gate,
        )
    }

    /// The same axial/angular loss composition under prepared material-point
    /// forces and couples, preserving their original application geometry.
    pub fn step_spatial_axially_damped_candidate(
        &self,
        dt: f64,
        axial: &PreparedRodAxialLoss,
        angular: RodAngularLosses<'_>,
        loads: &PreparedRodLoads,
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodMotionStep), RodMotionError> {
        self.step_with_material_losses(
            dt,
            StageLoads::Spatial(loads),
            angular,
            Some(axial),
            controls,
            gate,
        )
    }
}
