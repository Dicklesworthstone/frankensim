//! Intrinsic circular-shaft viscosity on the native finite-rotation trajectory.
//!
//! The original section law supplies c = sum(eta_s J)/L. For
//! phi = Log(R_i^-1 R_j), the twist-rate coordinate is d(phi_x)/dt, with
//! g_i = -J_left(phi)^-T e_x and g_j = J_right(phi)^-T e_x. These are the
//! SAME log-strain derivatives used by the native Cosserat energy; subtracting
//! axial spin components would be wrong under finite bend.
//!
//! A frozen-pose pair has rate v = sum(g dot omega), mobility
//! m = sum_free(g dot I_com^-1 g), and exact relaxation v(h)=exp(-c m h)v(0).
//! Applying its impulse through those same gradients gives nonnegative physical
//! heat and equal/opposite world moments. A palindromic pair sweep and symmetric
//! loss/native-motion/loss composition retain second-order temporal accuracy.
//! The elastic wave step still needs resolution; viscosity adds no explicit
//! stability restriction. Heat is reported, not deposited in an unnamed store.

use fs_exec::CancelGate;
use fs_ga::{So3, Vec3, so3_left_jacobian_inverse, so3_right_jacobian_inverse};
use fs_material::section::{ConcentricShearSection, torsion::ConcentricTorsionalSection};
use fs_solid::rod::MAX_MATERIAL_ROD_SEGMENTS;

use super::{
    KickAccounts, LieStepError, PreparedRodLoads, PreparedRodMotion, RodMotion, RodMotionControls,
    RodMotionError, RodMotionInertia, RodMotionStep, RodNodalWrench, StageLoads, a, add,
    checkpoint, finite, rotate,
};

/// Explicit fixed-state material shear loss, one source section per segment.
#[derive(Clone, Debug)]
pub struct PreparedRodTorsionalLoss {
    sections: Vec<ConcentricTorsionalSection>,
    resistance_nm_s: Vec<f64>,
    segment_length_m: f64,
}

/// Physical Kelvin shear heat in original segment order; never numerical damping.
#[derive(Clone, Debug, PartialEq)]
pub struct RodTorsionalHeat {
    /// Sum of intrinsic shear heat over this accepted step [J].
    pub total_heat_j: f64,
    /// Nonnegative physical heat in original source-segment order [J].
    pub segment_heat_j: Vec<f64>,
}

impl PreparedRodTorsionalLoss {
    /// Bind existing source-resolved torsional sections to the actual rod.
    /// Geometry, cards, selected mechanical properties and the complete query
    /// must match. Extra explicitly resolved viscosity is permitted, so a rod's
    /// mechanical source need not already have requested a loss property.
    ///
    /// Transients require constant coefficients at the held state and explicit
    /// source support over [0, f64::MAX] rad/s. A finite-band coupon cannot
    /// authorize arbitrary transients. No viscosity or damping ratio is inferred.
    pub fn new<I>(
        model: &PreparedRodMotion<I>,
        sections: Vec<ConcentricTorsionalSection>,
        gate: &CancelGate,
    ) -> Result<Self, RodMotionError> {
        checkpoint(gate)?;
        if sections.is_empty()
            || sections.len() > MAX_MATERIAL_ROD_SEGMENTS
            || sections.len() != model.source().sections().len()
        {
            return Err(RodMotionError::Invalid(
                "torsional loss needs matching bounded source sections",
            ));
        }
        let segment_length_m = model.source().reference_mass().segment_length_m();
        if !segment_length_m.is_finite() || segment_length_m <= 0.0 {
            return Err(RodMotionError::Invalid(
                "torsional loss requires positive finite reference length",
            ));
        }
        let mut resistance_nm_s = Vec::with_capacity(sections.len());
        for section in &sections {
            checkpoint(gate)?;
            if section.omega_band_rad_s() != (0.0, f64::MAX) {
                return Err(RodMotionError::Invalid(
                    "transient torsional viscosity requires explicit broadband source support",
                ));
            }
            let value = section.viscous_stiffness_n_m2_s() / segment_length_m;
            if !value.is_finite()
                || value < 0.0
                || (section.viscous_stiffness_n_m2_s() > 0.0 && value == 0.0)
            {
                return Err(RodMotionError::Invalid(
                    "unrepresentable torsional resistance",
                ));
            }
            resistance_nm_s.push(value);
        }
        let loss = Self {
            sections,
            resistance_nm_s,
            segment_length_m,
        };
        loss.admit(model, gate)?;
        Ok(loss)
    }

    /// Frozen annular section laws and their original source receipts.
    pub fn sections(&self) -> &[ConcentricTorsionalSection] {
        &self.sections
    }
    /// Per-segment resistance sum(eta_s J)/L, multiplying log-twist rate [N m s].
    pub fn resistance_nm_s(&self) -> &[f64] {
        &self.resistance_nm_s
    }

    pub(super) fn admit<I>(
        &self,
        model: &PreparedRodMotion<I>,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        let mechanical = model.source().sections();
        if self.sections.is_empty()
            || self.sections.len() > MAX_MATERIAL_ROD_SEGMENTS
            || self.sections.len() != mechanical.len()
            || self.segment_length_m != model.source().reference_mass().segment_length_m()
        {
            return Err(RodMotionError::Invalid(
                "torsional loss must match the reference rod segments",
            ));
        }
        for (source, current) in self.sections.iter().zip(mechanical) {
            checkpoint(gate)?;
            if !same_mechanics(current, source.shear()) {
                return Err(RodMotionError::Invalid(
                    "torsional loss source is stale or differs from the actual rod material state",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn heat(&self) -> RodTorsionalHeat {
        RodTorsionalHeat {
            total_heat_j: 0.0,
            segment_heat_j: vec![0.0; self.sections.len()],
        }
    }

    /// Exact individual dissipative flows; symmetric ordered composition when
    /// adjacent segments share inertia. Pose remains fixed throughout this map.
    pub(super) fn flow<I: RodMotionInertia>(
        &self,
        state: &mut RodMotion<I>,
        dt: f64,
        accounts: &mut KickAccounts,
        heat: &mut RodTorsionalHeat,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        // Bound scratch by the original rod ceiling; no matrix or per-edge
        // allocation, repeated source resolution or inertia factorization.
        let mut modes = [TwistMode::ZERO; MAX_MATERIAL_ROD_SEGMENTS];
        let n = self.sections.len();
        for (segment, mode) in modes[..n].iter_mut().enumerate() {
            checkpoint(gate)?;
            if self.resistance_nm_s[segment] == 0.0 {
                continue;
            }
            let directions =
                twist_gradients(state.pose.frames[segment], state.pose.frames[segment + 1])?;
            let mut mobility = 0.0;
            let mut free = false;
            for (side, gradient) in directions.iter().enumerate() {
                let node = segment + side;
                if !state.model.fixed_frames()[node] {
                    free = true;
                    mobility += dot(*gradient, state.kinetics[node].omega(*gradient)?);
                }
            }
            finite(&[mobility])?;
            if mobility < 0.0 || (free && mobility == 0.0) {
                return Err(RodMotionError::Invalid(
                    "unrepresentable positive torsional mobility",
                ));
            }
            *mode = TwistMode {
                directions,
                mobility,
            };
        }
        for segment in 0..n - 1 {
            self.relax(
                state,
                segment,
                modes[segment],
                0.5 * dt,
                accounts,
                heat,
                gate,
            )?;
        }
        self.relax(state, n - 1, modes[n - 1], dt, accounts, heat, gate)?;
        for segment in (0..n - 1).rev() {
            self.relax(
                state,
                segment,
                modes[segment],
                0.5 * dt,
                accounts,
                heat,
                gate,
            )?;
        }
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn relax<I: RodMotionInertia>(
        &self,
        state: &mut RodMotion<I>,
        segment: usize,
        mode: TwistMode,
        dt: f64,
        accounts: &mut KickAccounts,
        heat: &mut RodTorsionalHeat,
        gate: &CancelGate,
    ) -> Result<(), RodMotionError> {
        checkpoint(gate)?;
        let coefficient = self.resistance_nm_s[segment];
        if coefficient == 0.0 || mode.mobility == 0.0 {
            return Ok(());
        }
        let mut rate = 0.0;
        for side in 0..2 {
            let node = segment + side;
            rate += dot(
                mode.directions[side],
                state.kinetics[node].omega(state.spin[node])?,
            );
        }
        finite(&[rate])?;
        if rate == 0.0 {
            return Ok(());
        }
        let exposure = (coefficient * dt) * mode.mobility;
        if !exposure.is_finite() || exposure <= 0.0 {
            return Err(RodMotionError::Invalid(
                "unrepresentable torsional relaxation interval",
            ));
        }
        let fraction = -fs_math::det::expm1(-exposure);
        let impulse = rate * (fraction / mode.mobility);
        // Positive product form avoids subtracting nearly equal kinetic energies.
        let dissipated = (0.5 * impulse) * rate * (2.0 - fraction);
        finite(&[impulse, dissipated])?;
        if dissipated <= 0.0 || impulse == 0.0 {
            return Err(RodMotionError::Invalid(
                "unrepresentable positive torsional material heat",
            ));
        }
        for side in 0..2 {
            let node = segment + side;
            let gradient = mode.directions[side];
            if state.model.fixed_frames()[node] {
                let reaction = rotate(state.pose.frames[node], gradient.map(|g| impulse * g))?;
                accounts.angular = add(accounts.angular, reaction);
                accounts.supports[node].angular_nm_s =
                    add(accounts.supports[node].angular_nm_s, reaction);
                finite(&accounts.supports[node].angular_nm_s)?;
            } else {
                for c in 0..3 {
                    state.spin[node][c] = (-impulse).mul_add(gradient[c], state.spin[node][c]);
                }
                finite(&state.spin[node])?;
            }
        }
        heat.segment_heat_j[segment] += dissipated;
        heat.total_heat_j += dissipated;
        finite(&[heat.segment_heat_j[segment], heat.total_heat_j])?;
        Ok(())
    }
}

impl<I: RodMotionInertia> RodMotion<I> {
    /// Symmetric intrinsic-loss/native-motion/loss step under held world wrenches.
    /// Returned torsional heat belongs to the original segments, with no implied
    /// reservoir deposition or finite-temperature feedback. Changed source states
    /// must explicitly prepare a matching new loss law before continuing.
    pub fn step_torsional_candidate(
        &self,
        dt: f64,
        loss: &PreparedRodTorsionalLoss,
        loads: &[RodNodalWrench],
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodMotionStep), RodMotionError> {
        self.step_with_torsion_loads(dt, StageLoads::World(loads), Some(loss), controls, gate)
    }

    /// The same source-bound loss step with rotating material-point forces and
    /// body/world-held couples, evaluated at the native endpoint frames.
    pub fn step_spatial_torsional_candidate(
        &self,
        dt: f64,
        loss: &PreparedRodTorsionalLoss,
        loads: &PreparedRodLoads,
        controls: RodMotionControls,
        gate: &CancelGate,
    ) -> Result<(Self, RodMotionStep), RodMotionError> {
        self.step_with_torsion_loads(dt, StageLoads::Spatial(loads), Some(loss), controls, gate)
    }
}

#[derive(Clone, Copy)]
struct TwistMode {
    directions: [[f64; 3]; 2],
    mobility: f64,
}
impl TwistMode {
    const ZERO: Self = Self {
        directions: [[0.0; 3]; 2],
        mobility: 0.0,
    };
}

fn twist_gradients(left: So3, right: So3) -> Result<[[f64; 3]; 2], RodMotionError> {
    let relative = right.body_minus(left).map_err(LieStepError::from)?;
    if relative.angular.norm() >= core::f64::consts::PI * (1.0 - 64.0 * f64::EPSILON) {
        return Err(RodMotionError::Invalid(
            "torsional strain lies on the SO(3) logarithm cut",
        ));
    }
    let axis = Vec3::new(1.0, 0.0, 0.0);
    let left = so3_left_jacobian_inverse(relative)
        .map_err(LieStepError::from)?
        .transpose()
        .apply(axis);
    let right = so3_right_jacobian_inverse(relative)
        .map_err(LieStepError::from)?
        .transpose()
        .apply(axis);
    Ok([a(left.scale(-1.0)), a(right)])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|c| a[c] * b[c]).sum()
}

// The viscosity bundle may add the explicitly selected loss property, but may
// neither replace an original coefficient nor change its source/query meaning.
fn same_mechanics(original: &ConcentricShearSection, enriched: &ConcentricShearSection) -> bool {
    original.shear_property() == enriched.shear_property()
        && original.shear_correction() == enriched.shear_correction()
        && original.elastic().inner_radius_m() == enriched.elastic().inner_radius_m()
        && original.elastic().layers().len() == enriched.elastic().layers().len()
        && original
            .elastic()
            .layers()
            .iter()
            .zip(enriched.elastic().layers())
            .all(|(a, b)| {
                a.outer_radius_m == b.outer_radius_m
                    && a.material.card_identity() == b.material.card_identity()
                    && a.material.query_point() == b.material.query_point()
                    && a.material.axis_quantities() == b.material.axis_quantities()
                    && a.material
                        .properties()
                        .iter()
                        .all(|p| b.material.property(p.requirement().name()) == Some(p))
            })
}

#[cfg(test)]
mod geometry_tests {
    use super::*;
    use fs_ga::So3Tangent;

    #[test]
    fn finite_bend_log_gradient_matches_independent_rotation_probes_and_world_moments() {
        let left = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.1))).unwrap();
        let right = left
            .compose(So3::exp(So3Tangent::new(Vec3::new(0.3, 0.5, -0.4))).unwrap())
            .unwrap();
        let frames = [left, right];
        let gradients = twist_gradients(left, right).unwrap();
        let h = 2e-6;
        for side in 0..2 {
            for axis in 0..3 {
                let mut plus = frames;
                let mut minus = frames;
                let mut direction = [0.0; 3];
                direction[axis] = h;
                let tangent = So3Tangent::new(super::super::v(direction));
                plus[side] = frames[side].body_plus(tangent).unwrap();
                minus[side] = frames[side].body_plus(tangent.scale(-1.0)).unwrap();
                let derivative = (plus[1].body_minus(plus[0]).unwrap().angular.x
                    - minus[1].body_minus(minus[0]).unwrap().angular.x)
                    / (2.0 * h);
                assert!((derivative - gradients[side][axis]).abs() < 2e-9);
            }
        }
        let world_sum = add(
            rotate(left, gradients[0]).unwrap(),
            rotate(right, gradients[1]).unwrap(),
        );
        assert!(world_sum.iter().all(|v| v.abs() < 2e-14));
        let world_spin = [0.3, -0.7, 0.2];
        let rate = dot(gradients[0], rotate(left.inverse(), world_spin).unwrap())
            + dot(gradients[1], rotate(right.inverse(), world_spin).unwrap());
        assert!(
            rate.abs() < 2e-14,
            "a common rigid spin must not become material loss"
        );
    }
}
