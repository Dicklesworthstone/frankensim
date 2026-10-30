//! Two independently moving circular footprints on the original 3-D bridge.
//!
//! Unlike Young angles, hysteresis thresholds and line frictions retain their
//! original ordered material sources. The actual constant-mean-curvature
//! meridian is recovered from fixed volume, gap and two physical radii. The
//! original first-integral quadrature supplies both halves of that meridian;
//! no two-dimensional extrusion or replacement equilibrium source is used.
//!
//! In x=(h,b0,b1), the exact surface first variations are
//! E_h=pi*gamma*a*(2+q), E_bi=2*pi*gamma*bi*(cos(theta_i)-cos(theta_Yi)).
//! An average-vector-field discrete gradient integrates all three components
//! and their analytic Jacobian together along the straight physical path.
//! Each circular line has its own positive resistance 2*pi*zeta_i*mean(bi)/dt
//! and sourced advancing/receding pinning threshold. A coupled active-set
//! solve preserves cross-footprint compliance and separate physical losses.
//!
//! The branch is quasistatic, gravity-screened, acute-angle, connected and
//! line-friction dominated, with an interior neck and nonpositive pressure.
//! A conservative fixed-footprint axisymmetric stability screen and positive
//! local footprint Hessian restrict admission. No full 3-D stability,
//! bulk-flow inertia/loss, evaporation, or pinch-off model is claimed. Failed
//! shape continuation does not manufacture rupture or dissipated energy.
//! Quadrature and inverse-geometry uncertainties are estimates, not enclosures.
//!
//! The geometry and contact-line references are those of `contact_line`; the
//! joint physical-coordinate discrete gradient and perimeter integration are
//! derived here. They do not establish calibration for an arbitrary material.

pub use super::contact_line::AxisymmetricContactLineConfig as AsymmetricCircularContactLineConfig;
use super::{AxisymmetricError, MeridianPoint, SourceAxisymmetricCapillaryBridge};
use crate::capillary_bridge::contact_line::load_line::{
    DynamicContactLineLoadLine, DynamicContactLineLoadLineStep, DynamicContactLineModel,
    DynamicContactLineStepView, sealed,
};
use crate::capillary_bridge::contact_line::{
    ADVANCING_CONTACT_ANGLE_PROPERTY, CONTACT_LINE_FRICTION_PROPERTY,
    RECEDING_CONTACT_ANGLE_PROPERTY,
};
pub use crate::capillary_bridge::contact_line::{ContactLineError, ContactLineRegime};
use core::f64::consts::{FRAC_PI_2, PI};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_material::liquid::DYNAMIC_VISCOSITY_DIMS;
use fs_material::liquid::surface::require_frozen_coefficient_support;
use fs_math::det;
use fs_qty::{Dims, QuantitySpec};

mod evolution;
mod integration;
mod shape;

/// Original 3-D liquid geometry and work derivatives at physical coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsymmetricCircularContactLineProfile {
    pub gap_m: f64,
    pub wetted_radii_m: [f64; 2],
    pub contact_angles_radians: [f64; 2],
    pub suction_parameter: f64,
    pub neck_radius_m: f64,
    pub neck_height_m: f64,
    pub liquid_volume_m3: f64,
    pub liquid_gas_area_m2: f64,
    pub wetted_areas_m2: [f64; 2],
    pub liquid_minus_gas_pressure_pa: f64,
    pub surface_free_energy_j: f64,
    pub normal_force_n: f64,
    pub footprint_energy_gradients_n: [f64; 2],
    /// Complete Hessian in (gap, lower radius, upper radius), at fixed volume.
    pub energy_hessian_n_per_m: [[f64; 3]; 3],
    pub axisymmetric_stability_margin: f64,
    /// Canonical motion is separate from the numerical shape discrepancy.
    pub geometric_gap_residual_m: f64,
    pub geometric_radius_residuals_m: [f64; 2],
    pub quadrature_relative_error_estimate: f64,
    pub estimated_volume_error_m3: f64,
    pub estimated_surface_energy_error_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsymmetricCircularContactLineState {
    identity: ContentHash,
    evaluation: shape::Evaluation,
}
impl AsymmetricCircularContactLineState {
    pub const fn gap_m(self) -> f64 {
        self.evaluation.profile.gap_m
    }
    pub const fn wetted_radii_m(self) -> [f64; 2] {
        self.evaluation.profile.wetted_radii_m
    }
    pub const fn contact_angles_radians(self) -> [f64; 2] {
        self.evaluation.profile.contact_angles_radians
    }
    pub const fn profile(&self) -> &AsymmetricCircularContactLineProfile {
        &self.evaluation.profile
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AsymmetricCircularContactLineStep {
    state: AsymmetricCircularContactLineState,
    pub profile: AsymmetricCircularContactLineProfile,
    pub regimes: [ContactLineRegime; 2],
    pub mean_normal_force_n: f64,
    pub mean_force_gap_derivative_n_per_m: f64,
    pub surface_energy_change_j: f64,
    pub contact_line_dissipation_j: f64,
    pub viscous_line_dissipation_j: f64,
    pub pinning_dissipation_j: f64,
    pub per_line_viscous_dissipation_j: [f64; 2],
    pub per_line_pinning_dissipation_j: [f64; 2],
    pub work_on_liquid_j: f64,
    pub energy_balance_residual_j: f64,
    pub footprint_force_residuals_n: [f64; 2],
    pub estimated_force_error_n: f64,
    pub estimated_energy_error_j: f64,
    pub iterations: usize,
    /// Number of original half-meridian quadratures, including inverse solves.
    pub shape_evaluations: usize,
    pub quadrature_evaluations: u64,
}
impl AsymmetricCircularContactLineStep {
    pub const fn state(self) -> AsymmetricCircularContactLineState {
        self.state
    }
}

#[derive(Debug)]
pub struct SourceAsymmetricCircularContactLineBridge {
    original: SourceAxisymmetricCapillaryBridge,
    config: AsymmetricCircularContactLineConfig,
    identity: ContentHash,
    zeta: [f64; 2],
    advancing_per_radius_n_m: [f64; 2],
    receding_per_radius_n_m: [f64; 2],
}
impl SourceAsymmetricCircularContactLineBridge {
    pub fn new(
        original: SourceAxisymmetricCapillaryBridge,
        config: AsymmetricCircularContactLineConfig,
    ) -> Result<Self, ContactLineError> {
        config.validate()?;
        let mut values = [[0.0; 3]; 2];
        let mut advancing = [0.0; 2];
        let mut receding = [0.0; 2];
        for (i, wetting) in original.wetting().iter().enumerate() {
            for (j, (name, dims)) in [
                (CONTACT_LINE_FRICTION_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
                (ADVANCING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
                (RECEDING_CONTACT_ANGLE_PROPERTY, Dims::NONE),
            ]
            .into_iter()
            .enumerate()
            {
                let property =
                    wetting
                        .resolved()
                        .property(name)
                        .ok_or(ContactLineError::Invalid(
                            "missing original circular-line source",
                        ))?;
                if property.requirement().quantity() != QuantitySpec::dimensional(dims) {
                    return Err(ContactLineError::Invalid("circular-line source quantity"));
                }
                require_frozen_coefficient_support(property)?;
                values[i][j] = property.value_si();
            }
            let theta = wetting.contact_angle_radians();
            if values[i].iter().any(|v| !v.is_finite() || *v <= 0.0)
                || !(config.minimum_angle_radians <= values[i][2]
                    && values[i][2] <= theta
                    && theta <= values[i][1]
                    && values[i][1] <= config.maximum_angle_radians)
            {
                return Err(ContactLineError::Invalid(
                    "positive original line friction and ordered acute wetting angles",
                ));
            }
            let scale = 2.0 * PI * original.sigma;
            advancing[i] = if values[i][1] == theta {
                0.0
            } else {
                scale * (original.cosine[i] - det::cos(values[i][1]))
            };
            receding[i] = if values[i][2] == theta {
                0.0
            } else {
                scale * (det::cos(values[i][2]) - original.cosine[i])
            };
            if !advancing[i].is_finite()
                || advancing[i] < 0.0
                || !receding[i].is_finite()
                || receding[i] < 0.0
            {
                return Err(ContactLineError::Invalid(
                    "representable original circular pinning thresholds",
                ));
            }
        }
        let mut hash =
            DomainHasher::new("org.frankensim.fs-contact.asymmetric-circular-contact-line.v1");
        hash.update(original.identity().as_bytes());
        for v in config
            .numbers()
            .into_iter()
            .chain(values.into_iter().flatten())
        {
            hash.update(&v.to_bits().to_le_bytes());
        }
        for n in [
            config.maximum_iterations as u64,
            config.maximum_backtracks as u64,
            config.maximum_shape_evaluations as u64,
            u64::from(config.shape.maximum_root_iterations),
            u64::from(config.shape.maximum_quadrature_panels),
        ] {
            hash.update(&n.to_le_bytes());
        }
        Ok(Self {
            original,
            config,
            identity: hash.finalize(),
            zeta: [values[0][0], values[1][0]],
            advancing_per_radius_n_m: advancing,
            receding_per_radius_n_m: receding,
        })
    }
    pub const fn original(&self) -> &SourceAxisymmetricCapillaryBridge {
        &self.original
    }
    pub const fn config(&self) -> AsymmetricCircularContactLineConfig {
        self.config
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn line_frictions_pa_s(&self) -> [f64; 2] {
        self.zeta
    }

    /// Dynamic endpoint angles describe a physical initial shape, not sources.
    pub fn initialize_at_suction(
        &self,
        q: f64,
        angles_radians: [f64; 2],
        gate: &CancelGate,
    ) -> Result<AsymmetricCircularContactLineState, ContactLineError> {
        let evaluation = self.at_parameters(
            q,
            angles_radians.map(|a| FRAC_PI_2 - a),
            &mut shape::Work::default(),
            gate,
        )?;
        Ok(self.state_from(evaluation))
    }
    pub fn initialize(
        &self,
        gap_m: f64,
        angles_radians: [f64; 2],
        gate: &CancelGate,
    ) -> Result<AsymmetricCircularContactLineState, ContactLineError> {
        let evaluation =
            self.at_gap_angles(gap_m, angles_radians, &mut shape::Work::default(), gate)?;
        Ok(self.state_from(evaluation))
    }
    pub fn initialize_at_radii(
        &self,
        gap_m: f64,
        radii_m: [f64; 2],
        gate: &CancelGate,
    ) -> Result<AsymmetricCircularContactLineState, ContactLineError> {
        let mut work = shape::Work::default();
        let seed = self.at_gap_angles(
            gap_m,
            self.original
                .wetting()
                .each_ref()
                .map(|w| w.contact_angle_radians()),
            &mut work,
            gate,
        )?;
        let evaluation = self.at([gap_m, radii_m[0], radii_m[1]], seed, &mut work, gate)?;
        Ok(self.state_from(evaluation))
    }
    /// Authenticated constant-work read of the accepted shape.
    pub fn evaluate(
        &self,
        state: AsymmetricCircularContactLineState,
    ) -> Result<AsymmetricCircularContactLineProfile, ContactLineError> {
        self.require_state(state)?;
        Ok(state.evaluation.profile)
    }
    /// Fraction traverses the complete tangent-angle interval from lower plate
    /// to upper plate; the neck need not be at the midpoint.
    pub fn meridian_point(
        &self,
        state: AsymmetricCircularContactLineState,
        fraction: f64,
        gate: &CancelGate,
    ) -> Result<MeridianPoint, ContactLineError> {
        self.require_state(state)?;
        poll(gate)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(ContactLineError::Invalid("axisymmetric meridian fraction"));
        }
        let e = state.evaluation;
        let beta = -e.alpha[0] + fraction * (e.alpha[0] + e.alpha[1]);
        let (x, _, _) = super::numerics::radius(beta.abs(), e.q);
        let z = if beta == 0.0 {
            0.0
        } else {
            let integral = super::numerics::integrate(beta.abs(), e.q, self.config.shape, gate)
                .map_err(original_error)?;
            beta.signum() * e.profile.neck_radius_m * integral.values[0]
        };
        let (s, c) = det::sin_cos(beta);
        let radius = e.profile.neck_radius_m * x;
        let azimuthal = c / radius;
        Ok(MeridianPoint {
            height_m: e.profile.neck_height_m + z,
            radius_m: radius,
            radius_axial_derivative: s / c,
            azimuthal_curvature_per_m: azimuthal,
            meridional_curvature_per_m: -e.q / e.profile.neck_radius_m - azimuthal,
        })
    }
    fn state_from(&self, evaluation: shape::Evaluation) -> AsymmetricCircularContactLineState {
        AsymmetricCircularContactLineState {
            identity: self.identity,
            evaluation,
        }
    }
    fn require_state(
        &self,
        state: AsymmetricCircularContactLineState,
    ) -> Result<(), ContactLineError> {
        if state.identity == self.identity {
            Ok(())
        } else {
            Err(ContactLineError::StateMismatch)
        }
    }
}

pub type AsymmetricCircularContactLineLoadLine<'a> =
    DynamicContactLineLoadLine<'a, SourceAsymmetricCircularContactLineBridge>;
pub type AsymmetricCircularContactLineLoadLineStep =
    DynamicContactLineLoadLineStep<AsymmetricCircularContactLineStep>;
impl sealed::Sealed for SourceAsymmetricCircularContactLineBridge {}
impl DynamicContactLineModel for SourceAsymmetricCircularContactLineBridge {
    type State = AsymmetricCircularContactLineState;
    type Step = AsymmetricCircularContactLineStep;
    fn state_gap(state: Self::State) -> f64 {
        state.gap_m()
    }
    fn advance(
        &self,
        previous: Self::State,
        gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<Self::Step, ContactLineError> {
        SourceAsymmetricCircularContactLineBridge::advance(self, previous, gap_m, dt_s, gate)
    }
    fn step_view(step: &Self::Step) -> DynamicContactLineStepView {
        DynamicContactLineStepView {
            gap_m: step.state.gap_m(),
            mean_normal_force_n: step.mean_normal_force_n,
            mean_force_gap_derivative_n_per_m: step.mean_force_gap_derivative_n_per_m,
            estimated_force_error_n: step.estimated_force_error_n,
        }
    }
}
fn original_error(error: AxisymmetricError) -> ContactLineError {
    match error {
        AxisymmetricError::Physical(e) => ContactLineError::Original(e),
        AxisymmetricError::QuadratureBudget | AxisymmetricError::RootBudget => {
            ContactLineError::Budget
        }
        AxisymmetricError::StateMismatch => ContactLineError::StateMismatch,
        _ => ContactLineError::Domain("original axisymmetric shape branch"),
    }
}
fn poll(gate: &CancelGate) -> Result<(), ContactLineError> {
    if gate.is_requested() {
        Err(ContactLineError::Cancelled)
    } else {
        Ok(())
    }
}
