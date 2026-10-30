//! Finite-rate wetting on a genuine three-dimensional Young--Laplace bridge.
//!
//! This first dynamic branch is mirror-symmetric: the two original wetting
//! bundles must have identical Young angles, line friction and hysteresis
//! thresholds, and both circular footprints have the same evolving radius b.
//! Their original identities remain distinct. The complete original liquid
//! owner is retained; a different equilibrium angle is never fabricated as a
//! material source. Dynamic angle, neck and pressure follow from fixed volume,
//! actual gap and the independently evolved physical contact radius.
//!
//! The existing axisymmetric first integral and angular quadrature supply the
//! actual constant-mean-curvature meridian. Its first variations are
//! E_h=pi*gamma*a*(2+q) and E_b=4*pi*gamma*b*(cos(theta)-cos(theta_Y)).
//! Two circular lines dissipate 4*pi*zeta*b*b_dot^2. A symmetric discrete
//! gradient in (h,b), with the exact mean perimeter of the linear radius path,
//! accounts for normal work, signed surface storage, line friction and sourced
//! advancing/receding pinning separately. Analytic shape and implicit-step
//! derivatives supply the mechanical tangent. No runtime allocation is needed.
//!
//! The shape is quasistatic and line-friction dominated. It has an interior
//! neck, acute angles and nonpositive pressure on the original declared q
//! interval. A conservative fixed-footprint axisymmetric stability screen and
//! positive local footprint Hessian restrict admission; neither certifies full
//! 3-D or dynamic stability. Bulk-flow loss/inertia, asymmetric circular
//! footprints, evaporation and topology change are outside this owner. A shape
//! refusal does not create rupture, release, or heat. Numerical uncertainty is
//! estimated from quadrature, inverse geometry and arithmetic, not enclosed.
//!
//! Geometry: Fortes, https://doi.org/10.1016/0021-9797(82)90263-6.
//! Low-speed line law: https://doi.org/10.1038/s41526-021-00186-0, equation (2).
//! Pinning closure: https://doi.org/10.1021/acs.langmuir.3c01550.
//! The finite-step geometric resistance and work integration are derived here;
//! these references do not validate arbitrary materials or transient bridges.

use super::{
    AxisymmetricError, AxisymmetricNumerics, MeridianPoint, SourceAxisymmetricCapillaryBridge,
};
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

/// Explicit shape, motion and work bounds for the symmetric circular branch.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricContactLineConfig {
    pub shape: AxisymmetricNumerics,
    pub minimum_angle_radians: f64,
    pub maximum_angle_radians: f64,
    pub minimum_neck_radius_m: f64,
    pub maximum_gap_increment_m: f64,
    pub maximum_radius_increment_m: f64,
    pub maximum_step_s: f64,
    pub maximum_contact_line_speed_m_s: f64,
    pub maximum_angle_departure_radians: f64,
    /// Neglected bulk-flow applicability screens, not numerical error bounds.
    pub maximum_capillary_number: f64,
    pub maximum_weber_number: f64,
    pub maximum_iterations: usize,
    pub maximum_backtracks: usize,
    /// ALL original first-integral quadratures within one public initialization
    /// or advance, including inverse solves, backtracks and work nodes. Each
    /// quadrature additionally retains shape.maximum_quadrature_panels. Native
    /// mechanical continuation separately bounds its number of advance calls.
    pub maximum_shape_evaluations: usize,
    pub geometry_absolute_tolerance_m: f64,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
    pub energy_absolute_tolerance_j: f64,
    pub energy_relative_tolerance: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricContactLineProfile {
    /// Canonical physical coordinate. The converged shape residual is separate.
    pub gap_m: f64,
    pub wetted_radius_m: f64,
    pub contact_angle_radians: f64,
    pub suction_parameter: f64,
    pub neck_radius_m: f64,
    pub neck_height_m: f64,
    pub liquid_volume_m3: f64,
    pub liquid_gas_area_m2: f64,
    pub wetted_areas_m2: [f64; 2],
    pub liquid_minus_gas_pressure_pa: f64,
    pub surface_free_energy_j: f64,
    pub normal_force_n: f64,
    pub footprint_energy_gradient_n: f64,
    /// Derivatives at fixed physical footprint or fixed gap, respectively.
    pub normal_force_gap_derivative_n_per_m: f64,
    pub normal_force_radius_derivative_n_per_m: f64,
    pub footprint_gradient_radius_derivative_n_per_m: f64,
    pub axisymmetric_stability_margin: f64,
    pub geometric_gap_residual_m: f64,
    pub geometric_radius_residual_m: f64,
    pub quadrature_relative_error_estimate: f64,
    pub estimated_volume_error_m3: f64,
    pub estimated_surface_energy_error_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricContactLineState {
    identity: ContentHash,
    evaluation: shape::Evaluation,
}
impl AxisymmetricContactLineState {
    pub const fn gap_m(self) -> f64 {
        self.evaluation.profile.gap_m
    }
    pub const fn wetted_radius_m(self) -> f64 {
        self.evaluation.profile.wetted_radius_m
    }
    pub const fn contact_angle_radians(self) -> f64 {
        self.evaluation.profile.contact_angle_radians
    }
    pub const fn profile(&self) -> &AxisymmetricContactLineProfile {
        &self.evaluation.profile
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisymmetricContactLineStep {
    state: AxisymmetricContactLineState,
    pub profile: AxisymmetricContactLineProfile,
    pub regime: ContactLineRegime,
    pub mean_normal_force_n: f64,
    pub mean_force_gap_derivative_n_per_m: f64,
    pub surface_energy_change_j: f64,
    pub contact_line_dissipation_j: f64,
    pub viscous_line_dissipation_j: f64,
    pub pinning_dissipation_j: f64,
    pub work_on_liquid_j: f64,
    pub energy_balance_residual_j: f64,
    pub footprint_force_residual_n: f64,
    pub estimated_force_error_n: f64,
    pub estimated_energy_error_j: f64,
    pub iterations: usize,
    pub shape_evaluations: usize,
    pub quadrature_evaluations: u64,
}
impl AxisymmetricContactLineStep {
    pub const fn state(self) -> AxisymmetricContactLineState {
        self.state
    }
}

#[derive(Debug)]
pub struct SourceAxisymmetricContactLineBridge {
    original: SourceAxisymmetricCapillaryBridge,
    config: AxisymmetricContactLineConfig,
    identity: ContentHash,
    zeta: f64,
    /// Pinning resistance divided by common radius, for BOTH circular lines.
    advancing_per_radius_n_m: f64,
    receding_per_radius_n_m: f64,
}

impl SourceAxisymmetricContactLineBridge {
    pub fn new(
        original: SourceAxisymmetricCapillaryBridge,
        config: AxisymmetricContactLineConfig,
    ) -> Result<Self, ContactLineError> {
        config.validate()?;
        let mut values = [[0.0; 3]; 2];
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
        }
        let theta = original.wetting()[0].contact_angle_radians();
        if original.cosine[0] != original.cosine[1]
            || values[0] != values[1]
            || values[0].iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !(config.minimum_angle_radians <= values[0][2]
                && values[0][2] <= theta
                && theta <= values[0][1]
                && values[0][1] <= config.maximum_angle_radians)
        {
            return Err(ContactLineError::Invalid(
                "symmetric circular footprints require equal Young angle/friction/thresholds and ordered acute angles",
            ));
        }
        let scale = 4.0 * PI * original.sigma;
        let advancing = if values[0][1] == theta {
            0.0
        } else {
            scale * (original.cosine[0] - det::cos(values[0][1]))
        };
        let receding = if values[0][2] == theta {
            0.0
        } else {
            scale * (det::cos(values[0][2]) - original.cosine[0])
        };
        if advancing < 0.0 || receding < 0.0 || !advancing.is_finite() || !receding.is_finite() {
            return Err(ContactLineError::Invalid(
                "representable circular-line pinning",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.axisymmetric-contact-line.v1");
        hash.update(original.identity().as_bytes());
        for v in config.numbers().into_iter().chain(values[0]) {
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
            zeta: values[0][0],
            advancing_per_radius_n_m: advancing,
            receding_per_radius_n_m: receding,
        })
    }

    pub const fn original(&self) -> &SourceAxisymmetricCapillaryBridge {
        &self.original
    }
    pub const fn config(&self) -> AxisymmetricContactLineConfig {
        self.config
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn line_friction_pa_s(&self) -> f64 {
        self.zeta
    }

    /// An initial dynamic angle is a physical shape condition, not a new source.
    pub fn initialize_at_suction(
        &self,
        q: f64,
        angle_radians: f64,
        gate: &CancelGate,
    ) -> Result<AxisymmetricContactLineState, ContactLineError> {
        let mut work = shape::Work::default();
        let evaluation = self.at_parameters(q, FRAC_PI_2 - angle_radians, &mut work, gate)?;
        Ok(self.state_from(evaluation))
    }

    pub fn initialize(
        &self,
        gap_m: f64,
        angle_radians: f64,
        gate: &CancelGate,
    ) -> Result<AxisymmetricContactLineState, ContactLineError> {
        let mut work = shape::Work::default();
        let evaluation = self.at_gap_angle(gap_m, angle_radians, &mut work, gate)?;
        Ok(self.state_from(evaluation))
    }

    /// Independently initialize a physical footprint; inverse geometry may
    /// refuse when the stated local connected branch cannot support it.
    pub fn initialize_at_radius(
        &self,
        gap_m: f64,
        radius_m: f64,
        gate: &CancelGate,
    ) -> Result<AxisymmetricContactLineState, ContactLineError> {
        let mut work = shape::Work::default();
        let seed = self.at_gap_angle(
            gap_m,
            self.original.wetting()[0].contact_angle_radians(),
            &mut work,
            gate,
        )?;
        let evaluation = self.at(gap_m, radius_m, seed, &mut work, gate)?;
        Ok(self.state_from(evaluation))
    }

    /// Authenticated O(1) read of accepted shape and surface storage.
    pub fn evaluate(
        &self,
        state: AxisymmetricContactLineState,
    ) -> Result<AxisymmetricContactLineProfile, ContactLineError> {
        self.require_state(state)?;
        Ok(state.evaluation.profile)
    }

    /// Sample the actual Young--Laplace meridian, not a fitted circle/sphere.
    /// The canonical-gap residual in the profile reports the shape's geometric
    /// discrepancy. Fraction parameterizes tangent angle, not axial height.
    pub fn meridian_point(
        &self,
        state: AxisymmetricContactLineState,
        fraction: f64,
        gate: &CancelGate,
    ) -> Result<MeridianPoint, ContactLineError> {
        self.require_state(state)?;
        poll(gate)?;
        if !fraction.is_finite() || !(0.0..=1.0).contains(&fraction) {
            return Err(ContactLineError::Invalid("axisymmetric meridian fraction"));
        }
        let e = state.evaluation;
        let beta = (2.0 * fraction - 1.0) * e.alpha;
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

    fn state_from(&self, evaluation: shape::Evaluation) -> AxisymmetricContactLineState {
        AxisymmetricContactLineState {
            identity: self.identity,
            evaluation,
        }
    }
    fn require_state(&self, state: AxisymmetricContactLineState) -> Result<(), ContactLineError> {
        if state.identity == self.identity {
            Ok(())
        } else {
            Err(ContactLineError::StateMismatch)
        }
    }
}

pub type AxisymmetricContactLineLoadLine<'a> =
    DynamicContactLineLoadLine<'a, SourceAxisymmetricContactLineBridge>;
pub type AxisymmetricContactLineLoadLineStep =
    DynamicContactLineLoadLineStep<AxisymmetricContactLineStep>;
impl sealed::Sealed for SourceAxisymmetricContactLineBridge {}
impl DynamicContactLineModel for SourceAxisymmetricContactLineBridge {
    type State = AxisymmetricContactLineState;
    type Step = AxisymmetricContactLineStep;
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
        SourceAxisymmetricContactLineBridge::advance(self, previous, gap_m, dt_s, gate)
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
        AxisymmetricError::Physical(error) => ContactLineError::Original(error),
        AxisymmetricError::QuadratureBudget | AxisymmetricError::RootBudget => {
            ContactLineError::Budget
        }
        AxisymmetricError::StateMismatch => ContactLineError::StateMismatch,
        AxisymmetricError::StabilityUnestablished { .. } => {
            ContactLineError::Domain("axisymmetric fixed-footprint stability screen")
        }
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
impl AxisymmetricContactLineConfig {
    pub(super) fn numbers(self) -> [f64; 17] {
        [
            self.minimum_angle_radians,
            self.maximum_angle_radians,
            self.minimum_neck_radius_m,
            self.maximum_gap_increment_m,
            self.maximum_radius_increment_m,
            self.maximum_step_s,
            self.maximum_contact_line_speed_m_s,
            self.maximum_angle_departure_radians,
            self.maximum_capillary_number,
            self.maximum_weber_number,
            self.geometry_absolute_tolerance_m,
            self.force_absolute_tolerance_n,
            self.force_relative_tolerance,
            self.energy_absolute_tolerance_j,
            self.energy_relative_tolerance,
            self.shape.quadrature_relative_tolerance,
            self.shape.gap_relative_tolerance,
        ]
    }
    pub(super) fn validate(self) -> Result<(), ContactLineError> {
        self.shape.validate().map_err(original_error)?;
        if self.numbers().iter().any(|v| !v.is_finite() || *v <= 0.0)
            || self.minimum_angle_radians >= self.maximum_angle_radians
            || self.maximum_angle_radians >= FRAC_PI_2
            || self.maximum_capillary_number > 0.01
            || self.maximum_weber_number > 0.01
            || self.force_relative_tolerance >= 1.0
            || self.energy_relative_tolerance >= 1.0
            || self.maximum_iterations == 0
            || self.maximum_backtracks == 0
            || self.maximum_shape_evaluations == 0
        {
            return Err(ContactLineError::Invalid(
                "axisymmetric dynamic-wetting bounds",
            ));
        }
        Ok(())
    }
}
