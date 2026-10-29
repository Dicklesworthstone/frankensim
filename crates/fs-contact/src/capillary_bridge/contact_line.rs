//! Rate-dependent wetting and pinning of the original planar liquid bridge.
//!
//! The circular menisci equilibrate at fixed volume while their footprint is
//! an independent dissipative coordinate. Four identical contact lines move
//! at half the full-span speed. Their linear low-speed law is
//! zeta*u = gamma*(cos(theta_Y)-cos(theta)), supplemented by separately sourced
//! advancing/receding dry-friction thresholds. Thus the full-span resistance is
//! zeta*width*span_rate and its power is always nonnegative.
//!
//! This is the symmetric, gravity-free, two-dimensional, contact-line-dominated
//! model. It does not infer line friction from bulk viscosity, resolve the
//! liquid's inertia/bulk flow, or turn an inadmissible neck into pinch-off.
//! Linear line friction is a bounded near-equilibrium phenomenology; see
//! https://doi.org/10.1038/s41526-021-00186-0, equation (2). Pinning is a declared
//! rate-independent dissipation model using the supplied advancing/receding
//! angles, consistent with contact-line resistance as in
//! https://doi.org/10.1021/acs.langmuir.3c01550. Neither source validates this
//! restricted geometry for an arbitrary real bridge.
//!
//! A symmetric coordinate-increment discrete gradient transfers surface work
//! to the actual gap and footprint. Bounded Gauss integration evaluates the
//! gradients without subtracting nearly equal energies. Analytic Hessians give
//! the implicit mean-force tangent. Quadrature/root estimates and the original
//! energy budget must pass; they are numerical estimates, not certificates.

use super::{BridgeError, SourcePlanarCapillaryBridge};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_matdb::InterfaceSystemCard;
use fs_material::liquid::DYNAMIC_VISCOSITY_DIMS;
use fs_material::liquid::surface::adhesion::{
    EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, EquilibriumAngleOrigin,
};
use fs_material::liquid::surface::{
    EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, ResolvedLiquidSurfaceState, ResolvedLiquidWettingState,
    SURFACE_TENSION_DIMS, require_frozen_coefficient_support,
};
use fs_material::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_math::det;
use fs_qty::{Dims, QuantitySpec};

mod geometry;

/// Work-conjugate coupling of dynamic wetting to original mechanical compliance.
pub mod load_line;

pub const CONTACT_LINE_FRICTION_PROPERTY: &str = "contact-line-friction";
pub const ADVANCING_CONTACT_ANGLE_PROPERTY: &str = "advancing-contact-angle";
pub const RECEDING_CONTACT_ANGLE_PROPERTY: &str = "receding-contact-angle";

/// One original ordered wetting bundle, with all active line-law receipts.
pub fn resolve_contact_line_wetting(
    card: &InterfaceSystemCard,
    liquid: &ResolvedLiquidSurfaceState,
    origin: EquilibriumAngleOrigin,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidWettingState, MaterialStatePointError> {
    let (name, dims) = match origin {
        EquilibriumAngleOrigin::SourcedAngle => (EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, Dims::NONE),
        EquilibriumAngleOrigin::YoungDupre => {
            (EQUILIBRIUM_WORK_OF_ADHESION_PROPERTY, SURFACE_TENSION_DIMS)
        }
    };
    let requirements = [
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)?,
        ScalarPropertyRequirement::try_new(
            CONTACT_LINE_FRICTION_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            ADVANCING_CONTACT_ANGLE_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            RECEDING_CONTACT_ANGLE_PROPERTY,
            Dims::NONE,
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    let resolved =
        resolve_interface_state_point(card, &liquid.material().query(), &requirements, selection)?;
    let wetting = match origin {
        EquilibriumAngleOrigin::SourcedAngle => {
            ResolvedLiquidWettingState::try_from_resolved(card, resolved)
        }
        EquilibriumAngleOrigin::YoungDupre => {
            ResolvedLiquidWettingState::try_from_adhesion_work(card, resolved, liquid)
        }
    }?;
    wetting.require_liquid_state(liquid)?;
    Ok(wetting)
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineConfig {
    pub minimum_angle_radians: f64,
    pub maximum_angle_radians: f64,
    pub minimum_neck_width_m: f64,
    pub maximum_gap_increment_m: f64,
    pub maximum_span_increment_m: f64,
    pub maximum_step_s: f64,
    pub maximum_contact_line_speed_m_s: f64,
    pub maximum_angle_departure_radians: f64,
    pub maximum_iterations: usize,
    pub maximum_geometry_iterations: usize,
    /// Inverse-geometry/profile solves within ONE `advance`, including every
    /// quadrature node. Each solve retains `maximum_geometry_iterations`.
    /// A mechanical wrapper separately bounds its number of `advance` calls.
    pub maximum_evaluations: usize,
    pub geometry_absolute_tolerance_m: f64,
    pub geometry_relative_tolerance: f64,
    pub force_absolute_tolerance_n: f64,
    pub force_relative_tolerance: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ContactLineError {
    Original(BridgeError),
    Material(MaterialStatePointError),
    Invalid(&'static str),
    StateMismatch,
    Cancelled,
    Budget,
    /// The original connected acute-angle geometry or caller's regime ended.
    Domain(&'static str),
    NotConverged,
    Accuracy {
        estimated_error_n: f64,
        tolerance_n: f64,
    },
    Energy {
        residual_j: f64,
        tolerance_j: f64,
    },
}
impl core::fmt::Display for ContactLineError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "dynamic planar contact line refused: {self:?}")
    }
}
impl std::error::Error for ContactLineError {}
impl From<BridgeError> for ContactLineError {
    fn from(value: BridgeError) -> Self {
        Self::Original(value)
    }
}
impl From<MaterialStatePointError> for ContactLineError {
    fn from(value: MaterialStatePointError) -> Self {
        Self::Material(value)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContactLineRegime {
    Pinned,
    Advancing,
    Receding,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineProfile {
    pub gap_m: f64,
    /// Actual, independently evolved full wetted span on each plate.
    pub wetted_span_m: f64,
    pub contact_angle_radians: f64,
    pub neck_width_m: f64,
    pub liquid_volume_m3: f64,
    pub liquid_gas_area_m2: f64,
    pub wetted_area_m2: f64,
    pub liquid_minus_gas_pressure_pa: f64,
    pub surface_free_energy_j: f64,
    /// Normal force at FIXED footprint, along increasing plate separation.
    pub normal_force_n: f64,
    pub footprint_energy_gradient_n: f64,
    pub geometric_span_residual_m: f64,
    pub estimated_volume_error_m3: f64,
    pub estimated_surface_energy_error_j: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineState {
    identity: ContentHash,
    profile: ContactLineProfile,
}
impl ContactLineState {
    pub const fn gap_m(self) -> f64 {
        self.profile.gap_m
    }
    pub const fn wetted_span_m(self) -> f64 {
        self.profile.wetted_span_m
    }
    pub const fn contact_angle_radians(self) -> f64 {
        self.profile.contact_angle_radians
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContactLineStep {
    state: ContactLineState,
    pub profile: ContactLineProfile,
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
    pub evaluations: usize,
}
impl ContactLineStep {
    pub const fn state(self) -> ContactLineState {
        self.state
    }
}

#[derive(Debug)]
pub struct SourceContactLineBridge {
    original: SourcePlanarCapillaryBridge,
    config: ContactLineConfig,
    identity: ContentHash,
    zeta: f64,
    advancing_threshold_n: f64,
    receding_threshold_n: f64,
}

impl SourceContactLineBridge {
    pub fn new(
        original: SourcePlanarCapillaryBridge,
        config: ContactLineConfig,
    ) -> Result<Self, ContactLineError> {
        config.validate()?;
        let equilibrium = original.wetting()[0].contact_angle_radians();
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
                            "missing original contact-line source property",
                        ))?;
                if property.requirement().quantity() != QuantitySpec::dimensional(dims) {
                    return Err(ContactLineError::Invalid("contact-line source quantity"));
                }
                require_frozen_coefficient_support(property)?;
                values[i][j] = property.value_si();
            }
        }
        if values[0] != values[1]
            || values[0].iter().any(|v| !v.is_finite() || *v <= 0.0)
            || !(config.minimum_angle_radians <= values[0][2]
                && values[0][2] <= equilibrium
                && equilibrium <= values[0][1]
                && values[0][1] <= config.maximum_angle_radians)
        {
            return Err(ContactLineError::Invalid(
                "symmetric positive line friction and ordered receding/equilibrium/advancing angles required",
            ));
        }
        let scale = 2.0 * original.geometry.sigma_width;
        let advancing_threshold_n = if values[0][1] == equilibrium {
            0.0
        } else {
            scale * (original.geometry.cosine - det::cos(values[0][1]))
        };
        let receding_threshold_n = if values[0][2] == equilibrium {
            0.0
        } else {
            scale * (det::cos(values[0][2]) - original.geometry.cosine)
        };
        if advancing_threshold_n < 0.0 || receding_threshold_n < 0.0 {
            return Err(ContactLineError::Invalid(
                "unrepresentable pinning thresholds",
            ));
        }
        let mut hash = DomainHasher::new("org.frankensim.fs-contact.contact-line-bridge.v1");
        hash.update(original.identity().as_bytes());
        for value in config.numbers().into_iter().chain(values[0]) {
            hash.update(&value.to_bits().to_le_bytes());
        }
        for n in [
            config.maximum_iterations,
            config.maximum_geometry_iterations,
            config.maximum_evaluations,
        ] {
            hash.update(&(n as u64).to_le_bytes());
        }
        Ok(Self {
            original,
            config,
            identity: hash.finalize(),
            zeta: values[0][0],
            advancing_threshold_n,
            receding_threshold_n,
        })
    }

    pub const fn original(&self) -> &SourcePlanarCapillaryBridge {
        &self.original
    }
    pub const fn config(&self) -> ContactLineConfig {
        self.config
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn line_friction_pa_s(&self) -> f64 {
        self.zeta
    }

    /// The supplied angle is an actual initial condition; its excess surface
    /// energy is retained, not confused with an equilibrium material property.
    pub fn initialize(
        &self,
        gap_m: f64,
        angle_radians: f64,
    ) -> Result<ContactLineState, ContactLineError> {
        let span = self.span_at_angle(gap_m, angle_radians)?;
        let evaluation = self.at_angle(gap_m, span, angle_radians)?;
        Ok(ContactLineState {
            identity: self.identity,
            profile: evaluation.profile,
        })
    }

    /// Authenticate an already admitted source state without a numerical solve.
    pub fn evaluate(
        &self,
        state: ContactLineState,
    ) -> Result<ContactLineProfile, ContactLineError> {
        if state.identity != self.identity {
            return Err(ContactLineError::StateMismatch);
        }
        Ok(state.profile)
    }

    /// One finite prescribed-gap proposal. Its state remains unpublished until
    /// the receiving mechanical owner accepts its own actual-motion work.
    pub fn advance(
        &self,
        previous: ContactLineState,
        gap_m: f64,
        dt_s: f64,
        gate: &CancelGate,
    ) -> Result<ContactLineStep, ContactLineError> {
        poll(gate)?;
        let p = self.evaluate(previous)?;
        if !dt_s.is_finite() || dt_s <= 0.0 || dt_s > self.config.maximum_step_s {
            return Err(ContactLineError::Domain(
                "positive admitted physical timestep",
            ));
        }
        if !gap_m.is_finite() || (gap_m - p.gap_m).abs() > self.config.maximum_gap_increment_m {
            return Err(ContactLineError::Domain("physical gap increment"));
        }
        self.original.check_gap(gap_m)?;
        let mut evaluations = 0;
        let b0 = p.wetted_span_m;
        // Both coordinate orders are real connected geometries. Refuse a step
        // whose rectangle leaves this branch rather than substituting a state.
        let mut lower = self
            .span_at_angle(p.gap_m, self.config.maximum_angle_radians)?
            .max(self.span_at_angle(gap_m, self.config.maximum_angle_radians)?)
            .max(b0 - self.config.maximum_span_increment_m);
        let mut upper = self
            .span_at_angle(p.gap_m, self.config.minimum_angle_radians)?
            .min(self.span_at_angle(gap_m, self.config.minimum_angle_radians)?)
            .min(b0 + self.config.maximum_span_increment_m);
        if !(lower <= b0 && b0 <= upper) {
            return Err(ContactLineError::Domain(
                "connected coordinate-increment rectangle",
            ));
        }
        let friction = self.zeta * self.original.geometry.width / dt_s;
        let initial = self.footprint_gradient(p.gap_m, gap_m, b0, b0, &mut evaluations, gate)?;
        let (regime, threshold) = if initial.value < -self.advancing_threshold_n {
            lower = b0;
            (ContactLineRegime::Advancing, self.advancing_threshold_n)
        } else if initial.value > self.receding_threshold_n {
            upper = b0;
            (ContactLineRegime::Receding, -self.receding_threshold_n)
        } else {
            (ContactLineRegime::Pinned, -initial.value)
        };
        let mut b = b0;
        let mut gradient = initial;
        let mut residual = initial.value + threshold;
        let mut iterations = 0;
        if regime != ContactLineRegime::Pinned {
            for iteration in 0..self.config.maximum_iterations {
                iterations = iteration;
                let tolerance =
                    self.force_tolerance(gradient.value.abs().max((friction * (b - b0)).abs()));
                if residual.abs() + gradient.error <= tolerance {
                    break;
                }
                if residual < 0.0 {
                    lower = b;
                } else {
                    upper = b;
                }
                let derivative = friction + gradient.db;
                if !derivative.is_finite() || derivative <= 0.0 {
                    return Err(ContactLineError::Invalid(
                        "positive contact-line residual derivative",
                    ));
                }
                let candidate = b - residual / derivative;
                let next = if candidate > lower && candidate < upper {
                    candidate
                } else {
                    lower + 0.5 * (upper - lower)
                };
                if next == b || !(next > lower && next < upper) {
                    return Err(ContactLineError::NotConverged);
                }
                b = next;
                gradient =
                    self.footprint_gradient(p.gap_m, gap_m, b0, b, &mut evaluations, gate)?;
                residual = friction * (b - b0) + gradient.value + threshold;
            }
            let tolerance =
                self.force_tolerance(gradient.value.abs().max((friction * (b - b0)).abs()));
            if residual.abs() + gradient.error > tolerance {
                return Err(ContactLineError::NotConverged);
            }
        } else if gradient.error > self.force_tolerance(gradient.value.abs()) {
            return Err(ContactLineError::Accuracy {
                estimated_error_n: gradient.error,
                tolerance_n: self.force_tolerance(gradient.value.abs()),
            });
        }
        let db = b - b0;
        let dh = gap_m - p.gap_m;
        let endpoint = self.at(gap_m, b, &mut evaluations, gate)?;
        let normal = self.gap_gradient(p.gap_m, gap_m, b0, b, &mut evaluations, gate)?;
        let db_dh = if regime == ContactLineRegime::Pinned {
            0.0
        } else {
            -gradient.dh / (friction + gradient.db)
        };
        let mean_force = -normal.value;
        let tangent = -(normal.dh + normal.db * db_dh);
        let footprint_error = if regime == ContactLineRegime::Pinned {
            0.0
        } else {
            (residual.abs() + gradient.error) / (friction + gradient.db)
        };
        let force_error = normal.error + normal.db.abs() * footprint_error;
        let force_tolerance = self.force_tolerance(mean_force.abs());
        if !mean_force.is_finite()
            || !tangent.is_finite()
            || !force_error.is_finite()
            || force_error > force_tolerance
        {
            return Err(ContactLineError::Accuracy {
                estimated_error_n: force_error,
                tolerance_n: force_tolerance,
            });
        }
        let viscous = friction * db * db;
        let pinning = threshold * db;
        let dissipation = viscous + pinning;
        let work = normal.value * dh;
        let change = endpoint.profile.surface_free_energy_j - p.surface_free_energy_j;
        let energy_residual = change + dissipation - work;
        let limits = self.original.limits();
        let energy_scale = p
            .surface_free_energy_j
            .abs()
            .max(endpoint.profile.surface_free_energy_j.abs())
            .max(work.abs())
            .max(dissipation);
        let energy_tolerance =
            limits.energy_absolute_tolerance_j + limits.energy_relative_tolerance * energy_scale;
        let energy_error = endpoint.profile.estimated_surface_energy_error_j
            + p.estimated_surface_energy_error_j
            + normal.error * dh.abs()
            + gradient.error * db.abs();
        if !dissipation.is_finite()
            || dissipation < 0.0
            || !work.is_finite()
            || !energy_error.is_finite()
            || energy_residual.abs() + energy_error > energy_tolerance
        {
            return Err(ContactLineError::Energy {
                residual_j: energy_residual,
                tolerance_j: energy_tolerance,
            });
        }
        super::limit(
            work.abs(),
            limits.maximum_abs_step_work_j,
            "dynamic wetting work",
        )?;
        super::limit(
            dissipation,
            limits.maximum_abs_step_work_j,
            "dynamic wetting loss",
        )?;
        super::limit(
            mean_force.abs(),
            limits.maximum_abs_force_n,
            "dynamic wetting mean force",
        )?;
        let speed = (0.5 * db / dt_s).abs();
        if speed > self.config.maximum_contact_line_speed_m_s {
            return Err(ContactLineError::Domain("contact-line speed"));
        }
        let physical_speed = speed.max((dh / dt_s).abs());
        let sigma = self.original.surface_law().surface_tension_n_m();
        let liquid = self.original.liquid();
        super::limit(
            liquid.dynamic_viscosity_pa_s() * physical_speed / sigma,
            limits.maximum_capillary_number,
            "dynamic wetting capillary number",
        )?;
        let length = p.gap_m.max(gap_m).max(b0).max(b);
        super::limit(
            liquid.density_kg_m3() * physical_speed * physical_speed * length / sigma,
            limits.maximum_weber_number,
            "dynamic wetting Weber number",
        )?;
        poll(gate)?;
        Ok(ContactLineStep {
            state: ContactLineState {
                identity: self.identity,
                profile: endpoint.profile,
            },
            profile: endpoint.profile,
            regime,
            mean_normal_force_n: mean_force,
            mean_force_gap_derivative_n_per_m: tangent,
            surface_energy_change_j: change,
            contact_line_dissipation_j: dissipation,
            viscous_line_dissipation_j: viscous,
            pinning_dissipation_j: pinning,
            work_on_liquid_j: work,
            energy_balance_residual_j: energy_residual,
            footprint_force_residual_n: residual,
            estimated_force_error_n: force_error,
            estimated_energy_error_j: energy_error,
            iterations,
            evaluations,
        })
    }

    fn force_tolerance(&self, scale: f64) -> f64 {
        self.config.force_absolute_tolerance_n + self.config.force_relative_tolerance * scale
    }
}

impl ContactLineConfig {
    fn numbers(self) -> [f64; 12] {
        [
            self.minimum_angle_radians,
            self.maximum_angle_radians,
            self.minimum_neck_width_m,
            self.maximum_gap_increment_m,
            self.maximum_span_increment_m,
            self.maximum_step_s,
            self.maximum_contact_line_speed_m_s,
            self.maximum_angle_departure_radians,
            self.geometry_absolute_tolerance_m,
            self.geometry_relative_tolerance,
            self.force_absolute_tolerance_n,
            self.force_relative_tolerance,
        ]
    }
    fn validate(self) -> Result<(), ContactLineError> {
        if self.numbers().iter().any(|x| !x.is_finite() || *x <= 0.0)
            || self.minimum_angle_radians >= self.maximum_angle_radians
            || self.maximum_angle_radians >= core::f64::consts::FRAC_PI_2
            || self.geometry_relative_tolerance >= 1.0
            || self.force_relative_tolerance >= 1.0
            || !(1..=128).contains(&self.maximum_iterations)
            || !(1..=128).contains(&self.maximum_geometry_iterations)
            || self.maximum_evaluations == 0
        {
            return Err(ContactLineError::Invalid(
                "finite positive contact-line and numerical budgets",
            ));
        }
        Ok(())
    }
}

fn poll(gate: &CancelGate) -> Result<(), ContactLineError> {
    if gate.is_requested() {
        Err(ContactLineError::Cancelled)
    } else {
        Ok(())
    }
}
