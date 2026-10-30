//! Material interface laws integrated over a rigid, rectangular bonded pad.
//!
//! Local x/y are the pad's centroidal principal axes; z is its normal. Three
//! original traction/separation laws supply k and c per area. Integrating
//! `u(r) = translation + rotation cross r` over the rectangle gives the native
//! bushing's three translational and three rotational coefficients. No bulk
//! modulus, loss factor, or rotational stiffness is guessed.
//!
//! This is a bilateral, massless, small-local-rotation layer approximation.
//! Common rigid motion can be large. There is no debonding, unilateral contact,
//! cross-axis constitutive coupling, heat store, or automatic thermal feedback.
//! The native bushing owns the work-conjugate finite-rotation force map and
//! explicit/implicit dynamics. Source rebinding reports its actual storage work.

use crate::articulated::ResolvedArticulated;
use fs_ga::{GaError, Se3, Vec3};
use fs_matdb::{InterfaceSystemCard, QueryPoint};
use fs_material::state_point::{
    MaterialStatePointError, ResolvedInterfaceStatePoint, ResolvedLinearInterface,
    is_fixed_external_coordinate,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{
    BushingForces, MaterialFrame, SpatialBushing, evaluate_bushing_forces,
    step_free_floating_bushings_implicit,
};
use fs_mbd::articulated::motion::spatial_loads::implicit::{
    ImplicitSpatialConfig, ImplicitSpatialError, ImplicitSpatialStep,
};
use fs_mbd::articulated::motion::spatial_loads::{MAX_SPATIAL_MOTION_LOADS, SpatialMotionLoads};
use fs_mbd::articulated::motion::{FreeFloatingMotionError, FreeFloatingMotionState};
use fs_mbd::articulated::{ArticulatedError, BaseState, forward_kinematics};
use fs_qty::semantic::QuantityKind;

/// Original attachment poses and the admitted physical layer approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularPad {
    /// Ordered surface frames, with aligned material axes in the relaxed state.
    pub frames: [MaterialFrame; 2],
    /// Relaxed B-origin position in A; never inferred from a strained state.
    pub rest_translation_a: Vec3,
    /// Full rectangular widths in local x and y [m].
    pub size_m: [f64; 2],
    /// Caller-declared bound on relative origin displacement from rest [m].
    pub maximum_displacement_m: f64,
    /// Caller-declared local rotation bound in (0, 0.1] rad. The upper cap
    /// limits this linearized patch geometry; it is not material qualification.
    pub maximum_rotation_rad: f64,
}

/// No failed construction, state evaluation or rebinding changes the input.
#[derive(Debug)]
pub enum AttachmentError {
    Invalid(&'static str),
    Material(MaterialStatePointError),
    Kinematics(ArticulatedError),
    Frame(GaError),
    Motion(FreeFloatingMotionError),
    Implicit(ImplicitSpatialError),
    Cancelled,
}
impl core::fmt::Display for AttachmentError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material attachment: {self:?}")
    }
}
impl std::error::Error for AttachmentError {}
impl From<MaterialStatePointError> for AttachmentError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Material(error)
    }
}

/// One native attachment and its original ordered material/geometry owners.
#[derive(Clone, Debug)]
pub struct MaterialBushing<'a> {
    mechanism: &'a ResolvedArticulated,
    regions: [String; 2],
    geometry: RectangularPad,
    axes: [ResolvedLinearInterface; 3],
    bushing: SpatialBushing,
}

/// A new binding at the same mechanical state, with explicit external work.
#[derive(Clone, Debug)]
pub struct AttachmentRebind<'a> {
    pub binding: MaterialBushing<'a>,
    /// U_new - U_old [J], positive for work supplied to elastic storage.
    /// It is not damping heat or a complete thermoelastic heat model.
    pub parameter_work_j: f64,
}

impl<'a> MaterialBushing<'a> {
    /// Bind x/y shear and z normal interface properties to the actual two
    /// material regions. Each `(stiffness, viscosity)` pair names claims in
    /// ONE resolved bundle; all source selections and query coordinates stay
    /// attached. A sourced zero permits a pure spring or damper on an axis.
    ///
    /// Patch frames/coverage and linear-law applicability are caller-declared.
    /// Source frequency support must cover the application's actual response;
    /// `require_angular_frequency_rad_s` checks that shared intersection.
    pub fn rectangular_pad(
        mechanism: &'a ResolvedArticulated,
        regions: [&str; 2],
        geometry: RectangularPad,
        source: ResolvedInterfaceStatePoint,
        properties: [(&str, &str); 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, AttachmentError> {
        poll(&mut cancelled)?;
        let [x, y] = geometry.size_m;
        let area = x * y;
        if ![
            x,
            y,
            area,
            geometry.maximum_displacement_m,
            geometry.maximum_rotation_rad,
        ]
        .into_iter()
        .all(positive)
            || geometry.maximum_rotation_rad > 0.1
            || !finite_vector(geometry.rest_translation_a)
            || geometry.frames[0].link == geometry.frames[1].link
        {
            return Err(AttachmentError::Invalid(
                "positive pad geometry, finite rest pose and distinct links required",
            ));
        }
        let mut materials = Vec::with_capacity(2);
        for (frame, key) in geometry.frames.into_iter().zip(regions) {
            let region = mechanism
                .bodies()
                .get(frame.link)
                .and_then(|body| body.regions().iter().find(|r| r.key() == key))
                .ok_or(AttachmentError::Invalid(
                    "pad frame must name an original link material region",
                ))?;
            materials.push(region.material().material());
        }
        let mut axes = Vec::with_capacity(3);
        for (stiffness, viscosity) in properties {
            poll(&mut cancelled)?;
            let axis = ResolvedLinearInterface::new(source.clone(), stiffness, viscosity, area)?;
            axis.require_surface_materials(materials[0], materials[1])?;
            axis.require_angular_frequency_rad_s(0.0)?;
            admit_fixed_context(&axis)?;
            if axis.constant_loss_reference_omega_rad_s().is_some() {
                return Err(AttachmentError::Invalid(
                    "native Kelvin-Voigt attachments do not realize constant loss",
                ));
            }
            axes.push(axis);
        }
        let axes: [ResolvedLinearInterface; 3] =
            axes.try_into().expect("three original interface axes");
        let stiffness = axes.each_ref().map(ResolvedLinearInterface::stiffness_n_m);
        let damping = axes.each_ref().map(ResolvedLinearInterface::damping_n_s_m);
        let bushing = SpatialBushing {
            a: geometry.frames[0],
            b: geometry.frames[1],
            rest_translation_a: geometry.rest_translation_a,
            translation_stiffness: vector(stiffness),
            translation_damping: vector(damping),
            rotation_stiffness: rotational(stiffness, geometry.size_m)?,
            rotation_damping: rotational(damping, geometry.size_m)?,
        };
        poll(&mut cancelled)?;
        Ok(Self {
            mechanism,
            regions: regions.map(str::to_owned),
            geometry,
            axes,
            bushing,
        })
    }

    pub const fn mechanism(&self) -> &'a ResolvedArticulated {
        self.mechanism
    }
    pub const fn geometry(&self) -> RectangularPad {
        self.geometry
    }
    /// Original x/y/z interface laws, including source identity and selections.
    pub const fn axes(&self) -> &[ResolvedLinearInterface; 3] {
        &self.axes
    }

    /// Check the frequency against every original axis, without inferring an
    /// application's modal/excitation bandwidth from its sample rate.
    pub fn require_angular_frequency_rad_s(&self, omega: f64) -> Result<(), AttachmentError> {
        for axis in &self.axes {
            axis.require_angular_frequency_rad_s(omega)?;
        }
        Ok(())
    }

    /// Export the original native bushing after checking this actual state.
    /// Use with the SAME mechanism's existing solver. [`step_bushings_implicit`]
    /// also admits the actual accepted midpoint and endpoint automatically.
    /// Direct native stepping needs those checks from its caller. Neither path
    /// encloses the continuous trajectory or infers its full source bandwidth.
    pub fn for_state(
        &self,
        state: &FreeFloatingMotionState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialBushing, AttachmentError> {
        poll(&mut cancelled)?;
        let kinematics = forward_kinematics(
            self.mechanism.model(),
            BaseState::stationary(Se3::identity()),
            &state.joint_position,
            &state.joint_velocity,
        )
        .map_err(AttachmentError::Kinematics)?;
        let [a, b] = self.geometry.frames.map(|frame| {
            kinematics.world_from_link[frame.link].compose(frame.link_from_attachment)
        });
        let a = a.map_err(AttachmentError::Frame)?;
        let b = b.map_err(AttachmentError::Frame)?;
        let relative = a
            .inverse()
            .map_err(AttachmentError::Frame)?
            .compose(b)
            .map_err(AttachmentError::Frame)?;
        let displacement = (relative.translation() - self.geometry.rest_translation_a).norm();
        let angle = relative.rotation().log().angular.norm();
        if !displacement.is_finite()
            || !angle.is_finite()
            || displacement > self.geometry.maximum_displacement_m
            || angle > self.geometry.maximum_rotation_rad
        {
            return Err(AttachmentError::Invalid(
                "state exceeds the declared local pad approximation",
            ));
        }
        poll(&mut cancelled)?;
        Ok(self.bushing)
    }

    /// Evaluate original reciprocal wrenches, elastic storage and physical loss.
    pub fn evaluate(
        &self,
        state: &FreeFloatingMotionState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BushingForces, AttachmentError> {
        let bushing = self.for_state(state, &mut cancelled)?;
        evaluate_bushing_forces(self.mechanism.model(), state, &[bushing], &mut cancelled)
            .map_err(AttachmentError::Motion)
    }

    /// Atomically refresh the original source/query at the ACTUAL mechanical
    /// state. Original geometry, ordered materials and selected claims persist.
    /// The signed native storage difference is explicit external parameter work;
    /// the caller must account for it before adopting the new binding. A source
    /// update alone does not transfer momentum, advance motion, or deposit heat.
    pub fn rebind_at_state(
        &self,
        card: &InterfaceSystemCard,
        point: &QueryPoint,
        state: &FreeFloatingMotionState,
        maximum_parameter_work_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AttachmentRebind<'a>, AttachmentError> {
        poll(&mut cancelled)?;
        if !maximum_parameter_work_j.is_finite() || maximum_parameter_work_j < 0.0 {
            return Err(AttachmentError::Invalid(
                "finite nonnegative parameter-work budget required",
            ));
        }
        let before = self.evaluate(state, &mut cancelled)?.stored_energy;
        let source = self.axes[0].resolved().at_query_point(card, point)?;
        let binding = Self::rectangular_pad(
            self.mechanism,
            self.regions.each_ref().map(String::as_str),
            self.geometry,
            source,
            self.axes
                .each_ref()
                .map(|a| (a.stiffness_property(), a.viscosity_property())),
            &mut cancelled,
        )?;
        let work = binding.evaluate(state, &mut cancelled)?.stored_energy - before;
        if !work.is_finite() || work.abs() > maximum_parameter_work_j {
            return Err(AttachmentError::Invalid(
                "source rebinding exceeds the elastic parameter-work budget",
            ));
        }
        poll(&mut cancelled)?;
        Ok(AttachmentRebind {
            binding,
            parameter_work_j: work,
        })
    }
}

/// Advance material-bound attachments together through the original implicit
/// owner. All bindings must borrow the SAME original mechanism. Their physical
/// envelopes are checked at the input, exact accepted native midpoint and final
/// state before the unchanged native report is returned. Any refusal publishes
/// no motion or material change. No sequential per-attachment solve is performed.
///
/// This admits the accepted samples, not the continuous path between them or
/// every nonlinear trial/stencil. The source parameters stay fixed during the
/// step; the application's complete response frequency band remains explicit.
pub fn step_bushings_implicit(
    bindings: &[&MaterialBushing<'_>],
    state: &FreeFloatingMotionState,
    loads: SpatialMotionLoads<'_>,
    duration: f64,
    config: ImplicitSpatialConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ImplicitSpatialStep, AttachmentError> {
    poll(&mut cancelled)?;
    let first = bindings.first().ok_or(AttachmentError::Invalid(
        "at least one sourced attachment required",
    ))?;
    if bindings.len() > MAX_SPATIAL_MOTION_LOADS {
        return Err(AttachmentError::Invalid(
            "source attachments exceed the native spatial load bound",
        ));
    }
    let mut native = Vec::with_capacity(bindings.len());
    for binding in bindings {
        if !core::ptr::eq(first.mechanism, binding.mechanism) {
            return Err(AttachmentError::Invalid(
                "attachments belong to different original mechanisms",
            ));
        }
        native.push(binding.for_state(state, &mut cancelled)?);
    }
    let step = step_free_floating_bushings_implicit(
        first.mechanism.model(),
        state,
        loads,
        &native,
        duration,
        config,
        &mut cancelled,
    )
    .map_err(AttachmentError::Implicit)?;
    for binding in bindings {
        binding.for_state(&step.midpoint_state, &mut cancelled)?;
        binding.for_state(&step.spatial.motion.state_after, &mut cancelled)?;
    }
    poll(&mut cancelled)?;
    Ok(step)
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}
fn finite_vector(v: Vec3) -> bool {
    [v.x, v.y, v.z].into_iter().all(f64::is_finite)
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AttachmentError> {
    if cancelled() {
        Err(AttachmentError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit_fixed_context(axis: &ResolvedLinearInterface) -> Result<(), AttachmentError> {
    let source = axis.resolved();
    for name in [axis.stiffness_property(), axis.viscosity_property()] {
        let property = source.property(name).expect("original interface property");
        for coordinate in property.answer().evidence.model.validity.bounds().keys() {
            let quantity = source.axis_quantities().get(coordinate).copied();
            // The original linear-interface owner has already admitted the
            // unique static frequency axis and its full applicability band.
            let frequency = coordinate == "omega"
                || quantity
                    .and_then(|q| q.semantic_type())
                    .is_some_and(|ty| matches!(ty.kind(), QuantityKind::Frequency(_)));
            if !frequency
                && !is_fixed_external_coordinate(
                    coordinate,
                    quantity,
                    source.state_coordinate(coordinate),
                )
            {
                return Err(AttachmentError::Invalid(
                    "pad coefficients cannot freeze unowned strain, gap, rate or field applicability",
                ));
            }
        }
    }
    Ok(())
}
fn rotational(coefficients: [f64; 3], size: [f64; 2]) -> Result<Vec3, AttachmentError> {
    // Integral x² dA / A and y² dA / A; first/product moments vanish at the
    // rectangle's centroid in its principal axes. Translation uses full area.
    let [xx, yy] = size.map(|s| s * (s / 12.0));
    if !positive(xx) || !positive(yy) {
        return Err(AttachmentError::Invalid(
            "pad area moments are unrepresentable",
        ));
    }
    let moment = |c: f64, arm: f64| {
        let value = c * arm;
        if !value.is_finite() || (c > 0.0 && value == 0.0) {
            Err(AttachmentError::Invalid(
                "rotational interface coefficient is unrepresentable",
            ))
        } else {
            Ok(value)
        }
    };
    let result = Vec3::new(
        moment(coefficients[2], yy)?,
        moment(coefficients[2], xx)?,
        moment(coefficients[0], yy)? + moment(coefficients[1], xx)?,
    );
    if !finite_vector(result) {
        return Err(AttachmentError::Invalid(
            "pad torsional coefficient overflows",
        ));
    }
    Ok(result)
}
