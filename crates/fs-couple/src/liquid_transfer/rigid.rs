//! Release a source-bound rigid liquid into the existing moving-liquid owner.
//!
//! The original tetrahedral matter, principal frame and rigid momenta determine
//! every world vertex and velocity. Native tet/cell overlap, per-piece momentum
//! correction, viscous initialization and caloric transfer remain in `fs-lbm`.
//! No inertia, velocity field, phase fraction or heat is guessed here.
//!
//! This removes an ideal kinematic constraint on material ALREADY fully liquid.
//! The represented rigid model has no stored elastic/plastic/history energy.
//! It cannot accept a deformable solid owner or authorize erasing its memory.
//! Melting kinetics, stress relaxation, density-change volume work, capillarity
//! and a pressure projection need their own admitted physical owners.

mod kinetic;
/// Selected-cell release with a source-bound rigid remainder outside the grid.
pub mod partial;

use fs_exec::{CancelGate, Cx};
use fs_lbm::material::free_surface::enthalpy::rheology::RheologySolve;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceEnthalpyControls, FreeSurfaceEnthalpyError, MaterialFreeSurfaceEnthalpy3,
};
use fs_lbm::material::free_surface::transfer::enthalpy::{
    TetEnthalpyTransferError, TetEnthalpyTransferLimits, TetEnthalpyTransferReport,
    TetrahedralLiquidEnthalpyState,
};
use fs_lbm::material::free_surface::transfer::{LiquidTransferMoments, TetrahedralLiquidState};
use fs_lbm::material::free_surface::{MaterialFreeSurfaceDomain3, MaterialFreeSurfaceScaling};
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::phase::EquilibriumEnthalpyPhaseCurve;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, RigidBodyState, Vec3};
use fs_mbd_material::density_field::TetrahedralDensityBody;
use fs_mbd_material::principal::{PrincipalFrame, PrincipalFrameError};

/// One accepted no-storage rigid constraint and its original material matter.
///
/// `frame` must be the actual principal frame used to evolve `motion`, derived
/// from this exact geometry's inertia. Geometry vertices and nodal enthalpies
/// retain their original order. The caloric state must already be fully liquid;
/// this type does not represent a solid's elastic, plastic or damage history.
#[derive(Clone, Copy, Debug)]
pub struct RigidTetrahedralLiquidState<'a> {
    pub geometry: &'a TetrahedralDensityBody,
    pub frame: &'a PrincipalFrame,
    pub motion: RigidBodyState,
    /// Specific enthalpy [J/kg] on original reference vertices, affine per tet.
    pub specific_enthalpy_j_kg: &'a [f64],
    /// Accepted physical clock [s], preserved by the unadvanced receiver.
    pub time_s: f64,
}

/// Original receiver budgets and independent rigid-to-P1 identity gates.
#[derive(Clone, Copy, Debug)]
pub struct RigidLiquidReleaseControls {
    /// Shared overlap, mass, momenta, geometric projection and caloric budgets.
    /// Mass/P/L/h gates also compare the actual receiver with original rigid
    /// matter, so intermediate roundoff cannot spend these budgets twice.
    pub transfer: TetEnthalpyTransferLimits,
    pub thermal: FreeSurfaceEnthalpyControls,
    /// Absolute native rigid vs integrated P1 kinetic mismatch [J]. This
    /// checks the SAME motion representation, not an allowed physical loss.
    pub source_kinetic_tolerance_j: f64,
    /// Native rigid vs integrated P1 mass-centroid mismatch [m].
    pub source_centroid_tolerance_m: f64,
    /// Original generalized-Newtonian owner, required by an evolving shear law.
    pub rheology: Option<RheologySolve>,
    /// Original bounded source-density buoyancy; reference mass stays fixed.
    pub buoyancy: Option<BoussinesqConfig<3>>,
}

/// Which actual representation failed to reproduce the original rigid state.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RigidLiquidBalanceStage {
    /// Native selected and retained rigid states before any liquid projection.
    RigidPartition,
    TetrahedralKinematics,
    NativeLiquid,
}

/// No error changes the borrowed rigid matter, motion, h field or clock.
#[derive(Debug)]
pub enum RigidLiquidReleaseError {
    Invalid(&'static str),
    Frame(PrincipalFrameError),
    Mechanics(DynamicsError),
    Transfer(TetEnthalpyTransferError),
    Thermal(FreeSurfaceEnthalpyError),
    Balance {
        stage: RigidLiquidBalanceStage,
        mass_error_kg: f64,
        momentum_error_kg_m_s: [f64; 3],
        angular_momentum_error_kg_m2_s: [f64; 3],
        centroid_error_m: [f64; 3],
        kinetic_error_j: f64,
        gravitational_error_j: f64,
        enthalpy_error_j: f64,
    },
    Cancelled,
}

impl core::fmt::Display for RigidLiquidReleaseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid liquid release: {self:?}")
    }
}
impl std::error::Error for RigidLiquidReleaseError {}
impl From<PrincipalFrameError> for RigidLiquidReleaseError {
    fn from(error: PrincipalFrameError) -> Self {
        Self::Frame(error)
    }
}
impl From<DynamicsError> for RigidLiquidReleaseError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error)
    }
}
impl From<TetEnthalpyTransferError> for RigidLiquidReleaseError {
    fn from(error: TetEnthalpyTransferError) -> Self {
        Self::Transfer(error)
    }
}
impl From<FreeSurfaceEnthalpyError> for RigidLiquidReleaseError {
    fn from(error: FreeSurfaceEnthalpyError) -> Self {
        Self::Thermal(error)
    }
}

/// Source-native identity checks and actual liquid projection remain separate.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidLiquidReleaseReport {
    pub physical_time_s: f64,
    /// Native rigid M/P/L/K, about `transfer.mechanical.grid_origin_m`.
    pub rigid: LiquidTransferMoments,
    /// Original fixed cell masses times their affine nodal enthalpy [J].
    pub rigid_enthalpy_j: f64,
    /// Original overlap, per-piece lineage, source P1 and actual native receipt.
    pub transfer: TetEnthalpyTransferReport,
    /// Integrated P1 kinetic energy minus original rigid kinetic energy [J].
    pub source_kinetic_error_j: f64,
    /// Integrated P1 centroid minus original rigid COM [m].
    pub source_centroid_error_m: [f64; 3],
    /// Actual native receiver minus original rigid matter [J]. Both signed
    /// mechanical changes are numerical projection, never caloric heating.
    pub kinetic_projection_change_j: f64,
    pub gravitational_projection_change_j: f64,
    /// Actual receiver mass*h minus original rigid matter's mass*h [J].
    pub enthalpy_error_j: f64,
    kinetic_change: Dd,
}

impl RigidLiquidReleaseReport {
    /// Actual native liquid minus original rigid storage before rounding [J].
    /// Retains both the independently checked rigid-to-P1 discrepancy and the
    /// native remap increment, for composition with surviving rigid matter.
    pub const fn retained_kinetic_projection_change_j(&self) -> Dd {
        self.kinetic_change
    }
}

/// Return an ordinary thermal VOF candidate at the original rigid body's clock.
///
/// All nodal density bindings must retain the SAME immutable card, pinned
/// density claim, nonthermal context, quantity conventions and constant mass
/// density as the receiver. Different density fields cannot be projected away
/// by matching only their total mass. The supplied principal frame must retain
/// the exact original inertia, including products of inertia and COM offset.
/// Its native point-kinematics query derives `P/m + omega cross r` for each
/// vertex. No caller-authored velocity snapshot, copied nominal inertia or
/// zero-spin fallback can replace the accepted mechanical state.
///
/// Native transfer admits the complete P1 caloric interval as fully liquid and
/// retains its original thermal/rheological/Boussinesq sources. Each detached
/// face-connected piece keeps its own momentum and enthalpy. Receiver P/L/h
/// are additionally checked against the original rigid body, not only against
/// the intermediate tetrahedral field. Numerical kinetic and gravity changes
/// use the original projection budget and never enter h. The enthalpy datum is
/// retained; no conversion to internal energy or extra pressure work occurs.
///
/// The operation advances no time and applies no release impulse or work. It
/// removes only a rigid constraint with no stored mechanical energy. It does
/// not melt a stressed solid, detach a bonded assembly, or assign material
/// friction/heat to numerical remapping. On success the caller replaces the
/// rigid representation with this candidate; continuing both would duplicate
/// the same matter. Every refusal returns no receiver and leaves the source
/// available for a refined retry.
#[allow(clippy::too_many_arguments)]
pub fn release_rigid_liquid<'s, 'm>(
    source: RigidTetrahedralLiquidState<'_>,
    scaling: MaterialFreeSurfaceScaling<'s>,
    domain: MaterialFreeSurfaceDomain3<'_>,
    grid_origin_m: [f64; 3],
    law: &'m LiquidTemperatureLaw<'m>,
    caloric: &'m EquilibriumEnthalpyPhaseCurve,
    controls: RigidLiquidReleaseControls,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<
    (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        RigidLiquidReleaseReport,
    ),
    RigidLiquidReleaseError,
> {
    check(cx, gate)?;
    let maximum_vertices = controls
        .transfer
        .mechanical
        .max_tetrahedra
        .checked_mul(4)
        .ok_or(RigidLiquidReleaseError::Invalid(
            "source vertex budget overflows",
        ))?;
    if source.geometry.tetrahedra().len() > controls.transfer.mechanical.max_tetrahedra
        || source.geometry.vertices().len() > maximum_vertices
        || source.specific_enthalpy_j_kg.len() != source.geometry.vertices().len()
        || !source.time_s.is_finite()
        || source.time_s < 0.0
        || grid_origin_m.iter().any(|x| !x.is_finite())
        || [
            controls.source_kinetic_tolerance_j,
            controls.source_centroid_tolerance_m,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(RigidLiquidReleaseError::Invalid(
            "release needs bounded original vertices/tetrahedra, h, clock and identity tolerances",
        ));
    }
    if source.frame.inertia() != source.geometry.inertia() {
        return Err(RigidLiquidReleaseError::Invalid(
            "rigid motion must use the exact original matter's principal frame",
        ));
    }
    admit_reference_density(source.geometry, scaling, law, cx, gate)?;
    let properties = source.frame.properties();
    let rigid = rigid_moments(source.frame, source.motion, grid_origin_m)?;
    let mut positions = reserved(source.geometry.vertices().len())?;
    let mut velocities = reserved(source.geometry.vertices().len())?;
    for vertex in source.geometry.vertices() {
        check(cx, gate)?;
        let arm = source.frame.arm_for_body_point(vertex.position_m())?;
        let point = source.motion.point_kinematics(properties, arm)?;
        positions.push(array(point.point_world));
        velocities.push(array(point.point_velocity_world));
    }
    let mut tetrahedra = reserved(source.geometry.tetrahedra().len())?;
    let mut original_h = Dd::ZERO;
    for (tet, cell) in source
        .geometry
        .tetrahedra()
        .iter()
        .zip(source.geometry.cells())
    {
        check(cx, gate)?;
        let mut indices = [0; 4];
        let mut mean_h = Dd::ZERO;
        for (slot, &vertex) in tet.iter().enumerate() {
            indices[slot] = u32::try_from(vertex).map_err(|_| {
                RigidLiquidReleaseError::Invalid("original topology exceeds native vertex indices")
            })?;
            let h = source.specific_enthalpy_j_kg[vertex];
            if !h.is_finite() {
                return Err(RigidLiquidReleaseError::Invalid(
                    "original nodal enthalpy is nonfinite",
                ));
            }
            mean_h = mean_h + Dd::from_f64(h) * Dd::from_f64(0.25);
        }
        original_h = original_h + Dd::from_f64(cell.inertia().mass()) * mean_h;
        tetrahedra.push(indices);
    }
    let rigid_enthalpy_j = original_h.to_f64();
    if !rigid_enthalpy_j.is_finite() {
        return Err(RigidLiquidReleaseError::Invalid(
            "original caloric integral is unrepresentable",
        ));
    }
    let field = TetrahedralLiquidEnthalpyState {
        motion: TetrahedralLiquidState {
            positions_m: &positions,
            tetrahedra: &tetrahedra,
            velocities_m_s: &velocities,
        },
        specific_enthalpy_j_kg: source.specific_enthalpy_j_kg,
    };
    let (receiver, transfer) = match (controls.buoyancy, controls.rheology) {
        (Some(buoyancy), rheology) => scaling
            .free_surface3_from_tetrahedra_with_boussinesq_enthalpy(
                domain,
                grid_origin_m,
                field,
                law,
                caloric,
                controls.transfer,
                controls.thermal,
                buoyancy,
                rheology,
                cx,
                gate,
            ),
        (None, Some(rheology)) => scaling.free_surface3_from_tetrahedra_with_enthalpy_and_rheology(
            domain,
            grid_origin_m,
            field,
            law,
            caloric,
            controls.transfer,
            controls.thermal,
            rheology,
            cx,
            gate,
        ),
        (None, None) => scaling.free_surface3_from_tetrahedra_with_enthalpy(
            domain,
            grid_origin_m,
            field,
            law,
            caloric,
            controls.transfer,
            controls.thermal,
            cx,
            gate,
        ),
    }?;
    let reference = array(source.motion.center_of_mass_velocity_world(properties)?);
    let p1 = transfer
        .mechanical
        .source_kinetic_in_frame(reference)
        .map_err(TetEnthalpyTransferError::Transfer)?;
    let source_kinetic_change = kinetic::p1_change(
        rigid_storage(source.frame, source.motion),
        reference,
        p1.mass_kg,
        p1.relative_momentum_kg_m_s,
        p1.relative_kinetic_energy_j,
    )
    .ok_or(RigidLiquidReleaseError::Invalid(
        "original rigid-to-P1 kinetic change is unrepresentable",
    ))?;
    let source_balance = audit(
        rigid,
        rigid_enthalpy_j,
        transfer.mechanical.source,
        transfer.source_enthalpy_j,
        source_kinetic_change,
        domain.gravity_m_s2,
        controls,
        RigidLiquidBalanceStage::TetrahedralKinematics,
    )?;
    let target_balance = audit(
        rigid,
        rigid_enthalpy_j,
        transfer.mechanical.target,
        transfer.target_enthalpy_j,
        source_kinetic_change + transfer.mechanical.retained_kinetic_energy_change_j(),
        domain.gravity_m_s2,
        controls,
        RigidLiquidBalanceStage::NativeLiquid,
    )?;
    let receiver = receiver.with_time_origin(source.time_s)?;
    check(cx, gate)?;
    Ok((
        receiver,
        RigidLiquidReleaseReport {
            physical_time_s: source.time_s,
            rigid,
            rigid_enthalpy_j,
            transfer,
            source_kinetic_error_j: source_balance.kinetic_j.to_f64(),
            source_centroid_error_m: source_balance.centroid_m,
            kinetic_projection_change_j: target_balance.kinetic_j.to_f64(),
            gravitational_projection_change_j: target_balance.gravity_j,
            enthalpy_error_j: target_balance.enthalpy_j,
            kinetic_change: target_balance.kinetic_j,
        },
    ))
}

fn admit_reference_density(
    geometry: &TetrahedralDensityBody,
    scaling: MaterialFreeSurfaceScaling<'_>,
    law: &LiquidTemperatureLaw<'_>,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<(), RigidLiquidReleaseError> {
    let material = scaling.material().material();
    let density = material
        .property(DENSITY_PROPERTY)
        .ok_or(RigidLiquidReleaseError::Invalid(
            "receiver lost its original density source",
        ))?;
    for vertex in geometry.vertices() {
        check(cx, gate)?;
        let original = vertex.material();
        let original_density =
            original
                .property(DENSITY_PROPERTY)
                .ok_or(RigidLiquidReleaseError::Invalid(
                    "original vertex lost its density source",
                ))?;
        if original.card_identity() != material.card_identity()
            || original_density.answer().receipt.selected != density.answer().receipt.selected
            || vertex.density_kg_m3().to_bits() != scaling.material().density_kg_m3().to_bits()
            || original.axis_quantities() != material.axis_quantities()
            || !original
                .query_point()
                .iter()
                .filter(|(axis, _)| axis != law.temperature_axis())
                .eq(material
                    .query_point()
                    .iter()
                    .filter(|(axis, _)| axis != law.temperature_axis()))
        {
            return Err(RigidLiquidReleaseError::Invalid(
                "original matter must retain the receiver's density claim, context and uniform reference density",
            ));
        }
    }
    Ok(())
}

fn rigid_moments(
    frame: &PrincipalFrame,
    motion: RigidBodyState,
    grid_origin_m: [f64; 3],
) -> Result<LiquidTransferMoments, RigidLiquidReleaseError> {
    let properties = frame.properties();
    let center = motion.pose().position_world().sub(vector(grid_origin_m));
    let momentum = motion.linear_momentum_world();
    let spin = motion
        .pose()
        .orientation()
        .rotate_body_to_world(motion.angular_momentum_body());
    let moments = LiquidTransferMoments {
        mass_kg: properties.mass(),
        first_mass_moment_kg_m: array(center.scale(properties.mass())),
        momentum_kg_m_s: array(momentum),
        angular_momentum_kg_m2_s: array(center.cross(momentum).add(spin)),
        kinetic_energy_j: motion.kinetic_energy(properties)?,
    };
    if !finite_moments(moments) {
        return Err(RigidLiquidReleaseError::Invalid(
            "original rigid moments are unrepresentable",
        ));
    }
    Ok(moments)
}

fn rigid_storage(frame: &PrincipalFrame, motion: RigidBodyState) -> kinetic::RigidStorage {
    // PrincipalFrame lowers the original full tensor into these exact native
    // body axes. Do not derive a new inertia from the projected liquid grid.
    let properties = frame.properties();
    kinetic::RigidStorage {
        mass: properties.mass(),
        momentum_world: array(motion.linear_momentum_world()),
        principal_moments: array(properties.principal_inertia_body()),
        spin_principal: array(motion.angular_momentum_body()),
    }
}

struct Balance {
    centroid_m: [f64; 3],
    kinetic_j: Dd,
    gravity_j: f64,
    enthalpy_j: f64,
}

#[allow(clippy::too_many_arguments)]
fn audit(
    original: LiquidTransferMoments,
    original_h: f64,
    actual: LiquidTransferMoments,
    actual_h: f64,
    kinetic_change: Dd,
    gravity: [f64; 3],
    controls: RigidLiquidReleaseControls,
    stage: RigidLiquidBalanceStage,
) -> Result<Balance, RigidLiquidReleaseError> {
    let mass_error_kg = actual.mass_kg - original.mass_kg;
    let momentum_error_kg_m_s = difference(actual.momentum_kg_m_s, original.momentum_kg_m_s);
    let angular_momentum_error_kg_m2_s = difference(
        actual.angular_momentum_kg_m2_s,
        original.angular_momentum_kg_m2_s,
    );
    let centroid_error_m = core::array::from_fn(|axis| {
        actual.first_mass_moment_kg_m[axis] / actual.mass_kg
            - original.first_mass_moment_kg_m[axis] / original.mass_kg
    });
    let kinetic_error_j = kinetic_change.to_f64();
    let gravitational_error_j = (0..3)
        .fold(Dd::ZERO, |sum, axis| {
            sum - Dd::from_f64(gravity[axis])
                * (Dd::from_f64(actual.first_mass_moment_kg_m[axis])
                    - Dd::from_f64(original.first_mass_moment_kg_m[axis]))
        })
        .to_f64();
    let enthalpy_error_j = actual_h - original_h;
    let limits = controls.transfer.mechanical;
    let (centroid_limit, energy_ok) = match stage {
        RigidLiquidBalanceStage::RigidPartition
        | RigidLiquidBalanceStage::TetrahedralKinematics => (
            controls.source_centroid_tolerance_m,
            !Dd::from_f64(controls.source_kinetic_tolerance_j).lt(kinetic_change.abs()),
        ),
        RigidLiquidBalanceStage::NativeLiquid => (
            limits.maximum_centroid_shift_m,
            !Dd::from_f64(limits.maximum_projection_energy_change_j)
                .lt(kinetic_change.abs() + Dd::from_f64(gravitational_error_j.abs())),
        ),
    };
    if !finite_moments(actual)
        || !kinetic_change.lo.is_finite()
        || ![
            mass_error_kg,
            kinetic_error_j,
            gravitational_error_j,
            enthalpy_error_j,
        ]
        .iter()
        .all(|x| x.is_finite())
        || mass_error_kg.abs() > limits.mass_tolerance_kg
        || norm(momentum_error_kg_m_s) > limits.momentum_tolerance_kg_m_s
        || norm(angular_momentum_error_kg_m2_s) > limits.angular_momentum_tolerance_kg_m2_s
        || !norm(centroid_error_m).is_finite()
        || norm(centroid_error_m) > centroid_limit
        || enthalpy_error_j.abs() > controls.transfer.enthalpy_tolerance_j
        || !energy_ok
    {
        return Err(RigidLiquidReleaseError::Balance {
            stage,
            mass_error_kg,
            momentum_error_kg_m_s,
            angular_momentum_error_kg_m2_s,
            centroid_error_m,
            kinetic_error_j,
            gravitational_error_j,
            enthalpy_error_j,
        });
    }
    Ok(Balance {
        centroid_m: centroid_error_m,
        kinetic_j: kinetic_change,
        gravity_j: gravitational_error_j,
        enthalpy_j: enthalpy_error_j,
    })
}

fn reserved<T>(capacity: usize) -> Result<Vec<T>, RigidLiquidReleaseError> {
    let mut values = Vec::new();
    values.try_reserve_exact(capacity).map_err(|_| {
        RigidLiquidReleaseError::Invalid("rigid-to-P1 workspace allocation refused")
    })?;
    Ok(values)
}
fn check(cx: &Cx<'_>, gate: &CancelGate) -> Result<(), RigidLiquidReleaseError> {
    if gate.is_requested() || cx.checkpoint().is_err() {
        Err(RigidLiquidReleaseError::Cancelled)
    } else {
        Ok(())
    }
}
fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn finite_moments(m: LiquidTransferMoments) -> bool {
    m.mass_kg.is_finite()
        && m.mass_kg > 0.0
        && m.kinetic_energy_j.is_finite()
        && m.kinetic_energy_j >= 0.0
        && [
            m.first_mass_moment_kg_m,
            m.momentum_kg_m_s,
            m.angular_momentum_kg_m2_s,
        ]
        .iter()
        .flatten()
        .all(|x| x.is_finite())
}
