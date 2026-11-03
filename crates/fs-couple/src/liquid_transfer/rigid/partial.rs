//! Release selected original liquid cells and keep the remaining rigid matter.
//!
//! A conservative outside-grid test excludes retained material from the new
//! VOF domain. This is a simultaneous representation handoff, not a moving
//! cut-cell, fracture, melting, or future rigid/liquid interaction solver.

use super::{
    CancelGate, Cx, Dd, DynamicsError, EquilibriumEnthalpyPhaseCurve, FreeSurfaceEnthalpyError,
    LiquidTemperatureLaw, LiquidTransferMoments, MaterialFreeSurfaceDomain3,
    MaterialFreeSurfaceEnthalpy3, MaterialFreeSurfaceScaling, PrincipalFrame, PrincipalFrameError,
    RigidBodyState, RigidLiquidBalanceStage, RigidLiquidReleaseControls, RigidLiquidReleaseError,
    RigidLiquidReleaseReport, RigidTetrahedralLiquidState, TetrahedralDensityBody,
    admit_reference_density, array, audit, check, release_rigid_liquid, reserved, rigid_moments,
};
use fs_ga::Se3;
use fs_mbd::{Pose, UnitQuaternion};
use fs_mbd_material::MaterialInertiaError;
use fs_mbd_material::density_field::TetrahedralDensityPart;

/// One set of global transfer budgets, plus the new frame and separation gates.
#[derive(Clone, Copy, Debug)]
pub struct PartialRigidLiquidReleaseControls {
    pub release: RigidLiquidReleaseControls,
    /// Relative reconstruction tolerance for each actual full inertia, in
    /// (0, 1e-8]. Each piece uses the existing principal-frame owner.
    pub principal_relative_tolerance: f64,
    /// Expand the receiver's closed world box by this distance [m] before
    /// rejecting retained tetrahedra whose bounding boxes intersect it.
    pub minimum_retained_clearance_m: f64,
}

/// The returned native rigid body, with its exact original material receipts.
///
/// These are the surviving cells, not a nominal-density body with edited mass.
/// Original body coordinates, nodal h and accepted clock remain attached. Use
/// `frame.properties()` with the ordinary native rigid integrator; update
/// `motion` and `time_s` together when that owner accepts a subsequent step.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidLiquidRemainder {
    pub geometry: TetrahedralDensityBody,
    pub frame: PrincipalFrame,
    pub motion: RigidBodyState,
    pub specific_enthalpy_j_kg: Vec<f64>,
    pub time_s: f64,
}

impl RigidLiquidRemainder {
    /// Borrow the same accepted material state for another partial/full handoff.
    pub fn as_source(&self) -> RigidTetrahedralLiquidState<'_> {
        RigidTetrahedralLiquidState {
            geometry: &self.geometry,
            frame: &self.frame,
            motion: self.motion,
            specific_enthalpy_j_kg: &self.specific_enthalpy_j_kg,
            time_s: self.time_s,
        }
    }
}

/// Maps compact part indices to the input state's original geometry indices.
/// Compose these maps across successive releases to retain earlier ancestry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RigidLiquidPartLineage {
    pub original_vertex_indices: Vec<usize>,
    pub original_cell_indices: Vec<usize>,
}

/// Original rigid, split, and actual combined endpoint accounting about ONE
/// grid origin. The release receipt also retains native wet-cell lineage.
#[derive(Clone, Debug, PartialEq)]
pub struct PartialRigidLiquidReleaseReport {
    pub original: LiquidTransferMoments,
    pub original_enthalpy_j: f64,
    pub selected_lineage: RigidLiquidPartLineage,
    pub retained_lineage: RigidLiquidPartLineage,
    pub retained: LiquidTransferMoments,
    pub retained_enthalpy_j: f64,
    pub released: RigidLiquidReleaseReport,
    /// Both native child states minus the original body, before projection.
    pub rigid_partition_kinetic_error_j: f64,
    pub rigid_partition_centroid_error_m: [f64; 3],
    /// Actual native liquid plus the actual retained native rigid state.
    pub combined_target: LiquidTransferMoments,
    /// Numerical endpoint changes; neither is a release-work or heat source.
    pub kinetic_projection_change_j: f64,
    pub gravitational_projection_change_j: f64,
    pub enthalpy_error_j: f64,
    /// Whole-source phase checks plus the native selected-piece source queries,
    /// together bounded by the caller's one `max_source_queries` ceiling.
    pub source_queries: usize,
}

/// Every refusal leaves the original state available and returns neither child.
#[derive(Debug)]
pub enum PartialRigidLiquidReleaseError {
    Invalid(&'static str),
    Geometry(MaterialInertiaError),
    Release(RigidLiquidReleaseError),
    /// Conservative AABB overlap, in the original input cell numbering.
    RetainedOverlapsReceiver {
        original_cell: usize,
    },
}
impl core::fmt::Display for PartialRigidLiquidReleaseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "partial rigid liquid release: {self:?}")
    }
}
impl std::error::Error for PartialRigidLiquidReleaseError {}
impl From<MaterialInertiaError> for PartialRigidLiquidReleaseError {
    fn from(error: MaterialInertiaError) -> Self {
        Self::Geometry(error)
    }
}
impl From<RigidLiquidReleaseError> for PartialRigidLiquidReleaseError {
    fn from(error: RigidLiquidReleaseError) -> Self {
        Self::Release(error)
    }
}
impl From<PrincipalFrameError> for PartialRigidLiquidReleaseError {
    fn from(error: PrincipalFrameError) -> Self {
        Self::Release(error.into())
    }
}
impl From<DynamicsError> for PartialRigidLiquidReleaseError {
    fn from(error: DynamicsError) -> Self {
        Self::Release(error.into())
    }
}

/// Replace one no-storage rigid liquid constraint with a selected thermal VOF
/// liquid and the still-constrained original remainder, without release impulse.
///
/// `selected_cells` addresses `source.geometry.tetrahedra()`; both outputs must
/// be nonempty. Input order is immaterial and duplicates refuse. The geometry
/// owner partitions ORIGINAL cell mass, density receipts, full tensors and
/// topology. Child centre velocities and angular velocity come from the actual
/// parent's native rigid velocity field. No nominal tensor, velocity reset,
/// recoilless mass removal or inferred separation work substitutes for it.
///
/// The WHOLE source must already be liquid on the same constant reference
/// density caloric chart, with no stored elastic/plastic/cohesive history.
/// The complete nodal h span, including interior caloric knots, is admitted.
/// The selected part then uses the existing source-backed thermal/rheological/
/// Boussinesq transfer. The remainder keeps its original h values and clock.
/// This does not model a stressed solid losing a melted or fractured fragment.
///
/// Every retained tetrahedron's world AABB must lie strictly outside the entire
/// receiver grid plus caller clearance. Touching or ambiguous overlap refuses,
/// even in a dry/solid voxel: this API provides no moving occupancy closure.
/// The condition holds at handoff only. Subsequent interaction or grid entry
/// requires an admitted moving boundary/coupling owner before either advances.
///
/// Original M/P/L/h and zero-work rigid kinetic identity are checked before
/// projection. The ACTUAL liquid plus retained native body are checked again
/// against the original, with the same global budgets. Numerical kinetic and
/// gravitational changes stay separate from h. No time is advanced. On success
/// replace the original with BOTH returned states; also evolving the original
/// would duplicate matter. Any refusal leaves it unchanged and publishes none.
#[allow(clippy::too_many_arguments)]
pub fn release_rigid_liquid_cells<'s, 'm>(
    source: RigidTetrahedralLiquidState<'_>,
    selected_cells: &[usize],
    scaling: MaterialFreeSurfaceScaling<'s>,
    domain: MaterialFreeSurfaceDomain3<'_>,
    grid_origin_m: [f64; 3],
    law: &'m LiquidTemperatureLaw<'m>,
    caloric: &'m EquilibriumEnthalpyPhaseCurve,
    controls: PartialRigidLiquidReleaseControls,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<
    (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        RigidLiquidRemainder,
        PartialRigidLiquidReleaseReport,
    ),
    PartialRigidLiquidReleaseError,
> {
    check(cx, gate)?;
    let release = controls.release;
    let maximum_vertices = release
        .transfer
        .mechanical
        .max_tetrahedra
        .checked_mul(4)
        .ok_or(PartialRigidLiquidReleaseError::Invalid(
            "source vertex ceiling overflows",
        ))?;
    if source.geometry.vertices().len() > maximum_vertices
        || source.specific_enthalpy_j_kg.len() != source.geometry.vertices().len()
        || source.frame.inertia() != source.geometry.inertia()
        || !source.time_s.is_finite()
        || source.time_s < 0.0
        || !controls.principal_relative_tolerance.is_finite()
        || controls.principal_relative_tolerance <= 0.0
        || controls.principal_relative_tolerance > 1e-8
        || !controls.minimum_retained_clearance_m.is_finite()
        || controls.minimum_retained_clearance_m < 0.0
        || grid_origin_m.iter().any(|x| !x.is_finite())
        || domain.gravity_m_s2.iter().any(|x| !x.is_finite())
        || [
            release.source_kinetic_tolerance_j,
            release.source_centroid_tolerance_m,
            release.transfer.mechanical.mass_tolerance_kg,
            release.transfer.mechanical.momentum_tolerance_kg_m_s,
            release
                .transfer
                .mechanical
                .angular_momentum_tolerance_kg_m2_s,
            release.transfer.mechanical.maximum_centroid_shift_m,
            release
                .transfer
                .mechanical
                .maximum_projection_energy_change_j,
            release.transfer.enthalpy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(PartialRigidLiquidReleaseError::Invalid(
            "exact source frame, bounded h, clock, grid, clearance and numerical tolerances required",
        ));
    }
    let split = source.geometry.partition_cells(
        selected_cells,
        release.transfer.mechanical.max_tetrahedra,
        maximum_vertices,
        || gate.is_requested() || cx.checkpoint().is_err(),
    )?;
    admit_reference_density(source.geometry, scaling, law, cx, gate)?;
    let phase_queries = admit_liquid_source(source, scaling, law, caloric, release, cx, gate)?;
    let original = rigid_moments(source.frame, source.motion, grid_origin_m)?;
    let original_h = integrated_enthalpy(source, cx, gate)?;
    let world_from_body = source.frame.world_from_body(source.motion)?;
    let (selected, selected_lineage) = prepare_part(
        split.selected,
        source,
        world_from_body,
        controls.principal_relative_tolerance,
        cx,
        gate,
    )?;
    let (retained, retained_lineage) = prepare_part(
        split.retained,
        source,
        world_from_body,
        controls.principal_relative_tolerance,
        cx,
        gate,
    )?;
    admit_outside_receiver(
        &retained,
        &retained_lineage,
        scaling.cell_size_m(),
        domain.cells,
        grid_origin_m,
        controls.minimum_retained_clearance_m,
        cx,
        gate,
    )?;
    let selected_moments = rigid_moments(&selected.frame, selected.motion, grid_origin_m)?;
    let retained_moments = rigid_moments(&retained.frame, retained.motion, grid_origin_m)?;
    let selected_h = integrated_enthalpy(selected.as_source(), cx, gate)?;
    let retained_h = integrated_enthalpy(retained.as_source(), cx, gate)?;
    let partition = audit(
        original,
        original_h,
        add_moments(selected_moments, retained_moments),
        add(selected_h, retained_h),
        domain.gravity_m_s2,
        release,
        RigidLiquidBalanceStage::RigidPartition,
    )?;
    // Count the extra whole-source admission against the SAME work budget.
    let mut selected_controls = release;
    selected_controls.transfer.max_source_queries = release
        .transfer
        .max_source_queries
        .checked_sub(phase_queries)
        .filter(|remaining| *remaining > 0)
        .ok_or(PartialRigidLiquidReleaseError::Invalid(
            "source-query budget exhausted",
        ))?;
    let (liquid, released) = release_rigid_liquid(
        selected.as_source(),
        scaling,
        domain,
        grid_origin_m,
        law,
        caloric,
        selected_controls,
        cx,
        gate,
    )?;
    let combined_target = add_moments(released.transfer.mechanical.target, retained_moments);
    let endpoint = audit(
        original,
        original_h,
        combined_target,
        add(released.transfer.target_enthalpy_j, retained_h),
        domain.gravity_m_s2,
        release,
        RigidLiquidBalanceStage::NativeLiquid,
    )?;
    let source_queries = phase_queries
        .checked_add(released.transfer.source_queries)
        .filter(|count| *count <= release.transfer.max_source_queries)
        .ok_or(PartialRigidLiquidReleaseError::Invalid(
            "source-query accounting overflow",
        ))?;
    check(cx, gate)?;
    Ok((
        liquid,
        retained,
        PartialRigidLiquidReleaseReport {
            original,
            original_enthalpy_j: original_h,
            selected_lineage,
            retained_lineage,
            retained: retained_moments,
            retained_enthalpy_j: retained_h,
            released,
            rigid_partition_kinetic_error_j: partition.kinetic_j,
            rigid_partition_centroid_error_m: partition.centroid_m,
            combined_target,
            kinetic_projection_change_j: endpoint.kinetic_j,
            gravitational_projection_change_j: endpoint.gravity_j,
            enthalpy_error_j: endpoint.enthalpy_j,
            source_queries,
        },
    ))
}

fn prepare_part(
    part: TetrahedralDensityPart,
    source: RigidTetrahedralLiquidState<'_>,
    world_from_body: Se3,
    tolerance: f64,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<(RigidLiquidRemainder, RigidLiquidPartLineage), PartialRigidLiquidReleaseError> {
    let frame = part.body.principal_frame(tolerance, || {
        gate.is_requested() || cx.checkpoint().is_err()
    })?;
    // Sample at the child's actual COM. Reconstructing a twist at a distant
    // authored origin would introduce avoidable large-arm cancellation.
    let arm = source
        .frame
        .arm_for_body_point(frame.inertia().center_of_mass())?;
    let center = source
        .motion
        .point_kinematics(source.frame.properties(), arm)?;
    let placement = world_from_body
        .compose(frame.body_from_principal())
        .map_err(PrincipalFrameError::from)?;
    let q = *placement.rotation().as_quat();
    let orientation = UnitQuaternion::new(q.w, q.x, q.y, q.z)?;
    let angular_velocity = orientation.rotate_world_to_body(center.angular_velocity_world);
    let motion = RigidBodyState::new(
        Pose::new(center.point_world, orientation)?,
        center.point_velocity_world.scale(frame.properties().mass()),
        frame
            .properties()
            .angular_momentum_body_checked(angular_velocity)?,
    )?;
    let mut h = reserved(part.original_vertex_indices.len())?;
    for &original in &part.original_vertex_indices {
        check(cx, gate)?;
        h.push(source.specific_enthalpy_j_kg[original]);
    }
    Ok((
        RigidLiquidRemainder {
            geometry: part.body,
            frame,
            motion,
            specific_enthalpy_j_kg: h,
            time_s: source.time_s,
        },
        RigidLiquidPartLineage {
            original_vertex_indices: part.original_vertex_indices,
            original_cell_indices: part.original_cell_indices,
        },
    ))
}

#[allow(clippy::too_many_arguments)]
fn admit_outside_receiver(
    retained: &RigidLiquidRemainder,
    lineage: &RigidLiquidPartLineage,
    dx: f64,
    cells: [usize; 3],
    origin: [f64; 3],
    clearance: f64,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<(), PartialRigidLiquidReleaseError> {
    let low: [f64; 3] = core::array::from_fn(|a| origin[a] - clearance);
    let high: [f64; 3] = core::array::from_fn(|a| origin[a] + cells[a] as f64 * dx + clearance);
    if cells.contains(&0)
        || !dx.is_finite()
        || dx <= 0.0
        || (0..3).any(|a| !low[a].is_finite() || !high[a].is_finite() || high[a] <= low[a])
    {
        return Err(PartialRigidLiquidReleaseError::Invalid(
            "receiver world box is unrepresentable",
        ));
    }
    let mut positions = reserved(retained.geometry.vertices().len())?;
    for vertex in retained.geometry.vertices() {
        check(cx, gate)?;
        let arm = retained.frame.arm_for_body_point(vertex.position_m())?;
        positions.push(array(
            retained
                .motion
                .point_kinematics(retained.frame.properties(), arm)?
                .point_world,
        ));
    }
    for (tet, &original_cell) in retained
        .geometry
        .tetrahedra()
        .iter()
        .zip(&lineage.original_cell_indices)
    {
        check(cx, gate)?;
        let separated = (0..3).any(|a| {
            tet.iter().all(|&v| positions[v][a] < low[a])
                || tet.iter().all(|&v| positions[v][a] > high[a])
        });
        if !separated {
            return Err(PartialRigidLiquidReleaseError::RetainedOverlapsReceiver { original_cell });
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn admit_liquid_source(
    source: RigidTetrahedralLiquidState<'_>,
    scaling: MaterialFreeSurfaceScaling<'_>,
    law: &LiquidTemperatureLaw<'_>,
    caloric: &EquilibriumEnthalpyPhaseCurve,
    controls: RigidLiquidReleaseControls,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<usize, PartialRigidLiquidReleaseError> {
    let maximum = controls.transfer.max_source_queries;
    if caloric.knots().len() > maximum
        || source.specific_enthalpy_j_kg.len() > maximum
        || caloric.material_card_identity() != scaling.material().material().card_identity()
        || law.card().content_hash() != caloric.material_card_identity()
    {
        return Err(PartialRigidLiquidReleaseError::Invalid(
            "bounded original caloric source required",
        ));
    }
    let density = scaling.material().density_kg_m3();
    let (mut low, mut high) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut queries = 0usize;
    let mut admit = |h: f64| -> Result<(), PartialRigidLiquidReleaseError> {
        check(cx, gate)?;
        queries = queries.checked_add(1).filter(|n| *n <= maximum).ok_or(
            PartialRigidLiquidReleaseError::Invalid("source-query budget exhausted"),
        )?;
        // Admission needs the original chart's numeric state, not a freshly
        // hashed phase identity for every carried vertex.
        let phase = caloric
            .values_at_specific_enthalpy(h)
            .map_err(FreeSurfaceEnthalpyError::from)
            .map_err(RigidLiquidReleaseError::from)?;
        if phase.liquid_mass_fraction != 1.0
            || phase.bulk_density_kg_m3.to_bits() != density.to_bits()
        {
            return Err(PartialRigidLiquidReleaseError::Invalid(
                "whole original source must already be liquid at its reference density",
            ));
        }
        Ok(())
    };
    for &h in source.specific_enthalpy_j_kg {
        admit(h)?;
        low = low.min(h);
        high = high.max(h);
    }
    // P1 h is a convex combination. Admitting each interior chart knot also
    // excludes an unsupported phase/density interval hidden between vertices.
    for knot in caloric.knots() {
        check(cx, gate)?;
        if knot.specific_enthalpy_j_kg >= low && knot.specific_enthalpy_j_kg <= high {
            admit(knot.specific_enthalpy_j_kg)?;
        }
    }
    Ok(queries)
}

fn integrated_enthalpy(
    source: RigidTetrahedralLiquidState<'_>,
    cx: &Cx<'_>,
    gate: &CancelGate,
) -> Result<f64, PartialRigidLiquidReleaseError> {
    let mut sum = Dd::ZERO;
    for (tet, cell) in source
        .geometry
        .tetrahedra()
        .iter()
        .zip(source.geometry.cells())
    {
        check(cx, gate)?;
        let mean = tet.iter().fold(Dd::ZERO, |value, &v| {
            value + Dd::from_f64(source.specific_enthalpy_j_kg[v]) * Dd::from_f64(0.25)
        });
        sum = sum + Dd::from_f64(cell.inertia().mass()) * mean;
    }
    let h = sum.to_f64();
    if !h.is_finite() {
        return Err(PartialRigidLiquidReleaseError::Invalid(
            "original extensive enthalpy is unrepresentable",
        ));
    }
    Ok(h)
}

fn add(a: f64, b: f64) -> f64 {
    (Dd::from_f64(a) + Dd::from_f64(b)).to_f64()
}
fn add_moments(a: LiquidTransferMoments, b: LiquidTransferMoments) -> LiquidTransferMoments {
    LiquidTransferMoments {
        mass_kg: add(a.mass_kg, b.mass_kg),
        first_mass_moment_kg_m: core::array::from_fn(|i| {
            add(a.first_mass_moment_kg_m[i], b.first_mass_moment_kg_m[i])
        }),
        momentum_kg_m_s: core::array::from_fn(|i| add(a.momentum_kg_m_s[i], b.momentum_kg_m_s[i])),
        angular_momentum_kg_m2_s: core::array::from_fn(|i| {
            add(a.angular_momentum_kg_m2_s[i], b.angular_momentum_kg_m2_s[i])
        }),
        kinetic_energy_j: add(a.kinetic_energy_j, b.kinetic_energy_j),
    }
}
