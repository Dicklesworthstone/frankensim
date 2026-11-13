//! Reciprocal viscous momentum exchange between original liquid and rigid owners.
//!
//! The uniform-rate SOURCE is exact at fixed cell positions and rigid pose.
//! Each cell's actual retained liquid mass recoils against the original finite
//! body through one six-coordinate kinetic-metric exponential. The heterogeneous
//! source adds directional per-cell resistance through bounded symmetric
//! composition of exact local exchanges, with explicit temporal refinement.
//!
//! This operation advances neither owner clock, geometry, transport, pressure,
//! collision, free rotation nor thermodynamic state. Compose it with the native
//! drift/flow integrators under an explicitly chosen splitting. The model is
//! volumetric slip resistance on selected wet cells, not an impermeable wall,
//! excluded volume, inferred hydrodynamic added mass, or contact geometry.

/// Spatially varying directional drag with reciprocal finite-body recoil.
pub mod heterogeneous;
/// Same native source on the original material-bound liquid and its exact units.
pub mod material;
mod relaxation;

use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::{FreeSurface3, FreeSurfaceImpulseStencil3};
use fs_lbm::impulse::{ImpulseControls, ImpulseError, ImpulseReport};
use fs_math::dd::Dd;
use fs_mbd::{DynamicsError, MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3};

/// Original native/SI correspondence. Grid axes coincide with inertial axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidMomentumUnits3 {
    pub cell_size_m: f64,
    pub velocity_unit_m_s: f64,
    pub mass_unit_kg: f64,
    /// World location of the grid's lower corner, independent of the body COM.
    pub grid_origin_world_m: Vec3,
}

/// Physical and allocation gates for the complete reciprocal source operation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LiquidRigidDragControls {
    /// Original full domain cap, including inactive cells.
    pub maximum_cells: usize,
    /// Compact selected-cell capacity, checked before any stencil allocation.
    pub maximum_stencil_cells: usize,
    /// Both fluid and coupled rigid-point speeds, before and after, in m/s.
    pub maximum_speed_m_s: f64,
    pub maximum_momentum_defect_kg_m_s: f64,
    pub maximum_angular_momentum_defect_kg_m2_s: f64,
    pub maximum_work_defect_j: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum LiquidMomentumError {
    Invalid(&'static str),
    Fluid(ImpulseError),
    Mechanics(DynamicsError),
    Balance {
        account: &'static str,
        residual: f64,
        maximum: f64,
    },
    Cancelled,
}
impl core::fmt::Display for LiquidMomentumError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reciprocal liquid/rigid drag: {self:?}")
    }
}
impl std::error::Error for LiquidMomentumError {}
impl From<ImpulseError> for LiquidMomentumError {
    fn from(error: ImpulseError) -> Self {
        Self::Fluid(error)
    }
}
impl From<DynamicsError> for LiquidMomentumError {
    fn from(error: DynamicsError) -> Self {
        Self::Mechanics(error)
    }
}

/// Both owners' actual accepted momentum and kinetic accounts for one source.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct LiquidRigidDragReport {
    pub source_duration_s: f64,
    pub relaxation_rate_per_s: f64,
    pub selected_cells: usize,
    /// Original population receipt in NATIVE lattice units; units are retained.
    pub native_liquid: ImpulseReport<3>,
    pub units: LiquidMomentumUnits3,
    /// Original SI momentum owner, with the fixed before/after pose.
    pub body: MomentumImpulseReceipt,
    /// Common fixed COM origin for every angular account below.
    pub balance_origin_world_m: Vec3,
    pub liquid_impulse_kg_m_s: Vec3,
    pub liquid_angular_impulse_kg_m2_s: Vec3,
    pub liquid_momentum_change_kg_m_s: Vec3,
    pub liquid_angular_momentum_change_kg_m2_s: Vec3,
    pub liquid_kinetic_energy_change_j: f64,
    /// Physical kinetic loss of this declared drag source; no heat is deposited.
    pub drag_dissipation_j: f64,
    pub linear_residual_kg_m_s: Vec3,
    pub angular_residual_kg_m2_s: Vec3,
    pub kinetic_balance_residual_j: f64,
    /// Local readback, global balance and source-response defects cannot cancel.
    pub absolute_momentum_defect_kg_m_s: f64,
    pub absolute_angular_momentum_defect_kg_m2_s: f64,
    pub absolute_work_defect_j: f64,
}

/// Actual native momentum and kinetic accounts shared by reciprocal drag laws.
/// Source-law parameters and composition controls remain in their own report.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct LiquidRigidExchangeReport {
    pub selected_cells: usize,
    /// Original population receipt in NATIVE lattice units; units are retained.
    pub native_liquid: ImpulseReport<3>,
    pub units: LiquidMomentumUnits3,
    /// Original SI momentum owner, with the fixed before/after pose.
    pub body: MomentumImpulseReceipt,
    /// Common fixed COM origin for every angular account below.
    pub balance_origin_world_m: Vec3,
    pub liquid_impulse_kg_m_s: Vec3,
    pub liquid_angular_impulse_kg_m2_s: Vec3,
    pub liquid_momentum_change_kg_m_s: Vec3,
    pub liquid_angular_momentum_change_kg_m2_s: Vec3,
    pub liquid_kinetic_energy_change_j: f64,
    /// Physical kinetic loss of this declared drag source; no heat is deposited.
    pub drag_dissipation_j: f64,
    pub linear_residual_kg_m_s: Vec3,
    pub angular_residual_kg_m2_s: Vec3,
    pub kinetic_balance_residual_j: f64,
    /// Local readback, global balance and source-response defects cannot cancel.
    pub absolute_momentum_defect_kg_m_s: f64,
    pub absolute_angular_momentum_defect_kg_m2_s: f64,
    pub absolute_work_defect_j: f64,
}

impl LiquidRigidExchangeReport {
    fn with_uniform_rate(
        self,
        source_duration_s: f64,
        relaxation_rate_per_s: f64,
    ) -> LiquidRigidDragReport {
        LiquidRigidDragReport {
            source_duration_s,
            relaxation_rate_per_s,
            selected_cells: self.selected_cells,
            native_liquid: self.native_liquid,
            units: self.units,
            body: self.body,
            balance_origin_world_m: self.balance_origin_world_m,
            liquid_impulse_kg_m_s: self.liquid_impulse_kg_m_s,
            liquid_angular_impulse_kg_m2_s: self.liquid_angular_impulse_kg_m2_s,
            liquid_momentum_change_kg_m_s: self.liquid_momentum_change_kg_m_s,
            liquid_angular_momentum_change_kg_m2_s: self.liquid_angular_momentum_change_kg_m2_s,
            liquid_kinetic_energy_change_j: self.liquid_kinetic_energy_change_j,
            drag_dissipation_j: self.drag_dissipation_j,
            linear_residual_kg_m_s: self.linear_residual_kg_m_s,
            angular_residual_kg_m2_s: self.angular_residual_kg_m2_s,
            kinetic_balance_residual_j: self.kinetic_balance_residual_j,
            absolute_momentum_defect_kg_m_s: self.absolute_momentum_defect_kg_m_s,
            absolute_angular_momentum_defect_kg_m2_s: self.absolute_angular_momentum_defect_kg_m2_s,
            absolute_work_defect_j: self.absolute_work_defect_j,
        }
    }
}

#[derive(Clone, Copy)]
struct Units {
    impulse: f64,
    energy: f64,
}

fn dd(x: f64) -> Dd {
    Dd::from_f64(x)
}
fn vector(v: [Dd; 3]) -> Vec3 {
    Vec3::new(v[0].to_f64(), v[1].to_f64(), v[2].to_f64())
}
fn a(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn v(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn dot(a: Vec3, b: Vec3) -> Dd {
    dd(a.x) * dd(b.x) + dd(a.y) * dd(b.y) + dd(a.z) * dd(b.z)
}
fn cross(a: Vec3, b: Vec3) -> [Dd; 3] {
    [
        dd(a.y) * dd(b.z) - dd(a.z) * dd(b.y),
        dd(a.z) * dd(b.x) - dd(a.x) * dd(b.z),
        dd(a.x) * dd(b.y) - dd(a.y) * dd(b.x),
    ]
}
fn difference(a: Vec3, b: Vec3) -> Vec3 {
    vector([dd(a.x) - dd(b.x), dd(a.y) - dd(b.y), dd(a.z) - dd(b.z)])
}
fn l1(v: Vec3) -> f64 {
    v.x.abs() + v.y.abs() + v.z.abs()
}
fn invalid(reason: &'static str) -> LiquidMomentumError {
    LiquidMomentumError::Invalid(reason)
}
fn finite_vec(v: Vec3) -> Result<(), LiquidMomentumError> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(invalid("nonfinite liquid/rigid vector"))
    }
}
fn positive(v: f64) -> bool {
    v.is_finite() && v > 0.0
}
fn finite(v: f64) -> Result<f64, LiquidMomentumError> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(invalid("nonfinite liquid/rigid account"))
    }
}
fn poll(gate: &CancelGate) -> Result<(), LiquidMomentumError> {
    if gate.is_requested() {
        Err(LiquidMomentumError::Cancelled)
    } else {
        Ok(())
    }
}
fn check(account: &'static str, residual: f64, maximum: f64) -> Result<(), LiquidMomentumError> {
    if !residual.is_finite() || residual > maximum {
        Err(LiquidMomentumError::Balance {
            account,
            residual,
            maximum,
        })
    } else {
        Ok(())
    }
}
fn represented(input: f64, output: f64) -> Result<f64, LiquidMomentumError> {
    if !output.is_finite() || (input != 0.0 && output == 0.0) {
        Err(invalid("liquid/rigid unit conversion is unrepresentable"))
    } else {
        Ok(output)
    }
}

impl LiquidMomentumUnits3 {
    fn admit(self) -> Result<Units, LiquidMomentumError> {
        if !self.grid_origin_world_m.is_finite()
            || ![self.cell_size_m, self.velocity_unit_m_s, self.mass_unit_kg]
                .into_iter()
                .all(positive)
        {
            return Err(invalid(
                "finite positive original liquid scales and a finite origin are required",
            ));
        }
        let impulse = self.mass_unit_kg * self.velocity_unit_m_s;
        let energy = impulse * self.velocity_unit_m_s;
        if !positive(impulse) || !positive(energy) {
            return Err(invalid("liquid momentum/work units are unrepresentable"));
        }
        Ok(Units { impulse, energy })
    }
}
impl LiquidRigidDragControls {
    fn admit(self, len: usize) -> Result<(), LiquidMomentumError> {
        if len == 0
            || len > self.maximum_stencil_cells
            || self.maximum_stencil_cells > self.maximum_cells
            || !positive(self.maximum_speed_m_s)
            || [
                self.maximum_momentum_defect_kg_m_s,
                self.maximum_angular_momentum_defect_kg_m2_s,
                self.maximum_work_defect_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "bounded nonempty liquid stencil, speed and nonnegative physical tolerances required",
            ));
        }
        Ok(())
    }
    fn native(
        self,
        scale: LiquidMomentumUnits3,
        units: Units,
    ) -> Result<ImpulseControls, LiquidMomentumError> {
        Ok(ImpulseControls {
            max_cells: self.maximum_cells,
            maximum_speed: represented(
                self.maximum_speed_m_s,
                self.maximum_speed_m_s / scale.velocity_unit_m_s,
            )?,
            maximum_momentum_defect: represented(
                self.maximum_momentum_defect_kg_m_s,
                self.maximum_momentum_defect_kg_m_s / units.impulse,
            )?,
            maximum_work_defect: represented(
                self.maximum_work_defect_j,
                self.maximum_work_defect_j / units.energy,
            )?,
        })
    }
}
fn speed(v: Vec3, maximum: f64) -> Result<(), LiquidMomentumError> {
    let speed = v.x.hypot(v.y).hypot(v.z);
    if !speed.is_finite() || speed > maximum {
        Err(invalid(
            "liquid or coupled rigid point exceeds the physical speed bound",
        ))
    } else {
        Ok(())
    }
}

/// Apply exact reciprocal slip relaxation to an original native wet stencil.
///
/// The source law is du_i/dt = gamma (V + omega cross r_i - u_i), with each
/// cell's own retained liquid mass and equal/opposite finite-body force/torque.
/// Gamma is a homogeneous authored inverse time; it is not viscosity or an
/// inferred wall condition. A six-coordinate exponential resolves arbitrarily
/// stiff source relaxation without an iteration or stiff-step slip reversal.
/// Cell positions, body pose and both clocks stay fixed during this SOURCE.
/// Native pressure/transport/rotation evolution remains with its original owner.
///
/// Both states remain unchanged on any refusal or cancellation. A successful
/// source publishes only selected populations and the matching native momenta.
/// No whole-grid clone, dense impulse field, or per-cell factorization occurs.
#[allow(clippy::too_many_arguments)]
pub fn apply_rigid_drag_native(
    liquid: &mut FreeSurface3,
    body: &mut RigidBodyState,
    properties: MassProperties,
    indices: &[usize],
    units: LiquidMomentumUnits3,
    relaxation_rate_per_s: f64,
    source_duration_s: f64,
    controls: LiquidRigidDragControls,
    gate: &CancelGate,
) -> Result<LiquidRigidDragReport, LiquidMomentumError> {
    poll(gate)?;
    controls.admit(indices.len())?;
    let scales = units.admit()?;
    let stage = liquid.prepare_impulse_stencil(indices, controls.native(units, scales)?, gate)?;
    apply_staged(
        stage,
        body,
        properties,
        units,
        relaxation_rate_per_s,
        source_duration_s,
        controls,
        gate,
    )
}

#[allow(clippy::too_many_arguments)]
fn apply_staged(
    stage: FreeSurfaceImpulseStencil3<'_>,
    body: &mut RigidBodyState,
    properties: MassProperties,
    units: LiquidMomentumUnits3,
    relaxation_rate_per_s: f64,
    source_duration_s: f64,
    controls: LiquidRigidDragControls,
    gate: &CancelGate,
) -> Result<LiquidRigidDragReport, LiquidMomentumError> {
    poll(gate)?;
    if !relaxation_rate_per_s.is_finite()
        || relaxation_rate_per_s < 0.0
        || !source_duration_s.is_finite()
        || source_duration_s < 0.0
    {
        return Err(invalid(
            "finite nonnegative drag rate and source duration required",
        ));
    }
    let (cells, basis) = prepare_source(&stage, body, properties, units, controls, gate)?;
    // Positive overflow of gamma*duration is the resolved exponential-zero
    // limit. Both authored factors are finite; no infinity enters the states.
    let exposure = relaxation_rate_per_s * source_duration_s;
    let plan = relaxation::solve(&cells, basis, exposure, gate)?;
    let report = finish_source(
        stage, body, properties, units, controls, &cells, basis, plan, None, gate,
    )?;
    Ok(report.with_uniform_rate(source_duration_s, relaxation_rate_per_s))
}

fn prepare_source(
    stage: &FreeSurfaceImpulseStencil3<'_>,
    body: &RigidBodyState,
    properties: MassProperties,
    units: LiquidMomentumUnits3,
    controls: LiquidRigidDragControls,
    gate: &CancelGate,
) -> Result<(Vec<relaxation::Cell>, relaxation::Body), LiquidMomentumError> {
    let count = stage.cells().len();
    controls.admit(count)?;
    body.kinetic_energy(properties)?;
    let before = *body;
    let center = before.pose().position_world();
    let orientation = before.pose().orientation();
    let principal = properties.principal_axes_to_body();
    let moments = a(properties.principal_inertia_body());
    let basis = relaxation::Body {
        velocity: before.center_of_mass_velocity_world(properties)?,
        omega: orientation.rotate_body_to_world(
            properties.angular_velocity_body_checked(before.angular_momentum_body())?,
        ),
        inverse_sqrt_mass: 1.0 / properties.mass().sqrt(),
        inverse_sqrt_inertia: moments.map(|i| 1.0 / i.sqrt()),
        principal_axes: [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
        .map(|axis| orientation.rotate_body_to_world(principal.rotate_body_to_world(axis))),
    };
    speed(basis.velocity, controls.maximum_speed_m_s)?;
    let mut cells = Vec::new();
    cells
        .try_reserve_exact(count)
        .map_err(|_| invalid("physical wet stencil allocation refused"))?;
    for (index, sample) in stage.cells().enumerate() {
        if index % 64 == 0 {
            poll(gate)?;
        }
        let mass = represented(sample.mass, sample.mass * units.mass_unit_kg)?;
        let arm = vector(core::array::from_fn(|axis| {
            dd(a(units.grid_origin_world_m)[axis]) - dd(a(center)[axis])
                + dd(sample.position[axis]) * dd(units.cell_size_m)
        }));
        finite_vec(arm)?;
        let mut velocity_si = [0.0; 3];
        for axis in 0..3 {
            velocity_si[axis] = represented(
                sample.velocity[axis],
                sample.velocity[axis] * units.velocity_unit_m_s,
            )?;
        }
        let velocity = v(velocity_si);
        speed(velocity, controls.maximum_speed_m_s)?;
        speed(basis.point_velocity(arm), controls.maximum_speed_m_s)?;
        cells.push(relaxation::Cell {
            mass,
            arm,
            velocity,
        });
    }
    Ok((cells, basis))
}

#[allow(clippy::too_many_arguments)]
fn finish_source(
    stage: FreeSurfaceImpulseStencil3<'_>,
    body: &mut RigidBodyState,
    properties: MassProperties,
    units: LiquidMomentumUnits3,
    controls: LiquidRigidDragControls,
    cells: &[relaxation::Cell],
    basis: relaxation::Body,
    plan: relaxation::Plan,
    physical_dissipation_j: Option<f64>,
    gate: &CancelGate,
) -> Result<LiquidRigidExchangeReport, LiquidMomentumError> {
    let scales = units.admit()?;
    let count = cells.len();
    let before = *body;
    let center = before.pose().position_world();
    let orientation = before.pose().orientation();
    let mut applied = [Dd::ZERO; 3];
    let mut moment = [Dd::ZERO; 3];
    for (cell, &j) in cells.iter().zip(&plan.impulses) {
        let angular = cross(cell.arm, j);
        for axis in 0..3 {
            applied[axis] = applied[axis] + dd(a(j)[axis]);
            moment[axis] = moment[axis] + angular[axis];
        }
    }
    let applied = vector(applied);
    let moment = vector(moment);
    finite_vec(applied)?;
    finite_vec(moment)?;
    let source_linear = l1(applied.add(plan.body_linear));
    let source_angular = l1(moment.add(plan.body_angular));
    let mechanics = before.apply_momentum_impulse(
        properties,
        applied.scale(-1.0),
        orientation.rotate_world_to_body(moment.scale(-1.0)),
    )?;
    let after = mechanics.state_after;
    let velocity_after = after.center_of_mass_velocity_world(properties)?;
    let omega_after = orientation.rotate_body_to_world(
        properties.angular_velocity_body_checked(after.angular_momentum_body())?,
    );
    speed(velocity_after, controls.maximum_speed_m_s)?;
    let dv = difference(velocity_after, basis.velocity);
    let dw = difference(omega_after, basis.omega);
    let mut loss = Dd::ZERO;
    let mut native_impulses = Vec::new();
    native_impulses
        .try_reserve_exact(count)
        .map_err(|_| invalid("native compact impulse conversion allocation refused"))?;
    for (index, (cell, &j)) in cells.iter().zip(&plan.impulses).enumerate() {
        if index % 64 == 0 {
            poll(gate)?;
        }
        speed(
            velocity_after.add(omega_after.cross(cell.arm)),
            controls.maximum_speed_m_s,
        )?;
        let initial_slip = difference(cell.velocity, basis.point_velocity(cell.arm));
        let body_change = dv.add(dw.cross(cell.arm));
        let mean_relative = vector(core::array::from_fn(|axis| {
            dd(a(initial_slip)[axis])
                + dd(0.5) * (dd(a(j)[axis]) / dd(cell.mass) - dd(a(body_change)[axis]))
        }));
        loss = loss - dot(j, mean_relative);
        let mut converted = [0.0; 3];
        for axis in 0..3 {
            converted[axis] = represented(a(j)[axis], a(j)[axis] / scales.impulse)?;
        }
        native_impulses.push(converted);
    }
    let endpoint_loss = finite(loss.to_f64())?;
    if endpoint_loss < 0.0 {
        return Err(invalid(
            "liquid drag candidate increases relative kinetic storage",
        ));
    }
    let (loss, source_work_defect) = if let Some(integrated) = physical_dissipation_j {
        if !integrated.is_finite() || integrated < 0.0 {
            return Err(invalid(
                "local drag dissipation must be finite and nonnegative",
            ));
        }
        (
            integrated,
            (dd(integrated) - dd(endpoint_loss)).to_f64().abs(),
        )
    } else {
        (endpoint_loss, 0.0)
    };
    let candidate = stage.with_impulses(&native_impulses, gate)?;
    let mut change = [Dd::ZERO; 3];
    let mut angular_change = [Dd::ZERO; 3];
    let mut local_linear = Dd::ZERO;
    let mut local_angular = Dd::ZERO;
    let mut local_work = Dd::ZERO;
    for (index, ((read, cell), &command)) in candidate
        .readback()
        .zip(cells)
        .zip(&plan.impulses)
        .enumerate()
    {
        if index % 64 == 0 {
            poll(gate)?;
        }
        let observed = v(read.momentum_change.map(|x| x * scales.impulse));
        finite_vec(observed)?;
        let error = difference(observed, command);
        local_linear = local_linear + dd(l1(error));
        local_angular = local_angular + dd(l1(vector(cross(cell.arm, error))));
        let angular = cross(cell.arm, observed);
        for axis in 0..3 {
            change[axis] = change[axis] + dd(a(observed)[axis]);
            angular_change[axis] = angular_change[axis] + angular[axis];
        }
        let mean = vector(core::array::from_fn(|axis| {
            dd(0.5)
                * (dd(read.before.velocity[axis]) + dd(read.velocity_after[axis]))
                * dd(units.velocity_unit_m_s)
        }));
        local_work = local_work + dot(error, mean).abs();
    }
    let change = vector(change);
    let angular_change = vector(angular_change);
    let body_linear = mechanics.linear_increment_world;
    let body_angular = orientation.rotate_body_to_world(mechanics.angular_increment_body);
    let linear_residual = change.add(body_linear);
    let angular_residual = angular_change.add(body_angular);
    let linear_defect = finite(local_linear.to_f64() + l1(linear_residual) + source_linear)?;
    let angular_defect = finite(local_angular.to_f64() + l1(angular_residual) + source_angular)?;
    check(
        "linear momentum and reciprocal source response",
        linear_defect,
        controls.maximum_momentum_defect_kg_m_s,
    )?;
    check(
        "angular momentum and reciprocal source response",
        angular_defect,
        controls.maximum_angular_momentum_defect_kg_m2_s,
    )?;
    let native = candidate.report();
    let liquid_kinetic = finite(native.kinetic_energy_change * scales.energy)?;
    let kinetic_residual =
        (dd(liquid_kinetic) + dd(mechanics.work.kinetic_energy_change_j) + dd(loss)).to_f64();
    let work_defect = finite(
        kinetic_residual.abs()
            + local_work.to_f64()
            + native.absolute_work_defect * scales.energy
            + mechanics.work.work_energy_residual_j.abs()
            + source_work_defect,
    )?;
    check(
        "relative drag loss and actual kinetic storage",
        work_defect,
        controls.maximum_work_defect_j,
    )?;
    poll(gate)?;
    // All fallible physics and SI readback precede the native final commit.
    // There is no fallible work between its publication and the body assignment.
    let native_liquid = candidate.commit(gate)?;
    *body = after;
    Ok(LiquidRigidExchangeReport {
        selected_cells: count,
        native_liquid,
        units,
        body: mechanics,
        balance_origin_world_m: center,
        liquid_impulse_kg_m_s: applied,
        liquid_angular_impulse_kg_m2_s: moment,
        liquid_momentum_change_kg_m_s: change,
        liquid_angular_momentum_change_kg_m2_s: angular_change,
        liquid_kinetic_energy_change_j: liquid_kinetic,
        drag_dissipation_j: loss,
        linear_residual_kg_m_s: linear_residual,
        angular_residual_kg_m2_s: angular_residual,
        kinetic_balance_residual_j: kinetic_residual,
        absolute_momentum_defect_kg_m_s: linear_defect,
        absolute_angular_momentum_defect_kg_m2_s: angular_defect,
        absolute_work_defect_j: work_defect,
    })
}

#[cfg(test)]
mod tests;
