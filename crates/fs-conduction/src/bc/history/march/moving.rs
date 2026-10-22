//! Moving, power-normalized volumetric heating on the original P1 mesh.
//!
//! The declared input is power ABSORBED IN THIS DOMAIN, not incident laser
//! power. At each node the compact ellipsoidal weight is (1-r²)² for r<1.
//! Its P1 interpolant is normalized with the actual lumped nodal volumes, so
//! summing the existing consistent source load recovers the declared watts.
//! Boundary clipping therefore redistributes absorbed power; it does not model
//! optical capture, penetration or reflected/escaped power. P1 support extends
//! over elements incident to active nodes. No melt-pool or resolution claim is
//! inferred from power conservation or a minimum active-node count.

/// Footprint-relative motion limits for scan timestep selection.
pub mod motion;

/// Instantaneous normalized source derivatives for scan optimization.
pub mod sensitivity;

/// Exact segment rates of continuous absorbed-power schedules.
pub mod slew;

use super::{events, invalid, poll};
use crate::{ConductionError, ConductionMesh, ScalarField, COMPONENT_POWER_DIMS};
use fs_exec::Cx;
use fs_qty::{Dims, QtyAny};
use fs_scenario::signal::TimeSignal;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

/// Explicit spatial work caps and power-projection acceptance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MovingHeatBudget {
    /// Maximum mesh vertices visited by one sample.
    pub max_vertices: usize,
    /// Maximum tetrahedra visited by one sample.
    pub max_elements: usize,
    /// Minimum positively weighted, non-orphan nodes when power is nonzero.
    /// This detects a missed footprint, not interpolation accuracy.
    pub min_active_vertices: usize,
    /// Allowed relative error between requested and actually integrated watts.
    pub power_relative_tolerance: f64,
}

/// Immutable world-coordinate scan path, footprint and absorbed-power schedule.
#[derive(Debug, Clone, PartialEq)]
pub struct MovingHeatSource {
    center_m: [TimeSignal; 3],
    half_width_m: [f64; 3],
    absorbed_power_w: TimeSignal,
    budget: MovingHeatBudget,
}

/// One checked projection. Observations refer to this returned nodal field.
#[derive(Debug, Clone, PartialEq)]
pub struct MovingHeatSample {
    /// Volumetric generation [W/m³], directly accepted by ConductionProblem.
    pub source: ScalarField,
    /// Evaluated world-coordinate center [m].
    pub center_m: [f64; 3],
    /// Requested absorbed power [W].
    pub absorbed_power_w: f64,
    /// Integral of the returned P1 source [W], not a copy of the request.
    pub delivered_power_w: f64,
    /// Positively weighted non-orphan nodes; zero when the source is off.
    pub active_vertices: usize,
}

impl MovingHeatSource {
    /// Declare metre-valued center histories, positive world-axis half widths,
    /// a watt-valued absorbed-power history, and explicit sampling limits.
    /// All signal data are validated by their existing owner at the supplied
    /// reference time; later sampled power must also remain nonnegative.
    pub fn new(center_m: [TimeSignal; 3], half_width_m: [f64; 3],
        absorbed_power_w: TimeSignal, reference_time_s: f64, budget: MovingHeatBudget)
        -> Result<Self, ConductionError>
    {
        if budget.max_vertices == 0 || budget.max_elements == 0
            || budget.min_active_vertices == 0 || budget.min_active_vertices > budget.max_vertices
            || !budget.power_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&budget.power_relative_tolerance)
            || half_width_m.iter().any(|v| !v.is_finite() || *v <= 0.0)
        { return Err(invalid("moving heat requires positive spatial caps/half widths, a nonzero support floor and relative power tolerance in [0,1)")); }
        let source = Self { center_m, half_width_m, absorbed_power_w, budget };
        source.values(reference_time_s, false)?;
        Ok(source)
    }

    /// Original position and absorbed-power histories, for event scheduling.
    pub fn signals(&self) -> impl Iterator<Item = &TimeSignal> {
        self.center_m.iter().chain(std::iter::once(&self.absorbed_power_w))
    }

    /// Original sampling policy, unchanged by retries.
    #[must_use]
    pub const fn budget(&self) -> MovingHeatBudget { self.budget }

    /// Project at the original signals' right-continuous endpoint value.
    pub fn at(&self, cx: &Cx<'_>, mesh: &ConductionMesh, time_s: f64)
        -> Result<MovingHeatSample, ConductionError>
    { self.sample(cx, mesh, time_s, false) }

    /// Project with exact interval-end left limits for held path/power tables.
    /// Continuous signals retain endpoint interpolation; no epsilon shifts time.
    pub fn before(&self, cx: &Cx<'_>, mesh: &ConductionMesh, time_s: f64)
        -> Result<MovingHeatSample, ConductionError>
    { self.sample(cx, mesh, time_s, true) }

    fn values(&self, time: f64, before: bool) -> Result<([f64; 3], f64), ConductionError> {
        crate::require_finite("moving heat time", time)?;
        let value = |signal: &TimeSignal, dims: Dims| -> Result<f64, ConductionError> {
            let q: QtyAny = if before { events::value_before(signal, time)? } else {
                signal.eval(time).map_err(|e| invalid(&format!("moving heat signal: {e}")))?
            };
            if q.dims != dims { return Err(ConductionError::Dimensions {
                context: "moving heat source".into(), expected: dims.0, found: q.dims.0,
            }); }
            crate::require_finite("moving heat signal", q.value)?;
            Ok(q.value)
        };
        let center = [value(&self.center_m[0], LENGTH)?, value(&self.center_m[1], LENGTH)?,
            value(&self.center_m[2], LENGTH)?];
        let power = value(&self.absorbed_power_w, COMPONENT_POWER_DIMS)?;
        if power < 0.0 { return Err(invalid("absorbed moving-source power must be nonnegative watts")); }
        Ok((center, power))
    }

    fn sample(&self, cx: &Cx<'_>, mesh: &ConductionMesh, time: f64, before: bool)
        -> Result<MovingHeatSample, ConductionError>
    {
        poll(cx, 0)?;
        let n = mesh.vertex_count();
        if n > self.budget.max_vertices || mesh.element_count() > self.budget.max_elements {
            return Err(invalid("moving heat mesh exceeds its explicit sampling budget"));
        }
        let (center, power) = self.values(time, before)?;
        poll(cx, 0)?;
        if power == 0.0 { return Ok(MovingHeatSample { source: ScalarField::Uniform(0.0),
            center_m: center, absorbed_power_w: power, delivered_power_w: 0.0, active_vertices: 0 }); }
        let mut volumes = zeros(n)?;
        for (e, tet) in mesh.complex().tets.iter().enumerate() {
            if e % 256 == 0 { poll(cx, e)?; }
            let volume = mesh.element_volume(e) / 4.0;
            if !volume.is_finite() || volume <= 0.0 { return Err(invalid("moving heat nodal volume is not representable")); }
            for &v in tet { volumes[v as usize] += volume; }
        }
        let mut density = zeros(n)?;
        let mut maximum = 0.0_f64;
        let mut active = 0;
        for (i, point) in mesh.positions().iter().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            if !volumes[i].is_finite() { return Err(invalid("moving heat nodal volume overflow")); }
            if volumes[i] == 0.0 { continue; }
            let mut radius2 = 0.0;
            for axis in 0..3 {
                let distance = (point[axis] - center[axis]).abs();
                // An overflowing difference of finite coordinates is outside
                // any admitted finite half width, never a large in-domain weight.
                if distance >= self.half_width_m[axis] { radius2 = 1.0; break; }
                let r = distance / self.half_width_m[axis];
                radius2 += r * r;
            }
            let weight = (1.0_f64 - radius2).max(0.0).powi(2);
            density[i] = weight;
            maximum = maximum.max(weight);
            active += usize::from(weight > 0.0);
        }
        if active < self.budget.min_active_vertices || maximum == 0.0 {
            return Err(invalid("nonzero moving heat has insufficient resolved nodal support; refine the mesh or declare a wider footprint"));
        }
        // Scaling the largest weight to one avoids an overflowing P/integral
        // when the entire sampled footprint happens to lie near its edge.
        let mut volume = 0.0;
        for (i, (weight, &v)) in density.iter_mut().zip(&volumes).enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            *weight /= maximum;
            volume += *weight * v;
        }
        let scale = power / volume;
        if !volume.is_finite() || volume <= 0.0 || !scale.is_finite() || scale <= 0.0 {
            return Err(invalid("moving heat normalization is not representable"));
        }
        for (i, v) in density.iter_mut().enumerate() {
            if i % 256 == 0 { poll(cx, i)?; }
            *v *= scale;
        }
        // Independently integrate the returned field with the original P1
        // element load's row sum, retaining the actual mesh volumes and order.
        let mut delivered = 0.0;
        for (e, tet) in mesh.complex().tets.iter().enumerate() {
            if e % 256 == 0 { poll(cx, e)?; }
            let quarter = mesh.element_volume(e) / 4.0;
            for &v in tet { delivered += quarter * density[v as usize]; }
        }
        let discrepancy = (delivered - power).abs() / power;
        if !delivered.is_finite() || !discrepancy.is_finite()
            || discrepancy > self.budget.power_relative_tolerance
        { return Err(invalid("moving heat projection does not meet its declared power tolerance")); }
        poll(cx, n)?;
        Ok(MovingHeatSample { source: ScalarField::Nodal(density), center_m: center,
            absorbed_power_w: power, delivered_power_w: delivered, active_vertices: active })
    }
}

fn zeros(n: usize) -> Result<Vec<f64>, ConductionError> {
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| invalid("moving heat field allocation refused"))?;
    values.resize(n, 0.0);
    Ok(values)
}
