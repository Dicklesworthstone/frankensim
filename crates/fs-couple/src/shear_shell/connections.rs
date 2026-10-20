//! Source-owned bilateral compliance at physical shell footprints.
//!
//! Ordered interface cards supply stiffness and viscosity per contact area.
//! Exact offset force/moment rows enter the monolithic shell stiffness before
//! its eigenproblem, and retain one passive viscous channel per resultant.
//! These finite footprints do not assert continuum mortar, contact opening,
//! a massive support, material-state evolution or prestressed equilibrium.

use super::{invalid, owner, poll, ResolvedShearShell};
use crate::acoustic_realize::AcousticRealizeError;
use fs_exec::CancelGate;
use fs_matdb::MaterialStateId;
use fs_material::state_point::ResolvedInterfaceStatePoint;
use fs_plate::shell::mitc::connections::{
    ShellConnection, ShellConnectionBudget, ShellConnectionTerm,
};
use fs_qty::{Dims, DynViscosity, Length, Pressure};

/// A normalized quadrature weight at an explicit physical offset from a node.
#[derive(Clone, Debug)]
pub struct ShellInterfacePoint {
    /// Original physical midsurface node.
    pub node: usize,
    /// Positive displacement interpolation/force fraction; each side sums to 1.
    pub weight: f64,
    /// GLOBAL skin or attachment offset [m], supplying its true moment arm.
    pub offset_m: [f64; 3],
}

/// Each interface surface is a named material region or an explicit fixed side.
#[derive(Clone, Debug)]
pub enum ShellInterfaceSide {
    /// Complete weighted footprint belonging to one original material region.
    Region {
        /// Exact region name in the compiled shell source.
        name: String,
        /// Physical displacement/rotation interpolation.
        points: Vec<ShellInterfacePoint>,
    },
    /// Immobile counterpart with its explicitly named material state.
    /// This boundary condition supplies no inertia or exterior radiation.
    Fixed {
        /// Must match this ordered role of the original interface card.
        material: MaterialStateId,
    },
}

/// One bilateral interface resultant, with complete original data and geometry.
#[derive(Clone, Debug)]
pub struct MaterialShellConnection {
    /// Unique connection name.
    pub name: String,
    /// Original ordered interface, including fixed state and selected claims.
    pub interface: ResolvedInterfaceStatePoint,
    /// Explicit stiffness-per-area property [N/m³].
    pub stiffness_property: String,
    /// Explicit viscous-resistance-per-area property [N s/m³].
    pub viscosity_property: String,
    /// Actual associated contact area [m²].
    pub area_m2: f64,
    /// Unit GLOBAL direction of the scalar relative-displacement coordinate.
    pub direction: [f64; 3],
    /// Positive displacement side, matching interface surface A.
    pub left: ShellInterfaceSide,
    /// Negative displacement side, matching interface surface B.
    pub right: ShellInterfaceSide,
}

/// Immutable interface source and its geometry-resolved physical law.
#[derive(Clone, Debug)]
pub struct ResolvedShellConnection {
    input: MaterialShellConnection,
    law: ShellConnection,
    band: (f64, f64),
}
impl ResolvedShellConnection {
    /// Exact source card/state, property selections, area and both footprints.
    pub const fn input(&self) -> &MaterialShellConnection { &self.input }
    /// Complete signed nodal row and area-integrated stiffness/resistance.
    pub const fn law(&self) -> &ShellConnection { &self.law }
    /// Intersected original stiffness/viscosity applicability [rad/s].
    pub const fn omega_band_rad_s(&self) -> (f64, f64) { self.band }
}

/// Dimensions shared by material queries and the physical shell adapter.
pub fn stiffness_per_area_dims() -> Dims {
    Pressure::DIMS.checked_minus(Length::DIMS).expect("fixed SI interface exponents")
}
/// Dimensions shared by material queries and the physical shell adapter.
pub fn viscosity_per_area_dims() -> Dims {
    DynViscosity::DIMS.checked_minus(Length::DIMS).expect("fixed SI interface exponents")
}

fn coefficient(source: &ResolvedInterfaceStatePoint, key: &str, dims: Dims)
    -> Result<(f64, (f64, f64)), AcousticRealizeError>
{
    let p = source.property(key)
        .ok_or_else(|| invalid("missing or ambiguous original shell interface coefficient"))?;
    let value = p.value_si();
    if p.requirement().dims() != dims || !value.is_finite() || value < 0.0
        || !p.is_constant_at_fixed_temperature() {
        return Err(invalid("shell interface requires nonnegative fixed-state stiffness/viscosity in exact SI dimensions"));
    }
    let band = p.angular_frequency_band_rad_s().map_err(owner)?
        .ok_or_else(|| invalid("each shell interface coefficient requires explicit frequency applicability"))?;
    Ok((value, band))
}

fn point_count(side: &ShellInterfaceSide) -> usize {
    match side {
        ShellInterfaceSide::Region { points, .. } => points.len(),
        ShellInterfaceSide::Fixed { .. } => 0,
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]]
}

// Return the actual weighted physical application point for a moving side.
// Positive weights and this centroid preserve rigid translation/rotation work.
fn bind_side(shell: &ResolvedShearShell, side: &ShellInterfaceSide,
    expected: &MaterialStateId, direction: [f64; 3], sign: f64,
    terms: &mut Vec<ShellConnectionTerm>, gate: &CancelGate)
    -> Result<Option<[f64; 3]>, AcousticRealizeError>
{
    let (name, points) = match side {
        ShellInterfaceSide::Fixed { material } => {
            if material != expected {
                return Err(invalid("fixed shell interface counterpart does not match its ordered source role"));
            }
            return Ok(None);
        }
        ShellInterfaceSide::Region { name, points } => (name, points),
    };
    let region = shell.regions.iter().find(|r| &r.input.region.region.name == name)
        .ok_or_else(|| invalid("shell interface addresses an unknown material region"))?;
    if region.input.region.material.material() != expected || points.is_empty() {
        return Err(invalid("shell interface requires a nonempty footprint on its exact ordered material state"));
    }
    let mut sum = 0.0;
    let mut center = [0.0; 3];
    for p in points {
        poll(gate)?;
        if !p.weight.is_finite() || p.weight <= 0.0 || p.offset_m.iter().any(|x| !x.is_finite())
            || p.node >= shell.mesh.nodes.len()
            || !region.input.region.region.triangle_indices.iter()
                .any(|&e| shell.mesh.tris[e].contains(&p.node)) {
            return Err(invalid("shell interface point must have positive finite weight and belong to its material region"));
        }
        sum += p.weight;
        let translation = direction.map(|x| sign * p.weight * x);
        let rotation_m = cross(p.offset_m, direction).map(|x| sign * p.weight * x);
        for c in 0..3 { center[c] += p.weight * (shell.mesh.nodes[p.node][c] + p.offset_m[c]); }
        if translation.iter().chain(&rotation_m).any(|x| !x.is_finite()) {
            return Err(invalid("shell interface force or physical moment arm overflows"));
        }
        terms.push(ShellConnectionTerm { node: p.node, translation, rotation_m });
    }
    if !sum.is_finite() || (sum - 1.0).abs() > 1e-12
        || center.iter().any(|x| !x.is_finite()) {
        return Err(invalid("shell interface footprint needs finite application point and weights summing to one"));
    }
    Ok(Some(center))
}

impl ResolvedShearShell {
    /// Bind all physical interface resultants atomically to original regions.
    /// Coefficients are resolved from ordered interface data, never from the
    /// bulk material name. Missing viscosity does not imply zero or a rigid tie.
    ///
    /// Two moving sides must have the same weighted physical application point,
    /// so common rigid motion produces zero relative displacement. Offsets map
    /// rotations and forces by the exact d + theta cross r work-conjugate law.
    /// Different finite resultants may share nodes. Their area and placement are
    /// explicit caller data; continuum overlap/coverage is not inferred.
    pub fn with_interfaces(&self, inputs: Vec<MaterialShellConnection>,
        budget: ShellConnectionBudget, gate: &CancelGate)
        -> Result<Self, AcousticRealizeError>
    {
        poll(gate)?;
        if budget.max_connections > 256 || budget.max_terms > 4096
            || inputs.len() > budget.max_connections {
            return Err(invalid("material shell interfaces exceed declared count allowance"));
        }
        let count = inputs.iter().try_fold(0usize, |n, i|
            n.checked_add(point_count(&i.left))?.checked_add(point_count(&i.right)));
        if count.is_none_or(|n| n > budget.max_terms) {
            return Err(invalid("material shell interfaces exceed complete footprint allowance"));
        }
        let work = count.and_then(|n| self.mesh.tris.len().checked_mul(3)?
            .checked_add(36)?.checked_mul(n))
            .and_then(|w| inputs.len().checked_mul(self.mesh.nodes.len())?
                .checked_mul(3)?.checked_add(w));
        if work.is_none_or(|w| w > budget.max_work) {
            return Err(invalid("material shell interface binding exceeds complete work allowance"));
        }
        let mut names = std::collections::BTreeSet::new();
        let mut resolved = Vec::with_capacity(inputs.len());
        for input in inputs {
            poll(gate)?;
            if input.name.trim().is_empty() || !names.insert(input.name.clone())
                || input.stiffness_property.trim().is_empty() || input.viscosity_property.trim().is_empty()
                || !input.area_m2.is_finite() || input.area_m2 <= 0.0
                || input.direction.iter().any(|x| !x.is_finite())
                || (dot(input.direction, input.direction) - 1.0).abs() > 1e-12 {
                return Err(invalid("shell interface needs a unique name, explicit properties, positive area and unit direction"));
            }
            let (k, kb) = coefficient(&input.interface, &input.stiffness_property, stiffness_per_area_dims())?;
            let (c, cb) = coefficient(&input.interface, &input.viscosity_property, viscosity_per_area_dims())?;
            let band = (kb.0.max(cb.0), kb.1.min(cb.1));
            let stiffness_n_m = k * input.area_m2;
            let damping_n_s_m = c * input.area_m2;
            if !stiffness_n_m.is_finite() || !damping_n_s_m.is_finite()
                || (stiffness_n_m == 0.0 && damping_n_s_m == 0.0)
                || !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 || band.1 < band.0 {
                return Err(invalid("shell interface has no finite active law or intersected frequency support"));
            }
            let mut terms = Vec::with_capacity(point_count(&input.left) + point_count(&input.right));
            let a = bind_side(self, &input.left, &input.interface.surface_materials()[0],
                input.direction, 1.0, &mut terms, gate)?;
            let b = bind_side(self, &input.right, &input.interface.surface_materials()[1],
                input.direction, -1.0, &mut terms, gate)?;
            match (a, b) {
                (None, None) => return Err(invalid("shell interface cannot connect two immobile counterparts")),
                (Some(a), Some(b)) => {
                    let mut extent = 0.0_f64;
                    for node in &self.mesh.nodes {
                        poll(gate)?;
                        for c in 0..3 { extent = extent.max((node[c] - self.mesh.nodes[0][c]).abs()); }
                    }
                    let coordinate = a.iter().chain(&b).fold(0.0_f64, |s, x| s.max(x.abs()));
                    let tolerance = 1e-10 * extent.max(1e-6) + 64.0 * f64::EPSILON * coordinate;
                    if !(0..3).all(|i| (a[i] - b[i]).abs() <= tolerance) {
                        return Err(invalid("bonded shell footprint application points disagree"));
                    }
                }
                _ => {}
            }
            resolved.push(ResolvedShellConnection {
                input, law: ShellConnection { terms, stiffness_n_m, damping_n_s_m }, band,
            });
        }
        poll(gate)?;
        let mut out = self.clone();
        out.interfaces = resolved;
        out.interface_budget = budget;
        Ok(out)
    }

    /// Original ordered interface sources and their complete physical laws.
    pub fn interfaces(&self) -> &[ResolvedShellConnection] { &self.interfaces }

    /// Data-only interface rebind preserving geometry, area and property keys.
    /// A changed ordered surface identity must still match the original regions.
    pub fn with_interface_state(&self, name: &str, source: &ResolvedInterfaceStatePoint,
        gate: &CancelGate) -> Result<Self, AcousticRealizeError>
    {
        let index = self.interfaces.iter().position(|i| i.input.name == name)
            .ok_or_else(|| invalid("unknown material shell interface"))?;
        let mut inputs: Vec<_> = self.interfaces.iter().map(|i| i.input.clone()).collect();
        inputs[index].interface = source.clone();
        self.with_interfaces(inputs, self.interface_budget, gate)
    }
}
