//! Component power maps: per-component dissipation declared as typed
//! scenario objects, an EXACT total-power audit, and projection onto the
//! volumetric source the operator consumes.
//!
//! Component dissipation is the primary excitation for enclosure and board
//! cooling, and the commonest data-entry error is a power map whose parts do
//! not add up to the system power anyone believes. This module makes the
//! declared total a checked input rather than an assumption, and refuses a
//! mismatch with the per-component table in the message instead of solving a
//! problem nobody meant to pose.
//!
//! # The audit is exact, not approximate
//!
//! Assembly integrates a P1 nodal source with the exact tet mass matrix
//! `∫ λ_a λ_b dV = V(1+δ_ab)/20`. Summing that element load over `a` gives
//! `V/4 · Σ_b f_b`, so the power the discrete operator actually injects is
//!
//! ```text
//!   delivered = Σ_a f_a w_a        with   w_a = (1/4) Σ_{e ∋ a} V_e
//! ```
//!
//! `w_a` is the LUMPED NODAL VOLUME. Because delivered power is exactly that
//! weighted sum, a component holding `P` watts over a bound vertex set `S`
//! only has to take `f_a = P / Σ_{a∈S} w_a` for its delivered power to equal
//! `P` identically — not to discretization order, but to floating point. The
//! audit therefore compares three numbers that should agree exactly, and any
//! drift between them is a real defect rather than mesh resolution.
//!
//! Distributing by lumped nodal volume, rather than by a per-element density
//! smeared to nodes, is what buys that. The alternative loses power at
//! component boundaries where an element is shared, and the loss is invisible
//! precisely because it looks like discretization error.
//!
//! # Superposition
//!
//! Overlapping components sum their nodal densities. Delivered total stays
//! `Σ_c P_c` exactly, because the delivered functional is linear in `f`.
//! Overlap is therefore ADMITTED, not refused: two dies under one heat
//! spreader is a real configuration, and the audit still balances.

use crate::ConductionError;
use crate::field::ScalarField;
use crate::mesh::ConductionMesh;

/// SI exponents of component power, W.
pub const COMPONENT_POWER_DIMS: fs_qty::Dims = fs_qty::Dims([2, 1, -3, 0, 0, 0]);

/// Default relative tolerance for the declared-versus-summed total check.
///
/// Loose enough for a hand-entered datasheet total, tight enough that a
/// forgotten component or a transposed digit cannot pass.
pub const DEFAULT_POWER_TOLERANCE: f64 = 1e-9;

/// Maximum components admitted in one map.
pub const MAX_COMPONENTS: usize = 4_096;

/// Stated uncertainty on a declared component power.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PowerUncertainty {
    /// No band was supplied. An explicit unknown, never a zero-width band.
    Unstated,
    /// Symmetric half-width at a stated confidence.
    HalfWidth {
        /// Half-width in watts.
        half_width: f64,
        /// Confidence in `(0, 1)`.
        confidence: f64,
    },
}

/// One component's declared dissipation and the vertices it heats.
#[derive(Debug, Clone, PartialEq)]
pub struct ComponentPower {
    name: String,
    watts: f64,
    uncertainty: PowerUncertainty,
    vertices: Vec<usize>,
}

impl ComponentPower {
    /// Declare one component's dissipation over a bound vertex set.
    ///
    /// `vertices` is the component's footprint in the mesh. It is stored
    /// sorted and deduplicated: a vertex listed twice would otherwise take a
    /// double share of the component's power while the audit still balanced,
    /// which is the kind of quiet error this module exists to prevent.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for a blank name, a non-finite or
    /// negative power, or an empty vertex set.
    pub fn new(
        name: impl Into<String>,
        watts: f64,
        uncertainty: PowerUncertainty,
        vertices: Vec<usize>,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        admit_component_scalars(&name, watts, uncertainty)?;
        if vertices.is_empty() {
            return Err(power_error(
                &name,
                "component power row binds no vertices",
                "bind the component footprint to at least one mesh vertex",
            ));
        }
        let mut vertices = vertices;
        vertices.sort_unstable();
        vertices.dedup();
        Ok(Self {
            name,
            watts,
            uncertainty,
            vertices,
        })
    }

    /// Component name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Declared dissipation in watts.
    #[must_use]
    pub const fn watts(&self) -> f64 {
        self.watts
    }

    /// Stated uncertainty on the declared power.
    #[must_use]
    pub const fn uncertainty(&self) -> PowerUncertainty {
        self.uncertainty
    }

    /// Bound vertices, sorted and deduplicated.
    #[must_use]
    pub fn vertices(&self) -> &[usize] {
        &self.vertices
    }
}

/// One audited component row.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerAuditRow {
    name: String,
    declared_w: f64,
    bound_volume_m3: f64,
    density_w_per_m3: f64,
    delivered_w: f64,
}

impl PowerAuditRow {
    /// Component name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Declared dissipation, W.
    #[must_use]
    pub const fn declared_w(&self) -> f64 {
        self.declared_w
    }

    /// Lumped volume of the bound vertex set, m³.
    #[must_use]
    pub const fn bound_volume_m3(&self) -> f64 {
        self.bound_volume_m3
    }

    /// Volumetric density applied at the bound vertices, W/m³.
    #[must_use]
    pub const fn density_w_per_m3(&self) -> f64 {
        self.density_w_per_m3
    }

    /// Power the discrete operator actually injects for this component, W.
    #[must_use]
    pub const fn delivered_w(&self) -> f64 {
        self.delivered_w
    }
}

/// The retained total-power audit for one map on one mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerAudit {
    rows: Vec<PowerAuditRow>,
    declared_total_w: f64,
    component_total_w: f64,
    delivered_total_w: f64,
    total_half_width_w: Option<f64>,
}

impl PowerAudit {
    /// Per-component rows, in component order.
    #[must_use]
    pub fn rows(&self) -> &[PowerAuditRow] {
        &self.rows
    }

    /// The system total the caller declared.
    #[must_use]
    pub const fn declared_total_w(&self) -> f64 {
        self.declared_total_w
    }

    /// The sum of the declared component powers.
    #[must_use]
    pub const fn component_total_w(&self) -> f64 {
        self.component_total_w
    }

    /// The power the assembled operator will actually inject.
    ///
    /// Equal to [`PowerAudit::component_total_w`] to floating point for any
    /// admitted map; a difference is a defect, not discretization error.
    #[must_use]
    pub const fn delivered_total_w(&self) -> f64 {
        self.delivered_total_w
    }

    /// Conservative sum of the stated half-widths, or `None` when any
    /// component's uncertainty is [`PowerUncertainty::Unstated`].
    ///
    /// One unstated term makes the total unavailable rather than smaller:
    /// the same rule the series thermal-resistance budget uses. This is an
    /// arithmetic sum, not a probabilistic convolution, and no correlation
    /// among components is modelled.
    #[must_use]
    pub const fn total_half_width_w(&self) -> Option<f64> {
        self.total_half_width_w
    }

    /// Render the per-component table used in refusal messages and logs.
    #[must_use]
    pub fn render_table(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::new();
        let _ = writeln!(
            out,
            "component power map: declared {} W, components sum {} W, delivered {} W",
            self.declared_total_w, self.component_total_w, self.delivered_total_w
        );
        for row in &self.rows {
            let _ = writeln!(
                out,
                "  {}: declared {} W over {} m^3 -> {} W/m^3, delivered {} W",
                row.name,
                row.declared_w,
                row.bound_volume_m3,
                row.density_w_per_m3,
                row.delivered_w
            );
        }
        out
    }
}

/// A validated set of component power rows plus the declared system total.
#[derive(Debug, Clone, PartialEq)]
pub struct PowerMap {
    components: Vec<ComponentPower>,
    declared_total_w: f64,
}

impl PowerMap {
    /// Admit a power map.
    ///
    /// Components are sorted by name for deterministic downstream order, and
    /// a duplicated name refuses: two rows under one name is an ambiguity the
    /// audit cannot report usefully.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for an empty or oversized set, a
    /// duplicated component name, or a non-finite/negative declared total.
    pub fn new(
        components: Vec<ComponentPower>,
        declared_total_w: f64,
    ) -> Result<Self, ConductionError> {
        if components.is_empty() {
            return Err(power_error(
                "<power-map>",
                "power map declares no components",
                "declare at least one component dissipation row",
            ));
        }
        if components.len() > MAX_COMPONENTS {
            return Err(power_error(
                "<power-map>",
                format!(
                    "power map declares {} components, above the admitted maximum {MAX_COMPONENTS}",
                    components.len()
                ),
                "split the model or raise the declared component budget deliberately",
            ));
        }
        if !declared_total_w.is_finite() || declared_total_w < 0.0 {
            return Err(power_error(
                "<power-map>",
                format!("declared system power {declared_total_w} is not finite and non-negative"),
                "declare the system power the component rows are supposed to add up to",
            ));
        }
        let mut components = components;
        components.sort_by(|a, b| a.name.cmp(&b.name));
        for pair in components.windows(2) {
            if pair[0].name == pair[1].name {
                return Err(power_error(
                    &pair[0].name,
                    "component name is duplicated in the power map",
                    "declare each component exactly once; use distinct names for distinct dies",
                ));
            }
        }
        Ok(Self {
            components,
            declared_total_w,
        })
    }

    /// Components in deterministic name order.
    #[must_use]
    pub fn components(&self) -> &[ComponentPower] {
        &self.components
    }

    /// The declared system power.
    #[must_use]
    pub const fn declared_total_w(&self) -> f64 {
        self.declared_total_w
    }

    /// Project the map onto the nodal volumetric source and audit it.
    ///
    /// Returns the `W/m³` field the operator consumes together with the
    /// retained audit. The audit is computed from the SAME nodal field that
    /// is returned, so the reported delivered power is what the solve will
    /// actually receive rather than a parallel estimate of it.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for a vertex index outside the mesh,
    /// a component whose bound vertices carry no lumped volume, or a declared
    /// total that disagrees with the component sum beyond `relative_tolerance`.
    pub fn volumetric_source(
        &self,
        mesh: &ConductionMesh,
        relative_tolerance: f64,
    ) -> Result<(ScalarField, PowerAudit), ConductionError> {
        if !relative_tolerance.is_finite() || relative_tolerance < 0.0 {
            return Err(power_error(
                "<power-map>",
                "power tolerance is not finite and non-negative",
                "supply a non-negative relative tolerance, e.g. DEFAULT_POWER_TOLERANCE",
            ));
        }
        let vertex_count = mesh.vertex_count();
        let lumped = lumped_nodal_volumes(mesh);

        let mut density = vec![0.0f64; vertex_count];
        let mut rows = Vec::with_capacity(self.components.len());
        let mut component_total = 0.0f64;
        let mut half_width_total = Some(0.0f64);

        for component in &self.components {
            let mut bound_volume = 0.0f64;
            for &vertex in &component.vertices {
                if vertex >= vertex_count {
                    return Err(power_error(
                        &component.name,
                        format!(
                            "component binds vertex {vertex}, outside the {vertex_count}-vertex mesh"
                        ),
                        "bind the component footprint to vertices of the mesh being solved",
                    ));
                }
                bound_volume += lumped[vertex];
            }
            if !(bound_volume.is_finite() && bound_volume > 0.0) {
                return Err(power_error(
                    &component.name,
                    "component's bound vertices carry no lumped volume",
                    "bind the component to vertices that belong to at least one element",
                ));
            }
            // Delivered power is Σ_a f_a w_a, so this density makes the
            // component's delivered power exactly its declared power.
            let component_density = component.watts / bound_volume;
            for &vertex in &component.vertices {
                density[vertex] += component_density;
            }
            component_total += component.watts;
            half_width_total = match (half_width_total, component.uncertainty) {
                (Some(running), PowerUncertainty::HalfWidth { half_width, .. }) => {
                    Some(running + half_width)
                }
                _ => None,
            };
            rows.push(PowerAuditRow {
                name: component.name.clone(),
                declared_w: component.watts,
                bound_volume_m3: bound_volume,
                density_w_per_m3: component_density,
                delivered_w: component.watts,
            });
        }

        // Recompute delivered power from the ASSEMBLED field rather than from
        // the per-component bookkeeping, so overlap and any arithmetic drift
        // show up here instead of being assumed away.
        let mut delivered_total = 0.0f64;
        for (value, weight) in density.iter().zip(lumped.iter()) {
            delivered_total = value.mul_add(*weight, delivered_total);
        }
        for (value, field) in [
            (component_total, "component power sum"),
            (delivered_total, "delivered power"),
        ] {
            if !value.is_finite() {
                return Err(power_error(
                    "<power-map>",
                    format!("{field} is not finite"),
                    "check the declared component powers and the mesh volumes",
                ));
            }
        }

        let audit = PowerAudit {
            rows,
            declared_total_w: self.declared_total_w,
            component_total_w: component_total,
            delivered_total_w: delivered_total,
            total_half_width_w: half_width_total,
        };

        let allowed = relative_tolerance * self.declared_total_w.abs().max(component_total.abs());
        if (component_total - self.declared_total_w).abs() > allowed {
            return Err(power_error(
                "<power-map>",
                format!(
                    "declared system power {} W does not match the component sum {} W\n{}",
                    self.declared_total_w,
                    component_total,
                    audit.render_table()
                ),
                "correct the declared total or the component rows; a power map that does not add up is an entry error, not a modelling choice",
            ));
        }

        let field = ScalarField::nodal("component power source", vertex_count, density)?;
        Ok((field, audit))
    }
}

/// One component's declared dissipation across a set of BOUNDARY FACES.
///
/// The volumetric analogue heats a solid region; this heats a surface — a
/// heater film, a declared surface load, a component that couples through its
/// footprint rather than its bulk.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceComponentPower {
    name: String,
    watts: f64,
    uncertainty: PowerUncertainty,
    faces: Vec<usize>,
}

impl SurfaceComponentPower {
    /// Declare one surface component over boundary-face SLOTS.
    ///
    /// `faces` indexes [`ConductionMesh::boundary`]. Slots are stored sorted
    /// and deduplicated, for the same reason the volumetric path dedups
    /// vertices: a repeated face would take a double share while the totals
    /// still balanced.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for a blank name, a non-finite or
    /// negative power, a malformed uncertainty, or an empty face set.
    pub fn new(
        name: impl Into<String>,
        watts: f64,
        uncertainty: PowerUncertainty,
        faces: Vec<usize>,
    ) -> Result<Self, ConductionError> {
        let name = name.into();
        admit_component_scalars(&name, watts, uncertainty)?;
        if faces.is_empty() {
            return Err(power_error(
                &name,
                "surface component power row binds no boundary faces",
                "bind the component footprint to at least one boundary face slot",
            ));
        }
        let mut faces = faces;
        faces.sort_unstable();
        faces.dedup();
        Ok(Self {
            name,
            watts,
            uncertainty,
            faces,
        })
    }

    /// Component name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Declared dissipation in watts, delivered INTO the solid.
    #[must_use]
    pub const fn watts(&self) -> f64 {
        self.watts
    }

    /// Stated uncertainty on the declared power.
    #[must_use]
    pub const fn uncertainty(&self) -> PowerUncertainty {
        self.uncertainty
    }

    /// Bound boundary-face slots, sorted and deduplicated.
    #[must_use]
    pub fn faces(&self) -> &[usize] {
        &self.faces
    }
}

/// One audited surface-component row.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfacePowerAuditRow {
    name: String,
    declared_w: f64,
    face_area_m2: f64,
    lumped_area_m2: f64,
    outward_flux_w_per_m2: f64,
    delivered_w: f64,
}

impl SurfacePowerAuditRow {
    /// Component name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Declared dissipation, W.
    #[must_use]
    pub const fn declared_w(&self) -> f64 {
        self.declared_w
    }

    /// Geometric area of the component's own faces, m² — the intuitive
    /// footprint number.
    #[must_use]
    pub const fn face_area_m2(&self) -> f64 {
        self.face_area_m2
    }

    /// Lumped nodal area of the component's vertices over the WHOLE map, m².
    ///
    /// This, not [`SurfacePowerAuditRow::face_area_m2`], is the weight the
    /// distribution divides by. The two differ when a component shares
    /// vertices with a neighbouring footprint, because the shared vertex
    /// carries area from both; using the map-wide weight is exactly what
    /// keeps the delivered total exact in that case.
    #[must_use]
    pub const fn lumped_area_m2(&self) -> f64 {
        self.lumped_area_m2
    }

    /// The OUTWARD flux applied at this component's vertices, W/m².
    ///
    /// Negative: a surface heat source drives heat INTO the solid, and
    /// [`crate::bc::ThermalBc::Neumann`] is signed positive-leaving.
    #[must_use]
    pub const fn outward_flux_w_per_m2(&self) -> f64 {
        self.outward_flux_w_per_m2
    }

    /// Power the discrete operator injects for this component, W.
    #[must_use]
    pub const fn delivered_w(&self) -> f64 {
        self.delivered_w
    }
}

/// The retained audit for one surface map on one mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfacePowerAudit {
    rows: Vec<SurfacePowerAuditRow>,
    declared_total_w: f64,
    component_total_w: f64,
    delivered_total_w: f64,
    total_half_width_w: Option<f64>,
}

impl SurfacePowerAudit {
    /// Per-component rows, in component order.
    #[must_use]
    pub fn rows(&self) -> &[SurfacePowerAuditRow] {
        &self.rows
    }

    /// The system total the caller declared.
    #[must_use]
    pub const fn declared_total_w(&self) -> f64 {
        self.declared_total_w
    }

    /// The sum of the declared component powers.
    #[must_use]
    pub const fn component_total_w(&self) -> f64 {
        self.component_total_w
    }

    /// The power the assembled Neumann load will actually inject.
    #[must_use]
    pub const fn delivered_total_w(&self) -> f64 {
        self.delivered_total_w
    }

    /// Conservative sum of stated half-widths, `None` if any is unstated.
    #[must_use]
    pub const fn total_half_width_w(&self) -> Option<f64> {
        self.total_half_width_w
    }

    /// Render the per-component table used in refusals and logs.
    #[must_use]
    pub fn render_table(&self) -> String {
        use core::fmt::Write as _;
        let mut out = String::new();
        let _ = writeln!(
            out,
            "surface power map: declared {} W, components sum {} W, delivered {} W",
            self.declared_total_w, self.component_total_w, self.delivered_total_w
        );
        for row in &self.rows {
            let _ = writeln!(
                out,
                "  {}: declared {} W over {} m^2 (lumped {} m^2) -> outward flux {} W/m^2, delivered {} W",
                row.name,
                row.declared_w,
                row.face_area_m2,
                row.lumped_area_m2,
                row.outward_flux_w_per_m2,
                row.delivered_w
            );
        }
        out
    }
}

/// A validated set of surface component rows plus the declared total.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfacePowerMap {
    components: Vec<SurfaceComponentPower>,
    declared_total_w: f64,
}

impl SurfacePowerMap {
    /// Admit a surface power map.
    ///
    /// Components are sorted by name; a duplicated name refuses, and so does a
    /// boundary face claimed by two components — a face carries one surface
    /// load, and two components sharing it is a declaration error rather than
    /// a superposition the audit could report usefully.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for an empty or oversized set, a
    /// duplicated component name, a face claimed twice, or a non-finite or
    /// negative declared total.
    pub fn new(
        components: Vec<SurfaceComponentPower>,
        declared_total_w: f64,
    ) -> Result<Self, ConductionError> {
        if components.is_empty() {
            return Err(power_error(
                "<surface-power-map>",
                "surface power map declares no components",
                "declare at least one surface dissipation row",
            ));
        }
        if components.len() > MAX_COMPONENTS {
            return Err(power_error(
                "<surface-power-map>",
                format!(
                    "surface power map declares {} components, above the admitted maximum {MAX_COMPONENTS}",
                    components.len()
                ),
                "split the model or raise the declared component budget deliberately",
            ));
        }
        if !declared_total_w.is_finite() || declared_total_w < 0.0 {
            return Err(power_error(
                "<surface-power-map>",
                format!("declared surface power {declared_total_w} is not finite and non-negative"),
                "declare the surface power the component rows are supposed to add up to",
            ));
        }
        let mut components = components;
        components.sort_by(|a, b| a.name.cmp(&b.name));
        for pair in components.windows(2) {
            if pair[0].name == pair[1].name {
                return Err(power_error(
                    &pair[0].name,
                    "surface component name is duplicated in the map",
                    "declare each surface component exactly once",
                ));
            }
        }
        let mut claimed: Vec<(usize, &str)> = Vec::new();
        for component in &components {
            for &face in &component.faces {
                claimed.push((face, component.name.as_str()));
            }
        }
        claimed.sort_unstable();
        for pair in claimed.windows(2) {
            if pair[0].0 == pair[1].0 {
                return Err(power_error(
                    pair[1].1,
                    format!(
                        "boundary face {} is claimed by both {:?} and {:?}",
                        pair[0].0, pair[0].1, pair[1].1
                    ),
                    "a boundary face carries one surface load; merge the components or split the footprint",
                ));
            }
        }
        Ok(Self {
            components,
            declared_total_w,
        })
    }

    /// Components in deterministic name order.
    #[must_use]
    pub fn components(&self) -> &[SurfaceComponentPower] {
        &self.components
    }

    /// The declared surface total.
    #[must_use]
    pub const fn declared_total_w(&self) -> f64 {
        self.declared_total_w
    }

    /// Every boundary-face slot the map covers, sorted.
    #[must_use]
    pub fn covered_faces(&self) -> Vec<usize> {
        let mut faces: Vec<usize> = self
            .components
            .iter()
            .flat_map(|component| component.faces.iter().copied())
            .collect();
        faces.sort_unstable();
        faces
    }

    /// Project the map onto the nodal OUTWARD-flux field and audit it.
    ///
    /// The returned field is the `outward_flux` for a
    /// [`crate::bc::ThermalBc::Neumann`] region covering EXACTLY the faces
    /// this map claims. Use [`SurfacePowerMap::region_selector`] to declare
    /// that region.
    ///
    /// Values are NEGATIVE where the map applies: `ThermalBc::Neumann` is
    /// signed positive-leaving, and a surface source drives heat inward.
    ///
    /// # Errors
    /// [`ConductionError::ScenarioRow`] for a face slot outside the mesh
    /// boundary, a component whose faces carry no lumped area, or a declared
    /// total disagreeing with the component sum beyond `relative_tolerance`.
    pub fn neumann_flux(
        &self,
        mesh: &ConductionMesh,
        relative_tolerance: f64,
    ) -> Result<(ScalarField, SurfacePowerAudit), ConductionError> {
        if !relative_tolerance.is_finite() || relative_tolerance < 0.0 {
            return Err(power_error(
                "<surface-power-map>",
                "surface power tolerance is not finite and non-negative",
                "supply a non-negative relative tolerance, e.g. DEFAULT_POWER_TOLERANCE",
            ));
        }
        let boundary = mesh.boundary();
        self.admit_face_slots(boundary.len())?;

        // ONE map-wide lumped nodal area, used both to distribute and to
        // deliver. Anything else breaks the identity when two footprints share
        // a vertex: the assembled load integrates over the whole region, so
        // the distribution weight has to be the whole region's too.
        let vertex_count = mesh.vertex_count();
        let lumped = self.lumped_nodal_areas(mesh);

        let mut flux = vec![0.0f64; vertex_count];
        let mut rows = Vec::with_capacity(self.components.len());
        let mut component_total = 0.0f64;
        let mut half_width_total = Some(0.0f64);

        for component in &self.components {
            let mut vertices: Vec<usize> = component
                .faces
                .iter()
                .flat_map(|&slot| boundary[slot].vertices)
                .map(|vertex| vertex as usize)
                .collect();
            vertices.sort_unstable();
            vertices.dedup();

            let face_area: f64 = component
                .faces
                .iter()
                .map(|&slot| boundary[slot].area)
                .sum();
            let lumped_area: f64 = vertices.iter().map(|&vertex| lumped[vertex]).sum();
            if !(lumped_area.is_finite() && lumped_area > 0.0) {
                return Err(power_error(
                    &component.name,
                    "surface component's bound faces carry no lumped area",
                    "bind the component to boundary faces with positive area",
                ));
            }

            // Negative: positive Neumann flux LEAVES the domain.
            let outward = -component.watts / lumped_area;
            for &vertex in &vertices {
                flux[vertex] += outward;
            }
            component_total += component.watts;
            half_width_total = match (half_width_total, component.uncertainty) {
                (Some(running), PowerUncertainty::HalfWidth { half_width, .. }) => {
                    Some(running + half_width)
                }
                _ => None,
            };
            rows.push(SurfacePowerAuditRow {
                name: component.name.clone(),
                declared_w: component.watts,
                face_area_m2: face_area,
                lumped_area_m2: lumped_area,
                outward_flux_w_per_m2: outward,
                delivered_w: component.watts,
            });
        }

        // Delivered power INTO the solid is minus the outward load.
        let mut delivered_total = 0.0f64;
        for (value, weight) in flux.iter().zip(lumped.iter()) {
            delivered_total = (-value).mul_add(*weight, delivered_total);
        }
        for (value, field) in [
            (component_total, "surface component power sum"),
            (delivered_total, "delivered surface power"),
        ] {
            if !value.is_finite() {
                return Err(power_error(
                    "<surface-power-map>",
                    format!("{field} is not finite"),
                    "check the declared component powers and the boundary areas",
                ));
            }
        }

        let audit = SurfacePowerAudit {
            rows,
            declared_total_w: self.declared_total_w,
            component_total_w: component_total,
            delivered_total_w: delivered_total,
            total_half_width_w: half_width_total,
        };

        let allowed = relative_tolerance * self.declared_total_w.abs().max(component_total.abs());
        if (component_total - self.declared_total_w).abs() > allowed {
            return Err(power_error(
                "<surface-power-map>",
                format!(
                    "declared surface power {} W does not match the component sum {} W\n{}",
                    self.declared_total_w,
                    component_total,
                    audit.render_table()
                ),
                "correct the declared total or the component rows; a surface power map that does not add up is an entry error, not a modelling choice",
            ));
        }

        let field = ScalarField::nodal("surface power flux", vertex_count, flux)?;
        Ok((field, audit))
    }

    /// Refuse any bound face slot outside the mesh boundary.
    fn admit_face_slots(&self, boundary_len: usize) -> Result<(), ConductionError> {
        for component in &self.components {
            for &slot in &component.faces {
                if slot >= boundary_len {
                    return Err(power_error(
                        &component.name,
                        format!(
                            "component binds boundary face slot {slot}, outside the {boundary_len}-face boundary"
                        ),
                        "bind the footprint to boundary-face slots of the mesh being solved",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Map-wide lumped nodal areas `s_a = (1/3) Σ_{f ∋ a} A_f` over the faces
    /// this map claims.
    ///
    /// Summing the exact P1 face load `∫λ_a λ_b dA = A(1+δ_ab)/12` over `a`
    /// gives `A/3`, so these are the weights that turn a nodal outward flux
    /// into a delivered surface power. Slots are assumed already admitted by
    /// [`SurfacePowerMap::admit_face_slots`].
    fn lumped_nodal_areas(&self, mesh: &ConductionMesh) -> Vec<f64> {
        let boundary = mesh.boundary();
        let mut lumped = vec![0.0f64; mesh.vertex_count()];
        for &slot in &self.covered_faces() {
            let face = &boundary[slot];
            let third = face.area / 3.0;
            for vertex in face.vertices {
                lumped[vertex as usize] += third;
            }
        }
        lumped
    }

    /// A face predicate selecting exactly the faces this map claims.
    ///
    /// Pass it straight to [`crate::bc::ThermalBoundaryBuilder::region`]
    /// alongside the flux from [`SurfacePowerMap::neumann_flux`]. The two must
    /// describe the SAME face set: the lumping weight assumes the assembled
    /// Neumann load integrates over exactly these faces, and a region
    /// selecting more or fewer would deliver a different power than the audit
    /// reports — silently, because both halves would still look individually
    /// reasonable. Deriving the region from the map is what makes that
    /// mismatch unrepresentable rather than merely discouraged.
    ///
    /// # Two numbering spaces, and why this takes the mesh
    ///
    /// A component binds boundary-face SLOTS (indices into
    /// [`ConductionMesh::boundary`]), but the predicate is handed a
    /// [`crate::mesh::BoundaryFace`] whose `face` field is an index into the
    /// complex's CANONICAL FACE TABLE. Those are different numbering spaces
    /// and they are both plain `usize`, so comparing one to the other
    /// compiles, runs, and selects the wrong faces. Resolving slots through
    /// the mesh here is what keeps the two straight; the mesh argument is not
    /// ceremony, it is the translation.
    pub fn region_selector(
        &self,
        mesh: &ConductionMesh,
    ) -> impl Fn(&crate::mesh::BoundaryFace) -> bool + use<> {
        let boundary = mesh.boundary();
        let mut canonical: Vec<usize> = self
            .covered_faces()
            .into_iter()
            .filter_map(|slot| boundary.get(slot).map(|face| face.face))
            .collect();
        canonical.sort_unstable();
        move |face: &crate::mesh::BoundaryFace| canonical.binary_search(&face.face).is_ok()
    }
}

/// Shared scalar admission for a volumetric or surface component row.
fn admit_component_scalars(
    name: &str,
    watts: f64,
    uncertainty: PowerUncertainty,
) -> Result<(), ConductionError> {
    if name.trim().is_empty() {
        return Err(power_error(
            "<unnamed>",
            "component power row has a blank name",
            "name the component as the scenario declares it",
        ));
    }
    if !watts.is_finite() || watts < 0.0 {
        return Err(power_error(
            name,
            format!("component dissipation {watts} is not finite and non-negative"),
            "declare a non-negative watt value; a heat SINK is a boundary condition, not a negative source",
        ));
    }
    if let PowerUncertainty::HalfWidth {
        half_width,
        confidence,
    } = uncertainty
    {
        if !half_width.is_finite() || half_width < 0.0 {
            return Err(power_error(
                name,
                "power half-width is not finite and non-negative",
                "supply a non-negative half-width or PowerUncertainty::Unstated",
            ));
        }
        if !confidence.is_finite() || confidence <= 0.0 || confidence >= 1.0 {
            return Err(power_error(
                name,
                "power confidence is not strictly inside (0, 1)",
                "state the confidence the half-width was quoted at",
            ));
        }
    }
    Ok(())
}

/// Lumped nodal volumes `w_a = (1/4) Σ_{e ∋ a} V_e`.
///
/// These are exactly the weights that turn the assembled P1 source load into
/// a delivered power, which is why the projection uses them rather than a
/// geometric partition of the component footprint.
fn lumped_nodal_volumes(mesh: &ConductionMesh) -> Vec<f64> {
    let mut lumped = vec![0.0f64; mesh.vertex_count()];
    for element in 0..mesh.element_count() {
        let quarter = mesh.element_volume(element) / 4.0;
        for vertex in mesh.complex().tets[element] {
            lumped[vertex as usize] += quarter;
        }
    }
    lumped
}

fn power_error(
    component: impl Into<String>,
    what: impl Into<String>,
    fix: impl Into<String>,
) -> ConductionError {
    ConductionError::ScenarioRow {
        region: component.into(),
        what: what.into(),
        fix: fix.into(),
    }
}
