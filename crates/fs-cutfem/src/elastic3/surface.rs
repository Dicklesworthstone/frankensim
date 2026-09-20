//! Reference-configuration traction and pressure loads on retained interfaces.
//! The same implicit field builds both stiffness and oriented surface geometry.
//! Loads are integrated once, then reusable as independent dead-load vectors.
//! Density changes do not change their geometry or scale them implicitly.
use super::*;
use crate::quad3::surface::{surface_cell_rules3, SurfaceOptions3};

/// Integrated surface load and its unmasked physical totals. Resultant and
/// moment INCLUDE the parts applied at clamped nodes; `rhs` excludes those
/// entries since prescribed zero displacements do no virtual work there.
#[derive(Debug, Clone)]
pub struct SurfaceLoad3 {
    /// Nodal force in the operator's independent-displacement ordering.
    pub rhs: Vec<f64>,
    /// Integral of the applied traction vector (force units).
    pub resultant: [f64; 3],
    /// Integral of x cross traction about the coordinate origin (force*length).
    pub moment: [f64; 3],
    /// Numerical area of the entire retained zero-level surface, including
    /// portions on which the caller returned zero traction. Not a wet-area fit.
    pub area: f64,
}

impl CutElasticity3 {
    /// Build the original Cartesian operator and retain its matching oriented
    /// interface. Both traversals use this SAME pure implicit field and shared
    /// quadrature budget. A failed surface preparation returns no operator.
    /// The original bulk, stiffness, ghost and clamp arithmetic is untouched.
    /// Artificial box faces with phi<0 are not included in the loaded surface.
    #[allow(clippy::too_many_arguments)]
    pub fn build_with_surface(domain: HexCell, counts: [usize; 3], sdf: &dyn CutSdf3,
        material: &IsotropicElastic, clamp: &dyn Fn([f64; 3]) -> bool,
        options: ElasticityOptions3, surface_options: SurfaceOptions3,
        control: &mut QuadratureControl3<'_>) -> Result<Self, ElasticityError3> {
        surface_options.validate()?;
        let mut operator = Self::build(domain, counts, sdf, material, clamp, options, control)?;
        operator.integrate_surface(sdf, surface_options, control)?;
        Ok(operator)
    }

    // Called only on the not-yet-published result of either builder. Source SDF
    // identity is established by the construction call, not an unrelated field
    // supplied later to an already accepted physical operator.
    pub(super) fn integrate_surface(&mut self, sdf: &dyn CutSdf3, options: SurfaceOptions3,
        control: &mut QuadratureControl3<'_>) -> Result<(), ElasticityError3> {
        for cell in &mut self.cells {
            control.poll()?;
            let rules = surface_cell_rules3(sdf, cell.bounds, options, control)?;
            cell.rules.retain_surface(rules);
        }
        control.poll()?;
        Ok(())
    }

    /// Integrate the supplied traction per reference surface area. The callback
    /// receives global position and outward unit normal. Return zero to select
    /// only part of the interface; no unreported load normalization is applied.
    ///
    /// This is a DEAD load on the fixed reference surface, not follower pressure,
    /// shape differentiation, or embedded Dirichlet data. Legacy bulk-only
    /// construction refuses this operation instead of claiming zero loading.
    /// Each callback and each cell is bracketed by cancellation checkpoints;
    /// failure exposes neither a partial vector nor incomplete physical totals.
    pub fn surface_load(&self, traction: &dyn Fn([f64; 3], [f64; 3]) -> [f64; 3],
        mut checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<SurfaceLoad3, ElasticityError3> {
        surface_poll(&mut checkpoint)?;
        // Admit the entire stored rule family before executing any load law.
        for cell in &self.cells {
            surface_poll(&mut checkpoint)?;
            if cell.rules.surface().is_none() {
                return Err(ElasticityError3::Invalid("surface rules not retained; use build_with_surface"));
            }
        }
        let mut load = SurfaceLoad3 { rhs: vec![0.0; self.n()], resultant: [0.0; 3], moment: [0.0; 3], area: 0.0 };
        for cell in &self.cells {
            surface_poll(&mut checkpoint)?;
            for node in cell.rules.surface().expect("surface family admitted").points() {
                surface_poll(&mut checkpoint)?;
                let f = traction(node.position, node.normal);
                surface_poll(&mut checkpoint)?;
                if !f.iter().all(|v| v.is_finite()) { return Err(ElasticityError3::Invalid("nonfinite surface traction")); }
                let (values, _) = q1(cell.bounds, node.position);
                for a in 0..8 { if !self.fixed[cell.nodes[a]] { for c in 0..3 {
                    load.rhs[3*cell.nodes[a]+c] += node.weight*values[a]*f[c];
                } } }
                let force = f.map(|v| node.weight*v);
                for c in 0..3 {
                    load.resultant[c] += force[c];
                    load.moment[c] += node.position[(c+1)%3]*force[(c+2)%3]
                        - node.position[(c+2)%3]*force[(c+1)%3];
                }
                load.area += node.weight;
            }
        }
        if !load.area.is_finite() || !load.rhs.iter().chain(&load.resultant).chain(&load.moment).all(|v| v.is_finite()) {
            return Err(ElasticityError3::Invalid("surface load or physical total overflow"));
        }
        surface_poll(&mut checkpoint)?;
        Ok(load)
    }

    /// Positive pressure acts INTO the solid: traction = -pressure * normal.
    /// Signed pressure is allowed (negative means outward tension). Pressure is
    /// measured per reference area and is not updated as the solid deforms.
    pub fn pressure_load(&self, pressure: &dyn Fn([f64; 3]) -> f64,
        checkpoint: impl FnMut() -> ControlFlow<()>) -> Result<SurfaceLoad3, ElasticityError3> {
        self.surface_load(&|p, n| { let value = pressure(p); n.map(|v| -value*v) }, checkpoint)
    }
}
fn surface_poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), ElasticityError3> {
    if checkpoint().is_break() { Err(ElasticityError3::Cancelled) } else { Ok(()) }
}
