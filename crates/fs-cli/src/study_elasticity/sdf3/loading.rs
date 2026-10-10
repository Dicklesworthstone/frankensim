//! One declared reference load law feeds equilibrium AND enriched goal residuals.
//! Surface loads act on retained implicit boundaries, including cavity walls,
//! inside the declared x patch. They do not act on clipping-box faces.
use super::*;
use fs_cutfem::elastic3::surface::SurfaceForce3;
use fs_cutfem::elastic3::dirichlet::EmbeddedDirichletOptions3;

#[derive(Debug, Clone, Copy)]
pub(super) enum Force {
    Pressure(f64),
    Traction([f64; 3]),
}

#[derive(Debug, Clone, Copy)]
pub(super) struct SurfaceSpec {
    pub force: Force,
    /// Patch limits in normalized x coordinates of the declared box.
    pub x_fraction: [f64; 2],
}
impl SurfaceSpec {
    pub(super) fn validate(self, level: u32) -> Result<()> {
        if !(1..=2).contains(&level) {
            return Err(fail("cli-study-sdf3-input", "surface patch requires an admitted initial octree level"));
        }
        let finite_nonzero = match self.force {
            Force::Pressure(p) => p.is_finite() && p != 0.0 && p.abs() <= 1e12,
            Force::Traction(t) => t.iter().all(|v| v.is_finite() && v.abs() <= 1e12)
                && t.iter().any(|v| *v != 0.0),
        };
        let [a, b] = self.x_fraction;
        let cells = f64::from(1_u32 << level);
        if !finite_nonzero || !a.is_finite() || !b.is_finite()
            || !(0.0 <= a && a < b && b <= 1.0)
            || (a * cells).fract() != 0.0 || (b * cells).fract() != 0.0
        {
            return Err(fail("cli-study-sdf3-input",
                "surface pressure/traction must be finite and nonzero, with an ordered x-fraction patch aligned to the initial octree"));
        }
        Ok(())
    }

    pub(super) fn traction(self, p: [f64; 3], normal: [f64; 3],
        bounds: ([f64; 3], [f64; 3])) -> [f64; 3] {
        let x0 = bounds.0[0];
        let length = bounds.1[0] - x0;
        let left = x0 + self.x_fraction[0] * length;
        let right = x0 + self.x_fraction[1] * length;
        if p[0] < left || p[0] > right { return [0.0; 3]; }
        match self.force {
            Force::Pressure(value) => normal.map(|n| -value * n),
            Force::Traction(value) => value,
        }
    }
}

/// Surface loads OR embedded supports require surface rules. All rule families
/// use the same domain, support law and cumulative allowance on EVERY grid.
pub(super) fn build_operator(
    spec: &Spec, bounds: HexCell, tree: &Octree3, domain: &dyn CutSdf3,
    material: &IsotropicElastic, quadrature: &mut QuadratureControl3<'_>,
) -> std::result::Result<AdaptiveElasticity3, ElasticityError3> {
    // This is the shared construction seam for compliance, every goal-refined
    // background, replay, and minimum-volume stress studies. Never replace an
    // authored shape with the legacy height-field placeholder in one branch.
    let domain: &dyn CutSdf3 = match &spec.constructive {
        Some(shape) => {
            // Boolean switches can spoil a height direction away from the zero
            // set. Refine only those unresolved leaves; keep all global work
            // caps and the certified derivative/normal admission unchanged.
            quadrature.set_unresolved_refinement_limit(2)?;
            shape
        }
        None => domain,
    };
    let clamp = |p| spec.fixed.contains(p, spec.bounds);
    let options = ElasticityOptions3 {
        max_cells: spec.leaves, max_dofs: 50_000, ..Default::default()
    };
    if let Some(beta) = spec.fixed.embedded_penalty() {
        let conflict = std::cell::Cell::new(false);
        let patch = |p, n| {
            let selected = spec.fixed.embedded_contains(p, spec.bounds);
            // Diagnose the ACTUAL retained patch, including a surface on a
            // shared band endpoint. Selection itself remains pure and stable.
            // Never remove requested loads to manufacture lower compliance.
            if selected && spec.surfaces.iter().flatten().any(|law|
                law.traction(p, n, spec.bounds).iter().any(|v| *v != 0.0))
            {
                conflict.set(true);
            }
            selected
        };
        let operator = AdaptiveElasticity3::build_with_embedded_dirichlet(
            bounds, tree, domain, material, &clamp, &patch, options,
            EmbeddedDirichletOptions3 { beta }, Default::default(), quadrature,
        )?;
        if conflict.get() {
            return Err(ElasticityError3::Invalid(
                "surface load overlaps an embedded support; declare disjoint physical patches"));
        }
        Ok(operator)
    } else if spec.surfaces.iter().any(Option::is_some) {
        AdaptiveElasticity3::build_with_surface(bounds, tree, domain, material,
            &clamp, options, Default::default(), quadrature)
    } else {
        AdaptiveElasticity3::build(bounds, tree, domain, material,
            &clamp, options, quadrature)
    }
}

/// The lifetime of each law covers the complete numerical invocation. Body
/// and surface contributions within a case are summed; cases are not summed.
/// Pressure is inward-positive on the REFERENCE normal, never a follower load.
pub(super) fn with_laws<T>(spec: &Spec,
    apply: impl FnOnce(&[GoalReferenceLoad3<'_>]) -> T) -> T {
    let body: Vec<_> = spec.loads.iter()
        .map(|(value, _)| move |_: [f64; 3]| *value).collect();
    let surface: Vec<_> = spec.surfaces.iter().map(|law| {
        move |p, n| law.map_or([0.0; 3], |law| law.traction(p, n, spec.bounds))
    }).collect();
    let loads: Vec<_> = spec.loads.iter().enumerate().map(|(i, (force, weight))| {
        let body: Option<&dyn Fn([f64; 3]) -> [f64; 3]> =
            if force.iter().all(|v| *v == 0.0) { None } else { Some(&body[i]) };
        GoalReferenceLoad3 {
            load: ReferenceLoad3 {
                body,
                surface: spec.surfaces[i].map(|_| SurfaceForce3::Traction(&surface[i])),
            },
            weight: *weight,
        }
    }).collect();
    apply(&loads)
}

#[cfg(test)]
#[path = "loading_tests.rs"]
mod tests;
#[cfg(test)]
#[path = "embedded_tests.rs"]
mod embedded_tests;
