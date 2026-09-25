//! One declared reference load law feeds equilibrium AND enriched goal residuals.
//! Surface loads act only on the retained implicit graph, not the box faces.
use super::*;
use fs_cutfem::elastic3::surface::SurfaceForce3;

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

/// Only surface-bearing studies pay for surface rules. Both rule families use
/// the same domain, clamp and shared box/point allowance on EVERY background.
pub(super) fn build_operator(
    spec: &Spec, bounds: HexCell, tree: &Octree3, domain: &dyn CutSdf3,
    material: &IsotropicElastic, quadrature: &mut QuadratureControl3<'_>,
) -> std::result::Result<AdaptiveElasticity3, ElasticityError3> {
    let clamp = |p| spec.fixed.contains(p, spec.bounds);
    let options = ElasticityOptions3 {
        max_cells: spec.leaves, max_dofs: 50_000, ..Default::default()
    };
    if spec.surfaces.iter().any(Option::is_some) {
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
