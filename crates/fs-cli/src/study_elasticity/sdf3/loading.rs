//! One declared reference load law feeds equilibrium AND enriched goal residuals.
//! Surface loads act on retained implicit boundaries, including cavity walls,
//! inside the declared x patch. They do not act on clipping-box faces.
use super::*;
use fs_cutfem::elastic3::surface::SurfaceForce3;
use fs_cutfem::elastic3::dirichlet::EmbeddedDirichletOptions3;
use fs_cutfem::elastic3::ElasticMaterialField3;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_topopt::sdf3::Sdf3Elasticity;

/// A fixed reference acceleration, shared by every independent load case.
/// Material occupancy supplies the changing mass; the SIMP stiffness floor
/// never creates weight in a prescribed void.
#[derive(Debug, Clone, Copy)]
pub(super) struct MaterialAcceleration {
    pub density_kg_m3: f64,
    pub acceleration_m_s2: [f64; 3],
    pub max_terms: usize,
}

impl MaterialAcceleration {
    pub(super) fn validate(self) -> Result<()> {
        if !(1e-9..=1e8).contains(&self.density_kg_m3)
            || self.acceleration_m_s2.iter().any(|a| !a.is_finite() || a.abs() > 1e6)
            || self.acceleration_m_s2.iter().all(|a| *a == 0.0)
            || !(1..=4_000_000).contains(&self.max_terms)
        {
            return Err(fail("cli-study-sdf3-input",
                "material acceleration needs positive density in kg/m3, a finite nonzero vector in m/s2, and 1..=4000000 stencil terms"));
        }
        Ok(())
    }

    /// Covers simultaneous coarse/fine maps and sparse preparation workspace.
    /// This is a structural admission envelope, not measured peak RSS.
    pub(super) fn memory_envelope(self) -> usize { 128 * self.max_terms }
}

pub(super) fn material(spec: &Spec) -> Result<IsotropicElastic> {
    IsotropicElastic::new(spec.youngs, spec.poisson, 1.0)
        .map_err(|e| fail("cli-study-sdf3-material", e.to_string()))
}

/// Bind once on the actual native operator, before any solve or gradient gate.
/// Retaining the acceleration recipe lets the existing adaptive owner rebuild
/// its map on every enriched/refitted grid with the same reference law.
pub(super) fn bind_material_acceleration<O: Sdf3Elasticity>(
    study: CutDensityStudy3<O>, spec: &Spec, gate: &CancelGate,
    mut checkpoint: impl FnMut() -> ControlFlow<()>,
) -> Result<CutDensityStudy3<O>> {
    let Some(law) = spec.material_acceleration else { return Ok(study); };
    study.with_material_acceleration(move |_| law.acceleration_m_s2,
        law.max_terms, &mut checkpoint).map_err(|e| Failure {
            code: "cli-study-sdf3-material-acceleration",
            message: e.to_string(),
            exit: if gate.is_requested() { exit::CANCELLED }
                else if checkpoint().is_break() || matches!(e,
                    ElasticityError3::Invalid("material body-load term allowance exhausted")) {
                    exit::BUDGET
                } else { exit::REFUSED },
        })
}

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

/// An authored mechanical preload shared by the whole uncertain load ball.
/// It is fixed on the reference domain and independent of material occupancy;
/// the native material-acceleration recipe remains a separate changing load.
#[derive(Debug, Clone, Copy)]
pub(super) struct DeadLoad {
    pub body_n_m3: [f64; 3],
    pub surface: Option<SurfaceSpec>,
}
impl DeadLoad {
    pub(super) fn validate(self, level: u32) -> Result<()> {
        if self.body_n_m3.iter().any(|v| !v.is_finite() || v.abs() > 1e12) {
            return Err(fail("cli-study-sdf3-input",
                "affine dead-body load must be finite and bounded by 1e12 N/m3 per component"));
        }
        if let Some(surface) = self.surface {
            surface.validate(level)?;
        }
        Ok(())
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
    let domain: &dyn CutSdf3 = match spec.authored_domain() {
        Some(shape) => {
            // Boolean switches and sample seams can spoil a height direction.
            // Refine only those unresolved leaves; keep all global work
            // caps and the certified derivative/normal admission unchanged.
            quadrature.set_unresolved_refinement_limit(2)?;
            shape
        }
        None => domain,
    };
    let clamp = |p| spec.fixed.contains(p, spec.bounds);
    let options = ElasticityOptions3 {
        formulation: spec.formulation,
        max_cells: spec.leaves, max_dofs: 50_000, ..Default::default()
    };
    // The legacy E/nu constructor retains no mass density. Self-weight needs
    // the complete immutable native state, tied to this declared study source,
    // on bulk, surface-loaded and embedded-support grids alike.
    let native = spec.material_acceleration.map(|law|
        ElasticAcousticState::from_isotropic(material, law.density_kg_m3, spec.id)
            .map_err(|_| ElasticityError3::Invalid("invalid declared native elastic material")))
        .transpose()?;
    let field = native.as_ref().map(ElasticMaterialField3::Uniform);
    if let Some(beta) = spec.fixed.embedded_penalty() {
        let conflict = std::cell::Cell::new(false);
        let patch = |p, n| {
            let selected = spec.fixed.embedded_contains(p, spec.bounds);
            // Diagnose the ACTUAL retained patch, including a surface on a
            // shared band endpoint. Selection itself remains pure and stable.
            // Never remove requested loads to manufacture lower compliance.
            if selected && spec.surfaces.iter().flatten()
                .chain(spec.affine_dead_load.as_ref().and_then(|dead| dead.surface.as_ref()))
                .any(|law| law.traction(p, n, spec.bounds).iter().any(|v| *v != 0.0))
            {
                conflict.set(true);
            }
            selected
        };
        let operator = if let Some(field) = field {
            AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
                bounds, tree, domain, field, &clamp, &patch, options,
                EmbeddedDirichletOptions3 { beta }, Default::default(), quadrature,
            )
        } else {
            AdaptiveElasticity3::build_with_embedded_dirichlet(
                bounds, tree, domain, material, &clamp, &patch, options,
                EmbeddedDirichletOptions3 { beta }, Default::default(), quadrature,
            )
        }?;
        if conflict.get() {
            return Err(ElasticityError3::Invalid(
                "surface load overlaps an embedded support; declare disjoint physical patches"));
        }
        Ok(operator)
    } else if spec.surfaces.iter().any(Option::is_some)
        || spec.affine_dead_load.is_some_and(|dead| dead.surface.is_some())
    {
        if let Some(field) = field {
            AdaptiveElasticity3::build_with_material_field_and_surface(bounds, tree, domain, field,
                &clamp, options, Default::default(), quadrature)
        } else {
            AdaptiveElasticity3::build_with_surface(bounds, tree, domain, material,
                &clamp, options, Default::default(), quadrature)
        }
    } else {
        if let Some(field) = field {
            AdaptiveElasticity3::build_with_material_field(bounds, tree, domain, field,
                &clamp, options, quadrature)
        } else {
            AdaptiveElasticity3::build(bounds, tree, domain, material,
                &clamp, options, quadrature)
        }
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

/// Keep the known preload separate from every uncertain basis column. All
/// callbacks borrow reference laws for the complete numerical invocation, so
/// assessment reintegrates the same preload on each actual grid. In particular
/// pressure uses reference normals and is never a displacement-following load.
pub(super) fn with_affine_laws<T>(spec: &Spec,
    apply: impl FnOnce(ReferenceLoad3<'_>, &[GoalReferenceLoad3<'_>]) -> T) -> T {
    let body = |_: [f64; 3]| spec.affine_dead_load.map_or([0.0; 3], |dead| dead.body_n_m3);
    let surface = |p, n| spec.affine_dead_load.and_then(|dead| dead.surface)
        .map_or([0.0; 3], |law| law.traction(p, n, spec.bounds));
    let dead = ReferenceLoad3 {
        body: spec.affine_dead_load
            .filter(|dead| dead.body_n_m3.iter().any(|v| *v != 0.0))
            .map(|_| &body as &dyn Fn([f64; 3]) -> [f64; 3]),
        surface: spec.affine_dead_load.and_then(|dead| dead.surface)
            .map(|_| SurfaceForce3::Traction(&surface)),
    };
    with_laws(spec, |laws| apply(dead, laws))
}

#[cfg(test)]
#[path = "loading_tests.rs"]
mod tests;
#[cfg(test)]
#[path = "embedded_tests.rs"]
mod embedded_tests;
