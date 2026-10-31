//! Native material-density loads and the shared physical-occupancy chain rule.
use super::*;
use fs_cutfem::elastic3::ElasticityError3;
use std::ops::ControlFlow;
use std::sync::Arc;

#[derive(Clone)]
pub(super) struct MaterialAcceleration3 {
    law: Arc<dyn Fn([f64; 3]) -> [f64; 3] + Send + Sync>,
    max_terms: usize,
}

// Preserve the actual cumulative-control stop instead of collapsing it to the
// native callback's cancellation sentinel.
pub(crate) fn controlled<T>(
    control: &mut SolveControl<'_>,
    action: impl FnOnce(&mut dyn FnMut() -> ControlFlow<()>) -> Result<T, ElasticityError3>,
) -> Result<T, EvaluationStop> {
    control.checkpoint("sdf3-material-body-load")?;
    let mut stopped = None;
    let result = action(
        &mut || match control.checkpoint("sdf3-material-body-load") {
            Ok(()) => ControlFlow::Continue(()),
            Err(stop) => {
                stopped = Some(stop);
                ControlFlow::Break(())
            }
        },
    );
    if let Some(stop) = stopped {
        return Err(stop);
    }
    result.map_err(|_| failure("sdf3-material-body-load"))
}

// Both vectors already belong to the same reduced native displacement space.
// Kept separate from load integration so each case can reuse the common B*rho.
pub(crate) fn add_force(
    rhs: &mut [f64],
    body: Option<&[f64]>,
    control: &mut SolveControl<'_>,
) -> Result<(), EvaluationStop> {
    let Some(body) = body else {
        return Ok(());
    };
    if rhs.len() != body.len() {
        return Err(failure("sdf3-material-body-force-shape"));
    }
    for (i, (value, load)) in rhs.iter_mut().zip(body).enumerate() {
        if i % 192 == 0 {
            control.checkpoint("sdf3-material-body-force")?;
        }
        *value += load;
        if !value.is_finite() {
            return Err(failure("sdf3-material-body-force-overflow"));
        }
    }
    Ok(())
}

impl<O: Sdf3Elasticity> CutDensityStudy3<O> {
    /// Add one fixed reference-configuration acceleration to EVERY independent
    /// experiment in this study. Prepare it with the native operator's
    /// `prepare_material_body_load` before moving that operator into the study.
    /// Each experiment's existing nodal force remains an additional dead load.
    ///
    /// The force is B*rho_bar, where B contains the native source density and
    /// retained physical quadrature. It follows projected physical occupancy,
    /// including protected solid/void regions, independently of SIMP stiffness
    /// and its ersatz floor. An eroded realization uses its eroded occupancy.
    /// Gravity is one application; spatial acceleration must remain independent
    /// of design and deformation. This does not implement follower loading.
    ///
    /// The stencil is bound to this exact operator and displacement space.
    /// Explicitly refitted studies must prepare their own load on the new
    /// operator with the same acceleration law. Use `with_material_acceleration`
    /// to retain that law for DWR probes and automatic refinement. Replacing a
    /// prepared load here discards any previous acceleration recipe.
    pub fn with_material_body_load(
        mut self,
        load: MaterialBodyLoad3,
    ) -> Result<Self, ElasticityError3> {
        self.operator.validate_material_body_load(&load)?;
        self.material_body_load = Some(load);
        self.material_acceleration = None;
        Ok(self)
    }

    /// Prepare a common material acceleration and retain its pure reference-
    /// configuration law for refinement. Every new grid integrates its own
    /// native source density, retained cut quadrature, clamps and hanging-node
    /// constraints; neither nodal forces nor SIMP stiffness scales stand in
    /// for the physical body load.
    ///
    /// The law must remain independent of design, deformation and mutable
    /// external state for the entire study/campaign. Own any environmental
    /// parameters captured by the closure. It is evaluated during preparation
    /// once per grid, never during a density iteration or an adjoint solve.
    /// `max_terms` bounds each grid's native stencil, including refined grids.
    /// A failed or cancelled rebuild leaves the source study intact.
    ///
    /// This supports gravity and spatially varying fixed acceleration, not
    /// follower loads. Zero-centered spectral load envelopes still require a
    /// separate affine uncertainty formulation when this common load is present.
    pub fn with_material_acceleration(
        mut self,
        acceleration: impl Fn([f64; 3]) -> [f64; 3] + Send + Sync + 'static,
        max_terms: usize,
        checkpoint: impl FnMut() -> ControlFlow<()>,
    ) -> Result<Self, ElasticityError3> {
        let recipe = MaterialAcceleration3 {
            law: Arc::new(acceleration),
            max_terms,
        };
        let load =
            self.operator
                .prepare_material_body_load(recipe.law.as_ref(), max_terms, checkpoint)?;
        self.material_body_load = Some(load);
        self.material_acceleration = Some(recipe);
        Ok(self)
    }

    /// Whether this study can rebuild its current native body load on a new
    /// grid. A study without a body load needs no acceleration recipe.
    #[must_use]
    pub fn can_refine_material_body_load(&self) -> bool {
        self.material_body_load.is_none() || self.material_acceleration.is_some()
    }

    pub(crate) fn require_material_body_recipe(
        &self,
        stage: &'static str,
    ) -> Result<(), EvaluationStop> {
        if self.can_refine_material_body_load() {
            Ok(())
        } else {
            Err(failure(stage))
        }
    }

    // Build on the destination's exact native owner before preparing its
    // elasticity backend. Preserve the cumulative control's original stop.
    pub(crate) fn material_body_load_for<F: Sdf3Elasticity>(
        &self,
        operator: &F,
        control: &mut SolveControl<'_>,
    ) -> Result<Option<MaterialBodyLoad3>, EvaluationStop> {
        if self.material_body_load.is_none() {
            return Ok(None);
        }
        self.require_material_body_recipe("sdf3-material-body-refinement-requires-acceleration")?;
        let recipe = self
            .material_acceleration
            .as_ref()
            .expect("admitted acceleration recipe");
        controlled(control, |poll| {
            operator.prepare_material_body_load(recipe.law.as_ref(), recipe.max_terms, poll)
        })
        .map(Some)
    }

    // All automatic refitters use this same bounded rebuild. Install only
    // after complete preparation; the source and candidate are untouched on
    // failure. Physical occupancy remains the refitter's parent-map concern.
    pub(crate) fn inherit_material_acceleration<F: Sdf3Elasticity>(
        &self,
        candidate: &mut CutDensityStudy3<F>,
        control: &mut SolveControl<'_>,
    ) -> Result<(), EvaluationStop> {
        let load = self.material_body_load_for(&candidate.operator, control)?;
        candidate.material_body_load = load;
        candidate.material_acceleration = self.material_acceleration.clone();
        Ok(())
    }

    /// Whether every physical experiment includes a native material body load.
    #[must_use]
    pub const fn has_material_body_load(&self) -> bool {
        self.material_body_load.is_some()
    }

    /// Exact prepared native load, available for physical force inspection.
    #[must_use]
    pub const fn material_body_load(&self) -> Option<&MaterialBodyLoad3> {
        self.material_body_load.as_ref()
    }

    /// Stable digest of the prepared native force map, including its material
    /// attribution. Restarts bind this independently of the in-memory owner
    /// token; a freshly prepared equivalent operator can reproduce the digest.
    /// Other geometry, constitutive and callback bindings remain separate.
    #[must_use]
    pub fn material_body_load_fingerprint(&self) -> Option<fs_blake3::ContentHash> {
        self.material_body_load
            .as_ref()
            .map(MaterialBodyLoad3::fingerprint)
    }

    // Once per realized design, then shared by all independent cases.
    pub(crate) fn material_body_force(
        &self,
        projected: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Option<Vec<f64>>, EvaluationStop> {
        self.material_body_load
            .as_ref()
            .map(|load| controlled(control, |poll| load.load(projected, poll)))
            .transpose()
    }

    // Strong box clamps have zero admissible displacement and receive no
    // reduced load. Embedded support motion is added separately by its owner.
    pub(super) fn external_force(
        &self,
        force: &[f64],
        body: Option<&[f64]>,
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, EvaluationStop> {
        if force.len() != self.operator.n() || body.is_some_and(|b| b.len() != force.len()) {
            return Err(failure("sdf3-external-force-shape"));
        }
        let mut rhs = Vec::with_capacity(force.len());
        for (i, &f) in force.iter().enumerate() {
            if i % 192 == 0 {
                control.checkpoint("sdf3-external-force")?;
            }
            let value = if self.operator.fixed()[i / 3] {
                0.0
            } else {
                f + body.map_or(0.0, |b| b[i])
            };
            if !value.is_finite() {
                return Err(failure("sdf3-external-force-overflow"));
            }
            rhs.push(value);
        }
        Ok(rhs)
    }

    // Add dJ/df * df/drho_bar before the projection/filter transpose. The
    // caller supplies 2*w for compliance and 1 for a solved response adjoint.
    pub(super) fn add_material_body_gradient(
        &self,
        dual: &[f64],
        coefficient: f64,
        occupancy_gradient: &mut [f64],
        control: &mut SolveControl<'_>,
    ) -> Result<(), EvaluationStop> {
        let Some(load) = &self.material_body_load else {
            return Ok(());
        };
        if occupancy_gradient.len() != self.cells() || !coefficient.is_finite() {
            return Err(failure("sdf3-material-body-gradient"));
        }
        if coefficient == 0.0 {
            return Ok(());
        }
        let work = controlled(control, |poll| load.occupancy_pullback(dual, poll))?;
        for (i, (gradient, value)) in occupancy_gradient.iter_mut().zip(work).enumerate() {
            if i % 64 == 0 {
                control.checkpoint("sdf3-material-body-gradient")?;
            }
            *gradient += coefficient * value;
            if !gradient.is_finite() {
                return Err(failure("sdf3-material-body-gradient"));
            }
        }
        Ok(())
    }

    // One shared chain rule, including direct physical-occupancy effects.
    // No additional elasticity/adjoint or filter-transpose solve is introduced.
    pub(super) fn physical_gradient(
        &self,
        design: &Design3,
        scale_gradient: &[f64],
        occupancy_gradient: &[f64],
        control: &mut SolveControl<'_>,
    ) -> Result<Vec<f64>, EvaluationStop> {
        if scale_gradient.len() != self.cells() || occupancy_gradient.len() != self.cells() {
            return Err(failure("sdf3-physical-gradient-shape"));
        }
        let p = self.params;
        let mut local = Vec::with_capacity(self.cells());
        for (i, ((&scale, &occupancy), (&rho, &slope))) in scale_gradient
            .iter()
            .zip(occupancy_gradient)
            .zip(design.projected.iter().zip(&design.slope))
            .enumerate()
        {
            if i % 64 == 0 {
                control.checkpoint("sdf3-physical-gradient")?;
            }
            let power = if p.penal == 1.0 {
                1.0
            } else if rho == 0.0 {
                0.0
            } else {
                fs_math::det::pow(rho, p.penal - 1.0)
            };
            let dsimp = (1.0 - p.e_min) * p.penal * power;
            let value = (scale * dsimp + occupancy) * slope;
            if !value.is_finite() {
                return Err(failure("sdf3-physical-gradient"));
            }
            local.push(value);
        }
        let gradient = self.pullback(&local, control)?;
        if gradient.iter().any(|g| !g.is_finite()) {
            return Err(failure("sdf3-physical-gradient"));
        }
        Ok(gradient)
    }
}
