//! Shared-material layout and sizing for an explicitly supplied frame ground.
//!
//! This is the native pin-jointed, axial truss model: the objective is material
//! volume, not mass. Density supplies catalog mass and, when requested, the
//! body's equivalent nodal load. Bending frames and seismic histories remain
//! separate consumers; no source-backed status is inferred for them here.

use fs_exec::Cx;
use fs_truss::material::body_load::{BodyLoadDesignError, MaterialBodyLoadDesign};
use fs_truss::material::{MaterialLayout, TrussMaterial};
use fs_truss::{GroundStructure, LayoutCase, LayoutLimits, PdhgSettings, TrussConstructionError};

/// Explicit body acceleration and work bound for consistent catalog self-weight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BodyLoadControls {
    /// Gravity or equivalent inertial acceleration in the ground's axes [m/s²].
    pub acceleration_m_s2: [f64; 2],
    /// Maximum solves after the initial prescribed-load design.
    pub maximum_iterations: usize,
}

/// Existing solver/assembly controls and the caller's physical loading policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialLayoutControls {
    /// Native ground-to-LP resource limits.
    pub limits: LayoutLimits,
    /// Native PDHG convergence controls. Exhaustion never publishes a design.
    pub solver: PdhgSettings,
    /// Relative member-force threshold used by the original sizing owner.
    pub prune_fraction: f64,
    /// None means prescribed loads only; member mass is still reported.
    pub body_loads: Option<BodyLoadControls>,
}

impl Default for MaterialLayoutControls {
    fn default() -> Self {
        Self {
            limits: LayoutLimits::default(),
            solver: PdhgSettings::default(),
            prune_fraction: 1e-3,
            body_loads: None,
        }
    }
}

/// One design and the immutable geometry/material bindings that produced it.
pub struct MaterialLayoutReport<'a> {
    /// Original node/member ordering.
    pub ground: &'a GroundStructure,
    /// Original state-point bundles, including source claims and uncertainties.
    pub materials: &'a [TrussMaterial],
    /// One original material-table index per ground member.
    pub assignments: &'a [usize],
    /// Original external loads and supports, before adding body loads.
    pub prescribed_case: &'a LayoutCase,
    /// The exact requested numerical and body-load controls.
    pub controls: MaterialLayoutControls,
    /// Native accepted loads, areas, per-member mass, code checks and LP report.
    /// With prescribed loads only, body forces and iteration count are zero.
    pub design: MaterialBodyLoadDesign,
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), BodyLoadDesignError> {
    cx.checkpoint().map_err(|_| {
        TrussConstructionError::Cancelled {
            stage: "material frame layout",
        }
        .into()
    })
}

/// Lay out and size any admitted ground with per-member source material states.
///
/// The same immutable assignments supply yield costs, Euler stiffness, catalog
/// mass and optional self-weight. Self-weight delegates to the native bounded
/// catalog fixed-point solve: body forces are added once to the original load
/// case, and publication requires unchanged complete catalog sections. Without
/// body loads, a single native solve and sizing pass is sufficient.
///
/// Returned LP diagnostics describe a volume-minimization problem at its actual
/// applied loads. They are not certified optimum bounds, a mass optimum, or a
/// global optimum of the coupled self-weight problem. The original symmetric
/// yield, solid-square, pin-jointed model and nodal body-load approximation apply.
///
/// # Errors
/// Native material assignment, resource, numerical, catalog and body-load errors
/// propagate unchanged. Unconverged solves, failed sizing and cancellation return
/// no report; the caller's geometry, case and material states are never mutated.
pub fn layout_and_size<'a>(
    ground: &'a GroundStructure,
    case: &'a LayoutCase,
    materials: &'a [TrussMaterial],
    assignments: &'a [usize],
    catalog_m2: &[f64],
    controls: MaterialLayoutControls,
    cx: &Cx<'_>,
) -> Result<MaterialLayoutReport<'a>, BodyLoadDesignError> {
    checkpoint(cx)?;
    if catalog_m2.is_empty()
        || catalog_m2.iter().any(|a| !a.is_finite() || *a <= 0.0)
        || catalog_m2.windows(2).any(|a| a[0] >= a[1])
        || !controls.prune_fraction.is_finite()
        || !(0.0..=1.0).contains(&controls.prune_fraction)
    {
        return Err(TrussConstructionError::InvalidInput {
            field: "material frame catalog/pruning",
            requirement: "positive finite increasing areas and finite prune fraction in [0,1]",
        }
        .into());
    }
    let layout =
        MaterialLayout::try_new(ground, case, materials, assignments, controls.limits, cx)?;
    let design = if let Some(body) = controls.body_loads {
        layout.solve_with_body_loads(
            body.acceleration_m_s2,
            catalog_m2,
            controls.prune_fraction,
            controls.solver,
            body.maximum_iterations,
            cx,
        )?
    } else {
        layout.solve_and_size(catalog_m2, controls.prune_fraction, controls.solver, cx)?
    };
    checkpoint(cx)?;
    Ok(MaterialLayoutReport {
        ground,
        materials,
        assignments,
        prescribed_case: case,
        controls,
        design,
    })
}
