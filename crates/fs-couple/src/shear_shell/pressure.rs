//! Original midsurface ownership for exterior and interior fluid work.
//!
//! Mechanical members need not border either fluid. Both pressure rows project
//! onto the one assembled elastic basis; selecting a surface creates no extra
//! structural copy, cavity, or fitted gain. This is a signed midsurface trace,
//! not a finite-thickness skin extraction or closed-cavity certification.

use super::{AcousticRealizeError, ResolvedShearShell, owner, poll};
use fs_exec::CancelGate;
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace, ShellPressurePort};
use fs_plate::shell::profile::ProfileBudget;
use std::borrow::Cow;

/// Explicit, independent fluid-facing subsets of the original shell triangles.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ShellPressureSurfaces {
    /// Exterior normal-volume row used by the existing compact observation.
    pub exterior: Vec<OrientedShellFacet>,
    /// Outward-from-cavity volume row used for reciprocal compression work.
    pub cavity: Vec<OrientedShellFacet>,
}

#[derive(Clone, Debug)]
pub(super) struct ResolvedPressureSurfaces {
    pub(super) source: ShellPressureSurfaces,
    pub(super) exterior: ShellPressurePort,
    pub(super) cavity: ShellPressurePort,
}

impl ResolvedShearShell {
    fn pressure_geometry_budget(&self) -> ProfileBudget {
        ProfileBudget {
            max_nodes: self.budget.max_nodes,
            max_triangles: self.budget.max_triangles,
            max_feature_evaluations: 0,
        }
    }

    /// Bind two independent original facet selections before reduction. Empty
    /// means no pressure participation; repeated facets within a row refuse.
    /// The same facet can belong to both fluids with independently chosen sign.
    /// All triangles retain their original mechanical mass, stiffness and loss.
    /// A failed binding leaves this immutable source unchanged.
    pub fn with_pressure_surfaces(
        &self,
        surfaces: ShellPressureSurfaces,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        poll(gate)?;
        let budget = self.pressure_geometry_budget();
        let exterior =
            ShellPressurePort::new(&self.mesh, surfaces.exterior.clone(), budget, || {
                gate.is_requested()
            })
            .map_err(owner)?;
        let cavity = ShellPressurePort::new(&self.mesh, surfaces.cavity.clone(), budget, || {
            gate.is_requested()
        })
        .map_err(owner)?;
        let mut resolved = self.clone();
        resolved.pressure_surfaces = Some(ResolvedPressureSurfaces {
            source: surfaces,
            exterior,
            cavity,
        });
        poll(gate)?;
        Ok(resolved)
    }

    /// Original explicitly authored ownership. `None` preserves the historical
    /// full, positive midsurface trace for both fluids. `Some` with empty rows
    /// deliberately omits those fluid ports and retains all mechanical modes.
    pub fn pressure_surfaces(&self) -> Option<&ShellPressureSurfaces> {
        self.pressure_surfaces
            .as_ref()
            .map(|surfaces| &surfaces.source)
    }

    pub(super) fn pressure_ports(
        &self,
        gate: &CancelGate,
    ) -> Result<Cow<'_, ResolvedPressureSurfaces>, AcousticRealizeError> {
        poll(gate)?;
        if let Some(surfaces) = &self.pressure_surfaces {
            return Ok(Cow::Borrowed(surfaces));
        }
        let facets: Vec<_> = (0..self.mesh.tris.len())
            .map(|triangle| OrientedShellFacet {
                triangle,
                side: ShellFace::Positive,
            })
            .collect();
        let exterior = ShellPressurePort::new(
            &self.mesh,
            facets.clone(),
            self.pressure_geometry_budget(),
            || gate.is_requested(),
        )
        .map_err(owner)?;
        Ok(Cow::Owned(ResolvedPressureSurfaces {
            source: ShellPressureSurfaces {
                exterior: facets.clone(),
                cavity: facets,
            },
            cavity: exterior.clone(),
            exterior,
        }))
    }
}
