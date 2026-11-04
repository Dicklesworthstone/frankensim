//! Physical Kelvin--Voigt loss in the retained original director coordinates.

use super::*;
use crate::shell::mitc::viscosity::{ShearDissipation, ShearShellViscosity, projection_terms};

#[derive(Clone, Debug)]
pub(super) struct EnrichedViscosity {
    pub(super) matrix: Csr,
    dissipation: Vec<Option<ShearDissipation>>,
}

impl EnrichedMitcShell {
    /// Attach one explicit passive A/B/D/shear viscosity per original facet.
    /// Resultants use the facet strain frame and their original SI units. All
    /// nodal/internal cross terms survive in the fixed reference basis; neither
    /// numerical drilling stabilization nor inertia is damped by this law.
    ///
    /// Rebinding is immutable and changes no accepted state or stored energy.
    /// The law is frozen over each midpoint interval, including intervals with
    /// changing elastic stiffness. Refreshing elastic material retains this
    /// viscosity; temperature-dependent loss must be explicitly rebound by its
    /// source owner between intervals. No spectrum or source data is inferred.
    ///
    /// The existing shell-viscosity work estimate with zero projected modes
    /// bounds element assembly before allocation. This is an operation ceiling,
    /// not a wall-clock or discretization-error guarantee.
    pub fn with_viscosity(
        &self,
        laws: &[ShearShellViscosity],
        max_assembly_terms: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PlateError> {
        check(&mut cancelled)?;
        if laws.len() != self.reference_fields.len() {
            return Err(invalid(
                "one explicit viscosity per original shell facet required",
            ));
        }
        if projection_terms(self.model.mesh.nodes.len(), laws.len(), 0)
            .is_none_or(|terms| terms > max_assembly_terms)
        {
            return Err(invalid(
                "enriched shell viscosity exceeds its assembly work allowance",
            ));
        }
        let mut all_zero = true;
        for law in laws {
            check(&mut cancelled)?;
            all_zero &= law.is_zero();
        }
        let mut out = self.clone();
        if all_zero {
            check(&mut cancelled)?;
            out.viscosity = None;
            return Ok(out);
        }
        let mut matrix = Coo::new(self.dimension(), self.dimension());
        let mut dissipation = Vec::with_capacity(laws.len());
        for (element, (field, law)) in self.reference_fields.iter().zip(laws).enumerate() {
            check(&mut cancelled)?;
            if law.is_zero() {
                dissipation.push(None);
            } else {
                dissipation.push(Some(law.dissipation()?));
                let local = field.complete_viscosity(law, &mut cancelled)?;
                material::scatter(
                    &mut matrix,
                    &self.model,
                    element,
                    &local,
                    1.0,
                    &mut cancelled,
                )?;
            }
        }
        check(&mut cancelled)?;
        out.viscosity = Some(Arc::new(EnrichedViscosity {
            matrix: matrix.assemble(),
            dissipation,
        }));
        check(&mut cancelled)?;
        Ok(out)
    }

    /// Original spatial material operator C, in the same coordinates as M/K.
    /// No attached law (including explicitly all-zero laws) returns None.
    pub fn material_viscosity(&self) -> Option<&Csr> {
        self.viscosity.as_ref().map(|viscosity| &viscosity.matrix)
    }

    pub(super) fn same_viscosity(&self, other: &Self) -> bool {
        match (&self.viscosity, &other.viscosity) {
            (None, None) => true,
            (Some(a), Some(b)) => Arc::ptr_eq(a, b),
            _ => false,
        }
    }

    pub(super) fn material_dissipation(
        &self,
        midpoint_velocity: &[f64],
        dt_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(f64, Vec<f64>), PlateError> {
        let Some(viscosity) = &self.viscosity else {
            return Ok((0.0, Vec::new()));
        };
        let nodal = self.model.pencil.free;
        let velocity = self.model.expand(&midpoint_velocity[..nodal])?;
        let mut facets = Vec::with_capacity(self.reference_fields.len());
        let mut total = (0.0, 0.0);
        for (element, (field, law)) in self
            .reference_fields
            .iter()
            .zip(&viscosity.dissipation)
            .enumerate()
        {
            check(cancelled)?;
            let Some(law) = law else {
                facets.push(0.0);
                continue;
            };
            let triangle = self.model.mesh.tris[element];
            let local = core::array::from_fn(|i| velocity[6 * triangle[i / 6] + i % 6]);
            let first = nodal + 2 * element;
            let bubble = [midpoint_velocity[first], midpoint_velocity[first + 1]];
            let heat = finite(dt_s * field.viscous_power(&local, bubble, law, cancelled)?)?;
            facets.push(heat);
            total = dd_add(total, (heat, 0.0));
        }
        check(cancelled)?;
        Ok((finite(total.0 + total.1)?, facets))
    }
}
