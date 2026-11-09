//! Nodal thermal mismatch on one original elastic laminate.
//!
//! The three nodes carry independently integrated layer free strains in the
//! same original principal frames. P1 barycentric interpolation acts on those
//! strains, before constitutive energy is evaluated. It neither interpolates
//! temperature/CTE nor averages independently stored node energies.
//!
//! All nodes must share the SAME immutable section, including layer sources,
//! stiffness, axes, thickness and mass. Temperature-dependent stiffness needs
//! an explicitly resolved material discretization, not this fixed-section
//! interpolation. Evaluations reuse the uniform owner's original thickness
//! quadrature without allocating, rotating materials or looking up sources.

use super::{
    LaminateSection, LaminateThermalError, LaminateThermalResponse, LayerThermalStress,
    PreparedLaminateThermal, finite,
};
use std::sync::Arc;

/// P1 spatial interpolation of three admitted, layer-resolved thermal profiles.
#[derive(Clone, Debug, PartialEq)]
pub struct PreparedNodalLaminateThermal {
    nodes: [Arc<PreparedLaminateThermal>; 3],
}

impl PreparedNodalLaminateThermal {
    /// Bind profiles in the consumer's original triangle-node order.
    /// The same section allocation must own all three profiles; even equal
    /// independently authored sections cannot substitute their source identity.
    pub fn try_new(nodes: [Arc<PreparedLaminateThermal>; 3]) -> Result<Self, LaminateThermalError> {
        if nodes[1..]
            .iter()
            .any(|node| !Arc::ptr_eq(&nodes[0].section, &node.section))
        {
            return Err(LaminateThermalError::IncompatibleSection);
        }
        Ok(Self { nodes })
    }

    /// Original node profiles, retaining every layer's integrated input.
    pub const fn nodes(&self) -> &[Arc<PreparedLaminateThermal>; 3] {
        &self.nodes
    }

    pub fn section(&self) -> &LaminateSection {
        self.nodes[0].section()
    }

    /// Positive equivalent thermal membrane force and first moment at a point.
    /// The cached original integrals are linear in the interpolated free strain.
    pub fn thermal_resultants_at(
        &self,
        barycentric: [f64; 3],
    ) -> Result<([f64; 3], [f64; 3]), LaminateThermalError> {
        admit_barycentric(barycentric)?;
        let force = interpolate(barycentric, self.nodes.each_ref().map(|node| node.force));
        let moment = interpolate(barycentric, self.nodes.each_ref().map(|node| node.moment));
        finite(
            force.iter().chain(&moment),
            "nodal laminate thermal resultants",
        )?;
        Ok((force, moment))
    }

    /// Original layer-face engineering free strains, in section coordinates.
    pub fn layer_free_strain_in_section(
        &self,
        layer: usize,
        barycentric: [f64; 3],
    ) -> Result<[[f64; 3]; 2], LaminateThermalError> {
        admit_barycentric(barycentric)?;
        self.nodes[0].layer_free_strain_in_section(layer)?;
        let free = core::array::from_fn(|face| {
            interpolate(
                barycentric,
                self.nodes
                    .each_ref()
                    .map(|node| node.samples[layer].face_free[face]),
            )
        });
        finite(free.iter().flatten(), "nodal layer free strain")?;
        Ok(free)
    }

    /// Physical resultants and storage at one compatible section strain.
    /// Layer mismatch remains inside the original elastic-strain square.
    pub fn evaluate(
        &self,
        barycentric: [f64; 3],
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
        cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateThermalResponse, LaminateThermalError> {
        admit_barycentric(barycentric)?;
        self.nodes[0].evaluate_with_free(
            membrane,
            gradient_per_m,
            |layer, point| {
                interpolate(
                    barycentric,
                    self.nodes
                        .each_ref()
                        .map(|node| node.samples[layer].free[point]),
                )
            },
            cancelled,
        )
    }

    /// Actual stress at both faces of one original layer at this spatial point.
    pub fn layer_stress(
        &self,
        layer: usize,
        barycentric: [f64; 3],
        membrane: [f64; 3],
        gradient_per_m: [f64; 3],
    ) -> Result<LayerThermalStress, LaminateThermalError> {
        admit_barycentric(barycentric)?;
        self.nodes[0].layer_stress_with_free(layer, membrane, gradient_per_m, |layer| {
            core::array::from_fn(|face| {
                interpolate(
                    barycentric,
                    self.nodes
                        .each_ref()
                        .map(|node| node.samples[layer].face_free[face]),
                )
            })
        })
    }
}

fn admit_barycentric(lambda: [f64; 3]) -> Result<(), LaminateThermalError> {
    if lambda
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0 || *v > 1.0)
        || (lambda.iter().sum::<f64>() - 1.0).abs() > 64.0 * f64::EPSILON
    {
        return Err(LaminateThermalError::BarycentricCoordinates);
    }
    Ok(())
}

fn interpolate(lambda: [f64; 3], values: [[f64; 3]; 3]) -> [f64; 3] {
    if values[0] == values[1] && values[0] == values[2] {
        return values[0];
    }
    core::array::from_fn(|i| (0..3).map(|node| lambda[node] * values[node][i]).sum())
}
