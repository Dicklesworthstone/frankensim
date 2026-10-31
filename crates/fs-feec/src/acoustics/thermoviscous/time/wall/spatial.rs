//! Spatial heat spreading in a declared planar, transversely isotropic slab.
//! The original boundary gradient owner and original depth mass provide
//! K_t = T0 k_t K_surface tensor M_depth. No additional heat state is introduced.

use super::*;

/// In-plane conductivity [W/(m K)] of one homogeneous finite slab.
/// This is an EXPLICIT constitutive assertion: conductivity is isotropic in
/// the tangent plane and has no tangent/normal cross term. Normal conductivity
/// remains the independently supplied original wall value. A scalar measured
/// only along the normal does not supply this coefficient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TangentialWallConductivity {
    pub isotropic_w_m_k: f64,
}

impl PreparedThermoviscousTime<'_, '_> {
    /// Original declared in-plane transport, in conductive-region order.
    /// None denotes the original normal-only material model.
    pub fn tangential_wall_conductivity(&self) -> Option<&[TangentialWallConductivity]> {
        self.conductive_wall_fields.tangential.as_deref()
    }
}

impl<'c> ThermoviscousFe<'c> {
    /// Retain the complete P2 surface by P1 depth heat equation in every
    /// original planar slab, together with all original fluid/structure fields.
    /// Entries follow conductive-region order. Each patch has insulated outer
    /// sides; separate regions do not acquire an implicit tangential bond.
    /// The original front temperature and declared backing are unchanged.
    ///
    /// Actual boundary faces must be coplanar with the same outward normal.
    /// Thickness is physical, never a fitted penetration depth. The old local
    /// normal-profile approximation bounds do not constrain this spatial solve.
    /// Material coefficients and mean geometry remain fixed. Full mesh and time
    /// refinement are still required; the declared band alone proves no accuracy.
    pub fn prepare_time_with_spatial_wall_fields<'a>(
        &'a self,
        duration_s: f64,
        frozen_cp_j_kg_k: f64,
        reservoirs: &[ThermoviscousCaloricReservoir],
        walls: ThermoviscousWallDiscretization<'_>,
        tangential: &[TangentialWallConductivity],
        limits: ThermoviscousTimeLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedThermoviscousTime<'a, 'c>, ThermoviscousTimeError> {
        PreparedThermoviscousTime::new(
            self,
            duration_s,
            frozen_cp_j_kg_k,
            reservoirs,
            Some(walls),
            Some(tangential),
            limits,
            cancelled,
        )
    }
}

pub(super) fn admit(
    model: &ThermoviscousFe<'_>,
    transport: &[TangentialWallConductivity],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ThermoviscousTimeError> {
    let regions = model.conductive_walls();
    if transport.is_empty() || transport.len() != regions.len() {
        return Err(invalid(
            "complete nonempty tangential wall material set required",
        ));
    }
    for (region, tangential) in regions.iter().zip(transport) {
        poll(cancelled)?;
        let ThermalWallDepth::Slab { thickness_m, .. } = region.wall.depth else {
            return Err(invalid(
                "spatial wall diffusion requires the original finite thickness",
            ));
        };
        if [
            tangential.isotropic_w_m_k,
            thickness_m,
            region.wall.transport.conductivity_w_m_k,
            region.wall.transport.volumetric_heat_capacity_j_m3_k,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(invalid(
                "positive finite spatial wall transport and thickness required",
            ));
        }
        let (first, normal, _) = model
            .fe
            .boundary_geometry(region.faces[0])
            .map_err(ThermoviscousError::from)?;
        for &face in &region.faces {
            poll(cancelled)?;
            let (vertices, actual, _) = model
                .fe
                .boundary_geometry(face)
                .map_err(ThermoviscousError::from)?;
            if actual
                .iter()
                .zip(normal)
                .any(|(a, b)| (a - b).abs() > 2048.0 * f64::EPSILON)
            {
                return Err(invalid(
                    "spatial wall extrusion requires one actual outward normal",
                ));
            }
            for vertex in vertices {
                let difference: [f64; 3] = std::array::from_fn(|i| vertex[i] - first[0][i]);
                let distance: f64 = difference.iter().zip(normal).map(|(d, n)| d * n).sum();
                let scale = vertex
                    .iter()
                    .chain(&first[0])
                    .map(|v| v.abs())
                    .fold(thickness_m, f64::max);
                if !distance.is_finite() || distance.abs() > 2048.0 * f64::EPSILON * scale {
                    return Err(invalid("spatial wall faces must lie in one physical plane"));
                }
            }
        }
    }
    Ok(())
}

pub(super) fn assemble(
    fields: &WallFields,
    model: &ThermoviscousFe<'_>,
    transport: &[TangentialWallConductivity],
    conduction: &mut [f64],
    mixed: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), ThermoviscousTimeError> {
    for (region_index, ((region, layout), material)) in model
        .conductive_walls()
        .iter()
        .zip(&fields.layouts)
        .zip(transport)
        .enumerate()
    {
        poll(cancelled)?;
        let ThermalWallDepth::Slab { thickness_m, .. } = region.wall.depth else {
            unreachable!("admitted finite spatial slab")
        };
        let coefficient = depth_mass_coefficient(
            thickness_m,
            layout.normal_elements,
            material.isotropic_w_m_k,
            model.medium.temperature_k,
        )?;
        // Reuse the original hierarchical P2 surface gradients and quadrature.
        // The adapter preserves the mutable cancellation owner at every face.
        let surface = {
            let gate = std::cell::RefCell::new(&mut *cancelled);
            let stopped = std::cell::Cell::new(false);
            let result = model
                .fe
                .boundary_tangential_stiffness_controlled(&region.faces, &|| {
                    let stop = (*gate.borrow_mut())();
                    stopped.set(stopped.get() || stop);
                    stop
                });
            if stopped.get() {
                return Err(ThermoviscousError::Cancelled.into());
            }
            result.map_err(ThermoviscousError::from)?
        };
        for depth in 0..layout.normal_elements {
            poll(cancelled)?;
            for (i, &source_i) in layout.trace_dofs.iter().enumerate() {
                let (columns, values) = surface.row(source_i);
                for (&source_j, &value) in columns.iter().zip(values) {
                    let Ok(j) = layout.trace_dofs.binary_search(&source_j) else {
                        continue;
                    };
                    for a in 0..2 {
                        let Some(row) = fields.coordinate(model, region_index, depth + a, i) else {
                            continue;
                        };
                        for b in 0..2 {
                            let Some(col) = fields.coordinate(model, region_index, depth + b, j)
                            else {
                                continue;
                            };
                            conduction[row * mixed + col] +=
                                coefficient * value * if a == b { 2.0 } else { 1.0 };
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

// Exact P1 depth mass. The surface owner supplies the ORIGINAL dimensional
// gradient form, so this coefficient is k_t*T0*dz/6, never k_t*T0/dz.
fn depth_mass_coefficient(
    thickness_m: f64,
    elements: usize,
    conductivity: f64,
    temperature_k: f64,
) -> Result<f64, ThermoviscousTimeError> {
    let coefficient = conductivity * temperature_k * (thickness_m / elements as f64) / 6.0;
    if !coefficient.is_finite() || coefficient <= 0.0 {
        return Err(invalid(
            "finite positive tangential depth coefficient required",
        ));
    }
    Ok(coefficient)
}

#[cfg(test)]
mod tests;
