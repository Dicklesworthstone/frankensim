//! Shared adaptive Green rows with the original equivalent-fluid momentum.
use super::*;
use crate::helmholtz::molecular::MolecularNearFieldOptions;
use crate::near_field::{Geometry, VectorObservation, complex::Prepared};

/// Local heating from the same coherent interior field. These are macroscopic
/// rates per total body volume (gross volume for a porous body), not intrinsic
/// pore-fluid velocity or a gas/solid heat
/// partition. Quadrature discrepancies are retained on the paired fields;
/// they do not certify the spatial solve or relative accuracy near a null.
pub struct PorousAbsorptionObservation {
    pub fields: VectorObservation,
    pub viscous_loss_w_m3: Vec<f64>,
    /// Bulk-compliance loss, including pore heat diffusion and any coupled
    /// molecular reservoirs in the original material response.
    pub compression_loss_w_m3: Vec<f64>,
}

impl PorousTransmissionField<'_, '_> {
    /// Outgoing scattered p/u from the original retained triangles. The same
    /// adaptive row bank integrates pressure and all gradient components.
    pub fn scattered_pressure_velocity_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        self.observe_adaptive(
            false,
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            options,
            gate,
        )
    }

    /// Add the original analytic illumination after the outgoing field. Local
    /// intensity is Re(p*conj(u))/2; adding separate source powers would omit
    /// interference. Original source and panel work budgets remain distinct.
    pub fn total_pressure_velocity_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        molecular::admit_products(
            self.incident.source_count(),
            points.len(),
            self.maximum_source_point_products,
        )?;
        let mut result = self.scattered_pressure_velocity_adaptive_at(
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            options,
            gate,
        )?;
        for ((&point, pressure), velocity) in points
            .iter()
            .zip(&mut result.scalar.pressure)
            .zip(&mut result.particle_velocity_m_s)
        {
            poll(gate)?;
            let source = self.incident.sample_wave(
                self.interface.k_o,
                self.interface.exterior.density_kg_m3() * self.interface.omega,
                None,
                point,
                &mut || gate.is_requested(),
            )?;
            *pressure = *pressure + source.pressure_pa;
            for (component, incoming) in velocity.iter_mut().zip(source.velocity_m_s) {
                *component = *component + incoming;
            }
            if !finite(*pressure) || velocity.iter().any(|&value| !finite(value)) {
                return Err(bad("nonfinite porous transmission total receiver"));
            }
        }
        poll(gate)?;
        Ok(result)
    }

    /// Interior S_i*q_i-D_i*p and u=grad(p)/(i*omega*rho_eff). Both the boundary
    /// derivative and gradient conversion retain the same complex density.
    /// Superficial velocity is returned; porosity is not applied a second time.
    pub fn interior_pressure_velocity_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        self.observe_adaptive(
            true,
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            options,
            gate,
        )
    }

    /// Direct local dissipation from coherent p/u:
    /// omega/2 * [Im(rho_eff)|u|² + Im(1/K_eff)|p|²].
    /// All three velocity components participate. Neither a boundary-work
    /// residual nor 1-R-T supplies heat. The original row-storage budget
    /// excludes bounded receiver/output buffers, including these two vectors.
    pub fn interior_absorption_adaptive_at(
        &self,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<PorousAbsorptionObservation, HelmholtzError> {
        let fields = self.interior_pressure_velocity_adaptive_at(
            points,
            minimum_clearance_m,
            maximum_panel_point_products,
            options,
            gate,
        )?;
        absorption(
            fields,
            self.interface.density,
            self.interface.bulk,
            self.interface.omega,
            gate,
        )
    }

    fn observe_adaptive(
        &self,
        interior: bool,
        points: &[[f64; 3]],
        minimum_clearance_m: f64,
        maximum_panel_point_products: u64,
        options: MolecularNearFieldOptions,
        gate: &CancelGate,
    ) -> Result<VectorObservation, HelmholtzError> {
        poll(gate)?;
        let owner = self.interface;
        let panels = owner.surface.areas().len();
        molecular::admit_products(panels, points.len(), maximum_panel_point_products)?;
        options.admit(panels, points.len())?;
        let geometry = if interior {
            Geometry::new_interior(owner.surface, points, minimum_clearance_m)?
        } else {
            Geometry::new(owner.surface, points, minimum_clearance_m)?
        };
        let (wave, density, pressure, velocity) = if interior {
            (
                owner.k_i,
                owner.density,
                &self.total_pressure,
                &self.total_velocity,
            )
        } else {
            (
                owner.k_o,
                C64::from_re(owner.exterior.density_kg_m3()),
                &self.scattered_pressure,
                &self.scattered_velocity,
            )
        };
        let prepared = Prepared::new_equivalent_fluid(
            &geometry,
            wave,
            density.scale(owner.omega),
            options.quadrature,
            &mut || gate.is_requested(),
        )?;
        let mut result = prepared.evaluate(pressure, velocity, &mut || gate.is_requested())?;
        if interior {
            for (pressure, velocity) in result
                .scalar
                .pressure
                .iter_mut()
                .zip(&mut result.particle_velocity_m_s)
            {
                poll(gate)?;
                *pressure = -*pressure;
                for component in velocity {
                    *component = -*component;
                }
            }
        }
        poll(gate)?;
        Ok(result)
    }
}

// Shared direct work contraction for single and mutually coupled interiors.
pub(super) fn absorption(
    fields: VectorObservation,
    density: C64,
    bulk: C64,
    omega: f64,
    gate: &CancelGate,
) -> Result<PorousAbsorptionObservation, HelmholtzError> {
    // Apply the positive loss coefficient BEFORE squaring amplitudes.
    // This retains weak loss when an unweighted |p|² or |u|² would overflow.
    let common = det::sqrt(omega) * core::f64::consts::FRAC_1_SQRT_2;
    let viscous_root = common * det::sqrt(density.im);
    let compression_root = common * det::sqrt(bulk.recip().im);
    let mut viscous_loss_w_m3 = Vec::with_capacity(fields.scalar.pressure.len());
    let mut compression_loss_w_m3 = Vec::with_capacity(fields.scalar.pressure.len());
    for (&pressure, velocity) in fields
        .scalar
        .pressure
        .iter()
        .zip(&fields.particle_velocity_m_s)
    {
        poll(gate)?;
        let viscous = velocity
            .iter()
            .map(|value| value.scale(viscous_root).norm_sq())
            .sum();
        let compression = pressure.scale(compression_root).norm_sq();
        if !f64::is_finite(viscous) || !compression.is_finite() {
            return Err(bad("unrepresentable local porous absorption"));
        }
        viscous_loss_w_m3.push(viscous);
        compression_loss_w_m3.push(compression);
    }
    poll(gate)?;
    Ok(PorousAbsorptionObservation {
        fields,
        viscous_loss_w_m3,
        compression_loss_w_m3,
    })
}
