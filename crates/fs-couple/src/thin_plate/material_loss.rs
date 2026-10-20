//! Source-bound regional Kelvin-Voigt bending, projected on the existing modes.
//! The DKT element owner supplies every strain-energy contribution. No averaged
//! viscosity, damping of conservative prestress, or second modal solve.
//!
//! This is the same diagonal, weak-loss modal approximation used by the uniform
//! plate. Heterogeneous loss is generally nonproportional: off-diagonal modal
//! damping and strongly damped complex modes are not represented by this image.

use super::*;

impl ResolvedPlateChart {
    /// Select isotropic Kelvin-Voigt bending in EVERY region, from its retained
    /// material state. Viscosity, density, E and nu use the exact same admission
    /// and intersected frequency-domain policy as a uniform material specimen.
    /// Zero viscosity remains an explicit law with a checked applicability band.
    /// Replacing a region re-resolves this selection; it cannot silently drop it.
    ///
    /// The loss tensor is proportional to each region's elastic bending tensor.
    /// This does not supply anisotropic viscosity or interface/joint damping.
    /// Modal realization requires zero authored damping to avoid double counting.
    pub fn with_kelvin_voigt_bending_loss(mut self) -> Result<Self, AcousticRealizeError> {
        let mut viscosities = Vec::with_capacity(self.regions.len());
        for region in &self.regions {
            if region.input.model != PlateMaterialModel::Isotropic {
                return Err(plate_assignment_error(Some(&region.input.region.name), None,
                    "regional bending viscosity requires isotropic material elasticity"));
            }
            let law = resolve_viscosity(&region.input.material)?;
            let young = plate_property(&region.input.material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?;
            viscous_time(law, young)?;
            viscosities.push(law);
        }
        self.bending_viscosity = Some(viscosities);
        Ok(self)
    }

    /// Selected loss and source identities in the original region order.
    /// `None` retains the legacy authored-damping path. Geometry/specimen identity
    /// alone is not a complete scenario identity or a claim that loss is selected.
    #[must_use]
    pub fn bending_viscosity(&self) -> Option<&[IsotropicPlateBendingViscosity]> {
        self.bending_viscosity.as_deref()
    }

    /// Project selected material loss on an existing mode of this exact chart.
    /// The caller owns the assembly/eigensolve and must supply its matching DOF
    /// map and mode. This keeps valve, body and force ports on the same basis.
    pub(crate) fn modal_bending_damping_ratio(
        &self, model: &PlateModel, phi: &[f64], omega: f64,
    ) -> Result<Option<f64>, AcousticRealizeError> {
        self.bending_viscosity().map(|laws| {
            RegionalBendingLoss::new(self, laws, model)?.modal_zeta(model, phi, omega)
        }).transpose()
    }
}

/// Derive body modes, losses and observation ports from one regional specimen.
///
/// With source loss selected, each element contributes `(eta/E) * K_bending`;
/// its contribution is weighted by its ACTUAL modal curvature energy. Neither
/// area nor mass averages decide decay. The denominator includes conservative
/// geometric stiffness, so tension dilutes loss rather than acquiring it.
/// All retained modes must lie inside every region's source-law applicability.
///
/// Single- and multi-port source-bound bodies share the same reduction path.
/// Without selected source loss the original authored chart convention remains.
/// Source bindings remain with the specimen, not the numeric output. This
/// diagonal weak-loss image does not solve full nonproportional damping.
///
/// # Errors
/// Invalid law/controls, nonzero authored damping alongside material loss,
/// source-domain refusal, geometry/modal failure or unrepresentable projection.
pub fn certified_material_chart_radiators(
    specimen: &ResolvedPlateChart, options: &PlateChartRadiation,
) -> Result<Vec<CompactBody>, AcousticRealizeError> {
    ports::reduce(specimen, options, &[], 0).map(|(radiators, _)| radiators)
}

mod ports;

fn refuse(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}

// One material-property admission path for uniform and regional plates.
pub(super) fn resolve_viscosity(material: &ResolvedMaterialStatePoint)
    -> Result<IsotropicPlateBendingViscosity, AcousticRealizeError>
{
    let property = material
        .property(KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY)
        .ok_or_else(|| refuse("material plate loss needs kelvin_voigt_bending_viscosity"))?;
    let viscosity_pa_s = plate_property(
        material,
        KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
        DynViscosity::DIMS,
    )?;
    let mut band = property
        .angular_frequency_band_rad_s()
        .map_err(|_| refuse("material plate viscosity has invalid frequency-axis semantics or bounds"))?
        .ok_or_else(|| {
            refuse("material plate viscosity needs an explicit omega band in rad/s or a typed frequency band")
        })?;
    for key in [
        DENSITY_PROPERTY,
        YOUNG_MODULUS_PROPERTY,
        POISSON_RATIO_PROPERTY,
        KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY,
    ] {
        let coefficient = material
            .property(key)
            .ok_or_else(|| refuse("material plate loss needs density, E, nu and viscosity"))?;
        if !coefficient.is_constant_at_fixed_temperature() {
            return Err(refuse(
                "Kelvin-Voigt plate coefficients must be validity-wide scalar constants or temperature-only curves at fixed positive absolute T",
            ));
        }
        if let Some((lo, hi)) = coefficient.angular_frequency_band_rad_s().map_err(|_| {
            refuse("material plate coefficient has invalid frequency-axis semantics or bounds")
        })? {
            band = (band.0.max(lo), band.1.min(hi));
        }
    }
    Ok(IsotropicPlateBendingViscosity {
        viscosity_pa_s,
        omega_band_rad_s: band,
        material_state_identity: Some(material.identity()),
    })
}

// Same numeric refusal as the original uniform PlateBendingLaw.
pub(super) fn viscous_time(viscous: IsotropicPlateBendingViscosity, young: f64)
    -> Result<f64, AcousticRealizeError>
{
    let (lo, hi) = viscous.omega_band_rad_s;
    let viscous_time_s = viscous.viscosity_pa_s / young;
    if !viscous.viscosity_pa_s.is_finite()
        || viscous.viscosity_pa_s < 0.0
        || !viscous_time_s.is_finite()
        || (viscous.viscosity_pa_s > 0.0 && viscous_time_s == 0.0)
        || !lo.is_finite()
        || !hi.is_finite()
        || lo < 0.0
        || hi < lo
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "plate viscosity and angular-frequency applicability must be finite and nonnegative with an ordered band",
        });
    }
    Ok(viscous_time_s)
}

struct BendingElement {
    dofs: [Option<usize>; 9],
    stiffness: [f64; 81],
    region: usize,
}
struct RegionalBendingLoss {
    laws: Vec<PlateBendingLaw>,
    elements: Vec<BendingElement>,
}
impl RegionalBendingLoss {
    fn new(specimen: &ResolvedPlateChart, viscosities: &[IsotropicPlateBendingViscosity],
        model: &PlateModel) -> Result<Self, AcousticRealizeError>
    {
        let mut laws = Vec::with_capacity(viscosities.len());
        for (region, &viscous) in specimen.regions.iter().zip(viscosities) {
            let young = plate_property(&region.input.material, YOUNG_MODULUS_PROPERTY, Pressure::DIMS)?;
            laws.push(PlateBendingLaw { thermal: None, viscous: Some(viscous),
                viscous_time_s: viscous_time(viscous, young)?, thickness: region.section.thickness });
        }
        let mesh = &specimen.chart.mesh;
        let mut elements = Vec::with_capacity(mesh.tris.len());
        for (element, tri) in mesh.tris.iter().enumerate() {
            let region = specimen.region_for_element[element];
            let x = tri.map(|node| mesh.nodes[node].0);
            let y = tri.map(|node| mesh.nodes[node].1);
            // Original DKT element, same section, geometry and boundary map as
            // the modal pencil. No geometric stiffness enters viscous storage.
            let (stiffness, _) = fs_plate::dkt_stiffness(&x, &y,
                &specimen.regions[region].section.d, element).map_err(map_plate)?;
            let dofs = std::array::from_fn(|i| model.dof_map[3 * tri[i / 3] + i % 3]);
            elements.push(BendingElement { dofs, stiffness, region });
        }
        Ok(Self { laws, elements })
    }
}
impl ModalBendingLoss for RegionalBendingLoss {
    fn modal_zeta(&self, total: &PlateModel, phi: &[f64], omega: f64)
        -> Result<f64, AcousticRealizeError>
    {
        // Check ALL source bands, including silent/zero-viscosity regions.
        let zetas: Vec<f64> = self.laws.iter().map(|law| law.zeta(omega)).collect::<Result<_, _>>()?;
        let mut k_phi = vec![0.0; phi.len()];
        total.k.spmv(phi, &mut k_phi);
        let total_energy = phi.iter().zip(k_phi).map(|(p, k)| p * k).sum();
        let mut weighted_bending = 0.0;
        // Fixed element order, not region-declaration order or hash order.
        for element in &self.elements {
            let q: [f64; 9] = element.dofs.map(|dof| dof.map_or(0.0, |i| phi[i]));
            let mut energy = 0.0;
            for r in 0..9 {
                for c in 0..9 {
                    energy += q[r] * element.stiffness[9 * r + c] * q[c];
                }
            }
            weighted_bending += zetas[element.region] * energy;
        }
        bending_energy_fraction(weighted_bending, total_energy)
    }
}

#[cfg(test)]
mod tests;
