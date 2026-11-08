//! Smooth-wall SGS traction and optional heat on the original species/energy owner.

use super::*;
use fs_lbm::conjugate::FluidProperties;

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Add explicit equilibrium smooth-wall stress to one already configured
    /// adiabatic no-slip wall. Original molecular source coefficients and
    /// physical wall velocity are resampled on every stage and retry.
    pub fn with_equilibrium_sgs_wall(
        mut self,
        side: BoundarySide,
        config: EquilibriumWallConfig,
    ) -> Result<Self, DissipationError> {
        let index = side_index(side);
        config
            .admit_wall(self.walls[index])
            .map_err(DissipationError::Invalid)?;
        let sgs = self.sgs.as_mut().ok_or(DissipationError::Invalid(
            "equilibrium wall stress requires an explicit SGS operator",
        ))?;
        sgs.wall_models[index] = Some(config.into());
        Ok(self)
    }

    /// Add explicit equilibrium stress and thermal resistance to an already
    /// configured isothermal no-slip wall. The actual current composition and
    /// source heat capacities accompany the original face transport query;
    /// species remain impermeable and no molecular source query is repeated.
    pub fn with_equilibrium_sgs_thermal_wall(
        mut self,
        side: BoundarySide,
        config: EquilibriumThermalWallConfig,
    ) -> Result<Self, DissipationError> {
        let index = side_index(side);
        let model = EquilibriumWallModel::from(config);
        model.admit_wall(self.walls[index]).map_err(DissipationError::Invalid)?;
        let sgs = self.sgs.as_mut().ok_or(DissipationError::Invalid(
            "equilibrium thermal wall requires an explicit SGS operator",
        ))?;
        sgs.wall_models[index] = Some(model);
        Ok(self)
    }

    /// Exact wall controls, independent of local/dynamic momentum or scalars.
    pub fn equilibrium_sgs_wall(&self, side: BoundarySide) -> Option<EquilibriumWallConfig> {
        Some(self.sgs.as_ref()?.wall_models[side_index(side)]?.momentum)
    }

    pub fn equilibrium_sgs_thermal_wall(&self, side: BoundarySide) -> Option<EquilibriumThermalWallConfig> {
        self.sgs.as_ref()?.wall_models[side_index(side)]?.thermal_config()
    }
}

impl SubgridTransport {
    pub(in crate::compressible::spatial::mixture::molecular) fn has_equilibrium_wall(
        &self,
    ) -> bool {
        self.wall_models.iter().any(Option::is_some)
    }

    #[allow(clippy::too_many_arguments)]
    pub(in crate::compressible::spatial::mixture::molecular) fn add_wall(
        &self,
        gas: &MixtureGas,
        state: &[f64],
        side: BoundarySide,
        wall: MolecularWall3d,
        primitive: Primitive,
        distance: f64,
        mut receipt: TransportFaceReport,
        out: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TransportFaceReport, DissipationError> {
        let Some(config) = self.wall_models[side_index(side)] else {
            return Ok(receipt);
        };
        config.admit_wall(wall).map_err(DissipationError::Invalid)?;
        let molecular = receipt.dynamic_viscosity_pa_s;
        let viscosity = config
            .momentum
            .viscosity(
                wall,
                side.axis(),
                primitive.velocity,
                primitive.density,
                primitive.sound_speed,
                distance,
                molecular,
            )
            .map_err(DissipationError::Invalid)?;
        let MolecularWall3d::NoSlip { velocity, .. } = wall else {
            unreachable!();
        };
        let species = gas.species_count();
        // The source has already supplied impermeability, normal stress and
        // molecular wall work. Add only the missing tangential traction and
        // its work at the SAME physical wall velocity, in normal/tangent order.
        let extra = viscosity.secant - molecular;
        for (component, axis) in tangents(side.axis()).into_iter().enumerate() {
            let momentum =
                extra * side.sign() * (primitive.velocity[axis] - velocity[axis]) / distance;
            out[species + 1 + component] += momentum;
            out[species + 3] += momentum * velocity[axis];
        }
        // Preserve the source bound and add the nonlinear wall Jacobian's
        // excess over its molecular shear contribution. The owner subsequently
        // supplies the half-cell and three-axis inverse-length factors.
        receipt.diffusivity_m2_s += 4.0 * (viscosity.tangent - molecular) / primitive.density;
        if config.heat.is_some() {
            let mut cp = [0.0; MAX_SPECIES];
            for (k, model) in gas.source().models().iter().enumerate() {
                if k % 16 == 0 {
                    poll(cancelled)?;
                }
                if state[k] != 0.0 {
                    cp[k] = (state[k] / primitive.density) / model.molar_mass().value()
                        * model
                            .evaluate_properties(Temperature::new(receipt.enthalpy_temperature_k))
                            .map_err(FlowError::from)?
                            .cp()
                            .value();
                }
            }
            let conductivity = config.conductivity(
                wall,
                side.axis(),
                primitive.velocity,
                distance,
                &FluidProperties {
                    density_kg_m3: primitive.density,
                    specific_heat_j_kg_k: sum(cp[..species].iter().copied()),
                    conductivity_w_m_k: receipt.thermal_conductivity_w_m_k,
                    kinematic_viscosity_m2_s: molecular / primitive.density,
                },
                cancelled,
            )?;
            // The impermeable source receipt already carries k/(rho cv) at
            // the actual cell state. Scale that heat-only bound, preserving
            // source mechanical admission and avoiding another caloric solve.
            let ratio = conductivity / receipt.thermal_conductivity_w_m_k;
            let thermal_scale = ratio * receipt.nonmechanical_diffusivity_m2_s;
            if !positive(ratio) || !positive(thermal_scale) {
                return Err(DissipationError::Invalid(
                    "unrepresentable equilibrium SGS thermal wall scale",
                ));
            }
            let MolecularWall3d::NoSlip {
                thermal: WallThermalCondition::Isothermal(temperature),
                ..
            } = wall else {
                unreachable!();
            };
            let heat = -(conductivity - receipt.thermal_conductivity_w_m_k) * side.sign()
                * (temperature.value() - primitive.temperature) / distance;
            out[species + 3] += heat;
            receipt.reduced_heat_flux_w_m2 += heat;
            receipt.nonmechanical_energy_flux_w_m2 += heat;
            receipt.nonmechanical_diffusivity_m2_s = thermal_scale;
            receipt.diffusivity_m2_s = receipt.diffusivity_m2_s.max(thermal_scale);
        }
        if !out.iter().all(|v| v.is_finite())
            || !receipt.diffusivity_m2_s.is_finite()
            || !receipt.reduced_heat_flux_w_m2.is_finite()
            || !receipt.nonmechanical_energy_flux_w_m2.is_finite()
        {
            return Err(DissipationError::Invalid(
                "unrepresentable equilibrium SGS wall flux",
            ));
        }
        Ok(receipt)
    }
}
