//! Conservative runtime rebinding of material transport for live thermal LBM.
//!
//! This is the safe subset of evolving-property physics that does not require
//! remapping conserved state: viscosity, conductivity and buoyancy may change
//! while density, heat capacity and lattice units remain fixed.

use fs_material::fluid_transport::FluidTransport;

use super::{MaterialLbmError, MaterialScaling, MaterialSlab, forcing};

impl<F: FluidTransport + ?Sized> MaterialScaling<'_, F> {
    /// Update an already-evolving 2-D thermal slab from a new material snapshot.
    ///
    /// Momentum and temperature populations are preserved exactly. Density and
    /// specific heat must remain bitwise unchanged: changing either would
    /// reinterpret conserved mass/enthalpy and therefore requires an explicit
    /// conservative transfer operator rather than an in-place coefficient swap.
    /// The current lattice must still carry the previous plan's uniform
    /// relaxation values, so non-Newtonian or manually edited fields are never
    /// silently overwritten.
    ///
    /// # Errors
    /// Refuses dx/dt changes, density/cp changes, independently modified
    /// relaxation fields, and invalid buoyancy inputs before mutation.
    pub fn rebind_thermal_slab<G: FluidTransport + ?Sized>(
        &self,
        slab: &mut MaterialSlab,
        previous: &MaterialScaling<'_, G>,
        gravity_m_s2: f64,
        expansion_per_k: f64,
    ) -> Result<(), MaterialLbmError> {
        if self.cell_size_m.to_bits() != previous.cell_size_m.to_bits()
            || self.timestep_s.to_bits() != previous.timestep_s.to_bits()
        {
            return Err(MaterialLbmError::Invalid(
                "runtime material rebind cannot change the slab's physical dx/dt mapping",
            ));
        }
        if self.transport.density_kg_m3().to_bits()
            != previous.transport.density_kg_m3().to_bits()
            || self.transport.specific_heat_j_kg_k().to_bits()
                != previous.transport.specific_heat_j_kg_k().to_bits()
        {
            return Err(MaterialLbmError::Invalid(
                "runtime density or heat-capacity changes require conservative mass/enthalpy transfer",
            ));
        }
        if slab.lattice.tau_g.to_bits() != previous.thermal_relaxation_time.to_bits()
            || slab
                .lattice
                .grid
                .tau
                .iter()
                .any(|tau| tau.to_bits() != previous.flow.tau.to_bits())
        {
            return Err(MaterialLbmError::Invalid(
                "runtime material rebind refuses independently modified relaxation fields",
            ));
        }

        let bottom_temperature_k = slab.top_temperature_k + slab.temperature_difference_k;
        let (_, gbeta) = self.slab_boundary(
            bottom_temperature_k,
            slab.top_temperature_k,
            gravity_m_s2,
            expansion_per_k,
        )?;

        slab.lattice.grid.tau.fill(self.flow.tau);
        slab.lattice.tau_g = self.thermal_relaxation_time;
        slab.lattice.gbeta = gbeta;
        slab.force_units = forcing::ForceUnits::from_scale(self);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use fs_material::fluid_transport::FluidTransportState;

    use super::super::{MaterialLbmError, MaterialThermalSlab, plan_scaling_from_material};

    fn slab_input() -> MaterialThermalSlab {
        MaterialThermalSlab {
            nx: 4,
            ny: 8,
            bottom_temperature_k: 310.0,
            top_temperature_k: 300.0,
            gravity_m_s2: 9.81,
            expansion_per_k: 2.0e-4,
            max_cells: 40,
        }
    }

    #[test]
    fn rebind_updates_both_diffusivities_and_buoyancy_without_resetting_fields() {
        let first = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        let second = FluidTransportState::try_new(1000.0, 4000.0, 1.0, 0.002).unwrap();
        let old = plan_scaling_from_material(&first, 1.0e-3, 1.0e-3, 0.01).unwrap();
        let new = plan_scaling_from_material(&second, 1.0e-3, 1.0e-3, 0.01).unwrap();
        let mut slab = old.thermal_slab(slab_input()).unwrap();

        let cell = slab.lattice.grid.idx(1, 4);
        slab.lattice.grid.f[cell][1] += 1.0e-4;
        slab.lattice.grid.f[cell][3] -= 1.0e-4;
        slab.lattice.gpop[cell][0] += 2.0e-3;
        let flow_before = slab.lattice.grid.f.clone();
        let heat_before = slab.lattice.gpop.clone();

        new.rebind_thermal_slab(&mut slab, &old, 9.81, 2.0e-4).unwrap();

        assert_eq!(slab.lattice.grid.f, flow_before);
        assert_eq!(slab.lattice.gpop, heat_before);
        assert!(slab.lattice.grid.tau.iter().all(|&tau| tau == new.flow().tau));
        assert_eq!(slab.lattice.tau_g, new.thermal_relaxation_time());
        assert_eq!(
            slab.lattice.gbeta,
            new.acceleration_lattice(9.81 * 2.0e-4 * 10.0).unwrap()
        );
    }

    #[test]
    fn density_or_heat_capacity_change_refuses_without_mutation() {
        let first = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        let old = plan_scaling_from_material(&first, 1.0e-3, 1.0e-3, 0.01).unwrap();
        let mut slab = old.thermal_slab(slab_input()).unwrap();
        let tau_before = slab.lattice.grid.tau.clone();
        let tau_g_before = slab.lattice.tau_g;
        let gbeta_before = slab.lattice.gbeta;

        for changed in [
            FluidTransportState::try_new(900.0, 4000.0, 0.5, 0.001).unwrap(),
            FluidTransportState::try_new(1000.0, 4200.0, 0.5, 0.001).unwrap(),
        ] {
            let next = plan_scaling_from_material(&changed, 1.0e-3, 1.0e-3, 0.01).unwrap();
            let err = next
                .rebind_thermal_slab(&mut slab, &old, 9.81, 2.0e-4)
                .unwrap_err();
            assert!(matches!(
                err,
                MaterialLbmError::Invalid(
                    "runtime density or heat-capacity changes require conservative mass/enthalpy transfer"
                )
            ));
            assert_eq!(slab.lattice.grid.tau, tau_before);
            assert_eq!(slab.lattice.tau_g, tau_g_before);
            assert_eq!(slab.lattice.gbeta, gbeta_before);
        }
    }

    #[test]
    fn rebind_refuses_independent_relaxation_or_unit_changes() {
        let first = FluidTransportState::try_new(1000.0, 4000.0, 0.5, 0.001).unwrap();
        let second = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.0012).unwrap();
        let old = plan_scaling_from_material(&first, 1.0e-3, 1.0e-3, 0.01).unwrap();
        let new = plan_scaling_from_material(&second, 1.0e-3, 1.0e-3, 0.01).unwrap();
        let other_units = plan_scaling_from_material(&second, 2.0e-3, 1.0e-3, 0.01).unwrap();
        let mut slab = old.thermal_slab(slab_input()).unwrap();

        let cell = slab.lattice.grid.idx(0, 3);
        slab.lattice.grid.tau[cell] += 0.01;
        let tau_before = slab.lattice.grid.tau.clone();
        assert!(new.rebind_thermal_slab(&mut slab, &old, 9.81, 2.0e-4).is_err());
        assert_eq!(slab.lattice.grid.tau, tau_before);

        slab.lattice.grid.tau.fill(old.flow().tau);
        assert!(other_units
            .rebind_thermal_slab(&mut slab, &old, 9.81, 2.0e-4)
            .is_err());
    }
}
