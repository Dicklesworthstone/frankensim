//! Explicit material-state conductivity ingress for the existing energy kernel.

use super::{
    ChtError, EnergyConfig, EnergySolution, FlowField, FluidProperties,
    SolidMaterial, ThermalSetup, VoxelDomain, solve_energy_field,
};
use fs_exec::CancelGate;

impl EnergyConfig {
    /// Solve with one molecular conductivity `[kx, ky, kz]` in W/(m K) per
    /// cell, using the ordinary finite-volume assembly, boundary conditions,
    /// contact resistances, compact-component links and energy accounting.
    ///
    /// The coefficients are a frozen resolved state. They replace the bulk
    /// conductivities in `fluid`/`solids`; those tables still supply their
    /// validated reference properties and material-index/contact topology.
    /// No voxel is relabelled and no per-cell pseudo-material is created.
    /// Fluid heat-capacity transport stays `fluid.rho * fluid.cp`: this API
    /// does not admit variable-density enthalpy advection. The optional eddy
    /// conductivity in `setup` is added only to fluid cells, exactly once.
    ///
    /// Independent positive grid-aligned principal conductivities are supported.
    /// Off-diagonal tensors and coordinate rotations are not inferred. The caller
    /// owns source receipts, full state coordinates and nonlinear iteration;
    /// this method does not clamp, extrapolate, or average material data. Supplying
    /// raw values is not evidence of source qualification.
    ///
    /// # Errors
    /// Refuses an incomplete/non-positive/non-finite field, simultaneous authored
    /// conductivity tables (two conflicting material owners), ordinary solver
    /// failures, or cancellation. No field is returned after an observed cancel.
    #[allow(clippy::too_many_arguments)] // the existing physical inputs plus resolved coefficients
    pub fn solve_with_conductivity(
        &self,
        domain: &VoxelDomain,
        fluid: &FluidProperties,
        solids: &[SolidMaterial],
        flow: &FlowField,
        setup: &ThermalSetup,
        molecular_w_m_k: &[[f64; 3]],
        gate: &CancelGate,
    ) -> Result<EnergySolution, ChtError> {
        super::super::poll(gate)?;
        if molecular_w_m_k.len() != domain.cell_count() {
            return Err(ChtError::InvalidInput {
                field: "thermal.molecular_conductivity",
                reason: "needs exactly one grid-aligned conductivity triple per cell".into(),
            });
        }
        for (cell, values) in molecular_w_m_k.iter().enumerate() {
            if cell.is_multiple_of(4096) {
                super::super::poll(gate)?;
            }
            if values.iter().any(|k| !k.is_finite() || *k <= 0.0) {
                return Err(ChtError::InvalidInput {
                    field: "thermal.molecular_conductivity",
                    reason: format!("cell {cell}: all conductivities must be finite and positive"),
                });
            }
        }
        if solids.iter().any(|s| !s.conductivity_table.is_empty()) {
            return Err(ChtError::InvalidInput {
                field: "thermal.molecular_conductivity",
                reason: "explicit resolved field cannot also use authored conductivity tables".into(),
            });
        }
        let solution = solve_energy_field(
            domain, fluid, solids, flow, setup, self, None, None,
            Some(molecular_w_m_k), gate,
        )?;
        super::super::poll(gate)?;
        Ok(solution)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conjugate::{ContactResistance, ThermalFace, Voxel, solve_energy};

    fn fluid() -> FluidProperties {
        FluidProperties {
            density_kg_m3: 1.0,
            specific_heat_j_kg_k: 1.0,
            conductivity_w_m_k: 1.0,
            kinematic_viscosity_m2_s: 1.0,
        }
    }

    fn boundary() -> ThermalSetup {
        let mut faces = [ThermalFace::Adiabatic; 6];
        faces[0] = ThermalFace::HeatFlux(4.0);
        faces[1] = ThermalFace::Temperature(300.0);
        ThermalSetup::new(faces)
    }

    #[test]
    fn g1_local_material_field_changes_temperature_without_relabelling_voxels() {
        let domain = VoxelDomain::from_fn(4, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let before = domain.clone();
        let k = [1.0, 2.0, 4.0, 8.0];
        let field = k.map(|x| [x, 3.0 * x, 7.0 * x]);
        let solution = EnergyConfig::default().solve_with_conductivity(
            &domain, &fluid(), &[SolidMaterial::new("reference only", 100.0)],
            &FlowField::quiescent(&domain), &boundary(), &field, &CancelGate::new(),
        ).unwrap();
        for c in 0..4 {
            let resistance = 0.5 / k[c] + k[c + 1..].iter().map(|x| 1.0 / x).sum::<f64>();
            assert!((solution.temperature[c] - (300.0 + 4.0 * resistance)).abs() < 1e-8);
        }
        assert_eq!(solution.conductivity, field.to_vec());
        assert!(solution.report.balance.relative_residual < 1e-9);
        assert_eq!(domain, before);
    }

    #[test]
    fn g5_resolved_and_reference_paths_share_contacts_and_identical_arithmetic() {
        let domain = VoxelDomain::from_fn(4, 1, 1, 1.0, |p| {
            Voxel::Solid(u16::from(p[0] > 2.0))
        }).unwrap();
        let solids = [SolidMaterial::new("left", 2.0), SolidMaterial::new("right", 8.0)];
        let mut setup = boundary();
        setup.contacts.push(ContactResistance { materials: (0, 1), resistance_m2_k_w: 0.7 });
        let flow = FlowField::quiescent(&domain);
        let gate = CancelGate::new();
        let config = EnergyConfig::default();
        let original = solve_energy(&domain, &fluid(), &solids, &flow, &setup, &config, &gate).unwrap();
        let resolved = config.solve_with_conductivity(
            &domain, &fluid(), &solids, &flow, &setup,
            &[[2.0; 3], [2.0; 3], [8.0; 3], [8.0; 3]], &gate,
        ).unwrap();
        assert_eq!(original, resolved);
    }

    #[test]
    fn g0_fluid_eddy_conductivity_is_added_once_to_the_resolved_field() {
        let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
        let mut setup = boundary();
        setup.eddy_conductivity_w_m_k = vec![1.0; 3];
        let solution = EnergyConfig::default().solve_with_conductivity(
            &domain, &fluid(), &[], &FlowField::quiescent(&domain), &setup,
            &[[2.0, 3.0, 4.0]; 3], &CancelGate::new(),
        ).unwrap();
        assert_eq!(solution.conductivity, vec![[3.0, 4.0, 5.0]; 3]);
        assert!(solution.report.balance.relative_residual < 1e-9);
    }

    #[test]
    fn g0_shape_values_and_competing_conductivity_owners_refuse() {
        let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let solids = [SolidMaterial::new("reference", 1.0)];
        let check = |field: &[[f64; 3]], solids: &[SolidMaterial]| {
            EnergyConfig::default().solve_with_conductivity(
                &domain, &fluid(), solids, &FlowField::quiescent(&domain),
                &boundary(), field, &CancelGate::new(),
            )
        };
        for n in [0, 2, 4] {
            assert!(matches!(check(&vec![[1.0; 3]; n], &solids),
                Err(ChtError::InvalidInput { field: "thermal.molecular_conductivity", .. })));
        }
        for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(check(&[[1.0, bad, 1.0]; 3], &solids).is_err());
        }
        let competing = [solids[0].clone().with_conductivity_table(&[(300.0, 1.0), (400.0, 2.0)])];
        assert!(check(&[[1.0; 3]; 3], &competing).is_err());
    }
}
