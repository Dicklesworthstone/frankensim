//! Explicit material-state conductivity ingress for the existing energy kernel.

use super::{
    AdvectiveEnthalpyField, ChtError, EnergyConfig, EnergySolution, FlowField, FluidProperties,
    PseudoStep, SolidMaterial, ThermalSetup, VoxelDomain, solve_energy_field,
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
        self.solve_with_fields(
            domain,
            fluid,
            solids,
            flow,
            setup,
            molecular_w_m_k,
            None,
            None,
            gate,
        )
    }

    /// Native transient storage and radiation reuse the same admitted material
    /// field assembly, without translating storage into an external heat sink.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_step_with_conductivity(
        &self,
        domain: &VoxelDomain,
        fluid: &FluidProperties,
        solids: &[SolidMaterial],
        flow: &FlowField,
        setup: &ThermalSetup,
        molecular_w_m_k: &[[f64; 3]],
        step: &PseudoStep<'_>,
        gate: &CancelGate,
    ) -> Result<EnergySolution, ChtError> {
        self.solve_with_fields(
            domain,
            fluid,
            solids,
            flow,
            setup,
            molecular_w_m_k,
            None,
            Some(step),
            gate,
        )
    }

    /// The same native storage diagonal with source-caloric face transport.
    /// Nonlinear storage and caloric convergence remain the caller's owner.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve_step_with_enthalpy(
        &self,
        domain: &VoxelDomain,
        fluid: &FluidProperties,
        solids: &[SolidMaterial],
        flow: &FlowField,
        setup: &ThermalSetup,
        molecular_w_m_k: &[[f64; 3]],
        enthalpy: AdvectiveEnthalpyField<'_>,
        step: &PseudoStep<'_>,
        gate: &CancelGate,
    ) -> Result<EnergySolution, ChtError> {
        self.solve_with_fields(domain, fluid, solids, flow, setup, molecular_w_m_k,
            Some(enthalpy), Some(step), gate)
    }

    /// Steady conservative transport of a frozen source-caloric linearization.
    ///
    /// Each interior face carries one donor's `rho * Q * h(T)` on both sides;
    /// the reference Cp in `fluid` is not substituted for that enthalpy. The
    /// same original molecular/contact/boundary heat operator is used. All
    /// caloric entries and inlet/backflow values must share one datum.
    ///
    /// The caller owns nonlinear material iteration and must check the true
    /// chart's local face flux at the returned temperatures. This result audits
    /// the supplied affine law, not convergence of the nonlinear chart.
    ///
    /// # Errors
    /// Requires first-order upwind, positive finite caloric tangents and local
    /// volume conservation, in addition to `solve_with_conductivity` admission.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_enthalpy(
        &self,
        domain: &VoxelDomain,
        fluid: &FluidProperties,
        solids: &[SolidMaterial],
        flow: &FlowField,
        setup: &ThermalSetup,
        molecular_w_m_k: &[[f64; 3]],
        enthalpy: AdvectiveEnthalpyField<'_>,
        gate: &CancelGate,
    ) -> Result<EnergySolution, ChtError> {
        self.solve_with_fields(domain, fluid, solids, flow, setup, molecular_w_m_k,
            Some(enthalpy), None, gate)
    }

    #[allow(clippy::too_many_arguments)]
    fn solve_with_fields(
        &self,
        domain: &VoxelDomain,
        fluid: &FluidProperties,
        solids: &[SolidMaterial],
        flow: &FlowField,
        setup: &ThermalSetup,
        molecular_w_m_k: &[[f64; 3]],
        enthalpy: Option<AdvectiveEnthalpyField<'_>>,
        step: Option<&PseudoStep<'_>>,
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
            domain, fluid, solids, flow, setup, self, step, None,
            Some(molecular_w_m_k), enthalpy, gate,
        )?;
        super::super::poll(gate)?;
        Ok(solution)
    }

    /// One backward-Euler step with explicit material conductivity and cell
    /// heat capacity [J/K]. `previous_k` is the accepted old temperature, not
    /// an iterate. Capacities are frozen during this linear solve; a nonlinear
    /// owner may supply secants `m * delta(h) / delta(T)` and iterate to exact
    /// storage on its admitted single-phase enthalpy chart.
    ///
    /// Storage `C/dt * (T_new - T_old)` uses the same physical storage rows
    /// as the native transient solver. Compact aliases must have zero cell
    /// capacity; each component instead supplies one positive authored
    /// `heat_capacity_j_k`. Its aliases must share one old temperature.
    /// No second flux discretization or external storage reservoir is created.
    /// Physical storage also anchors an insulated domain.
    ///
    /// The returned `EnergySolution` reports physical sources, reservoirs and
    /// sinks, excluding storage. The accompanying record separately reports
    /// `stored + physical_outflow - source`, including any prescribed-cell
    /// temperature jump supplied by its reservoir. Its local `time_s` is
    /// `time_step_s`; a march supplies absolute time. Nonlinear enthalpy
    /// consumers must additionally audit actual `m * delta(h)`; the frozen
    /// linearization's `C * delta(T)` is not their final storage evidence.
    ///
    /// # Errors
    /// Refuses wrong field lengths, nonpositive ordinary/junction capacities,
    /// nonzero capacity on collapsed aliases, invalid durations/temperatures,
    /// inconsistent junction aliases, unrepresentable storage, or existing
    /// solver/cancellation errors. Fixed-temperature cells, domain-face
    /// temperatures, authored compact capacities and physical sinks are supported.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub fn step_with_conductivity(
        &self,
        domain: &VoxelDomain,
        fluid: &FluidProperties,
        solids: &[SolidMaterial],
        flow: &FlowField,
        setup: &ThermalSetup,
        molecular_w_m_k: &[[f64; 3]],
        capacity_j_k: &[f64],
        previous_k: &[f64],
        time_step_s: f64,
        gate: &CancelGate,
    ) -> Result<(EnergySolution, crate::conjugate::TransientRecord), ChtError> {
        self.step_with_fields(domain, fluid, solids, flow, setup, molecular_w_m_k,
            None, capacity_j_k, previous_k, time_step_s, gate)
    }

    /// Backward-Euler material storage and conservative source-enthalpy
    /// advection on the SAME finite-volume heat operator. Conductivity, contacts,
    /// thermal boundaries and physical sinks retain their ordinary owners.
    ///
    /// The advective field is a frozen caloric linearization: the caller owns
    /// nonlinear iteration and must check the true chart's local face fluxes
    /// and `m * delta(h)` at the accepted state. Interior transport uses the
    /// donor's complete affine enthalpy, not local `Cp * T`. Boundary values use
    /// that same chart/datum at inlet or backflow temperature. The returned
    /// report audits this frozen law; it is not nonlinear closure evidence.
    ///
    /// # Errors
    /// In addition to [`Self::step_with_conductivity`] refusals, requires
    /// [`super::ConvectionScheme::Upwind`], constant density, locally conservative
    /// volume flow, a finite positive-Cp tangent on every fluid cell and exact
    /// supplied enthalpy on each open boundary. No state is mutated on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_enthalpy(
        &self, domain: &VoxelDomain, fluid: &FluidProperties, solids: &[SolidMaterial],
        flow: &FlowField, setup: &ThermalSetup, molecular_w_m_k: &[[f64; 3]],
        enthalpy: AdvectiveEnthalpyField<'_>, capacity_j_k: &[f64], previous_k: &[f64],
        time_step_s: f64, gate: &CancelGate,
    ) -> Result<(EnergySolution, crate::conjugate::TransientRecord), ChtError> {
        self.step_with_fields(domain, fluid, solids, flow, setup, molecular_w_m_k,
            Some(enthalpy), capacity_j_k, previous_k, time_step_s, gate)
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn step_with_fields(
        &self, domain: &VoxelDomain, fluid: &FluidProperties, solids: &[SolidMaterial],
        flow: &FlowField, setup: &ThermalSetup, molecular_w_m_k: &[[f64; 3]],
        enthalpy: Option<AdvectiveEnthalpyField<'_>>, capacity_j_k: &[f64], previous_k: &[f64],
        time_step_s: f64, gate: &CancelGate,
    ) -> Result<(EnergySolution, crate::conjugate::TransientRecord), ChtError> {
        use super::super::transient::{junction_capacity, storage_record};
        let invalid = |reason: &str| ChtError::InvalidInput {
            field: "thermal.material_storage", reason: reason.into(),
        };
        super::super::poll(gate)?;
        let cells = domain.cell_count();
        if capacity_j_k.len() != cells || previous_k.len() != cells {
            return Err(invalid("one capacity and previous temperature per cell are required"));
        }
        if !time_step_s.is_finite() || time_step_s <= 0.0 {
            return Err(invalid("time step must be finite and positive"));
        }
        let owners = super::compact_owners(domain, &setup.compact_components, gate)?;
        let junctions = junction_capacity(domain, setup, &owners, previous_k, gate)?;
        let mut capacity = Vec::with_capacity(cells + junctions.len());
        let mut storage = Vec::with_capacity(cells + junctions.len());
        for (cell, (&value, &temperature)) in capacity_j_k.iter().zip(previous_k).enumerate() {
            if cell.is_multiple_of(4096) { super::super::poll(gate)?; }
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(invalid("previous temperatures must be finite positive kelvin"));
            }
            if owners.get(cell).is_some_and(|&owner| owner != usize::MAX) {
                if value != 0.0 {
                    return Err(invalid("compact aliases require zero cell capacity; the junction owns storage"));
                }
                capacity.push(0.0);
                storage.push(0.0);
                continue;
            }
            let rate = value / time_step_s;
            if !value.is_finite() || value <= 0.0 || !rate.is_finite() || rate <= 0.0
                || !(rate * temperature).is_finite()
            {
                return Err(invalid("ordinary capacity and storage rate must be representable and positive"));
            }
            capacity.push(value);
            storage.push(rate);
        }
        for (part, value) in setup.compact_components.iter().zip(junctions) {
            let old = previous_k[domain.index(part.lo[0], part.lo[1], part.lo[2])];
            let rate = value / time_step_s;
            if !rate.is_finite() || rate <= 0.0 || !(rate * old).is_finite() {
                return Err(invalid("junction storage rate must be representable and positive"));
            }
            capacity.push(value);
            storage.push(rate);
        }
        let step = PseudoStep {
            coefficient_w_k: &storage[..cells],
            previous: previous_k,
            junction_coefficient_w_k: &storage[cells..],
        };
        let solution = self.solve_with_fields(
            domain, fluid, solids, flow, setup, molecular_w_m_k, enthalpy, Some(&step), gate,
        )?;
        for (cell, &temperature) in solution.temperature.iter().enumerate() {
            if cell.is_multiple_of(4096) { super::super::poll(gate)?; }
            if !temperature.is_finite() || temperature <= 0.0 {
                return Err(invalid("solved temperature is not finite positive kelvin"));
            }
        }
        let solids_present = (0..cells).any(|cell| !domain.is_fluid(cell));
        let record = storage_record(
            domain, setup, &capacity, previous_k, &solution,
            (time_step_s, time_step_s, solids_present), (0.0, 0.0), gate,
        )?;
        super::super::poll(gate)?;
        Ok((solution, record))
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

    #[test]
    fn g1_storage_heats_an_insulated_cell_without_an_artificial_anchor() {
        let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
        setup.power_w = vec![6.0];
        let before = setup.clone();
        for capacity in [2.0, 4.0] {
            let (solution, record) = EnergyConfig::default().step_with_conductivity(
                &domain, &fluid(), &[SolidMaterial::new("synthetic", 1.0)],
                &FlowField::quiescent(&domain), &setup, &[[1.0; 3]],
                &[capacity], &[300.0], 0.25, &CancelGate::new(),
            ).unwrap();
            assert!((solution.temperature[0] - 300.0 - 1.5 / capacity).abs() < 1e-10);
            assert!((record.stored_energy_change_j - 1.5).abs() < 1e-10);
            assert_eq!(record.boundary_outflow_j, 0.0);
            assert!(record.closure_j.abs() < 1e-10);
        }
        assert_eq!(setup, before);
    }

    #[test]
    fn g1_g3_storage_cooling_is_backward_euler_and_first_order() {
        use crate::conjugate::CellSink;
        let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
        setup.cell_sinks.push(CellSink { cell: 0, conductance_w_k: 2.0, temperature_k: 300.0 });
        let mut errors = Vec::new();
        for steps in [8, 16, 32] {
            let dt = 1.0 / f64::from(steps);
            let mut temperature = 310.0;
            for _ in 0..steps {
                let (solution, record) = EnergyConfig::default().step_with_conductivity(
                    &domain, &fluid(), &[SolidMaterial::new("synthetic", 1.0)],
                    &FlowField::quiescent(&domain), &setup, &[[1.0; 3]],
                    &[2.0], &[temperature], dt, &CancelGate::new(),
                ).unwrap();
                let exact_step = 300.0 + (temperature - 300.0) / (1.0 + dt);
                assert!((solution.temperature[0] - exact_step).abs() < 1e-9);
                assert!(record.closure_j.abs() < 1e-8);
                assert!(record.boundary_outflow_j > 0.0);
                temperature = solution.temperature[0];
            }
            errors.push((temperature - (300.0 + 10.0 / std::f64::consts::E)).abs());
        }
        assert!(errors[0] / errors[1] > 1.8 && errors[1] / errors[2] > 1.8);
    }

    #[test]
    fn g1_g3_contact_and_spatial_capacity_share_the_existing_energy_operator() {
        let domain = VoxelDomain::from_fn(2, 1, 1, 1.0, |p| Voxel::Solid(u16::from(p[0] > 1.0))).unwrap();
        let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
        setup.contacts.push(ContactResistance { materials: (0, 1), resistance_m2_k_w: 0.5 });
        let (solution, record) = EnergyConfig::default().step_with_conductivity(
            &domain, &fluid(), &[SolidMaterial::new("a", 2.0), SolidMaterial::new("b", 4.0)],
            &FlowField::quiescent(&domain), &setup, &[[2.0; 3], [4.0; 3]],
            &[2.0, 4.0], &[310.0, 300.0], 0.5, &CancelGate::new(),
        ).unwrap();
        let g = 1.0 / (0.5 / 2.0 + 0.5 + 0.5 / 4.0);
        let drop = 10.0 / (1.0 + 0.5 * g * (0.5 + 0.25));
        assert!((solution.temperature[0] - solution.temperature[1] - drop).abs() < 1e-9);
        assert!((2.0 * solution.temperature[0] + 4.0 * solution.temperature[1] - 1820.0).abs() < 1e-8);
        assert!(record.closure_j.abs() < 1e-8);
    }

    #[test]
    fn g0_storage_rejects_invalid_capacity_shapes_and_cancelled_calls() {
        let domain = VoxelDomain::from_fn(1, 1, 1, 1.0, |_| Voxel::Solid(0)).unwrap();
        let call = |capacity: &[f64], dt: f64, gate: &CancelGate| {
            EnergyConfig::default().step_with_conductivity(
                &domain, &fluid(), &[SolidMaterial::new("synthetic", 1.0)],
                &FlowField::quiescent(&domain), &ThermalSetup::new([ThermalFace::Adiabatic; 6]),
                &[[1.0; 3]], capacity, &[300.0], dt, gate,
            )
        };
        for capacity in [vec![], vec![1.0; 2], vec![0.0], vec![-1.0], vec![f64::NAN]] {
            assert!(call(&capacity, 1.0, &CancelGate::new()).is_err());
        }
        assert!(call(&[1.0], 0.0, &CancelGate::new()).is_err());
        let gate = CancelGate::new();
        gate.request();
        assert!(matches!(call(&[1.0], 1.0, &gate), Err(ChtError::Cancelled)));
    }
}
