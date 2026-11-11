//! Inlet superheat and actual Darcy dissipation melt the same original cells.
use super::*;

fn pressure_channel(sign: f64, gauge: f64) -> SimpleConfig {
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Pressure { pressure_pa: gauge + sign * 20.0 };
    faces[1] = FvBoundary::Pressure { pressure_pa: gauge };
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-13;
    config.pressure_tolerance = 1e-12;
    config.momentum_tolerance = 1e-10;
    config
}

fn input<'a>(ids: &'a [Option<CellId>], curve: &'a EquilibriumEnthalpyPhaseCurve) -> FlowEnthalpyInput<'a> {
    let mut input = flow_input(ids, curve);
    // These are thermal enthalpies at the source's reference pressure. The
    // mechanical pressure gauge is not added to either inlet's caloric datum.
    input.inlet_specific_h_j_kg[0] = Some(150.0);
    input.inlet_specific_h_j_kg[1] = Some(150.0);
    input
}

fn beta(h: f64) -> f64 {
    let liquid_fraction = h / 100.0;
    100.0 * (1.0 - liquid_fraction).powi(2) / liquid_fraction.powi(3)
}

#[test]
fn pressure_driven_liquid_superheat_and_physical_drag_melt_cells_without_double_flow_work() {
    let card = card(claims(1.0, 100.0, false));
    let liquid = liquid(&card);
    let curve = continuity::chart(&card);
    let phase = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.1).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 4, _>(|i| Some(CellId(i as u64)));
    let gate = CancelGate::new();
    let initial = super::super::initial(&domain, &curve, 0.0);
    for sign in [-1.0, 1.0] {
        for gauge in [0.0, 1e6] {
            let config = pressure_channel(sign, gauge);
            let base = MaterialFlowFeedback::try_new(
                &initial, &domain, input(&ids, &curve), &config, (), &gate, 20,
                |_, _: &DepositedCell, t| liquid.at_temperature(Temperature::new(t))
                    .map_err(|error| fs_material::MaterialError::State { what: error.to_string() }),
            ).unwrap().with_mushy_mobility(&phase).unwrap();
            assert!(base.clone().with_liquid_carrier().is_err());
            let mut observer = base.with_pressure_work_heating(heating()).unwrap().with_liquid_carrier().unwrap();
            let mut current = initial.clone();
            let mut h = [50.0; 4];
            let inlet_x = if sign > 0.0 { 0 } else { 1 };
            for target in [0.1, 0.2] {
                let drag = h.map(beta);
                // Transverse symmetry gives uniform superficial velocity in
                // each two-cell row. The native series resistance is beta0+
                // beta1, including the physical half cells at the two mouths.
                let speed = 20.0 / (drag[0] + drag[1]);
                let stream_w = 2.0 * speed * (150.0 - 100.0);
                let pressure_w = 20.0 * 2.0 * speed;
                let (next, accepted, audit) = current.advance_observed(&observer, target, thermal(), || Ok(())).unwrap();
                assert_eq!(accepted.momentum_solves(), 1);
                let epoch = accepted.applied().unwrap();
                for x in 0..=2 {
                    for y in 0..2 { near(epoch.flow.field.flux_x(x, y, 0), sign * speed); }
                }
                let power = epoch.viscous_power.as_ref().unwrap();
                near(power.molecular_heat_w, 0.0);
                near(power.darcy_heat_w, pressure_w);
                near(power.heat_w, pressure_w);
                near(power.supplied_pressure_power_w, pressure_w);
                for c in 0..4 {
                    let physical_heat = drag[c] * speed * speed;
                    let stream_heat = if c % 2 == inlet_x { speed * 50.0 } else { 0.0 };
                    h[c] += 0.1 * (physical_heat + stream_heat);
                    near(power.cell_heat_w[c], physical_heat);
                    near(next.cells()[&CellId(c as u64)].specific_h_j_kg, h[c]);
                    assert_eq!(next.temperature_k(CellId(c as u64)).unwrap(), 500.0);
                }
                let work = accepted.viscous_heating_audit().unwrap();
                near(work.generated_heat_j, 0.1 * pressure_w);
                near(work.supplied_pressure_work_j, work.generated_heat_j);
                near(accepted.audit().external_enthalpy_j, 0.1 * stream_w);
                near(accepted.audit().inlet_mass_kg, 0.1 * 2.0 * speed);
                near(accepted.audit().outlet_mass_kg, accepted.audit().inlet_mass_kg);
                near(audit.external_heat_j, 0.1 * (stream_w + pressure_w));
                near(audit.stored_energy_change_j, audit.external_heat_j);
                assert_eq!(next.mass_kg(), initial.mass_kg());
                current = next;
                observer = accepted;
            }
        }
    }
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
#[cfg(feature = "implicit-enthalpy-coupling")]
fn implicit_open_liquid_flow_uses_its_own_melt_dependent_pressure_and_outlet_enthalpy() {
    let card = card(claims(1.0, 100.0, false));
    let liquid = liquid(&card);
    let curve = continuity::chart(&card);
    let phase = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.1).unwrap();
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 4, _>(|i| Some(CellId(i as u64)));
    let config = pressure_channel(1.0, 1e6);
    let gate = CancelGate::new();
    let initial = super::super::initial(&domain, &curve, 0.0);
    let observer = MaterialFlowFeedback::try_new(
        &initial, &domain, input(&ids, &curve), &config, (), &gate, 2048,
        |_, _: &DepositedCell, t| liquid.at_temperature(Temperature::new(t))
            .map_err(|error| fs_material::MaterialError::State { what: error.to_string() }),
    ).unwrap().with_mushy_mobility(&phase).unwrap()
        .with_pressure_work_heating(heating()).unwrap().with_liquid_carrier().unwrap();
    let controls = ImplicitFlowControls {
        max_unknowns: 4, max_newton_iterations: 32, difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-7, total_energy_tolerance_j: 4e-7,
        newton: fs_solver::NewtonKrylovConfig { absolute_tolerance: 2e-8,
            relative_tolerance: 0.0, linear_restart: 4, ..Default::default() },
    };
    let dt = 0.1;
    let (next, report) = observer.step_implicit(&initial, dt, controls, &[0.0; 4], || Ok(())).unwrap();
    let h: [f64; 4] = std::array::from_fn(|c| next.cells()[&CellId(c as u64)].specific_h_j_kg);
    let drag = h.map(beta);
    let speed = 20.0 / (drag[0] + drag[1]);
    near(report.endpoint.flow.field.flux_x(0, 0, 0), speed);
    for c in 0..4 {
        let stream_w = if c % 2 == 0 { speed * 50.0 } else { 0.0 };
        assert!((h[c] - 50.0 - dt * (drag[c] * speed * speed + stream_w)).abs() < 1e-7);
        assert_eq!(next.temperature_k(CellId(c as u64)).unwrap(), 500.0);
    }
    near(report.mechanical.generated_heat_j, dt * 20.0 * 2.0 * speed);
    near(report.mechanical.supplied_pressure_work_j, report.mechanical.generated_heat_j);
    near(report.net_stream_enthalpy_j, dt * 2.0 * speed * 50.0);
    near(report.thermal.external_heat_j, report.mechanical.generated_heat_j + report.net_stream_enthalpy_j);
    assert!(report.combined_residual_j.abs() < 4e-7);
    assert!(observer.applied().is_none());
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
#[cfg(feature = "implicit-enthalpy-coupling")]
fn live_reconstruction_carries_liquid_h_without_transporting_stationary_latent_inventory() {
    use fs_conduction::transient::deposition::advection::reconstruction::{
        AxialCellStencil, LimitedEnthalpyReconstruction,
    };
    use std::sync::Arc;
    let card = card(claims(1.0, 100.0, false));
    let liquid = liquid(&card);
    let curve = continuity::chart(&card);
    let phase = MushyPhaseLaw::try_new(&curve, 1.0, 1.0, 1.0, 0.1).unwrap();
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let ids = std::array::from_fn::<_, 4, _>(|i| Some(CellId(i as u64)));
    let config = pressure_channel(1.0, 1e6);
    let gate = CancelGate::new();
    let initial = state(&domain, &curve);
    let stencils = [
        AxialCellStencil { previous: CellId(0), cell: CellId(1), next: CellId(2) },
        AxialCellStencil { previous: CellId(1), cell: CellId(2), next: CellId(3) },
    ];
    let plan = Arc::new(LimitedEnthalpyReconstruction::try_new(&initial, &stencils, 2, || Ok(())).unwrap());
    let mut binding = input(&ids, &curve);
    binding.inlet_specific_h_j_kg = [Some(100.0); 6];
    let base = MaterialFlowFeedback::try_new(
        &initial, &domain, binding, &config, (), &gate, 2048,
        |_, _: &DepositedCell, t| liquid.at_temperature(Temperature::new(t))
            .map_err(|error| fs_material::MaterialError::State { what: error.to_string() }),
    ).unwrap().with_mushy_mobility(&phase).unwrap().with_pressure_work_heating(heating()).unwrap();
    for carrier_first in [true, false] {
        let observer = if carrier_first {
            base.clone().with_liquid_carrier().unwrap().with_reconstruction(Arc::clone(&plan)).unwrap()
        } else {
            base.clone().with_reconstruction(Arc::clone(&plan)).unwrap().with_liquid_carrier().unwrap()
        };
        let dt = 0.1;
        let old_h = [20.0, 40.0, 60.0, 80.0];
        let drag = old_h.map(beta);
        let q = 20.0 / drag.iter().sum::<f64>();
        let (next, accepted, audit) = initial.advance_observed(&observer, dt, thermal(), || Ok(())).unwrap();
        assert_eq!(accepted.momentum_solves(), 1);
        near(accepted.applied().unwrap().flow.field.flux_x(2, 0, 0), q);
        for c in 0..4 {
            // The reconstructed transported field is flat h_liquid=100 even
            // though the mixture h has a nonzero limited axial slope.
            near(next.cells()[&CellId(c as u64)].specific_h_j_kg, old_h[c] + dt * drag[c] * q * q);
            assert_eq!(next.temperature_k(CellId(c as u64)).unwrap(), 500.0);
        }
        near(accepted.audit().external_enthalpy_j, 0.0);
        near(audit.external_heat_j, dt * 20.0 * q);
        near(audit.stored_energy_change_j, audit.external_heat_j);
        assert!(std::ptr::eq(accepted.reconstruction().unwrap(), plan.as_ref()));
        assert_eq!(accepted.liquid_carrier().unwrap().caloric_identity(), curve.caloric_identity().unwrap());

        let controls = ImplicitFlowControls {
            max_unknowns: 4, max_newton_iterations: 32, difference_step_j: 1e-4,
            cell_energy_tolerance_j: 1e-7, total_energy_tolerance_j: 4e-7,
            newton: fs_solver::NewtonKrylovConfig { absolute_tolerance: 2e-8,
                relative_tolerance: 0.0, linear_restart: 4, ..Default::default() },
        };
        let (implicit, report) = observer.step_implicit(&initial, dt, controls, &[0.0; 4], || Ok(())).unwrap();
        let h: [f64; 4] = std::array::from_fn(|c| implicit.cells()[&CellId(c as u64)].specific_h_j_kg);
        let drag = h.map(beta);
        let q = 20.0 / drag.iter().sum::<f64>();
        for c in 0..4 { assert!((h[c] - old_h[c] - dt * drag[c] * q * q).abs() < 1e-7); }
        near(report.net_stream_enthalpy_j, 0.0);
        near(report.mechanical.generated_heat_j, dt * 20.0 * q);
        assert!(observer.applied().is_none());
        assert_eq!(initial.time_s(), 0.0);

        let bulk_entropy = plan.as_ref().clone().with_thermal_entropy_limit();
        let unsupported = observer.clone().with_reconstruction(bulk_entropy).unwrap();
        assert!(initial.advance_observed(&unsupported, dt, thermal(), || Ok(())).is_err());
        assert!(unsupported.step_implicit(&initial, dt, controls, &[0.0; 4], || Ok(())).is_err());
        assert!(unsupported.applied().is_none());
        assert_eq!(unsupported.momentum_solves(), 0);
    }
}
