//! Native finite gas storage, source transport and live duct/wall memory.
use super::*;
use fs_couple::bernoulli_aperture::wall::material::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgLoss, FiniteGasDgSolve,
};
use fs_duct::dg::{DgPropagator, DgSidePort, DgState, PortId, PortInput};
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{GasNozzle, GasNozzleStepSolve};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const RADIUS: f64 = 0.001;
const LENGTH: f64 = 0.3;
fn volume() -> f64 { core::f64::consts::PI * RADIUS * RADIUS * LENGTH }
fn point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static)), temperature).unwrap()
        .with_quantity("pressure", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Pressure, ValueForm::Static)), pressure).unwrap()
        .with_quantity("material-lot", QuantitySpec::dimensional(Dims([0; 6])), 42.0).unwrap()
}
fn initial_field() -> AxialGasField {
    let catalog = gas_field(1.0, 300.0).unwrap();
    AxialGasField::new(vec![0.0, LENGTH], 101325.0,
        vec![PrimitiveGasState { temperature_k: 300.0, composition: composition(1.0) }; 2],
        catalog.model().clone(),
        QueryPoint::new().with_quantity("material-lot", QuantitySpec::dimensional(Dims([0; 6])), 42.0).unwrap(),
        ThermochemicalPhase::FrozenGas,
    ).unwrap()
}
fn parcel(fraction: f64, temperature: f64, pressure: f64) -> ClosedGasParcel {
    let field = initial_field();
    ClosedGasParcel::new(field.model().with_composition(composition(fraction)).unwrap(),
        point(temperature, pressure), volume()).unwrap()
}
fn energy_solve() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 900.0,
        energy_tolerance_j_kg: 1e-8, maximum_iterations: 128 }
}
fn controls(diffusion: bool) -> FiniteGasDgSolve {
    FiniteGasDgSolve {
        work_bracket_j: [-1e-3, 1e-3], work_tolerance_j: 1e-13,
        energy_tolerance_j: 1e-11, maximum_iterations: 80,
        maximum_temperature_increment_k: 80.0,
        loss: if diffusion { FiniteGasDgLoss::UniformViscothermal { maximum_relative_error: 0.001 } }
            else { FiniteGasDgLoss::Lossless },
    }
}
fn acoustic(field: &AxialGasField, diffusion: bool) -> DgDuct {
    let duct = DgDuct::new_with_gas_field(
        &fs_duct::Duct { segments: vec![fs_duct::Segment::Cylinder { radius: RADIUS, length: LENGTH }] },
        field, DgMesh { elements: 2, order: 1 }, DgBoundary::Closed, DgBoundary::Closed,
    ).unwrap().with_side_ports(&[DgSidePort { id: PortId(7), downstream_element: 1 }]).unwrap();
    if diffusion { duct.with_viscothermal(16, [50.0, 2000.0], 0.001).unwrap() } else { duct }
}
fn regional() -> (RegionalMaterialWallDg, DgState, DgPropagator, ClosedGasParcel) {
    let source = source_wall(1800.0, 4e8);
    let hereditary = LocalMaterialWall::new(source.material().clone(), 2e-5,
        source.axial_constraint(), source.loss().clone(), source.band_hz()).unwrap();
    let bound = RegionalMaterialWallDg::new(
        vec![kelvin_wall(1000.0, 1e8, 2.0, 2e-5), hereditary],
        acoustic(&initial_field(), true),
    ).unwrap();
    let propagator = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|x| 1000.0 + 300.0*x, |_| 2e-6);
    for _ in 0..12 { propagator.step(&mut state, &[]).unwrap(); }
    (bound, state, propagator, parcel(1.0, 300.0, 101325.0))
}
fn energy_change(old: &ClosedGasParcel, next: &ClosedGasParcel) -> f64 {
    next.mass_kg() * (next.state().specific_internal_energy_j_kg() - old.state().specific_internal_energy_j_kg())
        + (next.mass_kg()-old.mass_kg())*old.state().specific_internal_energy_j_kg()
}
fn same_motion(a: &fs_duct::dg::DgCompliantWallState, b: &fs_duct::dg::DgCompliantWallState) {
    close(a.stored_energy_j, b.stored_energy_j);
    close(a.relaxation_stored_energy_j, b.relaxation_stored_energy_j);
    for (x,y) in a.displacement_m.iter().chain(&a.velocity_m_s)
        .zip(b.displacement_m.iter().chain(&b.velocity_m_s)) { close(*x,*y); }
    for (x,y) in a.relaxation_pressure_pa.iter().flatten().zip(b.relaxation_pressure_pa.iter().flatten()) {
        close(*x,*y);
    }
}

#[test]
fn heating_pays_actual_compressive_work_and_preserves_the_native_clock_and_ports() {
    let field = initial_field();
    let mut bound = kelvin_wall(1500.0, 1e9, 0.0, 2e-5).bind_dg(acoustic(&field, false)).unwrap();
    let mut propagator = bound.duct().propagator(1e-6, 1e-13).unwrap();
    let mut state = bound.duct().project_initial_state(|x| 3000.0 + 1000.0*x, |_| 0.0);
    let mut gas = parcel(1.0, 300.0, 101325.0);
    let old_gas = gas.clone();
    let old_state = state.clone();
    let ports = [PortInput { id: PortId(7), effort: 1e-7 }];
    let mut marched = state.clone();
    let original_step = propagator.step(&mut marched, &ports).unwrap();
    propagator.sync_nodal(&mut marched).unwrap();
    let motion = bound.duct().compliant_wall_state(&marched).unwrap();
    let old_bulk = gas.state().state().density * gas.state().state().sound_speed.powi(2);
    // Independent exact integral of each order-one pressure polynomial squared.
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let pressure_energy: f64 = marched.pressure().chunks_exact(2).map(|p|
        area*(LENGTH/2.0)/(6.0*old_bulk)*(p[0]*p[0]+p[0]*p[1]+p[1]*p[1])
    ).sum();
    let heat = 0.02;
    let step = bound.step_with_finite_gas_heat_controlled(
        &mut state, &mut propagator, &mut gas, &ports, heat,
        energy_solve(), controls(false), &mut || false,
    ).unwrap();
    // Synthetic source Cp/R = 3.5 + 0.001*T: solve its independent polynomial
    // caloric law plus physical B(T) pressure storage, never the runtime residual.
    let (t0, mass, r) = (300.0, old_gas.mass_kg(), old_gas.state().state().specific_gas_constant);
    let bulk = |t: f64| old_gas.state().state().density*r*t*(3.5+0.001*t)/(2.5+0.001*t);
    let balance = |t: f64| mass*r*(2.5*(t-t0)+0.0005*(t*t-t0*t0))
        + pressure_energy*(bulk(t)/bulk(t0)-1.0)-heat;
    let (mut lo, mut hi) = (300.0, 380.0);
    for _ in 0..60 { let mid=(lo+hi)/2.0; if balance(mid)>0.0 { hi=mid; } else { lo=mid; } }
    assert!((gas.state().state().temperature-(lo+hi)/2.0).abs() < 2e-7);
    let expected_work = pressure_energy*(bulk(gas.state().state().temperature)/bulk(t0)-1.0);
    assert!((step.refresh.parameter_work_j-expected_work).abs() < 2e-12);
    assert!(step.refresh.parameter_work_j > 1e-7);
    let mut heat_only = old_gas.clone();
    heat_only.add_heat(heat, energy_solve(), || false).unwrap();
    assert!(heat_only.state().state().temperature-gas.state().state().temperature > 1e-4);
    assert_eq!(state.step_count(), old_state.step_count()+1);
    assert_eq!(state.initial_energy_j(), old_state.initial_energy_j());
    assert_eq!(propagator.dt(), 1e-6);
    same_motion(&bound.duct().compliant_wall_state(&state).unwrap(), &motion);
    assert_eq!(step.acoustic.port_outputs[0].delta_x, original_step.port_outputs[0].delta_x);
    assert_eq!(step.acoustic.ledger.port_work_step_j, original_step.ledger.port_work_step_j);
    assert_eq!(step.acoustic.port_outputs[0].end_conjugate,
        bound.duct().port_pressures(&state).unwrap()[0].1);
    assert_ne!(step.acoustic.port_outputs[0].end_conjugate, original_step.port_outputs[0].end_conjugate);
    assert_eq!(step.acoustic.ledger.stored_energy_j, state.stored_energy_j());
    assert_eq!(bound.duct().gas_field().unwrap().uniform_state(), Some(gas.state()));
    assert_eq!(gas.point().axis_quantities().get("material-lot"), old_gas.point().axis_quantities().get("material-lot"));
    assert_eq!(gas.model().source_packs(), old_gas.model().source_packs());
    assert!(step.combined_energy_residual_j.abs() <= controls(false).energy_tolerance_j);
    let later = propagator.step(&mut state, &[]).unwrap();
    assert_eq!(later.ledger.cumulative_parameter_work_j, step.refresh.parameter_work_j);
}

#[test]
fn native_nozzle_species_and_enthalpy_reach_continuing_diffusion_and_regional_wall_history() {
    let (mut bound, mut state, mut propagator, mut gas) = regional();
    let reservoir = parcel(0.0, 400.0, 180000.0);
    let (before, old_gas, old_wall, old_prop) = (state.clone(), gas.clone(), bound.walls().to_vec(), propagator.clone());
    let mut marched = state.clone();
    old_prop.step(&mut marched, &[]).unwrap();
    let motion = bound.duct().compliant_wall_state(&marched).unwrap();
    assert!(motion.relaxation_stored_energy_j > 0.0);
    let old_model = bound.duct().viscothermal_model().unwrap().clone();
    let flow_controls = GasNozzleStepSolve {
        nozzle: energy_solve(),
        transfer: GasMassTransferSolve { donor: energy_solve(), receiver: energy_solve(),
            mass_tolerance_kg: 1e-16, energy_tolerance_j: 1e-11 },
        maximum_donor_fraction: 0.05,
    };
    let step = bound.step_with_finite_gas_controlled(
        &mut state, &mut propagator, &mut gas, &[], controls(true),
        |old, duration, work, cancel| {
            let (proposal, flow) = reservoir.through_nozzle_to(old,
                GasNozzle { area_m2: 1e-6, discharge_coefficient: 0.8 },
                duration, flow_controls, &mut *cancel)?;
            let (next_reservoir, mut next_gas, transfer) = proposal.into_parts();
            let model = next_gas.model().clone();
            let v = next_gas.volume_m3();
            let work_step = next_gas.replace_model_with_heat_work_volume(
                model, v, 0.0, work, energy_solve(), &mut *cancel,
            )?;
            // The caller retains the original physical direction after its
            // work correction; the original nozzle owner chose flow and inventories.
            if next_reservoir.state().state().pressure < next_gas.state().state().pressure {
                return Err(GasEnergyError::Invalid("endpoint work crossed original nozzle ordering"));
            }
            Ok(FiniteGasDgBackground { gas: next_gas, supplied_heat_j: 0.0,
                transported_enthalpy_j: transfer.transferred_enthalpy_j,
                source: (next_reservoir, transfer, flow, work_step) })
        }, &mut || false,
    ).unwrap();
    let (next_reservoir, transfer, flow, _) = &step.background;
    assert!(flow.mass_flow_kg_s > 0.0);
    assert_eq!(transfer.transferred_mass_kg, flow.mass_flow_kg_s*propagator.dt());
    assert!(gas.mass_kg() > old_gas.mass_kg());
    assert_eq!(gas.state().components().len(), 2);
    assert_eq!(gas.state().source_packs().len(), 2);
    assert!(gas.state().state().density > old_gas.state().state().density);
    assert_ne!(gas.state().state().sound_speed, old_gas.state().state().sound_speed);
    assert_ne!(bound.duct().viscothermal_model().unwrap().series_terms()[0].rate_per_s,
        old_model.series_terms()[0].rate_per_s);
    assert_eq!(bound.walls(), old_wall.as_slice());
    same_motion(&bound.duct().compliant_wall_state(&state).unwrap(), &motion);
    assert_eq!(bound.duct().gas_field().unwrap().uniform_state(), Some(gas.state()));
    assert_eq!(bound.duct().side_port_source(PortId(7)).unwrap().gas, *gas.state().state());
    assert_eq!(state.step_count(), before.step_count()+1);
    let total = energy_change(&old_gas, &gas)+energy_change(&reservoir, next_reservoir)
        + state.stored_energy_j()-before.stored_energy_j()+step.acoustic_dissipation_j;
    assert!(total.abs() < 2e-11, "{total}");
    assert_eq!(step.supplied_heat_j, 0.0);
    assert_eq!(step.transported_enthalpy_j, transfer.transferred_enthalpy_j);
    let (mut unmodified, mut difference) = (marched, 0.0);
    for _ in 0..16 {
        let frame = propagator.step(&mut state, &[]).unwrap();
        old_prop.step(&mut unmodified, &[]).unwrap();
        propagator.sync_nodal(&mut state).unwrap();
        old_prop.sync_nodal(&mut unmodified).unwrap();
        difference += state.pressure().iter().zip(unmodified.pressure())
            .map(|(a,b)| (a-b).powi(2)).sum::<f64>();
        assert_eq!(frame.ledger.cumulative_parameter_work_j, step.refresh.parameter_work_j);
    }
    assert!(difference > 0.0);
}

#[test]
fn late_background_failure_work_refusal_and_cancellation_preserve_all_accepted_owners() {
    let (mut bound, mut state, mut propagator, mut gas) = regional();
    let (old_bound, old_state, old_prop, old_gas) = (bound.clone(), state.clone(), propagator.clone(), gas.clone());
    let mut stale_gas = parcel(1.0, 301.0, 101325.0);
    let stale_before = stale_gas.clone();
    assert!(bound.step_with_finite_gas_heat_controlled(
        &mut state, &mut propagator, &mut stale_gas, &[], 0.02,
        energy_solve(), controls(true), &mut || false,
    ).is_err());
    assert_eq!(stale_gas, stale_before);
    assert_eq!(state, old_state);
    let mut trials = 0;
    let result = bound.step_with_finite_gas_controlled(
        &mut state, &mut propagator, &mut gas, &[], controls(true),
        |old, _, work, cancel| {
            trials += 1;
            if trials > 1 { return Err(GasEnergyError::Invalid("late background source failure")); }
            let mut next = old.clone();
            let source = next.replace_model_with_heat_work_volume(
                old.model().clone(), old.volume_m3(), 0.02, work, energy_solve(), cancel,
            )?;
            Ok(FiniteGasDgBackground { gas: next, supplied_heat_j: 0.02,
                transported_enthalpy_j: 0.0, source })
        }, &mut || false,
    );
    assert!(matches!(result, Err(FiniteGasDgError::Background(_))));
    assert_eq!(trials, 2);
    assert_eq!(state, old_state);
    assert_eq!(gas, old_gas);
    let zero_work = FiniteGasDgSolve { work_bracket_j: [0.0,0.0], ..controls(true) };
    assert!(bound.step_with_finite_gas_heat_controlled(
        &mut state, &mut propagator, &mut gas, &[], 0.02,
        energy_solve(), zero_work, &mut || false,
    ).is_err());
    let mut polls = 0;
    old_bound.clone().step_with_finite_gas_heat_controlled(
        &mut old_state.clone(), &mut old_prop.clone(), &mut old_gas.clone(),
        &[], 0.02, energy_solve(), controls(true), &mut || { polls += 1; false },
    ).unwrap();
    let mut count = 0;
    assert!(matches!(bound.step_with_finite_gas_heat_controlled(
        &mut state, &mut propagator, &mut gas, &[], 0.02,
        energy_solve(), controls(true), &mut || { count += 1; count == polls },
    ), Err(FiniteGasDgError::Cancelled)));
    assert_eq!(state, old_state);
    assert_eq!(gas, old_gas);
    assert_eq!(bound.walls(), old_bound.walls());
    assert_eq!(bound.duct().gas_field(), old_bound.duct().gas_field());
    assert_eq!(bound.duct().operator_identity(), old_bound.duct().operator_identity());
    let mut reference = old_state;
    assert_eq!(propagator.step(&mut state, &[]).unwrap(), old_prop.step(&mut reference, &[]).unwrap());
}

