//! Original-owner composition and transaction checks, not flute qualification.
use super::*;
use crate::bernoulli_aperture::radiation::environment::multiport::source::ComputedBaffledInput;
use crate::nsf_dg::{NsfDgRadiationSolve, NsfDgReservoirRadiation};
use fs_lbm::nsf3::{NsfOpenState, NsfReservoirAperture};
use fs_material::gas_species::input::{GasRecord, read_gas_record};

pub(super) const ACOUSTIC_RATE: u32 = 8000;
pub(super) const SUBSTEPS: usize = (RATE / ACOUSTIC_RATE) as usize;
const CELL_COUNT: usize = 108;

fn gas_record(text: &str) -> GasRecord {
    read_gas_record(&mut text.split_ascii_whitespace(), None, None).unwrap()
}

pub(super) fn exterior() -> GasRecord {
    gas_record("mixture 300 100000 mole 3 N2 0.78 O2 0.21 Ar 0.01")
}

pub(super) fn dx(gas: &GasState) -> f64 {
    (3.0 * gas.specific_gas_constant * gas.temperature).sqrt() / f64::from(RATE)
}

fn fluid(medium: NsfGasMedium, gas: GasState) -> NsfGrid3<D3q27> {
    let reference = NsfOpenState {
        density_kg_m3: gas.density, velocity_m_s: [0.0; 3], temperature_k: gas.temperature,
    };
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: gas.pressure, reference_temperature_k: gas.temperature,
        reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
    };
    NsfGrid3::new(NsfSpec {
        dims: [12, 3, 3], periodic: [false; 3], gas: medium,
        temperature_k: gas.temperature, pressure_pa: gas.pressure, dx_m: dx(&gas),
        temperature_window_k: [290.0, 320.0], maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0, bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "manufactured complete-field radiation fixture",
        maximum_cells: CELL_COUNT, wall: WallThermal::Adiabatic, workers: 1,
        transport_points: 17,
    }, |_| false, |_| (1.0, [0.0; 3], gas.temperature)).unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0, faces: [face; 2], maximum_relative_perturbation: 0.05,
        }).unwrap()
        .with_reservoir_apertures(vec![
            NsfReservoirAperture { axis: 1, high: false, cells: vec![[5, 0, 1]],
                state: NsfOpenState {
                    density_kg_m3: (gas.pressure + 2.0)
                        / (gas.specific_gas_constant * gas.temperature),
                    ..reference
                } },
            NsfReservoirAperture { axis: 1, high: true, cells: vec![[5, 2, 1]],
                state: reference },
            NsfReservoirAperture { axis: 1, high: true, cells: vec![[9, 2, 1]],
                state: reference },
        ]).unwrap()
}

// Compute the original exterior once per test; fresh owners clone that cold
// source and its original propagator instead of repeating the circular bake.
pub(super) fn computed_primary(source: &ThermoviscousInterior<'_>, gas: &GasRecord, gate: &CancelGate)
    -> DgPortLoad
{
    let area = (3.0 * dx(source.gas())).powi(2);
    let duct = DgDuct::new(&Duct { segments: vec![Segment::Cylinder {
        radius: (area / std::f64::consts::PI).sqrt(), length: 0.03,
    }] }, source.gas(), DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(7)), DgBoundary::Port(PortId(1))).unwrap()
        .with_molecular_relaxation(source.relaxation().unwrap().clone()).unwrap();
    DgPortLoad::new_with_computed_baffled(duct, ACOUSTIC_RATE, 1e-11,
        MAX_DG_PROPAGATOR_DOFS, gas, ComputedBaffledInput::parse(
            "20 65 4 8 12 20 0.1 67108864 2 100 1e-8 64 1048576 257",
        ).unwrap(), gate).unwrap()
}

fn load<'a, 'c>(source: &'a ThermoviscousInterior<'c>, primary: &DgPortLoad,
    gate: &CancelGate) -> ThermoviscousDgLoad<'a, 'c>
{
    let dt = primary.propagator().dt();
    let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(),
        &mut || false).unwrap();
    let mut initial = prepared.zero_state();
    initial.molecular_temperature_k[0][0] = 0.0002;
    let cork = DgDuct::new(&Duct { segments: vec![Segment::Cylinder {
        radius: (primary.duct().inlet_area_m2() / std::f64::consts::PI).sqrt(),
        length: 0.012,
    }] }, source.gas(), DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(7)), DgBoundary::Closed).unwrap()
        .with_molecular_relaxation(source.relaxation().unwrap().clone()).unwrap();
    let cork = DgPortLoad::new(cork, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap();
    ThermoviscousDgLoad::new(ThermoviscousDg::new(source, prepared, primary.clone(),
        1, ThermoviscousDgLimits::default(), gate).unwrap()
        .with_additional_duct(0, cork, gate).unwrap()
        .with_initial_state(initial).unwrap(), 2, gate).unwrap()
}

pub(super) fn placement() -> NsfDgReservoirRadiation {
    NsfDgReservoirRadiation {
        origin_m: [-0.08, 0.0], cells: [16, 16], terminal_radial_rings: 8,
        terminal_angular_points: 32, pressure_bracket_pa: [-20.0, 20.0],
    }
}

fn solve() -> NsfDgRadiationSolve {
    NsfDgRadiationSolve {
        difference_step: 1e-5, maximum_iterations: 12,
        newton: fs_solver::NewtonKrylovConfig {
            absolute_tolerance: 1e-8, relative_tolerance: 0.0,
            linear_restart: 3, max_linear_cycles: 8,
            forcing_minimum: 1e-5, forcing_maximum: 0.5,
            forcing_gamma: 0.9, forcing_exponent: 1.5,
            globalization: fs_solver::Globalization::LineSearch(fs_solver::LineSearchConfig {
                armijo: 1e-4, contraction: 0.5, minimum_step: 1e-6,
            }),
        },
    }
}

pub(super) fn system<'a, 'c>(source: &'a ThermoviscousInterior<'c>, medium: NsfGasMedium,
    primary: &DgPortLoad, evaluations: usize, gate: &CancelGate)
    -> NsfThermoviscousDg<'a, 'c, D3q27>
{
    let mut budget = limits(SUBSTEPS);
    budget.maximum_cells = CELL_COUNT;
    budget.maximum_evaluations = evaluations;
    budget.maximum_cell_steps = CELL_COUNT * SUBSTEPS * evaluations;
    budget.incoming_bracket_m_s = [-0.5, 0.5];
    NsfThermoviscousDg::new(fluid(medium, *source.gas()), load(source, primary, gate),
        true, budget, gate).unwrap()
}

pub(super) fn attach<'a, 'c>(system: NsfThermoviscousDg<'a, 'c, D3q27>, gate: &CancelGate)
    -> NsfThermoviscousDg<'a, 'c, D3q27>
{
    // Deliberately retain a different order from the original fluid owners.
    system.with_finite_reservoir(0, 1e-5, finite_limits(), gate).unwrap()
        .with_reservoir_radiations(vec![(2, placement()), (1, placement())],
            solve(), gate).unwrap()
}

macro_rules! acoustic_state {
    ($load:expr) => {{
        let load = $load;
        let field = load.system();
        (field.state().clone(), field.load().state().clone(),
            field.load_at_port(0).unwrap().state().clone(),
            field.load().terminal().unwrap().pole_states().to_vec(),
            field.load().last_independent_mouths_frame().cloned(),
            field.accepted_steps(),
            load.independent_mouths_response(&CancelGate::new_clock_free()).unwrap())
    }};
}

macro_rules! original_state {
    ($owner:expr) => {{
        let owner = $owner;
        (
            (owner.zone().state_hash(), owner.zone().cumulative_ledger(),
                owner.zone().boundary_drives(), owner.zone().last_reservoir_flows().to_vec(),
                owner.zone().last_reservoir_exchanges().to_vec(),
                (0..3).map(|i| owner.zone().reservoir_aperture_at(i).unwrap().state)
                    .collect::<Vec<_>>()),
            acoustic_state!(owner.load()),
            (owner.finite_reservoir_inventory(), owner.accepted_steps()),
        )
    }};
}

#[test]
fn complete_affine_mouth_response_replays_original_fields_and_retains_cross_reaction() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let zone = fluid(medium, gas);
    let mouths: Vec<_> = [2, 1].map(|aperture|
        crate::nsf_dg::reservoir::radiation_surface(&zone, gas, &primary, aperture,
            placement(), &gate).unwrap().0).to_vec();
    let mut load = load(&field, &primary, &gate)
        .with_independent_mouths(&mouths, ACOUSTIC_RATE, &gate).unwrap();
    let prepared = field.prepare_time(load.duration_s(), ThermoviscousTimeLimits::default(),
        &mut || false).unwrap();
    assert_eq!(load.mouth_count(), 2);
    assert!(load.port_response(&gate).is_err());
    assert!(load.preview_step(1e-8, &gate).is_err());
    let cold = acoustic_state!(&load);
    {
        let _dropped = load.preview_step_with_mouth_flows(2e-8, &[1e-8, -0.7e-8], &gate).unwrap();
    }
    assert_eq!(acoustic_state!(&load), cold);
    let mut cross_reaction = 0.0_f64;
    for flows in [[2e-8, 1e-8, -0.7e-8], [-1e-8, 0.3e-8, 1e-8],
        [0.4e-8, -2e-8, 0.0], [1e-8, 0.0, -1e-8]]
    {
        let response = load.independent_mouths_response(&gate).unwrap();
        assert_eq!(response.history_pressure_pa.len(), 3);
        assert_eq!(response.impedance_pa_s_m3.len(), 9);
        let held_pressure = load.preview_step_with_mouth_flows(flows[0], &[0.0; 2],
            &gate).unwrap().mean_pressure_pa();
        let initial = load.system().state().clone();
        let mut bore = load.system().load().clone();
        let mut cork = load.system().load_at_port(0).unwrap().clone();
        let trial = load.preview_step_with_mouth_flows(flows[0], &flows[1..], &gate).unwrap();
        let receipt = trial.independent_mouths_frame().unwrap();
        let actual = [trial.mean_pressure_pa(), receipt.mouths[0].mean_pressure_pa,
            receipt.mouths[1].mean_pressure_pa];
        for row in 0..3 {
            let terms: Vec<_> = response.impedance_pa_s_m3[3 * row..3 * row + 3]
                .iter().zip(flows).map(|(z, q)| z * q).collect();
            let expected = response.history_pressure_pa[row] + terms.iter().sum::<f64>();
            let scale = response.history_pressure_pa[row].abs() + actual[row].abs()
                + terms.iter().map(|x| x.abs()).sum::<f64>();
            assert!((actual[row] - expected).abs() <= 1e-8 * scale.max(1e-20));
        }
        let mutual = response.impedance_pa_s_m3[1] * flows[1]
            + response.impedance_pa_s_m3[2] * flows[2];
        cross_reaction = cross_reaction.max(mutual.abs());
        assert!((actual[0] - held_pressure - mutual).abs()
            <= 1e-8 * (actual[0].abs() + held_pressure.abs() + mutual.abs()).max(1e-20));
        let frame = trial.frame();
        let replay = prepared.step(&initial, &[frame.additional_ducts[0].mean_pressure_pa,
            frame.mean_pressure_pa, actual[0]], &[], &mut || false).unwrap();
        assert_eq!(replay.state, frame.field.state);
        assert_eq!(replay.mean_inward_flow_m3_s, frame.field.mean_inward_flow_m3_s);
        let original = bore.preview_step_with_mouth_flows(
            -replay.mean_inward_flow_m3_s[1], &flows[1..]).unwrap();
        assert_eq!(original.output(), &frame.duct);
        assert_eq!(original.independent_mouths_frame(), Some(receipt));
        original.accept();
        assert_eq!(cork.preview_step(-replay.mean_inward_flow_m3_s[0]).unwrap().accept(),
            frame.additional_ducts[0].duct);
        near(trial.inlet_work_j(), actual[0] * replay.mean_inward_flow_m3_s[2]
            / f64::from(ACOUSTIC_RATE), 1e-12);
        // Mouth work belongs to the shared exterior, while this balance owns
        // only FE + finite DG storage and their actual circular terminal work.
        assert!(frame.balance_residual_j.abs() <= frame.balance_tolerance_j);
        trial.accept();
        assert_eq!(load.system().load().state(), bore.state());
        assert_eq!(load.system().load_at_port(0).unwrap().state(), cork.state());
        assert_eq!(load.system().load().terminal().unwrap().pole_states(),
            bore.terminal().unwrap().pole_states());
    }
    assert!(cross_reaction > 1e-12, "independent mouths must react through the complete field");
    assert_ne!(load.independent_mouths_response(&gate).unwrap().history_pressure_pa,
        cold.6.history_pressure_pa);
}

#[test]
fn finite_source_and_joint_escapes_replay_original_nsf_fe_and_each_finite_duct() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let mut owner = attach(system(&field, medium.clone(), &primary, 256, &gate), &gate);
    let mut original_fluid = fluid(medium, gas);
    let prepared = field.prepare_time(owner.load().duration_s(), ThermoviscousTimeLimits::default(),
        &mut || false).unwrap();
    assert_eq!(owner.reservoir_radiations().iter().map(|(i, _)| *i).collect::<Vec<_>>(), [2, 1]);
    assert_eq!(owner.reservoir_radiation_at(1), Some(&placement()));
    assert_eq!(owner.radiation_solve(), Some(&solve()));
    for (rate, heat) in [(0.0, 0.0), (-2e-7, 0.01), (1e-7, -0.005)] {
        let before = owner.finite_reservoir_inventory().unwrap();
        let initial = owner.load().system().state().clone();
        let mut bore = owner.load().system().load().clone();
        let mut cork = owner.load().system().load_at_port(0).unwrap().clone();
        let frame = owner.step_with_drive(NsfDgDrive::FiniteReservoirActuation {
            aperture: 0, volume_rate_m3_s: rate, heat_power_w: heat,
            maximum_relative_volume_change: 0.05,
        }, &gate).unwrap();
        let mut controls = vec![(0, original_fluid.reservoir_aperture_state_for(0,
            before.gas.pressure, before.gas.temperature).unwrap())];
        for mouth in &frame.reservoir_radiations {
            let aperture = mouth.resolved_exchange.aperture;
            controls.push((aperture, original_fluid.reservoir_aperture_state_for(aperture,
                gas.pressure + mouth.ghost_pressure_offset_pa, gas.temperature).unwrap()));
        }
        let mut drives = [NsfBoundaryDrive::IncomingCharacteristic(0.0); 2];
        drives[1] = NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s);
        let original = original_fluid.preview_boundary_interval_with_reservoirs(drives,
            None, &controls, SUBSTEPS, CELL_COUNT * SUBSTEPS, || false).unwrap();
        assert_eq!(original.boundary_exchange(), &frame.faces);
        assert_eq!(original.reservoir_exchanges(), frame.reservoir_exchanges);
        assert_eq!(original.reservoir_flows(), frame.reservoir_flows);
        assert_eq!(original.energy_ledger(), &frame.flow);
        original.accept();
        assert_eq!(original_fluid.state_hash(), owner.zone().state_hash());
        let complete = owner.load().system().last_frame().unwrap();
        let replay = prepared.step(&initial, &[complete.additional_ducts[0].mean_pressure_pa,
            complete.mean_pressure_pa, frame.inlet_pressure_pa], &[], &mut || false).unwrap();
        assert_eq!(replay.state, complete.field.state);
        assert_eq!(replay.mean_inward_flow_m3_s, complete.field.mean_inward_flow_m3_s);
        let flows: Vec<_> = frame.reservoir_radiations.iter().map(|x| x.outward_flow_m3_s).collect();
        assert_eq!(frame.reservoir_radiations.iter()
            .map(|x| x.resolved_exchange.aperture).collect::<Vec<_>>(), [2, 1]);
        let original = bore.preview_step_with_mouth_flows(-replay.mean_inward_flow_m3_s[1],
            &flows).unwrap();
        assert_eq!(original.output(), &complete.duct);
        for (source, receipt) in frame.reservoir_radiations.iter()
            .zip(&original.independent_mouths_frame().unwrap().mouths)
        {
            assert_eq!(source.radiation_pressure_pa, receipt.mean_pressure_pa);
            assert_eq!(source.radiation_work_j, receipt.load_work_j);
            assert_eq!(source.pressure_residual_pa,
                source.mean_acoustic_pressure_pa - source.radiation_pressure_pa);
            assert!(source.pressure_residual_pa.abs() < 1e-6);
        }
        original.accept();
        assert_eq!(cork.preview_step(-replay.mean_inward_flow_m3_s[0]).unwrap().accept(),
            complete.additional_ducts[0].duct);
        assert_eq!(owner.load().system().load().state(), bore.state());
        assert_eq!(owner.load().system().load_at_port(0).unwrap().state(), cork.state());
        assert_eq!(frame.inlet_work_j, complete.authored_external_work_j);
        assert!(complete.balance_residual_j.abs() <= complete.balance_tolerance_j);
        let after = owner.finite_reservoir_inventory().unwrap();
        let source = frame.finite_reservoir.as_ref().unwrap();
        assert_eq!(source.source_exchange, frame.reservoir_exchanges[0]);
        let motion = source.actuation.as_ref().unwrap();
        assert!((after.mass_kg - before.mass_kg - source.source_exchange.mass_out_kg).abs() <= 1e-16);
        assert!((after.sensible_energy_j - before.sensible_energy_j
            - source.source_exchange.total_energy_out_j - motion.step.supplied_heat_j
            + motion.step.work_by_gas_j).abs() <= 1e-9);
        assert_eq!(motion.step.supplied_heat_j, heat * owner.load().duration_s());
        assert!((after.volume_m3 - before.volume_m3 - rate * owner.load().duration_s()).abs()
            <= 4.0 * f64::EPSILON * before.volume_m3);
    }
    assert_eq!(owner.zone().steps(), 3 * SUBSTEPS);
    assert_eq!(owner.load().system().accepted_steps(), 3);
}

#[test]
fn numerical_budget_and_late_chamber_failure_preserve_all_owners_and_next_warm_start() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let make = |evaluations| attach(system(&field, medium.clone(), &primary,
        evaluations, &gate), &gate);
    let drive = NsfDgDrive::FiniteReservoirActuation {
        aperture: 0, volume_rate_m3_s: -2e-7, heat_power_w: 0.01,
        maximum_relative_volume_change: 0.05,
    };
    let mut exhausted = make(3);
    let cold = original_state!(&exhausted);
    assert!(matches!(exhausted.step_with_drive(drive, &gate),
        Err(NsfDgError::Convergence { evaluations: 2, .. })));
    assert_eq!(original_state!(&exhausted), cold);
    let mut owner = make(256);
    owner.step_with_drive(drive, &gate).unwrap();
    let accepted = original_state!(&owner);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(owner.step_with_drive(drive, &cancelled), Err(NsfDgError::Cancelled)));
    assert_eq!(original_state!(&owner), accepted);
    let impossible = NsfDgDrive::FiniteReservoirActuation {
        aperture: 0, volume_rate_m3_s: -1.0, heat_power_w: 0.01,
        maximum_relative_volume_change: 0.05,
    };
    assert!(owner.validate_drive(impossible).is_ok());
    assert!(matches!(owner.step_with_drive(impossible, &gate), Err(NsfDgError::Input(
        "finite reservoir actuation loses an increment or exceeds its positive bounded volume movement; refine the interval"))));
    assert_eq!(original_state!(&owner), accepted);
    let retry = owner.step_with_drive(drive, &gate).unwrap();
    let mut uninterrupted = make(256);
    uninterrupted.step_with_drive(drive, &gate).unwrap();
    let expected = uninterrupted.step_with_drive(drive, &gate).unwrap();
    assert_eq!(original_state!(&owner), original_state!(&uninterrupted));
    assert_eq!(retry.evaluations, expected.evaluations);
    assert_eq!(retry.incoming_characteristic_m_s.to_bits(), expected.incoming_characteristic_m_s.to_bits());
    assert_eq!(retry.reservoir_radiations, expected.reservoir_radiations);
}

#[test]
fn cold_admission_retains_exclusive_aperture_geometry_clock_and_material_owners() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let field = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let make = || system(&field, medium.clone(), &primary, 256, &gate);
    for inputs in [Vec::new(), vec![(1, placement()), (1, placement())],
        vec![(3, placement())], vec![(1, NsfDgReservoirRadiation {
            pressure_bracket_pa: [f64::NAN, 20.0], ..placement()
        })]]
    {
        assert!(make().with_reservoir_radiations(inputs, solve(), &gate).is_err());
    }
    assert!(make().with_finite_reservoir(0, 1e-5, finite_limits(), &gate).unwrap()
        .with_reservoir_radiation_at(0, placement(), solve(), &gate).is_err());
    let one = make().with_reservoir_radiation_at(1, placement(), solve(), &gate).unwrap();
    assert_eq!(one.reservoir_radiations(), &[(1, placement())]);
    assert!(one.with_finite_reservoir(1, 1e-5, finite_limits(), &gate).is_err());
    let mut overlap = placement();
    overlap.origin_m[1] -= 4.0 * dx(&gas);
    assert!(make().with_reservoir_radiations(vec![(1, placement()), (2, overlap)],
        solve(), &gate).is_err());
    let surface = crate::nsf_dg::reservoir::radiation_surface(&fluid(medium.clone(), gas),
        gas, &primary, 1, placement(), &gate).unwrap().0;
    assert!(load(&field, &primary, &gate).with_independent_mouths(&[surface],
        ACOUSTIC_RATE * 2, &gate).is_err());
    let changed = computed_primary(&field, &gas_record("mixture 300 100000 mole 1 N2 1"), &gate);
    assert!(system(&field, medium.clone(), &changed, 256, &gate)
        .with_reservoir_radiation_at(1, placement(), solve(), &gate).is_err());
    let mut advanced = attach(make(), &gate);
    advanced.step_with_drive(NsfDgDrive::FiniteReservoir { aperture: 0 }, &gate).unwrap();
    assert!(advanced.with_reservoir_radiation_at(1, placement(), solve(), &gate).is_err());
}
