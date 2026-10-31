//! Two original complete acoustic assemblies on distinct NSF interfaces.
//! Manufactured composition, work and transaction checks, not jet qualification.
use super::*;
use super::radiation::{ACOUSTIC_RATE, computed_primary, dx, exterior, placement};
use crate::nsf_dg::NsfDgRadiationSolve;

pub(super) fn solve() -> NsfDgRadiationSolve {
    NsfDgRadiationSolve {
        difference_step: 1e-5, maximum_iterations: 12,
        newton: fs_solver::NewtonKrylovConfig {
            absolute_tolerance: 1e-8, relative_tolerance: 0.0,
            linear_restart: 4, max_linear_cycles: 8,
            forcing_minimum: 1e-5, forcing_maximum: 0.5,
            forcing_gamma: 0.9, forcing_exponent: 1.5,
            globalization: fs_solver::Globalization::LineSearch(fs_solver::LineSearchConfig {
                armijo: 1e-4, contraction: 0.5, minimum_step: 1e-6,
            }),
        },
    }
}

fn plain_primary(source: &ThermoviscousInterior<'_>, dt: f64, length: f64,
    end: DgBoundary) -> DgPortLoad
{
    let area = source.inner().plane_wave_port_area_m2(1).unwrap();
    let duct = DgDuct::new(&Duct { segments: vec![Segment::Cylinder {
        radius: (area / std::f64::consts::PI).sqrt(), length,
    }] }, source.gas(), DgMesh { elements: 2, order: 2 },
        DgBoundary::Port(PortId(7)), end).unwrap()
        .with_molecular_relaxation(source.relaxation().unwrap().clone()).unwrap();
    DgPortLoad::new(duct, dt, 1e-11, MAX_DG_PROPAGATOR_DOFS).unwrap()
}

pub(super) fn field_load<'a, 'c>(source: &'a ThermoviscousInterior<'c>,
    primary: &DgPortLoad, limits: ThermoviscousDgLimits, gate: &CancelGate)
    -> ThermoviscousDgLoad<'a, 'c>
{
    let dt = primary.propagator().dt();
    let prepared = source.prepare_time(dt, ThermoviscousTimeLimits::default(),
        &mut || false).unwrap();
    let mut initial = prepared.zero_state();
    // Distinct original caloric history from the downstream fixture's 0.0002 K.
    initial.molecular_temperature_k[0][0] = 0.0003;
    let cork = plain_primary(source, dt, 0.018, DgBoundary::Closed);
    ThermoviscousDgLoad::new(ThermoviscousDg::new(source, prepared, primary.clone(),
        1, limits, gate).unwrap().with_additional_duct(0, cork, gate).unwrap()
        .with_initial_state(initial).unwrap(), 2, gate).unwrap()
}

#[allow(clippy::too_many_arguments)]
pub(super) fn coupled<'a, 'c>(downstream: &'a ThermoviscousInterior<'c>,
    upstream: &'a ThermoviscousInterior<'c>, medium: NsfGasMedium,
    downstream_primary: &DgPortLoad, upstream_primary: &DgPortLoad,
    upstream_limits: ThermoviscousDgLimits, gate: &CancelGate)
    -> NsfThermoviscousDg<'a, 'c, D3q27>
{
    super::radiation::system(downstream, medium, downstream_primary, 256, gate)
        .with_upstream_load(field_load(upstream, upstream_primary, upstream_limits, gate),
            solve(), gate).unwrap()
        .with_finite_reservoir(0, 1e-5, finite_limits(), gate).unwrap()
        .with_reservoir_radiations(vec![(2, placement()), (1, placement())],
            solve(), gate).unwrap()
}

macro_rules! field_state {
    ($load:expr) => {{
        let load = $load;
        let system = load.system();
        (system.state().clone(), system.load().state().clone(),
            system.load_at_port(0).unwrap().state().clone(),
            system.load().terminal().map(|terminal| terminal.pole_states().to_vec()),
            system.load().last_independent_mouths_frame().cloned(), system.accepted_steps(),
            (if load.mouth_count() == 0 { Some(load.port_response(&CancelGate::new_clock_free()).unwrap()) }
                else { None },
                if load.mouth_count() != 0 {
                    Some(load.independent_mouths_response(&CancelGate::new_clock_free()).unwrap())
                } else { None }))
    }};
}

macro_rules! all_owners {
    ($owner:expr) => {{
        let owner = $owner;
        ((owner.zone().state_hash(), owner.zone().cumulative_ledger(),
            owner.zone().boundary_drives(), owner.zone().last_reservoir_flows().to_vec(),
            owner.zone().last_reservoir_exchanges().to_vec(),
            (0..owner.zone().reservoir_aperture_count())
                .map(|i| owner.zone().reservoir_aperture_at(i).unwrap().state).collect::<Vec<_>>()),
            field_state!(owner.load()), field_state!(owner.upstream_load().unwrap()),
            (owner.finite_reservoir_inventory(), owner.accepted_steps()))
    }};
}

const PATCH: NsfCharacteristicPatch = NsfCharacteristicPatch { start: [0, 0], end: [1, 2] };
const BACKGROUND: [NsfBoundaryDrive; 2] = [
    NsfBoundaryDrive::PressureOffset(-0.25), NsfBoundaryDrive::PressureOffset(0.15),
];

fn replay_zone(gas: GasState, patched: bool) -> NsfGrid3<D3q27> {
    let face = NsfCharacteristicBoundary {
        reference_pressure_pa: gas.pressure, reference_temperature_k: gas.temperature,
        reference_velocity_m_s: [0.0; 3], incoming_characteristic_m_s: 0.0,
    };
    let mut fluid = NsfGrid3::new(NsfSpec {
        dims: [8, 2, 2], periodic: [false; 3], gas: GasSpec::dry_air_ussa1976(),
        temperature_k: gas.temperature, pressure_pa: gas.pressure, dx_m: dx(&gas),
        temperature_window_k: [250.0, 600.0], maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0, bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "manufactured two-field interface fixture",
        maximum_cells: CELLS, wall: WallThermal::Adiabatic, workers: 1, transport_points: 17,
    }, |_| false, |x| {
        let z = (x[0] as f64 - 3.0) / 4.0;
        let pulse = 1e-5 * (0.7 + 0.1 * x[1] as f64) * fs_math::det::exp(-z * z);
        (1.0 + pulse / gas.gamma,
            [(gas.gamma / 3.0).sqrt() * pulse / gas.gamma, 0.0, 0.0],
            gas.temperature * (1.0 + (gas.gamma - 1.0) / gas.gamma * pulse))
    }).unwrap().with_characteristic_axis(NsfCharacteristicAxis {
        axis: 0, faces: [face; 2], maximum_relative_perturbation: 0.05,
    }).unwrap();
    if patched {
        fluid.set_boundary_drives(BACKGROUND).unwrap();
        fluid = fluid.with_characteristic_patch(false, PATCH, BACKGROUND[0]).unwrap()
            .with_characteristic_patch(true, PATCH, BACKGROUND[1]).unwrap();
    }
    fluid
}

#[test]
fn both_full_fields_replay_original_owners_at_whole_and_selected_interfaces() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let gas = reference();
    let dt = 1.0 / f64::from(RATE);
    for patched in [false, true] {
        let area = if patched { 2.0 } else { 4.0 } * dx(&gas).powi(2);
        let downstream = field_source(&mesh, &unit, area.sqrt(), gas);
        let upstream = field_source(&mesh, &unit, area.sqrt(), gas);
        let down_primary = plain_primary(&downstream, dt, 0.03, DgBoundary::Matched);
        let up_primary = plain_primary(&upstream, dt, 0.05, DgBoundary::Closed);
        let down_prepared = downstream.prepare_time(dt, ThermoviscousTimeLimits::default(),
            &mut || false).unwrap();
        let up_prepared = upstream.prepare_time(dt, ThermoviscousTimeLimits::default(),
            &mut || false).unwrap();
        for high in [false, true] {
            let mut budget = limits(1);
            budget.maximum_evaluations = 256;
            budget.maximum_cell_steps = CELLS * 256;
            let load = field_load(&downstream, &down_primary, ThermoviscousDgLimits::default(), &gate);
            let mut owner = if patched {
                NsfThermoviscousDg::new_on_patches(replay_zone(gas, true), load, high, budget, &gate)
            } else { NsfThermoviscousDg::new(replay_zone(gas, false), load, high, budget, &gate) }
                .unwrap().with_upstream_load(field_load(&upstream, &up_primary,
                    ThermoviscousDgLimits::default(), &gate), solve(), &gate).unwrap();
            let mut fluid = replay_zone(gas, patched);
            let mut complement = 0.0_f64;
            for rise in [1.0, -0.5, 0.0] {
                let down_initial = owner.load().system().state().clone();
                let up_initial = owner.upstream_load().unwrap().system().state().clone();
                let mut ducts = [owner.load().system().load().clone(),
                    owner.load().system().load_at_port(0).unwrap().clone(),
                    owner.upstream_load().unwrap().system().load().clone(),
                    owner.upstream_load().unwrap().system().load_at_port(0).unwrap().clone()];
                let frame = owner.step_with_upstream_pressure_rise(rise, &gate).unwrap();
                let up = frame.upstream.as_ref().unwrap();
                assert_eq!(up.face, usize::from(!high));
                let mut controls = [NsfBoundaryDrive::IncomingCharacteristic(0.0); 2];
                controls[usize::from(high)] =
                    NsfBoundaryDrive::IncomingCharacteristic(frame.incoming_characteristic_m_s);
                controls[usize::from(!high)] =
                    NsfBoundaryDrive::IncomingCharacteristic(up.incoming_characteristic_m_s);
                let original = if patched {
                    fluid.preview_boundary_interval_with_patches(BACKGROUND,
                        controls.map(Some), 1, CELLS, || false)
                } else { fluid.preview_boundary_interval(controls, 1, CELLS, || false) }.unwrap();
                assert_eq!(original.boundary_exchange(), &frame.faces);
                original.accept();
                assert_eq!(fluid.state_hash(), owner.zone().state_hash());
                for (index, (load, initial, prepared, pressure, work)) in [
                    (owner.load(), &down_initial, &down_prepared, frame.inlet_pressure_pa, frame.inlet_work_j),
                    (owner.upstream_load().unwrap(), &up_initial, &up_prepared, up.inlet_pressure_pa, up.inlet_work_j),
                ].into_iter().enumerate() {
                    let complete = load.system().last_frame().unwrap();
                    let replay = prepared.step(initial, &[complete.additional_ducts[0].mean_pressure_pa,
                        complete.mean_pressure_pa, pressure], &[], &mut || false).unwrap();
                    assert_eq!(complete.field.state, replay.state);
                    assert_eq!(complete.field.mean_inward_flow_m3_s, replay.mean_inward_flow_m3_s);
                    for (local, port) in [1, 0].into_iter().enumerate() {
                        let duct = &mut ducts[2 * index + local];
                        let output = duct.preview_step(-replay.mean_inward_flow_m3_s[port]).unwrap().accept();
                        assert_eq!(if local == 0 { &complete.duct } else { &complete.additional_ducts[0].duct }, &output);
                        assert_eq!(load.system().load_at_port(port).unwrap().state(), duct.state());
                    }
                    assert_eq!(work, complete.authored_external_work_j);
                    assert!(complete.balance_residual_j.abs() <= complete.balance_tolerance_j);
                }
                assert_eq!(up.inlet_flow_m3_s, up.field_exchange.normal_volume_flow_m3_s);
                assert_eq!(up.source_work_j, -rise * (up.inlet_flow_m3_s * dt));
                assert_eq!(up.acoustic_work_residual_j,
                    up.field_exchange.acoustic_pressure_work_out_j + up.source_work_j - up.inlet_work_j);
                assert!((up.field_exchange.mean_acoustic_pressure_pa - up.inlet_pressure_pa - rise).abs() < 1e-6);
                if patched {
                    assert_eq!(frame.field_exchange, frame.coupled_patch.unwrap());
                    assert_eq!(up.field_exchange, frame.source_patch.unwrap());
                    assert_eq!(owner.zone().boundary_drives(), Some(BACKGROUND));
                    for (side, selected) in [(usize::from(high), frame.field_exchange),
                        (usize::from(!high), up.field_exchange)]
                    {
                        let face = frame.faces[side];
                        let flow = face.normal_volume_flow_m3_s - selected.normal_volume_flow_m3_s;
                        let work = face.acoustic_pressure_work_out_j - selected.acoustic_pressure_work_out_j;
                        let pressure = if side == 0 { -0.25 } else { 0.15 };
                        let expected = pressure * flow * dt;
                        let scale = face.acoustic_pressure_work_out_j.abs()
                            + selected.acoustic_pressure_work_out_j.abs() + expected.abs();
                        assert!((work - expected).abs() <= 256.0 * f64::EPSILON * scale.max(1e-30));
                        complement = complement.max(flow.abs());
                    }
                }
            }
            if patched { assert!(complement > 1e-15); }
            assert_eq!(owner.load().system().accepted_steps(), 3);
            assert_eq!(owner.upstream_load().unwrap().system().accepted_steps(), 3);
        }
    }
}

#[test]
fn finite_supply_and_each_escape_receive_the_original_upstream_field_return() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let downstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let upstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&downstream, &exterior(), &gate);
    let dt = 1.0 / f64::from(ACOUSTIC_RATE);
    let closed = plain_primary(&upstream, dt, 0.03, DgBoundary::Closed);
    let matched = plain_primary(&upstream, dt, 0.03, DgBoundary::Matched);
    let make = |source| coupled(&downstream, &upstream, medium.clone(), &primary, source,
        ThermoviscousDgLimits::default(), &gate);
    let mut reflected = make(&closed);
    let mut escaping = make(&matched);
    assert_ne!(reflected.upstream_load().unwrap().port_response(&gate).unwrap().impedance_pa_s_m3,
        escaping.upstream_load().unwrap().port_response(&gate).unwrap().impedance_pa_s_m3);
    let mut return_difference = 0.0_f64;
    for _ in 0..16 {
        let relation = reflected.upstream_load().unwrap().port_response(&gate).unwrap();
        let drive = NsfDgDrive::FiniteReservoirActuation {
            aperture: 0, volume_rate_m3_s: -2e-7, heat_power_w: 0.01,
            maximum_relative_volume_change: 0.01,
        };
        let a = reflected.step_with_drive(drive, &gate).unwrap();
        let b = escaping.step_with_drive(drive, &gate).unwrap();
        let upstream = a.upstream.as_ref().unwrap();
        let expected = relation.history_pressure_pa + relation.impedance_pa_s_m3 * upstream.inlet_flow_m3_s;
        let scale = expected.abs() + relation.history_pressure_pa.abs()
            + (relation.impedance_pa_s_m3 * upstream.inlet_flow_m3_s).abs();
        assert!((upstream.inlet_pressure_pa - expected).abs() <= 1e-8 * scale.max(1e-20));
        assert_eq!(upstream.source_work_j, 0.0);
        assert!(a.finite_reservoir.as_ref().unwrap().mass_residual_kg.abs() <= 1e-16);
        assert!(a.finite_reservoir.as_ref().unwrap().energy_residual_j.abs() <= 1e-9);
        assert_eq!(a.reservoir_radiations.iter().map(|x| x.resolved_exchange.aperture)
            .collect::<Vec<_>>(), [2, 1]);
        return_difference = return_difference.max((a.inlet_pressure_pa - b.inlet_pressure_pa).abs());
        for (a, b) in a.reservoir_radiations.iter().zip(&b.reservoir_radiations) {
            return_difference = return_difference.max((a.radiation_pressure_pa - b.radiation_pressure_pa).abs());
        }
    }
    assert!(return_difference > 1e-9, "the second field's termination must affect the common resolved solve");
    assert_ne!(reflected.zone().state_hash(), escaping.zone().state_hash());
    assert_ne!(reflected.load().system().state(), escaping.load().system().state());
    assert_ne!(reflected.load().system().load().state(), escaping.load().system().load().state());
    assert_ne!(reflected.finite_reservoir_inventory(), escaping.finite_reservoir_inventory());
}

#[test]
fn late_upstream_and_finite_source_refusals_preserve_both_fields_and_joint_retry() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let downstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let upstream = field_source(&mesh, &unit, 3.0 * dx(&gas), gas);
    let primary = computed_primary(&downstream, &exterior(), &gate);
    let source = computed_primary(&upstream, &exterior(), &gate);
    let make = |limits| coupled(&downstream, &upstream, medium.clone(), &primary, &source, limits, &gate);
    let drive = NsfDgDrive::FiniteReservoirActuation {
        aperture: 0, volume_rate_m3_s: -2e-7, heat_power_w: 0.01,
        maximum_relative_volume_change: 0.01,
    };
    let mut refused = make(ThermoviscousDgLimits {
        maximum_abs_pressure_pa: 1e-12, ..ThermoviscousDgLimits::default()
    });
    let cold = all_owners!(&refused);
    assert!(matches!(refused.step_with_drive(drive, &gate), Err(NsfDgError::Thermoviscous(_))));
    assert_eq!(all_owners!(&refused), cold);
    let mut owner = make(ThermoviscousDgLimits::default());
    owner.step_with_drive(drive, &gate).unwrap();
    let accepted = all_owners!(&owner);
    let impossible = NsfDgDrive::FiniteReservoirActuation {
        aperture: 0, volume_rate_m3_s: -1.0, heat_power_w: 0.01,
        maximum_relative_volume_change: 0.01,
    };
    assert!(owner.validate_drive(impossible).is_ok());
    assert!(matches!(owner.step_with_drive(impossible, &gate), Err(NsfDgError::Input(
        "finite reservoir actuation loses an increment or exceeds its positive bounded volume movement; refine the interval"))));
    assert_eq!(all_owners!(&owner), accepted);
    let cancelled = CancelGate::new_clock_free(); cancelled.request();
    assert!(matches!(owner.step_with_drive(drive, &cancelled), Err(NsfDgError::Cancelled)));
    assert_eq!(all_owners!(&owner), accepted);
    let retry = owner.step_with_drive(drive, &gate).unwrap();
    let mut uninterrupted = make(ThermoviscousDgLimits::default());
    uninterrupted.step_with_drive(drive, &gate).unwrap();
    let expected = uninterrupted.step_with_drive(drive, &gate).unwrap();
    assert_eq!(all_owners!(&owner), all_owners!(&uninterrupted));
    assert_eq!(retry.evaluations, expected.evaluations);
    assert_eq!(retry.incoming_characteristic_m_s.to_bits(), expected.incoming_characteristic_m_s.to_bits());
    assert_eq!(retry.upstream.as_ref().unwrap().incoming_characteristic_m_s.to_bits(),
        expected.upstream.as_ref().unwrap().incoming_characteristic_m_s.to_bits());
    assert_eq!(retry.reservoir_radiations, expected.reservoir_radiations);
}

#[test]
fn upstream_admission_keeps_one_policy_and_the_original_interface_ownership() {
    let gate = CancelGate::new_clock_free();
    let (mesh, unit) = kuhn_cube(1);
    let (medium, gas) = finite_medium(false);
    let side = 3.0 * dx(&gas);
    let field = field_source(&mesh, &unit, side, gas);
    let primary = computed_primary(&field, &exterior(), &gate);
    let dt = 1.0 / f64::from(ACOUSTIC_RATE);
    let source = plain_primary(&field, dt, 0.03, DgBoundary::Closed);
    let make = || super::radiation::system(&field, medium.clone(), &primary, 256, &gate);
    let upstream = || field_load(&field, &source, ThermoviscousDgLimits::default(), &gate);
    assert!(make().validate_drive(NsfDgDrive::UpstreamPressureRise(1.0)).is_err());
    let attached = make().with_upstream_load(upstream(), solve(), &gate).unwrap();
    assert_eq!(attached.joint_solve(), Some(&solve()));
    assert!(attached.validate_drive(NsfDgDrive::Boundary(NsfBoundaryDrive::PressureOffset(1.0))).is_err());
    assert!(attached.validate_drive(NsfDgDrive::UpstreamPressureRise(f64::NAN)).is_err());
    assert!(attached.with_upstream_load(upstream(), solve(), &gate).is_err());
    let mut different = solve(); different.difference_step *= 2.0;
    assert!(make().with_upstream_load(upstream(), solve(), &gate).unwrap()
        .with_reservoir_radiations(vec![(1, placement())], different, &gate).is_err());
    assert!(make().with_reservoir_radiations(vec![(1, placement())], solve(), &gate).unwrap()
        .with_upstream_load(upstream(), different, &gate).is_err());
    let finite = make().with_finite_reservoir(0, 1e-5, finite_limits(), &gate).unwrap()
        .with_upstream_load(upstream(), solve(), &gate).unwrap();
    assert!(finite.validate_drive(NsfDgDrive::UpstreamPressureRise(0.0)).is_err());
    let surface = crate::nsf_dg::reservoir::radiation_surface(make().zone(), gas,
        &primary, 1, placement(), &gate).unwrap().0;
    let unbound_mouth = field_load(&field, &primary, ThermoviscousDgLimits::default(), &gate)
        .with_independent_mouths(&[surface], ACOUSTIC_RATE, &gate).unwrap();
    assert!(make().with_upstream_load(unbound_mouth, solve(), &gate).is_err());
    let changed_gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 120000.0).unwrap();
    for (gas, side, clock) in [(changed_gas, side, dt), (gas, side * 1.01, dt),
        (gas, side, dt * 2.0)]
    {
        let other = field_source(&mesh, &unit, side, gas);
        let duct = plain_primary(&other, clock, 0.03, DgBoundary::Closed);
        assert!(make().with_upstream_load(field_load(&other, &duct,
            ThermoviscousDgLimits::default(), &gate), solve(), &gate).is_err());
    }
    let mut warm = upstream();
    warm.preview_step(1e-8, &gate).unwrap().accept();
    assert!(make().with_upstream_load(warm, solve(), &gate).is_err());
    let mut advanced = make().with_upstream_load(upstream(), solve(), &gate).unwrap();
    advanced.step_with_upstream_pressure_rise(0.0, &gate).unwrap();
    assert!(advanced.with_upstream_load(upstream(), solve(), &gate).is_err());
}
