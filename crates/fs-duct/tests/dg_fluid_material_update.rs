//! Continuing real-fluid states use the native DG state and energy owner.
//! The explicit material fixture below is manufactured constitutive data,
//! not a measured liquid or a substitute for the original IF97 regression.

use fs_duct::dg::{
    DgBoundary, DgDuct, DgFluidMaterial, DgMedium, DgMesh, DgSidePort, DgState, PortId, PortInput,
};
use fs_duct::{Duct, Segment};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::WallPin;
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::regions::If97FluidState;
use fs_thermochem::water::if97::{If97Region, If97State};
use fs_thermochem::water::transport::If97TransportState;

fn shape(radius: f64, length: f64) -> Duct {
    Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    }
}
fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 5e-9 * scale.max(1e-30),
        "{actual:e} != {expected:e}"
    );
}
fn energy(values: &[f64]) -> f64 {
    0.5 * values.iter().map(|x| x * x).sum::<f64>()
}
fn coordinates<M: DgMedium>(duct: &DgDuct<M>, state: &DgState) -> Vec<f64> {
    duct.continuous_operator()
        .energy_coordinates(state)
        .unwrap()
        .to_vec()
}
fn water(pressure: f64, temperature: f64) -> If97TransportState {
    let state = If97State::at_pressure_temperature(
        Pressure::new(pressure),
        Temperature::new(temperature),
        If97Region::Region1,
    )
    .unwrap();
    If97TransportState::from_state(If97FluidState::Gibbs(state)).unwrap()
}

#[test]
fn pressure_refresh_of_real_water_continues_the_original_galerkin_wave_and_work() {
    let old = water(3e6, 300.0);
    let new = water(40e6, 300.0);
    let previous = old.acoustic_state().unwrap();
    let replacement = new.acoustic_state().unwrap();
    let (radius, length) = (0.01, 0.5);
    let area = core::f64::consts::PI * radius * radius;
    let mut duct = DgDuct::from_material(
        &shape(radius, length),
        &old,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    let dt = 1e-6;
    let old_propagator = duct.propagator(dt, 1e-13).unwrap();
    let mut state = duct.project_initial_state(|x| 1000.0 * (2.0 * x / length - 1.0), |_| 1e-9);
    for _ in 0..10 {
        old_propagator.step(&mut state, &[]).unwrap();
    }
    old_propagator.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let old_p = before.pressure()[1];
    let old_q = before.volume_flow()[0];
    let old_bulk = previous.adiabatic_bulk_modulus_pa();
    let new_bulk = replacement.adiabatic_bulk_modulus_pa();
    let old_rho = previous.density_kg_m3();
    let new_rho = replacement.density_kg_m3();
    let work = 0.5
        * area
        * length
        * ((new_bulk - old_bulk) * (old_p / old_bulk).powi(2) / 3.0
            + (new_rho.recip() - old_rho.recip()) * (old_rho * old_q / area).powi(2));
    let receipt = duct
        .rebind_uniform_material(&mut state, &new, 1.0, &mut || false)
        .unwrap();
    assert_eq!(duct.material(), &new);
    assert_eq!(duct.material().thermodynamics().pressure().value(), 40e6);
    assert_eq!(duct.acoustic_material(), replacement);
    near(receipt.parameter_work_j, work, before.stored_energy_j());
    near(
        receipt.energy_after_j - receipt.energy_before_j,
        work,
        before.stored_energy_j(),
    );
    assert_eq!(state.step_count(), 10);
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(
        state.cumulative_dissipated_j(),
        before.cumulative_dissipated_j()
    );
    let p0 = old_p * new_bulk / old_bulk;
    let q0 = old_q * old_rho / new_rho;
    near(state.pressure()[0], -p0, 1000.0);
    near(state.pressure()[1], p0, 1000.0);
    for &q in state.volume_flow() {
        near(q, q0, 1e-8);
    }
    let retained = state.clone();
    assert!(old_propagator.step(&mut state, &[]).is_err());
    assert_eq!(state, retained);

    // Independent two-coordinate Galerkin equations for this linear element:
    // P'=6 K_s Q/(A L), Q'=-2 A P/(rho L), omega=sqrt(12)*c/L.
    let omega = 12.0_f64.sqrt() * replacement.sound_speed_m_s() / length;
    let a = 6.0 * new_bulk / (area * length);
    let b = 2.0 * area / (new_rho * length);
    let propagator = duct.propagator(dt, 1e-13).unwrap();
    for step in 1..=30 {
        let output = propagator.step(&mut state, &[]).unwrap();
        propagator.sync_nodal(&mut state).unwrap();
        let (sin, cos) = (omega * dt * step as f64).sin_cos();
        let expected_p = p0 * cos + a * q0 / omega * sin;
        let expected_q = q0 * cos - b * p0 / omega * sin;
        near(state.pressure()[0], -expected_p, 1000.0);
        near(state.pressure()[1], expected_p, 1000.0);
        for &q in state.volume_flow() {
            near(q, expected_q, 1e-8);
        }
        assert_eq!(
            output.ledger.cumulative_parameter_work_j,
            receipt.parameter_work_j
        );
        near(
            output.ledger.energy_balance_residual_j,
            0.0,
            before.initial_energy_j(),
        );
    }
}

// Complete original query and material identity stay on the generic owner.
// These labels mark test data, not source authority or a named material law.
#[derive(Clone, Debug, PartialEq)]
struct FixtureMaterial {
    identity: u32,
    query: [f64; 2],
    state: ThermoviscousFluidState,
    different_acoustic_view: Option<FluidAcousticState>,
}
impl FluidAcoustics for FixtureMaterial {
    fn acoustic_state(&self) -> Result<FluidAcousticState, MaterialError> {
        Ok(self
            .different_acoustic_view
            .unwrap_or(self.state.acoustic()))
    }
}
impl FluidTransport for FixtureMaterial {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.state.transport())
    }
}
impl ThermoviscousFluid for FixtureMaterial {
    fn thermoviscous_state(&self) -> Result<ThermoviscousFluidState, MaterialError> {
        Ok(self.state)
    }
}
fn material(
    identity: u32,
    rho: f64,
    speed: f64,
    mu: f64,
    k: f64,
    cp: f64,
    cv: f64,
) -> FixtureMaterial {
    FixtureMaterial {
        identity,
        query: [
            290.0 + 10.0 * f64::from(identity),
            1e5 * f64::from(identity),
        ],
        state: ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(rho, speed).unwrap(),
            FluidTransportState::try_new(rho, cp, k, mu).unwrap(),
            cv,
        )
        .unwrap(),
        different_acoustic_view: None,
    }
}

#[test]
fn generic_refresh_preserves_live_diffusion_wall_history_ports_and_source_identity() {
    let old = material(1, 1000.0, 1500.0, 0.001, 0.6, 4200.0, 4100.0);
    let new = material(2, 850.0, 1300.0, 0.0015, 0.45, 3800.0, 3600.0);
    let pins = [
        WallPin {
            surface_density: 20.0,
            stiffness_per_area: 2e13,
            resistance: 2e6,
        },
        WallPin {
            surface_density: 30.0,
            stiffness_per_area: 3e13,
            resistance: 3e6,
        },
    ];
    let laws = [
        GeneralizedMaxwell::new(1e6, vec![(2e5, 2e-6)]).unwrap(),
        GeneralizedMaxwell::new(4e6, vec![(5e5, 5e-6), (9e5, 1e-6)]).unwrap(),
    ];
    let side = DgSidePort {
        id: PortId(7),
        downstream_element: 1,
    };
    let build = |source: &FixtureMaterial| {
        DgDuct::from_material(
            &shape(1e-4, 0.05),
            source,
            DgMesh {
                elements: 2,
                order: 1,
            },
            DgBoundary::Matched,
            DgBoundary::Port(PortId(1)),
        )
        .unwrap()
        .with_side_ports(&[side])
        .unwrap()
        .with_viscothermal(4, [20.0, 60.0], 0.001)
        .unwrap()
        .with_compliant_wall_relaxation_regions(&pins, &[Some(&laws[0]), Some(&laws[1])])
        .unwrap()
    };
    let mut duct = build(&old);
    assert_eq!(duct.thermoviscous_material(), Some(old.state));
    let dt = 1e-7;
    let propagator = duct.propagator(dt, 1e-12).unwrap();
    let mut state = duct.project_initial_state(|x| 100.0 + 10.0 * x, |_| 1e-10);
    for _ in 0..20 {
        propagator.step(&mut state, &[]).unwrap();
    }
    propagator.sync_nodal(&mut state).unwrap();
    let before = state.clone();
    let original = coordinates(&duct, &state);
    let original_ports = duct.continuous_operator().ports().to_vec();
    let wall = duct.compliant_wall_state(&state).unwrap();
    let old_model = duct.viscothermal_model().unwrap().clone();
    let m = duct.field_dofs();
    let series_end = m * (2 + old_model.series_terms().len());
    let prefix = series_end + m * old_model.shunt_terms().len();
    assert!(energy(&original[2 * m..series_end]) > 0.0);
    assert!(energy(&original[series_end..prefix]) > 0.0);
    assert!(wall.relaxation_stored_energy_j > 0.0);
    let old_bulk = old.state.acoustic().adiabatic_bulk_modulus_pa();
    let new_bulk = new.state.acoustic().adiabatic_bulk_modulus_pa();
    let pressure_factor = new_bulk / old_bulk;
    let momentum_factor =
        old.state.acoustic().density_kg_m3() / new.state.acoustic().density_kg_m3();
    let thermal_factor = pressure_factor * (old.state.heat_capacity_ratio() - 1.0)
        / (new.state.heat_capacity_ratio() - 1.0);
    let expected_work = (pressure_factor - 1.0) * energy(&original[..m])
        + (momentum_factor - 1.0) * energy(&original[m..series_end])
        + (thermal_factor - 1.0) * energy(&original[series_end..prefix]);
    let receipt = duct
        .rebind_uniform_viscothermal_material(&mut state, &new, 0.001, 1.0, &mut || false)
        .unwrap();
    assert_eq!(duct.material(), &new);
    assert_eq!(duct.thermoviscous_material(), Some(new.state));
    assert_eq!(duct.side_ports(), &[side]);
    assert_eq!(duct.viscothermal_band_rad_s(), Some(old_model.band_rad_s()));
    assert_eq!(state.step_count(), before.step_count());
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(
        state.cumulative_dissipated_j(),
        before.cumulative_dissipated_j()
    );
    let updated = coordinates(&duct, &state);
    near(
        receipt.parameter_work_j,
        expected_work,
        before.stored_energy_j(),
    );
    near(
        receipt.energy_after_j - receipt.energy_before_j,
        expected_work,
        before.stored_energy_j(),
    );
    assert_eq!(&updated[prefix..], &original[prefix..]);
    // Physical arm momentum / compression, up to the unchanged geometry root.
    for (range, old_storage, new_storage) in [
        (
            m..series_end,
            old.state.acoustic().density_kg_m3(),
            new.state.acoustic().density_kg_m3(),
        ),
        (
            series_end..prefix,
            (old.state.heat_capacity_ratio() - 1.0) / old_bulk,
            (new.state.heat_capacity_ratio() - 1.0) / new_bulk,
        ),
    ] {
        for i in range {
            near(
                updated[i] * new_storage.sqrt(),
                original[i] * old_storage.sqrt(),
                original[i].abs() * old_storage.sqrt(),
            );
        }
    }
    let after_wall = duct.compliant_wall_state(&state).unwrap();
    assert_eq!(after_wall.stored_energy_j, wall.stored_energy_j);
    assert_eq!(after_wall.dissipated_power_w, wall.dissipated_power_w);
    for (a, b) in after_wall.displacement_m.iter().zip(&wall.displacement_m) {
        near(*a, *b, b.abs().max(1e-18));
    }
    for (a, b) in after_wall.velocity_m_s.iter().zip(&wall.velocity_m_s) {
        near(*a, *b, b.abs().max(1e-12));
    }
    for (a, b) in after_wall
        .relaxation_pressure_pa
        .iter()
        .flatten()
        .zip(wall.relaxation_pressure_pa.iter().flatten())
    {
        near(*a, *b, b.abs().max(1e-12));
    }
    for ((id, column), (old_id, old_column)) in duct
        .continuous_operator()
        .ports()
        .iter()
        .zip(&original_ports)
    {
        assert_eq!(id, old_id);
        let effort: f64 = column.iter().zip(&updated).map(|(b, x)| b * x).sum();
        let old_effort: f64 = old_column.iter().zip(&original).map(|(b, x)| b * x).sum();
        near(effort, pressure_factor * old_effort, 1000.0);
    }
    let cold = build(&new);
    for (a, b) in duct
        .continuous_operator()
        .generator_per_s()
        .iter()
        .zip(cold.continuous_operator().generator_per_s())
    {
        near(*a, *b, b.abs().max(1.0));
    }
    let retained = state.clone();
    assert!(propagator.step(&mut state, &[]).is_err());
    assert_eq!(state, retained);
    let continued = duct.propagator(dt, 1e-12).unwrap();
    let inputs = [PortInput {
        id: side.id,
        effort: 1e-11,
    }];
    for _ in 0..20 {
        let output = continued.step(&mut state, &inputs).unwrap();
        assert_eq!(
            output.ledger.cumulative_parameter_work_j,
            receipt.parameter_work_j
        );
        near(
            output.ledger.energy_balance_residual_j,
            0.0,
            before.initial_energy_j(),
        );
    }

    // An otherwise identical query/source refresh remains attached even when
    // the numeric operator is reusable and no parameter work is owed.
    let original_id = duct.continuous_operator().operator_identity();
    let frozen = coordinates(&duct, &state);
    let mut source_only = new.clone();
    source_only.identity = 3;
    source_only.query = [325.0, 3e5];
    let same = duct
        .rebind_uniform_viscothermal_material(&mut state, &source_only, 0.001, 0.0, &mut || false)
        .unwrap();
    assert_eq!(duct.material(), &source_only);
    assert_eq!(duct.continuous_operator().operator_identity(), original_id);
    assert_eq!(coordinates(&duct, &state), frozen);
    assert_eq!(same.parameter_work_j, 0.0);

    // Transport-only data changes retain all instantaneous storage and change
    // subsequent physical diffusion, rather than scaling a rendered waveform.
    let transport_only = material(4, 850.0, 1300.0, 0.003, 1.35, 3800.0, 3600.0);
    let model_before = duct.viscothermal_model().unwrap().clone();
    let frozen = coordinates(&duct, &state);
    let mut control = state.clone();
    let old_continuation = duct.propagator(dt, 1e-12).unwrap();
    let change = duct
        .rebind_uniform_viscothermal_material(&mut state, &transport_only, 0.001, 0.0, &mut || {
            false
        })
        .unwrap();
    assert_eq!(change.parameter_work_j, 0.0);
    assert_eq!(coordinates(&duct, &state), frozen);
    assert_eq!(duct.thermoviscous_material(), Some(transport_only.state));
    let model_after = duct.viscothermal_model().unwrap();
    near(
        model_after.series_terms()[0].rate_per_s,
        2.0 * model_before.series_terms()[0].rate_per_s,
        model_before.series_terms()[0].rate_per_s,
    );
    near(
        model_after.shunt_terms()[0].rate_per_s,
        3.0 * model_before.shunt_terms()[0].rate_per_s,
        model_before.shunt_terms()[0].rate_per_s,
    );
    let refreshed = duct.propagator(dt, 1e-12).unwrap();
    for _ in 0..50 {
        refreshed.step(&mut state, &[]).unwrap();
        old_continuation.step(&mut control, &[]).unwrap();
    }
    refreshed.sync_nodal(&mut state).unwrap();
    old_continuation.sync_nodal(&mut control).unwrap();
    assert!(
        state
            .pressure()
            .iter()
            .zip(control.pressure())
            .any(|(a, b)| (a - b).abs() > 1e-9)
    );
}

#[test]
fn failed_live_material_changes_preserve_source_operator_history_and_clock() {
    let old = material(1, 1000.0, 1500.0, 0.001, 0.6, 4200.0, 4100.0);
    let new = material(2, 850.0, 1300.0, 0.0015, 0.45, 3800.0, 3600.0);
    let mut duct = DgDuct::from_material(
        &shape(1e-4, 0.05),
        &old,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
    .with_viscothermal(4, [20.0, 60.0], 0.001)
    .unwrap();
    let mut state = duct.project_initial_state(|_| 100.0, |_| 1e-10);
    let propagator = duct.propagator(1e-7, 1e-12).unwrap();
    for _ in 0..4 {
        propagator.step(&mut state, &[]).unwrap();
    }
    let accepted = state.clone();
    let identity = duct.continuous_operator().operator_identity();
    let generator = duct.continuous_operator().generator_per_s().to_vec();
    let model = duct.viscothermal_model().unwrap().clone();
    let same = |duct: &DgDuct<DgFluidMaterial<FixtureMaterial>>, state: &DgState| {
        assert_eq!(state, &accepted);
        assert_eq!(duct.material(), &old);
        assert_eq!(duct.thermoviscous_material(), Some(old.state));
        assert_eq!(duct.continuous_operator().operator_identity(), identity);
        assert_eq!(
            duct.continuous_operator().generator_per_s(),
            generator.as_slice()
        );
        assert_eq!(duct.viscothermal_model(), Some(&model));
    };
    assert!(
        duct.rebind_uniform_material(&mut state, &new, 1.0, &mut || false)
            .is_err()
    );
    same(&duct, &state);
    assert!(
        duct.rebind_uniform_viscothermal_material(&mut state, &new, 0.001, 0.0, &mut || false,)
            .is_err()
    );
    same(&duct, &state);
    let mut cancelled_calls = 0;
    assert!(
        duct.rebind_uniform_viscothermal_material(&mut state, &new, 0.001, 1.0, &mut || {
            cancelled_calls += 1;
            cancelled_calls == 5
        },)
            .is_err()
    );
    assert_eq!(cancelled_calls, 5);
    same(&duct, &state);
    let no_thermal = material(2, 850.0, 1300.0, 0.0015, 0.45, 3800.0, 3800.0);
    assert!(
        duct.rebind_uniform_viscothermal_material(&mut state, &no_thermal, 0.001, 1.0, &mut || {
            false
        },)
            .is_err()
    );
    same(&duct, &state);
    let mut mismatched = new.clone();
    mismatched.different_acoustic_view = Some(old.state.acoustic());
    assert!(
        duct.rebind_uniform_viscothermal_material(&mut state, &mismatched, 0.001, 1.0, &mut || {
            false
        },)
            .is_err()
    );
    same(&duct, &state);
    assert!(
        duct.rebind_uniform_viscothermal_material(&mut state, &new, 0.0, 1.0, &mut || false,)
            .is_err()
    );
    same(&duct, &state);
    // The accepted old propagator still owns the original state after refusals.
    propagator.step(&mut state, &[]).unwrap();
    assert_eq!(state.step_count(), accepted.step_count() + 1);
}
