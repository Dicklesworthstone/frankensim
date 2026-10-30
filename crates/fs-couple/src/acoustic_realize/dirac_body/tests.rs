//! G0/G1/G3: complete histories and physical common-velocity reaction, checked
//! against an independently eliminated mechanical/fluid DAE and time refinement.
use super::*;
use crate::thin_plate::{CompactBody, PlateBank, VkBody};
use fs_material::gas::{GasSpec, GasState};
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_phs::{QuadraticStorage, RelaxationBranch, modal_bank_ports, step_descriptor};

const AREA: [f64; 3] = [0.03, -0.02, 0.045];
const DRIVE: [f64; 3] = [-0.8, 0.6, -0.9];
const STRING_PORT: f64 = 0.7;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap()
}

fn waveguide() -> PortHamiltonian {
    modal_bank_ports(&[12.0], &[0.02], &[&[STRING_PORT], &[1.1]]).unwrap()
}

fn cold_bank(nonlinear: bool) -> PlateBank {
    let mut compact = CompactBody::from_radiator(fs_scenario::RadiatingPlate {
        area_m2: AREA[0],
        mass_kg: 2.0,
        frequency_hz: 18.0 / core::f64::consts::TAU,
        damping_ratio: 0.04,
    })
    .unwrap();
    compact.drive_participation = DRIVE[0];
    let channels = if nonlinear {
        vec![StressChannel {
            coefficient: 2e5,
            coupling: vec![1.0, -0.2, -0.2, 0.6],
        }]
    } else {
        vec![]
    };
    let vk = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![25.0, 37.0],
            channels,
        },
        &[1.4, -0.35, -0.35, 0.8],
        &DRIVE[1..],
        AREA[1..].to_vec(),
    )
    .unwrap()
    .with_relaxation_branches(
        vec![RelaxationBranch {
            projection: vec![0.7, 0.0, -0.4, 0.0],
            stiffness: 16.0,
            relaxation_time_s: 0.08,
        }],
        1,
    )
    .unwrap();
    let mut bank = PlateBank::default();
    bank.linear.push(compact);
    bank.vk = Some(vk);
    bank
}

fn mixed_bank(nonlinear: bool) -> PlateBank {
    let mut bank = cold_bank(nonlinear);
    let gas = gas();
    for _ in 0..64 {
        bank.advance_controlled(0.7, 0.0001, gas.density, 1.0, || false)
            .unwrap();
    }
    bank
}

fn image(bank: &mut PlateBank) -> PlateBankImage {
    bank.physical_image(&gas(), 10_000, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap()
}

fn fluid() -> DiracBodyLoad {
    DiracBodyLoad {
        system: PortHamiltonian::new(
            2,
            2,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 25.0],
            vec![1.0, 1.0, 0.0, 0.0],
            Box::new(QuadraticStorage::new(vec![20_000.0, 0.0, 0.0, 1.0 / 500.0], 2).unwrap()),
        )
        .unwrap(),
        state: vec![0.0003, -0.02],
        face: 1,
        drive_port: Some(0),
    }
}

fn compatible_string(image: &PlateBankImage) -> Vec<f64> {
    vec![0.015, image.system().output(image.state())[0] / STRING_PORT]
}

fn physical(state: &[f64]) -> Vec<f64> {
    let mut x = state[..11].to_vec();
    x[2] /= 2.0_f64.sqrt();
    x[3] /= 2.0_f64.sqrt();
    x.extend_from_slice(&[0.0, 0.0]);
    x
}

fn rhs(x: &[f64], force: f64, source: f64) -> Vec<f64> {
    let pressure = 20_000.0 * x[9];
    let flow = x[10] / 500.0;
    let strain = 4.0 * (0.7 * x[4] - 0.4 * x[6]) - x[8];
    let a_string = -144.0 * x[0] - 0.48 * x[1] + 1.1 * force;
    let a_compact = -324.0 * x[2] - 1.44 * x[3] + AREA[0] * pressure / 2.0;
    let a_first = -625.0 * x[4] - 1.4 * x[5] + 0.35 * x[7] - 2.8 * strain + AREA[1] * pressure;
    let a_second = -1369.0 * x[6] + 0.35 * x[5] - 0.8 * x[7] + 1.6 * strain + AREA[2] * pressure;
    // Differentiate the one physical velocity constraint. This analytical
    // elimination is only the reference, not a production DAE solver.
    let reaction = (DRIVE[0] * a_compact + DRIVE[1] * a_first + DRIVE[2] * a_second
        - STRING_PORT * a_string)
        / (STRING_PORT.powi(2) + DRIVE[0].powi(2) / 2.0 + DRIVE[1].powi(2) + DRIVE[2].powi(2));
    let volume = AREA[0] * x[3] + AREA[1] * x[5] + AREA[2] * x[7];
    vec![
        x[1],
        a_string + STRING_PORT * reaction,
        x[3],
        a_compact - DRIVE[0] * reaction / 2.0,
        x[5],
        a_first - DRIVE[1] * reaction,
        x[7],
        a_second - DRIVE[2] * reaction,
        strain / 0.08,
        flow - volume + source,
        -pressure - 25.0 * flow,
        1.1 * force * x[1] + source * pressure,
        0.48 * x[1].powi(2) + 2.88 * x[3].powi(2) + 1.4 * x[5].powi(2) - 0.7 * x[5] * x[7]
            + 0.8 * x[7].powi(2)
            + strain.powi(2) / 0.08
            + 25.0 * flow.powi(2),
    ]
}

fn rk4(x: &mut [f64], dt: f64, force: f64, source: f64) {
    let a = rhs(x, force, source);
    let stage = |k: &[f64], fraction: f64| {
        x.iter()
            .zip(k)
            .map(|(x, k)| x + dt * fraction * k)
            .collect::<Vec<_>>()
    };
    let b = rhs(&stage(&a, 0.5), force, source);
    let c = rhs(&stage(&b, 0.5), force, source);
    let d = rhs(&stage(&c, 1.0), force, source);
    for i in 0..x.len() {
        x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
    }
}

#[test]
fn mixed_history_reaction_and_observation_match_the_physical_dae() {
    let mut bank = mixed_bank(false);
    let original = image(&mut bank);
    let accepted = original.state().to_vec();
    assert_eq!(accepted.len(), 7);
    assert!(accepted[6].abs() > 1e-8, "warm Maxwell history is present");
    let string = compatible_string(&original);
    let joined = join_complete_body(waveguide(), string.clone(), original, Some(fluid())).unwrap();
    assert_eq!(joined.system.state_dim(), 12);
    assert_eq!(joined.system.differential_dim(), 11);
    assert_eq!(joined.fluid_port, Some(1));
    assert_eq!(joined.system.port_dim(), 2);
    assert_eq!(&joined.state[..2], &string);
    assert_eq!(&joined.state[2..9], &accepted);
    assert_eq!(&joined.state[9..11], &[0.0003, -0.02]);
    let duration = 0.01;
    let (force, source) = (0.3, 2e-4);
    let mut reference = physical(&joined.state);
    for _ in 0..8192 {
        rk4(&mut reference, duration / 8192.0, force, source);
    }
    let energy = joined.system.hamiltonian(&joined.state);
    let mut errors = Vec::new();
    for steps in [20, 40, 80] {
        let dt = duration / steps as f64;
        let rate = fs_phs::linear_exact::descriptor_pole_rate_bound(&joined.system).unwrap();
        fs_phs::linear_exact::admit_midpoint_band(rate, dt).unwrap();
        let mut state = joined.state.clone();
        let (mut work, mut heat) = (0.0, 0.0);
        for _ in 0..steps {
            let record = step_descriptor(&joined.system, &state, &[force, source], dt).unwrap();
            let before = physical(&state);
            let after = physical(&record.x);
            let pressure = gas().density
                * (AREA[0] * (after[3] - before[3])
                    + AREA[1] * (after[5] - before[5])
                    + AREA[2] * (after[7] - before[7]))
                / (2.0 * core::f64::consts::PI * 1.3 * dt);
            let observed = joined
                .observation
                .pressure_between(&state, &record.x, dt, gas().density, 1.3)
                .unwrap();
            assert!((pressure - observed).abs() < 1e-11);
            let mismatch = STRING_PORT * after[1]
                - DRIVE[0] * after[3]
                - DRIVE[1] * after[5]
                - DRIVE[2] * after[7];
            assert!(mismatch.abs() < 1e-10);
            work += record.supplied;
            heat += record.dissipated;
            state = record.x;
        }
        let actual = physical(&state);
        let weights = [
            144.0, 1.0, 648.0, 2.0, 625.0, 1.0, 1369.0, 1.0, 1.0, 20_000.0, 0.002,
        ];
        let error = actual
            .iter()
            .zip(&reference)
            .zip(weights)
            .map(|((a, b), w)| w * (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        errors.push(error);
        assert!((joined.system.hamiltonian(&state) - energy + heat - work).abs() < 1e-11);
        if steps == 80 {
            assert!(error < 2e-7, "state error {error:e}");
            assert!((work - reference[11]).abs() < 2e-9);
            assert!((heat - reference[12]).abs() < 2e-9);
        }
    }
    assert!(errors[0] / errors[1] > 3.8);
    assert!(errors[1] / errors[2] > 3.8);
    assert_eq!(image(&mut bank).state(), accepted);
}

#[test]
fn nonlinear_history_and_open_surface_are_preserved_and_bad_initial_joins_refuse() {
    let mut bank = mixed_bank(true);
    let source = image(&mut bank);
    assert!(!source.system().storage().quadratic_part().unwrap().exact);
    let original = source.state().to_vec();
    let initial = compatible_string(&source);
    let energy = waveguide().hamiltonian(&initial) + source.system().hamiltonian(&original);
    let joined = join_complete_body(waveguide(), initial, source, None).unwrap();
    assert_eq!(joined.fluid_port, None);
    assert_eq!(joined.system.port_dim(), 2); // String force and unforced body pressure.
    assert_eq!(&joined.state[2..9], &original);
    assert_eq!(joined.system.hamiltonian(&joined.state), energy);
    assert!(!joined.system.storage().quadratic_part().unwrap().exact);
    let record = step_descriptor(&joined.system, &joined.state, &[0.0, 0.0], 0.0001).unwrap();
    assert!((record.delta_h + record.dissipated).abs() < 1e-12);
    assert!(
        joined
            .observation
            .pressure_between(&joined.state[..9], &record.x, 0.0001, gas().density, 1.0,)
            .is_err()
    );
    assert!(join_complete_body(waveguide(), vec![0.015, 0.0], image(&mut bank), None).is_err());
    let source = image(&mut bank);
    let initial = compatible_string(&source);
    let mut bad = fluid();
    bad.drive_port = Some(bad.face);
    assert!(join_complete_body(waveguide(), initial, source, Some(bad)).is_err());
    assert_eq!(image(&mut bank).state(), original);
}

#[test]
fn attached_radiation_history_stays_in_the_descriptor_before_its_multiplier() {
    let gas = gas();
    let mut bank = PlateBank::default();
    bank.linear.push(
        CompactBody::from_radiator(fs_scenario::RadiatingPlate {
            area_m2: core::f64::consts::PI * 0.01_f64.powi(2),
            mass_kg: 2e-5,
            frequency_hz: 600.0 / core::f64::consts::TAU,
            damping_ratio: 0.01,
        })
        .unwrap(),
    );
    bank.declare_radiation_loads(&gas, 48_000).unwrap();
    for _ in 0..12 {
        bank.advance_controlled(2e-4, 1.0 / 48_000.0, gas.density, 1.0, || false)
            .unwrap();
    }
    let source = bank
        .physical_image(&gas, 48_000, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap();
    let accepted = source.state().to_vec();
    assert!(accepted.len() > 2 && accepted[2..].iter().any(|x| *x != 0.0));
    let initial = compatible_string(&source);
    let energy = waveguide().hamiltonian(&initial) + source.system().hamiltonian(&accepted);
    let joined = join_complete_body(waveguide(), initial, source, None).unwrap();
    assert_eq!(&joined.state[2..2 + accepted.len()], &accepted);
    assert_eq!(joined.system.differential_dim(), 2 + accepted.len());
    assert_eq!(joined.state.len(), 3 + accepted.len());
    assert_eq!(joined.system.hamiltonian(&joined.state), energy);
    let rate = fs_phs::linear_exact::descriptor_pole_rate_bound(&joined.system).unwrap();
    assert!(
        fs_phs::linear_exact::admit_midpoint_band(rate, 1.0 / 48_000.0).is_err(),
        "retained stiff acoustic history must not bypass original temporal admission"
    );
    let mut exact = fs_phs::linear_exact::PreparedCommonFlowZoh::new_owned(
        joined.system,
        1.0 / 48_000.0,
        fs_phs::linear_exact::LinearSplitBudget {
            max_modes: joined.state.len(),
            max_omega_dt: core::f64::consts::PI,
        },
    )
    .unwrap();
    let mut state = joined.state;
    let mut next = vec![0.0; state.len()];
    let input = vec![0.0; exact.port_dim()];
    let mut output = vec![0.0; input.len()];
    let mut heat = 0.0;
    for _ in 0..128 {
        let record = exact
            .step_into(&state, &input, &mut next, &mut output)
            .unwrap();
        heat += record.dissipated;
        core::mem::swap(&mut state, &mut next);
    }
    assert_eq!(state.len(), 3 + accepted.len());
    assert!(state.iter().all(|value| value.is_finite()));
    assert!(heat > 0.0);
    assert!(exact.hamiltonian(&state) < energy);
    assert!((exact.hamiltonian(&state) - energy + heat).abs() < 2e-10 * energy);
    assert_eq!(
        bank.physical_image(&gas, 48_000, &fs_exec::CancelGate::default())
            .unwrap()
            .unwrap()
            .state(),
        accepted,
    );
}

fn string(frequency_hz: f64) -> fs_scenario::PrestressedString {
    fs_scenario::PrestressedString {
        length_m: 0.5,
        tension_n: (2.0 * frequency_hz).powi(2),
        lin_density_kg_m: 1.0,
        axial_stiffness_n: 0.0,
        width_m: 0.001,
        n_modes: 1,
        damping_ratio: 0.0,
        rayleigh: None,
        bending_stiffness_n_m2: 0.0,
        kelvin_voigt_bending: None,
        relaxation_bending: None,
        polarization_detune: 0.0,
        moving_end: true,
        tension_loss_factor: None,
    }
}

fn render(
    bank: &mut PlateBank,
    string: &fs_scenario::PrestressedString,
    rate: u32,
    samples: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    super::super::realize_dirac_join(
        string,
        Some(fs_scenario::Pluck {
            station_frac: 0.5,
            height_m: 1e-3,
        }),
        None,
        &[],
        None,
        bank,
        None,
        None,
        None,
        None,
        &gas(),
        1.0,
        rate,
        samples,
    )
}

#[test]
fn actual_dirac_renderer_uses_mixed_and_vk_only_banks_with_signed_surface_maps() {
    let string = string(4.0);
    let mixed = render(&mut cold_bank(false), &string, 48_000, 96).unwrap();
    let mut compact = cold_bank(false);
    compact.vk = None;
    let compact = render(&mut compact, &string, 48_000, 96).unwrap();
    let mut only_vk = cold_bank(false);
    only_vk.linear.clear();
    let only_vk = render(&mut only_vk, &string, 48_000, 96).unwrap();
    let mut reversed = cold_bank(false);
    reversed.vk.as_mut().unwrap().areas[1] *= -1.0;
    let reversed = render(&mut reversed, &string, 48_000, 96).unwrap();
    let difference = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
            / a.iter().map(|a| a * a).sum::<f64>().sqrt()
    };
    assert!(mixed.iter().chain(&only_vk).all(|p| p.is_finite()));
    assert!(only_vk.iter().any(|p| p.abs() > 1e-6));
    assert!(difference(&mixed, &compact) > 1e-3);
    assert!(difference(&mixed, &only_vk) > 1e-3);
    assert!(difference(&mixed, &reversed) > 1e-3);
    assert!(
        matches!(
            render(&mut mixed_bank(false), &string, 48_000, 1),
            Err(AcousticRealizeError::InvalidDescription { .. })
        ),
        "an incompatible moving source cannot be silently reset by renderer dispatch"
    );
    assert!(matches!(
        super::super::realize_dirac_join(
            &string,
            None,
            None,
            &[],
            None,
            &mut cold_bank(false),
            None,
            None,
            None,
            Some(fs_scenario::VolumeVelocityPulse {
                peak_m3_s: 1e-4,
                duration_s: 0.01,
            }),
            &gas(),
            1.0,
            48_000,
            1,
        ),
        Err(AcousticRealizeError::InvalidDescription { .. })
    ));
}

#[test]
fn actual_dirac_renderer_resolves_quadratic_clocks_and_keeps_nonlinear_admission() {
    for rate in [48_000, 96_000, 192_000] {
        assert!(render(&mut cold_bank(false), &string(1.0), rate, 1).is_ok());
        let high = render(&mut cold_bank(false), &string(10_000.0), rate, 8).unwrap();
        assert!(high.iter().all(|sample| sample.is_finite()));
        assert!(high.iter().any(|sample| sample.abs() > 1e-6));
        assert!(matches!(
            render(&mut cold_bank(true), &string(10_000.0), rate, 1),
            Err(AcousticRealizeError::TimeResolution(_))
        ));
    }
}

#[test]
fn actual_dirac_renderer_retains_declared_radiation_at_the_audio_clock() {
    let mut bank = PlateBank::default();
    bank.linear.push(
        CompactBody::from_radiator(fs_scenario::RadiatingPlate {
            area_m2: core::f64::consts::PI * 0.01_f64.powi(2),
            mass_kg: 2e-5,
            frequency_hz: 600.0 / core::f64::consts::TAU,
            damping_ratio: 0.01,
        })
        .unwrap(),
    );
    bank.declare_radiation_loads(&gas(), 48_000).unwrap();
    let before = bank
        .physical_image(&gas(), 48_000, &fs_exec::CancelGate::default())
        .unwrap()
        .unwrap()
        .state()
        .to_vec();
    assert!(
        before.len() > 2,
        "the physical exterior has internal memory"
    );
    let long = render(&mut bank, &string(110.0), 48_000, 96).unwrap();
    let short = render(&mut bank, &string(110.0), 48_000, 8).unwrap();
    assert_eq!(&long[..short.len()], short);
    assert!(long.iter().all(|sample| sample.is_finite()));
    assert!(long.iter().any(|sample| sample.abs() > 1e-8));
    assert_eq!(
        bank.physical_image(&gas(), 48_000, &fs_exec::CancelGate::default())
            .unwrap()
            .unwrap()
            .state(),
        before,
    );
}
