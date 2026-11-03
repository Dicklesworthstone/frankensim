//! G1/G3: continuous wave scattering through the actual complete modal line.
//! Synthetic physical sources; the axial ladder does not validate real materials.
use core::f64::consts::{PI, TAU};
use fs_couple::mm_line::continuous::{ContinuousModalLine, ModalAdmittanceLoad};
use fs_duct::RadiusProfile;
use fs_duct::modal::ModeCount;
use fs_duct::modal::time_domain::ModalEnd;
use fs_material::gas::{GasSpec, GasState};
use fs_math::{c64::C64, det};
use fs_phs::harmonic::HarmonicBudget;
use fs_phs::harmonic::scattering::{ScatteringBudget, ScatteringResponse};
use fs_phs::linear_exact::LinearSplitBudget;
use fs_phs::{PortHamiltonian, QuadraticStorage};

fn gas(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101_325.).unwrap()
}
fn time_budget() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 124,
        max_omega_dt: 1.,
    }
}
fn harmonic_budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 124,
        max_dense_bytes: 16 * 1024 * 1024,
        max_dense_work: 4_000_000,
        relative_residual_tolerance: 2e-10,
    }
}
fn scattering_budget() -> ScatteringBudget {
    ScatteringBudget {
        maximum_ports: 4,
        maximum_payload_bytes: 1024 * 1024,
        maximum_work: 1_000_000,
    }
}
fn query(line: &ContinuousModalLine, omega: f64, input: &[C64]) -> ScatteringResponse {
    line.prepare_harmonic_scattering(harmonic_budget(), scattering_budget(), || false)
        .unwrap()
        .factor(omega, || false)
        .unwrap()
        .response(input, || false)
        .unwrap()
}
fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}
fn snapshot(line: &ContinuousModalLine) -> Vec<u64> {
    let mut values = line.state().to_vec();
    values.extend_from_slice(line.incoming());
    values.extend_from_slice(line.mean_pressures());
    for p in line.port_accounts() {
        values.extend([
            p.boundary_work,
            p.incident_energy,
            p.outgoing_energy,
            p.outgoing_wave,
        ]);
    }
    values.into_iter().map(f64::to_bits).collect()
}

#[test]
fn g1_cylinder_reflection_and_transmission_converge_on_original_two_mode_spatial_ladder() {
    let medium = gas(293.15);
    let (radius, length) = (0.025, 0.12);
    let profile = RadiusProfile::exponential(radius, 0., length).unwrap();
    let omega = TAU * 600.;
    let phase = omega * length / medium.sound_speed;
    let travel = C64::new(det::cos(phase), det::sin(phase));
    let z0 = medium.density * medium.sound_speed / (PI * radius * radius);
    for end in [ModalEnd::Ports, ModalEnd::Rigid, ModalEnd::PressureRelease] {
        let mut errors = Vec::new();
        for elements in [6, 12, 24] {
            let line = ContinuousModalLine::new(
                &profile,
                &medium,
                ModeCount(2),
                elements,
                end,
                1. / 96_000.,
                time_budget(),
                || false,
            )
            .unwrap();
            let before = snapshot(&line);
            let model = line
                .prepare_harmonic_scattering(harmonic_budget(), scattering_budget(), || false)
                .unwrap();
            assert_eq!(model.dimension(), line.state().len());
            assert_eq!(model.wave_ports().len(), line.port_count());
            assert!(model.wave_ports().iter().all(|&(_, z)| z == z0));
            let factor = model.factor(omega, || false).unwrap();
            let mut input = vec![C64::ZERO; model.port_count()];
            input[0] = C64::ONE;
            let result = factor.response(&input, || false).unwrap();
            assert_eq!(result.physical.state.len(), line.state().len());
            assert_eq!(result.physical.dissipated_power_w, 0.);
            let power = 1. / (2. * z0);
            close(result.wave_power_residual_w, 0., power, 2e-9);
            close(result.physical.power_residual_w, 0., power, 2e-9);
            for wave in &result.waves {
                close(wave.power_residual_w, 0., power, 2e-13);
            }
            // A uniform cylinder has no radial-mode mixing. Both retained
            // mode families remain in the original solve; the planar drive
            // leaves the second family unexcited by the physical symmetry.
            assert!(result.waves[1].outgoing.abs() < 1e-11);
            let error = match end {
                ModalEnd::Ports => {
                    assert!(result.waves[3].outgoing.abs() < 1e-11);
                    // Reuse the identical factor for excitation from the other
                    // physical face; equal references preserve reciprocity.
                    input[0] = C64::ZERO;
                    input[2] = C64::ONE;
                    let reverse = factor.response(&input, || false).unwrap();
                    assert!((reverse.waves[0].outgoing - result.waves[2].outgoing).abs() < 2e-11);
                    result.waves[0].outgoing.abs() + (result.waves[2].outgoing - travel).abs()
                }
                ModalEnd::Rigid => (result.waves[0].outgoing - travel * travel).abs(),
                ModalEnd::PressureRelease => (result.waves[0].outgoing + travel * travel).abs(),
            };
            errors.push(error);
            assert_eq!(snapshot(&line), before);
        }
        assert!(errors[0] > 1e-6, "live axial error: {end:?}, {errors:?}");
        assert!(
            errors[1] < 0.36 * errors[0] && errors[2] < 0.36 * errors[1],
            "{end:?}, {errors:?}"
        );
        assert!(errors[2] < 0.004, "{end:?}, {errors:?}");
    }
}

fn dynamic_load(radius: f64, medium: GasState) -> ModalAdmittanceLoad {
    // Authored acoustic displacement/momentum oscillator with original
    // resistance and two work-conjugate, mutually coupled pressure inputs.
    let original = PortHamiltonian::new(
        2,
        2,
        vec![0., 1., -1., 0.],
        vec![0., 0., 0., 12_000.],
        vec![0., 0., 1., 0.6],
        Box::new(QuadraticStorage::new(vec![1.8e8, 0., 0., 1. / 37.], 2).unwrap()),
    )
    .unwrap();
    let z = medium.density * medium.sound_speed / (PI * radius * radius);
    ModalAdmittanceLoad::new(original, radius, ModeCount(2), medium)
        .unwrap()
        .with_direct_conductance(vec![0.15 / z, 0.08 / z, 0.08 / z, 0.2 / z])
        .unwrap()
}

#[test]
fn g1_g3_regional_dynamic_load_keeps_all_states_and_original_heat_independent_of_clock_and_history()
{
    let (cold, warm) = (gas(293.15), gas(600.));
    let gases = [cold, cold, warm, warm];
    let radius = 0.025;
    let profile = RadiusProfile::exponential(radius, 0., 0.04).unwrap();
    let omega = TAU * 650.;
    let wave = [C64::from_re(50.), C64::ZERO];
    let zm = warm.density * warm.sound_speed / (PI * radius * radius);
    let conductance = [0.15 / zm, 0.08 / zm, 0.08 / zm, 0.2 / zm];
    let mut reference = None;
    for dt in [1. / 24_000., 1. / 48_000., 1. / 96_000.] {
        let mut line = ContinuousModalLine::new_regional_loaded(
            &profile,
            &gases,
            ModeCount(2),
            dynamic_load(radius, warm),
            dt,
            time_budget(),
            || false,
        )
        .unwrap();
        assert_eq!(line.source().regional_gases(), Some(gases.as_slice()));
        let n = line.source().coordinates();
        assert_eq!(n, 10);
        assert_eq!(line.state().len(), 22);
        let before = snapshot(&line);
        let response = query(&line, omega, &wave);
        assert_eq!(snapshot(&line), before);
        assert_eq!(response.physical.state.len(), 22);
        assert_eq!(response.waves.len(), 2);
        assert!(
            response.physical.state[2 * n..]
                .iter()
                .all(|z| z.abs() > 1e-12)
        );
        assert!(
            response.waves[1].outgoing.abs() > 1e-5,
            "original mutual mouth load must scatter modes"
        );
        // Pressure is the original momentum-block co-state at each mouth
        // node. Evaluate the actual two loss mechanisms from their own laws.
        let mouth = &response.physical.effort[n + (gases.len() * 2)..n + (gases.len() * 2) + 2];
        let momentum = response.physical.state[2 * n + 1];
        let oscillator_heat = 0.5 * 12_000. * (momentum.scale(1. / 37.)).norm_sq();
        let mut direct_heat = 0.;
        for i in 0..2 {
            for j in 0..2 {
                direct_heat += 0.5 * (mouth[i].conj() * mouth[j]).re * conductance[2 * i + j];
            }
        }
        assert!(oscillator_heat > 1e-10 && direct_heat > 1e-10);
        let heat = oscillator_heat + direct_heat;
        close(response.physical.dissipated_power_w, heat, heat, 2e-11);
        let incident = response
            .waves
            .iter()
            .map(|p| p.incident_power_w)
            .sum::<f64>();
        close(response.wave_power_residual_w, 0., incident, 2e-9);
        close(response.physical.power_residual_w, 0., incident, 2e-9);
        for step in 0..12 {
            line.push(&[35., -12. + step as f64]).unwrap();
        }
        assert_ne!(snapshot(&line), before);
        let advanced = snapshot(&line);
        assert_eq!(
            query(&line, omega, &wave),
            response,
            "history must not enter the continuous frequency query"
        );
        assert_eq!(snapshot(&line), advanced);
        if let Some(prior) = &reference {
            assert_eq!(
                &response, prior,
                "runtime dt must not warp harmonic scattering"
            );
        }
        reference = Some(response);
    }
}

#[test]
fn g1_static_mouth_loss_keeps_ordinary_source_slots_and_adds_no_fake_memory() {
    let medium = gas(293.15);
    let radius = 0.025;
    let profile = RadiusProfile::exponential(radius, 0., 0.04).unwrap();
    let z = medium.density * medium.sound_speed / (PI * radius * radius);
    let load = ModalAdmittanceLoad::from_conductance(
        vec![1. / z, 0., 0., 1. / z],
        radius,
        ModeCount(2),
        medium,
    )
    .unwrap();
    let line = ContinuousModalLine::new_loaded(
        &profile,
        &medium,
        ModeCount(2),
        12,
        load,
        1. / 96_000.,
        time_budget(),
        || false,
    )
    .unwrap();
    let model = line
        .prepare_harmonic_scattering(harmonic_budget(), scattering_budget(), || false)
        .unwrap();
    assert_eq!(line.load_state().len(), 0);
    assert_eq!(model.dimension(), 52);
    assert_eq!(model.port_count(), 4);
    assert_eq!(model.wave_ports(), &[(0, z), (1, z)]);
    let factor = model.factor(TAU * 600., || false).unwrap();
    let result = factor
        .response(&[C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO], || false)
        .unwrap();
    assert_eq!(&result.physical_inputs[2..], &[C64::ZERO, C64::ZERO]);
    assert_eq!(result.physical.state.len(), 52);
    assert!(result.waves[0].outgoing.abs() < 0.001);
    let heat =
        (result.physical.output[2].norm_sq() + result.physical.output[3].norm_sq()) / (2. * z);
    close(result.physical.dissipated_power_w, heat, heat, 2e-11);
    close(result.wave_power_residual_w, 0., 1. / (2. * z), 2e-9);
    assert_eq!(snapshot(&line), vec![0_u64; snapshot(&line).len()]);
}
