use super::*;
use fs_material::{
    gas::GasSpec,
    relaxation::{
        InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
        RelaxingMedium, SpeciesRelaxation,
    },
};
use fs_thermochem::SpeciesId;

fn gas(pressure: f64) -> GasState {
    GasState::try_new(
        &GasSpec {
            gamma: 1.25,
            ..GasSpec::dry_air_ussa1976()
        },
        293.15,
        pressure,
    )
    .unwrap()
}
fn relaxation(gas: GasState, relaxing: bool) -> RelaxationState {
    let id = SpeciesId::new("fixture").unwrap();
    let processes = if relaxing {
        vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 202.65,
                t_ref_k: 293.15,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 400.0,
            },
        }]
    } else {
        vec![]
    };
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes,
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
    .unwrap()
}
fn source(gas: GasState, loss: CavityNeckLoss) -> HelmholtzModalCavity {
    HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: 1e-7,
            neck_radius_m: 1e-4,
            neck_length_m: 0.1,
            inner_end: fs_scenario::NeckEnd::Flanged,
            outer_end: fs_scenario::NeckEnd::Flanged,
        },
        gas,
        loss,
    )
    .unwrap()
}
fn budget() -> ModalCouplingConfig {
    ModalCouplingConfig {
        max_modes: 8,
        max_connections: 6,
        max_setup_terms: 20000,
        nyquist_guard_fraction: 0.9,
        maximum_total_energy_j: 1.0,
        maximum_abs_pressure_pa: 1e5,
        maximum_abs_connection_force_n: 100.0,
        solve_relative_tolerance: 1e-10,
        energy_absolute_tolerance_j: 1e-14,
        energy_relative_tolerance: 1e-8,
    }
}
fn surface(rate: u32) -> CavitySurface {
    CavitySurface {
        model: ModalAcousticTimeModel::try_new(
            rate,
            vec![ModalAcousticMode {
                angular_frequency_rad_s: 500.0,
                damping_ratio: 0.0,
                pressure_per_modal_velocity: C64::ZERO,
            }],
            ModalAcousticTimeBudget::audible_reference(),
        )
        .unwrap(),
        outward_volume_shapes: vec![1e-4],
    }
}
fn create(rate: u32, viscous: bool, history: Option<Vec<f64>>) -> SharedCavity {
    let gas = gas(101325.0);
    SharedCavity::new_with_molecular_relaxation(
        source(
            gas,
            if viscous {
                CavityNeckLoss::ExteriorMouth
            } else {
                CavityNeckLoss::Resistance(0.0)
            },
        ),
        vec![surface(rate)],
        rate,
        ModalAcousticState {
            displacement_m_sqrt_kg: 2e-7,
            velocity_m_sqrt_kg_per_s: 0.0,
        },
        ModalAcousticTimeBudget::audible_reference(),
        budget(),
        viscous.then_some(ViscousNeckOptions {
            poles: 4,
            band_rad_s: [10.0, 1000.0],
            maximum_relative_error: 1e-3,
        }),
        relaxation(gas, true),
        history,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + 2e-10 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

#[test]
fn molecular_cavity_matches_independent_internal_temperature_dynamics() {
    // G1: integrate physical gas temperatures directly, independently of the
    // Maxwell transform, while the same pressure drives neck and elastic wall.
    let gate = CancelGate::new_clock_free();
    let relaxed = create(48000, false, None);
    let mut actual = create(48000, false, Some(vec![0.0]));
    let source = actual.source().clone();
    let gas = source.gas();
    let law = actual
        .cavity_molecular_compression()
        .unwrap()
        .bulk()
        .source()
        .clone();
    let cv = law.frozen_cv_j_kg_k();
    let capacity = law.processes()[0].heat_capacity_j_kg_k;
    let tau = law.processes()[0].relaxation_time_s;
    let volume = source.geometry().volume_m3;
    let area = source.neck_area_m2();
    let inertia = source.neck_mass_kg() / (area * area);
    let ratio = actual
        .cavity_molecular_compression()
        .unwrap()
        .frozen_stiffness_n_m()
        / actual
            .cavity_molecular_compression()
            .unwrap()
            .equilibrium_stiffness_n_m();
    close(
        actual.observe().unwrap().pressure_increment_pa,
        relaxed.observe().unwrap().pressure_increment_pa * ratio,
        1e-10,
    );
    let mut x = vec![
        actual.observe().unwrap().volume_change_m3,
        0.0,
        0.0,
        0.0,
        0.0,
    ];
    let temperature = |x: &[f64]| {
        (-gas.specific_gas_constant * gas.temperature * x[0] / volume - capacity * x[4]) / cv
    };
    let pressure = |x: &[f64]| gas.pressure * (temperature(x) / gas.temperature - x[0] / volume);
    let derivative = |x: &[f64], external: f64| {
        let p = pressure(x);
        vec![
            x[1] + 1e-4 * x[3],
            (p - external) / inertia,
            x[3],
            -500.0_f64.powi(2) * x[2] + 1e-4 * p,
            (temperature(x) - x[4]) / tau,
        ]
    };
    let initial = actual.network().total_energy_j().unwrap();
    let (mut work, mut loss) = (0.0, 0.0);
    for tick in 0..128 {
        let external = if tick < 16 { 0.1 } else { 0.0 };
        let dt = 1.0 / (48000.0 * 16.0);
        for _ in 0..16 {
            let stage = |k: &[f64], a: f64| {
                x.iter()
                    .zip(k)
                    .map(|(x, k)| x + a * dt * k)
                    .collect::<Vec<_>>()
            };
            let a = derivative(&x, external);
            let b = derivative(&stage(&a, 0.5), external);
            let c = derivative(&stage(&b, 0.5), external);
            let d = derivative(&stage(&c, 1.0), external);
            for i in 0..x.len() {
                x[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
            }
        }
        let frame = actual.step(&[0.0], external, &gate).unwrap().clone();
        let observed = actual.observe().unwrap();
        close(observed.pressure_increment_pa, pressure(&x), 1e-8);
        close(observed.neck_outflow_m3_s, x[1], 1e-18);
        close(observed.volume_change_m3, x[0], 1e-21);
        let theta = temperature(&x);
        let strain = x[0] / volume;
        let availability = 0.5
            * volume
            * (gas.pressure * strain * strain
                + gas.density / gas.temperature * (cv * theta * theta + capacity * x[4] * x[4]));
        close(observed.compression_energy_j, availability, 1e-19);
        let wall_energy = 0.5 * (x[3] * x[3] + 500.0_f64.powi(2) * x[2] * x[2]);
        close(
            frame.total_energy_j(),
            availability + 0.5 * inertia * x[1] * x[1] + wall_energy,
            1e-19,
        );
        assert!(frame.material_dissipation_j >= 0.0);
        work += frame.external_work_j;
        loss += frame.total_dissipation_j();
    }
    close(
        actual.network().total_energy_j().unwrap() + loss,
        initial + work,
        1e-19,
    );
    assert!(loss > 0.001 * initial);
}

#[test]
fn molecular_cavity_harmonic_pressure_and_loss_share_the_material_law() {
    let gate = CancelGate::new_clock_free();
    let cavity = create(48000, false, None);
    let volume = cavity.source().geometry().volume_m3;
    let source = cavity
        .cavity_molecular_compression()
        .unwrap()
        .bulk()
        .source();
    for omega in [100.0, 400.0, 1000.0, 5000.0] {
        let response = cavity
            .harmonic_response(
                omega,
                &[C64::from_re(0.01)],
                C64::from_re(0.03),
                20000,
                &gate,
            )
            .unwrap();
        let bulk = source.caloric_response(omega).unwrap().bulk_modulus_pa;
        let pressure = (bulk * response.volume_change_m3).scale(-1.0 / volume);
        close((response.pressure_increment_pa - pressure).abs(), 0.0, 1e-8);
        let loss = -0.5 * omega * bulk.im * response.fractional_volume_change.norm_sq() * volume;
        close(response.network.material_dissipation_w, loss, 1e-16);
        assert!(loss > 0.0);
    }
    assert_eq!(cavity.network().samples_rendered(), 0);
}

#[test]
fn molecular_compression_composes_and_rejects_without_publishing_state() {
    use crate::render::plate::impact::radiation::Pole;
    let gate = CancelGate::new_clock_free();
    let radiation = RadiationModel {
        ports: 1,
        poles: vec![Pole {
            omega: 700.0,
            zeta: 0.1,
            coupling: vec![2e6],
        }],
    };
    let mut cavity = create(48000, true, Some(vec![0.0]))
        .with_exterior_radiation(&radiation, &[1], 1, &gate)
        .unwrap();
    for _ in 0..12 {
        cavity.step(&[0.0], 0.03, &gate).unwrap();
    }
    let observation = cavity.observe().unwrap();
    let states: Vec<_> = cavity
        .network()
        .components()
        .iter()
        .map(|m| m.states().to_vec())
        .collect();
    let memory = cavity.network().material_memory().to_vec();
    let exterior = cavity.network().radiation_state().to_vec();
    let last = cavity.network().last_frame().unwrap().clone();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        cavity.step(&[0.0], 0.0, &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    assert!(cavity.step(&[0.0], 1e30, &gate).is_err());
    assert_eq!(cavity.observe().unwrap(), observation);
    assert_eq!(cavity.network().material_memory(), memory);
    assert_eq!(cavity.network().radiation_state(), exterior);
    assert_eq!(cavity.network().last_frame().unwrap(), &last);
    for (model, expected) in cavity.network().components().iter().zip(states) {
        assert_eq!(model.states(), expected);
    }
    let create = |relaxation, history, gate: &CancelGate| {
        SharedCavity::new_with_molecular_relaxation(
            source(gas(101325.0), CavityNeckLoss::Resistance(0.0)),
            vec![surface(48000)],
            48000,
            ModalAcousticState::default(),
            ModalAcousticTimeBudget::audible_reference(),
            budget(),
            None,
            relaxation,
            history,
            gate,
        )
    };
    assert!(create(relaxation(gas(90000.0), true), None, &gate).is_err());
    assert!(create(relaxation(gas(101325.0), true), Some(vec![]), &gate).is_err());
    assert!(matches!(
        create(relaxation(gas(101325.0), true), None, &cancelled),
        Err(ModalCouplingError::Cancelled)
    ));
    let nonrelaxing = create(relaxation(gas(101325.0), false), None, &gate).unwrap();
    assert!(nonrelaxing.network().material_memory().is_empty());
}
