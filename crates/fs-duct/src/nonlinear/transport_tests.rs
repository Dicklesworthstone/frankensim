use super::super::transport::{DissipatingDuct, DissipatingDuctError, DissipatingStepControl};
use super::*;
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};

fn species(name: &str, r: f64, linear_cp: f64, formation: f64) -> Nasa9StandardStateModelV1 {
    let base = gas(linear_cp, r);
    let source = &base.source().components()[0];
    let mut card = source.regions()[0].card().clone();
    let formation_name = required_nasa9_parameter_dims()
        .into_iter()
        .nth(7)
        .unwrap()
        .0;
    card.parameters.get_mut(formation_name).unwrap().value = formation;
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        source.molar_mass(),
        source.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn source(models: Vec<Nasa9StandardStateModelV1>, fractions: Vec<f64>) -> FrozenGas {
    FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            models,
            Composition::new(CompositionBasis::MoleFraction, fractions).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(5_000.0)),
        || false,
    )
    .unwrap()
}
fn law(gas: &FrozenGas, viscosity: f64, conductivity: f64) -> ConstantTransport {
    let models = gas.source().components().to_vec();
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap();
    let mut diffusion = vec![0.01; ids.len() * ids.len()];
    for i in 0..ids.len() {
        diffusion[i * ids.len() + i] = 0.0;
    }
    ConstantTransport::new(
        bank,
        "manufactured-longitudinal-duct-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 5_000.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(viscosity),
            bulk_viscosity: DynViscosity::new(0.0),
            thermal_conductivity: ThermalConductivityQuantityV1::new(conductivity),
        },
        diffusion,
    )
    .unwrap()
}
fn acoustic_gas() -> FrozenGas {
    source(vec![species("acoustic", 1.0, 0.0, 0.0)], vec![1.0])
}
fn prepare(gas: FrozenGas, n: usize, boundaries: Boundaries) -> NonlinearDuct {
    NonlinearDuct::new(&cylinder(0.1, 1.0), gas, &[n], boundaries, |x| {
        Ok(InitialGasState {
            density_kg_m3: 1.0,
            velocity_m_s: 0.0,
            pressure_pa: if x < 0.5 { 1.8 } else { 1.0 },
        })
    })
    .unwrap()
}

#[test]
fn g0_molecular_attachment_preserves_original_geometry_species_and_physical_energy() {
    let g = source(
        vec![species("A", 1.0, 0.0, 1e7), species("B", 0.5, 0.2, -3e7)],
        vec![0.25, 0.75],
    );
    let initial = NonlinearDuct::new(&cylinder(0.1, 1.0), g.clone(), &[12], closed(), |x| {
        Ok(InitialGasState {
            density_kg_m3: 1.0 + 0.2 * x,
            velocity_m_s: 0.37,
            pressure_pa: 1.8,
        })
    })
    .unwrap();
    let initial_fields = initial.samples().unwrap();
    let initial_momentum: Vec<_> = initial
        .solver()
        .states()
        .iter()
        .map(|q| q.momentum)
        .collect();
    let initial_energy: Vec<_> = initial.solver().states().iter().map(|q| q.energy).collect();
    let physical_energy: Vec<_> = initial
        .solver()
        .states()
        .iter()
        .map(|&q| g.absolute_energy_density(q).unwrap())
        .collect();
    let geometry = initial.geometry().clone();
    let mesh = initial.solver().mesh().clone();
    let flow = initial
        .with_molecular_transport(law(&g, 0.03, 0.02), || false)
        .unwrap();
    assert_eq!(flow.geometry(), &geometry);
    assert_eq!(flow.initial_gas(), &g);
    assert_eq!(flow.solver().branches()[0].mesh(), &mesh);
    assert_eq!(flow.time().value(), 0.0);
    assert_eq!(
        flow.solver().gas().source().models(),
        g.source().components()
    );
    for (i, expected) in initial_fields.iter().enumerate() {
        let q = flow.solver().branches()[0].cell(i).unwrap();
        relative(q.partial_densities[0], expected.density_kg_m3 / 7.0, 1e-14);
        relative(
            q.partial_densities[1],
            6.0 * expected.density_kg_m3 / 7.0,
            1e-14,
        );
        assert_eq!(q.momentum, initial_momentum[i]);
        assert_eq!(q.energy, initial_energy[i]);
        relative(
            flow.solver().gas().absolute_energy_density(q).unwrap(),
            physical_energy[i],
            3e-15,
        );
        let actual = flow.sample(i).unwrap();
        relative(actual.pressure_pa, expected.pressure_pa, 3e-12);
        relative(actual.temperature_k, expected.temperature_k, 3e-12);
    }
}

#[test]
fn g3_zero_molecular_coefficients_retain_the_inviscid_shock_at_the_same_clock() {
    let g = acoustic_gas();
    let mut inviscid = prepare(g.clone(), 32, closed());
    let mut coupled = prepare(g.clone(), 32, closed())
        .with_molecular_transport(law(&g, 0.0, 0.0), || false)
        .unwrap();
    for _ in 0..20 {
        let step = coupled
            .step(Time::new(0.01), DissipatingStepControl::default(), || false)
            .unwrap();
        let reference = inviscid
            .step(step.elapsed, StepControl::default(), || false)
            .unwrap();
        assert_eq!(reference.elapsed, step.elapsed);
        assert_eq!(inviscid.time(), coupled.time());
        for (i, expected) in inviscid.solver().states().iter().enumerate() {
            let actual = coupled.solver().branches()[0].cell(i).unwrap();
            relative(actual.partial_densities[0], expected.density, 2e-11);
            assert!((actual.momentum - expected.momentum).abs() < 2e-11);
            relative(actual.energy, expected.energy, 2e-11);
        }
    }
}

fn harmonic_error(n: usize) -> f64 {
    // Linearized compressible NSF in a periodic cylinder, with rho0=c0=1.
    // Pressure/velocity Fourier coefficients solve the independent exact 2x2
    // system p_t=-ik u, u_t=-ik p-nu_l k²u. The finite-amplitude solver is
    // tested in its small-amplitude limit; no solver output defines the oracle.
    let epsilon = 1e-4;
    let mu = 0.02;
    let g = acoustic_gas();
    let k = 2.0 * core::f64::consts::PI;
    let mut flow = NonlinearDuct::new(
        &cylinder(0.1, 1.0),
        g.clone(),
        &[n],
        Boundaries::Periodic,
        |x| {
            let wave = epsilon * det::sin(k * x);
            Ok(InitialGasState {
                density_kg_m3: 1.0 + wave,
                velocity_m_s: wave,
                pressure_pa: 1.0 / 1.4 + wave,
            })
        },
    )
    .unwrap()
    .with_molecular_transport(law(&g, mu, 0.0), || false)
    .unwrap();
    let end = 0.2;
    let mut steps = 0;
    while flow.time().value() < end {
        let report = flow
            .step(
                Time::new(end - flow.time().value()),
                DissipatingStepControl::default(),
                || false,
            )
            .unwrap();
        assert!(report.mass_residual_kg.abs() < 1e-15);
        assert!(report.absolute_energy_residual_j.abs() < 2e-14);
        assert_eq!(
            report.molecular.first.elapsed.value() * 2.0,
            report.elapsed.value()
        );
        assert_eq!(
            report.molecular.second.elapsed.value() * 2.0,
            report.elapsed.value()
        );
        steps += 1;
        assert!(steps < 10_000);
    }
    let a = (4.0 / 3.0) * mu * k * k / 2.0;
    let omega = det::sqrt(k * k - a * a);
    let factor = epsilon * det::exp(-a * end);
    let s = det::sin(omega * end) / omega;
    let c = det::cos(omega * end);
    let expected = [
        -factor * k * s,
        factor * (c + a * s),
        -factor * k * s,
        factor * (c - a * s),
    ];
    let mut measured = [0.0; 4];
    for q in flow.samples().unwrap() {
        let c = det::cos(k * q.axial_position_m) * 2.0 / n as f64;
        let s = det::sin(k * q.axial_position_m) * 2.0 / n as f64;
        measured[0] += (q.pressure_pa - 1.0 / 1.4) * c;
        measured[1] += (q.pressure_pa - 1.0 / 1.4) * s;
        measured[2] += q.velocity_m_s * c;
        measured[3] += q.velocity_m_s * s;
    }
    det::sqrt(
        measured
            .iter()
            .zip(expected)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>(),
    ) / epsilon
}

#[test]
fn g1_longitudinal_sound_attenuation_and_phase_converge_to_the_nsf_mode() {
    let errors = [16, 32, 64].map(harmonic_error);
    eprintln!("nonlinear duct NSF Fourier-mode relative errors: {errors:?}");
    assert!(errors[0] < 0.08);
    assert!(errors[1] < 0.34 * errors[0]);
    assert!(errors[2] < 0.34 * errors[1]);
    assert!(errors[2] < 0.005);
}

#[test]
fn g0_finite_shock_molecular_transport_heats_without_losing_total_energy() {
    let g = acoustic_gas();
    let mut flow = prepare(g.clone(), 48, closed())
        .with_molecular_transport(law(&g, 0.02, 0.03), || false)
        .unwrap();
    let before = flow.solver().totals();
    for _ in 0..30 {
        let report = flow
            .step(Time::new(0.01), DissipatingStepControl::default(), || false)
            .unwrap();
        assert!(report.molecular.first.accepted_steps > 0);
        assert!(report.molecular.second.accepted_steps > 0);
        assert!(report.absolute_energy_residual_j.abs() < 3e-14);
        assert!(report.transport.maximum_cfl <= 0.40000000000001);
    }
    let after = flow.solver().totals();
    relative(before.mass(), after.mass(), 1e-13);
    relative(before.energy, after.energy, 1e-13);
    let wave = flow.sample(30).unwrap();
    assert!(
        wave.temperature_k > 1.0 && wave.pressure_pa > 1.0 && wave.velocity_m_s > 0.0,
        "{wave:?}"
    );
}

#[test]
fn g0_molecular_refusals_and_cancellation_preserve_accepted_wave_and_report() {
    let g = acoustic_gas();
    let mut evolved = prepare(g.clone(), 12, closed());
    evolved
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    assert!(matches!(
        evolved.with_molecular_transport(law(&g, 0.01, 0.01), || false),
        Err(DissipatingDuctError::Invalid(_))
    ));
    assert!(matches!(
        prepare(g.clone(), 12, closed())
            .with_molecular_transport(law(&gas(0.0, 1.0), 0.01, 0.01), || false),
        Err(DissipatingDuctError::Invalid(_))
    ));
    let cone = Duct {
        segments: vec![Segment::Cone {
            length: 1.0,
            inlet_radius: 0.1,
            outlet_radius: 0.2,
        }],
    };
    let tapered = NonlinearDuct::new(&cone, g.clone(), &[12], closed(), |_| {
        Ok(InitialGasState {
            density_kg_m3: 1.0,
            velocity_m_s: 0.0,
            pressure_pa: 1.0,
        })
    })
    .unwrap();
    assert!(
        tapered
            .with_molecular_transport(law(&g, 0.01, 0.01), || false)
            .is_err()
    );
    let mut flow: DissipatingDuct<_> = prepare(g.clone(), 12, closed())
        .with_molecular_transport(law(&g, 0.01, 0.01), || false)
        .unwrap();
    flow.step(Time::new(0.001), DissipatingStepControl::default(), || {
        false
    })
    .unwrap();
    let before = flow.samples().unwrap();
    let time = flow.time();
    let report = flow.solver().last_report();
    let mut calls = 0;
    assert!(
        flow.step(Time::new(0.01), DissipatingStepControl::default(), || {
            calls += 1;
            calls > 80
        })
        .is_err()
    );
    assert_eq!(flow.samples().unwrap(), before);
    assert_eq!(flow.time(), time);
    assert_eq!(flow.solver().last_report(), report);
}
