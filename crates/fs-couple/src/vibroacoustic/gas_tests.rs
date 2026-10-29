//! G1 constitutive/coupled analytic solutions and G3 modal-basis invariance.
//! Synthetic rotational rates test implementation, not measured gas accuracy.
use super::*;
use fs_material::gas::ConductivityModel;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

fn relaxation(pressure: f64, reservoirs: bool) -> RelaxationState {
    let gas = GasState::try_new(
        &GasSpec {
            molar_mass: 0.0289644,
            gamma: 1.25,
            sutherland_beta: 1.458e-6,
            sutherland_s: 110.4,
            conductivity: ConductivityModel::Eucken,
        },
        293.15,
        pressure,
    )
    .unwrap();
    let species = SpeciesId::new("A").unwrap();
    let law = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: if reservoirs {
            vec![RelaxationProcess {
                reservoir: InternalReservoir::Rotational {
                    degrees_of_freedom: 2,
                },
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * 0.025,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 400.0,
                },
            }]
        } else {
            vec![]
        },
    }])
    .unwrap();
    law.evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
        .unwrap()
}

fn model(
    law: &RelaxationState,
    ws: &[f64],
    wa: &[f64],
    lambda: &[f64],
    coupling: &[f64],
) -> VibroacousticModel {
    VibroacousticModel::try_new(
        &StructuralModes {
            omegas: ws.to_vec(),
            shapes: vec![vec![1.0]; ws.len()],
            loss_factor: 0.015,
        },
        &CavityModes {
            omegas: wa.to_vec(),
            lambdas: lambda.to_vec(),
            interface: vec![vec![1.0]; wa.len()],
            loss_factor: 0.0,
            rho0: law.gas().density,
            c0: law.gas().sound_speed,
        },
        coupling.to_vec(),
        None,
    )
    .unwrap()
}

fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-25),
        "{a:e} != {b:e}"
    );
}
fn complex_close(a: C64, b: C64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * a.abs().max(b.abs()).max(1e-25),
        "{a:?} != {b:?}"
    );
}

// Exact one-reservoir K, calculated independently using a real rational form.
fn single_bulk(law: &RelaxationState, omega: f64) -> C64 {
    let gas = law.gas();
    let cv0 = gas.specific_heat_cp - gas.specific_gas_constant;
    let cvi = cv0 - gas.specific_gas_constant;
    let cp0 = gas.specific_heat_cp;
    let cpi = cvi + gas.specific_gas_constant;
    let y = omega * (101325.0 * 0.025 / gas.pressure);
    let denominator = cv0 * cv0 + y * y * cvi * cvi;
    C64::new(
        gas.pressure * (cp0 * cv0 + y * y * cpi * cvi) / denominator,
        gas.pressure * y * (cp0 * cvi - cpi * cv0) / denominator,
    )
}

#[test]
fn sealed_uniform_mode_has_the_gas_stiffness_and_real_molecular_heat() {
    let law = relaxation(101325.0, true);
    let (area, volume, ws) = (0.003, 0.002, 30.0);
    let coupled = model(&law, &[ws], &[0.0], &[volume], &[area])
        .with_cavity_relaxation(law.clone())
        .unwrap();
    let force = C64::new(1.0, 0.2);
    for omega in [0.1, 10.0, 30.0, 40.0, 80.0, 2000.0] {
        let k = single_bulk(&law, omega);
        let stiffness =
            C64::new(ws * ws - omega * omega, -0.015 * ws * ws) + k.scale(area * area / volume);
        let displacement = force / stiffness;
        let pressure = -k.scale(area / volume) * displacement;
        let response = coupled.frf(omega, &[force]).unwrap();
        complex_close(response.b[0], displacement, 5e-13);
        complex_close(response.a[0], pressure, 5e-13);
        let strain = displacement.scale(area / volume);
        let heat = -0.5 * omega * k.im * volume * strain.norm_sq();
        close(response.power.cavity, heat, 5e-13);
        close(response.power.interface, heat, 5e-12);
        assert!(heat > 0.0, "uniform compression must relax too");
        // Far above resonance, input work subtracts nearly reactive terms.
        // Bound that rounding by apparent force/velocity power, while the
        // independent displacement, pressure and heat checks stay relative.
        let roundoff = 64.0 * f64::EPSILON * force.abs() * omega * displacement.abs();
        let imbalance = response.power.input
            - response.power.structural
            - response.power.interface
            - response.power.radiated;
        assert!(imbalance.abs() <= roundoff.max(5e-13 * response.power.input.abs()));
        assert!(response.power.interface_residual.abs() < 5e-12);
        assert_eq!(response.a, coupled.frf(omega, &[force]).unwrap().a);
    }
}

#[test]
fn nonuniform_modes_match_independent_schur_solution_and_integrated_cycle_work() {
    let law = relaxation(101325.0, true);
    let ws = [150.0, 230.0];
    let wa = [0.0, 190.0];
    let lambda = [0.003, 0.0015];
    let c = [0.003, 0.0008, -0.001, 0.002];
    let forces = [C64::new(1.0, 0.3), C64::new(-0.4, 0.2)];
    let coupled = model(&law, &ws, &wa, &lambda, &c)
        .with_cavity_relaxation(law.clone())
        .unwrap();
    for omega in [20.0, 145.0, 187.0, 240.0, 370.0] {
        let beta = single_bulk(&law, omega).recip();
        let gas = law.gas();
        let h = [0, 1].map(|q| {
            let spatial = wa[q] * wa[q] / (gas.sound_speed * gas.sound_speed);
            C64::from_re(gas.density * omega * omega / lambda[q])
                / (C64::from_re(spatial) - beta.scale(gas.density * omega * omega))
        });
        let mut d = [C64::ZERO; 4];
        for r in 0..2 {
            for s in 0..2 {
                if r == s {
                    d[r * 2 + s] = C64::new(ws[r] * ws[r] - omega * omega, -0.015 * ws[r] * ws[r]);
                }
                for q in 0..2 {
                    d[r * 2 + s] = d[r * 2 + s] - h[q].scale(c[r * 2 + q] * c[s * 2 + q]);
                }
            }
        }
        let determinant = d[0] * d[3] - d[1] * d[2];
        let b = [
            (d[3] * forces[0] - d[1] * forces[1]) / determinant,
            (d[0] * forces[1] - d[2] * forces[0]) / determinant,
        ];
        let a = [0, 1].map(|q| h[q] * (b[0].scale(c[q]) + b[1].scale(c[2 + q])));
        let result = coupled.frf(omega, &forces).unwrap();
        for i in 0..2 {
            complex_close(result.b[i], b[i], 8e-13);
            complex_close(result.a[i], a[i], 8e-13);
        }
        // Integrate -p(t) * volumetric-strain-rate(t) over one complete
        // cycle, independently of both interface force and |p|^2 formulas.
        let mut work = 0.0;
        for sample in 0..128 {
            let phase = core::f64::consts::TAU * sample as f64 / 128.0;
            let (sine, cosine) = phase.sin_cos();
            for q in 0..2 {
                let strain_rate = C64::new(0.0, omega) * beta * a[q];
                let p = a[q].re * cosine + a[q].im * sine;
                let rate = strain_rate.re * cosine + strain_rate.im * sine;
                work -= lambda[q] * p * rate / 128.0;
            }
        }
        close(result.power.cavity, work, 8e-13);
        assert!(result.power.cavity > 0.0);
        assert!(result.power.balance_residual.abs() < 8e-13);
        assert!(result.power.interface_residual.abs() < 8e-12);
        let drive0 = coupled.frf(omega, &[C64::ONE, C64::ZERO]).unwrap();
        let drive1 = coupled.frf(omega, &[C64::ZERO, C64::ONE]).unwrap();
        complex_close(drive0.b[1], drive1.b[0], 8e-13);
    }
}

#[test]
fn modal_rescaling_and_truncation_keep_the_same_gas_work() {
    let law = relaxation(101325.0, true);
    let ws = [150.0, 230.0];
    let wa = [0.0, 190.0];
    let lambda = [0.003, 0.0015];
    let coupling = [0.003, 0.0008, -0.001, 0.002];
    let full = model(&law, &ws, &wa, &lambda, &coupling)
        .with_cavity_relaxation(law.clone())
        .unwrap();
    let scales = [-3.0, 0.25];
    let scaled = model(
        &law,
        &ws,
        &wa,
        &[
            lambda[0] * scales[0] * scales[0],
            lambda[1] * scales[1] * scales[1],
        ],
        &coupling
            .iter()
            .enumerate()
            .map(|(i, c)| c * scales[i % 2])
            .collect::<Vec<_>>(),
    )
    .with_cavity_relaxation(law.clone())
    .unwrap();
    let forces = [C64::new(1.0, 0.3), C64::new(-0.4, 0.2)];
    let first = full.frf(187.0, &forces).unwrap();
    let second = scaled.frf(187.0, &forces).unwrap();
    for i in 0..2 {
        complex_close(first.b[i], second.b[i], 8e-13);
        complex_close(first.a[i], second.a[i].scale(scales[i]), 8e-13);
    }
    close(first.power.cavity, second.power.cavity, 8e-13);
    close(first.power.input, second.power.input, 8e-13);
    let reduced = model(&law, &ws[..1], &wa[..1], &lambda[..1], &coupling[..1])
        .with_cavity_relaxation(law)
        .unwrap();
    let truncated = full.frf_truncated(187.0, &forces, 1, 1).unwrap();
    let expected = reduced.frf(187.0, &forces[..1]).unwrap();
    complex_close(truncated.b[0], expected.b[0], 8e-13);
    close(truncated.power.cavity, expected.power.cavity, 8e-13);
}

#[test]
fn material_binding_rejects_stale_gas_and_duplicate_loss_and_recovers_lossless_limit() {
    let law = relaxation(101325.0, true);
    let bare = model(&law, &[30.0], &[0.0], &[0.002], &[0.003]);
    assert!(
        bare.clone()
            .with_cavity_relaxation(relaxation(202650.0, true))
            .is_err()
    );
    let mut lossy = bare.clone();
    lossy.eta_a = 0.01;
    assert!(lossy.with_cavity_relaxation(law.clone()).is_err());
    let mut same_bulk_different_medium = bare.clone();
    same_bulk_different_medium.medium_a.rho0 *= 4.0;
    same_bulk_different_medium.medium_a.c0 *= 0.5;
    assert!(
        same_bulk_different_medium
            .with_cavity_relaxation(law.clone())
            .is_err()
    );
    let mut incoherent = *law.gas();
    incoherent.sound_speed *= 2.0;
    assert!(AcousticMedium::try_from_gas(&incoherent).is_err());
    let lossless = bare
        .clone()
        .with_cavity_relaxation(relaxation(101325.0, false))
        .unwrap();
    let reference = bare.frf(40.0, &[C64::ONE]).unwrap();
    let explicit = lossless.frf(40.0, &[C64::ONE]).unwrap();
    complex_close(reference.b[0], explicit.b[0], 2e-14);
    complex_close(reference.a[0], explicit.a[0], 2e-14);
    assert_eq!(explicit.power.cavity, 0.0);
    assert_eq!(
        bare.undamped_natural_frequencies().unwrap(),
        lossless.undamped_natural_frequencies().unwrap()
    );
    for omega in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::MAX,
        f64::MIN_POSITIVE,
    ] {
        assert!(lossless.frf(omega, &[C64::ONE]).is_err());
    }
}

#[test]
fn explicit_humid_air_species_change_coupled_cavity_response() {
    use fs_material::gas_species::input::read_gas_record;
    use fs_material::relaxation::input::medium_from_record;
    let law = RelaxationAbsorption::iso9613_terrestrial_air().unwrap();
    let solve = |record: &str| {
        let gas = read_gas_record(&mut record.split_ascii_whitespace(), None, None).unwrap();
        let relaxation = law.evaluate(&medium_from_record(&gas).unwrap()).unwrap();
        let cavity = rectangular_cavity_modes(
            0.4,
            0.2,
            0.15,
            AcousticMedium::try_from_gas(gas.state()).unwrap(),
            0.0,
            6,
            &[[0.1, 0.05]],
        )
        .unwrap();
        let structure = StructuralModes {
            omegas: vec![900.0],
            shapes: vec![vec![1.0]],
            loss_factor: 0.01,
        };
        let c = assemble_coupling(&structure, &cavity, &[0.01]).unwrap();
        let model = VibroacousticModel::try_new(&structure, &cavity, c, None)
            .unwrap()
            .with_cavity_relaxation(relaxation)
            .unwrap();
        model.frf(850.0, &[C64::ONE]).unwrap()
    };
    let dry = solve("mixture 293.15 101325 mole 3 N2 0.781 O2 0.209 Ar 0.01");
    let wet = solve("mixture 293.15 101325 mole 4 N2 0.77319 O2 0.20691 Ar 0.0099 H2O 0.01");
    assert!(dry.power.cavity > 0.0 && wet.power.cavity > 0.0);
    assert!((dry.b[0] - wet.b[0]).abs() > 1e-6 * dry.b[0].abs());
    assert!((dry.power.cavity / wet.power.cavity - 1.0).abs() > 0.01);
    assert!(dry.power.interface_residual.abs() < 1e-6);
    assert!(wet.power.interface_residual.abs() < 1e-6);
}
