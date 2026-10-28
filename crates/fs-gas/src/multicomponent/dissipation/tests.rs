use super::*;
use crate::multicomponent::tests::{bank, source};
use crate::multicomponent::{Boundary, CaloricDomain, MixtureGas};
use fs_qty::{DynViscosity, Length, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

fn gas(count: usize, formation: f64, unequal: bool) -> MixtureGas {
    let mut models = vec![source("A", R, 3.5, 0.0, formation)];
    if count == 2 {
        models.push(source(
            "B",
            if unequal { 2.0 * R } else { R },
            if unequal { 4.0 } else { 3.5 },
            0.0,
            -2.0 * formation,
        ));
    }
    MixtureGas::new(
        bank(models),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}
fn field(
    gas: MixtureGas,
    n: usize,
    boundaries: Boundaries,
    initial: impl Fn(f64) -> (Vec<f64>, f64, f64),
) -> MixtureEuler1d {
    let states = (0..n)
        .map(|i| {
            let (species, u, p) = initial((i as f64 + 0.5) / n as f64);
            gas.conserved(&species, u, p).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / n as f64),
        0.7,
        states,
        boundaries,
        || false,
    )
    .unwrap()
}
fn law(flow: &MixtureEuler1d, mu: f64, bulk: f64, lambda: f64, d: f64) -> ConstantTransport {
    let n = flow.gas().species_count();
    let mut binary = vec![d; n * n];
    for i in 0..n {
        binary[i * n + i] = 0.0;
    }
    ConstantTransport::new(
        flow.gas().source().clone(),
        "manufactured-constant-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(mu),
            bulk_viscosity: DynViscosity::new(bulk),
            thermal_conductivity: ThermalConductivityQuantityV1::new(lambda),
        },
        binary,
    )
    .unwrap()
}
fn operator(
    flow: &MixtureEuler1d,
    mu: f64,
    bulk: f64,
    lambda: f64,
    d: f64,
) -> DissipationOperator<ConstantTransport> {
    DissipationOperator::new(flow, law(flow, mu, bulk, lambda, d), || false).unwrap()
}
fn heat_error(n: usize) -> f64 {
    let mut flow = field(gas(1, 0.0, false), n, Boundaries::Periodic, |x| {
        (
            vec![1.0],
            0.0,
            2.0 + 0.2 * (2.0 * core::f64::consts::PI * x).sin(),
        )
    });
    let mut op = operator(&flow, 0.0, 0.0, 0.05, 0.0);
    let before = flow.totals();
    let report = op
        .advance_for(
            &mut flow,
            Time::new(0.1),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.elapsed.value(), 0.1);
    assert_eq!(flow.time().value(), 0.0);
    assert!(report.accepted_steps > 1);
    assert!(report.maximum_parabolic_number <= 0.40000000000001);
    assert!((flow.totals().energy - before.energy).abs() < 2e-13);
    assert!(report.absolute_energy_residual_j.abs() < 2e-13);
    let amplitude = 0.2 * (-0.02 * (2.0 * core::f64::consts::PI).powi(2) * 0.1).exp();
    let mut error = 0.0;
    for i in 0..n {
        let x = (i as f64 + 0.5) / n as f64;
        let q = flow.cell(i).unwrap();
        let p = flow.gas().primitive(q, || false).unwrap();
        assert_eq!(q.partial_densities, [1.0]);
        error += (p.temperature - (2.0 + amplitude * (2.0 * core::f64::consts::PI * x).sin()))
            .abs()
            / n as f64;
    }
    error
}

#[test]
fn g1_fourier_heat_mode_refines_and_conserves_total_energy_without_advancing_clock() {
    let coarse = heat_error(40);
    let fine = heat_error(80);
    eprintln!("Fourier mode L1: 40={coarse:e},80={fine:e}");
    assert!(coarse < 5e-5);
    assert!(fine < coarse * 0.3);
}

#[test]
fn g0_nonuniform_mesh_uses_physical_distances_and_cell_volumes() {
    let g = gas(1, 0.0, false);
    let widths = [0.2, 0.3, 0.5];
    let temperature = [2.1, 2.35, 2.75];
    let states = temperature
        .iter()
        .map(|t| g.conserved(&[1.0], 0.0, *t).unwrap())
        .collect();
    let mesh = Mesh1d::new(
        vec![0.0, 0.2, 0.5, 1.0],
        vec![0.7; 4],
        widths.iter().map(|w| 0.7 * w).collect(),
    )
    .unwrap();
    let mut flow = MixtureEuler1d::new_on_mesh(
        g,
        mesh,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    let dt = 1e-3;
    let diffusivity = 0.05 / 2.5;
    // Initially linear T gives identical interior conductive fluxes. The two
    // insulated end cells change first; the second RK stage carries their
    // changes through the actual unequal centre-to-centre distances.
    let first = [diffusivity / widths[0], 0.0, -diffusivity / widths[2]];
    let curvature_left = diffusivity * (first[1] - first[0]) / 0.25;
    let curvature_right = diffusivity * (first[2] - first[1]) / 0.4;
    let second = [
        curvature_left / widths[0],
        (curvature_right - curvature_left) / widths[1],
        -curvature_right / widths[2],
    ];
    let mut op = operator(&flow, 0.0, 0.0, 0.05, 0.0);
    let report = op
        .advance_for(
            &mut flow,
            Time::new(dt),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(report.accepted_steps, 1);
    for i in 0..3 {
        let actual = flow
            .gas()
            .primitive(flow.cell(i).unwrap(), || false)
            .unwrap()
            .temperature;
        let expected = temperature[i] + dt * first[i] + 0.5 * dt * dt * second[i];
        assert!((actual - expected).abs() < 2e-12);
    }
    assert!(report.absolute_energy_residual_j.abs() < 1e-13);
}

#[test]
fn g1_maxwell_stefan_binary_mode_matches_analytic_decay_and_keeps_temperature() {
    let n = 80;
    let mut flow = field(gas(2, 0.0, false), n, Boundaries::Periodic, |x| {
        let a = 0.5 + 0.2 * (2.0 * core::f64::consts::PI * x).sin();
        (vec![a, 1.0 - a], 0.0, 2.0)
    });
    let mut op = operator(&flow, 0.0, 0.0, 0.0, 0.03);
    let report = op
        .advance_for(
            &mut flow,
            Time::new(0.1),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let amplitude = 0.2 * (-0.03 * (2.0 * core::f64::consts::PI).powi(2) * 0.1).exp();
    let mut error = 0.0;
    for i in 0..n {
        let x = (i as f64 + 0.5) / n as f64;
        let q = flow.cell(i).unwrap();
        error += (q.partial_densities[0]
            - (0.5 + amplitude * (2.0 * core::f64::consts::PI * x).sin()))
        .abs()
            / n as f64;
        let p = flow.gas().primitive(q, || false).unwrap();
        assert!((p.temperature - 2.0).abs() < 5e-11);
        assert!((p.density - 1.0).abs() < 5e-13);
    }
    assert!(error < 1.5e-5, "binary diffusion error={error:e}");
    assert!(
        report.conservation_residual.species_masses[..2]
            .iter()
            .all(|v| v.abs() < 5e-14)
    );
    assert!(report.maximum_mass_frame_residual_kg_m2_s < 1e-14);
}

#[test]
fn g1_viscous_mode_heats_the_gas_while_conserving_total_energy() {
    let n = 80;
    let mut flow = field(gas(1, 0.0, false), n, Boundaries::Periodic, |x| {
        (
            vec![1.0],
            0.2 * (2.0 * core::f64::consts::PI * x).sin(),
            2.0,
        )
    });
    let initial_kinetic = kinetic(&flow);
    let before = flow.totals();
    let mut op = operator(&flow, 0.03, 0.01, 0.0, 0.0);
    let report = op
        .advance_for(
            &mut flow,
            Time::new(0.1),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    let amplitude = 0.2 * (-0.05 * (2.0 * core::f64::consts::PI).powi(2) * 0.1).exp();
    let mut error = 0.0;
    let mut mean_t = 0.0;
    for i in 0..n {
        let q = flow.cell(i).unwrap();
        let p = flow.gas().primitive(q, || false).unwrap();
        let x = (i as f64 + 0.5) / n as f64;
        error +=
            (p.velocity - amplitude * (2.0 * core::f64::consts::PI * x).sin()).abs() / n as f64;
        mean_t += p.temperature / n as f64;
    }
    assert!(error < 2e-5, "viscous mode error={error:e}");
    assert!(mean_t > 2.0);
    let heat_gain = 0.7 * 2.5 * (mean_t - 2.0);
    assert!((heat_gain - (initial_kinetic - kinetic(&flow))).abs() < 2e-12);
    assert!((flow.totals().energy - before.energy).abs() < 2e-13);
    assert!(report.absolute_energy_residual_j.abs() < 2e-13);
}
fn kinetic(flow: &MixtureEuler1d) -> f64 {
    (0..flow.cell_count())
        .map(|i| {
            let q = flow.cell(i).unwrap();
            0.5 * q.momentum * q.momentum / q.density() * flow.mesh().volumes()[i]
        })
        .sum()
}

#[test]
fn g0_unequal_mass_diffusion_carries_enthalpy_and_is_reference_invariant() {
    let initial = |x: f64| {
        let a = 0.5 + 0.2 * (2.0 * core::f64::consts::PI * x).sin();
        (vec![a, 2.0 * (1.0 - a)], 0.0, 2.0)
    };
    let mut ordinary = field(gas(2, 0.0, true), 32, Boundaries::Periodic, initial);
    let mut shifted = field(gas(2, 1e18, true), 32, Boundaries::Periodic, initial);
    let mut first = operator(&ordinary, 0.0, 0.0, 0.02, 0.03);
    let mut second = operator(&shifted, 0.0, 0.0, 0.02, 0.03);
    let report = first
        .advance_for(
            &mut ordinary,
            Time::new(0.03),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    second
        .advance_for(
            &mut shifted,
            Time::new(0.03),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(ordinary.states, shifted.states);
    assert!(
        report.conservation_residual.species_masses[..2]
            .iter()
            .all(|x| x.abs() < 1e-13)
    );
    assert!(report.absolute_energy_residual_j.abs() < 2e-12);
    assert!(report.maximum_mass_frame_residual_kg_m2_s < 2e-14);
    let temperatures = (0..ordinary.cell_count())
        .map(|i| {
            ordinary
                .gas()
                .primitive(ordinary.cell(i).unwrap(), || false)
                .unwrap()
                .temperature
        })
        .collect::<Vec<_>>();
    assert!(
        temperatures.iter().any(|t| (t - 2.0).abs() > 1e-5),
        "unequal molecular transport must include enthalpy/thermal coupling"
    );
}

#[test]
fn g0_stationary_adiabatic_walls_exert_traction_with_zero_energy_work() {
    let mut flow = field(
        gas(1, 0.0, false),
        24,
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        |_| (vec![1.0], 0.3, 2.0),
    );
    let before = flow.totals();
    let mut op = operator(&flow, 0.03, 0.0, 0.0, 0.0);
    let report = op
        .advance_for(
            &mut flow,
            Time::new(0.03),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert!(report.boundary_change.momentum < 0.0);
    assert_eq!(report.boundary_change.energy, 0.0);
    assert_eq!(report.boundary_change.mass(), 0.0);
    assert!(flow.totals().momentum < before.momentum);
    assert!((flow.totals().energy - before.energy).abs() < 1e-13);
    assert!(report.conservation_residual.momentum.abs() < 1e-13);
}

#[test]
fn g4_late_cancel_work_limit_and_geometry_refusal_preserve_the_field() {
    let setup = || {
        field(gas(1, 0.0, false), 20, Boundaries::Periodic, |x| {
            (
                vec![1.0],
                0.0,
                2.0 + 0.2 * (2.0 * core::f64::consts::PI * x).sin(),
            )
        })
    };
    let mut reference = setup();
    let mut reference_op = operator(&reference, 0.0, 0.0, 0.05, 0.0);
    let mut polls = 0;
    reference_op
        .advance_for(
            &mut reference,
            Time::new(0.1),
            DissipationControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut flow = setup();
    let original = flow.states.clone();
    let mut op = operator(&flow, 0.0, 0.0, 0.05, 0.0);
    let mut count = 0;
    let error = op
        .advance_for(
            &mut flow,
            Time::new(0.1),
            DissipationControl::default(),
            || {
                count += 1;
                count >= polls - 2
            },
        )
        .unwrap_err();
    assert_eq!(error, DissipationError::Cancelled);
    assert_eq!(flow.states, original);
    assert_eq!(flow.time().value(), 0.0);
    let error = op
        .advance_for(
            &mut flow,
            Time::new(0.1),
            DissipationControl {
                maximum_steps: 1,
                ..DissipationControl::default()
            },
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, DissipationError::WorkLimit { .. }));
    assert_eq!(flow.states, original);
    op.advance_for(
        &mut flow,
        Time::new(0.1),
        DissipationControl::default(),
        || false,
    )
    .unwrap();
    for (a, b) in flow.states.iter().zip(&reference.states) {
        assert!((a - b).abs() < 1e-11);
    }
    let state = flow.states.clone();
    let noop = op
        .advance_for(
            &mut flow,
            Time::new(0.0),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(noop.accepted_steps, 0);
    assert_eq!(flow.states, state);
    let g = gas(1, 0.0, false);
    let states = vec![g.conserved(&[1.0], 0.0, 2.0).unwrap(); 2];
    let mesh = Mesh1d::new(vec![0.0, 0.5, 1.0], vec![1.0, 1.1, 1.2], vec![0.525, 0.575]).unwrap();
    let taper = MixtureEuler1d::new_on_mesh(
        g,
        mesh,
        states,
        Boundaries::Ends {
            left: Boundary::Transmissive,
            right: Boundary::Transmissive,
        },
        || false,
    )
    .unwrap();
    assert!(matches!(
        DissipationOperator::new(&taper, law(&taper, 0.0, 0.0, 0.05, 0.0), || false),
        Err(DissipationError::Invalid(_))
    ));
}
