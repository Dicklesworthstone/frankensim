use super::*;
use crate::compose::ComposeOptions;
use crate::porous::{PorousLayer, PressureCovariance};
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};

const AREA: f64 = 0.004;
fn materials() -> [PorousFluidLaw; 3] {
    std::array::from_fn(|i| {
        let mut card = PorousCard::melamine_uf_nguyen2024();
        card.name = "passive backing numerical fixture";
        card.provenance = "authored constitutive comparison, not measured qualification";
        card.flow_resistivity *= [0.6, 1.8, 1.1][i];
        card.porosity = [0.74, 0.92, 0.85][i];
        PorousFluidLaw::try_new(card, 1.81e-5).unwrap()
    })
}
fn gas(t: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), t, 101_325.0).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-25),
        "{a:e} != {b:e}"
    );
}
fn complex_close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1.0),
        "{a:?} != {b:?}"
    );
}
fn wave(q: C64, x: f64) -> C64 {
    C64::new((q.re * x).cos(), (q.re * x).sin()).scale((-q.im * x).exp())
}
fn product(a: [C64; 4], b: [C64; 4]) -> [C64; 4] {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}

// Independent pressure/normal-velocity transfer from ORIGINAL rho and K.
// The modest optical depths permit growing waves in this test-only oracle.
// Production uses only bounded scattering and positive analytic integrals.
fn reference(
    layers: &[PorousLayer<'_>],
    incident: C64,
    load: PassivePorousTermination,
) -> (C64, Vec<[f64; 3]>, [C64; 2]) {
    let omega = layers[0].material_response().angular_frequency_rad_s();
    let kt = layers[0].tangential_wavenumber_rad_m();
    let r = layers[0].two_port().reference_impedance;
    let mut total = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    let mut coefficients = Vec::new();
    for layer in layers {
        let rho = layer.material_response().dynamic_density_kg_m3().conj();
        let bulk = layer.material_response().dynamic_bulk_modulus_pa().conj();
        let q = ((rho / bulk).scale(omega * omega) - C64::from_re(kt * kt)).sqrt();
        let z = rho.scale(omega) / q;
        let e = wave(q, layer.thickness_m());
        let c = (e + e.recip()).scale(0.5);
        let h = (e.recip() - e).scale(0.5);
        total = product(total, [c, z * h, h / z, c]);
        coefficients.push((rho, bulk, q, z));
    }
    let terminal = load.impedance_pa_s_m3().map_or([C64::ONE, C64::ZERO], |z| {
        [z.scale(1.0 / r), C64::from_re(1.0 / (r * AREA))]
    });
    let start = [
        total[0] * terminal[0] + total[1] * terminal[1],
        total[2] * terminal[0] + total[3] * terminal[1],
    ];
    let amplitude = incident.scale(2.0) / (start[0] + start[1].scale(r * AREA));
    let mut p = start[0] * amplitude;
    let mut v = start[1] * amplitude;
    let reflected = p - incident;
    let mut work = Vec::new();
    for (layer, (rho, bulk, q, z)) in layers.iter().zip(coefficients) {
        let f = (p + z * v).scale(0.5);
        let b = (p - z * v).scale(0.5);
        let n = 2048;
        let dx = layer.thickness_m() / n as f64;
        let mut sum = [0.0; 3];
        for i in 0..=n {
            let x = layer.thickness_m() * i as f64 / n as f64;
            let forward = f * wave(q, x);
            let backward = b * wave(q, -x);
            p = forward + backward;
            v = (forward - backward) / z;
            let vt = p / rho * C64::from_re(kt / omega);
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            for (sum, value) in sum.iter_mut().zip([
                0.5 * omega * rho.im * v.norm_sq(),
                0.5 * omega * rho.im * vt.norm_sq(),
                0.5 * omega * bulk.recip().im * p.norm_sq(),
            ]) {
                *sum += weight * value;
            }
        }
        work.push(sum.map(|x| x * AREA * dx / 3.0));
    }
    complex_close(p, terminal[0] * amplitude, 2e-10);
    complex_close(v, terminal[1] * amplitude, 2e-9);
    (reflected, work, [p, v])
}

#[test]
fn g1_terminated_oblique_stack_matches_independent_density_bulk_volume_work() {
    let materials = materials();
    let exterior = gas(293.15);
    let states = [gas(293.15), gas(347.0), gas(412.0)];
    let incident = C64::new(170.0, -65.0);
    for omega in [4400.0, 13200.0] {
        for mu in [1.0_f64, 0.55, 0.08] {
            let kt = omega / exterior.sound_speed * ((1.0 - mu) * (1.0 + mu)).sqrt();
            let r = exterior.characteristic_impedance / (AREA * mu);
            let layers: [PorousLayer<'_>; 3] = std::array::from_fn(|i| {
                PorousLayer::try_new_oblique(
                    &materials[i],
                    &states[i],
                    omega,
                    [0.017, 0.023, 0.011][i],
                    AREA,
                    r,
                    kt,
                )
                .unwrap()
            });
            let stack =
                PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
            for load in [
                PassivePorousTermination::rigid(),
                PassivePorousTermination::try_new_impedance(C64::ZERO).unwrap(),
                PassivePorousTermination::try_new_impedance(C64::new(0.3 * r, 0.7 * r)).unwrap(),
                PassivePorousTermination::try_new_impedance(C64::new(0.0, 1.1 * r)).unwrap(),
                PassivePorousTermination::try_new_impedance(C64::from_re(r)).unwrap(),
            ] {
                let solved = stack
                    .prepare_terminated_excitation(incident, load, || false)
                    .unwrap();
                let (reflected, work, terminal) = reference(&layers, incident, load);
                complex_close(
                    solved.boundary_power().outgoing_pressure_pa[0],
                    reflected,
                    2e-11,
                );
                for (field, expected) in solved.layer_excitations().iter().zip(&work) {
                    let actual = field
                        .dissipation_between(0.0, field.layer().thickness_m())
                        .unwrap();
                    close(actual.viscous_power_w, expected[0] + expected[1], 2e-10);
                    close(actual.thermal_power_w, expected[2], 2e-10);
                }
                let account = solved.source_power();
                close(
                    account.incident_power_w,
                    incident.norm_sq() / (2.0 * r),
                    3e-15,
                );
                close(
                    account.outgoing_power_w,
                    reflected.norm_sq() / (2.0 * r),
                    2e-11,
                );
                let expected_load = load
                    .impedance_pa_s_m3()
                    .map_or(0.0, |z| 0.5 * z.re * terminal[1].scale(AREA).norm_sq());
                close(account.load_absorbed_power_w, expected_load, 2e-9);
                close(
                    account.incident_power_w,
                    account.outgoing_power_w
                        + account.load_absorbed_power_w
                        + account.porous_absorbed_power_w,
                    5e-13,
                );
                if mu < 0.6 {
                    assert!(
                        work.iter().map(|p| p[1]).sum::<f64>()
                            > 0.01 * account.porous_absorbed_power_w
                    );
                }
            }
        }
    }
}

#[test]
fn g3_weak_passive_loss_and_rigid_matched_limits_keep_source_power_distinct() {
    let materials = materials();
    let gas = gas(293.15);
    let r = gas.characteristic_impedance / AREA;
    let layers = [PorousLayer::try_new(&materials[0], &gas, 6800.0, 0.04, AREA, r).unwrap()];
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    let a = C64::new(200.0, 30.0);
    let open = stack.prepare_excitation([a, C64::ZERO], || false).unwrap();
    let matched = stack
        .prepare_terminated_excitation(
            a,
            PassivePorousTermination::try_new_impedance(C64::from_re(r)).unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(open.dissipation(), matched.dissipation());
    assert_eq!(
        open.at_depth(0.013).unwrap(),
        matched.at_depth(0.013).unwrap()
    );
    close(
        matched.source_power().load_absorbed_power_w,
        open.boundary_power().outgoing_pressure_pa[1].norm_sq() / (2.0 * r),
        4e-15,
    );
    let rigid = stack
        .prepare_terminated_excitation(a, PassivePorousTermination::rigid(), || false)
        .unwrap();
    let surface = layers[0]
        .material_response()
        .surface_impedance_rigid(0.04)
        .unwrap()
        .conj()
        .scale(1.0 / AREA);
    complex_close(
        rigid.boundary_power().outgoing_pressure_pa[0],
        (surface - C64::from_re(r)) / (surface + C64::from_re(r)) * a,
        2e-12,
    );
    assert!(rigid.boundary_power().incident_power_w > rigid.source_power().incident_power_w);
    assert_eq!(rigid.source_power().load_absorbed_power_w, 0.0);
    for real in [1e-20 * r, 1e-12 * r, 0.1 * r] {
        let load = PassivePorousTermination::try_new_impedance(C64::new(real, 0.7 * r)).unwrap();
        let excitation = stack
            .prepare_terminated_excitation(a, load, || false)
            .unwrap();
        let flux = excitation.boundary_power().inward_volume_flow_m3_s[1];
        let expected = 0.5 * real * flux.norm_sq();
        assert!(excitation.source_power().load_absorbed_power_w > 0.0);
        close(
            excitation.source_power().load_absorbed_power_w,
            expected,
            2e-14,
        );
    }
}

#[test]
fn g0_termination_refuses_active_and_ambiguous_sources_and_bounds_feedback() {
    for z in [
        C64::new(-1.0, 0.0),
        C64::new(f64::INFINITY, 0.0),
        C64::new(0.0, f64::NAN),
    ] {
        assert!(PassivePorousTermination::try_new_impedance(z).is_err());
    }
    let materials = materials();
    let gas = gas(293.15);
    let r = gas.characteristic_impedance / AREA;
    let layers = [PorousLayer::try_new(&materials[0], &gas, 6800.0, 0.04, AREA, r).unwrap()];
    let stack = PorousLayerStack::try_new(
        &layers,
        ComposeOptions {
            maximum_condition: 1.0,
            ..ComposeOptions::default()
        },
        || false,
    )
    .unwrap();
    assert!(matches!(
        stack.prepare_terminated_excitation(C64::ONE, PassivePorousTermination::rigid(), || false),
        Err(DuctError::Singular)
    ));
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    let erased_loss =
        PassivePorousTermination::try_new_impedance(C64::new(f64::from_bits(1), r)).unwrap();
    assert!(
        stack
            .prepare_terminated_excitation(C64::ONE, erased_loss, || false)
            .is_err()
    );
    assert!(matches!(
        stack.prepare_terminated_excitation(C64::ONE, PassivePorousTermination::rigid(), || true),
        Err(DuctError::Cancelled { .. })
    ));
    assert!(
        stack
            .prepare_terminated_covariance(
                PressureCovariance::from_peak_pressure([C64::ONE, C64::ONE]).unwrap(),
                PassivePorousTermination::rigid(),
                || false
            )
            .is_err()
    );
    assert_eq!(
        stack
            .prepare_terminated_excitation(C64::ZERO, PassivePorousTermination::rigid(), || false)
            .unwrap()
            .source_power()
            .incident_power_w,
        0.0
    );
    let layers = [PorousLayer::try_new(&materials[0], &gas, 6800.0, 1000.0, AREA, r).unwrap()];
    let stack = PorousLayerStack::try_new(&layers, ComposeOptions::default(), || false).unwrap();
    let excitation = stack
        .prepare_terminated_excitation(C64::ONE, PassivePorousTermination::rigid(), || false)
        .unwrap();
    assert!(excitation.dissipation().absorbed_power_w > 0.0);
    assert_eq!(excitation.boundary_power().pressure_pa[1], C64::ZERO);
}
