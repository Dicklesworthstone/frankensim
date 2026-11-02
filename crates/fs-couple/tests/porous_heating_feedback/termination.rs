use super::*;
use fs_couple::porous_heating::feedback::{
    covariance::{
        PorousCovarianceBin, PorousCrossSpectrumHeating,
        angular::{DiffuseAngularNode, DiffuseHemisphere, PorousAngularMode},
    },
    spectrum::{PorousSpectrumHeating, PorousSpectrumOptions},
};
use fs_duct::porous::{PassivePorousTermination, PressureCovariance};

const LENGTHS: [f64; 2] = [0.03, 0.05];
const INITIAL: [f64; 2] = [0.0, 8000.0];
const AUTHORED: [f64; 2] = [0.4, -0.2];
const CONTACT: f64 = 0.1;

fn source() -> (PorousHarmonicDrive, f64, PassivePorousTermination) {
    let mut d = drive();
    let mu = 0.63_f64;
    let exterior = GasState::try_new(&SPEC, 300.0, 101325.0).unwrap();
    d.reference_impedance_pa_s_m3 /= mu;
    d.incident_pressure_pa = [C64::new(2200.0, -400.0), C64::ZERO];
    let kt = d.angular_frequency_rad_s / exterior.sound_speed * ((1.0 - mu) * (1.0 + mu)).sqrt();
    let load = PassivePorousTermination::try_new_impedance(
        C64::new(0.8, -0.45).scale(d.reference_impedance_pa_s_m3),
    )
    .unwrap();
    (d, kt, load)
}
fn initial() -> DepositionSimulation<'static> {
    let mut e = events(INITIAL, AUTHORED);
    e.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 1,
            a: id(0),
            b: id(1),
            area_m2: AREA,
            conductance_w_k: CONTACT,
        }),
    });
    thermal(e, [&SENSIBLE; 2])
}
fn wave(q: C64, x: f64) -> C64 {
    C64::new((q.re * x).cos(), (q.re * x).sin()).scale((-q.im * x).exp())
}

// Independent [p,v_n] transfer and direct volume integration from rho and K.
// No native acoustic field, scattering, termination or work API supplies this
// reference. All growing-wave operations are test-only at modest optical depth.
fn physical_power(
    h: [f64; 2],
    materials: &[PorousFluidLaw; 2],
    d: PorousHarmonicDrive,
    kt: f64,
    load: PassivePorousTermination,
    panels: usize,
) -> [f64; 5] {
    let omega = d.angular_frequency_rad_s;
    let r = d.reference_impedance_pa_s_m3;
    let coeff: [(C64, C64, C64, C64); 2] = std::array::from_fn(|i| {
        let state = GasState::try_new(&SPEC, 300.0 + h[i] / CP, 101325.0).unwrap();
        let response = materials[i].effective_fluid(&state, omega).unwrap();
        let rho = response.dynamic_density_kg_m3().conj();
        let bulk = response.dynamic_bulk_modulus_pa().conj();
        let q = ((rho / bulk).scale(omega * omega) - C64::from_re(kt * kt)).sqrt();
        (rho, bulk, q, rho.scale(omega) / q)
    });
    let mut total = [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE];
    for i in 0..2 {
        let (_, _, q, z) = coeff[i];
        let e = wave(q, LENGTHS[i]);
        let c = (e + e.recip()).scale(0.5);
        let s = (e.recip() - e).scale(0.5);
        let t = [c, z * s, s / z, c];
        total = [
            total[0] * t[0] + total[1] * t[2],
            total[0] * t[1] + total[1] * t[3],
            total[2] * t[0] + total[3] * t[2],
            total[2] * t[1] + total[3] * t[3],
        ];
    }
    let terminal = load.impedance_pa_s_m3().map_or([C64::ONE, C64::ZERO], |z| {
        [z.scale(1.0 / r), C64::from_re(1.0 / (r * AREA))]
    });
    let inlet = [
        total[0] * terminal[0] + total[1] * terminal[1],
        total[2] * terminal[0] + total[3] * terminal[1],
    ];
    let scale = d.incident_pressure_pa[0].scale(2.0) / (inlet[0] + inlet[1].scale(r * AREA));
    let mut p = inlet[0] * scale;
    let mut v = inlet[1] * scale;
    let reflected = p - d.incident_pressure_pa[0];
    let mut heat = [0.0; 2];
    for i in 0..2 {
        let (rho, bulk, q, z) = coeff[i];
        let forward = (p + z * v).scale(0.5);
        let backward = (p - z * v).scale(0.5);
        let mut sum = 0.0;
        for j in 0..=panels {
            let x = LENGTHS[i] * j as f64 / panels as f64;
            let f = forward * wave(q, x);
            let b = backward * wave(q, -x);
            p = f + b;
            v = (f - b) / z;
            let vt = (p / rho).scale(kt / omega);
            let work = 0.5
                * omega
                * (rho.im * (v.norm_sq() + vt.norm_sq()) + bulk.recip().im * p.norm_sq());
            let weight = if j == 0 || j == panels {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight * work;
        }
        heat[i] = AREA * LENGTHS[i] * sum / (3.0 * panels as f64);
    }
    let load_w = load
        .impedance_pa_s_m3()
        .map_or(0.0, |z| 0.5 * z.re * v.scale(AREA).norm_sq());
    [
        d.incident_pressure_pa[0].norm_sq() / (2.0 * r),
        reflected.norm_sq() / (2.0 * r),
        load_w,
        heat[0],
        heat[1],
    ]
}
fn rates(
    y: [f64; 6],
    materials: &[PorousFluidLaw; 2],
    d: PorousHarmonicDrive,
    kt: f64,
    load: PassivePorousTermination,
) -> [f64; 6] {
    let p = physical_power([y[0], y[1]], materials, d, kt, load, 64);
    let exchange = CONTACT * (y[1] - y[0]) / CP;
    [
        (p[3] + AUTHORED[0] + exchange) / (100.0 * AREA * LENGTHS[0]),
        (p[4] + AUTHORED[1] - exchange) / (100.0 * AREA * LENGTHS[1]),
        p[1],
        p[2],
        p[3],
        p[4],
    ]
}
fn rk4(
    materials: &[PorousFluidLaw; 2],
    d: PorousHarmonicDrive,
    kt: f64,
    load: PassivePorousTermination,
    steps: usize,
) -> [f64; 6] {
    let dt = 0.5 / steps as f64;
    let mut y = [INITIAL[0], INITIAL[1], 0.0, 0.0, 0.0, 0.0];
    for _ in 0..steps {
        let a = rates(y, materials, d, kt, load);
        let b = rates(
            std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]),
            materials,
            d,
            kt,
            load,
        );
        let c = rates(
            std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]),
            materials,
            d,
            kt,
            load,
        );
        let e = rates(
            std::array::from_fn(|i| y[i] + dt * c[i]),
            materials,
            d,
            kt,
            load,
        );
        y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + e[i]) / 6.0);
    }
    y
}

#[test]
fn g1_terminated_native_feedback_uses_accepted_start_load_work_and_refines_independent_fields() {
    let materials = materials();
    let routes = cells(&materials);
    let initial = initial();
    let (d, kt, load) = source();
    let heating = TemperaturePorousHeating::try_new_oblique_terminated(
        &initial,
        &routes,
        gas(),
        d,
        kt,
        load,
        options(),
        || Ok(()),
    )
    .unwrap();
    let expected = physical_power(INITIAL, &materials, d, kt, load, 256);
    let start = heating.power();
    close(start.incident_power_w, expected[0], 1e-12);
    close(start.outgoing_power_w, expected[1], 2e-10);
    close(start.load_absorbed_power_w, expected[2], 2e-10);
    for i in 0..2 {
        close(
            heating.cells()[i].dissipation.absorbed_power_w,
            expected[3 + i],
            2e-9,
        );
    }
    assert!(heating.boundary_power().incident_power_w > start.incident_power_w);
    let dt = 0.0625;
    let step = heating
        .advance_to(&initial, dt, config(dt), || Ok(()))
        .unwrap();
    close(
        step.step.acoustic.incident_j,
        start.incident_power_w * dt,
        1e-12,
    );
    close(
        step.step.acoustic.outgoing_j,
        start.outgoing_power_w * dt,
        1e-12,
    );
    close(
        step.step.acoustic.load_absorbed_j,
        start.load_absorbed_power_w * dt,
        1e-12,
    );
    close(
        step.step.acoustic.deposited_j,
        start.dissipation.absorbed_power_w * dt,
        1e-12,
    );
    assert_eq!(step.work.gas_evaluations, 2);
    assert_eq!(step.work.stack_solves, 1);
    assert_ne!(step.heating.power(), start);
    // The frozen-stack adapter uses the same source/load account, too. One
    // native step must deposit the identical accepted-start pore work.
    let stack = PorousLayerStack::try_new(heating.acoustic_layers(), options(), || false).unwrap();
    let excitation = stack
        .prepare_terminated_excitation(d.incident_pressure_pa[0], load, || false)
        .unwrap();
    let partitions = [
        [PorousHeatCell {
            cell: id(0),
            end_depth_m: LENGTHS[0],
        }],
        [PorousHeatCell {
            cell: id(1),
            end_depth_m: LENGTHS[1],
        }],
    ];
    let frozen = PorousStackHeating::try_new(
        &excitation,
        &initial,
        &[&partitions[0], &partitions[1]],
        2,
        || Ok(()),
    )
    .unwrap();
    let frozen_step = frozen
        .advance_to(&initial, dt, config(dt), || Ok(()))
        .unwrap();
    assert_eq!(frozen_step.thermal.cells(), step.step.thermal.cells());
    assert_eq!(frozen_step.acoustic, step.step.acoustic);
    let endpoint_h = std::array::from_fn(|i| step.step.thermal.cells()[&id(i)].specific_h_j_kg);
    let endpoint = physical_power(endpoint_h, &materials, d, kt, load, 256);
    close(
        step.heating.power().load_absorbed_power_w,
        endpoint[2],
        2e-10,
    );
    for i in 0..2 {
        close(
            step.heating.cells()[i].dissipation.absorbed_power_w,
            endpoint[3 + i],
            2e-9,
        );
    }
    let reference = rk4(&materials, d, kt, load, 512);
    let refined = rk4(&materials, d, kt, load, 1024);
    for i in 0..6 {
        close(reference[i], refined[i], 2e-9);
    }
    let mut errors = Vec::new();
    for n in [8, 32, 128] {
        let result = heating
            .advance_to(&initial, 0.5, config(0.5 / n as f64), || Ok(()))
            .unwrap();
        let actual = [
            result.step.thermal.cells()[&id(0)].specific_h_j_kg,
            result.step.thermal.cells()[&id(1)].specific_h_j_kg,
            result.step.acoustic.outgoing_j,
            result.step.acoustic.load_absorbed_j,
            result.step.cell_energy_j[0].1,
            result.step.cell_energy_j[1].1,
        ];
        errors.push(
            (0..6)
                .map(|i| (actual[i] - reference[i]).abs() / reference[i].abs().max(1.0))
                .fold(0.0, f64::max),
        );
        close(result.step.acoustic.residual_j, 0.0, 3e-10);
        assert_eq!(result.work.gas_evaluations, 2 * result.work.stack_solves);
    }
    assert!(
        errors[1] < 0.3 * errors[0] && errors[2] < 0.3 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 3e-4, "{errors:?}");
    assert_eq!(initial.cells()[&id(0)].specific_h_j_kg, INITIAL[0]);
    eprintln!(
        "terminated native refinement errors={errors:?}, load_W={:.8}, prescribed_W={:.8}",
        start.load_absorbed_power_w, start.incident_power_w
    );
}

#[test]
fn g3_terminated_spectrum_covariance_and_angular_modes_share_original_native_account() {
    let materials = materials();
    let routes = cells(&materials);
    let initial = initial();
    let (d, kt, load) = source();
    let harmonic = TemperaturePorousHeating::try_new_oblique_terminated(
        &initial,
        &routes,
        gas(),
        d,
        kt,
        load,
        options(),
        || Ok(()),
    )
    .unwrap();
    let spectral = PorousSpectrumHeating::try_new_oblique_terminated(
        &initial,
        &routes,
        gas(),
        &[d],
        &[kt],
        &[load],
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let bin = PorousCovarianceBin {
        angular_frequency_rad_s: d.angular_frequency_rad_s,
        area_m2: AREA,
        reference_impedance_pa_s_m3: d.reference_impedance_pa_s_m3,
        covariance: PressureCovariance::from_peak_pressure(d.incident_pressure_pa).unwrap(),
    };
    let cov = PorousCrossSpectrumHeating::try_new_oblique_terminated(
        &initial,
        &routes,
        gas(),
        &[bin],
        &[kt],
        &[load],
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let mode = PorousAngularMode {
        bin,
        tangential_wavenumber_rad_m: kt,
    };
    let angular = PorousCrossSpectrumHeating::try_new_independent_modes_terminated(
        &initial,
        &routes,
        gas(),
        &[mode],
        &[load],
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(harmonic.power(), spectral.power());
    assert_eq!(harmonic.power(), cov.power());
    assert_eq!(harmonic.power(), angular.power());
    let h = harmonic
        .advance_to(&initial, 0.25, config(0.0625), || Ok(()))
        .unwrap();
    let s = spectral
        .advance_to(&initial, 0.25, config(0.0625), || Ok(()))
        .unwrap();
    let c = cov
        .advance_to(&initial, 0.25, config(0.0625), || Ok(()))
        .unwrap();
    let a = angular
        .advance_to(&initial, 0.25, config(0.0625), || Ok(()))
        .unwrap();
    assert_eq!(h.step.thermal.cells(), s.step.thermal.cells());
    assert_eq!(h.step.thermal.cells(), c.step.thermal.cells());
    assert_eq!(h.step.thermal.cells(), a.step.thermal.cells());
    assert_eq!(h.step.acoustic, s.step.acoustic);
    assert_eq!(h.step.acoustic, c.step.acoustic);
    assert_eq!(h.step.acoustic, a.step.acoustic);
    // Normal-incidence overloads retain the identical core.
    let normal = TemperaturePorousHeating::try_new_terminated(
        &initial,
        &routes,
        gas(),
        d,
        load,
        options(),
        || Ok(()),
    )
    .unwrap();
    let normal_s = PorousSpectrumHeating::try_new_terminated(
        &initial,
        &routes,
        gas(),
        &[d],
        &[load],
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let normal_c = PorousCrossSpectrumHeating::try_new_terminated(
        &initial,
        &routes,
        gas(),
        &[bin],
        &[load],
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(normal.power(), normal_s.power());
    assert_eq!(normal.power(), normal_c.power());
    // Different loads remain attached to their input mode across canonical
    // sorting, including two distinct tangential modes at one frequency.
    let normal_mode = PorousAngularMode {
        bin,
        tangential_wavenumber_rad_m: 0.0,
    };
    let rigid = PassivePorousTermination::rigid();
    let pair = PorousCrossSpectrumHeating::try_new_independent_modes_terminated(
        &initial,
        &routes,
        gas(),
        &[mode, normal_mode],
        &[load, rigid],
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    let reversed = PorousCrossSpectrumHeating::try_new_independent_modes_terminated(
        &initial,
        &routes,
        gas(),
        &[normal_mode, mode],
        &[rigid, load],
        PorousSpectrumOptions::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(pair.power(), reversed.power());
    assert_eq!(pair.cells(), reversed.cells());
    assert_eq!(pair.termination(0), Some(rigid));
    assert_eq!(pair.termination(1), Some(load));
    assert_eq!(pair.bin_power(1), Some(harmonic.power()));
    assert!(
        TemperaturePorousHeating::try_new_terminated(
            &initial,
            &routes,
            gas(),
            drive(),
            load,
            options(),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        PorousSpectrumHeating::try_new_terminated(
            &initial,
            &routes,
            gas(),
            &[d],
            &[],
            PorousSpectrumOptions::default(),
            || Ok(())
        )
        .is_err()
    );
    assert!(
        PorousCrossSpectrumHeating::try_new_independent_modes_terminated(
            &initial,
            &routes,
            gas(),
            &[mode],
            &[],
            PorousSpectrumOptions::default(),
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn g1_rigid_backed_diffuse_heat_refines_angle_and_preserves_shared_gas_work() {
    let materials = materials();
    let routes = cells(&materials);
    let initial = initial();
    let exterior = GasState::try_new(&SPEC, 300.0, 101325.0).unwrap();
    let hemisphere = DiffuseHemisphere::try_new(
        &exterior,
        core::f64::consts::TAU * 1100.0,
        AREA,
        [1200.0, 0.0],
    )
    .unwrap();
    let mut totals = Vec::new();
    for n in [4, 16, 64] {
        let nodes: Vec<_> = (0..n)
            .map(|i| DiffuseAngularNode {
                cosine: (i as f64 + 0.5) / n as f64,
                weight_mu: 1.0 / n as f64,
            })
            .collect();
        let modes = hemisphere.modes(&nodes, n, || Ok(())).unwrap();
        let loads = vec![PassivePorousTermination::rigid(); n];
        let source = PorousCrossSpectrumHeating::try_new_independent_modes_terminated(
            &initial,
            &routes,
            gas(),
            &modes,
            &loads,
            PorousSpectrumOptions::default(),
            || Ok(()),
        )
        .unwrap();
        let mut expected = [0.0; 2];
        for mode in &modes {
            let drive = PorousHarmonicDrive {
                angular_frequency_rad_s: mode.bin.angular_frequency_rad_s,
                area_m2: AREA,
                reference_impedance_pa_s_m3: mode.bin.reference_impedance_pa_s_m3,
                incident_pressure_pa: mode.bin.covariance.peak_components()[0],
            };
            let power = physical_power(
                INITIAL,
                &materials,
                drive,
                mode.tangential_wavenumber_rad_m,
                PassivePorousTermination::rigid(),
                256,
            );
            for i in 0..2 {
                expected[i] += power[3 + i];
            }
        }
        for i in 0..2 {
            close(
                source.cells()[i].dissipation.absorbed_power_w,
                expected[i],
                3e-9,
            );
        }
        close(source.power().incident_power_w, AREA * 1200.0, 2e-14);
        assert_eq!(source.power().load_absorbed_power_w, 0.0);
        totals.push(source.power().dissipation.absorbed_power_w);
        if n == 16 {
            let step = source
                .advance_to(&initial, 0.0625, config(0.0625), || Ok(()))
                .unwrap();
            assert_eq!(step.work.gas_evaluations, 2);
            assert_eq!(step.work.stack_solves, n);
            close(step.step.acoustic.incident_j, AREA * 1200.0 * 0.0625, 1e-13);
            assert_eq!(step.step.acoustic.load_absorbed_j, 0.0);
        }
    }
    assert!(
        (totals[2] - totals[1]).abs() < 0.15 * (totals[1] - totals[0]).abs(),
        "{totals:?}"
    );
    eprintln!("rigid-backed diffuse angular heat W={totals:?}");
}
