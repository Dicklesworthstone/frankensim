use super::*;
use crate::helmholtz::bc::{PlaneWave, incident::IncidentField};
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};

fn source(relaxing: bool) -> MolecularExterior {
    let text = if relaxing {
        "declared 1 synthetic 1 vibration 1 100 1 power 101.325 293.15 0 200 500"
    } else {
        "declared 1 synthetic 0"
    };
    let law = read_relaxation(&mut text.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [100.0, 2000.0]).unwrap()
}

fn options(surface: &SpherePanels) -> ComplexOptions {
    let n = surface.areas().len();
    ComplexOptions {
        max_panels: n,
        max_bytes: 32 * n * n,
        max_work: 512 * (n * n) as u64,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-30),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_original_caloric_heat_closes_coherent_wave_energy_balance() {
    let material = source(true);
    let omega = 1000.0;
    let rho_omega = material.source().gas().density * omega;
    let k = material.source().relaxing_wavenumber(omega).unwrap();
    let bulk = material
        .source()
        .caloric_response(omega)
        .unwrap()
        .bulk_modulus_pa;
    let root_loss = (0.5 * omega * bulk.recip().im).sqrt();
    for z in [0.0, 0.2, 0.8] {
        let phase = C64::new((k.re * z).cos(), (k.re * z).sin()).scale((-k.im * z).exp());
        let forward = C64::new(0.8, 0.3) * phase;
        let backward = C64::new(-0.4, 0.7) / phase;
        let p = forward + backward;
        let v = (forward - backward) * k.scale(1.0 / rho_omega);
        let dp = (forward - backward) * k * C64::new(0.0, 1.0);
        let dv = p * k * k * C64::new(0.0, 1.0 / rho_omega);
        let power = density(p, [C64::ZERO, C64::ZERO, v], root_loss).unwrap();
        // Independent conservation identity from the differential equations,
        // including standing-wave interference and reactive particle velocity.
        let minus_divergence = -0.5 * (dp * v.conj() + p * dv.conj()).re;
        let compression_heat = -0.5 * omega * bulk.im * (p / bulk).norm_sq();
        close(power.absorbed_power_density_w_m3, minus_divergence);
        close(power.absorbed_power_density_w_m3, compression_heat);
        assert!(power.absorbed_power_density_w_m3 > 0.0);
    }
}

#[test]
fn g1_spherical_shell_volume_heat_converges_to_lost_radial_flux() {
    let material = source(true);
    let omega = 1000.0;
    let rho_omega = material.source().gas().density * omega;
    let k = material.source().relaxing_wavenumber(omega).unwrap();
    let bulk = material
        .source()
        .caloric_response(omega)
        .unwrap()
        .bulk_modulus_pa;
    let root_loss = (0.5 * omega * bulk.recip().im).sqrt();
    let field = |r: f64| {
        let p = C64::new((k.re * r).cos(), (k.re * r).sin()).scale((-k.im * r).exp() / r);
        let v = p * C64::new(-1.0 / r - k.im, k.re) * C64::new(0.0, -1.0 / rho_omega);
        density(p, [v, C64::ZERO, C64::ZERO], root_loss).unwrap()
    };
    let area = |r: f64| 4.0 * std::f64::consts::PI * r * r;
    let (a, b) = (0.2, 2.0);
    let lost_flux =
        area(a) * field(a).active_intensity_w_m2[0] - area(b) * field(b).active_intensity_w_m2[0];
    let mut previous = f64::INFINITY;
    for n in [8, 16, 32] {
        let step = (b - a) / n as f64;
        let heat = (0..=n)
            .map(|i| {
                let r = a + i as f64 * step;
                let weight = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                weight * area(r) * field(r).absorbed_power_density_w_m3
            })
            .sum::<f64>()
            * step
            / 3.0;
        let error = (heat - lost_flux).abs();
        assert!(
            error < previous / 12.0,
            "Simpson refinement: {error} vs {previous}"
        );
        previous = error;
    }
    assert!(previous < 1e-7 * lost_flux);
}

#[test]
fn g0_g3_source_bound_observers_reuse_fields_and_preserve_lossless_limit() {
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let velocity = vec![C64::new(0.4, -0.2); surface.areas().len()];
    let points = [[0.0, 0.0, 0.8], [0.0, 0.0, 1.0]];
    for relaxing in [false, true] {
        let material = source(relaxing);
        let solution = material
            .solve_batch(
                &surface,
                1000.0,
                &[&velocity],
                Formulation::BurtonMiller,
                &options(&surface),
                || false,
            )
            .unwrap()
            .remove(0);
        let controls = MolecularNearFieldOptions::default();
        let result = solution
            .observe_power(&surface, &points, 0.1, controls, || false)
            .unwrap();
        let prepared = solution
            .prepare_near_field(&surface, &points, 0.1, controls, || false)
            .unwrap();
        let reused = prepared.evaluate_power(&solution, || false).unwrap();
        assert_eq!(
            result.acoustic.scalar.pressure,
            reused.acoustic.scalar.pressure
        );
        assert_eq!(
            result.acoustic.particle_velocity_m_s,
            reused.acoustic.particle_velocity_m_s
        );
        assert_eq!(result.power, reused.power);
        let bulk = material
            .source()
            .caloric_response(1000.0)
            .unwrap()
            .bulk_modulus_pa;
        for (i, row) in result.power.iter().enumerate() {
            let p = result.acoustic.scalar.pressure[i];
            close(
                row.absorbed_power_density_w_m3,
                -500.0 * bulk.im * (p / bulk).norm_sq(),
            );
            assert!(row.active_intensity_w_m2[2] > 0.0);
            assert_eq!(row.absorbed_power_density_w_m3 == 0.0, !relaxing);
        }
        assert!(prepared.evaluate_power(&solution, || true).is_err());
    }
}

#[test]
fn g3_scattering_heat_uses_total_coherent_pressure() {
    let material = source(true);
    let surface = SpherePanels::icosphere(0.2, 0).unwrap();
    let incident = IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::ONE,
    });
    let fields = material
        .solve_scattering_batch(
            &surface,
            1000.0,
            &[incident],
            Admittance::Rigid,
            Formulation::BurtonMiller,
            &MolecularScatteringOptions {
                boundary: options(&surface),
                minimum_source_clearance_m: 0.1,
                maximum_source_point_products: 10000,
            },
            || false,
        )
        .unwrap();
    let points = [[0.0, 0.0, 0.7]];
    let total = fields[0]
        .total_power_at(&points, 0.1, MolecularNearFieldOptions::default(), || false)
        .unwrap();
    let scattered = fields[0]
        .solution()
        .scattered
        .observe_power(
            &surface,
            &points,
            0.1,
            MolecularNearFieldOptions::default(),
            || false,
        )
        .unwrap();
    let incident = material
        .incident_pressure_velocity_at(1000.0, incident, points[0], || false)
        .unwrap();
    let compliance = material
        .source()
        .caloric_response(1000.0)
        .unwrap()
        .bulk_modulus_pa
        .recip();
    let coefficient = 500.0 * compliance.im;
    close(
        total.power[0].absorbed_power_density_w_m3,
        coefficient * total.acoustic.scalar.pressure[0].norm_sq(),
    );
    let incoherent = scattered.power[0].absorbed_power_density_w_m3
        + coefficient * incident.pressure_pa.norm_sq();
    assert!(
        (total.power[0].absorbed_power_density_w_m3 - incoherent).abs() > 1e-3 * incoherent,
        "fixture must resolve coherent heat cross terms"
    );
}

#[test]
fn g0_weak_absorption_preserves_representable_heat_and_rejects_overflow() {
    let power = density(C64::from_re(1e200), [C64::from_re(1e-200); 3], 1e-200).unwrap();
    close(power.absorbed_power_density_w_m3, 1.0);
    assert_eq!(power.active_intensity_w_m2, [0.5; 3]);
    assert!(density(C64::from_re(1e200), [C64::ZERO; 3], 1.0).is_err());
    assert!(density(C64::from_re(f64::NAN), [C64::ZERO; 3], 0.0).is_err());
}
