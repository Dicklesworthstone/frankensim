//! G1/G2: sound-driven breathing sphere, with independent incident/scattered waves.
use super::MolecularIncidentLimits;
use crate::bem_loaded_modes::{BemHarmonicOptions, BemLoadedModes, BemModesOptions};
use fs_bem::helmholtz::{
    Formulation, Medium,
    bc::{IncidentField, PlaneWave},
    laplace::complex::ComplexOptions,
    molecular::MolecularExterior,
};
use fs_bem::observation_transfers::{MolecularResponseLimits, ReceiverSpec};
use fs_bem::panel3d::SpherePanels;
use fs_material::{
    gas::{GasSpec, GasState},
    relaxation::{RelaxingMedium, input::read_relaxation},
};
use fs_math::{c64::C64, det};
use fs_phs::RelaxationBranch;
use std::f64::consts::PI;

const OMEGA: f64 = 1000.0;
const RADIUS: f64 = 0.3;
const MASS: f64 = 0.1;
const STIFFNESS: f64 = 600_000.0;
const RESISTANCE: f64 = 6.0;
const EXCESS: f64 = 50_000.0;
const TAU: f64 = 0.0007;

fn exterior() -> MolecularExterior {
    // Synthetic constitutive data exercise significant relaxation, not an
    // empirical atmospheric fit or a claimed sound-pressure calibration.
    let data = "declared 1 synthetic 1 vibration 1 100 1 power 101.325 293.15 0 200 500";
    let law = read_relaxation(&mut data.split_ascii_whitespace()).unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    let medium = RelaxingMedium::new(gas, vec![(law.species()[0].species.clone(), 1.0)]).unwrap();
    MolecularExterior::new(law.evaluate(&medium).unwrap(), [800.0, 1200.0]).unwrap()
}

fn body(level: usize, gas: &GasState) -> BemLoadedModes {
    let surface = SpherePanels::icosphere(RADIUS, level).unwrap();
    let p = surface.areas().len();
    let n = p + 2;
    let weight = 1.0 / det::sqrt(MASS);
    BemLoadedModes::build(
        surface,
        vec![vec![weight; p]],
        vec![STIFFNESS / MASS],
        vec![RESISTANCE / MASS],
        Medium::from(gas),
        BemModesOptions {
            formulation: Formulation::PlainCbie,
            boundary: ComplexOptions {
                max_panels: p,
                max_bytes: 64 * p * p,
                max_work: 512 * (p * p) as u64,
            },
            panel_refinements: 0,
            reference_angular_frequency_rad_s: OMEGA,
            physical_residual_tolerance: 1e-10,
            max_dimension: n,
            max_dense_bytes: 256 * n * n + 1024 * n,
            max_operator_evaluations: 4,
            max_projection_work: 32 * (n * n) as u64,
        },
        &mut || false,
    )
    .unwrap()
    .with_relaxation_branches(
        &[RelaxationBranch {
            projection: vec![weight],
            stiffness: EXCESS,
            relaxation_time_s: TAU,
        }],
        1,
        || false,
    )
    .unwrap()
}

fn exp_i(z: C64) -> C64 {
    C64::new(det::cos(z.re), det::sin(z.re)).scale(det::exp(-z.im))
}

// NIST DLMF 10.53.1 (https://dlmf.nist.gov/10.53.E1) and 10.51.1:
// independent small-complex-argument j_l series and outgoing h_l recurrence.
fn spherical(z: C64, count: usize) -> (Vec<C64>, Vec<C64>) {
    let mut j = Vec::with_capacity(count);
    let mut leading = C64::ONE;
    for n in 0..count {
        if n > 0 {
            leading = leading * z.scale(1.0 / (2 * n + 1) as f64);
        }
        let (mut term, mut sum) = (leading, leading);
        for m in 1..48 {
            term = term * (z * z).scale(-1.0 / (2 * m * (2 * n + 2 * m + 1)) as f64);
            sum = sum + term;
        }
        j.push(sum);
    }
    let ez = exp_i(z) / z;
    let i = C64::new(0.0, 1.0);
    let mut h = vec![-i * ez, -ez * (C64::ONE + i / z)];
    for n in 1..count - 1 {
        h.push(h[n].scale((2 * n + 1) as f64) / z - h[n - 1]);
    }
    (j, h)
}

fn derivative(f: &[C64], z: C64, n: usize) -> C64 {
    if n == 0 {
        -f[1]
    } else {
        f[n - 1] - f[n].scale((n + 1) as f64) / z
    }
}

// NIST DLMF 10.60.7 expands the incident plane wave. Enforce zero TOTAL
// radial derivative at a: a_l=-P(2l+1)i^l j'_l(ka)/h'_l(ka).
// Derivative of the outgoing receiver series supplies velocity through i*rho*w.
fn blocked_receiver(k: C64, rho_omega: f64, amplitude: C64, radius: f64, mu: f64) -> (C64, C64) {
    let ka = k.scale(RADIUS);
    let kr = k.scale(radius);
    let (j, h) = spherical(ka, 24);
    let (_, hr) = spherical(kr, 24);
    let (mut previous, mut legendre) = (1.0, 1.0);
    let (mut pressure, mut gradient) = (C64::ZERO, C64::ZERO);
    let mut phase = C64::ONE;
    for n in 0..24 {
        if n == 1 {
            legendre = mu;
        } else if n > 1 {
            let next = ((2 * n - 1) as f64 * mu * legendre - (n - 1) as f64 * previous) / n as f64;
            previous = legendre;
            legendre = next;
        }
        let a = -amplitude * phase.scale((2 * n + 1) as f64) * derivative(&j, ka, n)
            / derivative(&h, ka, n);
        pressure = pressure + (a * hr[n]).scale(legendre);
        gradient = gradient + (a * k * derivative(&hr, kr, n)).scale(legendre);
        phase = phase * C64::new(0.0, 1.0);
    }
    let incident = amplitude * exp_i(k.scale(radius * mu));
    (
        pressure + incident,
        gradient * C64::new(0.0, -1.0 / rho_omega) + (incident * k).scale(mu / rho_omega),
    )
}

#[test]
fn incident_sound_drives_mass_and_memory_with_total_field_energy_interference() {
    let gas = exterior();
    let k = gas.source().relaxing_wavenumber(OMEGA).unwrap();
    let rho_omega = gas.source().gas().density * OMEGA;
    let amplitude = C64::new(1.0, -0.3);
    let area = 4.0 * PI * RADIUS * RADIUS;
    let s = C64::new(0.0, -OMEGA);
    let ik = C64::new(-k.im, k.re);
    // The l=0 rigid-blocked pressure is P[j0-j0' h0/h0'], simplifying
    // exactly to P exp(-ika)/(1-ika). Bare incident pressure misses reflection.
    let blocked_mean = amplitude * exp_i(-k.scale(RADIUS)) / (C64::ONE - ik.scale(RADIUS));
    let blocked_force = -blocked_mean.scale(area);
    let radiation = C64::new(0.0, rho_omega) / (ik - C64::from_re(1.0 / RADIUS));
    let mechanical = C64::from_re(RESISTANCE)
        + s.scale(MASS)
        + C64::from_re(STIFFNESS) / s
        + C64::from_re(EXCESS * TAU) / (C64::ONE + s.scale(TAU));
    let velocity = blocked_force / (mechanical + radiation.scale(area));
    let memory = velocity.scale(det::sqrt(EXCESS) * TAU) / (C64::ONE + s.scale(TAU));
    let heat = 0.5 * RESISTANCE * velocity.norm_sq() + 0.5 * memory.norm_sq() / TAU;
    let self_work = 0.5 * area * radiation.re * velocity.norm_sq();
    let sound_work = 0.5 * (blocked_force * velocity.conj()).re;
    assert!(heat > 0.0 && self_work > 0.0 && sound_work > 0.0);
    assert!((sound_work - self_work - heat).abs() < 1e-13 * sound_work);
    assert!(k.im > 0.01 * k.re);
    let expected_total_reaction = -(blocked_mean + radiation * velocity).scale(area);
    let receivers = [
        ReceiverSpec {
            id: "forward".into(),
            position_m: [0.0, 0.0, 0.9],
        },
        ReceiverSpec {
            id: "transverse".into(),
            position_m: [0.0, 1.2, 0.0],
        },
    ];
    let expected = [(0.9, 1.0), (1.2, 0.0)].map(|(r, mu)| {
        let (p, u) = blocked_receiver(k, rho_omega, amplitude, r, mu);
        let moving_p = radiation * velocity * exp_i(k.scale(r - RADIUS)).scale(RADIUS / r);
        (
            p + moving_p,
            u + moving_p * (ik - C64::from_re(1.0 / r)) / C64::new(0.0, rho_omega),
        )
    });
    let illumination = IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: amplitude,
    });
    let mut errors = Vec::new();
    for level in 0..=2 {
        let mechanics = body(level, gas.source().gas());
        let p = mechanics.panel_count();
        let prepared = mechanics
            .prepare_molecular_harmonic(
                &gas,
                OMEGA,
                &receivers,
                0.05,
                None,
                MolecularResponseLimits {
                    maximum_retained_entries: 2 * p + 32,
                    maximum_panel_products: 3 * p as u64,
                },
                BemHarmonicOptions {
                    max_load_cases: 1,
                    max_response_bytes: 131_072,
                    max_work: 1_000_000,
                    power_relative_tolerance: 1e-10,
                },
                || false,
            )
            .unwrap();
        let result = prepared
            .harmonic_incident_response_batch(
                &[illumination],
                &[&[C64::ZERO]],
                MolecularIncidentLimits {
                    minimum_source_clearance_m: 1e-4,
                    maximum_source_point_products: 10_000_000,
                    maximum_live_bytes: 8_000_000,
                    maximum_panel_products: 10_000_000,
                },
                || false,
            )
            .unwrap()
            .remove(0);
        let actual_velocity = result.motion.modal_velocity[0].scale(1.0 / det::sqrt(MASS));
        let actual_force = result.blocked_acoustic_force[0].scale(det::sqrt(MASS));
        let power = result.power;
        assert_eq!(power.external_mechanical_input_w, 0.0);
        assert!(power.kelvin_w > 0.0 && power.relaxation_w > 0.0);
        assert!(
            power.total_exterior_w < 0.0,
            "the sound-driven body absorbs energy from its total exterior field"
        );
        let structural_heat = power.kelvin_w + power.relaxation_w;
        assert!(
            (power.total_exterior_w + structural_heat).abs()
                < 1e-10 * power.blocked_acoustic_input_w
        );
        assert!(
            (power.interference_exterior_w + power.blocked_acoustic_input_w).abs()
                < 1e-10 * power.blocked_acoustic_input_w
        );
        assert!(power.blocked_exterior_w.abs() < 1e-10 * power.blocked_acoustic_input_w);
        assert_eq!(power.total_wall_w, 0.0);
        assert!(power.relative_balance < 1e-10);
        // Direct total p/v surface work independently checks that pressure
        // interference is retained when converting field amplitudes to energy.
        let work = result
            .total_surface_pressure_pa
            .iter()
            .zip(&result.total_surface_normal_velocity_m_s)
            .zip(mechanics.surface().areas())
            .map(|((&p, &v), &a)| 0.5 * (p * v.conj()).re * a)
            .sum::<f64>();
        assert!((work - power.total_exterior_w).abs() < 1e-11 * power.blocked_acoustic_input_w);
        let (mut dp, mut pn, mut du, mut un) = (0.0, 0.0, 0.0, 0.0);
        for (i, &(pressure, radial)) in expected.iter().enumerate() {
            let field = result.total_receivers[i];
            dp += (field.pressure_pa - pressure).norm_sq();
            pn += pressure.norm_sq();
            du += (field.particle_velocity_m_s[if i == 0 { 2 } else { 1 }] - radial).norm_sq();
            un += radial.norm_sq();
        }
        errors.push([
            (actual_force - blocked_force).abs() / blocked_force.abs(),
            (actual_velocity - velocity).abs() / velocity.abs(),
            (dp / pn).sqrt(),
            (du / un).sqrt(),
            (structural_heat / heat - 1.0).abs(),
        ]);
        if level == 2 {
            let reaction = result.total_fluid_reaction[0].scale(det::sqrt(MASS));
            assert!(
                (reaction - expected_total_reaction).abs() < 0.08 * expected_total_reaction.abs()
            );
            assert!((power.blocked_acoustic_input_w / sound_work - 1.0).abs() < 0.12);
            assert!((power.motion_exterior_w / self_work - 1.0).abs() < 0.12);
        }
    }
    eprintln!("incident sphere [blocked force, motion, total p/u, heat] errors: {errors:?}");
    for quantity in 0..5 {
        assert!(
            errors[1][quantity] < 0.85 * errors[0][quantity],
            "{errors:?}"
        );
        assert!(
            errors[2][quantity] < 0.85 * errors[1][quantity],
            "{errors:?}"
        );
        assert!(errors[2][quantity] < 0.08, "{errors:?}");
    }
}
