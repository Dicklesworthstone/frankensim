//! Radiation-limited Q of the authored Gt-S1 soundboard (bead pexfp.283).
//!
//! The renders' body lines ring because the material scenario is nearly all the
//! loss they have. fs-bem's complete modal loading (Z_rs = integral phi_r p_s dS,
//! planar baffled Rayleigh route) supplies the missing radiation resistance. This
//! test only MEASURES it on the declared plate; it does not wire it into a render.
//! Disclosed limits: an infinite rigid baffle, the top plate alone (no back, ribs,
//! braces, sound hole or joints), and the Navier modes of a specially orthotropic
//! simply supported Kirchhoff plate, grain along the length.
//!
//! Prediction and tolerances were written before the first fs-bem call (Agent Mail,
//! thread irp-instruments-swarm, sha256 prefix 47df76168f2c9642).

use fs_bem::helmholtz::Medium;
use fs_bem::helmholtz::rayleigh_modal::AreaSample;
use fs_bem::modal_loading::{LoadingRoute, SurfaceModeTraces, modal_radiation_loading};
use std::f64::consts::PI;
use std::time::Instant;

const CEDAR: &str = include_str!("../../../data/instruments/descriptions/guitar/cedar-top.tsv");
const PERFORMANCE: &str =
    include_str!("../../../data/instruments/descriptions/guitar/gt-s1-solo-g3.material-strings");
const MODES: usize = 14;
const RUNGS: [(usize, usize); 3] = [(24, 18), (48, 36), (96, 72)];
/// Carcagno et al. 2018, JASA 144(6) 3533, Table I, BR-rosewood guitar:
/// (frequency Hz, Q). Ledger row carcagno-2018-guitar.
const MEASURED: [(f64, f64); 3] = [(97.0, 34.0), (177.0, 18.0), (336.0, 36.0)];

fn scalar(property: &str) -> f64 {
    CEDAR
        .lines()
        .find_map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            (fields.first() == Some(&"scalar") && fields.get(3) == Some(&property))
                .then(|| fields[4].parse::<f64>().expect("source scalar"))
        })
        .unwrap_or_else(|| panic!("required source property missing: {property}"))
}

struct Plate {
    a: f64,
    b: f64,
    h: f64,
    rho: f64,
    d11: f64,
    d22: f64,
    d12: f64,
    d66: f64,
    /// Kelvin-Voigt loss time k: tan(delta) = omega k, from the pack's viscosity / stiffness.
    loss_time_s: f64,
}

fn plate() -> Plate {
    let body: Vec<&str> = PERFORMANCE
        .lines()
        .find(|line| line.starts_with("body "))
        .expect("body record")
        .split_whitespace()
        .collect();
    let (a, b, h): (f64, f64, f64) = (
        body[2].parse().unwrap(),
        body[3].parse().unwrap(),
        body[4].parse().unwrap(),
    );
    let (e1, e2) = (scalar("young_modulus_1"), scalar("young_modulus_2"));
    let (nu12, g12) = (scalar("poisson_ratio_12"), scalar("shear_modulus_12"));
    let nu21 = nu12 * e2 / e1;
    let den = 1.0 - nu12 * nu21;
    let c11 = e1 / den;
    let stiffness_to_d = h * h * h / 12.0;
    Plate {
        a,
        b,
        h,
        rho: scalar("density"),
        d11: c11 * stiffness_to_d,
        d22: e2 / den * stiffness_to_d,
        d12: nu12 * e2 / den * stiffness_to_d,
        d66: g12 * stiffness_to_d,
        loss_time_s: scalar("plane_stress_viscosity_11") / c11,
    }
}

/// (frequency Hz, m, n) of the lowest `MODES` Navier modes.
fn lowest_modes(p: &Plate) -> Vec<(f64, usize, usize)> {
    let mut all = Vec::new();
    for m in 1..=8usize {
        for n in 1..=8usize {
            let (x, y) = (m as f64 / p.a, n as f64 / p.b);
            let w2 = PI * PI * PI * PI / (p.rho * p.h)
                * (p.d11 * x * x * x * x
                    + 2.0 * (p.d12 + 2.0 * p.d66) * x * x * y * y
                    + p.d22 * y * y * y * y);
            all.push((w2.sqrt() / (2.0 * PI), m, n));
        }
    }
    all.sort_by(|u, v| u.0.total_cmp(&v.0));
    all.truncate(MODES);
    all
}

struct Row {
    omega: f64,
    resistance: f64,
    reactance: f64,
}

/// Diagonal modal radiation impedance at each mode's own frequency (kg/s per unit
/// modal velocity of the peak-normalized trace), on an nx-by-ny cell grid.
fn bake(p: &Plate, modes: &[(f64, usize, usize)], nx: usize, ny: usize) -> Vec<Row> {
    let (dx, dy) = (p.a / nx as f64, p.b / ny as f64);
    let mut cells = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            let (x, y) = ((i as f64 + 0.5) * dx, (j as f64 + 0.5) * dy);
            cells.push(AreaSample {
                position_m: [x, y],
                area_m2: dx * dy,
                mode_shape: modes
                    .iter()
                    .map(|&(_, m, n)| {
                        (m as f64 * PI * x / p.a).sin() * (n as f64 * PI * y / p.b).sin()
                    })
                    .collect(),
            });
        }
    }
    let omegas: Vec<f64> = modes.iter().map(|m| 2.0 * PI * m.0).collect();
    assert!(omegas.windows(2).all(|w| w[0] < w[1]));
    let provenance = format!(
        "cedar plate {} x {} x {} m; {nx}x{ny} cells; {MODES} Navier modes",
        p.a, p.b, p.h
    );
    let bake = modal_radiation_loading(
        SurfaceModeTraces::PlanarBaffled { samples: &cells },
        LoadingRoute::PlanarBaffledRayleigh,
        &omegas,
        Medium::air(),
        fs_blake3::hash_domain("gt-s1-plate-radiation", provenance.as_bytes()),
        MODES * MODES * omegas.len(),
        || false,
    )
    .expect("planar baffled modal loading");
    omegas
        .iter()
        .enumerate()
        .map(|(r, &omega)| {
            let z = bake.at_omega(omega).expect("sampled frequency").impedance()[r * MODES + r];
            Row {
                omega,
                resistance: z.re,
                reactance: z.im,
            }
        })
        .collect()
}

#[test]
fn authored_soundboard_radiation_q_meets_the_prediction_and_is_validated_against_carcagno() {
    let started = Instant::now();
    let p = plate();
    let modes = lowest_modes(&p);
    // Peak-normalized sin*sin trace: m_r = rho h a b / 4.
    let mass = p.rho * p.h * p.a * p.b / 4.0;
    let rungs: Vec<Vec<Row>> = RUNGS
        .iter()
        .map(|&(nx, ny)| bake(&p, &modes, nx, ny))
        .collect();
    let (coarse, finest) = (&rungs[1], &rungs[2]);
    let index = |m: usize, n: usize| {
        modes
            .iter()
            .position(|x| x.1 == m && x.2 == n)
            .expect("mode in the lowest set")
    };
    let q_rad = |r: usize| finest[r].omega * mass / finest[r].resistance;
    let q_mat = |r: usize| 1.0 / (finest[r].omega * p.loss_time_s);
    for (r, &(f, m, n)) in modes.iter().enumerate() {
        let change = (finest[r].resistance - coarse[r].resistance).abs() / finest[r].resistance;
        let added_mass = -finest[r].reactance / finest[r].omega;
        println!(
            "{{\"suite\":\"gt-s1-plate-radiation\",\"mode\":[{m},{n}],\"f_hz\":{f},\"resistance_kg_s\":{},\"q_radiation\":{},\"q_material\":{},\"q_total\":{},\"added_mass_kg\":{added_mass},\"added_mass_over_modal_mass\":{},\"first_order_frequency_shift_fraction\":{},\"resistance_change_last_rung\":{change},\"cells\":{:?}}}",
            finest[r].resistance,
            q_rad(r),
            q_mat(r),
            1.0 / (1.0 / q_rad(r) + 1.0 / q_mat(r)),
            added_mass / mass,
            -added_mass / (2.0 * mass),
            RUNGS[2],
        );
    }
    // Pre-written prediction bands (Agent Mail record).
    for (m, n, want, tol) in [
        (1, 1, 68.0, 0.25),
        (1, 3, 190.0, 0.35),
        (3, 1, 87.0, 0.35),
        (1, 5, 218.0, 0.35),
        (3, 3, 610.0, 0.35),
    ] {
        let q = q_rad(index(m, n));
        println!(
            "{{\"suite\":\"gt-s1-plate-radiation\",\"prediction\":[{m},{n}],\"predicted_q_radiation\":{want},\"band\":{tol},\"measured\":{q}}}"
        );
        assert!(
            (q / want - 1.0).abs() <= tol,
            "mode ({m},{n}): Q_rad {q} outside {want} +-{}%",
            tol * 100.0
        );
    }
    let reference = q_rad(index(1, 1));
    for (r, &(_, m, n)) in modes.iter().enumerate() {
        assert!(
            finest[r].resistance >= 0.0,
            "radiation resistance must be passive"
        );
        if m % 2 == 0 || n % 2 == 0 {
            assert!(
                q_rad(r) > 10.0 * reference,
                "even-index mode ({m},{n}) should be a weak radiator: {}",
                q_rad(r)
            );
        }
    }
    let r11 = index(1, 1);
    let converged =
        (finest[r11].resistance - coarse[r11].resistance).abs() / finest[r11].resistance;
    assert!(converged < 0.05, "R(1,1) cell ladder change {converged}");
    // Validation against the measured guitar: how much of the measured damping rate do material
    // loss plus top-plate radiation explain at the nearest plate mode? The pre-written gate asks
    // for at least 70 % at the nearest mode of EACH row; the prediction was that it FAILS.
    let mut all_pass = true;
    for &(f_meas, q_meas) in &MEASURED {
        let r = (0..modes.len())
            .min_by(|&u, &v| {
                (modes[u].0 - f_meas)
                    .abs()
                    .total_cmp(&(modes[v].0 - f_meas).abs())
            })
            .unwrap();
        let explained = (1.0 / q_mat(r) + 1.0 / q_rad(r)) / (1.0 / q_meas);
        all_pass &= explained >= 0.7;
        println!(
            "{{\"suite\":\"gt-s1-plate-radiation\",\"carcagno_row_hz\":{f_meas},\"q_measured\":{q_meas},\"nearest_mode\":[{},{}],\"mode_f_hz\":{},\"q_material\":{},\"q_radiation\":{},\"fraction_of_measured_damping_rate_explained\":{explained}}}",
            modes[r].1,
            modes[r].2,
            modes[r].0,
            q_mat(r),
            q_rad(r)
        );
    }
    println!(
        "{{\"suite\":\"gt-s1-plate-radiation\",\"gate_70_percent_at_each_row\":\"{}\",\"baffle\":\"infinite rigid\",\"elapsed_s\":{}}}",
        if all_pass { "pass" } else { "fail" },
        started.elapsed().as_secs_f64()
    );
}
