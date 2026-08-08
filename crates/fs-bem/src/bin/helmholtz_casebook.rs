//! End-to-end Helmholtz radiation casebook
//! (bead frankensim-fsim-helmholtz-bem-k1ryv): mesh -> Burton-Miller
//! solve -> radiation impedance matrix -> spherical-harmonic
//! directivity table -> radiation efficiency, emitting one JSON line
//! per stage.
//!
//! Output contract: every field is deterministic (result values are
//! hashed with the canonical FNV-1a-64 helper over little-endian f64
//! bytes) EXCEPT `elapsed_ms` and `estimated_dense_bytes`'s companion
//! `cap_headroom_panels`, which are run evidence: `elapsed_ms` is
//! wall-clock and varies per host, `estimated_dense_bytes` is a
//! deterministic model (3 dense complex matrices, 16 bytes per entry),
//! not a measured RSS. Rerunning on any host must reproduce every hash
//! bit-identically; a hash drift is a determinism regression, not
//! noise.

use std::time::Instant;

use fs_bem::helmholtz::{
    Formulation, MAX_DENSE_PANELS, Medium, RadiationSolution, directivity_sh_table,
    radiation_efficiency, radiation_impedance_matrix, solve_radiation,
};
use fs_bem::panel3d::SpherePanels;
use fs_casebook::fnv1a64;
use fs_math::c64::C64;

fn hash_c64s(values: &[C64]) -> u64 {
    let mut bytes = Vec::with_capacity(values.len() * 16);
    for v in values {
        bytes.extend_from_slice(&v.re.to_le_bytes());
        bytes.extend_from_slice(&v.im.to_le_bytes());
    }
    fnv1a64(&bytes)
}

fn solve_stage(
    surface: &SpherePanels,
    ka: f64,
    velocity: &[C64],
    label: &str,
) -> RadiationSolution {
    let n = surface.centroids().len();
    let start = Instant::now();
    let sol = solve_radiation(
        surface,
        ka,
        Medium::air(),
        velocity,
        Formulation::BurtonMiller,
    )
    .unwrap_or_else(|e| {
        // A refusal in the fixed casebook configuration is itself a
        // reportable outcome: emit the structured row and stop.
        println!(
            "{{\"suite\":\"fs-bem-helmholtz-casebook\",\"case\":\"solve\",\"mode\":\"{label}\",\
             \"ka\":{ka},\"verdict\":\"refused\",\"error\":\"{e}\"}}"
        );
        std::process::exit(1)
    });
    let elapsed_ms = start.elapsed().as_secs_f64() * 1e3;
    println!(
        "{{\"suite\":\"fs-bem-helmholtz-casebook\",\"case\":\"solve\",\"mode\":\"{label}\",\
         \"ka\":{ka},\"panels\":{n},\"ppw\":{:.2},\"condition_lower_bound\":{:.3e},\
         \"cap_utilization\":{:.4},\"cap_headroom_panels\":{},\"estimated_dense_bytes\":{},\
         \"radiated_power\":{:.6e},\"pressure_fnv\":\"{:016x}\",\"elapsed_ms\":{elapsed_ms:.1}}}",
        sol.panels_per_wavelength,
        sol.condition_lower_bound,
        sol.dense_cap_utilization,
        MAX_DENSE_PANELS - n,
        3 * 16 * n * n,
        sol.radiated_power,
        hash_c64s(&sol.pressure),
    );
    sol
}

fn main() {
    // Stage 1: mesh. Deterministic icosphere refinements.
    let coarse = SpherePanels::icosphere(1.0, 1).expect("icosphere-1");
    let surface = SpherePanels::icosphere(1.0, 2).expect("icosphere-2");
    println!(
        "{{\"suite\":\"fs-bem-helmholtz-casebook\",\"case\":\"mesh\",\"panels_coarse\":{},\
         \"panels\":{},\"dense_cap\":{MAX_DENSE_PANELS}}}",
        coarse.centroids().len(),
        surface.centroids().len(),
    );

    let n = surface.centroids().len();
    let uniform = vec![C64::ONE; n];
    let dipole: Vec<C64> = surface
        .normals()
        .iter()
        .map(|nrm| C64::from_re(nrm[2]))
        .collect();

    for &ka in &[0.5, 1.0, 2.0, core::f64::consts::PI] {
        let sol = solve_stage(&surface, ka, &uniform, "monopole");
        // Stage 3: spherical-harmonic directivity table.
        let start = Instant::now();
        let table =
            directivity_sh_table(&surface, &sol, Medium::air(), 8).expect("directivity table");
        let elapsed_ms = start.elapsed().as_secs_f64() * 1e3;
        println!(
            "{{\"suite\":\"fs-bem-helmholtz-casebook\",\"case\":\"sh-table\",\"mode\":\"monopole\",\
             \"ka\":{ka},\"l_max\":{},\"captured_fraction\":{:.6},\"coefficients_fnv\":\"{:016x}\",\
             \"elapsed_ms\":{elapsed_ms:.1}}}",
            table.l_max,
            table.captured_fraction,
            hash_c64s(&table.coefficients),
        );
        // Stage 4: radiation efficiency vs the analytic monopole oracle.
        let sigma = radiation_efficiency(&surface, &sol, Medium::air()).expect("efficiency");
        let oracle = ka * ka / (1.0 + ka * ka);
        println!(
            "{{\"suite\":\"fs-bem-helmholtz-casebook\",\"case\":\"efficiency\",\
             \"mode\":\"monopole\",\"ka\":{ka},\"sigma\":{sigma:.6},\"oracle\":{oracle:.6}}}"
        );
    }
    // One dipole arm at ka = 1 (the mode-shape route the vibroacoustic
    // bead consumes).
    let sol = solve_stage(&surface, 1.0, &dipole, "dipole");
    let table = directivity_sh_table(&surface, &sol, Medium::air(), 8).expect("directivity table");
    let spec = table.power_by_degree();
    let total: f64 = spec.iter().sum();
    println!(
        "{{\"suite\":\"fs-bem-helmholtz-casebook\",\"case\":\"sh-table\",\"mode\":\"dipole\",\
         \"ka\":1.0,\"l1_fraction\":{:.6},\"coefficients_fnv\":\"{:016x}\"}}",
        spec[1] / total,
        hash_c64s(&table.coefficients),
    );

    // Stage 2 (coarse mesh: the full n x n matrix is the deliverable):
    // radiation impedance matrix with the reciprocity residual.
    let nc = coarse.centroids().len();
    let start = Instant::now();
    let z = radiation_impedance_matrix(&coarse, 1.0, Medium::air(), Formulation::BurtonMiller)
        .expect("impedance matrix");
    let elapsed_ms = start.elapsed().as_secs_f64() * 1e3;
    let mut asym = 0.0f64;
    let mut scale = 0.0f64;
    for i in 0..nc {
        for j in 0..nc {
            let zij = z[i * nc + j].scale(coarse.areas()[i]);
            let zji = z[j * nc + i].scale(coarse.areas()[j]);
            asym = asym.max((zij - zji).abs());
            scale = scale.max(zij.abs());
        }
    }
    println!(
        "{{\"suite\":\"fs-bem-helmholtz-casebook\",\"case\":\"impedance-matrix\",\"ka\":1.0,\
         \"panels\":{nc},\"reciprocity_rel\":{:.4},\"z_fnv\":\"{:016x}\",\
         \"elapsed_ms\":{elapsed_ms:.1}}}",
        asym / scale,
        hash_c64s(&z),
    );
    println!(
        "{{\"suite\":\"fs-bem-helmholtz-casebook\",\"case\":\"summary\",\"verdict\":\"complete\"}}"
    );
}
