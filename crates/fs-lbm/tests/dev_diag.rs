//! D3Q19 deviation-map diagnostic (bead 84hv): prints the per-cell
//! relative deviation of the simulated duct section from the analytic
//! rectangular series — the tool that localized the halfway-bounce-back
//! wall-slip defect and proved it is the BGK τ-dependence, not
//! resolution: at τ = 0.8 the corner cells read −11.4% (12×12) and
//! −8.7% (32×32) — resolution-independent; at the magic
//! τ = ½ + √3/4 (as configured below, the (τ−½)² = 3/16 slip
//! cancellation) the corner drops to −3.7% at 12×12 and the wall rows
//! to −0.7%. Run with
//! `cargo test -p fs-lbm --release --test dev_diag -- --ignored --nocapture`.

use fs_lbm::{Duct, duct_analytic};
use std::fmt::Write as _;

#[test]
#[ignore = "diagnostic: prints the duct deviation map, asserts finiteness only"]
fn duct_deviation_map() {
    let (nx, ny, nz, tau, gz) = (12, 12, 8, 0.9330127018922193, 1e-6);
    let mut duct = Duct::new(nx, ny, nz, tau, gz);
    duct.run(6000);
    let nu = duct.viscosity();
    let section = duct.z_velocity_section();
    for y in 0..ny {
        let mut row = String::new();
        for x in 0..nx {
            let sim = section[y * nx + x];
            let ana = duct_analytic(gz, nu, nx, ny, x, y);
            assert!(sim.is_finite() && ana.is_finite());
            let _ = write!(row, "{:6.3} ", (sim - ana) / ana);
        }
        println!("{row}");
    }
    let c = section[(ny / 2) * nx + nx / 2];
    let ca = duct_analytic(gz, nu, nx, ny, nx / 2, ny / 2);
    println!("center sim {c:.6e} ana {ca:.6e} rel {:+.4}", (c - ca) / ca);
}
