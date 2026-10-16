//! Jet-labium collision-operator probe driver (bead 3zkcr): runs
//! the 9ok02 turbulent-probe rig family under a chosen
//! CollisionModel2 and prints per-run verdicts — the exploratory
//! sibling of the pinned `turbulent_regime_probe` test. Edit the
//! mode list to probe an operator; the pinned table lives in
//! tests/edgetone_staging.rs.

use fs_aeroac::jetlab::{JetLabiumConfig, run_jet_labium, transverse_force_peak};
use fs_lbm::core2::CollisionModel2;

fn main() {
    for mode in [CollisionModel2::CentralMoment] {
        for (re, tau) in [(2304.0, 0.5025), (9216.0, 0.500_625)] {
            let cfg = JetLabiumConfig {
                nx: 768,
                ny: 256,
                slot_half: 12.0,
                slot_smoothing: 4.8,
                u_jet: 0.08,
                tau,
                edge_distance: 240,
                plate_length: 200,
                fringe_width: 128,
                fringe_sigma: 0.3,
                steps_settle: 12_000,
                steps_record: 16_384,
                seed_amplitude: 0.005,
                nozzle_thickness: 4,
                collision: mode,
            };
            match run_jet_labium(&cfg) {
                Ok(run) => {
                    let d = &run.diagnostics;
                    let imb =
                        (d.flux_plate_plane - d.flux_fringe_plane).abs() / d.flux_plate_plane.abs();
                    let p = transverse_force_peak(&run.force_series, 12.0, 0.08, 6).unwrap();
                    println!(
                        "mode={mode:?} re={re} RAN st={:.5} bin={} prom={:.2e} rms={:.3e} imb={imb:.5} mach={:.3}",
                        p.strouhal, p.bin, p.prominence, p.force_rms, d.mach_max_lattice
                    );
                }
                Err(e) => println!("mode={mode:?} re={re} REFUSED: {e}"),
            }
        }
    }
}
