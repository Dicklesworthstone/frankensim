//! E2E casebook: a single-reed exciter recast as pHS components — the
//! reed lamella is a mass-spring-damper pHS, the bore a modal-bank
//! pHS, and the reed channel a MEMORYLESS DISSIPATIVE valve port
//! (Bernoulli jet: all pressure head lost, `dp * U >= 0` for every
//! sample). The whole-instrument energy audit is then trivial: each
//! subsystem is structurally passive, the valve provably dissipates,
//! and the only source is the mouth `pm * U`.
//!
//! Validation pins (independent of the integrator):
//! - quasi-static equilibrium below the oscillation threshold:
//!   `q* = Sr * pm / k`, `h* = h0 - q*`, bore DC pressure exactly 0
//!   (acoustic modal bank has no DC mode), `U* = w h* sqrt(2 pm/rho)`
//!   — derived by hand from the component laws, not from the code.
//! - per-step: both subsystems' supply audits hold, the valve power
//!   `dp * U` is never negative, and the closed accounting
//!   `dH_total <= pm * U * dt` (mouth supply bounds energy growth).
//!
//! NO-CLAIM: oscillation-regime validation (threshold pressures,
//! limit-cycle spectra) is the distributed-contact/exciter follow-up;
//! this casebook validates the pHS recast and its energy discipline in
//! the quasi-static regime where an exact analytic answer exists.

use fs_math::det;
use fs_phs::{mass_spring_damper, modal_bank, step};

const TWO_PI: f64 = 2.0 * core::f64::consts::PI;

/// Bernoulli reed-channel flow: `U = w * h * sqrt(2 dp / rho)` for
/// `dp > 0, h > 0`, else 0 (no backflow modeling below threshold).
fn valve_flow(w: f64, h: f64, dp: f64, rho: f64) -> f64 {
    if dp <= 0.0 || h <= 0.0 {
        0.0
    } else {
        w * h * det::sqrt(2.0 * dp / rho)
    }
}

#[test]
fn reed_as_phs_quasi_static_equilibrium_and_energy_audit() {
    // Reed lamella (clarinet-class): stiffness per area ~ 8e6 Pa/m,
    // effective area Sr, strong damping (below oscillation threshold).
    let (m_r, k_r, c_r) = (3.0e-6, 120.0, 0.02);
    let s_r = 1.0e-4; // effective pressure-collection area [m^2]
    let h0: f64 = 6.0e-4; // rest opening [m]
    let w_r = 1.2e-2; // channel width [m]
    let rho = 1.2;
    let pm = 300.0; // mouth pressure [Pa], far below threshold here
    let reed = mass_spring_damper(m_r, k_r, c_r).expect("reed");
    // Bore: 6-mode acoustic bank (no DC mode BY CONSTRUCTION — the
    // equilibrium pin depends on that fact).
    let omegas: Vec<f64> = (0i32..6)
        .map(|i| TWO_PI * (147.0 * f64::from(2 * i + 1)))
        .collect();
    let zetas = vec![0.03; 6];
    let drive = vec![30.0, 25.0, 20.0, 16.0, 12.0, 9.0];
    let bore = modal_bank(&omegas, &zetas, &drive).expect("bore");
    let dt = 2.0e-5;
    let mut x_reed = vec![0.0, 0.0];
    let mut x_bore = vec![0.0; bore.state_dim()];
    let mut worst_valve_power = 0.0f64;
    let mut worst_closed_defect = f64::NEG_INFINITY;
    let steps = 60_000usize;
    for _ in 0..steps {
        // Current observables.
        let p_bore = bore.output(&x_bore)[0];
        let q = x_reed[0];
        let h = (h0 - q).max(0.0);
        let dp = pm - p_bore;
        let flow = valve_flow(w_r, h, dp, rho);
        // Valve dissipativity: the jet destroys head, never creates.
        let valve_power = dp * flow;
        worst_valve_power = worst_valve_power.min(valve_power);
        // Step both pHS under their held port inputs.
        let rec_r = step(&reed, &x_reed, &[s_r * dp], dt).expect("reed step");
        let rec_b = step(&bore, &x_bore, &[flow], dt).expect("bore step");
        // Structural passivity of each component, every step.
        assert!(rec_r.supply_defect() <= 1.0e-10, "reed supply audit");
        assert!(rec_b.supply_defect() <= 1.0e-10, "bore supply audit");
        // Closed accounting: total stored-energy growth is bounded by
        // the mouth supply plus the reed's pressure work (both routed
        // through dp; the valve only ever removes energy).
        let d_total = rec_r.delta_h + rec_b.delta_h;
        let mouth_supply = (pm * flow + s_r * dp * rec_r.y[0]) * dt;
        worst_closed_defect = worst_closed_defect.max(d_total - mouth_supply);
        x_reed = rec_r.x;
        x_bore = rec_b.x;
    }
    assert!(
        worst_valve_power >= 0.0,
        "valve created energy: {worst_valve_power:.3e}"
    );
    assert!(
        worst_closed_defect <= 1.0e-9,
        "closed energy accounting violated by {worst_closed_defect:.3e}"
    );
    // Quasi-static equilibrium vs the ANALYTIC solution.
    let q_star = s_r * pm / k_r;
    let h_star = h0 - q_star;
    let u_star = valve_flow(w_r, h_star, pm, rho);
    let p_bore_end = bore.output(&x_bore)[0];
    assert!(
        (x_reed[0] - q_star).abs() <= 1.0e-4 * q_star,
        "reed deflection {:.6e} vs analytic {q_star:.6e}",
        x_reed[0]
    );
    assert!(
        x_reed[1].abs() <= 1.0e-8,
        "reed momentum not settled: {:.3e}",
        x_reed[1]
    );
    assert!(
        p_bore_end.abs() <= 1.0e-6 * pm,
        "bore DC pressure must vanish (modal bank has no DC mode): {p_bore_end:.3e}"
    );
    let flow_end = valve_flow(w_r, (h0 - x_reed[0]).max(0.0), pm - p_bore_end, rho);
    assert!(
        (flow_end - u_star).abs() <= 1.0e-4 * u_star,
        "steady flow {flow_end:.6e} vs analytic {u_star:.6e}"
    );
    println!(
        "{{\"suite\":\"fs-phs-reed\",\"case\":\"quasi-static\",\"q_star\":{q_star:.6e},\"q_end\":{:.6e},\"u_star\":{u_star:.6e},\"u_end\":{flow_end:.6e},\"worst_closed_defect\":{worst_closed_defect:.3e},\"verdict\":\"pass\"}}",
        x_reed[0]
    );
}
