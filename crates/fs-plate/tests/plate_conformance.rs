//! fs-plate CONFORMANCE battery (music bead
//! `frankensim-music-v8-root-3ez8g.13.2`): the reimplementation
//! contract through the public surface with independent oracles. The
//! inline `#[cfg(test)]` unit modules (including the FULL five-mode
//! Olson–Hazell literature case and the MAC-paired convergence
//! machinery) stay untouched — the conformance surface GROWS.
//!
//! Cases:
//! - pt-001: the simply supported isotropic plate lands on the exact
//!   Navier frequency `ω_mn = π²(m²/a² + n²/b²)√(D/ρh)` with mesh
//!   convergence.
//! - pt-002: the clamped square lands on Leissa's `λ = 35.992`
//!   (NASA SP-160).
//! - pt-003: the fs-qty front door is BIT-IDENTICAL to the plain-f64
//!   door (`isotropic_qty` vs `isotropic` through assemble + modes).
//! - pt-004: stiffener term isolation — a vanishing rib is a no-op
//!   and a real rib stiffens the fundamental.
//! - pt-005: Olson–Hazell 1977 fundamental-region pin through the
//!   public API (theory 718.1/751.4 Hz class; the corpus row
//!   `acoustic-olson-hazell-1977-mode3` retains 997.4 Hz; the full
//!   five-mode gate lives inline).
//! - pt-006: refusals by name (Poisson bound, non-physical section,
//!   out-of-range boundary node).

use fs_plate::{
    AssemblyOptions, EdgeSupport, PlateError, PlateMesh, PlateSection, Stiffener, assemble, modes,
};

const PI: f64 = std::f64::consts::PI;

fn steel() -> PlateSection {
    PlateSection::isotropic(200e9, 0.3, 0.002, 7800.0).expect("steel")
}

fn fundamental(
    sec: &PlateSection,
    a: f64,
    b: f64,
    nx: usize,
    ny: usize,
    support: EdgeSupport,
    braces: &[Stiffener],
    window: (f64, f64),
) -> f64 {
    let mesh = PlateMesh::rectangle(a, b, nx, ny);
    let boundary = PlateMesh::rectangle_boundary(nx, ny);
    let model = assemble(
        &mesh,
        sec,
        &boundary,
        braces,
        &AssemblyOptions {
            pretension: 0.0,
            support,
        },
    )
    .expect("assemble");
    let rep = modes(
        &model,
        (window.0 * window.0, window.1 * window.1),
        &fs_modal::SliceOptions::default(),
    )
    .expect("modes");
    assert!(!rep.modes.is_empty(), "no mode in window");
    rep.modes[0].lambda.sqrt()
}

#[test]
fn pt_001_simply_supported_matches_navier() {
    let sec = steel();
    let (a, b) = (0.5, 0.4);
    let d = sec.d[0];
    let rho_h = sec.density * sec.thickness;
    let w11 = PI * PI * (1.0 / (a * a) + 1.0 / (b * b)) * (d / rho_h).sqrt();
    let mut errs = Vec::new();
    for n in [10usize, 20] {
        let w = fundamental(
            &sec,
            a,
            b,
            n,
            n,
            EdgeSupport::SimplySupported,
            &[],
            (0.5 * w11, 1.8 * w11),
        );
        errs.push((w - w11).abs() / w11);
    }
    assert!(errs[1] < 0.01, "fine-mesh Navier error {:.3e}", errs[1]);
    assert!(errs[1] < errs[0], "convergence trend {errs:?}");
    println!(
        "{{\"suite\":\"fs-plate\",\"case\":\"pt-001-navier\",\"analytic_hz\":{:.2},\
         \"errors\":[{:.3e},{:.3e}],\"verdict\":\"pass\"}}",
        w11 / (2.0 * PI),
        errs[0],
        errs[1]
    );
}

#[test]
fn pt_002_clamped_square_matches_leissa() {
    let sec = steel();
    let a = 1.0;
    let d = sec.d[0];
    let rho_h = sec.density * sec.thickness;
    let w1 = 35.992 / (a * a) * (d / rho_h).sqrt();
    let w = fundamental(
        &sec,
        a,
        a,
        20,
        20,
        EdgeSupport::Clamped,
        &[],
        (0.7 * w1, 1.6 * w1),
    );
    let err = (w - w1).abs() / w1;
    assert!(err < 0.02, "Leissa clamped fundamental error {err:.3e}");
    println!(
        "{{\"suite\":\"fs-plate\",\"case\":\"pt-002-leissa\",\"lambda\":35.992,\
         \"error\":{err:.3e},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn pt_003_qty_front_door_is_bit_identical() {
    use fs_qty::{Density, Length, Pressure};
    let plain = PlateSection::isotropic(200e9, 0.3, 0.002, 7800.0).expect("plain");
    let qty = PlateSection::isotropic_qty(
        Pressure::new(200e9),
        0.3,
        Length::new(0.002),
        Density::new(7800.0),
    )
    .expect("qty");
    let (a, b) = (0.5, 0.4);
    let w_plain = fundamental(
        &plain,
        a,
        b,
        12,
        12,
        EdgeSupport::SimplySupported,
        &[],
        (100.0, 800.0),
    );
    let w_qty = fundamental(
        &qty,
        a,
        b,
        12,
        12,
        EdgeSupport::SimplySupported,
        &[],
        (100.0, 800.0),
    );
    assert!(
        w_plain.to_bits() == w_qty.to_bits(),
        "front doors diverged: {w_plain:e} vs {w_qty:e}"
    );
    println!(
        "{{\"suite\":\"fs-plate\",\"case\":\"pt-003-qty-bit-identity\",\"verdict\":\"pass\"}}"
    );
}

#[test]
fn pt_004_stiffener_terms_isolate() {
    let sec = steel();
    let (a, b, n) = (0.5, 0.4, 12usize);
    let rib = |scale: f64| -> Stiffener {
        let (bw, bd) = (6.0e-3 * scale, 10.0e-3 * scale);
        Stiffener {
            nodes: (0..=n).map(|i| (n / 3) * (n + 1) + i).collect(),
            e: 200e9,
            g: 200e9 / 2.6,
            area: bw * bd,
            inertia: bw * bd * bd * bd / 12.0,
            torsion: 0.229 * bd * bw * bw * bw,
            eccentricity: f64::midpoint(0.002, bd),
            density: 7800.0,
        }
    };
    let bare = fundamental(
        &sec,
        a,
        b,
        n,
        n,
        EdgeSupport::SimplySupported,
        &[],
        (100.0, 800.0),
    );
    // A vanishing rib is a no-op within numerical dust.
    let tiny = fundamental(
        &sec,
        a,
        b,
        n,
        n,
        EdgeSupport::SimplySupported,
        &[rib(1.0e-4)],
        (100.0, 800.0),
    );
    assert!(
        (tiny / bare - 1.0).abs() < 1.0e-6,
        "vanishing rib moved the fundamental: {bare:e} vs {tiny:e}"
    );
    // A real off-center rib stiffens it.
    let ribbed = fundamental(
        &sec,
        a,
        b,
        n,
        n,
        EdgeSupport::SimplySupported,
        &[rib(1.0)],
        (100.0, 1600.0),
    );
    assert!(
        ribbed > 1.05 * bare,
        "rib must stiffen: {bare:.1} -> {ribbed:.1} rad/s"
    );
    println!(
        "{{\"suite\":\"fs-plate\",\"case\":\"pt-004-stiffener-isolation\",\
         \"bare_hz\":{:.2},\"ribbed_hz\":{:.2},\"vanishing_rib_rel\":{:.2e},\
         \"verdict\":\"pass\"}}",
        bare / (2.0 * PI),
        ribbed / (2.0 * PI),
        (tiny / bare - 1.0).abs()
    );
}

#[test]
fn pt_005_olson_hazell_fundamental_region_pin() {
    // Public-API pin on the classical clamped rib-stiffened panel:
    // the first two theory rows are 718.1 and 751.4 Hz (Srivastava
    // 2004 Table 4; Thinh 2013 Table 1 — the corpus row
    // acoustic-olson-hazell-1977-mode3 retains the 997.4 Hz third
    // mode). The FULL five-mode/8% gate with MAC pairing lives in
    // the inline unit module; this conformance pin asserts the
    // fundamental lands in the same 8% band through the public
    // surface alone.
    let a = 0.203;
    let h = 1.37e-3;
    let sec = PlateSection::isotropic(68.7e9, 0.3, h, 2820.0).expect("aluminium");
    let (bw, bd) = (6.35e-3, 12.7e-3);
    let n = 16usize;
    let mid = n / 2;
    let rib = Stiffener {
        nodes: (0..=n).map(|i| mid * (n + 1) + i).collect(),
        e: 68.7e9,
        g: 68.7e9 / 2.6,
        area: bw * bd,
        inertia: bw * bd * bd * bd / 12.0,
        torsion: 0.229 * bd * bw * bw * bw,
        eccentricity: f64::midpoint(h, bd),
        density: 2820.0,
    };
    let w = fundamental(
        &sec,
        a,
        a,
        n,
        n,
        EdgeSupport::Clamped,
        &[rib],
        (2.0 * PI * 400.0, 2.0 * PI * 900.0),
    );
    let f = w / (2.0 * PI);
    let dev = (f / 718.1 - 1.0).abs();
    assert!(
        dev < 0.08,
        "Olson-Hazell fundamental {f:.1} Hz vs theory 718.1 Hz ({dev:.3})"
    );
    println!(
        "{{\"suite\":\"fs-plate\",\"case\":\"pt-005-olson-hazell-pin\",\"f_hz\":{f:.1},\
         \"theory_hz\":718.1,\"dev\":{dev:.4},\
         \"corpus\":\"acoustic-olson-hazell-1977-mode3 (997.4 Hz third mode)\",\
         \"verdict\":\"pass\"}}"
    );
}

#[test]
fn pt_006_refusals_fire_by_name() {
    assert!(matches!(
        PlateSection::isotropic(200e9, 0.5, 0.002, 7800.0),
        Err(PlateError::BadSection { .. })
    ));
    assert!(matches!(
        PlateSection::isotropic(200e9, 0.3, -0.002, 7800.0),
        Err(PlateError::BadSection { .. })
    ));
    let sec = steel();
    let mesh = PlateMesh::rectangle(0.5, 0.4, 8, 8);
    let bad_boundary = vec![mesh.node_count() + 5];
    // Found by THIS battery: the out-of-range boundary node used to
    // panic inside assemble; it now refuses by name.
    assert!(matches!(
        assemble(
            &mesh,
            &sec,
            &bad_boundary,
            &[],
            &AssemblyOptions {
                pretension: 0.0,
                support: EdgeSupport::SimplySupported,
            },
        ),
        Err(PlateError::BadBoundary { .. })
    ));
    println!("{{\"suite\":\"fs-plate\",\"case\":\"pt-006-refusals\",\"verdict\":\"pass\"}}");
}
