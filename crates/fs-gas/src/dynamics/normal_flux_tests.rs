use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * (1.0 + expected.abs()),
        "actual={actual:e}, expected={expected:e}",
    );
}

#[test]
fn g0_planar_flux_retains_original_one_dimensional_results() {
    let g = gas(0.0, 1e16);
    for (left, right) in [
        ([1.0, 0.0, 1.0], [4.0, 0.0, 1.0]),
        ([1.0, 0.0, 1.0], [0.125, 0.0, 0.1]),
        ([1.0, -2.0, 0.4], [1.0, 2.0, 0.4]),
        ([1.2, 10.0, 1.7], [0.8, 8.0, 1.3]),
    ] {
        let l = g.conserved(left[0], left[1], left[2]).unwrap();
        let r = g.conserved(right[0], right[1], right[2]).unwrap();
        let axial = g.hllc_flux(l, r, || false).unwrap();
        let planar = g.hllc_normal_flux(l.into(), r.into(), || false).unwrap();
        assert_eq!(FaceFlux::from(planar), axial);
        assert_eq!(planar.transverse_momentum, 0.0);
    }
}

#[test]
fn g0_planar_shear_contact_advects_full_energy_without_pressure_work_duplication() {
    let g = gas(0.4, -1e16);
    for u in [-8.0, -0.2, 0.0, 0.7, 10.0] {
        let l = g.normal_conserved(1.0, u, 3.0, 1.0).unwrap();
        let r = g.normal_conserved(0.3, u, -4.0, 1.0).unwrap();
        let f = g.hllc_normal_flux(l, r, || false).unwrap();
        let upwind = if u >= 0.0 { l } else { r };
        close(f.mass, upwind.normal_momentum);
        close(f.normal_momentum, upwind.normal_momentum * u + 1.0);
        close(f.transverse_momentum, upwind.transverse_momentum * u);
        close(f.energy, (upwind.energy + 1.0) * u);
        assert!(!f.used_hll_fallback);
        let p = g.normal_primitive(upwind, || false).unwrap();
        close(p.velocity, u);
        close(p.pressure, 1.0);
    }
}

#[test]
fn g3_planar_hllc_and_hll_are_covariant_under_transverse_boosts() {
    let g = gas(0.0, 0.0);
    for (left, right, fallback) in [
        ([0.9, -0.1, 1.5, 2.0], [0.35, 0.2, -2.0, 0.2], false),
        ([1.0, -2.0, 0.5, 0.4], [1.0, 2.0, -0.7, 0.4], true),
    ] {
        let l = g
            .normal_conserved(left[0], left[1], left[2], left[3])
            .unwrap();
        let r = g
            .normal_conserved(right[0], right[1], right[2], right[3])
            .unwrap();
        let f = g.hllc_normal_flux(l, r, || false).unwrap();
        assert_eq!(f.used_hll_fallback, fallback);
        for boost in [-3.0, 2.0] {
            let shift = |q: NormalConserved| NormalConserved {
                transverse_momentum: q.transverse_momentum + boost * q.density,
                energy: q.energy + boost * q.transverse_momentum + 0.5 * boost * boost * q.density,
                ..q
            };
            let shifted = g.hllc_normal_flux(shift(l), shift(r), || false).unwrap();
            assert_eq!(shifted.used_hll_fallback, fallback);
            close(shifted.signal_speed, f.signal_speed);
            close(shifted.mass, f.mass);
            close(shifted.normal_momentum, f.normal_momentum);
            close(
                shifted.transverse_momentum,
                f.transverse_momentum + boost * f.mass,
            );
            close(
                shifted.energy,
                f.energy + boost * f.transverse_momentum + 0.5 * boost * boost * f.mass,
            );
        }
    }
}

#[test]
fn g0_planar_hll_accounts_for_shear_heating_in_source_domain_admission() {
    let g = gas(0.0, 0.0);
    // Both traces are cold and admissible. The near-vacuum expansion forces
    // HLL, whose mixed state converts opposing transverse KE to heat. Its
    // temperature exceeds this source's upper limit, even as waves expand.
    let left = g.normal_conserved(1.0, -2.0, 30.0, 0.4).unwrap();
    let right = g.normal_conserved(1.0, 2.0, -30.0, 0.4).unwrap();
    assert_eq!(
        g.hllc_normal_flux(left, right, || false),
        Err(Error::Invalid(
            "HLL intermediate state outside admitted source domain"
        )),
    );
    // Discarding transverse energy would incorrectly admit that same face.
    assert!(
        g.hllc_flux(
            g.conserved(1.0, -2.0, 0.4).unwrap(),
            g.conserved(1.0, 2.0, 0.4).unwrap(),
            || false,
        )
        .unwrap()
        .used_hll_fallback
    );
}

#[test]
fn g0_planar_slip_wall_retains_caloric_shock_pressure_and_zero_work() {
    let g = gas(0.4, 1e18);
    for normal in [-1.0, 1.0] {
        for velocity in [-0.4, 0.0, 2.0] {
            let axial = g
                .normal_conserved(1.2, normal * velocity, 0.0, 1.0)
                .unwrap();
            let shear = g
                .normal_conserved(1.2, normal * velocity, 5.0, 1.0)
                .unwrap();
            let reference = g
                .stationary_wall_normal_flux(axial, normal, || false)
                .unwrap();
            let flux = g
                .stationary_wall_normal_flux(shear, normal, || false)
                .unwrap();
            close(flux.normal_momentum, reference.normal_momentum);
            close(flux.signal_speed, reference.signal_speed);
            assert_eq!(flux.mass, 0.0);
            assert_eq!(flux.transverse_momentum, 0.0);
            assert_eq!(flux.energy, 0.0);
        }
    }
}

#[test]
fn g4_planar_flux_refuses_erased_heat_invalid_states_and_cancellation() {
    let g = gas(0.4, 0.0);
    assert!(g.normal_conserved(1.0, 0.0, 1e12, 1.0).is_err());
    assert!(g.normal_conserved(1.0, 0.0, f64::NAN, 1.0).is_err());
    let q = g.normal_conserved(1.0, 0.2, 0.3, 1.0).unwrap();
    assert!(
        g.normal_primitive(NormalConserved { density: 0.0, ..q }, || false)
            .is_err()
    );
    assert!(
        g.normal_primitive(NormalConserved { energy: -1.0, ..q }, || false)
            .is_err()
    );
    assert!(g.stationary_wall_normal_flux(q, 0.0, || false).is_err());
    assert_eq!(g.normal_primitive(q, || true), Err(Error::Cancelled));
    assert_eq!(g.hllc_normal_flux(q, q, || true), Err(Error::Cancelled));
    assert_eq!(
        g.stationary_wall_normal_flux(q, 1.0, || true),
        Err(Error::Cancelled)
    );
    let mut polls = 0;
    assert_eq!(
        g.hllc_normal_flux(q, q, || {
            polls += 1;
            polls == 4
        }),
        Err(Error::Cancelled)
    );
}
