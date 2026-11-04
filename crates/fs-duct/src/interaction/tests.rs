use super::*;
use fs_material::gas::GasSpec;

const TAU: f64 = core::f64::consts::TAU;

fn air() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air")
}

fn close(a: C64, b: C64, tolerance: f64) -> bool {
    (a - b).abs() <= tolerance * b.abs().max(a.abs())
}

fn interaction(epsilon: f64) -> OpenHoleInteraction {
    OpenHoleInteraction::new(epsilon).unwrap()
}

#[test]
fn the_radiation_space_factor_is_positive_and_finite() {
    assert_eq!(interaction(2.0).epsilon(), 2.0);
    assert_eq!(interaction(4.0).epsilon(), 4.0);
    for bad in [0.0, -2.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(OpenHoleInteraction::new(bad).is_err(), "{bad}");
    }
}

#[test]
fn the_mutual_impedance_is_the_conjugate_of_the_printed_monopole_form() {
    // Lefebvre, Scavone and Kergomard 2013, Eq. 16, for exp(+j w t):
    //   Z = j k rho c exp(-j k d) / (2 pi d) = s (sin kd + j cos kd),  s = k rho c / (2 pi d).
    // This crate's exp(-i w t) takes the conjugate: s (sin kd - i cos kd).
    let gas = air();
    let (omega, d) = (TAU * 1500.0, 0.04);
    let k = omega / gas.sound_speed;
    let s = k * gas.density * gas.sound_speed / (2.0 * core::f64::consts::PI * d);
    let expected = C64::new(s * (k * d).sin(), -s * (k * d).cos());
    let z = mutual_radiation_impedance(&gas, omega, d, interaction(2.0)).unwrap();
    assert!(close(z, expected, 1.0e-12), "{z:?} vs {expected:?}");
    // The full-space factor halves it.
    let full = mutual_radiation_impedance(&gas, omega, d, interaction(4.0)).unwrap();
    assert!(close(full, expected.scale(0.5), 1.0e-12));
    // Reciprocal in the two stations, positive radiation resistance below kd = pi.
    assert!(z.re > 0.0 && k * d < core::f64::consts::PI);
    // As kd -> 0 the real part is the half-space point-source radiation resistance rho w^2 / (2 pi c).
    let near = mutual_radiation_impedance(&gas, omega, 1.0e-6, interaction(2.0)).unwrap();
    let source = gas.density * omega * omega / (2.0 * core::f64::consts::PI * gas.sound_speed);
    assert!((near.re - source).abs() < 1.0e-6 * source, "{} vs {source}", near.re);
    for distance in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        assert!(mutual_radiation_impedance(&gas, omega, distance, interaction(2.0)).is_err());
    }
}

#[test]
fn zero_length_chimneys_give_the_inverse_of_the_two_opening_radiation_matrix() {
    // Hand solution. With no chimney (A = D = 1, B = C = 0) each mouth pressure is the junction
    // pressure, p = Z U with Z = [[Zs1, Zm], [Zm, Zs2]], so the flow into the openings is Z^-1 p:
    //   Y = [[Zs2, -Zm], [-Zm, Zs1]] / (Zs1 Zs2 - Zm^2).
    let gas = air();
    let omega = TAU * 800.0;
    let (zs1, zs2) = (C64::new(2.0e4, -9.0e4), C64::new(3.0e4, -5.0e4));
    let ports = [CoupledPort::open_end(0.0, zs1), CoupledPort::open_end(0.05, zs2)];
    let zm = mutual_radiation_impedance(&gas, omega, 0.05, interaction(2.0)).unwrap();
    let determinant = zs1 * zs2 - zm * zm;
    let expected = [zs2 / determinant, -zm / determinant, -zm / determinant, zs1 / determinant];
    let coupling = couple(&ports, &gas, omega, interaction(2.0)).unwrap();
    for (got, want) in coupling.admittance.iter().zip(expected) {
        assert!(close(*got, want, 1.0e-12), "{got:?} vs {want:?}");
    }
}

fn chimney_ports(positions: &[f64]) -> Vec<CoupledPort> {
    let gas = air();
    positions
        .iter()
        .map(|&position| {
            CoupledPort::open_hole(position, &gas, &gas, 0.0025,
                0.0065 + fs_phs::side_hole_inner_length(0.0025, 0.00735), TAU * 700.0,
                LossModel::Bessel, WallModel::Rigid).unwrap()
        })
        .collect()
}

#[test]
fn distant_openings_decouple_to_the_independent_chimney_admittance() {
    // The chimney matrix and the self radiation must close to exactly the independent open-hole
    // impedance of the original owner; far apart, the coupling is all that remains and it vanishes.
    let gas = air();
    let omega = TAU * 700.0;
    let (b, h, a) = (0.0025, 0.0065, 0.00735);
    let (open, _) = crate::chimney_open_response_with_wall(
        &gas, &gas, b, h, a, omega, LossModel::Bessel, crate::hole_lengths::HoleLengths::Engine,
        WallModel::Rigid).unwrap();
    let ports = chimney_ports(&[0.0, 1.0e7]);
    let coupling = couple(&ports, &gas, omega, interaction(2.0)).unwrap();
    let independent = open.impedance.recip();
    assert!(close(coupling.admittance[0], independent, 1.0e-12), "{:?} vs {independent:?}", coupling.admittance[0]);
    assert!(close(coupling.admittance[3], independent, 1.0e-12));
    for off in [coupling.admittance[1], coupling.admittance[2]] {
        assert!(off.abs() < 1.0e-6 * independent.abs(), "{off:?}");
    }
    // Zs recovered from the chimney: mouth pressure over mouth flow.
    let [pa, pb, pc, pd] = ports[0].chimney;
    let zs = ports[0].self_radiation;
    assert!(close((pa * zs + pb) / (pc * zs + pd), open.impedance, 1.0e-12));
    assert!(close(pa * pd - pb * pc, C64::ONE, 1.0e-9), "reciprocal chimney");
}

#[test]
fn the_coupled_admittance_is_reciprocal() {
    let gas = air();
    let omega = TAU * 900.0;
    let mut ports = chimney_ports(&[0.10, 0.13, 0.19]);
    ports.push(CoupledPort::open_end(0.4, C64::new(1.2e4, -3.0e4)));
    let q = ports.len();
    let coupling = couple(&ports, &gas, omega, interaction(2.0)).unwrap();
    let scale = coupling.admittance.iter().map(|y| y.abs()).fold(0.0, f64::max);
    for i in 0..q {
        for j in 0..q {
            let (ij, ji) = (coupling.admittance[i * q + j], coupling.admittance[j * q + i]);
            assert!((ij - ji).abs() < 1.0e-9 * scale, "({i},{j}) {ij:?} vs {ji:?}");
        }
    }
    // Coupling is real: the off-diagonal terms are nonzero for stations 30 mm apart.
    assert!(coupling.admittance[1].abs() > 1.0e-3 * coupling.admittance[0].abs());
}
