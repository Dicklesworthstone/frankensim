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

#[test]
fn coupled_radiation_survives_half_and_full_wave_bore_runs() {
    // G2: at integer half wavelengths each lossless run is an exact signed
    // through connection. Its isolated admittance does not exist, but finite
    // passive aperture loads make the complete instrument perfectly regular.
    // A 0.4 m wavelength places this fixture near 858 Hz in room-temperature air.
    let gas = air();
    let wavelength = 0.4;
    let omega = TAU * gas.sound_speed / wavelength;
    let reference = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.007_f64.powi(2));
    let z1 = C64::new(0.7, -0.2).scale(reference);
    let z2 = C64::new(0.4, -0.1).scale(reference);
    let through = TwoPort::through(reference).unwrap();
    let signed = |phase: f64| TwoPort::new(
        [C64::ZERO, C64::from_re(phase), C64::from_re(phase), C64::ZERO], reference).unwrap();
    for inlet_phase in [-1.0, 1.0] {
        for bell_phase in [-1.0, 1.0] {
            let first_length = wavelength * if inlet_phase < 0.0 { 0.5 } else { 1.0 };
            let last_length = wavelength * if bell_phase < 0.0 { 0.5 } else { 1.0 };
            let chain = PreparedChain {
                blocks: vec![signed(inlet_phase), through, signed(bell_phase)],
                reference_impedance: reference,
                min_shear_number: f64::INFINITY,
            };
            assert!(chain.blocks[0].admittance(1.0e12).is_err(),
                "this fixture must exercise the isolated-component pole");
            let hole = PreparedOpenHole {
                block_index: 1,
                axial_position_m: first_length,
                effective_radius_m: 0.0025,
                chimney_length_m: 0.0,
                sigma: 1.0,
                half_series_impedance: C64::ZERO,
                open_chimney_impedance: z1,
                closed_return_coordinates: [C64::ONE, C64::ZERO],
                chimney_u_mouth_over_u_in: C64::ONE,
                shunt_port_reference_pa_s_m3: reference,
                orifice_quadratic_coeff_pa_s2_m6: 0.0,
            };
            let ports = [CoupledPort::open_end(first_length, z1)];
            let bell = CoupledPort::open_end(first_length + last_length, z2);
            let zm = mutual_radiation_impedance(&gas, omega, last_length, interaction(2.0)).unwrap();
            // Mouth pressures have relative sign bell_phase. For v=[1,t],
            // Zin=1/(v^T Zrad^-1 v); the inlet sign squares to one in Zin but
            // remains audible in the two radiating mouth-flow phasors.
            let denominator = z1 + z2 - zm.scale(2.0 * bell_phase);
            let expected_z = (z1 * z2 - zm * zm) / denominator;
            let expected_hole = (z2 - zm.scale(bell_phase)).scale(inlet_phase) / denominator;
            let expected_bell = (z1.scale(bell_phase) - zm).scale(inlet_phase) / denominator;
            let observed = solve(&chain, &[hole], &ports, Some(bell), &gas,
                omega, interaction(2.0), 1.0e12).unwrap();
            assert!(close(observed.input_impedance, expected_z, 1e-11));
            assert!(close(observed.bell_pressure, expected_z.scale(inlet_phase * bell_phase), 1e-11));
            assert!(close(observed.hole_mouth_flows[0], expected_hole, 1e-11));
            assert!(close(observed.hole_chimney_inflows[0], expected_hole, 1e-11));
            assert!(close(observed.bell_flow, expected_bell, 1e-11));
            assert!(close(observed.hole_chimney_inflows[0]
                + observed.bell_flow.scale(bell_phase), C64::from_re(inlet_phase), 1e-11));
            assert!(observed.condition.is_finite() && observed.condition < 1e12);
            assert!(observed.input_impedance.re > 0.0);
            let closed = solve(&chain, &[hole], &ports, None, &gas,
                omega, interaction(2.0), 1.0e12).unwrap();
            assert!(close(closed.input_impedance, z1, 1e-11));
            assert!(close(closed.bell_pressure, z1.scale(inlet_phase * bell_phase), 1e-11));
            assert!(close(closed.hole_mouth_flows[0], C64::from_re(inlet_phase), 1e-11));
            assert_eq!(closed.bell_flow, C64::ZERO);
        }
    }
}

fn zero_chimney_port(position_m: f64, radius_m: f64, load: C64) -> CoupledPort {
    CoupledPort { position_m, radius_m, chimney: [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE], self_radiation: load }
}

#[test]
fn a_pair_is_coupled_at_its_axial_distance_unless_the_mouths_would_overlap() {
    // Separate holes of a real bore (axial spacing above the radii) keep the original distance bit for bit.
    let load = C64::new(1.2e4, -3.0e4);
    let (a, b) = (zero_chimney_port(0.10, 0.0025, load), zero_chimney_port(0.13, 0.0025, load));
    assert_eq!(pair_distance(&a, &b), (0.10_f64 - 0.13).abs());
    assert_eq!(pair_distance(&b, &a), pair_distance(&a, &b));
    // The bore's open end has radius 0: only the axial distance counts.
    let end = zero_chimney_port(0.4, 0.0, load);
    assert_eq!(pair_distance(&a, &end), (0.10_f64 - 0.4).abs());
    // Mouths 0.5 mm apart on the axis cannot overlap: they are at least the sum of the radii apart.
    let close_by = zero_chimney_port(0.1005, 0.0025, load);
    assert_eq!(pair_distance(&a, &close_by), 0.005);
    // One station: the touching distance, with no spacing invented.
    let same = zero_chimney_port(0.10, 0.0021, load);
    let twin = zero_chimney_port(0.10, 0.0021, load);
    assert_eq!(pair_distance(&same, &twin), 0.0042);
}

#[test]
fn two_holes_at_one_station_couple_like_two_holes_a_touching_distance_apart() {
    let gas = air();
    let omega = TAU * 900.0;
    let load = C64::new(2.0e4, -9.0e4);
    let radius = 0.0021;
    let same = [zero_chimney_port(0.230, radius, load), zero_chimney_port(0.230, radius, load)];
    let apart = [zero_chimney_port(0.230, radius, load), zero_chimney_port(0.230 + 2.0 * radius, radius, load)];
    let (a, b) = (couple(&same, &gas, omega, interaction(2.0)).unwrap(), couple(&apart, &gas, omega, interaction(2.0)).unwrap());
    for (x, y) in a.admittance.iter().zip(&b.admittance) {
        assert!(close(*x, *y, 1.0e-12), "{x:?} vs {y:?}");
    }
    // Hand solution for two identical zero-chimney mouths: Z = [[z, zm], [zm, z]], so Y = [[z, -zm], [-zm, z]] / (z^2 - zm^2),
    // with zm the monopole impedance at the touching distance 2 r.
    let zm = mutual_radiation_impedance(&gas, omega, 2.0 * radius, interaction(2.0)).unwrap();
    let determinant = load * load - zm * zm;
    let expected = [load / determinant, -zm / determinant, -zm / determinant, load / determinant];
    for (got, want) in a.admittance.iter().zip(expected) {
        assert!(close(*got, want, 1.0e-12), "{got:?} vs {want:?}");
    }
    // Two zero-radius openings at one station still refuse: no distance exists.
    let bells = [zero_chimney_port(0.2, 0.0, load), zero_chimney_port(0.2, 0.0, load)];
    assert!(couple(&bells, &gas, omega, interaction(2.0)).is_err());
}

#[test]
fn a_doubled_hole_solves_and_shares_the_flow_between_its_two_mouths() {
    // The Bb4 oboe card has two holes at one station. With ideal through bore blocks around them and no chimney, both mouths see the
    // junction pressure p, so the input impedance of the closed bore is 1 / (1^T Z^-1 1) = (z + zm) / 2 with Z = [[z, zm], [zm, z]],
    // zm at the touching distance, and each mouth carries half of the unit flow.
    let gas = air();
    let omega = TAU * 900.0;
    let reference = gas.density * gas.sound_speed / (core::f64::consts::PI * 0.007_f64.powi(2));
    let load = C64::new(0.6, -0.4).scale(reference);
    let through = TwoPort::through(reference).unwrap();
    let chain = PreparedChain { blocks: vec![through; 4], reference_impedance: reference, min_shear_number: f64::INFINITY };
    let radius = 0.0021;
    let hole = |block_index: usize| PreparedOpenHole {
        block_index,
        axial_position_m: 0.230,
        effective_radius_m: radius,
        chimney_length_m: 0.0,
        sigma: 1.0,
        half_series_impedance: C64::ZERO,
        open_chimney_impedance: load,
        closed_return_coordinates: [C64::ONE, C64::ZERO],
        chimney_u_mouth_over_u_in: C64::ONE,
        shunt_port_reference_pa_s_m3: reference,
        orifice_quadratic_coeff_pa_s2_m6: 0.0,
    };
    let holes = [hole(1), hole(2)];
    let ports = [zero_chimney_port(0.230, radius, load), zero_chimney_port(0.230, radius, load)];
    let zm = mutual_radiation_impedance(&gas, omega, 2.0 * radius, interaction(2.0)).unwrap();
    let solved = solve(&chain, &holes, &ports, None, &gas, omega, interaction(2.0), 1.0e12).unwrap();
    let expected = (load + zm).scale(0.5);
    assert!(close(solved.input_impedance, expected, 1.0e-10), "{:?} vs {expected:?}", solved.input_impedance);
    for flow in &solved.hole_mouth_flows {
        assert!(close(*flow, C64::from_re(0.5), 1.0e-10), "{flow:?}");
    }
    assert!(solved.input_impedance.re > 0.0);
    // Uncoupled (eps -> infinity) it is the two independent holes in parallel, z / 2.
    let free = solve(&chain, &holes, &ports, None, &gas, omega, interaction(1.0e15), 1.0e12).unwrap();
    assert!(close(free.input_impedance, load.scale(0.5), 1.0e-9));
}
