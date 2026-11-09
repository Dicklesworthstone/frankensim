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

#[test]
fn mutual_loading_regularizes_an_isolated_chimney_pole_without_cancellation() {
    // G2: a lossless reactive limiting fixture, with a series compliance cancelling
    // the self-load inertance. At kd = pi the mutual load is reactive to rounding.
    // The isolated denominator is r, but the complete mouth system is
    // [[r, Zm], [Zm, r]], regular even at r = 0. Its closed-form inverse gives
    // finite mouth flows through the resonance; no artificial damping is needed.
    let gas = air();
    let distance = 0.05;
    let omega = core::f64::consts::PI * gas.sound_speed / distance;
    let zm = mutual_radiation_impedance(&gas, omega, distance, interaction(2.0)).unwrap();
    for relative_resistance in [1.0e-2, 1.0e-6, 1.0e-10, 1.0e-14, 0.0] {
        let r = C64::from_re(zm.abs() * relative_resistance);
        let make = |position_m| CoupledPort {
            position_m,
            radius_m: 0.002,
            chimney: [C64::ONE, r + C64::new(0.0, 1.0e5), C64::ZERO, C64::ONE],
            self_radiation: C64::new(0.0, -1.0e5),
            mouth_radius_m: 0.0,
        };
        let ports = [make(0.0), make(distance)];
        let determinant = r * r - zm * zm;
        let expected = [r / determinant, -zm / determinant, -zm / determinant, r / determinant];
        let coupling = couple(&ports, &gas, omega, interaction(2.0)).unwrap();
        assert!(coupling.condition < 1.1, "{}", coupling.condition);
        for map in [&coupling.admittance, &coupling.mouth_admittance] {
            for (got, want) in map.iter().zip(expected) {
                assert!((*got - want).abs() * zm.abs() < 2.0e-12,
                    "r/|Zm|={relative_resistance}: {got:?} vs {want:?}");
            }
        }
        // A genuine uncoupled pole still refuses; the fix does not regularize it.
        if relative_resistance == 0.0 {
            assert!(matches!(couple(&ports[..1], &gas, omega, interaction(2.0)), Err(DuctError::Singular)));
        }
    }
}

#[test]
fn unequal_lossless_chimneys_preserve_the_full_boundary_equations_and_acoustic_power() {
    // Independent checks at the original pressure/flow boundaries. Lossless
    // lines under exp(-i wt) have negative imaginary off-diagonal entries and
    // preserve Re(p conj(U)) even when their mouths exchange radiation power.
    let gas = air();
    let omega = TAU * 500.0;
    let ports: Vec<CoupledPort> = [(0.0, 0.2, 1.2e7), (0.04, 0.9, 2.0e7), (0.11, 1.7, 1.0e7)]
        .into_iter()
        .map(|(position_m, phase, reference): (f64, f64, f64)| {
            let (s, c) = phase.sin_cos();
            CoupledPort {
                position_m,
                radius_m: 0.002,
                chimney: [C64::from_re(c), C64::new(0.0, -reference * s),
                    C64::new(0.0, -s / reference), C64::from_re(c)],
                self_radiation: C64::new(2.0e5, -3.0e5),
                mouth_radius_m: 0.0,
            }
        })
        .collect();
    let junction = [C64::new(1.0, 0.4), C64::new(-0.5, 0.2), C64::new(0.3, -0.7)];
    let coupling = couple(&ports, &gas, omega, interaction(2.0)).unwrap();
    let response = |map: &[C64], i: usize| (0..3).fold(C64::ZERO, |sum, j| sum + map[i * 3 + j] * junction[j]);
    let mouth: Vec<_> = (0..3).map(|i| response(&coupling.mouth_admittance, i)).collect();
    let inflow: Vec<_> = (0..3).map(|i| response(&coupling.admittance, i)).collect();
    let power = |pressure: C64, flow: C64| pressure.re * flow.re + pressure.im * flow.im;
    let (mut bore_power, mut mouth_power) = (0.0, 0.0);
    for (i, port) in ports.iter().enumerate() {
        let mut pressure = port.self_radiation * mouth[i];
        for j in 0..3 {
            if i != j {
                pressure = pressure + mutual_radiation_impedance(
                    &gas, omega, pair_distance(port, &ports[j]), interaction(2.0)).unwrap() * mouth[j];
            }
        }
        let [a, b, c, d] = port.chimney;
        assert!(close(a * pressure + b * mouth[i], junction[i], 2.0e-12));
        assert!(close(c * pressure + d * mouth[i], inflow[i], 2.0e-12));
        bore_power += power(junction[i], inflow[i]);
        mouth_power += power(pressure, mouth[i]);
    }
    assert!(mouth_power > 0.0);
    assert!((bore_power - mouth_power).abs() <= 2.0e-11 * mouth_power,
        "bore={bore_power}, mouths={mouth_power}");
}

#[test]
fn the_joint_chimney_condition_respects_the_callers_limit() {
    let gas = air();
    let ports = [CoupledPort::open_end(0.0, C64::ONE), CoupledPort::open_end(0.1, C64::from_re(100.0))];
    let weak = interaction(1.0e9);
    assert!(matches!(couple_with_limit(&ports, &gas, 1.0, weak, 50.0), Err(DuctError::Singular)));
    let coupling = couple_with_limit(&ports, &gas, 1.0, weak, 200.0).unwrap();
    assert!(coupling.condition > 99.0 && coupling.condition < 101.0, "{}", coupling.condition);
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
    CoupledPort { position_m, radius_m, chimney: [C64::ONE, C64::ZERO, C64::ZERO, C64::ONE], self_radiation: load, mouth_radius_m: 0.0 }
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

#[test]
fn an_active_coupled_window_falls_back_to_independent_holes_and_stays_passive() {
    // Bead frankensim-1xe5t: the TMMI mutual impedance can make the loaded system active for large adjacent holes. Planted: the doubled-hole setup of the test
    // above (two zero-chimney mouths at one station, closed bore) at 55 kHz, where k d = 4.23 puts Re zm = s sin(kd) far below -Re z, so the input resistance
    // (z + zm) / 2 is negative. The coupled solve alone names the refusal; `solve` falls back to the independent holes (z / 2), marks it, and the result is passive.
    let gas = air();
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
    let omega = TAU * 55_000.0;
    let zm = mutual_radiation_impedance(&gas, omega, 2.0 * radius, interaction(2.0)).unwrap();
    assert!(load.re + zm.re < 0.0, "the planted window must be active: {load:?} {zm:?}");
    let coupled = solve_coupled(&chain, &holes, &ports, None, &gas, omega, interaction(2.0), 1.0e12);
    assert!(matches!(coupled, Err(DuctError::BadParameter { what: ACTIVE })), "the coupled solve alone must name the refusal");
    let solved = solve(&chain, &holes, &ports, None, &gas, omega, interaction(2.0), 1.0e12).unwrap();
    assert!(solved.independent_fallback);
    assert!(solved.input_impedance.re > 0.0, "{:?}", solved.input_impedance);
    // The independent holes in parallel, exactly the network with the mutual impedance off.
    assert!(close(solved.input_impedance, load.scale(0.5), 1.0e-9));
    for flow in &solved.hole_mouth_flows {
        assert!(close(*flow, C64::from_re(0.5), 1.0e-9), "{flow:?}");
    }
    let free = solve_coupled(&chain, &holes, &ports, None, &gas, omega, OpenHoleInteraction::independent(), 1.0e12).unwrap();
    assert_eq!(solved.input_impedance.re.to_bits(), free.input_impedance.re.to_bits());
    assert_eq!(solved.input_impedance.im.to_bits(), free.input_impedance.im.to_bits());
    // Where no window is active nothing changes: the same holes at 900 Hz take the coupled solve (the value the test above checks), unmarked, bit for bit.
    let calm = TAU * 900.0;
    let direct = solve_coupled(&chain, &holes, &ports, None, &gas, calm, interaction(2.0), 1.0e12).unwrap();
    let wrapped = solve(&chain, &holes, &ports, None, &gas, calm, interaction(2.0), 1.0e12).unwrap();
    assert!(!wrapped.independent_fallback && !direct.independent_fallback);
    assert_eq!(wrapped.input_impedance.re.to_bits(), direct.input_impedance.re.to_bits());
    assert_eq!(wrapped.input_impedance.im.to_bits(), direct.input_impedance.im.to_bits());
    // The independent network itself is never a fallback of a fallback.
    assert!(solve(&chain, &holes, &ports, None, &gas, omega, OpenHoleInteraction::independent(), 1.0e12).is_ok());
}

fn prepared_coupled_hole(block_index: usize, port: CoupledPort, reference: f64) -> PreparedOpenHole {
    let [a, b, c, d] = port.chimney;
    let inward = c * port.self_radiation + d;
    PreparedOpenHole {
        block_index,
        axial_position_m: port.position_m,
        effective_radius_m: port.radius_m,
        chimney_length_m: 0.0,
        sigma: 1.0,
        half_series_impedance: C64::ZERO,
        open_chimney_impedance: (a * port.self_radiation + b) / inward,
        closed_return_coordinates: [C64::ONE, C64::ZERO],
        chimney_u_mouth_over_u_in: inward.recip(),
        shunt_port_reference_pa_s_m3: reference,
        orifice_quadratic_coeff_pa_s2_m6: 0.0,
    }
}

fn series_resistance(resistance: f64, reference: f64) -> TwoPort {
    TwoPort::from_abcd([C64::ONE, C64::from_re(resistance), C64::ZERO, C64::ONE], reference).unwrap()
}

#[test]
fn the_complete_bore_regularizes_joint_chimney_poles_and_preserves_physical_power() {
    // G2 whole-instrument oracle. The passive two-mouth load is
    // Z = z*11^T + (i*Zc*tan(theta) + epsilon)*I, Re(z)>0, epsilon>=0.
    // Identical lossless chimneys cancel that diagonal reactance, so A Z+B
    // is rank one at epsilon=0. A real bore resistance between the openings
    // fixes the dark mouth-flow mode; their isolated admittance cannot do so.
    // The closed end carries no current. Summing/subtracting the ORIGINAL
    // chimney and bore equations gives the analytic sum/difference below.
    let gas = air();
    let omega = TAU * 700.0;
    let positions = [0.15, 0.20];
    let z = mutual_radiation_impedance(&gas, omega, positions[1] - positions[0], interaction(2.0)).unwrap();
    assert!(z.re > 0.0);
    let reference = z.abs();
    let middle_resistance = 0.9 * reference;
    let through = TwoPort::through(reference).unwrap();
    let check = |actual: C64, expected: C64, scale: f64| {
        assert!((actual - expected).abs() <= 2.0e-11 * scale.max(expected.abs()),
            "{actual:?} vs {expected:?}");
    };
    let power = |p: C64, u: C64| 0.5 * (p.re * u.re + p.im * u.im);
    for cosine in [1.0_f64, 0.5] {
        let sine = (1.0 - cosine * cosine).sqrt();
        let characteristic = 1.3 * reference;
        let a = C64::from_re(cosine);
        let b = C64::new(0.0, -characteristic * sine);
        let c = C64::new(0.0, -sine / characteristic);
        let diagonal_reactance = -b.scale(1.0 / cosine);
        for relative_resistance in [1.0e-3, 1.0e-7, 1.0e-11, 0.0] {
            let epsilon = relative_resistance * reference;
            let ports = positions.map(|position_m| CoupledPort {
                position_m,
                radius_m: 0.002,
                chimney: [a, b, c, a],
                self_radiation: z + diagonal_reactance + C64::from_re(epsilon),
                mouth_radius_m: 0.0,
            });
            let holes = [prepared_coupled_hole(1, ports[0], reference),
                prepared_coupled_hole(3, ports[1], reference)];
            let bound = 1.0e6;
            if relative_resistance <= 1.0e-7 {
                assert!(matches!(couple_with_limit(&ports, &gas, omega, interaction(2.0), bound),
                    Err(DuctError::Singular)), "the isolated chimney chart must refuse");
            }
            // A resistive inlet exercises regular bore admittances. Signed
            // through runs also put the bore itself exactly at a half/full wave.
            for inlet_phase in [None, Some(1.0), Some(-1.0)] {
                let phase = inlet_phase.unwrap_or(1.0);
                let inlet_resistance = if inlet_phase.is_some() { 0.0 } else { 0.7 * reference };
                let inlet = if inlet_phase.is_some() {
                    TwoPort::new([C64::ZERO, C64::from_re(phase), C64::from_re(phase), C64::ZERO], reference).unwrap()
                } else {
                    series_resistance(inlet_resistance, reference)
                };
                let chain = PreparedChain {
                    blocks: vec![inlet, through, series_resistance(middle_resistance, reference),
                        through, series_resistance(0.4 * reference, reference)],
                    reference_impedance: reference,
                    min_shear_number: f64::INFINITY,
                };
                let solved = solve(&chain, &holes, &ports, None, &gas, omega, interaction(2.0), bound).unwrap();
                let local = C64::from_re(1.0 / cosine) + c.scale(epsilon);
                let sum = (local + (c * z).scale(2.0)).recip().scale(phase);
                let difference = (local + C64::from_re(2.0 * cosine * epsilon / middle_resistance)).recip().scale(phase);
                let expected_mouth = [(sum + difference).scale(0.5), (sum - difference).scale(0.5)];
                let expected_pressure = expected_mouth.map(|u| (z * sum + u.scale(epsilon)).scale(cosine));
                check(solved.input_impedance, C64::from_re(inlet_resistance) + expected_pressure[0].scale(phase), reference);
                check(solved.bell_pressure, expected_pressure[1], reference);
                for i in 0..2 {
                    check(solved.hole_mouth_flows[i], expected_mouth[i], 1.0);
                }
                assert_eq!(solved.bell_flow, C64::ZERO);
                assert!(solved.condition <= bound);
                if relative_resistance <= 1.0e-7 {
                    assert!(solved.condition < 1000.0, "rejected chart condition leaked: {}", solved.condition);
                }

                // Primitive p/U checks and the whole acoustic power ledger;
                // no admittance from the implementation is used by this oracle.
                let mouth = &solved.hole_mouth_flows;
                let inflow = &solved.hole_chimney_inflows;
                let mouth_pressure: [C64; 2] = core::array::from_fn(|i|
                    ports[i].self_radiation * mouth[i] + z * mouth[1 - i]);
                let junction: [C64; 2] = core::array::from_fn(|i| a * mouth_pressure[i] + b * mouth[i]);
                for i in 0..2 {
                    check(inflow[i], c * mouth_pressure[i] + a * mouth[i], 1.0);
                    assert!((power(junction[i], inflow[i]) - power(mouth_pressure[i], mouth[i])).abs()
                        <= 2.0e-11 * reference);
                }
                check(inflow[0] + inflow[1], C64::from_re(phase), 1.0);
                check(junction[0] - junction[1], inflow[1].scale(middle_resistance), reference);
                check(solved.input_impedance, junction[0].scale(phase) + C64::from_re(inlet_resistance), reference);
                check(solved.bell_pressure, junction[1], reference);
                let radiation = power(mouth_pressure[0], mouth[0]) + power(mouth_pressure[1], mouth[1]);
                let bore_loss = 0.5 * (inlet_resistance + middle_resistance * inflow[1].abs().powi(2));
                assert!(radiation >= -2.0e-11 * reference);
                assert!((0.5 * solved.input_impedance.re - radiation - bore_loss).abs() <= 2.0e-11 * reference);
            }
        }
    }
}

#[test]
fn a_genuine_undamped_mouth_flow_mode_still_refuses() {
    // With a through connection between identical zero-chimney openings,
    // Z=z*11^T leaves an arbitrary circulating (U,-U) mode. The complete
    // physical system is now singular: it must not acquire artificial damping.
    let gas = air();
    let omega = TAU * 700.0;
    let positions = [0.15, 0.20];
    let z = mutual_radiation_impedance(&gas, omega, positions[1] - positions[0], interaction(2.0)).unwrap();
    let reference = z.abs();
    let ports = positions.map(|position| zero_chimney_port(position, 0.002, z));
    let holes = [prepared_coupled_hole(1, ports[0], reference), prepared_coupled_hole(3, ports[1], reference)];
    let through = TwoPort::through(reference).unwrap();
    let chain = PreparedChain {
        blocks: vec![series_resistance(reference, reference), through, through, through,
            series_resistance(reference, reference)],
        reference_impedance: reference,
        min_shear_number: f64::INFINITY,
    };
    assert!(matches!(solve(&chain, &holes, &ports, None, &gas, omega, interaction(2.0), 1.0e8),
        Err(DuctError::Singular)));
    assert!(matches!(port_admittance(&chain, &holes, &ports, None, &gas, omega, interaction(2.0), 1.0e8, 2),
        Err(DuctError::Singular)), "a matched inlet cannot fix the circulating dark mode");
}

#[test]
fn the_ordinary_loaded_nodal_system_must_respect_the_callers_condition_bound() {
    let gas = air();
    let reference = 1.0e6;
    let port = zero_chimney_port(0.1, 0.002, C64::from_re(100.0 * reference));
    let hole = prepared_coupled_hole(1, port, reference);
    let run = series_resistance(reference, reference);
    let chain = PreparedChain {
        blocks: vec![run, TwoPort::through(reference).unwrap(), run],
        reference_impedance: reference,
        min_shear_number: f64::INFINITY,
    };
    let bound = 20.0;
    assert!(run.admittance(bound).is_ok());
    assert!(couple_with_limit(&[port], &gas, TAU * 700.0, interaction(2.0), bound).is_ok());
    // Every component chart passes; only the complete system exceeds the bound.
    assert!(matches!(solve(&chain, &[hole], &[port], None, &gas, TAU * 700.0, interaction(2.0), bound),
        Err(DuctError::Singular)));
    let admitted = solve(&chain, &[hole], &[port], None, &gas, TAU * 700.0, interaction(2.0), 1.0e8).unwrap();
    assert!(close(admitted.input_impedance, C64::from_re(101.0 * reference), 2.0e-12));
    assert!(admitted.condition > bound && admitted.condition < 1.0e8);
}

#[test]
fn coupled_orifices_keep_mutual_radiation_at_exact_bore_resonances() {
    // Integer half-wave runs impose p_i=s_i*p_in, s_i=+/-1. With a
    // matched inlet, KCL gives p_in=-R*s^T*U, so the ORIGINAL mouth
    // equations reduce to (Zrad+R*s*s^T)*U=-dp. The two-by-two inverse
    // below is an analytic oracle, independent of the nodal/retained solve.
    let gas = air();
    let wavelength = 0.4;
    let omega = TAU * gas.sound_speed / wavelength;
    let reference = 1.0e5;
    let loads = [C64::new(0.7, -0.2), C64::new(0.4, -0.1), C64::new(0.9, -0.3)]
        .map(|z| z.scale(reference));
    let through = TwoPort::through(reference).unwrap();
    let signed = |phase: f64| TwoPort::new(
        [C64::ZERO, C64::from_re(phase), C64::from_re(phase), C64::ZERO], reference).unwrap();
    for middle_phase in [-1.0, 1.0] {
        for end_phase in [-1.0, 1.0] {
            let middle_length = wavelength * if middle_phase < 0.0 { 0.5 } else { 1.0 };
            let end_length = wavelength * if end_phase < 0.0 { 0.5 } else { 1.0 };
            let positions = [0.2, 0.2 + middle_length, 0.2 + middle_length + end_length];
            let ports = [zero_chimney_port(positions[0], 0.002, loads[0]),
                zero_chimney_port(positions[1], 0.002, loads[1])];
            let holes = [prepared_coupled_hole(1, ports[0], reference),
                prepared_coupled_hole(3, ports[1], reference)];
            let chain = PreparedChain {
                blocks: vec![signed(-1.0), through, signed(middle_phase), through, signed(end_phase)],
                reference_impedance: reference,
                min_shear_number: f64::INFINITY,
            };
            assert!(assemble_bore(&chain, &holes, 1.0e10).unwrap().keep_waves);
            let mutual = |i: usize, j: usize| mutual_radiation_impedance(
                &gas, omega, positions[j] - positions[i], interaction(2.0)).unwrap();
            for radiating_bell in [false, true] {
                let bell = radiating_bell.then(|| CoupledPort::open_end(positions[2], loads[2]));
                let mut a = loads[0] + C64::from_re(reference);
                let mut b = mutual(0, 1) + C64::from_re(reference * middle_phase);
                let mut d = loads[1] + C64::from_re(reference);
                if radiating_bell {
                    // The bell is linear (dp_bell=0). Eliminate its physical
                    // mouth current from the three-by-three analytic system.
                    let v0 = mutual(0, 2) + C64::from_re(reference * middle_phase * end_phase);
                    let v1 = mutual(1, 2) + C64::from_re(reference * end_phase);
                    let bell_diagonal = loads[2] + C64::from_re(reference);
                    a = a - v0 * v0 / bell_diagonal;
                    b = b - v0 * v1 / bell_diagonal;
                    d = d - v1 * v1 / bell_diagonal;
                }
                let determinant = a * d - b * b;
                let expected = [d / determinant, -b / determinant, -b / determinant, a / determinant];
                let observed = port_admittance(&chain, &holes, &ports, bell, &gas,
                    omega, interaction(2.0), 1.0e10, 2).unwrap();
                for (actual, expected) in observed.iter().zip(expected) {
                    assert!(close(*actual, expected, 2.0e-11), "{actual:?} vs {expected:?}");
                }
                assert!(close(observed[1], observed[2], 2.0e-12), "reciprocity");
                // Asking for one nonlinear hole still retains every linear
                // opening and its mutual radiation in that hole's response.
                let first = port_admittance(&chain, &holes, &ports, bell, &gas,
                    omega, interaction(2.0), 1.0e10, 1).unwrap();
                assert!(close(first[0], expected[0], 2.0e-11));
                let independent = port_admittance(&chain, &holes, &ports, bell, &gas,
                    omega, OpenHoleInteraction::independent(), 1.0e10, 2).unwrap();
                assert!((observed[0] - independent[0]).abs() > 1.0e-5 * observed[0].abs(),
                    "fixture must distinguish the former independent fallback");
                // Every complex orifice drive delivers nonnegative power
                // into the matched inlet and passive mouth loads.
                for drive in [[C64::ONE, C64::ONE], [C64::ONE, -C64::ONE],
                    [C64::ONE, C64::new(0.0, 1.0)]] {
                    let power = (0..2).map(|i| {
                        let inward = (0..2).fold(C64::ZERO, |sum, j| sum + observed[i * 2 + j] * drive[j]);
                        0.5 * (drive[i].re * inward.re + drive[i].im * inward.im)
                    }).sum::<f64>();
                    assert!(power > 0.0, "{power}");
                }
            }
        }
    }
}

#[test]
fn coupled_orifices_survive_joint_chimney_poles_without_changing_the_load() {
    // The individual chimney chart M=A*Z+B loses its differential mode,
    // while the real bore resistance between mouths uniquely fixes it.
    // With K=[[r,r],[r,r+rb]], the matched inlet/bore equations give
    // p=-K*U. Thus (M+K*N)*w=-dp and U=N*w, N=C*Z+D.
    let gas = air();
    let omega = TAU * 700.0;
    let positions = [0.15, 0.20];
    let z = mutual_radiation_impedance(&gas, omega, 0.05, interaction(2.0)).unwrap();
    let reference = z.abs();
    let r = 1.7 * reference; // Matched source plus the inlet bore resistance.
    let rb = 0.9 * reference;
    let through = TwoPort::through(reference).unwrap();
    let chain = PreparedChain {
        blocks: vec![series_resistance(0.7 * reference, reference), through,
            series_resistance(rb, reference), through, series_resistance(0.4 * reference, reference)],
        reference_impedance: reference,
        min_shear_number: f64::INFINITY,
    };
    for cosine in [1.0_f64, 0.5] {
        let sine = (1.0 - cosine * cosine).sqrt();
        let characteristic = 1.3 * reference;
        let c = C64::new(0.0, -sine / characteristic);
        for relative_resistance in [1.0e-3, 1.0e-7, 1.0e-11, 0.0] {
            let epsilon = relative_resistance * reference;
            let ports = positions.map(|position_m| CoupledPort {
                position_m,
                radius_m: 0.002,
                chimney: [C64::from_re(cosine), C64::new(0.0, -characteristic * sine), c, C64::from_re(cosine)],
                self_radiation: z + C64::new(epsilon, characteristic * sine / cosine),
                mouth_radius_m: 0.0,
            });
            let holes = [prepared_coupled_hole(1, ports[0], reference),
                prepared_coupled_hole(3, ports[1], reference)];
            let bound = 1.0e6;
            assert!(!assemble_bore(&chain, &holes, bound).unwrap().keep_waves,
                "this fixture must exercise the joint chimney chart, not bore poles");
            if relative_resistance <= 1.0e-7 {
                assert!(couple_with_limit(&ports, &gas, omega, interaction(2.0), bound).is_err());
            }
            let m0 = (z + C64::from_re(epsilon)).scale(cosine);
            let m1 = z.scale(cosine);
            let n1 = c * z;
            let n0 = n1 + C64::from_re(1.0 / cosine) + c.scale(epsilon);
            let h00 = m0 + (n0 + n1).scale(r);
            let h01 = m1 + (n0 + n1).scale(r);
            let h10 = m1 + n0.scale(r) + n1.scale(r + rb);
            let h11 = m0 + n1.scale(r) + n0.scale(r + rb);
            let determinant = h00 * h11 - h01 * h10;
            let expected = [(n0 * h11 - n1 * h10) / determinant,
                (n1 * h00 - n0 * h01) / determinant,
                (n1 * h11 - n0 * h10) / determinant,
                (n0 * h00 - n1 * h01) / determinant];
            let observed = port_admittance(&chain, &holes, &ports, None, &gas,
                omega, interaction(2.0), bound, 2).unwrap();
            for (actual, expected) in observed.iter().zip(expected) {
                assert!(close(*actual, expected, 2.0e-10), "eps={epsilon}, {actual:?} vs {expected:?}");
            }
            assert!(close(observed[1], observed[2], 2.0e-10));
        }
    }
}
