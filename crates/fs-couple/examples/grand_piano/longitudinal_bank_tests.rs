use super::super::super::geometry;
use super::super::BoardMode;
use super::*;

fn course() -> Course {
    Course {
        length_m: 1.0,
        linear_density_kg_m: 0.01,
        tension_n: 400.0,
        flexural_rigidity_nm2: 0.0,
        unison: 1,
        duplex_length_m: 0.0,
        detune_cents: 0.0,
        ..geometry::demonstration_scale().unwrap()[48]
    }
}

fn board() -> [BoardMode; 2] {
    [
        BoardMode {
            frequency_hz: 85.0,
            damping_ratio: 0.0,
            bridge: [0.8; 88],
            volume: 0.1,
        },
        BoardMode {
            frequency_hz: 240.0,
            damping_ratio: 0.0,
            bridge: [-0.45; 88],
            volume: -0.05,
        },
    ]
}

fn specification(courses: &[Course], band: f64, cap: usize, damping: f64) -> Specification {
    let mut text = format!("{HEADER}\n");
    for c in courses {
        text.push_str(&format!(
            "waves,{},100000,0.2,{band},{cap},{damping},0.05\n",
            c.midi
        ));
    }
    Specification::read(&text, courses).unwrap().with_fixture_axial_rows(&board())
}

fn advance(bank: &mut Bank, contacts: &[f64]) -> Result<(), &'static str> {
    bank.begin_string_stretching_step();
    for _ in 0..MAX_TRIALS {
        bank.predict();
        bank.finish(contacts);
        if bank.correct_string_stretching_step()? {
            return Ok(());
        }
    }
    Err("longitudinal string trial budget exceeded")
}

#[test]
fn longitudinal_cards_admit_the_complete_scale_and_enforce_band_and_coordinate_budget() {
    let c = course();
    let second = Course { midi: 70, ..c };
    let good = format!("{HEADER}\nwaves,69,100000,0.2,3300,2,0.005,0.05\nlinear,70\n");
    let spec = Specification::read(&good, &[c, second]).unwrap();
    assert_eq!(spec.waves.len(), 1);
    for row in [
        "waves,69,100000,0.2,3300,0,0.005,0.05",
        "waves,69,100000,0.2,3300,33,0.005,0.05",
        "waves,69,100000,0.2,NaN,2,0.005,0.05",
        "waves,69,100000,0.2,3300,2,-0.005,0.05",
        "waves,69,100000,0.2,3300,2,NaN,0.05",
        "waves,69,100000,0.2,3300,2,0.005,0.11",
        "waves,69,100000,0.2,3300,2,0.005,0",
        "waves,69,100000,0.2,3300,2,0.005",
    ] {
        assert!(
            Specification::read(&format!("{HEADER}\n{row}"), &[c]).is_err(),
            "{row}"
        );
    }
    assert!(Specification::read(&good, &[c]).is_err());
    assert!(
        Specification::read(
            &format!("{HEADER}\nwaves,69,100000,0.2,3300,2,0.005,0.05"),
            &[c, second]
        )
        .is_err()
    );
    let build = || Bank::new(&[c], &board(), 48_000, 3000.0, 4, false).unwrap();
    let mut admitted = build();
    admitted
        .configure_string_stretching(&[c], &specification(&[c], 3300.0, 2, 0.0))
        .unwrap();
    assert_eq!(admitted.longitudinal_mode_count(), 2);
    assert_eq!(admitted.string_longitudinal_state(0).unwrap().0.len(), 2);
    assert!(
        build()
            .configure_string_stretching(&[c], &specification(&[c], 5000.0, 2, 0.0))
            .is_err()
    );
    assert!(
        build()
            .configure_string_stretching(&[c], &specification(&[c], 30_000.0, 32, 0.0))
            .is_err()
    );
    let mut condensed = build();
    condensed
        .configure_string_stretching(&[c], &specification(&[c], 1000.0, 2, 0.0))
        .unwrap();
    assert_eq!(condensed.longitudinal_mode_count(), 0);
    assert!(condensed.has_string_stretching());
    assert!(condensed.string_longitudinal_state(0).unwrap().0.is_empty());
}

#[test]
fn longitudinal_force_and_work_remain_reciprocal_for_two_planes_unisons_and_duplexes() {
    let c = Course {
        unison: 3,
        duplex_length_m: 0.4,
        detune_cents: 3.0,
        ..course()
    };
    let mut bank = Bank::new_with_transverse_bridge(
        &[c],
        &board(),
        192_000,
        6000.0,
        4,
        false,
        Some(&[vec![0.3, 0.6]]),
    )
    .unwrap();
    bank.configure_string_stretching(&[c], &specification(&[c], 6500.0, 4, 0.0))
        .unwrap();
    assert_eq!(bank.strings.len(), 12);
    assert_eq!(bank.longitudinal_mode_count(), 15);
    for (j, q) in bank.q.iter_mut().enumerate() {
        *q = 1e-6 * (0.31 * j as f64 + 0.2).cos();
    }
    let state = bank.stretching.as_mut().unwrap();
    assert_eq!(state.channels.len(), 6);
    for (ci, channel) in state.channels.iter_mut().enumerate() {
        let wave = channel.waves.as_mut().unwrap();
        for (j, q) in wave.q.iter_mut().enumerate() {
            *q = 1e-8 * (j as f64 + ci as f64 + 0.4).sin();
        }
    }
    // Fixture preparation puts current and candidate axial coordinates at the
    // same nonzero state, so finite differences vary only actual bank motion.
    state.predict();
    state.commit();
    state
        .evaluate(&bank.q, &bank.q, &bank.strings, &bank.modes)
        .unwrap();
    let force = state.required.clone();
    for k in 0..bank.q.len() {
        let h = 1e-9;
        let mut a = bank.q.clone();
        let mut b = bank.q.clone();
        a[k] -= h;
        b[k] += h;
        let fd = (state.energy(&b, &bank.strings, &bank.modes)
            - state.energy(&a, &bank.strings, &bank.modes))
            / (2.0 * h);
        assert!(
            (force[k] + fd).abs() < 1e-7 * (1.0 + fd.abs()),
            "coordinate {k}: force={} fd={fd}",
            force[k]
        );
    }
    let mut end = bank.q.clone();
    for (j, q) in end.iter_mut().enumerate() {
        *q += 2e-7 * (j as f64 + 0.1).cos();
    }
    state
        .evaluate(&bank.q, &end, &bank.strings, &bank.modes)
        .unwrap();
    let work = state
        .required
        .iter()
        .zip(bank.q.iter().zip(&end))
        .map(|(force, (a, b))| force * (b - a))
        .sum::<f64>();
    let delta = state.energy(&end, &bank.strings, &bank.modes)
        - state.energy(&bank.q, &bank.strings, &bank.modes);
    assert!((delta + work).abs() < 1e-13);
    for channel in &state.channels {
        let other = channel
            .second
            .expect("both transverse directions share axial inertia");
        let first = channel.observe(&bank.q, &bank.strings, &bank.modes);
        assert!(first.minimum_tension_n > 0.0);
        assert_eq!(
            bank.strings[channel.string].member,
            bank.strings[other].member
        );
        assert_eq!(
            bank.strings[channel.string].duplex,
            bank.strings[other].duplex
        );
    }
}

/// A per-mode row that is NOT proportional to the vertical bridge row (the one thing the former scalar could not
/// express): the axial reaction on each board coordinate must still be the exact gradient of the stored energy,
/// and the work over a finite step must still close.
#[test]
fn a_per_mode_axial_row_keeps_force_and_work_reciprocal() {
    let c = Course { unison: 2, duplex_length_m: 0.4, detune_cents: 3.0, ..course() };
    let mut bank = Bank::new(&[c], &board(), 192_000, 6000.0, 4, false).unwrap();
    let mut spec = specification(&[c], 6500.0, 4, 0.0);
    // Vertical rows are 0.8 and -0.45; these are not a multiple of them (opposite sign ratio, different size).
    spec.set_bridge_axial_rows(BTreeMap::from([(c.midi, vec![0.35, 0.9])])).unwrap();
    bank.configure_string_stretching(&[c], &spec).unwrap();
    for (j, q) in bank.q.iter_mut().enumerate() {
        *q = 1e-6 * (0.31 * j as f64 + 0.2).cos();
    }
    let state = bank.stretching.as_mut().unwrap();
    for (ci, channel) in state.channels.iter_mut().enumerate() {
        for (j, q) in channel.waves.as_mut().unwrap().q.iter_mut().enumerate() {
            *q = 1e-8 * (j as f64 + ci as f64 + 0.4).sin();
        }
    }
    state.predict();
    state.commit();
    state.evaluate(&bank.q, &bank.q, &bank.strings, &bank.modes).unwrap();
    let force = state.required.clone();
    let mut strongest = 0.0f64;
    for k in 0..bank.q.len() {
        let h = 1e-9;
        let (mut a, mut b) = (bank.q.clone(), bank.q.clone());
        a[k] -= h;
        b[k] += h;
        let fd = (state.energy(&b, &bank.strings, &bank.modes) - state.energy(&a, &bank.strings, &bank.modes)) / (2.0 * h);
        strongest = strongest.max(fd.abs());
        assert!((force[k] + fd).abs() < 1e-7 * (1.0 + fd.abs()), "coordinate {k}: force={} fd={fd}", force[k]);
    }
    assert!(strongest > 0.0, "the check would be vacuous with no stored-energy gradient");
    let mut end = bank.q.clone();
    for (j, q) in end.iter_mut().enumerate() {
        *q += 2e-7 * (j as f64 + 0.1).cos();
    }
    state.evaluate(&bank.q, &end, &bank.strings, &bank.modes).unwrap();
    let work = state.required.iter().zip(bank.q.iter().zip(&end)).map(|(f, (a, b))| f * (b - a)).sum::<f64>();
    let delta = state.energy(&end, &bank.strings, &bank.modes) - state.energy(&bank.q, &bank.strings, &bank.modes);
    assert!((delta + work).abs() < 1e-13, "work {work} against stored energy change {delta}");
}

/// The axial coupling is geometry from the board, never defaulted: a `waves` key without its row, with a row of
/// the wrong length, or with a row for a key that has no `waves` row refuses, and leaves the bank unconfigured.
#[test]
fn waves_without_an_axial_participation_row_refuse() {
    let c = course();
    let text = format!("{HEADER}\nwaves,{},100000,0.2,3300,2,0.0,0.05\n", c.midi);
    let bare = Specification::read(&text, &[c]).unwrap();
    let mut bank = Bank::new(&[c], &board(), 96_000, 6000.0, 4, false).unwrap();
    let refusal = bank.configure_string_stretching(&[c], &bare).unwrap_err();
    assert!(refusal.contains("axial participation row") && refusal.contains("no default"), "{refusal}");
    assert!(!bank.has_string_stretching());
    // A row for the wrong number of board modes.
    let mut short = bare.clone();
    short.set_bridge_axial_rows(BTreeMap::from([(c.midi, vec![0.1])])).unwrap();
    assert!(bank.configure_string_stretching(&[c], &short).unwrap_err().contains("entries for 2 board modes"));
    assert!(!bank.has_string_stretching());
    // Rows that cannot be admitted at all: no waves key, empty, non-finite.
    let mut spec = bare.clone();
    assert!(spec.set_bridge_axial_rows(BTreeMap::from([(c.midi + 1, vec![0.1, 0.2])])).is_err());
    assert!(spec.set_bridge_axial_rows(BTreeMap::from([(c.midi, vec![])])).is_err());
    assert!(spec.set_bridge_axial_rows(BTreeMap::from([(c.midi, vec![0.1, f64::NAN])])).is_err());
    // With its row the same file configures.
    spec.set_bridge_axial_rows(BTreeMap::from([(c.midi, vec![0.1, -0.2])])).unwrap();
    bank.configure_string_stretching(&[c], &spec).unwrap();
    assert!(bank.has_string_stretching());
}

#[test]
fn longitudinal_energy_transfer_and_physical_loss_close_the_original_driven_work_ledger() {
    let c = Course {
        unison: 3,
        duplex_length_m: 0.4,
        detune_cents: 3.0,
        ..course()
    };
    let mut bank = Bank::new(&[c], &board(), 96_000, 6000.0, 4, false).unwrap();
    bank.configure_string_stretching(&[c], &specification(&[c], 6500.0, 4, 0.006))
        .unwrap();
    let initial = bank.energy();
    let mut work = 0.0;
    let mut loss = 0.0;
    let mut peak_axial_kinetic = 0.0_f64;
    for tick in 0..600 {
        let forces = vec![if tick < 120 { 0.2 } else { 0.0 }; bank.contact_strings.len()];
        advance(&mut bank, &forces).unwrap();
        work += (0..forces.len())
            .map(|j| {
                forces[j]
                    * (bank.contact_position(j, &bank.next_q) - bank.contact_position(j, &bank.q))
            })
            .sum::<f64>();
        loss += bank.last_modal_loss_j;
        assert!((bank.next_energy() - initial + loss - work).abs() < 1e-9 + 1e-8 * work.abs());
        bank.commit();
        peak_axial_kinetic = peak_axial_kinetic.max(
            bank.strings
                .iter()
                .enumerate()
                .filter(|(_, s)| s.polarization == 0)
                .map(|(i, _)| {
                    bank.string_stretching_observation(i)
                        .unwrap()
                        .longitudinal_kinetic_energy_j
                })
                .sum(),
        );
        assert!((bank.energy() - initial + loss - work).abs() < 1e-9 + 1e-8 * work.abs());
    }
    assert!(work > 0.0 && loss > 0.0 && peak_axial_kinetic > 0.0);
    assert!(bank.longitudinal_energy_j() > 0.0);
}

#[test]
fn longitudinal_trajectory_converges_to_independent_continuous_axial_transverse_equations() {
    let c = course();
    let ea = 100_000.0;
    let amplitude = 0.003 * det::sqrt(c.modal_mass_kg());
    let omega_t =
        core::f64::consts::PI / c.length_m * det::sqrt(c.tension_n / c.linear_density_kg_m);
    let slope_scale = core::f64::consts::PI / c.length_m / det::sqrt(c.modal_mass_kg());
    let g = 0.25 * slope_scale * slope_scale;
    let stiffness = 0.5 * ea * c.length_m;
    let d = 2.0 * slope_scale;
    let omega_l = d * det::sqrt(stiffness);
    // Independent continuous four-state vector field from the span's analytic
    // single-sine strain. This does not call the production energy or force.
    let derivative = |x: [f64; 4]| {
        let [q, v, u, w] = x;
        let strain = d * u + g * q * q;
        [
            v,
            -omega_t * omega_t * q
                - 4.0 * stiffness * g * g * q * q * q
                - 2.0 * stiffness * g * q * strain,
            w,
            -stiffness * d * strain,
        ]
    };
    let mut reference = [amplitude, 0.0, 0.0, 0.0];
    let reference_rate = 1_536_000;
    let dt = 1.0 / reference_rate as f64;
    for _ in 0..reference_rate / 100 {
        let a = derivative(reference);
        let b = derivative(std::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
        let c = derivative(std::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
        let d = derivative(std::array::from_fn(|i| reference[i] + dt * c[i]));
        for j in 0..4 {
            reference[j] += dt * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
        }
    }
    let run = |rate| {
        let fixed = [BoardMode {
            frequency_hz: 85.0,
            damping_ratio: 0.0,
            bridge: [0.0; 88],
            volume: 0.0,
        }];
        let mut bank = Bank::new(&[c], &fixed, rate, 3000.0, 1, false).unwrap();
        // This board has one mode (all-zero vertical row), so the fixture row has one (zero) entry too.
        bank.configure_string_stretching(
            &[c],
            &specification(&[c], 3300.0, 2, 0.0).with_fixture_axial_rows(&fixed),
        )
        .unwrap();
        bank.q[0] = amplitude;
        let initial = bank.energy();
        for _ in 0..rate / 100 {
            advance(&mut bank, &[0.0]).unwrap();
            assert!((bank.next_energy() - initial).abs() < 1e-10);
            bank.commit();
        }
        let (u, v) = bank.string_longitudinal_state(0).unwrap();
        assert_eq!(u.len(), 2);
        assert_eq!(u[0], 0.0); // no chord or odd geometric-strain harmonic
        assert_eq!(v[0], 0.0);
        assert!(u[1] != 0.0 && v[1] != 0.0);
        [bank.q[0], bank.v[0], u[1], v[1]]
    };
    let axial_scale = g * amplitude * amplitude / d;
    let scales = [
        amplitude,
        omega_t * amplitude,
        axial_scale,
        omega_l * axial_scale,
    ];
    let error = |actual: [f64; 4]| {
        det::sqrt(
            actual
                .iter()
                .zip(reference)
                .zip(scales)
                .map(|((a, b), scale)| ((a - b) / scale).powi(2))
                .sum(),
        )
    };
    let coarse_error = error(run(48_000));
    let fine_error = error(run(96_000));
    eprintln!(
        "longitudinal independent trajectory: coarse={coarse_error:.8e}, fine={fine_error:.8e}"
    );
    assert!(
        fine_error < 0.5 * coarse_error,
        "coarse={coarse_error}, fine={fine_error}"
    );
    assert!(
        fine_error < 0.01,
        "normalized full phase-space error={fine_error}"
    );
}

#[test]
fn bridge_port_without_retained_axial_modes_stores_the_work_it_does() {
    // A band below the first axial mode (1581 Hz here) retains none, yet the
    // channel keeps its quasi-static boundary extension eta*b and books that
    // port's bridge reaction: the potential must include it, or the force does
    // work no energy sees (the reference refused every such key at forte).
    let c = course();
    let mut bank = Bank::new(&[c], &board(), 96_000, 3000.0, 2, false).unwrap();
    bank.configure_string_stretching(&[c], &specification(&[c], 1000.0, 2, 0.0))
        .unwrap();
    assert_eq!(bank.longitudinal_mode_count(), 0);
    let board_offset = bank.modes.len();
    bank.q[0] = 1e-4;
    bank.q[board_offset] = 1e-5;
    let initial = bank.energy();
    assert!(initial > 0.0);
    for _ in 0..120 {
        advance(&mut bank, &[0.0]).unwrap();
        assert!((bank.next_energy() - initial).abs() < 1e-12 * initial);
        bank.commit();
    }
}

#[test]
fn axial_strain_drives_bridge_at_zero_initial_bridge_displacement_and_closes_work() {
    let c = course();
    let mut bank = Bank::new(&[c], &board(), 96_000, 3000.0, 2, false).unwrap();
    bank.configure_string_stretching(&[c], &specification(&[c], 3300.0, 2, 0.0))
        .unwrap();
    // Excite only an internal longitudinal wave coordinate while all transverse
    // string and soundboard coordinates start at rest (bridge b = 0).
    {
        let state = bank.stretching.as_mut().unwrap();
        let wave = state.channels[0].waves.as_mut().unwrap();
        wave.q[0] = 2e-7;
        state.predict();
        state.commit();
        state
            .evaluate(&bank.q, &bank.q, &bank.strings, &bank.modes)
            .unwrap();
        let board_offset = bank.modes.len();
        assert!(
            state.required[board_offset..].iter().any(|f| f.abs() > 1e-3),
            "axial strain must exert direct force on the bridge even at zero bridge displacement"
        );
    }
    let initial = bank.energy();
    assert!(initial > 0.0);
    let board_offset = bank.modes.len();
    for _ in 0..120 {
        advance(&mut bank, &[0.0]).unwrap();
        assert!((bank.next_energy() - initial).abs() < 1e-12 * initial);
        bank.commit();
    }
    let board_kinetic: f64 = bank.v[board_offset..].iter().map(|v| 0.5 * v * v).sum();
    assert!(
        board_kinetic > 0.0,
        "axial wave must transfer energy into soundboard motion through the bridge port"
    );
}

