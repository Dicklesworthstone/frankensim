//! Independent force reference and actual loaded-bank contact trajectories.
use super::super::super::linear::BoardMode;
use super::*;

fn course() -> Course {
    Course {
        midi: 69,
        unison: 1,
        length_m: 1.,
        linear_density_kg_m: 0.02,
        tension_n: 1000.,
        flexural_rigidity_nm2: 0.,
        bonded_section: None,
        strike_fraction: 0.12,
        hammer_mass_kg: 0.008,
        felt_thickness_m: 0.008,
        felt_area_m2: 0.0001,
        duplex_length_m: 0.,
        detune_cents: 0.,
    }
}

fn text(stations: usize, normal: [f64; 2], elastic: bool, closed: f64) -> String {
    let material = if elastic {
        "felt,69,400000,0.2,2,2,0,0.8,2500000\n"
    } else {
        "felt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\nbranch,69,2000000,0.0002\nbranch,69,500000,0.004\n"
    };
    format!(
        "{HEADER}\npad,69,0.31,0.39,0.0001,0.008,{closed},0.002,0.5,{},{},{stations}\n{material}",
        normal[0], normal[1]
    )
}

fn prepared(
    c: Course,
    stations: usize,
    normal: [f64; 2],
    elastic: bool,
    closed: f64,
    coupled: bool,
) -> (Bank, Prepared) {
    let board = [BoardMode {
        frequency_hz: 173.,
        damping_ratio: 0.,
        bridge: [if coupled { 0.7 } else { 0. }; 88],
        volume: 0.1,
    }];
    let secondary = (normal[1] != 0.).then(|| vec![vec![0.4]]);
    let mut bank = Bank::new_with_transverse_bridge(
        &[c],
        &board,
        48_000,
        21_600.,
        if coupled { 8 } else { 1 },
        false,
        secondary.as_deref(),
    )
    .unwrap();
    let spec = Specification::read(&text(stations, normal, elastic, closed), &[c]).unwrap();
    let pads = Prepared::new(&spec, &[c], &bank).unwrap();
    bank.prepare_external_force();
    (bank, pads)
}

fn solve(bank: &mut Bank, pads: &mut Prepared) -> usize {
    let hammer = vec![0.; bank.contact_strings.len()];
    for trial in 0..64 {
        bank.set_external_force(pads.generalized_force()).unwrap();
        bank.predict();
        bank.finish(&hammer);
        if pads.correct(&bank.next_q).unwrap() {
            return trial;
        }
    }
    panic!("felt damper fixture exhausted the actual coupled trial budget");
}

fn tick(bank: &mut Bank, pads: &mut Prepared, sustain: f64, held: bool) -> (f64, f64, f64) {
    let before = bank.energy() + pads.energy();
    pads.begin_step(sustain, |_| held).unwrap();
    solve(bank, pads);
    let result = pads.accept(&bank.next_q).unwrap();
    let modal = bank.last_modal_loss_j;
    bank.commit();
    let defect = bank.energy() + pads.energy() - before + modal + result.0 - result.2;
    assert!(
        defect.abs() < 1e-10,
        "combined damper work defect {defect:e}"
    );
    result
}

#[test]
fn station_normals_project_both_fields_and_the_actual_moving_bridge_once() {
    let c = Course {
        unison: 3,
        duplex_length_m: 0.12,
        ..course()
    };
    let (mut bank, mut pads) = prepared(c, 4, [0.6, 0.8], false, 0., true);
    assert_eq!(pads.resolution(), (3, 12));
    let q: Vec<_> = (0..bank.q.len())
        .map(|k| 1e-5 * (0.73 * k as f64).sin())
        .collect();
    for site in &pads.sites {
        let fraction = 0.31 + 0.08 * (site.station as f64 + 0.5) / 4.;
        let mut expected = 0.;
        for plane in 0..2 {
            let s = bank
                .strings
                .iter()
                .find(|s| {
                    s.course == 0 && s.member == site.member && !s.duplex && s.polarization == plane
                })
                .unwrap();
            let bridge = s
                .bridge
                .iter()
                .zip(&q[bank.modes.len()..])
                .map(|(g, q)| g * q)
                .sum::<f64>();
            // Reconstruct fixed-end modal displacement plus the physical
            // endpoint interpolation independently of the prepared projection.
            let mut y = fraction * bridge;
            for (n, k) in s.modes.clone().enumerate() {
                let phi = ((n + 1) as f64 * std::f64::consts::PI * fraction).sin()
                    / c.modal_mass_kg().sqrt();
                y += phi * (q[k] - bank.modes[k].beta * bridge);
            }
            expected += [0.6, 0.8][plane] * y;
        }
        assert!((site.projection.position(&q) - expected).abs() < 1e-16);
        for (&index, &value) in site.projection.indices.iter().zip(&site.projection.values) {
            assert!(value.is_finite());
            if index < bank.modes.len() {
                let s = &bank.strings[bank.modes[index].string];
                assert!(!s.duplex && s.member == site.member);
            }
        }
    }
    for (i, force) in pads.force.iter_mut().enumerate() {
        *force = 0.1 + i as f64 * 0.03;
    }
    pads.regather().unwrap();
    let work: f64 = pads
        .generalized_force()
        .iter()
        .zip(&q)
        .map(|(f, q)| f * q)
        .sum();
    let contact: f64 = pads
        .sites
        .iter()
        .zip(&pads.force)
        .map(|(s, f)| f * s.projection.position(&q))
        .sum();
    assert!((work + contact).abs() < 1e-17);

    bank.predict();
    let hammer = vec![0.; bank.contact_strings.len()];
    bank.finish(&hammer);
    let free = bank.next_q.clone();
    bank.set_external_force(pads.generalized_force()).unwrap();
    bank.predict();
    bank.finish(&hammer);
    let group = &pads.groups[0];
    for row in 0..12 {
        let expected = -(0..12)
            .map(|col| group.compliance[row * 12 + col] * pads.force[col])
            .sum::<f64>();
        let actual = pads.sites[row].projection.position(&bank.next_q)
            - pads.sites[row].projection.position(&free);
        assert!((actual - expected).abs() < 1e-16);
        for col in 0..12 {
            assert!(
                (group.compliance[row * 12 + col] - group.compliance[col * 12 + row]).abs() < 1e-18
            );
        }
    }
    assert!(bank.q.iter().chain(&bank.v).all(|x| *x == 0.));
}

#[test]
fn conservative_first_contact_matches_an_independent_quadratic_root_and_surface_work() {
    let (mut bank, mut pads) = prepared(course(), 1, [1., 0.], true, 0., false);
    bank.v[0] = 0.5;
    let before = bank.energy();
    pads.begin_step(0., |_| false).unwrap();
    bank.predict();
    bank.finish(&[0.]);
    let free = pads.sites[0].projection.position(&bank.next_q);
    let compliance = pads.groups[0].compliance[0];
    // F(e)=K*e^2, and first-contact work averaging from e0=0 gives
    // Fbar=K*e1^2/3. Solve e1+H*K*e1^2/3=e_free independently.
    let k = 0.0001 * 400_000. / (0.008_f64 * 0.2).powi(2);
    let exact_end = 2. * free / (1. + (1. + 4. * compliance * k * free / 3.).sqrt());
    let exact_force = k * exact_end.powi(2) / 3.;
    assert!(exact_force > 1e-3);
    let saved_q = bank.q.clone();
    let saved_contact = pads.contact_state(0);
    let trials = solve(&mut bank, &mut pads);
    assert!(trials > 0 && trials < 64);
    assert_eq!(bank.q, saved_q);
    assert_eq!(
        pads.contact_state(0),
        saved_contact,
        "trials must retain accepted material history"
    );
    assert!((pads.force[0] - exact_force).abs() < 2e-5);
    assert!((pads.sites[0].projection.position(&bank.next_q) - exact_end).abs() < 1e-10);
    let (loss, viscous, work) = pads.accept(&bank.next_q).unwrap();
    let modal = bank.last_modal_loss_j;
    bank.commit();
    assert_eq!(viscous, 0.);
    assert_eq!(work, 0.);
    assert!(loss.abs() < 1e-9);
    assert!((bank.energy() + pads.energy() - before + modal + loss).abs() < 1e-12);
    assert!(pads.energy() > 0.);
    // Raising a still-loaded prescribed surface extracts work. It must not be
    // called material loss or omitted from the original mechanics ledger.
    let result = tick(&mut bank, &mut pads, 1., false);
    assert!(result.2 < 0.);
    assert_eq!(pads.pad_height_m(0), Some(0.5 / 48_000.));
    assert!(
        pads.accept(&bank.q).is_err(),
        "one tick cannot be accepted twice"
    );
}

#[test]
fn half_pedal_has_a_real_amplitude_threshold_and_recontact_retains_felt_memory() {
    let (mut quiet_bank, mut quiet) = prepared(course(), 1, [1., 0.], false, 0., false);
    let (mut loud_bank, mut loud) = prepared(course(), 1, [1., 0.], false, 0., false);
    for _ in 0..120 {
        assert_eq!(tick(&mut quiet_bank, &mut quiet, 0.5, false), (0., 0., 0.));
        assert_eq!(tick(&mut loud_bank, &mut loud, 0.5, false), (0., 0., 0.));
    }
    assert_eq!(quiet.pad_height_m(0), Some(0.001));
    let shape = quiet.sites[0].projection.values[0];
    let omega = quiet_bank.modes[0].omega;
    quiet_bank.v[0] = 0.0003 * omega / shape;
    loud_bank.v[0] = 0.0018 * omega / shape;
    let mut loss = 0.;
    let mut relaxation_loss = 0.;
    let mut peak_force = 0.0_f64;
    for _ in 0..900 {
        assert_eq!(tick(&mut quiet_bank, &mut quiet, 0.5, false), (0., 0., 0.));
        let out = tick(&mut loud_bank, &mut loud, 0.5, false);
        loss += out.0;
        relaxation_loss += out.1;
        peak_force = peak_force.max(loud.contact_state(0).unwrap().force_n);
    }
    assert_eq!(quiet.energy(), 0.);
    assert_eq!(quiet.contact_state(0).unwrap().maximum_felt_strain, 0.);
    assert!(loss > 0. && relaxation_loss > 0. && peak_force > 0.1);
    let crushed = loud.contact_state(0).unwrap().maximum_felt_strain;
    assert!(crushed > 0.);
    for _ in 0..400 {
        tick(&mut loud_bank, &mut loud, 1., false);
    }
    assert_eq!(loud.pad_height_m(0), Some(0.002));
    assert_eq!(loud.contact_state(0).unwrap().force_n, 0.);
    let lifted_crush = loud.contact_state(0).unwrap().maximum_felt_strain;
    assert!(lifted_crush >= crushed);
    let memory_before = loud.contact_state(0).unwrap().creep_coordinates_m;
    for _ in 0..120 {
        tick(&mut loud_bank, &mut loud, 0., true);
    }
    assert_eq!(
        loud.pad_height_m(0),
        Some(0.002),
        "key hold overrides released sustain"
    );
    assert_eq!(
        loud.contact_state(0).unwrap().maximum_felt_strain,
        lifted_crush
    );
    let memory_after = loud.contact_state(0).unwrap().creep_coordinates_m;
    assert!(memory_after.iter().zip(memory_before).all(|(a, b)| *a <= b));
    let mut contacted_again = false;
    for _ in 0..1000 {
        tick(&mut loud_bank, &mut loud, 0., false);
        contacted_again |= loud.contact_state(0).unwrap().force_n > 1e-4;
    }
    assert!(contacted_again);
    assert!(loud.contact_state(0).unwrap().maximum_felt_strain >= lifted_crush);
}

#[test]
fn accepted_travel_and_material_state_restore_after_later_trial_refusal() {
    let (mut bank, mut pads) = prepared(course(), 1, [1., 0.], false, 0., false);
    bank.v[0] = 0.05;
    for _ in 0..100 {
        tick(&mut bank, &mut pads, 0., false);
    }
    let q = bank.q.clone();
    let v = bank.v.clone();
    let contact = pads.contact_state(0);
    let energy = pads.energy();
    let height = pads.pad_height_m(0);
    pads.checkpoint();
    // An output frame may have accepted earlier substeps before its later
    // failure. Restore all of those accepted pad histories and surface travel.
    for _ in 0..3 {
        tick(&mut bank, &mut pads, 1., false);
    }
    assert_ne!(pads.pad_height_m(0), height);
    pads.begin_step(1., |_| false).unwrap();
    assert!(pads.correct(&vec![f64::NAN; bank.q.len()]).is_err());
    pads.restore();
    bank.q.copy_from_slice(&q);
    bank.v.copy_from_slice(&v);
    assert_eq!(pads.contact_state(0), contact);
    assert_eq!(pads.pad_height_m(0), height);
    assert_eq!(pads.energy(), energy);
    let first = tick(&mut bank, &mut pads, 0.75, false);
    let expected_q = bank.q.clone();
    let expected_v = bank.v.clone();
    let expected_contact = pads.contact_state(0);
    let expected_height = pads.pad_height_m(0);
    pads.restore();
    bank.q.copy_from_slice(&q);
    bank.v.copy_from_slice(&v);
    let retry = tick(&mut bank, &mut pads, 0.75, false);
    assert_eq!(retry, first);
    assert_eq!(bank.q, expected_q);
    assert_eq!(bank.v, expected_v);
    assert_eq!(pads.contact_state(0), expected_contact);
    assert_eq!(pads.pad_height_m(0), expected_height);
}

#[test]
fn incomplete_materials_unsupported_normals_and_nonphysical_clearances_are_refused() {
    let c = course();
    let good = text(1, [1., 0.], false, 0.);
    for bad in [
        good.replace("felt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\n", ""),
        good.replace("0.008,0,0.002", "0.008,-0.001,0.002"),
        good.replace("0.5,1,0,1", "0.5,1,1,1"),
        good.replace("0.5,1,0,1", "0.5,1,0,5"),
        good.replace("0.5,1,0,1", "0,1,0,1"),
        format!("{good}free,69\n"),
        good.replace("branch,69,2000000,0.0002", "branch,60,2000000,0.0002"),
    ] {
        assert!(Specification::read(&bad, &[c]).is_err(), "accepted {bad}");
    }
    let (bank, _) = prepared(c, 1, [1., 0.], false, 0., false);
    let lateral = Specification::read(&text(1, [0., 1.], false, 0.), &[c]).unwrap();
    assert!(Prepared::new(&lateral, &[c], &bank).is_err());
    let free = Specification::read(&format!("{HEADER}\nfree,69\n"), &[c]).unwrap();
    let mut pads = Prepared::new(&free, &[c], &bank).unwrap();
    assert_eq!(pads.resolution(), (0, 0));
    pads.begin_step(0.5, |_| false).unwrap();
    assert!(pads.correct(&bank.q).unwrap());
    assert_eq!(pads.accept(&bank.q).unwrap(), (0., 0., 0.));
    assert!(
        Specification::read(
            &format!("{HEADER}\nfree,69\nfelt,69,400000,0.2,2,2,0,0.8,2500000\n"),
            &[c]
        )
        .is_err()
    );
}
