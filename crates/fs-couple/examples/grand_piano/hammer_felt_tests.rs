//! The felt pack reproduces its source's printed values, refuses what it cannot
//! answer, states its slow-loading error truthfully, and a harder card shortens
//! the contact and brightens the sound through the contact law alone.
use super::*;
use fs_matdb::{ClaimSet, PropertyValue};
use fs_material::Uniaxial;

const COEFFICIENTS: [&str; 12] = [
    "hammer_instantaneous_force_at_1mm_regression_prefactor",
    "hammer_instantaneous_force_regression_exponential_rate",
    "hammer_compression_exponent_regression_constant",
    "hammer_compression_exponent_regression_linear_coefficient",
    "hammer_hereditary_amplitude_regression_constant",
    "hammer_hereditary_amplitude_regression_linear_coefficient",
    "hammer_relaxation_time_regression_constant",
    "hammer_relaxation_time_regression_linear_coefficient_magnitude_negative",
    "hammer_relaxation_time_regression_quadratic_coefficient",
    "hammer_mass_regression_constant",
    "hammer_mass_regression_linear_coefficient_magnitude_negative",
    "hammer_mass_regression_quadratic_coefficient",
];

/// A copy of `source` under `id`, with its observations, and with every scalar
/// claim passed through `edit`; `None` drops the claim.
fn edited(
    source: &NormalizedPack,
    id: &str,
    edit: impl Fn(&str, f64) -> Option<f64>,
) -> NormalizedPack {
    let mut claims = ClaimSet::new();
    for id in source.claims().observation_ids() {
        let dataset = source.claims().observation(id).unwrap().clone();
        assert_eq!(claims.register_observation(dataset).unwrap(), id);
    }
    for (_, claim) in source.claims().claims_ordered() {
        let mut claim = claim.clone();
        if let PropertyValue::Scalar { value, .. } = &mut claim.value {
            match edit(claim.key.name(), *value) {
                Some(edited) => *value = edited,
                None => continue,
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    NormalizedPack::new(
        id,
        "hammer-felt-test-v1",
        source.source_artifact(),
        source.redistribution_terms(),
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

#[test]
fn the_set_reproduces_its_printed_values() {
    let card = stulov_abel().unwrap();
    assert_eq!(card.pack_id(), "hammer-felt-abel-stulov2005");
    let h = |midi| card.hammer(midi).unwrap();
    let per_mm = |x: Hammer| x.instantaneous_n_m_p / det::pow(1.0e3, x.exponent);
    // Thesis p.16: hammer 1 weighs 11.0 g and hammer 88 5.3 g.
    assert!(
        (h(21).mass_kg - 11.0e-3).abs() < 0.05e-3,
        "{}",
        h(21).mass_kg
    );
    assert!(
        (h(108).mass_kg - 5.3e-3).abs() < 0.05e-3,
        "{}",
        h(108).mass_kg
    );
    // Acta Mech. 225 eqs (4)-(7): the printed ranges, to one unit of their last
    // printed digit (their rounding is not uniform: 0.997144 is printed 0.9972
    // and 1.65696 microseconds 1.65).
    assert!((per_mm(h(21)) - 16_440.0).abs() < 5.0, "{}", per_mm(h(21)));
    assert!(
        (per_mm(h(108)) - 2_787_300.0).abs() < 50.0,
        "{}",
        per_mm(h(108))
    );
    assert!((h(21).exponent - 3.72).abs() < 1e-2);
    assert!((h(21).hereditary_amplitude - 0.9895).abs() < 1e-4);
    assert!((h(108).hereditary_amplitude - 0.9972).abs() < 1e-4);
    assert!((h(21).relaxation_time_s - 2.70e-6).abs() < 5e-9);
    assert!((h(108).relaxation_time_s - 1.65e-6).abs() < 1e-8);
    // Thesis p.17, hammer 85: p = 4.975 and alpha = 0.5312 ms. The printed
    // alpha(n) polynomial contradicts that value and is not used; the
    // four-parameter set gives gamma tau0 / (1 - gamma) within 0.5 % of it.
    let a7 = h(105);
    assert!((a7.exponent - 4.975).abs() < 1e-12);
    let alpha = a7.hereditary_amplitude * a7.relaxation_time_s / (1.0 - a7.hereditary_amplitude);
    assert!((alpha / 0.5312e-3 - 1.0).abs() < 0.01, "{alpha}");
    // Toward the treble the hammers get lighter and stiffer.
    for midi in 22..=108 {
        assert!(h(midi).mass_kg < h(midi - 1).mass_kg);
        assert!(h(midi).instantaneous_n_m_p > h(midi - 1).instantaneous_n_m_p);
    }
    // The course card spreads the whole hammer's K u^p over its patch.
    let c = super::super::steinway_scale::courses().unwrap()[48];
    let (law, elastic) = h(69).material(&c).unwrap();
    for strain in [0.01, 0.05, 0.1] {
        let force = c.felt_area_m2 * law.stress(strain, &law.initial_state());
        let expected = h(69).static_n_m_p() * det::pow(strain * c.felt_thickness_m, h(69).exponent);
        assert!((force / expected - 1.0).abs() < 1e-12);
    }
    assert!(elastic.terms.is_empty());
    assert!(
        (h(69).rate_per_string(&c) * c.unison as f64 / h(69).rate_n_s_m_p() - 1.0).abs() < 1e-15
    );
    assert!(h(70).material(&c).is_err());
}

#[test]
fn a_missing_coefficient_key_or_inadmissible_law_refuses() {
    assert!(load_pack(&[]).is_err());
    assert!(load_pack(&STULOV_ABEL[..STULOV_ABEL.len() / 2]).is_err());
    let source = load_pack(STULOV_ABEL).unwrap();
    for missing in COEFFICIENTS {
        let partial = edited(&source, "felt-without-coefficient", |key, value| {
            (key != missing).then_some(value)
        });
        let error = FeltCard::from_pack(&partial)
            .unwrap()
            .hammer(69)
            .unwrap_err();
        assert!(
            error.contains("felt-without-coefficient") && error.contains(missing),
            "{missing}: {error}"
        );
    }
    let card = stulov_abel().unwrap();
    assert!(card.hammer(20).is_err() && card.hammer(109).is_err());
    // A hereditary amplitude reaching 1, as SK2009's printed epsilon(n) does
    // from hammer 11, has no relaxed stiffness: that hammer refuses by name.
    let past_one = edited(&source, "felt-gamma-past-one", |key, value| {
        Some(
            if key == "hammer_hereditary_amplitude_regression_linear_coefficient" {
                2.0e-4
            } else {
                value
            },
        )
    });
    let card = FeltCard::from_pack(&past_one).unwrap();
    assert!(card.hammer(21).is_ok());
    let error = card.hammer(108).unwrap_err();
    assert!(error.contains("outside the hereditary model"), "{error}");
}

/// The module states how far the slow-loading form is from the exact
/// hereditary law. Half-sine compression u = U sin(pi t / T), U = 0.5 mm; the
/// one-term memory m (tau0 m' = u^p - m) advances exactly for u^p linear on a
/// step, and F = F0 (u^p - gamma m).
#[test]
fn the_slow_loading_error_is_what_the_module_states() {
    let card = stulov_abel().unwrap();
    for (midi, period_s, bound) in [
        (21, 4.0e-3, 2.0e-3),
        (60, 2.0e-3, 6.0e-3),
        (105, 0.8e-3, 2.2e-2),
    ] {
        let h = card.hammer(midi).unwrap();
        let (p, gamma, tau) = (h.exponent, h.hereditary_amplitude, h.relaxation_time_s);
        let amplitude = 0.5e-3;
        let steps = 400_000;
        let dt = period_s / f64::from(steps);
        let decay = det::exp(-dt / tau);
        // Clamped: the last step's sine can round a hair below zero.
        let compression =
            |t: f64| amplitude * det::sin(std::f64::consts::PI * t / period_s).max(0.0);
        let g = |t: f64| det::pow(compression(t), p);
        let (mut memory, mut peak, mut worst) = (0.0f64, 0.0f64, 0.0f64);
        for k in 1..=steps {
            let t = f64::from(k) * dt;
            let (g0, g1) = (g(t - dt), g(t));
            memory =
                decay * memory + (1.0 - decay) * g0 + (g1 - g0) * (1.0 - tau / dt * (1.0 - decay));
            let exact = h.instantaneous_n_m_p * (g1 - gamma * memory);
            let phase = std::f64::consts::PI * t / period_s;
            let rate = p
                * det::pow(compression(t), p - 1.0)
                * amplitude
                * det::cos(phase)
                * std::f64::consts::PI
                / period_s;
            let slow = h.static_n_m_p() * g1 + h.rate_n_s_m_p() * rate;
            peak = peak.max(exact.abs());
            worst = worst.max((exact - slow).abs());
        }
        eprintln!(
            "hammer {}: T = {period_s} s, slow-loading error {:.3e} of the peak force",
            h.number,
            worst / peak
        );
        assert!(
            worst / peak < bound,
            "hammer {}: {}",
            h.number,
            worst / peak
        );
    }
}

/// One strike on the demonstration board: contact frames (the felt's rate
/// loss grows only while the hammer touches), the output's spectral centroid
/// and its energy fraction from the eighth partial up (8 x the nominal f1).
fn strike(card: &FeltCard, midi: u8, velocity: f64) -> (usize, f64, f64) {
    let mut c = super::super::steinway_scale::courses().unwrap()[usize::from(midi - 21)];
    let h = card.hammer(midi).unwrap();
    c.hammer_mass_kg = h.mass_kg;
    let mut p = super::super::engine::Instrument::new_with_contact_geometry(
        vec![c],
        &super::super::board::demonstration(),
        48_000,
        4,
        12,
        true,
        vec![h.material(&c).unwrap()],
        Some(super::super::engine::ShankGeometry::published()),
        None,
    )
    .unwrap();
    p.configure_source_hammer_dissipation(&[h.rate_per_string(&c)])
        .unwrap();
    p.note_on(midi, velocity).unwrap();
    let (mut losses, mut output) = (Vec::new(), Vec::new());
    for _ in 0..8192 {
        output.push(p.step().unwrap());
        losses.push(p.accounting.felt_relaxation_loss_j);
    }
    let total = *losses.last().unwrap();
    assert!(total > 0.0);
    let contact = losses
        .windows(2)
        .filter(|w| w[1] - w[0] > 1e-9 * total)
        .count();
    let n = output.len();
    let twiddle: Vec<(f64, f64)> = (0..n)
        .map(|j| det::sin_cos(std::f64::consts::TAU * j as f64 / n as f64))
        .collect();
    let windowed: Vec<f64> = output
        .iter()
        .enumerate()
        .map(|(i, x)| x * (0.5 - 0.5 * twiddle[i].1))
        .collect();
    let eighth_partial_hz = 8.0 * 440.0 * det::pow(2.0, (f64::from(midi) - 69.0) / 12.0);
    let (mut weighted, mut power, mut high) = (0.0, 0.0, 0.0);
    for k in 1..n / 2 {
        let (mut re, mut im) = (0.0, 0.0);
        for (i, x) in windowed.iter().enumerate() {
            let (s, co) = twiddle[k * i % n];
            re += x * co;
            im -= x * s;
        }
        let e = re * re + im * im;
        let hz = k as f64 * 48_000.0 / n as f64;
        weighted += e * hz;
        power += e;
        if hz >= eighth_partial_hz {
            high += e;
        }
    }
    (contact, weighted / power, high / power)
}

#[test]
fn a_harder_felt_card_shortens_contact_and_brightens_the_sound() {
    let source = load_pack(STULOV_ABEL).unwrap();
    let measured = FeltCard::from_pack(&source).unwrap();
    let harder = FeltCard::from_pack(&edited(&source, "felt-four-times-stiffer", |key, value| {
        Some(
            if key == "hammer_instantaneous_force_at_1mm_regression_prefactor" {
                4.0 * value
            } else {
                value
            },
        )
    }))
    .unwrap();
    for midi in [48, 69] {
        let (soft_contact, soft_centroid, soft_high) = strike(&measured, midi, 2.5);
        let (hard_contact, hard_centroid, hard_high) = strike(&harder, midi, 2.5);
        eprintln!(
            "key {midi} at 2.5 m/s: contact {soft_contact} -> {hard_contact} frames, centroid {soft_centroid:.1} -> {hard_centroid:.1} Hz, energy from the 8th partial up {soft_high:.3e} -> {hard_high:.3e}"
        );
        assert!(
            soft_contact > 0 && hard_contact < soft_contact,
            "{soft_contact} -> {hard_contact}"
        );
        assert!(
            hard_centroid > soft_centroid && hard_high > soft_high,
            "centroid {soft_centroid} -> {hard_centroid}, from the 8th partial {soft_high} -> {hard_high}"
        );
    }
}
