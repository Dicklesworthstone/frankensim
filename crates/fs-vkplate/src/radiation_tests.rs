//! Head radiation: the shapes against independent references, the shipped
//! table against the analytic low-ka dipole, refusals, and the stepper's
//! energy account with implicit air mass and centered resistance, including
//! the shared response used during contact and friction.

use super::*;
use crate::radiation::{ModeShape, RadiationTable, SHIPPED_TABLE, bessel_j, bessel_zero};

#[test]
fn bessel_zeros_and_j1_match_independent_references() {
    // Abramowitz and Stegun Table 9.5 / DLMF 10.21.
    for (m, n, zero) in [
        (0, 1, 2.404_825_557_695_773),
        (1, 1, 3.831_705_970_207_512),
        (1, 2, 7.015_586_669_815_619),
        (2, 1, 5.135_622_301_840_683),
        (5, 1, 8.771_483_815_959_954),
        (10, 1, 14.475_500_686_554_54),
    ] {
        let found = bessel_zero(m, n);
        assert!(
            (found - zero).abs() < 1e-12 * zero,
            "j_{m},{n} = {found}, not {zero}"
        );
    }
    // fs-math's separately certified J1.
    for x in [0.01, 0.7, 3.0, 9.5, 17.0, 31.0, 45.0] {
        let (ours, reference) = (bessel_j(1, x), fs_math::special::bessel_j1(x));
        assert!(
            (ours - reference).abs() < 1e-13,
            "J1({x}): {ours} against {reference}"
        );
    }
}

#[test]
fn mode_shapes_have_unit_mean_square_and_are_orthogonal() {
    let shapes = [
        ModeShape::new(1, 1),
        ModeShape::new(1, 2),
        ModeShape::new(3, 1),
        ModeShape::new(3, 2),
    ];
    let h = 1.0 / 400.0;
    let mut gram = [[0.0; 8]; 8];
    let nodes = (-400..=400).flat_map(|j| (-400..=400).map(move |i| (i as f64 * h, j as f64 * h)));
    for (x, y) in nodes {
        let values: Vec<f64> = shapes
            .iter()
            .flat_map(|s| [s.value(x, y, false), s.value(x, y, true)])
            .collect();
        for a in 0..8 {
            for b in 0..8 {
                gram[a][b] += values[a] * values[b] * h * h;
            }
        }
    }
    for (a, row) in gram.iter().enumerate() {
        for (b, g) in row.iter().enumerate() {
            let expected = if a == b { PI } else { 0.0 };
            assert!((g - expected).abs() < 2e-3 * PI, "Gram[{a}][{b}] = {g}");
        }
    }
}

#[test]
fn the_shipped_table_meets_the_analytic_baffled_dipole_at_low_ka() {
    // An m = 1 shape radiates as an in-plane dipole of moment d q', with
    // d = integral phi x dA = sqrt(2) pi / j, so Re zhat -> ka^4 d^2 / (6 pi)
    // = ka^4 pi / (3 j^2) in an infinite baffle.
    let table = RadiationTable::shipped().unwrap();
    for (r, shape) in table
        .modes()
        .iter()
        .enumerate()
        .filter(|(_, s)| s.m == 1 && s.n <= 3)
    {
        for ka in [0.05, 0.1] {
            let z = table.zhat(r, r, ka).unwrap();
            let dipole = ka.powi(4) * PI / (3.0 * shape.zero * shape.zero);
            let resistance = z.resistance.expect("resolved at this ka");
            assert!(
                (resistance / dipole - 1.0).abs() < 1e-2,
                "(1,{}) at ka {ka}: {resistance} against {dipole}",
                shape.n
            );
            assert!(z.reactance < 0.0, "mass-like");
        }
    }
    let (lo, hi) = table.ka_range();
    assert!(table.zhat(0, 0, 0.5 * lo).is_err() && table.zhat(0, 0, 2.0 * hi).is_err());
}

#[test]
fn a_table_for_another_owner_or_with_damaged_rows_refuses() {
    assert!(RadiationTable::parse(SHIPPED_TABLE).is_ok());
    let damage = |from: &str, to: &str| SHIPPED_TABLE.replacen(from, to, 1);
    for text in [
        damage(
            "frankensim.vkplate-head-radiation.v1",
            "frankensim.vkplate-head-radiation.v2",
        ),
        damage("route\tplanar-baffled-rayleigh", "route\tclosed-body-bem"),
        damage(
            "geometry\tunit-disk-in-infinite-rigid-baffle",
            "geometry\tunit-square",
        ),
        damage("mode\t1\t1\t3.8317", "mode\t1\t1\t3.9317"),
        damage("resistance-resolved-from\t1\t", "note\t"),
        damage("\nka\t", "\nka\t-"),
        // A shape evaluated differently by the bake.
        damage("probe\t1\t", "probe\t2\t"),
    ] {
        assert!(RadiationTable::parse(&text).is_err());
    }
}

/// A clamped tensioned head without film loss, carrying the shipped load.
fn loaded_head(grid: f64, nonlinear: Nonlinear) -> Body {
    let (thickness, rest) = (|_: f64| 1.9e-4, |_: f64| 0.0);
    let mut head = Body::plate_without_surface_loss(
        "radiating head",
        &PlateSpec {
            radius: 0.32,
            hole: 0.0,
            grid,
            edge: Edge::Clamped,
            thickness: &thickness,
            rest_height: &rest,
            young: 4.9e9,
            poisson: 0.38,
            density: 1390.0,
            tension: 5360.0,
            bending: true,
            nonlinear,
            t60: [(100.0, 1.0), (1000.0, 1.0)],
            support: None,
            muffle: None,
        },
    )
    .unwrap();
    let table = RadiationTable::shipped().unwrap();
    head.attach_head_radiation(&table, 0.32, &standard_gas(), None)
        .unwrap();
    head
}

#[test]
fn g0_head_air_mass_and_resistance_close_the_discrete_energy() {
    // Berger tension modulation shares the Woodbury solve with the air mass.
    let mut head = loaded_head(0.02, Nonlinear::Berger);
    assert_eq!(head.loss(), (0.0, 0.0));
    let shapes = [ModeShape::new(1, 1), ModeShape::new(2, 1)];
    for (i, p) in head.xy.clone().iter().enumerate() {
        let (x, y) = (p[0] / 0.32, p[1] / 0.32);
        head.w[i] = 1e-3 * (shapes[0].value(x, y, false) + 0.5 * shapes[1].value(x, y, true));
        head.wp[i] = head.w[i];
    }
    head.sync();
    let mut ins = Instrument::new(48_000.0);
    ins.add(head);
    ins.validate_head_radiation().unwrap();
    assert!(ins.stability().0 < 1.0);
    let k = ins.k();
    let initial = ins.bodies[0].discrete_energy(k);
    let mut previous = initial;
    for _ in 0..2400 {
        ins.step(&mut []);
        let body = &ins.bodies[0];
        let energy = body.discrete_energy(k);
        assert!(
            energy <= previous + 1e-12 * initial,
            "radiation only removes energy"
        );
        assert!((energy + body.head_radiated_j() - initial).abs() < 1e-9 * initial);
        previous = energy;
    }
    assert!(ins.bodies[0].head_radiated_j() > 1e-3 * initial);
    let mut copy = ins.clone();
    for _ in 0..64 {
        ins.step(&mut []);
        copy.step(&mut []);
    }
    assert_eq!(ins.bodies[0].w, copy.bodies[0].w);
}

#[test]
fn g0_one_contact_on_a_radiating_head_closes_structural_air_and_contact_energy() {
    let mut head = loaded_head(0.04, Nonlinear::Berger);
    head.air_mass = 0.025;
    let node = head.node_at(0.12, 0.08).unwrap();
    let mut ins = Instrument::new(48_000.0);
    let k = ins.k();
    ins.add(head);
    let mut striker = Body::point("manufactured striker", 0.004);
    striker.w[0] = 1e-5;
    striker.wp[0] = striker.w[0] + 0.8 * k;
    ins.add(striker);
    let stiffness = 2e7;
    ins.contacts.push(Contact::new(
        vec![(0, node, 1.0)],
        vec![(1, 0, 1.0)],
        0.0,
        Law::Power {
            k: stiffness,
            alpha: 1.5,
            mu: 0.0,
        },
    ));
    ins.validate_head_radiation().unwrap();
    let energy = |ins: &Instrument| {
        let d = (ins.bodies[0].w[node] - ins.bodies[1].w[0]).max(0.0);
        let dp = (ins.bodies[0].wp[node] - ins.bodies[1].wp[0]).max(0.0);
        // The contact's two-step discrete gradient stores the average of
        // its two endpoint potentials at the same half step as the bodies.
        ins.bodies.iter().map(|b| b.discrete_energy(k)).sum::<f64>()
            + 0.5 * stiffness / 2.5 * (d * d * d.sqrt() + dp * dp * dp.sqrt())
    };
    let initial = energy(&ins);
    let mut sav_seen = false;
    for _ in 0..1200 {
        ins.step(&mut []);
        sav_seen |= ins.bodies[0].sav;
        let closed = energy(&ins) + ins.bodies[0].head_radiated_j();
        assert!(
            (closed - initial).abs() < 1e-8 * initial,
            "one contact must use the loaded body's inertia: {closed} vs {initial}"
        );
    }
    assert!(sav_seen && ins.contacts[0].peak_force > 0.0);
    assert!(ins.bodies[0].head_radiated_j() > 0.0);
}

/// Independent multiplication by the physical step operator; no prepared
/// inverse or Woodbury solve enters this reference.
fn response_force(body: &Body, displacement: &[f64], k: f64) -> Vec<f64> {
    let mut force: Vec<_> = displacement
        .iter()
        .zip(&body.cf)
        .map(|(q, c)| q / c[0])
        .collect();
    if body.sav {
        let a = 0.25 * dot(&body.g, displacement);
        for (f, g) in force.iter_mut().zip(&body.g) {
            *f += a * g;
        }
    }
    if body.air_mass > 0.0 {
        let a = body.air_mass * body.h.powi(4) / (body.area * body.area * k * k)
            * displacement.iter().sum::<f64>();
        for f in &mut force {
            *f += a;
        }
    }
    if let Some(load) = body.head_radiation.as_deref() {
        let p: Vec<_> = load.cols.iter().map(|u| dot(u, displacement)).collect();
        for (r, u) in load.cols.iter().enumerate() {
            let a: f64 = p
                .iter()
                .enumerate()
                .map(|(s, q)| {
                    (load.mass[r * p.len() + s] / (k * k)
                        + load.damping[r * p.len() + s] / (2.0 * k))
                        * q
                })
                .sum();
            for (f, u) in force.iter_mut().zip(u) {
                *f += a * u;
            }
        }
    }
    force
}

#[test]
fn g0_friction_on_a_loaded_head_uses_the_same_nonlocal_force_response() {
    let k = 1.0 / 48_000.0;
    let mut head = loaded_head(0.04, Nonlinear::Berger);
    head.air_mass = 0.025;
    let shape = ModeShape::new(1, 1);
    for (i, p) in head.xy.iter().enumerate() {
        let a = shape.value(p[0] / 0.32, p[1] / 0.32, false);
        head.w[i] = 3e-4 * a;
        head.wp[i] = head.w[i] - 0.03 * k * a;
    }
    head.sync();
    head.linear_forces();
    head.predict(k);
    assert!(head.sav);
    let a = head.node_at(0.12, 0.08).unwrap();
    let b = head.node_at(0.16, 0.08).unwrap();
    let c = head.node_at(0.12, 0.12).unwrap();
    for cap in [1e-5, 1e5] {
        let mut bodies = vec![head.clone()];
        let mut normal = Contact::new(
            vec![],
            vec![],
            0.0,
            Law::Power {
                k: 1.0,
                alpha: 1.5,
                mu: 0.0,
            },
        );
        normal.force = 1.0; // prescribed normal preload for this slip test
        let mut friction = [Friction {
            contact: 0,
            sx: vec![(0, a, 1.0), (0, b, -1.0)],
            sy: vec![(0, a, 1.0), (0, c, -1.0)],
            mu: cap,
            force: 0.0,
            stuck: false,
        }];
        let before = bodies[0].wn.clone();
        solve_frictions(&mut bodies, &[normal], &mut friction, k);
        let body = &bodies[0];
        let change: Vec<_> = body.wn.iter().zip(&before).map(|(x, y)| x - y).collect();
        let force = response_force(body, &change, k);
        let (fx, fy) = (-force[b], -force[c]);
        let scale = friction[0].force;
        assert!(scale > 0.0);
        for (i, f) in force.iter().enumerate() {
            let expected = if i == a {
                fx + fy
            } else if i == b {
                -fx
            } else if i == c {
                -fy
            } else {
                0.0
            };
            assert!(
                (f - expected).abs() < 1e-8 * scale,
                "a local friction force must satisfy the full loaded operator at node {i}: {f} vs {expected}"
            );
        }
        assert!((fx.hypot(fy) - scale).abs() < 1e-8 * scale);
        assert!(
            change
                .iter()
                .enumerate()
                .any(|(i, x)| i != a && i != b && i != c && *x != 0.0)
        );
        assert_eq!(friction[0].stuck, cap > 1.0);
        if friction[0].stuck {
            let slip = |j| (body.wn[a] - body.w[a]) - (body.wn[j] - body.w[j]);
            assert!(slip(b).hypot(slip(c)) < 1e-12);
        } else {
            assert!((scale / cap - 1.0).abs() < 1e-10);
        }
    }
}

#[test]
fn g1_a_single_shape_decays_at_its_narrowband_radiation_rate() {
    let mut head = loaded_head(0.016, Nonlinear::None);
    let load = head.head_radiation().unwrap().clone();
    // The receipt names each shape's loaded frequency and radiation T60.
    let line = load
        .receipt
        .iter()
        .find(|l| l.starts_with("(1,1) cos"))
        .unwrap();
    let t60: f64 = line
        .split("radiation T60 ")
        .nth(1)
        .and_then(|rest| rest.split(" s;").next())
        .unwrap()
        .parse()
        .unwrap();
    let mac: f64 = line
        .split("; MAC ")
        .nth(1)
        .and_then(|rest| rest.split(';').next())
        .unwrap()
        .parse()
        .unwrap();
    assert!(
        mac > 0.95,
        "the (1,1) shape is close to the FD head's own mode: {line}"
    );
    let shape = ModeShape::new(1, 1);
    for (i, p) in head.xy.clone().iter().enumerate() {
        head.w[i] = 1e-4 * shape.value(p[0] / 0.32, p[1] / 0.32, false);
        head.wp[i] = head.w[i];
    }
    head.sync();
    let mut ins = Instrument::new(48_000.0);
    ins.add(head);
    ins.validate_head_radiation().unwrap();
    let k = ins.k();
    let initial = ins.bodies[0].discrete_energy(k);
    let span = (0.1 * t60 / k) as usize;
    for _ in 0..span {
        ins.step(&mut []);
    }
    // Energy decays at 2 sigma = 13.8 / T60 for one radiating shape.
    let decay = ins.bodies[0].discrete_energy(k) / initial;
    let expected = (-2.0 * 6.907755 / t60 * span as f64 * k).exp();
    assert!(
        (decay.ln() / expected.ln() - 1.0).abs() < 0.05,
        "energy ratio {decay} against {expected}"
    );
}

#[test]
fn g0_a_shape_outside_the_table_is_refused_alone_unless_it_is_the_fundamental() {
    // Films without bending: the membrane wave speed sets every shape's ka.
    let attach = |tension: f64| {
        let (thickness, rest) = (|_: f64| 1.9e-4, |_: f64| 0.0);
        let mut head = Body::plate_without_surface_loss(
            "film head",
            &PlateSpec {
                radius: 0.32,
                hole: 0.0,
                grid: 0.02,
                edge: Edge::Clamped,
                thickness: &thickness,
                rest_height: &rest,
                young: 4.9e9,
                poisson: 0.38,
                density: 1390.0,
                tension,
                bending: false,
                nonlinear: Nonlinear::None,
                t60: [(100.0, 1.0), (1000.0, 1.0)],
                support: None,
                muffle: None,
            },
        )
        .unwrap();
        let table = RadiationTable::shipped().unwrap();
        head.attach_head_radiation(&table, 0.32, &standard_gas(), None)
    };
    // Slack: the (1,1) shape lies below the table's smallest ka, so the
    // receipt's main claim is gone and the drum refuses whole.
    let slack = attach(0.2);
    assert!(
        matches!(&slack, Err(error) if error.contains("(1,1)") && error.contains("the drum refuses")),
        "{slack:?}"
    );
    // Taut (wave speed about 275 m/s): (1,1) sits near ka 3, while the
    // highest resolved shapes pass ka 12 and are refused alone, never clamped.
    let receipt = attach(20_000.0).unwrap();
    assert!(
        receipt.iter().any(|line| line.starts_with("(1,1) cos: ")),
        "{receipt:#?}"
    );
    assert!(
        receipt
            .iter()
            .any(|line| line.starts_with("head radiation (")
                && line.contains("refused for this shape")
                && line.contains("outside the table")
                && line.contains("its T60 is an upper bound")),
        "{receipt:#?}"
    );
    assert!(
        receipt
            .iter()
            .any(|line| line.starts_with("head radiation: ")
                && line.contains("shapes refused by the table, the lowest")),
        "{receipt:#?}"
    );
}
