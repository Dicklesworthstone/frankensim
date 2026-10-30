//! The film's generalized-Maxwell law: its field recursion against
//! fs-material's `GeneralizedMaxwell::step`, its fit and certificate, the
//! glassy admission, and the stepped head's exact discrete energy account
//! with air mass, radiation and the Berger term.

use super::*;
use crate::film::{FilmField, FilmLaw};
use fs_material::visco::GeneralizedMaxwell;

/// A deterministic history in [-1, 1] (no clock or entropy).
fn history(i: usize, node: usize) -> f64 {
    let x = ((i * 7919 + node * 104_729) % 10_007) as f64 / 10_007.0;
    2.0 * x - 1.0
}

/// Three branches whose relaxation times straddle the clock: slow (dt/tau ~ 7e-4),
/// mid (0.13) and faster than the step (2.8).
fn three_branch_law() -> FilmLaw {
    let tau = |hz: f64| 1.0 / (2.0 * PI * hz);
    FilmLaw {
        law: GeneralizedMaxwell::new(
            0.9,
            vec![
                (0.02, tau(10.0)),
                (0.01, tau(2000.0)),
                (0.05, tau(43_000.0)),
            ],
        )
        .unwrap(),
        band_hz: (20.0, 20_000.0),
        anchor_hz: 632.0,
        ladder: 3,
        sup_rel_tan_delta_err: 0.0,
        verification_points: 0,
    }
}

#[test]
fn g0_the_field_recursion_is_fs_material_generalized_maxwell_step() {
    // With K = I the branch fields are the branch strains themselves: each node must follow
    // GeneralizedMaxwell::step bit for bit, and the force must be minus its stress.
    let law = three_branch_law();
    let (n, k) = (5, 1.0 / 96_000.0);
    let mut field = FilmField::new(law.clone(), n);
    let mut states: Vec<_> = (0..n).map(|_| law.law.state()).collect();
    let mut wp = vec![0.0; n];
    for i in 0..4000 {
        let w: Vec<f64> = (0..n).map(|node| 1e-4 * history(i, node)).collect();
        let mut f: Vec<f64> = w.iter().map(|x| -x).collect(); // linear_forces: -K w
        let prev = (i == 0).then_some(wp.as_slice());
        field.step_with(k, &w, prev, &w, &wp, &mut f);
        assert!(
            field.last_j >= -1e-30,
            "a step's heat is a square: {}",
            field.last_j
        );
        for node in 0..n {
            let stress = law.law.step(&mut states[node], w[node], k);
            for j in 0..3 {
                assert_eq!(
                    field.branch(j)[node],
                    states[node].q[j],
                    "step {i}, node {node}, branch {j}"
                );
            }
            assert!(
                (f[node] + stress).abs() <= 4e-15 * (w[node].abs() + stress.abs()),
                "step {i}: {} vs {stress}",
                f[node]
            );
        }
        wp = w;
    }
    // On a resolved history (200 Hz at 96 kHz) the midpoint heat is the exact ramp heat.
    let mut field = FilmField::new(law.clone(), 1);
    let mut state = law.law.state();
    let mut wp = 0.0;
    for i in 1..=1440 {
        let w = 1e-3 * (2.0 * PI * 200.0 * i as f64 * k).sin();
        let mut f = [-w];
        let prev = [wp];
        field.step_with(
            k,
            &[w],
            (i == 1).then_some(prev.as_slice()),
            &[w],
            &[wp],
            &mut f,
        );
        law.law.step(&mut state, w, k);
        wp = w;
    }
    assert!(
        (field.heat_j / state.dissipated - 1.0).abs() < 2e-4,
        "{} vs {}",
        field.heat_j,
        state.dissipated
    );
}

/// A smooth stand-in for the card's curve (manufactured, not measured): 0.016 at 20 Hz rising
/// linearly in log frequency to 0.037 at 20 kHz.
fn manufactured_curve(hz: f64) -> Result<Option<f64>, String> {
    Ok((20.0..=20_000.0)
        .contains(&hz)
        .then(|| 0.016 + 0.007 * (hz / 20.0).log10()))
}

#[test]
fn g0_the_fit_meets_its_certificate_and_its_anchor() {
    let law = crate::film::fit(
        &manufactured_curve,
        (20.0, 20_000.0),
        crate::film::matsuzawa_anchor_hz(),
    )
    .unwrap();
    assert!(
        law.sup_rel_tan_delta_err < 0.01,
        "{}",
        law.sup_rel_tan_delta_err
    );
    assert_eq!(law.ladder, 16);
    assert!(law.law.terms.len() >= 4 && law.law.terms.iter().all(|t| t.0 > 0.0));
    for hz in [20.0, 63.0, 632.0, 5000.0, 20_000.0] {
        let curve = manufactured_curve(hz).unwrap().unwrap();
        assert!((law.tan_delta(hz) / curve - 1.0).abs() <= law.sup_rel_tan_delta_err + 1e-12);
    }
    assert!((law.storage(law.anchor_hz) - 1.0).abs() < 1e-12);
    // Causality: the storage modulus rises through the band, from the relaxed to the glassy.
    assert!(law.law.e_inf < law.storage(20.0) && law.storage(20.0) < 1.0);
    assert!(1.0 < law.storage(20_000.0) && law.storage(20_000.0) < law.glassy());
    assert!(
        law.receipt().contains("never to a recording")
            && law.receipt().contains("Matsuzawa 1959's 100-4000 c/s")
    );
    // Refusals: an anchor outside the band, a curve that is silent inside it.
    assert!(crate::film::fit(&manufactured_curve, (20.0, 20_000.0), 30_000.0).is_err());
    assert!(crate::film::fit(&|_| Ok(None), (20.0, 20_000.0), 632.0).is_err());
}

fn film_head(grid: f64, tension: f64, nonlinear: Nonlinear) -> Body {
    let (thickness, rest) = (|_: f64| 1.9e-4, |_: f64| 0.0);
    Body::plate_without_surface_loss(
        "film head",
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
            tension,
            bending: true,
            nonlinear,
            t60: [(100.0, 1.0), (1000.0, 1.0)],
            support: None,
            muffle: None,
        },
    )
    .unwrap()
}

#[test]
fn g0_admission_takes_the_film_law_at_its_glassy_modulus() {
    // Without tension the stiffness is all bending, so the law's highest rate is the
    // elastic one times the glassy modulus.
    let law = crate::film::fit(
        &manufactured_curve,
        (20.0, 20_000.0),
        crate::film::matsuzawa_anchor_hz(),
    )
    .unwrap();
    let mut plate = film_head(0.02, 0.0, Nonlinear::None);
    let elastic = plate.max_rate2();
    plate.attach_film_maxwell(law.clone()).unwrap();
    let glassy = plate.max_rate2();
    assert!(
        (glassy / elastic / law.glassy() - 1.0).abs() < 1e-9,
        "{glassy} / {elastic} against {}",
        law.glassy()
    );
    // Refusals: a second law, a plate whose fitted T60 surface loss already lumps the film's,
    // a stiffness operator that also holds support springs, a body without bending.
    let second = plate.attach_film_maxwell(law.clone()).unwrap_err();
    assert!(second.contains("attaches once"), "{second}");
    let (thickness, rest) = (|_: f64| 1.9e-4, |_: f64| 0.0);
    let mut fitted = Body::plate(
        "fitted",
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
            tension: 100.0,
            bending: true,
            nonlinear: Nonlinear::None,
            t60: [(100.0, 1.0), (1000.0, 1.0)],
            support: None,
            muffle: None,
        },
    )
    .unwrap();
    let lumped = fitted.attach_film_maxwell(law.clone()).unwrap_err();
    assert!(lumped.contains("fitted T60"), "{lumped}");
    let mut supported = Body::plate_without_surface_loss(
        "supported",
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
            tension: 100.0,
            bending: true,
            nonlinear: Nonlinear::None,
            t60: [(100.0, 1.0), (1000.0, 1.0)],
            support: Some((0.28, 0.32, 1e6, 0.0)),
            muffle: None,
        },
    )
    .unwrap();
    assert!(
        supported
            .attach_film_maxwell(law.clone())
            .unwrap_err()
            .contains("support springs")
    );
    let mut wire = Body::string("wire", 8, 0.01, 1e-3, 10.0, 1.0, 0.0, 0.0);
    assert!(wire.attach_film_maxwell(law.clone()).is_err());
    // The law comes before the head radiation, so that the load is built knowing the head's
    // own loss is modelled and retains every shape the grid resolves.
    let mut loaded = film_head(0.02, 5360.0, Nonlinear::None);
    let table = crate::radiation::RadiationTable::shipped().unwrap();
    loaded
        .attach_head_radiation(&table, 0.32, &standard_gas(), None)
        .unwrap();
    let refused = loaded.attach_film_maxwell(law).unwrap_err();
    assert!(refused.contains("before the head radiation"), "{refused}");
}

#[test]
fn g0_the_film_law_heat_closes_the_discrete_energy_with_air_and_berger() {
    // Berger stretching (elastic, through the auxiliary variable), the implicit air mass and
    // the explicit radiation resistance share the step with the film law's branches.
    let law = crate::film::fit(
        &manufactured_curve,
        (20.0, 20_000.0),
        crate::film::matsuzawa_anchor_hz(),
    )
    .unwrap();
    let mut head = film_head(0.02, 5360.0, Nonlinear::Berger);
    let lines = head.attach_film_maxwell(law).unwrap();
    assert!(
        lines
            .iter()
            .any(|l| l.starts_with("film law on the bending operator: every grid mode carries it"))
    );
    let table = crate::radiation::RadiationTable::shipped().unwrap();
    // The film law alone (no viscous layer) is the head's own loss: the load retains every
    // shape the grid resolves, and its scope says the film is not in it.
    let receipt = head
        .attach_head_radiation(&table, 0.32, &standard_gas(), None)
        .unwrap();
    assert!(
        receipt
            .iter()
            .any(|l| l.contains("extended to every shape the grid resolves"))
    );
    let scope = receipt
        .iter()
        .find(|l| l.starts_with("head radiation scope:"))
        .unwrap();
    assert!(
        scope
            .contains("the film law is this head's own, on its bending operator, not in this load"),
        "{scope}"
    );
    let shapes = [
        crate::radiation::ModeShape::new(1, 1),
        crate::radiation::ModeShape::new(4, 3),
    ];
    for (i, p) in head.xy.clone().iter().enumerate() {
        let (x, y) = (p[0] / 0.32, p[1] / 0.32);
        head.w[i] = 1e-3 * (shapes[0].value(x, y, false) + 0.5 * shapes[1].value(x, y, true));
        head.wp[i] = 0.9 * head.w[i];
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
            "the film and the air only remove energy"
        );
        assert!(
            (energy + body.head_radiated_j() + body.film_heat_j() - initial).abs() < 1e-9 * initial,
            "{energy} + {} + {} against {initial}",
            body.head_radiated_j(),
            body.film_heat_j()
        );
        previous = energy;
    }
    assert!(ins.bodies[0].film_heat_j() > 1e-3 * initial);
    // A clone steps exactly as the original (the branches are state).
    let mut copy = ins.clone();
    for _ in 0..64 {
        ins.step(&mut []);
        copy.step(&mut []);
    }
    assert_eq!(ins.bodies[0].w, copy.bodies[0].w);
    assert_eq!(ins.bodies[0].film_heat_j(), copy.bodies[0].film_heat_j());
}
