//! Original Maxwell source to passive series-creep realization.
use fs_material::fiber::WoolFelt;
use fs_material::fiber::felt_rate::{MAX_FELT_RATE_BRANCHES, RateFelt};
fn envelope() -> WoolFelt {
    WoolFelt::new(30000.0, 0.2, 2.2, 3.0, 0.15, 0.7).unwrap()
}
fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-100),
        "{a:e} != {b:e}, relative tolerance {relative:e}"
    );
}
#[test]
fn g1_prony_single_arm_matches_analytic_retardation_including_weak_loss() {
    use fs_material::visco::GeneralizedMaxwell;
    for (equilibrium, branch, tau) in [
        (1e6, 2e6, 0.02),
        (1e6, 1e-12, 0.2),
        (2e4, 4e5, 0.0001),
        (1.0, 1e20, 1.0),
    ] {
        let card = GeneralizedMaxwell::new(equilibrium, vec![(branch, tau)]).unwrap();
        let law = RateFelt::from_prony(envelope(), &card).unwrap();
        assert_eq!(law.creep().len(), 1);
        let instant = equilibrium + branch;
        close(
            law.creep()[0].modulus_pa,
            equilibrium * (instant / branch),
            2e-12,
        );
        close(
            law.creep()[0].viscosity_pa_s,
            instant * (instant / branch) * tau,
            2e-12,
        );
        close(
            law.creep()[0].viscosity_pa_s / law.creep()[0].modulus_pa,
            tau * (instant / equilibrium),
            2e-12,
        );
    }
}

#[test]
fn g1_prony_multi_arm_preserves_original_complex_compliance() {
    use fs_material::visco::GeneralizedMaxwell;
    let cards = [
        GeneralizedMaxwell::new(
            2e6,
            vec![(8e5, 0.0001), (3e6, 0.02), (1e6, 0.7), (4e5, 4.0)],
        )
        .unwrap(),
        GeneralizedMaxwell::new(2e6, vec![(1e6, 0.03), (2e6, 0.03), (4e5, 0.8)]).unwrap(),
        GeneralizedMaxwell::new(2e6, vec![(5e-8, 0.0004), (2e4, 0.03)]).unwrap(),
    ];
    for card in cards {
        let law = RateFelt::from_prony(envelope(), &card).unwrap();
        let instant = card.e_inf + card.terms.iter().map(|x| x.0).sum::<f64>();
        for omega in [0.0, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 1e6] {
            // Independent reciprocal of the ORIGINAL authored Maxwell
            // modulus, compared to the realized spring-plus-Kelvin network.
            let (storage, loss) = card.modulus(omega);
            let denominator = storage * storage + loss * loss;
            let real = storage / denominator;
            let imag = -loss / denominator;
            let mut realized_real = 1.0 / instant;
            let mut realized_imag = 0.0;
            for arm in law.creep() {
                let wt = omega * (arm.viscosity_pa_s / arm.modulus_pa);
                let response = 1.0 / (arm.modulus_pa * (1.0 + wt * wt));
                realized_real += response;
                realized_imag -= response * wt;
            }
            close(realized_real, real, 1e-12);
            if omega > 0.0 {
                close(realized_imag, imag, 1e-12);
            } else {
                assert_eq!(realized_imag, 0.0);
            }
        }
    }
}

#[test]
fn g0_prony_admission_keeps_empty_and_zero_arms_and_refuses_invalid_cards() {
    use fs_material::visco::GeneralizedMaxwell;
    for terms in [vec![], vec![(0.0, 0.01)]] {
        let card = GeneralizedMaxwell::new(2e6, terms).unwrap();
        assert!(
            RateFelt::from_prony(envelope(), &card)
                .unwrap()
                .creep()
                .is_empty()
        );
    }
    let too_many =
        GeneralizedMaxwell::new(2e6, vec![(1e5, 0.01); MAX_FELT_RATE_BRANCHES + 1]).unwrap();
    assert!(RateFelt::from_prony(envelope(), &too_many).is_err());
    let invalid = GeneralizedMaxwell {
        e_inf: 2e6,
        terms: vec![(-1.0, 0.1)],
    };
    assert!(RateFelt::from_prony(envelope(), &invalid).is_err());
}

#[test]
fn g1_prony_extremely_weak_mixed_arm_is_accurate_or_explicitly_unresolved() {
    use fs_material::visco::GeneralizedMaxwell;
    let card = GeneralizedMaxwell::new(1.0, vec![(1.0, 1.0), (1e-40, 4.0)]).unwrap();
    match RateFelt::from_prony(envelope(), &card) {
        Ok(law) => {
            // Original rational compliance residue at the weak material pole.
            // As E_weak -> 0, J_weak/E_weak = 9/4 for this source card.
            // O(E_weak^2) corrections are far below the asserted relative gate.
            let weak = law
                .creep()
                .iter()
                .min_by(|a, b| (1.0 / a.modulus_pa).total_cmp(&(1.0 / b.modulus_pa)))
                .unwrap();
            close(1.0 / weak.modulus_pa, 2.25e-40, 1e-9);
            close(weak.viscosity_pa_s / weak.modulus_pa, 4.0, 1e-9);
        }
        Err(error) => {
            // A source-relative numerical refusal is acceptable. A silently
            // inaccurate retained weight is not a passive-source realization.
            let reason = error.to_string();
            assert!(
                reason.contains("unresolved")
                    || reason.contains("component")
                    || reason.contains("relative"),
                "{reason}"
            );
        }
    }
}
