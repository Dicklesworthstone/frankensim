use super::*;
use crate::kinetics::KineticWorkspaceV1;

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn complete_source_axis_original_calorics_and_channel_identity() {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let models = mechanism.thermodynamics().models();
    assert_eq!(models.len(), 10);
    assert_eq!(mechanism.reactions().len(), 29);
    for (i, model) in models.iter().enumerate() {
        assert_eq!(
            model.species(),
            &SpeciesId::new(CANONICAL_SPECIES[i]).unwrap()
        );
        assert_eq!(model.molar_mass().value(), MOLAR_MASSES_KG_PER_MOL[i]);
        assert!(model.caloric_reference().is_none());
        for (case, t) in [700.0, 1800.0].into_iter().enumerate() {
            let p = model.evaluate_properties(Temperature::new(t)).unwrap();
            for (actual, expected) in [p.cp().value(), p.h().value(), p.s().value()]
                .into_iter()
                .zip(THERMO[case][i])
            {
                close(actual, expected, 5e-13, 1e-9);
            }
        }
    }
    assert!(
        models[0]
            .evaluate_properties(Temperature::new(299.0))
            .is_err()
    );
    assert!(
        models[2]
            .evaluate_properties(Temperature::new(3501.0))
            .is_err()
    );
    for (i, reaction) in mechanism.reactions().iter().enumerate() {
        assert_eq!(
            reaction.reaction(),
            &ReactionId::new(format!("cantera-h2o2-r{:02}", i + 1)).unwrap()
        );
    }
    // Source duplicate channels remain separate kinetic contributions.
    for (a, b) in [(23, 28), (24, 25), (26, 27)] {
        assert_eq!(
            mechanism.reactions()[a].reactants(),
            mechanism.reactions()[b].reactants()
        );
        assert_eq!(
            mechanism.reactions()[a].products(),
            mechanism.reactions()[b].products()
        );
        assert_ne!(
            mechanism.reactions()[a].forward_rate(),
            mechanism.reactions()[b].forward_rate()
        );
    }
    // Explicit collision partners are present on both sides, including an
    // extra O2 molecule in channel 7; no second generic-M multiplier is added.
    assert_eq!(
        mechanism.reactions()[6]
            .reactants()
            .iter()
            .find(|(s, _)| s == &SpeciesId::new("O2").unwrap())
            .unwrap()
            .1,
        2
    );
    for i in [6, 7, 8, 9, 12, 13] {
        assert_eq!(
            mechanism.reactions()[i].pressure_dependence(),
            &PressureDependenceV1::Independent
        );
    }
}

#[test]
fn all_29_rates_match_independent_original_yaml_nasa7_and_cgs_algebra() {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let mut work = KineticWorkspaceV1::new(&mechanism).unwrap();
    let c = [1.0, 0.02, 4.0, 0.5, 0.05, 0.03, 6.0, 0.01, 2.0, 0.04];
    for (case, t) in [700.0, 1800.0].into_iter().enumerate() {
        let result = mechanism
            .evaluate(Temperature::new(t), &c, &mut work, || false)
            .unwrap();
        for (j, (&actual, &expected)) in work
            .reaction_rates()
            .iter()
            .zip(&NET_RATES[case])
            .enumerate()
        {
            assert!(
                (actual - expected).abs() <= 2e-11 * expected.abs() + 2e-10,
                "source channel {}, T={t}: actual={actual:.17e}, expected={expected:.17e}",
                j + 1
            );
        }
        close(
            result.thermodynamics.pressure.value(),
            c.iter().sum::<f64>() * R * t,
            1e-14,
            1e-9,
        );
        assert!(result.entropy_production_w_m3_k >= 0.0);
        let mass_rate: f64 = work
            .species_rates()
            .iter()
            .zip(MOLAR_MASSES_KG_PER_MOL)
            .map(|(rate, mass)| rate * mass)
            .sum();
        let absolute_mass_rate: f64 = work
            .species_rates()
            .iter()
            .zip(MOLAR_MASSES_KG_PER_MOL)
            .map(|(rate, mass)| (rate * mass).abs())
            .sum();
        assert!(mass_rate.abs() <= 2e-14 * absolute_mass_rate);
        let elements = mechanism.thermodynamics().elemental_matrix();
        for row in 0..elements.row_count() {
            let residual: f64 = work
                .species_rates()
                .iter()
                .enumerate()
                .map(|(s, rate)| elements.get(row, s).unwrap() as f64 * rate)
                .sum();
            let scale: f64 = work
                .species_rates()
                .iter()
                .enumerate()
                .map(|(s, rate)| (elements.get(row, s).unwrap() as f64 * rate).abs())
                .sum();
            assert!(residual.abs() <= 2e-14 * scale);
        }
    }
}

#[test]
fn continuity_is_explicit_bounded_and_preserves_original_low_temperature_references() {
    let raw = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let adjusted = cantera_h2o2(HydrogenCalorics::ContinuousAt300K).unwrap();
    assert!(
        raw.admit_temperature_interval(Temperature::new(300.0), Temperature::new(3500.0), || false)
            .is_err()
    );
    adjusted
        .admit_temperature_interval(Temperature::new(300.0), Temperature::new(3500.0), || false)
        .unwrap();
    for (original, continuous) in raw
        .thermodynamics()
        .models()
        .iter()
        .zip(adjusted.thermodynamics().models())
    {
        let reference = continuous.caloric_reference().unwrap();
        assert_eq!(reference.anchor_temperature().value(), 300.0);
        assert_eq!(reference.regions()[0].enthalpy_offset().value(), 0.0);
        assert_eq!(reference.regions()[0].entropy_offset().value(), 0.0);
        assert!(reference.regions()[1].enthalpy_offset().value().abs() <= 0.01);
        assert!(reference.regions()[1].entropy_offset().value().abs() <= 2e-5);
        for t in [300.0, 700.0, 1800.0, 3500.0] {
            let a = original.evaluate_properties(Temperature::new(t)).unwrap();
            let b = continuous.evaluate_properties(Temperature::new(t)).unwrap();
            assert_eq!(a.cp().value().to_bits(), b.cp().value().to_bits());
            if t < 1000.0 {
                assert_eq!(a, b);
            }
        }
    }
    let n2 = adjusted.thermodynamics().models()[6]
        .caloric_reference()
        .unwrap();
    close(
        n2.regions()[1].enthalpy_offset().value(),
        0.005198202023166232,
        0.0,
        1e-9,
    );
    close(
        n2.regions()[1].entropy_offset().value(),
        -1.5006413178753064e-5,
        0.0,
        1e-12,
    );
}

#[test]
fn source_equilibrium_has_detailed_balance_in_every_duplicate_and_pressure_channel() {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let t = Temperature::new(1500.0);
    let bank = mechanism.thermodynamics();
    // An independent chemical-potential equilibrium construction: mu_i/RT is
    // a linear combination of the four conserved element counts. Every exact
    // elemental stoichiometric column therefore has zero reaction affinity.
    let c: Vec<f64> = bank
        .models()
        .iter()
        .enumerate()
        .map(|(s, model)| {
            let atoms: i128 = (0..bank.elemental_matrix().row_count())
                .map(|e| bank.elemental_matrix().get(e, s).unwrap())
                .sum();
            let g = model.evaluate_properties(t).unwrap().g().value();
            101325.0 / (R * t.value()) * (-g / (R * t.value()) - 20.0 * atoms as f64).exp()
        })
        .collect();
    let mut work = KineticWorkspaceV1::new(&mechanism).unwrap();
    let evaluation = mechanism.evaluate(t, &c, &mut work, || false).unwrap();
    for affinity in work.dimensionless_affinities() {
        assert!(affinity.unwrap().abs() < 2e-13);
    }
    assert!(evaluation.entropy_production_w_m3_k >= 0.0);
}

#[test]
fn pure_inert_and_cancelled_reuse_do_not_create_radicals_or_drop_channels() {
    let mechanism = cantera_h2o2(HydrogenCalorics::ContinuousAt300K).unwrap();
    let mut work = KineticWorkspaceV1::new(&mechanism).unwrap();
    let c = [40.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert_eq!(
        mechanism.evaluate(Temperature::new(300.0), &c, &mut work, || true),
        Err(KineticsErrorV1::Cancelled)
    );
    let report = mechanism
        .evaluate(Temperature::new(300.0), &c, &mut work, || false)
        .unwrap();
    assert_eq!(work.reaction_rates(), &[0.0; 29]);
    assert_eq!(work.species_rates(), &[0.0; 10]);
    assert_eq!(report.entropy_production_w_m3_k, 0.0);
    let first = work.reaction_rates().to_vec();
    mechanism
        .evaluate(Temperature::new(300.0), &c, &mut work, || false)
        .unwrap();
    assert_eq!(first, work.reaction_rates());
}

// Independent development oracle: parse the pinned source YAML, evaluate its
// NASA7 equations directly with Python binary64/libm, then use its original
// cgs Arrhenius laws, stated collider efficiencies and standard Troe expression.
// Kc = exp(-sum(nu*g)/(RT)) * (101325/(RT))^sum(nu), c in mol/m3.
// This is source-algebra conformance, not a Cantera runtime or experiment claim.
#[rustfmt::skip]
const NET_RATES: [[f64; 29]; 2] = [
    [4.50514285714285700e-1, 2.83571428571428630, 7.95991490887316786e2, 5.99999999999995271e3, 1.33032444453160707e2,
     1.66142906885145203e3, 4.93437046002690693e2, 1.54986501664284492e3, 1.85038892251008997e3, 1.48278983145155536e2,
     -5.64312126992563026e3, 5.01714285714285779, 2.82685588275704403, 3.33279572882393804, 4.76816326530612287e2,
     1.47045018737716237e3, 1.24735859504202344e4, 3.19285045283912077e4, 1.41059899006771190e2, 7.51712697457856393e2,
     5.80492753593640882e4, 4.55771039330533767e2, 1.75275941922677430e3, 2.49260322728095343e4, 2.94270858475868954e3,
     2.23587356833318429, 3.77645740549834386e2, 6.77638522016085716e1, 2.33125571262920310e1],
    [1.74849139431163869e-1, 1.10262947865424210, 1.62997443057739583e5, 5.99928455893260252e3, 5.06397197414314087e3,
     -9.42380293161550071e1, -1.95483383028249875e1, -9.66174815690789046e1, -7.33062686355885944e1, -8.90091098730545127,
     5.53079573170053118e4, 1.93415622304675794, 1.59004178110828231, 1.01461239121542368, 7.21085275701627637e1,
     1.97456728526645497e3, 1.97878979328672394e4, 4.22016153704374155e4, 5.09516638872088515e3, 3.65518245536259246e3,
     1.08590153017486073e6, -1.98033423964696209e5, 4.90992757355752838e3, 2.00097897602739322e4, 3.54340395897729513e3,
     9.11786877104279236e5, 1.81336658688996039e2, 1.29693051464013006e4, 4.72095941622185346e4],
];

// Cp [J/(mol K)], h [J/mol], s [J/(mol K)] in the canonical species order.
#[rustfmt::skip]
const THERMO: [[[f64; 3]; 10]; 2] = [
    [[20.7861565453831, 8352.917007762198, 172.47272517645],
     [20.786156546066245, 226350.10256829578, 132.45785350391796],
     [29.344387908470154, 11750.979185037055, 155.62195271118554],
     [37.508301539904316, -227633.0043352295, 218.73271751908655],
     [56.250668612580434, -115856.62184438464, 276.1953573408057],
     [43.644636946364244, 28408.76692500034, 262.23524199901607],
     [30.684496851242105, 11943.85136142621, 216.78555941427757],
     [21.061931709503046, 257742.64133031195, 179.30723472156166],
     [32.99224827897792, 12499.679024836601, 231.47078425047835],
     [29.684113135655444, 51234.06107202984, 208.99814453887325]],
    [[20.7861565453831, 31217.689207683605, 192.10445202892035],
     [20.786156532153157, 249214.87477998526, 152.08958027712606],
     [33.503822066261804, 46164.23110437554, 184.85849374010894],
     [50.148446887194275, -178981.87749643147, 259.5470596180469],
     [72.35232482979043, -43603.36688918289, 337.03382720320184],
     [54.450097795224536, 83132.66641318776, 308.4266570946896],
     [35.59505654629914, 48972.38607571014, 248.21390508566842],
     [20.828155218676077, 280722.5120234869, 199.05232620989327],
     [37.335379424268915, 51691.50402752842, 264.8124604367768],
     [34.091715693385694, 86252.5968532739, 238.7208904626411]],
];
