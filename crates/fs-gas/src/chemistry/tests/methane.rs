use super::*;
use fs_thermochem::kinetics::gri30::{
    CANONICAL_SPECIES, Gri30Calorics, MOLAR_MASSES_KG_PER_MOL, cantera_gri30,
};

// Independent Cantera 3.2.0 CVODES integration of the complete 53-species,
// 325-channel source, including the 16 originally irreversible channels.
// Reproduce with data/reference/thermochem/gri30_cantera32.py. ContinuousAt300K
// changes only upper NASA h/s integration constants, with each species' own
// seam. The oracle applies exactly that declared source treatment independently.
// CH4:O2:N2=1:4:15.04, phi=0.5; 101325 Pa, one litre, NO radical seeds.
// Cantera rtol=1e-11, atol=1e-24, max_dt=duration/100, max_steps=200000.
struct Reference {
    initial_temperature: f64,
    duration: f64,
    temperature: f64,
    pressure: f64,
    energy: f64,
    entropy_change: f64,
    amounts: [f64; 53],
}

const INDUCTION: Reference = Reference {
    initial_temperature: 1300.0,
    duration: 0.002,
    temperature: 1300.5483131534243,
    pressure: 101368.2069938082,
    energy: 177.0971670785941,
    entropy_change: 0.00026760492537425407,
    amounts: [
        0.0,
        3.43115769788306732e-26,
        1.06140310996990990e-20,
        1.42174740898263780e-12,
        5.04991480829283950e-15,
        8.76511920824774114e-09,
        2.94152714996550214e-11,
        1.36629660965665878e-07,
        6.89783482530889528e-16,
        1.18572028214342473e-11,
        8.54303700807712062e-20,
        1.94095727619463079e-13,
        3.00880254629825631e-14,
        1.49003540758122192e-14,
        9.66372990776719041e-12,
        5.47923334353633314e-07,
        1.15776629354559805e-14,
        2.41159576727754582e-07,
        3.83015798298790955e-13,
        6.01723355315999326e-11,
        2.34379715892824635e-09,
        4.66666297607757144e-04,
        1.42737014289160528e-27,
        3.03109341887122088e-08,
        6.80790497505722825e-10,
        2.11343560307675804e-10,
        8.43276119036807679e-08,
        1.28313958737337556e-22,
        7.89829192485988918e-07,
        8.83602852049171541e-10,
        7.87520862654774059e-18,
        9.89186391103302012e-18,
        4.99864454986670983e-19,
        2.82594961731290390e-22,
        1.36241115765274114e-23,
        2.41903894745121819e-12,
        5.80384502272779273e-23,
        2.34408310816820553e-20,
        1.85686425558819954e-08,
        1.53916126800979234e-24,
        2.52751888585911878e-20,
        7.03540647352913617e-03,
        1.04275331537916860e-12,
        5.18843531521263183e-23,
        2.47631597831873221e-20,
        3.28527488919869210e-23,
        2.43150317199862137e-25,
        1.81572061316724083e-16,
        7.82954174293349958e-17,
        1.68944933823977689e-18,
        2.87925132713912357e-10,
        1.87041287382875296e-03,
        7.03854429960079426e-10,
    ],
};

const IGNITION: Reference = Reference {
    initial_temperature: 1500.0,
    duration: 0.002,
    temperature: 2600.310562321808,
    pressure: 177652.8770327451,
    energy: 200.88857379341917,
    entropy_change: 0.16699916354559008,
    amounts: [
        8.08668507371600084e-36,
        1.54721424014234740e-17,
        3.71451441403976596e-25,
        6.44223297554377638e-24,
        7.76973809009003871e-29,
        1.09963801436724108e-29,
        1.31047108852590532e-34,
        1.86464445535167813e-36,
        2.70081209802058159e-50,
        3.55091357038654156e-52,
        8.82394956394892589e-19,
        3.52725041788705959e-19,
        3.09796749881699335e-20,
        1.98404402863889824e-27,
        3.10816835524801180e-22,
        7.83005904888764137e-14,
        5.12328294155041582e-19,
        3.10136874493185609e-19,
        1.37208875759446300e-29,
        1.03465628373742960e-20,
        9.58920623586726957e-21,
        2.78367731766664837e-20,
        3.12906372463854879e-15,
        5.29879997956849617e-05,
        3.52421043959703180e-04,
        1.09012066582410259e-05,
        1.70722097989625436e-05,
        5.67590054093007360e-20,
        7.34103263776898058e-04,
        2.61790783304460633e-09,
        1.60897583888051148e-21,
        9.47804950408468721e-25,
        8.23224298191408146e-14,
        2.32702800938227308e-22,
        3.04261675441333645e-18,
        1.53445332830716471e-11,
        1.37695658083673934e-12,
        6.68243893614902961e-10,
        9.10224761796315466e-08,
        8.97145982662897260e-15,
        2.24859623145309477e-09,
        6.08154285366771519e-03,
        8.91695934542712967e-09,
        3.85116179141695339e-13,
        1.38789634637244048e-10,
        1.03621277669272131e-11,
        4.52906841404481967e-12,
        4.51776751019591690e-11,
        3.15751131235300357e-05,
        2.27388144947452572e-08,
        4.47102555781983522e-05,
        7.83262239372465381e-04,
        1.08286928545102407e-04,
    ],
};

fn methane_controls(relative: f64) -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        method: IsochoricKineticsMethod::Sdirk2,
        minimum_temperature_k: 300.0,
        maximum_temperature_k: 3000.0,
        energy_tolerance_j: 1e-8,
        relative_tolerance: relative,
        absolute_amount_tolerance_mol: 1e-16,
        initial_step_s: 1e-9,
        maximum_step_s: 2e-4,
        maximum_steps: 2000,
        ..IsochoricKineticsSolve::default()
    }
}

fn methane_run(case: &Reference, relative: f64) -> IsochoricKineticsResult {
    let mechanism = cantera_gri30(Gri30Calorics::ContinuousAt300K).unwrap();
    let unit = 101325.0 * 0.001 / (R * case.initial_temperature * 20.04);
    advance_adiabatic(
        &mechanism,
        &[
            (id("CH4"), Amount::new(unit)),
            (id("O2"), Amount::new(4.0 * unit)),
            (id("N2"), Amount::new(15.04 * unit)),
        ],
        Volume::new(0.001),
        Temperature::new(case.initial_temperature),
        Time::new(case.duration),
        methane_controls(relative),
        || false,
    )
    .unwrap()
}

// Independent original source element counts; deliberately not read back from
// the mechanism's stoichiometric matrix. Canonical rows C,H,N,O; AR remains zero.
const ATOMS: [[u8; 53]; 4] = [
    [
        0, 1, 2, 2, 2, 2, 2, 2, 3, 3, 1, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0,
        2, 2, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    ],
    [
        0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 2, 3, 2, 2, 3, 3, 4, 3, 4, 4, 0, 0, 0, 1, 2, 2, 2, 2,
        1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 2, 3, 1, 0, 0, 0, 0, 1,
    ],
    [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 1, 2, 1, 0, 1, 1, 0, 1, 1, 2, 2, 1, 1, 1, 1, 2, 1, 1, 0, 0, 0,
    ],
    [
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 2, 0, 0, 0, 1, 2,
        1, 1, 0, 0, 1, 1, 1, 1, 2, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 2, 1, 2, 1,
    ],
];

fn check_endpoint(actual: &IsochoricKineticsResult, case: &Reference) {
    let t_error = (actual.temperature().value() - case.temperature).abs();
    assert!(t_error < 1.0, "T error={t_error}, {:?}", actual.balance());
    assert!((actual.pressure().value() / case.pressure - 1.0).abs() < 0.001);
    for ((name, a), b) in CANONICAL_SPECIES
        .iter()
        .zip(actual.amounts())
        .zip(case.amounts)
    {
        assert!(
            (a.value() - b).abs() <= 1e-12 + 0.005 * b.abs(),
            "{name}: actual={}, Cantera={b}, {:?}",
            a.value(),
            actual.balance(),
        );
    }
    assert_eq!(actual.elapsed_time(), Time::new(case.duration));
    assert!(actual.amounts().iter().all(|n| n.value() >= 0.0));
    assert_eq!(actual.amounts()[0].value(), 0.0);
    let unit = 101325.0 * 0.001 / (R * case.initial_temperature * 20.04);
    for (atoms, ratio) in ATOMS.iter().zip([1.0, 4.0, 30.08, 8.0]) {
        let inventory: f64 = actual
            .amounts()
            .iter()
            .zip(atoms)
            .map(|(n, count)| n.value() * f64::from(*count))
            .sum();
        assert!((inventory - ratio * unit).abs() < 1e-14);
    }
    let mass: f64 = actual
        .amounts()
        .iter()
        .zip(MOLAR_MASSES_KG_PER_MOL)
        .map(|(n, m)| n.value() * m)
        .sum();
    let initial_mass = unit * (0.016043 + 4.0 * 0.031998 + 15.04 * 0.028014);
    assert!((mass - initial_mass).abs() < 1e-16);
    assert!((actual.internal_energy().value() - case.energy).abs() < 2e-8);
    assert!(actual.balance().energy_residual_j.abs() <= methane_controls(2e-5).energy_tolerance_j);
    assert!(actual.balance().entropy_change_j_k > 0.0);
    assert!((actual.balance().entropy_change_j_k - case.entropy_change).abs() < 2e-6);
    assert!(actual.balance().maximum_estimated_local_error <= 1.0);
    eprintln!(
        "GRI methane initial={} T error={t_error:e}, {:?}",
        case.initial_temperature,
        actual.balance()
    );
}

#[test]
fn full_gri_methane_induction_and_ignition_match_independent_cantera() {
    for case in [&INDUCTION, &IGNITION] {
        let result = methane_run(case, 2e-5);
        check_endpoint(&result, case);
        if case.initial_temperature == 1500.0 {
            assert!(result.balance().direct_fallbacks > 0);
        }
    }
}
