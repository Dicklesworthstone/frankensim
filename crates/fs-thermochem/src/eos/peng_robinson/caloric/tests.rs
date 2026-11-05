use super::super::{PengRobinsonComponentV1, PengRobinsonPhaseV1};
use super::*;
use fs_qty::Pressure;

fn components() -> [PengRobinsonComponentV1; 2] {
    [
        PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6),
            acentric_factor: 0.01142,
            molar_mass: MolarMass::new(0.028),
        },
        PengRobinsonComponentV1 {
            critical_temperature: Temperature::new(150.687),
            critical_pressure: Pressure::new(4.863e6),
            acentric_factor: -0.00219,
            molar_mass: MolarMass::new(0.004),
        },
    ]
}
fn model() -> PengRobinsonMixtureV1 {
    PengRobinsonMixtureV1::new(
        components().to_vec(),
        vec![0.75, 0.25],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap()
}
const RHO: f64 = 47.889_003_167_633_395;

fn residuals(change: &PengRobinsonCaloricChangeV1) -> [f64; 3] {
    [
        change.residual_internal_energy_change_j_mol(),
        change.residual_enthalpy_change_j_mol(),
        change.residual_entropy_at_density_change_j_mol_k(),
    ]
}

#[test]
fn finite_calorics_match_independent_high_precision_changes_above_coordinate_resolution() {
    // Independent 85-digit Decimal evaluation of the ORIGINAL dimensional PR
    // equation and logarithmic attractive integral, using exact f64 inputs.
    // These pin differences in J/mol and J/(mol K), not endpoint subtraction
    // from the production evaluator. The density-only step spans 6740 rho
    // ULPs and the thermal step spans 175922 T ULPs. Rounded endpoint
    // subtraction measurably loses departure precision on these weak steps.
    let references = [
        (
            300.000_000_01,
            RHO,
            [
                5.232_029_349_116_007_4e-9,
                1.642_113_721_392_722_6e-8,
                1.744_009_783_009_602_5e-11,
            ],
        ),
        (
            300.0,
            47.889_003_167_681_29,
            [
                -5.330_510_331_974_845_2e-10,
                -6.901_644_689_819_363_8e-10,
                -1.118_946_750_134_369_3e-12,
            ],
        ),
        (
            300.000_000_01,
            47.889_003_172_422_3,
            [
                -4.807_133_404_069_961e-8,
                -5.259_305_691_970_866_3e-8,
                -9.445_092_483_308_131_6e-11,
            ],
        ),
        (
            360.0,
            80.0,
            [
                -298.286_923_042_961_8,
                -281.769_595_727_296_29,
                -0.596_352_174_402_559_05,
            ],
        ),
        // The two components have different alpha zeros; this endpoint
        // crosses one component's sign and keeps its original derivative jump.
        (
            2200.0,
            50.0,
            [
                537.554_802_989_760_56,
                1840.638_700_738_443_5,
                0.639_795_576_887_383_09,
            ],
        ),
    ];
    let model = model();
    for (temperature, density, expected) in references {
        let reference = model
            .prepare_density_caloric_change(
                Temperature::new(300.0),
                Density::new(RHO),
                Density::new(density),
            )
            .unwrap();
        let change = reference.evaluate(Temperature::new(temperature)).unwrap();
        assert_eq!(
            *change.state(),
            model
                .evaluate_density_caloric(Temperature::new(temperature), Density::new(density))
                .unwrap()
        );
        for (actual, expected) in residuals(&change).into_iter().zip(expected) {
            assert!(
                (actual / expected - 1.0).abs() < 8e-13,
                "T={temperature}, rho={density}, actual={actual:e}, expected={expected:e}"
            );
        }
    }
}

#[test]
fn caloric_changes_reverse_exactly_and_keep_the_full_component_bound() {
    let maximum = PengRobinsonMixtureV1::new(
        vec![components()[0]; MAX_PENG_ROBINSON_COMPONENTS_V1],
        vec![1.0 / MAX_PENG_ROBINSON_COMPONENTS_V1 as f64; MAX_PENG_ROBINSON_COMPONENTS_V1],
        vec![],
    )
    .unwrap();
    for model in [model(), maximum] {
        for (t1, rho1) in [
            (300.000_000_01, RHO),
            (300.0, RHO * (1.0 + 1e-12)),
            (360.0, 80.0),
            (2200.0, 50.0),
            (4000.0, 30.0),
        ] {
            let forward = model
                .prepare_density_caloric_change(
                    Temperature::new(300.0),
                    Density::new(RHO),
                    Density::new(rho1),
                )
                .unwrap()
                .evaluate(Temperature::new(t1))
                .unwrap();
            let reverse = model
                .prepare_density_caloric_change(
                    Temperature::new(t1),
                    Density::new(rho1),
                    Density::new(RHO),
                )
                .unwrap()
                .evaluate(Temperature::new(300.0))
                .unwrap();
            for (forward, reverse) in residuals(&forward).into_iter().zip(residuals(&reverse)) {
                assert_eq!(
                    forward.to_bits(),
                    (-reverse).to_bits(),
                    "T={t1}, rho={rho1}"
                );
            }
            assert_eq!(
                forward.ideal_entropy_density_correction_j_mol_k().to_bits(),
                (-reverse.ideal_entropy_density_correction_j_mol_k()).to_bits()
            );
        }
        let zero = model
            .prepare_density_caloric_change(
                Temperature::new(300.0),
                Density::new(RHO),
                Density::new(RHO),
            )
            .unwrap()
            .evaluate(Temperature::new(300.0))
            .unwrap();
        assert_eq!(residuals(&zero), [0.0; 3]);
        assert_eq!(zero.ideal_entropy_density_correction_j_mol_k(), 0.0);
    }
}

#[test]
fn combined_entropy_uses_the_density_reference_and_original_cv() {
    let model = model();
    let reference = model
        .prepare_density_caloric_change(
            Temperature::new(300.0),
            Density::new(RHO),
            Density::new(RHO),
        )
        .unwrap();
    let t1 = 300.000_000_01;
    let change = reference.evaluate(Temperature::new(t1)).unwrap();
    let dt = t1 - 300.0;
    // Independent source Cp/R=3.25+0.00075T. Source Cv is Cp-R, not Cp.
    let ideal_standard_entropy = R * (3.25 * log_ratio(t1, 300.0) + 0.00075 * dt);
    let total_entropy = (ideal_standard_entropy
        + change.ideal_entropy_density_correction_j_mol_k()
        + change.residual_entropy_at_density_change_j_mol_k())
        / model.molar_mass().value();
    let expected = 3.197_199_127_111_231_4e-8; // Independent 85-digit total entropy.
    assert!((total_entropy / expected - 1.0).abs() < 5e-13);
    let derivative = change.residual_internal_energy_change_j_mol() / dt;
    assert!((derivative / reference.initial_state().residual_cv_j_mol_k() - 1.0).abs() < 1e-10);
    let entropy_derivative = change.residual_entropy_at_density_change_j_mol_k() / dt;
    assert!((entropy_derivative * 300.0 / derivative - 1.0).abs() < 1e-10);
}

#[test]
fn preparation_and_trials_keep_original_source_and_homogeneous_refusals() {
    let model = PengRobinsonMixtureV1::new(vec![components()[0]], vec![1.0], vec![]).unwrap();
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .prepare_density_caloric_change(
                    Temperature::new(300.0),
                    Density::new(bad),
                    Density::new(RHO)
                )
                .is_err()
        );
        assert!(
            model
                .prepare_density_caloric_change(
                    Temperature::new(300.0),
                    Density::new(RHO),
                    Density::new(bad)
                )
                .is_err()
        );
        let reference = model
            .prepare_density_caloric_change(
                Temperature::new(300.0),
                Density::new(RHO),
                Density::new(RHO),
            )
            .unwrap();
        assert_eq!(
            reference.evaluate(Temperature::new(bad)).unwrap_err(),
            model
                .evaluate_density_caloric(Temperature::new(bad), Density::new(RHO))
                .unwrap_err()
        );
    }
    assert!(matches!(
        model.prepare_density_caloric_change(
            Temperature::new(300.0),
            Density::new(RHO),
            Density::new(2000.0)
        ),
        Err(PengRobinsonErrorV1::OutsideCovolume)
    ));
    let pt = model
        .evaluate_pt(Temperature::new(150.0), Pressure::new(1e6))
        .unwrap();
    let unstable = pt.phase_roots()[1].mass_density_kg_per_m3();
    let reference = model
        .prepare_density_caloric_change(
            Temperature::new(300.0),
            Density::new(RHO),
            Density::new(unstable),
        )
        .unwrap();
    assert!(matches!(
        reference.evaluate(Temperature::new(150.0)),
        Err(PengRobinsonErrorV1::MechanicallyUnstable)
    ));
    let liquid = pt.select(PengRobinsonPhaseV1::Liquid);
    let reference = model
        .prepare_density_caloric_change(
            Temperature::new(150.0),
            Density::new(liquid.mass_density_kg_per_m3()),
            Density::new(liquid.mass_density_kg_per_m3()),
        )
        .unwrap();
    assert_eq!(
        residuals(&reference.evaluate(Temperature::new(150.0)).unwrap()),
        [0.0; 3]
    );
}
