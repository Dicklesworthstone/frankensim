//! Simultaneous creep/plasticity: independent solutions and constitutive checks.
use fs_material::creep::plastic::{NortonJ2, NortonJ2State};
use fs_material::creep::{
    CreepError, CreepLaw, CreepTolerances, NortonCreep, NortonCreepParameters,
};
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::tensor::{contract, rotate, rotation, von_mises};
use fs_material::{SmallStrainLaw, Voigt};

fn elastic() -> IsotropicElastic {
    IsotropicElastic::new(1.2e9, 0.2, 0.05).unwrap()
}
fn creep(n: f64, dt: f64) -> NortonCreep {
    NortonCreep::new(
        elastic(),
        NortonCreepParameters {
            reference_stress_pa: 3.0e6,
            reference_rate_s: 2.0e-3,
            stress_exponent: n,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 500.0,
            temperature_range_k: [300.0, 900.0],
        },
        500.0,
        dt,
    )
    .unwrap()
}
fn law(n: f64, dt: f64, yield_pa: f64, hardening: f64) -> NortonJ2 {
    NortonJ2::new(
        creep(n, dt),
        J2Plasticity::new(elastic(), yield_pa, hardening).unwrap(),
    )
    .unwrap()
}
fn tolerances() -> CreepTolerances {
    CreepTolerances {
        stress_pa: 2.0e-6,
        energy_j_m3: 1.0e-7,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}
fn close_tensor(a: &Voigt, b: &Voigt, tolerance: f64) {
    for (x, y) in a.iter().zip(b) {
        close(*x, *y, tolerance);
    }
}
fn deviatoric(e: f64) -> Voigt {
    [e, -0.5 * e, -0.5 * e, 0.0, 0.0, 0.0]
}

#[test]
fn linear_simultaneous_solution_and_energy_include_perfect_plasticity() {
    let (_, g) = elastic().lame();
    let h = 3.0 * g;
    let strain = deviatoric(0.006);
    let q_trial = h * strain[0];
    let a = 0.25 * 2.0e-3 / 3.0e6; // dc = a q for the linear Norton law.
    for hardening in [0.0, 1.0e8] {
        let law = law(1.0, 0.25, 2.0e6, hardening);
        let old = law.initial_state();
        let next = law
            .creep_increment(&[0.0; 6], &old, &strain, tolerances())
            .unwrap();
        // Solve the two linear equations independently, without a return map.
        let q = (hardening * q_trial + h * 2.0e6) / (hardening + h + h * hardening * a);
        let dc = a * q;
        let dp = (q_trial - q) / h - dc;
        close(von_mises(&next.stress_pa), q, 2.0e-8);
        close(next.state.creep.equivalent_creep_strain, dc, 2.0e-17);
        close(next.state.plastic.alpha, dp, 2.0e-17);
        close_tensor(&next.state.creep.creep_strain, &deviatoric(dc), 2.0e-17);
        close_tensor(&next.state.plastic.plastic_strain, &deviatoric(dp), 2.0e-17);
        let stored = q * q / (2.0 * h) + 0.5 * hardening * dp * dp;
        close(next.energy.final_stored_j_m3, stored, 2.0e-8);
        close(next.energy.numerical_dissipation_j_m3, stored, 2.0e-8);
        close(next.energy.creep_dissipation_j_m3, q * dc, 2.0e-8);
        close(next.energy.plastic_dissipation_j_m3, 2.0e6 * dp, 2.0e-8);
        close(
            next.energy.physical_dissipation_j_m3,
            q * dc + 2.0e6 * dp,
            2.0e-8,
        );
        close(next.energy.mechanical_work_j_m3, q * strain[0], 2.0e-8);
        assert!(next.energy.residual_j_m3.abs() < 1.0e-8);
        let axial_direction = deviatoric(1.0);
        let tangent_direction: Voigt = core::array::from_fn(|i| {
            (0..6)
                .map(|j| next.tangent_pa[i][j] * axial_direction[j])
                .sum()
        });
        let radial = hardening / (hardening + h + h * hardening * a);
        close_tensor(&tangent_direction, &deviatoric(2.0 * g * radial), 2.0e-6);
        assert_eq!(old, law.initial_state());
    }
}

#[test]
fn creep_can_relax_an_above_yield_predictor_without_plastic_flow() {
    let law = law(1.0, 0.25, 2.0e6, 1.0e8);
    let old = law.initial_state();
    let strain = deviatoric(0.0016); // q_trial=2.4 MPa, relaxed q=1.92 MPa.
    let expected = law
        .creep()
        .advance_with_energy(&[0.0; 6], &old.creep, &strain, tolerances())
        .unwrap();
    let next = law
        .advance_with_energy(&[0.0; 6], &old, &strain, tolerances())
        .unwrap();
    assert_eq!(next.stress_pa, expected.stress_pa);
    assert_eq!(next.tangent_pa, expected.tangent_pa);
    assert_eq!(next.state.creep, expected.state);
    assert_eq!(next.state.plastic, old.plastic);
    assert_eq!(next.energy, expected.energy);
    assert_eq!(
        next,
        law.advance_with_energy(&[0.0; 6], &old, &strain, tolerances())
            .unwrap()
    );
    let zero = law.response(&[0.0; 6], &old).unwrap();
    assert_eq!(zero.stress_pa, [0.0; 6]);
    assert_eq!(zero.state, old);
    assert_eq!(zero.tangent_pa, law.creep().tangent(&[0.0; 6], &old.creep));
}

#[test]
fn nonlinear_coupled_tangent_matches_all_tensor_columns() {
    let law = law(3.0, 0.125, 2.0e6, 2.0e8);
    let strain = [0.005, -0.001, -0.002, 0.0003, -0.0005, 0.0002];
    let old = law.initial_state();
    let next = law
        .advance_with_energy(&[0.0; 6], &old, &strain, tolerances())
        .unwrap();
    assert!(next.state.plastic.alpha > 0.001);
    assert!(next.state.creep.equivalent_creep_strain > 0.0);
    let q = von_mises(&next.stress_pa);
    close(
        next.state.creep.equivalent_creep_strain,
        0.125 * 2.0e-3 * (q / 3.0e6).powi(3),
        2.0e-17,
    );
    close(q, 2.0e6 + 2.0e8 * next.state.plastic.alpha, 2.0e-8);
    let step = 2.0e-9;
    for j in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += step;
        minus[j] -= step;
        let sp = law.response(&plus, &old).unwrap().stress_pa;
        let sm = law.response(&minus, &old).unwrap().stress_pa;
        for i in 0..6 {
            let fd = (sp[i] - sm[i]) / (2.0 * step);
            close(
                fd,
                next.tangent_pa[i][j],
                2.0e-6 * next.tangent_pa[i][j].abs().max(1.0e6),
            );
            let wi = if i < 3 { 1.0 } else { 2.0 };
            let wj = if j < 3 { 1.0 } else { 2.0 };
            close(
                wi * next.tangent_pa[i][j],
                wj * next.tangent_pa[j][i],
                2.0e-6,
            );
        }
    }
}

#[test]
fn both_retained_histories_rotate_and_ignore_hydrostatic_shifts() {
    let law = law(3.0, 0.125, 2.0e6, 2.0e8);
    let first_strain = deviatoric(0.004);
    let first = law
        .advance_with_energy(&[0.0; 6], &law.initial_state(), &first_strain, tolerances())
        .unwrap();
    let strain = [0.007, -0.003, -0.004, 0.0004, -0.0003, 0.0007];
    let next = law
        .advance_with_energy(&first_strain, &first.state, &strain, tolerances())
        .unwrap();
    let frame = rotation([1.0 / 3.0_f64.sqrt(); 3], 0.7);
    let mut rotated_state = first.state.clone();
    rotated_state.creep.creep_strain = rotate(&rotated_state.creep.creep_strain, &frame);
    rotated_state.plastic.plastic_strain = rotate(&rotated_state.plastic.plastic_strain, &frame);
    let rotated = law
        .advance_with_energy(
            &rotate(&first_strain, &frame),
            &rotated_state,
            &rotate(&strain, &frame),
            tolerances(),
        )
        .unwrap();
    close_tensor(&rotated.stress_pa, &rotate(&next.stress_pa, &frame), 2.0e-7);
    close_tensor(
        &rotated.state.creep.creep_strain,
        &rotate(&next.state.creep.creep_strain, &frame),
        2.0e-17,
    );
    close_tensor(
        &rotated.state.plastic.plastic_strain,
        &rotate(&next.state.plastic.plastic_strain, &frame),
        2.0e-17,
    );
    close(
        rotated.state.plastic.alpha,
        next.state.plastic.alpha,
        2.0e-17,
    );
    close(
        rotated.energy.physical_dissipation_j_m3,
        next.energy.physical_dissipation_j_m3,
        2.0e-8,
    );
    let shifted = |e: Voigt| core::array::from_fn(|i| e[i] + if i < 3 { 0.002 } else { 0.0 });
    let hydro = law
        .advance_with_energy(
            &shifted(first_strain),
            &first.state,
            &shifted(strain),
            tolerances(),
        )
        .unwrap();
    close_tensor(
        &hydro.state.creep.creep_strain,
        &next.state.creep.creep_strain,
        2.0e-17,
    );
    close_tensor(
        &hydro.state.plastic.plastic_strain,
        &next.state.plastic.plastic_strain,
        2.0e-17,
    );
    close(
        hydro.energy.physical_dissipation_j_m3,
        next.energy.physical_dissipation_j_m3,
        2.0e-8,
    );
}

#[test]
fn unloading_keeps_both_permanent_strains_and_reversal_grows_both_paths() {
    let law = law(1.0, 0.25, 2.0e6, 1.0e8);
    let strain = deviatoric(0.006);
    let first = law
        .advance_with_energy(&[0.0; 6], &law.initial_state(), &strain, tolerances())
        .unwrap();
    let permanent: Voigt = core::array::from_fn(|i| {
        first.state.creep.creep_strain[i] + first.state.plastic.plastic_strain[i]
    });
    let unloaded = law
        .advance_with_energy(&strain, &first.state, &permanent, tolerances())
        .unwrap();
    close_tensor(&unloaded.stress_pa, &[0.0; 6], 2.0e-8);
    close(unloaded.energy.physical_dissipation_j_m3, 0.0, 1.0e-12);
    close(
        unloaded.state.plastic.alpha,
        first.state.plastic.alpha,
        2.0e-17,
    );
    close(
        unloaded.state.creep.equivalent_creep_strain,
        first.state.creep.equivalent_creep_strain,
        2.0e-17,
    );
    let reversed = law
        .advance_with_energy(
            &permanent,
            &unloaded.state,
            &deviatoric(-0.006),
            tolerances(),
        )
        .unwrap();
    assert!(reversed.state.creep.creep_strain[0] < 0.0);
    assert!(reversed.state.plastic.plastic_strain[0] < 0.0);
    assert!(
        reversed.state.creep.equivalent_creep_strain > first.state.creep.equivalent_creep_strain
    );
    assert!(reversed.state.plastic.alpha > first.state.plastic.alpha);
    assert!(reversed.energy.creep_dissipation_j_m3 > 0.0);
    assert!(reversed.energy.plastic_dissipation_j_m3 > 0.0);
}

#[test]
fn changed_yield_offset_causes_held_strain_flow_with_zero_external_work() {
    let hardening = 1.0e8;
    let strong = law(1.0, 0.25, 20.0e6, hardening);
    let soft = law(1.0, 0.25, 2.0e6, hardening);
    let strain = deviatoric(0.006);
    let first = strong
        .advance_with_energy(&[0.0; 6], &strong.initial_state(), &strain, tolerances())
        .unwrap();
    assert_eq!(first.state.plastic.alpha, 0.0);
    let next = soft
        .advance_with_energy(&strain, &first.state, &strain, tolerances())
        .unwrap();
    assert!(next.state.plastic.alpha > 0.0);
    assert!(next.state.creep.equivalent_creep_strain > first.state.creep.equivalent_creep_strain);
    close(next.energy.mechanical_work_j_m3, 0.0, 0.0);
    close(
        next.energy.initial_stored_j_m3,
        first.energy.final_stored_j_m3,
        2.0e-8,
    );
    close(
        next.energy.initial_stored_j_m3 - next.energy.final_stored_j_m3,
        next.energy.physical_dissipation_j_m3 + next.energy.numerical_dissipation_j_m3,
        2.0e-8,
    );
    // Hardening storage and integration loss must not be counted as plastic heat.
    close(
        next.energy.plastic_dissipation_j_m3,
        2.0e6 * next.state.plastic.alpha,
        2.0e-8,
    );
    let e: Voigt = core::array::from_fn(|i| {
        strain[i] - next.state.creep.creep_strain[i] - next.state.plastic.plastic_strain[i]
    });
    let elastic_storage = 0.5 * contract(&e, &elastic().stress(&e, &()));
    close(
        next.energy.final_stored_j_m3 - elastic_storage,
        0.5 * hardening * next.state.plastic.alpha.powi(2),
        2.0e-8,
    );
}

#[test]
fn unsupported_models_and_hidden_invalid_component_histories_refuse() {
    let mut invalid = J2Plasticity::new(elastic(), 2.0e6, 1.0e8).unwrap();
    invalid.hardening = -1.0;
    assert!(matches!(
        NortonJ2::new(creep(1.0, 0.25), invalid),
        Err(CreepError::Plastic(_))
    ));
    let mixed = J2Plasticity::new_mixed(elastic(), 2.0e6, 1.0e8, 1.0e8).unwrap();
    assert!(NortonJ2::new(creep(1.0, 0.25), mixed).is_err());
    let different = J2Plasticity::new(
        IsotropicElastic::new(1.3e9, 0.2, 0.05).unwrap(),
        2.0e6,
        1.0e8,
    )
    .unwrap();
    assert!(NortonJ2::new(creep(1.0, 0.25), different).is_err());
    let law = law(1.0, 0.25, 2.0e6, 1.0e8);
    let old = law.initial_state();
    let mut hidden: NortonJ2State = old.clone();
    hidden.creep.creep_strain = deviatoric(0.06);
    hidden.creep.equivalent_creep_strain = 0.06;
    hidden.plastic.plastic_strain = deviatoric(-0.06);
    hidden.plastic.alpha = 0.06;
    assert!(law.validate_history(&[0.0; 6], &hidden).is_err());
    let mut traceful = old.clone();
    traceful.plastic.plastic_strain[0] = 0.001;
    traceful.plastic.alpha = 0.001;
    assert!(matches!(
        law.validate_history(&[0.0; 6], &traceful),
        Err(CreepError::Plastic(_))
    ));
    assert!(law.response(&[f64::NAN; 6], &old).is_err());
    assert!(law.stress(&[f64::NAN; 6], &old).iter().all(|x| x.is_nan()));
    let invalid_tolerance = CreepTolerances {
        stress_pa: 0.0,
        ..tolerances()
    };
    assert!(
        law.advance_with_energy(&[0.0; 6], &old, &deviatoric(0.006), invalid_tolerance)
            .is_err()
    );
    assert_eq!(old, law.initial_state());
}

#[test]
fn exact_plastic_activation_and_adjacent_ulps_remain_admissible() {
    let mut parameters = creep(1.0, 1.0).parameters().clone();
    parameters.reference_stress_pa = 1.0e6;
    parameters.reference_rate_s = 1.0e-5;
    let creep = NortonCreep::new(elastic(), parameters, 500.0, 1.0).unwrap();
    let law = NortonJ2::new(creep, J2Plasticity::new(elastic(), 2.0e6, 4.0e8).unwrap()).unwrap();
    let old = law.initial_state();
    let h = 3.0 * elastic().lame().1;
    // At q=Y the exact linear rate gives dc=2e-5, hence q_trial=2.03 MPa.
    // log-space evaluation must not turn its roundoff into negative plasticity.
    let switch: f64 = 2.03e6;
    for offset in -8_i64..=8 {
        let q_trial = f64::from_bits((switch.to_bits() as i64 + offset) as u64);
        let strain = deviatoric(q_trial / h);
        let next = law
            .advance_with_energy(&[0.0; 6], &old, &strain, tolerances())
            .unwrap_or_else(|error| panic!("activation offset {offset}: {error}"));
        close(von_mises(&next.stress_pa), 2.0e6, 1.0e-7);
        close(next.state.creep.equivalent_creep_strain, 2.0e-5, 2.0e-18);
        assert!(next.state.plastic.alpha >= 0.0 && next.state.plastic.alpha < 2.0e-17);
        assert!(next.energy.plastic_dissipation_j_m3 >= 0.0);
        assert!(next.energy.residual_j_m3.abs() < tolerances().energy_j_m3);
    }
    assert_eq!(old, law.initial_state());
}
