//! Original source calorics on the existing implicit Fourier face/storage owner.
use super::*;
use fs_lbm::conjugate::{
    TensorNonlinearConfig, ThermalFace, ThermalSetup, ThermochemicalConductionControls,
    ThermochemicalConductionScheme,
};

const CONDUCTION_BUDGET_J: f64 = 1e-4;

fn conduction_controls(scheme: ThermochemicalConductionScheme) -> ThermochemicalConductionControls {
    let mut thermal = TensorNonlinearConfig::default();
    thermal.linear.tolerance = 1e-11;
    ThermochemicalConductionControls {
        thermal,
        scheme,
        maximum_caloric_queries: 100_000,
        maximum_enthalpy_defect_j: CONDUCTION_BUDGET_J,
    }
}

fn insulated() -> ThermalSetup {
    ThermalSetup::new([ThermalFace::Adiabatic; 6])
}

fn pure_h(bank: &ThermochemicalGasModel, t: f64) -> f64 {
    bank.evaluate(&source::point(t, 100_000.0))
        .unwrap()
        .specific_enthalpy_j_kg()
}

#[test]
fn gas_fourier_heat_flows_at_zero_species_diffusivity_with_exact_be_and_sdirk_storage() {
    let domain = VoxelDomain::graded([vec![0.5, 1.5], vec![1.0], vec![1.0]]).unwrap();
    let mass = [1.0, 3.0];
    let carrier = carrier(&domain, &mass);
    let before = carrier.clone();
    let species = [[mass[0]], [mass[1]]];
    let bank = source::constant_cp_model(-4000.0, 1.0);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let diffusivity = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let source =
        GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
            .unwrap();
    let initial_t = [480.0, 320.0];
    let h = initial_t.map(|t| pure_h(&bank, t));
    let initial = source
        .prepare(&domain, &carrier, &species, &h, &gate)
        .unwrap();
    let k = initial.transport(0).unwrap().thermal_conductivity_w_m_k();
    near(
        initial.transport(1).unwrap().thermal_conductivity_w_m_k(),
        k,
        1e-14,
    );
    let cp = 3.0 * R / MOLAR_MASS[0];
    near(
        initial.transport(0).unwrap().specific_heat_j_kg_k(),
        cp,
        1e-9,
    );
    let capacities = mass.map(|m| m * cp);
    let conductance = 1.0 / (0.25 / k + 0.75 / k);
    let duration = 10_000.0;
    let z = duration * conductance * (1.0 / capacities[0] + 1.0 / capacities[1]);
    let mean = (capacities[0] * initial_t[0] + capacities[1] * initial_t[1])
        / (capacities[0] + capacities[1]);
    for scheme in [
        ThermochemicalConductionScheme::BackwardEuler,
        ThermochemicalConductionScheme::Sdirk2,
    ] {
        let result = carrier
            .conduct_thermochemical_species(
                &domain,
                &species,
                &h,
                &insulated(),
                2.0,
                duration,
                &source,
                conduction_controls(scheme),
                &gate,
            )
            .unwrap();
        // Independent two-capacity linear mode. The native integrator must
        // retain its own BE or L-stable two-stage SDIRK stability function.
        let factor = match scheme {
            ThermochemicalConductionScheme::BackwardEuler => 1.0 / (1.0 + z),
            ThermochemicalConductionScheme::Sdirk2 => {
                let gamma = 1.0 - 1.0 / 2.0_f64.sqrt();
                (1.0 + (2.0 * gamma - 1.0) * z) / (1.0 + gamma * z).powi(2)
            }
        };
        let difference = (initial_t[0] - initial_t[1]) * factor;
        let expected = [
            mean + capacities[1] / (capacities[0] + capacities[1]) * difference,
            mean - capacities[0] / (capacities[0] + capacities[1]) * difference,
        ];
        for c in 0..2 {
            near(result.state.temperatures_k()[c], expected[c], 2e-7);
            near(
                result.specific_enthalpy_j_kg()[c] - h[c],
                cp * (expected[c] - initial_t[c]),
                3e-4,
            );
        }
        assert!(result.state.temperatures_k()[0] < initial_t[0] - 1.0);
        assert!(result.state.temperatures_k()[1] > initial_t[1] + 1.0);
        assert_eq!(result.state.diffusivity_m2_s(), diffusivity);
        assert!(result.thermal.record().closure_j.abs() <= CONDUCTION_BUDGET_J);
        near(
            result.thermal.record().stored_energy_change_j,
            0.0,
            CONDUCTION_BUDGET_J,
        );
        near(result.thermal.record().time_s, 2.0 + duration, 1e-10);
        assert_eq!(result.state.species(), &names);
        assert_eq!(result.state.source().source_packs(), bank.source_packs());
    }
    assert_eq!(carrier, before);
    assert_eq!(species, [[mass[0]], [mass[1]]]);
}

fn nonlinear_pair_heat(
    initial_t: [f64; 2],
    species: [[f64; 2]; 2],
    conductance: f64,
    duration: f64,
) -> f64 {
    let c = species.map(capacity);
    // Strip formation terms analytically: only the original sensible
    // integral 3T+.0005T² changes within each immutable composition.
    let g = |t: f64| 3.0 * t + 0.0005 * t * t;
    let total = [c[0] * g(initial_t[0]), c[1] * g(initial_t[1])];
    let joined = core::array::from_fn(|k| species[0][k] + species[1][k]);
    let equilibrium = temperature(total[0] + total[1], joined, 0.0);
    let mut lower = 0.0;
    let mut upper =
        c[0] * (initial_t[0] - equilibrium) * (3.0 + 0.0005 * (initial_t[0] + equilibrium));
    // Independent scalar monotone root of Q=G*dt*(T_left(Q)-T_right(Q)).
    // This does not invoke the tensor Newton solver or the source h inverse.
    for _ in 0..100 {
        let q = 0.5 * (lower + upper);
        let left = temperature(total[0] - q, species[0], 0.0);
        let right = temperature(total[1] + q, species[1], 0.0);
        if q > conductance * duration * (left - right) {
            upper = q;
        } else {
            lower = q;
        }
    }
    0.5 * (lower + upper)
}

#[test]
fn gas_fourier_graded_nonlinear_nasa_storage_is_conservative_and_reference_covariant() {
    let domain = VoxelDomain::graded([vec![0.5, 1.5], vec![1.0], vec![1.0]]).unwrap();
    let mass = [1.0, 2.0];
    let carrier = carrier(&domain, &mass);
    let species = [[0.75, 0.25], [0.25, 1.75]];
    let initial_t = [480.0, 320.0];
    let names = ids();
    let context = source::point(400.0, 100_000.0);
    let diffusivity = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let duration = 10_000.0;
    let mut previous: Option<[f64; 2]> = None;
    for datum in [0.0, -4000.0, 10000.0] {
        let bank = source::moist_model(0.5, datum);
        let source =
            GasSpeciesEnthalpySource::new(&bank, &names, &context, &diffusivity, energy(), &gate)
                .unwrap();
        let h: [f64; 2] = core::array::from_fn(|c| enthalpy(&bank, species[c], initial_t[c]));
        let opening = source
            .prepare(&domain, &carrier, &species, &h, &gate)
            .unwrap();
        let k = core::array::from_fn::<_, 2, _>(|c| {
            opening.transport(c).unwrap().thermal_conductivity_w_m_k()
        });
        let conductance = 1.0 / (0.25 / k[0] + 0.75 / k[1]);
        let q = nonlinear_pair_heat(initial_t, species, conductance, duration);
        let result = carrier
            .conduct_thermochemical_species(
                &domain,
                &species,
                &h,
                &insulated(),
                0.0,
                duration,
                &source,
                conduction_controls(ThermochemicalConductionScheme::BackwardEuler),
                &gate,
            )
            .unwrap();
        let expected = [
            temperature(mass[0] * h[0] - q, species[0], datum),
            temperature(mass[1] * h[1] + q, species[1], datum),
        ];
        for c in 0..2 {
            near(result.state.temperatures_k()[c], expected[c], 5e-7);
            let signed = if c == 0 { -q } else { q };
            near(
                mass[c] * (result.specific_enthalpy_j_kg()[c] - h[c]),
                signed,
                5e-4,
            );
        }
        let actual = [
            result.state.temperatures_k()[0],
            result.state.temperatures_k()[1],
        ];
        if let Some(old) = previous {
            for c in 0..2 {
                near(actual[c], old[c], 5e-7);
            }
        }
        previous = Some(actual);
        near(
            result.thermal.record().stored_energy_change_j,
            0.0,
            CONDUCTION_BUDGET_J,
        );
        assert!(result.thermal.record().closure_j.abs() <= CONDUCTION_BUDGET_J);
        assert!(result.state.absolute_caloric_residual_j() <= CONDUCTION_BUDGET_J);
        // The manufactured law is deliberately nonlinear enough that holding
        // opening Cp would fail the independent caloric check.
        let cp = core::array::from_fn::<_, 2, _>(|c| {
            capacity(species[c]) * (3.0 + 0.001 * initial_t[c])
        });
        let frozen_q = conductance * duration * (initial_t[0] - initial_t[1])
            / (1.0 + conductance * duration * (1.0 / cp[0] + 1.0 / cp[1]));
        assert!((q - frozen_q).abs() > 1.0);
    }
    assert_eq!(carrier.mass_kg, mass);
    assert_eq!(species, [[0.75, 0.25], [0.25, 1.75]]);
}

#[test]
fn gas_fourier_uniform_source_and_fixed_face_limits_keep_one_heat_account() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let carrier = carrier(&domain, &[1.0, 2.0]);
    let species = [[1.0], [2.0]];
    let bank = source::constant_cp_model(10000.0, 1.0);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let gate = CancelGate::new_clock_free();
    let d = [0.0; 2];
    let source =
        GasSpeciesEnthalpySource::new(&bank, &names, &context, &d, energy(), &gate).unwrap();
    let h = [pure_h(&bank, 400.0); 2];
    let uniform = carrier
        .conduct_thermochemical_species(
            &domain,
            &species,
            &h,
            &insulated(),
            0.0,
            10.0,
            &source,
            conduction_controls(ThermochemicalConductionScheme::BackwardEuler),
            &gate,
        )
        .unwrap();
    assert_eq!(uniform.specific_enthalpy_j_kg(), h);
    near(
        uniform.thermal.record().stored_energy_change_j,
        0.0,
        CONDUCTION_BUDGET_J,
    );

    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let carrier = super::carrier(&domain, &[2.0]);
    let species = [[2.0]];
    let d = [0.0];
    let source =
        GasSpeciesEnthalpySource::new(&bank, &names, &context, &d, energy(), &gate).unwrap();
    let h = [pure_h(&bank, 320.0)];
    let mut heating = insulated();
    heating.power_w = vec![2.0];
    let duration = 200.0;
    for scheme in [
        ThermochemicalConductionScheme::BackwardEuler,
        ThermochemicalConductionScheme::Sdirk2,
    ] {
        let result = carrier
            .conduct_thermochemical_species(
                &domain,
                &species,
                &h,
                &heating,
                0.0,
                duration,
                &source,
                conduction_controls(scheme),
                &gate,
            )
            .unwrap();
        near(
            2.0 * (result.specific_enthalpy_j_kg()[0] - h[0]),
            400.0,
            CONDUCTION_BUDGET_J,
        );
        near(result.thermal.record().source_j, 400.0, CONDUCTION_BUDGET_J);
        near(
            result.thermal.record().boundary_outflow_j,
            0.0,
            CONDUCTION_BUDGET_J,
        );
    }
    let opening = source
        .prepare(&domain, &carrier, &species, &h, &gate)
        .unwrap();
    let conductance = 2.0 * opening.transport(0).unwrap().thermal_conductivity_w_m_k();
    let capacity = 2.0 * 3.0 * R / MOLAR_MASS[0];
    let mut faces = insulated();
    faces.faces[1] = ThermalFace::Temperature(480.0);
    let result = carrier
        .conduct_thermochemical_species(
            &domain,
            &species,
            &h,
            &faces,
            0.0,
            duration,
            &source,
            conduction_controls(ThermochemicalConductionScheme::BackwardEuler),
            &gate,
        )
        .unwrap();
    let expected =
        (capacity * 320.0 + duration * conductance * 480.0) / (capacity + duration * conductance);
    let heat = capacity * (expected - 320.0);
    near(result.state.temperatures_k()[0], expected, 2e-7);
    near(
        result.thermal.record().boundary_outflow_j,
        -heat,
        CONDUCTION_BUDGET_J,
    );
    near(
        result.thermal.record().stored_energy_change_j,
        heat,
        CONDUCTION_BUDGET_J,
    );

    // An authored cell temperature obtains its unique total h from the same
    // original mixture chart; its reservoir supplies the actual storage jump.
    let mut pinned = insulated();
    pinned.fixed_temperature.push((0, 400.0));
    let fixed = carrier
        .conduct_thermochemical_species(
            &domain,
            &species,
            &h,
            &pinned,
            0.0,
            duration,
            &source,
            conduction_controls(ThermochemicalConductionScheme::BackwardEuler),
            &gate,
        )
        .unwrap();
    let prescribed_h = pure_h(&bank, 400.0);
    near(fixed.specific_enthalpy_j_kg()[0], prescribed_h, 1e-7);
    near(fixed.state.temperatures_k()[0], 400.0, 1e-9);
    let prescribed_heat = 2.0 * (prescribed_h - h[0]);
    near(
        fixed.thermal.record().stored_energy_change_j,
        prescribed_heat,
        CONDUCTION_BUDGET_J,
    );
    near(
        fixed.thermal.record().boundary_outflow_j,
        -prescribed_heat,
        CONDUCTION_BUDGET_J,
    );
}

#[test]
fn gas_fourier_source_domain_query_and_cancellation_refusals_are_atomic() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let carrier = carrier(&domain, &[1.0; 2]);
    let species = [[1.0]; 2];
    let bank = source::constant_cp_model(0.0, 1.0);
    let names = [SpeciesId::new("manufactured").unwrap()];
    let context = source::point(400.0, 100_000.0);
    let d = [0.0; 2];
    let gate = CancelGate::new_clock_free();
    let source =
        GasSpeciesEnthalpySource::new(&bank, &names, &context, &d, energy(), &gate).unwrap();
    let h = [pure_h(&bank, 320.0), pure_h(&bank, 480.0)];
    let before = carrier.clone();
    let mut limited = conduction_controls(ThermochemicalConductionScheme::BackwardEuler);
    // Four bracket evaluations and two original initial inversions fit; the
    // first native trial then exhausts the shared query owner.
    limited.maximum_caloric_queries = 6;
    let exhausted = carrier
        .conduct_thermochemical_species(
            &domain,
            &species,
            &h,
            &insulated(),
            0.0,
            1.0,
            &source,
            limited,
            &gate,
        )
        .unwrap_err();
    assert!(
        matches!(exhausted, ChtError::InvalidInput { ref reason, .. }
        if reason.contains("caloric-query budget"))
    );
    // Sources are initially valid, but requested heating leaves their original
    // admitted temperature branch. No partial hot endpoint may be returned.
    let mut heating = insulated();
    heating.power_w = vec![1e6; 2];
    assert!(
        carrier
            .conduct_thermochemical_species(
                &domain,
                &species,
                &h,
                &heating,
                0.0,
                10.0,
                &source,
                conduction_controls(ThermochemicalConductionScheme::BackwardEuler),
                &gate,
            )
            .is_err()
    );
    let mut fixed_cell = insulated();
    fixed_cell.fixed_temperature.push((0, 900.0));
    assert!(
        carrier
            .conduct_thermochemical_species(
                &domain,
                &species,
                &h,
                &fixed_cell,
                0.0,
                1.0,
                &source,
                conduction_controls(ThermochemicalConductionScheme::BackwardEuler),
                &gate,
            )
            .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        carrier.conduct_thermochemical_species(
            &domain,
            &species,
            &h,
            &insulated(),
            0.0,
            1.0,
            &source,
            conduction_controls(ThermochemicalConductionScheme::BackwardEuler),
            &cancelled,
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(carrier, before);
    assert_eq!(species, [[1.0]; 2]);
    assert_eq!(h, [pure_h(&bank, 320.0), pure_h(&bank, 480.0)]);
}
