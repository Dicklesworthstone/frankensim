//! Manufactured finite-gas/film laws, not empirical membrane accuracy claims.
//! Every numerical update is the original public material/gas/membrane owner.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{
    boundary::GasSpeciesBoundarySolve,
    network::transport::membrane::{
        MEMBRANE_DIFFUSIVITY_DIMS, MOLAR_HENRY_SOLUBILITY_DIMS, MembranePermeabilityMaterial,
        transient::{
            MEMBRANE_DIFFUSIVITY_PROPERTY, MOLAR_HENRY_SOLUBILITY_PROPERTY,
            MembraneDiffusionLayerSource, MembraneDiffusionSolve,
            mixture::{
                IsothermalGasMembraneMixture, MembraneMixtureChannelSource,
                adaptive::{MembraneMixtureAdvanceControls, MembraneMixtureAdvanceError},
            },
            reservoir::{IsothermalGasMembraneError, IsothermalGasMembraneSolve},
        },
    },
};
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, DynViscosity, MolarMass, Pressure, QuantitySpec};
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R, mixture::FrozenIdealGasMixtureModelV1,
    required_nasa9_parameter_dims,
};

const T: f64 = 300.0;
const AREA: f64 = 0.01;
const THICKNESS: f64 = 1e-4;
const S: f64 = 2e-5;
const D: f64 = 1e-10;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn typed(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), T)
        .unwrap()
        .with_quantity("pressure", typed(QuantityKind::Pressure), p)
        .unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "manufactured finite membrane oracle".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    }
}
fn claim(name: &str, value: f64, dims: Dims, pressures: [f64; 2]) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                typed(QuantityKind::AbsoluteTemperature),
                200.0,
                600.0,
            )
            .with_quantity(
                "pressure",
                typed(QuantityKind::Pressure),
                pressures[0],
                pressures[1],
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: provenance(),
    }
}
fn card(name: &str, values: &[(&str, f64, Dims)], pressures: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(key, value, dims) in values {
        claims
            .insert_claim(claim(key, value, dims, pressures))
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "synthetic".into(),
            process: "analytic test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn membrane(name: &str, d: f64, s: f64, maximum_pressure: f64) -> MaterialCard {
    card(
        name,
        &[
            (MEMBRANE_DIFFUSIVITY_PROPERTY, d, MEMBRANE_DIFFUSIVITY_DIMS),
            (
                MOLAR_HENRY_SOLUBILITY_PROPERTY,
                s,
                MOLAR_HENRY_SOLUBILITY_DIMS,
            ),
        ],
        [0.0, maximum_pressure],
    )
}
fn species(name: &str, mass: f64, reference: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = 3.5;
    coefficients[7] = mass * reference / R;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
        sources: vec![],
        provenance: provenance(),
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic-reference").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn model(fractions: &[f64], reference: f64) -> ThermochemicalGasModel {
    let names = ["A", "B"];
    let masses = [0.028, 0.044];
    let components = (0..fractions.len())
        .map(|i| species(names[i], masses[i], reference))
        .collect();
    let transport = (0..fractions.len())
        .map(|i| SpeciesTransportMaterial {
            species: id(names[i]),
            material: card(
                names[i],
                &[
                    ("dynamic_viscosity", 1.7e-5, DynViscosity::DIMS),
                    (
                        THERMAL_CONDUCTIVITY_PROPERTY,
                        0.025,
                        THERMAL_CONDUCTIVITY_DIMS,
                    ),
                ],
                [100.0, 1e7],
            ),
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            components,
            Composition::new(CompositionBasis::MoleFraction, fractions.to_vec()).unwrap(),
        )
        .unwrap(),
        transport,
        None,
    )
    .unwrap()
}
fn gas(fractions: &[f64], reference: f64, pressure: f64, volume: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        model(fractions, reference),
        point(pressure).with("specimen-lot", 17.0).unwrap(),
        volume,
    )
    .unwrap()
}
fn channel<'a>(
    name: &str,
    card: &'a MaterialCard,
    context: &'a QueryPoint,
    cells: usize,
    loading: f64,
) -> MembraneMixtureChannelSource<'a> {
    MembraneMixtureChannelSource {
        layers: vec![MembraneDiffusionLayerSource {
            thickness_m: THICKNESS,
            cells,
            material: MembranePermeabilityMaterial {
                species: id(name),
                material: card,
                context,
                selection: MaterialPropertySelection::SingleClaimOnly,
            },
        }],
        interfaces: vec![],
        initial_partial_pressures_pa: vec![loading; cells],
    }
}
fn native_controls() -> IsothermalGasMembraneSolve {
    IsothermalGasMembraneSolve {
        membrane: MembraneDiffusionSolve {
            absolute_amount_tolerance_mol: 1e-14,
            relative_amount_tolerance: 1e-9,
        },
        gas: [GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 600.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-16,
            energy_tolerance_j: 1e-8,
        }; 2],
        mass_tolerance_kg: 1e-16,
    }
}
fn controls(step: f64, relative: f64) -> MembraneMixtureAdvanceControls {
    MembraneMixtureAdvanceControls {
        initial_step_s: step,
        minimum_step_s: 1e-8,
        maximum_step_s: step,
        maximum_evaluations: 20000,
        maximum_total_cells: 1024,
        absolute_pressure_tolerance_pa: 0.01,
        absolute_amount_tolerance_mol: 1e-14,
        relative_tolerance: relative,
        species_balance_tolerance_mol: 1e-12,
        mass_balance_tolerance_kg: 1e-14,
        gas_energy_tolerance_j: 1e-10,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; tolerance {tolerance}"
    );
}
fn film<'a>(state: &'a IsothermalGasMembraneMixture<'_>, name: &str) -> &'a [f64] {
    state
        .finite_state(&id(name))
        .unwrap()
        .membrane()
        .partial_pressures_pa()
}
fn equal_capacity_volume() -> f64 {
    AREA * S * THICKNESS * R * T
}

#[test]
fn adaptive_pure_relaxation_converges_to_independent_three_capacity_mode() {
    let card = membrane("A", D, S, 4e5);
    let context = point(1e5);
    let sources = [channel("A", &card, &context, 1, 1e5)];
    let volume = equal_capacity_volume();
    let gases = [
        gas(&[1.0], 0.0, 1.4e5, volume),
        gas(&[1.0], 0.0, 0.6e5, volume),
    ];
    let state = IsothermalGasMembraneMixture::new(
        T,
        AREA,
        &sources,
        gases.clone(),
        native_controls(),
        || false,
    )
    .unwrap();
    // One film cell and equal gas/film pressure capacities. The film remains
    // at the mean. Independent ODE: d(pL-pR)/dt = -(2 D/L²)(pL-pR).
    let duration = 50.0;
    let delta = 40000.0 * (-2.0 * D * duration / THICKNESS.powi(2)).exp();
    let expected = [100000.0 + delta, 100000.0 - delta];
    let mut errors = Vec::new();
    for relative in [0.01, 0.0025] {
        let mut work = state.adaptive_workspace().unwrap();
        let result = state
            .advance_candidate(duration, controls(duration, relative), &mut work, || false)
            .unwrap();
        let report = result.report();
        assert!(report.accepted_intervals > 1);
        assert!(report.rejected_intervals > 0);
        assert_eq!(
            report.evaluations,
            3 * (report.accepted_intervals + report.rejected_intervals)
        );
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert_eq!(report.interval_s, [0.0, duration]);
        assert_eq!(result.state().elapsed_s(), duration);
        assert!(report.maximum_accepted_error_ratio <= 1.0);
        let actual = result
            .state()
            .gases()
            .each_ref()
            .map(|g| g.state().state().pressure);
        errors.push((actual[0] - expected[0]).abs());
        near(actual[0] + actual[1], 2e5, 1e-7);
        near(film(result.state(), "A")[0], 1e5, 1e-7);
        let row = &report.species[0];
        assert!(row.amount_into_gas_mol[0] < 0.0 && row.amount_into_gas_mol[1] > 0.0);
        for side in 0..2 {
            let initial = [1.4e5, 0.6e5][side];
            let dn = volume / (R * T) * (actual[side] - initial);
            near(row.amount_into_gas_mol[side], dn, 1e-17);
            near(row.thermostat_heat_into_gas_j[side], -R * T * dn, 1e-13);
            near(row.gas_boundary_residuals_mol[side], 0.0, 1e-17);
            assert_eq!(report.source_relative_gas_energy_residual_j[side], 0.0);
        }
        assert_eq!(row.gross_amount_into_gas_mol[0], 0.0);
        assert_eq!(row.gross_amount_out_of_gas_mol[1], 0.0);
        near(
            row.gross_amount_out_of_gas_mol[0],
            -row.amount_into_gas_mol[0],
            1e-20,
        );
        near(
            row.gross_amount_into_gas_mol[1],
            row.amount_into_gas_mol[1],
            1e-20,
        );
        near(row.closed_species_residual_mol, 0.0, 1e-17);
        near(report.total_mass_residual_kg, 0.0, 1e-18);
    }
    assert!(errors[1] < 0.4 * errors[0], "refinement errors: {errors:?}");
    assert!(errors[1] < 60.0, "tight mode error: {} Pa", errors[1]);
    assert_eq!(state.gases(), &gases);
    assert_eq!(film(&state, "A"), &[1e5]);
    assert_eq!(state.elapsed_s(), 0.0);
}

#[test]
fn decimal_terminal_ulp_is_solved_in_the_preceding_interval_not_clock_snapped() {
    let card = membrane("A", D, S, 4e5);
    let context = point(1e5);
    let sources = [channel("A", &card, &context, 1, 1e5)];
    let volume = equal_capacity_volume();
    let state = IsothermalGasMembraneMixture::new(
        T,
        AREA,
        &sources,
        [
            gas(&[1.0], 0.0, 1.4e5, volume),
            gas(&[1.0], 0.0, 0.6e5, volume),
        ],
        native_controls(),
        || false,
    )
    .unwrap();
    let duration = 0.005_f64;
    let mut work = state.adaptive_workspace().unwrap();
    // The nominal ceiling deliberately leaves an indivisible one-ULP tail.
    // Both actual physical halves must instead cover the full decimal target.
    let result = state
        .advance_candidate(
            duration,
            controls(duration.next_down(), 0.1),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(result.state().elapsed_s(), duration);
    assert_eq!(result.report().accepted_intervals, 1);
    assert_eq!(result.report().accepted_substeps, 2);
    assert_eq!(result.report().evaluations, 3);
    assert_eq!(result.report().minimum_accepted_substep_s, 0.5 * duration);
    assert_eq!(result.report().maximum_accepted_substep_s, 0.5 * duration);
    // Independent discrete antisymmetric mode for the retained two BE halves.
    let delta = 40000.0 / (1.0 + D * duration / THICKNESS.powi(2)).powi(2);
    near(
        result.state().gases()[0].state().state().pressure,
        100000.0 + delta,
        1e-8,
    );
    near(
        result.state().gases()[1].state().state().pressure,
        100000.0 - delta,
        1e-8,
    );
    // A whole requested window with no represented midpoint still refuses.
    let tiny = f64::from_bits(1);
    let mut c = controls(tiny, 0.1);
    c.minimum_step_s = tiny;
    assert!(
        matches!(state.advance_candidate(tiny, c, &mut work, || false),
        Err(MembraneMixtureAdvanceError::TimeResolution { time_s: 0.0, target_s })
            if target_s == tiny)
    );
    assert_eq!(state.elapsed_s(), 0.0);
}

#[test]
fn opposing_species_errors_are_controlled_when_total_bath_pressure_cancels() {
    let a = membrane("A", D, S, 4e5);
    let b = membrane("B", D, S, 4e5);
    let context = point(1e5);
    let sources = [
        channel("A", &a, &context, 1, 5e4),
        channel("B", &b, &context, 1, 5e4),
    ];
    let volume = equal_capacity_volume();
    let state = IsothermalGasMembraneMixture::new(
        T,
        AREA,
        &sources,
        [
            gas(&[0.9, 0.1], 0.0, 1e5, volume),
            gas(&[0.1, 0.9], 0.0, 1e5, volume),
        ],
        native_controls(),
        || false,
    )
    .unwrap();
    let mut native = state.workspace().unwrap();
    let coarse = state.step_candidate(50.0, &mut native, || false).unwrap();
    let half = state.step_candidate(25.0, &mut native, || false).unwrap();
    let fine = half
        .state()
        .step_candidate(25.0, &mut native, || false)
        .unwrap();
    for side in 0..2 {
        near(
            coarse.state().gases()[side].state().state().pressure,
            fine.state().gases()[side].state().state().pressure,
            1e-8,
        );
    }
    let partial = |s: &IsothermalGasMembraneMixture<'_>| {
        s.gases()[0].species_partial_pressure_pa(&id("A")).unwrap()
    };
    assert!((partial(coarse.state()) - partial(fine.state())).abs() > 1000.0);
    let mut work = state.adaptive_workspace().unwrap();
    let result = state
        .advance_candidate(50.0, controls(50.0, 0.005), &mut work, || false)
        .unwrap();
    assert!(result.report().rejected_intervals > 0);
    let expected = 50000.0 + 40000.0 * (-1.0_f64).exp();
    near(partial(result.state()), expected, 120.0);
    for side in 0..2 {
        near(
            result.state().gases()[side].state().state().pressure,
            1e5,
            1e-7,
        );
    }
    let rows = &result.report().species;
    assert!(rows[0].amount_into_gas_mol[0] < 0.0);
    assert!(rows[1].amount_into_gas_mol[0] > 0.0);
    near(
        rows[0].amount_into_gas_mol[0] + rows[1].amount_into_gas_mol[0],
        0.0,
        1e-17,
    );
    assert!(rows[0].gross_amount_out_of_gas_mol[0] > 1e-7);
    assert!(rows[1].gross_amount_into_gas_mol[0] > 1e-7);
}

#[test]
fn blocked_species_reference_offsets_and_channel_order_survive_complete_windows() {
    let a = membrane("A", D, S, 4e5);
    let b = membrane("blocked B", 0.0, 1e-5, 4e5);
    let context = point(1e5);
    let mut expected = None;
    for (reference, reversed) in [(0.0, false), (1e24, false), (1e24, true)] {
        let mut sources = vec![
            channel("A", &a, &context, 4, 5e4),
            channel("B", &b, &context, 4, 12345.0),
        ];
        if reversed {
            sources.reverse();
        }
        let state = IsothermalGasMembraneMixture::new(
            T,
            AREA,
            &sources,
            [
                gas(&[0.8, 0.2], reference, 1e5, 1e-7),
                gas(&[0.2, 0.8], reference, 1e5, 3e-7),
            ],
            native_controls(),
            || false,
        )
        .unwrap();
        let mut work = state.adaptive_workspace().unwrap();
        let result = state
            .advance_candidate(2.0, controls(1.0, 0.1), &mut work, || false)
            .unwrap();
        let report = result.report();
        let blocked = &report.species[1];
        assert_eq!(blocked.species, id("B"));
        assert_eq!(blocked.amount_into_gas_mol, [0.0; 2]);
        assert_eq!(blocked.gross_amount_into_gas_mol, [0.0; 2]);
        assert_eq!(blocked.gross_amount_out_of_gas_mol, [0.0; 2]);
        assert_eq!(blocked.actual_sorbed_amount_change_mol, 0.0);
        assert_eq!(film(result.state(), "B"), &[12345.0; 4]);
        assert_eq!(report.source_relative_gas_energy_residual_j, [0.0; 2]);
        for side in 0..2 {
            near(
                result.state().gases()[side]
                    .species_mass_kg(&id("B"))
                    .unwrap(),
                state.gases()[side].species_mass_kg(&id("B")).unwrap(),
                1e-20,
            );
            assert_eq!(
                result.state().gases()[side]
                    .point()
                    .axes()
                    .get("specimen-lot"),
                Some(&17.0)
            );
        }
        let values = (
            result.state().gases().each_ref().map(|g| g.mass_kg()),
            result
                .state()
                .gases()
                .each_ref()
                .map(|g| g.state().state().pressure),
            film(result.state(), "A").to_vec(),
            report.thermostat_heat_into_gas_j,
            report.evaluations,
        );
        if let Some(expected) = &expected {
            assert_eq!(&values, expected);
        } else {
            expected = Some(values);
        }
        if reference > 0.0 {
            assert!(report.enthalpy_into_gas_j[0].abs() > 1e12);
        }
    }
}

#[test]
fn accepted_private_prefix_work_exhaustion_cancellation_and_cell_bound_are_atomic() {
    let card = membrane("A", D, S, 4e5);
    let context = point(1e5);
    let sources = [channel("A", &card, &context, 1, 1e5)];
    let volume = equal_capacity_volume();
    let gases = [
        gas(&[1.0], 0.0, 1.4e5, volume),
        gas(&[1.0], 0.0, 0.6e5, volume),
    ];
    let state = IsothermalGasMembraneMixture::new(
        T,
        AREA,
        &sources,
        gases.clone(),
        native_controls(),
        || false,
    )
    .unwrap();
    let saved = state.finite_state(&id("A")).unwrap().clone();
    let mut work = state.adaptive_workspace().unwrap();
    let mut c = controls(0.5, 0.1);
    c.maximum_evaluations = 6;
    let error = state
        .advance_candidate(5.0, c, &mut work, || false)
        .unwrap_err();
    match error {
        MembraneMixtureAdvanceError::WorkLimit {
            time_s,
            evaluations,
        } => {
            assert!(
                time_s > state.elapsed_s(),
                "failure must follow accepted private intervals"
            );
            assert_eq!(evaluations, 6);
        }
        error => panic!("expected late work refusal, got {error}"),
    }
    let mut polls = 0;
    assert!(matches!(
        state.advance_candidate(5.0, controls(0.5, 0.1), &mut work, || {
            polls += 1;
            polls == 30
        }),
        Err(MembraneMixtureAdvanceError::Cancelled)
    ));
    let mut c = controls(0.5, 0.1);
    c.maximum_total_cells = 0;
    assert!(matches!(
        state.advance_candidate(5.0, c, &mut work, || false),
        Err(MembraneMixtureAdvanceError::Invalid(_))
    ));
    assert_eq!(state.gases(), &gases);
    assert_eq!(state.finite_state(&id("A")).unwrap(), &saved);
    assert_eq!(state.elapsed_s(), 0.0);
}

#[test]
fn selective_loading_hits_original_total_pressure_support_without_partial_publication() {
    // Both baths and initial cells total 100 kPa. A fills the film rapidly
    // while B desorbs slowly, so a later INTERIOR total leaves 110 kPa support.
    let a = membrane("fast A", 1e-9, S, 110000.0);
    let b = membrane("slow B", 1e-12, 1e-5, 110000.0);
    let context = point(1e5);
    let sources = [
        channel("A", &a, &context, 4, 0.0),
        channel("B", &b, &context, 4, 1e5),
    ];
    let state = IsothermalGasMembraneMixture::new(
        T,
        AREA,
        &sources,
        [
            gas(&[0.9, 0.1], 0.0, 1e5, 1e-6),
            gas(&[0.1, 0.9], 0.0, 1e5, 1e-6),
        ],
        native_controls(),
        || false,
    )
    .unwrap();
    let mut work = state.adaptive_workspace().unwrap();
    let warm = state
        .advance_candidate(0.005, controls(0.005, 0.1), &mut work, || false)
        .unwrap()
        .into_parts()
        .0;
    assert_eq!(warm.elapsed_s(), 0.005);
    assert!(film(&warm, "A").iter().any(|&p| p > 0.0));
    let saved_gases = warm.gases().clone();
    let saved_a = warm.finite_state(&id("A")).unwrap().clone();
    let saved_b = warm.finite_state(&id("B")).unwrap().clone();
    let error = warm
        .advance_candidate(1.0, controls(0.005, 0.1), &mut work, || false)
        .unwrap_err();
    assert!(
        matches!(
            error,
            MembraneMixtureAdvanceError::Native(IsothermalGasMembraneError::Source(_))
        ),
        "{error}"
    );
    assert_eq!(warm.gases(), &saved_gases);
    assert_eq!(warm.finite_state(&id("A")).unwrap(), &saved_a);
    assert_eq!(warm.finite_state(&id("B")).unwrap(), &saved_b);
    assert_eq!(warm.elapsed_s(), 0.005);
}
