//! Held-boundary heat exchange must retain the original caloric coordinate.
//!
//! G1 references integrate affine source Cp/Cv independently; G3 checks physical
//! amount scaling. These are synthetic material-card fixtures, not measured data.

mod support;

#[path = "lumped_caloric_exchange/forced.rs"]
mod forced;

#[path = "lumped_caloric_exchange/pair.rs"]
mod pair;

use fs_conduction::ConductionError;
use fs_conduction::lumped::{
    BiotGate,
    exchange::caloric::{
        LumpedCaloricExchangeConfig, solve_lumped_caloric_exchange as solve,
    },
};
use fs_evidence::ValidityDomain;
use fs_material::thermal::transport::{
    conduction::CaloricConductionLayer,
    isobaric::{
        IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery, conduction::IsobaricConductionLaw,
    },
    isochoric::{
        IsochoricHeatCapacityLaw, IsochoricHeatCapacityQuery, conduction::IsochoricConductionLaw,
    },
    storage::{
        HeatStorageAmount, HeatStorageComponentQuery, HeatStorageThermalLaw, PreparedHeatStorage,
        conduction::HeatStorageConductionLaw,
    },
};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, MaterialId, Reservoir, ThermalStorageBasis,
};
use support::{with_cancelled_cx, with_cx};

// SI dimensions for the synthetic source properties; no private material API.
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CONDUCTIVITY: Dims = Dims([1, 1, -3, -1, 0, 0]);
const RHO: f64 = 1000.0;
const GATE: BiotGate = BiotGate::corpus_default();
const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

// The same original material cards used by fs-material's integrated conduction
// tests: Cv=200+0.5T at fixed density, Cp=400+T at fixed pressure, each with
// independent conductivity knots. No replacement ThermalLaw is implemented here.
struct Fixture {
    card: MaterialCard,
    capacity: PropertyKey,
    conductivity: PropertyKey,
    lower: QueryPoint,
    upper: QueryPoint,
}

impl Fixture {
    fn new(isochoric: bool) -> Self {
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let (axis, quantity, fixed) = if isochoric {
            ("rho", QuantitySpec::dimensional(DENSITY), RHO)
        } else {
            ("p", QuantitySpec::dimensional(Pressure::DIMS), 100000.0)
        };
        let capacity = PropertyKey::new(
            if isochoric { "source_cv" } else { "source_cp" },
            SPECIFIC_HEAT,
        );
        let conductivity = PropertyKey::new("normal_conductivity", CONDUCTIVITY);
        let validity = ValidityDomain::unconstrained()
            .with_quantity("T", temperature, 300.0, 700.0)
            .with_quantity(axis, quantity, fixed, fixed);
        let capacity_temperatures = if isochoric {
            [300.0, 440.0, 700.0]
        } else {
            [300.0, 475.0, 700.0]
        };
        let conductivity_temperatures = if isochoric {
            [300.0, 450.0, 600.0, 700.0]
        } else {
            [300.0, 350.0, 550.0, 700.0]
        };
        let capacity_knots = capacity_temperatures
            .into_iter()
            .map(|t| {
                (
                    t,
                    if isochoric {
                        200.0 + 0.5 * t
                    } else {
                        400.0 + t
                    },
                )
            })
            .collect();
        let conductivity_knots = conductivity_temperatures
            .into_iter()
            .map(|t| (t, if isochoric { 0.04 * t } else { 0.02 * t }))
            .collect();
        let mut claims = ClaimSet::new();
        for (key, dims, knots) in [
            (capacity.clone(), SPECIFIC_HEAT, capacity_knots),
            (conductivity.clone(), CONDUCTIVITY, conductivity_knots),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key,
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots,
                        dims,
                    },
                    validity: validity.clone(),
                    interpolation: InterpolationPolicy::LinearInside,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic independent caloric and conductivity paths".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: if isochoric {
                    "Cv fixture".into()
                } else {
                    "Cp fixture".into()
                },
                phase: "single phase".into(),
                process: "integrated conduction numerical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let point = |t| {
            QueryPoint::new()
                .with_quantity("T", temperature, t)
                .unwrap()
                .with_quantity(axis, quantity, fixed)
                .unwrap()
        };
        Self {
            card,
            capacity,
            conductivity,
            lower: point(300.0),
            upper: point(700.0),
        }
    }

    fn cv(&self) -> IsochoricHeatCapacityLaw {
        IsochoricHeatCapacityLaw::try_new(
            &self.card,
            IsochoricHeatCapacityQuery {
                specific_heat_cv: &self.capacity,
                temperature_axis: "T",
                density_axis: "rho",
                lower: &self.lower,
                upper: &self.upper,
                selection: SELECTION,
            },
            300.0,
        )
        .unwrap()
    }

    fn cp(&self) -> IsobaricHeatCapacityLaw {
        IsobaricHeatCapacityLaw::try_new(
            &self.card,
            IsobaricHeatCapacityQuery {
                specific_heat_cp: &self.capacity,
                temperature_axis: "T",
                pressure_axis: "p",
                lower: &self.lower,
                upper: &self.upper,
                selection: SELECTION,
            },
            300.0,
        )
        .unwrap()
    }
}


fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

fn config(duration_s: f64) -> LumpedCaloricExchangeConfig {
    LumpedCaloricExchangeConfig {
        duration_s,
        characteristic_length_m: 0.001,
        energy_tolerance_j: 1e-7,
    }
}

fn cell(specific_energy_j_kg: f64, ambient_k: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(42),
        segment_id: 7,
        body_id: 3,
        material_id: MaterialId(11),
        mass_kg: 0.4,
        specific_h_j_kg: specific_energy_j_kg,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, 0.0],
            end_m: [0.02, 0.0, 0.0],
            volume_m3: 0.4 / RHO,
            surface_area_m2: 0.25,
        },
        boundary: Boundary {
            ambient_k,
            htc_w_m2_k: 8.0, // G=2 W/K before a support covers any area.
            reservoir: None,
            source_w: 0.0,
        },
    }
}

// Independent primitives of the admitted affine source laws, with zero energy
// at 300 K. The logarithmic formula is used only for resolved temperature gaps.
fn energy(isochoric: bool, temperature_k: f64) -> f64 {
    let x = temperature_k - 300.0;
    let factor = if isochoric { 1.0 } else { 2.0 };
    factor * (350.0 * x + 0.25 * x * x)
}

fn exposure(isochoric: bool, initial_k: f64, final_k: f64, bath_k: f64) -> f64 {
    let factor = if isochoric { 1.0 } else { 2.0 };
    factor
        * ((200.0 + 0.5 * bath_k) * ((bath_k - initial_k) / (bath_k - final_k)).ln()
            - 0.5 * (final_k - initial_k))
}

#[test]
fn g1_affine_cv_and_cp_use_native_energy_over_finite_heating_and_cooling() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_conduction =
        IsochoricConductionLaw::try_new(&cv_source.card, &cv, &cv_source.conductivity, SELECTION)
            .unwrap();
    let cp_conduction =
        IsobaricConductionLaw::try_new(&cp_source.card, &cp, &cp_source.conductivity, SELECTION)
            .unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_conduction).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_conduction).unwrap();

    with_cx(|cx| {
        for (layer, isochoric) in [(&cv_layer, true), (&cp_layer, false)] {
            for (initial_k, final_k, bath_k) in [(360.0, 600.0, 680.0), (650.0, 380.0, 320.0)] {
                let body = cell(energy(isochoric, initial_k), bath_k);
                let original = body.clone();
                let conductance = body.boundary.htc_w_m2_k * body.geometry.surface_area_m2;
                let duration = body.mass_kg
                    * exposure(isochoric, initial_k, final_k, bath_k)
                    / conductance;
                let result = solve(cx, &body, layer, GATE, config(duration)).unwrap();
                let expected_heat =
                    body.mass_kg * (energy(isochoric, final_k) - energy(isochoric, initial_k));
                assert_eq!(
                    result.caloric.storage_basis,
                    if isochoric {
                        ThermalStorageBasis::InternalEnergy
                    } else {
                        ThermalStorageBasis::Enthalpy
                    }
                );
                assert_eq!(result.cell_id, body.id.0);
                assert_eq!(
                    result.caloric_identity,
                    layer.caloric_law().caloric_identity().unwrap()
                );
                near(result.caloric.temperature_k, final_k, 2e-9);
                near(
                    result.caloric.final_specific_energy_j_kg,
                    energy(isochoric, final_k),
                    3e-7,
                );
                near(result.heat_into_body_j, expected_heat, 2e-7);
                near(result.stored_energy_change_j, expected_heat, 2e-7);
                near(result.ambient_heat_j, expected_heat, 2e-7);
                assert_eq!(result.reservoir_heat_j, 0.0);
                near(
                    result.caloric.mean_temperature_k,
                    bath_k - expected_heat / (conductance * duration),
                    2e-9,
                );
                assert!(result.heat_residual_j.abs() <= config(duration).energy_tolerance_j);
                near(
                    result.maximum_biot,
                    8.0 * 0.001 / (if isochoric { 12.0 } else { 6.0 }),
                    1e-15,
                );
                assert_eq!(body, original);

                // A representative-temperature exponential is measurably wrong.
                let initial_capacity = if isochoric {
                    200.0 + 0.5 * initial_k
                } else {
                    400.0 + initial_k
                };
                let frozen = bath_k
                    + (initial_k - bath_k)
                        * (-conductance * duration / (body.mass_kg * initial_capacity)).exp();
                assert!((frozen - final_k).abs() > 5.0);

                // Exact held-bath exchange has no time-step discretization error.
                let first = solve(cx, &body, layer, GATE, config(duration * 0.5)).unwrap();
                let mut midpoint = body.clone();
                midpoint.specific_h_j_kg = first.caloric.final_specific_energy_j_kg;
                let second = solve(cx, &midpoint, layer, GATE, config(duration * 0.5)).unwrap();
                near(second.caloric.temperature_k, result.caloric.temperature_k, 2e-9);
                near(
                    first.heat_into_body_j + second.heat_into_body_j,
                    result.heat_into_body_j,
                    2e-7,
                );
            }
        }
    });
}

#[test]
fn g1_two_baths_exclude_support_area_and_retain_opposing_heat_at_equilibrium() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let mut body = cell(energy(true, 650.0), 320.0);
    body.boundary.reservoir = Some(Reservoir {
        area_m2: 0.0625,
        conductance_w_k: 0.5,
        temperature_k: 680.0,
    });
    // Ambient G=8*(.25-.0625)=1.5 W/K, support G=.5 W/K, bath=410 K.
    let duration = body.mass_kg * exposure(true, 650.0, 500.0, 410.0) / 2.0;
    let heat = body.mass_kg * (energy(true, 500.0) - energy(true, 650.0));
    let mean = 410.0 - heat / (2.0 * duration);

    with_cx(|cx| {
        let result = solve(cx, &body, &layer, GATE, config(duration)).unwrap();
        assert_eq!(result.total_conductance_w_k, 2.0);
        assert_eq!(result.effective_reservoir_temperature_k, 410.0);
        near(result.caloric.temperature_k, 500.0, 2e-9);
        near(result.heat_into_body_j, heat, 1e-7);
        near(result.caloric.mean_temperature_k, mean, 2e-9);
        near(result.ambient_heat_j, 1.5 * duration * (320.0 - mean), 1e-7);
        near(result.reservoir_heat_j, 0.5 * duration * (680.0 - mean), 1e-7);
        assert!(result.ambient_heat_j < 0.0 && result.reservoir_heat_j > 0.0);
        near(
            result.ambient_heat_j + result.reservoir_heat_j,
            result.heat_into_body_j,
            config(duration).energy_tolerance_j,
        );
        assert!(result.heat_partition_residual_j.abs() <= config(duration).energy_tolerance_j);

        // The body can be stationary while finite heat passes between the baths.
        body.specific_h_j_kg = energy(true, 410.0);
        let steady = solve(cx, &body, &layer, GATE, config(80.0)).unwrap();
        assert_eq!(steady.caloric.temperature_k, 410.0);
        assert_eq!(steady.caloric.final_specific_energy_j_kg, body.specific_h_j_kg);
        assert_eq!(steady.heat_into_body_j, 0.0);
        assert_eq!(steady.stored_energy_change_j, 0.0);
        near(steady.ambient_heat_j, -10800.0, 1e-9);
        near(steady.reservoir_heat_j, 10800.0, 1e-9);
        assert_eq!(steady.heat_partition_residual_j, 0.0);
    });
}

#[test]
fn g1_native_energy_bath_and_duration_gradients_match_independent_perturbations() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let mut body = cell(energy(true, 600.0), 320.0);
    body.boundary.reservoir = Some(Reservoir {
        area_m2: 0.0625,
        conductance_w_k: 0.5,
        temperature_k: 680.0,
    });
    with_cx(|cx| {
        let initial = body.specific_h_j_kg;
        let duration = 100.0;
        let result = solve(cx, &body, &layer, GATE, config(duration)).unwrap();
        let heat = |e: f64, ambient: f64, reservoir: f64, time: f64| {
            let mut perturbed = body.clone();
            perturbed.specific_h_j_kg = e;
            perturbed.boundary.ambient_k = ambient;
            perturbed.boundary.reservoir.as_mut().unwrap().temperature_k = reservoir;
            solve(cx, &perturbed, &layer, GATE, config(time))
                .unwrap()
                .heat_into_body_j
        };
        let energy_step = 0.5;
        let temperature_step = 0.01;
        let time_step = 0.01;
        let differences = [
            (heat(initial + energy_step, 320.0, 680.0, duration)
                - heat(initial - energy_step, 320.0, 680.0, duration))
                / (2.0 * energy_step),
            (heat(initial, 320.0 + temperature_step, 680.0, duration)
                - heat(initial, 320.0 - temperature_step, 680.0, duration))
                / (2.0 * temperature_step),
            (heat(initial, 320.0, 680.0 + temperature_step, duration)
                - heat(initial, 320.0, 680.0 - temperature_step, duration))
                / (2.0 * temperature_step),
            (heat(initial, 320.0, 680.0, duration + time_step)
                - heat(initial, 320.0, 680.0, duration - time_step))
                / (2.0 * time_step),
        ];
        let analytic = [
            result.heat_initial_energy_derivative_kg,
            result.heat_ambient_temperature_derivative_j_k,
            result.heat_reservoir_temperature_derivative_j_k,
            result.heat_duration_derivative_w,
        ];
        for (derivative, measured) in analytic.into_iter().zip(differences) {
            near(derivative, measured, 5e-7 * measured.abs().max(1.0));
        }
        assert!(result.heat_initial_energy_derivative_kg < 0.0);
        assert!(result.heat_ambient_temperature_derivative_j_k > 0.0);
        assert!(result.heat_reservoir_temperature_derivative_j_k > 0.0);
        near(
            result.heat_ambient_temperature_derivative_j_k,
            3.0 * result.heat_reservoir_temperature_derivative_j_k,
            1e-9,
        );
        let final_k = result.caloric.temperature_k;
        near(
            result.heat_duration_derivative_w,
            1.5 * (320.0 - final_k) + 0.5 * (680.0 - final_k),
            1e-8,
        );
    });
}

#[test]
fn g3_assembled_capacity_keeps_reference_mass_and_scales_joules_once() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    // .25*Cv + .75*Cp = 1.75*Cv J/K for the one-kilogram assembly.
    // Scaling all physical amounts and boundary area leaves T(t) unchanged.
    for scale in [1.0, 4.0] {
        let storage = PreparedHeatStorage::try_new(&[
            HeatStorageComponentQuery {
                id: 1,
                card: &cv_source.card,
                capacity: &cv_source.capacity,
                amount: HeatStorageAmount::MassKg(0.25 * scale),
                temperature_axis: "T",
                lower: &cv_source.lower,
                upper: &cv_source.upper,
                selection: SELECTION,
            },
            HeatStorageComponentQuery {
                id: 2,
                card: &cp_source.card,
                capacity: &cp_source.capacity,
                amount: HeatStorageAmount::MassKg(0.75 * scale),
                temperature_axis: "T",
                lower: &cp_source.lower,
                upper: &cp_source.upper,
                selection: SELECTION,
            },
        ])
        .unwrap();
        let caloric = HeatStorageThermalLaw::try_new(&storage, scale, 300.0).unwrap();
        let conduction = HeatStorageConductionLaw::try_new(
            &cv_source.card,
            &caloric,
            &cv_source.conductivity,
            "T",
            &cv_source.lower,
            &cv_source.upper,
            SELECTION,
        )
        .unwrap();
        let layer = CaloricConductionLayer::from_heat_storage(&conduction).unwrap();
        let mut body = cell(1.75 * energy(true, 360.0), 680.0);
        body.mass_kg = scale;
        body.geometry.volume_m3 = scale / RHO;
        body.geometry.surface_area_m2 = 0.25 * scale;
        let duration = 1.75 * exposure(true, 360.0, 600.0, 680.0) / 2.0;
        with_cx(|cx| {
            let result = solve(cx, &body, &layer, GATE, config(duration)).unwrap();
            assert_eq!(result.caloric.storage_basis, ThermalStorageBasis::DeclaredSensibleStorage);
            near(result.caloric.temperature_k, 600.0, 2e-9);
            near(result.caloric.final_specific_energy_j_kg, 1.75 * energy(true, 600.0), 3e-7);
            near(
                result.heat_into_body_j,
                scale * 1.75 * (energy(true, 600.0) - energy(true, 360.0)),
                5e-7,
            );
            near(result.total_conductance_w_k, 2.0 * scale, 1e-14);
            let mut wrong_mass = body.clone();
            wrong_mass.mass_kg *= 1.01;
            assert!(solve(cx, &wrong_mass, &layer, GATE, config(duration)).is_err());
        });
    }
}

#[test]
fn g1_zero_exposure_and_equal_baths_preserve_the_original_state() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    with_cx(|cx| {
        let body = cell(energy(true, 500.0), 680.0);
        let zero_time = solve(cx, &body, &layer, GATE, config(0.0)).unwrap();
        assert_eq!(zero_time.caloric.final_specific_energy_j_kg, body.specific_h_j_kg);
        assert_eq!(zero_time.heat_into_body_j, 0.0);
        near(zero_time.heat_duration_derivative_w, 2.0 * (680.0 - 500.0), 1e-12);

        let mut covered = body.clone();
        covered.boundary.htc_w_m2_k = 1e6;
        covered.boundary.reservoir = Some(Reservoir {
            area_m2: covered.geometry.surface_area_m2,
            conductance_w_k: 0.0,
            temperature_k: 900.0,
        });
        let isolated = solve(cx, &covered, &layer, GATE, config(1e6)).unwrap();
        assert_eq!(isolated.maximum_biot, 0.0);
        assert_eq!(isolated.total_conductance_w_k, 0.0);
        assert_eq!(isolated.caloric.final_specific_energy_j_kg, body.specific_h_j_kg);
        assert_eq!(isolated.heat_into_body_j, 0.0);
        assert_eq!(isolated.ambient_heat_j, 0.0);
        assert_eq!(isolated.reservoir_heat_j, 0.0);

        let mut equilibrium = body;
        equilibrium.boundary.ambient_k = 500.0;
        equilibrium.boundary.reservoir = Some(Reservoir {
            area_m2: 0.0625,
            conductance_w_k: 0.7,
            temperature_k: 500.0,
        });
        let same = solve(cx, &equilibrium, &layer, GATE, config(1000.0)).unwrap();
        assert_eq!(same.effective_reservoir_temperature_k, 500.0);
        assert_eq!(same.caloric.final_specific_energy_j_kg, equilibrium.specific_h_j_kg);
        assert_eq!(same.heat_into_body_j, 0.0);
        assert_eq!(same.ambient_heat_j, 0.0);
        assert_eq!(same.reservoir_heat_j, 0.0);
    });
}

#[test]
fn g0_original_density_and_applied_source_constraints_are_not_bypassed() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let body = cell(energy(true, 500.0), 680.0);
    with_cx(|cx| {
        for case in 0..4 {
            let mut invalid = body.clone();
            match case {
                0 => invalid.mass_kg *= 1.1,
                1 => invalid.geometry.volume_m3 *= 1.1,
                2 => invalid.boundary.source_w = 1.0,
                _ => invalid.boundary.source_w = -1.0,
            }
            // Even a zero-duration proposal must honor the original owner.
            for duration in [0.0, 20.0] {
                assert!(solve(cx, &invalid, &layer, GATE, config(duration)).is_err());
            }
        }
        assert!(solve(cx, &body, &layer, GATE, config(20.0)).is_ok());
    });
}

#[test]
fn g0_biot_uses_original_minimum_conductivity_and_active_support_coefficient() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let body = cell(energy(true, 650.0), 600.0);
    with_cx(|cx| {
        // k(650)=26 would admit this body, but original minimum k(300)=12
        // gives Bi=.1333... and correctly refuses the declared uniform rung.
        let mut cfg = config(10.0);
        cfg.characteristic_length_m = 0.2;
        assert!(8.0 * cfg.characteristic_length_m / 26.0 < GATE.ceiling());
        assert!(8.0 * cfg.characteristic_length_m / 12.0 > GATE.ceiling());
        assert!(solve(cx, &body, &layer, GATE, cfg).is_err());

        let mut support = body.clone();
        support.boundary.htc_w_m2_k = 0.0;
        support.boundary.reservoir = Some(Reservoir {
            area_m2: 0.00025,
            conductance_w_k: 0.5,
            temperature_k: 600.0,
        });
        // The active support has h=2000 W/(m2 K) even with no ambient film.
        assert!(solve(cx, &support, &layer, GATE, config(10.0)).is_err());
    });
}

#[test]
fn g0_material_support_limits_the_traversed_state_not_the_external_bath() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    with_cx(|cx| {
        let body = cell(energy(true, 500.0), 900.0);
        let short = solve(cx, &body, &layer, GATE, config(1.0)).unwrap();
        assert!((500.0..700.0).contains(&short.caloric.temperature_k));
        assert!(solve(cx, &body, &layer, GATE, config(1000.0)).is_err());
        for invalid in [
            cell(energy(true, 700.0) + 1.0, 600.0),
            cell(energy(true, 700.0), 900.0),
            cell(energy(true, 300.0), 200.0),
        ] {
            assert!(solve(cx, &invalid, &layer, GATE, config(1.0)).is_err());
        }
    });
}

#[test]
fn g1_weak_heat_is_retained_separately_from_native_energy_rounding() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let body = cell(energy(true, 500.0), 680.0);
    let duration = 1e-9;
    // Independent second-order ODE expansion; its next relative term is below
    // 1e-21 here. Subtracting two native 80000 J/kg energies loses heat accuracy.
    let attenuation = 2.0 * duration / (body.mass_kg * 450.0);
    let expected = 2.0 * duration * 180.0 * (1.0 - 0.5 * attenuation);
    with_cx(|cx| {
        let result = solve(cx, &body, &layer, GATE, config(duration)).unwrap();
        near(result.heat_into_body_j, expected, 128.0 * f64::EPSILON * expected);
        assert!(result.heat_into_body_j > 0.0);
        assert_ne!(result.stored_energy_change_j, result.heat_into_body_j);
        assert!(result.heat_residual_j.abs() > 1e-13);
        near(
            result.stored_energy_change_j,
            body.mass_kg
                * (result.caloric.final_specific_energy_j_kg - body.specific_h_j_kg),
            1e-20,
        );
        // A declared budget below the representable endpoint mismatch refuses.
        let mut strict = config(duration);
        strict.energy_tolerance_j = 1e-13;
        assert!(solve(cx, &body, &layer, GATE, strict).is_err());
    });
}

#[test]
fn g1_rounded_bath_equilibrium_retains_original_boundary_heat() {
    let source = Fixture::new(true);
    // Use the original source with its datum at an interior temperature.
    // Small signed heat is representable in native energy on both sides of
    // zero even when the conductance-weighted bath temperature rounds to 500 K.
    let initial_k: f64 = 500.0;
    let cv = IsochoricHeatCapacityLaw::try_new(
        &source.card,
        IsochoricHeatCapacityQuery {
            specific_heat_cv: &source.capacity,
            temperature_axis: "T",
            density_axis: "rho",
            lower: &source.lower,
            upper: &source.upper,
            selection: SELECTION,
        },
        initial_k,
    )
    .unwrap();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    with_cx(|cx| {
        for support_k in [initial_k.next_up(), initial_k.next_down()] {
            let mut body = cell(0.0, initial_k);
            body.boundary.reservoir = Some(Reservoir {
                area_m2: 0.125,
                conductance_w_k: 1.0,
                temperature_k: support_k,
            });
            // Ambient G=8*(.25-.125)=1 W/K and support G=1 W/K.
            // The original boundary has signed nonzero power, although its
            // rounded equivalent bath is indistinguishable from the body.
            assert_eq!(0.5 * (initial_k + support_k), initial_k);
            let original_power_w = support_k - initial_k;
            assert_ne!(original_power_w, 0.0);
            let duration = 1.0;
            let attenuation = 2.0 * duration / (body.mass_kg * 450.0);
            let retained_fraction = -(-attenuation).exp_m1() / attenuation;
            let expected_heat = original_power_w * duration * retained_fraction;
            let result = solve(cx, &body, &layer, GATE, config(duration)).unwrap();
            assert_eq!(result.caloric.temperature_k, initial_k);
            assert_ne!(result.caloric.final_specific_energy_j_kg, 0.0);
            assert_eq!(
                result.heat_into_body_j.is_sign_positive(),
                original_power_w.is_sign_positive()
            );
            near(
                result.heat_into_body_j,
                expected_heat,
                256.0 * f64::EPSILON * expected_heat.abs(),
            );
            near(
                result.caloric.final_specific_energy_j_kg,
                expected_heat / body.mass_kg,
                256.0 * f64::EPSILON * (expected_heat / body.mass_kg).abs(),
            );
            near(
                result.heat_duration_derivative_w,
                original_power_w * (-attenuation).exp(),
                256.0 * f64::EPSILON * original_power_w.abs(),
            );

            let stopped = solve(cx, &body, &layer, GATE, config(0.0)).unwrap();
            assert_eq!(stopped.total_conductance_w_k, 2.0);
            assert_eq!(stopped.heat_duration_derivative_w, original_power_w);
            assert_eq!(stopped.effective_reservoir_temperature_k, initial_k);
            assert_eq!(stopped.caloric.final_specific_energy_j_kg, 0.0);
            assert_eq!(stopped.heat_into_body_j, 0.0);
            assert_eq!(stopped.stored_energy_change_j, 0.0);
            assert_eq!(stopped.ambient_heat_j, 0.0);
            assert_eq!(stopped.reservoir_heat_j, 0.0);
        }
    });
}

#[test]
fn g0_cancelled_proposal_leaves_the_original_cell_reusable() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let body = cell(energy(true, 500.0), 680.0);
    let original = body.clone();
    let before = with_cx(|cx| solve(cx, &body, &layer, GATE, config(20.0)).unwrap());
    with_cancelled_cx(|cx| {
        assert!(matches!(
            solve(cx, &body, &layer, GATE, config(20.0)),
            Err(ConductionError::Cancelled { .. })
        ));
    });
    assert_eq!(body, original);
    let after = with_cx(|cx| solve(cx, &body, &layer, GATE, config(20.0)).unwrap());
    assert_eq!(before, after);
}
