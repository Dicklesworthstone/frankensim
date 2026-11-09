use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::thermal::transport::isochoric::IsochoricHeatCapacityQuery;
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_qty::{Dims, QuantitySpec, Temperature};

const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CONDUCTIVITY: Dims = Dims([1, 1, -3, -1, 0, 0]);

const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Fixture {
    card: MaterialCard,
    cv: PropertyKey,
    k: PropertyKey,
    lower: QueryPoint,
    upper: QueryPoint,
}
impl Fixture {
    // Synthetic source, not measured liquid data: Cv(T)=1e-5*T, k=constant.
    fn new(conductivities: [f64; 3], density: f64) -> Self {
        let cv = PropertyKey::new("jet_test_cv", SPECIFIC_HEAT);
        let k = PropertyKey::new("jet_test_conductivity", CONDUCTIVITY);
        let temperature = QuantitySpec::dimensional(Temperature::DIMS);
        let rho = QuantitySpec::dimensional(DENSITY);
        let validity = ValidityDomain::unconstrained()
            .with_quantity("T", temperature, 280.0, 360.0)
            .with_quantity("rho", rho, density, density);
        let mut claims = ClaimSet::new();
        for (key, dims, values) in [
            (cv.clone(), SPECIFIC_HEAT, [0.0028, 0.0032, 0.0036]),
            (k.clone(), CONDUCTIVITY, conductivities),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key,
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: [280.0, 320.0, 360.0].into_iter().zip(values).collect(),
                        dims,
                    },
                    validity: validity.clone(),
                    interpolation: InterpolationPolicy::LinearInside,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic coupled-jet caloric oracle".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic numerical liquid".into(),
                phase: "single incompressible phase".into(),
                process: "fixed-density test".into(),
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
                .with_quantity("rho", rho, density)
                .unwrap()
        };
        Self {
            card,
            cv,
            k,
            lower: point(280.0),
            upper: point(360.0),
        }
    }
    fn caloric(&self, datum: f64) -> IsochoricHeatCapacityLaw {
        IsochoricHeatCapacityLaw::try_new(
            &self.card,
            IsochoricHeatCapacityQuery {
                specific_heat_cv: &self.cv,
                temperature_axis: "T",
                density_axis: "rho",
                lower: &self.lower,
                upper: &self.upper,
                selection: SELECTION,
            },
            datum,
        )
        .unwrap()
    }
    fn conduction<'a>(&self, caloric: &'a IsochoricHeatCapacityLaw) -> IsochoricConductionLaw<'a> {
        IsochoricConductionLaw::try_new(&self.card, caloric, &self.k, SELECTION).unwrap()
    }
}

fn material() -> CarreauSlenderJetParameters {
    CarreauSlenderJetParameters {
        density_kg_m3: 2.0,
        surface_tension_n_m: 0.0,
        temperature_k: 300.0,
        rheology: ThermalCarreauYasuda::new(
            CarreauYasuda::new(0.8, 0.02, 10.0, 0.5, 2.0, [0.0, 1000.0]).unwrap(),
            TemperatureShift::new(
                TemperatureShiftModel::WilliamsLandelFerry {
                    c1: 1.0,
                    c2_k: 100.0,
                },
                300.0,
                [280.0, 360.0],
            )
            .unwrap(),
            ViscosityShiftMode::VerticalAndRate,
        ),
    }
}
fn limits() -> SlenderJetLimits {
    SlenderJetLimits {
        minimum_radius_m: 0.1,
        maximum_abs_slope: 0.75,
        minimum_radius_cells: 2.0,
    }
}
fn flow(velocities: &[f64], gate: &CancelGate) -> PeriodicSlenderJet {
    PeriodicSlenderJet::new_carreau(material(), 4.0, &[1.0; 16], velocities, limits(), gate)
        .unwrap()
}
fn alternating(gate: &CancelGate) -> PeriodicSlenderJet {
    let velocities: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    flow(&velocities, gate)
}
fn temperatures() -> Vec<f64> {
    (0..16)
        .map(|i| if i % 2 == 0 { 300.0 } else { 330.0 })
        .collect()
}
fn budget(dt: f64) -> ThermalJetStepBudget {
    ThermalJetStepBudget {
        flow: SlenderJetStepBudget {
            maximum_step_s: dt,
            cfl: 0.2,
            maximum_substeps: 10000,
        },
        energy_tolerance_j: 1e-10,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} vs {b}, tolerance {tolerance}"
    );
}

// Independent reduced conservative ODE: even/odd cell masses change, but the
// staggered dual-cell momentum advection is zero. No production material,
// transport or heat routines enter this reference.
fn reference_rhs(y: [f64; 5], k: f64) -> [f64; 5] {
    let [b, v, ee, eo, _heat] = y;
    let dx = 0.25;
    let rho = 2.0;
    let factor = rho * core::f64::consts::PI * dx;
    let ae = 1.0 - b;
    let ao = 1.0 + b;
    let ue = ee / (factor * ae);
    let uo = eo / (factor * ao);
    let te = (300.0 * 300.0 + 2.0 * ue / 1e-5).sqrt();
    let to = (300.0 * 300.0 + 2.0 * uo / 1e-5).sqrt();
    let mu = |t: f64| {
        let d = t - 300.0;
        let shift = 10.0_f64.powf(-d / (100.0 + d));
        let rate = shift * 2.0 * 3.0_f64.sqrt() * v.abs() / dx;
        shift * (0.02 + 0.78 * (1.0 + (10.0 * rate).powi(2)).powf(-0.25))
    };
    let me = mu(te);
    let mo = mu(to);
    let pe = 12.0 * core::f64::consts::PI * ae * me * v * v / dx;
    let po = 12.0 * core::f64::consts::PI * ao * mo * v * v / dx;
    // Positive v sends mass from even cells through BOTH adjacent faces.
    assert!(v > 0.0);
    let advected = 2.0 * rho * core::f64::consts::PI * v * ue;
    let harmonic = 2.0 * k * ae * ao / (ae + ao);
    let conducted = 2.0 * core::f64::consts::PI * harmonic * (to - te) / dx;
    [
        2.0 * v / dx,
        -6.0 * (ae * me + ao * mo) * v / (rho * dx * dx),
        -advected + pe + conducted,
        advected + po - conducted,
        8.0 * (pe + po),
    ]
}
fn reference(k: f64) -> [f64; 5] {
    let factor = 2.0 * core::f64::consts::PI * 0.25;
    let mut y = [
        0.0,
        0.125,
        0.0,
        factor * 0.5e-5 * (330.0 * 330.0 - 300.0 * 300.0),
        0.0,
    ];
    let dt = 0.02 / 10000.0;
    for _ in 0..10000 {
        let k1 = reference_rhs(y, k);
        let k2 = reference_rhs(std::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]), k);
        let k3 = reference_rhs(std::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]), k);
        let k4 = reference_rhs(std::array::from_fn(|i| y[i] + dt * k3[i]), k);
        y = std::array::from_fn(|i| y[i] + dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}

#[test]
fn g1_coupled_mass_momentum_and_nonlinear_calorics_converge_at_third_order() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let conductivity = fixture.conduction(&caloric);
    for conduction in [None, Some(&conductivity)] {
        let expected = reference(if conduction.is_some() { 1e-5 } else { 0.0 });
        let mut errors = Vec::new();
        for dt in [0.001, 0.0005, 0.00025] {
            let mut jet = ThermalCarreauJet::new(
                alternating(&gate),
                &caloric,
                conduction,
                &temperatures(),
                100.0,
                &gate,
            )
            .unwrap();
            let report = jet.advance(0.02, budget(dt), &gate).unwrap();
            assert_eq!(report.flow.outcome, SlenderJetOutcome::Completed);
            assert_eq!(
                report.flow.integration_method,
                SlenderJetIntegrationMethod::ExplicitSspRk3
            );
            assert!(report.flow.temporal_error.is_none());
            let actual = [
                1.0 - jet.flow().radius_squared_m2()[0],
                jet.flow().velocities_m_s().next().unwrap(),
                jet.cell_internal_energy_from_reference_j()[0],
                jet.cell_internal_energy_from_reference_j()[1],
                report.flow.viscous_dissipation_j,
            ];
            errors.push(std::array::from_fn::<_, 5, _>(|i| {
                (actual[i] - expected[i]).abs()
            }));
            near(report.energy_balance_residual_j, 0.0, 2e-13);
            near(report.mass_balance_residual_kg, 0.0, 2e-13);
            near(
                report.internal_energy_change_j,
                report.flow.viscous_dissipation_j,
                2e-13,
            );
            assert_eq!(report.flow.upwind_dissipation_j, 0.0);
            let cell_heat = jet.flow().cell_viscous_heat().unwrap();
            near(
                cell_heat.cell_heat_j.iter().sum(),
                report.internal_energy_change_j,
                2e-13,
            );
            assert!(cell_heat.cell_heat_j.iter().all(|heat| *heat > 0.0));
            let t: Vec<_> = jet.temperatures_k().collect();
            assert!(t[0] > 301.0 && t[1] > 330.0);
            let sample = jet.flow().carreau_sample_at(0).unwrap().unwrap();
            assert_eq!(sample.temperature().temperature_k(), t[0]);
            assert_eq!(sample.law(), &material().rheology);
            near(
                jet.flow().parameters().kinematic_viscosity_m2_s,
                jet.flow().prepared_carreau_cells().unwrap()[0].zero_shear_viscosity_pa_s() / 2.0,
                0.0,
            );
        }
        for pair in errors.windows(2) {
            for i in 0..5 {
                assert!(pair[1][i] < 0.20 * pair[0][i], "errors={errors:?}");
            }
        }
    }
}

#[test]
fn g0_uniform_temperature_moves_with_changing_mass_and_source_datum_is_irrelevant() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let a = fixture.caloric(280.0);
    let b = fixture.caloric(360.0);
    let radii: Vec<_> = (0..16)
        .map(|i| 1.0 + 0.03 * fs_math::det::cos(core::f64::consts::TAU * i as f64 / 16.0))
        .collect();
    let make = || {
        PeriodicSlenderJet::new_carreau(material(), 4.0, &radii, &[0.125; 16], limits(), &gate)
            .unwrap()
    };
    let mut first = ThermalCarreauJet::new(make(), &a, None, &[317.0; 16], 100.0, &gate).unwrap();
    let mut second = ThermalCarreauJet::new(make(), &b, None, &[317.0; 16], 100.0, &gate).unwrap();
    let masses: Vec<_> = first.cell_masses_kg().collect();
    let report = first.advance(0.02, budget(0.0005), &gate).unwrap();
    second.advance(0.02, budget(0.0005), &gate).unwrap();
    assert!(
        first
            .cell_masses_kg()
            .zip(masses)
            .any(|(after, before)| after != before)
    );
    for temperature in first.temperatures_k() {
        assert_eq!(temperature, 317.0);
    }
    assert_eq!(
        first.flow().radius_squared_m2(),
        second.flow().radius_squared_m2()
    );
    assert_eq!(
        first.flow().radius_squared_flux_m3_s(),
        second.flow().radius_squared_flux_m3_s()
    );
    assert_eq!(
        first.cell_internal_energy_from_reference_j(),
        second.cell_internal_energy_from_reference_j()
    );
    assert_eq!(
        first.temperatures_k().collect::<Vec<_>>(),
        second.temperatures_k().collect::<Vec<_>>()
    );
    near(report.internal_energy_change_j, 0.0, 1e-24);
    near(report.mass_balance_residual_kg, 0.0, 2e-13);
}

#[test]
fn g1_actual_mass_flux_advects_temperature_even_when_native_dual_momentum_flux_vanishes() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(280.0);
    let mut jet = ThermalCarreauJet::new(
        alternating(&gate),
        &caloric,
        None,
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    let report = jet.advance(0.02, budget(0.00025), &gate).unwrap();
    let masses: Vec<_> = jet.cell_masses_kg().collect();
    assert!(masses[0] < 0.5 * core::f64::consts::PI && masses[1] > 0.5 * core::f64::consts::PI);
    let heat = jet.flow().cell_viscous_heat().unwrap();
    // Cold donor cells export some newly generated heat with their moving mass.
    assert!(jet.cell_internal_energy_from_reference_j()[0] < heat.cell_heat_j[0]);
    // A frozen mass heat deposit would predict a different temperature.
    let frozen = caloric
        .temperature_after_specific_heat(300.0, heat.cell_heat_j[0] / (0.5 * core::f64::consts::PI))
        .unwrap();
    assert!((jet.temperatures_k().next().unwrap() - frozen).abs() > 0.005);
    near(report.internal_energy_change_j, 0.10364679290553484, 1e-8);
}

#[test]
fn g0_zero_duration_and_source_identity_density_and_policy_admission() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let other = fixture.caloric(310.0);
    let conduction = fixture.conduction(&other);
    assert!(matches!(
        ThermalCarreauJet::new(
            alternating(&gate),
            &caloric,
            Some(&conduction),
            &temperatures(),
            100.0,
            &gate,
        ),
        Err(ThermalJetError::InvalidInput(_))
    ));
    let wrong_density = Fixture::new([1e-5; 3], 3.0).caloric(300.0);
    assert!(matches!(
        ThermalCarreauJet::new(
            alternating(&gate),
            &wrong_density,
            None,
            &temperatures(),
            100.0,
            &gate,
        ),
        Err(ThermalJetError::InvalidInput(_))
    ));
    let mut jet = ThermalCarreauJet::new(
        alternating(&gate),
        &caloric,
        None,
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    let before = jet.cell_internal_energy_from_reference_j().to_vec();
    let report = jet.advance(0.0, budget(0.001), &gate).unwrap();
    assert_eq!(report.flow.accepted_substeps, 0);
    assert_eq!(report.internal_energy_change_j, 0.0);
    assert_eq!(jet.cell_internal_energy_from_reference_j(), before);
    assert_eq!(
        jet.flow().cell_viscous_heat().unwrap().cell_heat_j,
        &[0.0; 16]
    );
}

fn snapshot(
    jet: &ThermalCarreauJet<'_>,
) -> (
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    Vec<f64>,
    f64,
    SlenderJetDiagnostics,
) {
    (
        jet.flow().radius_squared_m2().to_vec(),
        jet.flow().radius_squared_flux_m3_s().to_vec(),
        jet.cell_internal_energy_from_reference_j().to_vec(),
        jet.temperatures_k().collect(),
        jet.flow().time_s(),
        jet.flow().diagnostics(),
    )
}

#[test]
fn g4_late_work_and_caloric_source_refusals_keep_all_accepted_state_and_heat() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let mut jet = ThermalCarreauJet::new(
        alternating(&gate),
        &caloric,
        None,
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    jet.advance(0.001, budget(0.00025), &gate).unwrap();
    let before = snapshot(&jet);
    let heat = jet.flow().cell_viscous_heat().unwrap().cell_heat_j.to_vec();
    let mut limited = budget(0.00025);
    limited.flow.maximum_substeps = 2;
    assert!(matches!(
        jet.advance(0.02, limited, &gate),
        Err(ThermalJetError::Flow(NonlinearJetError::SubstepBudget))
    ));
    assert_eq!(snapshot(&jet), before);
    assert_eq!(jet.flow().cell_viscous_heat().unwrap().cell_heat_j, heat);
    // Source upper bound is 360 K. Heating a positive-strain hot cell beyond
    // it refuses at the original inverse; no endpoint or source is clipped.
    let mut hot = ThermalCarreauJet::new(
        alternating(&gate),
        &caloric,
        None,
        &[359.999; 16],
        100.0,
        &gate,
    )
    .unwrap();
    let before = snapshot(&hot);
    assert!(matches!(
        hot.advance(0.02, budget(0.00025), &gate),
        Err(ThermalJetError::Thermal { .. })
    ));
    assert_eq!(snapshot(&hot), before);
    assert_eq!(
        hot.flow().cell_viscous_heat().unwrap().cell_heat_j,
        &[0.0; 16]
    );
}

#[test]
fn g0_variable_conductivity_reuses_relative_temperature_without_source_datum_roundtrip() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5, 3e-5, 1.5e-5], 2.0);
    let a = fixture.caloric(280.0);
    let b = fixture.caloric(360.0);
    let ka = fixture.conduction(&a);
    let kb = fixture.conduction(&b);
    let mut first = ThermalCarreauJet::new(
        alternating(&gate),
        &a,
        Some(&ka),
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    let mut second = ThermalCarreauJet::new(
        alternating(&gate),
        &b,
        Some(&kb),
        &temperatures(),
        100.0,
        &gate,
    )
    .unwrap();
    first.advance(0.02, budget(0.0005), &gate).unwrap();
    second.advance(0.02, budget(0.0005), &gate).unwrap();
    assert_eq!(snapshot(&first), snapshot(&second));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let before = snapshot(&first);
    assert!(matches!(
        first.advance(0.001, budget(0.001), &cancelled),
        Err(ThermalJetError::Flow(NonlinearJetError::Cancelled))
    ));
    assert_eq!(snapshot(&first), before);
}

#[test]
fn g0_geometry_stop_publishes_only_complete_coupled_substeps_and_their_source() {
    let gate = CancelGate::new_clock_free();
    let fixture = Fixture::new([1e-5; 3], 2.0);
    let caloric = fixture.caloric(300.0);
    let velocities: Vec<_> = (0..16)
        .map(|i| if i % 2 == 0 { 0.125 } else { -0.125 })
        .collect();
    let limited = SlenderJetLimits {
        minimum_radius_m: 0.998,
        ..limits()
    };
    let flow =
        PeriodicSlenderJet::new_carreau(material(), 4.0, &[1.0; 16], &velocities, limited, &gate)
            .unwrap();
    let mut jet =
        ThermalCarreauJet::new(flow, &caloric, None, &temperatures(), 100.0, &gate).unwrap();
    let report = jet.advance(0.02, budget(0.0005), &gate).unwrap();
    assert!(matches!(
        report.flow.outcome,
        SlenderJetOutcome::DomainLimit(_)
    ));
    assert!(report.flow.accepted_substeps > 0 && report.flow.end_time_s < 0.02);
    for ((energy, mass), temperature) in jet
        .cell_internal_energy_from_reference_j()
        .iter()
        .zip(jet.cell_masses_kg())
        .zip(jet.temperatures_k())
    {
        assert_eq!(
            temperature,
            caloric
                .temperature_at_relative_internal_energy(
                    jet.reference_temperature_k(),
                    energy / mass,
                )
                .unwrap()
        );
    }
    near(
        report.internal_energy_change_j,
        report.flow.viscous_dissipation_j,
        2e-13,
    );
    assert_eq!(jet.flow().time_s(), report.flow.end_time_s);
}
