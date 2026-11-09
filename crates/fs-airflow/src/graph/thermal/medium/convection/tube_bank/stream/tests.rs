use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_convection::tube_bank::TubeBankLayout;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, Provenance, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_material::thermal::transport::isobaric::IsobaricHeatCapacityQuery;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 229,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    with_gate(&CancelGate::new_clock_free(), f)
}
fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), p)
        .unwrap()
}
fn card(varying: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, constant, changing) in [
        (DENSITY_PROPERTY, Density::DIMS, 1000.0, [1025.0, 925.0]),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            1000.0,
            [500.0, 2500.0],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            Dims([1, 1, -3, -1, 0, 0]),
            0.5,
            [0.45, 0.65],
        ),
        (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            0.001,
            [0.00225, 0.00125],
        ),
    ] {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        let values = if varying { changing } else { [constant; 2] };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(250.0, values[0]), (450.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::dimensional(Temperature::DIMS),
                        250.0,
                        450.0,
                    )
                    .with_quantity(
                        "p",
                        QuantitySpec::dimensional(Pressure::DIMS),
                        100_000.0,
                        200_000.0,
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured finite bank fluid, not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured bank coolant".into(),
            phase: "liquid".into(),
            process: "fixed isobar".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn caloric(card: &MaterialCard, pressure: f64, datum: f64) -> IsobaricHeatCapacityLaw {
    IsobaricHeatCapacityLaw::try_new(
        card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &PropertyKey::new(
                SPECIFIC_HEAT_CAPACITY_PROPERTY,
                SPECIFIC_HEAT_CAPACITY_DIMS,
            ),
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &point(250.0, pressure),
            upper: &point(450.0, pressure),
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        },
        datum,
    )
    .unwrap()
}
fn liquid(card: &MaterialCard) -> LiquidTemperatureLaw<'_> {
    LiquidTemperatureLaw::try_new(
        card,
        &point(300.0, 100_000.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn surface(columns: u32) -> TubeBankSurface {
    let geometry = TubeBankGeometry::try_new(
        TubeBankLayout::Inline,
        16,
        Length::new(0.01),
        Length::new(0.02),
        Length::new(0.02),
    )
    .unwrap();
    TubeBankSurface::new(geometry, columns, Length::new(1.0)).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}

#[test]
fn g1_finite_bank_heating_and_cooling_match_independent_ntu_and_robin_heat() {
    with_cx(|cx| {
        let card = card(false);
        let liquid = liquid(&card);
        let caloric = caloric(&card, 100_000.0, 300.0);
        let area = surface(4);
        near(area.frontal_area_m2(), 0.08, 1e-15);
        near(area.surface_area_m2(), 0.64 * core::f64::consts::PI, 1e-15);
        for (inlet, wall) in [(300.0, 400.0), (400.0, 300.0), (350.0, 350.0)] {
            let stream = LiquidTubeBankStream::new(&liquid, &caloric, area, inlet, 80.0).unwrap();
            let state = stream
                .solve(cx, wall, TubeBankStreamConfig::default())
                .unwrap();
            // Independent constant-Cp source formula, not the caloric kernel.
            let h = 0.27 * 20_000_f64.powf(0.63) * 2_f64.powf(0.36) * 50.0;
            let ntu = h * area.surface_area_m2() / 80_000.0;
            let outlet = inlet + (wall - inlet) * -(-ntu).exp_m1();
            near(state.exchange.temperature_k, outlet, 2e-11);
            near(state.approach_speed_m_s, 1.0, 1e-15);
            near(state.heat_into_fluid_w, 80_000.0 * (outlet - inlet), 2e-8);
            assert_eq!(state.heat_into_wall_w(), -state.heat_into_fluid_w);
            near(state.endpoint_energy_residual_w, 0.0, 2e-8);
            near(state.robin_energy_residual_w, 0.0, 2e-8);
            let robin = state.robin_boundary().unwrap();
            near(robin.coefficient().value.value(), h, 1e-8);
            near(
                state.conductance_w_k * (wall - state.exchange.mean_temperature_k),
                state.heat_into_fluid_w,
                2e-8,
            );
            assert!(state.exchange.mean_temperature_k >= inlet.min(outlet));
            assert!(state.exchange.mean_temperature_k <= inlet.max(outlet));
        }
    });
}

#[test]
fn g1_variable_source_calorics_and_mean_transport_recover_a_manufactured_outlet() {
    with_cx(|cx| {
        let card = card(true);
        let liquid = liquid(&card);
        let caloric = caloric(&card, 100_000.0, 250.0);
        let area = surface(4);
        // Cp=1000+10(T-300); choose the target outlet first. Analytically
        // integrate Cp/(400-T), then invert the source's m_dot^(.63-1) scaling.
        let exposure = 2000.0 * (100_f64 / 80.0).ln() - 200.0;
        let pr: f64 = 0.00195 * 1100.0 / 0.51;
        let re_per_mass: f64 = 0.01 * 2.0 / (0.08 * 0.00195);
        let h_per_mass =
            0.27 * re_per_mass.powf(0.63) * pr.powf(0.36) * (pr / 5.0).powf(0.25) * 51.0;
        let mass = (h_per_mass * area.surface_area_m2() / exposure).powf(1.0 / 0.37);
        assert!((1000.0..200_000.0).contains(&(mass * re_per_mass)));
        let stream = LiquidTubeBankStream::new(&liquid, &caloric, area, 300.0, mass).unwrap();
        let state = stream
            .solve(cx, 400.0, TubeBankStreamConfig::default())
            .unwrap();
        near(state.exchange.temperature_k, 320.0, 2e-9);
        near(state.property_temperature_k, 310.0, 2e-9);
        near(state.bulk_source.specific_heat_j_kg_k(), 1100.0, 3e-8);
        near(state.inlet_speed_m_s, mass / 80.0, 1e-14);
        near(state.approach_speed_m_s, mass / (0.08 * 995.0), 3e-12);
        near(
            state.convection.groups()["Re_max"],
            mass * re_per_mass,
            3e-7,
        );
        near(state.surface_source.prandtl(), 5.0, 1e-14);
        near(state.exchange.specific_enthalpy_change_j_kg, 22_000.0, 3e-6);
        near(state.heat_into_fluid_w, mass * 22_000.0, 2e-4);
        // Independent 60-digit evaluation of 400 - 22000/exposure.
        near(state.exchange.mean_temperature_k, 310.67335737352013, 2e-9);
    });
}

#[test]
fn g3_parallel_pitch_cells_scale_heat_without_changing_the_stream_state() {
    with_cx(|cx| {
        let card = card(true);
        let liquid = liquid(&card);
        let first_law = caloric(&card, 100_000.0, 250.0);
        let second_law = caloric(&card, 100_000.0, 400.0);
        let a = LiquidTubeBankStream::new(&liquid, &first_law, surface(4), 300.0, 80.0)
            .unwrap()
            .solve(cx, 400.0, TubeBankStreamConfig::default())
            .unwrap();
        let b = LiquidTubeBankStream::new(&liquid, &second_law, surface(8), 300.0, 160.0)
            .unwrap()
            .solve(cx, 400.0, TubeBankStreamConfig::default())
            .unwrap();
        assert_eq!(a.exchange, b.exchange);
        assert_eq!(a.bulk_source, b.bulk_source);
        assert_eq!(a.convection, b.convection);
        assert_eq!(2.0 * a.heat_into_fluid_w, b.heat_into_fluid_w);
    });
}

#[test]
fn g0_source_flow_and_work_refusals_leave_original_owners_reusable() {
    with_cx(|cx| {
        let card = card(true);
        let liquid = liquid(&card);
        let other_pressure = caloric(&card, 200_000.0, 300.0);
        let caloric = caloric(&card, 100_000.0, 300.0);
        assert!(
            LiquidTubeBankStream::new(&liquid, &other_pressure, surface(4), 300.0, 80.0).is_err()
        );
        assert!(LiquidTubeBankStream::new(&liquid, &caloric, surface(4), 300.0, 0.0).is_err());
        assert!(TubeBankSurface::new(surface(4).bank(), 0, Length::new(1.0)).is_err());
        let stream = LiquidTubeBankStream::new(&liquid, &caloric, surface(4), 300.0, 80.0).unwrap();
        assert!(
            stream
                .solve(cx, 500.0, TubeBankStreamConfig::default())
                .is_err()
        );
        assert!(matches!(
            stream.solve(
                cx,
                400.0,
                TubeBankStreamConfig {
                    maximum_iterations: 1,
                    ..TubeBankStreamConfig::default()
                }
            ),
            Err(TubeBankStreamError::NonConvergence { .. })
        ));
        let low_flow =
            LiquidTubeBankStream::new(&liquid, &caloric, surface(4), 300.0, 0.001).unwrap();
        assert!(matches!(
            low_flow.solve(cx, 400.0, TubeBankStreamConfig::default()),
            Err(TubeBankStreamError::Convection(_))
        ));
        let again = stream
            .solve(cx, 400.0, TubeBankStreamConfig::default())
            .unwrap();
        assert!(again.heat_into_fluid_w > 0.0);
        assert_eq!(liquid.reference_point(), &point(300.0, 100_000.0));
    });
}

#[test]
fn g4_cancelled_stream_never_returns_a_partial_exchanger() {
    let card = card(false);
    let liquid = liquid(&card);
    let caloric = caloric(&card, 100_000.0, 300.0);
    let stream = LiquidTubeBankStream::new(&liquid, &caloric, surface(4), 300.0, 80.0).unwrap();
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_gate(&gate, |cx| {
        assert!(matches!(
            stream.solve(cx, 400.0, TubeBankStreamConfig::default()),
            Err(TubeBankStreamError::Cancelled)
        ))
    });
    with_cx(|cx| {
        assert!(
            stream
                .solve(cx, 400.0, TubeBankStreamConfig::default())
                .is_ok()
        )
    });
}
