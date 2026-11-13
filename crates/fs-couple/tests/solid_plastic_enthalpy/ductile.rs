//! Original coupled heat/momentum/Chaboche owners gated at actual ductile
//! endpoints. All coefficients and source claims are synthetic test data.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;
use fs_material::ductile::{DuctileDomain, JohnsonCookInitiation, JohnsonCookParameters};
use fs_solid::ductile::{DuctileAdvanceError, DuctileHistory, DuctileLawField, DuctileSolidError};

fn jc_law(critical: f64) -> JohnsonCookInitiation {
    JohnsonCookInitiation::new(
        JohnsonCookParameters {
            d1: critical,
            d2: 0.0,
            d3: 0.0,
            d4: 0.02,
            d5: -0.5,
            reference_strain_rate_s_inv: 1.0,
            transition_temperature_k: 280.0,
            melting_temperature_k: 330.0,
        },
        DuctileDomain {
            triaxiality: [-10.0, 10.0],
            strain_rate_s_inv: [1e-8, 1e6],
            temperature_k: [280.0, 320.0],
        },
    )
    .unwrap()
}

// Independent scalar specialization of the authored law, with no call to
// ductile::trial, fracture_strain or native stress reconstruction.
fn utilization(delta_alpha: f64, duration: f64, temperature: f64, critical: f64) -> f64 {
    delta_alpha
        / (critical
            * (1.0 + 0.02 * (delta_alpha / duration).ln())
            * (1.0 - 0.5 * (temperature - 280.0) / 50.0))
}

#[derive(Clone, Default)]
struct CountSteps(usize);
impl DepositionObserver for CountSteps {
    type Error = DepositionError;
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.0 += 1;
        Ok(())
    }
}

#[test]
fn heated_native_endpoints_use_actual_temperature_and_interval_without_repartitioning_heat() {
    context(|cx| {
        let fixture = Fixture::new(100.0);
        let dt = 0.001;
        let thermal = initial(&[&fixture], &[], &[RHO * V * 100.0 * 2.0 / dt]);
        let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
        let law = jc_law(0.08);
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&law), 1, 1).unwrap();
        let force = [(3, 3.4e6 * V / L)];
        let native = solid
            .advance(
                &thermal,
                &(),
                dt,
                &ROLLERS,
                &force,
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let first = solid
            .advance_with_ductility(
                &history,
                &thermal,
                &(),
                dt,
                &ROLLERS,
                &force,
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(
            first.mechanics.solid().solution(),
            native.solid().solution()
        );
        assert_eq!(first.mechanics.thermal().cells(), native.thermal().cells());
        assert_eq!(
            first.mechanics.report().physical_heat_j,
            native.report().physical_heat_j
        );
        assert_eq!(
            first.mechanics.report().numerical_loss_j,
            native.report().numerical_loss_j
        );
        let alpha = first.mechanics.solid().solution().element_states[0].alpha;
        let temperature = first
            .mechanics
            .thermal()
            .temperature_k(CellId(100))
            .unwrap();
        assert!(alpha > 0.0 && temperature > 302.0);
        close(
            first.ductility.points()[0].utilization(),
            utilization(alpha, dt, temperature, 0.08),
            1e-13,
        );
        assert!(
            (first.ductility.points()[0].utilization() - utilization(alpha, dt, 300.0, 0.08)).abs()
                > 1e-3
        );
        close(
            first
                .mechanics
                .report()
                .thermal_audit
                .stored_energy_change_j,
            first.mechanics.report().physical_heat_j + first.mechanics.report().external_heat_j,
            1e-8,
        );
        assert_eq!(history.points()[0].utilization(), 0.0);
        assert_eq!(solid.solution().element_states[0].alpha, 0.0);
        assert_eq!(thermal.time_s(), 0.0);

        let second = first
            .mechanics
            .solid()
            .advance_with_ductility(
                &first.ductility,
                first.mechanics.thermal(),
                first.mechanics.observer(),
                4.0 * dt,
                &ROLLERS,
                &[(3, 3.6e6 * V / L)],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let next_alpha = second.mechanics.solid().solution().element_states[0].alpha;
        let next_t = second
            .mechanics
            .thermal()
            .temperature_k(CellId(100))
            .unwrap();
        assert!(next_alpha > alpha && next_t > temperature);
        close(
            second.ductility.points()[0].accumulated_plastic_strain(),
            next_alpha,
            0.0,
        );
        close(
            second.ductility.points()[0].utilization(),
            first.ductility.points()[0].utilization()
                + utilization(next_alpha - alpha, 3.0 * dt, next_t, 0.08),
            1e-13,
        );
    });
}

#[test]
fn initiation_and_stale_history_refuse_the_complete_coupled_proposal_and_allow_refinement() {
    context(|cx| {
        let fixture = Fixture::new(10.0);
        let thermal = initial(&[&fixture], &[], &[0.0]);
        let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
        let before = solid.solution().clone();
        let cells = thermal.cells().clone();
        let observer = CountSteps::default();
        let law = jc_law(0.005);
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&law), 1, 1).unwrap();
        let failed = solid.advance_with_ductility(
            &history,
            &thermal,
            &observer,
            0.001,
            &ROLLERS,
            &[(3, 3.4e6 * V / L)],
            IncrementSettings::default(),
            thermal_config(),
            config(),
            cx,
        );
        assert!(matches!(
            failed,
            Err(DuctileAdvanceError::Ductility(
                DuctileSolidError::Initiated { point: 0, .. }
            ))
        ));
        assert_eq!(solid.solution(), &before);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(observer.0, 0);
        assert_eq!(history.points()[0].utilization(), 0.0);
        let refined = solid
            .advance_with_ductility(
                &history,
                &thermal,
                &observer,
                0.0005,
                &ROLLERS,
                &[(3, 3.1e6 * V / L)],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert!(refined.mechanics.report().physical_heat_j > 0.0);
        assert!(refined.mechanics.observer().0 > 0);
        assert!(refined.ductility.points()[0].utilization() > 0.0);
        assert!(refined.ductility.points()[0].utilization() < 1.0);
        let reset = refined.mechanics.solid().advance_with_ductility(
            &history,
            refined.mechanics.thermal(),
            refined.mechanics.observer(),
            0.001,
            &ROLLERS,
            &[(3, 3.1e6 * V / L)],
            IncrementSettings::default(),
            thermal_config(),
            config(),
            cx,
        );
        assert!(matches!(
            reset,
            Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))
        ));
        assert_eq!(observer.0, 0);
        assert_eq!(solid.solution(), &before);
        assert_eq!(thermal.cells(), &cells);
    });
}

#[test]
fn native_momentum_ductility_and_material_heat_share_one_physical_interval() {
    context(|cx| {
        let fixture = Fixture::new(5.0);
        let thermal = initial(&[&fixture], &[], &[0.0]);
        let mut velocity = vec![[0.0; 3]; 4];
        velocity[1][0] = 30.0;
        let moving = ThermoplasticDynamics::try_new(
            bind(&thermal, &[&fixture], &NODES, &TETS, cx),
            velocity.clone(),
            cx,
        )
        .unwrap();
        let dt = 2e-6;
        let law = jc_law(0.08);
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&law), 1, 1).unwrap();
        let native = moving
            .advance(
                &thermal,
                &(),
                dt,
                &ROLLERS,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let gated = moving
            .advance_with_ductility(
                &history,
                &thermal,
                &(),
                dt,
                &ROLLERS,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(
            gated.mechanics.dynamics().solid().solution(),
            native.dynamics().solid().solution()
        );
        assert_eq!(
            gated.mechanics.dynamics().velocity_mps(),
            native.dynamics().velocity_mps()
        );
        assert_eq!(gated.mechanics.thermal().cells(), native.thermal().cells());
        let report = gated.mechanics.report();
        assert!(report.material.physical_heat_j > 0.0);
        assert!(report.inertia.kinetic_change_j < 0.0);
        assert!(report.inertia.numerical_loss_j > 0.0);
        assert_eq!(
            report.material.physical_heat_j,
            native.report().material.physical_heat_j
        );
        assert_eq!(
            report.inertia.numerical_loss_j,
            native.report().inertia.numerical_loss_j
        );
        close(
            report.material.thermal_audit.stored_energy_change_j,
            report.material.physical_heat_j,
            1e-8,
        );
        let alpha = gated.mechanics.dynamics().solid().solution().element_states[0].alpha;
        let temperature = gated
            .mechanics
            .thermal()
            .temperature_k(CellId(100))
            .unwrap();
        close(
            gated.ductility.points()[0].utilization(),
            utilization(alpha, dt, temperature, 0.08),
            1e-13,
        );
        let fragile = jc_law(1e-5);
        let fragile_history =
            DuctileHistory::virgin(DuctileLawField::Uniform(&fragile), 1, 1).unwrap();
        assert!(matches!(
            moving.advance_with_ductility(
                &fragile_history,
                &thermal,
                &(),
                dt,
                &ROLLERS,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx
            ),
            Err(DuctileAdvanceError::Ductility(
                DuctileSolidError::Initiated { .. }
            ))
        ));
        assert_eq!(moving.velocity_mps(), &velocity);
        assert_eq!(moving.solid().solution().element_states[0].alpha, 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        let gate = CancelGate::new();
        gate.request();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 1,
                    kernel_id: 2,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(matches!(
                moving.advance_with_ductility(
                    &history,
                    &thermal,
                    &(),
                    dt,
                    &ROLLERS,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    &cancelled
                ),
                Err(DuctileAdvanceError::Ductility(DuctileSolidError::Cancelled))
            ));
        });
        assert_eq!(moving.velocity_mps(), &velocity);
        assert_eq!(history.points()[0].utilization(), 0.0);
    });
}

#[test]
fn complete_native_chaboche_backstress_and_reverse_plastic_accumulation_survive_the_gate() {
    context(|cx| {
        let fixture = chaboche::fixture(None, false);
        let thermal = initial(&[&fixture], &[], &[0.0]);
        let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx)
            .with_voce_hardening(&[Some(MaterialPropertySelection::SingleClaimOnly)], cx)
            .unwrap()
            .with_chaboche_branches(
                &[(
                    MaterialPropertySelection::SingleClaimOnly,
                    chaboche::schema(),
                )],
                cx,
            )
            .unwrap();
        let law = jc_law(0.08);
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&law), 1, 1).unwrap();
        let shear = |gamma: f64| -> Vec<_> {
            (0..12)
                .map(|dof| {
                    (
                        dof,
                        if dof % 3 == 0 {
                            gamma * NODES[dof / 3][1]
                        } else {
                            0.0
                        },
                    )
                })
                .collect()
        };
        let first = solid
            .advance_with_ductility(
                &history,
                &thermal,
                &(),
                0.001,
                &shear(0.02),
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let state = &first.mechanics.solid().solution().element_states[0];
        let alpha = state.plastic.alpha;
        assert!(alpha > 0.0);
        assert!(
            state
                .backstress_pa
                .iter()
                .all(|branch| branch.iter().any(|value| *value != 0.0))
        );
        let native = first
            .mechanics
            .solid()
            .advance(
                first.mechanics.thermal(),
                &(),
                0.003,
                &shear(-0.02),
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        let reverse = first
            .mechanics
            .solid()
            .advance_with_ductility(
                &first.ductility,
                first.mechanics.thermal(),
                &(),
                0.003,
                &shear(-0.02),
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(
            reverse.mechanics.solid().solution(),
            native.solid().solution()
        );
        assert_eq!(
            reverse.mechanics.thermal().cells(),
            native.thermal().cells()
        );
        assert_eq!(
            reverse.mechanics.report().physical_heat_j,
            native.report().physical_heat_j
        );
        let next = &reverse.mechanics.solid().solution().element_states[0];
        assert!(next.plastic.alpha > alpha);
        assert_ne!(next.backstress_pa, state.backstress_pa);
        let temperature = reverse
            .mechanics
            .thermal()
            .temperature_k(CellId(100))
            .unwrap();
        close(
            reverse.ductility.points()[0].accumulated_plastic_strain(),
            next.plastic.alpha,
            0.0,
        );
        close(
            reverse.ductility.points()[0].utilization(),
            first.ductility.points()[0].utilization()
                + utilization(next.plastic.alpha - alpha, 0.002, temperature, 0.08),
            1e-13,
        );
        assert_eq!(solid.solution().element_states[0].plastic.alpha, 0.0);
        assert_eq!(
            solid.solution().element_states[0].backstress_pa,
            [[0.0; 6]; 2]
        );
    });
}

use fs_matdb::ClaimSelection;
use fs_material::ductile::source::{
    DuctileSourceRequest, FRACTURE_PLASTIC_STRAIN_PROPERTY, PLASTIC_STRAIN_RATE_AXIS,
    STRESS_TRIAXIALITY_AXIS, resolve_ductile_locus,
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn typed_temperature(value: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_spec(), value)
        .unwrap()
}
fn source_fixture(chemistry: &str) -> Fixture {
    let original = Fixture::new(100.0);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.card.claims().claims_ordered() {
        let mut claim = claim.clone();
        if claim.key.name() == YIELD_STRESS_PROPERTY {
            claim.validity = ValidityDomain::unconstrained().with_quantity(
                "T",
                temperature_spec(),
                280.0,
                320.0,
            );
        }
        claims.insert_claim(claim).unwrap();
    }
    let mut fracture = scalar(FRACTURE_PLASTIC_STRAIN_PROPERTY, Dims::NONE, 0.02);
    fracture.value = PropertyValue::Curve {
        abscissa: STRESS_TRIAXIALITY_AXIS.into(),
        abscissa_dims: Dims::NONE,
        knots: vec![(-2.0, 0.02), (2.0, 0.02)],
        dims: Dims::NONE,
    };
    fracture.interpolation = InterpolationPolicy::LinearInside;
    fracture.validity = ValidityDomain::unconstrained()
        .with_quantity("T", temperature_spec(), 280.0, 320.0)
        .with_quantity(
            STRESS_TRIAXIALITY_AXIS,
            QuantitySpec::dimensional(Dims::NONE),
            -2.0,
            2.0,
        )
        .with_quantity(
            PLASTIC_STRAIN_RATE_AXIS,
            QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
            1e-8,
            1e6,
        );
    claims.insert_claim(fracture).unwrap();
    let mut identity = original.card.id().clone();
    identity.chemistry = chemistry.into();
    let card = MaterialCard::assemble(identity, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &typed_temperature(280.0),
            upper: &typed_temperature(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -2000.0,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: false,
    }
}

#[test]
fn original_source_reference_and_fixed_temperature_calibration_cannot_be_bypassed() {
    context(|cx| {
        let fixture = source_fixture("synthetic ductile thermal source");
        let foreign = source_fixture("numerically identical foreign source");
        let resolve = |fixture: &Fixture, point: QueryPoint| {
            let point = point
                .with_quantity(
                    STRESS_TRIAXIALITY_AXIS,
                    QuantitySpec::dimensional(Dims::NONE),
                    -2.0,
                )
                .unwrap()
                .with_quantity(
                    PLASTIC_STRAIN_RATE_AXIS,
                    QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0])),
                    1.0,
                )
                .unwrap();
            resolve_ductile_locus(
                &fixture.card,
                &point,
                &DuctileSourceRequest {
                    temperature_axis: "T".into(),
                    maximum_points: 2,
                },
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .unwrap()
        };
        let locus = resolve(&fixture, typed_temperature(300.0));
        let foreign_locus = resolve(&foreign, typed_temperature(300.0));
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&locus), 1, 1).unwrap();
        let foreign_history =
            DuctileHistory::virgin(DuctileLawField::Uniform(&foreign_locus), 1, 1).unwrap();
        let fixed: Vec<_> = (0..12).map(|dof| (dof, 0.0)).collect();
        for heating_k in [0.0, 5.0] {
            let thermal = initial(&[&fixture], &[], &[RHO * V * 100.0 * heating_k / 0.001]);
            let mut query = fixture.query();
            query.context = typed_temperature(300.0);
            let solid = ThermoplasticSolid::try_new(
                &thermal,
                &NODES,
                &TETS,
                vec![ThermoplasticCell {
                    cell: CellId(100),
                    query,
                }],
                TetAssemblyBudget::standard(),
                cx,
            )
            .unwrap();
            let cells = thermal.cells().clone();
            let before = solid.solution().clone();
            let observer = CountSteps::default();
            assert!(matches!(
                solid.advance_with_ductility(
                    &foreign_history,
                    &thermal,
                    &observer,
                    0.001,
                    &fixed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx
                ),
                Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))
            ));
            let proposal = solid.advance_with_ductility(
                &history,
                &thermal,
                &observer,
                0.001,
                &fixed,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            );
            if heating_k == 0.0 {
                let accepted = proposal.unwrap();
                assert_eq!(accepted.ductility.points()[0].utilization(), 0.0);
                assert!(accepted.mechanics.observer().0 > 0);
            } else {
                // The table's validity box is not a temperature calibration:
                // the resolved locus remains bound to its original T=300 K.
                assert!(matches!(
                    proposal,
                    Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))
                ));
                let warmed = solid
                    .advance(
                        &thermal,
                        &observer,
                        0.001,
                        &fixed,
                        &[],
                        IncrementSettings::default(),
                        thermal_config(),
                        config(),
                        cx,
                    )
                    .unwrap();
                let rebound = resolve(&fixture, warmed.solid().materials()[0].material().query());
                assert!(
                    rebound
                        .source()
                        .matches_material(warmed.solid().materials()[0].material())
                );
                let rebound_history =
                    DuctileHistory::virgin(DuctileLawField::Uniform(&rebound), 1, 1).unwrap();
                // Matching only the latest material snapshot must not turn a
                // fresh table at a changed T into original source authority.
                assert!(matches!(
                    warmed.solid().advance_with_ductility(
                        &rebound_history,
                        warmed.thermal(),
                        warmed.observer(),
                        0.002,
                        &fixed,
                        &[],
                        IncrementSettings::default(),
                        thermal_config(),
                        config(),
                        cx
                    ),
                    Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(
                        "ductile calibration must match the original thermal/mechanical source reference"
                    )))
                ));
            }
            assert_eq!(solid.solution(), &before);
            assert_eq!(thermal.cells(), &cells);
            assert_eq!(observer.0, 0);
            assert_eq!(history.points()[0].utilization(), 0.0);
        }
    });
}

#[test]
fn original_johnson_cook_fit_tracks_heated_native_endpoint_without_repartitioning_heat() {
    use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter};
    use fs_material::ductile::source::{
        JOHNSON_COOK_INITIATION_LAW_ID, resolve_johnson_cook_initiation,
    };
    context(|cx| {
        let mut fixture = source_fixture("synthetic temperature-dependent JC source");
        let rate_dims = Dims([0, 0, -1, 0, 0, 0]);
        let rate = QuantitySpec::dimensional(rate_dims);
        let eta = QuantitySpec::dimensional(Dims::NONE);
        let artifact = fs_blake3::hash_bytes(b"synthetic complete thermal JC fit");
        let model = ConstitutiveModelCard {
            law: LawId(JOHNSON_COOK_INITIATION_LAW_ID.into()),
            law_version: 1,
            parameters: [
                ("d1", 0.08, Dims::NONE),
                ("d2", 0.0, Dims::NONE),
                ("d3", 0.0, Dims::NONE),
                ("d4", 0.02, Dims::NONE),
                ("d5", -0.5, Dims::NONE),
                ("reference_strain_rate_s_inv", 1.0, rate_dims),
                ("transition_temperature_k", 280.0, Temperature::DIMS),
                ("melting_temperature_k", 330.0, Temperature::DIMS),
            ]
            .into_iter()
            .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
            .collect(),
            state_schema_version: 1,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained()
                .with_quantity(STRESS_TRIAXIALITY_AXIS, eta, -10.0, 10.0)
                .with_quantity(PLASTIC_STRAIN_RATE_AXIS, rate, 1e-8, 1e6)
                .with_quantity("T", temperature_spec(), 280.0, 320.0),
            sources: vec![artifact],
            provenance: Provenance {
                source: "synthetic coupled JC regression".into(),
                license: "test-only".into(),
                artifact: Some(artifact),
            },
        };
        let pin = model.content_hash();
        fixture.card = MaterialCard::assemble(
            fixture.card.id().clone(),
            fixture.card.claims().clone(),
            vec![model.clone()],
        )
        .unwrap();
        fixture.caloric = ResolvedHeatCapacityCurve::try_new(
            &fixture.card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(CP, CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "T",
                lower: &typed_temperature(280.0),
                upper: &typed_temperature(320.0),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -2000.0,
            1e-5,
        )
        .unwrap();
        let point = typed_temperature(300.0)
            .with_quantity(STRESS_TRIAXIALITY_AXIS, eta, 1.0 / 3.0)
            .unwrap()
            .with_quantity(PLASTIC_STRAIN_RATE_AXIS, rate, 1.0)
            .unwrap();
        let law = resolve_johnson_cook_initiation(&fixture.card, pin, &point, "T").unwrap();
        // A lineage starts at revision 0, so the next revision of the same
        // material comes from supersede.
        let foreign_card =
            MaterialCard::supersede(&fixture.card, fixture.card.claims().clone(), vec![model])
                .unwrap();
        let foreign = resolve_johnson_cook_initiation(&foreign_card, pin, &point, "T").unwrap();
        let pressure =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        let pressure_point = point
            .clone()
            .with_quantity("ambient-pressure", pressure, 101325.0)
            .unwrap();
        let wrong_context =
            resolve_johnson_cook_initiation(&fixture.card, pin, &pressure_point, "T").unwrap();

        let dt = 0.001;
        let thermal = initial(&[&fixture], &[], &[RHO * V * 100.0 * 2.0 / dt]);
        let mut query = fixture.query();
        query.context = typed_temperature(300.0);
        let solid = ThermoplasticSolid::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![ThermoplasticCell {
                cell: CellId(100),
                query,
            }],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let before = solid.solution().clone();
        let cells = thermal.cells().clone();
        let observer = CountSteps::default();
        let force = [(3, 3.4e6 * V / L)];
        let solver = IncrementSettings::default();
        let heat = thermal_config();
        let coupled = config();
        let native = solid
            .advance(
                &thermal, &observer, dt, &ROLLERS, &force, solver, heat, coupled, cx,
            )
            .unwrap();
        let advance = |history| {
            solid.advance_with_ductility(
                history, &thermal, &observer, dt, &ROLLERS, &force, solver, heat, coupled, cx,
            )
        };
        let history = DuctileHistory::virgin(DuctileLawField::Uniform(&law), 1, 1).unwrap();
        let accepted = advance(&history).unwrap();
        let actual = &accepted.mechanics;
        assert_eq!(actual.solid().solution(), native.solid().solution());
        assert_eq!(actual.thermal().cells(), native.thermal().cells());
        let (report, baseline) = (actual.report(), native.report());
        assert_eq!(report.physical_heat_j, baseline.physical_heat_j);
        assert_eq!(report.numerical_loss_j, baseline.numerical_loss_j);
        assert_eq!(actual.observer().0, native.observer().0);
        let alpha = actual.solid().solution().element_states[0].alpha;
        let temperature = actual.thermal().temperature_k(CellId(100)).unwrap();
        assert!(alpha > 0.0 && temperature > 302.0 && temperature < 320.0);
        let endpoint = actual.solid().materials()[0].material();
        assert!(law.source().matches_material(endpoint));
        assert_eq!(endpoint.state_coordinate("T"), Some(temperature));
        let actual_utilization = accepted.ductility.points()[0].utilization();
        close(
            actual_utilization,
            utilization(alpha, dt, temperature, 0.08),
            1e-13,
        );
        assert!((actual_utilization - utilization(alpha, dt, 300.0, 0.08)).abs() > 1e-3);
        // Neither an identical fit on another card nor a different fixed
        // environment can authorize the same mechanical/thermal proposal.
        // Built before the loop: `advance` takes one inferred history borrow,
        // which must outlive every call, so a history dropped per iteration fails.
        let mismatched_histories = [&foreign, &wrong_context].map(|mismatched| {
            DuctileHistory::virgin(DuctileLawField::Uniform(mismatched), 1, 1).unwrap()
        });
        for bad in &mismatched_histories {
            assert!(matches!(
                advance(bad),
                Err(DuctileAdvanceError::Ductility(DuctileSolidError::Input(_)))
            ));
            assert_eq!(bad.points()[0].utilization(), 0.0);
        }
        assert_eq!(solid.solution(), &before);
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(observer.0, 0);
        assert_eq!(history.points()[0].utilization(), 0.0);
    });
}
