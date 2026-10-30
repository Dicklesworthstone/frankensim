//! G1/G3/G4: chemical DATA reaches finite atmospheric heat and actual VOF motion.
//! Manufactured isomers and liquid laws test coupling, not an empirical flame.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::cavity_phs::piston::reactive::{GasReaction, GasReactionSolve};
use fs_couple::liquid_gas::{ReactiveLiquidGasError, ReactiveLiquidGasSolve};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceCellHeat, FreeSurfaceEnthalpyControls, MaterialFreeSurfaceEnthalpy3,
    conduction::{FreeSurfaceConductionControls, FreeSurfaceGasFace},
};
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceInitialCell3, plan_free_surface_from_material,
};
use fs_matdb::*;
use fs_material::liquid::surface::{LiquidCapillarity, resolve_liquid_surface_state};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{GasEnergySolve, IsochoricGasParcel, parcel::ClosedGasParcel},
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_thermal::enthalpy::CellGeometry;
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::*;

#[allow(dead_code)]
#[path = "material_exchange/gas.rs"]
mod gas_data;

#[path = "reactive_liquid_gas/surface.rs"]
mod surface;

const EXPOSED: [usize; 3] = [3, 3, 3];
const VOLUME: f64 = 1e-8;
const DT: f64 = 0.001;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn reaction(rate: f64, release: f64, conductivity_factor: f64) -> GasReaction {
    let (a, ta) = gas_data::component("A", 0.028, 0.0);
    let (b, mut tb) = gas_data::component("B", 0.028, -release);
    let mut claims = ClaimSet::new();
    for (_, original) in tb.material.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == "thermal_conductivity" {
            claim.value = PropertyValue::Scalar {
                value: 0.03 * conductivity_factor,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            };
        }
        claims.insert_claim(claim).unwrap();
    }
    tb.material = MaterialCard::assemble(tb.material.id().clone(), claims, vec![]).unwrap();
    let species = vec![a.species().clone(), b.species().clone()];
    let id = ReactionId::new("manufactured-isomerization").unwrap();
    let matrix =
        StoichiometricMatrix::new(species.clone(), vec![id.clone()], vec![vec![-1], vec![1]])
            .unwrap();
    let elements = ElementalMatrix::new(
        vec![ElementId::new("X").unwrap()],
        species.clone(),
        vec![vec![1, 1]],
    )
    .unwrap();
    let charge = ChargeVector::new(species, vec![0, 0]).unwrap();
    let certificate = verify_conservation(&elements, &matrix, &charge).unwrap();
    let equilibrium =
        IdealGasReactionEquilibriumV1::new(matrix, certificate, id, vec![a.clone(), b.clone()])
            .unwrap();
    let kinetics = DeclaredStoichiometricMassActionV1::new(
        ThermodynamicReverseRateClosureV1::new(
            equilibrium,
            ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(rate)).unwrap(),
        ),
        MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients,
    );
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![a, b],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![ta, tb],
        None,
    )
    .unwrap();
    GasReaction::new(catalog, kinetics).unwrap()
}

fn gas(reaction: &GasReaction) -> IsochoricGasParcel {
    let model = reaction
        .catalog()
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            400.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
        .unwrap();
    IsochoricGasParcel::new(model, point, VOLUME).unwrap()
}

fn controls() -> ReactiveLiquidGasSolve {
    ReactiveLiquidGasSolve {
        reaction: GasReactionSolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 300.0,
                maximum_temperature_k: 690.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            maximum_consumed_fraction: 0.2,
            species_tolerance_mol: 1e-15,
            entropy_tolerance_j_k: 1e-12,
        },
        thermal: FreeSurfaceConductionControls {
            maximum_faces: 8,
            maximum_substeps: 1000,
            maximum_dt_s: DT / 4.0,
            energy_tolerance_j: 1e-9,
        },
        maximum_gas_temperature_increment_k: 200.0,
        total_energy_tolerance_j: 2e-9,
        with_viscous_heating: false,
    }
}

fn face() -> FreeSurfaceGasFace {
    FreeSurfaceGasFace {
        liquid: EXPOSED,
        gas: [4, 3, 3],
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    }
}

fn liquid_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), [1000.0; 2]),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            [2000.0; 2],
        ),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, [0.5; 2]),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), [0.1, 0.02]),
    ] {
        let key = if name == "thermal_conductivity" {
            PropertyKey::with_quantity(name, quantity(QuantityKind::ThermalConductivity))
        } else {
            PropertyKey::new(name, dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (500.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 500.0)
                    .with("pressure", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured reacting-liquid fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured liquid".into(),
            phase: "liquid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn with_liquid<R>(
    moving: bool,
    run: impl FnOnce(&MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = liquid_card();
    let point = QueryPoint::new()
        .with("T", 400.0)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap();
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Liquid,
        [300.0, 600.0]
            .map(|temperature_k| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: (temperature_k - 300.0) * 2000.0,
                temperature_k,
                liquid_mass_fraction: 1.0,
                bulk_density_kg_m3: 1000.0,
            })
            .to_vec(),
    )
    .unwrap();
    let mut cells = Vec::new();
    for z in 1..=6 {
        for y in 1..=6 {
            for x in 0..=3 {
                cells.push(MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    velocity_m_s: if moving {
                        [0.0, 0.003 * (x as f64 - 1.5), 0.0]
                    } else {
                        [0.0; 3]
                    },
                });
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let flow = plan_free_surface_from_material(&source, 0.001, DT, 0.1)
            .unwrap()
            .free_surface3_from_cells(
                MaterialFreeSurfaceDomain3 {
                    cells: [8; 3],
                    gravity_m_s2: [0.0; 3],
                    contact: ContactModel3::Neutral,
                    solid_tiles: &[],
                    max_cells: 512,
                },
                &cells,
                &cx,
            )
            .unwrap();
        let liquid = flow
            .with_enthalpy(
                &law,
                &caloric,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        run(&liquid, &gate)
    })
}

#[test]
fn g1_reaction_and_heat_match_independent_species_energy_and_pressure() {
    with_liquid(false, |liquid, gate| {
        let law = reaction(500.0, 20000.0, 8.0);
        let original = gas(&law);
        let bits = liquid.flow().lattice().grid().state_bits();
        let before = liquid.total_enthalpy_j().unwrap();
        let next = law
            .liquid_candidate(
                liquid,
                &original,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[face()],
                controls(),
                gate,
            )
            .unwrap();
        let report = next.report();
        // Equal molar masses and viscosities make WMS reduce to mole-weighting.
        // Both authored species have cp/R=3.5, hence cv/R=2.5. A starts at
        // unit activity and B at zero, so the first explicit extent is dt*V*rate.
        let extent = DT * VOLUME * 500.0;
        let amount = 100000.0 * VOLUME / (UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * 400.0);
        let capacity = 2.5 * 100000.0 * VOLUME / 400.0;
        let conductivity = 0.03 * (1.0 + extent / amount * 7.0);
        let conductance = 1e-6 / (0.0005 / 0.5 + 0.0005 / conductivity);
        let (mut tl, mut tg) = (400.0, 400.0 + 20000.0 * extent / capacity);
        let mut transferred = 0.0;
        for _ in 0..4 {
            let q = DT / 4.0 * conductance * (tg - tl);
            tl += q / 0.002;
            tg -= q / capacity;
            transferred += q;
        }
        near(report.duration_s, DT, 0.0);
        near(report.reaction.extent_mol, extent, 1e-20);
        near(report.reaction.species[0].final_mol, amount - extent, 1e-18);
        near(report.reaction.species[1].final_mol, extent, 1e-18);
        near(next.gas().state().state().temperature, tg, 2e-9);
        near(
            next.gas().state().state().pressure,
            100000.0 * tg / 400.0,
            1e-6,
        );
        near(
            next.gas().state().state().thermal_conductivity,
            conductivity,
            1e-14,
        );
        let environment = report.thermal.conduction.environment.as_ref().unwrap();
        near(
            environment.thermal_endpoint[0].liquid_temperature_k,
            tl,
            1e-10,
        );
        near(environment.heat_to_liquid_j, transferred, 1e-14);
        near(report.liquid_enthalpy_change_j, transferred, 2e-12);
        near(report.gas_internal_energy_change_j, -transferred, 2e-14);
        near(report.supplied_thermal_energy_j, 0.0, 0.0);
        assert!(report.combined_energy_residual_j.abs() < controls().total_energy_tolerance_j);
        assert_eq!(next.gas().mass_kg(), original.mass_kg());
        assert_eq!(next.gas().volume_m3(), original.volume_m3());
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), before);
        assert_eq!(original, gas(&law));
    });
}

#[test]
fn g3_rate_caloric_and_transport_data_reach_actual_liquid_momentum() {
    let run = |rate, release, conductivity| {
        with_liquid(true, |liquid, gate| {
            let law = reaction(rate, release, conductivity);
            let original = gas(&law);
            let next = law
                .liquid_candidate(
                    liquid,
                    &original,
                    &[],
                    &[],
                    &[],
                    &[],
                    &[],
                    &[face()],
                    controls(),
                    gate,
                )
                .unwrap();
            (
                next.liquid().flow().lattice().grid().state_bits(),
                next.liquid()
                    .material_at(EXPOSED)
                    .unwrap()
                    .dynamic_viscosity_pa_s(),
                next.report()
                    .thermal
                    .conduction
                    .environment
                    .as_ref()
                    .unwrap()
                    .heat_to_liquid_j,
            )
        })
    };
    let base = run(500.0, 20000.0, 1.0);
    for changed in [
        run(1000.0, 20000.0, 1.0),
        run(500.0, 40000.0, 1.0),
        run(500.0, 20000.0, 8.0),
    ] {
        assert!(
            changed.2 > base.2,
            "source change must increase delivered heat"
        );
        assert!(
            changed.1 < base.1,
            "same liquid mu(T) must see chemical heat"
        );
        assert_ne!(
            changed.0, base.0,
            "source changes must reach native flow populations"
        );
    }
}

#[test]
fn g1_continuation_keeps_cooled_inventory_and_counts_authored_and_viscous_heat_once() {
    with_liquid(true, |liquid, gate| {
        let law = reaction(500.0, 20000.0, 4.0);
        let original = gas(&law);
        let first = law
            .liquid_candidate(
                liquid,
                &original,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[face()],
                controls(),
                gate,
            )
            .unwrap();
        let expected = law
            .parcel_candidate(
                &ClosedGasParcel::from_isochoric(first.gas().clone()),
                DT,
                controls().reaction,
                || false,
            )
            .unwrap();
        let mut limits = controls();
        limits.with_viscous_heating = true;
        let heat = [FreeSurfaceCellHeat {
            coordinates: EXPOSED,
            energy_j: 1e-5,
        }];
        let next = law
            .liquid_candidate(
                first.liquid(),
                first.gas(),
                &[],
                &[],
                &heat,
                &[],
                &[],
                &[face()],
                limits,
                gate,
            )
            .unwrap();
        assert_eq!(next.report().reaction, *expected.step());
        let flow = &next.report().thermal.conduction.flow;
        assert!(flow.viscous_heat_j > 0.0);
        near(
            next.report().supplied_thermal_energy_j,
            1e-5 + flow.viscous_heat_j,
            1e-16,
        );
        assert!(next.report().combined_energy_residual_j.abs() < limits.total_energy_tolerance_j);
        assert!(next.report().reaction.species[1].initial_mol > 0.0);
        assert!(
            next.report().reaction.species[1].final_mol
                > next.report().reaction.species[1].initial_mol
        );
    });
}

#[test]
fn g1_finite_material_wall_shares_the_same_energy_transaction() {
    use fs_lbm::material::free_surface::enthalpy::conduction::{
        FreeSurfaceBodyContact, FreeSurfaceThermalBody,
    };
    let mut claims = ClaimSet::new();
    let conductivity = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    claims
        .insert_claim(PropertyClaim {
            key: conductivity.clone(),
            value: PropertyValue::Scalar {
                value: 0.5,
                dims: THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured finite solid".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured solid".into(),
            phase: "solid".into(),
            process: "authored".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        [300.0, 600.0]
            .map(|temperature_k| EnthalpyPhaseKnot {
                temperature_k,
                specific_enthalpy_j_kg: 500.0 * (temperature_k - 300.0),
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 2000.0,
            })
            .to_vec(),
    )
    .unwrap();
    let point = QueryPoint::new().with("T", 450.0).unwrap();
    let body = FreeSurfaceThermalBody::try_new(
        7,
        ThermalMaterialPoint {
            material: &card,
            context: &point,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &caloric,
        CellGeometry {
            start_m: [-0.001, 0.0035, 0.0035],
            end_m: [0.0, 0.0035, 0.0035],
            volume_m3: 1e-9,
            surface_area_m2: 6e-6,
        },
        75000.0,
        &conductivity,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let contact = FreeSurfaceBodyContact {
        body_id: 7,
        liquid: [0, 3, 3],
        wall: [-1, 3, 3],
        area_m2: 1e-6,
        distances_m: [0.0005; 2],
    };
    with_liquid(false, |liquid, gate| {
        let law = reaction(500.0, 20000.0, 1.0);
        let original = gas(&law);
        let initial_h = body.enthalpy_j();
        let bodies = [body.clone()];
        let next = law
            .liquid_candidate(
                liquid,
                &original,
                &bodies,
                &[],
                &[],
                &[contact],
                &[],
                &[face()],
                controls(),
                gate,
            )
            .unwrap();
        assert!(next.report().wall_enthalpy_change_j < 0.0);
        near(
            next.bodies()[0].enthalpy_j() - initial_h,
            next.report().wall_enthalpy_change_j,
            0.0,
        );
        let delivered = next.report().thermal.conduction.bodies.as_ref().unwrap();
        near(
            next.report().wall_enthalpy_change_j,
            -delivered.heat_to_liquid_j,
            1e-14,
        );
        near(
            next.report().liquid_enthalpy_change_j
                + next.report().wall_enthalpy_change_j
                + next.report().gas_internal_energy_change_j,
            0.0,
            2e-9,
        );
        assert_eq!(bodies[0].enthalpy_j(), initial_h);
        assert_eq!(body.enthalpy_j(), initial_h);
    });
}

#[test]
fn g4_late_energy_geometry_and_cancellation_refusals_preserve_every_original() {
    with_liquid(true, |liquid, gate| {
        let law = reaction(500.0, 20000.0, 1.0);
        let original = gas(&law);
        let bits = liquid.flow().lattice().grid().state_bits();
        let enthalpy = liquid.total_enthalpy_j().unwrap();
        let material = liquid.material_at(EXPOSED).unwrap().clone();
        let mut bad_face = face();
        bad_face.gas = [2, 3, 3]; // Occupied liquid is not an exposed gas face.
        let mut no_defect = controls();
        no_defect.total_energy_tolerance_j = 0.0;
        for (face, limits) in [(bad_face, controls()), (face(), no_defect)] {
            assert!(
                law.liquid_candidate(
                    liquid,
                    &original,
                    &[],
                    &[],
                    &[],
                    &[],
                    &[],
                    &[face],
                    limits,
                    gate
                )
                .is_err()
            );
            assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
            assert_eq!(liquid.total_enthalpy_j().unwrap(), enthalpy);
            assert_eq!(liquid.material_at(EXPOSED).unwrap(), &material);
            assert_eq!(original, gas(&law));
        }
        gate.request();
        assert!(matches!(
            law.liquid_candidate(
                liquid,
                &original,
                &[],
                &[],
                &[],
                &[],
                &[],
                &[face()],
                controls(),
                gate
            ),
            Err(ReactiveLiquidGasError::Cancelled)
        ));
        assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
    });
}
