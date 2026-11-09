//! Original-source analytic fixtures exercise the real mixed operator and heat owner.
use super::{DEPTH, binding, initial_thermal, near, open_boundaries, thermal_config, zero};
use fs_couple::stokes_motion::thermal_dilation::{
    StokesThermalDilationControls, StokesThermalDilationError,
    isobaric::{
        StokesIsobaricDilation, StokesIsobaricDilationBinding, StokesIsobaricDilationError,
        solve_at_thermal_state,
    },
};
use fs_evidence::ValidityDomain;
use fs_flux::ns::boundary::thermal_dilation::{
    CaloricDilationError,
    isobaric::{IsobaricDilationCell, IsobaricDilationError},
};
use fs_flux::{FluxParams, FluxSystem, StokesBoundary, StokesBoundaryError, TriMesh};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::SolidLiquidPhase;
use fs_material::thermal::transport::{
    expansion::{
        FreeExpansionSourceQuery, FreeExpansionStrainSources, PreparedFreeExpansionDensity,
    },
    isobaric::{
        IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
        expansion::{IsobaricFreeExpansion, liquid::LiquidIsobaricExpansion},
    },
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::Mesh2;
use fs_thermal::enthalpy::{CellId, DepositionSimulation, MaterialId, ThermalLaw};
use std::collections::BTreeMap;

const STRAIN: [(f64, f64); 2] = [(290.0, 0.0), (310.0, 0.02)];
const KINK: [(f64, f64); 3] = [(290.0, 0.0), (302.0, 0.024), (310.0, 0.016)];
const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64, p: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute(), t)
        .unwrap()
        .with_quantity("p", QuantitySpec::dimensional(Pressure::DIMS), p)
        .unwrap()
}
fn keys() -> [PropertyKey; 3] {
    [
        PropertyKey::new("cp", Dims([2, 0, -2, -1, 0, 0])),
        PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new("engineering-length", Dims::NONE),
    ]
}
fn sources(
    strains: &[(f64, f64)],
    datum: f64,
) -> (
    MaterialCard,
    IsobaricHeatCapacityLaw,
    PreparedFreeExpansionDensity,
) {
    let keys = keys();
    let mut claims = ClaimSet::new();
    for (i, knots) in [
        &[(290.0, 2.0), (300.0, 2.0), (310.0, 4.0)][..],
        &[(300.0, 1000.0), (300.0, 1000.0)][..],
        strains,
    ]
    .into_iter()
    .enumerate()
    {
        let value = if i == 1 {
            PropertyValue::Scalar {
                value: 1000.0,
                dims: keys[i].quantity().dims(),
            }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: knots.to_vec(),
                dims: keys[i].quantity().dims(),
            }
        };
        claims
            .insert_claim(PropertyClaim {
                key: keys[i].clone(),
                value,
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", absolute(), knots[0].0, knots[knots.len() - 1].0)
                    .with_quantity(
                        "p",
                        QuantitySpec::dimensional(Pressure::DIMS),
                        100_000.0,
                        200_000.0,
                    ),
                interpolation: if i == 1 {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored original Cp/free-expansion liquid fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "analytic liquid".into(),
            phase: "liquid".into(),
            process: "maintained-pressure free volume".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let (cp, density) = prepare(&card, datum, 100_000.0);
    (card, cp, density)
}
fn prepare(
    card: &MaterialCard,
    datum: f64,
    pressure: f64,
) -> (IsobaricHeatCapacityLaw, PreparedFreeExpansionDensity) {
    let keys = keys();
    let (lower, upper, reference) = (
        point(290.0, pressure),
        point(310.0, pressure),
        point(300.0, pressure),
    );
    let cp = IsobaricHeatCapacityLaw::try_new(
        card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &keys[0],
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &lower,
            upper: &upper,
            selection: SELECTION,
        },
        datum,
    )
    .unwrap();
    let density = PreparedFreeExpansionDensity::try_new(
        card,
        FreeExpansionSourceQuery {
            density: &keys[1],
            density_selection: SELECTION,
            reference: &reference,
            strains: FreeExpansionStrainSources::Isotropic {
                strain: &keys[2],
                selection: SELECTION,
            },
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
        },
    )
    .unwrap();
    (cp, density)
}
fn solve(
    system: &FluxSystem<'_>,
    thermal: &DepositionSimulation<'_>,
    ids: &[CellId],
    caloric: &LiquidIsobaricExpansion<'_>,
    boundaries: &[StokesBoundary<'_>],
) -> Result<StokesIsobaricDilation, StokesIsobaricDilationError> {
    let laws = BTreeMap::from([(MaterialId(9), caloric)]);
    solve_at_thermal_state(
        system,
        thermal,
        StokesIsobaricDilationBinding {
            cells: binding(ids),
            depth_m: DEPTH,
            caloric_by_material: &laws,
        },
        FluxParams::default(),
        &|_, _| 2.0,
        &zero,
        boundaries,
        StokesThermalDilationControls {
            maximum_dilation_defect_per_s: 2e-11,
            maximum_volume_rate_defect_m3_s: 2e-11,
        },
        || false,
    )
}

#[test]
fn native_measured_source_heat_and_directional_cooling_produce_exact_outlet_velocity() {
    let negative = [(290.0, 0.02), (310.0, 0.0)];
    let cases: [(&[(f64, f64)], f64, f64, f64); 6] = [
        (&STRAIN, 305.0, 6.0, 6.0 * 0.001 / 1.015),
        (&STRAIN, 305.0, -6.0, -6.0 * 0.001 / 1.015),
        (&negative, 305.0, 6.0, -6.0 * 0.001 / 1.005),
        (&KINK, 302.0, -2.4, -3.0 * 0.002 / 1.024),
        (&KINK, 302.0, 2.4, -3.0 * 0.001 / 1.024),
        (&STRAIN, 310.0, 0.0, 0.0),
    ];
    for n in [1, 2] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n));
        let system = FluxSystem::new(&mesh);
        let boundaries = open_boundaries(&mesh);
        for &(strains, temperature, rate, expected) in &cases {
            let (card, cp, density) = sources(strains, 300.0);
            let original = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
            let caloric =
                LiquidIsobaricExpansion::try_new(&original, SolidLiquidPhase::Liquid).unwrap();
            let h = cp.specific_enthalpy_j_kg(temperature).unwrap();
            let rho = density.at_temperature(temperature).unwrap().density_kg_m3();
            let (thermal, ids) = initial_thermal(&mesh, &original, false, |_, cell| {
                cell.specific_h_j_kg = h;
                cell.mass_kg = rho * cell.geometry.volume_m3;
                cell.boundary.source_w = cell.mass_kg * rate;
            });
            let before = thermal.cells().clone();
            let answer = solve(&system, &thermal, &ids, &caloric, &boundaries).unwrap();
            near(answer.expected_volume_rate_m3_s, DEPTH * expected, 1e-14);
            near(answer.boundary_volume_rate_m3_s, DEPTH * expected, 2e-11);
            for point in [[0.2, 0.3], [0.7, 0.4]] {
                let velocity = system.velocity_at(&answer.flow.flow.field.x, point);
                near(velocity[0], expected * point[0], 2e-11);
                near(velocity[1], 0.0, 2e-11);
            }
            for (cell, source) in answer.cells.iter().zip(&answer.flow.source) {
                assert_eq!(cell.caloric_identity, *original.identity().as_bytes());
                assert_eq!(source.caloric_identity, original.identity());
                near(source.dilation_per_s, expected, 1e-15);
                near(
                    cell.volume_rate_m3_s,
                    before[&cell.cell].geometry.volume_m3 * expected,
                    1e-14,
                );
            }
            assert_eq!(thermal.cells(), &before);
            assert_eq!(thermal.time_s(), 0.0);
        }
    }
}

#[test]
fn native_conduction_volume_uses_original_variable_cp_and_cubic_j_without_external_heat() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let boundaries = open_boundaries(&mesh);
    let mut first_field: Option<Vec<f64>> = None;
    for datum in [300.0, 295.0] {
        let (card, cp, density) = sources(&STRAIN, datum);
        let original = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
        let caloric =
            LiquidIsobaricExpansion::try_new(&original, SolidLiquidPhase::Liquid).unwrap();
        let (thermal, ids) = initial_thermal(&mesh, &original, true, |index, cell| {
            let t = if index == 0 { 295.0 } else { 305.0 };
            cell.specific_h_j_kg = cp.specific_enthalpy_j_kg(t).unwrap();
            cell.mass_kg =
                density.at_temperature(t).unwrap().density_kg_m3() * cell.geometry.volume_m3;
        });
        let answer = solve(&system, &thermal, &ids, &caloric, &boundaries).unwrap();
        // Native G=0.2 W/K transfers 2 W. For the original affine length law,
        // dv/dh = 3*(1+epsilon)^2*epsilon_T/(rho_ref*(1+epsilon_ref)^3*Cp).
        let cold_slope = 3.0 * 1.005_f64.powi(2) * 0.001 / (1000.0 * 1.01_f64.powi(3) * 2.0);
        let hot_slope = 3.0 * 1.015_f64.powi(2) * 0.001 / (1000.0 * 1.01_f64.powi(3) * 3.0);
        near(answer.heat.net_power_w[&ids[0]], 2.0, 1e-13);
        near(answer.heat.net_power_w[&ids[1]], -2.0, 1e-13);
        near(answer.heat.external_power_w, 0.0, 1e-13);
        near(answer.cells[0].volume_rate_m3_s, 2.0 * cold_slope, 1e-17);
        near(answer.cells[1].volume_rate_m3_s, -2.0 * hot_slope, 1e-17);
        let volume_rate = 2.0 * (cold_slope - hot_slope);
        assert!(volume_rate > 0.0);
        near(answer.expected_volume_rate_m3_s, volume_rate, 1e-17);
        near(answer.boundary_volume_rate_m3_s, volume_rate, 2e-11);
        if let Some(first) = &first_field {
            for (&actual, &expected) in answer.flow.flow.field.x.iter().zip(first) {
                near(actual, expected, 2e-11);
            }
        } else {
            first_field = Some(answer.flow.flow.field.x.clone());
        }
    }
}

#[test]
fn original_source_identity_current_volume_and_sealed_compatibility_are_enforced() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let boundaries = open_boundaries(&mesh);
    let (card, cp, density) = sources(&STRAIN, 300.0);
    let original = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let caloric = LiquidIsobaricExpansion::try_new(&original, SolidLiquidPhase::Liquid).unwrap();
    let h = cp.specific_enthalpy_j_kg(305.0).unwrap();
    let (thermal, ids) = initial_thermal(&mesh, &original, false, |_, cell| {
        cell.specific_h_j_kg = h;
        cell.mass_kg = original.bulk_density_kg_m3(h).unwrap().unwrap() * cell.geometry.volume_m3;
        cell.boundary.source_w = cell.mass_kg * 6.0;
    });
    let closed: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, _)| mesh.is_boundary(*edge))
        .map(|(edge, _)| StokesBoundary::Velocity {
            edge,
            value: &zero as &dyn Fn([f64; 2]) -> [f64; 2],
        })
        .collect();
    assert!(matches!(
        solve(&system, &thermal, &ids, &caloric, &closed),
        Err(StokesIsobaricDilationError::Flow(
            IsobaricDilationError::Flow(CaloricDilationError::Flow(
                StokesBoundaryError::IncompatibleFlux { .. }
            ))
        ))
    ));
    let (other_cp, other_density) = prepare(&card, 295.0, 100_000.0);
    let other = IsobaricFreeExpansion::try_new(&card, &other_cp, &other_density).unwrap();
    let other_liquid = LiquidIsobaricExpansion::try_new(&other, SolidLiquidPhase::Liquid).unwrap();
    assert!(matches!(
        solve(&system, &thermal, &ids, &other_liquid, &boundaries),
        Err(StokesIsobaricDilationError::Admission(
            StokesThermalDilationError::Binding(_)
        ))
    ));
    let warmed = thermal
        .advance_to(0.1, thermal_config(), || Ok(()))
        .unwrap()
        .0;
    assert!(matches!(
        solve(&system, &warmed, &ids, &caloric, &boundaries),
        Err(StokesIsobaricDilationError::Admission(
            StokesThermalDilationError::Binding(_)
        ))
    ));
}

#[test]
fn connected_sources_on_different_isobars_refuse_before_native_solve() {
    let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 1, 1));
    let system = FluxSystem::new(&mesh);
    let (card, cp, density) = sources(&STRAIN, 300.0);
    let first = IsobaricFreeExpansion::try_new(&card, &cp, &density).unwrap();
    let first_liquid = LiquidIsobaricExpansion::try_new(&first, SolidLiquidPhase::Liquid).unwrap();
    let (other_cp, other_density) = prepare(&card, 300.0, 200_000.0);
    let second = IsobaricFreeExpansion::try_new(&card, &other_cp, &other_density).unwrap();
    let second_liquid =
        LiquidIsobaricExpansion::try_new(&second, SolidLiquidPhase::Liquid).unwrap();
    let cells = [&first_liquid, &second_liquid].map(|caloric| IsobaricDilationCell {
        caloric,
        specific_enthalpy_j_kg: 0.0,
        material_specific_enthalpy_rate_w_kg: 0.0,
    });
    assert!(matches!(
        system.assemble_stokes_isobaric_dilation(
            FluxParams::default(),
            &|_, _| 2.0,
            &zero,
            &open_boundaries(&mesh),
            &cells,
        ),
        Err(IsobaricDilationError::Flow(CaloricDilationError::Invalid(
            _
        )))
    ));
}
