//! G1/G3/G4: sourced tensile limits on native rods, thermal layers and bundles.
//! Synthetic materials verify equations and integration, not measured strengths.

use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::section::strength::AxialTensileCapacity;
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::section::{
    ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, integrate_isotropic_thermal_expansion, resolve_material_state_point,
};
use fs_qty::{Density, Pressure, Temperature};
use fs_solid::rod::material_field::strength::RodTensileMember;
use fs_solid::rod::material_field::supports::RodSupport;
use fs_solid::rod::material_field::tension::{MaterialTensionLink, TensionLinkAttachment};
use fs_solid::rod::material_field::{SpatialMaterialRod, SpatialRodError};
use fs_solid::rod::{TipLoad, WorldTipLoad};

const RADIUS: f64 = 0.05;
const LIMIT: &str = "specified_tensile_limit";

// Each entry is (Young's modulus, tensile limit at 300 K, constant alpha).
// Layers have equal areas, so the thermal mismatch oracle is independent of
// the section owner's force aggregation and can be written in closed form.
fn section(
    tag: u8,
    temperature: f64,
    data: &[(f64, f64, f64)],
) -> (ConcentricShearSection, ConcentricThermoelasticSection) {
    let point = QueryPoint::new().with("T", temperature).unwrap();
    let reference = QueryPoint::new().with("T", 300.0).unwrap();
    let mut layers = Vec::new();
    let mut paths = Vec::new();
    for (index, &(young, strength, alpha)) in data.iter().enumerate() {
        let mut claims = ClaimSet::new();
        for (name, value, dims) in [
            ("density", 1000.0, Density::DIMS),
            ("young_modulus", young, Pressure::DIMS),
            ("shear_modulus", 0.4 * young, Pressure::DIMS),
            (LIMIT, strength, Pressure::DIMS),
            (
                "linear_thermal_expansion_coefficient",
                alpha,
                INVERSE_TEMPERATURE_DIMS,
            ),
        ] {
            let (value, interpolation) = if name == LIMIT {
                (
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (400.0, 0.5 * value)],
                        dims,
                    },
                    InterpolationPolicy::LinearInside,
                )
            } else {
                (
                    PropertyValue::Scalar { value, dims },
                    InterpolationPolicy::ConstantWithinValidity,
                )
            };
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value,
                    interpolation,
                    validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0),
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: format!("synthetic rod tensile source {tag}:{index}"),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("synthetic rod layer {tag}:{index}"),
                phase: "solid".into(),
                process: "test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let requirements = [
            ("density", Density::DIMS),
            ("young_modulus", Pressure::DIMS),
            ("shear_modulus", Pressure::DIMS),
            (LIMIT, Pressure::DIMS),
        ]
        .map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
                .unwrap()
        });
        let material = resolve_material_state_point(
            &card,
            &point,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        paths.push(
            integrate_isotropic_thermal_expansion(
                &card,
                &reference,
                &point,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap(),
        );
        layers.push(ConcentricMaterialLayer {
            outer_radius_m: RADIUS * ((index + 1) as f64 / data.len() as f64).sqrt(),
            material,
        });
    }
    let elastic = ConcentricElasticSection::try_new(0.0, layers).unwrap();
    let thermal = ConcentricThermoelasticSection::try_new(elastic.clone(), paths).unwrap();
    (
        elastic.with_shear_modulus("shear_modulus", 0.8).unwrap(),
        thermal,
    )
}

fn capacity(section: &ConcentricShearSection) -> AxialTensileCapacity {
    section.elastic().axial_tensile_capacity(LIMIT).unwrap()
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "actual={actual}, expected={expected}"
    );
}

#[test]
fn g1_signed_native_axial_stress_uses_each_segments_own_strength() {
    let gate = CancelGate::new();
    let weak = section(1, 300.0, &[(1e6, 40.0, 1e-5)]).0;
    let strong = section(2, 300.0, &[(1e6, 80.0, 1e-5)]).0;
    let capacities = [capacity(&weak), capacity(&strong)];
    let mut rod = SpatialMaterialRod::try_straight(1.0, vec![weak, strong]).unwrap();
    for force in [0.5, -0.5] {
        rod.solve_static(
            &TipLoad {
                force: [force, 0.0, 0.0],
                moment: [0.0; 3],
            },
            1,
            1e-9,
            &gate,
        )
        .unwrap();
        let report = rod
            .assess_tensile_strength(&capacities, &[], 2, &gate)
            .unwrap();
        let stress = force / (core::f64::consts::PI * RADIUS * RADIUS);
        assert_eq!(report.samples().len(), 2);
        for (index, limit) in [40.0, 80.0].into_iter().enumerate() {
            let sample = &report.samples()[index];
            assert_eq!(sample.member(), RodTensileMember::Segment(index));
            assert_eq!(sample.capacity(), &capacities[index]);
            near(sample.axial_force_n(), force, 2e-7);
            let layer = sample.assessment().layers()[0];
            near(layer.minimum_axial_stress_pa, stress, 2e-7);
            near(layer.maximum_axial_stress_pa, stress, 2e-7);
            near(layer.tensile_utilization, stress.max(0.0) / limit, 2e-7);
        }
        assert_eq!(report.reaches_limit(), force > 0.0);
        assert_eq!(
            report.governing_sample().member(),
            RodTensileMember::Segment(0)
        );
    }
}

#[test]
fn g1_g3_native_bending_detects_surface_tension_at_zero_axial_resultant() {
    let gate = CancelGate::new();
    let source = section(3, 300.0, &[(1e6, 1000.0, 1e-5)]).0;
    let capacities = [capacity(&source)];
    let inertia = core::f64::consts::PI * RADIUS.powi(4) / 4.0;
    let expected_stress = 0.1 * RADIUS / inertia;
    for bending in [[0.0, 0.1], [0.06, 0.08], [-0.08, -0.06]] {
        let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source.clone()]).unwrap();
        rod.solve_static_world(
            &WorldTipLoad {
                force: [0.0; 3],
                moment: [0.0, bending[0], bending[1]],
            },
            1,
            1e-9,
            &gate,
        )
        .unwrap();
        let report = rod
            .assess_tensile_strength(&capacities, &[], 1, &gate)
            .unwrap();
        let sample = report.governing_sample();
        near(sample.axial_force_n(), 0.0, 1e-8);
        let layer = sample.assessment().layers()[0];
        near(layer.minimum_axial_stress_pa, -expected_stress, 2e-6);
        near(layer.maximum_axial_stress_pa, expected_stress, 2e-6);
        near(report.utilization(), expected_stress / 1000.0, 2e-6);
        assert!(report.reaches_limit());
    }
}

#[test]
fn g1_thermal_mismatch_retains_tensile_layers_and_the_original_strain_bound() {
    let gate = CancelGate::new();
    let (source, thermal) = section(4, 400.0, &[(1e6, 1000.0, 1e-5), (1e6, 1000.0, 3e-5)]);
    let capacities = [capacity(&source)];
    let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source])
        .unwrap()
        .with_thermal_sections(vec![thermal], 0.004, &gate)
        .unwrap();
    rod.solve_static(&TipLoad::default(), 1, 1e-9, &gate)
        .unwrap();
    let report = rod
        .assess_tensile_strength(&capacities, &[], 2, &gate)
        .unwrap();
    let sample = &report.samples()[0];
    near(sample.axial_force_n(), 0.0, 1e-8);
    near(sample.total_axial_strain(), 0.002, 1e-8);
    near(
        sample.assessment().layers()[0].maximum_axial_stress_pa,
        1000.0,
        1e-8,
    );
    near(
        sample.assessment().layers()[1].maximum_axial_stress_pa,
        -1000.0,
        1e-8,
    );
    assert_eq!(sample.capacity().layer_limits_pa(), &[500.0, 500.0]);
    near(report.utilization(), 2.0, 1e-8);
    assert_eq!(sample.assessment().governing_layer(), 0);
    assert!(report.reaches_limit());

    // Native axial/thermal strains stay in the declared band, but adding the
    // actual bending fiber strain leaves that band. Assessment must refuse.
    rod.solve_static_world(
        &WorldTipLoad {
            force: [0.0; 3],
            moment: [0.0, 0.0, 0.3],
        },
        1,
        1e-9,
        &gate,
    )
    .unwrap();
    assert!(rod.thermal_axial_response(0).unwrap().is_some());
    let before = rod.positions_m().to_vec();
    let energy = rod.energy_j().unwrap();
    assert!(matches!(
        rod.assess_tensile_strength(&capacities, &[], 2, &gate),
        Err(SpatialRodError::Material(_))
    ));
    assert_eq!(rod.positions_m(), before);
    assert_eq!(rod.energy_j().unwrap(), energy);
}

#[test]
fn g3_native_tension_bundles_use_force_per_strand_and_slack_has_zero_demand() {
    let gate = CancelGate::new();
    let source = section(5, 300.0, &[(1e6, 1e9, 1e-5)]).0;
    let strand = section(6, 300.0, &[(1e6, 1e4, 1e-5)]).0;
    let capacities = [capacity(&source)];
    let link_capacities = [capacity(&strand), capacity(&strand)];
    let bare = SpatialMaterialRod::try_straight(1.0, vec![source])
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(1.0), Some(0.0), Some(0.0)],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap();
    for count in [1, 7] {
        let links = [0.99, 1.01].map(|natural_length| {
            MaterialTensionLink::try_new(
                strand.elastic().clone(),
                [
                    TensionLinkAttachment {
                        node: 0,
                        offset_body_m: [0.0; 3],
                    },
                    TensionLinkAttachment {
                        node: 1,
                        offset_body_m: [0.0; 3],
                    },
                ],
                natural_length,
                count,
                0.05,
            )
            .unwrap()
        });
        let rod = bare.with_tension_links(links.into(), &gate).unwrap();
        let native = rod.tension_link_states().unwrap();
        let report = rod
            .assess_tensile_strength(&capacities, &link_capacities, 3, &gate)
            .unwrap();
        let taut = &report.samples()[1];
        let strain = (1.0 - 0.99) / 0.99;
        assert_eq!(taut.member(), RodTensileMember::TensionLink(0));
        near(
            taut.axial_force_n() * f64::from(count),
            native[0].tension_n,
            1e-14,
        );
        near(taut.total_axial_strain(), strain, 1e-14);
        near(
            taut.assessment().layers()[0].maximum_axial_stress_pa,
            1e6 * strain,
            1e-14,
        );
        near(report.utilization(), 1e6 * strain / 1e4, 1e-14);
        assert_eq!(
            report.governing_sample().member(),
            RodTensileMember::TensionLink(0)
        );
        assert_eq!(report.samples()[2].assessment().utilization(), 0.0);
        assert_eq!(report.samples()[2].axial_force_n(), 0.0);
        assert!(
            rod.assess_tensile_strength(&capacities, &[], 3, &gate)
                .is_err()
        );
        assert!(
            rod.assess_tensile_strength(&capacities, &link_capacities, 2, &gate)
                .is_err()
        );
        assert!(
            rod.assess_tensile_strength(
                &capacities,
                &[capacities[0].clone(), capacities[0].clone()],
                3,
                &gate,
            )
            .is_err()
        );
    }
}

#[test]
fn g3_g4_source_state_budget_and_cancellation_refuse_without_mutating_the_rod() {
    let gate = CancelGate::new();
    let data = [(1e6, 1000.0, 1e-5)];
    let source = section(7, 300.0, &data).0;
    let capacities = [capacity(&source)];
    let foreign = [capacity(&section(8, 300.0, &data).0)];
    let refreshed = [capacity(&section(7, 400.0, &data).0)];
    let mut rod = SpatialMaterialRod::try_straight(1.0, vec![source]).unwrap();
    rod.solve_static(
        &TipLoad {
            force: [0.5, 0.0, 0.0],
            moment: [0.0; 3],
        },
        1,
        1e-9,
        &gate,
    )
    .unwrap();
    let before = rod.positions_m().to_vec();
    let frames = rod.frames().to_vec();
    let energy = rod.energy_j().unwrap();
    let first = rod
        .assess_tensile_strength(&capacities, &[], 1, &gate)
        .unwrap();
    assert!(
        rod.assess_tensile_strength(&foreign, &[], 1, &gate)
            .is_err()
    );
    assert!(
        rod.assess_tensile_strength(&refreshed, &[], 1, &gate)
            .is_err()
    );
    assert!(rod.assess_tensile_strength(&[], &[], 1, &gate).is_err());
    assert!(
        rod.assess_tensile_strength(&capacities, &[], 0, &gate)
            .is_err()
    );
    assert_eq!(
        first,
        rod.assess_tensile_strength(&capacities, &[], 1, &gate)
            .unwrap()
    );
    gate.request();
    assert_eq!(
        rod.assess_tensile_strength(&capacities, &[], 1, &gate),
        Err(SpatialRodError::Cancelled)
    );
    assert_eq!(rod.positions_m(), before);
    assert_eq!(rod.frames(), frames);
    assert_eq!(rod.energy_j().unwrap(), energy);
}
