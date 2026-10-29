//! G1/G3: existing radiation, native heat and source mechanics share one state.
use super::*;
use fs_conduction::radiation::{GrayDiffuseEnclosure, RadiationSurface, ViewFactorMatrix};
use fs_conduction::transient::enthalpy::heterogeneous::{EnthalpyHeatStep, EnthalpyRadiation};
use fs_conduction::{
    AmbientRadiationConfig, AmbientRadiationPatch, EMISSIVITY_DIMS, PointwiseAmbientRadiation,
    PointwiseEnclosureRadiation, PointwiseEnclosureTolerance, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ThermalBc,
};

const CP: f64 = 10.0;
const RHO: f64 = 6.0;
const EMISSIVITY: f64 = 0.8;
const DT: f64 = 0.002;

fn regular_vertices() -> [[f64; 3]; 4] {
    [
        [0.0; 3],
        [1.0, 0.0, 0.0],
        [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
        [0.5, 3.0_f64.sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()],
    ]
}
fn emissivity(name: &str, upper: f64) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: EMISSIVITY,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 290.0, upper),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured gray boundary for native coupling".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray".into(),
            phase: "solid".into(),
            process: "numerical oracle".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    SurfaceEmissivity::from_card(name, &card, 300.0, SelectionPolicy::SingleClaimOnly).unwrap()
}
fn radiation_config() -> AmbientRadiationConfig {
    AmbientRadiationConfig {
        max_iterations: 128,
        temperature_tolerance_k: 1e-10,
        balance_tolerance_w: 1e-9,
        balance_relative_tolerance: 1e-12,
        relaxation: 0.5,
    }
}
fn root(mut lower: f64, mut upper: f64, residual: impl Fn(f64) -> f64) -> f64 {
    for _ in 0..80 {
        let value = 0.5 * (lower + upper);
        if residual(value) > 0.0 {
            upper = value;
        } else {
            lower = value;
        }
    }
    0.5 * (lower + upper)
}

#[test]
fn g3_no_radiation_selector_preserves_original_heat_and_solid_bits() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2.0, 1000.0, 1e-5, CP, 500.0, 750.0);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), RHO, cx);
        let context = point(300.0);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let boundary = boundary(&mesh);
        let source = ScalarField::uniform("finite generation", 60.0).unwrap();
        let old = [1000.0, 1200.0, 800.0, 1000.0];
        let original = model
            .advance(
                &old,
                DT,
                &boundary,
                None,
                &source,
                thermal_config(),
                &loads,
                &FREE,
                solid_config(),
                cx,
            )
            .unwrap();
        let selected = model
            .advance_with_radiation(
                &old,
                DT,
                &boundary,
                None,
                &source,
                thermal_config(),
                EnthalpyRadiation::None,
                &loads,
                &FREE,
                solid_config(),
                cx,
            )
            .unwrap();
        assert!(matches!(&selected.thermal, EnthalpyHeatStep::Native(_)));
        assert_eq!(
            format!("{:?}", original.thermal),
            format!("{:?}", selected.thermal.conduction())
        );
        assert_eq!(original.state.materials(), selected.state.materials());
        assert_eq!(
            original.equilibrium.deformed_mesh.nodes_m(),
            selected.equilibrium.deformed_mesh.nodes_m()
        );
    });
}

#[test]
fn g1_reservoir_radiation_heats_native_mass_and_changes_compatible_thermal_strain() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2.0, 1000.0, 1e-5, CP, 500.0, 750.0);
        let vertices = regular_vertices();
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
            vertices.to_vec(),
        )
        .unwrap();
        let native = native(&mesh, caloric.curve(), RHO, cx);
        let context = point(300.0);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("skin", |_| true, ThermalBc::neumann(0.0).unwrap())
            .unwrap()
            .finish()
            .unwrap();
        let patches =
            [AmbientRadiationPatch::new("skin", emissivity("skin", 500.0), 450.0).unwrap()];
        let pointwise = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let zero = ScalarField::uniform("no volumetric generation", 0.0).unwrap();
        let capacity = RHO * CP * 2.0_f64.sqrt() / 12.0;
        let area = 3.0_f64.sqrt();
        let expected = root(300.0, 450.0, |t| {
            capacity * (t - 300.0)
                + DT * area
                    * EMISSIVITY
                    * STEFAN_BOLTZMANN_W_M2_K4
                    * (t.powi(4) - 450.0_f64.powi(4))
        });
        for radiation in [
            EnthalpyRadiation::Ambient {
                patches: &patches,
                config: radiation_config(),
            },
            EnthalpyRadiation::PointwiseAmbient(&pointwise),
        ] {
            let result = model
                .advance_with_radiation(
                    &[0.0; 4],
                    DT,
                    &boundary,
                    None,
                    &zero,
                    thermal_config(),
                    radiation,
                    &loads,
                    &FREE,
                    solid_config(),
                    cx,
                )
                .unwrap();
            let thermal = result.thermal.conduction();
            for &temperature in &thermal.temperature {
                near(temperature, expected, 2e-8);
            }
            near(
                thermal.stored_energy_change_j,
                capacity * (expected - 300.0),
                1e-8,
            );
            near(result.thermal.physical_energy_residual_j(), 0.0, 1e-8);
            let radiative_heat = match &result.thermal {
                EnthalpyHeatStep::Ambient(step) => {
                    assert!(step.convective_robin_fluxes.is_empty());
                    near(step.convective_out_w, 0.0, 0.0);
                    step.radiation.nonlinear_radiation_out_w
                }
                EnthalpyHeatStep::Native(step) => {
                    near(step.robin_out_w, 0.0, 0.0);
                    near(step.neumann_out_w, 0.0, 0.0);
                    step.pointwise_radiation.as_ref().unwrap().outward_heat_w
                }
                _ => panic!("ambient selection must retain an ambient heat report"),
            };
            assert!(radiative_heat < 0.0);
            near(thermal.stored_energy_change_j, -DT * radiative_heat, 1e-8);
            let rise = expected - 300.0;
            let strain = 1e-5 * (rise + 0.005 * rise * rise);
            for (actual, reference) in result
                .equilibrium
                .deformed_mesh
                .nodes_m()
                .iter()
                .zip(vertices)
            {
                for axis in 0..3 {
                    near(actual[axis], reference[axis] * (1.0 + strain), 2e-9);
                }
            }
            near(
                result.state.samples()[0].phase.temperature_k(),
                expected,
                2e-8,
            );
            near(result.state.materials()[0].density_kg_m3, RHO, 0.0);
        }

        // Heat succeeds but a later geometric displacement gate refuses. The
        // unchanged source/history can still produce the same admitted step.
        let old = [0.0_f64; 4];
        let mut constrained = solid_config();
        constrained.maximum_nodal_displacement_m = 1e-15;
        assert!(
            model
                .advance_with_radiation(
                    &old,
                    DT,
                    &boundary,
                    None,
                    &zero,
                    thermal_config(),
                    EnthalpyRadiation::PointwiseAmbient(&pointwise),
                    &loads,
                    &FREE,
                    constrained,
                    cx
                )
                .is_err()
        );
        assert_eq!(old.map(f64::to_bits), [0.0_f64.to_bits(); 4]);
        let narrow =
            [AmbientRadiationPatch::new("skin", emissivity("skin", 300.01), 450.0).unwrap()];
        let narrow = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &narrow).unwrap();
        assert!(
            model
                .advance_with_radiation(
                    &old,
                    DT,
                    &boundary,
                    None,
                    &zero,
                    thermal_config(),
                    EnthalpyRadiation::PointwiseAmbient(&narrow),
                    &loads,
                    &FREE,
                    solid_config(),
                    cx
                )
                .is_err()
        );
        assert_eq!(old.map(f64::to_bits), [0.0_f64.to_bits(); 4]);
    });
}

#[test]
fn g1_closed_exchange_changes_both_material_states_without_external_heat() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2.0, 1000.0, 1e-5, CP, 500.0, 750.0);
        let mut positions = regular_vertices().to_vec();
        positions.extend(regular_vertices().map(|[x, y, z]| [x + 3.0, y, z]));
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(8, vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
            positions,
        )
        .unwrap();
        let native = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[ReferenceEnthalpyMaterial {
                curve: caloric.curve(),
                reference_density_kg_m3: RHO,
            }],
            &[0, 0],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let context = point(300.0);
        let model = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "cold",
                |f| f.vertices.iter().all(|&v| v < 4),
                ThermalBc::neumann(0.0).unwrap(),
            )
            .unwrap()
            .region(
                "hot",
                |f| f.vertices.iter().all(|&v| v >= 4),
                ThermalBc::neumann(0.0).unwrap(),
            )
            .unwrap()
            .finish()
            .unwrap();
        let cold = RadiationSurface::new(
            &mesh,
            "cold",
            |f| f.vertices.iter().all(|&v| v < 4),
            emissivity("cold", 500.0),
        )
        .unwrap();
        let hot = RadiationSurface::new(
            &mesh,
            "hot",
            |f| f.vertices.iter().all(|&v| v >= 4),
            emissivity("hot", 500.0),
        )
        .unwrap();
        // Explicit synthetic two-patch enclosure; geometric visibility is not inferred.
        let factors = ViewFactorMatrix::infinite_parallel_plates(cold.area_m2()).unwrap();
        let enclosure = GrayDiffuseEnclosure::new(vec![cold, hot], factors).unwrap();
        let pointwise = PointwiseEnclosureRadiation::new(
            cx,
            &mesh,
            &boundary,
            &enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let fixed: Vec<_> = (0..24).collect();
        let zero = ScalarField::uniform("closed finite bodies", 0.0).unwrap();
        let old = [0.0, 0.0, 0.0, 0.0, 1000.0, 1000.0, 1000.0, 1000.0];
        let capacity = RHO * CP * 2.0_f64.sqrt() / 12.0;
        let conductance = 3.0_f64.sqrt() * STEFAN_BOLTZMANN_W_M2_K4 / (2.0 / EMISSIVITY - 1.0);
        let expected_cold = root(300.0, 350.0, |t| {
            capacity * (t - 300.0) - DT * conductance * ((700.0 - t).powi(4) - t.powi(4))
        });
        let expected_hot = 700.0 - expected_cold;
        for radiation in [
            EnthalpyRadiation::Enclosure {
                enclosure: &enclosure,
                config: radiation_config(),
            },
            EnthalpyRadiation::PointwiseEnclosure(&pointwise),
        ] {
            let mut config = thermal_config();
            config.newton.linear_restart = 8;
            let result = model
                .advance_with_radiation(
                    &old,
                    DT,
                    &boundary,
                    None,
                    &zero,
                    config,
                    radiation,
                    &loads,
                    &fixed,
                    solid_config(),
                    cx,
                )
                .unwrap();
            let thermal = result.thermal.conduction();
            for &temperature in &thermal.temperature[..4] {
                near(temperature, expected_cold, 2e-8);
            }
            for &temperature in &thermal.temperature[4..] {
                near(temperature, expected_hot, 2e-8);
            }
            near(thermal.stored_energy_change_j, 0.0, 1e-8);
            near(result.thermal.physical_energy_residual_j(), 0.0, 1e-8);
            match &result.thermal {
                EnthalpyHeatStep::Enclosure(step) => {
                    near(step.convective_out_w, 0.0, 0.0);
                    near(step.prescribed_neumann_out_w, 0.0, 0.0);
                }
                EnthalpyHeatStep::Native(step) => {
                    let report = step.pointwise_enclosure_radiation.as_ref().unwrap();
                    near(report.outward_heat_w, 0.0, 1e-8);
                    assert!(report.patches[0].outward_heat_w < 0.0);
                    assert!(report.patches[1].outward_heat_w > 0.0);
                }
                _ => panic!("closed selection must retain its enclosure report"),
            }
            for (index, expected) in [expected_cold, expected_hot].into_iter().enumerate() {
                near(
                    result.state.samples()[index].phase.temperature_k(),
                    expected,
                    2e-8,
                );
                let rise = expected - 300.0;
                let strain = 1e-5 * (rise + 0.005 * rise * rise);
                let young = 1000.0 * (1.0 - 0.0025 * rise);
                let stress = -young * strain / (1.0 - 2.0 * 0.25);
                for &actual in
                    &result.equilibrium.solution.stress.elements()[index].stress_mandel_pa()[..3]
                {
                    near(actual, stress, 1e-7);
                }
            }
        }
        assert_eq!(old, [0.0, 0.0, 0.0, 0.0, 1000.0, 1000.0, 1000.0, 1000.0]);
    });
}
