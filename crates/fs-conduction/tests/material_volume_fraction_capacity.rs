//! Constituent mass, graded composition and source domains in thermal storage.

mod support;

use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{
    HeatSourceHistory, ThermalHistoryCapacity, ThermalHistoryMethod, ThermalHistoryProblem,
    ThermalHistoryRun,
};
use fs_conduction::transient::backward_euler::{BackwardEuler, StepConfig};
use fs_conduction::transient::{
    HeatCapacity, MaterialHeatCapacity, MaterialHeatCapacityVolumeFraction, TransientConfig,
    march_with_capacity,
};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig,
    ScalarField, ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY,
};
use fs_qty::Density;
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use support::{with_cancelled_cx, with_cx};

fn material_capacity(rho: f64, cp: f64, low: f64, high: f64) -> MaterialHeatCapacity {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (DENSITY_PROPERTY, rho, Density::DIMS),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            cp,
            SPECIFIC_HEAT_CAPACITY_DIMS,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", low, high),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "manufactured constituent; not measured data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "composite storage regression".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    MaterialHeatCapacity::from_card(
        &card,
        &QueryPoint::new().with("T", 0.5 * (low + high)).unwrap(),
        low,
        high,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn mesh(two: bool) -> ConductionMesh {
    let mut points = vec![[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    let mut tets = vec![[0, 1, 2, 3]];
    if two {
        points.extend([[2., 0., 0.], [3., 0., 0.], [2., 1., 0.], [2., 0., 1.]]);
        tets.push([4, 5, 6, 7]);
    }
    ConductionMesh::new(TetComplex::from_tets(points.len(), tets), points).unwrap()
}

fn controls() -> StepConfig {
    StepConfig {
        linear: LinearConfig {
            tolerance: 1e-12,
            max_iterations: 512,
            restart: 16,
        },
        energy_tolerance_j: 2e-6,
    }
}

fn constituent<'a>(
    source: &'a MaterialHeatCapacity,
    fractions: &'a [f64],
) -> MaterialHeatCapacityVolumeFraction<'a> {
    MaterialHeatCapacityVolumeFraction {
        source,
        nodal_volume_fraction: fractions,
    }
}

#[test]
fn constituent_masses_set_the_heating_rate_and_preparation_is_immutable() {
    let mesh = mesh(false);
    let a = material_capacity(1000.0, 1000.0, 280.0, 400.0);
    let b = material_capacity(8000.0, 500.0, 280.0, 400.0);
    let half = [0.5; 4];
    let constituents = [constituent(&a, &half), constituent(&b, &half)];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let conductivity = ConductivityModel::isotropic_declared(20.0).unwrap();
    let source = ScalarField::uniform("volumetric heating", 100_000.0).unwrap();
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &conductivity,
        element_materials: None,
        source: &source,
    };
    // Per m3, the two fixed constituent masses are 500 kg and 4000 kg.
    // Their heat capacities add to 2.5e6 J/K. Averaging rho and Cp separately
    // would give 3.375e6 J/K and the wrong temperature response.
    let expected: f64 = 300.0 + 100_000.0 * 2.0 / (500.0 * 1000.0 + 4000.0 * 500.0);
    let wrong: f64 = 300.0 + 100_000.0 * 2.0 / (4500.0 * 750.0);
    assert!((expected - wrong).abs() > 0.02);
    with_cx(|cx| {
        for theta in [0.5, 1.0] {
            let config = TransientConfig::new(theta, 0.5, controls().linear).unwrap();
            let solved = march_with_capacity(
                cx,
                problem,
                HeatCapacity::MaterialVolumeFractions(&constituents),
                &config,
                &[300.; 4],
                4,
            )
            .unwrap();
            assert!(
                solved
                    .temperature
                    .iter()
                    .all(|&t| (t - expected).abs() < 1e-9)
            );
        }
        let mut fraction_a = half;
        let mut fraction_b = half;
        let stepper = {
            let input = [constituent(&a, &fraction_a), constituent(&b, &fraction_b)];
            BackwardEuler::material_volume_fractions(cx, &mesh, &input).unwrap()
        };
        // The original input can change after setup: no subsequent step may
        // resample composition or assemble a different reference inventory.
        fraction_a.fill(1.0);
        fraction_b.fill(0.0);
        assert_eq!(fraction_a[0] + fraction_b[0], 1.0);
        let mut temperature = vec![300.; 4];
        for _ in 0..2 {
            temperature = stepper
                .advance(cx, problem, None, &temperature, 1.0, controls())
                .unwrap()
                .temperature;
        }
        assert!(temperature.iter().all(|&t| (t - expected).abs() < 1e-9));
        assert_eq!(stepper.material_capacities(), &[a.clone(), b.clone()]);

        let history = ThermalBoundaryHistory::new(
            cx,
            ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder(),
            &Environment::earth_lab(),
            Vec::new(),
            0.0,
        )
        .unwrap();
        let mut run = ThermalHistoryRun::new_with_capacity(
            cx,
            ThermalHistoryProblem {
                boundary: &history,
                material: &conductivity,
                element_materials: None,
                interfaces: None,
                source: HeatSourceHistory::Fixed(source.clone()),
            },
            ThermalHistoryCapacity::MaterialVolumeFractions(&constituents),
            &[300.; 4],
            2,
            ThermalHistoryMethod::Linear,
        )
        .unwrap();
        run.advance_to(cx, 0.75, controls()).unwrap();
        let solved = run.advance_to(cx, 2.0, controls()).unwrap();
        assert!(
            solved
                .temperature
                .iter()
                .all(|&t| (t - expected).abs() < 1e-9)
        );
        assert_eq!(run.material_capacities(), &[a.clone(), b.clone()]);
    });
}

// Independent degree-two quadrature of sum(phi_i rho_i Cp_i T_h).
fn constituent_heat(mesh: &ConductionMesh, phi: &[f64], temperature: &[f64]) -> f64 {
    let mut total = 0.0;
    for (element, tet) in mesh.complex().tets.iter().enumerate() {
        for point in 0..4 {
            let shape: [f64; 4] = std::array::from_fn(|a| {
                if a == point {
                    0.585_410_196_624_968_5
                } else {
                    0.138_196_601_125_010_5
                }
            });
            let fraction: f64 = tet
                .iter()
                .zip(shape)
                .map(|(&v, n)| n * phi[v as usize])
                .sum();
            let t: f64 = tet
                .iter()
                .zip(shape)
                .map(|(&v, n)| n * temperature[v as usize])
                .sum();
            total += mesh.element_volume(element)
                * 0.25
                * (fraction * 1000.0 * 1000.0 + (1.0 - fraction) * 8000.0 * 500.0)
                * t;
        }
    }
    total
}

#[test]
fn graded_composition_preserves_constituent_energy_and_uniform_specific_heating() {
    let mesh = mesh(false);
    let a = material_capacity(1000.0, 1000.0, 280.0, 400.0);
    let b = material_capacity(8000.0, 500.0, 280.0, 400.0);
    let phi = [0.0, 1.0, 0.25, 0.75];
    let complement = phi.map(|value| 1.0 - value);
    let constituents = [constituent(&a, &phi), constituent(&b, &complement)];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let conductivity = ConductivityModel::isotropic_declared(10_000.0).unwrap();
    let zero = ScalarField::uniform("no generation", 0.0).unwrap();
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &conductivity,
        element_materials: None,
        source: &zero,
    };
    with_cx(|cx| {
        let stepper = BackwardEuler::material_volume_fractions(cx, &mesh, &constituents).unwrap();
        let initial = [300.0, 340.0, 310.0, 320.0];
        let solved = stepper
            .advance(cx, problem, None, &initial, 2.0, controls())
            .unwrap();
        assert!(
            (constituent_heat(&mesh, &phi, &solved.temperature)
                - constituent_heat(&mesh, &phi, &initial))
            .abs()
                < 2e-6
        );
        assert!(
            initial
                .iter()
                .zip(&solved.temperature)
                .any(|(old, new)| (old - new).abs() > 1.0)
        );
        assert!(solved.energy_residual_j.abs() <= controls().energy_tolerance_j);

        let heating = ScalarField::nodal(
            "composition-matched heating",
            4,
            phi.iter()
                .map(|fraction| 0.2 * (fraction * 1e6 + (1.0 - fraction) * 4e6))
                .collect(),
        )
        .unwrap();
        let heating_problem = ConductionProblem {
            source: &heating,
            ..problem
        };
        let config = TransientConfig::new(0.75, 0.5, controls().linear).unwrap();
        let solved = march_with_capacity(
            cx,
            heating_problem,
            HeatCapacity::MaterialVolumeFractions(&constituents),
            &config,
            &[300.; 4],
            4,
        )
        .unwrap();
        assert!(solved.temperature.iter().all(|&t| (t - 300.4).abs() < 1e-9));
        assert!(
            (constituent_heat(&mesh, &phi, &solved.temperature)
                - constituent_heat(&mesh, &phi, &[300.; 4])
                - constituent_heat(&mesh, &phi, &[0.4; 4]))
            .abs()
                < 2e-6
        );
    });
}

#[test]
fn absent_constituents_leave_domains_free_but_mixed_interiors_require_every_source() {
    let mesh = mesh(true);
    let a = material_capacity(1000.0, 1000.0, 280.0, 400.0);
    let b = material_capacity(8000.0, 500.0, 280.0, 310.0);
    let pure_a = [1., 1., 1., 1., 0., 0., 0., 0.];
    let pure_b = pure_a.map(|value| 1.0 - value);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let conductivity = ConductivityModel::isotropic_declared(20.0).unwrap();
    let source = ScalarField::uniform("no generation", 0.0).unwrap();
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &conductivity,
        element_materials: None,
        source: &source,
    };
    with_cx(|cx| {
        let initial = [300., 350., 350., 350., 300., 300., 300., 300.];
        let separated = [constituent(&a, &pure_a), constituent(&b, &pure_b)];
        let stepper = BackwardEuler::material_volume_fractions(cx, &mesh, &separated).unwrap();
        stepper
            .advance(cx, problem, None, &initial, 0.1, controls())
            .unwrap();

        let mut mixed_a = pure_a;
        let mut mixed_b = pure_b;
        mixed_a[0] = 0.9;
        mixed_b[0] = 0.1;
        let mixed = [constituent(&a, &mixed_a), constituent(&b, &mixed_b)];
        let mixed_stepper = BackwardEuler::material_volume_fractions(cx, &mesh, &mixed).unwrap();
        // B's sole nonzero corner is at 300 K; its fraction is still positive
        // inside the tet near the 350 K zero-fraction corners. Those must fail.
        assert!(matches!(
            mixed_stepper.advance(cx, problem, None, &initial, 0.1, controls()),
            Err(ConductionError::OutsideTemperatureSpan { .. })
        ));

        let all_a = [1.; 8];
        let absent_b = [0.; 8];
        let inactive = [constituent(&a, &all_a), constituent(&b, &absent_b)];
        let inactive_stepper =
            BackwardEuler::material_volume_fractions(cx, &mesh, &inactive).unwrap();
        inactive_stepper
            .advance(cx, problem, None, &[350.; 8], 0.1, controls())
            .unwrap();
        assert_eq!(
            inactive_stepper.material_capacities(),
            &[a.clone(), b.clone()]
        );
    });
}

#[test]
fn incomplete_nonphysical_or_cancelled_compositions_publish_no_storage() {
    let mesh = mesh(false);
    let source = material_capacity(1000.0, 1000.0, 280.0, 400.0);
    let half = [0.5; 4];
    with_cx(|cx| {
        assert!(BackwardEuler::material_volume_fractions(cx, &mesh, &[]).is_err());
        let short = [
            constituent(&source, &half[..3]),
            constituent(&source, &half),
        ];
        assert!(matches!(
            BackwardEuler::material_volume_fractions(cx, &mesh, &short),
            Err(ConductionError::FieldLength { .. })
        ));
        for bad in [[0.4; 4], [-0.1; 4], [1.1; 4], [f64::NAN; 4]] {
            let invalid = [constituent(&source, &bad), constituent(&source, &half)];
            assert!(BackwardEuler::material_volume_fractions(cx, &mesh, &invalid).is_err());
        }
    });
    with_cancelled_cx(|cx| {
        let valid = [constituent(&source, &half), constituent(&source, &half)];
        assert!(matches!(
            BackwardEuler::material_volume_fractions(cx, &mesh, &valid),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}
