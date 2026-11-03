//! Original source thermal potentials through native plane continuum solves.
//! Synthetic source data; no calibration or measured-material claim.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::thermal::{AffineThermoelastic, ThermoelasticSymmetry};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::state_point::*;
use fs_material::{SmallStrainLaw, tensor};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem,
};

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", temperature_quantity(), t).unwrap()
}

fn potential(isotropic: bool, axial_only: bool) -> AffineThermoelastic {
    let mut claims = ClaimSet::new();
    let mut insert = |name: &str, dims: Dims, base: f64, slope: f64| {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, base - 20.0 * slope), (360.0, base + 60.0 * slope)],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", temperature_quantity(), 280.0, 360.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic source thermal mechanics".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    };
    let inverse_t = Dims([0, 0, 0, -1, 0, 0]);
    let mut elastic_names = Vec::new();
    let mut expansion_names = Vec::new();
    if isotropic {
        insert(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 12e9, -3e7);
        insert(POISSON_RATIO_PROPERTY, Dims::NONE, 0.23, 0.0);
        insert(LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, inverse_t, 1e-5, 1e-8);
        elastic_names.extend([
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS),
            (POISSON_RATIO_PROPERTY, Dims::NONE),
        ]);
        expansion_names.push((LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, inverse_t));
    } else {
        for i in 0..3 {
            for (name, dims, value, slope) in [
                (ORTHOTROPIC_YOUNG_MODULUS_PROPERTIES[i], Pressure::DIMS,
                    [12e9, 5e9, 3e9][i], [-3e7, 1e7, 0.0][i]),
                (ORTHOTROPIC_POISSON_RATIO_PROPERTIES[i], Dims::NONE,
                    if axial_only { 0.0 } else { [0.1, 0.08, 0.15][i] }, 0.0),
                (ORTHOTROPIC_SHEAR_MODULUS_PROPERTIES[i], Pressure::DIMS,
                    [3e9, 2e9, 1e9][i], [0.0, -1e6, 0.0][i]),
            ] {
                insert(name, dims, value, slope);
                elastic_names.push((name, dims));
            }
            let alpha = if axial_only { [1e-5, 0.0, 0.0] } else { [1e-5, -0.4e-5, 2e-5] };
            insert(ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[i], inverse_t,
                alpha[i], [1e-8, 0.0, 0.0][i]);
            expansion_names.push((ORTHOTROPIC_THERMAL_EXPANSION_COEFFICIENT_PROPERTIES[i], inverse_t));
        }
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "native thermal plane test".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap();
    let resolve = |names: Vec<(&str, Dims)>| {
        let requirements: Vec<_> = names.into_iter().map(|(name, dims)| {
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
        }).collect();
        resolve_material_state_point(&card, &point(300.0), &requirements,
            MaterialPropertySelection::SingleClaimOnly).unwrap()
    };
    let elastic = resolve(elastic_names);
    let expansion = resolve(expansion_names);
    let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let (stiffness, expansion_symmetry) = if isotropic {
        (ThermoelasticSymmetry::Isotropic, ThermoelasticSymmetry::Isotropic)
    } else {
        (
            ThermoelasticSymmetry::OrientedOrthotropic {
                principal_to_target: if axial_only { identity } else {
                    tensor::rotation([1.0, 1.0, 2.0].map(|x| x / 6.0_f64.sqrt()), 0.41)
                },
            },
            ThermoelasticSymmetry::OrientedOrthotropic {
                principal_to_target: if axial_only { identity } else {
                    tensor::rotation([0.0, 1.0, 1.0].map(|x| x / 2.0_f64.sqrt()), 0.54)
                },
            },
        )
    };
    AffineThermoelastic::new(
        &card, &elastic, &expansion, stiffness, expansion_symmetry, [280.0, 360.0], 0.02,
    ).unwrap()
}

fn close(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() <= 2e-9 * scale.max(a.abs()).max(b.abs()), "{a:e} != {b:e}");
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey { seed: 123, kernel_id: 123, tile: 0, iteration: 0 },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

#[test]
fn g1_frozen_mechanics_retains_original_thermal_potential_and_tangent() {
    for isotropic in [true, false] {
        let potential = potential(isotropic, false);
        for t in [280.0, 300.0, 340.0, 360.0] {
            let source = potential.at_temperature(t).unwrap();
            let law = source.try_mechanical_law().unwrap();
            assert!(std::ptr::eq(law.source(), &source));
            let e = [0.002, -0.0007, 0.0004, 0.0003, -0.0002, 0.0001];
            let state = source.state(tensor::to_mandel(&e), 1.0).unwrap();
            let s = law.stress(&e, &());
            let c = law.tangent(&e, &());
            for (a, b) in tensor::to_mandel(&s).into_iter().zip(state.stress_mandel_pa) {
                close(a, b, 1e7);
            }
            for j in 0..6 {
                let h = 1e-6;
                let mut plus = e;
                let mut minus = e;
                plus[j] += h;
                minus[j] -= h;
                let sp = law.stress(&plus, &());
                let sm = law.stress(&minus, &());
                for i in 0..6 {
                    close((sp[i] - sm[i]) / (2.0 * h), c[i][j], 1e10);
                }
                let energy = |x| source.state(tensor::to_mandel(&x), 1.0).unwrap().free_energy_j;
                close((energy(plus) - energy(minus)) / (2.0 * h),
                    s[j] * if j < 3 { 1.0 } else { 2.0 }, 1e7);
            }
        }
    }
    let potential = potential(true, false);
    let hot = potential.at_temperature(360.0).unwrap();
    let law = hot.try_mechanical_law().unwrap();
    for bad in [
        [0.021, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.0195, 0.0, 0.0, 0.0, 0.0, 0.0],
    ] {
        assert!(hot.state(tensor::to_mandel(&bad), 1.0).is_err());
        assert!(law.stress(&bad, &()).iter().all(|s| s.is_nan()));
        assert!(law.tangent(&bad, &()).iter().flatten().all(|s| s.is_nan()));
    }
}

#[test]
fn g1_native_p1_q1_plane_stress_recovers_free_anisotropic_thermal_expansion() {
    let potential = potential(false, false);
    let hot = potential.at_temperature(350.0).unwrap();
    let law = hot.try_mechanical_law().unwrap();
    let free = tensor::from_mandel(&hot.free_strain_mandel());
    assert!(free[4].abs() > 1e-6 && free[5].abs() > 1e-6,
        "exercise transverse thermal shears");
    with_cx(|cx| {
        for mesh in [Mesh2::triangles(2.0, 1.0, 2, 2), Mesh2::quads(2.0, 1.0, 2, 2)] {
            let right = mesh.nodes.iter().position(|&p| p == [2.0, 0.0]).unwrap();
            let supports = [(0, 0.0), (1, 0.0), (2 * right + 1, 2.0 * free[3])];
            let problem = SmallStrainPlaneProblem {
                mesh: &mesh,
                laws: SmallStrainLawField::Uniform(&law),
                plane: PlaneKind::Stress,
                formulation: Formulation::Standard,
                thickness_m: 0.2,
                reduction: ReductionConfig {
                    max_iterations: 8, stress_tolerance_pa: 1e-6, strain_limit: 0.02,
                },
                rotation_limit_rad: 0.02,
                prescribed_m: &supports,
                nodal_forces_n: &[],
                load_path: &[1.0],
                budget: PlaneAssemblyBudget::default(),
                settings: IncrementSettings::default(),
            };
            let solved = problem.solve(cx).unwrap();
            for (&p, u) in mesh.nodes.iter().zip(&solved.displacement_m) {
                close(u[0], free[0] * p[0] + free[3] * p[1], 0.01);
                close(u[1], free[3] * p[0] + free[1] * p[1], 0.01);
            }
            for (state, stress) in solved.quadrature_states.iter()
                .zip(&solved.quadrature_stress_pa)
            {
                for i in 0..6 {
                    close(state.strain[i], free[i], 0.01);
                    close(stress[i], 0.0, 1e7);
                }
                close(hot.state(tensor::to_mandel(&state.strain), 1.0)
                    .unwrap().free_energy_j, 0.0, 1e-6);
            }
        }
    });
}

#[test]
fn g1_nonuniform_heating_and_dead_load_use_each_original_source_stiffness() {
    let potential = potential(false, true);
    let temperatures = [310.0, 350.0];
    let sources = temperatures.map(|t| potential.at_temperature(t).unwrap());
    let laws: Vec<_> = sources.iter().map(|s| s.try_mechanical_law().unwrap()).collect();
    let mesh = Mesh2::quads(2.0, 1.0, 2, 1);
    let mut supports = Vec::new();
    let mut forces = Vec::new();
    for (node, p) in mesh.nodes.iter().enumerate() {
        supports.push((2 * node + 1, 0.0));
        if p[0] == 0.0 { supports.push((2 * node, 0.0)); }
        if p[0] == 2.0 { forces.push((2 * node, 1e5)); }
    }
    let strain = temperatures.map(|t| {
        let dt = t - 300.0;
        1e-5 * dt + 0.5e-8 * dt * dt + 1e6 / (12e9 - 3e7 * dt)
    });
    with_cx(|cx| {
        let problem = SmallStrainPlaneProblem {
            mesh: &mesh,
            laws: SmallStrainLawField::PerElement(&laws),
            plane: PlaneKind::Stress,
            formulation: Formulation::Standard,
            thickness_m: 0.2,
            reduction: ReductionConfig {
                max_iterations: 8, stress_tolerance_pa: 1e-6, strain_limit: 0.02,
            },
            rotation_limit_rad: 0.02,
            prescribed_m: &supports,
            nodal_forces_n: &forces,
            load_path: &[1.0],
            budget: PlaneAssemblyBudget::default(),
            settings: IncrementSettings::default(),
        };
        let solved = problem.solve(cx).unwrap();
        for (&p, u) in mesh.nodes.iter().zip(&solved.displacement_m) {
            close(u[0], strain[0] * p[0].min(1.0) + strain[1] * (p[0] - 1.0).max(0.0), 0.01);
            close(u[1], 0.0, 0.01);
        }
        for stress in &solved.quadrature_stress_pa {
            close(stress[0], 1e6, 1e6);
            for &s in &stress[1..] { close(s, 0.0, 1e6); }
        }
    });
}
