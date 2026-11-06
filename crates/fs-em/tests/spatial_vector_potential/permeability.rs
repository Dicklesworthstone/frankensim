//! G1/G3: full tensor curl-curl, magnetic energy and retained induction.
//! All cards are manufactured equation oracles, not measured magnetic materials.
use super::*;
use fs_em::spatial::vector_potential::induction::{
    SpatialInduction, SpatialInductionControls, SpatialInductionLimits,
};
use fs_material::{
    conductor::{
        ELECTRICAL_CONDUCTIVITY_DIMS, ELECTRICAL_CONDUCTIVITY_PROPERTY, ResolvedConductorState,
    },
    field_storage::spatial::{PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES, SpatialMagneticLaw},
    state_point::ResolvedMaterialStatePoint,
};
use fs_qty::{QuantitySpec, Temperature};

const AXES: [[f64; 3]; 3] = [[0.6, 0.8, 0.0], [-0.8, 0.6, 0.0], [0.0, 0.0, 1.0]];
const PRINCIPAL: [f64; 3] = [2.0, 5.0, 11.0];

fn law() -> SpatialMagneticLaw {
    SpatialMagneticLaw::orthotropic(AXES).unwrap()
}

fn card(thermal: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value) in PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES
        .into_iter()
        .zip(PRINCIPAL)
    {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, MAGNETIC_PERMEABILITY_DIMS),
                value: if thermal {
                    PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: vec![(300.0, value), (400.0, 2.0 * value)],
                        dims: MAGNETIC_PERMEABILITY_DIMS,
                    }
                } else {
                    PropertyValue::Scalar {
                        value,
                        dims: MAGNETIC_PERMEABILITY_DIMS,
                    }
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    400.0,
                ),
                interpolation: if thermal {
                    InterpolationPolicy::LinearInside
                } else {
                    InterpolationPolicy::ConstantWithinValidity
                },
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured magnetic tensor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                ELECTRICAL_CONDUCTIVITY_PROPERTY,
                ELECTRICAL_CONDUCTIVITY_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: 5.0,
                dims: ELECTRICAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured conductor".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "fixed tensor axes".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn at(card: &MaterialCard, temperature: f64) -> ResolvedMaterialStatePoint {
    let mut requirements: Vec<_> = PRINCIPAL_MAGNETIC_PERMEABILITY_PROPERTIES
        .into_iter()
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                name,
                MAGNETIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()
        })
        .collect();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            ELECTRICAL_CONDUCTIVITY_PROPERTY,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    resolve_material_state_point(
        card,
        &QueryPoint::new().with("T", temperature).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn model(n: usize, material: &ResolvedMaterialStatePoint) -> SpatialVectorMagnetostatics {
    let (positions, tets) = grid(n, n, n, 1.0 / n as f64, false);
    let count = tets.len();
    SpatialVectorMagnetostatics::with_magnetic_laws(
        positions,
        tets,
        vec![SpatialMagneticMaterial::Linear(material.clone()); count],
        vec![law(); count],
        limits(),
        &mut || false,
    )
    .unwrap()
}

// Independently rotated reciprocal tensor for the authored 3/5,4/5 frame.
// This is NOT elementwise inversion of the rotated permeability.
fn inverse(b: [f64; 3]) -> [f64; 3] {
    [
        0.308 * b[0] + 0.144 * b[1],
        0.144 * b[0] + 0.392 * b[1],
        b[2] / 11.0,
    ]
}

#[test]
fn rotated_reluctivity_recovers_uniform_field_energy_and_gauge_invariance() {
    let source = at(&card(false), 300.0);
    let problem = model(2, &source);
    let b = [0.2, -0.3, 0.4];
    let h = inverse(b);
    for gauge in [false, true] {
        let boundary = trace(&problem, &uniform_trace(&problem, b, gauge));
        let field = problem
            .solve(
                &vec![[0.0; 3]; problem.mesh().tets.len()],
                &boundary,
                controls(),
                &mut || false,
            )
            .unwrap();
        assert_physical(&field);
        for (actual_b, actual_h) in field
            .flux_density_t()
            .iter()
            .zip(field.field_strength_a_m())
        {
            for i in 0..3 {
                near(actual_b[i], b[i]);
                near(actual_h[i], h[i]);
            }
        }
        near(field.stored_energy_j(), 0.5 * dot(b, h));
    }
    assert!(problem.magnetic_laws().iter().all(|l| *l == law()));
    assert!(
        problem
            .permeability()
            .iter()
            .all(|p| p.isotropic_h_m().is_none())
    );
}

fn induction_limits() -> SpatialInductionLimits {
    SpatialInductionLimits {
        maximum_unknowns: 20000,
        maximum_conducting_cells: 4000,
        maximum_assembly_terms: 1_000_000,
    }
}
fn induction_controls() -> SpatialInductionControls {
    SpatialInductionControls {
        relative_tolerance: 2e-11,
        absolute_residual_a: 1e-11,
        restart: 64,
        maximum_cycles: 100,
        maximum_krylov_scalars: 2_000_000,
        current_tolerance_a: 2e-8,
        flux_tolerance_wb: 2e-9,
        energy_tolerance_j: 2e-8,
        maximum_cell_magnetic_peclet: 100.0,
    }
}

#[test]
fn current_carrying_tensor_and_native_induction_match_six_tet_scalar_mode() {
    let source = at(&card(false), 300.0);
    let problem = model(1, &source);
    let zero_trace = trace(&problem, &vec![0.0; problem.mesh().edges.len()]);
    let initial = problem
        .solve(&[[4.0; 3]; 6], &zero_trace, controls(), &mut || false)
        .unwrap();
    assert_physical(&initial);
    // The sole interior edge has curls +/-2 e_i on pairs of volume-1/6 tets.
    // Therefore K=4 tr(mu^-1)/3, J-load=2 and M=3 sigma/40, independently
    // of the production Hodge assembly. Tangential traces fix all other edges.
    let k = (4.0 / 3.0) * (1.0 / 2.0 + 1.0 / 5.0 + 1.0 / 11.0);
    let a0 = 2.0 / k;
    let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
    near(initial.potential_wb()[edge], a0);
    near(initial.stored_energy_j(), 0.5 * k * a0 * a0);
    for (&b, &h) in initial
        .flux_density_t()
        .iter()
        .zip(initial.field_strength_a_m())
    {
        let expected = inverse(b);
        for j in 0..3 {
            near(h[j], expected[j]);
        }
    }
    let conductor = ResolvedConductorState::try_from_material(source).unwrap();
    let state = SpatialInduction::new(
        &initial,
        vec![Some(conductor); 6],
        induction_limits(),
        2e-9,
        &mut || false,
    )
    .unwrap();
    let dt = 0.03;
    let mass = 3.0 * 5.0 / 40.0;
    let expected = mass * a0 / (mass + 0.5 * dt * k);
    let candidate = state
        .step_candidate(
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            dt,
            induction_controls(),
            &mut || false,
        )
        .unwrap();
    near(candidate.state().potential_wb()[edge], expected);
    let step = candidate.step();
    near(step.stored_after_j, 0.5 * k * expected * expected);
    near(step.heat_j, mass * ((expected - a0) / dt).powi(2) * dt);
    near(step.impressed_work_j, expected - a0);
    near(step.energy_residual_j, 0.0);
    assert_eq!(step.parameter_work_j, 0.0);
    assert_eq!(state.time_s(), 0.0);
    let held = candidate
        .state()
        .step_candidate(
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            dt,
            induction_controls(),
            &mut || false,
        )
        .unwrap();
    near(
        held.state().potential_wb()[edge],
        expected * (mass - 0.5 * dt * k) / (mass + 0.5 * dt * k),
    );
    assert!(held.step().stored_after_j < step.stored_after_j);
    assert_eq!(held.step().parameter_work_j, 0.0);
}

#[test]
fn temperature_refresh_retains_tensor_axes_and_closes_signed_parameter_work() {
    let card = card(true);
    let source = at(&card, 300.0);
    let problem = model(1, &source);
    let initial = problem
        .solve(
            &[[4.0; 3]; 6],
            &trace(&problem, &vec![0.0; problem.mesh().edges.len()]),
            controls(),
            &mut || false,
        )
        .unwrap();
    let conductor = ResolvedConductorState::try_from_material(source).unwrap();
    let state = SpatialInduction::new(
        &initial,
        vec![Some(conductor); 6],
        induction_limits(),
        2e-9,
        &mut || false,
    )
    .unwrap();
    let dt = 0.03;
    let k0 = (4.0 / 3.0) * (1.0 / 2.0 + 1.0 / 5.0 + 1.0 / 11.0);
    let k1 = k0 / 1.5;
    let mass = 3.0 * 5.0 / 40.0;
    let a0 = 2.0 / k0;
    let a1 = (mass * a0 - 0.5 * dt * k0 * a0 + dt * 2.0) / (mass + 0.5 * dt * k1);
    let refreshed = at(&card, 350.0);
    let next_conductor = ResolvedConductorState::try_from_material(refreshed.clone()).unwrap();
    let candidate = state
        .step_with_materials_candidate(
            &[[4.0; 3]; 6],
            &[[0.0; 3]; 8],
            dt,
            vec![SpatialMagneticMaterial::Linear(refreshed); 6],
            vec![Some(next_conductor); 6],
            limits(),
            induction_controls(),
            &mut || false,
        )
        .unwrap();
    let edge = problem.mesh().edges.binary_search(&[0, 7]).unwrap();
    near(candidate.state().potential_wb()[edge], a1);
    near(candidate.step().stored_after_j, 0.5 * k1 * a1 * a1);
    near(candidate.step().parameter_work_j, 0.5 * a0 * a1 * (k1 - k0));
    near(candidate.step().energy_residual_j, 0.0);
    assert!(candidate.step().parameter_work_j < 0.0);
    assert_eq!(
        candidate.state().problem().magnetic_laws(),
        problem.magnetic_laws()
    );
    assert_eq!(state.time_s(), 0.0);
}

#[test]
fn distributed_windings_retain_tensor_pairing_and_material_scaling() {
    use fs_em::spatial::vector_potential::winding::{SpatialWindingBasis, SpatialWindingConfig};
    let card = card(true);
    let problem = model(2, &at(&card, 300.0));
    // Curl of an interior P1 hat times a constant vector is a closed, exactly
    // normal-continuous current. Two differently directed loops exercise mutual L.
    let origin = problem
        .positions_m()
        .iter()
        .position(|&p| p == [0.5; 3])
        .unwrap() as u32;
    let sub = |a: [f64; 3], b: [f64; 3]| std::array::from_fn(|i| a[i] - b[i]);
    let mut profiles = vec![Vec::new(), Vec::new()];
    for tet in &problem.mesh().tets {
        let gradient = if tet.contains(&origin) {
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != origin).collect();
            let p = problem.positions_m()[face[0] as usize];
            let normal = cross(
                sub(problem.positions_m()[face[1] as usize], p),
                sub(problem.positions_m()[face[2] as usize], p),
            );
            let height = dot(normal, sub(problem.positions_m()[origin as usize], p));
            normal.map(|v| v / height)
        } else {
            [0.0; 3]
        };
        profiles[0].push(cross(gradient, [0.0, 0.0, 2.0]));
        profiles[1].push(cross(gradient, [2.0, 0.0, 2.0]));
    }
    let config = SpatialWindingConfig {
        max_windings: 4,
        max_basis_values: 100_000,
        pairing_tolerance_h: 2e-9,
    };
    let basis = SpatialWindingBasis::new(
        &problem,
        profiles.clone(),
        trace(&problem, &vec![0.0; problem.mesh().edges.len()]),
        vec![0.0; 2],
        config,
        controls(),
        &mut || false,
    )
    .unwrap();
    let currents = [0.7, -0.2];
    let field = basis.field_for_currents(&currents, &mut || false).unwrap();
    let mut energy = 0.0;
    for (cell, (&b, &h)) in field
        .flux_density_t
        .iter()
        .zip(&field.field_strength_a_m)
        .enumerate()
    {
        let expected = inverse(b);
        for i in 0..3 {
            near(h[i], expected[i]);
        }
        energy += 0.5 * problem.cell_volume_m3(cell).unwrap() * dot(b, expected);
    }
    near(field.field_energy_j, energy);
    near(field.energy_residual_j, 0.0);
    near(basis.inductance_h()[1], basis.inductance_h()[2]);
    assert!(basis.maximum_pairing_defect_h() < config.pairing_tolerance_h);
    let changed = model(2, &at(&card, 350.0));
    let changed_basis = SpatialWindingBasis::new(
        &changed,
        profiles,
        trace(&changed, &vec![0.0; changed.mesh().edges.len()]),
        vec![0.0; 2],
        config,
        controls(),
        &mut || false,
    )
    .unwrap();
    for (&old, &new) in basis
        .inductance_h()
        .iter()
        .zip(changed_basis.inductance_h())
    {
        near(new, 1.5 * old);
    }
}
