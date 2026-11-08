//! Physical heat, history and consistent work on source-bound native Hex8 solids.
//! All coefficients are authored fixtures, not measured or calibrated material data.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::HillState;
use fs_material::plastic::hill::source::{
    HILL_COEFFICIENT_PROPERTIES, HILL_HARDENING_MODULUS_PROPERTY,
    HILL_OVERSTRESS_VISCOSITY_PROPERTY, HILL_REFERENCE_YIELD_STRESS_PROPERTY,
    HILL_VOCE_SATURATION_INCREMENT_PROPERTY, HILL_VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::plastic::hill::thermal::{ThermalHillSource, ThermalHillSourceQuery};
use fs_material::plastic::section::thermal::source::ThermalSectionSourceFamily;
use fs_material::state_point::{
    DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
    MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, Temperature};
use fs_solid::plastic3::hex::thermal::{
    ThermalHexError, ThermalHexHistory, ThermalHexProblem, ThermalHexStep, ThermalHexTolerances,
};
use fs_solid::plastic3::hex::{HexAssemblyBudget, SmallStrainHexProblem};
use fs_solid::{Formulation, IncrementSettings, SmallStrainLawField};

const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const CELLS: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const IDENTITY: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const E: f64 = 1000.0;
const NU: f64 = 0.3;
const Y: f64 = 1.0;
const H: f64 = 60.0;
const Q: f64 = 2.0;
const Q_SLOPE: f64 = -0.1;
const B: f64 = 500.0;
const RHO: f64 = 2.0;
const CP: f64 = 500.0;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const DT: f64 = 0.1;

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= allowance,
        "{actual:e} != {expected:e}, allowance {allowance:e}"
    );
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 11,
                kernel_id: 97,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn property(name: &str, dims: Dims, value: f64, slope: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: if slope == 0.0 {
            PropertyValue::Scalar { value, dims }
        } else {
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, value - 20.0 * slope), (320.0, value + 20.0 * slope)],
                dims,
            }
        },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        interpolation: if slope == 0.0 {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        },
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored thermal anisotropic solid".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
}
impl Fixture {
    fn new(expansion: f64) -> Self {
        let mut properties = vec![
            property(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E, 0.0),
            property(POISSON_RATIO_PROPERTY, Dims([0; 6]), NU, 0.0),
            property(HILL_REFERENCE_YIELD_STRESS_PROPERTY, Pressure::DIMS, Y, 0.0),
            property(HILL_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, H, 0.0),
            property(
                HILL_VOCE_SATURATION_INCREMENT_PROPERTY,
                Pressure::DIMS,
                Q,
                Q_SLOPE,
            ),
            property(HILL_VOCE_SATURATION_RATE_PROPERTY, Dims([0; 6]), B, 0.0),
            property("specific_heat_capacity", CP_DIMS, CP, 0.0),
            property(DENSITY_PROPERTY, Density::DIMS, RHO, 0.0),
            property(
                HILL_OVERSTRESS_VISCOSITY_PROPERTY,
                Dims([-1, 1, -1, 0, 0, 0]),
                20.0,
                0.0,
            ),
            property(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                expansion,
                0.0,
            ),
        ];
        properties.extend(
            HILL_COEFFICIENT_PROPERTIES
                .into_iter()
                .zip([0.35, 0.8, 0.45, 1.5, 1.1, 2.0])
                .map(|(name, value)| property(name, Dims([0; 6]), value, 0.0)),
        );
        let mut claims = ClaimSet::new();
        for claim in properties {
            claims.insert_claim(claim).unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "authored anisotropic metal".into(),
                phase: "solid".into(),
                process: "thermal solid regression".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let lower = QueryPoint::new().with("T", 280.0).unwrap();
        let upper = QueryPoint::new().with("T", 320.0).unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(
            &card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("specific_heat_capacity", CP_DIMS),
                density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -20.0 * CP,
            1e-8,
        )
        .unwrap();
        Self { card, caloric }
    }
    fn source(&self) -> ThermalHillSource<'_> {
        ThermalHillSource::try_new(
            ThermalHillSourceQuery {
                card: &self.card,
                caloric: &self.caloric,
                context: QueryPoint::new().with("T", 300.0).unwrap(),
                temperature_axis: "T".into(),
                selection: MaterialPropertySelection::SingleClaimOnly,
                thermal_expansion: Some(MaterialPropertySelection::SingleClaimOnly),
                strain_limit: 0.04,
                material_to_spatial: IDENTITY,
            },
            0.0,
        )
        .unwrap()
        .with_voce_isotropic_hardening(MaterialPropertySelection::SingleClaimOnly)
        .unwrap()
    }
}
fn h(temperature: f64) -> f64 {
    CP * (temperature - 300.0)
}
fn tolerances() -> ThermalHexTolerances {
    ThermalHexTolerances {
        material: EnergyTolerances {
            energy_j_m3: 1e-9,
            yield_pa: 1e-7,
        },
        energy_j: 1e-8,
    }
}
fn problem<'p, F>(
    nodes: &'p [[f64; 3]],
    families: &'p [F],
    formulation: Formulation,
    fixed: &'p [(usize, f64)],
    forces: &'p [(usize, f64)],
) -> ThermalHexProblem<'p, F> {
    ThermalHexProblem {
        nodes_m: nodes,
        hexahedra: &CELLS,
        families,
        formulation,
        strain_limit: 0.04,
        rotation_limit_rad: 0.04,
        prescribed_m: fixed,
        nodal_forces_n: forces,
        budget: HexAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}
fn symmetry() -> Vec<(usize, f64)> {
    (0..24)
        .filter(|d| NODES[d / 3][d % 3] == 0.0)
        .map(|d| (d, 0.0))
        .collect()
}
fn history<'h>(
    u: &'h [f64],
    states: &'h [HillState],
    enthalpy: &'h [f64],
) -> ThermalHexHistory<'h, HillState> {
    ThermalHexHistory {
        displacement_m: u,
        quadrature_states: states,
        specific_enthalpy_j_kg: enthalpy,
    }
}
fn closed(step: &ThermalHexStep<HillState>) {
    near(
        step.energy.equilibrium_work_residual_j,
        0.0,
        tolerances().energy_j,
    );
    near(step.energy.balance_residual_j, 0.0, tolerances().energy_j);
    for point in &step.quadrature_energy {
        assert!(point.physical_heat_j >= 0.0 && point.numerical_loss_j >= 0.0);
        near(point.creep_heat_j, 0.0, 1e-14);
    }
}

#[test]
fn free_and_restrained_heating_close_separate_storage_entropy_and_numerical_work() {
    with_cx(|cx| {
        let expansion = 1e-5;
        let fixture = Fixture::new(expansion);
        let families = [fixture.source()];
        let (u, states, old_h) = ([0.0; 24], vec![HillState::default(); 8], [0.0; 8]);
        let new_h = [h(320.0); 8];
        let thermal_strain = 20.0 * expansion;
        let bulk = E / (3.0 * (1.0 - 2.0 * NU));
        for formulation in [Formulation::Standard, Formulation::BBar] {
            for restrained in [false, true] {
                let fixed = if restrained {
                    (0..24).map(|d| (d, 0.0)).collect()
                } else {
                    symmetry()
                };
                let step = problem(&NODES, &families, formulation, &fixed, &[])
                    .solve_from(history(&u, &states, &old_h), &new_h, DT, tolerances(), cx)
                    .unwrap();
                closed(&step);
                near(
                    step.energy.applied_work_j + step.energy.support_work_j,
                    0.0,
                    1e-14,
                );
                near(step.energy.material.physical_heat_j, 0.0, 1e-14);
                for (node, displacement) in step.solution.displacement_m.iter().enumerate() {
                    for j in 0..3 {
                        near(
                            displacement[j],
                            if restrained {
                                0.0
                            } else {
                                thermal_strain * NODES[node][j]
                            },
                            2e-11,
                        );
                    }
                }
                for (state, stress) in step
                    .solution
                    .quadrature_states
                    .iter()
                    .zip(&step.solution.quadrature_stress_pa)
                {
                    assert_eq!(*state, HillState::default());
                    for (j, value) in stress.iter().enumerate() {
                        near(
                            *value,
                            if restrained && j < 3 {
                                -3.0 * bulk * thermal_strain
                            } else {
                                0.0
                            },
                            2e-8,
                        );
                    }
                }
                // psi=9K(a*dT)^2/2, entropy=a*tr(sigma).
                // Endpoint quadrature has a positive numerical elastic loss;
                // the reversible caloric transfer is -T_old*delta(entropy).
                let stored = if restrained {
                    4.5 * bulk * thermal_strain.powi(2)
                } else {
                    0.0
                };
                let entropy = if restrained {
                    -9.0 * bulk * expansion * thermal_strain
                } else {
                    0.0
                };
                near(step.energy.material.stored_change_j, stored, 2e-10);
                near(
                    step.energy.material.thermoelastic_storage_change_j,
                    320.0 * entropy,
                    2e-10,
                );
                near(
                    step.energy.material.reversible_heat_j,
                    -300.0 * entropy,
                    2e-10,
                );
                near(step.energy.material.numerical_loss_j, stored, 2e-10);
                near(
                    step.energy.material.coefficient_storage_change_j,
                    0.0,
                    1e-14,
                );
                assert_eq!(step.specific_enthalpy_j_kg, new_h);
                assert_eq!(step.quadrature_temperature_k, [320.0; 8]);
            }
        }
    });
}

// Independent uniaxial associated Hill flow: q=sqrt(G+H_hill)*sigma,
// ep=alpha*[sqrt(G+H_hill), -H_hill/sqrt(G+H_hill), -G/sqrt(G+H_hill)].
fn uniaxial_alpha(sigma: f64, temperature: f64) -> f64 {
    let q = 1.25_f64.sqrt() * sigma;
    let saturation = Q + Q_SLOPE * (temperature - 300.0);
    let (mut lo, mut hi) = (0.0, (q - Y) / H);
    for _ in 0..80 {
        let mid = 0.5 * lo + 0.5 * hi;
        if Y + H * mid + saturation * -(-B * mid).exp_m1() < q {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * lo + 0.5 * hi
}
fn plastic_uniaxial(alpha: f64) -> [f64; 3] {
    let r = 1.25_f64.sqrt();
    [r * alpha, -0.45 / r * alpha, -0.8 / r * alpha]
}
fn voce_shape(alpha: f64) -> f64 {
    alpha + (-B * alpha).exp_m1() / B
}

#[test]
fn loaded_thermal_softening_preserves_original_history_then_unloads_to_permanent_set() {
    with_cx(|cx| {
        let fixture = Fixture::new(0.0);
        let families = [fixture.source()];
        let fixed = symmetry();
        let sigma = 2.0;
        let forces: Vec<_> = (0..8)
            .filter(|n| NODES[*n][0] == 1.0)
            .map(|n| (3 * n, sigma / 4.0))
            .collect();
        let states = vec![HillState::default(); 8];
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let p = problem(&NODES, &families, formulation, &fixed, &forces);
            let cold = p
                .solve_from(
                    history(&[0.0; 24], &states, &[0.0; 8]),
                    &[0.0; 8],
                    DT,
                    tolerances(),
                    cx,
                )
                .unwrap();
            let saved = cold.solution.quadrature_states.clone();
            let cold_alpha = uniaxial_alpha(sigma, 300.0);
            let hot_alpha = uniaxial_alpha(sigma, 320.0);
            assert!(hot_alpha > 8.0 * cold_alpha);
            let hot = p
                .solve_from(cold.history(), &[h(320.0); 8], DT, tolerances(), cx)
                .unwrap();
            assert_eq!(cold.solution.quadrature_states, saved);
            for (step, alpha) in [(&cold, cold_alpha), (&hot, hot_alpha)] {
                closed(step);
                let plastic = plastic_uniaxial(alpha);
                for (node, displacement) in step.solution.displacement_m.iter().enumerate() {
                    for j in 0..3 {
                        let elastic = if j == 0 { sigma / E } else { -NU * sigma / E };
                        near(
                            displacement[j],
                            (plastic[j] + elastic) * NODES[node][j],
                            3e-10,
                        );
                    }
                }
                for (state, stress) in step
                    .solution
                    .quadrature_states
                    .iter()
                    .zip(&step.solution.quadrature_stress_pa)
                {
                    near(state.alpha, alpha, 3e-10);
                    for j in 0..3 {
                        near(state.plastic_strain[j], plastic[j], 3e-10);
                    }
                    for j in 0..6 {
                        near(stress[j], if j == 0 { sigma } else { 0.0 }, 3e-8);
                    }
                }
            }
            near(
                hot.energy.material.physical_heat_j,
                Y * (hot_alpha - cold_alpha),
                4e-10,
            );
            near(
                hot.energy.material.coefficient_storage_change_j,
                -Q * voce_shape(cold_alpha),
                4e-10,
            );
            let delta_ts =
                -Q_SLOPE * (320.0 * voce_shape(hot_alpha) - 300.0 * voce_shape(cold_alpha));
            near(
                hot.energy.material.thermoelastic_storage_change_j,
                delta_ts,
                4e-9,
            );
            near(
                hot.energy.material.reversible_heat_j,
                Q * voce_shape(cold_alpha) - delta_ts,
                4e-9,
            );
            near(
                hot.energy.material.numerical_loss_j,
                0.5 * H * (hot_alpha - cold_alpha).powi(2),
                4e-10,
            );
            let unloaded = problem(&NODES, &families, formulation, &fixed, &[])
                .solve_from(hot.history(), &[h(320.0); 8], DT, tolerances(), cx)
                .unwrap();
            closed(&unloaded);
            near(unloaded.energy.material.physical_heat_j, 0.0, 1e-12);
            let permanent = plastic_uniaxial(hot_alpha);
            for (node, displacement) in unloaded.solution.displacement_m.iter().enumerate() {
                for j in 0..3 {
                    near(displacement[j], permanent[j] * NODES[node][j], 3e-10);
                }
            }
            for (state, stress) in unloaded
                .solution
                .quadrature_states
                .iter()
                .zip(&unloaded.solution.quadrature_stress_pa)
            {
                near(state.alpha, hot_alpha, 3e-10);
                for value in stress {
                    near(*value, 0.0, 3e-8);
                }
            }
        }
    });
}

#[test]
fn dynamic_heating_retains_reference_mass_and_support_work_without_turning_kinetic_loss_into_heat()
{
    with_cx(|cx| {
        let fixture = Fixture::new(1e-5);
        let families = [fixture.source()];
        let strain = 20.0e-5;
        let fixed: Vec<_> = NODES
            .iter()
            .flat_map(|x| x.map(|v| strain * v))
            .enumerate()
            .collect();
        let states = vec![HillState::default(); 8];
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let step = problem(&NODES, &families, formulation, &fixed, &[])
                .advance_backward_euler(
                    history(&[0.0; 24], &states, &[0.0; 8]),
                    &[[0.0; 3]; 8],
                    &[h(320.0); 8],
                    DT,
                    tolerances(),
                    cx,
                )
                .unwrap();
            closed(&step);
            let inertia = step.inertia.unwrap();
            let kinetic = 0.5 * RHO * (strain / DT).powi(2);
            near(inertia.reference_mass_kg, RHO, 2e-13);
            near(inertia.kinetic_change_j, kinetic, 2e-13);
            near(inertia.numerical_loss_j, kinetic, 2e-13);
            near(inertia.inertial_work_j, 2.0 * kinetic, 2e-13);
            near(step.energy.support_work_j, 2.0 * kinetic, 2e-13);
            near(step.energy.material.physical_heat_j, 0.0, 1e-14);
            near(step.energy.material.numerical_loss_j, 0.0, 1e-14);
            // Consistent mass couples even the stationary x=0/y=0/z=0 supports.
            for (&(dof, _), &reaction) in fixed.iter().zip(&step.solution.increments[0].reactions_n)
            {
                let moment = if NODES[dof / 3][dof % 3] == 0.0 {
                    1.0 / 24.0
                } else {
                    1.0 / 12.0
                };
                near(reaction, RHO * strain / (DT * DT) * moment, 2e-11);
            }
            for (node, velocity) in step.velocity_mps.as_ref().unwrap().iter().enumerate() {
                for j in 0..3 {
                    near(velocity[j], strain * NODES[node][j] / DT, 2e-13);
                }
            }
        }
    });
}

#[test]
fn local_thermal_laws_keep_nonuniform_bbar_histories_and_the_frozen_history_tangent() {
    with_cx(|cx| {
        let fixture = Fixture::new(0.0);
        let families = [fixture.source()];
        let nodes = NODES.map(|[x, y, z]| [x + 0.18 * x * y, y + 0.10 * y * z, z + 0.12 * x * z]);
        let u: Vec<_> = nodes
            .iter()
            .flat_map(|x| {
                [
                    0.006 * x[0] + 0.002 * x[1] * x[2],
                    -0.002 * x[1] + 0.001 * x[0] * x[2],
                    -0.001 * x[2] + 0.0015 * x[0] * x[1],
                ]
            })
            .collect();
        let fixed: Vec<_> = u.iter().copied().enumerate().collect();
        let states = vec![HillState::default(); 8];
        let enthalpy = core::array::from_fn::<_, 8, _>(|p| h(288.0 + 4.0 * p as f64));
        let p = problem(&nodes, &families, Formulation::BBar, &fixed, &[]);
        let step = p
            .solve_from(
                history(&[0.0; 24], &states, &[0.0; 8]),
                &enthalpy,
                DT,
                tolerances(),
                cx,
            )
            .unwrap();
        closed(&step);
        for i in 0..8 {
            near(
                step.quadrature_temperature_k[i],
                288.0 + 4.0 * i as f64,
                1e-12,
            );
            let independent = families[0]
                .section_increment(
                    0.0,
                    enthalpy[i],
                    &[0.0; 6],
                    &states[i],
                    &step.solution.quadrature_strain[i],
                    tolerances().material,
                    DT,
                )
                .unwrap();
            assert_eq!(step.solution.quadrature_states[i], independent.state);
            near(
                step.quadrature_energy[i].physical_heat_j,
                Y * independent.state.alpha * step.solution.quadrature_volume_m3[i],
                1e-12,
            );
        }
        let alpha0 = step.solution.quadrature_states[0].alpha;
        assert!(
            step.solution
                .quadrature_states
                .iter()
                .all(|s| s.alpha > 0.0)
        );
        assert!(
            step.solution
                .quadrature_states
                .iter()
                .any(|s| (s.alpha - alpha0).abs() > 1e-4)
        );
        let trace0: f64 = step.solution.quadrature_strain[0][..3].iter().sum();
        for strain in &step.solution.quadrature_strain {
            near(strain[..3].iter().sum(), trace0, 2e-14);
        }
        // Resolve the same original point laws, then perturb displacement with
        // all OLD histories frozen. This independently checks B-bar columns and
        // the local thermal Hill-Voce consistent tangent at the actual endpoint.
        let laws: Vec<_> = enthalpy
            .iter()
            .map(|h| families[0].resolve_section_step(0.0, *h, DT).unwrap().law)
            .collect();
        let fields = [laws.as_slice()];
        let native = SmallStrainHexProblem {
            nodes_m: &nodes,
            hexahedra: &CELLS,
            laws: SmallStrainLawField::PerQuadrature(&fields),
            formulation: Formulation::BBar,
            strain_limit: 0.04,
            rotation_limit_rad: 0.04,
            prescribed_m: &fixed,
            nodal_forces_n: &[],
            load_path: &[1.0],
            budget: HexAssemblyBudget::default(),
            settings: IncrementSettings::default(),
        };
        let tangent = native.tangent_matrix(&u, &states, cx).unwrap();
        let (mut error, mut scale) = (0.0_f64, 0.0_f64);
        let delta = 1e-8;
        for j in [0, 4, 8, 12, 16, 20] {
            let (mut plus, mut minus) = (u.clone(), u.clone());
            plus[j] += delta;
            minus[j] -= delta;
            let fp = native.internal_force(&plus, &states, cx).unwrap();
            let fm = native.internal_force(&minus, &states, cx).unwrap();
            for i in 0..24 {
                let derivative = (fp[i] - fm[i]) / (2.0 * delta);
                error = error.hypot(tangent.get(i, j) - derivative);
                scale = scale.hypot(derivative);
            }
        }
        assert!(
            error / scale < 3e-6,
            "thermal B-bar tangent error {}",
            error / scale
        );
        // A late point-domain refusal cannot advance any previously admitted
        // point. The full accepted caloric and mechanical carrier is reusable.
        let saved = step.solution.clone();
        let mut outside = enthalpy;
        outside[7] = h(321.0);
        assert!(matches!(
            p.solve_from(step.history(), &outside, DT, tolerances(), cx),
            Err(ThermalHexError::Material { point: 7, .. })
        ));
        assert_eq!(step.solution, saved);
        assert_eq!(step.specific_enthalpy_j_kg, enthalpy);
    });
}

#[test]
fn source_bound_rate_load_and_hold_consume_each_physical_interval_once() {
    with_cx(|cx| {
        let fixture = Fixture::new(0.0);
        let families = [fixture
            .source()
            .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            .unwrap()];
        let gamma = 0.008;
        let fixed: Vec<_> = NODES
            .iter()
            .flat_map(|x| [gamma * x[1], 0.0, 0.0])
            .enumerate()
            .collect();
        let mu = E / (2.0 * (1.0 + NU));
        // Pure xy shear has q=sqrt(2N)*tau with N=2, and gamma_p=2*alpha.
        let r = 2.0;
        let eta = 20.0;
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let p = problem(&NODES, &families, formulation, &fixed, &[]);
            let (mut old_u, mut old_states, mut old_h) =
                (vec![0.0; 24], vec![HillState::default(); 8], vec![0.0; 8]);
            let mut old_alpha = 0.0;
            for (interval, (dt, temperature)) in [(0.1, 300.0), (0.4, 300.0), (0.05, 310.0)]
                .into_iter()
                .enumerate()
            {
                let saturation = Q + Q_SLOPE * (temperature - 300.0);
                let (mut lo, mut hi) = (
                    old_alpha,
                    (r * mu * gamma - Y + eta / dt * old_alpha) / (r * r * mu + H + eta / dt),
                );
                for _ in 0..80 {
                    let mid = 0.5 * lo + 0.5 * hi;
                    let residual = r * mu * gamma
                        - r * r * mu * mid
                        - Y
                        - H * mid
                        - saturation * -(-B * mid).exp_m1()
                        - eta / dt * (mid - old_alpha);
                    if residual > 0.0 {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let alpha = 0.5 * lo + 0.5 * hi;
                let delta_alpha = alpha - old_alpha;
                let tau = mu * (gamma - r * alpha);
                assert!(delta_alpha > 0.0);
                let step = p
                    .solve_from(
                        history(&old_u, &old_states, &old_h),
                        &[h(temperature); 8],
                        dt,
                        tolerances(),
                        cx,
                    )
                    .unwrap();
                closed(&step);
                for (state, stress) in step
                    .solution
                    .quadrature_states
                    .iter()
                    .zip(&step.solution.quadrature_stress_pa)
                {
                    near(state.alpha, alpha, 3e-11);
                    near(state.plastic_strain[3], alpha, 3e-11);
                    for j in 0..6 {
                        near(stress[j], if j == 3 { tau } else { 0.0 }, 3e-8);
                    }
                }
                near(step.duration_s, dt, 0.0);
                near(
                    step.energy.material.physical_heat_j,
                    Y * delta_alpha + eta / dt * delta_alpha * delta_alpha,
                    4e-10,
                );
                near(
                    step.energy.material.mechanical_work_j,
                    if interval == 0 { tau * gamma } else { 0.0 },
                    4e-10,
                );
                // The accepted stress must be from OLD history. Advancing
                // this timed law again at committed history relaxes a second
                // physical interval and cannot reproduce that endpoint.
                let endpoint = families[0]
                    .resolve_section_step(old_h[0], h(temperature), dt)
                    .unwrap();
                let twice = endpoint.law.stress(
                    &step.solution.quadrature_strain[0],
                    &step.solution.quadrature_states[0],
                );
                assert!((twice[3] - tau).abs() > 1e-4);
                old_alpha = alpha;
                old_u = step.solution.displacement_m.as_flattened().to_vec();
                old_states = step.solution.quadrature_states;
                old_h = step.specific_enthalpy_j_kg;
            }
        }
    });
}
