//! G1: measured-role K/G sources enter the original reversible thermal trajectory.
//! Numerical fixtures only; no measured-material qualification is claimed.
use super::*;
use fs_solid::linear3::source::ThermalStiffnessModel;
use fs_solid::loads3::transient::TransientConfig3;

const BULK: f64 = 1000.0;
const BULK_SLOPE: f64 = 2.0;
const SHEAR: f64 = 500.0;
const SHEAR_SLOPE: f64 = -1.0;
const ALPHA: f64 = 0.0003;
const ALPHA_SLOPE: f64 = 0.000001;
const CAPACITY: f64 = 2.0;

fn fixture() -> (MaterialCard, EquilibriumEnthalpyPhaseCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        ("bulk_modulus", Pressure::DIMS, BULK, BULK_SLOPE),
        ("shear_modulus", Pressure::DIMS, SHEAR, SHEAR_SLOPE),
        (
            "linear_thermal_expansion_coefficient",
            INVERSE_TEMPERATURE_DIMS,
            ALPHA,
            ALPHA_SLOPE,
        ),
        ("density", Density::DIMS, 6.0, 0.0),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            CAPACITY,
            0.0,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, value - 20.0 * slope), (360.0, value + 60.0 * slope)],
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    quantity(),
                    280.0,
                    360.0,
                ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "authored independent K/G solid".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "bulk/shear thermal fixture".into(),
            phase: "solid".into(),
            process: "numerical source".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        vec![
            EnthalpyPhaseKnot {
                temperature_k: 280.0,
                specific_enthalpy_j_kg: 0.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 6.0,
            },
            EnthalpyPhaseKnot {
                temperature_k: 360.0,
                specific_enthalpy_j_kg: 80.0 * CAPACITY,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 6.0,
            },
        ],
    )
    .unwrap();
    (card, curve)
}

// Independent isotropic invariants with only epsilon_xx=u. The thermal
// expansion changes the trace; the deviatoric norm is exactly 2u^2/3.
fn response(t: f64, u: f64) -> (f64, f64) {
    let dt = t - 300.0;
    let k = BULK + BULK_SLOPE * dt;
    let g = SHEAR + SHEAR_SLOPE * dt;
    let free = ALPHA * dt + 0.5 * ALPHA_SLOPE * dt * dt;
    let alpha = ALPHA + ALPHA_SLOPE * dt;
    let trace = u - 3.0 * free;
    let dev2 = 2.0 * u * u / 3.0;
    let psi = VOLUME * (0.5 * k * trace * trace + g * dev2);
    let entropy =
        VOLUME * (3.0 * k * alpha * trace - 0.5 * BULK_SLOPE * trace * trace - SHEAR_SLOPE * dev2);
    ((k + 4.0 * g / 3.0) * u - 3.0 * k * free, psi + t * entropy)
}

#[test]
fn g1_bulk_shear_cards_heat_freely_and_under_restraint_without_synthetic_young_poisson_claims() {
    with_cx(|cx| {
        for clamped in [false, true] {
            let (card, curve) = fixture();
            let heat = thermal(&[&curve], &[CAPACITY], &[20.0], 0.0);
            let context = context(300.0);
            let mut binding = source(30, &card, &curve, &context);
            binding.stiffness_model = ThermalStiffnessModel::IsotropicBulkShear;
            let source = SolidEnthalpy::try_new(
                &heat,
                &NODES,
                &TETS,
                vec![binding],
                TetAssemblyBudget::standard(),
                0.05,
                mechanical().mesh_frame,
                cx,
            )
            .unwrap();
            let loads = source
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let fixed = if clamped {
                (0..12).collect::<Vec<_>>()
            } else {
                FREE.to_vec()
            };
            let body = ThermoelasticSolid::try_new(
                &source,
                &heat,
                &loads,
                &fixed,
                &[[280.0, 360.0]],
                mechanical(),
                solve(),
                cx,
            )
            .unwrap();
            let next = body
                .advance(&heat, &(), 0.1, &loads, thermal_controls(), solve(), cx)
                .unwrap();
            let expected = if clamped {
                root(280.0, 359.0, |t| {
                    CAPACITY * (t - 300.0) + response(t, 0.0).1 - 2.0
                })
            } else {
                301.0
            };
            near(next.thermal.temperature_k(IDS[0]).unwrap(), expected, 3e-8);
            let dt = expected - 300.0;
            let free = ALPHA * dt + 0.5 * ALPHA_SLOPE * dt * dt;
            let actual = &next.solid.equilibrium().solution;
            for (node, x) in NODES.iter().enumerate() {
                for axis in 0..3 {
                    near(
                        actual.displacement.displacement_m()[node][axis],
                        if clamped { 0.0 } else { free * x[axis] },
                        2e-10,
                    );
                }
            }
            let expected_stress = if clamped {
                response(expected, 0.0).0
            } else {
                0.0
            };
            for &stress in &actual.stress.elements()[0].stress_mandel_pa()[..3] {
                near(stress, expected_stress, 2e-8);
            }
            let sample = &next.solid.state().samples()[0];
            assert!(sample.elastic.property("young_modulus").is_none());
            assert!(sample.elastic.property("poisson_ratio").is_none());
            assert_eq!(sample.elastic.card_identity(), card.content_hash());
            near(sample.reference_mass_kg, 1.0, 1e-14);
            near(
                next.thermal.energy_j() - heat.energy_j() + next.report.internal_energy_change_j,
                2.0,
                1e-9,
            );
            near(next.report.combined_residual_j, 0.0, 1e-9);
            if clamped {
                assert!(next.report.reversible_heat_j > 0.01);
            } else {
                near(next.report.reversible_heat_j, 0.0, 1e-10);
            }
            assert_eq!(body.time_s(), 0.0);
            assert_eq!(heat.time_s(), 0.0);
        }
    });
}

#[test]
fn g1_bulk_shear_heat_feedback_preserves_native_inertia_and_actual_work() {
    with_cx(|cx| {
        let (card, curve) = fixture();
        let heat = thermal(&[&curve], &[CAPACITY], &[200.0], 0.0);
        let context = context(300.0);
        let mut binding = source(30, &card, &curve, &context);
        binding.stiffness_model = ThermalStiffnessModel::IsotropicBulkShear;
        let source = SolidEnthalpy::try_new(
            &heat,
            &NODES,
            &TETS,
            vec![binding],
            TetAssemblyBudget::standard(),
            0.05,
            mechanical().mesh_frame,
            cx,
        )
        .unwrap();
        let loads = source
            .geometry()
            .assemble(&[], VolumeLoad3::None, &[(3, 0.6)], cx)
            .unwrap();
        let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).collect();
        let equilibrium = ThermoelasticSolid::try_new(
            &source,
            &heat,
            &loads,
            &fixed,
            &[[280.0, 360.0]],
            mechanical(),
            solve(),
            cx,
        )
        .unwrap();
        let u0 = 0.6 / (VOLUME * (BULK + 4.0 * SHEAR / 3.0));
        let v0 = 0.01;
        near(
            equilibrium
                .equilibrium()
                .solution
                .displacement
                .displacement_m()[1][0],
            u0,
            1e-12,
        );
        let velocity = [[0.0; 3], [v0, 0.0, 0.0], [0.0; 3], [0.0; 3]];
        let body = equilibrium
            .begin_transient(
                &loads,
                &velocity,
                TransientConfig3 {
                    maximum_steps: 10,
                    maximum_displacement_m: 0.05,
                    absolute_force_tolerance_n: 1e-11,
                    absolute_energy_tolerance_j: 1e-10,
                    relative_energy_tolerance: 1e-8,
                    ..TransientConfig3::default()
                },
                cx,
            )
            .unwrap();
        let dt = 0.01;
        let mass = 0.1;
        let evaluate = |t: f64| {
            let delta = t - 300.0;
            let k = BULK + BULK_SLOPE * delta;
            let g = SHEAR + SHEAR_SLOPE * delta;
            let free = ALPHA * delta + 0.5 * ALPHA_SLOPE * delta * delta;
            let k0 = VOLUME * (BULK + 4.0 * SHEAR / 3.0);
            let k1 = VOLUME * (k + 4.0 * g / 3.0);
            let thermal_force = VOLUME * 3.0 * k * free;
            let u = ((2.0 * mass / (dt * dt) - 0.5 * k0) * u0
                + 2.0 * mass / dt * v0
                + 0.6
                + 0.5 * thermal_force)
                / (2.0 * mass / (dt * dt) + 0.5 * k1);
            let v = 2.0 * (u - u0) / dt - v0;
            let old = response(300.0, u0);
            let new = response(t, u);
            let work = 0.5 * VOLUME * (old.0 + new.0) * (u - u0);
            let reversible = work - (new.1 - old.1);
            (
                CAPACITY * (t - 300.0) - 200.0 * dt - reversible,
                [u, v, reversible],
            )
        };
        let expected_t = root(280.0, 359.0, |t| evaluate(t).0);
        let [expected_u, expected_v, expected_heat] = evaluate(expected_t).1;
        let next = body
            .advance(&heat, &(), dt, &loads, thermal_controls(), solve(), cx)
            .unwrap();
        let motion = next.solid.nodal_motion(cx).unwrap();
        near(
            next.thermal.temperature_k(IDS[0]).unwrap(),
            expected_t,
            3e-8,
        );
        near(motion[0][1][0], expected_u, 1e-10);
        near(motion[1][1][0], expected_v, 1e-9);
        near(next.report.reversible_heat_j, expected_heat, 2e-9);
        near(
            next.report.mechanical_work_j,
            0.6 * (expected_u - u0),
            1e-10,
        );
        near(
            next.report.kinetic_energy_change_j,
            0.5 * mass * (expected_v * expected_v - v0 * v0),
            1e-10,
        );
        near(next.report.external_heat_j, 2.0, 1e-9);
        near(next.report.combined_residual_j, 0.0, 1e-9);
        near(
            next.solid.state().samples()[0].reference_mass_kg,
            1.0,
            1e-14,
        );
        assert_eq!(next.solid.step(), 1);
        assert_eq!(body.step(), 0);
        assert_eq!(heat.time_s(), 0.0);
    });
}
