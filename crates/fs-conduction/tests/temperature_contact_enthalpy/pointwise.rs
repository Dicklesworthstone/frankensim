//! Independent simplex moments for one coupled contact/radiation enthalpy solve.
#[path = "pointwise/surface.rs"]
mod surface;
#[path = "pointwise/nodal_heat.rs"]
mod nodal_heat;
use super::*;
use fs_conduction::radiation::{GrayDiffuseEnclosure, RadiationSurface, ViewFactorMatrix};
use fs_conduction::{
    AmbientRadiationPatch, ConductionError, EMISSIVITY_DIMS, PointwiseAmbientRadiation,
    PointwiseEnclosureRadiation, PointwiseEnclosureTolerance, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, ThermalBc, ThermalBoundary,
};

pub(super) const DT: f64 = 0.2;
pub(super) const TARGET_H: [f64; 8] = [120., 250., 90., 50., 400., 550., 700., 450.];
const TARGET_T: [f64; 8] = [310., 315., 309., 305., 340., 350., 360., 345.];
const COLD_FACE: [usize; 3] = [0, 1, 3];
const HOT_FACE: [usize; 3] = [4, 5, 7];

pub(super) fn emissivity(region: &str, value: f64) -> SurfaceEmissivity {
    emissivity_in(region, value, 280., 380.)
}

fn emissivity_in(region: &str, value: f64, lower: f64, upper: f64) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", lower, upper),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic coupled contact/radiation simplex reference".into(),
                license: "internal-test-use".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = fs_matdb::MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic gray surface".into(),
            phase: "solid".into(),
            process: "numerical reference".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    SurfaceEmissivity::from_card(
        region,
        &card,
        (lower + upper) / 2.,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap()
}

pub(super) fn ambient_boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "exterior",
            |face| face.vertices.iter().all(|v| [0, 1, 3].contains(v)),
            ThermalBc::robin(2., 300.).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

pub(super) fn enclosure_boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "left",
            |face| face.vertices.iter().all(|v| [0, 1, 3].contains(v)),
            ThermalBc::adiabatic(),
        )
        .unwrap()
        .region(
            "right",
            |face| face.vertices.iter().all(|v| [4, 5, 7].contains(v)),
            ThermalBc::adiabatic(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

pub(super) fn enclosure(mesh: &ConductionMesh, epsilon: [f64; 2]) -> GrayDiffuseEnclosure {
    let left = RadiationSurface::new(
        mesh,
        "left",
        |face| face.vertices.iter().all(|v| [0, 1, 3].contains(v)),
        emissivity("left", epsilon[0]),
    )
    .unwrap();
    let right = RadiationSurface::new(
        mesh,
        "right",
        |face| face.vertices.iter().all(|v| [4, 5, 7].contains(v)),
        emissivity("right", epsilon[1]),
    )
    .unwrap();
    // Prescribed two-patch operator for this numerical reference; these view
    // factors are not inferred from the displayed tetrahedral geometry.
    let factors = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).unwrap();
    GrayDiffuseEnclosure::new(vec![left, right], factors).unwrap()
}

// Exact triangle moment int(N_i T_h^4) = A/105 sum_{|alpha|=4}
// (alpha_i+1) t^alpha, independent of the production quadrature rule.
fn emission_moments(t: [f64; 3]) -> [f64; 3] {
    let mut moments = [0.; 3];
    for a in 0..=4 {
        for b in 0..=4 - a {
            let c = 4 - a - b;
            let powers = [a, b, c];
            let monomial = t[0].powi(a) * t[1].powi(b) * t[2].powi(c);
            for i in 0..3 {
                moments[i] += 0.5 / 105. * f64::from(powers[i] + 1) * monomial;
            }
        }
    }
    moments
}

fn ambient_load(ambient: f64) -> [f64; 8] {
    let mut load = [0.; 8];
    let moments = emission_moments(COLD_FACE.map(|v| TARGET_T[v]));
    for i in 0..3 {
        load[COLD_FACE[i]] =
            0.2 * STEFAN_BOLTZMANN_W_M2_K4 * (moments[i] - 0.5 / 3. * ambient.powi(4));
    }
    load
}

fn convection_load() -> [f64; 8] {
    let mut load = [0.; 8];
    let sum: f64 = COLD_FACE.iter().map(|&v| TARGET_T[v]).sum();
    for &v in &COLD_FACE {
        load[v] = 2. * 0.5 / 12. * (sum + TARGET_T[v] - 4. * 300.);
    }
    load
}

fn enclosure_load(epsilon: [f64; 2]) -> [f64; 8] {
    let cold = emission_moments(COLD_FACE.map(|v| TARGET_T[v]))
        .map(|moment| STEFAN_BOLTZMANN_W_M2_K4 * moment);
    let hot = emission_moments(HOT_FACE.map(|v| TARGET_T[v]))
        .map(|moment| STEFAN_BOLTZMANN_W_M2_K4 * moment);
    let eb = [
        cold.iter().sum::<f64>() / 0.5,
        hot.iter().sum::<f64>() / 0.5,
    ];
    // Solve the two gray radiosity equations analytically, retaining uniform
    // patch irradiation and exact local emission of the declared model.
    let j0 = (epsilon[0] * eb[0] + (1. - epsilon[0]) * epsilon[1] * eb[1])
        / (epsilon[0] + epsilon[1] - epsilon[0] * epsilon[1]);
    let j1 = epsilon[1] * eb[1] + (1. - epsilon[1]) * j0;
    let mut load = [0.; 8];
    for i in 0..3 {
        load[COLD_FACE[i]] = epsilon[0] * (cold[i] - 0.5 / 3. * j1);
        load[HOT_FACE[i]] = epsilon[1] * (hot[i] - 0.5 / 3. * j0);
    }
    load
}

fn manufactured_history(exterior: [f64; 8]) -> Vec<f64> {
    let mut load = exterior;
    // Full tetrahedral P1 stiffness for k=2 W/m/K, V=1/6 m^3.
    for start in [0, 4] {
        load[start] += (3. * TARGET_T[start]
            - TARGET_T[start + 1]
            - TARGET_T[start + 2]
            - TARGET_T[start + 3])
            / 3.;
        for local in 1..4 {
            load[start + local] += (TARGET_T[start + local] - TARGET_T[start]) / 3.;
        }
    }
    let mean = (0..3).map(|i| TARGET_T[i] + TARGET_T[i + 4]).sum::<f64>() / 6.;
    let resistance = 0.1 + 0.006 * (mean - 280.);
    let difference: [f64; 3] = std::array::from_fn(|i| TARGET_T[i] - TARGET_T[i + 4]);
    let sum: f64 = difference.iter().sum();
    for i in 0..3 {
        let heat = (sum + difference[i]) / (24. * resistance);
        load[i] += heat;
        load[i + 4] -= heat;
    }
    // m(h*-h_old)+dt*q(h*)=0. Manufacture history, not a thermal response;
    // the production solve must discover h* across the two latent plateaus.
    (0..8)
        .map(|v| TARGET_H[v] + DT * load[v] / if v < 4 { 1. / 24. } else { 1. / 12. })
        .collect()
}

fn shifted_chart(hot: bool, offset: f64) -> EquilibriumEnthalpyPhaseCurve {
    let original = chart(hot);
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([if hot { 32 } else { 31 }; 32]),
        original
            .knots()
            .iter()
            .map(|knot| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: knot.specific_enthalpy_j_kg + offset,
                ..*knot
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn g1_coupled_pointwise_ambient_recovers_nonuniform_latent_enthalpy_and_separate_heat_accounts() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ambient_boundary(&mesh);
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let convection = convection_load();
        for ambient in [300., 340.] {
            let radiation_load = ambient_load(ambient);
            let heat: f64 = radiation_load.iter().sum();
            let exterior = std::array::from_fn(|v| radiation_load[v] + convection[v]);
            let original_history = manufactured_history(exterior);
            assert!(original_history[0] < 100.);
            let patches =
                [
                    AmbientRadiationPatch::new("exterior", emissivity("exterior", 0.2), ambient)
                        .unwrap(),
                ];
            let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
            for offsets in [[0., 0.], [1200., -800.]] {
                let cold = shifted_chart(false, offsets[0]);
                let hot = shifted_chart(true, offsets[1]);
                let engine = HeterogeneousEnthalpyBackwardEuler::new(
                    cx,
                    &mesh,
                    &[
                        ReferenceEnthalpyMaterial {
                            curve: &cold,
                            reference_density_kg_m3: 1.,
                        },
                        ReferenceEnthalpyMaterial {
                            curve: &hot,
                            reference_density_kg_m3: 2.,
                        },
                    ],
                    &[0, 1],
                    EnthalpyBudget {
                        max_vertices: 8,
                        max_elements: 2,
                    },
                )
                .unwrap();
                let old: Vec<_> = original_history
                    .iter()
                    .enumerate()
                    .map(|(v, h)| h + offsets[v / 4])
                    .collect();
                let before = old.clone();
                let next = contacts
                    .advance_enthalpy_newton_with_pointwise_ambient_radiation(
                        cx,
                        &engine,
                        problem,
                        &old,
                        DT,
                        coupled_newton(),
                        &radiation,
                    )
                    .unwrap();
                assert_eq!(old, before);
                for v in 0..8 {
                    close(
                        next.conduction.specific_enthalpy_j_kg[v],
                        TARGET_H[v] + offsets[v / 4],
                        3e-7,
                    );
                    close(next.conduction.temperature[v], TARGET_T[v], 3e-8);
                }
                close(next.conduction.liquid_mass_fraction[0], 0.2, 3e-9);
                close(next.conduction.liquid_mass_fraction[5], 0.5, 3e-9);
                close(next.physical.convective_out_w, 10., 1e-8);
                close(next.physical.radiation_out_w, heat, 1e-8);
                close(
                    next.physical.stored_energy_change_j,
                    -DT * (10. + heat),
                    1e-8,
                );
                close(next.physical.energy_residual_j, 0., 1e-8);
                assert!(next.physical.residual_norm_j <= coupled_newton().residual_tolerance_j);
                assert_eq!(next.physical.convective_robin_fluxes.len(), 1);
                close(
                    next.conduction
                        .pointwise_radiation
                        .as_ref()
                        .unwrap()
                        .outward_heat_w,
                    heat,
                    1e-8,
                );
                assert!(next.conduction.pointwise_enclosure_radiation.is_none());
                assert!(next.radiation.is_none() && next.combined_boundary.is_none());
                assert_eq!((next.iterations, next.radiation_iterations), (1, 0));
                assert_eq!(
                    next.physical_contacts.states(),
                    next.applied_contacts.states()
                );
                close(
                    next.physical_contacts.states()[0]
                        .resistance
                        .value_m2_k_per_w(),
                    0.404,
                    1e-10,
                );
            }
            let mean = COLD_FACE.iter().map(|&v| TARGET_T[v]).sum::<f64>() / 3.;
            let mean_heat = 0.2 * STEFAN_BOLTZMANN_W_M2_K4 * 0.5 * (mean.powi(4) - ambient.powi(4));
            assert!(
                (heat - mean_heat).abs() > 0.01,
                "pointwise hot spots must not collapse to a patch mean"
            );
        }
    });
}

#[test]
fn g1_coupled_pointwise_enclosure_recovers_local_exchange_when_total_heat_cancels() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = enclosure_boundary(&mesh);
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &cold,
                    reference_density_kg_m3: 1.,
                },
                ReferenceEnthalpyMaterial {
                    curve: &hot,
                    reference_density_kg_m3: 2.,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let epsilon = [0.2, 0.7];
        let enclosure = enclosure(&mesh, epsilon);
        let radiation = PointwiseEnclosureRadiation::new(
            cx,
            &mesh,
            &boundary,
            &enclosure,
            PointwiseEnclosureTolerance::default(),
        )
        .unwrap();
        let load = enclosure_load(epsilon);
        let old = manufactured_history(load);
        let before = old.clone();
        let next = contacts
            .advance_enthalpy_newton_with_pointwise_enclosure_radiation(
                cx,
                &engine,
                problem,
                &old,
                DT,
                coupled_newton(),
                &radiation,
            )
            .unwrap();
        assert_eq!(old, before);
        for v in 0..8 {
            close(next.conduction.specific_enthalpy_j_kg[v], TARGET_H[v], 3e-7);
            close(next.conduction.temperature[v], TARGET_T[v], 3e-8);
        }
        close(next.physical.stored_energy_change_j, 0., 1e-8);
        close(next.physical.radiation_out_w, 0., 1e-8);
        close(next.physical.energy_residual_j, 0., 1e-8);
        assert!(next.physical.residual_norm_j <= coupled_newton().residual_tolerance_j);
        let report = next
            .conduction
            .pointwise_enclosure_radiation
            .as_ref()
            .unwrap();
        close(report.outward_heat_w, 0., 1e-8);
        let left: f64 = COLD_FACE.iter().map(|&v| load[v]).sum();
        assert!(
            left < -10.,
            "closed enclosure must retain internal heat exchange"
        );
        close(report.patches[0].outward_heat_w, left, 1e-8);
        close(report.patches[1].outward_heat_w, -left, 1e-8);
        assert!(next.conduction.pointwise_radiation.is_none());
        assert_eq!((next.iterations, next.radiation_iterations), (1, 0));
    });
}

#[test]
fn g4_coupled_pointwise_contact_refuses_wrong_binding_and_unsupported_surface_without_history_changes()
 {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let boundary = ambient_boundary(&mesh);
        let other_boundary = ambient_boundary(&mesh);
        let card = card();
        let contacts =
            TemperatureContacts::new(cx, &mesh, &boundary, vec![], vec![law(&mesh, &card)], 1)
                .unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx,
            &mesh,
            &[
                ReferenceEnthalpyMaterial {
                    curve: &cold,
                    reference_density_kg_m3: 1.,
                },
                ReferenceEnthalpyMaterial {
                    curve: &hot,
                    reference_density_kg_m3: 2.,
                },
            ],
            &[0, 1],
            EnthalpyBudget {
                max_vertices: 8,
                max_elements: 2,
            },
        )
        .unwrap();
        let material = ConductivityModel::isotropic_declared(2.).unwrap();
        let source = ScalarField::Uniform(0.);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let old = TARGET_H;
        let patches =
            [AmbientRadiationPatch::new("exterior", emissivity("exterior", 0.2), 300.).unwrap()];
        let other = PointwiseAmbientRadiation::new(cx, &mesh, &other_boundary, &patches).unwrap();
        assert!(matches!(
            contacts.advance_enthalpy_newton_with_pointwise_ambient_radiation(
                cx,
                &engine,
                problem,
                &old,
                DT,
                coupled_newton(),
                &other
            ),
            Err(EnthalpyError::InvalidInput(
                "pointwise radiation must retain the exact transport mesh and boundary"
            ))
        ));
        let unsupported = [AmbientRadiationPatch::new(
            "exterior",
            emissivity_in("exterior", 0.2, 306., 380.),
            300.,
        )
        .unwrap()];
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &unsupported).unwrap();
        assert!(
            matches!(contacts.advance_enthalpy_newton_with_pointwise_ambient_radiation(
            cx, &engine, problem, &old, DT, coupled_newton(), &radiation),
            Err(EnthalpyError::Conduction(ConductionError::MaterialQuery { property, .. }))
                if property == SURFACE_EMISSIVITY_PROPERTY)
        );
        let radiation = PointwiseAmbientRadiation::new(cx, &mesh, &boundary, &patches).unwrap();
        let gate = CancelGate::new_clock_free();
        gate.request();
        context(&gate, |cancelled| {
            assert!(matches!(
                contacts.advance_enthalpy_newton_with_pointwise_ambient_radiation(
                    cancelled,
                    &engine,
                    problem,
                    &old,
                    DT,
                    coupled_newton(),
                    &radiation
                ),
                Err(EnthalpyError::Conduction(ConductionError::Cancelled { .. }))
            ));
        });
        assert_eq!(old, TARGET_H);
    });
}

