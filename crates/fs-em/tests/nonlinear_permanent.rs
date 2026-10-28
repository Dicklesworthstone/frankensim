//! G1/G3/G4: fixed-remanence sources in nonlinear fields and actual induction.
//! These manufactured constitutive laws establish implementation behavior,
//! not experimental qualification of a magnetic material or machine.

mod common;
#[path = "nonlinear_eddy/common.rs"]
mod induction;
#[path = "nonlinear_winding_circuit/common.rs"]
mod winding;

use fs_em::planar::{
    PlanarMagnetostatics,
    eddy::{
        EddyConductor,
        nonlinear::{NonlinearEddyError, PlanarNonlinearEddy},
    },
    nonlinear::{IsotropicMagneticMaterial, PlanarNonlinearMagnetostatics},
    permanent::{MagneticEnergyReference, PlanarMagneticMaterial, PlanarRecoilMaterial},
};
use fs_em::winding::{
    nonlinear::{NonlinearWindingCircuit, NonlinearWindingInitial, NonlinearWindings},
    ohmic::WindingConductor,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS,
        permanent::{
            FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
        },
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_rep_mesh::{Metric2, TriComplex2};

fn recoil(
    mu: f64,
    br: f64,
    parallel: [f64; 2],
    transverse: f64,
    sigma: f64,
) -> (ResolvedRecoilMaterial, ResolvedConductorState) {
    let basis = ElectricalPropertyBasis::Conductivity;
    let entries = [
        (RECOIL_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, mu),
        (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, br),
        (
            RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            parallel[0],
        ),
        (
            RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            parallel[1],
        ),
        (
            RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY,
            FIELD_STRENGTH_DIMS,
            transverse,
        ),
        (basis.property_name(), basis.dimensions(), sigma),
    ];
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in entries {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured conductive permanent-magnet law".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured magnet".into(),
            phase: "solid".into(),
            process: "nonlinear-recoil-regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    (
        ResolvedRecoilMaterial::try_from_material(material.clone()).unwrap(),
        ResolvedConductorState::try_from_material_with_basis(material, basis).unwrap(),
    )
}

fn triangle() -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "nonlinear-natural-recoil-source",
        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        vec![[0, 1, 2]],
        Metric2::planar(0.2).unwrap(),
    )
    .unwrap()
}

#[test]
fn nonlinear_owner_retains_the_recoil_natural_load_and_actual_h_energy() {
    let mu = 2e-6;
    let br = 0.8;
    for sign in [-1.0, 1.0] {
        for j in [0.0, 3e5] {
            let (law, _) = recoil(mu, br, [-1e6, 1e6], 1e6, 18.0);
            let magnet = PlanarRecoilMaterial::new(law, [sign, 0.0]).unwrap();
            let problem = PlanarNonlinearMagnetostatics::new(
                triangle(),
                vec![IsotropicMagneticMaterial::Recoil(magnet.clone())],
                common::limits(),
                &mut || false,
            )
            .unwrap();
            let linear = PlanarMagnetostatics::new_with_material_models(
                triangle(),
                vec![PlanarMagneticMaterial::Recoil(magnet)],
                common::limits(),
                &mut || false,
            )
            .unwrap();
            let boundary = [(0, 0.0), (1, 0.0)];
            let field = problem
                .solve(
                    &[j],
                    &boundary,
                    &[0.0; 3],
                    induction::controls().nonlinear,
                    &mut || false,
                )
                .unwrap();
            let reference = linear
                .solve(&[j], &boundary, common::config(), &mut || false)
                .unwrap();
            // Only the bottom edge is prescribed. The exact P1 equation is
            // B_x = signed Br + mu J/3, H_x = J/3. An interior curl(Br)
            // source alone misses this natural-boundary contribution entirely.
            let b = sign * br + mu * j / 3.0;
            let h = j / 3.0;
            common::near(field.potential_wb_m()[2], b);
            common::near(field.magnetic_flux_density_t()[0][0], b);
            common::near(field.magnetic_field_a_m()[0][0], h);
            common::near(field.magnetic_field_a_m()[0][1], 0.0);
            let energy = 0.1 * 0.5 * mu * h * h;
            common::near(field.energy_j(), energy);
            common::near(field.coenergy_j(), 0.1 * b * h - energy);
            assert_eq!(
                field.energy_reference(),
                MagneticEnergyReference::FixedRemanence
            );
            common::near(
                field.recoil_operating_points()[0]
                    .unwrap()
                    .parallel_field_a_m,
                sign * h,
            );
            common::near(field.energy_j(), reference.energy_j());
            for (actual, expected) in field
                .potential_wb_m()
                .iter()
                .zip(reference.potential_wb_m())
            {
                common::near(*actual, *expected);
            }
            assert!(field.report().residual_norm_a_m <= field.report().residual_target_a_m);
        }
    }
}

#[test]
fn nonlinear_iron_and_permanent_magnets_share_the_source_derived_layered_solution() {
    let br = 2.5;
    for saturating in [false, true] {
        let knots = if saturating {
            vec![(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)]
        } else {
            vec![(0.0, 0.0), (4.0, 4.0)]
        };
        let core = induction::law(&induction::source(&knots, 18.0));
        for sign in [-1.0, 1.0] {
            let (law, _) = recoil(1.0, br, [-5.0, 5.0], 5.0, 18.0);
            let magnet = PlanarRecoilMaterial::new(law, [0.0, sign]).unwrap();
            let problem = PlanarNonlinearMagnetostatics::new(
                common::layered(1.0, 2.0, 0.3),
                vec![
                    IsotropicMagneticMaterial::Recoil(magnet.clone()),
                    IsotropicMagneticMaterial::Recoil(magnet),
                    IsotropicMagneticMaterial::Reversible(core.clone()),
                    IsotropicMagneticMaterial::Reversible(core.clone()),
                ],
                common::limits(),
                &mut || false,
            )
            .unwrap();
            // Equal layer widths and equal endpoint A require B_core=-B_magnet.
            // Tangential H continuity gives q + H_core(q)=Br. The two source
            // laws therefore give q=1.25 and q=1.1 independently of FE assembly.
            let q = if saturating { 1.1 } else { 1.25 };
            let h = -sign * (br - q);
            for gauge in [0.0, 7.0] {
                let boundary: Vec<_> = common::ends(0.0)
                    .into_iter()
                    .map(|(v, a)| (v, a + gauge))
                    .collect();
                let field = problem
                    .solve(
                        &[0.0; 4],
                        &boundary,
                        &[gauge; 6],
                        induction::controls().nonlinear,
                        &mut || false,
                    )
                    .unwrap();
                for cell in 0..4 {
                    common::near(field.magnetic_flux_density_t()[cell][0], 0.0);
                    common::near(
                        field.magnetic_flux_density_t()[cell][1],
                        if cell < 2 { sign * q } else { -sign * q },
                    );
                    common::near(field.magnetic_field_a_m()[cell][0], 0.0);
                    common::near(field.magnetic_field_a_m()[cell][1], h);
                }
                // The chosen zero is B=Br in magnet cells: its Legendre dual
                // can be negative. Here the total B.H integral vanishes, so
                // total coenergy is exactly minus total constitutive energy.
                assert!(field.cell_coenergy_j()[0] < 0.0);
                common::near(field.energy_j() + field.coenergy_j(), 0.0);
                let core_psi = if saturating {
                    0.5 + (q - 1.0) + 2.0 * (q - 1.0) * (q - 1.0)
                } else {
                    0.5 * q * q
                };
                common::near(field.energy_j(), 0.6 * (0.5 * h * h + core_psi));
            }
        }
    }
}

#[test]
fn conductive_recoil_history_relaxes_toward_its_nonzero_bias_with_separate_losses() {
    let depth = 0.3;
    let mu = 1.0;
    let br = 0.5;
    let sigma = 18.0;
    let (law, conductor) = recoil(mu, br, [-5.0, 5.0], 5.0, sigma);
    let directions = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]];
    let materials: Vec<_> = directions
        .into_iter()
        .map(|direction| {
            IsotropicMagneticMaterial::Recoil(
                PlanarRecoilMaterial::new(law.clone(), direction).unwrap(),
            )
        })
        .collect();
    let problem = PlanarNonlinearMagnetostatics::new(
        induction::mesh(depth),
        materials.clone(),
        common::limits(),
        &mut || false,
    )
    .unwrap();
    let boundary = induction::boundary(0.0);
    let initial = problem
        .solve(
            &[6.0; 4],
            &boundary,
            &[0.0; 5],
            induction::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let a0 = 0.75;
    let equilibrium = br / 2.0;
    common::near(initial.potential_wb_m()[4], a0);
    let conductors: Vec<_> = (0..4)
        .map(|triangle| EddyConductor {
            triangle,
            region: 7,
            conductor: conductor.clone(),
        })
        .collect();
    let original = PlanarNonlinearEddy::new(
        &initial,
        conductors,
        &[],
        &boundary,
        0.0,
        induction::controls(),
        &mut || false,
    )
    .unwrap();
    let c = sigma * depth / 18.0;
    let k = 4.0 * depth / mu;
    let dt = 0.1;
    let a1 = equilibrium + (a0 - equilibrium) * c / (c + dt * k);
    let candidate = original
        .step_to_candidate(
            dt,
            original.conductors().to_vec(),
            &[0.0; 4],
            induction::controls(),
            &mut || false,
        )
        .unwrap();
    let step = candidate.step();
    common::near(candidate.state().potential_wb_m()[4], a1);
    common::near(step.stored_before_j, 0.5 * k * (a0 - equilibrium).powi(2));
    common::near(step.stored_after_j, 0.5 * k * (a1 - equilibrium).powi(2));
    common::near(step.heat_j, c * (a1 - a0).powi(2) / dt);
    common::near(step.numerical_dissipation_j, 0.5 * k * (a1 - a0).powi(2));
    common::near(step.supplied_work_j, 0.0);
    common::near(step.parameter_work_j, 0.0);
    common::near(step.energy_residual_j, 0.0);
    assert!(step.heat_j > 0.0 && step.numerical_dissipation_j > 0.0);
    assert_eq!(step.region_currents_a.len(), 1);
    common::near(step.region_currents_a[0].1, 0.0);
    for cell in &step.magnetic_cells {
        let db: [f64; 2] = [0, 1].map(|i| cell.flux_density_t[1][i] - cell.flux_density_t[0][i]);
        common::near(
            cell.numerical_dissipation_j,
            cell.volume_m3 * (db[0] * db[0] + db[1] * db[1]) / (2.0 * mu),
        );
    }
    let continued = candidate
        .state()
        .step_to_candidate(
            2.0 * dt,
            candidate.state().conductors().to_vec(),
            &[0.0; 4],
            induction::controls(),
            &mut || false,
        )
        .unwrap();
    common::near(
        continued.state().potential_wb_m()[4],
        equilibrium + (a1 - equilibrium) * c / (c + dt * k),
    );
    assert!(continued.state().potential_wb_m()[4] > equilibrium);

    // A fixed-source continuation may not silently rotate remanence or debit
    // its changed energy as Ohmic heat. Failed trials retain the accepted state.
    let mut changed = materials;
    changed[0] =
        IsotropicMagneticMaterial::Recoil(PlanarRecoilMaterial::new(law, [-1.0, 0.0]).unwrap());
    assert!(
        original
            .step_with_materials_candidate(
                dt,
                changed,
                original.conductors().to_vec(),
                &[0.0; 4],
                common::limits(),
                induction::controls(),
                &mut || false,
            )
            .is_err()
    );
    assert!(matches!(
        original.step_to_candidate(
            dt,
            original.conductors().to_vec(),
            &[0.0; 4],
            induction::controls(),
            &mut || true,
        ),
        Err(NonlinearEddyError::Cancelled)
    ));
    common::near(original.time_s(), 0.0);
    common::near(original.potential_wb_m()[4], a0);
    common::near(original.energy_j(), step.stored_before_j);
}

#[test]
fn native_winding_circuit_preserves_fixed_magnet_flux_without_inventing_source_work() {
    let depth = 0.3;
    let br = 0.5;
    let mu = 1.0;
    let sigma = 18.0;
    let (law, conductor) = recoil(mu, br, [-5.0, 5.0], 5.0, sigma);
    let directions = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]];
    let materials = (0..8)
        .map(|cell| {
            let orientation = if cell < 4 { 1.0 } else { -1.0 };
            let direction = directions[cell % 4].map(|x| orientation * x);
            IsotropicMagneticMaterial::Recoil(
                PlanarRecoilMaterial::new(law.clone(), direction).unwrap(),
            )
        })
        .collect();
    let problem = PlanarNonlinearMagnetostatics::new(
        winding::mesh(depth),
        materials,
        common::limits(),
        &mut || false,
    )
    .unwrap();
    let profiles = winding::profiles(1);
    let conductors: Vec<_> = (0..8)
        .map(|triangle| WindingConductor {
            triangle,
            conductor: conductor.clone(),
        })
        .collect();
    let graph_resistance = 0.3;
    let original = NonlinearWindingCircuit::new(
        winding::initial_field(&problem, &profiles, &[0.0]),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (
                    2,
                    0,
                    Branch::Resistor {
                        ohms: graph_resistance,
                    },
                ),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: profiles,
            end_inductance_h: vec![0.0],
        },
        conductors.clone(),
        NonlinearWindingInitial {
            currents_a: vec![0.0],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        winding::config(),
        &mut || false,
    )
    .unwrap();
    let bias = 2.0 * depth * br;
    let inductance = 2.0 * depth * mu;
    let winding_resistance = 72.0 * depth / sigma;
    common::near(original.flux_linkages_wb()[0], bias);
    common::near(original.field().energy_j(), 0.0);
    let dt = 0.2;
    let voltage = 9.0;
    let current = voltage / (winding_resistance + graph_resistance + inductance / dt);
    let next = original
        .step_to_candidate(
            dt,
            vec![],
            conductors,
            &[voltage],
            winding::controls(),
            &mut || false,
        )
        .unwrap();
    winding::near(next.circuit().currents_a()[0], current);
    winding::near(
        next.circuit().flux_linkages_wb()[0],
        bias + inductance * current,
    );
    winding::near(
        next.circuit().field().potential_wb_m()[6],
        (br + mu * current) / 2.0,
    );
    winding::near(
        next.circuit().field().potential_wb_m()[7],
        -(br + mu * current) / 2.0,
    );
    let step = next.step();
    winding::near(step.stored_before_j, 0.0);
    winding::near(step.stored_after_j, 0.5 * inductance * current * current);
    winding::near(
        step.numerical_dissipation_j,
        0.5 * inductance * current * current,
    );
    winding::near(
        step.heat_j,
        (winding_resistance + graph_resistance) * current * current * dt,
    );
    winding::near(step.supplied_work_j, voltage * current * dt);
    winding::near(step.energy_residual_j, 0.0);
    // The unchanged bias is retained linkage, never an extra voltage source.
    common::near(original.currents_a()[0], 0.0);
    common::near(original.flux_linkages_wb()[0], bias);
}

#[test]
fn nonlinear_fields_cannot_publish_recoil_points_outside_the_original_envelope() {
    let (law, _) = recoil(2e-6, 0.8, [-1e5, 2e5], 1e4, 18.0);
    let problem = PlanarNonlinearMagnetostatics::new(
        triangle(),
        vec![IsotropicMagneticMaterial::Recoil(
            PlanarRecoilMaterial::new(law, [1.0, 0.0]).unwrap(),
        )],
        common::limits(),
        &mut || false,
    )
    .unwrap();
    for (at_x, at_y) in [(0.0, 0.4), (0.0, 1.4), (-0.04, 0.8)] {
        let potential = [0.0, at_x, at_y];
        assert!(
            problem
                .solve(
                    &[0.0],
                    &[(0, 0.0), (1, at_x), (2, at_y)],
                    &potential,
                    induction::controls().nonlinear,
                    &mut || false,
                )
                .is_err()
        );
    }
    // A valid remanent state remains admissible after all rejected candidates.
    let field = problem
        .solve(
            &[0.0],
            &[(0, 0.0), (1, 0.0), (2, 0.8)],
            &[0.0, 0.0, 0.8],
            induction::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    common::near(field.energy_j(), 0.0);
    common::near(field.magnetic_field_a_m()[0][0], 0.0);
}
