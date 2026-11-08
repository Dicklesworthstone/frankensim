//! G1/G3: real temporal and spatial quadrature on the original solved AC fields.
use super::*;
use fs_em::spatial::vector_potential::induction::{
    harmonic::{
        SpatialHarmonicInductionResponse,
        coherent_exposure::{SpatialHarmonicInductionExposureError as Error, SpatialHarmonicInductionExposureLimits as Limits},
    },
    mechanical::SpatialLorentzLoadLimits,
};
use fs_fft::harmonic_window::{HarmonicWindow, HarmonicWindowProfile};
use fs_material::conductor::{
    ElectricalPropertyBasis,
    spatial::{ResolvedOrthotropicConductor3, SpatialConductorMaterial},
};

fn limits() -> Limits {
    Limits {
        maximum_lines: 8,
        maximum_pair_products: 1_000_000,
        maximum_workspace_bytes: 1_000_000,
    }
}
const LOADS: SpatialLorentzLoadLimits = SpatialLorentzLoadLimits {
    max_vertices: 8,
    max_cells: 6,
};

// Independent composite Simpson time rule and four-point degree-two
// tetrahedron quadrature. Reconstruct REAL fields before multiplying; neither
// the production Gram kernel nor its P1 mass formula supplies the oracle.
fn quadrature<C: SpatialConductorMaterial>(
    responses: &[&SpatialHarmonicInductionResponse<'_, '_, C>],
    velocity: &[[f64; 3]],
    window: HarmonicWindow,
) -> (Vec<f64>, Vec<[f64; 3]>, f64) {
    let problem = responses[0].problem();
    let mut heat = vec![0.0; responses[0].cells().len()];
    let mut nodal = vec![[0.0; 3]; problem.mesh().vertex_count];
    let mut mechanical = 0.0;
    let high = (5.0 + 3.0 * 5.0_f64.sqrt()) / 20.0;
    let low = (5.0 - 5.0_f64.sqrt()) / 20.0;
    let steps = 2048;
    for step in 0..=steps {
        let time = window.start_s() + window.duration_s() * step as f64 / steps as f64;
        let temporal_weight = if step == 0 || step == steps { 1.0 } else if step % 2 == 0 { 2.0 } else { 4.0 } / (3 * steps) as f64;
        let phases: Vec<_> = responses.iter().map(|r| (r.omega_rad_s() * time).sin_cos()).collect();
        for (index, original) in responses[0].cells().iter().enumerate() {
            let tet = problem.mesh().tets[original.cell].map(|v| v as usize);
            for point in 0..4 {
                let shape: [f64; 4] = std::array::from_fn(|i| if i == point { high } else { low });
                let (mut j, mut e, mut b, mut v) = ([0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]);
                for (response, &(sine, cosine)) in responses.iter().zip(&phases) {
                    let real = |z: C64| z.re * cosine + z.im * sine;
                    let cell = &response.cells()[index];
                    for axis in 0..3 {
                        b[axis] += real(response.flux_density_t()[original.cell][axis]);
                        for vertex in 0..4 {
                            j[axis] += shape[vertex] * real(cell.current_density_a_m2[vertex][axis]);
                            e[axis] += shape[vertex] * real(cell.conductor_electric_field_v_m[vertex][axis]);
                        }
                    }
                }
                for axis in 0..3 {
                    for vertex in 0..4 { v[axis] += shape[vertex] * velocity[tet[vertex]][axis]; }
                }
                let force = cross(j, b);
                let weight = temporal_weight * original.volume_m3 / 4.0;
                for axis in 0..3 {
                    heat[index] += weight * j[axis] * e[axis];
                    mechanical += weight * force[axis] * v[axis];
                    for vertex in 0..4 {
                        nodal[tet[vertex]][axis] += weight * shape[vertex] * force[axis];
                    }
                }
            }
        }
    }
    (heat, nodal, mechanical)
}

fn check_coherent<C: SpatialConductorMaterial>(
    problem: &SpatialVectorMagnetostatics,
    sources: Vec<Option<C>>,
) {
    let initial = problem.solve(
        &[[0.0; 3]; 6], &trace(problem, [0.0; 3], false), static_controls(), &mut || false,
    ).unwrap();
    let state = SpatialInduction::new_with_conductors(
        &initial, sources, preparation(), 2e-9, &mut || false,
    ).unwrap();
    for alpha in [0.0, 0.35] {
        let velocity: Vec<_> = problem.positions_m().iter().map(|p| p.map(|x| alpha * x)).collect();
        let a = state.solve_harmonic(
            &[[C64::new(4.0, -1.2); 3]; 6], &velocity, 3.0, harmonic_controls(), &mut || false,
        ).unwrap();
        let b = state.solve_harmonic(
            &[[C64::new(2.5, 1.7); 3]; 6], &velocity, 4.0, harmonic_controls(), &mut || false,
        ).unwrap();
        let bank = state.prepare_harmonic_exposure(&[&b, &a], limits(), &mut || false).unwrap();
        assert_eq!(bank.distinct_frequencies(), 2);
        assert!(std::ptr::eq(bank.source(), &state));
        let window = HarmonicWindow::new(0.13, 0.47).unwrap();
        let exposure = bank.observe(window, &mut || false).unwrap();
        let (heat, nodal, mechanical) = quadrature(&[&a, &b], &velocity, window);
        near(exposure.mean_ohmic_power_w, heat.iter().sum(), 2e-9);
        near(exposure.ohmic_heat_j, window.duration_s() * heat.iter().sum::<f64>(), 2e-9);
        near(exposure.mean_mechanical_power_w, mechanical, 2e-9);
        near(exposure.mechanical_work_j, window.duration_s() * mechanical, 2e-9);
        for (record, expected) in exposure.cells().iter().zip(heat) {
            near(record.mean_ohmic_power_w, expected, 2e-9);
            near(record.ohmic_heat_j, expected * window.duration_s(), 2e-9);
            assert!(std::ptr::eq(record.source, state.conductors()[record.cell].as_ref().unwrap()));
        }
        let loads = exposure.lorentz_loads(&[5, 2, 3, 0, 1, 4], LOADS, &mut || false).unwrap();
        for (actual, expected) in loads.nodal_forces_n().iter().zip(nodal) {
            for axis in 0..3 { near(actual[axis], expected[axis], 2e-9); }
        }
        near(loads.power_w(&velocity, &mut || false).unwrap(), mechanical, 2e-9);
        let separate_means = a.balance().mean_ohmic_power_w + b.balance().mean_ohmic_power_w;
        assert!((exposure.mean_ohmic_power_w - separate_means).abs() > 0.02);
        // Unweighted deposited energy is additive over consecutive intervals.
        let left = bank.observe(HarmonicWindow::new(0.13, 0.2).unwrap(), &mut || false).unwrap();
        let right = bank.observe(HarmonicWindow::new(0.33, 0.27).unwrap(), &mut || false).unwrap();
        near(exposure.ohmic_heat_j, left.ohmic_heat_j + right.ohmic_heat_j, 2e-12);
        near(exposure.mechanical_work_j, left.mechanical_work_j + right.mechanical_work_j, 2e-12);
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(state.potential_wb(), initial.potential_wb());
    }
}

#[test]
fn coherent_three_dimensional_heat_and_forcing_match_real_space_time_quadrature() {
    let (problem, sources) = cube(1, None);
    check_coherent(&problem, sources);
}

#[test]
fn rotated_orthotropic_coherent_heat_retains_the_original_tensor_source() {
    let names = ["sigma-1", "sigma-2", "sigma-3"];
    let mut claims = ClaimSet::new();
    let mut requirements = vec![];
    let values = [
        (MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, 2.0),
        (names[0], ELECTRICAL_CONDUCTIVITY_DIMS, 2.0),
        (names[1], ELECTRICAL_CONDUCTIVITY_DIMS, 5.0),
        (names[2], ELECTRICAL_CONDUCTIVITY_DIMS, 11.0),
    ];
    for (name, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(), interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated, observations: vec![],
            provenance: Provenance { source: "manufactured coherent tensor conductor".into(), license: "test-only".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap());
    }
    let card = MaterialCard::assemble(
        MaterialStateId { chemistry: "synthetic tensor".into(), phase: "solid".into(), process: "fixed axes".into(), revision: 0 },
        claims, vec![],
    ).unwrap();
    let material = resolve_material_state_point(&card, &QueryPoint::new(), &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let c = 0.5_f64.sqrt();
    let conductor = ResolvedOrthotropicConductor3::try_from_material(
        material.clone(), ElectricalPropertyBasis::Conductivity,
        names.map(|name| PropertyKey::new(name, ELECTRICAL_CONDUCTIVITY_DIMS)),
        [[c, c, 0.0], [-c, c, 0.0], [0.0, 0.0, 1.0]],
    ).unwrap();
    let (positions, tets) = grid(1);
    let problem = SpatialVectorMagnetostatics::new(
        positions, tets, vec![SpatialMagneticMaterial::Linear(material); 6], assembly(), &mut || false,
    ).unwrap();
    check_coherent(&problem, vec![Some(conductor); 6]);
}

#[test]
fn full_cycles_recover_means_and_equal_frequency_fields_interfere_before_heating() {
    let (problem, sources) = cube(1, None);
    let initial = problem.solve(
        &[[0.0; 3]; 6], &trace(&problem, [0.0; 3], false), static_controls(), &mut || false,
    ).unwrap();
    let state = SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let velocity: Vec<_> = problem.positions_m().iter().map(|p| p.map(|x| 0.35 * x)).collect();
    let a = state.solve_harmonic(&[[C64::new(4.0, -2.0); 3]; 6], &velocity, 3.0, harmonic_controls(), &mut || false).unwrap();
    let opposite = state.solve_harmonic(&[[C64::new(-4.0, 2.0); 3]; 6], &velocity, 3.0, harmonic_controls(), &mut || false).unwrap();
    let close_frequency = state.solve_harmonic(&[[C64::new(4.0, -2.0); 3]; 6], &velocity, 3.0 + 1e-7, harmonic_controls(), &mut || false).unwrap();
    let period = std::f64::consts::TAU / 3.0;
    let window = HarmonicWindow::new(0.17, period).unwrap();
    let single_bank = state.prepare_harmonic_exposure(&[&a], limits(), &mut || false).unwrap();
    let single = single_bank.observe(window, &mut || false).unwrap();
    near(single.mean_ohmic_power_w, a.balance().mean_ohmic_power_w, 2e-9);
    near(single.mean_mechanical_power_w, a.balance().mean_mechanical_power_w, 2e-9);
    let original = a.lorentz_loads(&[0, 1, 2, 3, 4, 5], LOADS, &mut || false).unwrap();
    let loads = single.lorentz_loads(&[0, 1, 2, 3, 4, 5], LOADS, &mut || false).unwrap();
    for (actual, expected) in loads.nodal_forces_n().iter().zip(original.mean().nodal_forces_n()) {
        for axis in 0..3 { near(actual[axis], expected[axis], 2e-9); }
    }
    let doubled_bank = state.prepare_harmonic_exposure(&[&a, &a], limits(), &mut || false).unwrap();
    assert_eq!(doubled_bank.distinct_frequencies(), 1);
    let doubled = doubled_bank.observe(window, &mut || false).unwrap();
    near(doubled.ohmic_heat_j, 4.0 * single.ohmic_heat_j, 2e-12);
    let cancelled_bank = state.prepare_harmonic_exposure(&[&a, &opposite], limits(), &mut || false).unwrap();
    assert_eq!(cancelled_bank.distinct_frequencies(), 1);
    let cancelled = cancelled_bank.observe(HarmonicWindow::new(0.31, 0.013).unwrap(), &mut || false).unwrap();
    assert!(cancelled.ohmic_heat_j < 1e-24);
    assert!(cancelled.cells().iter().flat_map(|c| c.mean_lorentz_force_density_n_m3).flatten().all(|f| f.abs() < 1e-23));
    assert_eq!(state.prepare_harmonic_exposure(&[&a, &close_frequency], limits(), &mut || false).unwrap().distinct_frequencies(), 2);
}

#[test]
fn original_source_motion_budgets_and_atomic_cancellation_are_enforced() {
    let (problem, sources) = cube(1, None);
    let initial = problem.solve(
        &[[0.0; 3]; 6], &trace(&problem, [0.0; 3], false), static_controls(), &mut || false,
    ).unwrap();
    let state = SpatialInduction::new(&initial, sources, preparation(), 2e-9, &mut || false).unwrap();
    let other = state.clone();
    let current = [[C64::new(4.0, 1.0); 3]; 6];
    let a = state.solve_harmonic(&current, &[[0.0; 3]; 8], 3.0, harmonic_controls(), &mut || false).unwrap();
    let alien = other.solve_harmonic(&current, &[[0.0; 3]; 8], 3.0, harmonic_controls(), &mut || false).unwrap();
    let moved = state.solve_harmonic(&current, &[[0.1, 0.0, 0.0]; 8], 3.0, harmonic_controls(), &mut || false).unwrap();
    for responses in [&[&a, &alien][..], &[&a, &moved][..], &[][..]] {
        assert!(matches!(state.prepare_harmonic_exposure(responses, limits(), &mut || false), Err(Error::Invalid(_))));
    }
    for limit in [
        Limits { maximum_lines: 0, ..limits() },
        Limits { maximum_pair_products: 0, ..limits() },
        Limits { maximum_workspace_bytes: 0, ..limits() },
    ] {
        assert!(state.prepare_harmonic_exposure(&[&a], limit, &mut || false).is_err());
    }
    let exact = Limits { maximum_pair_products: 39 * 6 * 4, ..limits() };
    assert!(state.prepare_harmonic_exposure(&[&a, &a], exact, &mut || false).is_ok());
    assert!(state.prepare_harmonic_exposure(&[&a, &a], Limits {
        maximum_pair_products: exact.maximum_pair_products - 1, ..exact
    }, &mut || false).is_err());
    // Coherent grouping cannot bypass the authored-line hard ceiling.
    assert!(state.prepare_harmonic_exposure(&vec![&a; 257], Limits {
        maximum_lines: usize::MAX, maximum_pair_products: usize::MAX,
        maximum_workspace_bytes: usize::MAX,
    }, &mut || false).is_err());
    assert!(matches!(state.prepare_harmonic_exposure(&[&a], limits(), &mut || true), Err(Error::Cancelled)));
    let bank = state.prepare_harmonic_exposure(&[&a], limits(), &mut || false).unwrap();
    let window = HarmonicWindow::new(0.2, 0.7).unwrap();
    let first = bank.observe(window, &mut || false).unwrap();
    let mut calls = 0;
    assert!(matches!(bank.observe(window, &mut || { calls += 1; calls == 7 }), Err(Error::Cancelled)));
    assert_eq!(bank.observe(window, &mut || false).unwrap().ohmic_heat_j.to_bits(), first.ohmic_heat_j.to_bits());
    let shaped = HarmonicWindow::with_profile(0.2, 0.7, HarmonicWindowProfile::Triangular).unwrap();
    assert!(bank.observe(shaped, &mut || false).is_err());
    for selection in [&[][..], &[0, 0][..], &[6][..]] {
        assert!(first.lorentz_loads(selection, LOADS, &mut || false).is_err());
    }
    assert!(matches!(first.lorentz_loads(&[0], LOADS, &mut || true), Err(Error::Cancelled)));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.potential_wb(), initial.potential_wb());
}
