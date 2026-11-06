//! Local 3D flux derivatives, source domains and gauge/polarity invariance.
//! Manufactured material curves test the implementation, not measured coils.
mod common;
#[path = "nonlinear_eddy/common.rs"]
mod magnetic;

use common::{material, near};
use fs_em::spatial::{
    SpatialMagneticError, SpatialMagneticMaterial, SpatialRecoilMaterial,
    nonlinear::SpatialNonlinearMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
        VectorPotentialFace,
        nonlinear::{
            NonlinearVectorControls, NonlinearVectorError, SpatialNonlinearVectorMagnetostatics,
            winding::{
                SpatialIncrementalWindingBasis, SpatialIncrementalWindingConfig,
                SpatialIncrementalWindingError,
            },
        },
        winding::{SpatialWindingBasis, SpatialWindingConfig},
    },
};
use fs_material::field_storage::nonlinear::ReversibleMagneticMaterial;
use fs_solver::CheckedCgConfig;

fn controls() -> NonlinearVectorControls {
    let solve = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 3000,
        max_corrections: 2,
    };
    NonlinearVectorControls {
        field: VectorMagneticControls {
            field: solve,
            gauge: solve,
            maximum_gauge_product_terms: 2_000_000,
            current_tolerance_a: 2e-9,
            flux_tolerance_wb: 2e-9,
            energy_tolerance_j: 2e-8,
        },
        relative_residual_tolerance: 1e-10,
        maximum_newton_iterations: 40,
        maximum_backtracks: 35,
    }
}
fn config() -> SpatialIncrementalWindingConfig {
    SpatialIncrementalWindingConfig {
        winding: SpatialWindingConfig {
            max_windings: 4,
            max_basis_values: 100_000,
            pairing_tolerance_h: 2e-9,
        },
        minimum_smooth_field_margin_t: 1e-7,
    }
}
fn law() -> ReversibleMagneticMaterial {
    magnetic::law(&magnetic::source(&[(0.0, 0.0), (1.0, 1.0), (4.0, 13.0)], 2.0))
}
fn problem(material: SpatialNonlinearMaterial) -> SpatialNonlinearVectorMagnetostatics {
    let mut positions = Vec::new();
    for k in 0..=2 {
        for j in 0..=2 {
            for i in 0..=2 {
                positions.push([f64::from(i) / 2.0, f64::from(j) / 2.0, f64::from(k) / 2.0]);
            }
        }
    }
    let id = |i, j, k| ((k * 3 + j) * 3 + i) as u32;
    let mut tets = Vec::new();
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let (a, b, c, d) = (
                    id(i, j, k), id(i + 1, j, k),
                    id(i, j + 1, k), id(i + 1, j + 1, k),
                );
                let (e, f, g, h) = (
                    id(i, j, k + 1), id(i + 1, j, k + 1),
                    id(i, j + 1, k + 1), id(i + 1, j + 1, k + 1),
                );
                tets.extend([
                    [a, b, d, h], [a, d, c, h], [a, c, g, h],
                    [a, g, e, h], [a, e, f, h], [a, f, b, h],
                ]);
            }
        }
    }
    let materials = vec![material; tets.len()];
    SpatialNonlinearVectorMagnetostatics::new(
        positions, tets, materials, common::limits(), &mut || false,
    ).unwrap()
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1],
     a[2] * b[0] - a[0] * b[2],
     a[0] * b[1] - a[1] * b[0]]
}
// The same independently normalized closed loops as spatial_winding.rs:
// J0 = 2 curl(N e_z), J1 = 2 curl(N(e_x+e_z)), with the central P1 hat N.
// N vanishes on every exterior face, so current closes even on prescribed A.
fn profiles(p: &SpatialNonlinearVectorMagnetostatics) -> Vec<Vec<[f64; 3]>> {
    let origin = p.positions_m().iter().position(|&v| v == [0.5; 3]).unwrap() as u32;
    let mut result = vec![Vec::new(), Vec::new()];
    for tet in &p.mesh().tets {
        let gradient = if tet.contains(&origin) {
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != origin).collect();
            let x = p.positions_m()[face[0] as usize];
            let normal = cross(
                sub(p.positions_m()[face[1] as usize], x),
                sub(p.positions_m()[face[2] as usize], x),
            );
            let height = dot(normal, sub(p.positions_m()[origin as usize], x));
            normal.map(|v| v / height)
        } else { [0.0; 3] };
        result[0].push(cross(gradient, [0.0, 0.0, 2.0]));
        result[1].push(cross(gradient, [2.0, 0.0, 2.0]));
    }
    result
}
fn uniform(p: &SpatialNonlinearVectorMagnetostatics, b: [f64; 3], shifted: bool) -> Vec<f64> {
    p.mesh().edges.iter().map(|&[i, j]| {
        let (x, y) = (p.positions_m()[i as usize], p.positions_m()[j as usize]);
        let mid = std::array::from_fn(|k| 0.5 * (x[k] + y[k]));
        let psi = |v: u32| if shifted { f64::from((v * 7) % 11) * 0.17 } else { 0.0 };
        0.5 * dot(cross(b, mid), sub(y, x)) + psi(j) - psi(i)
    }).collect()
}
fn boundary(p: &SpatialNonlinearVectorMagnetostatics, a: &[f64]) -> VectorMagneticBoundary {
    VectorMagneticBoundary {
        potential: p.exterior_faces().into_iter().map(|face| VectorPotentialFace {
            face, edge_integrals_wb: p.face_edges(face).unwrap().map(|e| a[e]),
        }).collect(),
    }
}
fn basis<'a>(
    p: &'a SpatialNonlinearVectorMagnetostatics,
    profiles: Vec<Vec<[f64; 3]>>,
    initial: &[f64],
    currents: &[f64],
    ends: &[f64],
) -> SpatialIncrementalWindingBasis<'a> {
    SpatialIncrementalWindingBasis::new(
        p, profiles, boundary(p, initial), currents.to_vec(), ends.to_vec(),
        initial, config(), controls(), &mut || false,
    ).unwrap()
}
fn sum_current(profiles: &[Vec<[f64; 3]>], currents: &[f64]) -> Vec<[f64; 3]> {
    (0..profiles[0].len()).map(|cell| {
        std::array::from_fn(|k| profiles.iter().zip(currents)
            .map(|(p, i)| p[cell][k] * i).sum())
    }).collect()
}

#[test]
fn linear_limit_matches_native_3d_basis_and_keeps_bias_out_of_the_derivative() {
    let p = problem(SpatialNonlinearMaterial::Linear(material(2.0)));
    let profiles = profiles(&p);
    let zero = vec![0.0; p.mesh().edges.len()];
    let ends = [0.02, 0.03];
    let native = SpatialVectorMagnetostatics::new(
        p.positions_m().to_vec(), p.mesh().tets.clone(),
        vec![SpatialMagneticMaterial::Linear(material(2.0)); p.mesh().tets.len()],
        common::limits(), &mut || false,
    ).unwrap();
    let original = SpatialWindingBasis::new(
        &native, profiles.clone(), boundary(&p, &zero), ends.to_vec(),
        config().winding, controls().field, &mut || false,
    ).unwrap();
    let model = basis(&p, profiles.clone(), &zero, &[0.7, -0.2], &ends);
    for (&actual, &expected) in model.inductance_h().iter().zip(original.inductance_h()) {
        near(actual, expected);
    }
    for (&actual, expected) in model.operating_flux_linkages_wb().iter()
        .zip(original.flux_linkages_wb(&[0.7, -0.2]).unwrap())
    {
        near(actual, expected);
    }
    for response in model.unit_responses() {
        assert!(response.report().unwrap().converged_euclidean());
        assert!(response.maximum_equation_residual <= controls().field.current_tolerance_a);
        assert!(response.source_compatibility_defect <= controls().field.current_tolerance_a);
        assert!(response.maximum_flux_defect_h <= controls().field.flux_tolerance_wb);
        for (&fixed, &a) in model.operating_field().prescribed_edges().iter()
            .zip(response.potential_per_current_h())
        {
            if fixed { assert_eq!(a, 0.0); }
        }
    }
    let initial = uniform(&p, [0.3, 0.2, 0.4], false);
    let biased = basis(&p, profiles.clone(), &initial, &[0.0, 0.0], &ends);
    assert!(biased.operating_flux_linkages_wb().iter().any(|v| v.abs() > 1e-3));
    assert_eq!(biased.voltage_from_current_rate_v(&[0.0, 0.0]).unwrap(), vec![0.0; 2]);
    let shifted = basis(&p, profiles, &uniform(&p, [0.3, 0.2, 0.4], true), &[0.0, 0.0], &ends);
    for ((&actual, &base), &gauge) in biased.inductance_h().iter()
        .zip(model.inductance_h()).zip(shifted.inductance_h())
    {
        near(actual, base);
        near(actual, gauge);
    }
    for (&actual, &gauge) in biased.operating_flux_linkages_wb().iter()
        .zip(shifted.operating_flux_linkages_wb())
    {
        near(actual, gauge);
    }
}

#[test]
fn saturated_tangent_retains_both_longitudinal_and_transverse_reluctivity() {
    let p = problem(SpatialNonlinearMaterial::Reversible(law()));
    let profiles = profiles(&p);
    let zero = vec![0.0; p.mesh().edges.len()];
    let unbiased = basis(&p, profiles.clone(), &zero, &[0.0, 0.0], &[0.0, 0.0]);
    let initial = uniform(&p, [0.0, 0.0, 1.5], false);
    let biased = basis(&p, profiles, &initial, &[0.0, 0.0], &[0.0, 0.0]);
    // At |B|=1.5 T, H=3 A/m: transverse reluctivity is H/B=2,
    // longitudinal reluctivity is dH/dB=4. A secant scalar loses this tensor.
    for i in 0..2 {
        for j in 0..2 {
            let a = biased.unit_responses()[i].flux_density_per_current_t_a();
            let b = biased.unit_responses()[j].flux_density_per_current_t_a();
            let gram: f64 = a.iter().zip(b).enumerate().map(|(cell, (a, b))| {
                p.cell_volume_m3(cell).unwrap()
                    * (2.0 * a[0] * b[0] + 2.0 * a[1] * b[1] + 4.0 * a[2] * b[2])
            }).sum();
            near(biased.active_inductance_h()[i * 2 + j], gram);
        }
        let l0 = unbiased.active_inductance_h()[i * 2 + i];
        let l = biased.active_inductance_h()[i * 2 + i];
        assert!(l >= 0.25 * l0 - 2e-9 && l <= 0.5 * l0 + 2e-9);
        assert!((l - 0.5 * l0).abs() > 1e-3 * l0, "secant permeability erased saturation");
    }
    assert_eq!(biased.inductance_h()[1].to_bits(), biased.inductance_h()[2].to_bits());
    let delta = biased.linearized_flux_increment_wb(&[0.7, -0.2]).unwrap();
    assert!(0.7 * delta[0] - 0.2 * delta[1] > 0.0);
}

#[test]
fn actual_biased_flux_finite_differences_agree_with_all_signed_mutual_entries() {
    let p = problem(SpatialNonlinearMaterial::Reversible(law()));
    let profiles = profiles(&p);
    let initial = uniform(&p, [0.0, 0.0, 1.5], false);
    let currents = [0.2, -0.1];
    let ends = [0.02, 0.03];
    let model = basis(&p, profiles.clone(), &initial, &currents, &ends);
    let eps = 1e-4;
    for j in 0..2 {
        let mut plus = currents;
        let mut minus = currents;
        plus[j] += eps;
        minus[j] -= eps;
        let fp = p.solve(
            &sum_current(&profiles, &plus), &boundary(&p, &initial),
            model.operating_field().potential_wb(), controls(), &mut || false,
        ).unwrap();
        let fm = p.solve(
            &sum_current(&profiles, &minus), &boundary(&p, &initial),
            model.operating_field().potential_wb(), controls(), &mut || false,
        ).unwrap();
        let lp = model.pair_potential_wb(fp.potential_wb(), &mut || false).unwrap();
        let lm = model.pair_potential_wb(fm.potential_wb(), &mut || false).unwrap();
        for i in 0..2 {
            let finite_difference = (lp[i] - lm[i] + ends[i] * (plus[i] - minus[i])) / (2.0 * eps);
            let tangent = model.inductance_h()[i * 2 + j];
            assert!((finite_difference - tangent).abs() < 2e-5, "{finite_difference} != {tangent}");
        }
    }
    let mut reversed = profiles;
    for vector in &mut reversed[1] {
        *vector = vector.map(|v| -v);
    }
    let other = basis(&p, reversed, &initial, &[currents[0], -currents[1]], &ends);
    for i in 0..2 {
        let expected = if i == 1 {
            -model.operating_flux_linkages_wb()[i]
        } else {
            model.operating_flux_linkages_wb()[i]
        };
        near(other.operating_flux_linkages_wb()[i], expected);
        for j in 0..2 {
            let sign = if i == j { 1.0 } else { -1.0 };
            near(other.inductance_h()[i * 2 + j], sign * model.inductance_h()[i * 2 + j]);
        }
    }
}

#[test]
fn source_corners_and_envelope_edges_refuse_but_collinear_knots_remain_smooth() {
    let p = problem(SpatialNonlinearMaterial::Reversible(law()));
    for (b, corner) in [(1.0, true), (3.999, false)] {
        let initial = uniform(&p, [0.0, 0.0, b], false);
        let mut cfg = config();
        cfg.minimum_smooth_field_margin_t = 0.002;
        let result = SpatialIncrementalWindingBasis::new(
            &p, profiles(&p), boundary(&p, &initial), vec![0.0; 2], vec![0.0; 2],
            &initial, cfg, controls(), &mut || false,
        );
        if corner {
            assert!(matches!(result, Err(SpatialIncrementalWindingError::Nondifferentiable { .. })));
        } else {
            assert!(matches!(result, Err(SpatialIncrementalWindingError::SourceBoundary { .. })));
        }
    }
    let collinear = magnetic::law(&magnetic::source(
        &[(0.0, 0.0), (1.0, 1.0), (4.0, 4.0)], 2.0,
    ));
    let p = problem(SpatialNonlinearMaterial::Reversible(collinear));
    let initial = uniform(&p, [0.0, 0.0, 1.0], false);
    let model = basis(&p, profiles(&p), &initial, &[0.0; 2], &[0.0; 2]);
    assert!(model.inductance_h()[0] > 0.0);
}

#[test]
fn dependent_observation_modes_keep_physical_zero_leakage() {
    let p = problem(SpatialNonlinearMaterial::Linear(material(2.0)));
    let primary = profiles(&p).remove(0);
    let secondary: Vec<_> = primary.iter().map(|v| v.map(|x| -3.0 * x)).collect();
    let model = basis(&p, vec![primary, secondary], &vec![0.0; p.mesh().edges.len()],
        &[0.2, -0.1], &[0.0, 0.0]);
    let l = model.inductance_h();
    near(l[1], -3.0 * l[0]);
    near(l[2], -3.0 * l[0]);
    near(l[3], 9.0 * l[0]);
    for flux in model.linearized_flux_increment_wb(&[3.0, 1.0]).unwrap() {
        near(flux, 0.0);
    }
    assert_eq!(model.external_inductance_h(), &[0.0; 2]);
}

#[test]
fn bounds_current_closure_cancellation_and_solver_failure_publish_no_basis() {
    let p = problem(SpatialNonlinearMaterial::Linear(material(2.0)));
    let initial = vec![0.0; p.mesh().edges.len()];
    let run = |profiles, cfg, ctrl, cancelled: &mut dyn FnMut() -> bool| {
        SpatialIncrementalWindingBasis::new(
            &p, profiles, boundary(&p, &initial), vec![0.0; 2], vec![0.0; 2],
            &initial, cfg, ctrl, &mut || cancelled(),
        )
    };
    let mut cfg = config();
    cfg.winding.max_basis_values = 0;
    assert!(run(profiles(&p), cfg, controls(), &mut || false).is_err());
    let mut bad = profiles(&p);
    bad[0][0][0] += 1.0; // Operating current zero must not hide a broken unit profile.
    assert!(run(bad, config(), controls(), &mut || false).is_err());
    let mut ctrl = controls();
    ctrl.field.field.max_iterations = 0;
    ctrl.field.field.max_corrections = 0;
    assert!(run(profiles(&p), config(), ctrl, &mut || false).is_err());
    let mut total_polls = 0;
    let model = run(profiles(&p), config(), controls(), &mut || {
        total_polls += 1;
        false
    }).unwrap();
    for limit in [0, total_polls / 2, total_polls - 1] {
        let mut polls = 0;
        assert!(matches!(run(profiles(&p), config(), controls(), &mut || {
            polls += 1;
            polls > limit
        }), Err(SpatialIncrementalWindingError::Magnetic(
            NonlinearVectorError::Spatial(SpatialMagneticError::Cancelled)
        ))));
    }
    assert!(model.linearized_flux_increment_wb(&[0.0]).is_err());
    assert!(model.voltage_from_current_rate_v(&[f64::NAN, 0.0]).is_err());
}

fn recoil(mu: f64, br: f64, limit: f64) -> SpatialRecoilMaterial {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::{
        field_storage::{
            MAGNETIC_PERMEABILITY_DIMS,
            permanent::{
                FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS, RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY,
                RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY,
                RECOIL_PERMEABILITY_PROPERTY, REMANENT_FLUX_DENSITY_PROPERTY, ResolvedRecoilMaterial,
            },
        },
        state_point::MaterialPropertySelection,
    };
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (RECOIL_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS, mu),
        (REMANENT_FLUX_DENSITY_PROPERTY, FLUX_DENSITY_DIMS, br),
        (RECOIL_MINIMUM_PARALLEL_FIELD_PROPERTY, FIELD_STRENGTH_DIMS, -limit),
        (RECOIL_MAXIMUM_PARALLEL_FIELD_PROPERTY, FIELD_STRENGTH_DIMS, limit),
        (RECOIL_MAXIMUM_TRANSVERSE_FIELD_PROPERTY, FIELD_STRENGTH_DIMS, limit),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured differential recoil source".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic magnet".into(), phase: "solid".into(),
            process: "reversible test".into(), revision: 0,
        },
        claims,
        vec![],
    ).unwrap();
    let material = ResolvedRecoilMaterial::resolve(
        &card, &QueryPoint::new(), MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    SpatialRecoilMaterial::new(material, [0.0, 0.0, 1.0]).unwrap()
}

#[test]
fn recoil_bias_is_physical_but_unit_derivatives_need_not_fit_a_finite_amplitude_envelope() {
    let mu = 2.0;
    let envelope = 1e-4;
    let p = problem(SpatialNonlinearMaterial::Recoil(recoil(mu, 0.5, envelope)));
    let initial = uniform(&p, [0.0, 0.0, 0.5], false);
    let model = basis(&p, profiles(&p), &initial, &[0.0; 2], &[0.02, 0.03]);
    assert!(model.operating_flux_linkages_wb().iter().any(|v| v.abs() > 1e-3));
    for &h in model.operating_field().field_strength_a_m() {
        near(dot(h, h), 0.0);
    }
    assert!(model.unit_responses()[0].flux_density_per_current_t_a().iter()
        .any(|&b| dot(b, b).sqrt() / mu > envelope));
    let linear = problem(SpatialNonlinearMaterial::Linear(material(mu)));
    let reference = basis(
        &linear, profiles(&linear), &vec![0.0; linear.mesh().edges.len()],
        &[0.0; 2], &[0.02, 0.03],
    );
    for (&l, &expected) in model.inductance_h().iter().zip(reference.inductance_h()) {
        near(l, expected);
    }
    assert_eq!(model.voltage_from_current_rate_v(&[0.0; 2]).unwrap(), vec![0.0; 2]);
}
