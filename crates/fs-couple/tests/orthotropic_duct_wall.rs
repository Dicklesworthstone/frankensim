//! Manufactured principal-plane data, not measured wood or instrument fidelity.
use fs_couple::bernoulli_aperture::wall::material::{AxialConstraint, LocalMaterialWall, WallMaterialLoss};
use fs_duct::dg::{DgBoundary, DgDuct, DgMesh};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::gas::{GasSpec, GasState};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, DynViscosity, Pressure};

const BAND: [f64; 2] = [20.0, 4000.0];
const THICKNESS: f64 = 0.0002;
const RADIUS: f64 = 0.006;
const TIME: &str = "normal_plane_viscosity_time";

type Datum<'a> = (&'a str, Dims, f64);
fn data(e: [f64; 3], nu: [f64; 3], tau: f64) -> Vec<Datum<'static>> {
    vec![
        ("density", Density::DIMS, 1000.0),
        ("young_modulus_1", Pressure::DIMS, e[0]),
        ("young_modulus_2", Pressure::DIMS, e[1]),
        ("young_modulus_3", Pressure::DIMS, e[2]),
        ("poisson_ratio_12", Dims::NONE, nu[0]),
        ("poisson_ratio_13", Dims::NONE, nu[1]),
        ("poisson_ratio_23", Dims::NONE, nu[2]),
        (TIME, Dims([0, 0, 1, 0, 0, 0]), tau),
    ]
}
fn state(data: &[Datum<'_>], frequency_band: [f64; 2]) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in data {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0)
                .with("omega", frequency_band[0], frequency_band[1]),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance {
                source: "manufactured principal-plane membrane test, not measured data".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic-directional-wall".into(), phase: "solid".into(),
        process: "manufactured test".into(), revision: 0,
    }, claims, vec![]).unwrap();
    resolve_material_state_point(&card,
        &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn source(e: [f64; 3], nu: [f64; 3], tau: f64) -> ResolvedMaterialStatePoint {
    state(&data(e, nu, tau), [0.0, 1e6])
}
fn wall(source: ResolvedMaterialStatePoint, axes: [usize; 2], axial: AxialConstraint, viscous: bool) -> LocalMaterialWall {
    LocalMaterialWall::new_orthotropic(source, THICKNESS, axial, axes, viscous.then_some(TIME), BAND).unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 2e-13 * expected.abs().max(1e-20), "{actual:e} != {expected:e}");
}

#[test]
fn every_principal_plane_reuses_the_reciprocal_normal_compliance() {
    let e = [12e6, 3e6, 1.5e6];
    let nu = [0.25, -0.2, 0.3];
    let original = source(e, nu, 0.0);
    for axes in [[0, 1], [1, 0], [0, 2], [2, 0], [1, 2], [2, 1]] {
        let [h, a] = axes;
        let (i, j) = (h.min(a), h.max(a));
        let n = nu[match (i, j) { (0, 1) => 0, (0, 2) => 1, _ => 2 }];
        // Independent dimensional compliance and Schur complement. Reversing
        // the axes changes E_hoop, not the off-diagonal reciprocal compliance.
        let s_hh = 1.0 / e[h];
        let s_aa = 1.0 / e[a];
        let s_ha = -n / e[i];
        for constraint in [AxialConstraint::StressFree, AxialConstraint::StrainFixed] {
            let bound = wall(original.clone(), axes, constraint, false);
            let pin = bound.base_patch(RADIUS, 0.3).unwrap().wall;
            let modulus = match constraint {
                AxialConstraint::StressFree => s_hh.recip(),
                AxialConstraint::StrainFixed => (s_hh - s_ha * s_ha / s_aa).recip(),
            };
            close(pin.stiffness_per_area, modulus * THICKNESS / RADIUS.powi(2));
            close(pin.surface_density, 1000.0 * THICKNESS);
            assert_eq!(pin.resistance, 0.0);
            assert_eq!(bound.orthotropic_axes(), Some(axes));
            assert_eq!(bound.material(), &original);
        }
    }
    // A source Poisson ratio is not subject to the isotropic -1..0.5 bound.
    let anisotropic = source([16e6, 1e6, 1e6], [2.0, 0.0, 0.0], 0.0);
    close(wall(anisotropic, [0, 1], AxialConstraint::StrainFixed, false)
        .base_patch(RADIUS, 0.3).unwrap().wall.stiffness_per_area,
        (16e6 / 0.75) * THICKNESS / RADIUS.powi(2));
}

#[test]
fn isotropic_limit_and_proportional_viscosity_match_the_original_wall_owner() {
    let young = 3e6;
    let tau = 2e-5;
    let mut values = data([young; 3], [0.3; 3], tau);
    values.extend([
        ("young_modulus", Pressure::DIMS, young),
        ("poisson_ratio", Dims::NONE, 0.3),
        ("kelvin_voigt_young_viscosity", DynViscosity::DIMS, young * tau),
    ]);
    let state = state(&values, [0.0, 1e6]);
    for axial in [AxialConstraint::StressFree, AxialConstraint::StrainFixed] {
        let original = LocalMaterialWall::new(state.clone(), THICKNESS, axial, WallMaterialLoss::KelvinVoigt, BAND).unwrap();
        assert_eq!(original.orthotropic_axes(), None);
        for axes in [[0, 1], [1, 0], [2, 1]] {
            let orthotropic = wall(state.clone(), axes, axial, true);
            let a = original.base_patch(RADIUS, 0.3).unwrap();
            let b = orthotropic.base_patch(RADIUS, 0.3).unwrap();
            close(a.wall.stiffness_per_area, b.wall.stiffness_per_area);
            close(a.wall.resistance, b.wall.resistance);
            close(b.wall.resistance, tau * b.wall.stiffness_per_area);
            assert_eq!(a.wall.surface_density, b.wall.surface_density);
            assert_eq!(orthotropic.proportional_viscosity_time_property(), Some(TIME));
            assert_eq!(orthotropic.material(), &state);
        }
    }
}

#[test]
fn directional_requirements_do_not_invent_unused_radial_or_shear_data() {
    let needed = LocalMaterialWall::orthotropic_requirements([2, 0], Some(TIME)).unwrap();
    let names: Vec<_> = needed.iter().map(|p| p.name()).collect();
    assert_eq!(names, ["density", "young_modulus_3", "young_modulus_1", "poisson_ratio_13", TIME]);
    let values: Vec<_> = data([12e6, 3e6, 1.5e6], [0.25, 0.2, 0.3], 0.0)
        .into_iter().filter(|(name, _, _)| names.contains(name)).collect();
    let exact = state(&values, [0.0, 1e6]);
    assert!(LocalMaterialWall::new_orthotropic(exact.clone(), THICKNESS,
        AxialConstraint::StressFree, [2, 0], Some(TIME), BAND).is_ok());
    assert!(LocalMaterialWall::new_orthotropic(exact, THICKNESS,
        AxialConstraint::StressFree, [1, 0], None, BAND).is_err());
}

#[test]
fn incomplete_unstable_out_of_band_and_wrong_quantity_sources_refuse() {
    let baseline = data([4e6, 1e6, 1e6], [0.2; 3], 1e-5);
    for name in ["density", "young_modulus_1", "young_modulus_2", "poisson_ratio_12", TIME] {
        let partial: Vec<_> = baseline.iter().copied().filter(|(key, _, _)| *key != name).collect();
        assert!(LocalMaterialWall::new_orthotropic(state(&partial, [0.0, 1e6]), THICKNESS,
            AxialConstraint::StrainFixed, [0, 1], Some(TIME), BAND).is_err());
    }
    for (name, dims, value) in [
        ("density", Density::DIMS, 0.0),
        ("young_modulus_1", Pressure::DIMS, -1.0),
        ("poisson_ratio_12", Dims::NONE, 2.0), // singular at E2/E1=1/4
        ("poisson_ratio_12", Dims::NONE, -2.1),
        (TIME, Dims([0, 0, 1, 0, 0, 0]), -1e-5),
        (TIME, DynViscosity::DIMS, 1e-5),
    ] {
        let mut changed = baseline.clone();
        *changed.iter_mut().find(|p| p.0 == name).unwrap() = (name, dims, value);
        assert!(LocalMaterialWall::new_orthotropic(state(&changed, [0.0, 1e6]), THICKNESS,
            AxialConstraint::StrainFixed, [0, 1], Some(TIME), BAND).is_err());
    }
    let original = source([4e6, 1e6, 1e6], [0.2; 3], 1e-5);
    for axes in [[0, 0], [3, 0], [0, usize::MAX]] {
        assert!(LocalMaterialWall::orthotropic_requirements(axes, None).is_err());
        assert!(LocalMaterialWall::new_orthotropic(original.clone(), THICKNESS,
            AxialConstraint::StressFree, axes, None, BAND).is_err());
    }
    assert!(LocalMaterialWall::new_orthotropic(state(&baseline, [0.0, 2000.0]), THICKNESS,
        AxialConstraint::StressFree, [0, 1], None, BAND).is_err());
    for thickness in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(LocalMaterialWall::new_orthotropic(original.clone(), thickness,
            AxialConstraint::StressFree, [0, 1], None, BAND).is_err());
    }
    assert!(wall(original, [0, 1], AxialConstraint::StressFree, false)
        .base_patch(THICKNESS, 0.3).is_err());
}

#[test]
fn directional_source_and_gas_changes_reach_the_original_dg_pressure_and_energy() {
    let source = source([4e6, 1e6, 1e6], [0.2; 3], 0.0);
    let mut histories = Vec::new();
    for (axes, pressure) in [([0, 1], 101325.0), ([1, 0], 101325.0), ([0, 1], 80000.0)] {
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap();
        let geometry = fs_duct::Duct { segments: vec![fs_duct::Segment::Cylinder { radius: RADIUS, length: 0.3 }] };
        let duct = DgDuct::new(&geometry, &gas, DgMesh { elements: 2, order: 2 },
            DgBoundary::Closed, DgBoundary::Closed).unwrap();
        let bound = wall(source.clone(), axes, AxialConstraint::StrainFixed, false).bind_dg(duct).unwrap();
        let duct = bound.duct();
        let pin = duct.compliant_wall().unwrap();
        // Uniform closed-gas compression: p=1-beta*w; the material wall gives
        // sigma*w'' + (K+beta)*w=1. This is independent of DG implementation.
        let beta = 2.0 * gas.density * gas.sound_speed.powi(2) / RADIUS;
        let omega = ((pin.stiffness_per_area + beta) / pin.surface_density).sqrt();
        let dt = 1e-6;
        let propagator = duct.propagator(dt, 1e-13).unwrap();
        let mut current = duct.project_initial_state(|_| 1.0, |_| 0.0);
        let initial = current.stored_energy_j();
        let mut pressures = Vec::new();
        for step in 1..=100 {
            propagator.step(&mut current, &[]).unwrap();
            propagator.sync_nodal(&mut current).unwrap();
            let displacement = (1.0 - (omega * step as f64 * dt).cos()) / (pin.stiffness_per_area + beta);
            assert!((current.pressure()[0] - (1.0 - beta * displacement)).abs() < 1e-10);
            assert!((duct.compliant_wall_state(&current).unwrap().displacement_m[0] - displacement).abs() < 1e-17);
            assert!((current.stored_energy_j() - initial).abs() < 1e-10 * initial);
            pressures.push(current.pressure()[0]);
        }
        histories.push(pressures);
    }
    assert_ne!(histories[0], histories[1]);
    assert_ne!(histories[0], histories[2]);
    assert!(histories.iter().flatten().all(|p| p.is_finite()));
}
