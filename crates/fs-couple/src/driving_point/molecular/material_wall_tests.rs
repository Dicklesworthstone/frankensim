//! G1/G3: molecular caloric storage and source-wall motion reach the same port.
//! Synthetic material claims verify causality; they are not measured air/wall data.
use super::*;
use crate::bernoulli_aperture::wall::material::{
    AxialConstraint, WallMaterialLoss, relaxation::WallRelaxation,
};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::gas::GasSpec;
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use fs_thermochem::SpeciesId;

const RADIUS: f64 = 0.0073;
const LENGTH: f64 = 0.57;
const TAU: f64 = core::f64::consts::TAU;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn molecular(time_scale: f64) -> RelaxationState {
    let data = format!(
        "declared 1 synthetic 2 vibration 1 1300 1 power {} 293.15 0 200 500 \
         vibration 1 1900 1 power {} 293.15 0 200 500",
        101325.0 * 2e-5 * time_scale,
        101325.0 * 3.2e-4 * time_scale,
    );
    read_relaxation(&mut data.split_ascii_whitespace())
        .unwrap()
        .evaluate(
            &RelaxingMedium::new(gas(), vec![(SpeciesId::new("synthetic").unwrap(), 1.0)]).unwrap(),
        )
        .unwrap()
}

fn wall(time_s: f64, band: [f64; 2]) -> LocalMaterialWall {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, 3e7),
        ("poisson_ratio", Dims([0; 6]), 0.3),
        ("wall_delta_e", Pressure::DIMS, 2e7),
        ("wall_tau", Dims([0, 0, 1, 0, 0, 0]), time_s),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 290.0, 300.0)
                    .with("omega", TAU * band[0], TAU * band[1]),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic Maxwell wall for molecular port integration".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_with_quantity(
                key,
                QuantitySpec::dimensional(dims),
                ScalarAdmissibility::Finite,
            )
            .unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-wall".into(),
            phase: "solid".into(),
            process: "manufactured test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let spectrum = WallRelaxation::from_state_point(
        &material,
        &[("wall_delta_e".into(), "wall_tau".into())],
        band,
        0.5,
        0.5,
    )
    .unwrap();
    LocalMaterialWall::new(
        material,
        0.0002,
        AxialConstraint::StressFree,
        WallMaterialLoss::Maxwell(spectrum),
        band,
    )
    .unwrap()
}

fn geometry() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: RADIUS,
            length: LENGTH,
        }],
    }
}

fn point(gas_time: f64, wall_time: f64, diffusion: bool) -> MolecularDrivingPoint {
    MolecularDrivingPoint::new_with_material_wall(
        &geometry(),
        molecular(gas_time),
        diffusion.then_some(MolecularWallResolution {
            poles_per_field: 32,
            band_rad_s: [1e-4, 30000.0],
            maximum_relative_error: 0.05,
        }),
        wall(wall_time, [1e-5, 4500.0]),
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        gas(),
    )
    .unwrap()
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs().max(1e-20),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn molecular_material_wall_matches_independent_caloric_and_membrane_cylinder() {
    let point = point(1.0, 0.0003, false);
    let material = point.material().state();
    let surface = point.material_wall().unwrap();
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let z0 = gas().density * gas().sound_speed / area;
    for hz in [30.0, 100.0, 430.0, 700.0, 2000.0] {
        let omega = TAU * hz;
        // Independent Debye sum and physical pressure/volume-flow telegraph law.
        // Neither MolecularDuct's wave nor its scattering composition is used.
        let cv = material.processes().iter().fold(
            C64::from_re(material.frozen_cv_j_kg_k()),
            |sum, arm| {
                sum + C64::from_re(arm.heat_capacity_j_kg_k)
                    / C64::new(1.0, -omega * arm.relaxation_time_s)
            },
        );
        let cp = cv + C64::from_re(gas().specific_gas_constant);
        let gas_compliance = (cv / cp).scale(area / gas().pressure);
        let shunt = C64::new(0.0, -omega) * gas_compliance
            + surface
                .specific_impedance(RADIUS, omega)
                .unwrap()
                .recip()
                .scale(TAU * RADIUS);
        let series = C64::new(0.0, -omega * gas().density / area);
        let k = (-series * shunt).sqrt();
        let z = C64::from_re(omega * gas().density / area) / k;
        let phase = k.scale(LENGTH);
        let cosine = C64::new(
            phase.re.cos() * phase.im.cosh(),
            -phase.re.sin() * phase.im.sinh(),
        );
        let sine = C64::new(
            phase.re.sin() * phase.im.cosh(),
            phase.re.cos() * phase.im.sinh(),
        );
        let load = fs_duct::compact_radiation_impedance(
            Termination::FlangedOpen,
            point.exterior_gas(),
            RADIUS,
            omega,
        )
        .unwrap();
        let p = cosine * load + C64::new(0.0, -1.0) * z * sine;
        let u = cosine + C64::new(0.0, -1.0) * sine * load / z;
        close(
            point
                .reflectance(omega, Termination::FlangedOpen, z0)
                .unwrap(),
            (p - u.scale(z0)) / (p + u.scale(z0)),
        );
        close(
            point.flanged_outlet_flow(omega).unwrap(),
            C64::from_re(2.0) / (p + u.scale(z0)),
        );
    }
    assert_eq!(
        point.one_way_wavefront_s(),
        LENGTH / material.frozen_sound_speed_m_s()
    );
    assert_eq!(
        surface.material().identity(),
        point.clone().material_wall().unwrap().material().identity()
    );
}

#[test]
fn gas_and_wall_source_changes_both_reach_feedback_and_observed_pressure() {
    let baseline = point(1.0, 0.0003, true);
    let gas_changed = point(0.1, 0.0003, true);
    let wall_changed = point(1.0, 0.003, true);
    let z0 = gas().density * gas().sound_speed / (core::f64::consts::PI * RADIUS * RADIUS);
    let r = baseline
        .reflectance(TAU * 430.0, Termination::FlangedOpen, z0)
        .unwrap();
    for changed in [&gas_changed, &wall_changed] {
        let new = changed
            .reflectance(TAU * 430.0, Termination::FlangedOpen, z0)
            .unwrap();
        assert!(r.abs() <= 1.0 + 1e-12 && new.abs() <= 1.0 + 1e-12);
        assert!((new - r).abs() > 1e-4);
    }
    let receiver = CircularOutletReceiver {
        position_m: [0.0, 0.0, 1.0],
        radial_rings: 2,
        angular_points: 8,
        maximum_frequency_hz: 800.0,
    };
    let mut observers = [&baseline, &gas_changed, &wall_changed]
        .map(|p| p.baffled_outlet(8000, 512, receiver).unwrap());
    let first =
        observers[0].info().axial_delay_samples + observers[0].info().exterior_delay_samples.0;
    let mut energy = 0.0;
    let mut differences = [0.0; 2];
    for n in 0..1024 {
        let input = if n == 0 { 100.0 } else { 0.0 };
        let p = observers
            .each_mut()
            .map(|observer| observer.step(input).unwrap());
        assert!(p.iter().all(|value| value.is_finite()));
        if n < first {
            assert_eq!(p, [0.0; 3]);
        }
        energy += p[0] * p[0];
        for (difference, changed) in differences.iter_mut().zip(&p[1..]) {
            *difference += (changed - p[0]).powi(2);
        }
    }
    assert!(
        energy > 1e-12 && differences.iter().all(|d| *d > 1e-7 * energy),
        "energy={energy}, differences={differences:?}"
    );
}

#[test]
fn molecular_feedback_and_observation_refuse_missing_wall_source_coverage() {
    let point = MolecularDrivingPoint::new_with_material_wall(
        &geometry(),
        molecular(1.0),
        None,
        wall(0.0003, [1.0, 1000.0]),
        ConeSlices::Fixed(8),
        ComposeOptions::default(),
        gas(),
    )
    .unwrap();
    let z0 = gas().characteristic_impedance / (core::f64::consts::PI * RADIUS * RADIUS);
    assert!(
        point
            .reflectance(TAU * 430.0, Termination::FlangedOpen, z0)
            .is_ok()
    );
    assert!(
        point
            .reflectance(1e-3, Termination::FlangedOpen, z0)
            .is_err()
    );
    assert!(
        point
            .reflectance(TAU * 1200.0, Termination::FlangedOpen, z0)
            .is_err()
    );
    assert!(
        point
            .realize_line(Termination::FlangedOpen, 8000, z0, &LineLadder::default())
            .is_err()
    );
    assert!(
        point
            .baffled_outlet(
                8000,
                512,
                CircularOutletReceiver {
                    position_m: [0.0, 0.0, 1.0],
                    radial_rings: 2,
                    angular_points: 8,
                    maximum_frequency_hz: 800.0,
                }
            )
            .is_err()
    );
}
