//! Synthetic source cards test material propagation, not measured duct fidelity.
use super::*;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::bernoulli_aperture::wall::material::{
    AxialConstraint, MaterialWallAssignment, WallMaterialLoss, bind_material_walls,
    relaxation::WallRelaxation,
};
use crate::driving_point::regional::field::tests::{cylinder, field};
use fs_duct::{Duct, RadiusProfile, RiccatiSteps, Segment};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::{MaterialPropertySelection, ResolvedMaterialStatePoint,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

const BAND: [f64; 2] = [1.0, 4500.0];
const TAU: f64 = core::f64::consts::TAU;
const RADIUS: f64 = 0.0073;
const LENGTH: f64 = 0.57;
const THICKNESS: f64 = 0.0002;

// Changing a relaxation time changes only a source claim. Geometry, gas,
// equilibrium elasticity, density and every acoustic algorithm stay fixed.
fn material(time_s: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (key, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, 3e7),
        ("poisson_ratio", Dims([0; 6]), 0.3),
        ("wall_delta_e", Pressure::DIMS, 2e7),
        ("wall_tau", Dims([0, 0, 1, 0, 0, 0]), time_s),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0)
                .with("omega", TAU * BAND[0], TAU * BAND[1]),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic isotropic Maxwell wall; not measured material data".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_with_quantity(key,
            QuantitySpec::dimensional(dims), ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic-wall".into(), phase: "solid".into(),
        process: "manufactured test".into(), revision: 0,
    }, claims, vec![]).unwrap();
    resolve_material_state_point(&card,
        &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn wall(time_s: f64) -> LocalMaterialWall {
    let material = material(time_s);
    let loss = WallRelaxation::from_state_point(&material,
        &[("wall_delta_e".into(), "wall_tau".into())], BAND, 0.5, 0.5).unwrap();
    LocalMaterialWall::new(material, THICKNESS, AxialConstraint::StressFree,
        WallMaterialLoss::Maxwell(loss), BAND).unwrap()
}
fn receiver() -> CircularOutletReceiver {
    CircularOutletReceiver { position_m: [0.0, 0.0, 1.0], radial_rings: 2,
        // The single-rate circular receiver admits up to 0.1 x the 8 kHz rate.
        angular_points: 8, maximum_frequency_hz: 800.0 }
}
fn relative(actual: C64, expected: C64) -> f64 {
    (actual - expected).abs() / expected.abs().max(1e-30)
}

#[test]
fn tabulated_material_wall_refuses_when_static_wall_data_is_unavailable() {
    let gas = *field(300.0).state_at(0.0).unwrap().state();
    let bound = RegionalDuct::new(&cylinder(), &[gas]).unwrap().with_material_wall(wall(0.0003));
    // Every positive 8 kHz / 512-point FFT bin is supported by this wall.
    for omega in [TAU * 8000.0 / 512.0, TAU * 4000.0] {
        assert!(bound.wall().specific_impedance(RADIUS, omega).is_ok());
    }
    assert!(bound.characteristic_line(Termination::FlangedOpen, 8000, 4801, Some(512)).is_ok());
    let table = fs_duct::TabulatedLoad::try_new(
        vec![TAU * 10.0, TAU * 4000.0], vec![C64::ZERO; 2],
        "ideal-open static-domain control", 0.0,
    ).unwrap().with_mouth(fs_duct::TabulatedMouth {
        axial_position_m: LENGTH, radius_m: RADIUS, gas, dc_impedance: 0.0,
    }).unwrap();
    // Mouth metadata supplies only the terminal DC load, not a source-backed
    // static law for the compliant wall. Do not silently use rigid-wall DC.
    assert!(bound.wall().specific_impedance(RADIUS, 1e-3).is_err());
    assert!(matches!(bound.characteristic_line(MouthLoad::Tabulated(&table),
        8000, 4801, Some(512)), Err(DrivingPointError::Duct(DuctError::BadParameter { .. }))));
}

#[test]
fn sourced_relaxation_rejects_a_foreign_wall_even_when_equilibrium_modulus_matches() {
    let original = material(0.0003);
    let foreign = material(0.003);
    assert_eq!(original.property("young_modulus").unwrap().value_si(),
        foreign.property("young_modulus").unwrap().value_si());
    assert_ne!(original.identity(), foreign.identity());
    let spectrum = WallRelaxation::from_state_point(&original,
        &[("wall_delta_e".into(), "wall_tau".into())], BAND, 0.5, 0.5).unwrap();
    assert_eq!(spectrum.source_identity(), Some(original.identity()));
    assert!(LocalMaterialWall::new(original, THICKNESS, AxialConstraint::StressFree,
        WallMaterialLoss::Maxwell(spectrum.clone()), BAND).is_ok());
    assert!(LocalMaterialWall::new(foreign.clone(), THICKNESS, AxialConstraint::StressFree,
        WallMaterialLoss::Maxwell(spectrum.clone()), BAND).is_err());

    // The graph wall shares this binding; it cannot bypass source identity.
    let gas = *field(300.0).state_at(0.0).unwrap().state();
    let graph = TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
            NetworkNode::Termination { reflection: -0.9 }],
        sections: (0..2).map(|i| TubeSection { nodes: [i,i+1],
            length_m: LENGTH/2.0, radius_m: RADIUS, max_length_error_m: 0.005 }).collect(),
        sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
    };
    let assignment = MaterialWallAssignment { node: 1, material: foreign.clone(),
        model: Default::default(),
        thickness_m: THICKNESS, axial: AxialConstraint::StressFree,
        loss: WallMaterialLoss::Maxwell(spectrum.clone()) };
    assert!(bind_material_walls(graph, vec![assignment], BAND, 1,
        &CancelGate::new_clock_free()).is_err());

    // Explicit authored spectra retain their established unbound semantics.
    let authored = WallRelaxation::new(spectrum.material().clone(),
        "explicit authored numerical spectrum".into(), BAND, 0.5, 0.5).unwrap();
    assert_eq!(authored.source_identity(), None);
    assert!(LocalMaterialWall::new(foreign, THICKNESS, AxialConstraint::StressFree,
        WallMaterialLoss::Maxwell(authored), BAND).is_ok());
}

#[test]
fn continuous_material_wall_matches_original_graph_hereditary_impedance() {
    let law = wall(0.0003);
    let gas = *field(300.0).state_at(0.0).unwrap().state();
    let graph = TubeNetworkSpec {
        nodes: vec![NetworkNode::Inlet, NetworkNode::Junction,
            NetworkNode::Termination { reflection: -0.9 }],
        sections: (0..2).map(|i| TubeSection { nodes: [i, i + 1],
            length_m: LENGTH / 2.0, radius_m: RADIUS, max_length_error_m: 0.005 }).collect(),
        sound_speed_m_s: gas.sound_speed, max_wave_memory_bytes: 1_000_000,
    };
    let source = MaterialWallAssignment { node: 1, material: law.material().clone(),
        model: Default::default(),
        thickness_m: law.thickness_m(), axial: law.axial_constraint(), loss: law.loss().clone() };
    let (_, walls) = bind_material_walls(graph, vec![source.clone()], BAND, 1,
        &CancelGate::new_clock_free()).unwrap();
    let discrete = &walls[0];
    assert_eq!(discrete.source(), &source);
    assert_eq!(law.base_patch(RADIUS, discrete.length_m()).unwrap(), discrete.patch());

    // Independent evaluation of the graph's existing pressure/volume-flow RC
    // declaration. Its wetted area converts it to the continuum surface law.
    let load = discrete.load();
    let base = load.base();
    for frequency in [2.0, 70.0, 430.0, 1700.0, 4000.0] {
        let omega = TAU * frequency;
        let s = C64::new(0.0, -omega);
        let mut volume = C64::from_re(base.resistance_pa_s_m3)
            + s.scale(base.inertance_pa_s2_m3)
            + s.scale(base.compliance_m3_pa.unwrap()).recip();
        for term in load.compliance_terms() {
            volume = volume + (s + C64::from_re(term.rate_per_s)).recip()
                .scale(term.resistance_pa_s_m3 * term.rate_per_s);
        }
        let expected = volume.scale(discrete.patch().area_m2);
        let actual = law.specific_impedance(RADIUS, omega).unwrap();
        assert!(actual.re >= 0.0);
        assert!(relative(actual, expected) < 3e-13, "{frequency} Hz: {actual:?} != {expected:?}");
    }

    // The law uses each actual radius, not the graph's original strip radius.
    let narrow = law.base_patch(RADIUS, 1.0).unwrap().wall;
    let wide = law.base_patch(2.0 * RADIUS, 1.0).unwrap().wall;
    assert_eq!(wide.surface_density, narrow.surface_density);
    assert_eq!(wide.stiffness_per_area, narrow.stiffness_per_area / 4.0);
    assert!(law.specific_impedance(THICKNESS, TAU * 430.0).is_err());
    assert!(law.specific_impedance(RADIUS, TAU * 5000.0).is_err());
}

#[test]
fn elastic_material_recovers_the_independent_compliant_cylinder_limit() {
    let law = LocalMaterialWall::new(material(0.0003), THICKNESS,
        AxialConstraint::StressFree, WallMaterialLoss::Elastic, BAND).unwrap();
    let gas = *field(300.0).state_at(0.0).unwrap().state();
    let area = core::f64::consts::PI * RADIUS * RADIUS;
    let pin = law.base_patch(RADIUS, LENGTH).unwrap().wall;
    for geometry in [cylinder(), Duct { segments: vec![Segment::Profile {
        profile: std::sync::Arc::new(RadiusProfile::exponential(RADIUS, 0.0, LENGTH).unwrap()),
        steps: RiccatiSteps(32),
    }] }] {
        let regional = RegionalDuct::new(&geometry, &[gas]).unwrap();
        let bound = regional.clone().with_material_wall(law.clone());
        for frequency in [30.0, 100.0, 430.0, 700.0] {
            let omega = TAU * frequency;
            let compliance = area / (gas.density * gas.sound_speed * gas.sound_speed)
                + TAU * RADIUS / (pin.stiffness_per_area - omega * omega * pin.surface_density);
            let inertance = gas.density / area;
            let k = omega * fs_math::det::sqrt(inertance * compliance);
            let zc = fs_math::det::sqrt(inertance / compliance);
            let phase = k * LENGTH;
            let expected = C64::new(0.0, -zc * fs_math::det::sin(phase) / fs_math::det::cos(phase));
            let actual = bound.response(omega, LossModel::Lossless, Termination::IdealOpen).unwrap();
            assert!(relative(actual.impedance, expected) < 3e-11);
            let old = regional.response(omega, LossModel::Lossless,
                Termination::IdealOpen, Some(&pin)).unwrap();
            let explicit = regional.response_with_wall(omega, LossModel::Lossless,
                Termination::IdealOpen, WallModel::Pin(&pin)).unwrap();
            assert_eq!(old.impedance, explicit.impedance);
            assert_eq!(old.u_mouth_over_u_in, explicit.u_mouth_over_u_in);
            assert!(relative(actual.impedance, old.impedance) < 3e-11);
        }
    }
}

#[test]
fn simultaneous_gas_geometry_and_material_wall_refine_together() {
    let source = field(270.0);
    let geometry = Duct { segments: vec![Segment::Profile {
        profile: std::sync::Arc::new(RadiusProfile::exponential(RADIUS, 0.8, LENGTH).unwrap()),
        steps: RiccatiSteps(32),
    }] };
    let law = wall(0.0003);
    let evaluate = |steps| RegionalDuct::from_continuous_gas_field(
        &geometry, &source, RiccatiSteps(steps)).unwrap().with_material_wall(law.clone())
        .response(TAU * 430.0, LossModel::Bessel, Termination::FlangedOpen).unwrap();
    let reference = evaluate(256);
    let errors = [8, 16, 32].map(|steps| {
        let response = evaluate(steps);
        relative(response.impedance, reference.impedance)
            .max(relative(response.u_mouth_over_u_in, reference.u_mouth_over_u_in))
    });
    assert!(errors[0] > 1e-10, "joint spatial dependence must be resolved: {errors:?}");
    for pair in errors.windows(2) {
        let order = fs_math::det::ln(pair[0] / pair[1]) / fs_math::det::ln(2.0);
        assert!((order - 4.0).abs() < 0.35, "joint gas/material order={order}: {errors:?}");
    }
}

#[test]
fn changing_only_source_relaxation_reaches_feedback_and_physical_outlet_pressure() {
    let source = field(270.0);
    let duct = RegionalDuct::from_continuous_gas_field(
        &cylinder(), &source, RiccatiSteps(32)).unwrap();
    let a = duct.clone().with_material_wall(wall(0.0003));
    let b = duct.with_material_wall(wall(0.003));
    assert_ne!(a.wall().material().identity(), b.wall().material().identity());
    assert_eq!(a.wall().base_patch(RADIUS, LENGTH).unwrap(),
        b.wall().base_patch(RADIUS, LENGTH).unwrap());
    assert_eq!(a.duct().inlet_gas(), b.duct().inlet_gas());
    let omega = TAU * 430.0;
    let reflection = a.reflection(omega, LossModel::Bessel, Termination::FlangedOpen).unwrap();
    let response = a.response(omega, LossModel::Bessel, Termination::FlangedOpen).unwrap();
    let zc = C64::from_re(a.duct().inlet_impedance());
    assert!(relative(reflection, (response.impedance - zc) / (response.impedance + zc)) < 2e-12);
    assert!(reflection.abs() <= 1.0 + 1e-12);

    let mut feedback_a = a.characteristic_line(Termination::FlangedOpen, 8000, 4801, Some(512)).unwrap();
    let mut feedback_b = b.characteristic_line(Termination::FlangedOpen, 8000, 4801, Some(512)).unwrap();
    let mut mic_a = a.baffled_outlet(8000, 512, receiver()).unwrap();
    let mut mic_b = b.baffled_outlet(8000, 512, receiver()).unwrap();
    let first = mic_a.info().axial_delay_samples + mic_a.info().exterior_delay_samples.0;
    let (mut feedback_change, mut pressure_change, mut pressure_energy) = (0.0, 0.0, 0.0);
    for n in 0..1024 {
        let input = if n == 0 { 100.0 } else { 0.0 };
        let fa = feedback_a.push(input).unwrap();
        let fb = feedback_b.push(input).unwrap();
        let pa = mic_a.step(input).unwrap();
        let pb = mic_b.step(input).unwrap();
        assert!([fa, fb, pa, pb].iter().all(|x| x.is_finite()));
        if n < first { assert_eq!(pa, 0.0); assert_eq!(pb, 0.0); }
        feedback_change += (fa - fb) * (fa - fb);
        pressure_change += (pa - pb) * (pa - pb);
        pressure_energy += pa * pa;
    }
    assert!(feedback_change > 1e-6);
    assert!(pressure_energy > 1e-12 && pressure_change > 1e-6 * pressure_energy);
    // Both consumers refuse the same incomplete source-frequency coverage.
    assert!(a.characteristic_line(Termination::FlangedOpen, 10000, 4801, Some(512)).is_err());
    assert!(a.baffled_outlet(10000, 512, receiver()).is_err());
}

/// G1: solve two different local membranes and gases by independent analytical
/// pressure/volume-flow matrices. Neither wall nor gas is averaged at the join.
#[test]
fn regional_wall_sources_match_independent_two_material_transfer() {
    let gases = [*field(300.0).state_at(0.0).unwrap().state(),
        *field(270.0).state_at(LENGTH).unwrap().state()];
    let radii = [RADIUS, 1.3 * RADIUS];
    let lengths = [0.2, LENGTH - 0.2];
    let times = [0.0003, 0.003];
    let geometry = Duct { segments: (0..2).map(|i| Segment::Cylinder {
        radius: radii[i], length: lengths[i],
    }).collect() };
    let walls = times.map(wall);
    let bound = RegionalDuct::new(&geometry, &gases).unwrap()
        .with_material_wall_regions(&walls).unwrap();
    for frequency in [30.0, 100.0, 430.0, 700.0] {
        let omega = TAU * frequency;
        let blocks = [0, 1].map(|i| {
            let s = C64::new(0.0, -omega);
            let r = radii[i];
            let area = core::f64::consts::PI * r * r;
            // Original manufactured Maxwell source: E_inf=30 MPa,
            // delta E=20 MPa, rho=1000 kg/m^3, free axial stress.
            let wall_z = s.scale(1000.0 * THICKNESS)
                + s.recip().scale(3e7 * THICKNESS / (r * r))
                + (s + C64::from_re(1.0 / times[i])).recip()
                    .scale(2e7 * THICKNESS / (r * r));
            let z = s.scale(gases[i].density / area);
            let y = s.scale(area / (gases[i].density * gases[i].sound_speed.powi(2)))
                + wall_z.recip().scale(TAU * r);
            let gamma = (z * y).sqrt();
            let x = gamma.re * lengths[i];
            let phase = gamma.im * lengths[i];
            let cosine = C64::new(x.cosh() * phase.cos(), x.sinh() * phase.sin());
            let sine = C64::new(x.sinh() * phase.cos(), x.cosh() * phase.sin());
            [cosine, z / gamma * sine, y / gamma * sine, cosine]
        });
        let [a, b, c, d] = blocks[0];
        let [_, f, _, h] = blocks[1];
        let numerator = a * f + b * h;
        let denominator = c * f + d * h;
        let actual = bound.response(omega, LossModel::Lossless, Termination::IdealOpen, None).unwrap();
        assert!(relative(actual.impedance, numerator / denominator) < 5e-11);
        assert!(relative(actual.u_mouth_over_u_in, denominator.recip()) < 5e-11);
        assert!(actual.impedance.re >= 0.0);
        let reflection = bound.reflection_with_wall(omega, LossModel::Lossless,
            Termination::IdealOpen, WallModel::Rigid).unwrap();
        let zc = C64::from_re(bound.inlet_impedance());
        assert!(relative(reflection, (actual.impedance - zc) / (actual.impedance + zc)) < 2e-12);
        assert!(reflection.abs() <= 1.0 + 1e-12);
    }
}

/// G3/G4: regional subdivision preserves the local material, while incomplete
/// source partitions and external overrides cannot erase its constitutive law.
#[test]
fn regional_wall_partition_preserves_uniform_limit_and_source_refusals() {
    let gas = *field(300.0).state_at(0.0).unwrap().state();
    let law = wall(0.0003);
    let geometry = Duct { segments: vec![Segment::Cylinder {
        radius: RADIUS, length: LENGTH / 2.0,
    }; 2] };
    let original = RegionalDuct::new(&geometry, &[gas; 2]).unwrap();
    let bound = original.clone().with_material_wall_regions(&[law.clone(), law.clone()]).unwrap();
    let uniform = original.clone().with_material_wall(law.clone());
    for frequency in [30.0, 430.0, 1700.0] {
        let actual = bound.response(TAU * frequency, LossModel::Bessel,
            Termination::FlangedOpen, None).unwrap();
        let expected = uniform.response(TAU * frequency, LossModel::Bessel,
            Termination::FlangedOpen).unwrap();
        assert_eq!(actual.impedance, expected.impedance);
        assert_eq!(actual.u_mouth_over_u_in, expected.u_mouth_over_u_in);
    }
    let split_geometry = Duct { segments: vec![Segment::Cylinder {
        radius: RADIUS, length: LENGTH / 4.0,
    }; 4] };
    let split = RegionalDuct::new(&split_geometry, &[gas; 4]).unwrap()
        .with_material_wall_regions(&vec![law.clone(); 4]).unwrap();
    let omega = TAU * 430.0;
    let actual = bound.response(omega, LossModel::Bessel, Termination::FlangedOpen, None).unwrap();
    let refined = split.response(omega, LossModel::Bessel, Termination::FlangedOpen, None).unwrap();
    assert!(relative(actual.impedance, refined.impedance) < 2e-12);
    assert!(relative(actual.u_mouth_over_u_in, refined.u_mouth_over_u_in) < 2e-12);
    assert!(original.clone().with_material_wall_regions(&[law.clone()]).is_err());
    assert!(bound.clone().with_material_wall_regions(&[law.clone(), law.clone()]).is_err());
    assert!(bound.response_with_wall(omega, LossModel::Bessel, Termination::FlangedOpen,
        WallModel::Law(&law)).is_err());
    assert!(bound.response(TAU * 5000.0, LossModel::Bessel, Termination::FlangedOpen, None).is_err());
    let continuous = RegionalDuct::from_continuous_gas_field(&cylinder(), &field(270.0),
        RiccatiSteps(16)).unwrap();
    assert!(continuous.with_material_wall_regions(&[law]).is_err());
    let table = fs_duct::TabulatedLoad::try_new(
        vec![TAU * 10.0, TAU * 4000.0], vec![C64::ZERO; 2],
        "ideal-open source-domain test", 0.0,
    ).unwrap().with_mouth(fs_duct::TabulatedMouth {
        axial_position_m: LENGTH, radius_m: RADIUS, gas, dc_impedance: 0.0,
    }).unwrap();
    // The old constant-pin shortcut must not bypass the source's missing DC.
    assert!(bound.characteristic_line_with_wall(MouthLoad::Tabulated(&table),
        8000, 4801, WallModel::Rigid, Some(512)).is_err());
}
