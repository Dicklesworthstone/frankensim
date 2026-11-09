//! Original material walls, interior work and external sound share one solve.

use super::*;
use crate::bernoulli_aperture::wall::material::{
    AxialConstraint, LocalMaterialWall, WallMaterialLoss, relaxation::WallRelaxation,
};
use fs_duct::WallLaw;
use fs_duct::molecular::MolecularCompliantWall;
use fs_duct::network::field::{NetworkFieldBudget, NetworkFieldProbe, NetworkFieldSection};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};
use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};

mod exterior;

struct CountedWall {
    material: LocalMaterialWall,
    calls: AtomicUsize,
}
impl WallLaw for CountedWall {
    fn specific_impedance(&self, radius: f64, omega: f64) -> Result<C64, DuctError> {
        self.calls.fetch_add(1, Ordering::Relaxed);
        WallLaw::specific_impedance(&self.material, radius, omega)
    }
}

// The original resolved state and native Maxwell membrane are retained.
// These explicit synthetic claims are not measured material qualification.
fn wall() -> Arc<CountedWall> {
    let band = [1.0, 4500.0];
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 1000.0),
        ("young_modulus", Pressure::DIMS, 3e7),
        ("poisson_ratio", Dims([0; 6]), 0.3),
        ("wall_delta_e", Pressure::DIMS, 2e7),
        ("wall_tau", Dims([0, 0, 1, 0, 0, 0]), 0.0003),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 300.0)
                .with("omega", core::f64::consts::TAU * band[0], core::f64::consts::TAU * band[1]),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic resolved Maxwell wall for original-owner integration".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_with_quantity(
            name, QuantitySpec::dimensional(dims), ScalarAdmissibility::Finite,
        ).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic-wall".into(), phase: "solid".into(),
        process: "manufactured integration fixture".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1000.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let memory = WallRelaxation::from_state_point(
        &material, &[("wall_delta_e".into(), "wall_tau".into())], band, 0.5, 0.5,
    ).unwrap();
    Arc::new(CountedWall {
        material: LocalMaterialWall::new(
            material, 0.0002, AxialConstraint::StressFree, WallMaterialLoss::Maxwell(memory), band,
        ).unwrap(),
        calls: AtomicUsize::new(0),
    })
}

fn walls(model: &RadiatingNetwork, wall: &Arc<CountedWall>)
    -> Vec<(String, MolecularCompliantWall)> {
    model.driving_point().network().edges.iter().map(|edge| {
        let source: MolecularCompliantWall = wall.clone();
        (edge.id.clone(), source)
    }).collect()
}

fn probes(model: &RadiatingNetwork) -> (Vec<NetworkFieldProbe>, Vec<[usize; 2]>) {
    let mut probes = Vec::new();
    let mut cells = Vec::new();
    for edge in &model.driving_point().network().edges {
        let fs_duct::Segment::Cylinder { length, .. } = edge.duct.segments[0] else {
            panic!("original fixture consists of cylinders");
        };
        let start = probes.len();
        for offset_m in [0.0, length * 0.5, length] {
            probes.push(NetworkFieldProbe {
                edge_id: edge.id.clone(),
                section: NetworkFieldSection::Cylinder { segment: 0, offset_m },
            });
        }
        cells.push([start, start + 2]);
    }
    (probes, cells)
}

fn request<'a>(
    omega: f64,
    bank: &'a CoupledSourceBank,
    probes: &'a [NetworkFieldProbe],
) -> MolecularRadiatingFieldRequest<'a> {
    MolecularRadiatingFieldRequest {
        omega,
        source_clock: "original-molecular-BEM-clock",
        sources: &[],
        source_bank: Some(bank),
        scalar_loads: &[],
        probes,
        budget: NetworkFieldBudget::default(),
    }
}

// Independent complex transmission-line input. Only original material laws
// supply coefficients; no duct wave, scattering, field or network solver is used.
fn shorted_wall_line(state: &RelaxationState, material: &LocalMaterialWall, omega: f64) -> C64 {
    let radius = 0.006;
    let area = core::f64::consts::PI * radius * radius;
    let bulk = state.caloric_response(omega).unwrap().bulk_modulus_pa;
    let s = C64::new(0.0, -omega);
    let z = s.scale(state.gas().density / area);
    let y = s.scale(area) / bulk
        + material.specific_impedance(radius, omega).unwrap().recip()
            .scale(core::f64::consts::TAU * radius);
    let k = (-z * y).sqrt();
    let phase = k.scale(0.13);
    let (a, b) = phase.re.sin_cos();
    let sin = C64::new(a * phase.im.cosh(), b * phase.im.sinh());
    let cos = C64::new(b * phase.im.cosh(), -a * phase.im.sinh());
    z * sin / (k * cos)
}

#[test]
fn one_original_loaded_bank_supplies_material_wall_fields_work_and_microphones() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let source = wall();
    let bindings = walls(&model, &source);
    let retained = bindings[0].1.clone();
    let first_edge = bindings[0].0.clone();
    let molecular = model.bind_molecular_with_compliant_walls(
        edges(&model, 2, 4.0), bindings, 64, &gate,
    ).unwrap();
    assert!(std::ptr::eq(
        molecular.molecular_network().edge(&first_edge).unwrap().compliant_wall().unwrap(),
        retained.as_ref(),
    ));
    assert_eq!(source.calls.load(Ordering::Relaxed), 0);
    let (probes, pairs) = probes(&model);
    let omega = core::f64::consts::TAU * 1200.0;
    let bank = bank(omega);
    let prepared = molecular.prepare_fields(&request(omega, &bank, &probes), &gate).unwrap();
    let exterior = prepared.radiating();
    let interior = prepared.interior();
    let original = interior.fields().port_response();
    assert!(std::ptr::eq(interior.network(), molecular.molecular_network()));
    assert_eq!(original.source_coordinates(), exterior.source_coordinates());
    assert_eq!(original.default_excitations(), exterior.default_excitations());
    assert_eq!(source.calls.load(Ordering::Relaxed), model.driving_point().network().edges.len());
    let amplitudes = exterior.default_excitations();
    let mut outside = exterior.response_buffer();
    exterior.drive_into(amplitudes, &mut outside, &gate).unwrap();
    let mut native = original.response_buffer();
    original.drive_into(amplitudes, &mut native).unwrap();
    assert_eq!(outside.network().unwrap().fields.port_pressures, native.fields.port_pressures);
    assert_eq!(outside.network().unwrap().fields.port_inflows, native.fields.port_inflows);
    let local = interior.fields().sample(amplitudes, || false).unwrap();
    let absorption = interior.prepare_cylinder_absorption(&pairs, pairs.len(), || false).unwrap();
    let losses = absorption.sample(amplitudes, || false).unwrap();
    let internal = losses.iter().map(|p| p.total_removed_power_w()).sum::<f64>();
    let mechanical = losses.iter().map(|p| p.wall_work_w).sum::<f64>();
    assert!(internal > 0.0 && mechanical > 0.0);
    real_close(internal, native.fields.network_power_w);
    real_close(
        native.generator_power_w,
        native.source_dissipation_w + internal + outside.exterior_surface_power_w().unwrap(),
    );
    for (patch_index, patch) in model.patches().iter().enumerate() {
        let edge_index = model.driving_point().network().edges.iter()
            .position(|edge| &edge.id == &patch.port).unwrap();
        let sample = &local[3 * edge_index + 2];
        close(sample.volume_flow_m3_s, outside.outward_flows_m3_s().unwrap()[patch_index]);
        let pressure = outside.outward_flows_m3_s().unwrap().iter().enumerate()
            .fold(C64::ZERO, |sum, (j, &q)| sum + exterior.boundary().z[patch_index * 2 + j] * q);
        close(sample.pressure_pa, pressure);
    }
    for (a, b) in outside.receiver_fields().unwrap().iter()
        .zip(direct_fields(&model, omega, outside.outward_flows_m3_s().unwrap())) {
        close(a.pressure_pa, b.pressure_pa);
        for axis in 0..3 {
            close(a.particle_velocity_m_s[axis], b.particle_velocity_m_s[axis]);
        }
    }
    // Independent antisymmetric line: the junction is pressure release even
    // with the full mutual exterior, so wall feedback must reach the inlet.
    exterior.drive_into(&[C64::ONE, -C64::ONE], &mut outside, &gate).unwrap();
    let impedance = shorted_wall_line(
        molecular.molecular_network().edge("stem").unwrap().state(), &source.material, omega,
    );
    let expected = C64::ONE / (bank.impedance.z[0] - bank.impedance.z[1] + impedance);
    let aux = exterior.ports().iter().position(|p| p == "aux").unwrap();
    close(outside.network().unwrap().fields.port_inflows[aux], expected);
    real_close(outside.network().unwrap().fields.network_power_w, impedance.re * expected.norm_sq());
    assert_eq!(source.calls.load(Ordering::Relaxed), model.driving_point().network().edges.len());
}

#[test]
fn interior_loss_and_directional_radiation_keep_the_same_coherent_source_factors() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let source = wall();
    let molecular = model.bind_molecular_with_compliant_walls(
        edges(&model, 2, 4.0), walls(&model, &source), 64, &gate,
    ).unwrap();
    let (probes, pairs) = probes(&model);
    let omega = core::f64::consts::TAU * 1200.0;
    let bank = bank(omega);
    let directions = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let radiation = PreparedRadiationOptions {
        angular: helmholtz::DirectivityProjectionOptions::default(),
        initial_latitudes: 8,
        refinements: 3,
        maximum_kernel_evaluations: 4_000_000,
        maximum_source_entries: 4096,
    };
    let prepared = molecular.prepare_fields_with_radiation(
        &request(omega, &bank, &probes), &directions, radiation, &gate,
    ).unwrap();
    let exterior = prepared.radiating();
    let factors = [
        C64::new(0.7, -0.2), C64::new(0.3, 0.6),
        C64::new(-0.1, 0.4), C64::new(0.9, -0.3),
    ];
    let absorption = prepared.interior().prepare_cylinder_absorption(
        &pairs, pairs.len(), || false,
    ).unwrap();
    let cells = absorption.observe_coherence(
        None, &[&factors[..2], &factors[2..]], || false,
    ).unwrap();
    let mut statistics = exterior.statistics_buffer();
    exterior.coherence_into(
        &factors, 2, NetworkCoherenceBudget::default(), FIELD_BUDGET, &mut statistics, &gate,
    ).unwrap();
    let work = statistics.network().unwrap().work().unwrap();
    let internal = cells.iter().map(|cell| cell.total_removed_power_w()).sum::<f64>();
    real_close(internal, work.network_power_w);
    real_close(
        work.generator_power_w,
        work.source_dissipation_w + internal + statistics.exterior_surface_power_w().unwrap(),
    );
    let patterns = [
        FirstOrder::new(0.5, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.2, [0.0, 1.0, 0.0]).unwrap(),
    ];
    let mut microphone_covariance = [C64::ZERO; 4];
    let mut outside = exterior.response_buffer();
    let mut independent_cells = vec![0.0; pairs.len()];
    let far = exterior.far_field().unwrap();
    let angular_budget = HarmonicWorkBudget { maximum_products: 1_000_000 };
    let mut far_sum = [C64::ZERO; 2];
    for factor in factors.chunks_exact(2) {
        exterior.drive_into(factor, &mut outside, &gate).unwrap();
        let mic: Vec<_> = patterns.iter().zip(outside.receiver_fields().unwrap()).map(|(p, f)| {
            p.observe(f.pressure_pa, f.particle_velocity_m_s, Medium::from(&model.gas())).unwrap()
        }).collect();
        for i in 0..2 {
            for j in 0..2 {
                microphone_covariance[2 * i + j] = microphone_covariance[2 * i + j]
                    + (mic[i] * mic[j].conj()).scale(0.5);
            }
        }
        for (sum, cell) in independent_cells.iter_mut().zip(absorption.sample(factor, || false).unwrap()) {
            *sum += cell.total_removed_power_w();
        }
        let mut fields = [C64::ZERO; 2];
        let power = far.evaluate_coherent(factor, angular_budget, &mut fields, &gate).unwrap();
        assert!(power.total_power_w > 0.0);
        for i in 0..2 { far_sum[i] = far_sum[i] + fields[i]; }
    }
    for (sum, cell) in independent_cells.into_iter().zip(&cells) {
        real_close(sum, cell.total_removed_power_w());
    }
    let mut microphones = exterior.microphones(&patterns).unwrap();
    let covariance = microphones.coherence(&factors, 2, FIELD_BUDGET, &gate).unwrap();
    for (&a, b) in covariance.iter().zip(microphone_covariance) { close(a, b); }
    let coherent = [factors[0] + factors[2], factors[1] + factors[3]];
    let mut fields = [C64::ZERO; 2];
    far.observe_coherent(&coherent, angular_budget, &mut fields, &gate).unwrap();
    for i in 0..2 { close(fields[i], far_sum[i]); }
    assert_eq!(source.calls.load(Ordering::Relaxed), model.driving_point().network().edges.len());
}

#[test]
fn paired_preparation_preserves_original_domain_budgets_and_retry_after_refusal() {
    let gate = CancelGate::new_clock_free();
    let omega = core::f64::consts::TAU * 1200.0;
    let drive = bank(omega);
    let source = wall();
    let identity = source.material.material().identity();

    let mut tight = regional();
    let (n, p, r) = (tight.surface().areas().len(), tight.patches().len(), tight.receivers().len());
    let ports = tight.driving_point().network().nodes.iter()
        .filter(|node| matches!(node.boundary, NodeBoundary::Port { .. })).count();
    let count = drive.impedance.labels.len();
    // Exactly the old port/exterior allowance, without the extra retained
    // field-owner bank. Refusal occurs before any material-frequency query.
    tight.budget.max_field_entries = 5 * n * p + 8 * p * r + 4 * count * (ports + r);
    let constrained = tight.bind_molecular_with_compliant_walls(
        edges(&tight, 2, 4.0), walls(&tight, &source), 64, &gate,
    ).unwrap();
    let (points, _) = probes(&tight);
    assert!(matches!(
        constrained.prepare_fields(&request(omega, &drive, &points), &gate),
        Err(ExteriorNetworkError::Invalid(
            "complete prepared source/exterior basis exceeds its declared budget"
        ))
    ));
    assert_eq!(source.calls.load(Ordering::Relaxed), 0);

    let model = regional();
    let mut wrong = walls(&model, &source);
    wrong[0].0 = "foreign-edge".into();
    assert!(model.bind_molecular_with_compliant_walls(
        edges(&model, 2, 4.0), wrong, 64, &gate,
    ).is_err());
    let molecular = model.bind_molecular_with_compliant_walls(
        edges(&model, 2, 4.0), walls(&model, &source), 64, &gate,
    ).unwrap();
    let (points, _) = probes(&model);
    let mut request = request(omega, &drive, &points);
    request.budget.maximum_transfer_entries = 0;
    assert!(molecular.prepare_fields(&request, &gate).is_err());
    assert_eq!(source.calls.load(Ordering::Relaxed), 0);
    request.budget = NetworkFieldBudget::default();
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        molecular.prepare_fields(&request, &stopped),
        Err(ExteriorNetworkError::Cancelled)
    ));
    let accepted = molecular.prepare_fields(&request, &gate).unwrap();
    let original_pressure = accepted.interior().fields().pressure_transfer().to_vec();
    let original_flow = accepted.interior().fields().flow_transfer().to_vec();
    let unsupported_frequency = core::f64::consts::TAU * 4501.0;
    let unsupported_bank = bank(unsupported_frequency);
    let unsupported = MolecularRadiatingFieldRequest {
        omega: unsupported_frequency,
        source_clock: request.source_clock,
        sources: &[],
        source_bank: Some(&unsupported_bank),
        scalar_loads: &[],
        probes: &points,
        budget: request.budget,
    };
    assert!(matches!(
        molecular.prepare_fields(&unsupported, &gate),
        Err(ExteriorNetworkError::Molecular(_))
    ));
    let replay = molecular.prepare_fields(&request, &gate).unwrap();
    assert_eq!(replay.interior().fields().pressure_transfer(), original_pressure);
    assert_eq!(replay.interior().fields().flow_transfer(), original_flow);
    assert_eq!(accepted.radiating().boundary().z, replay.radiating().boundary().z);
    assert_eq!(source.material.material().identity(), identity);
}
