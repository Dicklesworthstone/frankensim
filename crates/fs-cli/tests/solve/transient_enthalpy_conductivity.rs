//! Native phase conductivity against an independent two-state tetrahedron.
//!
//! Three identical wall vertices share one solid state; the interior vertex
//! melts. Eliminating the wall equation leaves a bounded scalar energy root.
//! The oracle uses no material evaluator, finite-element assembler or solver.

use super::*;
use fs_project::{
    TransientEnthalpy, TransientEnthalpyKnot, TransientPhaseConductivity,
    TransientRegionEnthalpy, spec::dims,
};

#[allow(dead_code)]
#[path = "../../src/json_read.rs"]
mod json;
use json::JsonValue as J;

const AMBIENT_K: f64 = 300.0;
const MELTING_K: f64 = 301.0;
const POWER_W: f64 = 5.0;
const HTC: f64 = 10.0;
const HORIZON: f64 = 2.0;
const ENTHALPY_SOURCE: &str =
    "synthetic Cp=1 and 10 J/kg latent interval for a spatial conductivity oracle";
const CONDUCTIVITY_SOURCE: &str =
    "synthetic linear liquid-mass-fraction conductivity multiplier; not measured material data";

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

fn n(row: &J, key: &str) -> f64 {
    row.f64_field(key).unwrap()
}

fn values(row: &J, key: &str) -> Vec<f64> {
    row.get(key)
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|value| value.as_f64().unwrap())
        .collect()
}

fn face_area() -> f64 {
    3.0_f64.sqrt() / 2.0
}

fn cards(conductivity: f64) -> CardPackSet {
    CardPackSet::admit(vec![raw_pack(
        CardPackKind::Material,
        "fixtures/enthalpy-phase-conductivity.fsmcdpk",
        material_pack_bytes_with_property(
            "AA6061",
            &format!("enthalpy-phase-reference-k-{conductivity}"),
            280.0,
            360.0,
            fs_matdb::PropertyValue::Scalar {
                value: conductivity,
                dims: CONDUCTIVITY_DIMS,
            },
            fs_matdb::InterpolationPolicy::ConstantWithinValidity,
        ),
    )])
    .unwrap()
}

pub(super) fn project(cards: &CardPackSet, multipliers: Option<(f64, f64)>) -> ProjectSpec {
    let bytes = tetra_stl();
    let mut spec = transient_product::driven_face_project(&bytes);
    rebind_to(&mut spec, cards);
    let identity = fs_blake3::ContentHash::from_hex(
        &spec.materials.as_ref().unwrap()[0].card,
    )
    .unwrap();
    spec.power.as_mut().unwrap()[0].watts.value = POWER_W;
    spec.solver.as_mut().unwrap().tolerance_rel = 1e-12;
    let setup = spec
        .cooling
        .as_mut()
        .unwrap()
        .conduction
        .as_mut()
        .unwrap();
    setup.boundaries[0].condition = ThermalBoundaryCondition::Convection {
        coefficient: QtyAny::new(HTC, dims::HEAT_TRANSFER_COEFFICIENT),
        reference_temperature: QtyAny::new(AMBIENT_K, dims::TEMPERATURE),
    };
    let time = setup.transient.as_mut().unwrap();
    time.initial_temperature = fs_project::TransientInitialTemperature::Uniform(
        QtyAny::new(AMBIENT_K, dims::TEMPERATURE));
    time.horizon.value = HORIZON;
    time.max_step.value = 0.5;
    time.energy_tolerance.value = 1e-8;
    time.capacities.clear();
    time.enthalpy = Some(TransientEnthalpy {
        max_newton_iterations: 32,
        max_backtracks: 24,
        materials: vec![TransientRegionEnthalpy {
            region: "air".into(),
            material_card_identity: identity,
            source: ENTHALPY_SOURCE.into(),
            // Unit right tetrahedron V=1/6: rho0=6 gives total mass 1 kg.
            reference_density: QtyAny::new(6.0, dims::DENSITY),
            initial_specific_enthalpy: QtyAny::new(0.0, dims::SPECIFIC_ENTHALPY),
            phase_conductivity: multipliers.map(|(solid, liquid)| TransientPhaseConductivity {
                solid_multiplier: solid,
                liquid_multiplier: liquid,
                source: CONDUCTIVITY_SOURCE.into(),
            }),
            knots: [
                (-20.0, 280.0, 0.0),
                (1.0, MELTING_K, 0.0),
                (11.0, MELTING_K, 1.0),
                (40.0, 330.0, 1.0),
            ]
            .into_iter()
            .map(|(h, temperature, fraction)| TransientEnthalpyKnot {
                specific_enthalpy: QtyAny::new(h, dims::SPECIFIC_ENTHALPY),
                temperature: QtyAny::new(temperature, dims::TEMPERATURE),
                liquid_mass_fraction: fraction,
            })
            .collect(),
        }],
    });
    spec
}

struct ThermalRun {
    receipt: J,
    qoi: J,
    field: J,
}

impl ThermalRun {
    fn transient(&self) -> &J {
        self.receipt.get("transient").unwrap()
    }

    fn rows(&self, grid: &str) -> &[J] {
        self.transient().get(grid).unwrap().as_array().unwrap()
    }

    fn endpoint(&self) -> &J {
        self.transient().get("endpoint").unwrap()
    }
}

fn run(cards: &CardPackSet, multipliers: Option<(f64, f64)>) -> ThermalRun {
    let spec = project(cards, multipliers);
    let decoded = decode(&spec);
    assert!(decoded.findings().is_empty(), "{:?}", decoded.findings());
    let ledger = Ledger::open(":memory:").unwrap();
    import_fixture(&ledger, &spec, tetra_stl());
    let (_, qoi, receipt, field) = run_conjugate_to_completion(&ledger, &decoded, cards);
    ThermalRun {
        receipt: J::parse(&receipt).unwrap(),
        qoi: J::parse(&qoi).unwrap(),
        field: J::parse(&field).unwrap(),
    }
}

/// This oracle only needs the solid segment and the latent interval. It
/// deliberately refuses a candidate beyond the independently checked regime.
fn state(h: f64) -> (f64, f64) {
    assert!((0.0..=11.0).contains(&h));
    if h <= 1.0 {
        (AMBIENT_K + h, 0.0)
    } else {
        (MELTING_K, (h - 1.0) / 10.0)
    }
}

fn reference_step(
    old: [f64; 2],
    dt: f64,
    base_k: f64,
    solid: f64,
    liquid: f64,
) -> [f64; 2] {
    // Each vertex has 1/4 kg and receives P/4. The interior conduction
    // row is k*(Ta-Tb)/2. Each wall row is k*(Tb-Ta)/6 plus
    // HTC*A*(Tb-300)/3. Wall nodes remain solid with Tb=300+hb.
    let candidate = |ha: f64| {
        let (ta, fa) = state(ha);
        let k = base_k * (solid + (liquid - solid) * fa / 4.0);
        let hb = (0.25 * old[1] + dt * (k * (ta - AMBIENT_K) / 6.0 + POWER_W / 4.0))
            / (0.25 + dt * (k / 6.0 + HTC * face_area() / 3.0));
        let residual = 0.25 * (ha - old[0])
            + dt * (k * (ta - (AMBIENT_K + hb)) / 2.0 - POWER_W / 4.0);
        (residual, hb)
    };
    let (mut low, mut high) = (old[0], 11.0);
    assert!(candidate(low).0 <= 0.0 && candidate(high).0 >= 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (low + high);
        if candidate(middle).0 > 0.0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    let ha = 0.5 * (low + high);
    let hb = candidate(ha).1;
    assert!((0.0..1.0).contains(&hb), "the three wall nodes must remain solid");
    [ha, hb]
}

fn check_attribution(solved: &ThermalRun, multipliers: Option<(f64, f64)>) {
    let materials = solved
        .transient()
        .get("materials")
        .unwrap()
        .as_array()
        .unwrap();
    assert_eq!(materials.len(), 1);
    let material = &materials[0];
    assert_eq!(material.str_field("source"), Some(ENTHALPY_SOURCE));
    match multipliers {
        None => assert!(material.get("phase_conductivity").is_none()),
        Some((solid, liquid)) => {
            let law = material.get("phase_conductivity").unwrap();
            assert_eq!(law.str_field("law"), Some("linear-liquid-mass-fraction"));
            assert_eq!(
                law.str_field("coordinate"),
                Some("element-mean-nodal-liquid-mass-fraction")
            );
            assert_eq!(law.str_field("authority"), Some("Declared"));
            assert_eq!(law.str_field("source"), Some(CONDUCTIVITY_SOURCE));
            close(n(law, "solid_multiplier"), solid, 0.0);
            close(n(law, "liquid_multiplier"), liquid, 0.0);
        }
    }
    assert_eq!(solved.transient().str_field("authority"), Some("Estimated"));
}

fn check_reference(solved: &ThermalRun, base_k: f64, solid: f64, liquid: f64) {
    let history = solved.transient();
    close(n(history, "total_steps"), 12.0, 0.0);
    let mut final_h = [0.0; 2];
    let mut final_integral = 0.0;
    for (grid, count) in [("coarse", 4), ("fine", 8)] {
        let rows = solved.rows(grid);
        assert_eq!(rows.len(), count);
        let mut old = [0.0; 2];
        let mut elapsed = 0.0;
        let mut integrated = 0.0;
        for row in rows {
            let dt = n(row, "dt_s");
            let next = reference_step(old, dt, base_k, solid, liquid);
            let (ta, fa) = state(next[0]);
            let stored = 0.25 * (next[0] - old[0]) + 0.75 * (next[1] - old[1]);
            let net_w = POWER_W - HTC * face_area() * next[1];
            close(n(row, "time_s") - elapsed, dt, 1e-15);
            close(n(row, "final_region_max_k"), ta, 2e-6);
            close(n(row, "source_w"), POWER_W, 2e-12);
            close(n(row, "stored_energy_change_j"), stored, 2e-6);
            close(n(row, "net_input_w"), net_w, 2e-6);
            close(n(row, "mean_liquid_mass_fraction"), fa / 4.0, 2e-7);
            close(n(row, "liquid_mass_kg"), fa / 4.0, 2e-7);
            close(stored - dt * net_w, 0.0, 1e-10);
            assert!(n(row, "residual_j") <= n(row, "threshold_j"));
            assert!(n(row, "energy_residual_j").abs() <= n(history, "energy_tolerance_j"));
            integrated += dt * net_w;
            old = next;
            elapsed = n(row, "time_s");
        }
        close(elapsed, HORIZON, 0.0);
        close(0.25 * old[0] + 0.75 * old[1], integrated, 1e-10);
        if grid == "fine" {
            final_h = old;
            final_integral = integrated;
        }
    }

    let temperature = values(&solved.field, "temperature");
    let enthalpy = values(solved.endpoint(), "specific_enthalpy_j_kg");
    let fraction = values(solved.endpoint(), "liquid_mass_fraction");
    let mass = values(solved.endpoint(), "reference_nodal_mass_kg");
    assert_eq!(temperature.len(), 4, "the oracle requires the imported unit tetrahedron");
    assert_eq!(enthalpy.len(), 4);
    assert_eq!(fraction.len(), 4);
    assert_eq!(mass.len(), 4);
    let mut nodes: Vec<_> = (0..4)
        .map(|index| {
            close(mass[index], 0.25, 2e-12);
            (temperature[index], enthalpy[index], fraction[index])
        })
        .collect();
    nodes.sort_by(|left, right| left.0.total_cmp(&right.0));
    for &(t, h, f) in &nodes[..3] {
        close(t, AMBIENT_K + final_h[1], 2e-6);
        close(h, final_h[1], 2e-6);
        close(f, 0.0, 0.0);
    }
    let (ta, fa) = state(final_h[0]);
    close(nodes[3].0, ta, 2e-6);
    close(nodes[3].1, final_h[0], 2e-6);
    close(nodes[3].2, fa, 2e-7);
    assert!(nodes[3].1 > 1.0 && nodes[3].1 < 11.0);

    // Re-integrate the endpoint using only retained fields, known nodal
    // masses and the actual slanted-face trace, independently of receipt sums.
    let retained_storage: f64 = mass.iter().zip(&enthalpy).map(|(m, h)| m * h).sum();
    let wall_mean = nodes[..3].iter().map(|node| node.0).sum::<f64>() / 3.0;
    let boundary_w = HTC * face_area() * (wall_mean - AMBIENT_K);
    let energy = history.get("energy").unwrap();
    close(n(energy, "stored_change_j"), retained_storage, 2e-6);
    close(n(energy, "stored_change_j"), final_integral, 2e-6);
    close(n(energy, "integrated_net_input_j"), final_integral, 2e-6);
    close(
        n(solved.rows("fine").last().unwrap(), "net_input_w"),
        POWER_W - boundary_w,
        2e-6,
    );
    close(n(solved.receipt.get("energy").unwrap(), "robin_out_w"), boundary_w, 2e-6);
    let qoi = &solved.qoi.get("qoi").unwrap().as_array().unwrap()[0];
    close(n(qoi, "value"), nodes[3].0, 0.0);
    let scope = solved.qoi.get("thermal_time_scope").unwrap();
    assert_eq!(scope.str_field("kind"), Some("final"));
    close(n(scope, "time_s"), HORIZON, 0.0);
}

fn compare_fields(left: &ThermalRun, right: &ThermalRun) {
    for (left_row, right_row, keys) in [
        (&left.field, &right.field, &["temperature"][..]),
        (
            left.endpoint(),
            right.endpoint(),
            &["specific_enthalpy_j_kg", "liquid_mass_fraction", "reference_nodal_mass_kg"][..],
        ),
    ] {
        for key in keys {
            let a = values(left_row, key);
            let b = values(right_row, key);
            assert_eq!(a.len(), b.len());
            for (a, b) in a.into_iter().zip(b) {
                close(a, b, 2e-7);
            }
        }
    }
}

#[test]
fn g1_native_enthalpy_omitted_phase_conductivity_matches_explicit_identity() {
    let base = cards(2.0);
    let omitted = run(&base, None);
    let identity = run(&base, Some((1.0, 1.0)));
    check_attribution(&omitted, None);
    check_attribution(&identity, Some((1.0, 1.0)));
    check_reference(&omitted, 2.0, 1.0, 1.0);
    check_reference(&identity, 2.0, 1.0, 1.0);
    compare_fields(&omitted, &identity);
}

#[test]
fn g1_native_liquid_phase_conductivity_changes_spatial_heat_and_latent_storage() {
    let base = cards(2.0);
    let unchanged = run(&base, None);
    let coupled = run(&base, Some((1.0, 4.0)));
    check_reference(&unchanged, 2.0, 1.0, 1.0);
    check_reference(&coupled, 2.0, 1.0, 4.0);
    check_attribution(&coupled, Some((1.0, 4.0)));
    let wall = |run: &ThermalRun| {
        values(&run.field, "temperature")
            .into_iter()
            .fold(f64::INFINITY, f64::min)
    };
    assert!(wall(&coupled) > wall(&unchanged) + 0.005);
    let fraction = |run: &ThermalRun| {
        n(run.rows("fine").last().unwrap(), "mean_liquid_mass_fraction")
    };
    assert!(fraction(&coupled) < fraction(&unchanged) - 0.005);
    let storage = |run: &ThermalRun| {
        n(run.transient().get("energy").unwrap(), "stored_change_j")
    };
    assert!(storage(&coupled) < storage(&unchanged) - 0.05);
    // The plateau keeps the maximum fixed while the wall and phase change.
    for solved in [&unchanged, &coupled] {
        close(
            n(solved.rows("fine").last().unwrap(), "final_region_max_k"),
            MELTING_K,
            0.0,
        );
    }
}

#[test]
fn g1_native_constant_phase_multiplier_matches_independently_sourced_conductivity() {
    let base = cards(2.0);
    let doubled = cards(4.0);
    assert_ne!(base.materials()[0].card().to_hex(), doubled.materials()[0].card().to_hex());
    let scaled = run(&base, Some((2.0, 2.0)));
    let sourced = run(&doubled, None);
    check_attribution(&scaled, Some((2.0, 2.0)));
    check_attribution(&sourced, None);
    check_reference(&scaled, 2.0, 2.0, 2.0);
    check_reference(&sourced, 4.0, 1.0, 1.0);
    compare_fields(&scaled, &sourced);
}
