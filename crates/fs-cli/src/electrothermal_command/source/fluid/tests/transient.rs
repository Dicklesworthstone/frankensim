//! The real CLI material, current, face-flux and storage owners share one law.
use super::*;
use crate::electrothermal_command::source::transient::{load, march, render};
use fs_lbm::conjugate::ConvectionScheme;

pub(super) fn request_text() -> &'static str {
    r#"{"schema":"frankensim.electrothermal.transient.v1","time_step_s":0.1,"steps":8,
        "temperature_min_k":290,"temperature_max_k":350,"enthalpy_interpolation_error_j_kg":0.00001,
        "energy_absolute_tolerance_j":0.0000001,"energy_relative_tolerance":0.000000001,
        "schedule":[[1,1,1],[5,0,1]]}"#
}

fn request(text: &str, cells: usize) -> crate::electrothermal_command::source::transient::Request {
    use std::sync::atomic::{AtomicU64, Ordering};
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "frankensim-fluid-time-{}-{stamp}-{}.json",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::write(&path, text).unwrap();
    load(&path, cells).unwrap()
}

fn temperatures(value: &J) -> Vec<f64> {
    value
        .get("temperature_k")
        .unwrap()
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

// Independent three-control-volume backward Euler balance in temperature
// excess above the inlet. Unit density, volume, conductivity, velocity; Cp is
// shared by storage and mass enthalpy flux. Ideal electrode cells deposit the
// circuit's full V^2/R into the middle fluid cell. No native matrix is reused.
fn reference(old: [f64; 3], cp: f64, power: f64) -> ([f64; 3], f64, f64) {
    let rate = cp / 0.1;
    let mut diagonal = [rate + 3.0 + cp, rate + 2.0 + cp, rate + 1.0 + cp];
    let mut rhs = [rate * old[0], rate * old[1] + power, rate * old[2]];
    for row in 1..3 {
        let factor = -(1.0 + cp) / diagonal[row - 1];
        diagonal[row] += factor;
        rhs[row] -= factor * rhs[row - 1];
    }
    let mut next = [0.0; 3];
    next[2] = rhs[2] / diagonal[2];
    next[1] = (rhs[1] + next[2]) / diagonal[1];
    next[0] = (rhs[0] + next[1]) / diagonal[0];
    let stored = cp * next.iter().zip(old).map(|(a, b)| a - b).sum::<f64>();
    let outflow = 0.1 * (2.0 * next[0] + cp * next[2]);
    (next, stored, outflow)
}

#[test]
fn g1_g3_sourced_fluid_capacity_and_resistivity_reach_the_same_transient_energy_balance() {
    let mut final_peaks = Vec::new();
    for (cp, resistivity) in [(1.0, 2.0), (4.0, 2.0), (1.0, 4.0)] {
        let mut p = problem(cp, resistivity, 1.0, 1.0, true);
        p.energy.scheme = ConvectionScheme::Upwind;
        let initial = p.initial_temperature.clone();
        let bytes = p.materials[0].pack.to_bytes();
        let flow = p.fluid.as_ref().unwrap().flow.clone();
        let r = request(request_text(), 3);
        let result = march(&p, &r, &CancelGate::new()).unwrap();
        let output = render(&p, &result, true);
        let value = J::parse(&output).unwrap();
        let rows = value.get("records").unwrap().as_array().unwrap();
        let mut theta = [0.0; 3];
        for (index, row) in rows.iter().enumerate() {
            let power = if index < 4 { 9.0 / resistivity } else { 0.0 };
            let (next, stored, outflow) = reference(theta, cp, power);
            assert!((stored + outflow - 0.1 * power).abs() < 1e-14);
            assert!((row.f64_field("stored_energy_change_j").unwrap() - stored).abs() < 2e-7);
            assert!((row.f64_field("physical_outflow_j").unwrap() - outflow).abs() < 2e-7);
            assert!((row.f64_field("source_j").unwrap() - 0.1 * power).abs() < 1e-8);
            assert!(row.f64_field("closure_j").unwrap().abs() < 1e-7);
            theta = next;
        }
        let actual = temperatures(&value);
        for (cell, expected) in theta.iter().enumerate() {
            assert!((actual[cell] - 300.0 - expected).abs() < 2e-7);
        }
        assert_eq!(
            value
                .get("fluid_reference")
                .unwrap()
                .f64_field("specific_heat_j_kg_k"),
            Some(cp)
        );
        assert_eq!(p.initial_temperature, initial);
        assert_eq!(p.materials[0].pack.to_bytes(), bytes);
        assert_eq!(p.fluid.as_ref().unwrap().flow, flow);
        assert_eq!(
            output,
            render(&p, &march(&p, &r, &CancelGate::new()).unwrap(), true)
        );
        final_peaks.push(actual.iter().copied().fold(f64::NEG_INFINITY, f64::max));
    }
    assert!(final_peaks[1] < final_peaks[0]);
    assert!(final_peaks[2] < final_peaks[0]);
}

#[test]
fn g3_solid_capacity_changes_the_adjacent_advecting_fluid_through_shared_heat_storage() {
    let mut fluid_temperatures = Vec::new();
    for solid_cp in [3.0, 9.0] {
        let mut p = problem(1.0, 2.0, 1.0, 1.0, false);
        p.materials
            .push(material_source(solid_cp, 2.0, 1.0, false, "solid", |_| {}));
        p.material_indices = vec![0, 0, 0, 1, 1, 1];
        p.domain =
            VoxelDomain::from_fn(3, 2, 1, 1.0, |at| Voxel::Solid(u16::from(at[1] > 1.0))).unwrap();
        p.fluid = load_fluid(
            &input(1.0),
            &mut p.domain,
            &p.material_indices,
            &p.materials,
        )
        .unwrap();
        for electrode in &mut p.electric.electrodes {
            electrode.hi[1] = 2;
        }
        p.initial_temperature = vec![300.0; 6];
        let r = request(request_text(), 6);
        let result = march(&p, &r, &CancelGate::new()).unwrap();
        let value = J::parse(&render(&p, &result, true)).unwrap();
        let t = temperatures(&value);
        let stored: f64 = t
            .iter()
            .enumerate()
            .map(|(cell, t)| {
                if cell < 3 {
                    t - 300.0
                } else {
                    solid_cp * (t - 300.0)
                }
            })
            .sum();
        let mut net_input = 0.0;
        for row in value.get("records").unwrap().as_array().unwrap() {
            net_input +=
                row.f64_field("source_j").unwrap() - row.f64_field("physical_outflow_j").unwrap();
            assert!(row.f64_field("closure_j").unwrap().abs() < 1e-7);
        }
        assert!((stored - net_input).abs() < 8e-7);
        assert!((0..3).all(|cell| p.domain.is_fluid(cell)));
        assert!((3..6).all(|cell| !p.domain.is_fluid(cell)));
        fluid_temperatures.push(t[1]);
    }
    assert!(fluid_temperatures[1] < fluid_temperatures[0]);
}

#[test]
fn g0_g4_interior_property_variation_and_late_storage_escape_refuse_without_mutation() {
    for property in [DENSITY_PROPERTY] {
        let mut p = problem(1.0, 2.0, 1.0, 1.0, false);
        p.materials[0] = material_source(1.0, 2.0, 1.0, false, "liquid", |claim| {
            if claim.key.name() == property {
                claim.value = PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: fs_qty::Temperature::DIMS,
                    knots: vec![(290.0, 1.0), (300.0, 1.0), (320.0, 2.0), (350.0, 1.0)],
                    dims: claim.key.dims(),
                };
                claim.interpolation = InterpolationPolicy::LinearInside;
            }
        });
        p.fluid = load_fluid(
            &input(1.0),
            &mut p.domain,
            &p.material_indices,
            &p.materials,
        )
        .unwrap();
        let before = p.initial_temperature.clone();
        let error = march(&p, &request(request_text(), 3), &CancelGate::new())
            .err()
            .unwrap();
        assert!(error.to_string().contains("must remain constant"));
        assert_eq!(p.initial_temperature, before);
    }
    let mut p = problem(1.0, 2.0, 1.0, 1.0, false);
    let before = p.initial_temperature.clone();
    let tight = request_text().replace(
        "\"temperature_max_k\":350",
        "\"temperature_max_k\":300.00001",
    );
    assert!(march(&p, &request(&tight, 3), &CancelGate::new()).is_err());
    assert_eq!(p.initial_temperature, before);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        march(&p, &request(request_text(), 3), &gate),
        Err(ChtError::Cancelled)
    ));
    assert!(march(&p, &request(request_text(), 3), &CancelGate::new()).is_ok());
    p.thermal.faces[0] = ThermalFace::Inflow { temperature: 351.0 };
    assert!(
        march(&p, &request(request_text(), 3), &CancelGate::new())
            .err()
            .unwrap()
            .to_string()
            .contains("inflow and backflow temperatures")
    );
    assert_eq!(p.initial_temperature, before);
}

mod caloric;
