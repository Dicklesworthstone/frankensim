//! Authored compact storage alongside original source-caloric cell storage.

use super::*;
use fs_lbm::Face3;
use fs_lbm::conjugate::CompactComponent;

fn fixture(sourced_board: bool) -> (VoxelDomain, ThermalSetup, Vec<f64>) {
    let domain =
        VoxelDomain::graded_from_fn([vec![0.25, 0.75], vec![1.0], vec![1.0, 0.2, 1.0]], |p| {
            if p[2] > 1.0 && p[2] < 1.2 {
                // The alias index is intentionally beyond the source array.
                Voxel::Solid(u16::from(sourced_board))
            } else if sourced_board && p[2] < 1.0 {
                Voxel::Solid(0)
            } else {
                Voxel::Fluid
            }
        })
        .unwrap();
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.compact_components.push(CompactComponent {
        lo: [0, 0, 1],
        hi: [2, 1, 2],
        board_face: Face3::ZMin,
        power_w: 12.0,
        junction_to_case_k_w: 1.0,
        junction_to_board_k_w: 1.0,
        heat_capacity_j_k: Some(6.0),
    });
    let mut initial = vec![300.0; domain.cell_count()];
    for x in 0..2 {
        initial[domain.index(x, 0, 1)] = 320.0;
    }
    (domain, thermal, initial)
}

fn pulse(time: f64) -> f64 {
    if time <= 1.0 { 1.0 } else { 0.0 }
}

fn fluid_pair_step(tf: f64, tj: f64, power: f64) -> (f64, f64) {
    // Each outside layer has mass 3 kg, Cp=2+.002*(T-300), and G=2/3.
    // At dt=1, write dh=cp_old*dT+.001*dT². Conservation gives
    // Tj_new=Tj_old+P/6-dh; one fluid layer gives 3*dh=G*(Tj_new-Tf_new).
    // Eliminating Tj_new produces this independent positive quadratic root.
    let g = 2.0 / 3.0;
    let cp = 2.0 + 0.002 * (tf - 300.0);
    let a = 0.001 * (3.0 + g);
    let b = 3.0 * cp + g * (cp + 1.0);
    let drive = g * (tj - tf + power / 6.0);
    let delta = 2.0 * drive / (b + (b * b + 4.0 * a * drive).sqrt());
    let dh = cp * delta + 0.001 * delta * delta;
    (tf + delta, tj + power / 6.0 - dh)
}

#[test]
fn source_caloric_compact_pulse_counts_one_junction_and_schedules_its_power() {
    let fluid_card = card(0.002, 0.0);
    let source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &source, [290.0, 400.0]);
    let reference = transport(&source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let (domain, thermal, initial) = fixture(false);

    let run = bound
        .march_natural_convection_with_enthalpy(
            &domain,
            &[], // Compact aliases need no fictitious bulk material slot.
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &thermal,
            &energy(),
            &initial,
            bounds(),
            &caloric,
            controls(),
            &UnsteadyConfig::new(1.0, 3),
            |_| 1.0,
            pulse,
            None,
            &CancelGate::new(),
            |_, t, _| transport(&source, t),
        )
        .unwrap();

    let (mut tf, mut tj) = (300.0, 320.0);
    for heat in &run.convection.march.energy_records {
        let power = 12.0 * pulse(heat.time_s);
        (tf, tj) = fluid_pair_step(tf, tj, power);
        close(heat.max_solid_temperature_k, tj, 1e-10);
        close(heat.stored_energy_change_j, power, 1e-9);
        close(heat.source_j, power, 1e-12);
        close(heat.boundary_outflow_j, 0.0, 1e-12);
        close(heat.closure_j, 0.0, 2e-8);
    }
    let temperature = &run.convection.march.temperature;
    for (cell, &actual) in temperature.iter().enumerate() {
        close(actual, if domain.is_fluid(cell) { tf } else { tj }, 1e-10);
    }
    assert_eq!(
        temperature[domain.index(0, 0, 1)],
        temperature[domain.index(1, 0, 1)]
    );
    let fluid_energy: f64 = temperature
        .iter()
        .enumerate()
        .filter(|(cell, _)| domain.is_fluid(*cell))
        .map(|(cell, &t)| 3.0 * domain.volume(cell) * specific_energy(t, 0.002))
        .sum();
    // Initially 120 J resides in the junction; the pulse adds exactly 12 J.
    close(fluid_energy + 6.0 * (tj - 300.0), 132.0, 1e-9);
    assert_eq!(run.convection.endpoint.len(), 4);
    assert!(
        run.convection
            .endpoint
            .iter()
            .all(|sample| domain.is_fluid(sample.cell))
    );
}

#[test]
fn source_caloric_fixed_reservoir_jump_includes_junction_storage_once() {
    let fluid_card = card(0.002, 0.0);
    let source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &source, [290.0, 400.0]);
    let reference = transport(&source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let (domain, mut thermal, initial) = fixture(false);
    thermal.compact_components[0].power_w = 0.0;
    thermal.fixed_temperature = (0..domain.cell_count())
        .filter(|&cell| domain.is_fluid(cell))
        .map(|cell| (cell, 310.0))
        .collect();

    let run = bound
        .march_natural_convection_with_enthalpy(
            &domain,
            &[],
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &thermal,
            &energy(),
            &initial,
            bounds(),
            &caloric,
            controls(),
            &UnsteadyConfig::new(1.0, 1),
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| transport(&source, t),
        )
        .unwrap();

    // Both branches together give G=4/3 W/K. The thermostat supplies the
    // exact fluid enthalpy jump, while the junction cools by its own Cj law.
    let tj = (6.0 * 320.0 + (4.0 / 3.0) * 310.0) / (6.0 + 4.0 / 3.0);
    let fluid_storage = 6.0 * specific_energy(310.0, 0.002); // 120.6 J
    let junction_storage = 6.0 * (tj - 320.0);
    let stored = fluid_storage + junction_storage;
    for (cell, &actual) in run.convection.march.temperature.iter().enumerate() {
        close(
            actual,
            if domain.is_fluid(cell) { 310.0 } else { tj },
            1e-10,
        );
    }
    let heat = &run.convection.march.energy_records[0];
    close(heat.stored_energy_change_j, stored, 1e-9);
    close(heat.source_j, 0.0, 1e-12);
    close(heat.boundary_outflow_j, -stored, 1e-9);
    close(heat.closure_j, 0.0, 2e-8);
    assert!(run.records[0].maximum_storage_defect_w <= controls().maximum_storage_defect_w);
}

#[test]
fn compact_caloric_domain_refusal_and_mid_march_cancellation_publish_no_endpoint() {
    let fluid_card = card(0.002, 0.0);
    let source = law(&fluid_card, [290.0, 400.0]);
    let narrow = source_caloric(&fluid_card, &source, [290.0, 300.5]);
    let broad = source_caloric(&fluid_card, &source, [290.0, 400.0]);
    let reference = transport(&source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let (domain, thermal, initial) = fixture(false);
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);

    let refused = bound.march_natural_convection_with_enthalpy(
        &domain,
        &[],
        &flow,
        &thermal,
        &energy(),
        &initial,
        bounds(),
        &narrow,
        controls(),
        &UnsteadyConfig::new(1.0, 1),
        |_| 1.0,
        pulse,
        None,
        &CancelGate::new(),
        |_, t, _| transport(&source, t),
    );
    assert!(matches!(
        refused,
        Err(MaterialBuoyancyError::Caloric { temperature_k, .. })
            if temperature_k > 300.5
    ));

    let gate = CancelGate::new();
    let cancelled = bound.march_natural_convection_with_enthalpy(
        &domain,
        &[],
        &flow,
        &thermal,
        &energy(),
        &initial,
        bounds(),
        &broad,
        controls(),
        &UnsteadyConfig::new(1.0, 2),
        |_| 1.0,
        |time| {
            // The first physical step succeeds; cancellation arrives while
            // preparing the second step's scheduled component source.
            if time > 1.0 {
                gate.request();
            }
            pulse(time)
        },
        None,
        &gate,
        |_, t, _| transport(&source, t),
    );
    assert!(matches!(
        cancelled,
        Err(MaterialBuoyancyError::Solver(ChtError::Cancelled))
    ));

    let retry = bound
        .march_natural_convection_with_enthalpy(
            &domain,
            &[],
            &flow,
            &thermal,
            &energy(),
            &initial,
            bounds(),
            &broad,
            controls(),
            &UnsteadyConfig::new(1.0, 2),
            |_| 1.0,
            pulse,
            None,
            &CancelGate::new(),
            |_, t, _| transport(&source, t),
        )
        .unwrap();
    let (tf, tj) = fluid_pair_step(300.0, 320.0, 12.0);
    let (tf, tj) = fluid_pair_step(tf, tj, 0.0);
    for (cell, &actual) in retry.convection.march.temperature.iter().enumerate() {
        close(actual, if domain.is_fluid(cell) { tf } else { tj }, 1e-10);
    }
    assert_eq!(retry.convection.march.energy_records.len(), 2);
}

#[test]
fn joint_caloric_compact_aliases_have_no_bulk_source_or_source_solid_energy() {
    let fluid_card = card(0.002, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let reference = transport(&fluid_source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.0);
    let solid = law(&solid_card, [290.0, 400.0]);
    let (domain, thermal, initial) = fixture(true);

    let run = bound
        .march_natural_convection_with_enthalpy_and_material_solids(
            &domain,
            &[solid], // Ordinary material 0 only; compact material 1 is omitted.
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &thermal,
            &energy(),
            &initial,
            bounds(),
            &caloric,
            controls(),
            SolidThermalConfig::default(),
            &UnsteadyConfig::new(1.0, 1),
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| transport(&fluid_source, t),
        )
        .unwrap();

    // Each outside layer independently satisfies
    // 3*(2*theta + slope*theta²/2) = (2/3)*(theta_j-theta).
    // Substitute those exact quadratic roots into the junction BE equation.
    let g = 2.0 / 3.0;
    let layer_rise = |junction: f64, slope: f64| {
        let a = 1.5 * slope;
        let b = 6.0 + g;
        let drive = g * junction;
        2.0 * drive / (b + (b * b + 4.0 * a * drive).sqrt())
    };
    let (mut lo, mut hi) = (0.0, 22.0);
    for _ in 0..90 {
        let j = 0.5 * (lo + hi);
        let s = layer_rise(j, 0.2);
        let f = layer_rise(j, 0.002);
        if 6.0 * (j - 20.0) + g * (2.0 * j - s - f) > 12.0 {
            hi = j;
        } else {
            lo = j;
        }
    }
    let j = 0.5 * (lo + hi);
    let s = layer_rise(j, 0.2);
    let f = layer_rise(j, 0.002);
    let temperature = &run.fluid.convection.march.temperature;
    for (cell, &actual) in temperature.iter().enumerate() {
        close(actual, 300.0 + [s, j, f][domain.coords(cell)[2]], 1e-10);
    }

    let solid_storage = 3.0 * specific_energy(300.0 + s, 0.2);
    let fluid_storage = 3.0 * specific_energy(300.0 + f, 0.002);
    let junction_storage = 6.0 * (j - 20.0);
    close(
        solid_storage + fluid_storage + junction_storage,
        12.0,
        1e-12,
    );
    close(
        run.solid_records[0].solid_stored_energy_change_j,
        solid_storage,
        1e-9,
    );
    let heat = &run.fluid.convection.march.energy_records[0];
    close(heat.stored_energy_change_j, 12.0, 1e-9);
    close(heat.source_j, 12.0, 1e-12);
    close(heat.boundary_outflow_j, 0.0, 1e-12);
    close(heat.closure_j, 0.0, 2e-8);

    assert_eq!(run.solid_endpoint.len(), 2);
    assert!(
        run.solid_endpoint
            .iter()
            .all(|sample| sample.material_index == 0 && domain.coords(sample.cell)[2] == 0)
    );
    for x in 0..2 {
        assert_eq!(run.solid_energy_j[domain.index(x, 0, 1)], 0.0);
    }
    // Ordinary solid energy uses its declared 290 K datum. The component
    // contributes neither a source receipt nor a second bulk energy term.
    close(run.solid_energy_j.iter().sum(), 60.0 + solid_storage, 1e-9);
}

#[test]
fn source_solid_compact_storage_preserves_fixed_fluid_cp_and_solid_subtotal() {
    let fluid_card = card(0.0, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let reference = transport(&fluid_source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.0);
    let solid = law(&solid_card, [290.0, 400.0]);
    let (domain, thermal, initial) = fixture(true);

    let run = bound
        .march_natural_convection_with_material_solids(
            &domain,
            &[solid], // Ordinary material 0 only; compact material 1 is omitted.
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &thermal,
            &energy(),
            &initial,
            bounds(),
            &UnsteadyConfig::new(1.0, 1),
            SolidThermalConfig::default(),
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| transport(&fluid_source, t),
        )
        .unwrap();

    // Each outside layer has reference mass 3 kg and branch G=2/3 W/K.
    // Solid: .3*theta_s²+(6+G)*theta_s=G*theta_j.
    // Constant-Cp fluid: (6+G)*theta_f=G*theta_j, hence theta_f=theta_j/10.
    // Junction: 6*(theta_j-20)+G*(2*theta_j-theta_s-theta_f)=12.
    let g = 2.0 / 3.0;
    let solid_rise = |junction: f64| {
        let b = 6.0 + g;
        let drive = g * junction;
        2.0 * drive / (b + (b * b + 1.2 * drive).sqrt())
    };
    let (mut lo, mut hi) = (0.0, 22.0);
    for _ in 0..90 {
        let j = 0.5 * (lo + hi);
        if 6.0 * (j - 20.0) + g * (2.0 * j - solid_rise(j) - j / 10.0) > 12.0 {
            hi = j;
        } else {
            lo = j;
        }
    }
    let j = 0.5 * (lo + hi);
    let s = solid_rise(j);
    let f = j / 10.0;
    for (cell, &actual) in run.fluid.march.temperature.iter().enumerate() {
        close(actual, 300.0 + [s, j, f][domain.coords(cell)[2]], 1e-10);
    }

    let solid_storage = 3.0 * specific_energy(300.0 + s, 0.2);
    close(
        run.solid_records[0].solid_stored_energy_change_j,
        solid_storage,
        1e-9,
    );
    close(solid_storage + 6.0 * f + 6.0 * (j - 20.0), 12.0, 1e-12);
    let heat = &run.fluid.march.energy_records[0];
    close(heat.stored_energy_change_j, 12.0, 1e-9);
    close(heat.source_j, 12.0, 1e-12);
    close(heat.boundary_outflow_j, 0.0, 1e-12);
    close(heat.closure_j, 0.0, 2e-8);

    assert_eq!(run.solid_endpoint.len(), 2);
    assert!(
        run.solid_endpoint
            .iter()
            .all(|sample| sample.material_index == 0 && domain.coords(sample.cell)[2] == 0)
    );
    for x in 0..2 {
        assert_eq!(run.solid_energy_j[domain.index(x, 0, 1)], 0.0);
    }
}
