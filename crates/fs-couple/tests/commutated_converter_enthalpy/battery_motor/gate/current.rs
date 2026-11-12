//! G0/G1: typed original triangles and independent current/power/caloric identities.
//! G3: full source support, passive current, thermal authority and retained history.
//! Manufactured source fixtures establish no measured-device qualification.
use super::*;
use fs_couple::circuit_enthalpy::converter::commutation::{
    GateChannelMaterial, GateCurrentSurfaceBinding, GateCurrentSurfaceMaterial,
    GateCurrentThermalScale,
};
use fs_qty::{Current, Dims};

type CurrentDrive = MaterialGateBatteryMotorDrive<GateCurrentSurfaceMaterial>;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Defect {
    None,
    MissingZeroMesh,
    Hole,
    ActiveVertex,
    NonzeroZeroRow,
    NarrowDomain,
    ScaleDrainBias,
    ScaleGateBias,
}
struct CurrentFixture {
    gate: GateFixture,
    current: [PropertyKey; 4],
    scale: PropertyKey,
}
fn current_fixture(miller: f64, slope: f64, defect: Defect) -> CurrentFixture {
    let mut gate = gate_fixture(miller, false, 0.0, false, false);
    let current = core::array::from_fn(|i| {
        PropertyKey::new(format!("channel_current_surface_{i}"), Current::DIMS)
    });
    let scale = PropertyKey::new("channel_temperature_multiplier", Dims::NONE);
    let mut claims = ClaimSet::new();
    for (_, claim) in gate.base.card.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let thermal = || {
        ValidityDomain::unconstrained().with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            299.0,
            350.0,
        )
    };
    let mut add = |key: &PropertyKey, value, validity| {
        claims.insert_claim(PropertyClaim{key:key.clone(),value,validity,
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,observations:vec![],
            provenance:Provenance{source:"manufactured original triangular channel-only current; not qualified transistor data".into(),license:"CC0-1.0".into(),artifact:None}}).unwrap();
    };
    for (device, key) in current.iter().enumerate() {
        let drain: Vec<f64> = if defect == Defect::MissingZeroMesh {
            vec![-100.0, -1.0, 1.0, 100.0]
        } else {
            vec![-100.0, -1.0, 0.0, 1.0, 100.0]
        };
        let gate_v: [f64; 4] = [
            -4.0,
            0.0,
            1.0,
            if defect == Defect::NarrowDomain {
                3.0
            } else {
                4.0
            },
        ];
        let reversed = device % 2 == 1;
        let mut samples = Vec::new();
        for &vg in &gate_v {
            for &vd in &drain {
                let mut amps = vd.signum()
                    * vd.abs().min(1.0)
                    * vg.max(0.0)
                    * if vd < 0.0 { 0.4 } else { 1.0 };
                if defect == Defect::ActiveVertex && vd == 1.0 && vg == 4.0 {
                    amps = -1.0;
                }
                if defect == Defect::NonzeroZeroRow && vd == 0.0 && vg == 1.0 {
                    amps = 1e-6;
                }
                samples.push((if reversed { [vg, vd] } else { [vd, vg] }, amps));
            }
        }
        let mut triangles = Vec::new();
        for row in 0..gate_v.len() - 1 {
            for col in 0..drain.len() - 1 {
                // This interior hole retains all four vertices incident to neighboring cells.
                if defect == Defect::Hole && row == 1 && col == 2 {
                    continue;
                }
                let a = (row * drain.len() + col) as u32;
                let b = a + 1;
                let d = a + drain.len() as u32;
                let e = d + 1;
                for mut triangle in [[a, b, d], [b, e, d]] {
                    if reversed {
                        triangle.swap(1, 2);
                    }
                    triangles.push(triangle);
                }
            }
        }
        add(
            key,
            PropertyValue::Surface2D {
                axes: if reversed {
                    ["Vgs".into(), "Vds".into()]
                } else {
                    ["Vds".into(), "Vgs".into()]
                },
                axis_dims: [Voltage::DIMS; 2],
                samples,
                triangles,
                dims: Current::DIMS,
            },
            thermal()
                .with_quantity(
                    "Vds",
                    QuantitySpec::dimensional(Voltage::DIMS),
                    -100.0,
                    100.0,
                )
                .with_quantity(
                    "Vgs",
                    QuantitySpec::dimensional(Voltage::DIMS),
                    -4.0,
                    if defect == Defect::NarrowDomain {
                        3.0
                    } else {
                        4.0
                    },
                ),
        );
    }
    let mut domain = thermal();
    if defect == Defect::ScaleDrainBias {
        domain = domain.with_quantity("Vds", QuantitySpec::dimensional(Voltage::DIMS), 0.0, 0.0);
    }
    if defect == Defect::ScaleGateBias {
        domain = domain.with_quantity("Vgs", QuantitySpec::dimensional(Voltage::DIMS), -4.0, -4.0);
    }
    add(
        &scale,
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(299.0, 1.0 - slope), (350.0, 1.0 + 50.0 * slope)],
            dims: Dims::NONE,
        },
        domain,
    );
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-device".into(),
            phase: "solid".into(),
            process: "original signed triangular current and explicit thermal multiplier".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &common::point(299.0),
            upper: &common::point(350.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    gate.base.card = card;
    gate.base.caloric = caloric;
    CurrentFixture {
        gate,
        current,
        scale,
    }
}
impl CurrentFixture {
    fn bindings(&self) -> [GateCurrentSurfaceBinding<'_>; 4] {
        core::array::from_fn(|i| GateCurrentSurfaceBinding {
            cell: IDS[i + 1],
            key: &self.current[i],
            drain_voltage_axis: "Vds",
            gate_voltage_axis: "Vgs",
            drain_voltage_range_v: [-20.0, 20.0],
            gate_voltage_range_v: [-4.0, 4.0],
            thermal_scale: Some(GateCurrentThermalScale {
                key: &self.scale,
                selection: MaterialPropertySelection::SingleClaimOnly,
            }),
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
    }
}
fn material_bindings<'a>(
    q: &'a [DeviceChargeBinding<'a>; 4],
    cgs: &'a [DeviceCapacitanceBinding<'a>; 4],
    cgd: &'a [DeviceCapacitanceBinding<'a>; 4],
    current: &'a [GateCurrentSurfaceBinding<'a>; 4],
    driver: &'a [ConverterThermalChannel<'a>; 8],
) -> GateMaterialBindings<'a, GateCurrentSurfaceBinding<'a>> {
    GateMaterialBindings {
        intrinsic_drain_charge: q,
        gate_source_capacitance: cgs,
        gate_drain_capacitance: cgd,
        channel_conductance: current,
        driver_resistors: driver,
        drivers: drivers(),
    }
}
fn current_network<'a>(
    cells: &'a [MaterialCell<'a>],
    r: &'a [ConverterThermalChannel<'a>; 5],
    d: &'a [DeviceConductionBinding<'a>; 4],
    bindings: GateMaterialBindings<'a, GateCurrentSurfaceBinding<'a>>,
) -> GateBatteryMotorThermalNetwork<'a, GateCurrentSurfaceBinding<'a>> {
    GateBatteryMotorThermalNetwork {
        drive: CommutatedBatteryMotorThermalNetwork {
            cells,
            links: &[],
            resistors: r,
            devices: d,
            bearing_cell: IDS[0],
        },
        gate: bindings,
        maximum_cells_per_chord: 4096,
    }
}
fn current_initial(steady: &SteadyDrive, p: &GateCommutationParameters) -> [f64; 15] {
    let mut state = initial(steady, p);
    state[3..12].copy_from_slice(
        &p.state(0.3, [12.0, 6.0], [10.0, 7.5], [2.0, 0.2, 2.0, 0.2])
            .unwrap(),
    );
    state
}
fn command(steady: &SteadyDrive) -> BatteryMotorSwitchCommand {
    BatteryMotorSwitchCommand {
        auxiliary_load: None,
        ..steady.command
    }
}
fn approx(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-10, "{a:.17e} != {b:.17e}");
}

#[test]
fn original_triangles_retain_signed_saturation_cutoff_zero_mesh_and_named_partials() {
    let f = current_fixture(0.05, 0.2, Defect::None);
    let thermal = f.gate.thermal();
    let context = gate_context();
    let cells = f.gate.cells(&context);
    let r = drive_resistors(&f.gate.base);
    let d = f.gate.base.devices();
    let q = f.gate.charges();
    let cgs = f.gate.caps(false);
    let cgd = f.gate.caps(true);
    let current = f.bindings();
    let driver = f.gate.driver_resistors();
    let network = current_network(
        &cells,
        &r,
        &d,
        material_bindings(&q, &cgs, &cgd, &current, &driver),
    );
    let steady = steady_drive(false);
    let p = authored(0.05, false);
    let state = current_initial(&steady, &p);
    let bound: CurrentDrive = network
        .bind(&thermal, &steady.native, &p, state, || false)
        .unwrap();
    for (i, source) in bound
        .gate_material()
        .channel_conductance()
        .iter()
        .enumerate()
    {
        assert_eq!(source.source_claim().key, f.current[i]);
        assert_eq!(source.drain_voltage_axis(), "Vds");
        assert_eq!(source.gate_voltage_axis(), "Vgs");
        assert_eq!(
            source.response().axes()[0],
            if i % 2 == 0 { "Vds" } else { "Vgs" }
        );
        near(source.thermal_scale(), 1.0);
        assert!(source.thermal_scale_source().is_some());
        let law = source.native_response();
        for (vd, vg, expected) in [
            (2.0, 2.0, [2.0, 0.0, 1.0]),
            (-2.0, 2.0, [-0.8, 0.0, -0.4]),
            (0.75, 0.75, [0.5, 1.0, 1.0]),
            (2.0, -1.0, [0.0, 0.0, 0.0]),
        ] {
            let response = law.evaluate(vd, vg).unwrap();
            for (a, b) in [
                response.current_a,
                response.drain_derivative_s,
                response.gate_derivative_s,
            ]
            .into_iter()
            .zip(expected)
            {
                approx(a, b);
            }
            let e = 1e-5;
            approx(
                (law.evaluate(vd + e, vg).unwrap().current_a
                    - law.evaluate(vd - e, vg).unwrap().current_a)
                    / (2.0 * e),
                response.drain_derivative_s,
            );
            approx(
                (law.evaluate(vd, vg + e).unwrap().current_a
                    - law.evaluate(vd, vg - e).unwrap().current_a)
                    / (2.0 * e),
                response.gate_derivative_s,
            );
        }
        for vg in [-4.0, -0.2, 0.0, 0.7, 2.0, 4.0] {
            assert_eq!(law.evaluate(0.0, vg).unwrap().current_a, 0.0);
        }
        assert!(law.evaluate(20.01, 1.0).is_err());
        assert!(law.evaluate(-20.01, 1.0).is_err());
        assert!(law.evaluate(1.0, 4.01).is_err());
        let PropertyValue::Surface2D {
            samples, triangles, ..
        } = &source.source_claim().value
        else {
            panic!("original surface lost")
        };
        assert_eq!(source.response().samples(), samples.as_slice());
        assert_eq!(source.response().triangles(), triangles.as_slice());
    }
}

#[test]
fn accepted_temperature_scales_actual_channel_power_and_retains_original_source_history() {
    let f = current_fixture(0.05, 0.2, Defect::None);
    let thermal = f.gate.thermal();
    let context = gate_context();
    let cells = f.gate.cells(&context);
    let r = drive_resistors(&f.gate.base);
    let d = f.gate.base.devices();
    let q = f.gate.charges();
    let cgs = f.gate.caps(false);
    let cgd = f.gate.caps(true);
    let current = f.bindings();
    let driver = f.gate.driver_resistors();
    let bindings = material_bindings(&q, &cgs, &cgd, &current, &driver);
    let network = current_network(&cells, &r, &d, bindings);
    let steady = steady_drive(false);
    let p = authored(0.05, false);
    let state = current_initial(&steady, &p);
    let drive = network
        .bind(&thermal, &steady.native, &p, state, || false)
        .unwrap();
    let h = 1.0 / 1024.0;
    let deck = interval(h);
    let cmd = command(&steady);
    let next = network
        .step_candidate(
            &thermal,
            h,
            &drive,
            &deck,
            &[cmd],
            1,
            None,
            power_controls(),
            cycle_controls(),
            drive_controls(11),
            &[0.0; 11],
            || false,
        )
        .unwrap();
    let step = &next.report().electrical;
    let gate = step.converter.gate.as_ref().unwrap();
    let (d0, g0) = voltages(&state, 0.05);
    let (d1, g1) = voltages(next.drive().state(), 0.05);
    for i in 0..4 {
        let t = next.thermal().temperature_k(IDS[i + 1]).unwrap();
        let factor = 1.0 + 0.2 * (t - 300.0);
        let source = &next.drive().gate_material().channel_conductance()[i];
        near(source.thermal_scale(), factor);
        assert!(factor > 1.0);
        let vd = f64::midpoint(d0[i], d1[i]);
        let vg = f64::midpoint(g0[i], g1[i]);
        assert!(vd.abs() > 1.0 && vg > 0.0);
        let amps = vd.signum() * if vd < 0.0 { 0.4 } else { 1.0 } * vg * factor;
        near(
            source.native_response().evaluate(vd, vg).unwrap().current_a,
            amps,
        );
        near(step.converter.channel_heat_j[i], h * vd * amps);
        near(
            step.converter.diode_heat_j[i],
            h * (-vd) * ((-vd - 0.5) / 0.2).max(0.0),
        );
        near(
            t - 300.0,
            step.converter.channel_heat_j[i] + step.converter.diode_heat_j[i],
        );
        let old = &drive.gate_material().channel_conductance()[i];
        assert_eq!(source.source_claim(), old.source_claim());
        assert_eq!(source.response().samples(), old.response().samples());
        assert_eq!(source.response().triangles(), old.response().triangles());
        assert_eq!(
            source.source().card_identity(),
            old.source().card_identity()
        );
        near(
            source
                .source()
                .query_point()
                .iter()
                .find(|(axis, _)| axis == "T")
                .unwrap()
                .1,
            t,
        );
        near(
            source
                .thermal_scale_source()
                .unwrap()
                .query_point()
                .iter()
                .find(|(axis, _)| axis == "T")
                .unwrap()
                .1,
            t,
        );
        let driver = next.drive().gate_material().drivers()[i];
        let (u, rr) = if deck[0].gates.on[i] {
            (driver.on_voltage_v, driver.on_resistance_ohm)
        } else {
            (driver.off_voltage_v, driver.off_resistance_ohm)
        };
        let ig = (u - vg) / rr;
        near(gate.source_work_j[i], h * u * ig);
        near(gate.resistor_heat_j[i], h * rr * ig * ig);
        near(
            next.thermal().temperature_k(DRIVER[i]).unwrap() - 300.0,
            gate.resistor_heat_j[i],
        );
    }
    for (i, heat) in per_body(step.heat_joules(), 0.0).into_iter().enumerate() {
        near(next.thermal().temperature_k(IDS[i]).unwrap() - 300.0, heat);
    }
    let actual = next.drive().gate_parameters(4096);
    let mut replay = PreparedGateBatteryMotorCycle::new_gate_controlled(
        steady_drive(false).native,
        actual,
        &deck,
        1,
        power_controls(),
        cycle_controls(),
        || false,
    )
    .unwrap();
    let mut end = [0.0; 15];
    let raw = replay
        .step_interval_commands_controlled(&state, &[cmd], &mut end, || false)
        .unwrap();
    for (a, b) in next.drive().state().iter().zip(end) {
        near(*a, b);
    }
    for (a, b) in step
        .all_heat_joules()
        .into_iter()
        .zip(raw.all_heat_joules())
    {
        near(a, b);
    }
    near(
        next.report().deposited_joule_j,
        step.all_heat_joules().iter().sum(),
    );
    near(next.report().total_energy_residual_j, 0.0);
    let continued = network
        .step_candidate(
            next.thermal(),
            2.0 * h,
            next.drive(),
            &deck,
            &[cmd],
            1,
            None,
            power_controls(),
            cycle_controls(),
            drive_controls(11),
            &[0.0; 11],
            || false,
        )
        .unwrap();
    for (a, b) in next
        .drive()
        .gate_material()
        .channel_conductance()
        .iter()
        .zip(continued.drive().gate_material().channel_conductance())
    {
        assert_eq!(a.source_claim(), b.source_claim());
        assert_eq!(a.response().samples(), b.response().samples());
        assert_eq!(a.source().card_identity(), b.source().card_identity());
    }
    let changed = context.clone().with("lot", 8.0).unwrap();
    let changed_cells = f.gate.cells(&changed);
    let wrong = current_network(&changed_cells, &r, &d, bindings);
    assert!(
        wrong
            .step_candidate(
                next.thermal(),
                2.0 * h,
                next.drive(),
                &deck,
                &[cmd],
                1,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(11),
                &[0.0; 11],
                || false
            )
            .is_err()
    );
    assert_eq!(drive.state(), &state);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn missing_surface_support_active_vertices_and_bias_conditioned_multiplier_are_refused() {
    for defect in [
        Defect::MissingZeroMesh,
        Defect::Hole,
        Defect::ActiveVertex,
        Defect::NonzeroZeroRow,
        Defect::NarrowDomain,
        Defect::ScaleDrainBias,
        Defect::ScaleGateBias,
    ] {
        let f = current_fixture(0.0, 0.2, defect);
        let thermal = f.gate.thermal();
        let context = gate_context();
        let cells = f.gate.cells(&context);
        let r = drive_resistors(&f.gate.base);
        let d = f.gate.base.devices();
        let q = f.gate.charges();
        let cgs = f.gate.caps(false);
        let cgd = f.gate.caps(true);
        let current = f.bindings();
        let driver = f.gate.driver_resistors();
        let network = current_network(
            &cells,
            &r,
            &d,
            material_bindings(&q, &cgs, &cgd, &current, &driver),
        );
        let steady = steady_drive(false);
        let p = authored(0.0, false);
        let state = current_initial(&steady, &p);
        assert!(
            network
                .bind(&thermal, &steady.native, &p, state, || false)
                .is_err(),
            "defect {defect:?} admitted"
        );
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn current_surface_thermal_feedback_uses_the_original_generic_heat_assembly() {
    let f = current_fixture(0.05, 0.2, Defect::None);
    let thermal = f.gate.thermal();
    let context = gate_context();
    let cells = f.gate.cells(&context);
    let bulk = MaterialNetwork {
        cells: &cells,
        links: &[],
        resistors: &[],
    };
    let electrical = ElectricalContactNetwork {
        bulk: &bulk,
        contacts: &[],
        electrical_bases: &[],
    };
    let heat = MaterialHeatNetwork {
        electrical: &electrical,
        thermal_contacts: &[],
        radiation: None,
    };
    let cancel = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &cancel,
            arena,
            fs_exec::StreamKey {
                seed: 59,
                kernel_id: 103,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = heat
            .prepare_with_link_conductances(&cx, &thermal, &BTreeMap::new())
            .unwrap();
        let r = drive_resistors(&f.gate.base);
        let d = f.gate.base.devices();
        let q = f.gate.charges();
        let cgs = f.gate.caps(false);
        let cgd = f.gate.caps(true);
        let current = f.bindings();
        let driver = f.gate.driver_resistors();
        let network = GateBatteryMotorHeatNetwork {
            drive: CommutatedBatteryMotorHeatNetwork {
                assembly: &prepared,
                resistors: &r,
                devices: &d,
                bearing_cell: IDS[0],
            },
            gate: material_bindings(&q, &cgs, &cgd, &current, &driver),
            maximum_cells_per_chord: 4096,
        };
        let steady = steady_drive(false);
        let p = authored(0.05, false);
        let state = current_initial(&steady, &p);
        let drive = network
            .bind(&cx, &thermal, &steady.native, &p, state)
            .unwrap();
        let h = 1.0 / 1024.0;
        let next = network
            .step_candidate(
                &cx,
                &thermal,
                h,
                &drive,
                &interval(h),
                &[command(&steady)],
                1,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(11),
                &[0.0; 11],
            )
            .unwrap();
        let step = &next.report().electrical;
        for (i, q) in per_body(step.heat_joules(), 0.0).into_iter().enumerate() {
            near(next.thermal().temperature_k(IDS[i]).unwrap() - 300.0, q);
        }
        for (i, q) in step
            .converter
            .gate
            .as_ref()
            .unwrap()
            .resistor_heat_j
            .into_iter()
            .enumerate()
        {
            near(next.thermal().temperature_k(DRIVER[i]).unwrap() - 300.0, q);
        }
        for (i, source) in next
            .drive()
            .gate_material()
            .channel_conductance()
            .iter()
            .enumerate()
        {
            near(
                source.thermal_scale(),
                1.0 + 0.2 * (next.thermal().temperature_k(IDS[i + 1]).unwrap() - 300.0),
            );
        }
        near(
            next.thermal().energy_j() - thermal.energy_j(),
            step.all_heat_joules().iter().sum(),
        );
        near(
            next.report().deposited_joule_j,
            step.all_heat_joules().iter().sum(),
        );
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(next.report().source.drive.network.material.cells.len(), 11);
    });
}

// Planted negative on the accepted-state check (see the converter's): the drive stays bound at 300 K
// while the accepted heat state is another one at the same time, and every channel current carries an
// explicit thermal multiplier (slope 0.2 per K), so the current-surface material resolved at the
// accepted state differs in its multiplier bits from the bound one. The first attempt at the bound
// state is the positive: bound under SingleClaimOnly, stepped once, accepted.
#[test]
fn current_surface_drive_bound_at_one_state_refuses_a_step_from_another_accepted_state() {
    let f = current_fixture(0.05, 0.2, Defect::None);
    let thermal = f.gate.thermal();
    let hot = f.gate.thermal_at(320.0);
    assert!((hot.temperature_k(DRIVER[0]).unwrap() - 320.0).abs() < 1e-5);
    let context = gate_context();
    let cells = f.gate.cells(&context);
    let r = drive_resistors(&f.gate.base);
    let d = f.gate.base.devices();
    let q = f.gate.charges();
    let cgs = f.gate.caps(false);
    let cgd = f.gate.caps(true);
    let current = f.bindings();
    let driver = f.gate.driver_resistors();
    let bindings = material_bindings(&q, &cgs, &cgd, &current, &driver);
    let network = current_network(&cells, &r, &d, bindings);
    let steady = steady_drive(false);
    let p = authored(0.05, false);
    let state = current_initial(&steady, &p);
    let drive = network
        .bind(&thermal, &steady.native, &p, state, || false)
        .unwrap();
    let h = 1.0 / 1024.0;
    let deck = interval(h);
    let cmd = command(&steady);
    let attempt = |at: &DepositionSimulation<'_>| {
        network
            .step_candidate(
                at,
                h,
                &drive,
                &deck,
                &[cmd],
                1,
                None,
                power_controls(),
                cycle_controls(),
                drive_controls(11),
                &[0.0; 11],
                || false,
            )
            .map(|_| ())
            .map_err(|error| format!("{error:?}"))
    };
    assert_eq!(attempt(&thermal), Ok(()));
    let Err(text) = attempt(&hot) else {
        panic!("a step from another accepted state was taken");
    };
    assert!(
        text.contains(
            "commutated drive must retain every accepted resistance, conduction and capacitance source"
        ),
        "{text}"
    );
}
