//! Original finite masses exchange heat; native expansion pays its own work.
use super::*;
use fs_blake3::ContentHash;
use fs_lbm::material::free_surface::enthalpy::volume::heat::reservoir::{
    FreeSurfaceSlabReservoirHeatControls, FreeSurfaceSlabReservoirHeatError,
};
use fs_material::phase::exchange::pair::{FinitePhaseBody, FinitePhaseExchangeControls};

const DT: f64 = 0.001;
const HEIGHT: f64 = 0.002;
const CP: f64 = 2000.0;

fn limits(bracket: [f64; 2]) -> FreeSurfaceSlabReservoirHeatControls {
    FreeSurfaceSlabReservoirHeatControls {
        heat: heat_controls(bracket),
        exchange: FinitePhaseExchangeControls {
            absolute_heat_tolerance_j: 1e-12,
            relative_heat_tolerance: 1e-12,
            absolute_exposure_tolerance_j_k: 1e-14,
            relative_exposure_tolerance: 1e-12,
            ..FinitePhaseExchangeControls::default()
        },
        maximum_reservoir_energy_defect_j: 5e-11,
        maximum_total_energy_defect_j: 1e-10,
    }
}

fn reservoir_curve(cp: f64, datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([0x81; 32]), SolidLiquidPhase::Solid,
        [250.0, 550.0].map(|t| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: datum + cp * (t - 300.0),
            temperature_k: t, liquid_mass_fraction: 0.0, bulk_density_kg_m3: 1200.0,
        }).to_vec(),
    ).unwrap()
}

fn body(curve: &EquilibriumEnthalpyPhaseCurve, h: f64, mass: f64) -> FinitePhaseBody<'_> {
    FinitePhaseBody { curve, initial: curve.state_at_specific_enthalpy(h).unwrap(), mass_kg: mass }
}

#[test]
fn g1_finite_reservoir_pays_actual_expansion_work_and_retains_both_source_datums() {
    for hot_density in [980.0, 1020.0] {
        for cooling in [false, true] {
            with_case(Case { hot_density, start_hot: cooling,
                increment_j_kg: if cooling { -200_000.0 } else { 200_000.0 },
                datum: 1e8, ..Case::default() }, |old, _, _, gate, cx| {
                let mass = old.flow().ledger_mass_kg().unwrap();
                let t0 = if cooling { 500.0 } else { 300.0 };
                let tr = if cooling { 300.0 } else { 500.0 };
                let rho0 = if cooling { hot_density } else { 1000.0 };
                let rho1 = 0.5 * (1000.0 + hot_density);
                let u = ((rho0 - rho1) / rho1) * HEIGHT / DT;
                let work = mass * (1.0 / 3.0 + 1.0 / 24.0) * u * u;
                let heat = mass * CP * (400.0 - t0) + work;
                let c = mass * CP;
                let cr = 2.0 * c;
                let reduced = c * cr / (c + cr);
                // Independently invert the finite-two-body exponential so
                // its ONE Q supplies both the chosen thermal change and work.
                let g = -reduced * (-heat / (reduced * (tr - t0))).ln_1p() / DT;
                let chart = reservoir_curve(1000.0, 2e8);
                let original = body(&chart, 2e8 + 1000.0 * (tr - 300.0), 4.0 * mass);
                let bits = old.flow().lattice().grid().state_bits();
                let bracket = if cooling { [-400_000.0, 0.0] } else { [0.0, 400_000.0] };
                let plan = old.prepare_uniform_slab_reservoir_heat(
                    Face3::XMin, original, g, limits(bracket), gate, cx,
                ).unwrap();
                let step = plan.advance(gate, cx).unwrap();
                near(plan.native().endpoint().state().temperature_k(), 400.0, 3e-12);
                near(step.reservoir.initial.temperature_k(), tr - heat / cr, 3e-12);
                near(step.report.heat_to_liquid_j, heat, 3e-12);
                near(step.heat.report.pressure_correction_work_j, work, 3e-8);
                near(step.heat.volume.report.final_height_m, HEIGHT * rho0 / rho1, 3e-12);
                near(step.heat.volume.report.mass_after_kg, mass, 3e-13);
                // Each represented stage has its own 5e-11 J allowance; the
                // small temperature difference must respect those joule gates.
                assert!((c * step.report.mechanical_temperature_correction_k
                    + step.heat.report.pressure_correction_work_j).abs() <= 1e-10);
                assert!(std::ptr::eq(step.reservoir.curve, original.curve));
                assert_eq!(step.reservoir.mass_kg.to_bits(), original.mass_kg.to_bits());
                assert!(step.report.reservoir_energy_defect_j.abs() <= 5e-11);
                assert!(step.report.total_energy_defect_j.abs() <= 1e-10);
                assert!(step.thermal_stage.total_entropy_change_j_k > 0.0);
                assert_eq!(old.flow().lattice().grid().state_bits(), bits);
                assert_eq!(old.accepted_steps(), 0);
            });
        }
    }
}

#[test]
fn g1_reservoir_latent_inventory_and_zero_contact_use_original_states() {
    with_case(Case { increment_j_kg: 50_000.0, ..Case::default() }, |old, _, _, gate, cx| {
        let mass = old.flow().ledger_mass_kg().unwrap();
        let chart = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x82; 32]),
            [(0.0, 300.0, 0.0), (50_000.0, 350.0, 0.0),
                (150_000.0, 350.0, 1.0), (300_000.0, 500.0, 1.0)]
                .map(|(h, t, f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h,
                    temperature_k: t, liquid_mass_fraction: f,
                    bulk_density_kg_m3: 1200.0 - 100.0 * f }).to_vec()).unwrap();
        let original = body(&chart, 100_000.0, 4.0 * mass);
        let plan = old.prepare_uniform_slab_reservoir_heat(Face3::XMin, original,
            mass * CP * 2.0_f64.ln() / DT, limits([0.0, 100_000.0]), gate, cx).unwrap();
        let step = plan.advance(gate, cx).unwrap();
        near(step.thermal_stage.states[0].temperature_k(), 325.0, 1e-13);
        assert_eq!(step.reservoir.initial.temperature_k(), 350.0);
        assert_eq!(step.reservoir.initial.phase(), SolidLiquidPhase::SolidLiquid);
        near(step.reservoir.initial.liquid_mass_fraction(), 0.375, 1e-13);
        near(step.report.heat_to_liquid_j, mass * CP * 25.0, 2e-13);
        assert!(step.report.mechanical_temperature_correction_k < 0.0);
        for temperature in [300.0, 500.0] {
            let sensible = reservoir_curve(1000.0, 0.0);
            let reservoir = body(&sensible, 1000.0 * (temperature - 300.0), mass);
            let no_heat = old.prepare_uniform_slab_reservoir_heat(Face3::XMin, reservoir,
                if temperature == 300.0 { 1.0 } else { 0.0 },
                limits([0.0, 400_000.0]), gate, cx).unwrap();
            let stopped = no_heat.advance(gate, cx).unwrap();
            assert_eq!(stopped.report.heat_to_liquid_j, 0.0);
            assert_eq!(stopped.reservoir.initial, reservoir.initial);
            assert_eq!(stopped.heat.report.net_specific_enthalpy_increment_j_kg, 0.0);
            assert_eq!(stopped.heat.report.pressure_correction_work_j, 0.0);
        }
    });
}

// The SAME source-card/native fixture at authored dt and Cp. The small positive
// Cp used below is manufactured to expose mechanical feedback at modest speeds;
// this is an implementation reference, not an experimental liquid prediction.
fn with_parameters(dt: f64, hot_density: f64, cp: f64,
    run: impl FnOnce(&MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate, &Cx<'_>)) {
    let template = card(0.1, hot_density);
    let mut claims = ClaimSet::new();
    for (_, original) in template.claims().claims_ordered() {
        let mut claim = original.clone();
        if claim.key.name() == "specific_heat_capacity" {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                for (_, value) in knots { *value = cp; }
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(template.id().clone(), claims, vec![]).unwrap();
    let chart = ResolvedHeatCapacityCurve::try_new_exact(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::with_quantity("specific_heat_capacity",
            QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0]))),
        density: &PropertyKey::with_quantity("density",
            QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0]))),
        temperature_axis: "T", lower: &point(300.0), upper: &point(500.0),
        selection: fs_matdb::SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Liquid, 0.0).unwrap();
    let law = LiquidTemperatureLaw::try_new(&card, &point(300.0), "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let source = resolve_liquid_surface_state(&card, &point(300.0), LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected, MaterialPropertySelection::SingleClaimOnly).unwrap();
    let cells: Vec<_> = (0..8).flat_map(|z| (0..8).flat_map(move |y| (0..2).map(move |x|
        MaterialFreeSurfaceInitialCell3 { coordinates: [x,y,z], fill_fraction: 1.0,
            velocity_m_s: [0.0; 3] }))).collect();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan_free_surface_from_material(&source, 0.001, dt, 0.1).unwrap()
            .free_surface3_from_cells(domain(), &cells, cx).unwrap();
        let liquid = flow.with_enthalpy(&law, chart.curve(), FreeSurfaceEnthalpyControls::default(),
            |_| 0.0, &gate).unwrap();
        run(&liquid, &gate, cx);
    });
}

// Retain the actual source-owning plans across these bounded native intervals.
fn partitions(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>, reservoir: FinitePhaseBody<'_>,
    g: f64, cp: f64, count: usize, mut controls: FreeSurfaceSlabReservoirHeatControls,
    gate: &CancelGate, cx: &Cx<'_>) -> (f64, f64, f64, f64) {
    controls.heat.net_specific_enthalpy_bracket_j_kg =
        [0.0, cp * 200.0 - liquid.specific_enthalpy_at([0,0,0]).unwrap()];
    let plan = liquid.prepare_uniform_slab_reservoir_heat(Face3::XMin, reservoir, g,
        controls, gate, cx).unwrap();
    let step = plan.advance(gate, cx).unwrap();
    assert!(std::ptr::eq(step.reservoir.curve, reservoir.curve));
    assert_eq!(step.reservoir.mass_kg, reservoir.mass_kg);
    assert!(step.report.total_energy_defect_j.abs() <= controls.maximum_total_energy_defect_j);
    if count == 1 {
        (plan.native().endpoint().state().temperature_k(), step.reservoir.initial.temperature_k(),
            step.heat.volume.liquid.physical_time_s(), step.heat.report.pressure_correction_work_j)
    } else {
        let (t,r,time,w) = partitions(&step.heat.volume.liquid, step.reservoir, g, cp,
            count-1, controls, gate, cx);
        (t,r,time,w+step.heat.report.pressure_correction_work_j)
    }
}

#[test]
fn g3_finite_capacity_feedback_has_exact_zero_expansion_partition_limit() {
    for count in [1,2,4] {
        let dt = DT / count as f64;
        with_parameters(dt, 1000.0, CP, |old, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            let chart = reservoir_curve(CP, 0.0);
            let reservoir = body(&chart, CP*200.0, 2.0*mass);
            let mut controls = limits([0.0, CP*200.0]);
            controls.heat.volume.volume.projection.duration_s = dt;
            let (t,r,time,w) = partitions(old, reservoir, mass*CP/DT, CP, count, controls, gate, cx);
            let equilibrium = (300.0+2.0*500.0)/3.0;
            let attenuation = (-1.5_f64).exp();
            near(t, equilibrium+(300.0-equilibrium)*attenuation, 3e-12);
            near(r, equilibrium+(500.0-equilibrium)*attenuation, 3e-12);
            near(time, DT, 1e-14);
            assert_eq!(w, 0.0);
            assert!(t < 410.0 && r < 455.0);
        });
    }
}

// Independent continuous equations on the declared native face metric:
// H=M/(A rho(T)), u=H_T*Tdot; C*Tdot=G*(R-T)-B(H)*u*udot,
// Cr*Rdot=-G*(R-T). Integrating k=u²/2 avoids dividing by zero speed.
// This small test-only RK4 reference shares NO production thermal/root kernel.
fn continuous_reference(cp: f64, panels: usize) -> [f64; 3] {
    let mass = 128.0 * 1e-6;
    let capacity = mass * cp;
    let g = capacity / DT;
    let rhs = |[t,r,k]: [f64;3]| {
        let rho = 1000.0 - 0.1*(t-300.0);
        let height = HEIGHT * 1000.0/rho;
        let b = inertia(height, 0.001, rho*1e-9, 64.0);
        let tdot = (2.0*k.max(0.0)).sqrt()/(HEIGHT*1000.0*0.1/(rho*rho));
        let power = g*(r-t);
        [tdot, -power/(2.0*capacity), (power-capacity*tdot)/b]
    };
    let h = DT / panels as f64;
    let mut z = [300.0,500.0,0.0];
    for _ in 0..panels {
        let a = rhs(z);
        let b = rhs(core::array::from_fn(|i| z[i]+0.5*h*a[i]));
        let c = rhs(core::array::from_fn(|i| z[i]+0.5*h*b[i]));
        let d = rhs(core::array::from_fn(|i| z[i]+h*c[i]));
        z = core::array::from_fn(|i| z[i]+h*(a[i]+2.0*b[i]+2.0*c[i]+d[i])/6.0);
        assert!(z[2] >= 0.0);
    }
    z
}

#[test]
fn g3_nonzero_expansion_and_reservoir_feedback_refine_against_continuous_work_balance() {
    let cp = 2e-6;
    let reference = continuous_reference(cp, 32768);
    let coarser_reference = continuous_reference(cp, 16384);
    assert!((reference[0]-coarser_reference[0]).abs() < 1e-5);
    assert!((reference[1]-coarser_reference[1]).abs() < 1e-5);
    let mut errors = Vec::new();
    for count in [2,4,8] {
        let dt = DT/count as f64;
        with_parameters(dt, 980.0, cp, |old, gate, cx| {
            let mass = old.flow().ledger_mass_kg().unwrap();
            let chart = reservoir_curve(cp, 0.0);
            let reservoir = body(&chart, cp*200.0, 2.0*mass);
            let mut controls = limits([0.0, cp*200.0]);
            controls.heat.volume.volume.projection.duration_s = dt;
            controls.heat.maximum_heat_defect_j = 1e-16;
            controls.maximum_reservoir_energy_defect_j = 1e-16;
            controls.maximum_total_energy_defect_j = 2e-16;
            controls.exchange.absolute_heat_tolerance_j = 1e-18;
            controls.exchange.absolute_exposure_tolerance_j_k = 1e-20;
            let (t,r,time,w) = partitions(old, reservoir, mass*cp/DT, cp, count, controls, gate, cx);
            near(time, DT, 1e-14);
            assert!(w > 5e-9 && t > 360.0 && r < 450.0);
            // This account uses BOTH actual endpoints and the measured native
            // mechanical work, not the provisional isolated-pair entropy.
            let defect = mass*cp*(t-300.0)+2.0*mass*cp*(r-500.0)+w;
            assert!(defect.abs() < 2e-15);
            errors.push((t-reference[0]).abs()+(r-reference[1]).abs());
        });
    }
    assert!(errors[1] < errors[0] && errors[2] < errors[1]);
    assert!(errors[2] < 0.55*errors[0]);
}

#[test]
fn g4_source_limits_native_refusal_and_cancellation_preserve_both_participants() {
    with_case(Case { increment_j_kg: 200_000.0, ..Case::default() }, |old, _, _, gate, cx| {
        let mass = old.flow().ledger_mass_kg().unwrap();
        let chart = reservoir_curve(CP, 0.0);
        let reservoir = body(&chart, CP*200.0, 2.0*mass);
        let original_state = reservoir.initial;
        let bits = old.flow().lattice().grid().state_bits();
        let g = mass*CP/DT;
        assert!(matches!(old.prepare_uniform_slab_reservoir_heat(Face3::XMin, reservoir, g,
            limits([0.0,1.0]), gate, cx), Err(FreeSurfaceSlabReservoirHeatError::Native(
                FreeSurfaceSlabHeatError::Balance { .. }))));
        let mut short = limits([0.0,400_000.0]);
        short.exchange.maximum_chart_intervals = 1;
        assert!(old.prepare_uniform_slab_reservoir_heat(Face3::XMin, reservoir, g, short, gate, cx).is_err());
        let unsupported = body(&chart, -50.0*CP, 2.0*mass);
        assert!(old.prepare_uniform_slab_reservoir_heat(Face3::XMin, unsupported, g,
            limits([0.0,400_000.0]), gate, cx).is_err());
        let mut strict = limits([0.0,400_000.0]);
        strict.heat.volume.maximum_kinetic_mapping_change_j = 0.0;
        let rejected = old.prepare_uniform_slab_reservoir_heat(Face3::XMin, reservoir, g,
            strict, gate, cx).unwrap();
        assert!(rejected.advance(gate, cx).is_err());
        let plan = old.prepare_uniform_slab_reservoir_heat(Face3::XMin, reservoir, g,
            limits([0.0,400_000.0]), gate, cx).unwrap();
        gate.request();
        assert!(plan.advance(gate, cx).is_err());
        assert!(old.prepare_uniform_slab_reservoir_heat(Face3::XMin, reservoir, g,
            limits([0.0,400_000.0]), gate, cx).is_err());
        assert_eq!(reservoir.initial, original_state);
        assert_eq!(old.flow().lattice().grid().state_bits(), bits);
        assert_eq!(old.accepted_steps(), 0);
    });
}
