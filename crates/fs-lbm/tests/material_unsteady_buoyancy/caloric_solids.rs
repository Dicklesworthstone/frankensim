//! Original source Cp on native transient storage and donor enthalpy fluxes.
//! Synthetic scalar curves are numerical oracles, not measured fluid claims.

#[path = "compact_caloric.rs"]
mod compact;

use super::solid::{card, law, specific_energy};
use super::*;
use fs_lbm::conjugate::{ConvectionScheme, FvSgsConfig};
use fs_lbm::material::conjugate::buoyant::SolidThermalConfig;
use fs_lbm::material::conjugate::buoyant::enthalpy::{
    MaterialEnthalpyConfig, MaterialEnthalpyTransientConfig,
};
use fs_lbm::material::conjugate::solid::SolidThermalLaw;
use fs_lbm::sgs3::{SgsModel, SgsUse};
use fs_matdb::{ClaimSelection, MaterialCard};
use fs_material::thermal::transport::isobaric::{
    IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery,
};

fn source_caloric(
    card: &MaterialCard,
    source: &SolidThermalLaw<'_>,
    span: [f64; 2],
) -> IsobaricHeatCapacityLaw {
    let reference = source.conductivity().reference();
    let point = reference.query();
    let lower = point.clone().with("T", span[0]).unwrap();
    let upper = point.with("T", span[1]).unwrap();
    let selected = reference.property("cp").unwrap().answer().receipt.selected;
    let claim = card.claims().claim(selected).unwrap();
    IsobaricHeatCapacityLaw::try_new(
        card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &claim.key,
            temperature_axis: "T",
            pressure_axis: "pressure",
            lower: &lower,
            upper: &upper,
            selection: ClaimSelection::Pinned(selected),
        },
        300.0,
    )
    .unwrap()
}

fn transport(
    source: &SolidThermalLaw<'_>,
    temperature: f64,
) -> Result<FluidTransportState, MaterialError> {
    let sample = source
        .conductivity()
        .at_temperature(0, 0, temperature)
        .map_err(|error| MaterialError::State {
            what: error.to_string(),
        })?;
    let value = |name| sample.material.property(name).unwrap().value_si();
    FluidTransportState::try_new(value("rho"), value("cp"), value("k"), 1.0)
}

fn bounds() -> BoussinesqConfig<3> {
    let mut bounds = limits([0.0; 3]);
    bounds.max_relative_specific_heat_departure = 0.05;
    bounds
}

fn energy() -> EnergyConfig {
    EnergyConfig {
        scheme: ConvectionScheme::Upwind,
        tolerance: 1e-13,
        ..EnergyConfig::default()
    }
}

fn controls() -> MaterialEnthalpyTransientConfig {
    MaterialEnthalpyTransientConfig {
        caloric: MaterialEnthalpyConfig {
            reference_temperature_k: 300.0,
            maximum_flux_defect_w: 1e-9,
        },
        maximum_passes_per_step: 100,
        maximum_energy_solves: 1000,
        temperature_tolerance_k: 1e-9,
        maximum_storage_defect_w: 1e-8,
        maximum_energy_closure_j: 1e-8,
        temperature_relaxation: 1.0,
    }
}

fn sgs(cs: f64) -> FvSgsConfig {
    FvSgsConfig {
        model: SgsModel::Smagorinsky { cs },
        filter_width_factor: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn shear_cavity(with_solid: bool) -> (VoxelDomain, SimpleConfig, ThermalSetup, Vec<f64>, UnsteadyConfig) {
    let domain = VoxelDomain::from_fn(4, 5, 1, 0.25, |p| {
        if with_solid && p[1] < 0.25 { Voxel::Solid(0) } else { Voxel::Fluid }
    }).unwrap();
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall { velocity: [0.05, 0.0, 0.0] };
    faces[4] = FvBoundary::Symmetry;
    faces[5] = FvBoundary::Symmetry;
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.power_w = vec![0.01; domain.cell_count()];
    let initial = (0..domain.cell_count()).map(|c| {
        if domain.is_fluid(c) { 304.0 + domain.coords(c)[0] as f64 } else { 300.0 }
    }).collect();
    let mut unsteady = UnsteadyConfig::new(0.25, 4);
    unsteady.scheme = TimeScheme::BackwardEuler;
    unsteady.inner_iterations = 1000;
    unsteady.inner_tolerance = 1e-9;
    (domain, SimpleConfig::new(faces), thermal, initial, unsteady)
}

fn shear_transport(source: &SolidThermalLaw<'_>, t: f64) -> Result<FluidTransportState, MaterialError> {
    // Original Cp/k/rho queries plus an explicit synthetic molecular mu(T),
    // not a turbulent effective viscosity or an empirical source claim.
    let state = transport(source, t)?;
    FluidTransportState::try_new(state.density_kg_m3(), state.specific_heat_j_kg_k(),
        state.thermal_conductivity_w_m_k(), 0.03 / (1.0 + 0.15 * (t - 300.0)))
}

#[test]
fn sgs_joint_source_calorimetry_preserves_material_refresh_and_total_heat() {
    let fluid_card = card(0.004, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let reference = shear_transport(&fluid_source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.01);
    let solid = law(&solid_card, [290.0, 400.0]);
    let solids = std::slice::from_ref(&solid);
    let (domain, flow, thermal, initial, unsteady) = shear_cavity(true);
    let original = initial.clone();
    let gate = CancelGate::new();
    let run = |cs| bound.march_natural_convection_with_enthalpy_and_material_solids_and_sgs(
        &domain, solids, &flow, &thermal, &energy(), &initial, bounds(),
        &caloric, controls(), SolidThermalConfig::default(), &unsteady, sgs(cs),
        |_| 1.0, |_| 1.0, None, &gate, |_, t, _| shear_transport(&fluid_source, t),
    ).unwrap();
    let active = run(0.5);
    let zero = run(0.0);
    let laminar = bound.march_natural_convection_with_enthalpy_and_material_solids(
        &domain, solids, &flow, &thermal, &energy(), &initial, bounds(),
        &caloric, controls(), SolidThermalConfig::default(), &unsteady,
        |_| 1.0, |_| 1.0, None, &gate, |_, t, _| shear_transport(&fluid_source, t),
    ).unwrap();
    assert_eq!(zero.fluid.convection.march.flow, laminar.fluid.convection.march.flow);
    assert_eq!(zero.fluid.convection.march.temperature, laminar.fluid.convection.march.temperature);
    assert_eq!(zero.solid_energy_j, laminar.solid_energy_j);
    let march = &active.fluid.convection.march;
    assert!(march.flow.flow.eddy_viscosity_m2_s.iter().any(|&nu| nu > 1e-7));
    let velocity_difference = march.flow.flow.velocity_m_s.iter()
        .zip(&zero.fluid.convection.march.flow.flow.velocity_m_s)
        .flat_map(|(a, b)| a.iter().zip(b)).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
    assert!(velocity_difference > 1e-6, "SGS must reach actual material momentum");
    assert!(march.temperature.iter().zip(&zero.fluid.convection.march.temperature)
        .any(|(a, b)| (a - b).abs() > 1e-6), "SGS must reach the joint heat solve");

    // Independent original-law integrals, with fixed reference masses. Wall
    // motion supplies no invented heat in this explicitly non-heated entry.
    let mut solid_change = 0.0;
    let mut total_change = 0.0;
    for (cell, (&t, &old)) in march.temperature.iter().zip(&initial).enumerate() {
        let slope = if domain.is_fluid(cell) { 0.004 } else { 0.2 };
        let change = 3.0 * domain.volume(cell)
            * (specific_energy(t, slope) - specific_energy(old, slope));
        total_change += change;
        if !domain.is_fluid(cell) {
            solid_change += change;
            close(active.solid_energy_j[cell],
                3.0 * domain.volume(cell) * (specific_energy(t, 0.2) + 20.0), 1e-9);
        }
    }
    let supplied = unsteady.time_step_s * unsteady.steps as f64 * thermal.power_w.iter().sum::<f64>();
    close(total_change, supplied, 1e-8);
    close(active.solid_records.iter().map(|r| r.solid_stored_energy_change_j).sum(), solid_change, 1e-8);
    assert_eq!(active.fluid.records.len(), unsteady.steps);
    for (fluid_record, solid_record) in active.fluid.records.iter().zip(&active.solid_records) {
        assert_eq!(fluid_record.time_s, solid_record.time_s);
        assert_eq!(fluid_record.passes, solid_record.iterations);
        assert!(fluid_record.maximum_storage_defect_w <= controls().maximum_storage_defect_w);
    }
    for heat in &march.energy_records { assert!(heat.closure_j.abs() <= 1e-8); }
    for sample in active.fluid.convection.applied.iter().chain(&active.fluid.convection.endpoint) {
        close(sample.material.resolved_transport().dynamic_viscosity_pa_s(),
            0.03 / (1.0 + 0.15 * (sample.temperature_k - 300.0)), 1e-13);
    }

    // A refusal at the first evolved source endpoint returns no partially
    // accepted caloric/SGS result and leaves original physical inputs intact.
    let mut queries = 0usize;
    let refused = bound.march_natural_convection_with_enthalpy_and_material_solids_and_sgs(
        &domain, solids, &flow, &thermal, &energy(), &initial, bounds(),
        &caloric, controls(), SolidThermalConfig::default(), &unsteady, sgs(0.5),
        |_| 1.0, |_| 1.0, None, &gate, |_, t, _| {
            queries += 1;
            if queries > domain.fluid_count() {
                Err(MaterialError::State { what: "synthetic endpoint source refusal".into() })
            } else { shear_transport(&fluid_source, t) }
        },
    );
    assert!(matches!(refused, Err(MaterialBuoyancyError::Material { .. })));
    assert_eq!(queries, domain.fluid_count() + 1);
    assert_eq!(initial, original);
    assert_eq!(bound.transport(), &reference);
}

#[test]
fn sgs_source_caloric_and_solid_entries_recover_constant_fluid_heat_capacity() {
    let fluid_card = card(0.0, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let reference = shear_transport(&fluid_source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.01);
    let solid = law(&solid_card, [290.0, 400.0]);
    let gate = CancelGate::new();
    for with_solid in [false, true] {
        let (domain, flow, thermal, initial, unsteady) = shear_cavity(with_solid);
        let (caloric_march, fixed_cp_march) = if with_solid {
            let joint = bound.march_natural_convection_with_enthalpy_and_material_solids_and_sgs(
                &domain, std::slice::from_ref(&solid), &flow, &thermal, &energy(), &initial,
                bounds(), &caloric, controls(), SolidThermalConfig::default(), &unsteady, sgs(0.5),
                |_| 1.0, |_| 1.0, None, &gate, |_, t, _| shear_transport(&fluid_source, t),
            ).unwrap();
            let solid_only = bound.march_natural_convection_with_material_solids_and_sgs(
                &domain, std::slice::from_ref(&solid), &flow, &thermal, &energy(), &initial,
                bounds(), &unsteady, SolidThermalConfig::default(), sgs(0.5),
                |_| 1.0, |_| 1.0, None, &gate, |_, t, _| shear_transport(&fluid_source, t),
            ).unwrap();
            (joint.fluid.convection.march, solid_only.fluid.march)
        } else {
            let caloric_run = bound.march_natural_convection_with_enthalpy_and_sgs(
                &domain, &[], &flow, &thermal, &energy(), &initial, bounds(),
                &caloric, controls(), &unsteady, sgs(0.5), |_| 1.0, |_| 1.0,
                None, &gate, |_, t, _| shear_transport(&fluid_source, t),
            ).unwrap();
            let original = bound.march_natural_convection_with_sgs(
                &domain, &[], &flow, &thermal, &energy(), &initial, bounds(),
                &unsteady, sgs(0.5), |_| 1.0, |_| 1.0, None, &gate,
                |_, t, _| shear_transport(&fluid_source, t),
            ).unwrap();
            (caloric_run.convection.march, original.march)
        };
        assert!(caloric_march.flow.flow.eddy_viscosity_m2_s.iter().any(|&nu| nu > 1e-7));
        for (&t, &expected) in caloric_march.temperature.iter().zip(&fixed_cp_march.temperature) {
            assert!((t - expected).abs() < 1e-8, "{t} != {expected}");
        }
        for (u, expected) in caloric_march.flow.flow.velocity_m_s.iter()
            .zip(&fixed_cp_march.flow.flow.velocity_m_s) {
            for (&a, &b) in u.iter().zip(expected) { close(a, b, 1e-9); }
        }
    }
}

#[test]
fn g1_fluid_and_evolving_solid_share_one_exact_heat_balance_and_clock() {
    let fluid_card = card(0.002, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let reference = transport(&fluid_source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.01);
    let solid = law(&solid_card, [290.0, 400.0]);
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    let run = bound
        .march_natural_convection_with_enthalpy_and_material_solids(
            &domain,
            &[solid],
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &ThermalSetup::new([ThermalFace::Adiabatic; 6]),
            &energy(),
            &[300.0, 320.0],
            bounds(),
            &caloric,
            controls(),
            SolidThermalConfig::default(),
            &UnsteadyConfig::new(1.0, 3),
            |_| 1.0,
            |_| 1.0,
            None,
            &CancelGate::new(),
            |_, t, _| transport(&fluid_source, t),
        )
        .unwrap();
    let (mut ts, mut tf) = (300.0, 320.0);
    for (index, heat) in run.fluid.convection.march.energy_records.iter().enumerate() {
        let (mut lo, mut hi) = (ts, tf);
        for _ in 0..100 {
            let candidate = (lo + hi) * 0.5;
            let solid_delta = specific_energy(candidate, 0.2) - specific_energy(ts, 0.2);
            let fluid_h = specific_energy(tf, 0.002) - solid_delta;
            let (mut flo, mut fhi) = (300.0, tf);
            for _ in 0..100 {
                let mid = (flo + fhi) * 0.5;
                if specific_energy(mid, 0.002) > fluid_h {
                    fhi = mid;
                } else {
                    flo = mid;
                }
            }
            let next_f = (flo + fhi) * 0.5;
            let k = 1.0 + 0.01 * (candidate - 300.0);
            let conductance = 1.0 / (0.5 / k + 0.5);
            if 3.0 * solid_delta > conductance * (next_f - candidate) {
                hi = candidate;
            } else {
                lo = candidate;
            }
        }
        let next_s = (lo + hi) * 0.5;
        let solid_delta = 3.0 * (specific_energy(next_s, 0.2) - specific_energy(ts, 0.2));
        let fluid_h = specific_energy(tf, 0.002) - solid_delta / 3.0;
        let (mut lo, mut hi) = (300.0, tf);
        for _ in 0..100 {
            let mid = (lo + hi) * 0.5;
            if specific_energy(mid, 0.002) > fluid_h {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        ts = next_s;
        tf = (lo + hi) * 0.5;
        close(
            run.solid_records[index].solid_stored_energy_change_j,
            solid_delta,
            2e-9,
        );
        assert!(run.fluid.records[index].maximum_storage_defect_w <= 1e-8);
        assert_eq!(
            run.fluid.records[index].passes,
            run.solid_records[index].iterations
        );
        close(heat.stored_energy_change_j, 0.0, 2e-8);
        close(heat.closure_j, 0.0, 2e-8);
        assert_eq!(run.fluid.records[index].time_s, (index + 1) as f64);
        assert_eq!(
            run.solid_records[index].time_s,
            run.fluid.records[index].time_s
        );
    }
    close(run.fluid.convection.march.temperature[0], ts, 1e-10);
    close(run.fluid.convection.march.temperature[1], tf, 1e-10);
    close(
        run.solid_energy_j[0],
        3.0 * (20.0 + specific_energy(ts, 0.2)),
        1e-10,
    );
    assert_eq!(
        run.solid_endpoint[0].material.card_identity(),
        solid_card.content_hash()
    );
    assert!(core::ptr::eq(run.fluid.caloric, &caloric));
    let work: usize = run
        .fluid
        .records
        .iter()
        .map(|record| record.energy_solves)
        .sum();
    assert!(work > run.fluid.records.len() && work <= controls().maximum_energy_solves);
}

#[test]
fn g1_combined_caloric_fixed_reservoir_and_radiation_close_exact_storage_once() {
    let fluid_card = card(0.0001, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let reference = transport(&fluid_source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.01);
    let solid = law(&solid_card, [290.0, 400.0]);
    let dx = 0.01;
    let domain = VoxelDomain::from_fn(1, 1, 2, dx, |p| {
        if p[2] < dx {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.fixed_temperature = vec![(0, 330.0), (0, 350.0)];
    let mut surroundings = [Some(300.0); 6];
    surroundings[4] = None;
    let radiation = RadiationConfig::new(vec![0.9], surroundings);
    let run = bound
        .march_natural_convection_with_enthalpy_and_material_solids(
            &domain,
            &[solid],
            &SimpleConfig::new([FvBoundary::wall(); 6]),
            &thermal,
            &energy(),
            &[300.0; 2],
            bounds(),
            &caloric,
            controls(),
            SolidThermalConfig::default(),
            &UnsteadyConfig::new(1.0, 1),
            |_| 1.0,
            |_| 1.0,
            Some(&radiation),
            &CancelGate::new(),
            |_, t, _| transport(&fluid_source, t),
        )
        .unwrap();
    let solid_change = 3.0 * domain.volume(0) * specific_energy(350.0, 0.2);
    let fluid_change =
        3.0 * domain.volume(1) * specific_energy(run.fluid.convection.march.temperature[1], 0.0001);
    let heat = &run.fluid.convection.march.energy_records[0];
    close(
        run.solid_records[0].solid_stored_energy_change_j,
        solid_change,
        1e-13,
    );
    close(
        heat.stored_energy_change_j,
        solid_change + fluid_change,
        1e-12,
    );
    close(heat.boundary_outflow_j, -heat.stored_energy_change_j, 1e-8);
    close(heat.closure_j, 0.0, 1e-8);
    close(
        heat.radiated_w,
        0.9 * STEFAN_BOLTZMANN * dx * dx * (350.0_f64.powi(4) - 300.0_f64.powi(4)),
        1e-10,
    );
    assert!(run.fluid.records[0].energy_solves > run.fluid.records[0].passes);
}

#[test]
fn g0_combined_source_span_whole_march_budget_and_cancellation_preserve_inputs() {
    let fluid_card = card(0.002, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let reference = transport(&fluid_source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let solid_card = card(0.2, 0.01);
    let narrow = law(&solid_card, [300.0, 301.0]);
    let broad = law(&solid_card, [290.0, 400.0]);
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    let thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    let flow = SimpleConfig::new([FvBoundary::wall(); 6]);
    let unsteady = UnsteadyConfig::new(1.0, 1);
    let run = |source: &SolidThermalLaw<'_>, config, gate: &CancelGate| {
        bound.march_natural_convection_with_enthalpy_and_material_solids(
            &domain,
            std::slice::from_ref(source),
            &flow,
            &thermal,
            &energy(),
            &[300.0, 320.0],
            bounds(),
            &caloric,
            config,
            SolidThermalConfig::default(),
            &unsteady,
            |_| 1.0,
            |_| 1.0,
            None,
            gate,
            |_, t, _| transport(&fluid_source, t),
        )
    };
    assert!(matches!(
        run(&narrow, controls(), &CancelGate::new()),
        Err(MaterialBuoyancyError::SolidThermal { cell: 0, .. })
    ));
    let mut budget = controls();
    budget.maximum_energy_solves = 1;
    assert!(run(&broad, budget, &CancelGate::new()).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        run(&broad, controls(), &gate),
        Err(MaterialBuoyancyError::Solver(ChtError::Cancelled))
    ));
    let retry = run(&broad, controls(), &CancelGate::new()).unwrap();
    assert!(retry.fluid.convection.march.temperature[0] > 300.0);
    assert_eq!(broad.storage().temperature_span_k(), [290.0, 400.0]);
    assert_eq!(broad.conductivity().reference().query().axes()["T"], 300.0);
}

#[test]
fn g0_unrelated_large_reference_mass_cannot_mask_small_cell_caloric_defect() {
    use fs_matdb::{ClaimSet, MaterialStateId, PropertyValue};
    use fs_qty::Density;
    let fluid_card = card(0.0, 0.0);
    let fluid_source = law(&fluid_card, [290.0, 400.0]);
    let caloric = source_caloric(&fluid_card, &fluid_source, [290.0, 400.0]);
    let reference = transport(&fluid_source, 300.0).unwrap();
    let bound = MaterialChtFluid::try_new(&reference).unwrap();
    let small_card = card(0.2, 0.0);
    let small = law(&small_card, [290.0, 400.0]);
    // A distinct, unheated source specimen has huge fixed reference mass.
    // It has zero actual storage change and cannot supply the small cell's
    // omitted nonlinear Cp energy. No raw storage coefficient is injected.
    let mut claims = ClaimSet::new();
    for (_, source) in fluid_card.claims().claims_ordered() {
        let mut claim = source.clone();
        if let PropertyValue::Curve { dims, knots, .. } = &mut claim.value {
            if *dims == Density::DIMS {
                for (_, value) in knots {
                    *value *= 1e14;
                }
            }
        }
        claims.insert_claim(claim).unwrap();
    }
    let large_card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic reference-mass scale".into(),
            phase: "solid".into(),
            process: "local-energy regression".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let large = law(&large_card, [290.0, 400.0]);
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    domain.set(1, 0, 0, Voxel::Solid(1));
    let mut thermal = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    thermal.fixed_temperature = vec![(0, 310.0), (1, 300.0), (2, 300.0)];
    let mut fluid_controls = controls();
    fluid_controls.maximum_passes_per_step = 1;
    fluid_controls.temperature_tolerance_k = 20.0;
    fluid_controls.maximum_energy_closure_j = 1e6;
    let mut solid_controls = SolidThermalConfig::default();
    solid_controls.maximum_iterations = 1;
    solid_controls.temperature_tolerance_k = 20.0;
    // The small cell's exact 3*integral_300^310(2+0.2*(T-300))dT is 90 J;
    // its first tangent stores only 60 J. The large cell's roughly 5 kJ
    // roundoff scale must NOT waive this local 30 J constitutive defect.
    let result = bound.march_natural_convection_with_enthalpy_and_material_solids(
        &domain,
        &[small, large],
        &SimpleConfig::new([FvBoundary::wall(); 6]),
        &thermal,
        &energy(),
        &[300.0; 3],
        bounds(),
        &caloric,
        fluid_controls,
        solid_controls,
        &UnsteadyConfig::new(1.0, 1),
        |_| 1.0,
        |_| 1.0,
        None,
        &CancelGate::new(),
        |_, t, _| transport(&fluid_source, t),
    );
    assert!(
        matches!(result,Err(MaterialBuoyancyError::SolidNotConverged(report))
        if report.storage_defect_j>29.0 && report.storage_defect_j<31.0)
    );
}
