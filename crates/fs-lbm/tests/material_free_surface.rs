//! Data-only liquid swaps through the actual sparse VOF solver. Synthetic
//! liquid cards provide controlled constitutive inputs, not experimental
//! capillary or contact-angle validation.

use fs_evidence::ValidityDomain;
use fs_lbm::d3q19::freesurface3::{ContactModel3, FreeSurface3};
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceScaling, maximum_capillary_timestep_s,
    plan_free_surface_from_material,
};
use fs_lbm::material::plan_scaling_from_material;
use fs_lbm::CS2;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, SURFACE_TENSION_DIMS, resolve_liquid_surface_state,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::Dims;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "actual={actual}, expected={expected}");
}

fn liquid(rho: f64, mu: f64, sigma: Option<f64>) -> ResolvedLiquidSurfaceState {
    let mut rows = vec![
        ("density", Dims([-3, 1, 0, 0, 0, 0]), rho),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), mu),
    ];
    if let Some(value) = sigma {
        rows.push(("surface_tension", SURFACE_TENSION_DIMS, value));
    }
    let mut claims = ClaimSet::new();
    for (name, dims, value) in rows {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 300.0, 300.0)
                .with("surface_environment", 2.0, 2.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "synthetic liquid similitude fixture".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic liquid".into(), phase: "liquid".into(),
        process: "numerical fixture".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let point = QueryPoint::new().with("T", 300.0).unwrap()
        .with("surface_environment", 2.0).unwrap();
    resolve_liquid_surface_state(&card, &point, LiquidViscosityBasis::Dynamic,
        if sigma.is_some() { LiquidCapillarity::Sourced } else { LiquidCapillarity::Neglected },
        MaterialPropertySelection::SingleClaimOnly).unwrap()
}

fn plan(liquid: &ResolvedLiquidSurfaceState) -> MaterialFreeSurfaceScaling<'_> {
    plan_free_surface_from_material(liquid, 0.001, 0.001, 0.1).unwrap()
}
fn domain() -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [16; 3], gravity_m_s2: [0.0, 0.0, -0.1],
        contact: ContactModel3::Neutral, solid_tiles: &[], max_cells: 16 * 16 * 16,
    }
}
fn drop(x: usize, y: usize, z: usize) -> bool {
    let d = [x as f64 - 7.5, y as f64 - 7.5, z as f64 - 7.5];
    d.iter().map(|v| v * v).sum::<f64>() <= 3.5 * 3.5
}

fn with_cx<T>(gate: &fs_exec::CancelGate, f: impl FnOnce(&fs_exec::Cx<'_>) -> T) -> T {
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, Cx, ExecMode, StreamKey};
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(gate, arena, StreamKey { seed: 41, kernel_id: 27, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic))
    })
}

#[test]
fn transferred_fill_and_velocity_preserve_mass_and_momentum_with_material_substitution() {
    use fs_lbm::material::free_surface::MaterialFreeSurfaceInitialCell3;
    let initial = [
        MaterialFreeSurfaceInitialCell3 {
            coordinates: [6, 6, 6], fill_fraction: 1.0, velocity_m_s: [0.03, -0.02, 0.01],
        },
        MaterialFreeSurfaceInitialCell3 {
            coordinates: [7, 6, 6], fill_fraction: 0.25, velocity_m_s: [-0.01, 0.02, 0.0],
        },
        MaterialFreeSurfaceInitialCell3 {
            coordinates: [6, 7, 6], fill_fraction: 0.6, velocity_m_s: [0.01, 0.0, -0.02],
        },
    ];
    let gate = fs_exec::CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let mut bits = None;
        for rho in [1_000.0, 2_000.0] {
            let source = liquid(rho, rho * 0.0001, None);
            let mut flow = plan(&source).free_surface3_from_cells(domain(), &initial, cx).unwrap();
            let mass_unit = rho * 1e-9;
            near(flow.ledger_mass_kg().unwrap(), 1.85 * mass_unit, 2e-14);
            assert_eq!(flow.carry_mass_kg().unwrap(), 0.0);
            let mut momentum = [0.0; 3];
            let mut expected = [0.0; 3];
            for cell in initial {
                let [x, y, z] = cell.coordinates;
                let sample = flow.sample_at_si(x, y, z).unwrap().unwrap();
                near(sample.fill_fraction, cell.fill_fraction, 2e-14);
                for axis in 0..3 {
                    assert!((sample.velocity_m_s[axis] - cell.velocity_m_s[axis]).abs() < 2e-16);
                    momentum[axis] += sample.density_kg_m3 * 1e-9
                        * sample.fill_fraction * sample.velocity_m_s[axis];
                    expected[axis] += mass_unit * cell.fill_fraction * cell.velocity_m_s[axis];
                }
            }
            for axis in 0..3 { near(momentum[axis], expected[axis], 3e-13); }
            let mut reversed = initial;
            reversed.reverse();
            let reordered = plan(&source).free_surface3_from_cells(domain(), &reversed, cx).unwrap();
            assert_eq!(flow.lattice().grid().state_bits(), reordered.lattice().grid().state_bits());
            if let Some(previous) = &bits { assert_eq!(&flow.lattice().grid().state_bits(), previous); }
            bits = Some(flow.lattice().grid().state_bits());
            let before = flow.ledger_mass_kg().unwrap();
            flow.step().unwrap();
            near(flow.ledger_mass_kg().unwrap(), before, 2e-12);
        }
    });
}

#[test]
fn transferred_state_refuses_duplicate_wall_bad_fill_speed_and_cancelled_candidates() {
    use fs_lbm::d3q19::freesurface3::FreeSurfaceError3;
    use fs_lbm::material::MaterialLbmError;
    use fs_lbm::material::free_surface::MaterialFreeSurfaceInitialCell3;
    let source = liquid(1_000.0, 0.1, None);
    let cell = MaterialFreeSurfaceInitialCell3 {
        coordinates: [6, 6, 6], fill_fraction: 0.7, velocity_m_s: [0.03, 0.0, 0.0],
    };
    let gate = fs_exec::CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        assert!(plan(&source).free_surface3_from_cells(domain(), &[cell, cell], cx).is_err());
        let wall = MaterialFreeSurfaceDomain3 { solid_tiles: &[(1, 1, 1)], ..domain() };
        assert!(plan(&source).free_surface3_from_cells(wall, &[cell], cx).is_err());
        for fill_fraction in [0.0, -0.1, 1.1, f64::NAN] {
            let invalid = MaterialFreeSurfaceInitialCell3 { fill_fraction, ..cell };
            assert!(plan(&source).free_surface3_from_cells(domain(), &[invalid], cx).is_err());
        }
        let too_fast = MaterialFreeSurfaceInitialCell3 { velocity_m_s: [0.101, 0.0, 0.0], ..cell };
        assert!(plan(&source).free_surface3_from_cells(domain(), &[too_fast], cx).is_err());
    });
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(plan(&source).free_surface3_from_cells(domain(), &[cell], cx),
            Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))));
    });
}

#[test]
fn free_surface_shares_momentum_scaling_and_preserves_laplace_units() {
    let state = liquid(1_000.0, 0.1, Some(0.01));
    let caloric = FluidTransportState::try_new(1_000.0, 4_000.0, 0.5, 0.1).unwrap();
    for (dx, dt) in [(0.001, 0.001), (0.0005, 0.00025)] {
        let free = plan_free_surface_from_material(&state, dx, dt, 0.1).unwrap();
        let thermal = plan_scaling_from_material(&caloric, dx, dt, 0.1).unwrap();
        assert_eq!(free.flow(), thermal.flow());
        near(free.flow().viscosity * dx * dx / dt, 0.0001, 1e-14);
        near(free.mass_unit_kg(), 1_000.0 * dx * dx * dx, 1e-14);
        // Physical Young-Laplace pressure jump sigma*kappa is independent of
        // dx/dt bookkeeping; this is not a claim about native curvature accuracy.
        for physical_curvature in [-100.0, 100.0] {
            let lattice_curvature = physical_curvature * dx;
            let boundary_density = 1.0 + free.surface_tension_lattice() * lattice_curvature / CS2;
            near(free.pressure_perturbation_si(boundary_density).unwrap(),
                0.01 * physical_curvature, 2e-12);
        }
        assert!(free.capillary_density_excursion_bound() < 1.0);
    }
}

#[test]
fn proportional_rho_mu_sigma_swaps_preserve_trajectory_and_scale_physical_mass_and_pressure() {
    let original = liquid(1_000.0, 0.1, Some(0.005));
    let doubled = liquid(2_000.0, 0.2, Some(0.010));
    let mut a = plan(&original).free_surface3(domain(), drop).unwrap();
    let mut b = plan(&doubled).free_surface3(domain(), drop).unwrap();
    assert_ne!(a.scaling().material().material().card_identity(),
        b.scaling().material().material().card_identity());
    assert_eq!(a.scaling().flow(), b.scaling().flow());
    assert_eq!(a.scaling().surface_tension_lattice(), b.scaling().surface_tension_lattice());
    let initial_a = a.ledger_mass_kg().unwrap();
    let initial_b = b.ledger_mass_kg().unwrap();
    near(initial_b, 2.0 * initial_a, 1e-14);
    for _ in 0..12 {
        a.step().unwrap();
        b.step().unwrap();
    }
    assert_eq!(a.lattice().grid().state_bits(), b.lattice().grid().state_bits());
    near(a.ledger_mass_kg().unwrap(), initial_a, 1e-9);
    near(b.ledger_mass_kg().unwrap(), initial_b, 1e-9);
    let pa = a.sample_at_si(8, 8, 8).unwrap().unwrap();
    let pb = b.sample_at_si(8, 8, 8).unwrap().unwrap();
    assert_eq!(pa.velocity_m_s, pb.velocity_m_s);
    near(pb.density_kg_m3, 2.0 * pa.density_kg_m3, 1e-14);
    assert!(pa.pressure_perturbation_pa.abs() > 1e-6);
    near(pb.pressure_perturbation_pa, 2.0 * pa.pressure_perturbation_pa, 1e-14);
    near(a.region_mass_kg(|x, _, _| x < 8).unwrap()
        + a.region_mass_kg(|x, _, _| x >= 8).unwrap() + a.carry_mass_kg().unwrap(),
        a.ledger_mass_kg().unwrap(), 2e-14);
}

#[test]
fn changing_viscosity_or_surface_data_changes_actual_evolution() {
    let mut populations = Vec::new();
    let mut central_pressure = Vec::new();
    for (mu, sigma) in [(0.1, 0.005), (0.2, 0.005), (0.1, 0.010)] {
        let state = liquid(1_000.0, mu, Some(sigma));
        let mut flow = plan(&state).free_surface3(domain(), drop).unwrap();
        let mass = flow.ledger_mass_kg().unwrap();
        for _ in 0..12 { flow.step().unwrap(); }
        near(flow.ledger_mass_kg().unwrap(), mass, 1e-9);
        populations.push(flow.lattice().grid().state_bits());
        central_pressure.push(flow.sample_at_si(8, 8, 8).unwrap().unwrap().pressure_perturbation_pa);
    }
    assert_ne!(populations[0], populations[1], "source viscosity must reach collision");
    assert_ne!(populations[0], populations[2], "source sigma must reach interface reconstruction");
    assert!((central_pressure[0] - central_pressure[1]).abs() > 1e-6);
    assert!((central_pressure[0] - central_pressure[2]).abs() > 1e-6);
}

#[test]
fn physical_samples_use_guo_velocity_and_hide_gas_populations() {
    let source = liquid(1_000.0, 0.1, None);
    let free = plan(&source).free_surface3(domain(), drop).unwrap();
    assert_eq!(free.scaling().surface_tension_lattice(), 0.0);
    let sample = free.sample_at_si(8, 8, 8).unwrap().unwrap();
    near(sample.density_kg_m3, 1_000.0, 1e-14);
    near(sample.velocity_m_s[2], 0.5 * domain().gravity_m_s2[2] * 0.001, 1e-14);
    assert!(sample.pressure_perturbation_pa.abs() < 1e-10);
    assert!(free.sample_at_si(0, 0, 0).unwrap().is_none());
    assert!(free.sample_at_si(16, 8, 8).unwrap().is_none());
    assert!(free.lattice().sample_at(usize::MAX, 0, 0).is_none());
}

#[test]
fn capillary_reconstruction_and_allocation_domains_refuse_before_construction() {
    let source = liquid(1_000.0, 0.1, Some(0.01));
    // At this long timestep native sigma*kappa/CS2 can exceed one and
    // reconstruct nonpositive density. A zero speed does not make it valid.
    assert!(plan_free_surface_from_material(&source, 0.001, 0.1, 0.0).is_err());
    for (dx, dt, speed) in [(0.0, 0.001, 0.0), (0.001, -1.0, 0.0),
        (0.001, 0.001, 1.0), (f64::NAN, 0.001, 0.0), (1e100, 1e-100, 0.0)] {
        assert!(plan_free_surface_from_material(&source, dx, dt, speed).is_err());
    }
    let valid = plan(&source);
    for input in [
        MaterialFreeSurfaceDomain3 { cells: [15, 16, 16], ..domain() },
        MaterialFreeSurfaceDomain3 { max_cells: 4095, ..domain() },
        MaterialFreeSurfaceDomain3 { cells: [usize::MAX, 16, 16], ..domain() },
        MaterialFreeSurfaceDomain3 { gravity_m_s2: [f64::NAN, 0.0, 0.0], ..domain() },
        MaterialFreeSurfaceDomain3 { solid_tiles: &[(4, 0, 0)], ..domain() },
        MaterialFreeSurfaceDomain3 { solid_tiles: &[(1, 1, 1)], ..domain() },
    ] {
        assert!(valid.free_surface3(input, drop).is_err());
    }
    assert!(valid.free_surface3(domain(), |_, _, _| false).is_err());

    // Native clamp is the same owner used by the scaling admission.
    near(valid.capillary_density_excursion_bound(),
        valid.surface_tension_lattice() * FreeSurface3::CURVATURE_LIMIT / CS2, 1e-14);
}

#[test]
fn g1_capillary_timestep_resolves_the_shortest_wave_and_scales_with_material_and_mesh() {
    let source = liquid(1_000.0, 0.1, Some(0.01));
    let dx = 0.001;
    let limit = maximum_capillary_timestep_s(&source, dx).unwrap().unwrap();
    // Independent inviscid dispersion omega^2 = sigma*k^3/rho. Two waves
    // travelling in opposite directions cannot cover more than one cell in dt.
    let k = core::f64::consts::PI / dx;
    let omega = (0.01 * k.powi(3) / 1_000.0).sqrt();
    near(2.0 * omega / k * limit, dx, 3e-15);
    for (scale, expected) in [(0.25_f64, 0.125), (4.0, 8.0)] {
        let refined = maximum_capillary_timestep_s(&source, dx * scale).unwrap().unwrap();
        near(refined / limit, expected, 3e-15);
    }
    for (rho, sigma, expected) in [(4_000.0, 0.01, 2.0), (1_000.0, 0.04, 0.5),
        (2_000.0, 0.02, 1.0)] {
        let changed = liquid(rho, 0.1, Some(sigma));
        let changed_limit = maximum_capillary_timestep_s(&changed, dx).unwrap().unwrap();
        near(changed_limit / limit, expected, 3e-15);
    }
    let at_limit = plan_free_surface_from_material(&source, dx, limit, 0.0).unwrap();
    near(at_limit.capillary_timestep_ratio(), 1.0, 3e-15);
    let shorter = plan_free_surface_from_material(&source, dx, limit * 0.25, 0.0).unwrap();
    near(shorter.capillary_timestep_ratio(), 0.25, 3e-15);
}

#[test]
fn g0_zero_speed_and_positive_pressure_do_not_admit_underresolved_capillary_waves() {
    use fs_lbm::material::MaterialLbmError;
    let source = liquid(1_000.0, 0.1, Some(0.01));
    let dx = 0.001;
    let limit = maximum_capillary_timestep_s(&source, dx).unwrap().unwrap();
    let dt = 1.1 * limit;
    let sigma_lattice = 0.01 * dt.powi(2) / (1_000.0 * dx.powi(3));
    assert!(sigma_lattice * FreeSurface3::CURVATURE_LIMIT / CS2 < 1.0);
    assert!(matches!(
        plan_free_surface_from_material(&source, dx, dt, 0.0),
        Err(MaterialLbmError::Invalid(reason)) if reason.contains("capillary-wave")
    ));
    // Changing only source sigma tightens the bound for the same real timeline.
    let stronger = liquid(1_000.0, 0.1, Some(0.04));
    assert!(plan_free_surface_from_material(&source, dx, 0.75 * limit, 0.0).is_ok());
    assert!(plan_free_surface_from_material(&stronger, dx, 0.75 * limit, 0.0).is_err());
    let omitted = liquid(1_000.0, 0.1, None);
    assert_eq!(maximum_capillary_timestep_s(&omitted, dx).unwrap(), None);
    let no_capillarity = plan_free_surface_from_material(&omitted, dx, dt, 0.0).unwrap();
    assert_eq!(no_capillarity.capillary_timestep_ratio(), 0.0);
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(maximum_capillary_timestep_s(&source, invalid).is_err());
        assert!(maximum_capillary_timestep_s(&omitted, invalid).is_err());
    }
    // Native transferred-state ingress shares the same gate, independently of
    // the material constructor. The coefficient refuses before missing cells.
    let gate = fs_exec::CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        use fs_lbm::d3q19::freesurface3::FreeSurfaceError3;
        assert!(matches!(FreeSurface3::from_cells([8; 3], 0.8, [0.0; 3],
            sigma_lattice, ContactModel3::Neutral, &[], &[], cx),
            Err(FreeSurfaceError3::InvalidInitialState { reason, .. })
                if reason.contains("capillary-wave")));
    });
}
