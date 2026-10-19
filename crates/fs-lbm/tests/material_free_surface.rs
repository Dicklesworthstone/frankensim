//! Data-only liquid swaps through the actual sparse VOF solver. Synthetic
//! liquid cards provide controlled constitutive inputs, not experimental
//! capillary or contact-angle validation.

use fs_evidence::ValidityDomain;
use fs_lbm::d3q19::freesurface3::{ContactModel3, FreeSurface3};
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceScaling, plan_free_surface_from_material,
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
