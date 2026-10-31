//! G0/G1/G3: original material/interface owners drive capillary transport.
//! All cards below are manufactured numerical fixtures, not measured fluids.

use core::f64::consts::PI;
use std::cell::Cell;

use fs_evidence::ValidityDomain;
use fs_flux::CircularMeniscusEnd;
use fs_flux::capillary::material::{
    MaterialCapillaryError, MaterialCapillaryTube, MaterialMeniscusInput,
};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_material::liquid::surface::{
    EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, LiquidCapillarity, ResolvedLiquidSurfaceState,
    ResolvedLiquidWettingState, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
    resolve_liquid_surface_state, resolve_liquid_wetting_state,
};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{DENSITY_PROPERTY, MaterialPropertySelection};
use fs_qty::{Density, Dims, Frequency, QuantitySpec};

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
        .with("surface_environment", 7.0)
        .unwrap()
}

fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with("T", 280.0, 330.0)
            .with("pressure", 100_000.0, 100_000.0)
            .with("surface_environment", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured capillary constitutive fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn liquid_card(rho: f64, viscosity: f64, basis: LiquidViscosityBasis, sigma: f64) -> MaterialCard {
    let (viscosity_name, viscosity_dims) = match basis {
        LiquidViscosityBasis::Dynamic => (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
        LiquidViscosityBasis::Kinematic => (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS),
    };
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (DENSITY_PROPERTY, Density::DIMS, rho),
        (viscosity_name, viscosity_dims, viscosity),
        (SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS, sigma),
    ] {
        claims.insert_claim(claim(name, dims, value)).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic capillary liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn wetting(
    liquid: &ResolvedLiquidSurfaceState,
    point: &QueryPoint,
    angle: f64,
) -> ResolvedLiquidWettingState {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, Dims::NONE, angle))
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: MaterialStateId {
                chemistry: "synthetic tube wall".into(),
                phase: "solid".into(),
                process: "numerical fixture".into(),
                revision: 0,
            },
            texture_frame: "tube coordinates".into(),
        },
        SurfaceSpec {
            material: liquid.material().material().clone(),
            texture_frame: "liquid boundary".into(),
        },
        SystemContext {
            medium: "declared synthetic liquid/exterior pair".into(),
            third_body: None,
            environment: "surface_environment = 7".into(),
            history: "fixed equilibrium wetting fixture".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_liquid_wetting_state(&card, point, MaterialPropertySelection::SingleClaimOnly).unwrap()
}

fn fixture(
    rho: f64,
    mu: f64,
    sigma: f64,
    angle: f64,
) -> (ResolvedLiquidSurfaceState, ResolvedLiquidWettingState) {
    let card = liquid_card(rho, mu, LiquidViscosityBasis::Dynamic, sigma);
    let point = point(300.0);
    let liquid = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let wetting = wetting(&liquid, &point, angle);
    (liquid, wetting)
}

fn input() -> MaterialMeniscusInput {
    MaterialMeniscusInput {
        radius_m: 1e-4,
        tube_length_m: 0.1,
        wet_length_m: 0.0,
        applied_pressure_pa: 0.0,
        duration_s: 0.2,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-24),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_sourced_washburn_motion_closes_wall_work_and_signed_mass() {
    let (liquid, wetting) = fixture(880.0, 0.004, 0.03, PI / 3.0);
    let original_liquid = liquid.clone();
    let original_wetting = wetting.clone();
    let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
    let input = input();
    let step = tube.step_horizontal(input).unwrap();
    let expected_squared =
        input.radius_m * 0.03 * (PI / 3.0).cos() * input.duration_s / (2.0 * 0.004);
    close(step.meniscus().wet_length_m.powi(2), expected_squared);
    let volume = PI * input.radius_m.powi(2) * expected_squared.sqrt();
    close(step.meniscus().transferred_volume_m3, volume);
    close(step.transferred_mass_kg(), 880.0 * volume);
    close(
        step.meniscus().viscous_dissipation_j,
        -step.meniscus().wall_surface_energy_change_j,
    );
    assert_eq!(step.meniscus().end, CircularMeniscusEnd::TimeLimit);
    assert_eq!(step.discrete().gravitational_energy_change_j, 0.0);
    assert_eq!(step.input(), input);
    assert!(std::ptr::eq(step.tube().liquid(), &liquid));
    assert!(std::ptr::eq(step.tube().wetting(), &wetting));
    assert_eq!(liquid, original_liquid);
    assert_eq!(wetting, original_wetting);
}

#[test]
fn g3_material_swaps_change_physical_motion_without_solver_coefficients() {
    let (liquid, wetting) = fixture(880.0, 0.004, 0.03, PI / 3.0);
    let base = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .step_horizontal(input())
        .unwrap();
    for (mu, sigma, angle, ratio) in [
        (0.008, 0.03, PI / 3.0, 0.5),
        (0.004, 0.06, PI / 3.0, 2.0),
        (0.004, 0.03, PI / 6.0, (PI / 6.0).cos() / (PI / 3.0).cos()),
    ] {
        let (liquid, wetting) = fixture(880.0, mu, sigma, angle);
        let changed = MaterialCapillaryTube::new(&liquid, &wetting)
            .unwrap()
            .step_horizontal(input())
            .unwrap();
        close(
            changed.meniscus().wet_length_m.powi(2),
            ratio * base.meniscus().wet_length_m.powi(2),
        );
    }
}

#[test]
fn g1_same_source_density_sets_jurin_height_gravity_work_and_mass() {
    for rho in [500.0, 1000.0] {
        let (liquid, wetting) = fixture(rho, 0.004, 0.02, PI / 3.0);
        let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
        let input = MaterialMeniscusInput {
            wet_length_m: 0.001,
            duration_s: 1e6,
            ..input()
        };
        let step = tube.step_inclined(input, 9.81).unwrap();
        let height = 2.0 * 0.02 * (PI / 3.0).cos() / (input.radius_m * rho * 9.81);
        close(step.meniscus().wet_length_m, height);
        assert!(step.meniscus().wet_length_m <= height * (1.0 + 2e-15));
        close(
            step.transferred_mass_kg(),
            rho * step.meniscus().transferred_volume_m3,
        );
        close(
            step.meniscus().pressure_work_j - step.meniscus().wall_surface_energy_change_j,
            step.discrete().gravitational_energy_change_j + step.meniscus().viscous_dissipation_j,
        );
        assert_eq!(step.opposing_gravity_m_s2(), 9.81);
    }
}

#[test]
fn g1_kinematic_source_derives_viscosity_with_original_density() {
    let card = liquid_card(800.0, 5e-6, LiquidViscosityBasis::Kinematic, 0.03);
    let point = point(300.0);
    let liquid = resolve_liquid_surface_state(
        &card,
        &point,
        LiquidViscosityBasis::Kinematic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let wetting = wetting(&liquid, &point, PI / 3.0);
    let step = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .step_horizontal(input())
        .unwrap();
    close(liquid.dynamic_viscosity_pa_s(), 0.004);
    close(
        step.meniscus().wet_length_m.powi(2),
        input().radius_m * 0.03 * (PI / 3.0).cos() * input().duration_s / (2.0 * 0.004),
    );
    assert!(
        step.tube()
            .liquid()
            .material()
            .property(DYNAMIC_VISCOSITY_PROPERTY)
            .is_none()
    );
    assert!(
        step.tube()
            .liquid()
            .material()
            .property(KINEMATIC_VISCOSITY_PROPERTY)
            .is_some()
    );
}

#[test]
fn g0_mismatched_or_evolving_sources_refuse_before_motion() {
    let (liquid, wall) = fixture(880.0, 0.004, 0.03, PI / 3.0);
    let mismatched = wetting(&liquid, &point(310.0), PI / 3.0);
    assert!(matches!(
        MaterialCapillaryTube::new(&liquid, &mismatched),
        Err(MaterialCapillaryError::Source(_))
    ));

    let card = liquid_card(880.0, 0.004, LiquidViscosityBasis::Dynamic, 0.03);
    let neglected = resolve_liquid_surface_state(
        &card,
        &point(300.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(matches!(
        MaterialCapillaryTube::new(&neglected, &wall),
        Err(MaterialCapillaryError::Unsupported(_))
    ));

    let rate_point = point(300.0)
        .with_quantity(
            "shear_rate",
            QuantitySpec::dimensional(Frequency::DIMS),
            0.0,
        )
        .unwrap();
    let frozen = resolve_liquid_surface_state(
        &card,
        &rate_point,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let evolving = ResolvedLiquidSurfaceState::try_for_shear_rate_evolution_with_capillarity(
        frozen.material().clone(),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        "T",
        "shear_rate",
    )
    .unwrap();
    let evolving_wall = wetting(&evolving, &rate_point, PI / 3.0);
    assert!(matches!(
        MaterialCapillaryTube::new(&evolving, &evolving_wall),
        Err(MaterialCapillaryError::Unsupported(_))
    ));

    let invalid_input = MaterialMeniscusInput {
        radius_m: 0.0,
        ..input()
    };
    assert!(matches!(
        MaterialCapillaryTube::new(&liquid, &wall)
            .unwrap()
            .step_horizontal(invalid_input),
        Err(MaterialCapillaryError::Flow(_))
    ));
}

#[test]
fn g1_wetting_and_nonwetting_end_events_preserve_time_and_mass_sign() {
    for (angle, initial_length, end, direction) in [
        (PI / 3.0, 0.0, CircularMeniscusEnd::Outlet, 1.0),
        (2.0 * PI / 3.0, 0.1, CircularMeniscusEnd::Inlet, -1.0),
    ] {
        let (liquid, wetting) = fixture(880.0, 0.004, 0.03, angle);
        let input = MaterialMeniscusInput {
            wet_length_m: initial_length,
            duration_s: 1e6,
            ..input()
        };
        let step = MaterialCapillaryTube::new(&liquid, &wetting)
            .unwrap()
            .step_horizontal(input)
            .unwrap();
        assert_eq!(step.meniscus().end, end);
        assert!(step.meniscus().elapsed_time_s < input.duration_s);
        close(
            step.meniscus().elapsed_time_s,
            2.0 * 0.004 * input.tube_length_m.powi(2) / (input.radius_m * 0.03 * angle.cos().abs()),
        );
        close(
            step.transferred_mass_kg(),
            direction * 880.0 * PI * input.radius_m.powi(2) * input.tube_length_m,
        );
        assert!(step.meniscus().viscous_dissipation_j > 0.0);
    }
}

struct DisplacedFixture {
    label: &'static str,
    transport: FluidTransportState,
    reads: Cell<usize>,
    refuse: bool,
}

impl FluidTransport for DisplacedFixture {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        if self.refuse {
            return Err(MaterialError::Parameters {
                what: "outside the fixture's admitted state".into(),
            });
        }
        Ok(self.transport)
    }
}

#[test]
fn g1_two_fluid_binding_retains_original_owner_and_finite_startup_resistance() {
    let (liquid, wetting) = fixture(880.0, 0.004, 0.03, PI / 3.0);
    let displaced = DisplacedFixture {
        label: "original declared exterior fixture",
        transport: FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.004).unwrap(),
        reads: Cell::new(0),
        refuse: false,
    };
    let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
    let input = input();
    let step = tube.step_two_fluid_horizontal(input, &displaced).unwrap();
    let pressure = 2.0 * 0.03 * (PI / 3.0).cos() / input.radius_m;
    // Equal viscosities give constant full-tube resistance and linear motion.
    let expected_length =
        input.radius_m.powi(2) * pressure * input.duration_s / (8.0 * 0.004 * input.tube_length_m);
    close(step.discrete().meniscus.wet_length_m, expected_length);
    close(
        step.discrete().liquid_viscous_dissipation_j
            + step.discrete().displaced_viscous_dissipation_j,
        step.discrete().meniscus.viscous_dissipation_j,
    );
    let volume = PI * input.radius_m.powi(2) * expected_length;
    close(step.transferred_liquid_mass_kg(), 880.0 * volume);
    close(step.displaced_mass_kg(), 2.0 * volume);
    assert!(std::ptr::eq(step.displaced_fluid(), &displaced));
    assert_eq!(
        step.displaced_fluid().label,
        "original declared exterior fixture"
    );
    assert_eq!(step.displaced_transport(), &displaced.transport);
    assert_eq!(step.input(), input);
    assert!(std::ptr::eq(step.tube().liquid(), &liquid));
    assert_eq!(displaced.reads.get(), 1);
}

#[test]
fn g0_displaced_owner_refusal_propagates_without_publishing_a_step() {
    let (liquid, wetting) = fixture(880.0, 0.004, 0.03, PI / 3.0);
    let displaced = DisplacedFixture {
        label: "unavailable exterior state",
        transport: FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.004).unwrap(),
        reads: Cell::new(0),
        refuse: true,
    };
    let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
    assert!(matches!(
        tube.step_two_fluid_horizontal(input(), &displaced),
        Err(MaterialCapillaryError::Transport(_))
    ));
    assert_eq!(displaced.reads.get(), 1);
}

#[test]
fn g1_two_fluid_rise_uses_both_sourced_densities_and_closes_open_column_work() {
    let (liquid, wetting) = fixture(880.0, 0.004, 0.03, PI / 3.0);
    let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
    let input = MaterialMeniscusInput {
        wet_length_m: 0.02,
        applied_pressure_pa: 500.0,
        duration_s: 0.02,
        ..input()
    };
    let gravity = 9.81;
    for (density, direction) in [(400.0, 1.0), (900.0, -1.0)] {
        let displaced = DisplacedFixture {
            label: "original inclined displaced liquid",
            transport: FluidTransportState::try_new(density, 1000.0, 0.5, 0.002).unwrap(),
            reads: Cell::new(0),
            refuse: false,
        };
        let step = tube
            .step_two_fluid_inclined(input, &displaced, gravity)
            .unwrap();
        let inclined = step.inclined_discrete();
        let accepted = &step.discrete().meniscus;
        let volume = accepted.transferred_volume_m3;
        assert_eq!(volume.signum(), direction);
        let midpoint = 0.5 * (input.wet_length_m + accepted.wet_length_m);
        close(
            inclined.gravitational_storage_change_j,
            gravity * (880.0 - density) * midpoint * volume,
        );
        close(
            inclined.outlet_potential_energy_j,
            gravity * density * input.tube_length_m * volume,
        );
        close(
            accepted.pressure_work_j
                - accepted.wall_surface_energy_change_j
                - inclined.gravitational_storage_change_j
                - inclined.outlet_potential_energy_j,
            accepted.viscous_dissipation_j,
        );
        close(
            inclined.hydrostatic_pressure_pa,
            gravity
                * (880.0 * accepted.wet_length_m
                    + density * (input.tube_length_m - accepted.wet_length_m)),
        );
        close(step.transferred_liquid_mass_kg(), 880.0 * volume);
        close(step.displaced_mass_kg(), density * volume);
        assert!(accepted.viscous_dissipation_j > 0.0);
        assert_eq!(step.opposing_gravity_m_s2(), gravity);
        assert!(std::ptr::eq(step.tube().liquid(), &liquid));
        assert!(std::ptr::eq(step.displaced_fluid(), &displaced));
        assert_eq!(displaced.reads.get(), 1);
    }
}

#[test]
fn g3_equal_density_two_fluid_rise_retains_outlet_potential_transport() {
    let density = 880.0;
    let gravity = 9.81;
    let (liquid, wetting) = fixture(density, 0.004, 0.03, PI / 3.0);
    let tube = MaterialCapillaryTube::new(&liquid, &wetting).unwrap();
    let displaced = DisplacedFixture {
        label: "equal density second phase",
        transport: FluidTransportState::try_new(density, 1000.0, 0.5, 0.002).unwrap(),
        reads: Cell::new(0),
        refuse: false,
    };
    let input = MaterialMeniscusInput {
        wet_length_m: 0.02,
        applied_pressure_pa: 1200.0,
        duration_s: 0.02,
        ..input()
    };
    let inclined = tube
        .step_two_fluid_inclined(input, &displaced, gravity)
        .unwrap();
    let horizontal = tube
        .step_two_fluid_horizontal(
            MaterialMeniscusInput {
                applied_pressure_pa: input.applied_pressure_pa
                    - density * gravity * input.tube_length_m,
                ..input
            },
            &displaced,
        )
        .unwrap();
    close(
        inclined.discrete().meniscus.wet_length_m,
        horizontal.discrete().meniscus.wet_length_m,
    );
    close(
        inclined.discrete().meniscus.viscous_dissipation_j,
        horizontal.discrete().meniscus.viscous_dissipation_j,
    );
    assert_eq!(
        inclined.inclined_discrete().gravitational_storage_change_j,
        0.0
    );
    assert!(inclined.inclined_discrete().outlet_potential_energy_j > 0.0);
    let zero = tube
        .step_two_fluid_inclined(input, &displaced, 0.0)
        .unwrap();
    let original = tube.step_two_fluid_horizontal(input, &displaced).unwrap();
    assert_eq!(zero.discrete(), original.discrete());
    assert_eq!(zero.inclined_discrete().gravitational_storage_change_j, 0.0);
    assert_eq!(zero.inclined_discrete().outlet_potential_energy_j, 0.0);
}
