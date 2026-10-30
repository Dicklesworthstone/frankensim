//! Run `cargo run --release -p fs-flux --no-default-features --features experimental-compressible --example shock_box -- 32`.
//! CSV is written to stdout; accepted-step balances and support loads go to stderr.
//! This is an inviscid frozen-N2 engineering example, not experimental validation.

use fs_evidence::ValidityDomain;
use fs_flux::compressible::{
    AxisBoundaries, Boundaries2d, BoundaryCondition, BoundarySide, CartesianGrid, Euler2d,
    FrozenGas, State, StepControl,
};
use fs_gas::dynamics::CaloricDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure, Temperature, Time};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let n = args
        .next()
        .map(|s| s.parse::<usize>())
        .transpose()?
        .unwrap_or(32);
    if !(8..=128).contains(&n) || args.next().is_some() {
        return Err("usage: shock_box [cells_per_axis in 8..=128]".into());
    }
    let gas = nitrogen()?;
    let grid = CartesianGrid::new(n, n, 1.0 / n as f64, 1.0 / n as f64, 0.01)?;
    let ambient_pressure = gas.gas_constant() * 400.0; // rho=1 kg/m³, T=400 K.
    let ambient = gas.primitive(gas.conserved(1.0, 0.0, ambient_pressure)?, || false)?;
    // More than one box acoustic crossing allows the pulse to reach and reflect
    // from the walls. An off-centre pulse exposes a nonzero net support impulse.
    let duration = 1.05 / ambient.sound_speed;
    let states = (0..n * n)
        .map(|index| {
            let [x, y] = grid.cell_center(index).expect("valid cell index");
            let inside = (x - 0.35).powi(2) + (y - 0.5).powi(2) < 0.12_f64.powi(2);
            let pressure = ambient_pressure * if inside { 2.0 } else { 1.0 };
            State::from_primitive(&gas, 1.0, [0.0, 0.0], pressure)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let walls = AxisBoundaries::Ends {
        lower: BoundaryCondition::Reflecting,
        upper: BoundaryCondition::Reflecting,
    };
    let mut flow = Euler2d::new(
        gas,
        grid,
        states,
        Boundaries2d { x: walls, y: walls },
        || false,
    )?;
    let initial = flow.totals();
    let mut support_impulse = [0.0; 2];
    let mut steps = 0;
    let mut maximum_cfl: f64 = 0.0;
    while flow.time().value() < duration {
        if steps >= 20_000 {
            return Err("accepted-step budget exhausted".into());
        }
        let report = flow.step(
            Time::new(duration - flow.time().value()),
            StepControl::default(),
            || false,
        )?;
        maximum_cfl = maximum_cfl.max(report.maximum_cfl);
        for side in [
            BoundarySide::XMin,
            BoundarySide::XMax,
            BoundarySide::YMin,
            BoundarySide::YMax,
        ] {
            for impulse in flow.boundary_impulses(side) {
                for axis in 0..2 {
                    support_impulse[axis] += impulse.momentum[axis];
                }
            }
        }
        steps += 1;
    }
    let final_state = flow.totals();
    eprintln!(
        "N2 NASA9 reflected pressure pulse: grid={n}x{n}, accepted_steps={steps}, time_s={:.12e}, maximum_cfl={maximum_cfl:.6}",
        flow.time().value(),
    );
    eprintln!(
        "mass_residual_kg={:.6e}, referenced_energy_residual_j={:.6e}, support_impulse_kg_m_s=[{:.6e},{:.6e}], momentum_residual_kg_m_s=[{:.6e},{:.6e}]",
        final_state.mass - initial.mass,
        final_state.energy - initial.energy,
        support_impulse[0],
        support_impulse[1],
        final_state.momentum[0] - initial.momentum[0] + support_impulse[0],
        final_state.momentum[1] - initial.momentum[1] + support_impulse[1],
    );
    println!("x_m,y_m,rho_kg_m3,u_m_s,v_m_s,p_pa,T_K");
    for (index, &q) in flow.states().iter().enumerate() {
        let [x, y] = grid.cell_center(index).expect("valid cell index");
        let p = q.primitive(flow.gas(), || false)?;
        println!(
            "{x:.9},{y:.9},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
            p.density, p.velocity[0], p.velocity[1], p.pressure, p.temperature
        );
    }
    Ok(())
}

fn nitrogen() -> Result<FrozenGas, Box<dyn std::error::Error>> {
    // Same retained 200–1000 K NASA-9 coefficients as fs-gas's shock_tube example.
    let a = [
        2.210371497e4,
        -3.818461820e2,
        6.082738360,
        -8.530914410e-3,
        1.384646189e-5,
        -9.625793620e-9,
        2.519705809e-12,
        7.108460860e2,
        -1.076003744e1,
    ];
    let mut parameters = std::collections::BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims().into_iter().take(9).zip(a) {
        parameters.insert(name.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100_000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId("nasa9-standard-state".into()),
        law_version: 1,
        parameters,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 1000.0),
        sources: Vec::new(),
        provenance: Provenance {
            source: "NASA/TP-2002-211556 N2 example association".into(),
            license: "US government technical report".into(),
            artifact: None,
        },
    };
    let convention = StandardStateConventionV1::new(
        StandardStatePhaseV1::Gas,
        ReferenceEquationOfStateV1::IdealGas,
        Pressure::new(100_000.0),
        ElementalReferenceIdV1::new("nasa-glenn-2002-reference-elements")?,
    )?;
    let source = Nasa9StandardStateModelV1::new(
        SpeciesId::new("N2")?,
        MolarMass::new(0.0280134),
        convention,
        vec![Nasa9RegionV1::from_card(card)?],
    )?;
    let mixture = FrozenIdealGasMixtureModelV1::new(
        vec![source],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0])?,
    )?;
    Ok(FrozenGas::new(
        mixture,
        CaloricDomain::new(Temperature::new(200.0), Temperature::new(1000.0)),
        || false,
    )?)
}
