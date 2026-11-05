//! Nonlinear contact-loaded equilibrium, using the existing normal-reaction solve.
use super::*;
use super::super::contact::{ModalContact, ModalContactConfig, contact_column, force_tolerance};
use super::super::contact::multiple::MultiContactConfig;
use super::normal::{admit_setup, solve_joint};
use fs_dcontact::ContactStorage;
use fs_phs::Storage;

struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
}
struct Point {
    column: Vec<f64>,
    storage: ContactStorage,
    config: ModalContactConfig,
    gap: f64,
}
impl Point {
    fn force(&self, closure: f64) -> Result<f64, ModalCouplingError> {
        let mut gradient = [0.0; 2];
        self.storage.gradient(&[-closure, 0.0], &mut gradient);
        let force = finite(-gradient[0])?;
        if force < 0.0 { return Err(invalid("static contact reaction cannot be attractive")); }
        Ok(force)
    }
}

// Contacts condensed over omitted modes: their complete closure columns and the
// omitted modes' static compliance G among them, row-major [m/N]. Each closure
// reads x = N^T q - G R, and the omitted modes store R G R / 2.
pub(super) struct Condensed<'a> {
    pub(super) columns: &'a [Vec<f64>],
    pub(super) compliance: &'a [f64],
}

// Returns the stored energy and the solved reactions.
pub(super) fn initialize(
    network: &mut CoupledModalSystem,
    external: &[f64],
    contacts: &[(ModalContact, ModalContactConfig)],
    condensed: Option<&Condensed<'_>>,
    config: MultiContactConfig,
    gate: &CancelGate,
) -> Result<(f64, Vec<f64>), ModalCouplingError> {
    poll(Some(gate))?;
    network.require_static_initialization(external)?;
    admit_setup(network, contacts.len(), config)?;
    let mut points = Vec::with_capacity(contacts.len());
    for (index, (contact, c)) in contacts.iter().enumerate() {
        poll(Some(gate))?;
        let column = match condensed {
            Some(set) => set.columns[index].clone(),
            None => contact_column(network, contact, *c)?,
        };
        let storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![contact.law.clone()])
            .map_err(ModalCouplingError::ContactLaw)?;
        points.push(Point { column, storage, config: *c, gap: contact.law.gaps()[0] });
    }
    let response = StaticResponse::new(network, gate)?;
    // Contact can restrain this prediction; apply physical caps only to the
    // final contacted equilibrium, not a fictional unrestrained load state.
    let free_q = response.solve(network, external, true, false, gate)?;
    let mut free_x: Vec<f64> = points.iter().map(|p| dot(&p.column, &free_q)).collect::<Result<_,_>>()?;
    let p = points.len();
    let mut compliance = vec![0.0; p*p];
    for j in 0..p {
        poll(Some(gate))?;
        // Rest offsets belong in the free equilibrium, never in a derivative.
        let dq = response.solve(network, &points[j].column, false, false, gate)?;
        for i in 0..p { compliance[i*p+j] = dot(&points[i].column, &dq)?; }
        if compliance[j*p+j] <= 0.0 {
            return Err(invalid("contact attachment needs positive representable static compliance"));
        }
    }
    if let Some(set) = condensed {
        for (a, g) in compliance.iter_mut().zip(set.compliance) { *a = finite(*a+g)?; }
    }
    let mut reactions = vec![0.0; p];
    let mut sweeps = solve_joint(&compliance, &free_x, |i| points[i].config,
        |i,x| points[i].force(x), config.max_sweeps, &mut reactions, Some(gate))?;
    let mut corrected = false;
    let mut energy = loop {
        let mut forces = external.to_vec();
        for (point, reaction) in points.iter().zip(&reactions) {
            if *reaction != 0.0 {
                for (f,b) in forces.iter_mut().zip(&point.column) { *f = finite(*f-b*reaction)?; }
            }
        }
        let q = response.solve(network, &forces, true, true, gate)?;
        let staged = network.stage_static_equilibrium(&q, &forces, gate)?;
        let Some(set) = condensed else { break staged };
        if corrected { break staged; }
        // The solve's closure `free_x - A R` and the actual `N^T q - G R` each cancel
        // terms far larger than the closure, along different sums. A stiff stop turns
        // their rounded difference into a force past its tolerance. Where the gate
        // below would refuse for that, carry the measured difference into `free_x`
        // and solve once more. The law, load, tolerances and gate don't change.
        let mut measured = free_x.clone();
        let mut misses = false;
        for index in 0..p {
            let point = &points[index];
            let lag = dot(&set.compliance[index*p..(index+1)*p], &reactions)?;
            let actual = finite(extension(&network.candidates, &point.column, 0.0)?-lag)?;
            let expected = point.force(actual)?;
            let tolerance = force_tolerance(reactions[index], expected, point.config)?;
            misses |= finite(reactions[index]-expected)?.abs() > tolerance;
            let mut solved = free_x[index];
            for j in 0..p { solved = finite(solved-compliance[index*p+j]*reactions[j])?; }
            measured[index] = finite(free_x[index]+(actual-solved))?;
        }
        if !misses { break staged; }
        free_x = measured;
        corrected = true;
        sweeps = solve_joint(&compliance, &free_x, |i| points[i].config,
            |i,x| points[i].force(x), config.max_sweeps, &mut reactions, Some(gate))?;
    };
    // Recompute contact law on ACTUAL rounded states. In addition to the
    // bilateral force-balance gate above, report no equilibrium if one of
    // these contact equations fails. Neither candidate states nor reactions
    // have been published. There is no relaxation-time or dissipation term
    // in this stationary solve, and no discarded startup vibration.
    let mut actual_forces = external.to_vec();
    let mut actual_scale: Vec<f64> = external.iter().map(|f| f.abs()).collect();
    let mut contact_allowance = vec![0.0; network.mode_count()];
    for (index, (point, reaction)) in points.iter().zip(&reactions).enumerate() {
        poll(Some(gate))?;
        let mut x = extension(&network.candidates, &point.column, 0.0)?;
        if let Some(set) = condensed {
            let lag = dot(&set.compliance[index*p..(index+1)*p], &reactions)?;
            x = finite(x-lag)?;
            energy = finite(energy+0.5*reaction*lag)?;
        }
        let expected = point.force(x)?;
        let residual = finite(reaction-expected)?;
        let tolerance = force_tolerance(*reaction, expected, point.config)?;
        if residual.abs() > tolerance {
            return Err(ModalCouplingError::ContactSolve { residual_n: residual, tolerance_n: tolerance, iterations: sweeps });
        }
        limit("static contact force", expected, point.config.maximum_force_n)?;
        let penetration = finite(x-point.gap)?.max(0.0);
        limit("static contact penetration", penetration, point.config.maximum_penetration_m)?;
        let stored = finite(point.storage.hamiltonian(&[-x, 0.0]))?;
        if stored < 0.0 { return Err(invalid("static contact potential is negative")); }
        energy = finite(energy+stored)?;
        for j in 0..actual_forces.len() {
            let force = finite(point.column[j]*expected)?;
            actual_forces[j] = finite(actual_forces[j]-force)?;
            actual_scale[j] = finite(actual_scale[j]+force.abs())?;
            contact_allowance[j] = finite(contact_allowance[j]+point.column[j].abs()*tolerance)?;
        }
    }
    for (column, link) in network.columns.iter().zip(&network.connections) {
        poll(Some(gate))?;
        let reaction = finite(-link.stiffness_n_m*extension(&network.candidates, column, link.rest_extension_m)?)?;
        // The same rounding scale stage_static_equilibrium gives this link for these states: the extension is a
        // cancelling sum, so a link that ends exactly slack carries a rounded-zero reaction that is balanced, not
        // a 100% residual of its own magnitude.
        let reaction_scale = finite(link.stiffness_n_m*extension_terms(&network.candidates, column, link.rest_extension_m)?)?;
        for j in 0..actual_forces.len() {
            let force = finite(column[j]*reaction)?;
            actual_forces[j] = finite(actual_forces[j]+force)?;
            actual_scale[j] = finite(actual_scale[j]+column[j].abs()*reaction_scale)?;
        }
    }
    let mut j = 0;
    for model in &network.candidates {
        poll(Some(gate))?;
        for (mode, state) in model.modes().iter().zip(model.states()) {
            let actual = finite(mode.angular_frequency_rad_s*mode.angular_frequency_rad_s*state.displacement_m_sqrt_kg)?;
            let scale = finite(actual_scale[j]+actual.abs())?;
            // Contact force tolerances have units N and map through the
            // attachment shapes to N/sqrt(kg). Do not compare these with
            // a dimensionless matrix residual or silently loosen either.
            let tolerance = finite(network.config.solve_relative_tolerance*scale+contact_allowance[j])?;
            if finite(actual-actual_forces[j])?.abs() > tolerance {
                return Err(invalid("actual static contact/network force balance exceeds the mapped force tolerances"));
            }
            j += 1;
        }
    }
    limit("static total energy including contacts", energy, network.config.maximum_total_energy_j)?;
    poll(Some(gate))?;
    std::mem::swap(&mut network.models, &mut network.candidates);
    std::mem::swap(&mut network.material_memory, &mut network.candidate_material_memory);
    Ok((energy, reactions))
}
