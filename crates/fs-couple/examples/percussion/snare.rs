//! A bank of distinct tensioned strands attached to the EXISTING resonant head.
//! Commercial anchor: 20 steel-coil strands on a nominal 14-inch PureSound
//! Custom Pro. All detailed coil, installation and constitutive values below
//! are editable estimates, not a reverse-engineered or calibrated product.
use super::Error;
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::ImpactBody;
use fs_couple::render::plate::impact::linear::wire::{HelicalWire, LineContact, WireSpan, film_shapes};
use fs_dcontact::Obstacle;
use fs_plate::{ModePair, shell::head::TensionedDisk};
use std::ops::Range;

#[derive(Clone, Copy, Debug)]
pub struct SnareSet {
    pub strands: usize,
    pub modes_per_strand: usize,
    pub contact_cells: usize,
    pub length_m: f64,
    pub width_m: f64,
    pub coil: HelicalWire,
    pub tension_per_strand_n: f64,
    pub bending_per_strand_n_m2: f64,
    pub damping_per_s: f64,
    pub clearance_m: f64,
    pub contact_stiffness_per_length: f64,
    pub contact_exponent: f64,
    pub contact_internal_loss_s_m: f64,
}
impl SnareSet {
    pub fn reference(disengaged: bool) -> Self {
        Self { strands: 20, modes_per_strand: 8, contact_cells: 12,
            length_m: 0.30, width_m: 0.04,
            coil: HelicalWire { wire_radius_m: 0.00015, coil_radius_m: 0.00055,
                pitch_m: 0.00085, density_kg_m3: 7800.0 },
            tension_per_strand_n: 0.7, bending_per_strand_n_m2: 1e-6, damping_per_s: 4.0,
            clearance_m: if disengaged {0.003} else {0.00002},
            contact_stiffness_per_length: 5e8, contact_exponent: 1.5, contact_internal_loss_s_m: 0.05 }
    }
    pub fn mode_count(self) -> Result<usize, Error> {
        if !(1..=24).contains(&self.strands) || !(1..=16).contains(&self.modes_per_strand)
            || self.contact_cells < self.modes_per_strand+1 || self.contact_cells > 20
            || !self.width_m.is_finite() || self.width_m < 0.0
            || !self.length_m.is_finite() || self.length_m <= 0.0
        { return Err("snare reference requires 1..24 strands, 1..16 modes and explicit bounded contact sampling".into()); }
        Ok(self.strands*self.modes_per_strand)
    }
    /// Every strand is its own body, every station has its own force, and all
    /// of them share the SAME accepted drumhead coordinates in one joint solve.
    pub fn assemble(self, film: &TensionedDisk, modes: &[ModePair], receiver: Range<usize>,
        first_wire: usize, total_modes: usize,
    ) -> Result<(Vec<ImpactBody>, Vec<Obstacle>), Error> {
        let extra = self.mode_count()?;
        if first_wire.checked_add(extra) != Some(total_modes) || receiver.end > first_wire
            || receiver.start >= receiver.end || receiver.len() != modes.len()
        { return Err("snare body and receiver coordinate ranges are inconsistent".into()); }
        let mu = self.coil.linear_density_kg_m()?;
        let mut bodies = Vec::with_capacity(self.strands);
        let mut contacts = Vec::with_capacity(self.strands);
        for strand in 0..self.strands {
            let y = if self.strands == 1 {0.0} else {
                self.width_m*(strand as f64/(self.strands-1) as f64-0.5)
            };
            let wire = WireSpan { endpoints_m: [[-0.5*self.length_m,y],[0.5*self.length_m,y]],
                linear_density_kg_m: mu, tension_n: self.tension_per_strand_n,
                bending_stiffness_n_m2: self.bending_per_strand_n_m2,
                damping_per_s: vec![self.damping_per_s;self.modes_per_strand] };
            let line = LineContact::uniform(wire.length_m(), self.contact_cells, self.clearance_m,
                self.contact_stiffness_per_length, self.contact_exponent, self.contact_internal_loss_s_m,
                format!("estimated homogenized steel-coil snare strand {strand}; per-metre law, not measured commercial contact"))?;
            let positions = wire.positions(&line)?;
            let shapes = film_shapes(film, modes, &positions)?;
            let start = first_wire+strand*self.modes_per_strand;
            contacts.push(wire.contact(&line, &shapes, receiver.clone(),
                start..start+self.modes_per_strand, total_modes)?);
            bodies.push(wire.body(vec![ModalAcousticState::default();self.modes_per_strand])?);
        }
        eprintln!("snare reference: strands={}, modal_coordinates={}, contact_points={}, line_mass_kg_m={}, total_wire_mass_kg={}, clearance_m={}; coil/tension/loss are estimates; no direct wire radiation", self.strands, extra,
            self.strands*self.contact_cells, mu, mu*self.length_m*self.strands as f64, self.clearance_m);
        Ok((bodies,contacts))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reference_keeps_twenty_distinct_strands_and_explicit_line_sampling() {
        let spec=SnareSet::reference(false);
        assert_eq!(spec.strands,20);assert_eq!(spec.mode_count().unwrap(),160);
        assert_eq!(spec.strands*spec.contact_cells,240);
        assert!(SnareSet::reference(true).clearance_m > spec.clearance_m);
        assert!(SnareSet{contact_cells:2,..spec}.mode_count().is_err());
    }
    #[test]
    fn actual_head_and_wire_assembly_retains_all_contact_rows_and_unexcited_wire_states() {
        // Uses the original film assembly/eigensolve and prepared modal host.
        let experiment=super::super::drum_with_wires(8,2e-6,false,true,
            Some(SnareSet::reference(false))).unwrap();
        let super::super::mechanics::Mechanics::Prepared(system)=&experiment.system else {panic!()};
        assert_eq!(system.contact_count(),241);
        let n=experiment.force.len();
        assert!(system.state()[2*(n-160)..].iter().all(|v|*v==0.0));
        assert!(experiment.pressure.as_ref().unwrap().areas[n-160..].iter().all(|a|*a==0.0));
    }
}
