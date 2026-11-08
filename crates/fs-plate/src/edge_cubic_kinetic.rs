//! Consistent rotary inertia of the EXISTING edge-cubic displacement field.
//! Through-thickness velocity is (-z*w_dot,x, -z*w_dot,y, w_dot), so its
//! kinetic form is rho*h*int(w_dot^2) + rho*h^3/12*int(|grad w_dot|^2).
//! This is a Rayleigh/Kirchhoff kinetic model, not transverse-shear stiffness.
use super::{edge_cubic_displacement, EDGE_CUBIC_EXP, EDGE_CUBIC_MULT};
use crate::{AssemblyOptions, PlateChart, PlateError, PlateModel, Stiffener, StiffenerMass, TransverseMass};
use fs_sparse::Coo;

/// Exact polynomial area integral of the analytic cubic rotations. The scale
/// is rho*h^3*area/12; geometry/section admission belongs to the assembly owner.
fn rotary_mass(x: &[f64; 3], y: &[f64; 3], scale: f64) -> Result<[f64; 81], PlateError> {
    rotary_mass_inner(x, y, scale, None)
}

fn rotary_mass_inner(x: &[f64; 3], y: &[f64; 3], scale: f64,
    density: Option<[f64; 3]>) -> Result<[f64; 81], PlateError> {
    let bad = || PlateError::BadSection { what: "consistent cubic rotary inertia is not representable" };
    let a2 = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    if !a2.is_finite() || a2 <= 0.0 || !scale.is_finite() || scale <= 0.0 { return Err(bad()); }
    let gradients = [
        [(y[1]-y[2])/a2, (y[2]-y[0])/a2, (y[0]-y[1])/a2],
        [(x[2]-x[1])/a2, (x[0]-x[2])/a2, (x[1]-x[0])/a2],
    ];
    const POWERS: [[i32;3];6] = [[2,0,0],[0,2,0],[0,0,2],[1,1,0],[0,1,1],[1,0,1]];
    const FACT: [f64;5] = [1.0,1.0,2.0,6.0,24.0];
    let control = edge_cubic_displacement(x,y);
    // Quadratic MONOMIAL coefficients, not values or interpolated nodal slopes.
    let mut derivative = [[[0.0;9];6];2];
    for i in 0..10 {
        for axis in 0..3 {
            let power = EDGE_CUBIC_EXP[i][axis];
            if power == 0 { continue; }
            let mut exponent = EDGE_CUBIC_EXP[i]; exponent[axis] -= 1;
            let term = POWERS.iter().position(|p| *p == exponent).expect("quadratic derivative");
            for c in 0..2 {
                let factor = EDGE_CUBIC_MULT[i]*f64::from(power)*gradients[c][axis];
                for dof in 0..9 { derivative[c][term][dof] += factor*control[i][dof]; }
            }
        }
    }
    if derivative.iter().flatten().flatten().any(|v| !v.is_finite()) { return Err(bad()); }
    let mut mass = [0.0;81];
    for row in 0..9 { for col in row..9 {
        let mut value = 0.0;
        for i in 0..6 { for j in 0..6 {
            // int_triangle lambda^alpha / area = 2*prod(alpha!)/(degree+2)!.
            let gram = if let Some(rho) = density {
                crate::reference_density::moment(rho, core::array::from_fn(|a| {
                    (POWERS[i][a]+POWERS[j][a]) as usize
                }))
            } else {
                (2.0/720.0)*(0..3)
                    .map(|a| FACT[(POWERS[i][a]+POWERS[j][a]) as usize]).product::<f64>()
            };
            value += gram*(derivative[0][i][row]*derivative[0][j][col]
                + derivative[1][i][row]*derivative[1][j][col]);
        }}
        value *= scale;
        if !value.is_finite() || (row == col && value <= 0.0) { return Err(bad()); }
        mass[9*row+col] = value; mass[9*col+row] = value;
    }}
    Ok(mass)
}

impl PlateChart {
    /// Assemble the existing cubic transverse mass and replace lumped panel
    /// slope inertia by the exact area integral of its analytic rotations.
    /// The same field can supply pressure, bridge forces and physical motion.
    /// DKT stiffness, prestress, supports and the chosen beam inertia are kept.
    /// This opt-in image does not change any existing assembly method/default.
    ///
    /// The shared assembler owns geometry, section, support and beam admission.
    /// Only its identified panel rotary diagonal is replaced; no beam inertia
    /// is subtracted and no mass, stiffness or eigenvalue is fitted or clamped.
    ///
    /// # Errors
    /// Returns [`PlateError`] on invalid inputs, unrepresentable inertia, or
    /// unresolved subtraction of the original lumped panel rotary diagonal.
    pub fn assemble_edge_cubic_kinetic_mass(&self, stiffeners: &[Stiffener],
        opts: &AssemblyOptions, stiffener_mass: StiffenerMass) -> Result<PlateModel, PlateError> {
        let mut model = self.assemble_with_mass(stiffeners,opts,TransverseMass::EdgeCubic,stiffener_mass)?;
        let sections = self.section_field();
        let mut removed = vec![0.0;model.free];
        let mut consistent = Coo::new(model.free,model.free);
        for (element,tri) in self.mesh.tris.iter().enumerate() {
            let s = sections.at(element);
            let x = tri.map(|node| self.mesh.nodes[node].0);
            let y = tri.map(|node| self.mesh.nodes[node].1);
            let area = 0.5*((x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]));
            let h = s.thickness;
            let density = self.reference_density().map_or([s.density;3], |rho| rho[element]);
            let (lumped, rotary) = if let Some(moments) = self.reference_section_mass() {
                // The independent I2 field owns BOTH the removed slope row
                // sums and their consistent replacement. I0 owns translation.
                let (scale, weights) = crate::reference_density::normalized(
                    moments[element].map(|m| m.second_moment_kg()),
                );
                let total = scale * area;
                (crate::reference_density::lumped(total, weights),
                    rotary_mass_inner(&x, &y, total, Some(weights))?)
            } else if density[0] == density[1] && density[1] == density[2] {
                // Same operation order as the original assembler's slope mass.
                let total = density[0]*h*h*h/12.0*area;
                ([total/3.0;3], rotary_mass(&x,&y,total)?)
            } else {
                let (scale,rho) = crate::reference_density::normalized(density);
                let total = scale*h*h*h/12.0*area;
                (crate::reference_density::lumped(total,rho),
                    rotary_mass_inner(&x,&y,total,Some(rho))?)
            };
            for local in 0..9 {
                if local%3 != 0 {
                    if let Some(index) = model.dof_map[3*tri[local/3]+local%3] {
                        removed[index] += lumped[local/3];
                    }
                }
            }
            for row in 0..9 { for col in 0..9 {
                if let (Some(r),Some(c)) = (model.dof_map[3*tri[row/3]+row%3],
                    model.dof_map[3*tri[col/3]+col%3]) {
                    consistent.push(r,c,rotary[9*row+col]);
                }
            }}
        }
        for row in 0..model.free {
            let (columns,values) = model.m.row(row);
            for (&col,&original) in columns.iter().zip(values) {
                let value = original - if row == col { removed[row] } else { 0.0 };
                if !value.is_finite() || (row == col && value <= 0.0) {
                    return Err(PlateError::BadSection {
                        what: "cubic transverse/beam inertia cannot be resolved after removing lumped panel rotation",
                    });
                }
                consistent.push(row,col,value);
            }
        }
        model.m = consistent.assemble();
        if (0..model.free).any(|r| model.m.row(r).1.iter().any(|v| !v.is_finite())) {
            return Err(PlateError::BadSection { what: "assembled consistent cubic kinetic mass overflow" });
        }
        Ok(model)
    }
}

#[cfg(test)]
#[path = "edge_cubic_kinetic_tests.rs"]
mod tests;
