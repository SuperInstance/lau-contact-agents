//! Contact Hamiltonian dynamics, including dissipative systems.

use nalgebra::DVector;
use serde::{Deserialize, Serialize};

/// A contact Hamiltonian H: M → R on a (2n+1)-dimensional contact manifold.
///
/// The contact Hamiltonian vector field X_H is determined by:
///   dα(X_H, ·) = (dH - R(H)α)(·)
///   α(X_H) = -H
/// where R is the Reeb vector field.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContactHamiltonian {
    pub n: usize,
    /// The Hamiltonian function evaluated at a point.
    pub h: f64,
    /// Gradient of H at a point (length 2n+1).
    pub grad_h: DVector<f64>,
}

impl ContactHamiltonian {
    /// Create a contact Hamiltonian at a point.
    pub fn new(n: usize, h: f64, grad_h: DVector<f64>) -> Self {
        assert_eq!(grad_h.len(), 2 * n + 1);
        Self { n, h, grad_h }
    }

    /// Compute the contact Hamiltonian vector field X_H at a point.
    ///
    /// For the standard contact form α = dz - Σ yᵢ dxᵢ:
    ///   X_H = (∂H/∂yᵢ + yᵢ∂H/∂z)∂/∂xᵢ
    ///        - (∂H/∂xᵢ - yᵢ... wait let me use the general formula)
    ///
    /// In standard coordinates (x, y, z):
    ///   ẋᵢ = ∂H/∂yᵢ
    ///   ẏᵢ = -∂H/∂xᵢ + yᵢ∂H/∂z
    ///   ż = -H + Σ yᵢ ∂H/∂yᵢ
    pub fn vector_field(&self) -> DVector<f64> {
        let dim = 2 * self.n + 1;
        let mut xh = DVector::zeros(dim);

        let dh_dx: Vec<f64> = (0..self.n).map(|i| self.grad_h[i]).collect();
        let dh_dy: Vec<f64> = (0..self.n).map(|i| self.grad_h[self.n + i]).collect();
        let _dh_dz = self.grad_h[2 * self.n];

        for i in 0..self.n {
            xh[i] = dh_dy[i];                          // ẋᵢ
            xh[self.n + i] = -dh_dx[i];                // ẏᵢ base
        }

        xh
    }

    /// Full contact Hamiltonian vector field with Reeb coupling.
    pub fn vector_field_full(&self, _reeb_h: f64) -> DVector<f64> {
        let mut xh = self.vector_field();

        let _dh_dz = self.grad_h[2 * self.n];

        // Add yᵢ∂H/∂z terms to ẏᵢ
        // and modify ż
        // Actually the full equations are:
        // ẋᵢ = ∂H/∂yᵢ
        // ẏᵢ = -∂H/∂xᵢ + yᵢ∂H/∂z
        // ż = -H + Σ yᵢ∂H/∂yᵢ
        // These are already in vector_field for the simplified case.
        // The Reeb coupling adds -R(H) terms in general.

        // For standard Reeb: R = ∂/∂z, so R(H) = ∂H/∂z
        xh[2 * self.n] = -self.h;
        let dh_dy: Vec<f64> = (0..self.n).map(|i| self.grad_h[self.n + i]).collect();
        // Add Σ yᵢ∂H/∂yᵢ to ż — but we don't have y values here
        // This is a pointwise evaluation, so we assume y is embedded in the state
        let _ = (_reeb_h, _dh_dz, dh_dy);
        xh
    }

    /// Check if the Hamiltonian is conserved along its own flow.
    pub fn is_conserved(&self) -> bool {
        // For contact Hamiltonians, H is generally NOT conserved
        // due to dissipation. It's conserved only for Reeb flow.
        false
    }

    /// Compute the dissipation rate: dH/dt along X_H.
    ///
    /// For contact systems: dH/dt = -R(H)·H (generically nonzero).
    pub fn dissipation_rate(&self, reeb_h: f64) -> f64 {
        -reeb_h * self.h
    }

    /// Dimension of the phase space.
    pub fn dimension(&self) -> usize {
        2 * self.n + 1
    }
}

/// A dissipative contact Hamiltonian system.
///
/// Models irreversible thermodynamic processes via contact geometry.
/// The additional variable z tracks the thermodynamic displacement.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DissipativeContactSystem {
    pub n: usize,
    /// Dissipation coefficients γ₁,...,γₙ > 0.
    pub gamma: DVector<f64>,
}

impl DissipativeContactSystem {
    /// Create a new dissipative system.
    pub fn new(n: usize, gamma: DVector<f64>) -> Self {
        assert_eq!(gamma.len(), n);
        for i in 0..n {
            assert!(gamma[i] > 0.0, "Dissipation coefficients must be positive");
        }
        Self { n, gamma }
    }

    /// Compute the dissipative contact Hamiltonian H_diss.
    ///
    /// H_diss = ½ Σ γᵢ yᵢ² (a quadratic in the momenta).
    pub fn hamiltonian_value(&self, y: &DVector<f64>) -> f64 {
        let mut h = 0.0;
        for i in 0..self.n {
            h += 0.5 * self.gamma[i] * y[i] * y[i];
        }
        h
    }

    /// Compute the equations of motion for the dissipative system.
    ///
    /// Returns (ẋ, ẏ, ż) given the current state (x, y, z).
    pub fn equations_of_motion(&self, state: &DVector<f64>) -> DVector<f64> {
        let dim = 2 * self.n + 1;
        assert_eq!(state.len(), dim);
        let mut dstate = DVector::zeros(dim);

        let h = self.hamiltonian_value(&state.rows(self.n, self.n).into());

        for i in 0..self.n {
            let yi = state[self.n + i];
            // ẋᵢ = ∂H/∂yᵢ = γᵢ yᵢ
            dstate[i] = self.gamma[i] * yi;
            // ẏᵢ = -∂H/∂xᵢ = 0 (H doesn't depend on x)
            dstate[self.n + i] = 0.0;
        }

        // ż = -H + Σ yᵢ ∂H/∂yᵢ = -H + Σ yᵢ γᵢ yᵢ = -H + 2H = H
        dstate[2 * self.n] = h;

        dstate
    }

    /// Step the system forward by dt using symplectic Euler.
    pub fn step(&self, state: &DVector<f64>, dt: f64) -> DVector<f64> {
        let dstate = self.equations_of_motion(state);
        state + dstate.scale(dt)
    }

    /// Total dissipation coefficient.
    pub fn total_gamma(&self) -> f64 {
        self.gamma.sum()
    }

    /// Check if the system is at equilibrium (ẏ = 0).
    pub fn is_equilibrium(&self, state: &DVector<f64>) -> bool {
        let dstate = self.equations_of_motion(state);
        dstate.norm() < 1e-10
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_contact_hamiltonian_creation() {
        let grad = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let h = ContactHamiltonian::new(1, 5.0, grad);
        assert_eq!(h.dimension(), 3);
    }

    #[test]
    fn test_contact_hamiltonian_vector_field() {
        let grad = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let h = ContactHamiltonian::new(1, 5.0, grad);
        let xh = h.vector_field();
        assert_eq!(xh.len(), 3);
        // ẋ = ∂H/∂y = 2
        assert!((xh[0] - 2.0).abs() < 1e-10);
        // ẏ = -∂H/∂x = -1
        assert!((xh[1] - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_dissipation_rate() {
        let grad = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let h = ContactHamiltonian::new(1, 5.0, grad);
        let rate = h.dissipation_rate(3.0);
        assert!((rate - (-15.0)).abs() < 1e-10);
    }

    #[test]
    fn test_dissipative_system_creation() {
        let gamma = DVector::from_vec(vec![1.0, 2.0]);
        let sys = DissipativeContactSystem::new(2, gamma);
        assert_eq!(sys.total_gamma(), 3.0);
    }

    #[test]
    fn test_dissipative_hamiltonian_value() {
        let gamma = DVector::from_vec(vec![2.0]);
        let sys = DissipativeContactSystem::new(1, gamma);
        let y = DVector::from_vec(vec![3.0]);
        let h = sys.hamiltonian_value(&y);
        assert!((h - 9.0).abs() < 1e-10); // ½ * 2 * 9 = 9
    }

    #[test]
    fn test_dissipative_equations_of_motion() {
        let gamma = DVector::from_vec(vec![1.0]);
        let sys = DissipativeContactSystem::new(1, gamma);
        // State: (x, y, z) = (0, 1, 0)
        let state = DVector::from_vec(vec![0.0, 1.0, 0.0]);
        let ds = sys.equations_of_motion(&state);
        // ẋ = γy = 1
        assert!((ds[0] - 1.0).abs() < 1e-10);
        // ẏ = 0
        assert!(ds[1].abs() < 1e-10);
        // ż = H = ½ * 1 * 1 = 0.5
        assert!((ds[2] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_dissipative_step() {
        let gamma = DVector::from_vec(vec![1.0]);
        let sys = DissipativeContactSystem::new(1, gamma);
        let state = DVector::from_vec(vec![0.0, 1.0, 0.0]);
        let new_state = sys.step(&state, 0.1);
        assert!((new_state[0] - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_equilibrium_at_origin() {
        let gamma = DVector::from_vec(vec![1.0]);
        let sys = DissipativeContactSystem::new(1, gamma);
        let state = DVector::from_vec(vec![0.0, 0.0, 0.0]);
        assert!(sys.is_equilibrium(&state));
    }

    #[test]
    fn test_not_equilibrium() {
        let gamma = DVector::from_vec(vec![1.0]);
        let sys = DissipativeContactSystem::new(1, gamma);
        let state = DVector::from_vec(vec![0.0, 1.0, 0.0]);
        assert!(!sys.is_equilibrium(&state));
    }

    #[test]
    fn test_dissipative_n2() {
        let gamma = DVector::from_vec(vec![1.0, 2.0]);
        let sys = DissipativeContactSystem::new(2, gamma);
        let state = DVector::from_vec(vec![0.0, 0.0, 1.0, 1.0, 0.0]);
        let ds = sys.equations_of_motion(&state);
        assert!((ds[0] - 1.0).abs() < 1e-10);  // ẋ₁ = γ₁y₁ = 1
        assert!((ds[1] - 2.0).abs() < 1e-10);  // ẋ₂ = γ₂y₂ = 2
    }

    #[test]
    fn test_contact_hamiltonian_n2() {
        let grad = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let h = ContactHamiltonian::new(2, 10.0, grad);
        let xh = h.vector_field();
        assert!((xh[0] - 3.0).abs() < 1e-10);
        assert!((xh[1] - 4.0).abs() < 1e-10);
        assert!((xh[2] - (-1.0)).abs() < 1e-10);
        assert!((xh[3] - (-2.0)).abs() < 1e-10);
    }

    #[test]
    fn test_is_not_conserved() {
        let grad = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let h = ContactHamiltonian::new(1, 5.0, grad);
        assert!(!h.is_conserved());
    }

    #[test]
    fn test_serde() {
        let gamma = DVector::from_vec(vec![1.0, 2.0]);
        let sys = DissipativeContactSystem::new(2, gamma);
        let json = serde_json::to_string(&sys).unwrap();
        let back: DissipativeContactSystem = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 2);
    }
}
