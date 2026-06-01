//! Thermodynamic agent dynamics using contact geometry.
//!
//! Models thermodynamic systems (agents) on contact manifolds where:
//! - Position variables x represent extensive quantities
//! - Momentum variables y represent intensive quantities
//! - z tracks the thermodynamic displacement (entropy production)

use nalgebra::{DVector, DMatrix};
use serde::{Deserialize, Serialize};

/// A thermodynamic agent living on a contact manifold.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermoAgent {
    pub n: usize,
    /// Agent state (x₁,...,xₙ, y₁,...,yₙ, z) in the contact manifold.
    pub state: DVector<f64>,
    /// Internal energy function coefficients.
    pub energy_coeffs: DVector<f64>,
}

impl ThermoAgent {
    /// Create a new thermodynamic agent.
    pub fn new(n: usize, state: DVector<f64>) -> Self {
        assert_eq!(state.len(), 2 * n + 1);
        let mut energy_coeffs = DVector::zeros(n);
        for i in 0..n {
            energy_coeffs[i] = 1.0;
        }
        Self { n, state, energy_coeffs }
    }

    /// Create with custom energy coefficients.
    pub fn with_energy(n: usize, state: DVector<f64>, energy_coeffs: DVector<f64>) -> Self {
        assert_eq!(state.len(), 2 * n + 1);
        assert_eq!(energy_coeffs.len(), n);
        Self { n, state, energy_coeffs }
    }

    /// Get the extensive variables x.
    pub fn extensive(&self) -> DVector<f64> {
        self.state.rows(0, self.n).into()
    }

    /// Get the intensive variables y.
    pub fn intensive(&self) -> DVector<f64> {
        self.state.rows(self.n, self.n).into()
    }

    /// Get the entropy-like variable z.
    pub fn entropy(&self) -> f64 {
        self.state[2 * self.n]
    }

    /// Compute the internal energy E(x) = Σ cᵢ xᵢ².
    pub fn internal_energy(&self) -> f64 {
        let x = self.extensive();
        let mut e = 0.0;
        for i in 0..self.n {
            e += self.energy_coeffs[i] * x[i] * x[i];
        }
        e
    }

    /// Compute the free energy F = E - TS where T is temperature-like.
    pub fn free_energy(&self, temperature: f64) -> f64 {
        self.internal_energy() - temperature * self.entropy()
    }

    /// Compute entropy production rate.
    pub fn entropy_production(&self, gamma: &DVector<f64>) -> f64 {
        let y = self.intensive();
        let mut s_dot = 0.0;
        for i in 0..self.n {
            s_dot += gamma[i] * y[i] * y[i];
        }
        s_dot
    }

    /// Evolve the agent state by one step using contact Hamiltonian dynamics.
    pub fn evolve(&mut self, dt: f64, gamma: &DVector<f64>) {
        assert_eq!(gamma.len(), self.n);
        let y = self.intensive();
        let x = self.extensive();
        let z = self.entropy();

        // Contact Hamiltonian equations:
        // ẋᵢ = γᵢ yᵢ (dissipative flow)
        // ẏᵢ = -∂E/∂xᵢ = -2cᵢxᵢ
        // ż = H - Σ yᵢ ẋᵢ

        let mut new_state = DVector::zeros(2 * self.n + 1);

        for i in 0..self.n {
            new_state[i] = self.state[i] + gamma[i] * y[i] * dt;
            new_state[self.n + i] = self.state[self.n + i]
                - 2.0 * self.energy_coeffs[i] * x[i] * dt;
        }

        let h: f64 = (0..self.n)
            .map(|i| 0.5 * gamma[i] * y[i] * y[i])
            .sum();
        let y_dot_x: f64 = (0..self.n)
            .map(|i| y[i] * gamma[i] * y[i])
            .sum();
        new_state[2 * self.n] = z + (h - y_dot_x) * dt;

        self.state = new_state;
    }

    /// Check if the agent is at thermal equilibrium.
    pub fn is_equilibrium(&self) -> bool {
        let y = self.intensive();
        y.norm() < 1e-10
    }

    /// Compute the thermodynamic length (metric on state space).
    pub fn thermodynamic_length(&self, other: &ThermoAgent) -> f64 {
        let diff = &self.state - &other.state;
        // Use the Riemannian metric g_ij = δ_ij on the contact manifold
        diff.norm()
    }

    /// Total dimension of the contact manifold.
    pub fn dimension(&self) -> usize {
        2 * self.n + 1
    }
}

/// A system of interacting thermodynamic agents.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThermoMultiAgent {
    pub agents: Vec<ThermoAgent>,
    /// Coupling matrix between agents.
    pub coupling: DMatrix<f64>,
}

impl ThermoMultiAgent {
    /// Create a multi-agent system.
    pub fn new(agents: Vec<ThermoAgent>, coupling: DMatrix<f64>) -> Self {
        Self { agents, coupling }
    }

    /// Create with no coupling.
    pub fn decoupled(agents: Vec<ThermoAgent>) -> Self {
        let n_agents = agents.len();
        Self {
            agents,
            coupling: DMatrix::zeros(n_agents, n_agents),
        }
    }

    /// Number of agents.
    pub fn num_agents(&self) -> usize {
        self.agents.len()
    }

    /// Total entropy production across all agents.
    pub fn total_entropy_production(&self, gamma: &DVector<f64>) -> f64 {
        self.agents.iter().map(|a| a.entropy_production(gamma)).sum()
    }

    /// Evolve all agents for one step.
    pub fn evolve(&mut self, dt: f64, gamma: &DVector<f64>) {
        // Apply coupling forces first
        let n = self.agents.len();
        if n > 1 {
            let n_state = self.agents[0].dimension();
            // Compute mean state for coupling
            let mut mean_state = DVector::zeros(n_state);
            for a in &self.agents {
                mean_state += &a.state;
            }
            mean_state /= n as f64;

            // Apply coupling perturbation
            for (i, a) in self.agents.iter_mut().enumerate() {
                let coupling_strength = self.coupling.row(i).sum();
                let perturbation = (&mean_state - &a.state).scale(coupling_strength * dt * 0.01);
                a.state += &perturbation;
            }
        }

        // Evolve each agent
        for a in &mut self.agents {
            a.evolve(dt, gamma);
        }
    }

    /// Check if all agents are at equilibrium.
    pub fn all_equilibrium(&self) -> bool {
        self.agents.iter().all(|a| a.is_equilibrium())
    }

    /// Total internal energy.
    pub fn total_energy(&self) -> f64 {
        self.agents.iter().map(|a| a.internal_energy()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thermo_agent_creation() {
        let state = DVector::from_vec(vec![1.0, 0.0, 0.0]);
        let agent = ThermoAgent::new(1, state);
        assert_eq!(agent.dimension(), 3);
    }

    #[test]
    fn test_extensive_intensive() {
        let state = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let agent = ThermoAgent::new(2, state);
        let x = agent.extensive();
        let y = agent.intensive();
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 2.0).abs() < 1e-10);
        assert!((y[0] - 3.0).abs() < 1e-10);
        assert!((y[1] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_entropy() {
        let state = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let agent = ThermoAgent::new(1, state);
        assert!((agent.entropy() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_internal_energy() {
        let state = DVector::from_vec(vec![2.0, 0.0, 0.0]);
        let agent = ThermoAgent::new(1, state);
        assert!((agent.internal_energy() - 4.0).abs() < 1e-10); // c₁ * x₁² = 1 * 4
    }

    #[test]
    fn test_free_energy() {
        let state = DVector::from_vec(vec![2.0, 0.0, 5.0]);
        let agent = ThermoAgent::new(1, state);
        let f = agent.free_energy(1.0);
        assert!((f - (-1.0)).abs() < 1e-10); // 4 - 1*5 = -1
    }

    #[test]
    fn test_entropy_production() {
        let state = DVector::from_vec(vec![0.0, 2.0, 0.0]);
        let agent = ThermoAgent::new(1, state);
        let gamma = DVector::from_vec(vec![1.0]);
        let s = agent.entropy_production(&gamma);
        assert!((s - 4.0).abs() < 1e-10); // γ * y² = 1 * 4
    }

    #[test]
    fn test_evolve() {
        let state = DVector::from_vec(vec![0.0, 1.0, 0.0]);
        let mut agent = ThermoAgent::new(1, state);
        let gamma = DVector::from_vec(vec![1.0]);
        agent.evolve(0.1, &gamma);
        // ẋ = γy = 1, so x → 0.1
        assert!((agent.state[0] - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_equilibrium() {
        let state = DVector::from_vec(vec![1.0, 0.0, 0.0]);
        let agent = ThermoAgent::new(1, state);
        assert!(agent.is_equilibrium());
    }

    #[test]
    fn test_not_equilibrium() {
        let state = DVector::from_vec(vec![1.0, 1.0, 0.0]);
        let agent = ThermoAgent::new(1, state);
        assert!(!agent.is_equilibrium());
    }

    #[test]
    fn test_thermodynamic_length() {
        let s1 = DVector::from_vec(vec![0.0, 0.0, 0.0]);
        let s2 = DVector::from_vec(vec![3.0, 4.0, 0.0]);
        let a1 = ThermoAgent::new(1, s1);
        let a2 = ThermoAgent::new(1, s2);
        assert!((a1.thermodynamic_length(&a2) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_multi_agent_creation() {
        let a1 = ThermoAgent::new(1, DVector::from_vec(vec![0.0, 0.0, 0.0]));
        let a2 = ThermoAgent::new(1, DVector::from_vec(vec![1.0, 1.0, 0.0]));
        let sys = ThermoMultiAgent::decoupled(vec![a1, a2]);
        assert_eq!(sys.num_agents(), 2);
    }

    #[test]
    fn test_multi_agent_total_energy() {
        let a1 = ThermoAgent::new(1, DVector::from_vec(vec![2.0, 0.0, 0.0]));
        let a2 = ThermoAgent::new(1, DVector::from_vec(vec![3.0, 0.0, 0.0]));
        let sys = ThermoMultiAgent::decoupled(vec![a1, a2]);
        assert!((sys.total_energy() - 13.0).abs() < 1e-10); // 4 + 9
    }

    #[test]
    fn test_multi_agent_evolve() {
        let a1 = ThermoAgent::new(1, DVector::from_vec(vec![0.0, 1.0, 0.0]));
        let a2 = ThermoAgent::new(1, DVector::from_vec(vec![1.0, 0.0, 0.0]));
        let mut sys = ThermoMultiAgent::decoupled(vec![a1, a2]);
        let gamma = DVector::from_vec(vec![1.0]);
        sys.evolve(0.1, &gamma);
        assert_eq!(sys.num_agents(), 2);
    }

    #[test]
    fn test_total_entropy_production() {
        let a1 = ThermoAgent::new(1, DVector::from_vec(vec![0.0, 1.0, 0.0]));
        let a2 = ThermoAgent::new(1, DVector::from_vec(vec![0.0, 2.0, 0.0]));
        let sys = ThermoMultiAgent::decoupled(vec![a1, a2]);
        let gamma = DVector::from_vec(vec![1.0]);
        let s = sys.total_entropy_production(&gamma);
        assert!((s - 5.0).abs() < 1e-10); // 1 + 4
    }

    #[test]
    fn test_serde() {
        let state = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let agent = ThermoAgent::new(1, state);
        let json = serde_json::to_string(&agent).unwrap();
        let back: ThermoAgent = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 1);
    }
}
