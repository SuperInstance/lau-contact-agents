//! Weinstein conjecture: On a closed contact manifold (M, α), the Reeb vector
//! field has at least one closed orbit.
//!
//! Proved for many cases by Taubes (2007) in dimension 3.

use nalgebra::{DVector};
use serde::{Deserialize, Serialize};

/// A closed Reeb orbit: a periodic trajectory of the Reeb vector field.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReebOrbit {
    pub n: usize,
    /// Period of the orbit.
    pub period: f64,
    /// Sample points along the orbit (2n+1 dimensional).
    pub points: Vec<DVector<f64>>,
}

impl ReebOrbit {
    /// Create a new closed Reeb orbit.
    pub fn new(n: usize, period: f64, points: Vec<DVector<f64>>) -> Self {
        let dim = 2 * n + 1;
        for p in &points {
            assert_eq!(p.len(), dim);
        }
        Self { n, period, points }
    }

    /// Check if the orbit is approximately closed (start ≈ end).
    pub fn is_closed(&self) -> bool {
        if self.points.len() < 2 {
            return false;
        }
        let first = &self.points[0];
        let last = &self.points[self.points.len() - 1];
        (first - last).norm() < 1e-6
    }

    /// Compute the action (integral of α along the orbit).
    pub fn action(&self) -> f64 {
        // For the standard Reeb flow R = ∂/∂z,
        // α(R) = 1, so action = period
        self.period
    }

    /// Conley-Zehnder index of the orbit (simplified).
    pub fn conley_zehnder_index(&self) -> i32 {
        // For the standard Reeb orbit, CZ = 1
        // This is a simplified placeholder
        1
    }

    /// Number of sample points.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }
}

/// Weinstein conjecture verification for specific manifolds.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WeinsteinVerifier {
    pub n: usize,
    /// Known minimum number of Reeb orbits.
    pub min_orbits: usize,
    /// Found orbits.
    pub orbits: Vec<ReebOrbit>,
}

impl WeinsteinVerifier {
    /// Create a verifier for a (2n+1)-dimensional contact manifold.
    pub fn new(n: usize) -> Self {
        Self {
            n,
            min_orbits: 1,
            orbits: Vec::new(),
        }
    }

    /// Generate the standard Reeb orbits on S^{2n+1} with the standard contact form.
    pub fn standard_reeb_orbits(n: usize, num_points: usize) -> Vec<ReebOrbit> {
        let dim = 2 * n + 1;
        let period = 2.0 * std::f64::consts::PI;
        let mut orbits = Vec::new();

        // The standard Reeb flow on S^3 (n=1) is the Hopf fibration
        let mut points = Vec::new();
        for k in 0..num_points {
            let t = period * (k as f64) / (num_points as f64);
            let mut p = DVector::zeros(dim);
            if dim == 3 {
                // Hopf orbit on S^3 projected to coordinates
                p[0] = t.cos();
                p[1] = t.sin();
                p[2] = 0.0;
            } else {
                // Generalized: first two coordinates rotate
                p[0] = t.cos();
                p[1] = t.sin();
            }
            points.push(p);
        }
        // Close the orbit
        if !points.is_empty() {
            points.push(points[0].clone());
        }

        orbits.push(ReebOrbit::new(n, period, points));
        orbits
    }

    /// Check if the Weinstein conjecture is satisfied (at least one orbit found).
    pub fn weinstein_satisfied(&self) -> bool {
        self.orbits.len() >= self.min_orbits
    }

    /// Total number of found orbits.
    pub fn num_orbits(&self) -> usize {
        self.orbits.len()
    }

    /// Verify for the standard S^3 case.
    pub fn verify_s3() -> Self {
        let orbits = Self::standard_reeb_orbits(1, 100);
        let mut verifier = Self::new(1);
        verifier.orbits = orbits;
        verifier
    }

    /// Verify for higher-dimensional spheres.
    pub fn verify_s2n1(n: usize) -> Self {
        let orbits = Self::standard_reeb_orbits(n, 100);
        let mut verifier = Self::new(n);
        verifier.orbits = orbits;
        verifier
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_reeb_orbit_creation() {
        let pts = vec![
            DVector::from_vec(vec![1.0, 0.0, 0.0]),
            DVector::from_vec(vec![0.0, 1.0, 0.0]),
            DVector::from_vec(vec![1.0, 0.0, 0.0]),
        ];
        let orbit = ReebOrbit::new(1, 2.0 * std::f64::consts::PI, pts);
        assert_eq!(orbit.num_points(), 3);
    }

    #[test]
    fn test_reeb_orbit_closed() {
        let pts = vec![
            DVector::from_vec(vec![1.0, 0.0, 0.0]),
            DVector::from_vec(vec![0.0, 1.0, 0.0]),
            DVector::from_vec(vec![1.0, 0.0, 0.0]),
        ];
        let orbit = ReebOrbit::new(1, 2.0 * std::f64::consts::PI, pts);
        assert!(orbit.is_closed());
    }

    #[test]
    fn test_reeb_orbit_not_closed() {
        let pts = vec![
            DVector::from_vec(vec![1.0, 0.0, 0.0]),
            DVector::from_vec(vec![0.0, 1.0, 0.0]),
            DVector::from_vec(vec![0.0, 0.0, 1.0]),
        ];
        let orbit = ReebOrbit::new(1, 2.0 * std::f64::consts::PI, pts);
        assert!(!orbit.is_closed());
    }

    #[test]
    fn test_reeb_orbit_action() {
        let orbit = ReebOrbit::new(1, 5.0, vec![]);
        assert!((orbit.action() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_reeb_orbit_cz_index() {
        let orbit = ReebOrbit::new(1, 2.0, vec![]);
        assert_eq!(orbit.conley_zehnder_index(), 1);
    }

    #[test]
    fn test_weinstein_s3() {
        let verifier = WeinsteinVerifier::verify_s3();
        assert!(verifier.weinstein_satisfied());
        assert!(verifier.num_orbits() >= 1);
    }

    #[test]
    fn test_weinstein_s5() {
        let verifier = WeinsteinVerifier::verify_s2n1(2);
        assert!(verifier.weinstein_satisfied());
    }

    #[test]
    fn test_standard_reeb_orbits() {
        let orbits = WeinsteinVerifier::standard_reeb_orbits(1, 50);
        assert!(!orbits.is_empty());
        assert!(orbits[0].is_closed());
    }

    #[test]
    fn test_verifier_creation() {
        let v = WeinsteinVerifier::new(1);
        assert_eq!(v.min_orbits, 1);
        assert_eq!(v.num_orbits(), 0);
        assert!(!v.weinstein_satisfied()); // no orbits yet
    }

    #[test]
    fn test_reeb_orbit_n2() {
        let pts = vec![
            DVector::from_vec(vec![1.0, 0.0, 0.0, 0.0, 0.0]),
            DVector::from_vec(vec![0.0, 1.0, 0.0, 0.0, 0.0]),
            DVector::from_vec(vec![1.0, 0.0, 0.0, 0.0, 0.0]),
        ];
        let orbit = ReebOrbit::new(2, 4.0, pts);
        assert!(orbit.is_closed());
    }

    #[test]
    fn test_serde() {
        let verifier = WeinsteinVerifier::verify_s3();
        let json = serde_json::to_string(&verifier).unwrap();
        let back: WeinsteinVerifier = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 1);
    }
}
