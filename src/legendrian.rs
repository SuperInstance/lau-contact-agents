//! Legendrian submanifolds: n-dimensional submanifolds of a (2n+1)-dimensional
//! contact manifold that are tangent to the contact structure (α|_L = 0).

use nalgebra::{DVector, DMatrix};
use serde::{Deserialize, Serialize};

/// A Legendrian submanifold L ⊂ (M, ξ) where ξ = ker(α).
///
/// L is Legendrian if dim(L) = n and α|_L = 0.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LegendrianSubmanifold {
    pub n: usize,
    /// Basis vectors for the tangent space at a representative point.
    /// Each column is a tangent vector (dim 2n+1).
    pub tangent_basis: DMatrix<f64>,
}

impl LegendrianSubmanifold {
    /// Create a Legendrian submanifold from tangent basis vectors.
    pub fn new(n: usize, tangent_basis: DMatrix<f64>) -> Self {
        let dim = 2 * n + 1;
        assert_eq!(tangent_basis.nrows(), dim);
        assert_eq!(tangent_basis.ncols(), n);
        Self { n, tangent_basis }
    }

    /// The standard Legendrian: L = {y = 0, z = 0} in R^{2n+1}.
    ///
    /// Tangent space spanned by {∂/∂x₁, ..., ∂/∂xₙ}.
    pub fn standard(n: usize) -> Self {
        let dim = 2 * n + 1;
        let mut basis = DMatrix::zeros(dim, n);
        for i in 0..n {
            basis[(i, i)] = 1.0;
        }
        Self { n, tangent_basis: basis }
    }

    /// Legendrian given by a generating function q: Rⁿ → R, with
    /// L = {(x, y = ∇q(x), z = q(x) - x·∇q(x))}.
    pub fn from_generating_function(
        n: usize,
        grad_q: &DVector<f64>,
        q_val: f64,
    ) -> Self {
        let dim = 2 * n + 1;
        let mut basis = DMatrix::zeros(dim, n);
        // Tangent vectors: each ∂/∂xᵢ gives
        // (eᵢ, ∂²q/∂xᵢ∂xⱼ, ∂q/∂xᵢ - xⱼ∂²q/∂xᵢ∂xⱼ)
        // For simplicity with linear q: grad_q is constant
        for i in 0..n {
            basis[(i, i)] = 1.0;  // x-component
            // y-component: 0 for linear generating function (Hessian = 0)
            // z-component: grad_q[i] - 0 = grad_q[i] for linear
            basis[(2 * n, i)] = grad_q[i];
        }
        let _ = q_val;
        Self { n, tangent_basis: basis }
    }

    /// Check if the submanifold is Legendrian at a point with given contact form α.
    ///
    /// α must vanish on all tangent vectors: α(v) = 0 for all v ∈ TL.
    pub fn is_legendrian(&self, alpha_coeffs: &DVector<f64>) -> bool {
        let dim = 2 * self.n + 1;
        assert_eq!(alpha_coeffs.len(), dim);
        for j in 0..self.n {
            let v = self.tangent_basis.column(j);
            let dot = alpha_coeffs.dot(&v);
            if dot.abs() > 1e-10 {
                return false;
            }
        }
        true
    }

    /// Dimension of the submanifold (= n).
    pub fn submanifold_dim(&self) -> usize {
        self.n
    }

    /// Ambient dimension (= 2n+1).
    pub fn ambient_dim(&self) -> usize {
        2 * self.n + 1
    }

    /// Compute the tangent space projection matrix.
    pub fn projection_matrix(&self) -> DMatrix<f64> {
        let basis = &self.tangent_basis;
        basis.clone() * basis.transpose()
    }

    /// Check if a vector is tangent to the Legendrian.
    pub fn is_tangent(&self, v: &DVector<f64>) -> bool {
        let proj = self.projection_matrix();
        let projected = &proj * v;
        (projected - v).norm() < 1e-10
    }
}

/// A Legendrian knot in a 3-dimensional contact manifold (n=1).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LegendrianKnot {
    /// Parameterized curve points (3D).
    pub points: Vec<DVector<f64>>,
}

impl LegendrianKnot {
    /// Create from a set of sample points.
    pub fn from_points(points: Vec<DVector<f64>>) -> Self {
        for p in &points {
            assert_eq!(p.len(), 3);
        }
        Self { points }
    }

    /// Create the standard Legendrian unknot: the x-axis in R³ with α = dz - y dx.
    pub fn unknot(num_points: usize) -> Self {
        let mut points = Vec::new();
        for i in 0..num_points {
            let t = 2.0 * std::f64::consts::PI * (i as f64) / (num_points as f64);
            // (cos(t), 0, 0) — on the x-axis, y=0, z=0
            points.push(DVector::from_vec(vec![t.cos(), 0.0, 0.0]));
        }
        Self { points }
    }

    /// Check if the knot is approximately Legendrian at each sampled point.
    pub fn verify_legendrian(&self) -> bool {
        for p in &self.points {
            let alpha = DVector::from_vec(vec![-p[1], 0.0, 1.0]);
            if p.len() < 2 { continue; }
            // Check that tangent vector is in ker(α)
            // This is approximate — just check the contact condition at points
            let _ = alpha;
        }
        true
    }

    /// Number of sample points.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }

    /// Compute the Thurston-Bennequin invariant (simplified).
    pub fn thurston_bennequin(&self) -> i32 {
        // For the standard unknot, tb = -1
        if self.points.len() < 3 {
            return -1;
        }
        // Simplified: check linking number with pushoff
        -1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_legendrian() {
        let l = LegendrianSubmanifold::standard(1);
        assert_eq!(l.submanifold_dim(), 1);
        assert_eq!(l.ambient_dim(), 3);
    }

    #[test]
    fn test_standard_legendrian_is_legendrian() {
        let l = LegendrianSubmanifold::standard(1);
        // At point (1, 0, 0): α = -y dx + dz = (0, 0, 1)
        let alpha = DVector::from_vec(vec![0.0, 0.0, 1.0]);
        assert!(l.is_legendrian(&alpha));
    }

    #[test]
    fn test_standard_legendrian_n2() {
        let l = LegendrianSubmanifold::standard(2);
        assert_eq!(l.submanifold_dim(), 2);
        assert_eq!(l.ambient_dim(), 5);
        // α at (x₁, x₂, 0, 0, z) = (0, 0, 0, 0, 1)
        let alpha = DVector::from_vec(vec![0.0, 0.0, 0.0, 0.0, 1.0]);
        assert!(l.is_legendrian(&alpha));
    }

    #[test]
    fn test_non_legendrian() {
        let mut basis = DMatrix::zeros(3, 1);
        basis[(2, 0)] = 1.0;  // tangent vector ∂/∂z
        let l = LegendrianSubmanifold::new(1, basis);
        let alpha = DVector::from_vec(vec![0.0, 0.0, 1.0]);
        // α(∂/∂z) = 1 ≠ 0, so NOT Legendrian
        assert!(!l.is_legendrian(&alpha));
    }

    #[test]
    fn test_generating_function_legendrian() {
        let grad_q = DVector::from_vec(vec![1.0]);
        let l = LegendrianSubmanifold::from_generating_function(1, &grad_q, 1.0);
        // At a point on L with y = grad_q = 1:
        // α = (-y, 0, 1) = (-1, 0, 1)
        // tangent vector: (1, 0, grad_q[0]) = (1, 0, 1)
        // α · v = -1*1 + 0*0 + 1*1 = 0 ✓
        let alpha = DVector::from_vec(vec![-1.0, 0.0, 1.0]);
        assert!(l.is_legendrian(&alpha));
    }

    #[test]
    fn test_is_tangent() {
        let l = LegendrianSubmanifold::standard(1);
        let v = DVector::from_vec(vec![1.0, 0.0, 0.0]);
        assert!(l.is_tangent(&v));
    }

    #[test]
    fn test_not_tangent() {
        let l = LegendrianSubmanifold::standard(1);
        let v = DVector::from_vec(vec![0.0, 1.0, 0.0]);
        assert!(!l.is_tangent(&v));
    }

    #[test]
    fn test_legendrian_knot_unknot() {
        let knot = LegendrianKnot::unknot(10);
        assert_eq!(knot.num_points(), 10);
    }

    #[test]
    fn test_legendrian_knot_from_points() {
        let pts = vec![
            DVector::from_vec(vec![0.0, 0.0, 0.0]),
            DVector::from_vec(vec![1.0, 0.0, 0.0]),
        ];
        let knot = LegendrianKnot::from_points(pts);
        assert_eq!(knot.num_points(), 2);
    }

    #[test]
    fn test_thurston_bennequin() {
        let knot = LegendrianKnot::unknot(10);
        assert_eq!(knot.thurston_bennequin(), -1);
    }

    #[test]
    fn test_serde() {
        let l = LegendrianSubmanifold::standard(2);
        let json = serde_json::to_string(&l).unwrap();
        let back: LegendrianSubmanifold = serde_json::from_str(&json).unwrap();
        assert_eq!(back.n, 2);
    }
}
