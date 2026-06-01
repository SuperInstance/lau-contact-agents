# lau-contact-agents

> Contact geometry for agents: Reeb flows, Legendrian submanifolds, contact Hamiltonians, and thermodynamic dynamics

## What This Does

This crate implements contact geometry — the odd-dimensional cousin of symplectic geometry — applied to agent dynamics. It provides contact forms, Reeb vector fields, contactomorphisms, Legendrian submanifolds, contact Hamiltonian dynamics, Darboux's theorem, Gray stability, the Weinstein conjecture, contact homology, and a thermodynamic agent model where extensive/intensive quantities live on a contact manifold with entropy production.

## The Key Idea

Symplectic geometry is for conservative systems (energy preserved). Contact geometry is for dissipative systems — energy can flow in and out. A contact manifold has an odd dimension 2n+1, and the extra dimension tracks a "thermodynamic displacement" (like entropy). The contact form α ∧ (dα)ⁿ ≠ 0 ensures the geometry never degenerates. For agents, this means: position = extensive variables (energy, volume), momentum = intensive variables (temperature, pressure), the extra coordinate z = entropy production. Contact Hamiltonian dynamics naturally model dissipation.

## Install

```toml
[dependencies]
lau-contact-agents = { git = "https://github.com/SuperInstance/lau-contact-agents" }
```

## Quick Start

```rust
use lau_contact_agents::*;
use nalgebra::DVector;

// Create standard contact form on R^3 (n=1, so dim=3)
let alpha = ContactForm::standard(1);

// Evaluate at a point p = (x, y, z)
let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);
let alpha_p = alpha.evaluate_at(&p);
println!("α at p: {:?}", alpha_p); // [−y, 0, 1] = [−2, 0, 1]

// Compute dα
let dalpha = alpha.exterior_derivative(&p);
println!("dα:\n{:?}", dalpha); // [[0, 1, 0], [−1, 0, 0], [0, 0, 0]]

// Reeb vector field (unique vector R such that α(R)=1, dα(R,·)=0)
let reeb = alpha.reeb_vector_field(&p);
println!("Reeb vector: {:?}", reeb); // [0, 0, 1] for standard form

// Contact Hamiltonian dynamics
let hamiltonian = ContactHamiltonian::new(
    1,
    5.0,                                        // H(p) = 5.0
    DVector::from_vec(vec![1.0, 2.0, 0.5]),     // ∇H
);
let xh = hamiltonian.vector_field();
println!("Contact Hamiltonian VF: {:?}", xh);

// Legendrian submanifold (n-dim, tangent to contact structure)
let legendre = LegendrianSubmanifold::standard(1);
assert!(legendre.is_legendrian(&alpha, &p));

// Darboux chart (local normal form)
let chart = DarbouxChart::new(1, DVector::from_vec(vec![0.0, 0.0, 0.0]));
let darboux = chart.to_darboux(&p);

// Gray stability: check if two contact forms are isotopic
let stability = GrayStability::check_path(
    1, &alpha0_coeffs, &alpha1_coeffs,
    &dalpha0, &dalpha1, 100,
);
println!("Isotopic: {}", stability.is_isotopic);

// Thermodynamic agent
let mut agent = ThermoAgent::new(1, DVector::from_vec(vec![1.0, 0.0, 0.0]));
agent.evolve(0.1); // Evolve along contact Hamiltonian flow
println!("Agent state: {:?}", agent.state);
```

## API Reference

### `contact_form`

| Type | Description |
|------|-------------|
| `ContactForm::standard(n)` | Standard form α = dz − Σ yᵢdxᵢ on R^{2n+1}. |
| `evaluate_at(p)` | Covector α_p at point p. |
| `exterior_derivative(p)` | dα as skew-symmetric matrix. |
| `reeb_vector_field(p)` | Reeb vector R: α(R)=1, dα(R,·)=0. |
| `verify_contact_condition(p)` | Check α ∧ (dα)ⁿ ≠ 0. |

### `contactomorphism`

| Type | Description |
|------|-------------|
| `Contactomorphism::identity(n)` | Identity map. |
| `from_jacobian(n, J, strict)` | From Jacobian matrix. |
| `compose(other)` | Compose two contactomorphisms. |
| `inverse()` | Inverse map. |
| `preserves_contact(alpha, p)` | Check φ*α = f·α. |

### `legendrian`

| Type | Description |
|------|-------------|
| `LegendrianSubmanifold::standard(n)` | {y=0, z=0} in R^{2n+1}. |
| `from_generating_function(n, grad_q, q_val)` | Legendrian from generating function q: Rⁿ→R. |
| `is_legendrian(alpha, p)` | Verify α|_L = 0. |

### `contact_hamiltonian`

| Type | Description |
|------|-------------|
| `ContactHamiltonian::new(n, h, grad_h)` | Contact Hamiltonian at a point. |
| `vector_field()` | X_H: dα(X_H,·) = (dH−R(H)α)(·), α(X_H)=−H. |
| `vector_field_full(reeb_h)` | Full vector field with Reeb coupling. |

### `darboux`

| Type | Description |
|------|-------------|
| `DarbouxChart::new(n, center)` | Local chart where α takes standard form. |
| `to_darboux(point)` | Transform to Darboux coordinates. |
| `from_darboux(point)` | Transform back. |

### `gray_stability`

| Function | Description |
|----------|-------------|
| `GrayStability::check_path(n, α₀, α₁, dα₀, dα₁, steps)` | Check isotopy along linear path. |

### `weinstein`

Weinstein conjecture: every contact form on a closed manifold has a closed Reeb orbit.

### `contact_homology`

| Type | Description |
|------|-------------|
| `HomologyGenerator::new(id, grading, action)` | Closed Reeb orbit as chain generator. |
| `Differential::add_edge(source, target)` | Holomorphic curve between orbits. |
| `compute_homology(generators, diff)` | Contact homology groups. |

### `thermo_dynamics`

| Type | Description |
|------|-------------|
| `ThermoAgent::new(n, state)` | Agent on contact manifold. x=extensive, y=intensive, z=entropy. |
| `evolve(dt)` | Flow along contact Hamiltonian vector field. |
| `entropy_production()` | Rate of entropy change ż. |
| `energy()` | Internal energy at current state. |

## How It Works

1. **Contact Form**: A 1-form α on M^{2n+1} with α ∧ (dα)ⁿ ≠ 0. This "maximal non-integrability" means there's no hypersurface tangent to ker(α).
2. **Reeb Flow**: The unique vector field R preserving α. Its closed orbits encode the contact topology.
3. **Legendrian Submanifolds**: n-dimensional submanifolds L with α|_L = 0 (maximally tangent to ξ = ker(α)).
4. **Contact Hamiltonian**: Generates dissipative dynamics. Unlike symplectic Hamiltonians, energy is not conserved.
5. **Thermodynamic Interpretation**: (x, y, z) = (extensive, intensive, entropy). Contact flow satisfies the first and second laws simultaneously.

## The Math

- **Contact Form**: α on M^{2n+1} with α ∧ (dα)ⁿ ≠ 0.
- **Reeb Vector Field**: ι_R dα = 0, ι_R α = 1.
- **Darboux Theorem**: Locally, every contact form is α₀ = dz − Σ yᵢ dxᵢ.
- **Gray Stability**: If αₜ = (1−t)α₀ + tα₁ remains contact for t ∈ [0,1], the structures are isotopic.
- **Contact Hamiltonian**: X_H satisfies dα(X_H, ·) = dH − R(H)α, α(X_H) = −H.

## Testing

102 tests covering:
- Contact form evaluation and exterior derivative
- Reeb vector field computation
- Contact condition verification
- Legendrian submanifold verification
- Contact Hamiltonian vector field
- Darboux chart normal form
- Gray stability path checking
- Contact homology differential ∂² = 0
- Thermodynamic agent evolution
- Entropy production positivity

## License

MIT
