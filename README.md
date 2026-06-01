# lau-contact-agents

> Contact geometry for agents — contact manifolds, Reeb fields, contactomorphisms, Legendrian submanifolds, and thermodynamic agent dynamics.

A Rust crate implementing the structures and theorems of **contact geometry**: contact forms and their non-degeneracy condition, Reeb vector fields, contact Hamiltonian dynamics (inherently dissipative), Darboux charts, Gray stability, the Weinstein conjecture on closed Reeb orbits, contact homology, Legendrian submanifolds, and a thermodynamic agent model where extensive/intensive variables live on a contact manifold.

Every result is verified by **102 unit tests**.

---

## What This Does

This crate provides concrete, finite-dimensional implementations of contact-geometric objects on ℝ^(2n+1):

- **Contact forms** α with the non-degeneracy condition α ∧ (dα)ⁿ ≠ 0
- **Reeb vector fields** R satisfying α(R) = 1, dα(R, ·) = 0
- **Contact Hamiltonians** and their dissipative dynamics (unlike symplectic H, contact H is not conserved)
- **Darboux charts** — local coordinates where α takes the standard form dz − Σ yᵢ dxᵢ
- **Contactomorphisms** — diffeomorphisms preserving the contact structure (φ*α = f·α)
- **Legendrian submanifolds** — maximal submanifolds tangent to the contact distribution
- **Gray stability** — deformation of contact structures through isotopies
- **Weinstein conjecture** — closed Reeb orbits on compact contact manifolds
- **Contact homology** — Floer-type invariant built from Reeb orbits
- **Thermodynamic agents** — physical systems on contact manifolds with entropy production

All structures serialize via `serde`.

---

## Key Idea

**Contact geometry is the odd-dimensional cousin of symplectic geometry.** Where symplectic manifolds model conservative Hamiltonian mechanics on even-dimensional phase spaces, contact manifolds model *dissipative* systems on odd-dimensional spaces (2n+1).

The extra dimension is tracked by the variable z, and the contact form α = dz − Σ yᵢ dxᵢ forces z to evolve in a way that captures **entropy production**. Unlike symplectic systems where energy is conserved, contact Hamiltonians generically dissipate: dH/dt = −R(H)·H ≠ 0.

This crate models **thermodynamic agents** whose states live on contact manifolds:
- **Position variables x** → extensive quantities (volume, particle number)
- **Momentum variables y** → intensive quantities (pressure, chemical potential)
- **The z coordinate** → entropy production (cumulative dissipation)

The contact structure enforces the laws of thermodynamics geometrically.

---

## Install

```toml
[dependencies]
lau-contact-agents = "0.1"
```

Or clone directly:

```bash
git clone https://github.com/SuperInstance/lau-contact-agents.git
cargo build
```

### Dependencies

| Crate | Purpose |
|-------|---------|
| `nalgebra` 0.33 | Linear algebra (vectors, matrices, determinants, SVD) |
| `serde` 1 | Serialization of all structures |
| `serde_json` 1 | JSON serialization for tests and storage |

---

## Quick Start

### Contact form and Reeb field

```rust
use lau_contact_agents::{ContactForm, ReebVectorField};
use nalgebra::DVector;

let form = ContactForm::standard(1); // α = dz − y dx on ℝ³
let p = DVector::from_vec(vec![1.0, 2.0, 3.0]);

assert!(form.is_contact_at(&p)); // α ∧ (dα)¹ ≠ 0

let reeb = ReebVectorField::new(1);
assert!(reeb.verify_alpha_condition(&p, &form)); // α(R) = 1
assert!(reeb.verify_dalpha_condition(&p, &form)); // dα(R, ·) = 0
```

### Contact Hamiltonian dynamics

```rust
use lau_contact_agents::{ContactHamiltonian, DissipativeContactSystem};
use nalgebra::DVector;

// Contact Hamiltonian with ∇H at a point
let grad = DVector::from_vec(vec![1.0, 2.0, 3.0]);
let h = ContactHamiltonian::new(1, 5.0, grad);
let xh = h.vector_field();
// ẋ = ∂H/∂y = 2, ẏ = −∂H/∂x = −1

assert!(!h.is_conserved()); // Contact H is dissipative!
let rate = h.dissipation_rate(3.0); // dH/dt = −R(H)·H = −15

// Dissipative system with friction coefficients
let gamma = DVector::from_vec(vec![1.0]);
let sys = DissipativeContactSystem::new(1, gamma);
let state = DVector::from_vec(vec![0.0, 1.0, 0.0]);
let new_state = sys.step(&state, 0.1);
```

### Legendrian submanifold

```rust
use lau_contact_agents::LegendrianSubmanifold;

// Standard Legendrian: {y = 0, z = 0} with tangent space spanned by ∂/∂x
let l = LegendrianSubmanifold::standard(1);
let alpha = DVector::from_vec(vec![0.0, 0.0, 1.0]); // α at y=0
assert!(l.is_legendrian(&alpha)); // α|_L = 0

// From a generating function q(x)
let grad_q = DVector::from_vec(vec![1.0]);
let l2 = LegendrianSubmanifold::from_generating_function(1, &grad_q, 1.0);
```

### Thermodynamic agent

```rust
use lau_contact_agents::ThermoAgent;
use nalgebra::DVector;

let mut agent = ThermoAgent::new(1, DVector::from_vec(vec![2.0, 0.0, 0.0]));
let energy = agent.internal_energy(); // E = c₁·x₁² = 4.0
let free = agent.free_energy(1.0);    // F = E − T·S = 4.0

let gamma = DVector::from_vec(vec![1.0]);
agent.evolve(0.1, &gamma); // Step forward in contact Hamiltonian flow
```

---

## API Reference

### `contact_form` — Contact Forms and Reeb Vector Fields

| Type / Function | Description |
|----------------|-------------|
| `ContactForm` | A 1-form α on ℝ^(2n+1). Methods: `standard(n)`, `from_coefficients(n, c)`, `evaluate_at(p)`, `exterior_derivative(p)`, `volume_form(p)`, `is_contact_at(p)`, `dimension` |
| `ReebVectorField` | The Reeb vector R with α(R)=1, dα(R,·)=0. Methods: `new(n)`, `evaluate(p)`, `verify_alpha_condition(p, form)`, `verify_dalpha_condition(p, form)` |

The standard contact form is α = dz − Σ yᵢ dxᵢ. Its exterior derivative dα = Σ dxᵢ ∧ dyᵢ is the standard symplectic form. The volume form α ∧ (dα)ⁿ is computed as a determinant and must be nonzero everywhere.

### `contact_hamiltonian` — Dissipative Hamiltonian Dynamics

| Type / Function | Description |
|----------------|-------------|
| `ContactHamiltonian` | H on a contact manifold with gradient. Methods: `new(n, h, grad_h)`, `vector_field`, `vector_field_full(reeb_h)`, `is_conserved`, `dissipation_rate(reeb_h)`, `dimension` |
| `DissipativeContactSystem` | Friction model with coefficients γᵢ > 0. Methods: `new(n, gamma)`, `hamiltonian_value(y)`, `equations_of_motion(state)`, `step(state, dt)`, `total_gamma`, `is_equilibrium(state)` |

Contact Hamiltonian equations: ẋᵢ = ∂H/∂yᵢ, ẏᵢ = −∂H/∂xᵢ, ż = −H + Σ yᵢ ∂H/∂yᵢ. Energy dissipation: dH/dt = −R(H)·H.

### `darboux` — Darboux's Theorem

| Type / Function | Description |
|----------------|-------------|
| `DarbouxChart` | Local coordinates where α = dz − Σ yᵢ dxᵢ. Methods: `new(n, center)`, `with_transform(n, center, T)`, `to_darboux(point)`, `from_darboux(darboux_point)`, `verify_standard_form(alpha_coeffs)`, `dimension` |
| `construct_darboux_chart(n, point, α, dα)` | Build a chart by finding a symplectic basis for dα restricted to ker(α) |

### `legendrian` — Legendrian Submanifolds and Knots

| Type / Function | Description |
|----------------|-------------|
| `LegendrianSubmanifold` | n-dimensional submanifold with α\|_L = 0. Methods: `new(n, tangent_basis)`, `standard(n)`, `from_generating_function(n, grad_q, q_val)`, `is_legendrian(alpha)`, `submanifold_dim`, `ambient_dim`, `projection_matrix`, `is_tangent(v)` |
| `LegendrianKnot` | Legendrian curve in a 3D contact manifold. Methods: `from_points(pts)`, `unknot(n_pts)`, `verify_legendrian`, `num_points`, `thurston_bennequin` |

The standard Legendrian is L = {y = 0, z = 0} with tangent space spanned by {∂/∂x₁, …, ∂/∂xₙ}. Legendrian knots are 1D Legendrian submanifolds of 3D contact manifolds, classified by the Thurston-Bennequin invariant.

### `contactomorphism` — Structure-Preserving Diffeomorphisms

| Type / Function | Description |
|----------------|-------------|
| `Contactomorphism` | φ: M → M with φ*α = f·α. Methods: `identity(n)`, `from_jacobian(n, J, strict)`, `compose`, `inverse`, `preserves_contact_structure(α)`, `dimension`, `z_translation(n, dz)`, `xy_rotation(n, i, θ)` |

A *strict* contactomorphism has φ*α = α (f = 1). Rotations in the (xᵢ, yᵢ) planes preserve the symplectic form dα and hence the contact structure.

### `gray_stability` — Gray's Stability Theorem

| Type / Function | Description |
|----------------|-------------|
| `GrayStability` | Checks isotopy of contact structures. Methods: `check_path(n, α₀, α₁, dα₀, dα₁, steps)`, `moser_vector_field(αₜ, α̇ₜ, R)`, `verify(n)` |

Given two contact forms α₀, α₁, checks that the linear interpolation αₜ = (1−t)α₀ + tα₁ remains contact for all t ∈ [0,1]. The Moser vector field generates the isotopy.

### `weinstein` — Weinstein Conjecture

| Type / Function | Description |
|----------------|-------------|
| `ReebOrbit` | A closed periodic orbit of the Reeb field. Methods: `new(n, period, points)`, `is_closed`, `action`, `conley_zehnder_index`, `num_points` |
| `WeinsteinVerifier` | Verifies the conjecture on specific manifolds. Methods: `new(n)`, `standard_reeb_orbits(n, num_pts)`, `weinstein_satisfied`, `num_orbits`, `verify_s3`, `verify_s2n1(n)` |

The Weinstein conjecture (proved by Taubes in dimension 3): every Reeb field on a closed contact manifold has at least one closed orbit.

### `contact_homology` — Floer-Type Invariant

| Type / Function | Description |
|----------------|-------------|
| `HomologyGenerator` | A Reeb orbit with grading and action. Methods: `new(id, grading, action)`, `degree`, `parity` |
| `Differential` | Counts pseudoholomorphic curves between orbits (mod 2). Methods: `zero`, `add_edge(src, tgt)`, `apply(generators)`, `check_d_squared_zero` |
| `ContactHomologyComplex` | Chain complex (generators + differential). Methods: `new(gens, diff)`, `compute_homology`, `betti_numbers`, `euler_characteristic` |
| `ContactHomology` | Result of homology computation. Methods: `zero`, `rank`, `is_trivial` |

Contact homology is generated by closed Reeb orbits graded by the Conley-Zehnder index, with the differential counting rigid J-holomorphic curves in the symplectization. It satisfies ∂² = 0.

### `thermo_dynamics` — Thermodynamic Agent Dynamics

| Type / Function | Description |
|----------------|-------------|
| `ThermoAgent` | Agent on a contact manifold with state (x, y, z). Methods: `new(n, state)`, `with_energy(n, state, coeffs)`, `extensive`, `intensive`, `entropy`, `internal_energy`, `free_energy(T)`, `entropy_production(γ)`, `evolve(dt, γ)`, `is_equilibrium`, `thermodynamic_length(other)`, `dimension` |
| `ThermoMultiAgent` | System of interacting agents. Methods: `new(agents, coupling)`, `decoupled(agents)`, `num_agents`, `total_entropy_production(γ)`, `evolve(dt, γ)`, `all_equilibrium`, `total_energy` |

Internal energy E = Σ cᵢ xᵢ². Free energy F = E − T·S. Entropy production rate = Σ γᵢ yᵢ² ≥ 0. The contact structure ensures the second law of thermodynamics.

---

## How It Works

### Architecture

```
contact_form ──→ contact_hamiltonian ──→ thermo_dynamics
     │                    │
     ├── darboux          ├── legendrian
     ├── contactomorphism ├── gray_stability
     └── reeb             ├── weinstein (Reeb orbits)
                          └── contact_homology
```

### Coordinates

The crate works in local coordinates (x₁,…,xₙ, y₁,…,yₙ, z) on ℝ^(2n+1), where:

- x ∈ ℝⁿ: position/extensive variables
- y ∈ ℝⁿ: momentum/intensive variables
- z ∈ ℝ: the "extra" coordinate (entropy production / thermodynamic displacement)

The standard contact form α = dz − Σ yᵢ dxᵢ defines the contact structure ξ = ker(α).

### Non-Degeneracy Check

The defining condition α ∧ (dα)ⁿ ≠ 0 is checked by computing the determinant of a (2n+1)×(2n+1) matrix whose first row is the coefficients of α and remaining rows encode dα. This is the algebraic manifestation of the contact condition.

### Dissipation via Contact Hamiltonians

Unlike symplectic Hamiltonians (conserved by the flow), contact Hamiltonians satisfy:

```
dH/dt = −R(H)·H
```

This is always nonzero (generically), making contact geometry the natural setting for **irreversible** processes. The z-coordinate accumulates this dissipation as entropy production.

---

## The Math

### Contact Manifolds

A **contact manifold** (M, ξ) is a (2n+1)-dimensional manifold with a maximally non-integrable hyperplane distribution ξ. Locally, ξ = ker(α) for a **contact form** α satisfying:

```
α ∧ (dα)ⁿ ≠ 0   (everywhere)
```

This is an odd-dimensional analog of a symplectic structure. The standard example: α = dz − y dx on ℝ³.

### Reeb Vector Field

The **Reeb vector field** R is uniquely determined by:

```
α(R) = 1    and    dα(R, ·) = 0
```

For the standard contact form, R = ∂/∂z. The Reeb flow is a special geodesic-like flow on the contact manifold.

### Darboux's Theorem

Every contact form is locally equivalent to the standard form: there exist local coordinates where α = dz − Σ yᵢ dxᵢ. This is the contact analog of Darboux's theorem in symplectic geometry.

Contact geometry has **no local invariants** — all contact structures look the same locally. Global invariants (Gray stability, contact homology) distinguish different contact structures.

### Contact Hamiltonian Dynamics

Given a contact form α and a function H (the contact Hamiltonian), the **contact Hamiltonian vector field** X_H satisfies:

```
dα(X_H, ·) = dH − R(H)α
α(X_H) = −H
```

In coordinates: ẋᵢ = ∂H/∂yᵢ, ẏᵢ = −∂H/∂xᵢ + yᵢ ∂H/∂z, ż = −H + Σ yᵢ ∂H/∂yᵢ.

Key difference from symplectic dynamics: **H is not conserved**. The dissipation rate dH/dt = −R(H)·H makes this ideal for modeling thermodynamic systems.

### Legendrian Submanifolds

A **Legendrian submanifold** L ⊂ (M, ξ) has maximal dimension (n) among submanifolds tangent to ξ: dim(L) = n and α|_L = 0.

Legendrian submanifolds can be constructed from **generating functions** q: ℝⁿ → ℝ via:

```
L = {(x, ∇q(x), q(x) − x·∇q(x))}
```

In 3D contact manifolds, Legendrian knots are classified by the **Thurston-Bennequin invariant** tb(K) and the **rotation number** rot(K).

### Gray Stability

Gray's theorem: if {αₜ}ₜ∈[0,1] is a smooth family of contact forms on a closed manifold, then there exists an isotopy {φₜ} with φₜ*α₀ = fₜ·α₁ₜ.

This means contact structures cannot be destroyed by smooth deformation — they are **stable**. The proof uses the **Moser trick**: find a vector field Xₜ whose flow implements the isotopy.

### Weinstein Conjecture

On any closed contact manifold (M, α), the Reeb vector field has at least one closed orbit.

Proved by Taubes (2007) for all 3-dimensional contact manifolds, using Seiberg-Witten theory. For higher dimensions, it remains open in general but is known for many cases (standard contact sphere, toric contact manifolds, etc.).

### Contact Homology

A **Floer homology** for contact structures, generated by closed Reeb orbits and graded by the Conley-Zehnder index. The differential counts rigid J-holomorphic curves in the symplectization ℝ × M.

The resulting homology HC(M, ξ) is an invariant of the contact structure. It satisfies ∂² = 0 (gluing theorem for holomorphic curves) and its Euler characteristic is related to classical invariants.

### Thermodynamic Contact Geometry

The Gibbs contact form on the thermodynamic phase space is:

```
α = dS − (1/T)dE − (p/T)dV − (μ/T)dN
```

where S = entropy, T = temperature, E = energy, V = volume, p = pressure, μ = chemical potential, N = particle count. The contact structure encodes the equation of state and the first and second laws of thermodynamics.

---

## License

MIT
