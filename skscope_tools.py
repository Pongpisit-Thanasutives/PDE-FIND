import jax.numpy as jnp
from skscope import ScopeSolver

def best_subset_solution(X_pre, y_pre, sparsity, p=1):
    objective_function = lambda coefs: jnp.linalg.norm(y_pre-X_pre@coefs)**p
    scope_solver = ScopeSolver(dimensionality=X_pre.shape[-1], sparsity=sparsity)
    return scope_solver.solve(objective_function)
